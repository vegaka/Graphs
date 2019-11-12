#include <algorithm/MaxFlow.h>
#include <queue>
#include <iostream>
#include <thread>
#include <cmath>
#include <atomic>

#define DEFAULT_NUM_THREADS 4

using data_vec = std::vector<long>;

static long BFSColoring(Graph &g, data_vec &heights, data_vec &residuals, data_vec &wave, data_vec &reverse,
                        data_vec &color, long startVertex, long startLevel, long currentWave) {
    long coloredVertices = 0;
    long currentLevel = startLevel;
    std::queue<long> queue;
    queue.push(startVertex);

    while(!queue.empty()) {
        long curNode = queue.front();
        queue.pop();

        currentLevel++;
        data_vec neighbours = g.getNeighbourListFor(curNode);
        long edgeId;
        for (auto &n : neighbours) {
            edgeId = g.getIdFromSrcDst(curNode, n);
            edgeId = edgeId == -1 ? reverse[g.getIdFromSrcDst(n, curNode)] : edgeId;

            if (residuals[edgeId] <= 0) continue;

            if (color[n] == 0) {
                color[n] = 1;
                coloredVertices++;

                if (heights[n] < currentLevel) {
                    heights[n] = currentLevel;
                }
                wave[n] = currentWave;
                queue.push(n);
            }
        }
    }

    return coloredVertices;
}

static void globalRelabel(Graph &g, data_vec &heights, data_vec &residuals, data_vec &wave, data_vec &reverse,
                          const long s, const long t, long &currentWave) {
    currentWave++;
    data_vec color(g.getNE(), 0);
    std::queue<long> queue;
    long coloredVertices = BFSColoring(g, heights, residuals, wave, reverse, color, t, 0, currentWave);
    if (coloredVertices < g.getNV()) {
        BFSColoring(g, heights, residuals, wave, reverse, color, s, g.getNV(), currentWave);
    }
}

static void processWithRelabel(long node, Graph &g, std::queue<long> &queue, data_vec &heights,
                               std::vector<std::atomic_long> &excess, data_vec &capacities, data_vec &reverse,
                               std::vector<std::atomic_long> &residuals, data_vec &wave, std::vector<bool> &isActive,
                               const long s, const long t) {
    //std::cout << "Processing node: " << node << std::endl;
    while (excess[node] > 0) {
        long e = excess[node];
        long h = std::numeric_limits<long>::max();
        long nextV = -1;
        long edgeId;
        data_vec neighbours = g.getNeighbourListFor(node);

        if (neighbours.empty()) return;

        for (auto &n: neighbours) {
            edgeId = g.getIdFromSrcDst(node, n);
            edgeId = edgeId == -1 ? reverse[g.getIdFromSrcDst(n, node)] : edgeId;

            if (residuals[edgeId].load() <= 0) continue;

            if (heights[n] < h) {
                nextV = n;
                h = heights[n];
            }
        }

        edgeId = g.getIdFromSrcDst(node, nextV);
        edgeId = edgeId == -1 ? reverse[g.getIdFromSrcDst(nextV, node)] : edgeId;

        if (heights[node] > h) {
            long delta = std::min(e, residuals[edgeId].load());
            if (wave[node] <= wave[nextV] && heights[node] > heights[nextV]) {
                residuals[edgeId].fetch_sub(delta);
                residuals[reverse[edgeId]].fetch_add(delta);
                excess[node].fetch_sub(delta);
                excess[nextV].fetch_add(delta);

                if (nextV != t && nextV != s && !isActive[nextV]) {
                    queue.push(nextV);
                    isActive[nextV] = true;
                }
            }
        } else if (heights[node] < h + 1){
            heights[node] = h + 1;
        }
    }
    queue.pop();
    isActive[node] = false;
}

static void execute(int threadId, Graph &g, std::queue<long> &queue, data_vec &heights,
                    std::vector<std::atomic_long> &excess, data_vec &capacities, data_vec &reverse,
                    std::vector<std::atomic_long> &residuals, data_vec &wave, std::vector<bool> &isActive,
                    const long s, const long t) {
    std::cout << "Thread id: " << threadId << std::endl;
    if (queue.empty()) return;
    while (excess[s] + excess[t] < 0 && !queue.empty()) {
        long curNode = queue.front();

        processWithRelabel(curNode, g, queue, heights, excess, capacities, reverse, residuals, wave, isActive, s, t);
    }
}

static void fillReverseAndCapacityVectors(Graph &g, std::vector<long>& reverse, std::vector<long>& capacities) {
    long newEdges = 0;
    long numEdges = 0;

    for (int i = 0; i < g.getNV(); ++i) {
        std::vector<long> neighbours = g.getNeighbourListFor(i);

        for (auto &n : neighbours) {
            long edgeId = g.getIdFromSrcDst(i, n);
            long revEdgeId = g.getIdFromSrcDst(n, i);

            // Don't process and edge if the reverse has already been processed
            if (reverse[std::max(edgeId, revEdgeId)] != -1) continue;

            if (edgeId == -1) {
                // This edge was created when undirecting the graph
                // The reverse edge has to exist so revEdgeId != -1
                reverse[revEdgeId] = g.getNE() + newEdges;
                reverse[g.getNE() + newEdges] = revEdgeId;

                capacities[g.getNE() + newEdges] = g.getWeightFromId(revEdgeId);

                newEdges++;
            } else if (revEdgeId == -1) {
                // The reverse edge was created when undirecting the graph
                reverse[edgeId] = g.getNE() + newEdges;
                reverse[g.getNE() + newEdges] = edgeId;

                capacities[g.getNE() + newEdges] = g.getWeightFromId(edgeId);

                newEdges++;
            } else {
                // Both edges exists in the original graph
                reverse[edgeId] = revEdgeId;
                reverse[revEdgeId] = edgeId;

                capacities[edgeId] += g.getWeightFromId(revEdgeId);
                capacities[revEdgeId] += g.getWeightFromId(edgeId);
            }

            numEdges++;
            if (numEdges % 50000 == 0) {
                std::cout << numEdges << std::endl;
            }
        }
    }
}

// A parallel lock free max flow algorithm
std::pair<data_vec, data_vec> PLFFlow(Graph &g, const long s, const long t) {
    std::vector<bool> isActive(g.getNV(), false);

    unsigned int numThreads = std::thread::hardware_concurrency();
    if (numThreads == 0) numThreads = DEFAULT_NUM_THREADS;
    std::vector<std::thread> threads(numThreads);
    std::vector<std::queue<long>> queues(numThreads);
    std::vector<std::atomic_bool> exchange(numThreads);
    for (auto &e : exchange) {
        e.store(false);
    }

    std::cout << "Creating neighbour lists" << std::endl;
    g.createNeighbourList(true);

    data_vec capacities = g.getWeights();
    capacities.resize(g.getUndirectedNumEdges(), 0);
    data_vec reverse(g.getUndirectedNumEdges(), -1);
    std::cout << "Filling reverse and capacity vectors" << std::endl;
    fillReverseAndCapacityVectors(g, reverse, capacities);

    std::vector<std::atomic_long> residuals(capacities.begin(), capacities.end());

    data_vec heights(g.getNV(), 0);
    std::vector<std::atomic_long> excess(g.getNV());
    for (auto & e : excess) {
        e.store(0);
    }

    data_vec wave(g.getNV(), 0);
    long currentWave = 0;

    heights[s] = g.getNV();

    std::cout << "Creating initial preflow" << std::endl;
    // Initial preflow
    std::vector<long> neighbours = g.getNeighbours(s);
    for (int i = 0; i < neighbours.size(); ++i) {
        long n = neighbours[i];
        long edgeId = g.getIdFromSrcDst(s, n);
        residuals[edgeId].store(0);
        residuals[reverse[edgeId]].store(capacities[edgeId]);
        excess[n].store(g.getWeightFromId(edgeId));
        excess[s].store(excess[s] - g.getWeightFromId(edgeId));

        queues[i % numThreads].push(n);
        isActive[n] = true;
    }

    const int itersBetweenGlobalRelabel = std::floor(g.getNV() / 2);
    int itersSinceGlobalRelabel = 0;
    unsigned long totalIters = 0;

    for (int i = 0; i < threads.size(); ++i) {
        threads[i] = std::thread(execute, i, std::ref(g), std::ref(queues[i]), std::ref(heights), std::ref(excess),
                                 std::ref(capacities), std::ref(reverse), std::ref(residuals), std::ref(wave), std::ref(isActive), s, t);
    }

    for (auto &thread : threads) {
        thread.join();
    }

    /*while (excess[s] + excess[t] < 0) {
        long curnode = activeNodes.front();

        //std::cout << "Processing node: " << curnode << ", numActive: " << activeNodes.size() << std::endl;

        //processWithRelabel(curnode, g, heights, excess, capacities, reverse, residuals, wave, activeNodes, isActive, s, t);

        itersSinceGlobalRelabel++;
        if (itersSinceGlobalRelabel >= itersBetweenGlobalRelabel) {
            //globalRelabel(g, heights, residuals, wave, reverse, s, t, currentWave);
            //std::cout << "Finished global relabel" << std::endl;
            itersSinceGlobalRelabel = 0;
        }

        totalIters++;
        if (totalIters % 10000000 == 0) {
            std::cout << "Total iterations: " << totalIters << std::endl;
            break;
        }
    }*/

    return std::make_pair(capacities, heights);
}