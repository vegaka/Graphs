#include <algorithm/MaxFlow.h>
#include <queue>
#include <iostream>
#include <thread>
#include <cmath>
#include <atomic>

#define DEFAULT_NUM_THREADS 4

using data_vec = std::vector<long>;

// Simple BFS from the sink to create the initial distance labels
static void createDistanceLabels(Graph &g, data_vec &heights, long t) {
    long level = 0;
    std::queue<long> queue;
    queue.push(t);
    heights[t] = level;

    while(!queue.empty()) {
        long curNode = queue.front();
        queue.pop();

        data_vec neighbours = g.getNeighbourListFor(curNode);
        for (auto &n : neighbours) {
            if (heights[n] == -1) {
                heights[n] = heights[curNode] + 1;
                queue.push(n);
            }
        }
    }
}

static long BFSColoring(Graph &g, data_vec &heights, data_vec &residuals, data_vec &wave, data_vec &reverse,
                        data_vec &color, long startVertex, long startLevel, long currentWave) {
    long coloredVertices = 0;
    std::queue<long> queue;
    queue.push(startVertex);
    color[startVertex] = 1;
    heights[startVertex] = startLevel;
    wave[startVertex] = currentWave;

    while(!queue.empty()) {
        long curNode = queue.front();
        queue.pop();

        data_vec neighbours = g.getNeighbourListFor(curNode);
        long edgeId;
        for (auto &n : neighbours) {
            edgeId = g.getIdFromSrcDst(n, curNode);
            edgeId = edgeId == -1 ? reverse[g.getIdFromSrcDst(curNode, n)] : edgeId;

            if (residuals[edgeId] <= 0) continue;

            if (color[n] == 0) {
                color[n] = 1;
                coloredVertices++;

                heights[n] = heights[curNode] + 1;

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
    data_vec color(g.getNV(), 0);
    std::queue<long> queue;
    long coloredVertices = BFSColoring(g, heights, residuals, wave, reverse, color, t, 0, currentWave);
    if (coloredVertices < g.getNV()) {
        BFSColoring(g, heights, residuals, wave, reverse, color, s, g.getNV(), currentWave);
    }
}

static void processWithRelabel(long node, Graph &g, std::queue<long> &queue, data_vec &heights,
                               std::vector<std::atomic_long> &excess, data_vec &reverse,
                               std::atomic_long *residuals, data_vec &wave, std::vector<std::atomic_bool> &isActive,
                               const long s, const long t) {

    if (!g.isValidNode(node)) {
        queue.pop();
        return;
    }

    while (excess[node].load() > 0) {
        long h = std::numeric_limits<long>::max();
        long nextV = -1;
        long edgeId, lowEdgeId = -1;
        data_vec neighbours = g.getNeighbourListFor(node);

        if (neighbours.empty() || node == s || node == t) {
            queue.pop();
            isActive[node].store(false);
            return;
        }

        for (auto &n: neighbours) {
            edgeId = g.getIdFromSrcDst(node, n);
            edgeId = edgeId == -1 ? reverse[g.getIdFromSrcDst(n, node)] : edgeId;

            if (residuals[edgeId].load() <= 0) continue;

            if (heights[n] < h) {
                nextV = n;
                h = heights[n];
                lowEdgeId = edgeId;
            }
        }

        if (heights[node] > h) {
            long delta = std::min(excess[node].load(), residuals[lowEdgeId].load());
            if (wave[node] <= wave[nextV] && heights[node] > heights[nextV]) {
                residuals[lowEdgeId].fetch_sub(delta);
                residuals[reverse[lowEdgeId]].fetch_add(delta);
                excess[node].fetch_sub(delta);
                excess[nextV].fetch_add(delta);

                if (!isActive[nextV].load()) {
                    isActive[nextV].store(true);
                    queue.emplace(nextV);
                }
            }
        } else if (heights[node] < h + 1){
            heights[node] = h + 1;
        }
    }

    queue.pop();
    isActive[node].store(false);
}

struct ExchangeFlag {
    std::atomic_int flag;
    int receiver;
};

static void sendActiveNodes(int from, std::vector<std::queue<long>> &queues, std::vector<ExchangeFlag> &exchange) {
    int receiver = exchange[from].receiver;
    int nodesToSend;// = std::floor(queues[from].size() / 2);
    nodesToSend = 1;
    for (int i = 0; i < nodesToSend; ++i) {
        queues[receiver].emplace(queues[from].front());
        queues[from].pop();
    }

    exchange[receiver].flag.store(0);
    exchange[from].flag.store(0);
}

static bool requestActiveNodes(int threadId, std::vector<std::queue<long>> &queues, std::vector<ExchangeFlag> &exchange) {
    std::vector<int> possibleQueues;
    bool nodesAvailable = false;

    for (size_t i = 0; i < queues.size(); ++i) {
        if (queues[i].size() >= 2 && !exchange[i].flag.load()) {
            possibleQueues.emplace_back(i);
        }

        nodesAvailable = !queues[i].empty() ? true : nodesAvailable;
    }

    if (!possibleQueues.empty()) {
        int queue = possibleQueues.front();
        int ref = 0;
        if (exchange[queue].flag.compare_exchange_strong(ref, 1)) {
            if (exchange[threadId].flag.compare_exchange_strong(ref, 1)) {
                exchange[queue].receiver = threadId;
            }
        } else {
            exchange[queue].flag.store(0);
        }
    }

    return nodesAvailable;
}

static void execute(int threadId, Graph &g, std::vector<std::queue<long>> &queues, data_vec &heights,
                    std::vector<std::atomic_long> &excess, data_vec &reverse, std::atomic_long *residuals,
                    data_vec &wave, std::vector<ExchangeFlag> &exchange, std::vector<std::atomic_bool> &isActive,
                    const long s, const long t) {

    //std::cout << "Thread[" << threadId << "]: Queue size is " << queues[threadId].size() << std::endl;
    long numIters = 0;

    while (excess[s].load() + excess[t].load() < 0) {
        if (queues[threadId].empty()) {
            bool nodesAvailable = requestActiveNodes(threadId, queues, exchange);
            if (!nodesAvailable) break;
            while (exchange[threadId].flag.load() && nodesAvailable && queues[threadId].empty()) {
                nodesAvailable = false;
                for (auto & q: queues) {
                    nodesAvailable = !q.empty() ? true : nodesAvailable;
                }
            }
        } else if (exchange[threadId].flag.load() && queues[threadId].size() > 1) {
            sendActiveNodes(threadId, queues, exchange);
        }

        if (!queues[threadId].empty()) {
            long curNode = queues[threadId].front();
            processWithRelabel(curNode, g, queues[threadId], heights, excess, reverse, residuals, wave, isActive, s, t);
            numIters++;
        }
    }

    // Cleanup

    if (exchange[threadId].flag.load()) {
        int receiver = exchange[threadId].receiver;
        exchange[receiver].flag.store(0);
        exchange[threadId].flag.store(0);
    }

    while (!queues[threadId].empty()) {
        queues[threadId].pop();
    }

    //std::cout << "Finished executing thread " << threadId << std::endl;
}

static void fillReverseVector(Graph &g, std::vector<long>& reverse) {
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

                newEdges++;
            } else if (revEdgeId == -1) {
                // The reverse edge was created when undirecting the graph
                reverse[edgeId] = g.getNE() + newEdges;
                reverse[g.getNE() + newEdges] = edgeId;

                newEdges++;
            } else {
                // Both edges exists in the original graph
                reverse[edgeId] = revEdgeId;
                reverse[revEdgeId] = edgeId;
            }

            numEdges++;
            if (numEdges % 50000 == 0) {
                std::cout << numEdges << std::endl;
            }
        }
    }
}

// A parallel lock free max flow algorithm
long PLFFlow(Graph &g, const long s, const long t) {
    std::vector<std::atomic_bool> isActive(g.getNV());
    for (auto & a: isActive) {
        a.store(false);
    }

    unsigned int numThreads = std::thread::hardware_concurrency();
    if (numThreads == 0) numThreads = DEFAULT_NUM_THREADS;
    std::vector<std::thread> threads(numThreads);
    std::vector<std::queue<long>> queues(numThreads);
    std::vector<ExchangeFlag> exchange(numThreads);
    for (auto &e : exchange) {
        e.flag.store(0);
    }

    std::cout << "Creating neighbour lists" << std::endl;
    if (!g.isNeighbourListAvailable()) {
        g.createNeighbourList(true);
    } else if (!g.isNeighbourListUndirected()) {
        std::cerr << "Graph with neighbour list supplied, but it is not undirected." << std::endl;
        return -1;
    }

    data_vec capacities = g.getWeights();
    capacities.resize(g.getUndirectedNumEdges(), 0);
    data_vec reverse(g.getUndirectedNumEdges(), -1);
    std::cout << "Filling reverse vector" << std::endl;
    fillReverseVector(g, reverse);

    auto residuals = new std::atomic_long[g.getUndirectedNumEdges()];
    for (size_t i = 0; i < capacities.size(); ++i) {
        residuals[i].store(capacities[i]);
    }

    std::vector<std::atomic_long> excess(g.getNV());
    for (auto & e : excess) {
        e.store(0);
    }

    data_vec wave(g.getNV(), 0);
    long currentWave = 0;

    data_vec heights(g.getNV(), -1);
    createDistanceLabels(g, heights, t);
    heights[s] = g.getNV();

    std::cout << "Creating initial preflow" << std::endl;
    // Initial preflow
    std::vector<long> neighbours = g.getNeighbours(s);
    for (size_t i = 0; i < neighbours.size(); ++i) {
        long n = neighbours[i];
        long edgeId = g.getIdFromSrcDst(s, n);
        residuals[edgeId].store(0);
        residuals[reverse[edgeId]].store(capacities[reverse[edgeId]] + capacities[edgeId]);
        excess[n].store(g.getWeightFromId(edgeId));
        excess[s].store(excess[s] - g.getWeightFromId(edgeId));

        queues[i % numThreads].emplace(n);
        isActive[n].store(true);
    }

    const int itersBetweenGlobalRelabel = std::floor(g.getNV() / 2);
    int itersSinceGlobalRelabel = 0;
    unsigned long totalIters = 0;

    for (int i = 0; i < threads.size(); ++i) {
        threads[i] = std::thread(execute, i, std::ref(g), std::ref(queues), std::ref(heights), std::ref(excess),
                                 std::ref(reverse), residuals, std::ref(wave), std::ref(exchange),
                                 std::ref(isActive), s, t);
    }

    for (auto &thread : threads) {
        thread.join();
    }

    long maxFlow = 0;
    neighbours = g.getNeighbours(s);
    for (auto &n : neighbours) {
        long edgeId = g.getIdFromSrcDst(s, n);
        maxFlow += (capacities[edgeId] - residuals[edgeId].load());
    }

    delete[] residuals;

    return maxFlow;
}