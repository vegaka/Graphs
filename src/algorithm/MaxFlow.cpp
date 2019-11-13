#include <queue>
#include <iostream>
#include <cmath>
#include "algorithm/MaxFlow.h"
#include "algorithm/BFS.h"

static void preflow(Graph &g, std::vector<long> &flows, std::vector<long> &excess, const long s, std::queue<long> &activeNodes) {
    auto neighbours = g.getNeighbours(s);
    for (auto &n : neighbours) {
        long id = g.getIdFromSrcDst(s, n);
        flows[id] = g.getWeightFromId(id);
        excess[n] += flows[id];
        activeNodes.push(n);
    }
}

static void push(Graph &g, std::vector<long> &flows, std::vector<long> &excess, std::vector<long> &distances,
                 std::queue<long> &activeNodes, const long curNode, const long s, const long t) {
    auto neighbours = g.getNeighbours(curNode);
    for (auto &n : neighbours) {
        if (distances[curNode] == distances[n] + 1) {

            // Check that we have excess flow to push
            if (excess[curNode] == 0) break;

            long id = g.getIdFromSrcDst(curNode, n);
            long res = g.getWeightFromId(id) - flows[id];
            long delta = std::min(excess[curNode], res);
            excess[curNode] -= delta;
            excess[n] += delta;
            flows[id] += delta;
            if (n != s && n!= t) activeNodes.push(n);
        }
    }

    if (excess[curNode] > 0) {
        // Relabel
        long newDist = std::numeric_limits<long>::max();
        for (auto &n : neighbours) {
            long id = g.getIdFromSrcDst(curNode, n);
            if (g.getWeightFromId(id) - flows[id] > 0) {
                newDist = std::min(newDist, distances[n] + 1);
            }
        }

        distances[curNode] = newDist;
        activeNodes.push(curNode);
    }
}

// Implementation of the Preflow-Push algorithm
std::vector<long> MaxFlow(Graph &g, const long s, const long t) {
    std::cout << "Init" << std::endl;
    std::vector<long> flows(g.getNE(), 0);
    std::vector<long> distances = BFS(g, t);
    std::vector<long> excess(g.getNV(), 0);
    std::queue<long> activeNodes;


    preflow(g, flows, excess, s, activeNodes);

    std::cout << "Preflow:" << std::endl;
    for (int i = 0; i < flows.size(); ++i) {
        std::cout << i << " : " << flows[i] << std::endl;
    }

    distances[s] = g.getNV();

    std::cout << "Distances:" << std::endl;
    for (int i = 0; i < distances.size(); ++i) {
        std::cout << i << " : " << distances[i] << std::endl;
    }

    while (!activeNodes.empty()) {
        long curNode = activeNodes.front();
        activeNodes.pop();

        push(g, flows, excess, distances, activeNodes, curNode, s, t);
    }

    return flows;
}

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

static void processWithRelabel(long node, Graph &g, data_vec &heights, data_vec &excess,
        data_vec &reverse, data_vec &residuals, data_vec &wave, std::queue<long> &activeNodes,
        std::vector<bool> &isActive, const long s, const long t) {
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

            if (residuals[edgeId] <= 0) continue;

            if (heights[n] < h) {
                nextV = n;
                h = heights[n];
            }
        }

        edgeId = g.getIdFromSrcDst(node, nextV);
        edgeId = edgeId == -1 ? reverse[g.getIdFromSrcDst(nextV, node)] : edgeId;

        if (heights[node] > h) {
            long delta = std::min(e, residuals[edgeId]);
            if (wave[node] <= wave[nextV] && heights[node] > heights[nextV]) {
                residuals[edgeId] -= delta;
                residuals[reverse[edgeId]] += delta;
                excess[node] -= delta;
                excess[nextV] += delta;

                if (nextV != t && nextV != s && !isActive[nextV]) {
                    activeNodes.push(nextV);
                    isActive[nextV] = true;
                }
            }
        } else if (heights[node] < h + 1){
            heights[node] = h + 1;
        }
    }
    activeNodes.pop();
    isActive[node] = false;
}

static void process(long node, Graph &g, data_vec &heights, data_vec &excess, data_vec &reverse, data_vec &residuals,
                    std::queue<long> &activeNodes, std::vector<bool> &isActive, const long s, const long t) {
    //std::cout << "Processing node: " << node << std::endl;
    while (excess[node] > 0) {
        long e = excess[node];
        long h = std::numeric_limits<long>::max();
        long nextV = -1;
        long edgeId = -1;
        data_vec neighbours = g.getNeighbourListFor(node);

        if (neighbours.empty()) return;

        for (auto &n: neighbours) {
            edgeId = g.getIdFromSrcDst(node, n);
            edgeId = edgeId == -1 ? reverse[g.getIdFromSrcDst(n, node)] : edgeId;

            if (residuals[edgeId] <= 0) continue;

            if (heights[n] < h) {
                nextV = n;
                h = heights[n];
            }
        }

        edgeId = g.getIdFromSrcDst(node, nextV);
        edgeId = edgeId == -1 ? reverse[g.getIdFromSrcDst(nextV, node)] : edgeId;

        if (heights[node] > h) {
            long delta = std::min(e, residuals[edgeId]);
            residuals[edgeId] -= delta;
            residuals[reverse[edgeId]] += delta;
            excess[node] -= delta;
            excess[nextV] += delta;

            if (nextV != t && nextV != s && !isActive[nextV]) {
                activeNodes.push(nextV);
                isActive[nextV] = true;
            }
        } else {
            heights[node] = h + 1;
        }
    }
    activeNodes.pop();
    isActive[node] = false;
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

// A lock free max flow algorithm
long LFFlow(Graph &g, const long s, const long t, bool globalRelabeling) {
    std::queue<long> activeNodes;
    std::vector<bool> isActive(g.getNV(), false);

    std::cout << "Creating neighbour lists" << std::endl;
    g.createNeighbourList(true);

    data_vec capacities = g.getWeights();
    capacities.resize(g.getUndirectedNumEdges(), 0);
    data_vec reverse(g.getUndirectedNumEdges(), -1);
    std::cout << "Filling reverse vector" << std::endl;
    fillReverseVector(g, reverse);

    data_vec residuals = capacities;

    data_vec excess(g.getNV(), 0);

    data_vec wave(g.getNV(), 0);
    long currentWave = 0;

    data_vec heights(g.getNV(), -1);
    createDistanceLabels(g, heights, t);
    heights[s] = g.getNV();

    std::cout << "Creating initial preflow" << std::endl;
    // Initial preflow
    std::vector<long> neighbours = g.getNeighbours(s);
    for (auto &n : neighbours) {
        long edgeId = g.getIdFromSrcDst(s, n);
        residuals[edgeId] = 0;
        residuals[reverse[edgeId]] = capacities[reverse[edgeId]] + capacities[edgeId];
        excess[n] = g.getWeightFromId(edgeId);
        excess[s] = excess[s] - g.getWeightFromId(edgeId);

        activeNodes.push(n);
        isActive[n] = true;
    }

    std::cout << "Active initial nodes: " << activeNodes.size() << std::endl;

    const int itersBetweenGlobalRelabel = std::floor(g.getNV() / 2);
    int itersSinceGlobalRelabel = 0;
    unsigned long totalIters = 0;

    while (excess[s] + excess[t] < 0) {
        long curnode = activeNodes.front();

        //std::cout << "Processing node: " << curnode << ", numActive: " << activeNodes.size() << std::endl;

        if (globalRelabeling)
            processWithRelabel(curnode, g, heights, excess, reverse, residuals, wave, activeNodes, isActive, s, t);
        else
            process(curnode, g, heights, excess, reverse, residuals, activeNodes, isActive, s, t);

        itersSinceGlobalRelabel++;
        if (globalRelabeling && itersSinceGlobalRelabel >= itersBetweenGlobalRelabel) {
            globalRelabel(g, heights, residuals, wave, reverse, s, t, currentWave);
            //std::cout << "Finished global relabel" << std::endl;
            itersSinceGlobalRelabel = 0;
        }

        totalIters++;
        if (totalIters % 10000000 == 0) {
            std::cout << "Total iterations: " << totalIters << std::endl;
        }
    }

    long maxFlow = 0;
    neighbours = g.getNeighbours(s);
    for (auto &n : neighbours) {
        long edgeId = g.getIdFromSrcDst(s, n);
        maxFlow += (capacities[edgeId] - residuals[edgeId]);
    }

    return maxFlow;
}