#include <queue>
#include <iostream>
#include "algorithm/MaxFlow.h"
#include "algorithm/BFS.h"

static void preflow(Graph &g, std::vector<long> &flows, std::vector<long> &excess, const long s, std::queue<long> &activeNodes) {
    auto neighbours = g.getNeighbours(s, false);
    for (auto &n : neighbours) {
        long id = g.getIdFromSrcDst(s, n);
        flows[id] = g.getWeightFromId(id);
        excess[n] += flows[id];
        activeNodes.push(n);
    }
}

static void push(Graph &g, std::vector<long> &flows, std::vector<long> &excess, std::vector<long> &distances,
                 std::queue<long> &activeNodes, const long curNode, const long s, const long t) {
    auto neighbours = g.getNeighbours(curNode, false);
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

using data_vec = std::vector<long>;

static void process(long node, Graph &g, data_vec &heights, data_vec &excess, data_vec &capacities, data_vec &reverse,
        data_vec &residuals, std::queue<long> &activeNodes, std::vector<bool> &isActive, const long s, const long t) {
    std::cout << "Processing node: " << node << std::endl;
    if (excess[node] > 0) {
        long e = excess[node];
        long h = std::numeric_limits<long>::max();
        long nextV = -1;
        long edgeId = -1;
        data_vec neighbours = g.getNeighbours(node, true);

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

            if (nextV != t && !isActive[nextV]) {
                activeNodes.push(nextV);
                isActive[nextV] = true;
            }
        } else {
            heights[node] = h + 1;
        }
    } else {
        activeNodes.pop();
        isActive[node] = false;
    }
}

static void fillReverseAndCapacityVectors(Graph &g, std::vector<long>& reverse, std::vector<long>& capacities) {
    long newEdges = 0;

    for (long i = 0; i < g.getNE(); ++i) {
        if (reverse[i] != -1) continue;

        auto srcdst = g.getSrcDstFromId(i);
        long src = srcdst.first;
        long dst = srcdst.second;

        if (!g.hasEdge(dst, src)) {
            reverse[i] = g.getNE() + newEdges;
            reverse[g.getNE() + newEdges] = i;

            capacities[g.getNE() + newEdges] = g.getWeightFromId(i);

            newEdges++;
        } else {
            long id = g.getIdFromSrcDst(dst, src);
            reverse[id] = i;
            reverse[i] = id;

            capacities[i] += g.getWeightFromId(id);
            capacities[id] += g.getWeightFromId(i);
        }
    }
}

// A lock free max flow algorithm
std::pair<std::vector<long>, std::vector<long>> LFFlow(Graph &g, const long s, const long t) {
    std::queue<long> activeNodes;
    std::vector<bool> isActive(g.getNV(), false);

    g.createNeighbourList(true);

    data_vec capacities = g.getWeights();
    capacities.resize(g.getUndirectedNumEdges(), 0);
    data_vec reverse(g.getUndirectedNumEdges(), -1);
    fillReverseAndCapacityVectors(g, reverse, capacities);

    data_vec residuals = g.getWeights();
    residuals.resize(g.getUndirectedNumEdges(), 0);

    data_vec heights(g.getNV(), 0);
    data_vec excess(g.getNV(), 0);

    heights[s] = g.getNV();

    // Initial preflow
    std::vector<long> neighbours = g.getNeighbours(s, true);
    for (auto &n : neighbours) {
        long edgeId = g.getIdFromSrcDst(s, n);
        residuals[edgeId] = 0;
        residuals[reverse[edgeId]] = capacities[edgeId];
        excess[n] = g.getWeightFromId(edgeId);
        excess[s] = excess[s] - g.getWeightFromId(edgeId);

        activeNodes.push(n);
        isActive[n] = true;
    }

    std::cout << activeNodes.size() << std::endl;

    while (excess[s] + excess[t] < 0) {
        long curnode = activeNodes.front();

        process(curnode, g, heights, excess, capacities, reverse, residuals, activeNodes, isActive, s, t);
    }

    return std::make_pair(capacities, residuals);
}