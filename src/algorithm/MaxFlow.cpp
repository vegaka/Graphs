#include <queue>
#include <iostream>
#include <graph/CSR_Graph.h>
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

using data_vec = std::vector<long>;

static void process(long node, Graph &g, data_vec &heights, data_vec &excess, flow_vec &flow,
        std::queue<long> &activeNodes, const long s, const long t, std::vector<bool> &isActive) {
    std::cout << "Processing node: " << node << std::endl;
    if (excess[node] > 0) {
        long e = excess[node];
        long h = std::numeric_limits<long>::max();
        long nextV = -1;
        data_vec neighbours = g.getNeighbours(node);

        if (neighbours.empty()) return;

        for (auto &n: neighbours) {
            if (g.getEdgeWeight(node, n) - flow[node][n] <= 0) continue;

            if (heights[n] < h) {
                nextV = n;
                h = heights[n];
            }
        }

        if (heights[node] > h) {
            long delta = std::min(e, g.getEdgeWeight(node, nextV) - flow[node][nextV]);
            flow[node][nextV] += delta;
            flow[nextV][node] -= delta;
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

static CSR_Graph createResidualGraph(Graph &g) {
    std::cout << "Creating residual graph." << std::endl;
    std::unordered_map<std::pair<long, long>, long> edges;

    for (int i = 0; i < g.getNE(); ++i) {
        auto srcdst = g.getSrcDstFromId(i);

        if (edges.find(srcdst) != edges.end()) {
            edges[srcdst] += g.getWeightFromId(i);
        } else {
            edges.insert(srcdst, g.getWeightFromId(i));
        }

        long src = srcdst.first;
        long dst = srcdst.second;
        auto revedge = std::make_pair(dst, src);

        if (edges.find(revedge) != edges.end()) {
            edges[revedge] += g.getWeightFromId(i);
        } else {
            edges.insert(revedge, g.getWeightFromId(i));
        }
    }

    return CSR_Graph(edges);
}

// A lock free max flow algorithm
flow_vec LFFlow(Graph &g, const long s, const long t) {
    std::queue<long> activeNodes;
    std::vector<bool> isActive(g.getNV(), false);

    auto residualGraph = createResidualGraph(g);

    data_vec capacities = residualGraph.getWeights();
    data_vec residuals(residualGraph.getNE(), 0);

    data_vec heights(g.getNV(), 0);
    data_vec excess(g.getNV(), 0);
    flow_vec flow(g.getNE(), std::vector<long>(g.getNE(), 0));

    heights[s] = g.getNV();

    // Initial preflow
    std::vector<long> neighbours = g.getNeighbours(s);
    for (auto &n : neighbours) {
        flow[s][n] = g.getEdgeWeight(s, n);
        flow[n][s] = -g.getEdgeWeight(s, n);
        excess[n] = g.getEdgeWeight(s, n);
        excess[s] = excess[s] - g.getEdgeWeight(s, n);

        activeNodes.push(n);
        isActive[n] = true;
    }

    std::cout << activeNodes.size() << std::endl;

    while (excess[s] + excess[t] < 0) {
        long curnode = activeNodes.front();

        process(curnode, g, heights, excess, flow, activeNodes, s, t, isActive);
    }

    return flow;
}