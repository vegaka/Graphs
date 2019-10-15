#include <queue>
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

static void push(Graph &g, std::vector<long> &flows, std::vector<long> excess, std::vector<long> &distances, const long curNode) {
    auto neighbours = g.getNeighbours(curNode);
    for (auto &n : neighbours) {
        if (distances[curNode] == distances[n] + 1) {
            // Push
            long id = g.getIdFromSrcDst(curNode, n);
            long res = g.getWeightFromId(id) - flows[id];
            long delta = std::min(excess[curNode], res);
            excess[curNode] -= delta;
            excess[n] += delta;
            flows[id] += delta;
        } else {
            // Relabel

        }
    }
}

// Implementation of the Preflow-Push algorithm
void MaxFlow(Graph &g, const long s, const long t) {
    std::vector<long> flows(g.getNE(), 0);
    std::vector<long> distances = BFS(g, t);
    std::vector<long> excess(g.getNV(), 0);
    std::queue<long> activeNodes;

    preflow(g, flows, excess, s, activeNodes);
    distances[s] = g.getNV();

    while (!activeNodes.empty()) {
        long curNode = activeNodes.front();
        activeNodes.pop();

        push(g, flows, excess, distances, curNode);
    }
}