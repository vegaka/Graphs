#include <queue>
#include <iostream>
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