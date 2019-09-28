#include <queue>
#include "algorithm/MaxFlow.h"
#include "algorithm/BFS.h"

static void preflow(Graph &g, std::vector<long> &flows, const long s) {
    auto neighbours = g.getNeighbours(s);
    for (auto &n : neighbours) {
        long id = g.getIdFromSrcDst(s, n);
        flows[id] = g.getWeightFromId(id);
    }
}

static void push() {

}

static void relabel() {

}

// Implementation of the Preflow-Push algorithm
void MaxFlow(Graph &g, const long s, const long t) {
    std::vector<long> flows(g.getNE(), 0);
    std::vector<long> distances = BFS(g, t);
    std::queue<long> activeNodes;

    preflow(g, flows, s);
    distances[s] = g.getNV();




}