#include "algorithm/MaxFlow.h"
#include "algorithm/BFS.h"

static void preflow(std::vector<long> &flows, const long s) {

}

// Implementation of the Preflow-Push algorithm
void MaxFlow(Graph &g, const long s, const long t) {
    std::vector<long> flows(g.getNE(), 0);
    std::vector<long> distances = BFS(g, t);


    preflow(flows, s);

}