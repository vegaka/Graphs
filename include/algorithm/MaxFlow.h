#ifndef GRAPHS_MAXFLOW_H
#define GRAPHS_MAXFLOW_H

#include "graph/Graph.h"

std::vector<long> MaxFlow(Graph &g, long s, long t);
std::pair<std::vector<long>, std::vector<long>> LFFlow(Graph &g, long s, long t);

#endif //GRAPHS_MAXFLOW_H
