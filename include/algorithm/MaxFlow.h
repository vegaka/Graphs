#ifndef GRAPHS_MAXFLOW_H
#define GRAPHS_MAXFLOW_H

#include "graph/Graph.h"

using flow_vec = std::vector<std::vector<long>>;

std::vector<long> MaxFlow(Graph &g, long s, long t);
flow_vec LFFlow(Graph &g, long s, long t);

#endif //GRAPHS_MAXFLOW_H
