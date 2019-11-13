#ifndef GRAPHS_MAXFLOW_H
#define GRAPHS_MAXFLOW_H

#include <atomic>
#include "graph/Graph.h"

using data_vec = std::vector<long>;

std::vector<long> MaxFlow(Graph &g, long s, long t);
long LFFlow(Graph &g, long s, long t, bool globalRelabeling);
long PLFFlow(Graph &g, long s, long t);

#endif //GRAPHS_MAXFLOW_H
