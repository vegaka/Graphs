#ifndef GRAPHS_MAXFLOW_H
#define GRAPHS_MAXFLOW_H

#include "graph/Graph.h"

using data_vec = std::vector<long>;

std::vector<long> MaxFlow(Graph &g, long s, long t);
std::pair<data_vec, data_vec> LFFlow(Graph &g, long s, long t, bool globalRelabeling);
std::pair<data_vec, data_vec> PLFFlow(Graph &g, long s, long t);

#endif //GRAPHS_MAXFLOW_H
