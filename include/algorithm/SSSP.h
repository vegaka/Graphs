#ifndef GRAPHS_SSSP_H
#define GRAPHS_SSSP_H

#include <graph/Graph.h>

std::pair<std::vector<long>, std::vector<long>> Dijkstra(Graph &g, long start);

std::pair<std::vector<long>, std::vector<long>> BellmanFord(Graph &g, long start);

#endif //GRAPHS_SSSP_H
