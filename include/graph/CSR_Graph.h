#ifndef GRAPHS_CSR_GRAPH_H
#define GRAPHS_CSR_GRAPH_H

#include <string>
#include <vector>
#include <tuple>
#include "Graph.h"

class CSR_Graph : public Graph {
private:
    typedef std::tuple<long, long, long> Edge;

    std::vector<long> w;
    std::vector<long> rowPtr;
    std::vector<long> colIdx;

    void edgeListToCSR(std::vector<Edge>);

public:
    // Path points to a matrix market file
    explicit CSR_Graph(std::string path);

    long getEdgeWeight(long src, long dst) override;

    std::vector<long> getNeighbours(long src) override;

};


#endif //GRAPHS_GRAPH_H
