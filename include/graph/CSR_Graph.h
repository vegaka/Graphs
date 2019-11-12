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

    void edgeListToCSR(const std::vector<Edge>&);

public:

    // Path points to a matrix market file
    explicit CSR_Graph(const std::string& path);

    long getEdgeWeight(long src, long dst) override;

    std::vector<long> getNeighbours(long src) override;

    std::pair<long, long> getSrcDstFromId(long id) override;

    long getIdFromSrcDst(long src, long dst) override;

    bool hasEdge(long src, long dst) override ;

    long getWeightFromId(long id) override;
    
    std::vector<long> getWeights() override;

    void createNeighbourList(bool undirected) override;

    friend std::ostream& operator<<(std::ostream& os, CSR_Graph g);
};


#endif //GRAPHS_GRAPH_H
