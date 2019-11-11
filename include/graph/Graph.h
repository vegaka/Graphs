#ifndef GRAPHS_GRAPH_H
#define GRAPHS_GRAPH_H

#include <vector>
#include <algorithm>

class Graph {
protected:
    std::vector<long> degrees;

    bool symmetric;
    unsigned long nv, ne;

    // Neighbour lists for Max Flow
    bool hasNeighbourList = false;
    bool undirectedNeighbourList = false;
    unsigned long undirectedNE;
    std::vector<std::vector<long>> neighbourList;

public:
    long getMaxDegreeNode() {
        if (degrees.empty()) return -1;
        auto maxElement = std::max_element(degrees.begin(), degrees.end());
        return std::distance(degrees.begin(), maxElement);
    }

    unsigned long getNV() {
        return nv;
    }

    unsigned long getNE() {
        return ne;
    }

    long getDegree(unsigned long idx) {
        return degrees.at(idx);
    }

    virtual long getEdgeWeight(long src, long dst)= 0;

    virtual std::vector<long> getNeighbours(long src, bool useNeighbourList)= 0;

    virtual std::pair<long, long> getSrcDstFromId(long id)= 0;

    virtual long getIdFromSrcDst(long src, long dst)= 0;

    virtual bool hasEdge(long src, long dst)= 0;

    virtual long getWeightFromId(long id)= 0;

    virtual std::vector<long> getWeights()= 0;

    virtual void createNeighbourList(bool undirected)= 0;

    bool isNeighbourListUndirected() {
        return undirectedNeighbourList;
    }

    bool isNeighbourListAvailable()  {
        return hasNeighbourList;
    }

    unsigned long getUndirectedNumEdges() {
        return undirectedNE;
    }

    std::vector<long> getNeighbourListFor(long node) {
        if (hasNeighbourList) {
            return neighbourList.at(node);
        }

        return std::vector<long>(0);
    }
};


#endif //GRAPHS_GRAPH_H
