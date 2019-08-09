#ifndef GRAPHS_GRAPH_H
#define GRAPHS_GRAPH_H

#include <vector>

class Graph {
protected:
    std::vector<long> degrees;

    bool symmetric;
    unsigned long nv, ne;

public:
    long getMaxDegreeNode() {
        auto maxElement = std::max_element(degrees.begin(), degrees.end());
        return std::distance(degrees.begin(), maxElement);
    }

    unsigned long getNV() {
        return nv;
    }

    long getDegree(unsigned long idx) {
        return degrees.at(idx);
    }

    virtual long getEdgeWeight(long src, long dst)= 0;

    virtual std::vector<long> getNeighbours(long src)= 0;
};


#endif //GRAPHS_GRAPH_H
