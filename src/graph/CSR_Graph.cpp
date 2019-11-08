#include <fstream>
#include <iostream>
#include "graph/CSR_Graph.h"

CSR_Graph::CSR_Graph(const long nv, const long ne, const bool sym) {
    this->ne = ne;
    this->nv = nv;
    symmetric = sym;
}

CSR_Graph::CSR_Graph(const std::unordered_map<std::pair<long, long>, long>& edges) {
    std::vector<Edge> edgeList;
    nv = 0;

    for (auto const& e : edges) {
        long src = e.first.first;
        long dst = e.first.second;
        long weight = e.second;

        edgeList.emplace_back(src, dst, weight);

        nv = nv > src ? nv : src;
        nv = nv > dst ? nv : dst;
    }

    symmetric = true;
    ne = edgeList.size();

    std::sort(edgeList.begin(), edgeList.end());
    edgeListToCSR(edgeList);
}

CSR_Graph::CSR_Graph(const std::string& path) {
    std::ifstream ifs {path};
    std::string line;
    if (ifs) {
        // Header
        getline(ifs, line);
        symmetric = line.find("symmetric") != std::string::npos;

        // Skip comments
        while (getline(ifs, line)) {
            if (line.front() == '%') {
                continue;
            } else {
                // Line now contains the dimensions
                break;
            }
        }

        nv = std::stoul(line);
        std::vector<long> degs(nv);
        std::vector<Edge> edges;

        while (getline(ifs, line)) {
            // The MM file is assumed to be 1-indexed
            long src = std::stol(line) - 1;
            auto spacepos = line.find(' ');
            long dst = std::stol(line.substr(spacepos)) - 1;
            auto lspos = line.find_last_of(' ');
            long data = std::stol(line.substr(lspos));
            edges.emplace_back(src, dst, data);
            degs[src]++;

            if (symmetric) {
                edges.emplace_back(dst, src, data);
                degs[dst]++;
            }
        }

        ne = edges.size();
        degrees = degs;

        std::sort(edges.begin(), edges.end());

        edgeListToCSR(edges);

    } else {
        std::cerr << "Could not open file: " << path << "!" << std::endl;
    }
}

void CSR_Graph::edgeListToCSR(const std::vector<Edge>& edges) {
    rowPtr.emplace_back(0);

    for (Edge e: edges) {
        w.emplace_back(std::get<2>(e));
        colIdx.emplace_back(std::get<1>(e));
    }

    for (int i = 1; i <= degrees.size(); ++i) {
        rowPtr.emplace_back(rowPtr[i - 1] + degrees[i - 1]);
    }
}

long CSR_Graph::getEdgeWeight(long src, long dst) {
    long id = getIdFromSrcDst(src, dst);
    if (id != -1) return getWeightFromId(id);
    else return std::numeric_limits<int>::max();
}

std::vector<long> CSR_Graph::getNeighbours(long src) {
    std::vector<long> neighbours;
    for (long i = rowPtr[src]; i < rowPtr[src + 1]; ++i) {
        neighbours.push_back(colIdx[i]);
    }

    return neighbours;
}

std::pair<long, long> CSR_Graph::getSrcDstFromId(long id) {
    // Find src from id
    long src = 0;
    for (auto i = rowPtr.size() - 1; i >= 0; --i) {
        if (rowPtr[i] > id) {
            continue;
        } else {
            src = i;
            break;
        }
    }

    long dst = colIdx[id];
    return std::make_pair(src, dst);
}

long CSR_Graph::getIdFromSrcDst(long src, long dst) {
    for (long i = rowPtr[src]; i < rowPtr[src + 1]; ++i) {
        if (colIdx[i] == dst) {
            return i;
        }
    }

    return -1;
}

long CSR_Graph::getWeightFromId(long id) {
    return w.at(id);
}

std::vector<long> CSR_Graph::getWeights() {
    return w;
}
