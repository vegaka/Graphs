#include <fstream>
#include <iostream>
#include "graph/CSR_Graph.h"

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
    // src = row, dst = col
    //std::cout << "src: " << src << ", dst: " << dst << std::endl;
    for (long i = rowPtr[src]; i < rowPtr[src + 1]; ++i) {
        //std::cout << "Weight: " << w[i] << ", colIdx: " << colIdx[i] << std::endl;
        if (colIdx[i] == dst) {
            return w[i];
        }
    }

    return std::numeric_limits<int>::max();
}

std::vector<long> CSR_Graph::getNeighbours(long src) {
    std::vector<long> neighbours;
    for (long i = rowPtr[src]; i < rowPtr[src + 1]; ++i) {
        neighbours.push_back(colIdx[i]);
    }

    return neighbours;
}

