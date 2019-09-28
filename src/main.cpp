#include <iostream>
#include <queue>
#include <chrono>
#include <algorithm/VertexCover.h>
#include <sstream>
#include "algorithm/SSSP.h"
#include "graph/CSR_Graph.h"

static void testSSSP(int argc, char* argv[]) {
    if (argc == 2) {
        CSR_Graph g {argv[1]};
        long startNode = g.getMaxDegreeNode();

        auto startTime = std::chrono::high_resolution_clock::now();
        auto result = Dijkstra(g, startNode);
        auto endTime = std::chrono::high_resolution_clock::now();

        auto distances = result.first;
        auto parent = result.second;
        for (int i = 0; i < distances.size(); ++i) {
            if (distances[i] == std::numeric_limits<long>::max()) continue;
            std::cout << "Distance to node " << i << ": " << distances[i] << std::endl;
        }

        for (int i = 0; i < parent.size(); ++i) {
            if (i == startNode || parent[i] == std::numeric_limits<long>::max()) continue;
            long curNode = parent[i];
            std::cout << i << "->" << startNode << ": ";
            while (curNode != startNode) {
                std::cout << curNode << ", ";
                curNode = parent[curNode];
            }
            std::cout << curNode << std::endl;
        }

        std::cout << "Time spent (ms): "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()
                  << std::endl;
    } else {
        for (int i = 0; i < argc; i++) {
            std::cout << argv[i] << std::endl;
        }
    }
}

static void testVertexCover(int argc, char* argv[]) {
    if (argc == 2) {
        CSR_Graph g {argv[1]};

        auto startTime = std::chrono::high_resolution_clock::now();
        auto result = VertexCover(g);
        auto endTime = std::chrono::high_resolution_clock::now();

        std::cout << "Vertices in cover: ";
        std::ostringstream s;
        for (auto &n : result) {
            s << n << ", ";
        }
        std::cout << s.str().substr(0, s.str().size() - 2) << std::endl;

        std::cout << "Time spent (ms): "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()
                  << std::endl;
    }
}

int main(int argc, char* argv[]) {

    //testSSSP(argc, argv);
    testVertexCover(argc, argv);

    return 0;
}