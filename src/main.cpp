#include <iostream>
#include <queue>
#include <chrono>
#include <algorithm/VertexCover.h>
#include <sstream>
#include <algorithm/MaxFlow.h>
#include <algorithm/BFS.h>
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

static void testMF(int argc, char* argv[]) {
    if (argc == 2) {
        CSR_Graph g {argv[1]};
        long src = 0;
        long sink = 7999;

        int iterations = 10;
        std::vector<long> flows(iterations, 0);
        std::vector<long> times(iterations);
        for (size_t i = 0; i < flows.size(); ++i) {
            auto pstartTime = std::chrono::high_resolution_clock::now();
            std::cout << "Parallel Max Flow run: " << i << std::endl;
            flows[i] = PLFFlow(g, src, sink);
            auto pendTime = std::chrono::high_resolution_clock::now();
            times[i] = std::chrono::duration_cast<std::chrono::milliseconds>(pendTime - pstartTime).count();
        }

        auto startTime = std::chrono::high_resolution_clock::now();
        auto maxFlow = LFFlow(g, src, sink, true);
        auto endTime = std::chrono::high_resolution_clock::now();

        long errors = 0;
        for (auto & f: flows) {
            if (f != maxFlow) {
                errors++;
            }
        }

        std::cout << "=======================================" << std::endl;
        std::cout << "(Parallel) Flows: [";
        for (size_t i = 0; i < flows.size() - 1; ++i) {
            std::cout << flows[i] << ", ";
        }
        std::cout << flows[flows.size() - 1] << "]" << std::endl;
        std::cout << "(Parallel) Errors: " << errors << std::endl;

        std::cout << "(Parallel) Times: [";
        for (size_t i = 0; i < times.size() - 1; ++i) {
            std::cout << times[i] << ", ";
        }
        std::cout << times[times.size() - 1] << "]" << std::endl;

        //std::cout << "(Parallel) The maximum flow is: " << flows[0] << "." << std::endl;
        //std::cout << "(Parallel) Time spent (ms): "
        //          << std::chrono::duration_cast<std::chrono::milliseconds>(pendTime - pstartTime).count()
        //          << std::endl;

        std::cout << "(Single) The maximum flow is: " << maxFlow << "." << std::endl;
        std::cout << "(Single) Time spent (ms): "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()
                  << std::endl;
    }
}

static void testBFS(int argc, char* argv[]) {
    if (argc == 2) {
        CSR_Graph g {argv[1]};
        long src = 0;

        auto startTime = std::chrono::high_resolution_clock::now();
        auto distances = BFS(g, src);
        auto endTime = std::chrono::high_resolution_clock::now();

        std::cout << "Distances:" << std::endl;
        for (int i = 0; i < distances.size(); ++i) {
            std::cout << i << " : " << distances[i] << std::endl;
        }

        std::cout << "Time spent (ms): "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()
                  << std::endl;
    }
}

int main(int argc, char* argv[]) {

    //testSSSP(argc, argv);
    //testVertexCover(argc, argv);
    testMF(argc, argv);
    //testBFS(argc, argv);

    return 0;
}
