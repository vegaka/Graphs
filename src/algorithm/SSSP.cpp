#include <queue>
#include <iostream>
#include "algorithm/SSSP.h"

static bool relax(std::vector<long> &dist, std::vector<long> &par, const long src, const long dst, const long w) {
    if (dist[src] == std::numeric_limits<long>::max()) return false;
    else if (dist[dst] > dist[src] + w) {
        dist[dst] = dist[src] + w;
        par[dst] = src;
        return true;
    } else return false;
}

std::pair<std::vector<long>, std::vector<long>> Dijkstra(Graph &g, const long start) {
    //std::cout << "Start node: " << start << std::endl;
    std::vector<long> distances(g.getNV(), std::numeric_limits<long>::max());
    distances[start] = 0;

    std::vector<long> parent(g.getNV(), std::numeric_limits<long>::max());

    auto cmp = [&distances](long n1, long n2) { return distances[n1] > distances[n2]; };
    std::priority_queue<long, std::vector<long>, decltype(cmp)> queue(cmp);

    queue.push(start);

    while (!queue.empty()) {
        const long curSrc = queue.top();
        queue.pop();

        std::vector<long> neighbours = g.getNeighbours(curSrc, false);

        for (long dst: neighbours) {
            long w = g.getEdgeWeight(curSrc, dst);
            if (relax(distances, parent, curSrc, dst, w)) {
                queue.push(dst);
            }
        }
    }

    return std::make_pair(distances, parent);
}

std::pair<std::vector<long>, std::vector<long>> BellmanFord(Graph &g, const long start) {
    std::vector<long> distances(g.getNV(), std::numeric_limits<long>::max());
    distances[start] = 0;

    std::vector<long> parent(g.getNV(), std::numeric_limits<long>::max());

    for (long i = 0; i < g.getNV(); ++i) {
        for (long j = 0; j < g.getNV(); ++j) {
            std::vector<long> neighbours = g.getNeighbours(j, false);

            for (long dst: neighbours) {
                long w = g.getEdgeWeight(j, dst);
                relax(distances, parent, j, dst, w);
            }
        }
    }

    return std::make_pair(distances, parent);
}