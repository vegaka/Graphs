#include <queue>
#include "algorithm/BFS.h"

std::vector<long> BFS(Graph &g, const long src) {
    std::queue<long> q;
    std::vector<long> distances(g.getNV(), 0);

    q.push(src);

    while (!q.empty()) {
        long curNode = q.front();
        q.pop();

        auto neighbours = g.getNeighbours(curNode);
        for (auto &n : neighbours) {
            if (distances[n] != 0) {
                distances[n] = distances[curNode] + 1;
                q.push(n);
            }
        }
    }

    return distances;
}