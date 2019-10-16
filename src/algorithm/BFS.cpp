#include <queue>
#include "algorithm/BFS.h"

std::vector<long> BFS(Graph &g, const long src) {
    std::queue<long> q;
    std::vector<long> distances(g.getNV(), -1);

    q.push(src);
    distances[src] = 0;

    while (!q.empty()) {
        long curNode = q.front();
        q.pop();

        auto neighbours = g.getNeighbours(curNode);
        for (auto &n : neighbours) {
            if (distances[n] == -1) {
                distances[n] = distances[curNode] + 1;
                q.push(n);
            }
        }
    }

    return distances;
}