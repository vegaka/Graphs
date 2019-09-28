#include "algorithm/VertexCover.h"

// 2-approximation for Vertex Cover
std::unordered_set<long> VertexCover(Graph &g) {
    std::unordered_set<long> cover;

    for (long i = 0; i < g.getNE(); i++) {
        auto vertices = g.getSrcDstFromId(i);
        auto src = vertices.first;
        auto dst = vertices.second;

        if (cover.find(src) == cover.end() && cover.find(dst) == cover.end()) {
            cover.insert(src);
            cover.insert(dst);
        }
    }

    return cover;
}
