#include <algorithm/MaxFlow.h>
#include <queue>
#include <iostream>
#include <thread>

#define DEFAULT_NUM_THREADS 4

using data_vec = std::vector<long>;

static void process(long node, Graph &g, data_vec &heights, data_vec &excess, data_vec &capacities, data_vec &reverse,
                    data_vec &residuals, std::queue<long> &activeNodes, std::vector<bool> &isActive, const long s, const long t) {
    std::cout << "Processing node: " << node << std::endl;
    if (excess[node] > 0) {
        long e = excess[node];
        long h = std::numeric_limits<long>::max();
        long nextV = -1;
        long edgeId = -1;
        data_vec neighbours = g.getNeighbourListFor(node);

        if (neighbours.empty()) return;

        for (auto &n: neighbours) {
            edgeId = g.getIdFromSrcDst(node, n);
            edgeId = edgeId == -1 ? reverse[g.getIdFromSrcDst(n, node)] : edgeId;

            if (residuals[edgeId] <= 0) continue;

            if (heights[n] < h) {
                nextV = n;
                h = heights[n];
            }
        }

        edgeId = g.getIdFromSrcDst(node, nextV);
        edgeId = edgeId == -1 ? reverse[g.getIdFromSrcDst(nextV, node)] : edgeId;

        if (heights[node] > h) {
            long delta = std::min(e, residuals[edgeId]);
            residuals[edgeId] -= delta;
            residuals[reverse[edgeId]] += delta;
            excess[node] -= delta;
            excess[nextV] += delta;

            if (nextV != t && !isActive[nextV]) {
                activeNodes.push(nextV);
                isActive[nextV] = true;
            }
        } else {
            heights[node] = h + 1;
        }
    } else {
        activeNodes.pop();
        isActive[node] = false;
    }
}

static void fillReverseAndCapacityVectors(Graph &g, std::vector<long>& reverse, std::vector<long>& capacities) {
    long newEdges = 0;

    for (long i = 0; i < g.getNE(); ++i) {
        if (reverse[i] != -1) continue;

        auto srcdst = g.getSrcDstFromId(i);
        long src = srcdst.first;
        long dst = srcdst.second;

        if (!g.hasEdge(dst, src)) {
            reverse[i] = g.getNE() + newEdges;
            reverse[g.getNE() + newEdges] = i;

            capacities[g.getNE() + newEdges] = g.getWeightFromId(i);

            newEdges++;
        } else {
            long id = g.getIdFromSrcDst(dst, src);
            reverse[id] = i;
            reverse[i] = id;

            capacities[i] += g.getWeightFromId(id);
            capacities[id] += g.getWeightFromId(i);
        }
    }
}

void execute(int tid) {
    std::cout << "Thread num: " << tid << std::endl;
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
}

// A parallel lock free max flow algorithm
std::pair<std::vector<long>, std::vector<long>> PLFFlow(Graph &g, const long s, const long t) {
    unsigned int numThreads = std::thread::hardware_concurrency();
    if (numThreads == 0) numThreads = DEFAULT_NUM_THREADS;

    std::vector<std::thread> threads(numThreads);
    for (int i = 0; i < numThreads; ++i) {
        threads[i] = std::thread(execute, i);
    }

    for (int i = 0; i < numThreads; ++i) {
        threads[i].join();
    }

    std::queue<long> activeNodes;
    std::vector<bool> isActive(g.getNV(), false);

    g.createNeighbourList(true);

    data_vec capacities = g.getWeights();
    capacities.resize(g.getUndirectedNumEdges(), 0);
    data_vec reverse(g.getUndirectedNumEdges(), -1);
    fillReverseAndCapacityVectors(g, reverse, capacities);

    data_vec residuals = g.getWeights();
    residuals.resize(g.getUndirectedNumEdges(), 0);

    data_vec heights(g.getNV(), 0);
    data_vec excess(g.getNV(), 0);

    heights[s] = g.getNV();

    // Initial preflow
    std::vector<long> neighbours = g.getNeighbours(s);
    for (auto &n : neighbours) {
        long edgeId = g.getIdFromSrcDst(s, n);
        residuals[edgeId] = 0;
        residuals[reverse[edgeId]] = capacities[edgeId];
        excess[n] = g.getWeightFromId(edgeId);
        excess[s] = excess[s] - g.getWeightFromId(edgeId);

        activeNodes.push(n);
        isActive[n] = true;
    }

    std::cout << activeNodes.size() << std::endl;

    while (excess[s] + excess[t] < 0) {
        long curnode = activeNodes.front();

        process(curnode, g, heights, excess, capacities, reverse, residuals, activeNodes, isActive, s, t);
    }

    return std::make_pair(capacities, residuals);
}