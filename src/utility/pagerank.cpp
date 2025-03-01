#include <vector>
#include <numeric>
#include <cstdlib>
#include "graph.h"

std::vector<std::vector<std::pair<int, double>>> constructWalkingMatrix(const DiGraph& G, double alpha) {
    int n = G.n;
    std::vector<std::vector<std::pair<int, double>>> W(G.adjList);
    std::vector<double> sums(n);

    for(int i = 0; i < n; ++i) {
        for(const auto& j : G.adjList[i]) {
            sums[i] += j.second;
        }
    }

    // scale such that W is walking matrix
    for(int i = 0; i < n; ++i) {
        for(auto& j : W[i]) {
            j.second *= (alpha/sums[i]);
        }
    }

    return W;
}

std::vector<double> pageRank(const std::vector<std::vector<std::pair<int, double>>>& W, std::vector<double>& r, int n, double eps, int matIter) {
    std::vector<double> p(r);

    for(int i = 0; i < matIter; ++i) {
        std::vector<double> newP(n, 0);
        for(int j = 0; j < n; ++j) {
            for(const auto& k : W[j]) {
                newP[k.first] += k.second * p[j];
            }
        }

        double sum = std::accumulate(newP.begin(), newP.end(), 0.0);
        for(int j = 0; j < n; ++j) {
            newP[j] += (1-sum) * r[j];
        }

        double diff = 0;
        for(int j = 0; j < n; ++j) {
            diff += std::abs(newP[j] - p[j]);
        }
        if (diff < eps) return newP;

        std::swap(newP, p);
    }

    std::cout << "stopped because max iteration" << std::endl;

    return p;
}

// compute (personalized) pagerank using power iteration
std::vector<double> pageRank(DiGraph& G, std::vector<double> r, double alpha, double eps, int matIter) {
    auto W = constructWalkingMatrix(G, alpha);
    return pageRank(W, r, G.n, eps, matIter);
}

double comparePPVs(GEdge& G, Sparsifier& H, double alpha) {
    std::vector<double> pH(H.H.n-1);
    std::vector<double> pG(G.n);

    DiGraph dirG(G.n);
    for(auto edge : G.edges){
        dirG.add_edge(edge.u, edge.v, edge.weight);
        dirG.add_edge(edge.v, edge.u, edge.weight);
    }

    double diff = 0;

    for(int i = 0; i < H.H.n-1; ++i) {
        std::cout << "PPV: " << i << std::endl;
        // compute pagerank of g and h for p_i
        // need correct mapping for p_i for G
        std::vector<double> rH(H.H.n, 0);
        rH[i] = 1;
        pH = pageRank(H.H, rH, alpha, 1e-12, 10000);

        std::vector<double> rG(G.n, 0);
        rG[H.mapHToG[i]] = 1;
        pG = pageRank(dirG, rG, alpha, 1e-12, 10000);
        // compare common p values and normalize
        double sumG = 0.0;
        double sumH = 0.0;
        for(int j = 0; j < H.H.n-1; ++j) {
            sumG += pG[H.mapHToG[j]];
            sumH += pH[j];
        }
        for(int j = 0; j < H.H.n-1; ++j) {
            auto localDiff = std::abs((pH[j]/sumH)-(pG[H.mapHToG[j]]/sumG));
            auto localDiff2 = (pH[j]/sumH)/(pG[H.mapHToG[j]]/sumG);
            std::cout << "diff " << pH[j]/sumH << " and " << pG[H.mapHToG[j]]/sumG << std::endl;
            if(localDiff > 0.000002) {
                std::cout << "high difference detected for node " << j << ":" << pH[j]/sumH << " and " << pG[H.mapHToG[j]]/sumG << std::endl;
            }
            diff += localDiff;
        }
    }

    return diff;
}