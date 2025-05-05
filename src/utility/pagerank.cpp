#include <vector>
#include <numeric>
#include <cstdlib>
#include <cmath>
#include <iostream>
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
            if(sums[i] == 0) {
                j.second == 0;
            } else if (j.second < 0) {
                if (std::abs(j.second) < 1e-7) {
                    j.second *= -(alpha/sums[i]);
                } else {
                    std::cout << "error" << std::endl;
                }
            } else {
                j.second *= (alpha/sums[i]);
            }
            //j.second *= (alpha/sums[i]);
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

double frobeniusNorm(std::vector<std::vector<double>>& PPR1, std::vector<std::vector<double>>& PPR2, int n) {
    double sum1 = 0.0;
    double sum2 = 0.0;
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            sum1 += PPR1[i][j] * PPR1[i][j];
            sum2 += (PPR1[i][j] - PPR2[i][j]) * (PPR1[i][j] - PPR2[i][j]);
        }
    }

    return std::sqrt(sum2/sum1);
}

std::vector<std::vector<double>> precomputePPRMatrix(GEdge& G, const std::vector<int>& K, double alpha) {
    int k = K.size();
    std::vector<std::vector<double>> PPR_G(k, std::vector<double>(k));

    DiGraph dirG(G.n);
    for (const auto& edge : G.edges) {
        dirG.add_edge(edge.u, edge.v, edge.weight);
        dirG.add_edge(edge.v, edge.u, edge.weight);
    }

    // Compute personalized PageRank for each terminal
    for (int i = 0; i < k; ++i) {
        std::vector<double> personalization(G.n, 0.0);
        personalization[K[i]] = 1.0;

        auto ppr = pageRank(dirG, personalization, alpha, 1e-12, 10000);

        // Normalize over terminals
        double norm_factor = 0.0;
        for (int node : K) {
            norm_factor += ppr[node];
        }

        for (int j = 0; j < k; ++j) {
            PPR_G[i][j] = ppr[K[j]] / norm_factor;
        }
    }

    return PPR_G;
}

double frobeniusInduced(GEdge& G, GEdge& induced, std::vector<int> mapHtoG, double alpha) {
    DiGraph dirG(G.n);
    for (const auto& edge : G.edges) {
        dirG.add_edge(edge.u, edge.v, edge.weight);
        dirG.add_edge(edge.v, edge.u, edge.weight);
    }
    std::vector<std::vector<double>> PPR(induced.n, std::vector<double> (induced.n));

    for(int i = 0; i < induced.n; ++i) {
        std::vector<double> rG(G.n, 0);
        rG[mapHtoG[i]] = 1;
        auto pG = pageRank(dirG, rG, alpha, 1e-12, 10000);
        // compare common p values and normalize
        double sumG = 0.0;
        for(int j = 0; j < induced.n; ++j) {
            sumG += pG[mapHtoG[j]];
        }
        for(int j = 0; j < induced.n; ++j) {
            PPR[i][j] = pG[mapHtoG[j]]/sumG;
        }
    }

    DiGraph dirInd(induced.n);
    for (const auto& edge : induced.edges) {
        dirInd.add_edge(edge.u, edge.v, edge.weight);
        dirInd.add_edge(edge.v, edge.u, edge.weight);
    }
    std::vector<std::vector<double>> PPR_ind(induced.n, std::vector<double> (induced.n));
    for(int i = 0; i < induced.n; ++i) {
        std::vector<double> rG(induced.n, 0);
        rG[i] = 1;
        PPR_ind[i] = pageRank(dirInd, rG, alpha, 1e-12, 10000);
    }

    return frobeniusNorm(PPR, PPR_ind, induced.n);
}


double comparePPVs(GEdge& G, Sparsifier& H, double alpha) {

    std::vector<double> pH(H.H.n-1);
    std::vector<double> pG(G.n);

    std::vector<std::vector<double>> PPR_G(H.H.n-1, std::vector<double> (H.H.n-1));
    std::vector<std::vector<double>> PPR_H(H.H.n-1, std::vector<double> (H.H.n-1));

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
            PPR_G[i][j] = pG[H.mapHToG[j]]/sumG;
            PPR_H[i][j] = pH[j]/sumH;
        }

        for(int j = 0; j < H.H.n-1; ++j) {
            auto localDiff = std::abs((pH[j]/sumH)-(pG[H.mapHToG[j]]/sumG));
            auto localDiff2 = (pH[j]/sumH)/(pG[H.mapHToG[j]]/sumG);
            //std::cout << "diff " << pH[j]/sumH << " and " << pG[H.mapHToG[j]]/sumG << std::endl;
            if(localDiff > 0.1) {
                std::cout << "high difference detected for node " << j << ":" << pH[j]/sumH << " and " << pG[H.mapHToG[j]]/sumG << std::endl;
            }
            diff += localDiff;
        }
    }

    auto norm = frobeniusNorm(PPR_G, PPR_H, H.H.n-1);

    return norm;
}

double comparePPVs(std::vector<std::vector<double>>& PPR_G, Sparsifier& H, double alpha) {

    std::vector<double> pH(H.H.n-1);
    std::vector<std::vector<double>> PPR_H(H.H.n-1, std::vector<double> (H.H.n-1));

    double diff = 0;

    for(int i = 0; i < H.H.n-1; ++i) {
        // std::cout << "PPV: " << i << std::endl;
        // compute pagerank of g and h for p_i
        // need correct mapping for p_i for G
        std::vector<double> rH(H.H.n, 0);
        rH[i] = 1;
        pH = pageRank(H.H, rH, alpha, 1e-12, 10000);

        // compare common p values and normalize
        double sumG = 0.0;
        double sumH = 0.0;
        for(int j = 0; j < H.H.n-1; ++j) {
            sumH += pH[j];
        }
        for(int j = 0; j < H.H.n-1; ++j) {
            PPR_H[i][j] = pH[j]/sumH;
        }

        for(int j = 0; j < H.H.n-1; ++j) {
            auto localDiff = std::abs(PPR_H[i][j]-PPR_G[i][j]);
            // auto localDiff = std::abs((pH[j]/sumH)-(pG[H.mapHToG[j]]/sumG));
            // auto localDiff2 = (pH[j]/sumH)/(pG[H.mapHToG[j]]/sumG);
            //std::cout << "diff " << pH[j]/sumH << " and " << pG[H.mapHToG[j]]/sumG << std::endl;
            if(localDiff > 0.01) {
                //std::cout << "high difference detected for node " << j << ":" << PPR_H[i][j] << " and " << PPR_G[i][j] << std::endl;
            }
            diff += localDiff;
        }
    }

    auto norm = frobeniusNorm(PPR_G, PPR_H, H.H.n-1);

    return norm;
}