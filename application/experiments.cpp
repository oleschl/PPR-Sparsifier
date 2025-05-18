#include <set>
#include <iostream>
#include <fstream>
#include "graph.h"
#include "utility/graph_utility.h"
#include "utility/pagerank.h"
#include "utility/utility.h"
#include "exact_ppr_sparsifier/SC_Adj_List.h"
#include "approximate_ppr_sparsifier/approximate_pp_sparsifier.h"

// code for creating UNWEIGHTED version of sparsifiers
std::vector<std::pair<int, int>> getUnweightedEdges(DiGraph& G, int num_edges) {
    // our sparsifiers are not guaranteed to return a normalized graph
    // to be consistent with [VCG11] we normalize the outdegree here
    std::vector<double> weighted_degree(G.n);
    for(int i = 0; i < G.n-1; ++i){
        for(auto edge : G.adjList[i]){
            weighted_degree[i] += edge.second;
        }
    }

    std::vector<Edge> edges;

    double count = 0;
    for(int i = 0; i < G.n-1; ++i){
        double sum = 0;
        for(int j = 0; j < G.adjList[i].size(); ++j){
            if(G.adjList[i][j].first != G.n -1 && i != G.adjList[i][j].first){
                edges.emplace_back(i, G.adjList[i][j].first, G.adjList[i][j].second/weighted_degree[i]);
            }
        }
    }

    std::sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) {
        if (a.weight != b.weight) {
            return a.weight > b.weight;
        }
        if (a.u != b.u) {
            return a.u < b.u;
        }
        return a.v < b.v;
    });

    std::vector<std::pair<int, int>> unweighted_edges(num_edges);
    for(int i = 0; i < num_edges; ++ i) {
        unweighted_edges[i] = {edges[i].u, edges[i].v};
    }

    return unweighted_edges;
}

// code used for comparing UNWEIGHTED to our approximate sparisfier by computing edge overlap using Jaccard Index
int experiment_unweighted(GEdge& G, DiGraph& exactSparsifier, DiGraph& approximateSparsifier, std::vector<int>& K) {
    auto induced = inducedSubgraph(G, K);
    std::cout << "undirected edges in induced: " << induced.m << std::endl;
    auto edges_exact = getUnweightedEdges(exactSparsifier, 2*induced.m);
    auto edges_approx = getUnweightedEdges(approximateSparsifier, 2*induced.m);

    std::vector<std::pair<int, int>> induced_edges;
    for(auto edge : induced.edges) {
        induced_edges.emplace_back(edge.u, edge.v);
        induced_edges.emplace_back(edge.v, edge.u);
    }
    auto set_induced = std::set(induced_edges.begin(), induced_edges.end());
    auto set_exact = std::set(edges_exact.begin(), edges_exact.end());
    auto set_approx = std::set(edges_approx.begin(), edges_approx.end());

    std::vector<std::pair<int, int>> intersection;
    std::set_intersection(set_exact.begin(), set_exact.end(),
                          set_approx.begin(), set_approx.end(),
                          std::back_inserter(intersection));

    std::vector<std::pair<int, int>> uni;
    std::set_union(set_exact.begin(), set_exact.end(),
                   set_approx.begin(), set_approx.end(),
                   std::back_inserter(uni));

    std::vector<std::pair<int, int>> intersection2;
    std::set_intersection(set_exact.begin(), set_exact.end(),
                          set_induced.begin(), set_induced.end(),
                          std::back_inserter(intersection2));

    std::vector<std::pair<int, int>> uni2;
    std::set_union(set_exact.begin(), set_exact.end(),
                   set_induced.begin(), set_induced.end(),
                   std::back_inserter(uni2));

    if (uni.empty()) return 0.0;  // avoid division by zero

    std::cout << "Jaccard Index Exact and Approx: " << static_cast<double>(intersection.size()) / uni.size() << std::endl;
    std::cout << "Jaccard Index Exact and Induced: " << static_cast<double>(intersection2.size()) / uni2.size() << std::endl;

    return static_cast<double>(intersection.size()) / uni.size();
}

// code used for comparing WEIGHTED to our approximate sparsifier
// returns the error norm between WEIGHTED subgraph PPR values and PPR values of input graph G
double experiment_weighted(std::vector<std::vector<double>> PPR_G, DiGraph& exactSparsifier, DiGraph& approximateSparsifier, std::vector<int>& K, double alpha) {
    std::vector<Edge> edges;
    for(int i = 0; i < exactSparsifier.n; ++i){
        for(auto j : exactSparsifier.adjList[i]){
            edges.emplace_back(i, j.first, j.second);
        }
    }

    std::sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) {
        if (a.weight != b.weight) {
            return a.weight > b.weight;  // Primary: descending weight
        }
        // Secondary: sort by (u, v) ascending
        if (a.u != b.u) {
            return a.u < b.u;
        }
        return a.v < b.v;
    });

    std::vector<Edge> new_edges(approximateSparsifier.m);
    for(int i = 0; i < approximateSparsifier.m; ++ i) {
        new_edges[i] = edges[i];
    }
    DiGraph new_spar(exactSparsifier.n);
    for(auto edge : new_edges){
        new_spar.add_edge(edge.u, edge.v, edge.weight);
    }

    Sparsifier sp = {new_spar, K};

    return comparePPVs(PPR_G, sp, 1-alpha);
}

// computes the error between the input and induced subgraph
double norm_induced(GEdge& G, std::vector<int>& K, double alpha){
    auto induced = inducedSubgraph(G, K);
    return frobeniusInduced(G, induced, K, 1-alpha);
}

void writeMatrixToFile(const std::vector<std::vector<double>>& matrix, const std::string& filename) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Failed to open file: " << filename << "\n";
        return;
    }

    for (const auto& row : matrix) {
        for (size_t j = 0; j < row.size(); ++j) {
            out << row[j];
            if (j < row.size() - 1) out << " ";  // space-separated
        }
        out << "\n";
    }

    out.close();
}

void writeMatrixToFile2(const std::vector<std::vector<int>>& matrix, const std::string& filename) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Failed to open file: " << filename << "\n";
        return;
    }

    for (const auto& row : matrix) {
        for (size_t j = 0; j < row.size(); ++j) {
            out << row[j];
            if (j < row.size() - 1) out << " ";  // space-separated
        }
        out << "\n";
    }

    out.close();
}

// code for creating PPR matrices for clustering experiment in python notebook
void create_PPR_matrices(GEdge& G, double alpha, std::string datasetName) {
    int seed = 42;
    int k = 1000;
    int split = 10;
    int merge = 10;
    auto K = getRandomTerminals(G.n, k , seed);
    std::string fileName = datasetName + "_exact_RN";
    // std::string fileName = datasetName + "_spiel_10_10_RN.txt";

    std::string version = "dynamic";
    auto exact = SC_Adj_List::constructPPRSparsifier(G, K, alpha, version);
    auto PPR = getPPRMatrixSparsifier(exact, 1-alpha);
    writeMatrixToFile(PPR, fileName);

    fileName = datasetName + "_spiel_10_10_RN.txt";
    auto spiel_10_10_RN = ApproximateSparsifier::constructPPRSparsifier(G, K, alpha, split, merge, seed, "elim_star");
    PPR = getPPRMatrixSparsifier(spiel_10_10_RN, 1-alpha);
    writeMatrixToFile(PPR, fileName);

    auto kyng_10_10_RN = ApproximateSparsifier::constructPPRSparsifier(G, K, alpha, split, merge, seed, "random_clique");
    fileName = datasetName + "_kyng_10_10_RN.txt";
    PPR = getPPRMatrixSparsifier(kyng_10_10_RN, 1-alpha);
    writeMatrixToFile(PPR, fileName);

    auto induced_RN = inducedSubgraph(G, K);
    fileName = datasetName + "_induced_RN.txt";
    PPR = getPPRMatrixG(induced_RN, 1-alpha);
    writeMatrixToFile(PPR, fileName);

    fileName = datasetName + "_map.txt";
    writeMatrixToFile2({K}, fileName);
}
