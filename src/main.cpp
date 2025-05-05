#include <iostream>
#include <chrono>
#include <functional>
#include <unistd.h>

#include <julia.h>
#include <set>
#include "graph.h"
#include "./utility/utility.h"
#include "./utility/graph_utility.h"
#include "./utility/pagerank.h"
#include "exact_ppr_sparsifier/NodeRemoval.h"
#include "exact_ppr_sparsifier/PPR_Matrix_Inv.h"
#include "exact_ppr_sparsifier/SC_Block_Elimination.h"
#include "exact_ppr_sparsifier/SC_Adj_List.h"
#include "exact_ppr_sparsifier/SC_Hash_Table.h"
#include "approximate_ppr_sparsifier/approximate_pp_sparsifier.h"

template <typename Func, typename... Args>
auto measureTime(Func func, Args&&... args) {
    auto start = std::chrono::high_resolution_clock::now();
    auto result = func(std::forward<Args>(args)...);
    auto end = std::chrono::high_resolution_clock::now();

    double elapsed_time = std::chrono::duration<double, std::milli>(end - start).count(); // Time in milliseconds
    return std::make_pair(elapsed_time, result);
}

void run_exact_from_file(const std::string& algorithm, const std::string& version, bool debug, int k, double alpha, bool o, const std::string& filename) {
    GEdge G;
    parseEdgeList(filename, G);
    std::cout << "num nodes: " << G.n << " , num edges: " << G.m << std::endl;
    auto K = getRandomTerminals(G.n, k, 42);
    // auto K = get_RNE_terminals(G, k, 42);

    std::pair<double, DiGraph> result;
    if (algorithm == "node_removal") {
        if (o) {
            result = measureTime(NodeRemoval::constructPPRSparsifier2, G, K, alpha);
        } else {
            result = measureTime(NodeRemoval::constructPPRSparsifier, G, K, alpha, version);
        }
    } else if (algorithm == "sc_hash_table") {
        if (o) {
            result = measureTime(SC_Hash_Table::constructPPRSparsifier2, G, K, alpha);
        } else {
            result = measureTime(SC_Hash_Table::constructPPRSparsifier, G, K, alpha, version);
        }
    } else if (algorithm == "sc_adj_list") {
        result = measureTime(SC_Adj_List::constructPPRSparsifier, G, K, alpha, version);
    } else if (algorithm == "block_elimination") {
        result = measureTime(SC_BlockElimination::constructPPRSparsifier, G, K, alpha);
    } else if (algorithm == "ppr_matrix_inv") {
        result = measureTime(PPR_Matrix_Inv::constructPPRSparsifier, G, K, alpha);
    } else {
        std::cerr << "Unknown algorithm: " << algorithm << std::endl;
    }

    std::cout << "runtime: " << result.first << " ms\n";

    Sparsifier sp = {result.second, K};
    if (debug) comparePPVs(G, sp, 1 - alpha);
}

void run_exact_chimera(const std::string& algorithm, const std::string& version, bool debug, int k, float alpha, bool o, int n, int runs) {
    std::vector<double> runtimes;
    for(int i = 0; i < runs; ++i) {
        GEdge G = chimera(n, i+1, false);
        std::cout << "num nodes: " << G.n << " , num edges: " << G.m << std::endl;
        auto K = getRandomTerminals(G.n, k, 42+i);

        std::pair<double, DiGraph> result;
        if (algorithm == "node_removal") {
            if (o) {
                result = measureTime(NodeRemoval::constructPPRSparsifier2, G, K, alpha);
            } else {
                result = measureTime(NodeRemoval::constructPPRSparsifier, G, K, alpha, version);
            }
        } else if (algorithm == "sc_hash_table") {
            if (o) {
                result = measureTime(SC_Hash_Table::constructPPRSparsifier2, G, K, alpha);
            } else {
                result = measureTime(SC_Hash_Table::constructPPRSparsifier, G, K, alpha, version);
            }
        } else if (algorithm == "sc_adj_list") {
            result = measureTime(SC_Adj_List::constructPPRSparsifier, G, K, alpha, version);
        } else if (algorithm == "block_elimination") {
            result = measureTime(SC_BlockElimination::constructPPRSparsifier, G, K, alpha);
        } else if (algorithm == "ppr_matrix_inv") {
            result = measureTime(PPR_Matrix_Inv::constructPPRSparsifier, G, K, alpha);
        } else if (algorithm == "approximate") {
            int split = 10;
            int merge = 10;
            result = measureTime(ApproximateSparsifier::constructPPRSparsifier, G, K, alpha, split, merge, 42, version);
        } else {
            std::cerr << "Unknown algorithm: " << algorithm << std::endl;
        }

        runtimes.push_back(result.first);
        std::cout << "runtime: " << result.first << " ms\n";

        Sparsifier sp = {result.second, K};
        if (debug) comparePPVs(G, sp, 1 - alpha);
    }

    std::cout << "median: " << compute_percentile(runtimes, 0.5) << std::endl;
    std::cout << "75th percentile: " << compute_percentile(runtimes, 0.75) << std::endl;
    std::cout << "worst case: " << *std::max_element(runtimes.begin(), runtimes.end()) << std::endl;
}

void run_approximate_from_file(const std::string& version, bool debug, int k, float alpha, bool o, const std::string& filename, int split, int merge, int runs) {
    GEdge G;
    parseEdgeList(filename, G);
    std::cout << "num nodes: " << G.n << " , num edges: " << G.m << std::endl;
    auto K = getRandomTerminals(G.n, k, 42);
    // auto K = get_RNE_terminals(G, k, 42);
    std::vector<std::vector<double>> PPR_G;
    if(debug) {
        PPR_G = precomputePPRMatrix(G, K, 1-alpha);
    }
    std::vector<int> params = {2, 3, 5, 10};
    for(auto p : params) {
        std::cout << "param value: " << p << std::endl;
        std::vector<double> runtimes;
        std::vector<double> norms;
        std::vector<double> num_edges;
        for (int i = 0; i < runs; ++i) {
            auto result = measureTime(ApproximateSparsifier::constructPPRSparsifier, G, K, alpha, p, p, 42 + i,
                                      version);
            std::cout << "runtime: " << result.first << " ms\n";
            std::cout << "edges: " << result.second.m << "\n";
            runtimes.push_back(result.first);
            if (debug) {
                Sparsifier sp = {result.second, K};
                auto norm = comparePPVs(PPR_G, sp, 1 - alpha);
                norms.push_back(norm);
                num_edges.push_back(result.second.m);
            }
        }
        std::cout << "mean runtimes: " << mean(runtimes) << "sample variance:  " << sample_stddev(runtimes) << std::endl;
        std::cout << "mean norms: " << mean(norms) << "sample variance:  " << sample_stddev(norms) << std::endl;
        std::cout << "mean edges: " << mean(num_edges) << "sample variance:  " << sample_stddev(num_edges) << std::endl;
    }

    std::string v = "random_clique";
    for(auto p : params) {
        std::cout << "param value: " << p << std::endl;
        std::vector<double> runtimes;
        std::vector<double> norms;
        std::vector<double> num_edges;
        for (int i = 0; i < runs; ++i) {
            auto result = measureTime(ApproximateSparsifier::constructPPRSparsifier, G, K, alpha, p, p, 42 + i,
                                      v);
            std::cout << "runtime: " << result.first << " ms\n";
            std::cout << "edges: " << result.second.m << "\n";
            runtimes.push_back(result.first);
            if (debug) {
                Sparsifier sp = {result.second, K};
                auto norm = comparePPVs(PPR_G, sp, 1 - alpha);
                norms.push_back(norm);
                num_edges.push_back(result.second.m);
            }
        }
        std::cout << "mean runtimes: " << mean(runtimes) << "sample variance:  " << sample_stddev(runtimes) << std::endl;
        std::cout << "mean norms: " << mean(norms) << "sample variance:  " << sample_stddev(norms) << std::endl;
        std::cout << "mean edges: " << mean(num_edges) << "sample variance:  " << sample_stddev(num_edges) << std::endl;
    }
}

void run_approximate_chimera(const std::string& version, bool debug, int k, float alpha, bool o, int n, int runs, int split, int merge) {
    std::vector<double> runtimes;
    std::vector<double> norms;
    std::vector<double> avg_edges;
    for (int i = 0; i < runs; ++i) {
        GEdge G = chimera(n, i+1, false);
        std::cout << "num nodes: " << G.n << " , num edges: " << G.m << std::endl;
        auto K = getRandomTerminals(G.n, k, 42);
        std::vector<std::vector<double>> PPR_G;
        if(debug) {
            PPR_G = precomputePPRMatrix(G, K, 1-alpha);
        }
        std::vector<double> tmp_runtimes;
        std::vector<double> tmp_norms;
        std::vector<double> num_edges;
        for (int j = 0; j < runs; ++j) {
            auto result = measureTime(ApproximateSparsifier::constructPPRSparsifier, G, K, alpha, split, merge, 42 + j,
                                      version);
            std::cout << "runtime: " << result.first << " ms\n";
            tmp_runtimes.push_back(result.first);
            if(debug){
                Sparsifier sp = {result.second, K};
                auto norm = comparePPVs(PPR_G, sp, 1 - alpha);
                norms.push_back(norm);
                tmp_norms.push_back(norm);
                norms.push_back(norm);
                num_edges.push_back(result.second.m);
                std::cout << "edges: " << result.second.m << std::endl;
            }
        }
        auto mean_runtimes = mean(tmp_runtimes);
        auto variance_runtimes = sample_stddev(tmp_runtimes);
        auto mean_norms = mean(tmp_norms);
        auto variance_norms = sample_stddev(tmp_norms);
        auto mean_edge = mean(num_edges);
        auto variance_edge = sample_stddev(num_edges);
        runtimes.push_back(mean_runtimes);
        std::cout << "runtime mean: " << mean_runtimes << "sample variance:  " << variance_runtimes << std::endl;
        std::cout << "norm mean: " << mean_norms << "sample variance:  " << variance_norms << std::endl;
        std::cout << "edges mean: " << mean_edge << "sample variance:  " << variance_edge << std::endl;
    }
    std::cout << "median: " << compute_percentile(runtimes, 0.5) << std::endl;
    std::cout << "75th percentile: " << compute_percentile(runtimes, 0.75) << std::endl;
    std::cout << "worst case: " << *std::max_element(runtimes.begin(), runtimes.end()) << std::endl;

    std::cout << "mean norm: " << mean(norms) << std::endl;
    std::cout << "variance norm: " << sample_stddev(norms) << std::endl;
}

void run_approximate_sachdeva_star(const std::string& version, bool debug, int k, float alpha, bool o, int runs,  int j, int l, int split, int merge) {
    GEdge G = sachdevaStar(j, l);
    std::cout << "num nodes: " << G.n << " , num edges: " << G.m << std::endl;
    auto K = getRandomTerminals(G.n, k, 42);
    for (int i = 0; i < runs; ++i) {
        auto result = measureTime(ApproximateSparsifier::constructPPRSparsifier, G, K, alpha, split, merge, 42+i, version);
        std::cout << "runtime: " << result.first << " ms\n";
        // TODO code for measuring error (frobenius norm, ...)
        Sparsifier sp = {result.second, K};
        comparePPVs(G, sp, 1 - alpha);
    }
}

std::vector<std::pair<int, int>> getUnweightedEdges(DiGraph& G, int num_edges) {
    std::vector<double> weighted_degree(G.n);
    for(int i = 0; i < G.n-1; ++i){
        for(auto edge : G.adjList[i]){
            weighted_degree[i] += edge.second;
        }
    }

    std::vector<Edge> edges;
    // TODO: digraphs should be that all edges except to sink have equal weight and so are undirected
    // true for ours but pagerank sparsifier not unique but because of normalization can have scaling
    // maybe better to normalize so that unique
    for(int i = 0; i < G.n-1; ++i){
        for(int j = 0; j < G.adjList[i].size(); ++j){
            if(G.adjList[i][j].first != G.n -1 && i != G.adjList[i][j].first)
                edges.emplace_back(i, G.adjList[i][j].first, G.adjList[i][j].second/weighted_degree[i]);
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

    std::vector<std::pair<int, int>> unweighted_edges(num_edges);
    for(int i = 0; i < num_edges; ++ i) {
        unweighted_edges[i] = {edges[i].u, edges[i].v};
    }

    return unweighted_edges;
}

double experiment_weighted(GEdge& G, DiGraph& exactSparsifier, DiGraph& approximateSparsifier, std::vector<int>& K, double alpha) {
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

    return comparePPVs(G, sp, 1 - alpha);
}

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

    std::cout << "Jaccard: " << static_cast<double>(intersection.size()) / uni.size() << std::endl;
    std::cout << "Jaccard ind: " << static_cast<double>(intersection2.size()) / uni2.size() << std::endl;

    return static_cast<double>(intersection.size()) / uni.size();
}

double norm_induced(GEdge& G, std::vector<int>& K, double alpha){
    auto induced = inducedSubgraph(G, K);
    return frobeniusInduced(G, induced, K, 1-alpha);
}


int main(int argc, char *argv[]) {
    // init julia here
    // jl_init();
    //jl_eval_string("using Laplacians");

    std::string mode;
    std::string algorithm;
    std::string version;
    std::string filename;
    bool debug = false;
    bool o_flag = false;
    int k = 1;
    float alpha = 0.85;
    int split = 2, merge = 2, runs = 1;
    int n = 1000, j = 10, l = 10;

    int opt;
    while ((opt = getopt(argc, argv, "m:a:v:dk:x:of:s:g:r:n:j:l:")) != -1) {
        switch (opt) {
            case 'm': mode = optarg; break;
            case 'a': algorithm = optarg; break;
            case 'v': version = optarg; break;
            case 'd': debug = true; break;
            case 'k': k = std::stoi(optarg); break;
            case 'x': alpha = std::stof(optarg); break;
            case 'o': o_flag = true; break;
            case 'f': filename = optarg; break;
            case 's': split = std::stoi(optarg); break;
            case 'g': merge = std::stoi(optarg); break;
            case 'r': runs = std::stoi(optarg); break;
            case 'n': n = std::stoi(optarg); break;
            case 'j': j = std::stoi(optarg); break;
            case 'l': l = std::stoi(optarg); break;
            default:
                std::cerr << "Unknown option: " << opt << "\n";
                return 1;
        }
    }

    if (mode == "exact_from_file") {
        if (filename.empty()) {
            std::cerr << "Error: -f (filename) is required for exact_from_file\n";
            return 1;
        }
        run_exact_from_file(algorithm, version, debug, k, alpha, o_flag, filename);
    } else if (mode == "exact_chimera") {
        run_exact_chimera(algorithm, version, debug, k, alpha, o_flag, n, runs);
    } else if (mode == "approximate_from_file") {
        if (filename.empty()) {
            std::cerr << "Error: -f (filename) is required for approximate_from_file\n";
            return 1;
        }
        run_approximate_from_file(version, debug, k, alpha, o_flag, filename, split, merge, runs);
    } else if (mode == "approximate_chimera") {
        run_approximate_chimera(version, debug, k, alpha, o_flag, n, runs, split, merge);
    } else if (mode == "approximate_sachdeva_star") {
        run_approximate_sachdeva_star(version, debug, k, alpha, o_flag, runs, j, l, split, merge);
    } else {
        std::cerr << "Error: unknown mode \"" << mode << "\"\n";
        return 1;
    }

    jl_atexit_hook(0);
    return 0;
}
