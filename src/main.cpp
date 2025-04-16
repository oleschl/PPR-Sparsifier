#include <iostream>
#include <chrono>
#include <functional>
#include <unistd.h>

#include <julia.h>
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
    std::vector<double> runtimes;
    std::vector<double> norms;
    auto K = getRandomTerminals(G.n, k, 42);
    for (int i = 0; i < runs; ++i) {
        auto result = measureTime(ApproximateSparsifier::constructPPRSparsifier, G, K, alpha, split, merge, 42+i, version);
        std::cout << "runtime: " << result.first << " ms\n";
        Sparsifier sp = {result.second, K};
        auto norm = comparePPVs(G, sp, 1 - alpha);
        runtimes.push_back(result.first);
        norms.push_back(norm);
    }
    std::cout << "mean: " << mean(runtimes) << "sample variance:  " << variance_sample(runtimes) << std::endl;
    std::cout << "mean: " << mean(norms) << "sample variance:  " << variance_sample(norms) << std::endl;
}

void run_approximate_chimera(const std::string& version, bool debug, int k, float alpha, bool o, int n, int runs, int split, int merge) {
    for (int i = 0; i < runs; ++i) {
        GEdge G = chimera(n, i+1, false);
        std::cout << "num nodes: " << G.n << " , num edges: " << G.m << std::endl;
        auto K = getRandomTerminals(G.n, k, 42+i);
        auto result = measureTime(ApproximateSparsifier::constructPPRSparsifier, G, K, alpha, split, merge, 42+i, version);
        std::cout << "runtime: " << result.first << " ms\n";
        // TODO code for measuring error (frobenius norm, ...)
        Sparsifier sp = {result.second, K};
        comparePPVs(G, sp, 1 - alpha);
    }
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

int main(int argc, char *argv[]) {
    // init julia here
    jl_init();
    jl_eval_string("using Laplacians");

    std::string mode;
    std::string algorithm;
    std::string version;
    std::string filename;
    bool debug = false;
    bool o_flag = false;
    int k = 1;
    float alpha = 0.85;
    int split = 1, merge = 1, runs = 1;
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
                std::cerr << "Unknown or malformed option: " << opt << "\n";
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
