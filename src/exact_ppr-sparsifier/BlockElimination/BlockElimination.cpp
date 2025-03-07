#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <julia.h>

#include "exact_ppr_sparsifier/SC_Block_Elimination.h"
#include "utility/utility.h"

JULIA_DEFINE_FAST_TLS

namespace SC_BlockElimination {

    struct CSCMatrix {
        int n, m;
        std::vector<int> col_ind;
        std::vector<int> row_ind;
        std::vector<double> v;
    };

    /*
    *  [A   B
    *   B^T C]
    */
    struct BlockMatrix {
        GEdge A, B, C;
    };

    // converts a n x m matrix to a sparse csr matrix
    // Warning: indices start at 1 and not 0 as julia uses 1-based-indexing!
    // row indices need to be sorted!
    // (think about introducing to() and from() julia methods)
    CSCMatrix GEdge_to_csc_matrix(GEdge &G) {
        // sort edges by column and row
        auto compareEdges = [](const Edge &a, const Edge &b) {
            if (a.v == b.v) {
                return a.u < b.u;
            }
            return a.v < b.v;
        };
        std::sort(G.edges.begin(), G.edges.end(), compareEdges);

        std::vector<int> col_ind(G.m + 1, 0);
        col_ind[0] = 1;
        std::vector<int> row_ind(G.edges.size());
        std::vector<double> v(G.edges.size());
        int c = 0;
        for (auto e: G.edges) {
            v[c] = e.weight;
            row_ind[c] = e.u + 1;
            ++col_ind[e.v + 1];
            ++c;
        }

        for (int i = 1; i <= G.m; i++) {
            col_ind[i] += col_ind[i - 1];
        }

        return {G.n, G.m, col_ind, row_ind, v};
    }

    // takes an n x n matrix M and splits into four block matrices
    BlockMatrix G_to_block_matrix(const GEdge &G, std::vector<int> &inv_perm, std::vector<double> &weighted_degree, int k) {
        BlockMatrix M;

        M.A.n = M.A.m = G.n - k;
        M.C.n = M.C.m = k;
        M.B.n = G.n - k;
        M.B.m = k;

        // add diagonal entries
        for (int i = 0; i < G.n - k; ++i) {
            M.A.edges.emplace_back(i, i, weighted_degree[i]);
        }
        for (int i = 0; i < k; ++i) {
            M.C.edges.emplace_back(i, i, weighted_degree[i + G.n - k]);
        }

        // add off diagonal entries
        for (auto e: G.edges) {
            e.weight = -e.weight;
            if (inv_perm[e.u] < G.n - k && inv_perm[e.v] < G.n - k) {
                e.u = inv_perm[e.u];
                e.v = inv_perm[e.v];
                M.A.edges.push_back(e);
                M.A.edges.emplace_back(e.v, e.u, e.weight);
            } else if (inv_perm[e.u] >= G.n - k && inv_perm[e.v] >= G.n - k) {
                e.u = inv_perm[e.u] - G.n + k;
                e.v = inv_perm[e.v] - G.n + k;
                M.C.edges.push_back(e);
                M.C.edges.emplace_back(e.v, e.u, e.weight);
            } else {
                if (inv_perm[e.u] >= G.n - k) {
                    e.u = inv_perm[e.u] - G.n + k;
                    e.v = inv_perm[e.v];
                    std::swap(e.u, e.v);
                } else {
                    e.v = inv_perm[e.v] - G.n + k;
                    e.u = inv_perm[e.u];
                }
                M.B.edges.push_back(e);
            }
        }

        return M;
    }

    jl_array_t *matrix_solve(CSCMatrix &M, CSCMatrix &B) {
        // std::cout << "started matrix solve" << std::endl;
        std::string full_path = __FILE__;  // Full path of main.cpp at compile-time
        std::string sourceDir = full_path.substr(0, full_path.find_last_of("/\\"));
        std::string command = "include(\"" + sourceDir + "/system_solver.jl\")";
        // call laplacian solver in julia
        jl_init();
        jl_eval_string(command.c_str());

        // matrix dimensions
        jl_value_t *n1 = jl_box_int32(M.n);
        jl_value_t *m1 = jl_box_int32(M.m);
        jl_value_t *n2 = jl_box_int32(B.n);
        jl_value_t *m2 = jl_box_int32(B.m);

        // lists of indices and values for csr matrices
        jl_value_t *array_type = jl_apply_array_type((jl_value_t *) jl_int32_type, 1);
        jl_array_t *C1 = jl_ptr_to_array_1d(array_type, M.col_ind.data(), M.col_ind.size(), 0);
        jl_array_t *R1 = jl_ptr_to_array_1d(array_type, M.row_ind.data(), M.row_ind.size(), 0);
        jl_array_t *C2 = jl_ptr_to_array_1d(array_type, B.col_ind.data(), B.col_ind.size(), 0);
        jl_array_t *R2 = jl_ptr_to_array_1d(array_type, B.row_ind.data(), B.row_ind.size(), 0);
        array_type = jl_apply_array_type((jl_value_t *) jl_float64_type, 1);
        jl_array_t *V1 = jl_ptr_to_array_1d(array_type, M.v.data(), M.v.size(), 0);
        jl_array_t *V2 = jl_ptr_to_array_1d(array_type, B.v.data(), B.v.size(), 0);

        jl_function_t *func = jl_get_function(jl_main_module, "solve_systems");
        jl_value_t *args[10] = {n1, n2, m1, m2, (jl_value_t *) C1, (jl_value_t *) R1, (jl_value_t *) V1,
                                (jl_value_t *) C2, (jl_value_t *) R2, (jl_value_t *) V2};

        // there should be nicer ways to transfer arguments?
        jl_array_t *x = (jl_array_t *) jl_call(func, args, 10);
        jl_atexit_hook(0);
        // std::cout << "finished matrix solve" << std::endl;
        return x;
    }

    std::vector<std::vector<double>> schur_complement(BlockMatrix &M) {
        // 1. compute X = L_SS_inv * L_SB
        auto L_SS = GEdge_to_csc_matrix(M.A);
        auto L_SB = GEdge_to_csc_matrix(M.B);
        auto X = matrix_solve(L_SS, L_SB);
        // X_T is stored in column major order (default of julia)
        double *X_T = (double *) jl_array_data(X);
        // 2. compute Prod = -L_BS * X
        std::vector<std::vector<double>> SC(L_SB.m, std::vector<double>(L_SB.m, 0));
        for (int i = 0; i < L_SB.m; ++i) {
            for (int j = 0; j < L_SB.m; ++j) {
                for (int k = L_SB.col_ind[i] - 1; k < L_SB.col_ind[i + 1] - 1; ++k) {
                    SC[i][j] -= L_SB.v[k] * X_T[j * L_SB.n + L_SB.row_ind[k] - 1];
                }
            }
        }
        // 3. compute L_BB + X_prod
        for (auto edge: M.C.edges) {
            SC[edge.u][edge.v] += edge.weight;
        }

        return SC;
    }

    DiGraph constructPPRSparsifier(const GEdge &G, std::vector<int>& K, double alpha) {
        std::vector<int> inv_K = getNonTerminals(G.n, K);
        // mapping to divide nodes in inv_K from nodes in K
        std::vector<int> inv_perm(G.n);
        for (std::size_t i = 0; i < inv_K.size(); ++i) {
            inv_perm[inv_K[i]] = i;
        }
        for (int i = 0; i < K.size(); ++i) {
            inv_perm[K[i]] = inv_K.size() + i;
        }

        std::vector<double> weighted_degree(G.n, 0);
        for (auto edge: G.edges) {
            weighted_degree[inv_perm[edge.u]] += edge.weight;
            weighted_degree[inv_perm[edge.v]] += edge.weight;
        }
        for (int i = 0; i < G.n; ++i) {
            weighted_degree[i] += (alpha / (1 - alpha)) * weighted_degree[i];
        }

        BlockMatrix M = G_to_block_matrix(G, inv_perm, weighted_degree, K.size());

        auto S = schur_complement(M);
        std::vector<double> new_weighted_degree(K.size(), 0);
        DiGraph sparsifier(K.size()+1);
        for (int i = 0; i < S.size(); ++i) {
            for (int j = 0; j < S[i].size(); ++j) {
                if (i == j) {
                    sparsifier.add_edge(i, i, weighted_degree[inv_K.size() + i] - S[i][i]);
                } else {
                    sparsifier.add_edge(i, j, -S[i][j]);
                    new_weighted_degree[i] -= S[i][j];
                }
            }
        }

        for (int i = 0; i < K.size(); ++i) {
            sparsifier.add_edge(i, K.size(),
                                S[i][i] - new_weighted_degree[i] - alpha * weighted_degree[inv_K.size() + i]);
        }

        return sparsifier;
    }
}