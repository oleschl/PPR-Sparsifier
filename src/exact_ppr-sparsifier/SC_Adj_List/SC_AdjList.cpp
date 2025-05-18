#include "graph.h"
#include "utility/utility.h"
#include "utility/orderings.h"
#include "symbolic_factorization.cpp"
#include "exact_ppr_sparsifier/SC_Adj_List.h"

namespace SC_Adj_List {

    DiGraph constructPPRSparsifier(const GEdge &G, std::vector<int> &K, double alpha, const std::string &order) {
        // get set of non-terminal nodes
        auto inv_K = getNonTerminals(G.n, K);
        std::vector<double> weighted_degree(G.n, 0);
        for (auto edge: G.edges) {
            weighted_degree[edge.u] += edge.weight;
            weighted_degree[edge.v] += edge.weight;
        }

        for (int i = 0; i < G.n; ++i) {
            weighted_degree[i] += (alpha / (1 - alpha)) * weighted_degree[i];
        }
        std::vector<double> initial_diag(weighted_degree);

        // create compressed graph structure for computing minimum degree ordering
        std::vector<int> xadj(G.n + 1);
        std::vector<int> adj(2 * G.m);
        std::vector<double> weights(2 * G.m);
        for (auto edge: G.edges) {
            ++xadj[edge.u + 1];
            ++xadj[edge.v + 1];
        }
        for (int i = 1; i <= G.n; ++i) {
            xadj[i] += xadj[i - 1];
        }
        std::vector<int> temp_pos(G.n, 0);
        for (auto edge: G.edges) {
            adj[xadj[edge.u] + temp_pos[edge.u]] = edge.v;
            weights[xadj[edge.u] + temp_pos[edge.u]] = -edge.weight;
            ++temp_pos[edge.u];
            adj[xadj[edge.v] + temp_pos[edge.v]] = edge.u;
            weights[xadj[edge.v] + temp_pos[edge.v]] = -edge.weight;
            ++temp_pos[edge.v];
        }

        // get ordering
        std::pair<std::vector<int>, std::vector<int>> ordering;
        if (order == "random") {
            ordering = getRandomOrdering(G, K, inv_K);
        } else if (order == "static") {
            ordering = getStaticMinDegOrdering(G, K, inv_K);
        } else if (order == "dynamic") {
            ordering = getDynamicMinDegOrdering(G.n, G.m, xadj, adj, K, inv_K);
        }
        auto [perm, inv_perm] = ordering;

        auto M = symbolic_factorization(G.n, inv_K.size(), adj, xadj, perm, inv_perm);

        std::vector<int> link(G.n, -1);
        std::vector<int> first(G.n);
        std::vector<double> working_row(G.n, 0);

        // for each row in inv_K
        for (int i = 0; i < inv_K.size(); ++i) {
            // eliminate e
            int e = perm[i];
            //std::cout << "eliminating node " << i << " with id " << e << std::endl;
            // insert value of row[e] into tmp
            for (int j = xadj[e]; j < xadj[e + 1]; ++j) {
                if (inv_perm[adj[j]] > i) {
                    working_row[inv_perm[adj[j]]] = weights[j];
                }
            }
            // get rows that effect row e
            int linked_row = link[i];
            while (linked_row != -1) {
                double factor = M.lnz[M.xlnz[linked_row] + first[linked_row]] / weighted_degree[perm[linked_row]];
                weighted_degree[e] -= M.lnz[M.xlnz[linked_row] + first[linked_row]] * factor;
                int c = 1;
                for (int k = M.xlnz[linked_row] + first[linked_row] + 1; k < M.xlnz[linked_row + 1]; ++k) {
                    int col = M.nzsub[M.xnzsub[linked_row] + first[linked_row] + c];
                    working_row[col] -= M.lnz[k] * factor;
                    ++c;
                }

                int llinked_row = linked_row;
                linked_row = link[linked_row];

                ++first[llinked_row];
                if (first[llinked_row] < M.xlnz[llinked_row + 1] - M.xlnz[llinked_row]) {
                    int new_linked_row = M.nzsub[M.xnzsub[llinked_row] + first[llinked_row]];
                    int old_linked_row = link[new_linked_row];
                    link[new_linked_row] = llinked_row;
                    link[llinked_row] = old_linked_row;
                }
            }

            first[i] = 0;
            int new_linked_row = M.nzsub[M.xnzsub[i]];
            // std::cout << "linking with row " << new_linked_row << std::endl;
            int old_linked_row = link[new_linked_row];
            link[new_linked_row] = i;
            link[i] = old_linked_row;

            // store into matrix structure + reset tmp vector
            int c = 0;
            for (int j = M.xlnz[i]; j < M.xlnz[i + 1]; ++j) {
                int col = M.nzsub[M.xnzsub[i] + c];
                //std::cout << "using col " << col << std::endl;
                if (j == M.lnz.size()) {
                    std::cout << j << std::endl;
                }
                M.lnz[j] = working_row[col];
                //std::cout << "updating " << i << "  " << M.nzsub[M.xnzsub[i]+c] << " value " << working_row[col] << std::endl;
                working_row[col] = 0;
                ++c;
            }
        }

        // add code for applying elimination to nodes in K (triangular part)
        DiGraph sparsifier(K.size() + 1);
        for (int i = 0; i < K.size(); ++i) {
            int e = K[i];

            // insert value of row[e] into tmp
            for (int j = xadj[e]; j < xadj[e + 1]; ++j) {
                if (inv_perm[adj[j]] > inv_K.size() + i) {
                    working_row[inv_perm[adj[j]]] = weights[j];
                    //std::cout << "setting col " << inv_perm[adj[j]] << " with weight " << weights[j] << std::endl;
                }
            }

            int linked_row = link[inv_K.size() + i];
            while (linked_row != -1) {
                double factor = M.lnz[M.xlnz[linked_row] + first[linked_row]] / weighted_degree[perm[linked_row]];
                weighted_degree[e] -= M.lnz[M.xlnz[linked_row] + first[linked_row]] * factor;
                int c = 1;
                for (int k = M.xlnz[linked_row] + first[linked_row] + 1; k < M.xlnz[linked_row + 1]; ++k) {
                    int col = M.nzsub[M.xnzsub[linked_row] + first[linked_row] + c];
                    working_row[col] -= M.lnz[k] * factor;
                    //std::cout << "setting col " << col << " new weight " << working_row[col] << std::endl;
                    ++c;
                }

                int llinked_row = linked_row;
                linked_row = link[linked_row];

                ++first[llinked_row];
                if (first[llinked_row] < M.xlnz[llinked_row + 1] - M.xlnz[llinked_row]) {
                    int new_linked_row = M.nzsub[M.xnzsub[llinked_row] + first[llinked_row]];
                    int old_linked_row = link[new_linked_row];
                    link[new_linked_row] = llinked_row;
                    link[llinked_row] = old_linked_row;
                }
            }

            // store into matrix structure + reset tmp vector
            // will be dense anyways so just insert and reset whole temporary vector!
            for (int j = i + 1; j < K.size(); ++j) {
                int col = inv_K.size() + j;
                if (working_row[col] != 0) {
                    sparsifier.add_edge(i, j, -working_row[col]);
                    sparsifier.add_edge(j, i, -working_row[col]);
                    working_row[col] = 0;
                }
            }
        }
        for (int i = 0; i < K.size(); ++i) {
            double sum = 0;
            for (auto &j: sparsifier.adjList[i]) {
                sum += j.second;
            }
            sparsifier.add_edge(i, i, initial_diag[K[i]] - weighted_degree[K[i]]);
            double sink_weight = (weighted_degree[K[i]] - sum) - alpha * initial_diag[K[i]];
            sparsifier.add_edge(i, K.size(), sink_weight);
        }

        return sparsifier;
    }

}