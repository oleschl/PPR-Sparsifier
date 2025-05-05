#include <vector>
#include <random>
#include <iostream>
#include <algorithm>

#include "approximate_ppr_sparsifier/approximate_pp_sparsifier.h"
#include "utility/utility.h"
#include "utility/min_degree_pq.h"

namespace ApproximateSparsifier {

    //std::random_device rd;
    //std::mt19937 gen(rd());
    std::mt19937 gen;

    // linked-list based graph structure
    struct col_el {
        int row;
        double v;
        col_el *next;
        col_el *rev;
    };

    struct EliminationGraph {
        int n;
        std::vector<col_el *> cols;
    };

    bool compareByValueThenRow(const col_el *a, const col_el *b) {
        return a->v > b->v || (b->v == a->v && a->row < b->row);
    }

    bool compareByRow(const col_el *a, const col_el *b) {
        return a->row < b->row;
    }

    int compressAndAverageColumn(std::vector<col_el *> &column, int merge, MinDegreePQ &pq, int numNonTerminals) {
        // if there is no or just one edge we do not need to compress it
        if (column.size() <= 1) {
            return column.size();
        }

        std::sort(column.begin(), column.end(), compareByRow);

        int prev_edge_new = 0;
        int prev_edge = 0;
        int count = 1;
        double sum = column[0]->v;
        for (int i = 1; i < column.size(); ++i) {
            if (column[prev_edge]->row == column[i]->row) {
                ++count;
                sum += column[i]->v;
            } else {
                if (count <= merge) {
                    double newWeight = sum / count;
                    for (int j = prev_edge; j < prev_edge + count; ++j) {
                        column[j]->v = newWeight;
                        column[j]->rev->v = newWeight;
                        column[prev_edge_new] = column[j];
                        ++prev_edge_new;
                    }
                } else {
                    double newWeight = sum / merge;
                    for (int j = prev_edge; j < prev_edge + merge; ++j) {
                        column[j]->v = newWeight;
                        column[j]->rev->v = newWeight;
                        column[prev_edge_new] = column[j];
                        ++prev_edge_new;
                    }
                    // decrease degree in pq of rev
                    for (int j = prev_edge + merge; j < prev_edge + count; ++j) {
                        column[j]->rev->v = 0;
                        if (column[j]->row < numNonTerminals) pq.decrease(column[j]->row);
                    }
                }
                count = 1;
                sum = column[i]->v;
                prev_edge = i;
            }
        }

        // add last edges of last row
        if (count <= merge) {
            double newWeight = sum / count;
            for (int j = prev_edge; j < prev_edge + count; ++j) {
                column[j]->v = newWeight;
                column[j]->rev->v = newWeight;
                column[prev_edge_new] = column[j];
                ++prev_edge_new;
            }
        } else {
            double newWeight = sum / merge;
            for (int j = prev_edge; j < prev_edge + merge; ++j) {
                column[j]->v = newWeight;
                column[j]->rev->v = newWeight;
                column[prev_edge_new] = column[j];
                ++prev_edge_new;
            }
            // decrease degree in pq of rev
            for (int j = prev_edge + merge; j < prev_edge + count; ++j) {
                column[j]->rev->v = 0;
                if (column[j]->row < numNonTerminals) pq.decrease(column[j]->row);
            }
        }

        std::sort(column.begin(), column.begin() + prev_edge_new, compareByValueThenRow);
        return prev_edge_new;
    }

    // clique sampling procedure based on ...
    void sample_random_clique(EliminationGraph& G, int v, MinDegreePQ& pq, int merge, int numNonTerminals) {
        // load column
        std::vector<col_el*> column = {};
        auto cur = G.cols[v];
        while(cur != nullptr){
            if(cur->v != 0) {
                column.push_back(cur);
            }
            cur = cur->next;
        }
        int length = column.size();
        if(length == 0) return;

        std::vector<double> weights(length);
        for(int i = 0; i < length; ++i) {
            weights[i] = column[i]->v;
        }

        std::discrete_distribution<int> weighted_dist(weights.begin(), weights.end());
        std::uniform_int_distribution<> uniform_dist(0, length - 1);
        // sample edges and add to graph
        for(int i = 0; i < length; ++i) {
            auto j_idx = weighted_dist(gen);
            auto k_idx = uniform_dist(gen);

            if (column[j_idx]->row != column[k_idx]->row) {
                auto new_weight = column[j_idx]->v * column[k_idx]->v / (column[j_idx]->v + column[k_idx]->v);
                // Create new edges as it is difficult to reuse existing ones
                // (requires doubly linked list or careful mapping from existing to new edges (both appears too expensive))
                auto* jk = new col_el{column[k_idx]->row, new_weight, G.cols[column[j_idx]->row], nullptr};
                G.cols[column[j_idx]->row] = jk;
                if (column[j_idx]->row < numNonTerminals) pq.increase(column[j_idx]->row);

                auto* kj = new col_el{column[j_idx]->row, new_weight, G.cols[column[k_idx]->row], jk};
                G.cols[column[k_idx]->row] = kj;
                if (column[k_idx]->row < numNonTerminals) pq.increase(column[k_idx]->row);

                jk->rev = kj;
            }
        }
        // remove column i and its reverse edges
        for(int i = 0; i < length; ++i){
            // we indicate deletion by setting the weight to zero
            // deletion would be more memory friendly but requires using a double linked list
            // problem becomes simpler when elimination order is known advance as we could then store only one direction
            // per edge and thus do not have to deal with deletion of reverse edges
            column[i]->rev->v = 0;
            if (column[i]->row < numNonTerminals) pq.decrease(column[i]->row);
            delete column[i];
        }
    }

    // clique sampling procedure based on ...
    void sample_elim_star(EliminationGraph &G, int v, MinDegreePQ &pq, int merge, int numNonTerminals) {
        // load column
        std::vector<col_el *> column = {};
        auto cur = G.cols[v];
        while (cur != nullptr) {
            if (cur->v != 0) {
                column.push_back(cur);
            }
            cur = cur->next;
        }

        int length = compressAndAverageColumn(column, merge, pq, numNonTerminals);

        if (length == 0) return;

        std::vector<double> weights(length);
        double S = 0;
        for (int i = 0; i < length; ++i) {
            S += column[i]->v;
            weights[i] = column[i]->v;
        }
        double sum = S;

        std::vector<double> csumRev(length+1, 0);
        double csum = 0;
        for (int i = length - 1; i >= 0; --i) {
            csum += column[i]->v;
            csumRev[i] = csum;
        }

        int next = 1;
        for (int i = 0; i < length - 1; ++i) {
            S -= column[i]->v;
            // avoid sampling multi edges
            while (next < length && column[i]->row == column[next]->row) {
                ++next;
            }
            // take care if only multi edges of same edge are left by removing all multi edges and decreasing degree in pq
            if (next == length) {
                while (i < length - 1) {
                    column[i]->rev->v = 0;
                    if (column[i]->row < numNonTerminals) pq.decrease(column[i]->row);
                    ++i;
                }
                break;
            }
            // sample other node
            std::uniform_real_distribution<double> dist(0.0, csumRev[next]);
            double r = dist(gen);
            auto it = std::upper_bound(csumRev.begin() + next, csumRev.end(), r,std::greater<>());
            int sample = std::distance(csumRev.begin(), it) - 1;
            // reassign edge j -> i to j -> k; degree stays the same
            column[i]->rev->v = (csumRev[next] * column[i]->v) / sum;
            column[i]->rev->row = column[sample]->row;
            // add reverse edge k -> j by reusing edge i->j and increment degree
            column[i]->v = (csumRev[next] * column[i]->v) / sum;
            column[i]->next = G.cols[column[sample]->row];
            G.cols[column[sample]->row] = column[i];
            // if node in K
            if (column[sample]->row < numNonTerminals) pq.increase(column[sample]->row);
        }

        // for last element no new edge will be added -> decrease degree
        column[length - 1]->rev->v = 0;
        if (column[length - 1]->row < numNonTerminals) pq.decrease(column[length - 1]->row);
    }

    // computes an approximate personalized PageRank sparsifier
    DiGraph constructPPRSparsifier(const GEdge &G, std::vector<int> &K, double alpha, int split, int merge, int seed, const std::string& sampling) {
        gen.seed(seed);
        auto inv_K = getNonTerminals(G.n, K);
        // mapping from node to elimination order, here it is used to divide nodes in terminal nodes K from
        // non-terminal nodes inv_K
        std::vector<int> inv_perm(G.n + 1);
        for (int i = 0; i < inv_K.size(); ++i) {
            inv_perm[inv_K[i]] = i;
        }
        for (int i = 0; i < K.size(); ++i) {
            inv_perm[K[i]] = inv_K.size() + i;
        }
        inv_perm[G.n] = G.n;
        // compute degrees
        std::vector<int> degrees(G.n, split);
        std::vector<double> weighted_degree(G.n, 0);
        for (const auto &edge: G.edges) {
            degrees[edge.u] += split;
            degrees[edge.v] += split;
            weighted_degree[edge.u] += edge.weight;
            weighted_degree[edge.v] += edge.weight;
        }
        // create elimination graph and split edges
        EliminationGraph H = {G.n + 1, std::vector<col_el *>(G.n + 1, nullptr)};
        for (const auto &edge: G.edges) {
            for (int i = 0; i < split; ++i) {
                col_el *prev_u = H.cols[inv_perm[edge.u]];
                col_el *prev_v = H.cols[inv_perm[edge.v]];
                auto *uv = new col_el{inv_perm[edge.v], edge.weight / split, prev_u, nullptr};
                auto *vu = new col_el{inv_perm[edge.u], edge.weight / split, prev_v, uv};
                H.cols[inv_perm[edge.u]] = uv;
                H.cols[inv_perm[edge.v]] = vu;
                uv->rev = vu;
            }
        }
        // add lifted sink node
        for (int i = 0; i < G.n; ++i) {
            double weight = (alpha / (1 - alpha)) * weighted_degree[i] / split;
            for (int j = 0; j < split; ++j) {
                col_el *prev_u = H.cols[inv_perm[i]];
                col_el *prev_v = H.cols[G.n];
                auto *uv = new col_el{G.n, weight, prev_u, nullptr};
                auto *vu = new col_el{inv_perm[i], weight, prev_v, uv};
                H.cols[inv_perm[i]] = uv;
                H.cols[G.n] = vu;
                uv->rev = vu;
            }
            weighted_degree[i] += weight;
        }

        MinDegreePQ pq(G.n, inv_K.size());
        for (int i = 0; i < inv_K.size(); ++i) {
            pq.add(i, degrees[inv_K[i]]);
        }

        auto sample = (sampling == "random_clique")
                      ? sample_random_clique
                      : sample_elim_star;

        for (int i = 0; i < inv_K.size(); ++i) {
            auto eliminationVertex = pq.pop();
            // std::cout << "Eliminating vertex " << vertex << " at step " << i << std::endl;
            sample(H, eliminationVertex, pq, merge, inv_K.size());
        }

        // transform multi graph to sparsifier
        DiGraph sparsifier(K.size() + 1);
        for (int i = 0; i < K.size(); ++i) {
            auto columnPtr = H.cols[inv_K.size() + i];
            std::vector<col_el *> column = {};
            while (columnPtr) {
                if(columnPtr->row >= inv_K.size())
                    column.push_back(columnPtr);
                columnPtr = columnPtr->next;
            }
            if (column.empty()) continue;
            std::sort(column.begin(), column.end(), compareByRow);

            int prevRow = column[0]->row;
            double new_weight = 0;
            double new_weighted_degree = 0;
            for (const auto &j: column) {
                // Consider only terminal-to-terminal edges
                if (j->row >= inv_K.size()) {
                    if (prevRow != j->row) {
                        sparsifier.add_edge(i, prevRow - inv_K.size(), new_weight);
                        new_weighted_degree += new_weight;
                        prevRow = j->row;
                        new_weight = 0.0;
                    }
                    new_weight += j->v;
                }
            }

            // add edge to source node (after sorting last edge always points to source node)
            // and remove initial weight of lifted node
            sparsifier.add_edge(i, prevRow - inv_K.size(), new_weight - alpha * weighted_degree[K[i]]);

            // add self loop if weighted degree has decreased
            if (weighted_degree[K[i]] - new_weighted_degree > 0) {
                sparsifier.add_edge(i, i, weighted_degree[K[i]] - new_weighted_degree);
            }
        }

        return sparsifier;
    }
}