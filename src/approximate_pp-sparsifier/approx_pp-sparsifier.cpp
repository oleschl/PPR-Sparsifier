#include <vector>
#include <random>
#include <iostream>
#include <algorithm>

#include "approximate_ppr_sparsifier/approximate_pp_sparsifier.h"
#include "utility/utility.h"
#include "utility/min_degree_pq.h"

/*
https://student.cs.uwaterloo.ca/~cs466/notes/Notes_5.1.1_DynamicMinDegree.pdf
*/
namespace ApproximateSparsifier {

    //std::random_device rd;
    //std::mt19937 gen(rd());
    std::mt19937 gen(42);

    struct col_el {
        int row;
        double v;
        col_el *next;
        col_el *rev;
    };

    struct col_el2 {
        int row;
        double v;
        col_el2 *next;
        col_el2 *prev;
        col_el2 *rev;
    };

    struct EliminationGraph {
        int n;
        std::vector<col_el *> cols;
    };

    struct new_edge {
        int row, col;
        double v;
    };

    bool compareByValueThenRow(const col_el *a, const col_el *b) {
        return a->v > b->v || (b->v == a->v && a->row < b->row);
    }

    bool compareByRow(const col_el *a, const col_el *b) {
        return a->row < b->row;
    }

    int compressAndAverageColumn(std::vector<col_el *> &column, int merge, MinDegreePQ &pq, int numNonTerminals) {
        // if there is no or just one edge we do not need to order compress it
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

    // consider doing differnt elimination order (as described in paper)
    void sample_clique2(EliminationGraph& G, int v, MinDegreePQ& pq, int numNonTerminals) {
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
        //std::uniform_int_distribution<> weighted_dist(0, length - 1);
        std::uniform_int_distribution<> uniform_dist(0, length - 1);
        // sample edges and add to graph
        for(int i = 0; i < length; ++i) {
            auto j_idx = weighted_dist(gen);
            auto k_idx = uniform_dist(gen);

            if (column[j_idx]->row != column[k_idx]->row) {
                auto new_weight = column[j_idx]->v * column[k_idx]->v / (column[j_idx]->v + column[k_idx]->v);
                // Create new edges as it is difficult to reuse existing ones
                // (requires doubly linked list or careful mapping from existing to new edges (both appears expensive))
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
            // per edge and thus dont have to deal with deletion of reverse edges
            column[i]->rev->v = 0;
            //column[i]->rev->rev = nullptr;
            if (column[i]->row < numNonTerminals) pq.decrease(column[i]->row);
            delete column[i];
        }
    }

    void sample_clique(EliminationGraph &G, int v, MinDegreePQ &pq, int merge, int numNonTerminals) {
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
            auto it = std::upper_bound(csumRev.begin() + next, csumRev.end(), r,std::greater<double>());
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

/*
// K - id of nodes that should be eliminated in sorted order
Graph approx_chol(Graph& G, std::vector<int> K, int split, int merge){
    EliminationGraph H = {G.n, std::vector<col_el*>(G.n, nullptr)};
    std::vector<int> invK;
    std::vector<int> degrees(G.n, 0);
    int cur_k = 0;
    // Replace each edge e by split parallel edges each with weight w_e/split
    for (int i = 0; i < G.n; ++i) {
        for (int j = 0; j < G.adjList[i].size(); ++j) {
            for(int k = 0; k < split; ++k){
                col_el* prev_i = H.cols[i];
                col_el* prev_j = H.cols[G.adjList[i][j].first];
                auto* ij = new col_el{G.adjList[i][j].first, G.adjList[i][j].second/split, prev_i, nullptr};
                auto* ji = new col_el{i, G.adjList[i][j].second/split, prev_j, ij};
                H.cols[i] = ij;
                H.cols[G.adjList[i][j].first] = ji;
                H.cols[i]->rev = ji;
            }
            degrees[i] += split;
            degrees[G.adjList[i][j].first] += split;
        }
        if(cur_k < K.size() && i == K[cur_k]){
            ++cur_k;
        } else {
            invK.push_back(i);
        }
    }

    min_degree_pq pq(invK, degrees, G.n, split);

    for(int i = 0; i < invK.size(); ++i){
        auto node_i = pq.popMinDegree();
        std::cout<< "deleted K[i]: "<< i << std::endl;
        sample_clique(H, node_i, pq, merge);
    }

    // TODO: convert multi graph to normal graph
    Graph newG(K.size());
    // need mapping from G to nodes in sparsifier and backwards
    std::vector<int> mapping(G.n);
    int curK = 0;
    for(int i = 0; i < G.n; ++i) {
        if(i == K[curK]) {
            mapping[i] = curK;
            ++curK;
        } else {
            mapping[i] = -1;
        }
    }

    for (int i = 0; i < K.size(); ++i) {
        int columnIndex = K[i];
        auto columnPtr = H.cols[columnIndex];
        std::vector<col_el*> column = {};
        while(columnPtr){
            column.push_back(columnPtr);
            columnPtr = columnPtr->next;
        }
        std::sort(column.begin(), column.end(), compareByRow);

        int prevRow = -1;
        for(auto & j : column) {
            if(i < mapping[j->row]) {
                if(prevRow != j->row){
                    newG.adjList[i].emplace_back(mapping[j->row], j->v);
                    prevRow = j->row;
                } else {
                    newG.adjList[i].back().second += j->v;
                }
            }
        }
    }

    return newG;
}
 */
    DiGraph constructPPRSparsifier(const GEdge &G, std::vector<int> &K, double alpha, int split, int merge) {
        auto inv_K = getNonTerminals(G.n, K);
        // mapping from node to elimination order, here it is used to divide nodes in inv_K from nodes in K
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

        for (int i = 0; i < inv_K.size(); ++i) {
            auto eliminationVertex = pq.pop();
            std::cout << "deleted K[i]: " << i << "vertexid: " << eliminationVertex << std::endl;
            //sample_clique(H, eliminationVertex, pq, merge, inv_K.size());
            sample_clique2(H, eliminationVertex, pq, inv_K.size());
        }

        // TODO: convert multi graph to normal graph
        DiGraph sparsifier(K.size() + 1);
        // TODO add code for reduce weight to sink node!

        for (int i = 0; i < K.size(); ++i) {
            auto columnPtr = H.cols[inv_K.size() + i];
            std::vector<col_el *> column = {};
            while (columnPtr) {
                column.push_back(columnPtr);
                columnPtr = columnPtr->next;
            }
            std::sort(column.begin(), column.end(), compareByRow);
            // TODO test if it contains rows that are non terminals

            if (column.empty()) continue;
            int prevRow = column[0]->row;
            double new_weight = 0;
            double new_weighted_degree = 0;
            for (const auto &j: column) {
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
            // add self loop if degree has decreased
            if (weighted_degree[K[i]] - new_weighted_degree > 0) {
                sparsifier.add_edge(i, i, weighted_degree[K[i]] - new_weighted_degree);
            }
        }

        return sparsifier;
    }

/*
Sparsifier getApproximateSparsifier(GEdge& G, std::vector<int> K, int split, int merge, double alpha) {
    Graph H(G);
    std::vector<double> weightedDegree(G.n + 1, 0);

    // for mapping from K to G, K can be used if it is ordered!

    // compute (1-alpha) A
    for(int i = 0; i < G.n; ++i) {
        for(int j = 0; j < G.adjList[i].size(); ++j) {
            weightedDegree[i] += G.adjList[i][j].second;
            weightedDegree[G.adjList[i][j].first] += G.adjList[i][j].second;
            H.adjList[i][j].second *= alpha;
        }
    }

    auto max = *std::min_element(weightedDegree.begin(), weightedDegree.end());
    std::cout << *std::max(weightedDegree.begin(), weightedDegree.end()) << std::endl;

    // add new sink node with weights alpha * outdegree
    int sinkID = G.n;
    for(int i = 0; i < G.n; ++i) {
        H.add_edge(i, sinkID, (1-alpha) * weightedDegree[i]);
        weightedDegree[sinkID] += (1-alpha) * weightedDegree[i];
    }
    // call approx cholesky
    // maybe we must add sink node to K?!
    std::vector<int> newK(K);
    newK.push_back(sinkID);
    H = approx_chol(H, newK, split, merge);
    DiGraph sparsifier;
    // H must be a directed graph!
    // scale graph by D^-1
    std::vector<double> newWeightedDegrees(H.n, 0);
    for(int i = 0; i < H.n; ++i) {
        for(int j = 0; j < H.adjList[i].size(); ++j) {
            sparsifier.add_edge(i, H.adjList[i][j].first, H.adjList[i][j].second/weightedDegree[K[i]]);
            sparsifier.add_edge(H.adjList[i][j].first, i, H.adjList[i][j].second/weightedDegree[newK[H.adjList[i][j].first]]);
            newWeightedDegrees[i] += H.adjList[i][j].second/weightedDegree[K[i]];
            newWeightedDegrees[H.adjList[i][j].first] += H.adjList[i][j].second/weightedDegree[newK[H.adjList[i][j].first]];
        }
    }
    // add self loops
    for(int i = 0; i < H.n; ++i) {
        if(1-newWeightedDegrees[i] > 0) {
            sparsifier.adjList[i].emplace_back(i, 1-newWeightedDegrees[i]);
        }
    }
    // reduce edge to sink by (1-alpha)
    for(int i = 0; i < H.n; ++i) {
        for(int j = 0; j < sparsifier.adjList[i].size()-1; ++j){
            if(sparsifier.adjList[i][j].first == K.size()) {
                sparsifier.adjList[i][j].second -= (1-alpha);
            }
        }
    }
    // remove outoging edges of sink node
    sparsifier.adjList.back() = {};

    return {sparsifier, K};
}*/
}