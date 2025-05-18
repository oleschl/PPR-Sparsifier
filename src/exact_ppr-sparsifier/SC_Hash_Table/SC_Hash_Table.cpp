#include <vector>
#include <unordered_map>
#include <map>
#include "exact_ppr_sparsifier/SC_Hash_Table.h"
#include "utility/utility.h"
#include "utility/orderings.h"
#include "utility/min_degree_pq.h"

// implementation for computing exact PPR sparsifiers using a matrix augmentation and the Schur complement
// this version computes the schur complement by repeatedly eliminating vertices representing the elimination
// graph as a list of hash tables
namespace SC_Hash_Table {
    DiGraph constructPPRSparsifier(const GEdge &G, std::vector<int>& K, double alpha, const std::string& order) {
        // get V\K, the set of non-terminals
        auto inv_K = getNonTerminals(G.n, K);
        // compute static minimum degree ordering
        std::pair<std::vector<int>, std::vector<int>> ordering;
        if (order == "random"){
            ordering = getRandomOrdering(G, K, inv_K);
        } else if (order == "static") {
            ordering = getStaticMinDegOrdering(G, K, inv_K);
        } else if (order == "dynamic") {
            ordering = getDynamicMinDegOrdering(G, K, inv_K);
        }
        auto [perm, invPerm] = ordering;
        // create graph structure for elimination
        std::vector<std::unordered_map<int, double>> H(G.n);
        std::vector<double> weighted_degree(G.n, 0);
        std::vector<double> initial_weighted_degrees(G.n, 0);
        for (auto edge: G.edges) {
            if (invPerm[edge.u] < invPerm[edge.v]) {
                H[invPerm[edge.u]][invPerm[edge.v]] = edge.weight;
                weighted_degree[invPerm[edge.u]] += edge.weight;
            } else {
                H[invPerm[edge.v]][invPerm[edge.u]] = edge.weight;
                weighted_degree[invPerm[edge.v]] += edge.weight;
            }
            initial_weighted_degrees[invPerm[edge.u]] += edge.weight;
            initial_weighted_degrees[invPerm[edge.v]] += edge.weight;
        }
        // perform elimination on D-(1-alpha)A by lifting with additional sink node
        // instead of scaling A, do elimination on (alpha/(1-alpha))D-A
        perm.push_back(G.n);
        invPerm.push_back(G.n);
        for (int i = 0; i < G.n; ++i) {
            double edge_weight = (alpha / (1 - alpha)) * initial_weighted_degrees[i];
            H[i][G.n] = edge_weight;
            weighted_degree[i] += edge_weight;
            initial_weighted_degrees[i] += edge_weight;
        }

        // schur complement / elimination
        for (int i = 0; i < inv_K.size(); ++i) {
            // std::cout << "eliminating node " << i << " with id " << perm[i] << std::endl;
            for (auto e1 = H[i].begin(); e1 != H[i].end(); ++e1) {
                for (auto e2 = std::next(e1); e2 != H[i].end(); ++e2) {
                    auto newEdgeWeight = e1->second * e2->second / weighted_degree[i];
                    if (e2->first > e1->first) {
                        H[e1->first][e2->first] += newEdgeWeight;
                        weighted_degree[e1->first] += newEdgeWeight;
                    } else {
                        H[e2->first][e1->first] += newEdgeWeight;
                        weighted_degree[e2->first] += newEdgeWeight;
                    }
                }
            }
        }

        // compute new degree for each node
        // diff to initial degrees corresponds to self loops in graph
        std::vector<double> final_weighted_degree(K.size() + 1, 0);
        for (int i = 0; i < K.size(); ++i) {
            for (auto e: H[inv_K.size() + i]) {
                final_weighted_degree[i] += e.second;
                final_weighted_degree[e.first - inv_K.size()] += e.second;
            }
        }

        DiGraph sparsifier(K.size()+1);
        for (int i = 0; i < K.size(); ++i) {
            int node = inv_K.size() + i;
            for (auto e: H[node]) {
                if (e.first != G.n) {
                    if(e.second == 0) std::cout << "zero weight" << std::endl;
                    sparsifier.add_edge(i, e.first - inv_K.size(), e.second);
                    sparsifier.add_edge(e.first - inv_K.size(), i, e.second);
                } else {
                    // remove initial weight of lifted node
                    double new_weight = e.second - alpha * initial_weighted_degrees[node];
                    sparsifier.add_edge(i, K.size(), new_weight);
                }
            }
            // add self loops
            sparsifier.add_edge(i, i, initial_weighted_degrees[node] - final_weighted_degree[i]);
        }

        return sparsifier;
    }

    DiGraph constructPPRSparsifier2(const GEdge& G, std::vector<int>& K, double alpha){

        // get V\K, the set of non-terminals
        auto inv_K = getNonTerminals(G.n, K);
        // compute initial degrees
        std::vector<int> degrees(G.n, 0);
        for(auto edge : G.edges){
            ++degrees[edge.u];
            ++degrees[edge.v];
        }
        // create minimum degree data structure
        MinDegreePQ pq(G.n, inv_K.size());
        for(int i = 0; i < inv_K.size(); ++i){
            pq.add(i, degrees[inv_K[i]]);
        }
        // mapping from node to elimination order, here it is used to divide nodes in inv_K from nodes in K
        std::vector<int> inv_perm(G.n+1);
        for(int i = 0; i < inv_K.size(); ++i){
            inv_perm[inv_K[i]] = i;
        }
        for(int i = 0; i < K.size(); ++i){
            inv_perm[K[i]] = inv_K.size() + i;
        }

        // next create elimination graph H as list of maps
        std::vector<std::unordered_map<int, double>> H(G.n+1);
        std::vector<double> weighted_degree(G.n+1,0);
        // since we do not know elimination in advance, we need to store the whole matrix
        for(auto edge : G.edges){
            H[inv_perm[edge.u]][inv_perm[edge.v]] = edge.weight;
            weighted_degree[inv_perm[edge.u]] += edge.weight;
            H[inv_perm[edge.v]][inv_perm[edge.u]] = edge.weight;
            weighted_degree[inv_perm[edge.v]] += edge.weight;
        }

        // perform elimination on D-(1-alpha)A by lifting with additional sink node
        // instead of scaling A, do elimination on (alpha/(1-alpha))D-A
        for(int i = 0; i < G.n; ++i){
            H[i][G.n] = (alpha/(1-alpha)) * weighted_degree[i];
            weighted_degree[i] += H[i][G.n];
        }
        std::vector<double> initial_weighted_degrees = weighted_degree;

        // schur complement computation of undirected graph
        for(int i = 0; i < inv_K.size(); ++i){
            //std::cout << "removing node " << i << std::endl;
            // get vertex with minimum degree
            int eliminationVertex = pq.pop();
            // create clique between all pairs of adjacent nodes
            for(auto e1 = H[eliminationVertex].begin(); e1 != H[eliminationVertex].end(); ++e1){
                for(auto e2 = std::next(e1); e2 != H[eliminationVertex].end(); ++e2){
                    auto newEdgeWeight = e1->second * e2->second / weighted_degree[eliminationVertex];
                    H[e1->first][e2->first] += newEdgeWeight;
                    weighted_degree[e1->first] += newEdgeWeight;
                    H[e2->first][e1->first] += newEdgeWeight;
                    weighted_degree[e2->first] += newEdgeWeight;
                }
            }

            // remove star around eliminationVertex
            for (auto e : H[eliminationVertex]) {
                weighted_degree[e.first] -= e.second;
                H[e.first].erase(eliminationVertex);
            }

            for(auto e : H[eliminationVertex]){
                // if e.first not in K
                if(e.first < inv_K.size()){
                    pq.update(e.first, H[e.first].size());
                }
            }
        }

        // remove initial weight of lifted
        for(int i = inv_K.size(); i < G.n; ++i){
            H[i][G.n] -= alpha * initial_weighted_degrees[i];
        }
        // transform to output format
        DiGraph sparsifier(K.size()+1);
        for(int i = 0; i < K.size(); ++i){
            for(auto e : H[inv_K.size()+i]){
                sparsifier.add_edge(i, e.first-inv_K.size(), e.second);
            }
            sparsifier.add_edge(i, i , initial_weighted_degrees[inv_K.size()+i] - weighted_degree[inv_K.size()+i]);
        }

        return sparsifier;
    }
}