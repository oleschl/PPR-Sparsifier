#include <vector>
#include <unordered_map>

#include "exact_ppr_sparsifier/NodeRemoval.h"
#include "utility/utility.h"
#include "utility/orderings.h"
#include "utility/min_degree_pq.h"

namespace NodeRemoval {

    // version for undirected input graphs
    // even though normalized graph is non symmetric, the non zero elements stay symmetric
    // incoming and outgoing adjacent nodes are the same
    DiGraph constructPPRSparsifier(const GEdge &G, std::vector<int>& K, double alpha, const std::string& order) {
        // get V\K, the set of non-terminals
        auto inv_K = getNonTerminals(G.n, K);
        // precompute elimination order
        std::pair<std::vector<int>, std::vector<int>> ordering;
        if (order == "random"){
            ordering = getRandomOrdering(G, K, inv_K);
        } else if (order == "static") {
            ordering = getStaticMinDegOrdering(G, K, inv_K);
        } else if (order == "dynamic") {
            ordering = getDynamicMinDegOrdering(G, K, inv_K);
        }
        auto [perm, invPerm] = ordering;
        // compute weighted degrees
        std::vector<double> weighted_degrees(G.n, 0);
        for (auto edge: G.edges) {
            weighted_degrees[edge.u] += edge.weight;
            weighted_degrees[edge.v] += edge.weight;
        }
        // create normalized elimination graph H with additional sink node
        std::vector<std::unordered_map<int, double>> H(G.n + 1);
        for (auto edge: G.edges) {
            H[edge.u][edge.v] = edge.weight / weighted_degrees[edge.u];
            H[edge.v][edge.u] = edge.weight / weighted_degrees[edge.v];
        }
        H[G.n][G.n] = 1;
        invPerm[G.n] = G.n;
        // node removal: replace j -> i -> k by j -> k
        for (int i = 0; i < inv_K.size(); ++i) {
            int e_node = perm[i];
            //std::cout << "removing node " << i << std::endl;
            // for all j --> i (using the fact that the unweighted graph is symmetric)
            for (auto j: H[e_node]) {
                if (j.first == G.n || j.first == e_node) continue;
                double ji = H[j.first][e_node];
                double new_weight = 0;
                // and for i --> k
                for (auto k: H[e_node]) {
                    if (k.first == e_node) continue;
                    // compute new edge weight, account for moving from x --> z --> y, x --> z --> z --> --> y, ...
                    // (1-alpha) * w_G(x,z) * w_G(z,y) * sum_(t=0)^∞ ((1 - a) * w_G(z,z))^t
                    // sum_(t=0)^∞ ((1 - a) * w_G(z,z))^t = 1/(a * w_G(z,z) - w_G(z,z) + 1) when abs((1 - a) * w_G(z,z))<1
                    // do not remove next line! If H[x][x] does not exist, it will be inserted and invalidate our loop iterators
                    double selfLoop = H[e_node].contains(e_node) ? H[e_node][e_node] : 0;
                    double new_edge_weight = (1 - alpha) * ji * k.second * (1.0 / (alpha * selfLoop - selfLoop + 1.0));
                    new_weight += new_edge_weight;
                    // create new edge x --> y
                    H[j.first][k.first] += new_edge_weight;
                }
                // add missing weight (accounts for moving from x to z and restarting at z)
                H[j.first][G.n] += (ji - new_weight);
                H[j.first].erase(e_node);
            }
        }

        DiGraph sparsifier(K.size()+1);
        for (int i = 0; i < K.size(); ++i) {
            for (auto e: H[K[i]]) {
                sparsifier.add_edge(i, invPerm[e.first] - inv_K.size(), e.second);
            }
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
        // compute weighted degrees
        std::vector<double> weighted_degrees(G.n, 0);
        for(auto edge : G.edges){
            weighted_degrees[edge.u] += edge.weight;
            weighted_degrees[edge.v] += edge.weight;
        }
        // create normalized elimination graph H with additional sink node
        std::vector<std::unordered_map<int, double>> H(G.n+1);
        for(auto edge : G.edges){
            H[edge.u][edge.v] = edge.weight/weighted_degrees[edge.u];
            H[edge.v][edge.u] = edge.weight/weighted_degrees[edge.v];
        }
        H[G.n][G.n] = 1;

        inv_perm.push_back(G.n);
        // node removal: replace j -> i -> k by j -> k
        for(int i = 0; i < inv_K.size(); ++i){
            int e_node = inv_K[pq.pop()];
            //std::cout << "removing node " << i << std::endl;
            // for all j --> i (using the fact that the unweighted graph is symmetric)
            for(auto j : H[e_node]){
                if(j.first == G.n || j.first == e_node) continue;
                double ji = H[j.first][e_node];
                double new_weight = 0;
                // and for i --> k
                for(auto k : H[e_node]) {
                    if(k.first == e_node) continue;
                    // compute new edge weight, account for moving from x --> z --> y, x --> z --> z --> --> y, ...
                    // (1-alpha) * w_G(x,z) * w_G(z,y) * sum_(t=0)^∞ ((1 - a) * w_G(z,z))^t
                    // sum_(t=0)^∞ ((1 - a) * w_G(z,z))^t = 1/(a * w_G(z,z) - w_G(z,z) + 1) when abs((1 - a) * w_G(z,z))<1
                    // do not remove next line! If H[x][x] does not exist, it will be inserted and invalidate our loop iterators
                    double selfLoop = H[e_node].contains(e_node)? H[e_node][e_node] : 0;
                    double new_edge_weight = (1-alpha) * ji * k.second * (1.0/(alpha * selfLoop - selfLoop + 1.0));
                    new_weight += new_edge_weight;
                    // create new edge x --> y
                    H[j.first][k.first] += new_edge_weight;
                }
                // add missing weight (accounts for moving from x to z and restarting at z)
                H[j.first][G.n] += (ji - new_weight);
                H[j.first].erase(e_node);
                // if j is non-terminal
                if(inv_perm[j.first] < inv_K.size()){
                    // then update degree
                    pq.update(inv_perm[j.first], H[j.first].size());
                }
            }
        }

        DiGraph sparsifier(K.size()+1);
        for(int i = 0; i < K.size(); ++i){
            for(auto e : H[K[i]]){
                sparsifier.add_edge(i, inv_perm[e.first]-inv_K.size(), e.second);
            }
        }

        return sparsifier;
    }
}
