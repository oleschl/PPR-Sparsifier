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

/*
// version for undirected input graphs
// even though normalized graph in non symmetric, the non zero elements stay symmetric
// incoming and outgoing adjacent nodes are the same
Sparsifier nodeRemoval(GEdge& G, std::vector<int> K, double alpha){

    // get V\K, the set of non-terminals
    auto inv_K = getNonTerminals(G.n, K);
    // compute initial degrees + weighted degrees
    auto degrees = getDegrees2(G);
    std::sort(inv_K.begin(), inv_K.end(), [&degrees](int x, int y) {
        return degrees[x] < degrees[y];
    });
    auto weighted_degrees = getWeightedDegrees(G);

    // create normalized elimination graph H with additional sink node
    std::vector<std::unordered_map<int, double>> H(G.n+1);
    // since we do not know elimination in advance, we need to store the whole matrix
    // know partial order: terminal nodes after non-terminals
    for(int i = 0; i < G.n; ++i){
        for(auto e : G.adjList[i]) {
            H[i][e.first] = e.second/weighted_degrees[i];
            H[e.first][i] = e.second/weighted_degrees[e.first];
        }
    }

    H[G.n][G.n] = 1;

    int count = 0;
    // schur complement: replace j -> i -> k by j -> k
    for(int i : inv_K){
        std::cout << "removing node " << count << std::endl;
        ++count;
        // for all j --> i (using the fact that the unweighted graph is symmetric)
        for(auto j : H[i]){
            if(j.first == G.n || j.first == i) continue;
            double ji = H[j.first][i];
            double new_weight = 0;
            // and for i --> k
            for(auto k : H[i]) {
                if(k.first == i) continue;
                // compute new edge weight, account for moving from x --> z --> y, x --> z --> z --> --> y, ...
                // (1-alpha) * w_G(x,z) * w_G(z,y) * sum_(t=0)^∞ ((1 - a) * w_G(z,z))^t
                // sum_(t=0)^∞ ((1 - a) * w_G(z,z))^t = 1/(a * w_G(z,z) - w_G(z,z) + 1) when abs((1 - a) * w_G(z,z))<1
                // do not remove next line! If H[x][x] does not exist, it will be inserted and invalidate our loop iterators
                double selfLoop = H[i].contains(i)? H[i][i] : 0;
                double new_edge_weight = (1-alpha) * ji * k.second * (1.0/(alpha * selfLoop - selfLoop + 1.0));
                new_weight += new_edge_weight;
                // create new edge x --> y
                H[j.first][k.first] += new_edge_weight;
            }
            // add missing weight (accounts for moving from x to z and restarting at z)
            H[j.first][G.n] += (ji - new_weight);
            H[j.first].erase(i);
        }
    }

    // transform to output format
    std::vector<int> map_G_to_K(G.n+1, -1);
    for(int i = 0; i < K.size(); ++i){
        map_G_to_K[K[i]] = i;
    }
    map_G_to_K[G.n] = K.size();
    DiGraph sparsifier;
    for(int i = 0; i < K.size(); ++i){
        for(auto e : H[K[i]]){
            sparsifier.add_edge(i, map_G_to_K[e.first], e.second);
        }
    }

    return {sparsifier, K};
}
 */

/*
Sparsifier nodeRemoval2(Graph& G, std::vector<int> K, double alpha){

    // get V\K, the set of non-terminals
    auto inv_K = getNonTerminals2(G.n, K);
    // map ids in K and inv_K to range 0 to size-1
    std::vector<int> map_G_to_K(G.n+1, -1);
    for(int i = 0; i < K.size(); ++i){
        map_G_to_K[K[i]] = i;
    }
    map_G_to_K[G.n] = K.size();
    std::vector<int> map_G_to_inv_K(G.n+1);
    for(int i = 0; i < inv_K.size(); ++i){
        map_G_to_inv_K[inv_K[i]] = i;
    }

    // compute initial degrees + weighted degrees
    auto degrees = getDegrees2(G);
    auto weighted_degrees = getWeightedDegrees(G);

    // create normalized elimination graph H with additional sink node
    std::vector<std::unordered_map<int, double>> H(G.n+1);
    // since we do not know elimination in advance, we need to store the whole matrix
    // know partial order: terminal nodes after non-terminals
    for(int i = 0; i < G.n; ++i){
        for(auto e : G.adjList[i]) {
            H[i][e.first] = e.second/weighted_degrees[i];
            H[e.first][i] = e.second/weighted_degrees[e.first];
        }
    }

    H[G.n][G.n] = 1;

    std::vector<std::pair<int, int>> degreeInvK(inv_K.size());
    // create priority queue
    for(int i = 0; i < inv_K.size(); ++i){
        degreeInvK[i] = {i, degrees[inv_K[i]]};
    }
    MinDegreePQ pq(G.n, inv_K.size());

    int count = 0;
    // schur complement: replace j -> i -> k by j -> k
    for(int i = 0; i < inv_K.size(); ++i){
        int eliminationVertex = inv_K[pq.pop()];
        std::cout << "removing node " << count << std::endl;
        ++count;
        // for all j --> i (using the fact that the unweighted graph is symmetric)
        for(auto j : H[eliminationVertex]){
            if(j.first == G.n || j.first == eliminationVertex) continue;
            double ji = H[j.first][eliminationVertex];
            double new_weight = 0;
            // and for i --> k
            for(auto k : H[eliminationVertex]) {
                if(k.first == eliminationVertex) continue;
                // compute new edge weight, account for moving from x --> z --> y, x --> z --> z --> --> y, ...
                // (1-alpha) * w_G(x,z) * w_G(z,y) * sum_(t=0)^∞ ((1 - a) * w_G(z,z))^t
                // sum_(t=0)^∞ ((1 - a) * w_G(z,z))^t = 1/(a * w_G(z,z) - w_G(z,z) + 1) when abs((1 - a) * w_G(z,z))<1
                // do not remove next line! If H[x][x] does not exist, it will be inserted and invalidate our loop iterators
                double selfLoop = H[eliminationVertex].contains(eliminationVertex)? H[eliminationVertex][eliminationVertex] : 0;
                double new_edge_weight = (1-alpha) * ji * k.second * (1.0/(alpha * selfLoop - selfLoop + 1.0));
                new_weight += new_edge_weight;
                // create new edge x --> y
                H[j.first][k.first] += new_edge_weight;
            }
            // add missing weight (accounts for moving from x to z and restarting at z)
            H[j.first][G.n] += (ji - new_weight);
            H[j.first].erase(eliminationVertex);
            // if j is non-terminal
            if(map_G_to_K[j.first] == -1){
                // then update degree
                pq.update(map_G_to_inv_K[j.first], H[j.first].size());
            }
        }
    }

    DiGraph sparsifier;
    for(int i = 0; i < K.size(); ++i){
        for(auto e : H[K[i]]){
            sparsifier.add_edge(i, map_G_to_K[e.first], e.second);
            if(map_G_to_K[e.first] == -1){
                std::cout << e.first << std::endl;
            }
        }
    }

    return {sparsifier, K};
}
 */

/*

struct RowEl3{
    int col;
    RowEl3* next;
    RowEl3* rev;
    double v;

    RowEl3(int col, RowEl3* next, RowEl3* rev, double v) : col(col), next(next), rev(rev), v(v) {}
};

Sparsifier nodeRemoval3(Graph& G, std::vector<int> K, double alpha){

    // get V\K, the set of non-terminals
    auto inv_K = getNonTerminals2(G.n, K);
    // compute initial degrees + weighted degrees
    auto degrees = getDegrees2(G);
    std::sort(inv_K.begin(), inv_K.end(), [&degrees](int x, int y) {
        return degrees[x] < degrees[y];
    });
    auto weighted_degrees = getWeightedDegrees(G);

    std::vector<RowEl3*> H(G.n+1);
    // since we do not know elimination in advance, we need to store the whole matrix
    // know partial order: terminal nodes after non-terminals
    for(int i = 0; i < G.n; ++i){
        for(auto e : G.adjList[i]) {
            H[i] = new RowEl3(e.first, H[i], nullptr, e.second/weighted_degrees[i]);
            H[e.first] = new RowEl3(e.first, H[e.first], H[i], e.second/weighted_degrees[e.first]);
            H[i]->rev = H[e.first];
        }
    }

    // have vector of pointers to self loop and to sink node
    H[G.n] = new RowEl3(G.n, H[G.n], nullptr, 1);
    std::vector<double> selfLoops(G.n, 0);

    int count = 0;
    // schur complement: replace j -> i -> k by j -> k
    for(int i : inv_K){
        std::cout << "removing node " << count << std::endl;
        ++count;
        // for all j --> i (using the fact that the unweighted graph is symmetric)
        auto ij = H[i];
        while(ij){
            if(ij->col == G.n || ij->col == i) continue;
            double ji = ij->rev->v;
            double new_weight = 0;
            // and for i --> k
            auto ik = H[i];
            while(ik) {
                if(ik->col == i) continue;
                // compute new edge weight, account for moving from x --> z --> y, x --> z --> z --> --> y, ...
                // (1-alpha) * w_G(x,z) * w_G(z,y) * sum_(t=0)^∞ ((1 - a) * w_G(z,z))^t
                // sum_(t=0)^∞ ((1 - a) * w_G(z,z))^t = 1/(a * w_G(z,z) - w_G(z,z) + 1) when abs((1 - a) * w_G(z,z))<1
                // do not remove next line! If H[x][x] does not exist, it will be inserted and invalidate our loop iterators
                double new_edge_weight = (1-alpha) * ji * ik->v * (1.0/(alpha * selfLoops[i] - selfLoops[i] + 1.0));
                new_weight += new_edge_weight;
                // create new edge x --> y
                // H[ij->col][ik->col]. += new_edge_weight;
                // also remember to add self loop;
            }
            // add missing weight (accounts for moving from x to z and restarting at z)
            H[j.first][G.n] += (ji - new_weight);
            H[j.first].erase(i);
        }
    }

    // transform to output format
    std::vector<int> map_G_to_K(G.n+1, -1);
    for(int i = 0; i < K.size(); ++i){
        map_G_to_K[K[i]] = i;
    }
    map_G_to_K[G.n] = K.size();
    DiGraph sparsifier;
    for(int i = 0; i < K.size(); ++i){
        for(auto e : H[K[i]]){
            sparsifier.add_edge(i, map_G_to_K[e.first], e.second);
        }
    }

    return {sparsifier, K};
}*/