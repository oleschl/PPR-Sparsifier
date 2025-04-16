#include <vector>
#include <unordered_map>
#include <map>
#include <algorithm>
#include "exact_ppr_sparsifier/SC_Hash_Table.h"
#include "utility/utility.h"
#include "utility/orderings.h"
#include "utility/min_degree_pq.h"

// TODO: Refactor remove linked lists version (maybe test performance once)
// test if it is better to precompute dynamic min degree or not
// test maps vs unordered maps

// this is an implementation for computing exact pagerank sparsifiers for undirected graphs
// using exact schur complement/gaussian elimination
// almost all operations are on a symmetric matrix and therefore should half the work
// this version uses a map for storing and modifying the triangular matrix
// TODO: refactor: make code clean!
// TODO: test on larger graph
// TODO: start version with linked list instead of list of maps: done for static version, also for adaptive but might not be clean

// static min degree version
// again assume that K is sorted
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
            //std::cout << "eliminating node " << i << " with id " << perm[i] << std::endl;
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


// TODO: do version with dynamic min degree, one with precomputing order, one without

// this versions stores the whole matrix even though the input is symmetric
// this because we dont know the order of eliminations
// one option is to split into two steps: analyze and factorize
Sparsifier exact_schur_complement2(GEdge& G, std::vector<int> K, double alpha){

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

    return {sparsifier, K};
}

Sparsifier exact_schur_complement3(GEdge& G, std::vector<int> K, double alpha){

    // get V\K, the set of non-terminals
    auto inv_K = getNonTerminals(G.n, K);

    auto [perm, invPerm] = getDynamicMinDegOrdering(G, K, inv_K);

    // create graph structure for elimination
    std::vector<std::unordered_map<int, double>> H(G.n);
    std::vector<double> weighted_degree(G.n,0);
    std::vector<double> initial_weighted_degrees(G.n, 0);
    for(auto edge : G.edges){
        if(invPerm[edge.u] < invPerm[edge.v]){
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
    for(int i = 0; i < G.n; ++i){
        double edge_weight = (alpha/(1-alpha)) * initial_weighted_degrees[i];
        H[i][G.n] = edge_weight;
        weighted_degree[i] += edge_weight;
        initial_weighted_degrees[i] += edge_weight;
    }

    // schur complement / elimination
    for(int i = 0; i < inv_K.size(); ++i){
        //std::cout << "eliminating node " << i << " with id " << perm[i] << std::endl;
        for(auto e1 = H[i].begin(); e1 != H[i].end(); ++e1){
            for(auto e2 = std::next(e1); e2 != H[i].end(); ++e2){
                auto newEdgeWeight = e1->second * e2->second / weighted_degree[i];
                if(e2->first > e1->first) {
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
    std::vector<double> final_weighted_degree(K.size()+1, 0);
    for(int i = 0; i < K.size(); ++i){
        for(auto e : H[inv_K.size() + i]){
            final_weighted_degree[i] += e.second;
            final_weighted_degree[e.first-inv_K.size()] += e.second;
        }
    }

    DiGraph sparsifier(K.size()+1);
    for(int i = 0; i < K.size(); ++i){
        int node = inv_K.size()+i;
        for(auto e : H[node]){
            if(e.first != G.n){
                sparsifier.add_edge(i, e.first-inv_K.size(), e.second);
                sparsifier.add_edge(e.first-inv_K.size(), i, e.second);
            } else {
                // remove initial weight of lifted node
                double new_weight = e.second - alpha * initial_weighted_degrees[node];
                sparsifier.add_edge(i, K.size(), new_weight);
            }
        }
        // add self loops
        sparsifier.add_edge(i, i, initial_weighted_degrees[node] - final_weighted_degree[i]);
    }

    return {sparsifier, K};
}
/*

struct RowEl{
    int col;
    RowEl* next;
    double v;

    RowEl(int col, RowEl* next, double v) : col(col), next(next), v(v) {}
};

std::vector<double> compressionRow;

void compressRow(RowEl* row) {
    auto e = row;
    std::vector<int> indices;
    while(e){
        if(compressionRow[e->col] == 0){
            indices.push_back(e->col);
        }
        compressionRow[e->col] += e->v;
        e = e->next;
    }

    if(indices.empty()){
        return;
    }
    // create new edges
    e = row;
    int i = 0;
    for(; i < indices.size()-1; ++i){
        e->col = indices[i];
        e->v = compressionRow[indices[i]];
        e = e->next;
    }
    e->col = indices[i];
    e->v = compressionRow[indices[i]];
    auto next_e = e->next;
    e->next = nullptr;
    while(next_e){
        auto tmp = next_e->next;
        delete(next_e);
        next_e = tmp;

    }
    // reset compressionRow
    for(auto i : indices){
        compressionRow[i] = 0;
    }
}

// version with static min degree and linked list data structure (instead of maps)
Sparsifier exact_schur_complement3(Graph& G, std::vector<int> K, double alpha){

    compressionRow.resize(G.n+1, 0);
    // get V\K, the set of non-terminals
    auto inv_K = getNonTerminals(G.n, K);
    // compute initial degrees
    auto degrees = getDegrees(G);
    // sort inv_k by degree (number of adjacent nodes), used for static min degree ordering
    // could use bucket sort here since degree is in [0, n]
    std::sort(inv_K.begin(), inv_K.end(), [&degrees](int x, int y) {
        return degrees[x] < degrees[y];
    });
    // create permutation vector
    std::vector<int> permutation(G.n+1);
    for(int i = 0; i < inv_K.size(); ++i) {
        permutation[inv_K[i]] = i;
    }
    for(int i = 0; i < K.size(); ++i){
        permutation[K[i]] = inv_K.size() + i;
    }
    permutation[G.n] = G.n;

    // create graph structure for elimination
    std::vector<RowEl*> H(G.n);
    std::vector<double> weighted_degree(G.n,0);
    std::vector<double> initial_weighted_degrees(G.n, 0);
    // until here code is identical to exact_schur_complement1
    for(int i = 0; i < G.n; ++i){
        for(auto e : G.adjList[i]) {
            if(permutation[i] < permutation[e.first]) {
                // or map to permutation?
                // could use a vector and insert
                H[i] = new RowEl(e.first, H[i], e.second);
                weighted_degree[i] += e.second;
            } else {
                H[e.first] = new RowEl(i, H[e.first], e.second);
                weighted_degree[e.first] += e.second;
            }
            initial_weighted_degrees[i] += e.second;
            initial_weighted_degrees[e.first] += e.second;
        }
    }
    // perform elimination on D-(1-alpha)A by lifting with additional sink node
    // instead of scaling A, do elimination on (alpha/(1-alpha))D-A
    for(int i = 0; i < G.n; ++i){
        double v = (alpha/(1-alpha)) * initial_weighted_degrees[i];
        H[i] = new RowEl(G.n, H[i], v);
        weighted_degree[i] += v;
        initial_weighted_degrees[i] += v;
    }

    int count = 0;
    // schur complement computation of undirected graph
    for(int i : inv_K){
        std::cout << "removing node " << count << std::endl;
        ++count;
        auto e1 = H[i];
        // each elimination simply adds new edges and H becomes a multi-graph
        // to reduce work, multi-edges are merged before elimination
        // does not work, too much memory
        //compressRow(e1);
        while(e1){
            auto e2 = e1->next;
            while(e2){
                auto newEdgeWeight = e1->v * e2->v / weighted_degree[i];
                // could remove if here if we would store in ascending order of permutation index
                // this could be useful as it would allow us to directly modify compression row here
                // and not create new objects for each edge
                // however it is difficult to keep row sorted during eliminations because of fill ins
                if(permutation[e2->col] > permutation[e1->col]) {
                    H[e1->col] = new RowEl(e2->col, H[e1->col], newEdgeWeight);
                    weighted_degree[e1->col] += newEdgeWeight;
                } else {
                    H[e2->col] = new RowEl(e1->col, H[e2->col], newEdgeWeight);
                    weighted_degree[e2->col] += newEdgeWeight;
                }
                e2 = e2->next;
            }
            e1 = e1->next;
        }
        e1 = H[i];
        while(e1) {
            compressRow(H[e1->col]);
            e1 = e1->next;
        }
    }
    // compute new degree for each node
    // diff to initial degrees corresponds to self loops in graph
    std::vector<double> final_weighted_degree(G.n + 1, 0);
    for(auto i : K){
        auto e = H[i];
        compressRow(e);
        while(e) {
            final_weighted_degree[i] += e->v;
            final_weighted_degree[e->col] += e->v;
            if(e->col == G.n) {
                e->v -= alpha * initial_weighted_degrees[i];
            }
            e = e->next;
        }
        H[i] = new RowEl(i, H[i], initial_weighted_degrees[i] - final_weighted_degree[i]);
    }

    // transform to sparsifier
    std::vector<int> mapGToK(G.n, -1);
    for(int i = 0; i < K.size(); ++i){
        mapGToK[K[i]] = i;
    }
    DiGraph sparsifier;
    for(int i = 0; i < K.size(); ++i){
        auto e = H[K[i]];
        while(e){
            if(e->col != G.n){
                if(mapGToK[e->col] != i){
                    sparsifier.add_edge(i, mapGToK[e->col], e->v);
                    sparsifier.add_edge(mapGToK[e->col], i, e->v);
                } else {
                    sparsifier.add_edge(i, i, e->v);
                }
            } else {
                sparsifier.add_edge(i, K.size(), e->v);
            }
            e = e->next;
        }
    }

    return {sparsifier, K};
}

int compressRowAndEliminate(RowEl* row, int eliminationVertex) {
    auto e = row;
    std::vector<int> indices;
    while(e){
        if(e->col != eliminationVertex) {
            if (compressionRow[e->col] == 0) {
                indices.push_back(e->col);
            }
            compressionRow[e->col] += e->v;
        }
        e = e->next;
    }

    if(indices.empty()){
        return 0;
    }
    // create new edges
    e = row;
    int ii = 0;
    for(; ii+1 < indices.size(); ++ii){
        e->col = indices[ii];
        e->v = compressionRow[indices[ii]];
        e = e->next;
    }
    e->col = indices[ii];
    e->v = compressionRow[indices[ii]];
    //auto next_e = e->next;
    e->next = nullptr;

    while(next_e){
        auto tmp = next_e->next;
        delete(next_e);
        next_e = tmp;

    }
    // reset compressionRow
    for(auto i : indices){
        compressionRow[i] = 0;
    }
    return indices.size();
}

void compressTerminalRow(RowEl* row, std::vector<int>& map_G_to_K){
    auto e = row;
    std::vector<int> indices;
    while(e){
        if(map_G_to_K[e->col] != -1) {
            if (compressionRow[e->col] == 0) {
                indices.push_back(e->col);
            }
            compressionRow[e->col] += e->v;
        }
        e = e->next;
    }

    // create new edges
    e = row;
    int i = 0;
    for(; i < indices.size()-1; ++i){
        e->col = indices[i];
        e->v = compressionRow[indices[i]];
        e = e->next;
    }
    e->col = indices[i];
    e->v = compressionRow[indices[i]];
    e->next = nullptr;
    // reset compressionRow
    for(auto i : indices){
        compressionRow[i] = 0;
    }
}

// version with dynamic min degree and linked list data structure (instead of map)
Sparsifier exact_schur_complement4(Graph& G, std::vector<int> K, double alpha) {

    compressionRow.resize(G.n+1, 0);
    // get V\K, the set of non-terminals
    auto inv_K = getNonTerminals(G.n, K);
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
    // compute initial degrees
    auto degrees = getDegrees(G);
    // create min degree pq
    std::vector<std::pair<int, int>> degreeInvK(inv_K.size());
    for(int i = 0; i < inv_K.size(); ++i){
        degreeInvK[i] = {i, degrees[inv_K[i]]};
    }
    MinDegreePQ pq(G.n, inv_K.size());

    // next create graph H as list of pointers to head of linked lists
    std::vector<RowEl*> H(G.n+1);
    std::vector<double> weighted_degree(G.n+1,0);
    // since we do not know elimination in advance, we need to store the whole
    // actually know partial order: terminal nodes after non-terminals
    for(int i = 0; i < G.n; ++i){
        for(auto e : G.adjList[i]) {
            H[i] = new RowEl(e.first, H[i], e.second);
            weighted_degree[i] += e.second;
            H[e.first] = new RowEl(i, H[e.first], e.second);
            weighted_degree[e.first] += e.second;
        }
    }

    // perform elimination on D-(1-alpha)A by lifting with additional sink node
    // instead of scaling A, do elimination on (alpha/(1-alpha))D-A
    for(int i = 0; i < G.n; ++i){
        double v = (alpha/(1-alpha)) * weighted_degree[i];
        H[i] = new RowEl(G.n, H[i], v);
        H[G.n] = new RowEl(i, H[G.n], v);
        weighted_degree[i] += v;
    }

    std::vector<double> initial_weighted_degrees = weighted_degree;

    int count = 0;
    // schur complement computation of undirected graph
    for(int i = 0; i < inv_K.size(); ++i){
        std::cout << "removing node " << count << std::endl;
        ++count;
        // get vertex with minimum degree
        int eliminationVertex = inv_K[pq.pop()];
        auto e1 = H[eliminationVertex];
        // create clique between all pairs of adjacent nodes
        while(e1){
            auto e2 = e1->next;
            while(e2){
                auto newEdgeWeight = e1->v * e2->v / weighted_degree[eliminationVertex];
                H[e1->col] = new RowEl(e2->col, H[e1->col], newEdgeWeight);
                weighted_degree[e1->col] += newEdgeWeight;
                H[e2->col] = new RowEl(e1->col, H[e2->col], newEdgeWeight);
                weighted_degree[e2->col] += newEdgeWeight;
                e2 = e2->next;
            }
            e1 = e1->next;
        }

        // now compress row and remove elimination vertex and update degrees
        e1 = H[eliminationVertex];
        while(e1) {
            // update degree
            int degree = compressRowAndEliminate(H[e1->col], eliminationVertex);
            if(map_G_to_K[e1->col] == -1){
                //int degree = compressRowAndEliminate(H[e1->col], eliminationVertex);
                pq.update(map_G_to_inv_K[e1->col], degree);
            }
            weighted_degree[e1->col] -= e1->v;
            e1= e1->next;
        }

        // adaptive version more complex
        // to determine degree we need to insert elements at correct positions
        // however linked list only allow efficient insertion in front position
        // only option appears to be to compress all rows here??? (except row in K)
    }

    // compute new degree for each node
    // diff to initial degrees corresponds to self loops in graph
    for(auto i : K){
        auto e = H[i];
        // dont need degree here
        compressTerminalRow(e,map_G_to_K);
        while(e) {
            if(e->col == G.n) {
                e->v -= alpha * initial_weighted_degrees[i];
            }
            e = e->next;
        }
        H[i] = new RowEl(i, H[i], initial_weighted_degrees[i] - weighted_degree[i]);
    }

    // transform to output format
    DiGraph sparsifier;

    for(int i = 0; i < K.size(); ++i){
        auto e = H[K[i]];
        while(e){
            if(map_G_to_K[e->col] != -1){
                sparsifier.add_edge(i, map_G_to_K[e->col], e->v);
            }
            e = e->next;
        }
    }

    return {sparsifier, K};
}
*/