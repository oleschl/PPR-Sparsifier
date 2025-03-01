#include <iostream>
#include <iomanip>
#include <chrono>
#include <functional>

#include "graph.h"
#include "graph_io.h"
#include "./utility/pagerank.h"
#include "exact_ppr_sparsifier/NodeRemoval.h"
#include "exact_ppr_sparsifier/PPR_Matrix_Inv.h"
#include "exact_ppr_sparsifier/SC_Block_Elimination.h"
#include "exact_ppr_sparsifier/SC_Adj_List.h"
#include "exact_ppr_sparsifier/SC_Hash_Table.h"

enum Implementation {
    NodeRemovalImplementation,
    PPR_MatrixInvImplementation,
    SC_AdjListImplementation,
    SC_HashTableImplementation,
    SC_BlockEliminationImplementation
};

void printAdjMatrix(const std::vector<std::vector<std::pair<int, double>>>& adjList, int n) {
    // Initialize adjacency matrix with zeros
    std::vector<std::vector<double>> adjMatrix(n, std::vector<double>(n, 0.0));

    // Fill adjacency matrix using adjacency list
    for (int u = 0; u < n; ++u) {
        for (const auto& [v, weight] : adjList[u]) {
            adjMatrix[u][v] = weight;
        }
    }

    // Print the adjacency matrix
    std::cout << std::fixed << std::setprecision(2);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cout << adjMatrix[i][j] << " ";
        }
        std::cout << "\n";
    }
}

void printEdgesAboveThreshold(const std::vector<std::vector<std::pair<int, double>>>& adjList, double threshold) {
    std::cout << std::fixed << std::setprecision(2);
    for (int u = 0; u < adjList.size(); ++u) {
        for (const auto& [v, weight] : adjList[u]) {
            if (weight > threshold) {
                std::cout << "Edge (" << u << " -> " << v << ") with weight " << weight << "\n";
            }
        }
    }
}

template <typename Func, typename... Args>
auto measureTime(Func func, Args&&... args) {
    auto start = std::chrono::high_resolution_clock::now();
    auto result = func(std::forward<Args>(args)...);
    auto end = std::chrono::high_resolution_clock::now();

    double elapsed_time = std::chrono::duration<double, std::milli>(end - start).count(); // Time in milliseconds
    return std::make_pair(elapsed_time, result);
}

int main() {

    Implementation version = NodeRemovalImplementation;

    GEdge G;
    std::vector<int> K;
    parseFile3("in.txt", G, K, 100);

    std::pair<double, DiGraph> result;
    switch (version) {
        case Implementation::NodeRemovalImplementation:
            result = measureTime(NodeRemoval::constructPPRSparsifier, G, K, 0.15);
            break;
        case Implementation::PPR_MatrixInvImplementation:
            result = measureTime(PPR_Matrix_Inv::constructPPRSparsifier, G, K, 0.15);
            break;
        case Implementation::SC_AdjListImplementation:
            result = measureTime(SC_Adj_List::constructPPRSparsifier, G, K, 0.15);
            break;
        case Implementation::SC_HashTableImplementation:
            result = measureTime(SC_Hash_Table::constructPPRSparsifier, G, K, 0.15);
            break;
        case Implementation::SC_BlockEliminationImplementation:
            result = measureTime(SC_BlockElimination::constructPPRSparsifier, G, K, 0.15);
            break;
    }

    std::cout << "Function 1 time: " << result.first << " ms\n";
    Sparsifier sp = {result.second, K};
    comparePPVs(G, sp, 0.85);

    //auto sparsifier = getApproximateSparsifier(G, K, 10, 10, 0.85);
    // comparePPVs(G, sparsifier, 0.85);

   // auto newG = approx_chol(G, K, 2, 2);

    //printAdjMatrix(sparsifier.H.adjList, sparsifier.H.n);

    //printEdgesAboveThreshold(sparsifier.H.adjList, 0.1);

   //sparsifier.H.printAdjList();

    //newG.printAdjList();

    return 0;
}


/*

#include <bits/stdc++.h>
using namespace std;

static float alpha = 0.15;

// implementation of nodeRemoval algorithm
// invS is the set of nodes that will be removed!
void nodeRemoval(vector<unordered_map<int, double>>& G, const vector<int>& invS){

    // keep track of removed nodes
    vector<bool> removed(G.size(), false);

    // add sink node
    int sink = G.size();
    G.push_back({{sink, 1.0}});

    for(auto z : invS){

        cout << "removing node " << z << endl;
        removed[z] = true;

        // for all x --> z
        for(int x = 0; x < G.size(); ++x){
            if(removed[x] || !G[x].contains(z)) continue;

            double xz = G[x][z];
            double newEdgeSum = 0.0;

            // and z --> y
            for(auto y : G[z]){
                if(removed[y.first]) continue;

                // compute new edge weight, account for moving from x --> z --> y, x --> z --> z --> --> y, ...
                // (1-alpha) * w_G(x,z) * w_G(z,y) * sum_(t=0)^∞ ((1 - a) * w_G(z,z))^t
                // sum_(t=0)^∞ ((1 - a) * w_G(z,z))^t = 1/(a * w_G(z,z) - w_G(z,z) + 1) when abs((1 - a) * w_G(z,z))<1
                auto newEdgeWeight = (1-alpha) * xz * y.second * (1.0/(alpha * G[z][z] - G[z][z] + 1.0));
                newEdgeSum += newEdgeWeight;

                // create new edge x --> y
                G[x][y.first] += newEdgeWeight;
            }

            // add missing weight (accounts for moving from x to z and restarting at z)
            G[x][sink] += (xz - newEdgeSum);

            // can remove edge x --> z here or need to test if z has been removed above
            //G[x].erase(z);
        }
    }
}

// reads a weighted, directed graph G and subset S of nodes (the nodes that will be in the compressed graph) from file
// first line of file should contain a comma separated list of node ids for S
// rest of file should be edge list representation of G: 'sourceId targetId weight'
// node ids must be in range 0, ..., |V|-1
// out degrees must be normalized: degree of node v must sum up to 1 for all v
void readGraph(const string& filename, vector<unordered_map<int, double>>& G, set<int>& S) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    string firstLine;
    getline(file, firstLine);

    istringstream iss(firstLine);

    int value;
    while (iss >> value) {
        S.insert(value);
        if (iss.peek() == ',')
            iss.ignore();
    }

    int u, v;
    double w;
    while (file >> u >> v >> w) {
        while(G.size() <= u){
            G.emplace_back();
        }
        G[u][v] = w;
    }

    file.close();
}

// writes graph G as edge list to file
void saveGraph(vector<unordered_map<int, double>>& G, vector<int>& S, const string& filename){
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    auto tmp = set(S.begin(), S.end());

    for (int u = 0; u < G.size(); ++u) {
        if(tmp.contains(u)) continue;
        for (const auto& neighbor : G[u]) {
            int v = neighbor.first;
            double w = neighbor.second;
            if (!tmp.contains((v))){
                file << u << " " << v << " " << w << endl;
            }
        }
    }

    file.close();
}

int main(){

    vector<unordered_map<int, double>> G = {};
    set<int> S = {};

    readGraph("in.txt", G, S);

    vector<int> invS = {};
    vector<int> degree(G.size());

    for(int i = 0; i < G.size(); ++i){
        if(S.contains(i)) continue;
        invS.push_back(i);
        degree[i] = G[i].size();
    }

    // sort nodes by degree (min degree heuristic)
    std::sort(invS.begin(), invS.end(), [&degree](int v, int u) {
        return degree[v] < degree[u];
    });

    auto start = std::chrono::high_resolution_clock::now();
    nodeRemoval(G, invS);
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed_time = std::chrono::duration<double, std::milli>(end - start).count(); // Time in milliseconds
    std::cout << elapsed_time << std::endl;

    saveGraph(G, invS, "out.txt");
}

*/