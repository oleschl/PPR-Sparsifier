#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <iostream>

struct Edge {
    int u, v;
    double weight;
};

struct GEdge {
    int n, m;
    std::vector<Edge> edges;

    void addEdge(int u, int v, double w) {
        edges.push_back({u, v, w});
    }
};

/*
struct Graph {
    int n, m;
    bool isDirected = false;
    std::vector<std::vector<std::pair<int, double>>> adjList;

    Graph(): n(0), m(0), adjList() {};
    Graph(int n): n(n), adjList(n), m(0) {};
    Graph(int n, bool isDirected): n(n), adjList(n), isDirected(isDirected), m(0) {};

    void printAdjList() {
        for (int i = 0; i < n; ++i) {
            std::cout << "Node " << i << ": ";
            for (const auto& edge : adjList[i]) {
                // edge.first is the target node, edge.second is the weight
                std::cout << "(" << edge.first << ", " << edge.second << ") ";
            }
            std::cout << std::endl;
        }
    }

    void add_edge(int u, int v, double w) {
        int uu = std::min(u, v);
        int vv = std::max(u, v);

        while(adjList.size() <= vv){
            adjList.emplace_back();
            ++n;
        }
        adjList[uu].emplace_back(vv, w);
        ++m;
    }
};*/

struct DiGraph {
    int n, m;
    std::vector<std::vector<std::pair<int, double>>> adjList;

    DiGraph() = default;
    DiGraph(int n): n(n), m(0), adjList(n) {};

    void printAdjList() {
        for (int i = 0; i < n; ++i) {
            std::cout << "Node " << i << ": ";
            for (const auto& edge : adjList[i]) {
                // edge.first is the target node, edge.second is the weight
                std::cout << "(" << edge.first << ", " << edge.second << ") ";
            }
            std::cout << std::endl;
        }
    }

    void add_edge(int u, int v, double w) {
        adjList[u].emplace_back(v, w);
        ++m;
    }
};

struct Sparsifier {
    DiGraph H;
    std::vector<int> mapHToG;
};

#endif