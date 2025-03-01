#include "graph_io.h"
#include <algorithm>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <random>
/*
void parseFile(const std::string& filename, Graph& G, std::vector<int>& indices) {
    std::ifstream infile(filename);
    if (!infile) {
        throw std::runtime_error("Error opening file: " + filename);
    }

    // read list of node indices from first line of file
    std::string firstLine;
    std::getline(infile, firstLine);
    std::istringstream iss(firstLine);

    int index;
    while (iss >> index) {
        indices.push_back(index);
    }

    std::sort(indices.begin(), indices.end());

    // read in graph
    int u, v;
    int n = 0;
    int m = 0;
    double w;
    while (infile >> u >> v >> w) {
        G.add_edge(u, v, w);
    }

    infile.close();
}

void parseFile2(const std::string& filename, Graph& G, std::vector<int>& indices, int k) {
    std::ifstream infile(filename);
    if (!infile) {
        throw std::runtime_error("Error opening file: " + filename);
    }

    // read in graph
    int u, v;
    while (infile >> u >> v) {
        G.add_edge(u, v, 1.0);
    }

    std::random_device rd;
    std::mt19937 g(12345);

    std::vector<int> K(G.n);
    std::iota(K.begin(), K.end(), 0);
    std::shuffle(K.begin(), K.end(), g);

    for(int i = 0; i < k; ++i) {
        indices.push_back(K[i]);
    }

    std::sort(indices.begin(), indices.end());

    infile.close();
}*/

void parseFile3(const std::string& filename, GEdge& G, std::vector<int>& indices, int k) {
    std::ifstream infile(filename);
    if (!infile) {
        throw std::runtime_error("Error opening file: " + filename);
    }

    // read in edges
    int u, v;
    G.n = 0, G.m = 0;
    while (infile >> u >> v) {
        ++G.m;
        if(G.n < std::max(u, v)+1){
            G.n = std::max(u,v)+1;
        }
        G.addEdge(u, v, 1.0);
    }

    std::random_device rd;
    std::mt19937 g(12345);

    std::vector<int> K(G.n);
    std::iota(K.begin(), K.end(), 0);
    std::shuffle(K.begin(), K.end(), g);

    for(int i = 0; i < k; ++i) {
        indices.push_back(K[i]);
    }

    std::sort(indices.begin(), indices.end());

    infile.close();
}
