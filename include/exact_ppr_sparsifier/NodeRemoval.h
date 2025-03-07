#ifndef UNTITLED2_NODEREMOVAL_H
#define UNTITLED2_NODEREMOVAL_H

#include "graph.h"

/**
 * @brief Computes a Personalized PageRank (PPR) sparsifier.
 *
 * Based on pseudo code presented in:
 * Andrea Vattani, Deepayan Chakrabarti, and Maxim Gurevich. 2011.
 * Preserving personalized pagerank in subgraphs.
 *
 * @param G The undirected weighted input graph.
 * @param K The set of terminal nodes
 * @param alpha The restart probability
 * @return A PPR-sparsifier of G.
 */
namespace NodeRemoval {
    DiGraph constructPPRSparsifier2(const GEdge &G, std::vector<int>& K, double alpha);
    DiGraph constructPPRSparsifier(const GEdge &G, std::vector<int>& K, double alpha, std::string &order);
}

#endif //UNTITLED2_NODEREMOVAL_H
