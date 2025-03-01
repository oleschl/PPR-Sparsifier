#ifndef UNTITLED2_PPR_MATRIX_INV_H
#define UNTITLED2_PPR_MATRIX_INV_H

#include "graph.h"

/**
 * @brief Computes a Personalized PageRank (PPR) sparsifier.
 *
 * Based on alternative algorithm description in:
 * Andrea Vattani, Deepayan Chakrabarti, and Maxim Gurevich. 2011.
 * Preserving personalized pagerank in subgraphs.
 *
 * @param G The undirected weighted input graph.
 * @param K The set of terminal nodes
 * @param alpha The restart probability
 * @return A PPR-sparsifier of G.
 */
namespace PPR_Matrix_Inv {
    DiGraph constructPPRSparsifier(const GEdge& G, const std::vector<int>& K, double alpha);
}

#endif //UNTITLED2_PPR_MATRIX_INV_H
