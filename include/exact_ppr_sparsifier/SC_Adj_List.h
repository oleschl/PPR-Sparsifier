#ifndef UNTITLED2_SC_ADJ_LIST_H
#define UNTITLED2_SC_ADJ_LIST_H

#include "graph.h"

/**
 * @brief Computes a Personalized PageRank (PPR) sparsifier.
 *
 * This implementation leverages the fact that construction of PPR sparsifiers
 * can be reduced to computing and modifying the schur complement. It computes
 * the Schur complement by applying partial gaussian elimination and represents
 * the (elimination) graph using an adjacency list.
 *
 * @param G The undirected weighted input graph.
 * @param K The set of terminal nodes
 * @param alpha The restart probability
 * @return A PPR-sparsifier of G.
 */
namespace SC_Adj_List {
    DiGraph constructPPRSparsifier(const GEdge &G, std::vector<int>& K, double alpha);
}

#endif //UNTITLED2_SC_ADJ_LIST_H
