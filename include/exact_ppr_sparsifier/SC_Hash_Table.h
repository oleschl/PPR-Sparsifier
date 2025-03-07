#ifndef UNTITLED2_SC_HASH_TABLE_H
#define UNTITLED2_SC_HASH_TABLE_H

#include "graph.h"

/**
 * @brief Computes a Personalized PageRank (PPR) sparsifier.
 *
 * This implementation leverages the fact that construction of PPR sparsifiers
 * can be reduced to computing and modifying the schur complement. It computes
 * the Schur complement by applying partial gaussian elimination and represents
 * the (elimination) graph using a list of hash tables.
 *
 * @param G The undirected weighted input graph.
 * @param K The set of terminal nodes
 * @param alpha The restart probability
 * @return A PPR-sparsifier of G.
 */
namespace SC_Hash_Table {
    DiGraph constructPPRSparsifier2(const GEdge &G, std::vector<int>& K, double alpha);
    DiGraph constructPPRSparsifier(const GEdge &G, std::vector<int>& K, double alpha, std::string& order);
}

#endif //UNTITLED2_SC_HASH_TABLE_H
