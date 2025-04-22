#ifndef PPR_SPARSIFIER_GRAPH_UTILITY_H
#define PPR_SPARSIFIER_GRAPH_UTILITY_H

#include "graph.h"
#include <string>

void parseEdgeList(const std::string& filename, GEdge& G);
void writeEdgeList(const std::string& filename, const GEdge& G);
GEdge chimera(int n, int k, bool weighted);
GEdge sachdevaStar(int l, int k);
GEdge inducedSubgraph(const GEdge& G, const std::vector<int>& terminals);

#endif //PPR_SPARSIFIER_GRAPH_UTILITY_H
