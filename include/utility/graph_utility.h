#ifndef PPR_SPARSIFIER_GRAPH_UTILITY_H
#define PPR_SPARSIFIER_GRAPH_UTILITY_H

#include "graph.h"
#include <string>

void parseEdgeList(const std::string& filename, GEdge& G);
GEdge chimera(int n, int k, bool weighted);
GEdge sachdevaStar(int l, int k);

#endif //PPR_SPARSIFIER_GRAPH_UTILITY_H
