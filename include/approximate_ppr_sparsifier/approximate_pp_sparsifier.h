#ifndef PPR_SPARSIFIER_APPROXIMATE_PP_SPARSIFIER_H
#define PPR_SPARSIFIER_APPROXIMATE_PP_SPARSIFIER_H

#include "graph.h"

namespace ApproximateSparsifier {
    DiGraph constructPPRSparsifier(const GEdge &G, std::vector<int> &K, double alpha, int split, int merge);
}

#endif //PPR_SPARSIFIER_APPROXIMATE_PP_SPARSIFIER_H
