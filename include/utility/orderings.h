#ifndef UNTITLED2_ORDERINGS_H
#define UNTITLED2_ORDERINGS_H

#include "../graph.h"

std::pair<std::vector<int>, std::vector<int>> getRandomOrdering(const GEdge &G, std::vector<int> &K, std::vector<int> &invK);
std::pair<std::vector<int>, std::vector<int>> getStaticMinDegOrdering(const GEdge &G, std::vector<int> &K, std::vector<int> &invK);
std::pair<std::vector<int>, std::vector<int>> getDynamicMinDegOrdering(int n, int m, std::vector<int> &xadj, std::vector<int> &adj, std::vector<int>& K, std::vector<int> & inv_K);
std::pair<std::vector<int>, std::vector<int>> getDynamicMinDegOrdering(const GEdge &G, std::vector<int> &K, std::vector<int> &invK);

#endif //UNTITLED2_ORDERINGS_H
