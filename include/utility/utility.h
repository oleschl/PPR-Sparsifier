#ifndef UNTITLED2_UTILITY_H
#define UNTITLED2_UTILITY_H

#include "../graph.h"
#include <vector>

std::vector<int> getRandomTerminals(int n, int k, int seed);
std::vector<int> getNonTerminals(int n, const std::vector<int>& K);

#endif //UNTITLED2_UTILITY_H
