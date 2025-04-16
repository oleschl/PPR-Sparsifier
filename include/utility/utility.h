#ifndef UNTITLED2_UTILITY_H
#define UNTITLED2_UTILITY_H

#include "../graph.h"
#include <vector>

std::vector<int> getRandomTerminals(int n, int k, int seed);
double mean(const std::vector<double>& data);
double variance_sample(const std::vector<double>& data);
double compute_percentile(std::vector<double> data, double p);
std::vector<int> getNonTerminals(int n, const std::vector<int>& K);

#endif //UNTITLED2_UTILITY_H
