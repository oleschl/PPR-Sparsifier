#ifndef UNTITLED2_UTILITY_H
#define UNTITLED2_UTILITY_H

#include "../graph.h"
#include <vector>

std::vector<int> getRandomTerminals(int n, int k, int seed);
std::vector<int> get_RNE_terminals(const GEdge& G, int m, int seed);
double mean(const std::vector<double>& data);
double sample_variance(const std::vector<double>& data);
double sample_stddev(const std::vector<double>& data);
double compute_percentile(std::vector<double> data, double p);
std::vector<int> getNonTerminals(int n, const std::vector<int>& K);

#endif //UNTITLED2_UTILITY_H
