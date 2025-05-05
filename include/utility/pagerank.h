#ifndef PAGERANK_H
#define PAGERANK_H

# include "graph.h"

std::vector<std::vector<std::pair<int, double>>> constructWalkingMatrix(const DiGraph& G, double alpha);
std::vector<double> pageRank(const std::vector<std::vector<std::pair<int, double>>>& G, std::vector<double>& r, int n, double eps, int matIter);
std::vector<double> pageRank(const DiGraph& G, std::vector<double>& r, double alpha, double eps, int matIter);
std::vector<std::vector<double>> precomputePPRMatrix(GEdge& G, const std::vector<int>& K, double alpha);

double frobeniusInduced(GEdge& G, GEdge& induced, std::vector<int> mapHtoG, double alpha);

double comparePPVs(GEdge& G, Sparsifier& H, double alpha);
double comparePPVs(std::vector<std::vector<double>>& PPR_G, Sparsifier& H, double alpha);

#endif //PAGERANK_H
