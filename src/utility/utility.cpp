#include <random>
#include <unordered_set>
#include "utility/utility.h"

std::vector<int> getRandomTerminals(int n, int k, int seed) {
    if (k > n) throw std::invalid_argument("k cannot be greater than n");

    std::mt19937 rng(seed);
    std::vector<int> nodes(n);
    std::iota(nodes.begin(), nodes.end(), 0);  // Fill with 0, 1, ..., n-1
    std::shuffle(nodes.begin(), nodes.end(), rng);

    std::vector<int> terminals(nodes.begin(), nodes.begin() + k);
    std::sort(terminals.begin(), terminals.end());

    return terminals;
}

// RNE sampling for m unique nodes
std::vector<int> get_RNE_terminals(const GEdge& G, int m, int seed) {

    int n = G.n;
    std::vector<std::vector<int>> adj(G.n, std::vector<int>());
    for(auto edge : G.edges){
        adj[edge.v].push_back(edge.u);
        adj[edge.u].push_back(edge.v);
    }
    std::unordered_set<int> sampled_nodes;

    std::mt19937 gen(seed);
    std::uniform_int_distribution<> node_dist(0, n - 1);

    while ((int)sampled_nodes.size() < m) {
        int u = node_dist(gen);

        // Skip nodes with no neighbors
        if (adj[u].empty()) continue;

        std::uniform_int_distribution<> neigh_dist(0, adj[u].size() - 1);
        int v = adj[u][neigh_dist(gen)];

        sampled_nodes.insert(u);
        sampled_nodes.insert(v);
    }

    std::vector<int> result(sampled_nodes.begin(), sampled_nodes.end());
    std::sort(result.begin(), result.end());
    return result;
}

double mean(const std::vector<double>& data) {
    if (data.empty()) return 0.0;
    double sum = std::accumulate(data.begin(), data.end(), 0.0);
    return sum / data.size();
}

// Computes sample variance
double sample_variance(const std::vector<double>& data) {
    if (data.size() < 2) return 0.0;
    double mu = mean(data);
    double var = 0.0;
    for (double x : data) {
        var += (x - mu) * (x - mu);
    }
    return var / (data.size() - 1);
}

double sample_stddev(const std::vector<double>& data) {
    return std::sqrt(sample_variance(data));
}

double compute_percentile(std::vector<double> data, double p) {
    std::sort(data.begin(), data.end());

    double idx = p * (data.size() - 1);
    auto lower = static_cast<size_t>(std::floor(idx));
    auto upper = static_cast<size_t>(std::ceil(idx));

    if (lower == upper) {
        return data[lower];
    } else {
        // Linear interpolation between nearest ranks
        return data[lower] + (idx - lower) * (data[upper] - data[lower]);
    }
}

std::vector<int> getNonTerminals(int n, const std::vector<int>& K) {
    std::vector<int> non_terminals;
    non_terminals.reserve(n - K.size());

    size_t k_index = 0;
    for (int i = 0; i < n; ++i) {
        if (k_index < K.size() && i == K[k_index]) {
            ++k_index;
        } else {
            non_terminals.push_back(i);
        }
    }

    return non_terminals;
}