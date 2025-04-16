#include <random>
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

double mean(const std::vector<double>& data) {
    if (data.empty()) return 0.0;
    double sum = std::accumulate(data.begin(), data.end(), 0.0);
    return sum / data.size();
}

// Computes sample variance (denominator: N-1)
double variance_sample(const std::vector<double>& data) {
    if (data.size() < 2) return 0.0;
    double mu = mean(data);
    double var = 0.0;
    for (double x : data) {
        var += (x - mu) * (x - mu);
    }
    return var / (data.size() - 1);
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