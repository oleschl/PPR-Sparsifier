#include "utility/utility.h"

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