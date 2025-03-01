#include <vector>
#include <armadillo>

#include "exact_ppr_sparsifier/PPR_Matrix_Inv.h"
#include "utility/pagerank.h"

namespace PPR_Matrix_Inv {
    DiGraph constructPPRSparsifier(const GEdge& G, const std::vector<int>& K, double alpha) {
        int k = K.size();
        // transform input graph to directed graph and construct walking (transition) matrix W for computing PageRank
        DiGraph H;
        for (auto edge: G.edges) {
            H.add_edge(edge.u, edge.v, edge.weight);
            H.add_edge(edge.v, edge.u, edge.weight);
        }
        auto W_G = constructWalkingMatrix(H, 1-alpha);
        // construct PageRank matrix PPR
        arma::mat PPR(K.size() + 1, K.size() + 1, arma::fill::zeros);
        // set row of sink node
        PPR(K.size(), K.size()) = 1;
        // personalization vector r
        std::vector<double> r(G.n, 0);
        for (int i = 0; i < k; ++i) {
            // set personalization vector
            r[K[i]] = 1;
            auto ppr_G = pageRank(W_G, r, H.n, 1e-8, 1000);
            double cum_sum = 0;
            for (int j = 0; j < k; ++j) {
                PPR(i, j) = ppr_G[K[j]];
                cum_sum += ppr_G[K[j]];
            }
            // set column entry for sink node
            PPR(i, K.size()) = (1 - cum_sum);
            // reset personalization vector
            r[K[i]] = 0;
        }
        // solve W_T = (1/(1-alpha)) (I-alpha P_inv)
        arma::mat PPR_inv = arma::inv(PPR);
        arma::mat I = arma::eye<arma::mat>(K.size() + 1, K.size() + 1);
        arma::mat W = (1 / (1 - alpha)) * (I - alpha * PPR_inv);
        // transform new walking matrix W to graph
        DiGraph sparsifier(K.size()+1);
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < k + 1; ++j) {
                sparsifier.add_edge(i, j, W(i, j));
            }
        }

        return sparsifier;
    }
}
