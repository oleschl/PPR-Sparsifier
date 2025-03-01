#include <vector>
#include "graph.h"
#include "symbolic_factorization.cpp"

Sparsifier gaussian_elimination(Graph &G, std::vector<int> &K, double alpha){

    auto inv_K = getNonTerminals(G.n, K);

    // transform triangular part of symmetric graph to complete matrix
    Graph H(G);
    std::vector<double> diag(G.n, 0);
    for(int i = 0; i < G.n; ++i){
        for(auto e : G.adjList[i]){
            H.adjList[e.first].emplace_back(i, e.second);
            diag[i] += e.second;
            diag[e.first] += e.second;
        }
    }

    for(int i = 0; i < G.n; ++i){
        diag[i] += (alpha/(1-alpha)) * diag[i];
    }

    std::vector<double> initial_diag(diag);

    // create sparse matrix representation of input graph G
    std::vector<int> xadj(G.n+1);
    std::vector<int> adj(G.m*2);
    std::vector<double> weights(G.m*2);
    int count = 0;
    for(int i = 0; i < H.n; ++i){
        for(auto e : H.adjList[i]){
            adj[count] = e.first;
            weights[count] = -e.second;
            ++count;
        }
        xadj[i+1] = count;
    }

    // get minimum degree ordering
    auto [perm, inv_perm] = minimum_degree_ordering(G.n, G.m, xadj, adj, K, inv_K);
    // precompute fill-in structure
    auto M = symbolic_factorization(G.n, inv_K.size(), adj, xadj, perm, inv_perm);

    std::vector<int> link(G.n, -1);
    std::vector<int> first(G.n);
    std::vector<double> working_row(G.n, 0);

    // for each row in inv_K
    for(int i = 0; i < inv_K.size(); ++i){
        // eliminate e
        // TODO: adapt code
        // symbolic factorization maps node ids to elimination order
        // eg. 7619 (first node to eliminate) is mapped to 0 i think
        int e = perm[i];
        std::cout << "eliminating node " << i << " with id " << e << std::endl;

        // insert value of row[e] into tmp
        for(int j = xadj[e]; j < xadj[e+1]; ++j){
            if(inv_perm[adj[j]] > i){
                working_row[inv_perm[adj[j]]] = weights[j];
                //std::cout << "setting col " << inv_perm[adj[j]] << " with weight " << weights[j] << std::endl;
            }
        }
        // get rows that effect row
        int linked_row = link[i];
        while(linked_row != -1){
            // maybe error with diag change indexing!
            double factor = M.lnz[M.xlnz[linked_row]+first[linked_row]] / diag[perm[linked_row]];
            diag[e] -= M.lnz[M.xlnz[linked_row]+first[linked_row]] * factor;
            int c = 1;
            for(int k = M.xlnz[linked_row]+first[linked_row]+1; k < M.xlnz[linked_row+1]; ++k){
                int col = M.nzsub[M.xnzsub[linked_row]+first[linked_row]+c];
                if(M.lnz[k] == 0 ){
                    int abcefd = 0;
                }
                working_row[col] -= M.lnz[k] * factor;
                //std::cout << "setting col " << col << " new weight " << working_row[col] << std::endl;
                ++c;
            }

            int llinked_row = linked_row;
            linked_row = link[linked_row];

            ++first[llinked_row];
            if(first[llinked_row] < M.xlnz[llinked_row+1]-M.xlnz[llinked_row]){
                int new_linked_row = M.nzsub[M.xnzsub[llinked_row]+first[llinked_row]];
                int old_linked_row = link[new_linked_row];
                link[new_linked_row] = llinked_row;
                link[llinked_row] = old_linked_row;
            }
        }

        first[i] = 0;
        int new_linked_row = M.nzsub[M.xnzsub[i]];
        // std::cout << "linking with row " << new_linked_row << std::endl;
        int old_linked_row = link[new_linked_row];
        link[new_linked_row] = i;
        link[i] = old_linked_row;

        // store into matrix structure + reset tmp vector
        int c = 0;
        for(int j = M.xlnz[i]; j < M.xlnz[i+1]; ++j){
            int col = M.nzsub[M.xnzsub[i]+c];
            //std::cout << "using col " << col << std::endl;
            M.lnz[j] = working_row[col];
            //std::cout << "updating " << i << "  " << M.nzsub[M.xnzsub[i]+c] << " value " << working_row[col] << std::endl;
            working_row[col] = 0;
            if(c > 1){
                int asfd = 5;
            }
            ++c;
        }
    }

    // add code for applying elimination to nodes in K (triangular part)
    DiGraph sparsifier;
    for(int i = 0; i < K.size(); ++i){
        int e = K[i];

        // insert value of row[e] into tmp
        for(int j = xadj[e]; j < xadj[e+1]; ++j){
            if(inv_perm[adj[j]] > inv_K.size()+i){
                working_row[inv_perm[adj[j]]] = weights[j];
                //std::cout << "setting col " << inv_perm[adj[j]] << " with weight " << weights[j] << std::endl;
            }
        }

        int linked_row = link[inv_K.size()+i];
        while(linked_row != -1){
            // maybe error with diag change indexing!
            double factor = M.lnz[M.xlnz[linked_row]+first[linked_row]] / diag[perm[linked_row]];
            diag[e] -= M.lnz[M.xlnz[linked_row]+first[linked_row]] * factor;
            int c = 1;
            for(int k = M.xlnz[linked_row]+first[linked_row]+1; k < M.xlnz[linked_row+1]; ++k){
                int col = M.nzsub[M.xnzsub[linked_row]+first[linked_row]+c];
                if(M.lnz[k] == 0 ){
                    int abcefd = 0;
                }
                working_row[col] -= M.lnz[k] * factor;
                //std::cout << "setting col " << col << " new weight " << working_row[col] << std::endl;
                ++c;
            }

            int llinked_row = linked_row;
            linked_row = link[linked_row];

            ++first[llinked_row];
            if(first[llinked_row] < M.xlnz[llinked_row+1]-M.xlnz[llinked_row]){
                int new_linked_row = M.nzsub[M.xnzsub[llinked_row]+first[llinked_row]];
                int old_linked_row = link[new_linked_row];
                link[new_linked_row] = llinked_row;
                link[llinked_row] = old_linked_row;
            }
        }

        // store into matrix structure + reset tmp vector
        // will be dense anyways so just insert and reset whole temporary vector!
        for(int j = i+1; j < K.size(); ++j){
            int col = inv_K.size() + j;
            sparsifier.add_edge(i, j, -working_row[col]);
            sparsifier.add_edge(j, i, -working_row[col]);
            working_row[col] = 0;

        }
    }
    for(int i = 0; i < K.size(); ++i){
        double sum = 0;
        for(auto & j : sparsifier.adjList[i]){
            sum += j.second;
        }
        sparsifier.add_edge(i, i, initial_diag[K[i]]-diag[K[i]]);
        double sink_weight = (diag[K[i]]-sum) - alpha * initial_diag[K[i]];
        sparsifier.add_edge(i, K.size(), sink_weight);
    }
    // test result
    return {sparsifier, K};
}
