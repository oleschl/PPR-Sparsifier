#include <vector>

struct SparseMatrix{
    int nnz;
    std::vector<double> lnz;
    std::vector<int> xlnz;
    std::vector<int> xnzsub;
    std::vector<int> nzsub;
};

SparseMatrix symbolic_factorization(int n, int k, std::vector<int> &adj, std::vector<int> &xadj, std::vector<int> &perm, std::vector<int> &inv_perm){

    int nz_beg = 0;
    int nz_end = -1;

    std::vector<int> merge_link(n, -1);
    std::vector<int> marker(n, -1);
    std::vector<int> reach_link(n);

    std::vector<int> xlnz(n+1);
    xlnz[0] = 0;
    std::vector<int> xnzsub(n+1);
    xnzsub[0] = 0;
    std::vector<int> nzsub;

    for(int i = 0; i < k; ++i){
        //std::cout << "processing row " << i << " or " << perm[i] << std::endl;
        // number of non-zeros of i
        int i_nz = 0;
        int merge_i = merge_link[i];
        bool merge_flag = false;
        if(merge_i != -1){
            marker[i] = marker[merge_i];
        }
        xnzsub[i] = nz_end;
        int node = perm[i];

        // find adj[i]-s_i-1
        reach_link[i] = n;
        for(int j = xadj[node]; j < xadj[node+1]; ++j){
            int neighbor = inv_perm[adj[j]];
            if(neighbor < i) continue;

            int cur = i;
            while(reach_link[cur] < neighbor) {
                cur = reach_link[cur];
            }
            ++i_nz;
            reach_link[neighbor] = reach_link[cur];
            reach_link[cur] = neighbor;

            if(marker[neighbor] != marker[i]) merge_flag = true;
        }

        int lmax = 0;

        if(merge_flag || merge_i == -1 || merge_link[merge_i] != -1){
            // iterate through each column j that effects i
            int j = merge_link[i];
            while(j != -1){
                int j_nz = xlnz[j+1] - xlnz[j] - 1;
                if(lmax <= j_nz) {
                    lmax = j_nz;
                    xnzsub[i] = xnzsub[j] + 1;
                }

                int cur = i;
                for(int jj = xnzsub[j]+1; jj <= xnzsub[j] + j_nz; ++jj){
                    int neighbor = nzsub[jj];
                    while(reach_link[cur] < neighbor){
                        cur = reach_link[cur];
                    }
                    if(reach_link[cur] == neighbor) continue;
                    ++i_nz;
                    reach_link[neighbor] = reach_link[cur];
                    reach_link[cur] = neighbor;
                    cur = neighbor;
                }
                j = merge_link[j];
            }

            if(i_nz != lmax){
                bool exit_early = false;
                if(nz_beg <= nz_end){
                    int jj = reach_link[i];
                    int j_start = nz_beg;
                    while(j_start <= nz_end && j_start < jj){
                        ++j_start;
                    }

                    if(j_start == jj){
                        bool update_end = true;
                        xnzsub[i] = j_start;

                        for(int jjj = j_start; jjj <= nz_end; ++jjj){
                            if(nzsub[jjj] != jj){
                                // goto 1200
                                update_end = false;
                                break;
                            }
                            jj = reach_link[jj];
                            if(jj > n){
                                // goto 1400
                                exit_early = true;
                                break;
                            }
                        }
                        if(update_end)
                            nz_end = j_start - 1;
                    }
                }

                if(!exit_early) {
                    // copy structure from reach_link to data structure
                    nz_beg = nz_end + 1;
                    nz_end += i_nz;
                    int cur = i;
                    for (int jj = nz_beg; jj <= nz_end; ++jj) {
                        cur = reach_link[cur];
                        nzsub.push_back(cur);
                        marker[cur] = i;
                    }
                    xnzsub[i] = nz_beg;
                    marker[i] = i;
                }
            }
        } else {
            xnzsub[i] = xnzsub[merge_i] + 1;
            i_nz = xlnz[merge_i+1] - xlnz[merge_i] - 1;
        }

        if(i_nz > 0){
            int j = nzsub[xnzsub[i]];
            //std::cout << "first non-zero " << j << std::endl;
            merge_link[i] = merge_link[j];
            merge_link[j] = i;
        }
        xlnz[i+1] = xlnz[i] + i_nz;
    }
    int nnz = xlnz[k];

    return {nnz, std::vector<double> (nnz), xlnz, xnzsub, nzsub};
}