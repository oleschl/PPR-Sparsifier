#include <vector>
#include <algorithm>
#include <climits>
#include "graph.h"
#include "utility/min_degree_pq.h"


std::pair<std::vector<int>, std::vector<int>> getStaticMinDegOrdering(const GEdge &G, std::vector<int> &K, std::vector<int> &invK){
    // compute initial degrees
    std::vector<int> degrees(G.n, 0);
    for(auto edge : G.edges){
        ++degrees[edge.u];
        ++degrees[edge.v];
    }
    // sort invK by degree (number of adjacent nodes)
    // could use bucket sort here since degree is in [0, n]
    std::sort(invK.begin(), invK.end(), [&degrees](int x, int y) {
        return degrees[x] < degrees[y];
    });

    // create permutation vectors
    std::vector<int> perm(G.n);
    std::vector<int> invPerm(G.n);
    for(int i = 0; i < invK.size(); ++i) {
        perm[i] = invK[i];
        invPerm[invK[i]] = i;
    }
    // include nodes in K in arbitrary order
    for(int i = invK.size(); i < G.n; ++i) {
        perm[i] = K[i-invK.size()];
        invPerm[K[i-invK.size()]] = i;
    }

    return {perm, invPerm};
}

void eliminate(MinDegreePQ &pq, std::vector<int> &map_G_to_invK, int v, int tag, std::vector<int> &marker, std::vector<int> &needs_update, std::vector<int> &inv_perm, std::vector<int> &super_size, std::vector<int> &merge_parent, std::vector<int> &elim_next, std::vector<int> &xadj, std::vector<int> &adj, std::vector<int> &map_G_to_K){
    marker[v] = tag;

    int elmnt = 0;
    int rloc = xadj[v];
    int rlmt = xadj[v+1]-1;
    // iterate over adjacent nodes and find all super and eliminated nodes
    for(int i = xadj[v]; i < xadj[v+1]; ++i){
        int neighbor = adj[i];
        if(neighbor == 0) break;
        // why do i need this test??? to avoid merge super nodes!
        if(marker[neighbor] < tag){
            marker[neighbor] = tag;
            // if neighbor has not been eliminated
            if(inv_perm[neighbor] == 0) {
                adj[rloc] = neighbor;
                ++rloc;
            } else {
                elim_next[neighbor] = elmnt;
                elmnt = neighbor;
            }
        }
    }

    // build clique from super (already inserted above) and eliminated (done here) nodes and store in G.xadj[v]
    while(elmnt != 0) {
        adj[rlmt] = -elmnt;
        int j = xadj[elmnt];
        int jstop = xadj[elmnt+1];
        int node = adj[j];
        while (node != 0){
            if(node < 0) {
                j = xadj[-node];
                jstop = xadj[-node + 1];
            } else {
                // original code also test the next node with same degree but i dont know why??? & degnext[node] >= 0
                if(marker[node] < tag) {
                    marker[node] = tag;

                    while(rloc >= rlmt){
                        int link = -adj[rlmt];
                        rloc = xadj[link];
                        rlmt = xadj[link+1]-1;
                    }
                    adj[rloc] = node;
                    ++rloc;
                }
                ++j;
            }
            if(j >= jstop) break;
            node = adj[j];
        }
        elmnt = elim_next[elmnt];
    }

    if (rloc <= rlmt) adj[rloc] = 0;

    // TODO: somewhere here we also have to test if the node is in K or invK!
    // maybe enough to set needupdate[node] = 0???
    // find indistinguishable nodes and nodes where we have to update degrees
    int i = xadj[v];
    int istop = xadj[v+1];
    // iterate over reachable nodes
    int r_node = adj[i];
    while (r_node != 0){
        if(r_node < 0) {
            i = xadj[-r_node];
            istop = xadj[-r_node + 1];
        } else if(map_G_to_K[r_node-1] == -1){
            // i think this step is only necessary here because we do multiple elimination
            // and want to find an independent set of nodes to eliminate
            // if node is connected, it is not independent and we need to remove it
            // TODO: remove rnode

            // this step should only be applied to nodes in invK i think because here we modify the
            // adjacency list of the adjacent nodes and we dont have to update nodes in K
            // we somehow have to avoid updates in update function
            // if rnode in K continue (dont forget to increment i)

            int xqnbr = xadj[r_node];
            for(int j = xadj[r_node]; j < xadj[r_node+1]; ++j){
                int neighbor = adj[j];
                if(neighbor == 0) break;
                // count number of differnt neihbours
                if(marker[neighbor] < tag){
                    adj[xqnbr] = neighbor;
                    ++xqnbr;
                }
            }

            int nqnbrs = xqnbr - xadj[r_node];
            // r_node is after elimination equal to v -> merge
            if(nqnbrs <= 0){
                //std::cout << "merging node " << r_node << " with size: " << super_size[r_node] << std::endl;
                super_size[v] += super_size[r_node];
                super_size[r_node] = 0;
                merge_parent[r_node] = v;
                marker[r_node] = INT_MAX;
                // we need to remove node from degree structure; maybe in update
                pq.eliminate(map_G_to_invK[r_node - 1]);
            } else {
                needs_update[r_node] = nqnbrs + 1;
                adj[xqnbr] = v;
                ++xqnbr;
                if(xqnbr < xadj[r_node+1]){
                    adj[xqnbr] = 0;
                }
            }

            ++i;
        } else {
            ++i;
        }
        if(i >= istop) break;
        r_node = adj[i];
    }
}

int update(MinDegreePQ &pq, int tag, int min_deg, std::vector<int> & xadj, std::vector<int> &adj, int v, std::vector<int> &super_size, std::vector<int> &marker, std::vector<int> &needs_update, std::vector<int> &tmp_next, std::vector<int> &inv_perm, std::vector<int> &merge_parent, std::vector<int> &map_G_to_invK){

    int max_tag = tag + min_deg;

    int q2_head = 0;
    int qx_head = 0;
    int elim_size = 0;

    // iterate over all reachable nodes and divide them into group with two neighbors and group with more than two
    int i = xadj[v];
    int i_stop = xadj[v+1];
    int r_node = adj[i];
    while(r_node != 0){
        if(r_node < 0){
            i = xadj[-r_node];
            i_stop = xadj[-r_node + 1];
        } else {
            // still need to this for all nodes as nodes in K also contribute to elimsize
            if(super_size[r_node] != 0) {
                elim_size += super_size[r_node];
                marker[r_node] = max_tag;

                // we dont need this, so setting needupdate = 0 for nodes in K should be enough?
                // i think we dont need this test since we dont do multple eliminations but can use to avoid eliminations of nodes in K
                // actually we might want to merch nodes in K to reduce number of updates; we dont do updates since we dont eliminate them; does not work
                // not so easy because only want to merge a node of K with another node of K
                // because ohter node will be eliminated so maybe dont merch haha
                // i think its not worth the effort
                if(needs_update[r_node] > 0){

                    if(needs_update[r_node] != 2){
                        tmp_next[r_node] = qx_head;
                        qx_head = r_node;
                    } else {
                        tmp_next[r_node] = q2_head;
                        q2_head = r_node;
                    }
                    //tmp_next[r_node] = qx_head;
                    //qx_head = r_node;
                }
            }
            ++i;
        }
        if(i >= i_stop) break;
        r_node = adj[i];
    }

    // for each reachable node with two neighbors (enode form current elimination and some other)
    r_node = q2_head;
    while(r_node > 0){
        if(needs_update[r_node] > 0){
            ++tag;
            int deg = elim_size;

            // identify other adjacent node
            int j = xadj[r_node];
            int neighbor = adj[j];
            if(neighbor == v) neighbor = adj[j+1];

            // has not been eliminated
            if(inv_perm[neighbor] == 0){
                deg += super_size[neighbor];
            } else { // has been eliminated so need to look at elimination clique
                int k = xadj[neighbor];
                int k_stop = xadj[neighbor+1];
                int clique_node = adj[k];
                while(clique_node != 0){
                    if(clique_node < 0){
                        k = xadj[-clique_node];
                        k_stop = xadj[-clique_node + 1];
                    } else {

                        if(clique_node != r_node && super_size[clique_node] != 0) {

                            if (marker[clique_node] < tag) {
                                marker[clique_node] = tag;
                                deg += super_size[clique_node];
                            } else if (needs_update[clique_node] > 0) {
                                // r_node and clique_node are indistingusialbe and merge them
                                if (needs_update[clique_node] == 2) {
                                    super_size[r_node] += super_size[clique_node];
                                    super_size[clique_node] = 0;
                                    marker[clique_node] = INT_MAX;
                                    merge_parent[clique_node] = r_node;
                                    pq.eliminate(map_G_to_invK[clique_node - 1]);
                                    //std::cout << "merging node (update)" << clique_node << " with: " << r_node
                                              //<< std::endl;
                                } else {
                                    //std::cout << "delaying update node " << clique_node << " with" << r_node
                                    //        << std::endl;
                                }
                                // or enode is a subset of snode. in this case we can delay the degree update
                                // until enode has been eliminated
                                needs_update[clique_node] = 0;
                                // i have to see how if i still can use my pq for this; maybe enough to just not update?
                                // dont need to do this; is only important for removing it in elimination but we dont do this
                                //degprev[snode] = 0;
                                pq.update(map_G_to_invK[clique_node-1], INT_MAX-1);
                            }
                        }
                        ++k;
                    }
                    if(k >= k_stop) break;
                    clique_node = adj[k];
                }
            }
            needs_update[r_node] = 0;
            pq.update(map_G_to_invK[r_node-1], deg-super_size[r_node]);
        }
        r_node = tmp_next[r_node];
    }

    r_node = qx_head;
    while(r_node != 0){
        if(needs_update[r_node] > 0){
            ++tag;
            int deg = elim_size;

            for(int ii = xadj[r_node]; ii < xadj[r_node+1]; ++ii){
                int neighbor = adj[ii];

                if(neighbor == 0) break;
                if(marker[neighbor] < tag){
                    marker[neighbor] = tag;

                    if(inv_perm[neighbor] == 0){
                        deg += super_size[neighbor];
                    } else {
                        int j = xadj[neighbor];
                        int jstop = xadj[neighbor + 1];
                        int node = adj[j];
                        while(node != 0){
                            if(node < 0){
                                j = xadj[-node];
                                jstop = xadj[-node + 1];
                            } else {
                                if(marker[node] < tag){
                                    marker[node] = tag;
                                    deg += super_size[node];
                                }
                                ++j;
                            }
                            if(j >= jstop) break;
                            node = adj[j];
                        }
                    }
                }
            }
            needs_update[r_node] = 0;
            pq.update(map_G_to_invK[r_node-1], deg-super_size[r_node]);
        }
        r_node = tmp_next[r_node];
    }

    return max_tag;
}

void set_permutation(int n, int k, std::vector<int> &perm, std::vector<int> &inv_perm, std::vector<int> &merge_parent) {

    std::vector<int> merge_last_num(n, 0);
    for(int i = 0; i < n; ++i){
        if(merge_parent[i] == 0){
            merge_last_num[i] = inv_perm[i];
        }
    }

    for(int i = 1; i < n; ++i){
        int parent = merge_parent[i];
        if(parent > 0){
            int root = 0;
            while(parent > 0){
                root = parent;
                parent = merge_parent[parent];
            }

            int num = merge_last_num[root] + 1;
            inv_perm[i] = num;
            merge_last_num[root] = num;

            int j = i;
            while(j != root){
                parent = merge_parent[j];
                merge_parent[j] = root;
                j = parent;
            }
        }
    }

    for(int i = 1; i < n; ++i){
        if(inv_perm[i] != 0 ){
            perm[inv_perm[i]-1] = i-1;
        }
    }
    for(int i = 0; i < n; ++i){
        inv_perm[i] = inv_perm[i+1]-1;
    }
}


// if v is a positive value means v is adjacent, if v is negative means all of v neighbors are adjacent (pointer to clique)
// implementation without multiple elimination
std::pair<std::vector<int>, std::vector<int>> getDynamicMinDegOrdering(int n, int m, std::vector<int> &xadj, std::vector<int> &adj, std::vector<int>& K, std::vector<int> & inv_K){
    int size_inv_K = inv_K.size();
    // map node ids to range [1, 2, ..., n]
    std::vector<int> md_xadj(n+2);
    std::vector<int> md_adj(m*2);
    for(int i = 0; i < adj.size(); ++i){
        md_adj[i] = adj[i]+1;
    }
    for(int i = 0; i < xadj.size(); ++i){
        md_xadj[i+1] = xadj[i];
    }

    // either need to map nodes in inv_K to range 1, ..., size(inv_K) or create vector of size G.n+1
    std::vector<int> perm(inv_K.size());
    std::vector<int> inv_perm(n+1, 0);
    std::vector<int> marker(n+1, 0);
    std::vector<int> super_size(n+1, 1);
    std::vector<int> merge_parent(n+1, 0);
    std::vector<int> needs_update(n+1, 0);
    // linked list to keep track of temporary node relations
    std::vector<int> tmp_next(n+1, 0);
    // number of eliminations
    int num = 1;
    int tag = 0;

    MinDegreePQ pq(n, size_inv_K);

    std::vector<std::pair<int, int>> degreeInvK(size_inv_K);
    // create priority queue
    for(int i = 0; i < size_inv_K; ++i){
        int degree = md_xadj[inv_K[i]+2] - md_xadj[inv_K[i]+1];
        // eliminate isolated nodes; acutally we dont have any as all nodes are connected to source vertex
        if(degree == 0){
            marker[i+1] = INT_MAX;
            inv_perm[i+1] = num;
            ++num;
        } else {
            pq.add(i, degree);
        }
    }

    std::vector<int> map_G_to_K(n+1, -1);
    for(int i = 0; i < K.size(); ++i){
        map_G_to_K[K[i]] = i;
    }

    std::vector<int> map_G_to_inv_K(n+1);
    for(int i = 0; i < size_inv_K; ++i){
        map_G_to_inv_K[inv_K[i]] = i;
    }


    while(num <= size_inv_K){
        int eliminationVertex = inv_K[pq.pop()]+1;
        int mindeg = pq.get_min_deg();
        inv_perm[eliminationVertex] = num;
        //std::cout << "eliminating: " << eliminationVertex << " degree: " << mindeg << " num: " << num << std::endl;

        if(num + super_size[eliminationVertex] > size_inv_K){
            break;
        }
        ++tag;
        eliminate(pq, map_G_to_inv_K, eliminationVertex, tag, marker, needs_update, inv_perm, super_size, merge_parent, tmp_next, md_xadj, md_adj, map_G_to_K);
        num += super_size[eliminationVertex];
        if(num > size_inv_K) break;
        tag = update(pq, tag, mindeg, md_xadj, md_adj, eliminationVertex, super_size, marker, needs_update, tmp_next, inv_perm, merge_parent, map_G_to_inv_K);
    }

    set_permutation(n+1, K.size(), perm, inv_perm, merge_parent);

    for(int i = 0; i < K.size(); ++i){
        inv_perm[K[i]] = size_inv_K + i;
    }

    return {perm, inv_perm};
}

std::pair<std::vector<int>, std::vector<int>> getDynamicMinDegOrdering(const GEdge &G, std::vector<int> &K, std::vector<int> &invK){
    std::vector<int> xadj(G.n+1);
    std::vector<int> adj(2*G.m);
    for(auto edge : G.edges){
        ++xadj[edge.u+1];
        ++xadj[edge.v+1];
    }
    for(int i = 1; i <= G.n; ++i){
        xadj[i] += xadj[i-1];
    }
    std::vector<int> temp_pos(G.n, 0);
    for(auto edge : G.edges){
        adj[xadj[edge.u]+temp_pos[edge.u]] = edge.v;
        ++temp_pos[edge.u];
        adj[xadj[edge.v]+temp_pos[edge.v]] = edge.u;
        ++temp_pos[edge.v];
    }

    return getDynamicMinDegOrdering(G.n, G.m, xadj, adj, K, invK);
}
