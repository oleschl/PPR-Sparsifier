#include "utility/graph_utility.h"
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <julia.h>

GEdge sachdevaStar(int l, int k){
    int n = l * k + 1;
    int m = l * (k*(k-1)/2) + l;
    std::vector<Edge> edges;

    int num_vertices = 1;
    for(int i = 0; i < l; ++i){
        for(int j = 0; j < k; ++j){
            for(int jj = j+1; jj < k; ++jj){
                edges.emplace_back(num_vertices+j, num_vertices+jj, 1);
            }
        }
        edges.emplace_back(0, num_vertices, 1);
        num_vertices += k;
    }

    return {n, m, edges};
}

GEdge chimera(int n, int k, bool weighted) {
    jl_value_t* n_value = jl_box_int64(n);
    jl_value_t* k_value = jl_box_int64(k);
    jl_function_t* func;
    if(!weighted) {
        func = jl_get_function(jl_main_module, "uni_chimera");
    } else {
        func = jl_get_function(jl_main_module, "wted_chimera");
    }
    jl_value_t* ijv_obj =  jl_call2(func, n_value, k_value);

     if (jl_exception_occurred()) {
         jl_call2(jl_get_function(jl_base_module, "showerror"),
                  (jl_value_t*)jl_stderr_obj(),
                  jl_exception_occurred());
         jl_printf(jl_stderr_stream(), "\n");
     }

     if (ijv_obj == nullptr) {
         std::cerr << "Error: ijv_obj is null!" << std::endl;
     }

    jl_value_t* col_field = jl_get_field(ijv_obj, "colptr");
    jl_value_t* row_field = jl_get_field(ijv_obj, "rowval");
    jl_value_t* v_field = jl_get_field(ijv_obj, "nzval");

    int64_t* rows = (int64_t*) jl_array_data(row_field);
    int64_t* cols = (int64_t*) jl_array_data(col_field);
    double* v = (double*) jl_array_data(v_field);

    std::vector<Edge> edges;
    int m = 0;
    for(int i = 0; i < n; i++){
        for(int j = cols[i]-1; j < cols[i+1]-1; ++j){
            if(i < rows[j]-1){
                edges.emplace_back(i, rows[j]-1, v[j]);
                ++m;
            }
        }
    }

    return {n, m, edges};
}

GEdge inducedSubgraph(const GEdge& G, const std::vector<int>& terminals) {
    std::vector<int> inv(G.n, -1);
    for(int i = 0; i < terminals.size(); ++i){
        inv[terminals[i]] = i;
    }

    std::vector<Edge> edges;

    for (const auto& edge : G.edges) {
        if (inv[edge.u] != -1 && inv[edge.v] != -1) {
            edges.emplace_back(inv[edge.u], inv[edge.v], edge.weight);
        }
    }

    int n = terminals.size();
    int m = edges.size();

    return {n, m, std::move(edges)};
}

// parse weighted + unweighted edge list
void parseEdgeList(const std::string& filename, GEdge& G) {
    std::ifstream infile(filename);
    if (!infile) {
        throw std::runtime_error("Error opening file: " + filename);
    }
    G.n = 0, G.m = 0;
    // read in edges
    std::string line;
    int u, v;
    double w;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        if (!(iss >> u >> v)) continue;
        if (!(iss >> w)) {
            w = 1.0; // only two values -> unweighted
        }

        // currently we only support reading undirected graphs where each edge appears only once
        //if (u < v){
        G.n = std::max(G.n, std::max(u, v) + 1);
        ++G.m;
        G.addEdge(u, v, w);
        //}
    }

    infile.close();
}

void writeEdgeList(const std::string& filename, const GEdge& G) {
    std::ofstream outfile(filename);
    if (!outfile) {
        throw std::runtime_error("Error opening file for writing: " + filename);
    }

    for (const auto& e : G.edges) {
        outfile << e.u << " " << e.v << " " << e.weight << "\n";
    }

    outfile.close();
}

