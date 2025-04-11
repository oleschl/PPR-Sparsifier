#include "utility/graph_utility.h"
#include <string>
#include <vector>
#include <fstream>
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
    jl_value_t* arg1 = jl_box_int64(n);
    jl_value_t* arg2 = jl_box_int64(k);
    jl_value_t* args[5] = {arg1, arg2};
    jl_eval_string("using Laplacians");
    jl_module_t* lap_mod = (jl_module_t*)jl_eval_string("Laplacians");
    jl_function_t* func;
    if(!weighted) {
        func = jl_get_function(lap_mod, "uni_chimera");
    } else {
        func = jl_get_function(lap_mod, "wted_chimera");
    }
    jl_value_t* ijv_obj = jl_call(func, args, 2);

    if (jl_exception_occurred())
        jl_printf(jl_stderr_stream(), "Error: %s\n", jl_typeof_str(jl_exception_occurred()));

    //jl_value_t* nnz_field = jl_get_field(ijv_obj, "nnz");
    jl_value_t* col_field = jl_get_field(ijv_obj, "colptr");
    jl_value_t* row_field = jl_get_field(ijv_obj, "rowval");
    jl_value_t* v_field = jl_get_field(ijv_obj, "nzval");

    int64_t* rows = (int64_t*) jl_array_data(row_field);
    int64_t* cols = (int64_t*) jl_array_data(col_field);
    double* v = (double*) jl_array_data(v_field);

    std::cout << v[1] << std::endl;

    std::vector<Edge> edges;
    int m = 0;
    for(int i = 0; i < n; i++){
        for(int j = cols[i]-1; j < cols[i+1]-2; ++j){
            if(i < rows[j]-1){
                edges.emplace_back(i, rows[j]-1, v[j]);
                ++m;
            }
        }
    }

    return {n, m, edges};
}

void parseEdgeList(const std::string& filename, GEdge& G) {
    std::ifstream infile(filename);
    if (!infile) {
        throw std::runtime_error("Error opening file: " + filename);
    }

    // read in edges
    int u, v;
    G.n = 0, G.m = 0;
    while (infile >> u >> v) {
        ++G.m;
        if(G.n < std::max(u, v)+1){
            G.n = std::max(u,v)+1;
        }
        G.addEdge(u, v, 1.0);
    }

    infile.close();
}
