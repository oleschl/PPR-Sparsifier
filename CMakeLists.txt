cmake_minimum_required(VERSION 3.22.1)
project(PPR-Sparsifier)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-O0")
set(CMAKE_C_FLAGS "-O0")

find_package(Armadillo REQUIRED)

set(JULIA_DIR "/home/oles/.julia/juliaup/julia-1.10.4+0.x64.linux.gnu")

# Include Julia headers
include_directories(${CMAKE_SOURCE_DIR}/include "${JULIA_DIR}/include/julia")

# Link Julia library (Ensure libjulia.dll.a is used)
link_directories(${JULIA_DIR}/lib)

add_executable(PPR-Sparsifier
        src/main.cpp
        src/utility/pagerank.cpp
        src/utility/utility.cpp
        src/utility/orderings.cpp
        src/utility/min_degree_pq.cpp
        src/utility/graph_utility.cpp
        src/exact_ppr-sparsifier/NodeRemoval/NodeRemoval.cpp
        src/exact_ppr-sparsifier/PPR_Matrix_Inv/PPR_Matrix_Inv.cpp
        src/exact_ppr-sparsifier/SC_Adj_List/SC_AdjList.cpp
        src/exact_ppr-sparsifier/SC_Hash_Table/SC_Hash_Table.cpp
        src/exact_ppr-sparsifier/BlockElimination/BlockElimination.cpp
        src/approximate_pp-sparsifier/approx_pp-sparsifier.cpp
)

# Link with Julia
target_link_libraries(PPR-Sparsifier armadillo "${JULIA_DIR}/lib/libjulia.so")

target_compile_options(PPR-Sparsifier PRIVATE
        -Wall       # Enable most common warnings
        -Wextra     # Enable extra warnings
        -Wpedantic  # Enforce standard compliance
)

