# PPR-Sparsifier: Graph Compression Techniques for Personalized PageRank

This repository contains the source code for the implementations and experiments conducted as part of my Master's thesis.  
We study vertex sparsifiers that preserve **Personalized PageRank (PPR)** values.

The repository provides a framework for comparing the **runtime performance** and **quality** of different implementation approaches.

It includes five implementations for computing **exact** sparsifiers:

- **Two implementations** are based on the pseudocode and descriptions from the paper  
  *Andrea Vattani, Deepayan Chakrabarti, and Maxim Gurevich.  
  [Preserving Personalized PageRank in Subgraphs](https://icml.cc/Conferences/2011/papers/434_icmlpaper.pdf).  
  In *Proceedings of the 28th International Conference on Machine Learning (ICML)*, pages 793–800, 2011.*
- The **other three implementations** are novel contributions developed as part of this thesis. These approaches are similar but differ primarily in how they compute the **Schur complement**.

In addition, the repository includes two implementations for computing **approximate** sparsifiers.

## How to Run / Requirements

### Prerequisites

This project embeds Julia into C++, so you will need a working Julia installation along with a few specific Julia packages. Julia is used primarily for solving Laplacian systems and generating Chimera graphs. The project also depends on the Armadillo linear algebra library.

- [Julia](https://julialang.org/) (tested with version ≥ 1.10.4)
- Julia packages:
    - `SparseArrays`
    - `Laplacians`
- [Armadillo](https://arma.sourceforge.net/) (C++ linear algebra library)

Note: You may need to set the JULIA_DIR environment variable (for locating Julia headers and libraries), and a similar configuration may be required for Armadillo.

### Building with CMake

Ensure you have CMake and a C++ compiler installed (e.g., GCC or Clang). Then build the project as follows:
```
mkdir build
cd build
cmake ..
make
```
This will create the PPR-Sparsifier executable in the build/ directory.

### Running the Executable

After building the project, you can run it with:

```bash
./PPR-Sparsifier -m <mode> -a <algorithm> -v <version> -f <input.txt> -k <number_of_terminal_nodes> -x <alpha> -d 
```

### Command-Line Arguments

- `-m <mode>`: Choose the mode. Currently supported:
    - exact_from_file: Exact sparsifier using a graph from file
    - approximate_from_file: Approximate sparsifier using a graph from file
    - exact_chimera: Exact sparsifier on a randomly generated Chimera graph
    - approximate_chimera: Approximate sparsifier on a Chimera graph
- `-a <algorithm>`: (Exact modes only) Choose one of the five algorithms:
  - `node_removal`, `ppr_matrix_inv`, `block_elimination`, `sc_adj_list`, `sc_hash_table`
- `v <version>`: Version of the algorithm (depends on mode):
  - For exact: choose an elimination order - `random`, `static`, `dynamic`
  - For approximate: choose a sampling strategy - 'random_clique', 'elim_star'
- `-f <input.txt>`: Path to the input graph file in edge list format (for _from_file modes). 
- `-k <number_of_terminal_nodes>`: Number of terminal nodes to retain, selected uniformly at random from the input graph.
- `-x <alpha>`: The restart probability for Personalized PageRank (typically between 0.1 and 0.3).
- `-d`: If provided, the tool will compute PPR values in both the input and sparsifier graphs, report whether differences are detected, and calculate error norms.

Additional Parameters

For approximate modes:
- `-s <split>`: Splits each edge into s multi-edges.
- `-g <merge>`: Limits the number of multi-edges between any two nodes and merges them when number in exceeded.
- `-r <runs>`: Repeats the randomized sparsification r times and reports mean performance.

For chimera modes:
- `-n <number_of_nodes>`: Number of nodes in the generated Chimera graph.
- `-r <runs>`: (used again) Number of Chimera instances to generate; reports median performance over these runs.

### Other Experiments

This repository provides a pipeline for evaluating both the runtime performance and quality of PPR sparsifiers using error norms.
For completeness, the `experiments/` folder contains additional code snippets used in other experiments conducted as part of my Master's thesis.