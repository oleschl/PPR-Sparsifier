#ifndef GRAPH_IO_H
#define GRAPH_IO_H

#include "graph.h"
#include <string>

// Function declarations
//void parseFile(const std::string& filename, Graph& g, std::vector<int>& indices);
//void parseFile2(const std::string& filename, Graph& g, std::vector<int>& indices, int k);
void parseFile3(const std::string& filename, GEdge& g, std::vector<int>& indices, int k);

#endif  // GRAPH_IO_H