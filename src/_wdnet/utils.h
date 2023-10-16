#pragma once
#include <vector>

void node_strength_cpp(
    const std::vector<int> &snode,
    const std::vector<int> &tnode,
    const std::vector<double> &edgeweight,
    const int nnode,
    std::vector<double> &outs,
    std::vector<double> &ins);