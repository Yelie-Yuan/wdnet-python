#pragma once

#include <vector>

std::pair<std::vector<float>, std::vector<float>> node_strength(
    const std::vector<int>& snode,
    const std::vector<int>& tnode,
    const std::vector<float>& edgeweight,
    int nnode
);