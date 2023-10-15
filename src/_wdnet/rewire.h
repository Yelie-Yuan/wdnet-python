#pragma once

#include <vector>
#include <cstdint>

void dprewire_directed_cpp(
    const int iteration, 
    const int nattempts, 
    std::vector<int>& tnode,
    const std::vector<double>& sout,
    const std::vector<double>& sin,
    std::vector<double>& tout,
    std::vector<double>& tin,
    const std::vector<int>& index_s,
    std::vector<int>& index_t,
    const std::vector<std::vector<double>>& eta,
    const bool history,
    std::vector<std::vector<int>>& rewire_history,
    const uint32_t random_seed,
    std::vector<double>& outout,
    std::vector<double>& outin,
    std::vector<double>& inout,
    std::vector<double>& inin);
