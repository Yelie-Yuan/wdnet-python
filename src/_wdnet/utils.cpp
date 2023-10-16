#include <iostream>
#include <algorithm>

/**
 * @brief Compute the out-strengths and in-strengths for each node in the network.
 *
 * @param snode The list of source nodes.
 * @param tnode The list of target nodes.
 * @param edgeweight The list of edge weights corresponding to the edges in edgelist.
 * @param nnode The number of nodes in the network.
 * @param outs The vector to store the out-strengths.
 * @param ins The vector to store the in-strengths.
 * 
 * @return A pair of vectors containing out-strengths and in-strengths for each node.
 */
void node_strength_cpp(
    const std::vector<int>& snode,
    const std::vector<int>& tnode,
    const std::vector<double>& edgeweight,
    const int nnode,
    std::vector<double>& outs,
    std::vector<double>& ins)
{

    int src, tar;
    double weight;
    for (decltype(snode.size()) i = 0; i < snode.size(); i++)
    {
        src = snode[i];
        tar = tnode[i];
        weight = edgeweight[i];
        outs[src] += weight;
        ins[tar] += weight;
    }
}
