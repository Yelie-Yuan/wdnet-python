#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>

/**
 * Calculate node strengths for a directed network.
 *
 * This function calculates the out-strengths and in-strengths of nodes in a
 * directed network based on the provided edgelist and edgeweight.
 *
 * @param snode The list of source nodes.
 * @param tnode The list of target nodes.
 * @param edgeweight The list of edge weights corresponding to the edges in edgelist.
 * @param nnode The number of nodes in the network.
 * 
 * @return A pair of vectors containing out-strengths and in-strengths for each node.
 */
std::pair<std::vector<float>, std::vector<float>> node_strength(
    const std::vector<int>& snode,
    const std::vector<int>& tnode,
    const std::vector<float>& edgeweight,
    int nnode)
{
    std::vector<float> ins(nnode, 0.0);
    std::vector<float> outs(nnode, 0.0);

    int src, tar;
    float weight;
    for (decltype(snode.size()) i = 0; i < snode.size(); i++)
    {
        src = snode[i];
        tar = tnode[i];
        weight = edgeweight[i];
        outs[src] += weight;
        ins[tar] += weight;
    }

    return std::make_pair(outs, ins);
}
