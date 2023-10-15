#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <cstdint>

/**
 * @brief Compute the correlation between two vectors.
 * 
 * @param x Vector 1.
 * @param y Vector 2.
 * @param xsum Sum of vector 1.
 * @param ysum Sum of vector 2.
 * @param x2sum Sum of squares of vector 1.
 * @param y2sum Sum of squares of vector 2.
 * @return double 
 */
double compute_correlation(
    const std::vector<double> &x,
    const std::vector<double> &y,
    double xsum,
    double ysum,
    double x2sum,
    double y2sum)
{
    double xysum = 0;
    int n = x.size();

    for (int i = 0; i < n; ++i)
    {
        xysum += x[i] * y[i];
    }

    double numerator = n * xysum - xsum * ysum;
    double denominator = std::sqrt((n * x2sum - xsum * xsum) * (n * y2sum - ysum * ysum));

    return numerator / denominator;
}

/**
 * @brief Rewire a directed network towards a given eta.
 *
 * @param iteration Number of iterations. Each iteration consists of
 * nattempts attempts.
 * @param nattempts Number of attempts per iteration.
 * @param tnode Vector of target nodes; its order will be modified in
 * the rewiring process.
 * @param sout Vector of out-strength of source nodes.
 * @param sin Vector of in-strength of source nodes.
 * @param tout Vector of out-strength of target nodes.
 * @param tin Vector of in-strength of target nodes.
 * @param index_s Vector of indices of source nodes.
 * @param index_t Vector of indices of target nodes.
 * @param eta Matrix of eta values.
 * @param history Whether to record the rewiring history.
 * @param rewire_history Rewiring history.
 * @param outout Out-out assortativity coefficient after each
 * iteration.
 * @param outin Out-in assortativity coefficient after each iteration.
 * @param inout In-out assortativity coefficient after each iteration.
 * @param inin In-in assortativity coefficient after each iteration.
*/
void dprewire_directed_cpp(
    const int iteration,
    const int nattempts,
    std::vector<int> &tnode,
    const std::vector<double> &sout,
    const std::vector<double> &sin,
    std::vector<double> &tout,
    std::vector<double> &tin,
    const std::vector<int> &index_s,
    std::vector<int> &index_t,
    const std::vector<std::vector<double>> &eta,
    const bool history,
    std::vector<std::vector<int>> &rewire_history,
    const uint32_t random_seed,
    std::vector<double> &outout,
    std::vector<double> &outin,
    std::vector<double> &inout,
    std::vector<double> &inin)
{
    double soutsum = std::accumulate(sout.begin(), sout.end(), 0.0);
    double sinsum = std::accumulate(sin.begin(), sin.end(), 0.0);
    double toutsum = std::accumulate(tout.begin(), tout.end(), 0.0);
    double tinsum = std::accumulate(tin.begin(), tin.end(), 0.0);

    double sout2sum = 0, sin2sum = 0, tout2sum = 0, tin2sum = 0;
    double ratio;

    for (double val : sout)
        sout2sum += val * val;
    for (double val : sin)
        sin2sum += val * val;
    for (double val : tout)
        tout2sum += val * val;
    for (double val : tin)
        tin2sum += val * val;

    int nedge = tnode.size();

    // std::random_device rd;
    std::mt19937 gen(random_seed);
    std::uniform_real_distribution<> dis_nedge(0, nedge);
    std::uniform_real_distribution<> dis_unit(0, 1);

    int count = 0;
    for (int n = 0; n < iteration; ++n)
    {
        for (int i = 0; i < nattempts; ++i)
        {
            int e1 = std::floor(dis_nedge(gen));
            int e2 = std::floor(dis_nedge(gen));

            while (e1 == e2)
            {
                e2 = std::floor(dis_nedge(gen));
            }

            if (history)
            {
                rewire_history[count][0] = e1;
                rewire_history[count][1] = e2;
            }

            int s1 = index_s[e1];
            int s2 = index_s[e2];
            int t1 = index_t[e1];
            int t2 = index_t[e2];

            ratio = 1.0;
            if (eta[s1][t2] * eta[s2][t1] < eta[s1][t1] * eta[s2][t2])
            {
                ratio = eta[s1][t2] * eta[s2][t1] / (eta[s1][t1] * eta[s2][t2]);
            }

            if (dis_unit(gen) <= ratio)
            {
                std::swap(index_t[e1], index_t[e2]);
                std::swap(tnode[e1], tnode[e2]);
                std::swap(tout[e1], tout[e2]);
                std::swap(tin[e1], tin[e2]);

                if (history)
                {
                    rewire_history[count][2] = 1;
                }
            }

            ++count;
        }

        outout[n] = compute_correlation(sout, tout, soutsum, toutsum, sout2sum, tout2sum);
        outin[n] = compute_correlation(sout, tin, soutsum, tinsum, sout2sum, tin2sum);
        inout[n] = compute_correlation(sin, tout, sinsum, toutsum, sin2sum, tout2sum);
        inin[n] = compute_correlation(sin, tin, sinsum, tinsum, sin2sum, tin2sum);
    }
}
