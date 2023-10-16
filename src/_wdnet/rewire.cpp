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
 * @return double Correlation coefficient.
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
 * @param tnode Target nodes; its order will be modified in the
 * rewiring process.
 * @param sout Out-strength of source nodes.
 * @param sin In-strength of source nodes.
 * @param tout Out-strength of target nodes.
 * @param tin In-strength of target nodes.
 * @param index_s Indices of source nodes. The indices are used to
 * locate the values in eta.
 * @param index_t Indices of target nodes.
 * @param eta Target network structure.
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
    int n = 0, count = 0, i = 0;
    int e1, e2, s1, s2, t1, t2;

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

    for (n = 0; n < iteration; ++n)
    {
        for (i = 0; i < nattempts; ++i)
        {
            e1 = std::floor(dis_nedge(gen));
            e2 = std::floor(dis_nedge(gen));

            while (e1 == e2)
            {
                e2 = std::floor(dis_nedge(gen));
            }

            if (history)
            {
                rewire_history[count][0] = e1;
                rewire_history[count][1] = e2;
            }

            s1 = index_s[e1];
            s2 = index_s[e2];
            t1 = index_t[e1];
            t2 = index_t[e2];

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

/**
 * @brief Compute the mean of a vector.
 *
 * @param x Vector.
 * @return double Mean value.
 */
double mean(const std::vector<double> &x)
{
    double sum = std::accumulate(x.begin(), x.end(), 0.0);
    return sum / x.size();
}

/**
 * @brief Compute the correlation between two vectors.
 *
 * @param x Vector 1.
 * @param y Vector 2.
 * @return double Correlation coefficient.
 */
double corr(const std::vector<double> &x,
            const std::vector<double> &y)
{
    double mean_x = mean(x);
    double mean_y = mean(y);
    double numerator = 0, sum_x2 = 0.0, sum_y2 = 0.0;
    int n = x.size();

    for (int i = 0; i < n; ++i)
    {
        numerator += (x[i] - mean_x) * (y[i] - mean_y);
        sum_x2 += std::pow(x[i] - mean_x, 2);
        sum_y2 += std::pow(y[i] - mean_x, 2);
    }
    return numerator / std::sqrt(sum_x2 * sum_y2);
}

/**
 * @brief Rewire an undirected network towards a given eta.
 *
 * @param iteration Number of iterations. Each iteration consists of
 * nattempts rewiring attempts.
 * @param nattempts Number of attempts per iteration.
 * @param node1 Nodes in the first column.
 * @param node2 Nodes in the second column.
 * @param degree1 Degrees of nodes in the first column.
 * @param degree2 Degrees of nodes in the second column.
 * @param index1 Indices of nodes in the first column. The indices are
 * used to locate the values in eta.
 * @param index2 Indices of nodes in the second column.
 * @param eta Target network structure.
 * @param history Whether to record the rewiring history.
 * @param rewire_history Rewiring history.
 * @param random_seed Random seed.
 * @param rho Assortativity coefficient after each iteration.
 */
void dprewire_undirected_cpp(
    const int iteration,
    const int nattempts,
    std::vector<int> &node1,
    std::vector<int> &node2,
    std::vector<double> &degree1,
    std::vector<double> &degree2,
    std::vector<int> &index1,
    std::vector<int> &index2,
    const std::vector<std::vector<double>> &eta,
    const bool history,
    std::vector<std::vector<int>> &rewire_history,
    const uint32_t random_seed,
    std::vector<double> &rho)
{
    int nedge = node1.size();
    int n = 0, count = 0, i = 0;
    int e1, e2, s1, s2, t1, t2;
    double ratio;

    std::mt19937 gen(random_seed);
    std::uniform_real_distribution<> dis_nedge(0, nedge);
    std::uniform_real_distribution<> dis_unit(0, 1);

    for (n = 0; n < iteration; n++)
    {
        for (i = 0; i < nattempts; i++)
        {
            e1 = std::floor(dis_nedge(gen));
            e2 = std::floor(dis_nedge(gen));
            while (e1 == e2)
            {
                e2 = std::floor(dis_nedge(gen));
            }

            if (history)
            {
                rewire_history[count][0] = e1;
                rewire_history[count][1] = e2;
            }
            s1 = index1[e1];
            s2 = index1[e2];
            t1 = index2[e1];
            t2 = index2[e2];

            ratio = 1.0;
            if (dis_unit(gen) < 0.5)
            {
                if (eta[s1][t2] * eta[s2][t1] < eta[s1][t1] * eta[s2][t2])
                {
                    ratio = eta[s1][t2] * eta[s2][t1] / (eta[s1][t1] * eta[s2][t2]);
                }

                if (dis_unit(gen) <= ratio)
                {
                    std::swap(index2[e1], index2[e2]);
                    std::swap(node2[e1], node2[e2]);
                    std::swap(degree2[e1], degree2[e2]);
                    std::swap(degree1[e1 + nedge], degree1[e2 + nedge]);
                    if (history)
                    {
                        rewire_history[count][3] = 1;
                    }
                }
            }
            else
            {
                if (history)
                {
                    rewire_history[count][2] = 1;
                }
                if (eta[s1][s2] * eta[t1][t2] < eta[s1][t1] * eta[s2][t2])
                {
                    ratio = eta[s1][s2] * eta[t1][t2] / (eta[s1][t1] * eta[s2][t2]);
                }

                if (dis_unit(gen) <= ratio)
                {
                    std::swap(index2[e1], index1[e2]);
                    std::swap(node2[e1], node1[e2]);
                    std::swap(degree2[e1], degree1[e2]);
                    std::swap(degree1[e1 + nedge], degree2[e2 + nedge]);
                    if (history)
                    {
                        rewire_history[count][3] = 1;
                    }
                }
            }
            count++;
        }
        rho[n] = corr(degree1, degree2);
    }
}
