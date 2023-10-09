# Import necessary libraries
import numpy as np
cimport numpy as cnp
from libc.math cimport sqrt
from libcpp cimport bool
from libcpp.algorithm cimport swap
cnp.import_array()


# Define the function in Cython
cpdef cnp.float64_t compute_correlation(
        cnp.ndarray[cnp.float64_t, ndim=1] x,
        cnp.ndarray[cnp.float64_t, ndim=1] y,
        const cnp.float64_t xsum,
        const cnp.float64_t ysum,
        const cnp.float64_t x2sum,
        const cnp.float64_t y2sum):
    """
    Compute the correlation between arraries x and y.

    Parameters
    ----------
    x : numpy.ndarray
        The first array.
    y : numpy.ndarray
        The second array.
    xsum : double
        The sum of the first array.
    ysum : double
        The sum of the second array.
    x2sum : double
        The sum of the square of the first array.
    y2sum : double
        The sum of the square of the second array.
    
    Returns
    -------
    double
        The correlation between x and y.
    """
    cdef:
        cnp.float64_t xysum = 0
        int n = x.shape[0]
        int i
        cnp.float64_t numerator
        cnp.float64_t denominator

    for i in range(n):
        xysum += x[i] * y[i]

    numerator = n * xysum - xsum * ysum
    denominator = sqrt((n * x2sum - xsum * xsum) * (n * y2sum - ysum * ysum))

    return numerator / denominator

def dprewire_directed_cy(
        const int iteration,
        const int nattempts,
        cnp.ndarray[cnp.int_t, ndim=1] tnode,
        cnp.ndarray[cnp.float64_t, ndim=1] sout,
        cnp.ndarray[cnp.float64_t, ndim=1] sin,
        cnp.ndarray[cnp.float64_t, ndim=1] tout,
        cnp.ndarray[cnp.float64_t, ndim=1] tin,
        cnp.ndarray[cnp.int_t, ndim=1] index_s,
        cnp.ndarray[cnp.int_t, ndim=1] index_t,
        cnp.ndarray[cnp.float64_t, ndim=2] eta,
        const bool history):
    """
    Rewire the directed network towards the target eta.

    Parameters
    ----------
    iteration : int
        The number of iterations.
    nattempts : int
        The number of attempts in each iteration.
    tnode : vector[int]
        The target node of the edges.
    sout : vector[double]
        The out-strength of the source nodes.
    sin : vector[double]
        The in-strength of the source nodes.
    tout : vector[double]
        The out-strength of the target nodes.
    tin : vector[double]
        The in-strength of the target nodes.
    index_s : vector[int]
        The index of the source nodes, for looking up elements in eta.
    index_t : vector[int]
        The index of the target nodes, for looking up elements in eta.
    eta : vector[vector[double]]
        The target eta; matrix like.
    history : bool
        Whether to record the rewire_history of rewiring attempts.
    """
    cdef:
        cnp.ndarray[cnp.float64_t, ndim=1] outout = np.zeros(iteration, dtype=np.float64)
        cnp.ndarray[cnp.float64_t, ndim=1] outin = np.zeros(iteration, dtype=np.float64)
        cnp.ndarray[cnp.float64_t, ndim=1] inout = np.zeros(iteration, dtype=np.float64)
        cnp.ndarray[cnp.float64_t, ndim=1] inin = np.zeros(iteration, dtype=np.float64)
        cnp.float64_t soutsum = 0, sinsum = 0, toutsum = 0, tinsum = 0
        cnp.float64_t sout2sum = 0, sin2sum = 0, tout2sum = 0, tin2sum = 0
        cnp.int_t nedge = tnode.shape[0]
        cnp.int_t e1, e2, count = 0
        cnp.int_t s1, s2, t1, t2, hist_row
        cnp.float64_t u, ratio
        cnp.int_t i, n

    soutsum = np.sum(sout)
    sinsum = np.sum(sin)
    toutsum = np.sum(tout)
    tinsum = np.sum(tin)
    
    for i in range(sout.shape[0]):
        sout2sum += sout[i] * sout[i]
        sin2sum += sin[i] * sin[i]
        tout2sum += tout[i] * tout[i]
        tin2sum += tin[i] * tin[i]
    if history:
        hist_row = iteration * nattempts
    else:
        hist_row = 1

    cdef int[:, :] rewire_history = np.zeros((hist_row, 3), dtype=np.intc)

    for n in range(iteration):
        for i in range(nattempts):
            e1 = int(np.floor(np.random.uniform(0, nedge)))
            e2 = int(np.floor(np.random.uniform(0, nedge)))
            while e1 == e2:
                e2 = int(np.floor(np.random.uniform(0, nedge)))

            if history:
                rewire_history[count][0] = e1
                rewire_history[count][1] = e2

            s1 = index_s[e1]
            s2 = index_s[e2]
            t1 = index_t[e1]
            t2 = index_t[e2]

            if eta[s1][t2] * eta[s2][t1] < eta[s1][t1] * eta[s2][t2]:
                ratio = eta[s1][t2] * eta[s2][t1] / (eta[s1][t1] * eta[s2][t2])
            else:
                ratio = 1

            u = np.random.uniform()
            if u <= ratio:
                swap(index_t[e1], index_t[e2])
                swap(tnode[e1], tnode[e2])
                swap(tout[e1], tout[e2])
                swap(tin[e1], tin[e2])

                if history:
                    rewire_history[count][2] = 1
            count += 1

        outout[n] = compute_correlation(sout, tout, soutsum, toutsum, sout2sum, tout2sum)
        outin[n] = compute_correlation(sout, tin, soutsum, tinsum, sout2sum, tin2sum)
        inout[n] = compute_correlation(sin, tout, sinsum, toutsum, sin2sum, tout2sum)
        inin[n] = compute_correlation(sin, tin, sinsum, tinsum, sin2sum, tin2sum)

    return tnode, outout, outin, inout, inin, rewire_history
