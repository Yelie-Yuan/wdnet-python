# Import the necessary Cython and NumPy modules
import numpy as np
cimport numpy as cnp
cnp.import_array()

# Define the Cython function with appropriate type declarations
cpdef node_strength_cy(
        cnp.ndarray[cnp.int_t, ndim=2] edgelist, 
        cnp.ndarray[cnp.float64_t, ndim=1] edgeweight,
        int nnode,
        int nedge
    ):
    cdef cnp.ndarray[cnp.float64_t, ndim=1] outs = np.zeros(nnode, dtype=np.float64)
    cdef cnp.ndarray[cnp.float64_t, ndim=1] ins = np.zeros(nnode, dtype=np.float64)
    cdef cnp.int_t src, tar
    cdef cnp.float64_t weight
    cdef int i

    for i in range(nedge):
        src = edgelist[i, 0]
        tar = edgelist[i, 1]
        weight = edgeweight[i]
        outs[src] += weight
        ins[tar] += weight

    return outs, ins
