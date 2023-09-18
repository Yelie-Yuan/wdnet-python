# Import C++ standard library types using cimport
from libcpp.vector cimport vector
from libcpp.pair cimport pair

# Import the C++ function declaration from the header file
cdef extern from "../_wdnet/utils.h":
    pair[vector[float], vector[float]] node_strength(
        vector[int], vector[int], vector[float], int)

# Define the Python wrapper function
def node_strength_py(edgelist, edgeweight):
    # Convert Python lists to C++ vectors
    cdef vector[int] snode = [x[0] for x in edgelist]
    cdef vector[int] tnode = [x[1] for x in edgelist]
    cdef vector[float] edgeweight_cpp = edgeweight   
    cdef int nnode = max(max(snode), max(tnode)) + 1

    # Call the C++ function
    cdef pair[vector[float], vector[float]] ret = node_strength(
        snode, tnode, edgeweight_cpp, nnode)
    
    # Convert the C++ vectors to Python lists
    outs = list(ret.first)
    ins = list(ret.second)

    # Return the results
    return outs, ins
