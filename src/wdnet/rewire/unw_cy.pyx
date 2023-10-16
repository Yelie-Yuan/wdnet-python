# Import necessary libraries
import numpy as np
cimport numpy as cnp
from libc.math cimport sqrt
from libcpp cimport bool
from libcpp.algorithm cimport swap
from libcpp.vector cimport vector
from libc.stdint cimport uint32_t
cnp.import_array()

# Declare the C++ function
cdef extern from "../../_wdnet/rewire.h":
    void dprewire_directed_cpp(
        const int iteration, 
        const int nattempts, 
        vector[int]& tnode,
        const vector[double]& sout,
        const vector[double]& sin,
        vector[double]& tout,
        vector[double]& tin,
        const vector[int]& index_s,
        vector[int]& index_t,
        const vector[vector[double]]& eta,
        const bool history,
        vector[vector[int]]& rewire_history,
        const uint32_t random_seed,
        vector[double]& outout,
        vector[double]& outin,
        vector[double]& inout,
        vector[double]& inin)

def dprewire_directed_cpp_wrapper(
        int iteration,
        int nattempts,
        cnp.ndarray[cnp.int_t, ndim=1] tnode,
        cnp.ndarray[cnp.float64_t, ndim=1] sout,
        cnp.ndarray[cnp.float64_t, ndim=1] sin,
        cnp.ndarray[cnp.float64_t, ndim=1] tout,
        cnp.ndarray[cnp.float64_t, ndim=1] tin,
        cnp.ndarray[cnp.int_t, ndim=1] index_s,
        cnp.ndarray[cnp.int_t, ndim=1] index_t,
        cnp.ndarray[cnp.float64_t, ndim=2] eta,
        bool history,
        cnp.uint32_t random_seed):
    cdef:
        vector[int] c_tnode = tnode.tolist()
        vector[double] c_sout = sout.tolist()
        vector[double] c_sin = sin.tolist()
        vector[double] c_tout = tout.tolist()
        vector[double] c_tin = tin.tolist()
        vector[int] c_index_s = index_s.tolist()
        vector[int] c_index_t = index_t.tolist()
        vector[vector[double]] c_eta = [[d for d in inner] for inner in eta]
        vector[vector[int]] c_rewire_history = vector[vector[int]](iteration * nattempts if history else 1, vector[int](3, 0))
        vector[double] c_outout = vector[double](iteration, 0)
        vector[double] c_outin = vector[double](iteration, 0)
        vector[double] c_inout = vector[double](iteration, 0)
        vector[double] c_inin = vector[double](iteration, 0)

    # Call the C++ function
    dprewire_directed_cpp(
        iteration, 
        nattempts, 
        c_tnode, 
        c_sout, 
        c_sin, 
        c_tout, 
        c_tin, 
        c_index_s, 
        c_index_t, 
        c_eta, 
        history,
        c_rewire_history, 
        random_seed,
        c_outout, 
        c_outin, 
        c_inout, 
        c_inin)

    # Convert C++ vectors back to Python ndarrays
    py_tnode = np.array(c_tnode, dtype=np.int_)
    py_outout = np.array(c_outout, dtype=np.float64)
    py_outin = np.array(c_outin, dtype=np.float64)
    py_inout = np.array(c_inout, dtype=np.float64)
    py_inin = np.array(c_inin, dtype=np.float64)
    py_rewire_history = np.array(c_rewire_history, dtype=np.int_)

    return py_tnode, py_outout, py_outin, py_inout, py_inin, py_rewire_history

cdef extern from "../../_wdnet/rewire.h":
    void dprewire_undirected_cpp(
        const int iteration,
        const int nattempts,
        vector[int]& node1,
        vector[int]& node2,
        vector[double]& degree1,
        vector[double]& degree2,
        vector[int]& index1,
        vector[int]& index2,
        vector[vector[double]]& eta,
        const bool history,
        vector[vector[int]]& rewire_history,
        const uint32_t random_seed,
        vector[double]& rho)

def dprewire_undirected_cpp_wrapper(
        int iteration,
        int nattempts,
        cnp.ndarray[cnp.int_t, ndim=1] node1,
        cnp.ndarray[cnp.int_t, ndim=1] node2,
        cnp.ndarray[cnp.float64_t, ndim=1] degree1,
        cnp.ndarray[cnp.float64_t, ndim=1] degree2,
        cnp.ndarray[cnp.int_t, ndim=1] index1,
        cnp.ndarray[cnp.int_t, ndim=1] index2,
        cnp.ndarray[cnp.float64_t, ndim=2] eta,
        bool history,
        cnp.uint32_t random_seed):
    cdef:
        vector[int] c_node1 = node1.tolist()
        vector[int] c_node2 = node2.tolist()
        vector[double] c_degree1 = degree1.tolist()
        vector[double] c_degree2 = degree2.tolist()
        vector[int] c_index1 = index1.tolist()
        vector[int] c_index2 = index2.tolist()
        vector[vector[double]] c_eta = [[d for d in inner] for inner in eta]
        vector[vector[int]] c_rewire_history = vector[vector[int]](iteration * nattempts if history else 1, vector[int](4, 0))
        vector[double] c_rho = vector[double](iteration, 0)
    
    # Call the C++ function
    dprewire_undirected_cpp(
        iteration,
        nattempts,
        c_node1,
        c_node2,
        c_degree1,
        c_degree2,
        c_index1,
        c_index2,
        c_eta,
        history,
        c_rewire_history,
        random_seed,
        c_rho
    )

    # Convert C++ vectors back to Python ndarrays
    py_node1 = np.array(c_node1, dtype=np.int_)
    py_node2 = np.array(c_node2, dtype=np.int_)
    py_rho = np.array(c_rho, dtype=np.float64)
    py_rewire_history = np.array(c_rewire_history, dtype=np.int_)
    
    return py_node1, py_node2, py_rho, py_rewire_history
