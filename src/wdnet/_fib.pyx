cdef extern from "../_wdnet/fib.h":
    double cfib(int n)

def fib(n):
    ''' Returns the nth Fibonacci number.'''
    return cfib(n)