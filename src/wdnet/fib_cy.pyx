cdef extern from "../_wdnet/fib.h":
    void cfib(int n, int&ret)

def fib(n):
    ''' Returns the nth Fibonacci number.'''
    cdef int ret
    cfib(n, ret)
    return "result is " + str(ret)