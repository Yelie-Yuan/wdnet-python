#include "fib.h"

void cfib(int n, int& ret) {
    int i;
    double a=0.0, b=1.0, tmp;
    for (i=0; i<n; ++i) {
        tmp = a; a = a + b; b = tmp;
    }
    ret = (int)a; 
}
