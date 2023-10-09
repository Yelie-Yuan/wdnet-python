# Random Number Generation

In `R`, the `Rcpp` library is utilized to synchronize the Random
Number Generator (RNG) states between `R` and `C++`. However, I
haven't found a method to achieve this in `Python`. Therefore,
`numpy.random` is employed to generate random numbers in this `Python`
package. The `C++` (`Rcpp`) functions involving random number
generation have now been rewritten in `Cython`.
