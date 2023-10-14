# Random Number Generation

In `R`, the `Rcpp` library is utilized to synchronize the Random
Number Generator (RNG) states between `R` and `C++`. However, I
haven't found a method to achieve this in `Python`. Therefore,
`numpy.random` is employed to generate random numbers in this `Python`
package. The `C++` (`Rcpp`) functions involving random number
generation have now been rewritten in `Cython`.

## Possible alternatives

### Prepare random numbers in Python

+ Generate random numbers in `Python` and pass them to `Cython` or
  `C++` as arguments.
+ Easier to implement.
+ Control the RNG state via `numpy.random.seed()`.
+ The number of random numbers to be generated is not known in
  advance. To avoid out-of-bounds access, perhaps, the second approach
  is better.
+ In `Python`, the default `random.seed()` and `numpy.random.seed()`
  are independent of each other. If the `wdnet` package uses
  `numpy.random.seed()` to control the RNG state, users still need to
  set random seeds for their own functions if they do not use
  `numpy.random` for random number generation in their own functions.

### Pass random seed as an integer argument to `C++`

+ Add a new `int` arguments (`seed`) in `Python` functions and pass
  the arguments to `C++` functions.
  - Is this appropriate?
+ Sample all random numbers in `C++` functions.
+ Not convenient for reproducibility when using user-defined functions
  for sampling edge weights and number of edges per step (in `rpanet`,
  `RPACtrl.newedge`, and `RPACtrl.edgeweight`).
  - Users will need to set random seeds for their own functions and
  for `rpanet` separately.
  - Need to add a note in the documentation.
