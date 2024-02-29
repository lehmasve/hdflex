# Submission notes

## Purpose
* Fixed a bug in the computation of the time-varying coefficients in the first step of the `stsc()` method.
* The forgetting factor `delta` in the second step of the `stsc()` method now already applies to the most recent predictive likelihood score in t-1, as stated in Equation (13) in Adaemmer et al. (2023). Previously, the score in t-1 was given a weight of 1.0
* Added new argument to `stsc()` to decide whether the subset combinations in the second step of the method should be combined with equal weights (as proposed in Adaemmer et al. (2023)) or with weights derived from the predictive log-likelihood scores.

## Test environments
* R version 4.3.2 (2023-10-31)
* Platform: aarch64-apple-darwin20 (64-bit)
* R was compiled by
    Apple clang version 14.0.0 (clang-1400.0.29.202)
    GNU Fortran (GCC) 12.2.0
* Running under: macOS 14.3.1

## R CMD check results
0 errors ✔ | 0 warnings ✔ | 0 notes ✔