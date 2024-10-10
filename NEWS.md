# hdflex 0.3.0
* Enhanced parallelization using RcppThreads for the `stsc()` function.
* Improved (computational) performance 
* Added S3 class method for `stsc` and `dsc` objects: `summary.stsc_obj` and `summary.dsc_obj` for generating plots showing the evolution of the tuning parameter, as well as standard accuracy metrics such as Mean-Squared-Error, Continuous-Ranked-Probability-Score, and Predictive-Log-Likelihood-Score.
* Introduction of the new argument `bias` for `stsc()` and `tvc`, allowing users to decide whether bias correction should be applied to the F-Signals in the TVC-models.
* Addition of the new argument `incl` for `stsc()` and `dsc`, enabling users to specify whether certain signals are required to be included in the subsets.
* Improved internal structure and performance for `dsc()`.
* Renamed the argument `burn_in_tvc` to `burn_in` and `sample_length` to `init`.
* Consolidated the arguments `risk_aversion`, `min_weight`, and `max_weight` into `portfolio_arguments`.

# hdflex 0.2.1
* Fixed a bug in the computation of the time-varying coefficients in the first step of the `stsc()` method.
* The forgetting factor `delta` in the second step of the `stsc()` method now already applies to the most recent predictive likelihood score in t-1, as stated in Equation (13) in Adaemmer et al. (2023). Previously, the score in t-1 was given a weight of 1.0
* Added new argument to `stsc()` to decide whether the subset combinations in the second step of the method should be combined with equal weights (as proposed in Adaemmer et al. (2023)) or with weights derived from the predictive log-likelihood scores.

# hdflex 0.2.0
* Added the function `stsc()` to directly apply the STSC-algorithm from Adaemmer, Lehmann and Schuessler (2023). This function is faster and more memory efficient than subsequently applying `tvc()` and `dsc()` as it is now completely written in Rcpp.
* Fixed the package overview help file.
* Updated documentation
* Updated example

# hdflex 0.1.0
* Added a `NEWS.md` file to track changes to the package.
