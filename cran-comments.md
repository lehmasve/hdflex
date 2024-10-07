# Submission notes
### Reverse Dependency results
* checked with tools::dependsOnPkgs() and tools::check_packages_in_dir()
* no reverse dependencies

### Test environments
* Using R version 4.4.1 (2024-06-14)
* Using platform: aarch64-apple-darwin20
* R was compiled by
  Apple clang version 14.0.0 (clang-1400.0.29.202)
  GNU Fortran (GCC) 12.2.0
* Running under: macOS 15.0.1

### R CMD check results
Specified C++11: please drop specification unless essential
* 0 errors ✔ | 0 warnings ✔ | 1 note ✖
* Comment: C++11 specification needed for parallelization over RcppThread
