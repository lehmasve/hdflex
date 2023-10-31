# Submission notes

## Purpose
* Added the function `stsc()` to directly apply the STSC-algorithm from Adaemmer, Lehmann and Schuessler (2023). This function is faster and more memory efficient than subsequently applying `tvc()` and `dsc()` as it is now completely written in Rcpp.
* Fixed the package overview help file.
* Updated example. 
* Updated documentation. 

## Test environments
* R version 4.3.0 (2023-04-21)
* Platform: aarch64-apple-darwin20 (64-bit)
* Running under: macOS 14.1

## R CMD check results
0 errors ✔ | 0 warnings ✔ | 0 notes ✔