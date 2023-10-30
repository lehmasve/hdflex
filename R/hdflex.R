#' hdflex: High-Dimensional Density Forecasts
#'
#' hdflex contains the forecasting algorithm STSC 
#' developed by Adämmer, Lehmann and Schüssler (2023) <doi:10.2139/ssrn.4342487>.
#' STSC is a novel time series forecasting method designed to handle very large
#' sets of predictive signals, many of which are irrelevant or have only
#' short-lived predictive power.
#' Please cite the paper when using the package.
#'
#' @docType package
#' @author Philipp Adämmer, Sven Lehmann, Rainer Schüssler
#' @importFrom Rcpp evalCpp
#' @useDynLib hdflex, .registration = TRUE
#' @name hdflex
#' @aliases hdflex-package
NULL