% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-documentation.R
\docType{data}
\name{benchmark_ar2}
\alias{benchmark_ar2}
\title{AR(2) benchmark forecasts for quarterly U.S. inflation}
\format{
A \link{matrix} with 123 quarterly observations (rows) and 4 benchmarks (columns):
\describe{

   \item{GDPCTPI}{OOS-AR2-benchmark forecast for quarterly GDP deflator (GDPCTPI).}

   \item{PCECTPI}{OOS-AR2-benchmark forecast for quarterly PCE deflator (PCECTPI).}

   \item{CPIAUCSL}{OOS-AR2-benchmark forecast for quarterly Total CPI (CPIAUCSL).}

   \item{CPILFESL}{OOS-AR2-benchmark forecast for quarterly Core CPI (CPILFESL).}
}
}
\source{
<https://doi.org/10.1111/iere.12623>
}
\usage{
benchmark_ar2
}
\description{
Out-of-sample one-step-ahead AR(2) benchmark forecasts for the period
             from 1991-Q2 to 2021-Q4. The AR(2) models are estimated with OLS and intercept.
}
\references{
Koop, G. and Korobilis, D. (2023) "Bayesian dynamic variable selection in high dimensions."
\emph{International Economic Review}.
}
\keyword{datasets}
