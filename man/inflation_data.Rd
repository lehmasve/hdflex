% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-documentation.R
\docType{data}
\name{inflation_data}
\alias{inflation_data}
\title{Quarterly U.S. Inflation Dataset (Total CPI)}
\format{
A \link{matrix} with 245 quarterly observations (rows) and 462 signals (columns):
\describe{
\item{Column 1}{Transformed target variable: Total CPI (CPIAUCSL)}
\item{Columns 2-3}{First and second lag of the target variable}
\item{Columns 4-442}{Lagged and transformed signals from the sources listed above}
\item{Columns 443-462}{External point forecasts available from 1976-Q1 to 2021-Q4
for quarterly Total CPI (CPIAUCSL), including:
\describe{
\item{First 12 forecasts}{Generated using regression trees,
ridge regressions, and elastic nets
over expanding and rolling windows}
\item{Remaining 8 forecasts}{Based on models discussed in Koop and Korobilis (2023)
such as Gaussian process regressions (GPR_FAC5),
Unobserved Component Stochastic Volatility (UCSV),
and Variational Bayes Dynamic Variable Selection (VBDVS_X)}
}
}
}
}
\source{
\doi{10.1111/iere.12623}
}
\usage{
inflation_data
}
\description{
A high-dimensional dataset created by \emph{Koop and Korobilis (2023)}
that integrates predictive signals from various macroeconomic and financial sources.
}
\details{
The dataset includes data from the following sources:
\itemize{
\item \strong{FRED-QD dataset} (McCracken and Ng, 2020)
\item \strong{Portfolio data} (Jurado et al., 2015)
\item \strong{Stock market predictors} (Welch and Goyal, 2008)
\item \strong{University of Michigan consumer surveys}
\item \strong{World Bank’s Pink Sheet commodity prices}
\item \strong{Key macroeconomic indicators} from the Federal Reserve Economic Data for Canada, Germany, Japan, and the United Kingdom
}

The dataset is pre-processed for one-step-ahead forecasts and includes external point forecasts.
It spans from 1960-Q3 to 2021-Q4.
}
\references{
Jurado, K., Ludvigson, S. C., and Ng, S. (2015) "Measuring uncertainty." \emph{American Economic Review}, 105 (3): 1177–1216.

Koop, G. and Korobilis, D. (2023) "Bayesian dynamic variable selection in high dimensions." \emph{International Economic Review}.

McCracken, M., and S. Ng (2020) “FRED-QD: A Quarterly Database for Macroeconomic Research” \emph{National Bureau of Economic Research}, Working Paper 26872.

Welch, I. and Goyal, A. (2008) "A comprehensive look at the empirical performance of equity premium prediction." \emph{The Review of Financial Studies}, 21 (4): 1455–1508.
}
\keyword{datasets}
