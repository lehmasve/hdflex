% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dsc.R
\name{dsc}
\alias{dsc}
\title{Generate dynamic subset forecast combinations}
\usage{
dsc(
  y,
  point_forecasts,
  variance_forecasts,
  gamma_grid,
  psi_grid,
  delta,
  burn_in,
  burn_in_dsc,
  metric,
  equal_weight,
  incl,
  portfolio_params = NULL
)
}
\arguments{
\item{y}{A matrix of dimension \code{T * 1} or numeric vector of length \code{T}
containing the observations of the target variable.}

\item{point_forecasts}{A matrix with \code{T} rows containing
the first moments of (conditionally) normal distributed
predictive densities in each column.}

\item{variance_forecasts}{A matrix with \code{T} rows containing
the second moments of (conditionally) normal distributed
predictive densities in each column.}

\item{gamma_grid}{A numeric vector containing potential discount factors
between 0 and 1 to exponentially down-weight the past predictive performance
of the candidate forecasting models. The values of this tuning parameter
are chosen in a procedure that amounts to leave-one-out cross-validation,
taking into account the time series structure of the data.
For details, \emph{see Adaemmer et al. (2023)}.}

\item{psi_grid}{An integer vector that controls
the (possible) sizes of the subsets. The values of this tuning parameter
are chosen in a procedure that amounts to leave-one-out cross-validation,
taking taking into account the time series structure of the data.
For details, \emph{see Adaemmer et al. (2023)}.}

\item{delta}{A numeric value between 0 and 1 denoting the discount factor
applied to down-weight the past predictive performance of the
aggregate predictive densities.}

\item{burn_in}{An integer value \verb{>= 1} that denotes the number of
observations used to 'initialize' the rankings.
After 'burn_in' observations, the rankings for both,
the candidate forecasting models and aggregate predictive densities
are reset. \code{burn_in = 1} means no burn-in period is applied.}

\item{burn_in_dsc}{An integer value \verb{>= 1} that denotes the number of
observations used to 'initialize' the rankings.
After 'burn_in_dsc' observations, only the ranking of the
aggregate predictive densities is reset.
\code{burn_in_dsc = 1} means no burn-in period is applied.}

\item{metric}{An integer from the set \verb{1, 2, 3, 4, 5} representing
the metric used to rank the candidate forecasting models (TV-C models)
and subset combinations based on their predictive performance.
The default value is \code{metric = 5} which ranks them according to the
sum of (discounted) Continuous-Ranked-Probability-Scores (CRPS).
\code{metric = 1} uses discounted Predictive Log-Likelihoods,
\code{metric = 2} uses discounted Squared-Errors,
\code{metric = 3} uses discounted Absolute-Errors,
\code{metric = 4} uses discounted Compounded-Returns
(in this case the target variable y has to be a time series of
financial returns).}

\item{equal_weight}{A boolean that denotes whether equal weights are used to
combine the candidate forecasts within a subset. If \code{FALSE}, the weights are
calculated applying the softmax function on the ranking scores of
the candidate forecasting models. The method proposed in
Adaemmer et al. (2023) uses equal weights to combine the
candidate forecasting models.}

\item{incl}{An optional integer vector that denotes signals that
must be included in the subset combinations. For example, \code{incl = c(1, 3)}
includes all candidate forecasting models generated by
the first and third signals. If \code{NULL}, no signal is forced to be included.}

\item{portfolio_params}{A numeric vector of length 3
containing the following elements:
\describe{
\item{risk_aversion}{
A non-negative double representing the investor's
risk aversion. Higher values indicate more risk-averse behavior.
}
\item{min_weight}{
A double specifying the minimum weight allocated to the market.
A non-negative lower bound effectively rules out short sales.
}
\item{max_weight}{
A double specifying the maximum weight allocated to the market.
For example, a value of 2 allows for a maximum leverage ratio of two.
}
}
This parameter is only required if \code{metric = 4}.}
}
\value{
A list containing:
\describe{
\item{Forecasts}{A list containing:
\describe{
\item{Realization}{
A vector with the actual values of the target variable.
}
\item{Point_Forecasts}{
A vector with the first moments of the aggregate predictive densities
of the DSC model.
}
\item{Variance_Prediction}{
A vector with the second moments of the aggregate predictive densities
of the DSC model.
}
}
}
\item{Tuning_Parameters}{A list containing:
\describe{
\item{Gamma}{
A vector containing the selected values for the tuning parameter gamma.
}
\item{Psi}{
A vector containing the selected values for the tuning parameter psi.
}
\item{CFM}{
A matrix containing the selected candidate forecasting models.
}
}
}
\item{Model}{A list containing:
\describe{
\item{Gamma_grid}{The grid of gamma values used in the model.}
\item{Psi_grid}{The grid of psi values used in the model.}
\item{Delta}{The delta value used in the model.}
\item{Burn_in}{The burn-in period used in the model.}
\item{Burn_in_dsc}{The burn-in period used in the model.}
\item{Metric}{The ranking metric used in the model.}
\item{Equal_weight}{A boolean indicating if equal weighting was used.}
\item{Incl}{Additional included parameters.}
}
}
}
}
\description{
The \code{dsc()} function generates
dynamic forecast combinations from a set of
candidate density forecasts. For each period,
it selects and combines a subset of predictive densities
with the highest ranks regarding local predictive accuracy.
The identities of the candidate forecasting models and
the subset sizes used for building the aggregate predictive density
may vary over time based on the data.
If only one candidate forecast is picked,
the approach temporarily collapses to pure model selection.
}
\examples{
\donttest{

# See example for tvc().

}
}
\references{
Beckmann, J., Koop, G., Korobilis, D., and Schüssler, R. A. (2020)
"Exchange rate predictability and dynamic bayesian learning."
\emph{Journal of Applied Econometrics}, 35 (4): 410–421.

Dangl, T. and Halling, M. (2012)
"Predictive regressions with time-varying coefficients."
\emph{Journal of Financial Economics}, 106 (1): 157–181.

Del Negro, M., Hasegawa, R. B., and Schorfheide, F. (2016)
"Dynamic prediction pools:
An investigation of financial frictions and forecasting performance."
\emph{Journal of Econometrics}, 192 (2): 391–405.

Koop, G. and Korobilis, D. (2012)
"Forecasting inflation using dynamic model averaging."
\emph{International Economic Review}, 53 (3): 867–886.

Koop, G. and Korobilis, D. (2023)
"Bayesian dynamic variable selection in high dimensions."
\emph{International Economic Review}.

Raftery, A. E., Kárn`y, M., and Ettler, P. (2010)
"Online prediction under model uncertainty via dynamic model averaging:
Application to a cold rolling mill."
\emph{Technometrics}, 52 (1): 52–66.

West, M. and Harrison, J. (1997)
"Bayesian forecasting and dynamic models"
\emph{Springer}, 2nd edn.
}
\seealso{
\url{https://github.com/lehmasve/hdflex#readme}
}
\author{
Philipp Adämmer, Sven Lehmann, Rainer Schüssler
}
