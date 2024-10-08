#' @name dsc
#' @title Generate dynamic subset forecast combinations
#' @description The `dsc()` function generates
#' dynamic forecast combinations from a set of
#' candidate density forecasts. For each period,
#' it selects and combines a subset of predictive densities
#' with the highest ranks regarding local predictive accuracy.
#' The identities of the candidate forecasting models and
#' the subset sizes used for building the aggregate predictive density
#' may vary over time based on the data.
#' If only one candidate forecast is picked,
#' the approach temporarily collapses to pure model selection.
#' @param y A matrix of dimension `T * 1` or numeric vector of length `T`
#' containing the observations of the target variable.
#' @param point_forecasts A matrix with `T` rows containing
#' the first moments of (conditionally) normal distributed
#' predictive densities in each column.
#' @param variance_forecasts A matrix with `T` rows containing
#' the second moments of (conditionally) normal distributed
#' predictive densities in each column.
#' @param gamma_grid A numeric vector containing potential discount factors
#' between 0 and 1 to exponentially down-weight the past predictive performance
#' of the candidate forecasting models. The values of this tuning parameter
#' are chosen in a procedure that amounts to leave-one-out cross-validation,
#' taking into account the time series structure of the data.
#' For details, \emph{see Adaemmer et al. (2023)}.
#' @param psi_grid An integer vector that controls
#' the (possible) sizes of the subsets. The values of this tuning parameter
#' are chosen in a procedure that amounts to leave-one-out cross-validation,
#' taking taking into account the time series structure of the data.
#' For details, \emph{see Adaemmer et al. (2023)}.
#' @param delta A numeric value between 0 and 1 denoting the discount factor
#' applied to down-weight the past predictive performance of the
#' aggregate predictive densities.
#' @param burn_in An integer value `>= 1` that denotes the number of
#' observations used to 'initialize' the rankings.
#' After 'burn_in' observations, the rankings for both,
#' the candidate forecasting models and aggregate predictive densities
#' are reset. `burn_in = 1` means no burn-in period is applied.
#' @param burn_in_dsc An integer value `>= 1` that denotes the number of
#' observations used to 'initialize' the rankings.
#' After 'burn_in_dsc' observations, only the ranking of the
#' aggregate predictive densities is reset.
#' `burn_in_dsc = 1` means no burn-in period is applied.
#' @param metric An integer from the set `1, 2, 3, 4, 5` representing
#' the metric used to rank the candidate forecasting models (TV-C models)
#' and subset combinations based on their predictive performance.
#' The default value is `metric = 5` which ranks them according to the
#' sum of (discounted) Continuous-Ranked-Probability-Scores (CRPS).
#' `metric = 1` uses discounted Predictive Log-Likelihoods,
#' `metric = 2` uses discounted Squared-Errors,
#' `metric = 3` uses discounted Absolute-Errors,
#' `metric = 4` uses discounted Compounded-Returns
#' (in this case the target variable y has to be a time series of
#' financial returns).
#' @param equal_weight A boolean that denotes whether equal weights are used to
#' combine the candidate forecasts within a subset. If `FALSE`, the weights are
#' calculated applying the softmax function on the ranking scores of
#' the candidate forecasting models. The method proposed in
#' Adaemmer et al. (2023) uses equal weights to combine the
#' candidate forecasting models.
#' @param incl An optional integer vector that denotes signals that
#' must be included in the subset combinations. For example, `incl = c(1, 3)`
#' includes all candidate forecasting models generated by
#' the first and third signals. If `NULL`, no signal is forced to be included.
#' @param portfolio_params A numeric vector of length 3
#' containing the following elements:
#'  \describe{
#'   \item{risk_aversion}{
#'   A non-negative double representing the investor's
#'   risk aversion. Higher values indicate more risk-averse behavior.
#'   }
#'   \item{min_weight}{
#'   A double specifying the minimum weight allocated to the market.
#'   A non-negative lower bound effectively rules out short sales.
#'   }
#'   \item{max_weight}{
#'   A double specifying the maximum weight allocated to the market.
#'   For example, a value of 2 allows for a maximum leverage ratio of two.
#'   }
#' }
#' This parameter is only required if `metric = 4`.
#' @return A list containing:
#' \describe{
#'   \item{Forecasts}{A list containing:
#'     \describe{
#'       \item{Realization}{
#'       A vector with the actual values of the target variable.
#'       }
#'       \item{Point_Forecasts}{
#'       A vector with the first moments of the aggregate predictive densities
#'       of the DSC model.
#'       }
#'       \item{Variance_Prediction}{
#'       A vector with the second moments of the aggregate predictive densities
#'       of the DSC model.
#'       }
#'     }
#'   }
#'   \item{Tuning_Parameters}{A list containing:
#'     \describe{
#'       \item{Gamma}{
#'       A vector containing the selected values for the tuning parameter gamma.
#'       }
#'       \item{Psi}{
#'       A vector containing the selected values for the tuning parameter psi.
#'       }
#'       \item{CFM}{
#'       A matrix containing the selected candidate forecasting models.
#'       }
#'     }
#'   }
#'   \item{Model}{A list containing:
#'     \describe{
#'       \item{Gamma_grid}{The grid of gamma values used in the model.}
#'       \item{Psi_grid}{The grid of psi values used in the model.}
#'       \item{Delta}{The delta value used in the model.}
#'       \item{Burn_in}{The burn-in period used in the model.}
#'       \item{Burn_in_dsc}{The burn-in period used in the model.}
#'       \item{Metric}{The ranking metric used in the model.}
#'       \item{Equal_weight}{A boolean indicating if equal weighting was used.}
#'       \item{Incl}{Additional included parameters.}
#'     }
#'   }
#' }
#' @export
#' @seealso \url{https://github.com/lehmasve/hdflex#readme}
#' @author Philipp Adämmer, Sven Lehmann, Rainer Schüssler
#' @references
#' Beckmann, J., Koop, G., Korobilis, D., and Schüssler, R. A. (2020)
#' "Exchange rate predictability and dynamic bayesian learning."
#' \emph{Journal of Applied Econometrics}, 35 (4): 410–421.
#'
#' Dangl, T. and Halling, M. (2012)
#' "Predictive regressions with time-varying coefficients."
#' \emph{Journal of Financial Economics}, 106 (1): 157–181.
#'
#' Del Negro, M., Hasegawa, R. B., and Schorfheide, F. (2016)
#' "Dynamic prediction pools:
#' An investigation of financial frictions and forecasting performance."
#' \emph{Journal of Econometrics}, 192 (2): 391–405.
#'
#' Koop, G. and Korobilis, D. (2012)
#' "Forecasting inflation using dynamic model averaging."
#' \emph{International Economic Review}, 53 (3): 867–886.
#'
#' Koop, G. and Korobilis, D. (2023)
#' "Bayesian dynamic variable selection in high dimensions."
#' \emph{International Economic Review}.
#'
#' Raftery, A. E., Kárn`y, M., and Ettler, P. (2010)
#' "Online prediction under model uncertainty via dynamic model averaging:
#' Application to a cold rolling mill."
#' \emph{Technometrics}, 52 (1): 52–66.
#'
#' West, M. and Harrison, J. (1997)
#' "Bayesian forecasting and dynamic models"
#' \emph{Springer}, 2nd edn.
#' @import checkmate
#' @examples
#' \donttest{
#'
#' # See example for tvc().
#'
#' }

### Dynamic Subset Combination
dsc <- function(y,
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
                portfolio_params = NULL) {


  # Convert y to matrix
  y <- as.matrix(y, ncol = 1)

  ########################################################
  ### Checkmate
  # Check if numeric vector without missing or infinite values
  assertNumeric(y, min.len = 1, any.missing = FALSE, finite = TRUE)

  # Check if numeric matrix and same length as y
  assertMatrix(point_forecasts, mode = "double", any.missing = FALSE, nrows = length(y))
  assertMatrix(variance_forecasts, mode = "double", any.missing = FALSE, nrows = length(y), ncols = ncol(point_forecasts))
  qassert(variance_forecasts, c("M+[0,]"))

  # Check if numeric vectors with values between exp(-10) and 1
  assertNumeric(gamma_grid, lower = exp(-10), upper = 1, min.len = 1, any.missing = FALSE, finite = TRUE)

  # Check if integer vector with values between 1 and ncol(point_forecasts)
  assertIntegerish(psi_grid, lower = 1, upper = ncol(point_forecasts), min.len = 1, any.missing = FALSE)

  # Check if numeric value between exp(-10) and 1
  assertNumber(delta, lower = exp(-10), upper = 1, finite = TRUE, na.ok = FALSE, null.ok = FALSE)

  # Check if integers between 1 and the length of y
  assertInt(burn_in, lower = 1, upper = nrow(y), na.ok = FALSE, null.ok = FALSE)
  assertInt(burn_in_dsc, lower = 1, upper = nrow(y), na.ok = FALSE, null.ok = FALSE)

  # Check if element of the set {1, 2, 3, 4, 5}
  assertChoice(metric, c(1, 2, 3, 4, 5), null.ok = FALSE)

  # Check if boolean
  assertLogical(equal_weight, len = 1, any.missing = FALSE)

  # Additional checks if metric == 4
  if (metric == 4) {

    # Check if "returns"
    assertNumeric(y, lower = -1)

    # Check if numeric vector of length 3
    assertNumeric(portfolio_params, len = 3, any.missing = FALSE)

    # Extract values from portfolio_params
    risk_aversion <- portfolio_params[1]
    min_weight <- portfolio_params[2]
    max_weight <- portfolio_params[3]

    # Check if numeric value and at least 0.0
    assertNumber(risk_aversion, lower = 0, upper = Inf)

    # Check if min_weight is a number smaller than max_weight
    assertNumber(min_weight, lower = -Inf, upper = max_weight)

    # Check if max_weight is a number greater than min_weight
    assertNumber(max_weight, lower = min_weight, upper = Inf)
  }

  # Check if nullable integer vector
  assertIntegerish(incl, lower = 1, upper = ncol(point_forecasts), null.ok = TRUE)

  # Check if minimum psi matches the keep argument
  if (!is.null(incl)) {
    assertTRUE(min(psi_grid) >= length(incl))
  }

  ########################################################
  # Apply Rcpp-Function
  dsc_results <- dsc_(y,
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
                      portfolio_params)

  # Assign Results
  dsc_forecast <- as.numeric(dsc_results[[1]])
  dsc_variance <- as.numeric(dsc_results[[2]])
  dsc_comb_mod <- as.integer(dsc_results[[3]])
  dsc_cand_mod <- dsc_results[[4]]

  # Get Values for Gamma & Psi
  para_grid <- expand.grid(psi_grid, gamma_grid)
  chosen_psi <- para_grid[dsc_comb_mod + 1, 1]
  chosen_gamma <- para_grid[dsc_comb_mod + 1, 2]

  # Candidate Forecast Model Names
  if (!is.null(colnames(point_forecasts))) {
    cfm_names <- colnames(point_forecasts)
  } else {
    cfm_names <- paste0("CFM", seq_len(ncol(point_forecasts)))
  }

  # Set up matrix for selected CFM
  chosen_cfm <- matrix(0, nrow = nrow(y),
                       ncol = length(cfm_names),
                       dimnames = list(NULL, cfm_names))

  # Fill matrices
  for (t in seq(max(burn_in, burn_in_dsc) + 1, nrow(y))) {
    col_names <- cfm_names[dsc_cand_mod[[t]] + 1]
    chosen_cfm[t, col_names] <- 1
  }

  # Return Results
  return(
    structure(
      list(
        Forecasts = list(
          Realization = y,
          Point_Forecasts = dsc_forecast,
          Variance_Forecasts = dsc_variance
        ),
        Tuning_Parameters = list(
          Gamma = chosen_gamma,
          Psi = chosen_psi,
          CFM = chosen_cfm
        ),
        Model = list(
          Gamma_grid = gamma_grid,
          Psi_grid = psi_grid,
          Delta = delta,
          Burn_in = burn_in,
          Burn_in_dsc = burn_in_dsc,
          Metric = metric,
          Equal_weight = equal_weight,
          Incl = incl
        )
      ), class = "dsc_obj"
    )
  )
}
