#' @name stsc
#' @title Signal-Transform-Subset-Combination (STSC)
#' @description `stsc()` is a time series forecasting method designed to handle
#' vast sets of predictive signals, many of which are irrelevant or short-lived.
#' The method transforms heterogeneous scalar-valued signals into
#' candidate density forecasts via time-varying coefficient models (TV-C),
#' and subsequently, combines them into a final density forecast
#' via dynamic subset combination (DSC).
#' @param y A matrix of dimension `T * 1` or numeric vector of length `T`
#' containing the observations of the target variable.
#' @param X A matrix with `T` rows containing
#' the lagged 'simple' signals in each column.
#' Use NULL if no 'simple' signal shall be included.
#' @param Ext_F A matrix with `T` rows containing
#' point forecasts of y in each column.
#' Use NULL if no point forecasts shall be included.
#' @param sample_length An integer that denotes the number of observations used
#' to initialize the observational variance and the coefficients' variance
#' in the TV-C models.
#' @param lambda_grid A numeric vector with values between 0 and 1 denoting the
#' discount factor(s) that control the dynamics of the time-varying
#' coefficients. Each signal in combination with each value of
#' lambda provides a separate candidate forecast.
#' Constant coefficients are nested for the case `lambda = 1`.
#' @param kappa_grid A numeric vector between 0 and 1 to accommodate
#' time-varying volatility in the TV-C models. The observational variance
#' is estimated via Exponentially Weighted Moving Average.
#' Constant variance is nested for the case `kappa = 1`.
#' Each signal in combination with each value of
#' kappa provides a separate forecast.
#' @param burn_in_tvc An integer value `>= 1` that denotes the number of
#' observations used to 'initialize' the TV-C models.
#' After 'burn_in_tvc' observations the generated sum of discounted
#' predictive log-likelihoods (DPLLs) of each Candidate Model (TV-C model)
#' and Subset Combination (combination of gamma and psi) is resetted.
#' `burn_in_tvc = 1` means no burn-in period is applied.
#' @param gamma_grid A numerical vector that contains discount factors
#' between 0 and 1 to exponentially down-weight the past predictive performance
#' of the candidate forecasts.
#' @param psi_grid An integer vector that controls
#' the (possible) sizes of the active subsets.
#' @param delta A numeric value between 0 and 1 denoting the discount factor
#' used to down-weight the past predictive performance of the
#' subset combinations.
#' @param burn_in_dsc An integer value `>= 1` that denotes the number of
#' observations used to 'initialize' the Dynamic Subset Combinations.
#' After 'burn_in_dsc' observations the generated sum of discounted
#' predictive log-likelihoods (DPLLs) of each Subset Combination
#' (combination of gamma and psi) is resetted.
#' `burn_in_dsc = 1` means no burn-in period is applied.
#' @param method An integer of the set `1, 2, 3, 4` that denotes
#' the method used to rank the Candidate Models (TV-C models)
#' and Subset Combinations according to their performance.
#' Default is `method = 1` which ranks according to their
#' generated sum of discounted predictive log-likelihoods (DPLLs).
#' `method = 2` uses Squared-Error (SE) instead of DPLLs.
#' `method = 3` uses Absolute-Error (AE) and
#' `method = 4` uses Compounded-Returns
#' (in this case the target variable y has to be a time series of
#' financial returns).
#' @param equal_weight A boolean that denotes whether equal weights are used to
#' combine the candidate forecasts within a subset. If `FALSE`, the weights are
#' calculated using the softmax-function on the predictive log-scores of
#' the candidate models. The method proposed in Adaemmer et al (2023) uses
#' equal weights to combine the candidate forecasts.
#' @param risk_aversion A double `>= 0` that denotes the risk aversion
#' of an investor. A higher value indicates a risk avoiding behaviour.
#' @param min_weight A double that denotes the lower bound
#' for the weight placed on the market.
#' A non-negative value rules out short sales.
#' @param max_weight A double that denotes the upper bound
#' for the weight placed on the market.
#' A value of e.g. 2 allows for a maximum leverage ratio of two.
#' @return A list that contains:
#' * (1) a vector with the first moments (point forecasts) of the STSC-Model,
#' * (2) a vector with the second moments (variance) of the STSC-Model,
#' * (3) a vector that contains the selected values for gamma,
#' * (4) a vector that contains the selected values for psi and
#' * (5) a matrix that indicates the selected signals for every point in time.
#'
#' @seealso \url{https://github.com/lehmasve/hdflex#readme}
#' @author Philipp Adämmer, Sven Lehmann, Rainer Schüssler
#' @references
#' Beckmann, J., Koop, G., Korobilis, D., and Schüssler, R. A. (2020) "Exchange rate predictability and dynamic bayesian learning."
#' \emph{Journal of Applied Econometrics}, 35 (4): 410–421.
#'
#' Dangl, T. and Halling, M. (2012) "Predictive regressions with time-varying coefficients."
#' \emph{Journal of Financial Economics}, 106 (1): 157–181.
#'
#' Del Negro, M., Hasegawa, R. B., and Schorfheide, F. (2016) "Dynamic prediction pools: An investigation of financial frictions and forecasting performance."
#' \emph{Journal of Econometrics}, 192 (2): 391–405.
#'
#' Koop, G. and Korobilis, D. (2012) "Forecasting inflation using dynamic model averaging."
#' \emph{International Economic Review}, 53 (3): 867–886.
#'
#' Koop, G. and Korobilis, D. (2023) "Bayesian dynamic variable selection in high dimensions."
#' \emph{International Economic Review}.
#'
#' Raftery, A. E., Kárn`y, M., and Ettler, P. (2010) "Online prediction under model uncertainty via dynamic model averaging: Application to a cold rolling mill."
#' \emph{Technometrics}, 52 (1): 52–66.
#'
#' West, M. and Harrison, J. (1997) "Bayesian forecasting and dynamic models"
#' \emph{Springer}, 2nd edn.
#'
#' @export
#' @import checkmate
#' @examples
#' \donttest{
#'
#'    #########################################################
#'    ######### Forecasting quarterly U.S. inflation ##########
#'    #### Please see Koop & Korobilis (2023) for further  ####
#'    #### details regarding the data & external forecasts ####
#'    #########################################################
#'
#'    # Packages
#'    library("hdflex")
#'
#'    ########## Get Data ##########
#'    # Load Data
#'    inflation_data   <-  inflation_data
#'    benchmark_ar2    <-  benchmark_ar2
#'
#'    # Set Index for Target Variable
#'    i  <-  1   # (1 -> GDPCTPI; 2 -> PCECTPI; 3 -> CPIAUCSL; 4 -> CPILFESL)
#'
#'    # Subset Data (keep only data relevant for target variable i)
#'    dataset  <-  inflation_data[, c(1+(i-1),                          # Target Variable
#'                                    5+(i-1),                          # Lag 1
#'                                    9+(i-1),                          # Lag 2
#'                                    (13:16)[-i],                      # Remaining Price Series
#'                                    17:452,                           # Exogenous Predictor Variables
#'                                    seq(453+(i-1)*16,468+(i-1)*16))]  # Ext. Point Forecasts
#'
#'    ########## STSC ##########
#'    # Set Target Variable
#'    y  <-  dataset[,  1, drop = FALSE]
#'
#'    # Set 'Simple' Signals
#'    X  <-  dataset[, 2:442, drop = FALSE]
#'
#'    # Set External Point Forecasts (Koop & Korobilis 2023)
#'    Ext_F  <-  dataset[, 443:458, drop = FALSE]
#'
#'    # Set Dates
#'    dates  <-  rownames(dataset)
#'
#'    # Set TV-C-Parameter
#'    sample_length  <-  4 * 5
#'    lambda_grid    <-  c(0.90, 0.95, 1)
#'    kappa_grid     <-  0.98
#'
#'    # Set DSC-Parameter
#'    gamma_grid  <-  c(0.40, 0.50, 0.60, 0.70, 0.80, 0.90,
#'                      0.91, 0.92, 0.93, 0.94, 0.95, 0.96,
#'                      0.97, 0.98, 0.99, 1.00)
#'    psi_grid    <-  c(1:100)
#'    delta       <-  0.95
#'
#'    # Apply STSC-Function
#'    results  <-  hdflex::stsc(y,
#'                              X,
#'                              Ext_F,
#'                              sample_length,
#'                              lambda_grid,
#'                              kappa_grid,
#'                              burn_in_tvc = 79,
#'                              gamma_grid,
#'                              psi_grid,
#'                              delta,
#'                              burn_in_dsc = 1,
#'                              method = 1,
#'                              equal_weight = TRUE,
#'                              risk_aversion = NULL,
#'                              min_weight = NULL,
#'                              max_weight = NULL)
#'
#'    # Assign DSC-Results
#'    forecast_stsc    <-  results[[1]]
#'    variance_stsc    <-  results[[2]]
#'    chosen_gamma     <-  results[[3]]
#'    chosen_psi       <-  results[[4]]
#'    chosen_signals   <-  results[[5]]
#'
#'    # Define Evaluation Period
#'    eval_date_start      <-  "1991-01-01"
#'    eval_date_end        <-  "2021-12-31"
#'    eval_period_idx      <-  which(dates > eval_date_start & dates <= eval_date_end)
#'
#'    # Trim Objects
#'    oos_y                <-  y[eval_period_idx, ]
#'    oos_forecast_stsc    <-  forecast_stsc[eval_period_idx]
#'    oos_variance_stsc    <-  variance_stsc[eval_period_idx]
#'    oos_chosen_gamma     <-  chosen_gamma[eval_period_idx]
#'    oos_chosen_psi       <-  chosen_psi[eval_period_idx]
#'    oos_chosen_signals   <-  chosen_signals[eval_period_idx, , drop = FALSE]
#'    oos_dates            <-  dates[eval_period_idx]
#'
#'    # Add Dates
#'    names(oos_forecast_stsc)     <-  oos_dates
#'    names(oos_variance_stsc)     <-  oos_dates
#'    names(oos_chosen_gamma)      <-  oos_dates
#'    names(oos_chosen_psi)        <-  oos_dates
#'    rownames(oos_chosen_signals) <-  oos_dates
#'
#'    ### Part 2: Evaluation ###
#'    # Apply Summary-Function
#'    summary_results  <-  summary_stsc(oos_y,
#'                                      benchmark_ar2[, i],
#'                                      oos_forecast_stsc)
#'
#'    # Assign Summary-Results
#'    cssed  <-  summary_results[[3]]
#'    mse    <-  summary_results[[4]]
#'
#'    ########## Results ##########
#'    # Relative MSE
#'    print(paste("Relative MSE:", round(mse[[1]] / mse[[2]], 4)))
#'
#'    # Plot CSSED
#'    plot(x    = as.Date(oos_dates),
#'         y    = cssed,
#'         ylim = c(-0.0008, 0.0008),
#'         main = "Cumulated squared error differences",
#'         type = "l",
#'         lwd  = 1.5,
#'         xlab = "Date",
#'         ylab = "CSSED") + abline(h = 0, lty = 2, col = "darkgray")
#'
#'    # Plot Predictive Signals
#'    vec  <-  seq_len(dim(oos_chosen_signals)[2])
#'    mat  <-  oos_chosen_signals %*% diag(vec)
#'    mat[mat == 0]  <- NA
#'    matplot(x    = as.Date(oos_dates),
#'            y    = mat,
#'            cex  = 0.4,
#'            pch  = 20,
#'            type = "p",
#'            main = "Evolution of selected signal(s)",
#'            xlab = "Date",
#'            ylab = "Predictive Signal")
#'
#'    # Plot Psi
#'    plot(x    = as.Date(oos_dates),
#'         y    = oos_chosen_psi,
#'         ylim = c(1, 100),
#'         main = "Evolution of the subset size",
#'         type = "p",
#'         cex  = 0.75,
#'         pch  = 20,
#'         xlab = "Date",
#'         ylab = "Psi")
#'  }

### New STSC-Function
stsc  <-  function(y,
                   X,
                   Ext_F,
                   sample_length,
                   lambda_grid,
                   kappa_grid,
                   burn_in_tvc,
                   gamma_grid,
                   psi_grid,
                   delta,
                   burn_in_dsc,
                   method,
                   equal_weight,
                   risk_aversion = NULL,
                   min_weight = NULL,
                   max_weight = NULL) {

  ########################################################
  ### Checkmate
  # Check if y is numeric vector without missing / infinite values
  checkmate::assertNumeric(y,
                           min.len = 1,
                           any.missing = FALSE,
                           finite = TRUE)

  # Either x or f must not be null
  checkmate::assert(checkmate::checkMatrix(X),
                    checkmate::checkMatrix(Ext_F),
                    combine = "or")

  # Check if x is numeric matrix and has the same number of observations as y
  checkmate::assertMatrix(X,
                          mode = "numeric",
                          nrow = length(y),
                          null.ok = TRUE)

  # Check if Ext_F is numeric matrix and has the same number of observations as y
  checkmate::assertMatrix(Ext_F,
                          mode = "numeric",
                          nrow = length(y),
                          null.ok = TRUE)

  # Check if sample_length is Integer between 2 and N
  checkmate::assertInt(sample_length,
                       lower = 2,
                       upper = length(y))

  # Check if lambda_grid is numeric vector with values between 0 and 1
  checkmate::assertNumeric(lambda_grid,
                           lower = 0,
                           upper = 1,
                           min.len = 1,
                           any.missing = FALSE,
                           finite = TRUE)

  # Check if kappa_grid is numeric vector with values between 0 and 1
  checkmate::assertNumeric(kappa_grid,
                           lower = 0,
                           upper = 1,
                           min.len = 1,
                           any.missing = FALSE,
                           finite = TRUE)

  # Check if gamma_grid is numeric vector with values between 0 and 1
  checkmate::assertNumeric(gamma_grid,
                           lower = 0,
                           upper = 1,
                           min.len = 1,
                           any.missing = FALSE,
                           finite = TRUE)

  # Check if psi_grid is Integer vector with values between 1 and ncol(mu_mat)
  checkmate::assertIntegerish(psi_grid,
                              lower = 1,
                              upper = (ncol(cbind(if (exists("X")) X,
                                                  if (exists("Ext_F")) Ext_F)) *
                                      length(lambda_grid) * length(kappa_grid)),
                              min.len = 1,
                              any.missing = FALSE,
                              null.ok = FALSE)

  # Check if delta is Numeric between 0 and 1
  checkmate::assertNumber(delta,
                          lower = exp(-10),
                          upper = 1,
                          finite = TRUE,
                          na.ok = FALSE,
                          null.ok = FALSE)

  # Check if burn_in_tvc is Integer between 1 and N
  checkmate::assertInt(burn_in_tvc,
                       lower = 1,
                       upper = length(y),
                       na.ok = FALSE,
                       null.ok = FALSE)

  # Check if burn_in_dsc is Integer between 1 and N
  checkmate::assertInt(burn_in_dsc,
                       lower = 1,
                       upper = length(y),
                       na.ok = FALSE,
                       null.ok = FALSE)

  # Check if method is element of set {1, 2, 3, 4}
  checkmate::assertChoice(method,
                          c(1, 2, 3, 4),
                          null.ok = FALSE)

  # Check if equal_weight is Boolean
  checkmate::assertLogical(equal_weight,
                           len = 1,
                           any.missing = FALSE)

  # Check if method == 4: risk_aversion, min_weight & max_weight are given ...
  # ... & y are returns
  if (checkmate::testChoice(method, c(4))) {

    # Check if returns
    checkmate::assertNumeric(y, lower = -1)

    # Check if not NULL
    checkmate::assert(checkmate::checkNumber(risk_aversion, na.ok = FALSE),
                      checkmate::checkNumber(min_weight, na.ok = FALSE),
                      checkmate::checkNumber(max_weight, na.ok = FALSE),
                      combine = "and")
  }

  # Check if risk_aversion is Numeric between 0 and 1
  checkmate::assertNumber(risk_aversion,
                          lower = 0,
                          upper = Inf,
                          null.ok = TRUE,
                          na.ok = FALSE)

  # Check if min_weight Number smaller than max_weight
  checkmate::assertNumber(min_weight,
                          lower = -Inf,
                          upper = max_weight,
                          null.ok = TRUE,
                          na.ok = FALSE)

  # Check if max_weight is Number greater than min_weight
  checkmate::assertNumber(max_weight,
                          lower = min_weight,
                          upper = Inf,
                          null.ok = TRUE,
                          na.ok = FALSE)

  # Check if there are only NA values in any row
  all_na <- any(apply(cbind(if (exists("X")) X,
                            if (exists("Ext_F")) Ext_F),
                      1, function(x) sum(is.na(x)) == length(x)))
  checkmate::assertFALSE(all_na)

  # Check if there are any Na-values not from the start
  non_consec_na  <-  any(apply(is.na(cbind(if (exists("X")) X,
                                           if (exists("Ext_F")) Ext_F)),
                               2, function(x) {
                                    sum(abs(diff(x))) > 1 |
                                    sum(diff(x)) == 1 }))
  checkmate::assertFALSE(non_consec_na)
  ########################################################

  # Convert y to matrix
  y  <-  as.matrix(y)

  # Apply Rcpp-Function
  stsc_results  <-  stsc_loop(y,
                              X,
                              Ext_F,
                              sample_length,
                              lambda_grid,
                              kappa_grid,
                              burn_in_tvc,
                              gamma_grid,
                              psi_grid,
                              delta,
                              burn_in_dsc,
                              method,
                              equal_weight,
                              risk_aversion,
                              min_weight,
                              max_weight)

  # Assign Results
  stsc_forecast  <-  stsc_results[[1]]
  stsc_variance  <-  stsc_results[[2]]
  stsc_comb_mod  <-  stsc_results[[3]]
  stsc_cand_mod  <-  stsc_results[[4]]

  # Get Values for Gamma & Psi
  para_grid     <-  expand.grid(psi_grid, gamma_grid)
  chosen_psi    <-  para_grid[stsc_comb_mod + 1, 1]
  chosen_gamma  <-  para_grid[stsc_comb_mod + 1, 2]

  # Create / Get Raw-Signal Names
  if (!is.null(X)) {
    if (!is.null(colnames(X))) {
      x_names  <-  colnames(X)
    } else {
      x_names  <-  paste0("X", as.character(seq_len(ncol(X))))
    }
  }

  # Create / Get Point-Forecast Names
  if (!is.null(Ext_F)) {
    if (!is.null(colnames(Ext_F))) {
      f_names  <-  colnames(Ext_F)
    } else {
      f_names  <-  paste0("Ext_F", as.character(seq_len(ncol(Ext_F))))
    }
  }

  # Combine Signal Names
  signal_names  <-  c(if (exists("x_names")) x_names,
                      if (exists("f_names")) f_names)

  # Create Signal-Parameter-Grid
  signal_grid  <-  expand.grid(signal_names,
                               kappa_grid,
                               lambda_grid,
                               stringsAsFactors = FALSE)

  # Set up matrix for selected signals
  mat  <-  matrix(0,
                  nrow = nrow(y),
                  ncol = length(signal_names),
                  dimnames = list(NULL, signal_names))

  # Fill matrix with selected signals
  for (t in seq(max(burn_in_tvc, burn_in_dsc) + 1, nrow(y))) {
    col_names  <-  signal_grid[stsc_cand_mod[[t]] + 1, 1]
    mat[t, col_names]  <- 1
  }

  # Return Results
  return(list(stsc_forecast,
              stsc_variance,
              chosen_gamma,
              chosen_psi,
              mat))
}