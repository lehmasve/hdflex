#' @name tvc
#' @title Compute density forecasts based on univariate time-varying
#' coefficient (TV-C) models in state-space form
#' @description `tvc()` can be used to generate density forecasts based on
#' univariate time-varying coefficient models. In each forecasting model,
#' we include an intercept and one predictive signal. The predictive signal
#' either represents the value of a 'simple' signal
#' or the the value of an external point forecast.
#' All models are estimated independently from each other and
#' estimation and forecasting are carried out recursively.
#' @param y A matrix of dimension `T * 1` or numeric vector of length `T`
#' containing the observations of the target variable.
#' @param X A matrix with `T` rows containing
#' the lagged 'simple' signals in each column.
#' Use NULL if no 'simple' signal shall be included.
#' @param Ext_F A matrix with `T` rows containing
#' point forecasts of y in each column.
#' Use NULL if no point forecasts shall be included.
#' @param sample_length An integer that denotes the number of observations used
#' to initialize the observational variance and the coefficients' variance.
#' @param lambda_grid A numeric vector denoting the discount factor(s)
#' that control the dynamics of the coefficients.
#' Each signal in combination with each value of
#' lambda provides a separate candidate forecast.
#' Constant coefficients are nested for the case `lambda = 1`.
#' @param kappa_grid A numeric vector to accommodate time-varying volatility.
#' The observational variance is estimated via
#' Exponentially Weighted Moving Average.
#' Constant variance is nested for the case `kappa = 1`.
#' Each signal in combination with each value of
#' kappa provides a separate forecast.
#' @param bias A boolean to indicate whether the TVC-Models should
#' perform a bias correction to external point forecasts
#' (TRUE -> time-varying intercept) or
#' take it 'as is' (FALSE -> constant intercept of 0.0).
#' @return A list that contains:
#'
#' * (1) a matrix with the first moments (point forecasts)
#' of the conditionally normal predictive distributions and
#'
#' * (2) a matrix with the second moments (variance)
#' of the conditionally normal predictive distributions.
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
#' @export
#' @import parallel
#' @import checkmate
#' @importFrom stats na.omit
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
#'    inflation_data <- inflation_data
#'    benchmark_ar2  <- benchmark_ar2
#'
#'    # Set Index for Target Variable
#'    i <- 1   # (1 -> GDPCTPI; 2 -> PCECTPI; 3 -> CPIAUCSL; 4 -> CPILFESL)
#'
#'    # Subset Data (keep only data relevant for target variable i)
#'    dataset <- inflation_data[, c(1+(i-1),                          # Target Variable
#'                                  5+(i-1),                          # Lag 1
#'                                  9+(i-1),                          # Lag 2
#'                                  (13:16)[-i],                      # Remaining Price Series
#'                                  17:452,                           # Exogenous Predictor Variables
#'                                  seq(453+(i-1)*16,468+(i-1)*16))]  # Ext. Point Forecasts
#'
#'    ########## STSC ##########
#'    ### Part 1: TV-C Model ###
#'    # Set Target Variable
#'    y <- dataset[,  1, drop = FALSE]
#'
#'    # Set 'Simple' Signals
#'    X <- dataset[, 2:442, drop = FALSE]
#'
#'    # Set External Point Forecasts (Koop & Korobilis 2023)
#'    Ext_F <- dataset[, 443:458, drop = FALSE]
#'
#'    # Set TV-C-Parameter
#'    sample_length <- 4 * 5
#'    lambda_grid   <- c(0.90, 0.95, 1)
#'    kappa_grid    <- 0.98
#'    bias          <- TRUE
#'
#'    # Apply TV-C-Function
#'    results  <-  hdflex::tvc(y,
#'                             X,
#'                             Ext_F,
#'                             sample_length,
#'                             lambda_grid,
#'                             kappa_grid,
#'                             bias)
#'
#'    # Assign TV-C-Results
#'    forecast_tvc      <-  results[[1]]
#'    variance_tvc      <-  results[[2]]
#'
#'    # Define Burn-In Period
#'    sample_period_idx  <-  80:nrow(dataset)
#'    sub_forecast_tvc   <-  forecast_tvc[sample_period_idx, , drop = FALSE]
#'    sub_variance_tvc   <-  variance_tvc[sample_period_idx, , drop = FALSE]
#'    sub_y              <-  y[sample_period_idx, , drop = FALSE]
#'    sub_dates          <-  rownames(dataset)[sample_period_idx]
#'
#'    ### Part 2: Dynamic Subset Combination ###
#'    # Set DSC-Parameter
#'    nr_mods     <-  ncol(sub_forecast_tvc)
#'    gamma_grid  <-  c(0.40, 0.50, 0.60, 0.70, 0.80, 0.90,
#'                      0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00)
#'    psi_grid    <-  c(1:100)
#'    delta       <-  0.95
#'    n_cores     <-  1
#'
#'    # Apply DSC-Function
#'    results  <-  hdflex::dsc(gamma_grid,
#'                             psi_grid,
#'                             sub_y,
#'                             sub_forecast_tvc,
#'                             sub_variance_tvc,
#'                             delta,
#'                             n_cores)
#'
#'    # Assign DSC-Results
#'    sub_forecast_stsc    <-  results[[1]]
#'    sub_variance_stsc    <-  results[[2]]
#'    sub_chosen_gamma     <-  results[[3]]
#'    sub_chosen_psi       <-  results[[4]]
#'    sub_chosen_signals   <-  results[[5]]
#'
#'    # Define Evaluation Period
#'    eval_date_start      <-  "1991-01-01"
#'    eval_date_end        <-  "2021-12-31"
#'    eval_period_idx      <-  which(sub_dates > eval_date_start & sub_dates <= eval_date_end)
#'
#'    # Trim Objects
#'    oos_y                <-  sub_y[eval_period_idx, ]
#'    oos_forecast_stsc    <-  sub_forecast_stsc[eval_period_idx]
#'    oos_variance_stsc    <-  sub_variance_stsc[eval_period_idx]
#'    oos_chosen_gamma     <-  sub_chosen_gamma[eval_period_idx]
#'    oos_chosen_psi       <-  sub_chosen_psi[eval_period_idx]
#'    oos_chosen_signals   <-  sub_chosen_signals[eval_period_idx, , drop = FALSE]
#'    oos_dates            <-  sub_dates[eval_period_idx]
#'
#'    # Add Dates
#'    names(oos_forecast_stsc)     <-  oos_dates
#'    names(oos_variance_stsc)     <-  oos_dates
#'    names(oos_chosen_gamma)      <-  oos_dates
#'    names(oos_chosen_psi)        <-  oos_dates
#'    rownames(oos_chosen_signals) <-  oos_dates
#'
#'    ### Part 3: Evaluation ###
#'    # Apply Summary-Function
#'    summary_results  <-  summary_stsc(oos_y,
#'                                      benchmark_ar2[, i],
#'                                      oos_forecast_stsc)
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

### Time-Varying Coefficient Model
tvc  <- function(y,
                 X,
                 Ext_F,
                 sample_length,
                 lambda_grid,
                 kappa_grid,
                 bias) {

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

  # Check if F is numeric matrix and has the same number of observations as y
  checkmate::assertMatrix(Ext_F,
                          mode = "numeric",
                          nrow = length(y),
                          null.ok = TRUE)

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

  # Check if bias is Boolean
  checkmate::assertLogical(bias,
                           len = 1,
                           any.missing = FALSE)

  # Check if sample_length is Integer between 2 and N
  checkmate::assertInt(sample_length,
                       lower = 2,
                       upper = length(y))

  # Check if any column in X is constant for the first observations
  if (!is.null(X)) {
    if (any(apply(X, 2, function(x) length(unique(na.omit(x)[1:sample_length])) == 1))) {
      print("One or more columns in X are constant for the first 1:sample_length observations.")
    }
  }

  # Check if any column in Ext_F is constant for the first observations
  if (!is.null(Ext_F)) {
    if (any(apply(Ext_F, 2, function(x) length(unique(na.omit(x)[1:sample_length])) == 1))) {
      print("One or more columns in Ext_F are constant for the first 1:sample_length observations.")
    }
  }
  ########################################################

  ### Apply Rcpp-Function
  tvc_results <- tvc_(y,
                      X,
                      Ext_F,
                      sample_length,
                      lambda_grid,
                      kappa_grid,
                      bias)

  ### Assign Results
  forecast_tvc <- tvc_results[[1]]
  variance_tvc <- tvc_results[[2]]

  ### Remove
  rm(list = c("tvc_results"))

  ### Create / Get Raw-Signal Names
  x_names <- if (!is.null(X)) {
    if (!is.null(colnames(X))) {
      colnames(X)
    } else {
      paste0("X", as.character(seq_len(ncol(X))))
    }
  }

  if (!is.null(X)) {
    
    # Preallocate TVC-Model Names
    tvc_x_name <- vector("character", length(lambda_grid) * length(kappa_grid) * ncol(X))

    # Create TVC-Model Names
    i <- 1
    for (l in lambda_grid) {
      for (k in kappa_grid) {
        for (j in seq_len(ncol(X))) {
          # Create TVC-Model Name
          tvc_x_name[i] <- paste(x_names[j], l, k, sep = "_")
          i <- i + 1
        }
      }
    }
  }

  ### Create / Get Point-Forecast Names
  f_names <- if (!is.null(Ext_F)) {
    if (!is.null(colnames(Ext_F))) {
      colnames(Ext_F)
    } else {
      paste0("Ext_F", as.character(seq_len(ncol(Ext_F))))
    }
  }

  if (!is.null(Ext_F)) {

    # Preallocate TVC-Model Names
    tvc_f_name <- vector("character", length(lambda_grid) * length(kappa_grid) * ncol(Ext_F))

    # Create TVC-Model Names
    i <- 1
    for (l in lambda_grid) {
      for (k in kappa_grid) {
        for (j in seq_len(ncol(Ext_F))) {
          # Create TVC-Model Name
          tvc_f_name[i] <- paste(f_names[j], l, k, sep = "_")
          i <- i + 1
        }
      }
    }
  }

  # Combine Signal Names
  model_names_tvc  <-  c(if (exists("tvc_x_name")) tvc_x_name,
                         if (exists("tvc_f_name")) tvc_f_name)

  # Assign Model Names (-> Column Names)
  colnames(forecast_tvc)  <-  model_names_tvc
  colnames(variance_tvc)  <-  model_names_tvc

  # Return Results
  return(list(forecast_tvc, variance_tvc))
}