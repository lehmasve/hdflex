#' @name tvc
#' @title Compute Univariate Time-Varying Coefficient (TV-C) Regressions
#' @description `tvc()` can be used to generate density forecasts based on
#' univariate time-varying coefficient regressions
#' (nesting constant coefficients as a special case)
#' with time-varying volatility. In addition to “raw” predictors,
#' `tvc()` accommodates point forecast as predictors.
#' @param y A matrix of dimension `T * 1` or Numeric Vector of length `T`
#' containing the observations of the target variable.
#' @param x A matrix with `T` rows containing
#' the lagged raw predictors in each column.
#' Use `NULL` if no "raw" predictors shall be included.
#' @param f A matrix with `T` rows containing
#' point forecasts of y in each column.
#' Use NULL if no point forecasts shall be included.
#' @param lambda_grid A numeric vector denoting the discount factor(s)
#' that control the dynamics of the coefficients.
#' Each predictor in combination with each value of
#' lambda defines a separate forecasting model.
#' Constant coefficients are nested for the case `lambda = 1`.
#' @param kappa_grid A numeric vector to accomodate time-varying volatility.
#' The observational variance is estimated via
#' Exponentially Weighted Moving Average.
#' Constant variance is nested for the case `kappa = 1`.
#' Each predictor in combination with each value of
#' kappa defines a separate forecasting model.
#' @param init_length Integer. Denotes the number of observations used
#' to initialize (estimate) the the TV-C models
#' (e.g. the observational variance).
#' @param x_names Character vector with all 'raw' predictor names.
#' Use NULL if no 'raw' predictors are included.
#' @param f_names Character vector with all point forecast names.
#' Use NULL if no point forecasts are included.
#' @param n_cores Integer. The number of cores to use for the estimation.
#' @return List that contains (1) a matrix with all forecasts,
#' (2) a matrix with all variances and (3) a vector with all model names.
#' @export
#' @import parallel
#' @import checkmate
#' @importFrom stats na.omit
#' @example
#' \donttest{
#
#   ### Simulate Data
#   set.seed(123)
#
#   # Set Dimensions
#   numb_obs   <-  500
#   numb_pred  <-  50
#   numb_forc  <-  10
#
#   # Create Random Target-Variable
#   equity_premium  <-  rnorm(n = numb_obs, mean = 0, sd = 1)
#
#   # Create Random Text-Variables
#   raw_preds            <-  replicate(numb_pred, sample(0:10, numb_obs, rep = TRUE), )
#   raw_names            <-  paste0("X", as.character(seq_len(ncol(raw_preds))))
#   colnames(raw_preds)  <-  raw_names
#
#   # Create Random Point Forecasts
#   f_preds            <-  replicate(10, rnorm(n = numb_obs, mean = 0, sd = 0.5), )
#   f_names            <-  paste0("F", as.character(seq_len(ncol(f_preds))))
#   colnames(f_preds)  <-  f_names
#
#   # Create Benchmark
#   benchmark  <-  dplyr::lag(roll::roll_mean(equity_premium,
#                                             width = length(equity_premium),
#                                             min_obs = 1), n = 1)
#
#   # Specify TV-C-Parameter
#   sample_length  <-  floor(numb_obs / 10)
#   lambda_grid    <-  c(0.9995, 0.9999, 1.0000)
#   kappa_grid     <-  c(0.94)
#   n_cores        <-  1
#
#   # Apply TV-C-Function
#   results  <-  hdflex::tvc(equity_premium,
#                            raw_preds,
#                            f_preds,
#                            lambda_grid,
#                            kappa_grid,
#                            sample_length,
#                            raw_names,
#                            f_names,
#                            n_cores)
#
#   # Assign Results
#   forecast_tvc      <-  results[[1]]
#   variance_tvc      <-  results[[2]]
#   model_names_tvc   <-  results[[3]]
#
#   # Cut Initialization-Period
#   nr_drp              <-  dim(forecast_tvc)[1] - dim(na.omit(forecast_tvc))[1] + 1
#   sample_period_idx   <-  (nr_drp + sample_length):numb_obs
#
#   # Trim Objects
#   sub_forecast_tvc    <-  forecast_tvc[sample_period_idx, , drop = FALSE]
#   sub_variance_tvc    <-  variance_tvc[sample_period_idx, , drop = FALSE]
#   sub_benchmark       <-  benchmark[sample_period_idx]
#   sub_equity_premium  <-  equity_premium[sample_period_idx]
#
#   ##### Dynamic Subset Combination #####
#   # Set DSC Parameter
#   nr_mods     <-  length(model_names_tvc)
#   gamma_grid  <-  c(0.9, 0.95, 0.99, 1)
#   psi_grid    <-  c(1, 2, 3, 4, 5)
#   ew_grid     <-  1
#   delta       <-  0.9992
#   n_cores     <-  1
#
#   # Apply DSC-Function
#   results  <-  hdflex::dsc(gamma_grid,
#                            psi_grid,
#                            ew_grid,
#                            sub_equity_premium,
#                            sub_forecast_tvc,
#                            sub_variance_tvc,
#                            delta,
#                            n_cores)
#
#   # Assign Results
#   sub_forecast_dsc    <-  results[[1]]
#   sub_variance_dsc    <-  results[[2]]
#   model_names_comb      <-  results[[3]]
#   sub_chosen_parameter  <-  results[[4]]
#   sub_models_idx        <-  results[[5]]
#
#   # Define Evaluation Period
#   eval_period_idx      <-  50:length(sub_equity_premium)
#
#   # Trim Objects
#   oos_equity_premium   <-  sub_equity_premium[eval_period_idx]
#   oos_benchmark        <-  sub_benchmark[eval_period_idx]
#   oos_forecast_dsc   <-  sub_forecast_dsc[eval_period_idx]
#   oos_variance_dsc   <-  sub_variance_dsc[eval_period_idx]
#   oos_models_idx       <-  lapply(seq_along(model_names_comb), function(i) {
#                                       sub_models_idx[[i]][eval_period_idx]})
#   oos_chosen_parameter <-  sub_chosen_parameter[eval_period_idx, , drop = FALSE]
#   oos_dates            <-  seq(as.Date("1989-12-01"), by = "day", length.out = length(eval_period_idx))
#
#   # Assign Names
#   names(oos_forecast_dsc)  <-  oos_dates
#   names(oos_variance_dsc)  <-  oos_dates
#
#   # Apply Statistial-Evaluation-Function
#   eval_results  <-  eval_dsc(oos_equity_premium,
#                              oos_benchmark,
#                              oos_forecast_dsc,
#                              oos_dates,
#                              oos_chosen_parameter,
#                              model_names_tvc,
#                              oos_models_idx)
#
#   # Assign Results
#   cw_t          <-  eval_results[[1]]
#   oos_r2        <-  eval_results[[2]]
#   csed          <-  eval_results[[3]]
#   pred_pockets  <-  eval_results[[4]]
#' }

### Time-Varying Coefficient Model
tvc  <- function(y,
                 x,
                 f,
                 lambda_grid,
                 kappa_grid,
                 init_length,
                 x_names,
                 f_names,
                 n_cores) {

  ### Checkmate
  # Check if y is numeric vector without missing / infinite values
  checkmate::assertNumeric(y,
                           min.len = 1,
                           any.missing = FALSE,
                           finite = TRUE)

  # Either x or f must not be null
  checkmate::assert(checkmate::checkMatrix(x),
                    checkmate::checkMatrix(f),
                    combine = "or")

  # If x is provided, x_names must be provided and vice versa
  checkmate::assertTRUE(
              checkmate::checkNull(x) == checkmate::checkNull(x_names))

  # If f is provided, f_names must be provided and vice versa
  checkmate::assertTRUE(
              checkmate::checkNull(f) == checkmate::checkNull(f_names))

  # Check if x is numeric matrix and has the same number of observations as y
  checkmate::assertMatrix(x,
                          mode = "numeric",
                          nrow = length(y),
                          null.ok = TRUE)

  # Check if f is numeric matrix and has the same number of observations as y
  checkmate::assertMatrix(f,
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

  # Check if init_length is Integer between 2 and N
  checkmate::assertInt(init_length,
                       lower = 2,
                       upper = length(y))

  # Check if x_names is character vector with an entry for every column of x
  checkmate::assertCharacter(x_names,
                             len = ncol(x),
                             unique = TRUE,
                             null.ok = TRUE)

  # Check if f_names is character vector with an entry for every column of f
  checkmate::assertCharacter(f_names,
                             len = ncol(f),
                             unique = TRUE,
                             null.ok = TRUE)

  # Check if n_cores is Integer bigger or equal to 1
  checkmate::assertInt(n_cores,
                       lower = 1)

  ### 1) TV-C-Model for 'raw' Predictors
  if (!is.null(x)) {

    # Set Variables and Indices to Subset Matrices
    nr_preds      <-  ncol(x)
    start_cols    <-  1
    end_cols      <-  nr_preds
    mu_tmp_seq    <-  seq(1, 2 * nr_preds, 2)
    var_tmp_seq   <-  seq(2, 2 * nr_preds, 2)

    # Set Column-Index-Grid
    col_grid  <-  seq(1, ncol(x))

    # Set Number Models
    nr_mods   <- length(lambda_grid) * length(kappa_grid) * nr_preds

    # Set up Result Matrices (for Predictive Density)
    max_length        <-  length(y)
    mu_mat_raw        <-  matrix(NA, ncol = nr_mods, nrow = max_length)
    variance_mat_raw  <-  matrix(NA, ncol = nr_mods, nrow = max_length)

    # Loop over Lambda- and Kappa-Grid
    for (i in seq_along(lambda_grid)) {
      for (j in seq_along(kappa_grid)) {

        # Set Lambda and Kappa
        lam  <-  lambda_grid[i]
        kap   <-  kappa_grid[j]

        # Set up Backend for Parallel Processing
        cores   <-  n_cores
        cl      <-  parallel::makeCluster(cores, type = "PSOCK")
        parallel::clusterExport(cl = cl, varlist = c("y",
                                                     "x",
                                                     "max_length",
                                                     "init_length",
                                                     "lam",
                                                     "kap"),
                                  envir = environment())
        parallel::clusterEvalQ(cl, library("hdflex"))

        # Parallelize with parLapply
        mu_var_tmp     <-  do.call("cbind", parallel::parLapply(cl, col_grid, function(j) {

          # Get Predictor Matrices and Set Lengths
          y_x          <-  cbind(y, x[, j])
          ts_y_x       <-  stats::na.omit(y_x)
          drop_length  <-  max_length - nrow(ts_y_x)
          ts_length    <-  nrow(ts_y_x) - 1

          # Select Y-Variable
          y_var  <-  ts_y_x[, 1]

          # Select X-Variable
          x_var  <-  ts_y_x[, -1]

          # Initialize TV-C Model
          init_results  <-  init_tvc(y_var, x_var, init_length)

          # Assign Init-Results
          theta    <-  init_results[[1]]
          cov_mat  <-  init_results[[2]]
          h        <-  init_results[[3]]

          # Update & Predict: Recursively Apply TV-C-Model-Function
          tvc_results  <-  tvc_model_loop(y_var,
                                          x_var,
                                          lam,
                                          kap,
                                          theta,
                                          cov_mat,
                                          h,
                                          ts_length,
                                          drop_length,
                                          max_length)

          # Return Predictive Densities (Mu & Variance) for Predictor j
          return(cbind(tvc_results[[1]], tvc_results[[2]]))
        }))

        # Stop Cluster
        parallel::stopCluster(cl)

        # Fill Result Matrices mu_mat and variance_mat
        mu_mat_raw[, start_cols:end_cols]        <-  mu_var_tmp[, mu_tmp_seq]
        variance_mat_raw[, start_cols:end_cols]  <-  mu_var_tmp[, var_tmp_seq]

        # Increase Column Indices
        start_cols  <-  start_cols + nr_preds
        end_cols    <-  start_cols + nr_preds - 1
      }
    }

    ### Get Model Names
    # Set up Vector
    model_names_raw  <-  rep(NA, nr_mods)
                  i  <-  1

    # Loop over lambda, kappa and col-grids
    for (lam in lambda_grid) {
      for (kap in kappa_grid)  {
        for (col_ind in col_grid) {

          # Set Col-Name
          col_name    <-  x_names[col_ind]

          # Append
          model_names_raw[i] <-  paste(col_name, lam, kap, sep = "-")
                          i  <-  i + 1
        }
      }
    }

    # Remove Objects
    rm(list = c("mu_var_tmp", "x", "mu_tmp_seq", "var_tmp_seq"))
  }

  ### 2.) TV-C-Model for Point Forecasts
  if (!is.null(f)) {

    # Set Variables and Indices to Subset Matrices mu_mat and variance_mat
    nr_preds      <-  ncol(f)
    start_cols    <-  1
    end_cols      <-  nr_preds
    mu_tmp_seq    <-  seq(1, 2 * nr_preds, 2)
    var_tmp_seq   <-  seq(2, 2 * nr_preds, 2)

    # Set Column-Index-Grid
    col_grid  <-  seq(1, ncol(f))

    # Number Models
    nr_mods  <- length(kappa_grid) * nr_preds

    # Set up Matrices
    max_length      <-  length(y)
    mu_mat_f        <-  matrix(NA, ncol = nr_mods, nrow = max_length)
    variance_mat_f  <-  matrix(NA, ncol = nr_mods, nrow = max_length)

    # Loop over Kappa Grid
    for (j in seq_along(kappa_grid)) {

      # Set Kappa
      kap   <-  kappa_grid[j]

      # Set Up Backend for Parallel Processing
      cores   <-  n_cores
      cl      <-  parallel::makeCluster(cores, type = "PSOCK")
      parallel::clusterExport(cl = cl, varlist = c("y",
                                                   "f",
                                                   "max_length",
                                                   "init_length",
                                                   "kap"),
                              envir = environment())
      parallel::clusterEvalQ(cl, library("hdflex"))

      # Parallelize with parLapply
      mu_var_tmp  <-  do.call("cbind", parallel::parLapply(cl, col_grid, function(j) {

        # Get Predicor Matrices (->  Point Forecasts) and Set Lengths
        y_x          <-  cbind(y, f[, j])
        ts_y_x       <-  na.omit(y_x)
        drop_length  <-  max_length - nrow(ts_y_x)
        ts_length    <-  nrow(ts_y_x) - 1

        # Select Y-Variable
        y_var  <-  ts_y_x[, 1]

        # Select X-Variable
        x_var  <-  ts_y_x[, -1]

        # Initialize TV-C-Model
        init_results  <-  init_tvc_forecast(y_var, init_length)

        # Assign Init-Results
        theta    <-  init_results[[1]]
        cov_mat  <-  init_results[[2]]
        h        <-  init_results[[3]]

        # Update & Predict: Recursively Apply TV-C-Model-Function
        tvc_results  <-  tvc_model_loop_forecasts(y_var,
                                                  x_var,
                                                  kap,
                                                  theta,
                                                  cov_mat,
                                                  h,
                                                  ts_length,
                                                  drop_length,
                                                  max_length)

        # Return Predictive Densities (Mu & Variance) for Predictor j
        return(cbind(tvc_results[[1]], tvc_results[[2]]))
      }))

      # Stop Cluster
      parallel::stopCluster(cl)

      # Fill Result Matrices mu_mat and variance_mat
      mu_mat_f[, start_cols:end_cols]        <-  mu_var_tmp[, mu_tmp_seq]
      variance_mat_f[, start_cols:end_cols]  <-  mu_var_tmp[, var_tmp_seq]

      # Increase Column Indexes
      start_cols  <-  start_cols + nr_preds
      end_cols    <-  start_cols + nr_preds - 1
    }

    ### Get Model Names
    # Set up Vector
    model_names_f  <-  rep(NA, nr_mods)
          i        <-  1

    # Loop over Kappa Grid
    for (kap in kappa_grid) {
      for (col_ind in col_grid) {

        # Set Colname
        col_name  <-  f_names[col_ind]

        # Append
        model_names_f[i]  <-  paste(col_name, "1", kap, sep = "-")
                i         <-  i + 1
      }
    }

    # Remove Objects
    rm(list = c("mu_var_tmp", "f", "mu_tmp_seq", "var_tmp_seq"))
  }

  ### 3) Combine Results
  # Combine Forecasts
  forecast_tvc     <-  cbind(if (exists("mu_mat_raw")) mu_mat_raw,
                             if (exists("mu_mat_f")) mu_mat_f)

  # Combine Variances
  variance_tvc     <-  cbind(if (exists("variance_mat_raw")) variance_mat_raw,
                             if (exists("variance_mat_f")) variance_mat_f)

  # Get / Save TVc Model Names
  model_names_tvc  <-  c(if (exists("model_names_raw")) model_names_raw,
                         if (exists("model_names_f")) model_names_f)

  # Return Results
  return(list(forecast_tvc, variance_tvc, model_names_tvc))
}