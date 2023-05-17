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
#' @param F A matrix with `T` rows containing
#' point forecasts of y in each column.
#' Use NULL if no point forecasts shall be included.
#' @param lambda_grid A numeric vector denoting the discount factor(s)
#' that control the dynamics of the coefficients.
#' Each signal in combination with each value of
#' lambda provides a separate candidate forecast.
#' Constant coefficients are nested for the case `lambda = 1`.
#' @param kappa_grid A numeric vector to accomodate time-varying volatility.
#' The observational variance is estimated via
#' Exponentially Weighted Moving Average.
#' Constant variance is nested for the case `kappa = 1`.
#' Each signal in combination with each value of
#' kappa provides a separate forecast.
#' @param init_length An integer that denotes the number of observations used
#' to initialize the observational variance and the coefficients' variance.
#' @param n_cores An integer that denotes the number of CPU-cores used
#' for the computation.
#' @return A list that contains
#' (1) a matrix with the first moments (point forecasts)
#' of the conditionally normal predictive distributions and
#' (2) a matrix with the second moments (variance)
#' of the conditionally normal predictive distributions.
#' @export
#' @import parallel
#' @import checkmate
#' @import tidyverse
#' @importFrom stats na.omit
#' @examples
#' \donttest{
#'
#'  # Packages
#'  library("tidyverse")
#'  library("hdflex")
#'
#'  # Set Target-Variables
#'  target_var_names  <- c("GDPCTPI", "PCECTPI", "CPIAUCSL", "CPILFESL")
#'
#'  # Load Data
#'  data  <-  inflation_data
#'
#'  # Loop over Target Variables
#'  results <-  do.call("rbind", lapply(X = seq_along(target_var_names), FUN = function(p) {
#'
#'      # Y-Column-Name
#'      y_target    <-  paste(target_var_names[p], "h_1", sep = "_")
#'      y_signal    <-  target_var_names[p]
#'      not_target  <-  setdiff(paste(target_var_names, "h_1", sep = "_"),
#'                              y_target)
#'
#'      # Create Forecast-Dataset
#'      dataset  <-  data                                                    %>%
#'                    select(-any_of(c(y_signal, not_target)))               %>%
#'                    mutate(across(all_of(y_target), .fns = list("lag_1" =
#'                           ~dplyr::lag(., 1))), .after = 2)                %>%
#'                    mutate(across(all_of(y_target), .fns = list("lag_2" =
#'                           ~dplyr::lag(., 2))), .after = 3)                %>%
#'                    mutate(across(-c(1, 2, 3, 4), dplyr::lag))             %>%
#'                    slice(-c(1:3))                                         %>%
#'                    column_to_rownames("Date")                             %>%
#'                    as.matrix()
#'
#'      # Get Dates & Length
#'      full_dates   <-  rownames(dataset)
#'      full_length  <-  length(full_dates)
#'
#'      # Create Time-Sequence for loop
#'      T_full      <-  full_length -1
#'      T_sequence  <-  122:T_full
#'
#'      ### Benchmark Model ###
#'      # Create Result Matrices for Predictions
#'      preds_ar2  <-  matrix(NA, ncol = 1, nrow = full_length,
#'                            dimnames = list(full_dates, "AR"))
#'
#'      # Create Result Matrices for Squared Errors
#'      se_ar2     <-  matrix(NA, ncol = 1, nrow = full_length,
#'                            dimnames = list(full_dates, "AR"))
#'
#'      # Loop over t
#'      for(t in T_sequence) {
#'
#'          ### Pre-Process Data ###
#'          # Train Data
#'          x_train     <-  scale(dataset[1:t, -1, drop = FALSE])
#'          y_train     <-        dataset[1:t,  1, drop = FALSE]
#'
#'          # Predict Data
#'          x_pred      <-  scale(dataset[1:(t+1), -1, drop = FALSE])[(t+1), , drop = FALSE]
#'          y_pred      <-        dataset[t+1, 1]
#'
#'          ### Model 1: AR(2) ###
#'          # Train Data
#'          x_train_ar  <-  cbind(int = 1, x_train[, c(1:2), drop = FALSE])
#'
#'          # Predict Data
#'          x_pred_ar   <-  cbind(int = 1,  x_pred[, c(1:2), drop = FALSE])
#'
#'          # Fit Regressions
#'          model_ar    <-  .lm.fit(x_train_ar, y_train)
#'
#'          # Predict & Combine
#'          pred                  <-  model_ar$coefficients %*% x_pred_ar[,]
#'          preds_ar2[t+1, "AR"]  <-  pred
#'          se_ar2[t+1, "AR"]     <-  (y_pred - pred) ** 2
#'      }
#'
#'      ##### TV-C Models #####
#'      # Set Target Variable
#'      Y  <-  dataset[,  1, drop = FALSE]
#'
#'      # Set 'Simple' Signals
#'      X  <-  dataset[, -1, drop = FALSE]
#'
#'      # Load External Point Forecasts (Koop & Korobilis 2023)
#'      F  <-  get(paste0("Ext_PF_", target_var_names[p]))
#'
#'      # Set TV-C-Parameter
#'      sample_length  <-  4 * 3
#'      lambda_grid    <-  c(0.90, 0.95, 0.99, 0.999, 1)
#'      kappa_grid     <-  0.98
#'      n_cores        <-  1
#'
#'      # Apply TV-C-Function
#'      results  <-  hdflex::tvc(Y,
#'                               X,
#'                               F,
#'                               lambda_grid,
#'                               kappa_grid,
#'                               sample_length,
#'                               n_cores)
#'
#'      # Assign Results
#'      forecast_tvc      <-  results[[1]]
#'      variance_tvc      <-  results[[2]]
#'      model_names_tvc   <-  colnames(forecast_tvc)
#'
#'      # Define Cut Length and Trim Objects
#'      sample_period_idx  <-  80:full_length
#'      sub_forecast_tvc   <-  forecast_tvc[sample_period_idx, , drop = FALSE]
#'      sub_variance_tvc   <-  variance_tvc[sample_period_idx, , drop = FALSE]
#'      sub_Y              <-  Y[sample_period_idx, , drop = FALSE]
#'      sub_dates          <-  full_dates[sample_period_idx]
#'      sub_length         <-  length(sub_dates)
#'
#'      ##### Dynamic Subset Combination #####
#'      # Set DSC-Parameter
#'      nr_mods     <-  ncol(sub_forecast_tvc)
#'      gamma_grid  <-  c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
#'                        0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99,
#'                        0.999, 1.00)
#'      psi_grid    <-  c(1:100)
#'      delta       <-  0.95
#'      n_cores     <-  1
#'
#'      # Apply DSC-Function
#'      results  <-  hdflex::dsc(gamma_grid,
#'                               psi_grid,
#'                               sub_Y,
#'                               sub_forecast_tvc,
#'                               sub_variance_tvc,
#'                               delta,
#'                               n_cores)
#'
#'      # Assign Results
#'      sub_forecast_stsc    <-  results[[1]]
#'      sub_variance_stsc    <-  results[[2]]
#'      sub_chosen_gamma     <-  results[[3]]
#'      sub_chosen_psi       <-  results[[4]]
#'      sub_pred_pockets     <-  results[[5]]
#'
#'      # Define Evaluation Period
#'      eval_date_start      <-  "1991-01-01"
#'      eval_date_end        <-  "2021-12-31"
#'      eval_period_idx      <-  which(sub_dates > eval_date_start &
#'                                     sub_dates <= eval_date_end)
#'
#'      # Trim Objects
#'      oos_Y                <-  sub_Y[eval_period_idx, ]
#'      oos_benchmark        <-  preds_ar2[rownames(preds_ar2) > eval_date_start, "AR"]
#'      oos_forecast_stsc    <-  sub_forecast_stsc[eval_period_idx]
#'      oos_variance_stsc    <-  sub_variance_stsc[eval_period_idx]
#'      oos_chosen_gamma     <-  sub_chosen_gamma[eval_period_idx]
#'      oos_chosen_psi       <-  sub_chosen_psi[eval_period_idx]
#'      oos_pred_pockets     <-  sub_pred_pockets[eval_period_idx, , drop = FALSE]
#'      oos_length           <-  length(eval_period_idx)
#'      oos_dates            <-  sub_dates[eval_period_idx]
#'
#'      # Add Dates
#'      names(oos_forecast_stsc)   <-  oos_dates
#'      names(oos_variance_stsc)   <-  oos_dates
#'      names(oos_chosen_gamma)    <-  oos_dates
#'      names(oos_chosen_psi)      <-  oos_dates
#'      rownames(oos_pred_pockets) <-  oos_dates
#'
#'      ##### Evaluate #####
#'      # Apply Statistial-Evaluation-Function
#'      summary_results  <-  summary_stsc(oos_Y,
#'                                        oos_benchmark,
#'                                        oos_forecast_stsc)
#'       # Assign MSE-Results
#'      mse  <-  summary_results[[4]]
#'
#'      # Return
#'      return(c(mse[[2]], mse[[1]]))
#'      }))
#'
#'  # MSE-Results
#'  dimnames(results)  <-  list(target_var_names, c("AR", "STSC"))
#'  round(results / results[, 1], 4)
#'  }

### Time-Varying Coefficient Model
tvc  <- function(y,
                 X,
                 F,
                 lambda_grid,
                 kappa_grid,
                 init_length,
                 n_cores) {

  ### Checkmate
  # Check if y is numeric vector without missing / infinite values
  checkmate::assertNumeric(y,
                           min.len = 1,
                           any.missing = FALSE,
                           finite = TRUE)

  # Either x or f must not be null
  checkmate::assert(checkmate::checkMatrix(X),
                    checkmate::checkMatrix(F),
                    combine = "or")

  # Check if x is numeric matrix and has the same number of observations as y
  checkmate::assertMatrix(X,
                          mode = "numeric",
                          nrow = length(y),
                          null.ok = TRUE)

  # Check if F is numeric matrix and has the same number of observations as y
  checkmate::assertMatrix(F,
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

  # Check if n_cores is integer bigger or equal to 1
  checkmate::assertInt(n_cores,
                       lower = 1)

  ### 1) TV-C-Model for 'simple' Signals
  if (!is.null(X)) {

    # Set Variables and Indices to Subset Matrices
    nr_preds      <-  ncol(X)
    start_cols    <-  1
    end_cols      <-  nr_preds
    mu_tmp_seq    <-  seq(1, 2 * nr_preds, 2)
    var_tmp_seq   <-  seq(2, 2 * nr_preds, 2)

    # Set Column-Index-Grid
    col_grid  <-  seq(1, nr_preds)

    # Set number of predictive densities to create
    nr_mods   <- length(lambda_grid) * length(kappa_grid) * nr_preds

    # Set up result matrices (for Predictive Density)
    max_length        <-  length(y)
    mu_mat_raw        <-  matrix(NA, ncol = nr_mods, nrow = max_length)
    variance_mat_raw  <-  matrix(NA, ncol = nr_mods, nrow = max_length)

    # Loop over Lambda- and Kappa-Grid
    for (i in seq_along(lambda_grid)) {
      for (j in seq_along(kappa_grid)) {

        # Set Lambda and Kappa
        lambda  <-  lambda_grid[i]
        kappa   <-  kappa_grid[j]

        # Set up Backend for Parallel Processing
        cores   <-  n_cores
        cl      <-  parallel::makeCluster(cores, type = "PSOCK")
        parallel::clusterExport(cl = cl, varlist = c("y",
                                                     "X",
                                                     "max_length",
                                                     "init_length",
                                                     "lambda",
                                                     "kappa"),
                                  envir = environment())
        parallel::clusterEvalQ(cl, library("hdflex"))

        # Parallelize with parLapply
        mu_var_tmp     <-  do.call("cbind", parallel::parLapply(cl, col_grid, function(j) {

          # Select signal-column and align lenght with target variable
          y_x          <-  cbind(y, X[, j])
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

          # Update & Predict: Recursively apply TV-C-Model-Function
          tvc_results  <-  tvc_model_loop(y_var,
                                          x_var,
                                          lambda,
                                          kappa,
                                          theta,
                                          cov_mat,
                                          h,
                                          ts_length,
                                          drop_length,
                                          max_length)

          # Return Predictive Densities (Mu & Variance) for Signal j
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

    ### Create Candidate Forecast Names
    # Get / Create Signal Names
    if (!is.null(colnames(X))) {
      x_names  <-  colnames(X)
    } else {
      x_names  <-  paste0("X", as.character(seq_len(nr_preds)))
    }

    # Set up Vector
    model_names_raw  <-  rep(NA, nr_mods)
                  i  <-  1

    # Loop over lambda-, kappa- and col-grids
    for (lambda in lambda_grid) {
      for (kappa in kappa_grid)  {
        for (col_ind in col_grid) {

          # Set Col-Name
          col_name    <-  x_names[col_ind]

          # Append
          model_names_raw[i] <-  paste(col_name, lambda, kappa, sep = "-")
                          i  <-  i + 1
        }
      }
    }

    # Remove Objects
    rm(list = c("mu_var_tmp", "X", "mu_tmp_seq", "var_tmp_seq"))
  }

  ### 2.) TV-C-Model for Point Forecasts
  if (!is.null(F)) {

    # Set Variables and Indices to Subset Matrices mu_mat and variance_mat
    nr_preds      <-  ncol(F)
    start_cols    <-  1
    end_cols      <-  nr_preds
    mu_tmp_seq    <-  seq(1, 2 * nr_preds, 2)
    var_tmp_seq   <-  seq(2, 2 * nr_preds, 2)

    # Set Column-Index-Grid
    col_grid  <-  seq(1, nr_preds)

    # Number Models
    nr_mods  <- length(lambda_grid) * length(kappa_grid) * nr_preds

    # Set up Matrices
    max_length      <-  length(y)
    mu_mat_f        <-  matrix(NA, ncol = nr_mods, nrow = max_length)
    variance_mat_f  <-  matrix(NA, ncol = nr_mods, nrow = max_length)

    # Loop over Lambda- and Kappa-Grid
    for (i in seq_along(lambda_grid)) {
      for (j in seq_along(kappa_grid)) {

        # Set Lambda and Kappa
        lambda  <-  lambda_grid[i]
        kappa   <-  kappa_grid[j]

        # Set Up Backend for Parallel Processing
        cores   <-  n_cores
        cl      <-  parallel::makeCluster(cores, type = "PSOCK")
        parallel::clusterExport(cl = cl, varlist = c("y",
                                                     "F",
                                                     "max_length",
                                                     "init_length",
                                                     "lambda",
                                                     "kappa"),
                                envir = environment())
        parallel::clusterEvalQ(cl, library("hdflex"))

        # Parallelize with parLapply
        mu_var_tmp  <-  do.call("cbind", parallel::parLapply(cl, col_grid, function(j) {

          # Get Point-Forecasts Column and align lengths
          y_x          <-  cbind(y, F[, j])
          ts_y_x       <-  na.omit(y_x)
          drop_length  <-  max_length - nrow(ts_y_x)
          ts_length    <-  nrow(ts_y_x) - 1

          # Select Y-Variable
          y_var  <-  ts_y_x[,  1]

          # Select X-Variable
          x_var  <-  ts_y_x[, -1]

          # Initialize TV-C-Model
          init_results  <-  init_tvc_forecast(y_var, x_var, init_length)

          # Assign Init-Results
          theta    <-  init_results[[1]]
          cov_mat  <-  init_results[[2]]
          h        <-  init_results[[3]]

          # Update & Predict: Recursively Apply TV-C-Model-Function
          tvc_results  <-  tvc_model_loop(y_var,
                                          x_var,
                                          lambda,
                                          kappa,
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
    }

    ### Get Candidate Forecast Names
    # Get / Create Point Forecast Names
    if (!is.null(colnames(F))) {
      f_names  <-  colnames(F)
    } else {
      f_names  <-  paste0("F", as.character(seq_len(nr_preds)))
    }

    # Set up Vector
    model_names_f  <-  rep(NA, nr_mods)
          i        <-  1

    # Loop over Kappa Grid
    for (lambda in lambda_grid) {
      for (kappa in kappa_grid) {
        for (col_ind in col_grid) {

          # Set Colname
          col_name  <-  f_names[col_ind]

          # Append
          model_names_f[i]  <-  paste(col_name, lambda, kappa, sep = "-")
                  i         <-  i + 1
        }
      }
    }

    # Remove Objects
    rm(list = c("mu_var_tmp", "F", "mu_tmp_seq", "var_tmp_seq"))
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

  # Assign Model Names (-> Column Names)
  colnames(forecast_tvc)  <-  model_names_tvc
  colnames(variance_tvc)  <-  model_names_tvc

  # Return Results
  return(list(forecast_tvc, variance_tvc))
}