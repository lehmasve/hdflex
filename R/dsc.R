#' @name dsc
#' @title Generate Dynamic Subset Forecast Combinations
#' @description `dsc()` can be used to generate forecast combinations
#' from a set of candidate density forecasts. For each period,
#' `dsc()` selects a subset of predictive densities with highest ranks
#' regarding local predictive accuracy.
#' By using a discount factor gamma,
#' past predictive performance is exponentially down-weighted.
#' Both the identities of the forecasting models
#' that are used for building the combined forecast and
#' the subset sizes may vary over time based on the data.
#' If only one forecasting model is picked, the approach (temporarily)
#' collapses to pure model selection.
#' @param gamma_grid Numerical Vector.
#' Tuning parameter to exponentially down-weight past predictive performance.
#' @param psi_grid Integer Vector.
#' Tuning parameter that controls the size of the active subsets
#' (number of predictive densities that are combined).
#' @param y A matrix of dimension `T * 1` or Numeric Vector of length `T`
#' containing the observations of the target variable.
#' @param mu_mat A matrix with `T` rows containing
#' the point forecasts of each predictive density in each column.
#' @param var_mat A matrix with `T` rows containing
#' the variance of each predictive density in each column.
#' @param delta Double. Value denoting the discount factor used
#' to down-weight past predictive performance.
#' @param n_cores Integer. The number of cores to use for the estimation.
#' @return List that contains
#' (1) the point forecast (2) the variance of the Dynamic Subset Combination
#' and (3) a list that contains the column-index of the
#' selected density forecasts.
#' @export
#' @import parallel
#' @import checkmate
#' @examples
#' \donttest{
#'
#'   ### Simulate Data
#'   set.seed(123)
#'
#'   # Set Dimensions
#'   numb_obs   <-  500
#'   numb_pred  <-  50
#'   numb_forc  <-  10
#'
#'   # Create Random Target-Variable
#'   equity_premium  <-  rnorm(n = numb_obs, mean = 0, sd = 1)
#'
#'   # Create Random Text-Variables
#'   raw_preds            <-  replicate(numb_pred, sample(0:10,
#'                                                        numb_obs,
#'                                                        rep = TRUE), )
#'   raw_names            <-  paste0("X", as.character(seq_len(numb_pred)))
#'   colnames(raw_preds)  <-  raw_names
#'
#'   # Create Random Point Forecasts
#'   f_preds            <-  replicate(10, rnorm(n    = numb_obs,
#'                                              mean = 0,
#'                                              sd   = 0.5), )
#'   f_names            <-  paste0("F", as.character(seq_len(numb_forc)))
#'   colnames(f_preds)  <-  f_names
#'
#'   # Create Benchmark
#'   benchmark  <-  dplyr::lag(roll::roll_mean(equity_premium,
#'                                             width = length(equity_premium),
#'                                             min_obs = 1), n = 1)
#'
#'   # Specify TV-C-Parameter
#'   sample_length  <-  floor(numb_obs / 10)
#'   lambda_grid    <-  c(0.9995, 0.9999, 1.0000)
#'   kappa_grid     <-  c(0.94)
#'   n_cores        <-  1
#'
#'   # Apply TV-C-Function
#'   results  <-  hdflex::tvc(equity_premium,
#'                            raw_preds,
#'                            f_preds,
#'                            lambda_grid,
#'                            kappa_grid,
#'                            sample_length,
#'                            n_cores)
#'
#'   # Assign Results
#'   forecast_tvc      <-  results[[1]]
#'   variance_tvc      <-  results[[2]]
#'   model_names_tvc   <-  results[[3]]
#'
#'   # Cut Initialization-Period
#'   nr_drp              <-  dim(forecast_tvc)[1] -
#'                           dim(na.omit(forecast_tvc))[1] + 1
#'   sample_period_idx   <-  (nr_drp + sample_length):numb_obs
#'
#'   # Trim Objects
#'   sub_forecast_tvc    <-  forecast_tvc[sample_period_idx, , drop = FALSE]
#'   sub_variance_tvc    <-  variance_tvc[sample_period_idx, , drop = FALSE]
#'   sub_benchmark       <-  benchmark[sample_period_idx]
#'   sub_equity_premium  <-  equity_premium[sample_period_idx]
#'
#'   ##### Dynamic Subset Combination #####
#'   # Set DSC Parameter
#'   nr_mods     <-  length(model_names_tvc)
#'   gamma_grid  <-  c(0.9, 0.95, 0.99, 1)
#'   psi_grid    <-  c(1, 2, 3, 4, 5)
#'   delta       <-  0.9992
#'   n_cores     <-  1
#'
#'   # Apply DSC-Function
#'   results  <-  hdflex::dsc(gamma_grid,
#'                            psi_grid,
#'                            sub_equity_premium,
#'                            sub_forecast_tvc,
#'                            sub_variance_tvc,
#'                            delta,
#'                            n_cores)
#'
#'   # Assign Results
#'   sub_forecast_dsc    <-  results[[1]]
#'   sub_variance_dsc    <-  results[[2]]
#'   model_names_comb    <-  results[[3]]
#'   sub_chosen_para     <-  results[[4]]
#'   sub_models_idx      <-  results[[5]]
#'
#'   # Define Evaluation Period
#'   eval_period_idx     <-  50:length(sub_equity_premium)
#'
#'   # Trim Objects
#'   oos_equity_premium  <-  sub_equity_premium[eval_period_idx]
#'   oos_benchmark       <-  sub_benchmark[eval_period_idx]
#'   oos_forecast_dsc    <-  sub_forecast_dsc[eval_period_idx]
#'   oos_variance_dsc    <-  sub_variance_dsc[eval_period_idx]
#'   oos_models_idx      <-  lapply(seq_along(model_names_comb), function(i) {
#'                                       sub_models_idx[[i]][eval_period_idx]})
#'   oos_chosen_para     <-  sub_chosen_para[eval_period_idx, , drop = FALSE]
#'   oos_dates           <-  seq(as.Date("1989-12-01"),
#'                               by = "day",
#'                               length.out = length(eval_period_idx))
#'
#'   # Assign Names
#'   names(oos_forecast_dsc)  <-  oos_dates
#'   names(oos_variance_dsc)  <-  oos_dates
#'
#'   # Apply Statistial-Evaluation-Function
#'   eval_results  <-  eval_dsc(oos_equity_premium,
#'                              oos_benchmark,
#'                              oos_forecast_dsc,
#'                              oos_dates,
#'                              oos_chosen_para,
#'                              model_names_tvc,
#'                              oos_models_idx)
#'
#'   # Assign Results
#'   cw_t          <-  eval_results[[1]]
#'   oos_r2        <-  eval_results[[2]]
#'   csed          <-  eval_results[[3]]
#'   pred_pockets  <-  eval_results[[4]]
#' }

### Dynamic Subset Combination
dsc  <-  function(gamma_grid,
                  psi_grid,
                  y,
                  mu_mat,
                  var_mat,
                  delta,
                  n_cores) {

  ### Checkmate
  # Check if y is numeric vector without missing / infinite values
  checkmate::assertNumeric(y,
                           min.len = 1,
                           any.missing = FALSE,
                           finite = TRUE)

  # Check if mu_mat is numeric matrix and has the same number of observations as y
  checkmate::assertMatrix(mu_mat,
                          mode = "double",
                          any.missing = FALSE,
                          nrows = length(y))

  # Check if var_mat is numeric matrix and has the same number of observations as y
  checkmate::assertMatrix(var_mat,
                          mode = "double",
                          any.missing = FALSE,
                          nrows = length(y),
                          ncols = ncol(mu_mat))

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
                              upper = ncol(mu_mat),
                              min.len = 1,
                              any.missing = FALSE)

  # Check if var_mat has only non-negativ entries
  checkmate::qassert(var_mat,
                     c("M+[0,]"))

  # Check if delta is Numeric between 0 and 1
  checkmate::assertNumber(delta,
                          lower = exp(-10),
                          upper = 1,
                          na.ok = FALSE)

  # Check if n_cores is Integer bigger or equal to 1
  checkmate::assertInt(n_cores,
                       lower = 1)

  ### Parameter Grid
  # Create List with Parameter Grids to use as ID
  n_combs         <-  length(gamma_grid) * length(psi_grid)
  parameter_grid  <-  vector(mode = "list", length = n_combs)
  index           <-  1L
  for (i in seq_along(gamma_grid)) {
    for (j in seq_along(psi_grid))  {
      parameter_grid[[index]]  <-  c(i, j)
                      index    <-  index + 1
    }
  }

  ### 1) Forecast Combinations
  # Set Size
  n_models  <-  ncol(mu_mat)
  len       <-  nrow(mu_mat)

  # Set Up Backend for Parallel Processing
  cores  <-  n_cores
  cl     <-  parallel::makeCluster(cores, type = "PSOCK")
  parallel::clusterExport(cl = cl, varlist = c("parameter_grid",
                                               "n_models",
                                               "n_combs",
                                               "len",
                                               "y",
                                               "mu_mat",
                                               "var_mat"),
                            envir = environment())

  # Parallelize with parLapply
  idx      <-  seq_along(parameter_grid)
  dsc_tmp  <-  parallel::parLapply(cl, idx, function(i) {

    # Set Gamma and Psi
    gam    <-  gamma_grid[parameter_grid[[i]][1]]
    psi    <-  psi_grid[parameter_grid[[i]][2]]

    # Create Initial Weights
    weights  <-  init_dsc(n_models)

    # Loop DSC
    dsc_results  <-  dsc_loop(weights,
                              gam,
                              psi,
                              y,
                              mu_mat,
                              var_mat,
                              n_combs,
                              len,
                              n_models)
    # Return Results
    return(list(cbind(dsc_results[[1]],
                      dsc_results[[2]],
                      dsc_results[[3]]),
                      dsc_results[[4]]))})

  # Stop Cluster
  parallel::stopCluster(cl)

  # Get Results
  tmp             <-  lapply(dsc_tmp, "[[", 1)
  forecasts_comb  <-  sapply(tmp, function(x) x[, 1])
  variances_comb  <-  sapply(tmp, function(x) x[, 2])
  ln_scores       <-  sapply(tmp, function(x) x[, 3])
  models_idx      <-  lapply(dsc_tmp, "[[", 2)

  # Remove Objects
  rm(list = c("dsc_tmp", "tmp"))

  ### 2) Get Model Names
  # Set up Vector
  model_names  <-  rep(NA, n_combs)

  # Loop over Grid
  for (i in as.integer(seq_along(parameter_grid))) {

        # Set Gamma, Psi and Equal-Weight
        gam  <-  gamma_grid[parameter_grid[[i]][1]]
        psi  <-  psi_grid[parameter_grid[[i]][2]]

        # Create Model Name
        mod_name  <-  paste("gamma", gam, "psi", psi, sep = "_")

        # Assign Name
        model_names[i]  <-  mod_name
  }

  ### 3) Tune Parameter and choose Subset Combination
  # Exponentially Discounted Sum of Predictive Log-Likelihoods
  weights            <-  delta^(seq(1:len) - 1)
  cum_ln_scores_lag  <-  dplyr::lag(roll::roll_sum(ln_scores,
                                                    weights = rev(weights),
                                                    width = len, min_obs = 1),
                                    n = 1L)

  # Tune Gamma, Psi and EW
  chosen_parameters  <-  matrix(FALSE, ncol = n_combs, nrow = len)
  chosen_parameters[cbind(seq_len(len), max.col(cum_ln_scores_lag, "first"))]  <-  TRUE 

  # Get chosen Dynamic Subset Combination Forecast
  forecast_dsc  <-  rowSums(forecasts_comb * chosen_parameters)
  variance_dsc  <-  rowSums(variances_comb * chosen_parameters)

  # Assign Model Names
  colnames(chosen_parameters)  <-  model_names

  # Return Results
  return(list(forecast_dsc,
              variance_dsc,
              model_names,
              chosen_parameters,
              models_idx))
}