#' @name summary_stsc
#' @title Statistical summary of the STSC-results
#' @description `summary_stsc()` returns a statistical summary
#' of the results from dsc(). It provides statistical measures
#' such as Clark-West-Statistic, OOS-R2, Mean-Squared-Error and
#' Cumulated Sum of Squared-Error-Differences.
#' @param oos_y A matrix of dimension `T * 1` or numeric vector of length `T`
#' containing the out-of-sample observations of the target variable.
#' @param oos_benchmark A matrix of dimension `T * 1` or
#' numeric vector of length `T` containing the
#' out-of-sample forecasts of an arbitrary benchmark
#' (i.e. prevailing historical mean).
#' @param oos_forecast_stsc A matrix of dimension `T * 1`
#' or numeric vector of length `T` containing the
#' out-of-sample forecasts of dsc().
#' @export
#' @return List that contains
#' (1) the Clark-West-Statistic,
#' (2) the Out-of-Sample R2,
#' (3) a vector with the CSSED between the STSC-Forecast and the benchmark and
#' (4) a list with the MSE of the STSC-Model and the benchmark.
#' @import checkmate
#' @importFrom stats t.test
#' @examples
#' \donttest{
#'
#'  ### Set Seed
#'  set.seed(123)
#'
#'  ### Simulate Coefficients
#'  # Set Dimensions
#'  numb_obs      <-  500
#'  numb_noise    <-  499
#'  numb_signals  <-  numb_noise + 1
#'  numb_forc     <-  10
#'  full_dates    <-  as.character(seq(as.Date("1989-12-01"),
#'                                     by = "day",
#'                                     length.out = numb_obs))
#'
#'  # Set up Coefficient-Matrix
#'  theta  <-  matrix(NA, nrow = numb_obs, ncol = numb_signals)
#'
#'  # Loop over Time
#'  for(t in seq_len(numb_obs)) {
#'
#'    ### Theta: Abrupt Change
#'    theta[t, 1]  <-  ifelse((t > 200 & t < 450), -0.40, 0.50)
#'
#'    ### Noise Predictors
#'    theta[t, -1] <-  0
#'  }
#'
#'  ### Draw Simple Signals
#'  raw_signals            <-  replicate(numb_signals, rnorm(numb_obs, 0, 1))
#'  colnames(raw_signals)  <-  paste0("X", as.character(seq_len(numb_signals)))
#'
#'  ### Draw Noise
#'  eps          <-  rnorm(numb_obs, 0, 0.25)
#'
#'  ### Compute Target Variable
#'  y            <-  as.matrix(rowSums(cbind(theta * raw_signals, eps)))
#'  colnames(y)  <-  "y"
#'
#'  # Create 'External' Forecasts
#'  f_signals            <-  y[,] + replicate(numb_forc, rnorm(numb_obs, 0, 1.0))
#'  colnames(f_signals)  <-  paste0("F", as.character(seq_len(numb_forc)))
#'
#'  # Create Benchmark (Prevailing Historical Mean)
#'  benchmark  <-  dplyr::lag(roll::roll_mean(y,
#'                                            width = length(y),
#'                                            min_obs = 1), n = 1)
#'  ### TV-C-Forecasts
#'  # Set TV-C-Parameter
#'  sample_length  <-  floor(numb_obs / 10)
#'  lambda_grid    <-  c(0.98, 0.99, 1.00)
#'  kappa_grid     <-  c(0.97)
#'  n_cores        <-  1
#'
#'  # Apply TV-C-Function
#'  results  <-  hdflex::tvc(y,
#'                           raw_signals,
#'                           f_signals,
#'                           lambda_grid,
#'                           kappa_grid,
#'                           sample_length,
#'                           n_cores)
#'
#'  # Assign Results
#'  forecast_tvc      <-  results[[1]]
#'  variance_tvc      <-  results[[2]]
#'  model_names_tvc   <-  colnames(forecast_tvc)
#'
#'  # Cut Initialization-Period
#'  sample_period_idx   <-  (sample_length+1):numb_obs
#'  sub_forecast_tvc    <-  forecast_tvc[sample_period_idx, , drop = FALSE]
#'  sub_variance_tvc    <-  variance_tvc[sample_period_idx, , drop = FALSE]
#'  sub_benchmark       <-  benchmark[sample_period_idx]
#'  sub_y               <-  y[sample_period_idx]
#'  sub_dates           <-  full_dates[sample_period_idx]
#'
#'  ##### Dynamic Subset Combination #####
#'  # Set DSC-Parameter
#'  nr_mods     <-  ncol(forecast_tvc)
#'  gamma_grid  <-  c(0.8, 0.9, 0.95, 0.99, 1)
#'  psi_grid    <-  c(1:25)
#'  delta       <-  0.99
#'  n_cores     <-  1
#'
#'  # Apply DSC-Function
#'  results  <-  hdflex::dsc(gamma_grid,
#'                           psi_grid,
#'                           sub_y,
#'                           sub_forecast_tvc,
#'                           sub_variance_tvc,
#'                           delta,
#'                           n_cores)
#'
#'   # Assign Results
#'   sub_forecast_stsc    <-  results[[1]]
#'   sub_variance_stsc    <-  results[[2]]
#'   sub_chosen_gamma     <-  results[[3]]
#'   sub_chosen_psi       <-  results[[4]]
#'   sub_pred_pockets     <-  results[[5]]
#'
#'   # Define Evaluation Period
#'   eval_period_idx      <-  51:length(sub_y)
#'
#'   # Keep only OOS-Period
#'   oos_y                <-  sub_y[eval_period_idx]
#'   oos_benchmark        <-  sub_benchmark[eval_period_idx]
#'   oos_forecast_stsc    <-  sub_forecast_stsc[eval_period_idx]
#'   oos_variance_stsc    <-  sub_variance_stsc[eval_period_idx]
#'   oos_chosen_gamma     <-  sub_chosen_gamma[eval_period_idx]
#'   oos_chosen_psi       <-  sub_chosen_psi[eval_period_idx]
#'   oos_pred_pockets     <-  sub_pred_pockets[eval_period_idx, , drop = FALSE]
#'   oos_length           <-  length(eval_period_idx)
#'   oos_dates            <-  sub_dates[eval_period_idx]
#'
#'   # Assign Names
#'   names(oos_forecast_stsc)   <-  oos_dates
#'   names(oos_variance_stsc)   <-  oos_dates
#'   names(oos_chosen_gamma)    <-  oos_dates
#'   names(oos_chosen_psi)      <-  oos_dates
#'   rownames(oos_pred_pockets) <-  oos_dates
#'
#'  # Apply Statistial-Evaluation-Function
#'  summary_results  <-  summary_stsc(oos_y,
#'                                    oos_benchmark,
#'                                    oos_forecast_stsc)
#'
#'  # Assign Results
#'  cw_t           <-  summary_results[[1]]
#'  oos_r2         <-  summary_results[[2]]
#'  cssed          <-  summary_results[[3]]
#'  mse            <-  summary_results[[4]]
#' }

### Evaluation
summary_stsc  <-  function(oos_y,
                          oos_benchmark,
                          oos_forecast_stsc) {

    ### Checkmate
    # Check if oos_y is numeric vector without missing / infinite values
    checkmate::assertNumeric(oos_y,
                             min.len = 1,
                             any.missing = FALSE,
                             finite = TRUE)

    # Check if oos_benchmark is numeric vector without missing / infinite values
    checkmate::assertNumeric(oos_benchmark,
                             len = length(oos_y),
                             any.missing = FALSE,
                             finite = TRUE)

    # Check if oos_forecast_stsc is numeric vector without missing / infinite values
    checkmate::assertNumeric(oos_forecast_stsc,
                             len = length(oos_y),
                             any.missing = FALSE,
                             finite = TRUE)

    ### 1) Clark-West and OOS-R2
    # Squared Error (SE) Target-Variable vs. Benchmark
    se_benchmark      <-  (oos_y - oos_benchmark) ** 2

    # SE Target-Variable vs. Dynamic Subset Combination
    se_stsc            <-  (oos_y - oos_forecast_stsc) ** 2

    # SE Benchmark vs. Dynamic Subset Combination
    se_benchmark_stsc  <-  (oos_benchmark - oos_forecast_stsc) ** 2

    # SED Benchmark vs. Dynamic Subset Combination
    sed_stsc  <-  se_benchmark - se_stsc

    # Cumulated SED
    cssed    <-  cumsum(sed_stsc)

    # Clark-West-Statistic
    cw_statistic  <-  se_benchmark - se_stsc + se_benchmark_stsc
    cw_t          <-  stats::t.test(cw_statistic,
                                    mu = 0,
                                    alternative = "greater")$statistic

    # Out-of-Sample R2
    oos_r2  <-  1 - sum(se_stsc) / sum(se_benchmark)

    # Mean-Squared Error
    mse_stsc       <-  mean(se_stsc)
    mse_benchmark  <-  mean(se_benchmark)
    mse            <-  list(STSC = mse_stsc, Benchmark = mse_benchmark)

    # Return Results
    return(list(cw_t,
                oos_r2,
                cssed,
                mse))
}