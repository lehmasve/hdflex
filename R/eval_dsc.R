#' @name eval_dsc
#' @title Evaluate Dynamic Subset Combination Results
#' @description `eval_dsc()` can be used to evaluate the results from [dsc()].
#' It provides statistical measures such as Clark-West-Statistic or the OOS-R2,
#' visualizes for the predictive performance as well as the selected parameters
#' and returns a matrix which indicates the selected predictors
#' at every point in time.
#' @param oos_y y A matrix of dimension `T * 1` or Numeric Vector of length `T`
#' containing the out-of-sample observations of the target variable.
#' @param oos_benchmark A matrix of dimension `T * 1` or
#' Numeric Vector of length `T` containing the
#' out-of-sample forecasts of a arbitrary benchmark.
#' @param oos_forecast_dsc A matrix of dimension `T * 1`
#' or Numeric Vector of length `T` containing the
#' out-of-sample forecasts of the Dynamic Subset Combination.
#' @param oos_dates A vector that contains all ouf-of-sample dates.
#' @param oos_chosen_parameter Matrix from 'dsc' that contains the
#' chosen (tuned) parameter from the Dynamic Subset Combination.
#' @param model_names_tvc Character Vector from [tvp()] that contains
#' all model names from the density forecasts.
#' @param oos_models_idx List from [dsc()] that contains
#' the column-indices of the selected density forecasts.
#' @export
#' @return List that contains
#' (1) the Clark-West-Statistic,
#' (2) the Out-of-Sample R2,
#' (3) the CSED between the DSC-Forecast and the benchmark,
#' (4) a boolean matrix that indicates when a predictor was chosen from the DSC
#' (5) and a CSED Plot.
#' @import checkmate
#' @importFrom graphics matplot
#' @importFrom stringr str_split
#' @importFrom stats t.test
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

### Evaluation
eval_dsc  <-  function(oos_y,
                       oos_benchmark,
                       oos_forecast_dsc,
                       oos_dates,
                       oos_chosen_parameter,
                       model_names_tvc,
                       oos_models_idx) {

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

    # Check if oos_forecast_dsc is numeric vector without missing / infinite values
    checkmate::assertNumeric(oos_forecast_dsc,
                             len = length(oos_y),
                             any.missing = FALSE,
                             finite = TRUE)

    # Check if oos_dates is numeric vector without missing / infinite values
    checkmate::assertDate(oos_dates,
                          len = length(oos_y),
                          any.missing = FALSE)

    # Check if oos_chosen_parameter is Logical Matrix
    checkmate::assertMatrix(oos_chosen_parameter,
                            mode = "logical",
                            any.missing = FALSE,
                            nrows = length(oos_y))

    # Check if model_names_tvc is Character Vector
    checkmate::assertCharacter(model_names_tvc,
                               unique = TRUE)

    # Check if oos_models_idx is a List
    checkmate::assertList(oos_models_idx,
                          types = "list")

    ### 1) Clark-West and OOS-R2
    # Squared Error (SE) Target-Variable vs. Benchmark
    se_benchmark      <-  (oos_y - oos_benchmark) ** 2

    # SE Target-Variable vs. Dynamic Subset Combination
    se_dsc            <-  (oos_y - oos_forecast_dsc) ** 2

    # SE Benchmark vs. Dynamic Subset Combination
    se_benchmark_dsc  <-  (oos_benchmark - oos_forecast_dsc) ** 2

    # SED Benchmark vs. Dynamic Subset Combination
    sed_dsc  <-  se_benchmark - se_dsc

    # Cumulated SED
    csed     <-  cumsum(sed_dsc)

    # Clark-West-Statistic
    cw_statistic  <-  se_benchmark - se_dsc + se_benchmark_dsc
    cw_t          <-  stats::t.test(cw_statistic,
                                    mu = 0,
                                    alternative = "greater")$statistic

    # Out-of-Sample R2
    oos_r2  <-  1 - sum(se_dsc) / sum(se_benchmark)

    # Performance Plot
    # plot_performance  <-  graphics::matplot(x = as.Date(oos_dates),
    #                                         y = csed,
    #                                         type = "l",
    #                                         xlab = "Time",
    #                                         ylab = "Cumulated SED",
    #                                         main = "Model Performance",
    #                                         verbose = FALSE)

    ### 2) Get the selected Predictors
    # Assign Dates
    rownames(oos_chosen_parameter) <-  as.character(oos_dates)

    # Get Index of the selected Forecast Combination at every point in time
    ind  <-  which(oos_chosen_parameter, arr.ind = TRUE, useNames = TRUE)
    ind  <-  ind[order(row.names(ind)), ]

    # Get the "active" Univariate Models for every point in time
    oos_uni_models_chr  <-  lapply(seq_len(nrow(ind)), function(x) model_names_tvc[oos_models_idx[[ind[x, 2]]][[ind[x, 1]]] + 1]) #nolint

    # Transform to 0 / 1 - Matrix
    oos_uni_models  <-  matrix(0, ncol = length(model_names_tvc),
                               nrow = nrow(oos_chosen_parameter),
                               dimnames = list(oos_dates, model_names_tvc))
    for (i in seq_len(nrow(oos_chosen_parameter))) {
        col_idx  <-  match(oos_uni_models_chr[[i]], model_names_tvc)
        oos_uni_models[i, col_idx]  <-  1
    }

    # Clean Column Names
    cols_clean  <-  sapply(stringr::str_split(model_names_tvc, "-"), `[`, 1)
    preds       <-  unique(cols_clean)
    n_preds     <-  length(preds)

    # Count how often a predictor was selected
    pred_pockets  <-  sapply(1:n_preds, function(x) rowSums(oos_uni_models[, which(preds[x] == cols_clean), drop = FALSE])) #nolint

    # Change to Selected / not Selected -> 0 vs. 1
    pred_pockets[pred_pockets > 1]  <-  1

    # Change Row and Column Names
    dimnames(pred_pockets)  <-  list(as.character(oos_dates), preds)

    # Return Results
    return(list(cw_t,
                oos_r2,
                csed,
                pred_pockets))
}
