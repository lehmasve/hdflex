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
#' @return List that contains:
#' * (1) the Clark-West-Statistic,
#' * (2) the Out-of-Sample R2,
#' * (3) a vector with the CSSED between the STSC-Forecast and the benchmark and
#' * (4) a list with the MSE of the STSC-Model and the benchmark.
#' @seealso \url{https://github.com/lehmasve/hdflex#readme}
#' @author Philipp Adämmer, Sven Lehmann, Rainer Schüssler
#' @references
#' Clark, T. E. and West, K. D. (2007) "Approximately normal tests for equal predictive accuracy in nested models."
#' \emph{Journal of Econometrics}, 138 (1): 291–311.
#'
#' @import checkmate
#' @importFrom stats t.test
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
#'    ### Part 1: TV-C Model ###
#'    # Set Target Variable
#'    y  <-  dataset[,  1, drop = FALSE]
#'
#'    # Set 'Simple' Signals
#'    X  <-  dataset[, 2:442, drop = FALSE]
#'
#'    # Set External Point Forecasts (Koop & Korobilis 2023)
#'    F  <-  dataset[, 443:458, drop = FALSE]
#'
#'    # Set TV-C-Parameter
#'    sample_length  <-  4 * 5
#'    lambda_grid    <-  c(0.90, 0.95, 1)
#'    kappa_grid     <-  0.98
#'    n_cores        <-  1
#'
#'    # Apply TV-C-Function
#'    results  <-  hdflex::tvc(y,
#'                             X,
#'                             F,
#'                             lambda_grid,
#'                             kappa_grid,
#'                             sample_length,
#'                             n_cores)
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