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
#'\donttest{
#'
#' # See example for tvc().
#'
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