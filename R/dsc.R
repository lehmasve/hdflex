#' @name dsc
#' @title Generate dynamic subset forecast combinations
#' @description `dsc()` can be used to generate forecast combinations
#' from a set of candidate density forecasts. For each period,
#' `dsc()` selects a subset of predictive densities with highest ranks
#' regarding (local) predictive accuracy.
#' Both the identities of the candidate forecasts
#' that are used for building the combined forecast and
#' the subset sizes may vary over time based on the data.
#' If only one candidate forecast is picked, the approach (temporarily)
#' collapses to pure model selection.
#' @param gamma_grid A numerical vector that contains discount factors
#' to exponentially down-weight the past predictive performance
#' of the candidate forecasts.
#' @param psi_grid An integer vector that controls
#' the (possible) sizes of the active subsets.
#' @param y A matrix of dimension `T * 1` or numeric vector of length `T`
#' containing the observations of the target variable.
#' @param mu_mat A matrix with `T` rows containing
#' the first moment of each predictive density in each column.
#' @param var_mat A matrix with `T` rows containing
#' the second moment of each predictive density in each column.
#' @param delta A numeric value denoting the discount factor used
#' to down-weight the past predictive performance of the subset combinations.
#' @param n_cores An integer that denotes the number of CPU-cores
#' used for the computational estimation.
#' @return A list that contains:
#' * (1) a vector with the first moments (point forecasts) of the STSC-Model,
#' * (2) a vector with the the second moments (variance) of the STSC-Model,
#' * (3) a vector that contains the selected values for gamma,
#' * (4) a vector that contains the selected values for psi and
#' * (5) a matrix that indicates the selected signals for every point in time.
#' @export
#' @seealso \url{https://github.com/lehmasve/hdflex#readme}
#' @author Philipp Adämmer, Sven Lehmann, Rainer Schüssler
#' @references
#' Beckmann, J., Koop, G., Korobilis, D., and Schüssler, R. A. (2020) "Exchange rate predictability and dynamic bayesian learning."
#' \emph{Journal of Applied Econometrics}, 35 (4): 410–421.
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
#' Del Negro, M., Hasegawa, R. B., and Schorfheide, F. (2016) "Dynamic prediction pools: An investigation of financial frictions and forecasting performance."
#' \emph{Journal of Econometrics}, 192 (2): 391–405.
#'
#' West, M. and Harrison, J. (1997) "Bayesian forecasting and dynamic models"
#' \emph{Springer}, 2nd edn.
#' @import parallel
#' @import checkmate
#' @importFrom stringr str_split
#' @importFrom dplyr lag
#' @importFrom roll roll_sum
#' @examples
#' \donttest{
#'
#' # See example for tvc().
#'
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

  ### 1) Compute (all) Subset Combinations
  # Set Size
  n_models  <-  ncol(mu_mat)
  len       <-  nrow(mu_mat)

  # Set up Backend for Parallel Processing
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
  idx        <-  seq_along(parameter_grid)
  dsc_tmp    <-  parallel::parLapply(cl, idx, function(i) {

    # Set Gamma and Psi
    gamma    <-  gamma_grid[parameter_grid[[i]][1]]
    psi      <-  psi_grid[parameter_grid[[i]][2]]

    # Create Initial Weights
    weights  <-  init_dsc(n_models)

    # Loop over DSC-Function
    dsc_results  <-  dsc_loop(weights,
                              gamma,
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

  # Assign Results
  tmp             <-  lapply(dsc_tmp, "[[", 1)
  forecasts_comb  <-  sapply(tmp, function(x) x[, 1])
  variances_comb  <-  sapply(tmp, function(x) x[, 2])
  ln_scores       <-  sapply(tmp, function(x) x[, 3])
  models_idx      <-  lapply(dsc_tmp, "[[", 2)

  # Remove Objects
  rm(list = c("dsc_tmp", "tmp"))

  ### 2) Create Combination Names
  #  Set up Vector & Loop over Grid
  model_names_comb  <-  rep(NA, n_combs)
  for (i in as.integer(seq_along(parameter_grid))) {

        # Set Gamma and Psi
        gamma  <-  gamma_grid[parameter_grid[[i]][1]]
        psi    <-  psi_grid[parameter_grid[[i]][2]]

        # Create Model Name
        mod_name  <-  paste("gamma", gamma, "psi", psi, sep = "_")

        # Assign Name
        model_names_comb[i]  <-  mod_name
  }

  ### 3) Compute Dynamic Subset Combination
  ### -> select subset combination for each point in time
  # Compute exponentially discounted sum of predictive log-likelihoods (DPLL)
  weights            <-  delta^(seq_len(len)) #  - 1)
  cum_ln_scores_lag  <-  dplyr::lag(roll::roll_sum(ln_scores,
                                                   weights = rev(weights),
                                                   width = len, min_obs = 1), n = 1L) #nolint

  # Select highest DPLL for each point in time
  chosen_parameter  <-  matrix(FALSE, ncol = n_combs, nrow = len, dimnames = list(NULL, model_names_comb))  #nolint
  chosen_parameter[cbind(seq_len(len), max.col(cum_ln_scores_lag, "first"))]  <-  TRUE #nolint

  # Set first Subset Combination deterministically
  chosen_parameter[1, n_combs]  <-  TRUE

  # Compute DSC-Forecast and DSC-Variance
  forecast_dsc  <-  rowSums(forecasts_comb * chosen_parameter)
  variance_dsc  <-  rowSums(variances_comb * chosen_parameter)


  ### 4) Diagnosis: Get selected values for gamma & psi
  ### and the names of the candidate forecasts
  # Get selected gamma & psi value for each point in time
  gamma_psi  <-  apply(chosen_parameter == TRUE, 1, function(x) model_names_comb[x])                     #nolint
  val_gamma  <-  unname(sapply(gamma_psi, function(x) as.numeric(stringr::str_split(x, "_")[[1]][2])))   #nolint
  val_psi    <-  unname(sapply(gamma_psi, function(x) as.numeric(stringr::str_split(x, "_")[[1]][4])))   #nolint

  # Get / Set the Candidate Forecast Names
  if (!is.null(colnames(mu_mat))) {
      model_names_tvc  <-  colnames(mu_mat)
  } else {
      model_names_tvc  <-  paste0("C_", as.character(seq_len(n_models)))
  }


  # Get index of the selected Subset Combination at every point in time
  ind  <-  which(chosen_parameter, arr.ind = TRUE, useNames = TRUE)
  ind  <-  ind[order(ind[, "row"]), ]

  # Get the Candidate Forecast(s) in the Subset Combination
  tmp  <-  lapply(seq_len(nrow(ind)), function(x) model_names_tvc[models_idx[[ind[x, 2]]][[ind[x, 1]]] + 1]) #nolint

  # Transform to 0 / 1 - Matrix
  cand_models  <-  matrix(0, ncol = length(model_names_tvc),
                             nrow = nrow(chosen_parameter),
                             dimnames = list(NULL, model_names_tvc))
  for (i in seq_len(nrow(chosen_parameter))) {
      col_idx  <-  match(tmp[[i]], model_names_tvc)
      cand_models[i, col_idx]  <-  1
  }

  # Clean Column Names
  cols_clean  <-  sapply(stringr::str_split(model_names_tvc, "-"), `[`, 1)
  preds       <-  unique(cols_clean)
  n_preds     <-  length(preds)

  # Count how often a signal was selected
  chosen_signals  <-  sapply(1:n_preds, function(x) rowSums(cand_models[, which(preds[x] == cols_clean), drop = FALSE])) #nolint

  # Change to Selected / not Selected -> 0 vs. 1
  chosen_signals[chosen_signals > 1]  <-  1

  # Change Row and Column Names
  colnames(chosen_signals)  <-  preds

  # Return Results
  return(list(forecast_dsc,
              variance_dsc,
              val_gamma,
              val_psi,
              chosen_signals))
}