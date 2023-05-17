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
#' @return A list that contains
#' (1) a vector with the first moments (point forecasts) of the STSC-Model,
#' (2) a vector with the the second moments (variance) of the STSC-Model,
#' (3) a vector that contains the selected values for gamma,
#' (4) a vector that contains the selected values for psi and
#' (5) a matrix that indicates the selected signals for every point in time.
#' @export
#' @import parallel
#' @import checkmate
#' @import tidyverse
#' @importFrom stringr str_split
#' @importFrom dplyr lag
#' @importFrom roll roll_sum
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
  weights            <-  delta^(seq_len(len) - 1)
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

  # Count how often a predictor was selected
  pred_pockets  <-  sapply(1:n_preds, function(x) rowSums(cand_models[, which(preds[x] == cols_clean), drop = FALSE])) #nolint

  # Change to Selected / not Selected -> 0 vs. 1
  pred_pockets[pred_pockets > 1]  <-  1

  # Change Row and Column Names
  colnames(pred_pockets)  <-  preds

  # Return Results
  return(list(forecast_dsc,
              variance_dsc,
              val_gamma,
              val_psi,
              pred_pockets))
}