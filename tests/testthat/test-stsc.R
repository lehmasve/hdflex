###########################################################
### Simulate Data
# Set Dimensions
n_obs <- 500
n_sigs <- 90

### Simulate Data
# Generate Covariates
X <- matrix(rnorm(n_obs * n_sigs), nrow = n_obs, ncol = n_sigs)

# Generate Beta-Coefficients
n_relevant <- 10
beta <- runif(n_relevant, -1.0, 1.0)

# Compute f(x)
f_x <- X[, seq(n_relevant)] %*% beta

# Generate Error-Term
eps <- rnorm(n_obs)

# Calculate Response
y <- as.matrix(f_x + eps, ncol = 1)
returns <- as.matrix(exp(f_x + eps), ncol = 1)

# F-Signals
Ext_F <- matrix(rep(y, 10), nrow = n_obs, ncol = 10) + rnorm(n_obs * 10)

# Add Names
colnames(X) <- paste0("X", seq_len(n_sigs))
colnames(y) <- "response"
colnames(Ext_F) <- paste0("F", seq_len(10))

###########################################################
### STSC Parameter
# TV-C-Parameter
init <- 10
lambda_grid <- c(0.95, 1.00)
kappa_grid <- c(0.95, 0.97)
bias <- TRUE

# Set DSC-Parameter
gamma_grid <- c(0.9, 0.95, 1)
psi_grid <- c(1:10)
delta <- 0.95
burn_in <- 5
burn_in_dsc <- 10
metric <- 5
equal_weight <- TRUE
incl <- NULL

# Parallel-Parameter
parallel <- FALSE
n_threads <- 1

# Set Portfolio-Parameter
portfolio_params <- c(3, 0, 2)

###########################################################
### Test STSC
test_that("Test whether the STSC-Function works", {

  apply_stsc <- function(y, metric) {
    stsc(y,
         X,
         Ext_F,
         init,
         lambda_grid,
         kappa_grid,
         bias,
         gamma_grid,
         psi_grid,
         delta,
         burn_in,
         burn_in_dsc,
         metric,
         equal_weight,
         incl,
         parallel,
         n_threads,
         portfolio_params)
  }

  check_results <- function(results, y) {
    # List Contains three Elements
    expect_equal(length(results), 3)

    # Forecasts List Contains three Elements
    expect_equal(length(results$Forecasts), 3)

    # Point Forecasts
    expect_numeric(results$Forecasts$Point_Forecasts, len = n_obs, finite = TRUE)

    # Variance Forecasts
    expect_numeric(results$Forecasts$Variance_Forecasts, len = n_obs, lower = 0, finite = TRUE)

    # Realization
    expect_equal(results$Forecasts$Realization, y)

    # Tuning Parameters List Contains five Elements
    expect_equal(length(results$Tuning_Parameters), 5)

    # Gamma-Vector
    expect_numeric(results$Tuning_Parameters$Gamma, len = n_obs, lower = min(gamma_grid), upper = max(gamma_grid), finite = TRUE)

    # Psi-Vector
    expect_numeric(results$Tuning_Parameters$Psi, len = n_obs, lower = min(psi_grid), upper = max(psi_grid), finite = TRUE)

    # Signals
    expect_matrix(results$Tuning_Parameters$Signals, mode = "integerish", nrows = n_obs, ncols = (ncol(X) + ncol(Ext_F)))

    # Lambda-Vector
    expect_matrix(results$Tuning_Parameters$Lambda, mode = "integerish", nrows = n_obs, ncols = length(lambda_grid))

    # Kappa-Vector
    expect_matrix(results$Tuning_Parameters$Kappa, mode = "integerish", nrows = n_obs, ncols = length(kappa_grid))

    # Model List Contains 12 Elements
    expect_equal(length(results$Model), 12)

    # Lambda Grid
    expect_equal(results$Model$Lambda_grid, lambda_grid)

    # Kappa Grid
    expect_equal(results$Model$Kappa_grid, kappa_grid)

    # Gamma Grid
    expect_equal(results$Model$Gamma_grid, gamma_grid)

    # Psi Grid
    expect_equal(results$Model$Psi_grid, psi_grid)

    # Delta
    expect_equal(results$Model$Delta, delta)

    # Init
    expect_equal(results$Model$Init, init)

    # Burn-in
    expect_equal(results$Model$Burn_in, burn_in)

    # Burn-in DSC
    expect_equal(results$Model$Burn_in_dsc, burn_in_dsc)

    # Metric
    expect_numeric(results$Model$Metric, len = 1, lower = 1, upper = 5)

    # Equal Weight
    expect_equal(results$Model$Equal_weight, equal_weight)

    # Bias
    expect_equal(results$Model$Bias, bias)

    # Incl
    expect_equal(results$Model$Incl, incl)
  }

  # Apply STSC-Function
  results1 <- apply_stsc(y, 1)
  results2 <- apply_stsc(y, 2)
  results3 <- apply_stsc(y, 3)
  results4 <- apply_stsc(returns, 4)
  results5 <- apply_stsc(y, 5)

  # Check results
  check_results(results1, y)
  check_results(results2, y)
  check_results(results3, y)
  check_results(results4, returns)
  check_results(results5, y)
})

###########################################################
### Test STSC-Parallel
test_that("Test whether the STSC-Parallel-Function works", {

  apply_stsc <- function(y, metric) {
    stsc(y,
         X,
         Ext_F,
         init,
         lambda_grid,
         kappa_grid,
         bias,
         gamma_grid,
         psi_grid,
         delta,
         burn_in,
         burn_in_dsc,
         metric,
         equal_weight,
         incl,
         TRUE,
         1,
         portfolio_params)
  }

  check_results <- function(results, y) {
    # List Contains three Elements
    expect_equal(length(results), 3)

    # Forecasts List Contains three Elements
    expect_equal(length(results$Forecasts), 3)

    # Point Forecasts
    expect_numeric(results$Forecasts$Point_Forecasts, len = n_obs, finite = TRUE)

    # Variance Forecasts
    expect_numeric(results$Forecasts$Variance_Forecasts, len = n_obs, lower = 0, finite = TRUE)

    # Realization
    expect_equal(results$Forecasts$Realization, y)

    # Tuning Parameters List Contains five Elements
    expect_equal(length(results$Tuning_Parameters), 5)

    # Gamma-Vector
    expect_numeric(results$Tuning_Parameters$Gamma, len = n_obs, lower = min(gamma_grid), upper = max(gamma_grid), finite = TRUE)

    # Psi-Vector
    expect_numeric(results$Tuning_Parameters$Psi, len = n_obs, lower = min(psi_grid), upper = max(psi_grid), finite = TRUE)

    # Signals
    expect_matrix(results$Tuning_Parameters$Signals, mode = "integerish", nrows = n_obs, ncols = (ncol(X) + ncol(Ext_F)))

    # Lambda-Vector
    expect_matrix(results$Tuning_Parameters$Lambda, mode = "integerish", nrows = n_obs, ncols = length(lambda_grid))

    # Kappa-Vector
    expect_matrix(results$Tuning_Parameters$Kappa, mode = "integerish", nrows = n_obs, ncols = length(kappa_grid))

    # Model List Contains 12 Elements
    expect_equal(length(results$Model), 12)

    # Lambda Grid
    expect_equal(results$Model$Lambda_grid, lambda_grid)

    # Kappa Grid
    expect_equal(results$Model$Kappa_grid, kappa_grid)

    # Gamma Grid
    expect_equal(results$Model$Gamma_grid, gamma_grid)

    # Psi Grid
    expect_equal(results$Model$Psi_grid, psi_grid)

    # Delta
    expect_equal(results$Model$Delta, delta)

    # Init
    expect_equal(results$Model$Init, init)

    # Burn-in
    expect_equal(results$Model$Burn_in, burn_in)

    # Burn-in DSC
    expect_equal(results$Model$Burn_in_dsc, burn_in_dsc)

    # Metric
    expect_numeric(results$Model$Metric, len = 1, lower = 1, upper = 5)

    # Equal Weight
    expect_equal(results$Model$Equal_weight, equal_weight)

    # Bias
    expect_equal(results$Model$Bias, bias)

    # Incl
    expect_equal(results$Model$Incl, incl)
  }

  # Apply STSC-Function
  results1 <- apply_stsc(y, 1)
  results2 <- apply_stsc(y, 2)
  results3 <- apply_stsc(y, 3)
  results4 <- apply_stsc(returns, 4)
  results5 <- apply_stsc(y, 5)

  # Check results
  check_results(results1, y)
  check_results(results2, y)
  check_results(results3, y)
  check_results(results4, returns)
  check_results(results5, y)
})

#############################################################
### Test for same results between STSC and STSC-Parallel for different metrics
test_that("Test whether the STSC-Function and STSC-Parallel-Function return the same results", {

  for (m in seq(5)) {

    # Use returns instead of y if m == 4
    y_input <- if (m == 4) returns else y

    # Apply STSC-Function
    results <- stsc(y_input,
                    X,
                    Ext_F,
                    init,
                    lambda_grid,
                    kappa_grid,
                    bias,
                    gamma_grid,
                    psi_grid,
                    delta,
                    burn_in,
                    burn_in_dsc,
                    m,
                    equal_weight,
                    incl,
                    parallel,
                    n_threads,
                    portfolio_params)

    # Apply STSC-Parallel-Function
    results_par <- stsc(y_input,
                        X,
                        Ext_F,
                        init,
                        lambda_grid,
                        kappa_grid,
                        bias,
                        gamma_grid,
                        psi_grid,
                        delta,
                        burn_in,
                        burn_in_dsc,
                        m,
                        equal_weight,
                        incl,
                        TRUE,
                        n_threads,
                        portfolio_params)

    # Forecasts
    expect_equal(results$Forecasts$Point_Forecasts,
                 results_par$Forecasts$Point_Forecasts,
                 info = paste("Mismatch in Point_Forecasts for m =", m))

    expect_equal(results$Forecasts$Variance_Forecasts,
                 results_par$Forecasts$Variance_Forecasts,
                 info = paste("Mismatch in Variance_Forecasts for m =", m))

    expect_equal(results$Forecasts$Realization,
                 results_par$Forecasts$Realization,
                 info = paste("Mismatch in Realization for m =", m))

    # Tuning Parameters
    expect_equal(results$Tuning_Parameters$Gamma,
                 results_par$Tuning_Parameters$Gamma,
                 info = paste("Mismatch in Gamma for m =", m))

    expect_equal(results$Tuning_Parameters$Psi,
                 results_par$Tuning_Parameters$Psi,
                 info = paste("Mismatch in Psi for m =", m))

    expect_equal(results$Tuning_Parameters$Signals,
                 results_par$Tuning_Parameters$Signals,
                 info = paste("Mismatch in Signals for m =", m))

    expect_equal(results$Tuning_Parameters$Lambda,
                 results_par$Tuning_Parameters$Lambda,
                 info = paste("Mismatch in Lambda for m =", m))

    expect_equal(results$Tuning_Parameters$Kappa,
                 results_par$Tuning_Parameters$Kappa,
                 info = paste("Mismatch in Kappa for m =", m))
  }
})

###########################################################
### Test same results between STSC and TVC/DSC
test_that("Test whether the STSC-Function and TVC/DSC-Function return the same results", {

  # Apply STSC-Function
  results <- stsc(y,
                  X,
                  Ext_F,
                  init,
                  lambda_grid,
                  kappa_grid,
                  bias,
                  gamma_grid,
                  psi_grid,
                  delta,
                  burn_in,
                  burn_in_dsc,
                  metric,
                  equal_weight,
                  incl,
                  parallel,
                  n_threads,
                  portfolio_params)

  # Apply STSC-Function
  results_par <- stsc(y,
                      X,
                      Ext_F,
                      init,
                      lambda_grid,
                      kappa_grid,
                      bias,
                      gamma_grid,
                      psi_grid,
                      delta,
                      burn_in,
                      burn_in_dsc,
                      metric,
                      equal_weight,
                      incl,
                      TRUE,
                      n_threads,
                      portfolio_params)

  # Apply TVC-Function
  tvc_results <- hdflex::tvc(y,
                             X,
                             Ext_F,
                             init,
                             lambda_grid,
                             kappa_grid,
                             bias)

  # Assign TVC-Results
  forecast_tvc <- tvc_results$Forecasts$Point_Forecasts
  variance_tvc <- tvc_results$Forecasts$Variance_Forecasts

  # First complete forecast period (no missing values)
  sub_period <- seq(which(complete.cases(forecast_tvc))[1], nrow(y))

  ### Part 2: DSC-Function
  # Apply DSC-Function
  dsc_results <- hdflex::dsc(y[sub_period],
                             forecast_tvc[sub_period, , drop = FALSE],
                             variance_tvc[sub_period, , drop = FALSE],
                             gamma_grid,
                             psi_grid,
                             delta,
                             burn_in,
                             burn_in_dsc,
                             metric,
                             equal_weight,
                             incl,
                             NULL)

  # Forecasts
  expect_true(
    all(
      all.equal(
        na.omit(results$Forecasts$Point_Forecasts),
        na.omit(results_par$Forecasts$Point_Forecasts),
        check.attributes = FALSE,
      ) == TRUE,
      all.equal(
        na.omit(results$Forecasts$Point_Forecasts),
        na.omit(dsc_results$Forecasts$Point_Forecasts),
        check.attributes = FALSE,
      ) == TRUE
    )
  )

  expect_true(
    all(
      all.equal(
        na.omit(results$Forecasts$Variance_Forecasts),
        na.omit(results_par$Forecasts$Variance_Forecasts),
        check.attributes = FALSE,
      ) == TRUE,
      all.equal(
        na.omit(results$Forecasts$Variance_Forecasts),
        na.omit(dsc_results$Forecasts$Variance_Forecasts),
        check.attributes = FALSE,
      ) == TRUE
    )
  )

  # Tuning Parameters
  expect_true(
    all(
      all.equal(
        na.omit(results$Tuning_Parameters$Gamma),
        na.omit(dsc_results$Tuning_Parameters$Gamma),
        check.attributes = FALSE
      ) == TRUE,
      all.equal(
        na.omit(results$Tuning_Parameters$Gamma),
        na.omit(results_par$Tuning_Parameters$Gamma),
        check.attributes = FALSE
      ) == TRUE
    )
  )

  expect_true(
    all(
      all.equal(
        na.omit(results$Tuning_Parameters$Psi),
        na.omit(dsc_results$Tuning_Parameters$Psi),
        check.attributes = FALSE
      ) == TRUE,
      all.equal(
        na.omit(results$Tuning_Parameters$Psi),
        na.omit(results_par$Tuning_Parameters$Psi),
        check.attributes = FALSE
      ) == TRUE
    )
  )
})

###########################################################
### Test STSC with missing values
test_that("Test whether the STSC-Function works with missing values", {

  # Set Missing Values
  X[1:20, 1] <- NA
  Ext_F[1:15, 1] <- NA

  # Apply STSC-Function
  results <- stsc(y,
                  X,
                  Ext_F,
                  init,
                  lambda_grid,
                  kappa_grid,
                  bias,
                  gamma_grid,
                  psi_grid,
                  delta,
                  burn_in,
                  burn_in_dsc,
                  metric,
                  equal_weight,
                  incl,
                  parallel,
                  n_threads,
                  portfolio_params)

  # Apply STSC-Function
  results_par <- stsc(y,
                      X,
                      Ext_F,
                      init,
                      lambda_grid,
                      kappa_grid,
                      bias,
                      gamma_grid,
                      psi_grid,
                      delta,
                      burn_in,
                      burn_in_dsc,
                      metric,
                      equal_weight,
                      incl,
                      TRUE,
                      NULL,
                      portfolio_params)

  # Compare Forecasts
  expect_equal(results$Forecasts$Point_Forecasts,
                         results_par$Forecasts$Point_Forecasts)

  # List Contains three Elements
  expect_equal(length(results), 3)

  # Forecasts List Contains three Elements
  expect_equal(length(results$Forecasts), 3)

  # Point Forecasts
  expect_numeric(results$Forecasts$Point_Forecasts,
                            len = n_obs,
                            finite = TRUE)

  # Variance Forecasts
  expect_numeric(results$Forecasts$Variance_Forecasts,
                            len = n_obs,
                            lower = 0,
                            finite = TRUE)
})

###########################################################
### Test STSC with inclusion
test_that("Test whether the STSC-Function works with inclusion", {

  # Set Inclusion
  incl <- c(1, 2)
  psi_grid <- c(8:20)

  # Apply STSC-Function
  results <- stsc(y,
                  X,
                  Ext_F,
                  init,
                  lambda_grid,
                  kappa_grid,
                  bias,
                  gamma_grid,
                  psi_grid,
                  delta,
                  burn_in,
                  burn_in_dsc,
                  metric,
                  equal_weight,
                  incl,
                  parallel,
                  n_threads,
                  portfolio_params)

  # Apply STSC-Function
  results_par <- stsc(y,
                      X,
                      Ext_F,
                      init,
                      lambda_grid,
                      kappa_grid,
                      bias,
                      gamma_grid,
                      psi_grid,
                      delta,
                      burn_in,
                      burn_in_dsc,
                      metric,
                      equal_weight,
                      incl,
                      TRUE,
                      n_threads,
                      portfolio_params)

  # Compare Forecasts
  expect_equal(results$Forecasts$Point_Forecasts,
                         results_par$Forecasts$Point_Forecasts)

  # Cut-Off
  cut_off <- seq(max(burn_in, burn_in_dsc))

  # List Contains three Elements
  expect_equal(length(results), 3)

  # Forecasts List Contains three Elements
  expect_equal(length(results$Forecasts), 3)

  # Point Forecasts
  expect_numeric(results$Forecasts$Point_Forecasts,
                            len = n_obs,
                            finite = TRUE)

  # Variance Forecasts
  expect_numeric(results$Forecasts$Variance_Forecasts,
                            len = n_obs,
                            lower = 0,
                            finite = TRUE)

  # Tuning Parameters List Contains five Elements
  expect_equal(length(results$Tuning_Parameters), 5)

  # Psi-Vector
  expect_numeric(results$Tuning_Parameters$Psi,
                            len = n_obs,
                            lower = min(psi_grid),
                            upper = max(psi_grid),
                            finite = TRUE)

  # Signals
  expect_matrix(results$Tuning_Parameters$Signals,
                           mode = "integerish",
                           nrows = n_obs,
                           ncols = (ncol(X) + ncol(Ext_F)))

  # Check if the Signals in incl were really selected
  for (i in incl) {
    expect_true(all(results$Tuning_Parameters$Signals[-cut_off, i] > 0),
                info = paste("Column", i, "contains zeros"))
  }

  # Lambda-Vector
  expect_matrix(results$Tuning_Parameters$Lambda,
                           mode = "integerish",
                           nrows = n_obs,
                           ncols = length(lambda_grid))

  # Check that the Lambda matrix does not contain any zeros
  expect_true(all(results$Tuning_Parameters$Lambda[-cut_off, ] > 0),
              info = "Lambda matrix contains zeros")

  # Kappa-Vector
  expect_matrix(results$Tuning_Parameters$Kappa,
                           mode = "integerish",
                           nrows = n_obs,
                           ncols = length(kappa_grid))

  # Check that the Kappa matrix does not contain any zeros
  expect_true(all(results$Tuning_Parameters$Kappa[-cut_off, ] > 0),
              info = "Kappa matrix contains zeros")
})

###########################################################
### Test STSC with X
test_that("Test whether the STSC-Function works with X", {

  # Apply STSC-Function
  results <- stsc(y,
                  X,
                  NULL,
                  init,
                  lambda_grid,
                  kappa_grid,
                  bias,
                  gamma_grid,
                  psi_grid,
                  delta,
                  burn_in,
                  burn_in_dsc,
                  metric,
                  equal_weight,
                  incl,
                  parallel,
                  n_threads,
                  portfolio_params)

  # Apply STSC-Function
  results_par <- stsc(y,
                      X,
                      NULL,
                      init,
                      lambda_grid,
                      kappa_grid,
                      bias,
                      gamma_grid,
                      psi_grid,
                      delta,
                      burn_in,
                      burn_in_dsc,
                      metric,
                      equal_weight,
                      incl,
                      TRUE,
                      n_threads,
                      portfolio_params)

  # Compare Forecasts
  expect_equal(results$Forecasts$Point_Forecasts,
                         results_par$Forecasts$Point_Forecasts)

  # List Contains three Elements
  expect_equal(length(results), 3)

  # Forecasts List Contains three Elements
  expect_equal(length(results$Forecasts), 3)

  # Point Forecasts
  expect_numeric(results$Forecasts$Point_Forecasts,
                            len = n_obs,
                            finite = TRUE)

  # Variance Forecasts
  expect_numeric(results$Forecasts$Variance_Forecasts,
                            len = n_obs,
                            lower = 0,
                            finite = TRUE)

  # Tuning Parameters List Contains five Elements
  expect_equal(length(results$Tuning_Parameters), 5)

  # Gamma-Vector
  expect_numeric(results$Tuning_Parameters$Gamma,
                            len = n_obs,
                            lower = min(gamma_grid),
                            upper = max(gamma_grid),
                            finite = TRUE)

  # Psi-Vector
  expect_numeric(results$Tuning_Parameters$Psi,
                            len = n_obs,
                            lower = min(psi_grid),
                            upper = max(psi_grid),
                            finite = TRUE)

  # Signals
  expect_matrix(results$Tuning_Parameters$Signals,
                           mode = "integerish",
                           nrows = n_obs,
                           ncols = ncol(X))

  # Lambda-Vector
  expect_matrix(results$Tuning_Parameters$Lambda,
                           nrows = n_obs,
                           mode = "integerish",
                           ncols = length(lambda_grid))

  # Kappa-Vector
  expect_matrix(results$Tuning_Parameters$Kappa,
                           nrows = n_obs,
                           mode = "integerish",
                           ncols = length(kappa_grid))
})

### Test STSC with Ext_F
test_that("Test whether the STSC-Function works with Ext_F", {

  # Apply STSC-Function
  results <- stsc(y,
                  NULL,
                  Ext_F,
                  init,
                  lambda_grid,
                  kappa_grid,
                  bias,
                  gamma_grid,
                  psi_grid,
                  delta,
                  burn_in,
                  burn_in_dsc,
                  metric,
                  equal_weight,
                  incl,
                  parallel,
                  n_threads,
                  portfolio_params)

  # Apply STSC-Function
  results_par <- stsc(y,
                      NULL,
                      Ext_F,
                      init,
                      lambda_grid,
                      kappa_grid,
                      bias,
                      gamma_grid,
                      psi_grid,
                      delta,
                      burn_in,
                      burn_in_dsc,
                      metric,
                      equal_weight,
                      incl,
                      TRUE,
                      n_threads,
                      portfolio_params)

  # Compare Forecasts
  expect_equal(results$Forecasts$Point_Forecasts,
                         results_par$Forecasts$Point_Forecasts)

  # List Contains three Elements
  expect_equal(length(results), 3)

  # Forecasts List Contains three Elements
  expect_equal(length(results$Forecasts), 3)

  # Point Forecasts
  expect_numeric(results$Forecasts$Point_Forecasts,
                            len = n_obs,
                            finite = TRUE)

  # Variance Forecasts
  expect_numeric(results$Forecasts$Variance_Forecasts,
                            len = n_obs,
                            lower = 0,
                            finite = TRUE)

  # Tuning Parameters List Contains five Elements
  expect_equal(length(results$Tuning_Parameters), 5)

  # Gamma-Vector
  expect_numeric(results$Tuning_Parameters$Gamma,
                            len = n_obs,
                            lower = min(gamma_grid),
                            upper = max(gamma_grid),
                            finite = TRUE)

  # Psi-Vector
  expect_numeric(results$Tuning_Parameters$Psi,
                            len = n_obs,
                            lower = min(psi_grid),
                            upper = max(psi_grid),
                            finite = TRUE)

  # Signals
  expect_matrix(results$Tuning_Parameters$Signals,
                           mode = "integerish",
                           nrows = n_obs,
                           ncols = ncol(Ext_F))

  # Lambda-Vector
  expect_matrix(results$Tuning_Parameters$Lambda,
                           nrows = n_obs,
                           ncols = length(lambda_grid))

  # Kappa-Vector
  expect_matrix(results$Tuning_Parameters$Kappa,
                           nrows = n_obs,
                           ncols = length(kappa_grid))
})

###########################################################
### Test STSC without Bias Correction
test_that("Test whether the STSC-Function works without Bias", {

  # Set Signal constant for init periods
  X[1:10, 1] <- 0
  Ext_F[1:10, 1] <- 0

  # Apply STSC-Function
  results <- stsc(y,
                  X,
                  Ext_F,
                  init,
                  lambda_grid,
                  kappa_grid,
                  FALSE,
                  gamma_grid,
                  psi_grid,
                  delta,
                  burn_in,
                  burn_in_dsc,
                  metric,
                  equal_weight,
                  incl,
                  parallel,
                  n_threads,
                  portfolio_params)

  # Apply STSC-Function
  results_par <- stsc(y,
                      X,
                      Ext_F,
                      init,
                      lambda_grid,
                      kappa_grid,
                      FALSE,
                      gamma_grid,
                      psi_grid,
                      delta,
                      burn_in,
                      burn_in_dsc,
                      metric,
                      equal_weight,
                      incl,
                      TRUE,
                      n_threads,
                      portfolio_params)

  # Compare Forecasts
  expect_equal(results$Forecasts$Point_Forecasts,
                         results_par$Forecasts$Point_Forecasts)

  # List Contains three Elements
  expect_equal(length(results), 3)

  # Forecasts List Contains three Elements
  expect_equal(length(results$Forecasts), 3)

  # Point Forecasts
  expect_numeric(results$Forecasts$Point_Forecasts,
                            len = n_obs,
                            finite = TRUE)

  # Variance Forecasts
  expect_numeric(results$Forecasts$Variance_Forecasts,
                            len = n_obs,
                            lower = 0,
                            finite = TRUE)

  # Tuning Parameters List Contains five Elements
  expect_equal(length(results$Tuning_Parameters), 5)

  # Gamma-Vector
  expect_numeric(results$Tuning_Parameters$Gamma,
                            len = n_obs,
                            lower = min(gamma_grid),
                            upper = max(gamma_grid),
                            finite = TRUE)

  # Psi-Vector
  expect_numeric(results$Tuning_Parameters$Psi,
                            len = n_obs,
                            lower = min(psi_grid),
                            upper = max(psi_grid),
                            finite = TRUE)

  # Signals
  expect_matrix(results$Tuning_Parameters$Signals,
                           mode = "integerish",
                           nrows = n_obs,
                           ncols = (ncol(X) + ncol(Ext_F)))

  # Lambda-Vector
  expect_matrix(results$Tuning_Parameters$Lambda,
                           mode = "integerish",
                           nrows = n_obs,
                           ncols = length(lambda_grid))

  # Kappa-Vector
  expect_matrix(results$Tuning_Parameters$Kappa,
                           mode = "integerish",
                           nrows = n_obs,
                           ncols = length(kappa_grid))
})

###########################################################
### Test STSC with equal weight option
test_that("Test whether the STSC-Function works with equal weight option", {

  # Apply STSC-Function
  results <- stsc(y,
                  X,
                  Ext_F,
                  init,
                  lambda_grid,
                  kappa_grid,
                  bias,
                  gamma_grid,
                  psi_grid,
                  delta,
                  burn_in,
                  burn_in_dsc,
                  metric,
                  FALSE,
                  incl,
                  parallel,
                  n_threads,
                  portfolio_params)

  # Apply STSC-Function
  results_par <- stsc(y,
                      X,
                      Ext_F,
                      init,
                      lambda_grid,
                      kappa_grid,
                      bias,
                      gamma_grid,
                      psi_grid,
                      delta,
                      burn_in,
                      burn_in_dsc,
                      metric,
                      FALSE,
                      incl,
                      TRUE,
                      n_threads,
                      portfolio_params)

  # Compare Forecasts
  expect_equal(results$Forecasts$Point_Forecasts,
                         results_par$Forecasts$Point_Forecasts)

  # List Contains three Elements
  expect_equal(length(results), 3)

  # Forecasts List Contains three Elements
  expect_equal(length(results$Forecasts), 3)

  # Point Forecasts
  expect_numeric(results$Forecasts$Point_Forecasts,
                            len = n_obs,
                            finite = TRUE)

  # Variance Forecasts
  expect_numeric(results$Forecasts$Variance_Forecasts,
                            len = n_obs,
                            lower = 0,
                            finite = TRUE)

  # Tuning Parameters List Contains five Elements
  expect_equal(length(results$Tuning_Parameters), 5)

  # Gamma-Vector
  expect_numeric(results$Tuning_Parameters$Gamma,
                            len = n_obs,
                            lower = min(gamma_grid),
                            upper = max(gamma_grid),
                            finite = TRUE)

  # Psi-Vector
  expect_numeric(results$Tuning_Parameters$Psi,
                            len = n_obs,
                            lower = min(psi_grid),
                            upper = max(psi_grid),
                            finite = TRUE)

  # Signals
  expect_matrix(results$Tuning_Parameters$Signals,
                           mode = "integerish",
                           nrows = n_obs,
                           ncols = (ncol(X) + ncol(Ext_F)))

  # Lambda-Vector
  expect_matrix(results$Tuning_Parameters$Lambda,
                           mode = "integerish",
                           nrows = n_obs,
                           ncols = length(lambda_grid))

  # Kappa-Vector
  expect_matrix(results$Tuning_Parameters$Kappa,
                           mode = "integerish",
                           nrows = n_obs,
                           ncols = length(kappa_grid))
})
