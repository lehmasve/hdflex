###########################################################
### Simulate Data
# Set Seed
set.seed(123)

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
n_threads <- NULL

# Set Portfolio-Parameter
portfolio_params <- c(3, 0, 2)

###########################################################
### Create Density Forecast
# Apply TVC-Function
tvc_results <- tvc(y,
                           X,
                           Ext_F,
                           init,
                           lambda_grid,
                           kappa_grid,
                           bias)

# Assign TVC-Results
forecast_tvc <- tvc_results$Forecasts$Point_Forecasts
variance_tvc <- tvc_results$Forecasts$Variance_Forecasts

# Remove NAs
y <- y[-1, , drop = FALSE]
returns <- returns[-1, , drop = FALSE]
forecast_tvc <- forecast_tvc[-1, ]
variance_tvc <- variance_tvc[-1, ]

###########################################################
### Test DSC (with test_that)
test_that("DSC-Function works correctly", {

  # Apply DSC-Function
  results <- dsc(y,
                 forecast_tvc,
                 variance_tvc,
                 gamma_grid,
                 psi_grid,
                 delta,
                 burn_in,
                 burn_in_dsc,
                 1,
                 equal_weight,
                 incl,
                 portfolio_params)

  # Apply DSC-Function
  results <- dsc(y,
                 forecast_tvc,
                 variance_tvc,
                 gamma_grid,
                 psi_grid,
                 delta,
                 burn_in,
                 burn_in_dsc,
                 2,
                 equal_weight,
                 incl,
                 portfolio_params)

  # Apply DSC-Function
  results <- dsc(y,
                 forecast_tvc,
                 variance_tvc,
                 gamma_grid,
                 psi_grid,
                 delta,
                 burn_in,
                 burn_in_dsc,
                 3,
                 equal_weight,
                 incl,
                 portfolio_params)

  # Apply DSC-Function
  results <- dsc(returns,
                 forecast_tvc,
                 variance_tvc,
                 gamma_grid,
                 psi_grid,
                 delta,
                 burn_in,
                 burn_in_dsc,
                 4,
                 equal_weight,
                 incl,
                 portfolio_params)

  # Apply DSC-Function
  results <- dsc(y,
                 forecast_tvc,
                 variance_tvc,
                 gamma_grid,
                 psi_grid,
                 delta,
                 burn_in,
                 burn_in_dsc,
                 metric,
                 equal_weight,
                 incl,
                 portfolio_params)

  # List Contains three Elements
  expect_equal(length(results), 3)

  # Forecasts List Contains three Elements
  expect_equal(length(results$Forecasts), 3)

  # Point Forecasts
  expect_numeric(results$Forecasts$Point_Forecasts,
                            len = n_obs - 1,
                            finite = TRUE)

  # Variance Forecasts
  expect_numeric(results$Forecasts$Variance_Forecasts,
                            len = n_obs - 1,
                            lower = 0,
                            finite = TRUE)

  # Realization
  expect_numeric(results$Forecasts$Realization,
                            len = n_obs - 1,
                            finite = TRUE)

  # Tuning Parameters List Contains three Elements
  expect_equal(length(results$Tuning_Parameters), 3)

  # Gamma-Vector
  expect_numeric(results$Tuning_Parameters$Gamma,
                            len = n_obs - 1,
                            lower = min(gamma_grid),
                            upper = max(gamma_grid),
                            finite = TRUE)

  # Psi-Vector
  expect_numeric(results$Tuning_Parameters$Psi,
                            len = n_obs - 1,
                            lower = min(psi_grid),
                            upper = max(psi_grid),
                            finite = TRUE)

  # CFM
  expect_matrix(results$Tuning_Parameters$CFM,
                           mode = "integerish",
                           nrows = n_obs - 1,
                           ncols = ncol(forecast_tvc))

  # Model List Contains 8 Elements
  expect_equal(length(results$Model), 8)

  # Gamma Grid
  expect_numeric(results$Model$Gamma_grid,
                            len = length(gamma_grid),
                            finite = TRUE)

  # Psi Grid
  expect_numeric(results$Model$Psi_grid,
                            len = length(psi_grid),
                            finite = TRUE)

  # Delta
  expect_numeric(results$Model$Delta,
                            len = 1,
                            finite = TRUE)

  # Burn-in
  expect_numeric(results$Model$Burn_in,
                            len = 1,
                            finite = TRUE)

  # Burn-in DSC
  expect_numeric(results$Model$Burn_in_dsc,
                            len = 1,
                            finite = TRUE)

  # Metric
  expect_numeric(results$Model$Metric,
                            len = 1,
                            finite = TRUE)

  # Equal Weight
  expect_equal(results$Model$Equal_weight, equal_weight)

  # Incl
  expect_equal(results$Model$Incl, incl)
})

###########################################################
### Test DSC with inclusion
###########################################################
### Test STSC with inclusion
test_that("Test whether the STSC-Function works with inclusion", {

  # Set Inclusion
  incl <- c(1, 2)
  psi_grid <- c(8:20)

  # Apply DSC-Function
  results <- dsc(y,
                         forecast_tvc,
                         variance_tvc,
                         gamma_grid,
                         psi_grid,
                         delta,
                         burn_in,
                         burn_in_dsc,
                         metric,
                         equal_weight,
                         incl,
                         portfolio_params)

  # Cut-Off
  cut_off <- seq(max(burn_in, burn_in_dsc))

  # List Contains three Elements
  expect_equal(length(results), 3)

  # Forecasts List Contains three Elements
  expect_equal(length(results$Forecasts), 3)

  # Point Forecasts
  expect_numeric(results$Forecasts$Point_Forecasts,
                            len = n_obs - 1,
                            finite = TRUE)

  # Variance Forecasts
  expect_numeric(results$Forecasts$Variance_Forecasts,
                            len = n_obs - 1,
                            lower = 0,
                            finite = TRUE)

  # Tuning Parameters List Contains three Elements
  expect_equal(length(results$Tuning_Parameters), 3)

  # Psi-Vector
  expect_numeric(results$Tuning_Parameters$Psi,
                            len = n_obs - 1,
                            lower = min(psi_grid),
                            upper = max(psi_grid),
                            finite = TRUE)

  # CFM
  expect_matrix(results$Tuning_Parameters$CFM,
                           mode = "integerish",
                           nrows = n_obs - 1,
                           ncols = ncol(forecast_tvc))

  # Check if the CFMs in incl were really selected
  for (i in incl) {
    expect_true(all(results$Tuning_Parameters$CFM[-cut_off, i] > 0),
                info = paste("Column", i, "contains zeros"))
  }
})

###########################################################
### Test STSC with equal weight option
test_that("Test whether the STSC-Function works with equal weight option", {

  # Apply DSC-Function
  results <- dsc(y,
                         forecast_tvc,
                         variance_tvc,
                         gamma_grid,
                         psi_grid,
                         delta,
                         burn_in,
                         burn_in_dsc,
                         metric,
                         FALSE,
                         incl,
                         portfolio_params)

  # List Contains three Elements
  expect_equal(length(results), 3)

  # Forecasts List Contains three Elements
  expect_equal(length(results$Forecasts), 3)

  # Point Forecasts
  expect_numeric(results$Forecasts$Point_Forecasts,
                            len = n_obs - 1,
                            finite = TRUE)

  # Variance Forecasts
  expect_numeric(results$Forecasts$Variance_Forecasts,
                            len = n_obs - 1,
                            lower = 0,
                            finite = TRUE)

  # Tuning Parameters List Contains three Elements
  expect_equal(length(results$Tuning_Parameters), 3)

  # Gamma-Vector
  expect_numeric(results$Tuning_Parameters$Gamma,
                            len = n_obs - 1,
                            lower = min(gamma_grid),
                            upper = max(gamma_grid),
                            finite = TRUE)

  # Psi-Vector
  expect_numeric(results$Tuning_Parameters$Psi,
                            len = n_obs - 1,
                            lower = min(psi_grid),
                            upper = max(psi_grid),
                            finite = TRUE)

  # CFM
  expect_matrix(results$Tuning_Parameters$CFM,
                           mode = "integerish",
                           nrows = n_obs - 1,
                           ncols = ncol(forecast_tvc))
})
