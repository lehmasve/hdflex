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
portfolio_params <- NULL

###########################################################
### Create Density Forecast
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

###########################################################
### Create STSC and DSC-Objects
# Apply DSC-Function
dsc_results <- dsc(y[-1, , drop = FALSE],
                   forecast_tvc[-1, ],
                   variance_tvc[-1, ],
                   gamma_grid,
                   psi_grid,
                   delta,
                   burn_in,
                   burn_in_dsc,
                   metric,
                   equal_weight,
                   incl,
                   portfolio_params)

# Apply STSC-Function
stsc_results <- stsc(y,
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

###########################################################
### STSC - Object
# Test Metrics
test_that("summary calculates metrics correctly", {
  result <- summary(stsc_results, eval_period = 50:500)

  eval_length <- length(50:500)

  expect_true(is.list(result))

  expect_true("MSE" %in% names(result))
  expect_numeric(result$MSE[[1]], lower = 0, len = 1, any.missing = FALSE)
  expect_numeric(result$MSE[[2]], lower = 0, len = eval_length, any.missing = FALSE)

  expect_true("ACRPS" %in% names(result))
  expect_numeric(result$ACRPS[[1]], lower = 0, len = 1, any.missing = FALSE)
  expect_numeric(result$ACRPS[[2]], lower = 0, len = eval_length, any.missing = FALSE)

  expect_true("APLL" %in% names(result))
  expect_numeric(result$APLL[[1]], len = 1, any.missing = FALSE)
  expect_numeric(result$APLL[[2]], len = eval_length, any.missing = FALSE)
})

# Test Plots
test_that("summary generates plots", {
  result <- summary(stsc_results)
  expect_true("Plots" %in% names(result))
  expect_true(is.ggplot(result$Plots$Gamma))
  expect_true(is.ggplot(result$Plots$Psi))
  expect_true(is.ggplot(result$Plots$Signals))
  expect_true(is.ggplot(result$Plots$Lambda))
  expect_true(is.ggplot(result$Plots$Kappa))
})

### DSC - Object
# Test Metrics
test_that("summary calculates metrics correctly", {
  result <- summary(dsc_results, eval_period = 50:499)

  eval_length <- length(50:499)

  expect_true(is.list(result))

  expect_true("MSE" %in% names(result))
  expect_numeric(result$MSE[[1]], lower = 0, len = 1, any.missing = FALSE)
  expect_numeric(result$MSE[[2]], lower = 0, len = eval_length, any.missing = FALSE)

  expect_true("ACRPS" %in% names(result))
  expect_numeric(result$ACRPS[[1]], lower = 0, len = 1, any.missing = FALSE)
  expect_numeric(result$ACRPS[[2]], lower = 0, len = eval_length, any.missing = FALSE)

  expect_true("APLL" %in% names(result))
  expect_numeric(result$APLL[[1]], len = 1, any.missing = FALSE)
  expect_numeric(result$APLL[[2]], len = eval_length, any.missing = FALSE)
})

# Test Plots
test_that("summary generates plots", {
  result <- summary(dsc_results)
  expect_true("Plots" %in% names(result))
  expect_true(is.ggplot(result$Plots$Gamma))
  expect_true(is.ggplot(result$Plots$Psi))
  expect_true(is.ggplot(result$Plots$CFM))
})