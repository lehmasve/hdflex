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

###########################################################
### Test TVC
test_that("Test TVC", {

  # Apply TVC-Function
  results <- hdflex::tvc(y,
                         X,
                         Ext_F,
                         init,
                         lambda_grid,
                         kappa_grid,
                         bias)

  # List Contains three Elements
  expect_equal(length(results), 2)

  # Forecasts List Contains three Elements
  expect_equal(length(results$Forecasts), 3)

  # Point Forecasts
  expect_matrix(results$Forecasts$Point_Forecasts,
                           nrow = n_obs,
                           ncol = (ncol(X) + ncol(Ext_F)) * length(lambda_grid) * length(kappa_grid))

  # Variance Forecasts
  expect_matrix(results$Forecasts$Variance_Forecasts,
                           nrow = n_obs,
                           ncol = (ncol(X) + ncol(Ext_F)) * length(lambda_grid) * length(kappa_grid))

  # Ensure minimum value is 0
  expect_true(all(results$Forecasts$Variance_Forecasts[-1, ] >= 0))

  # Realization
  expect_numeric(results$Forecasts$Realization,
                            len = n_obs,
                            finite = TRUE)

  # Model List Contains 4 Elements
  expect_equal(length(results$Model), 4)

  # Lambda Grid
  expect_numeric(results$Model$Lambda_grid,
                            len = length(lambda_grid),
                            finite = TRUE)

  # Kappa Grid
  expect_numeric(results$Model$Kappa_grid,
                            len = length(kappa_grid),
                            finite = TRUE)

  # Init
  expect_numeric(results$Model$Init,
                            len = 1,
                            finite = TRUE)

  # Bias
  expect_logical(results$Model$Bias,
                            len = 1,
                            any.missing = FALSE)
})

###########################################################
### Test STSC with X
test_that("Test whether the STSC-Function works with X", {

  # Apply TVC-Function
  results <- hdflex::tvc(y,
                         X,
                         NULL,
                         init,
                         lambda_grid,
                         kappa_grid,
                         bias)

  # List Contains three Elements
  expect_equal(length(results), 2)

  # Forecasts List Contains three Elements
  expect_equal(length(results$Forecasts), 3)

  # Point Forecasts
  expect_matrix(results$Forecasts$Point_Forecasts,
                           nrow = n_obs,
                           ncol = ncol(X) * length(lambda_grid) * length(kappa_grid))

  # Variance Forecasts
  expect_matrix(results$Forecasts$Variance_Forecasts,
                           nrow = n_obs,
                           ncol = ncol(X) * length(lambda_grid) * length(kappa_grid))

  # Ensure minimum value is 0
  expect_true(all(results$Forecasts$Variance_Forecasts[-1, ] >= 0))
})

###########################################################
### Test STSC with Ext_F
test_that("Test whether the STSC-Function works with Ext_F", {

  # Apply TVC-Function
  results <- hdflex::tvc(y,
                         NULL,
                         Ext_F,
                         init,
                         lambda_grid,
                         kappa_grid,
                         bias)

  # List Contains three Elements
  expect_equal(length(results), 2)

  # Forecasts List Contains three Elements
  expect_equal(length(results$Forecasts), 3)

  # Point Forecasts
  expect_matrix(results$Forecasts$Point_Forecasts,
                           nrow = n_obs,
                           ncol =  ncol(Ext_F) * length(lambda_grid) * length(kappa_grid))

  # Variance Forecasts
  expect_matrix(results$Forecasts$Variance_Forecasts,
                           nrow = n_obs,
                           ncol = ncol(Ext_F) * length(lambda_grid) * length(kappa_grid))

  # Ensure minimum value is 0
  expect_true(all(results$Forecasts$Variance_Forecasts[-1, ] >= 0))
})

###########################################################
### Test TVC without Bias Correction
test_that("Test TVC without Bias Correction", {

  # Apply TVC-Function
  results <- hdflex::tvc(y,
                         X,
                         Ext_F,
                         init,
                         lambda_grid,
                         kappa_grid,
                         FALSE)

  # List Contains three Elements
  expect_equal(length(results), 2)

  # Forecasts List Contains three Elements
  expect_equal(length(results$Forecasts), 3)

  # Point Forecasts
  expect_matrix(results$Forecasts$Point_Forecasts,
                           nrow = n_obs,
                           ncol = (ncol(X) + ncol(Ext_F)) * length(lambda_grid) * length(kappa_grid))

  # Variance Forecasts
  expect_matrix(results$Forecasts$Variance_Forecasts,
                           nrow = n_obs,
                           ncol = (ncol(X) + ncol(Ext_F)) * length(lambda_grid) * length(kappa_grid))

  # Ensure minimum value is 0
  expect_true(all(results$Forecasts$Variance_Forecasts[-1, ] >= 0))
})

###########################################################
### Test TVC with Missing Values
test_that("Test TVC with Missing / Constant Values", {

  # Add Missing Values
  X[1, 1:10] <- NA
  Ext_F[1, 1:10] <- NA

  X[2, 1:10] <- 1
  Ext_F[2, 1:10] <- 1

  # Apply TVC-Function
  results <- hdflex::tvc(y,
                         X,
                         Ext_F,
                         init,
                         lambda_grid,
                         kappa_grid,
                         bias)

  # List Contains three Elements
  expect_equal(length(results), 2)

  # Forecasts List Contains three Elements
  expect_equal(length(results$Forecasts), 3)

  # Point Forecasts
  expect_matrix(results$Forecasts$Point_Forecasts,
                           nrow = n_obs,
                           ncol = (ncol(X) + ncol(Ext_F)) * length(lambda_grid) * length(kappa_grid))

  # Variance Forecasts
  expect_matrix(results$Forecasts$Variance_Forecasts,
                           nrow = n_obs,
                           ncol = (ncol(X) + ncol(Ext_F)) * length(lambda_grid) * length(kappa_grid))

  # Ensure minimum value is 0
  expect_true(all(results$Forecasts$Variance_Forecasts[-c(1:11), ] >= 0))
})

###########################################################
### Test TVC without providing colnames
test_that("Test TVC without providing colnames", {

  # Remove Colnames
  colnames(X) <- NULL
  colnames(y) <- NULL

  # Apply TVC-Function
  results <- hdflex::tvc(y,
                         X,
                         Ext_F,
                         init,
                         lambda_grid,
                         kappa_grid,
                         bias)

  # List Contains three Elements
  expect_equal(length(results), 2)

  # Forecasts List Contains three Elements
  expect_equal(length(results$Forecasts), 3)

  # Point Forecasts
  expect_matrix(results$Forecasts$Point_Forecasts,
                           nrow = n_obs,
                           ncol = (ncol(X) + ncol(Ext_F)) * length(lambda_grid) * length(kappa_grid))

  # Variance Forecasts
  expect_matrix(results$Forecasts$Variance_Forecasts,
                           nrow = n_obs,
                           ncol = (ncol(X) + ncol(Ext_F)) * length(lambda_grid) * length(kappa_grid))

  # Ensure minimum value is 0
  expect_true(all(results$Forecasts$Variance_Forecasts[-1, ] >= 0))
})