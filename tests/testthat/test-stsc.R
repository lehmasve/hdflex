###########################################################
### Simulate Data
# Set Seed
set.seed(123)

# Set Dimensions
numb_obs      <-  500
numb_noise    <-  49
numb_signals  <-  numb_noise + 1
numb_forc     <-  50

# Set up Coefficient-Matrix
theta  <-  matrix(NA, nrow = numb_obs, ncol = numb_signals)

# Loop over Time
for (t in seq_len(numb_obs)) {

   ### Theta: Abrupt Change
  theta[t, 1]  <-  ifelse((t > 200 & t < 450), -0.30, 0.30)

   ### Noise Predictors
  theta[t, -1] <-  0
}

### Draw Simple Signals
raw_signals            <-  replicate(numb_signals, rnorm(numb_obs, 0, 0.5))
colnames(raw_signals)  <-  paste0("X", as.character(seq_len(numb_signals)))

### Draw Noise
eps  <-  rnorm(numb_obs, 0, 0.1)

### Compute Target Variable
y    <-  rowSums(cbind(theta * raw_signals, eps))

# Create 'External' Forecasts
f_signals            <-  y + replicate(numb_forc, rnorm(numb_obs, 0, 0.5))
colnames(f_signals)  <-  paste0("F", as.character(seq_len(numb_forc)))
###########################################################

###########################################################
### STSC Parameter
# TV-C-Parameter
sample_length  <-  100
lambda_grid    <-  c(0.95, 1.00)
kappa_grid     <-  0.97

# Set DSC-Parameter
gamma_grid  <-  c(0.9, 0.95, 1)
psi_grid    <-  c(1:10)
delta       <-  0.95

# Set Method-Parameter
burn_in_tvc   <-  5
burn_in_dsc   <-  5
method        <-  1

# Set Portfolio-Parameter
risk_aversion <-  3
min_weight    <-  0
max_weight    <-  2

###########################################################

###########################################################
### Tests on Y
test_that("Test whether y is Numeric Vector", {

  y  <-  as.data.frame(y)
  testthat::expect_error(stsc(y,
                              raw_signals,
                              f_signals,
                              sample_length,
                              lambda_grid,
                              kappa_grid,
                              burn_in_tvc,
                              gamma_grid,
                              psi_grid,
                              delta,
                              burn_in_dsc,
                              method,
                              risk_aversion,
                              min_weight,
                              max_weight),
              "Must be of type 'numeric', not 'data.frame'.", fixed = TRUE)
})

test_that("Test whether y is not NULL", {

  y  <-  NULL
  testthat::expect_error(stsc(y,
                              raw_signals,
                              f_signals,
                              sample_length,
                              lambda_grid,
                              kappa_grid,
                              burn_in_tvc,
                              gamma_grid,
                              psi_grid,
                              delta,
                              burn_in_dsc,
                              method,
                              risk_aversion,
                              min_weight,
                              max_weight),
              "Must be of type 'numeric', not 'NULL'.", fixed = TRUE)
})

test_that("Test whether y has only numeric values", {

  y[10]  <-  "test"
  testthat::expect_error(stsc(y,
                              raw_signals,
                              f_signals,
                              sample_length,
                              lambda_grid,
                              kappa_grid,
                              burn_in_tvc,
                              gamma_grid,
                              psi_grid,
                              delta,
                              burn_in_dsc,
                              method,
                              risk_aversion,
                              min_weight,
                              max_weight),
              "Must be of type 'numeric', not 'character'.", fixed = TRUE)
})

test_that("Test whether y has no NA-Values", {

  y[10]  <-  NA
  testthat::expect_error(stsc(y,
                              raw_signals,
                              f_signals,
                              sample_length,
                              lambda_grid,
                              kappa_grid,
                              burn_in_tvc,
                              gamma_grid,
                              psi_grid,
                              delta,
                              burn_in_dsc,
                              method,
                              risk_aversion,
                              min_weight,
                              max_weight),
              "Contains missing values", fixed = TRUE)
})

### Tests on X
test_that("Test whether x is matrix", {

  raw_signals  <-  as.data.frame(raw_signals)
  testthat::expect_error(stsc(y,
                              raw_signals,
                              f_signals,
                              sample_length,
                              lambda_grid,
                              kappa_grid,
                              burn_in_tvc,
                              gamma_grid,
                              psi_grid,
                              delta,
                              burn_in_dsc,
                              method,
                              risk_aversion,
                              min_weight,
                              max_weight),
              "Must be of type 'matrix' (or 'NULL')", fixed = TRUE)
})

test_that("Test whether x has the same number of observations as y", {

  raw_signals  <-  raw_signals[1:10, ]
  testthat::expect_error(stsc(y,
                              raw_signals,
                              f_signals,
                              sample_length,
                              lambda_grid,
                              kappa_grid,
                              burn_in_tvc,
                              gamma_grid,
                              psi_grid,
                              delta,
                              burn_in_dsc,
                              method,
                              risk_aversion,
                              min_weight,
                              max_weight),
                        "Must have exactly", fixed = TRUE)
})

test_that("Test whether exception works when cov_mat cannot be initialised", {

  raw_signals[1:100, 10]  <-  0
  testthat::expect_no_error(stsc(y,
                                 raw_signals,
                                 f_signals,
                                 sample_length,
                                 lambda_grid,
                                 kappa_grid,
                                 burn_in_tvc,
                                 gamma_grid,
                                 psi_grid,
                                 delta,
                                 burn_in_dsc,
                                 method,
                                 risk_aversion,
                                 min_weight,
                                 max_weight))
})

### Tests on f
test_that("Test whether f is matrix", {

  f_signals  <-  as.data.frame(f_signals)
  testthat::expect_error(stsc(y,
                              raw_signals,
                              f_signals,
                              sample_length,
                              lambda_grid,
                              kappa_grid,
                              burn_in_tvc,
                              gamma_grid,
                              psi_grid,
                              delta,
                              burn_in_dsc,
                              method,
                              risk_aversion,
                              min_weight,
                              max_weight),
                        "Must be of type 'matrix' (or 'NULL')", fixed = TRUE)
})

test_that("Test whether f has the same number of observations as y", {

  f_signals  <-  f_signals[1:10, ]
  testthat::expect_error(stsc(y,
                              raw_signals,
                              f_signals,
                              sample_length,
                              lambda_grid,
                              kappa_grid,
                              burn_in_tvc,
                              gamma_grid,
                              psi_grid,
                              delta,
                              burn_in_dsc,
                              method,
                              risk_aversion,
                              min_weight,
                              max_weight),
                        "Must have exactly", fixed = TRUE)
})

### Tests on x and f
test_that("Test whether either x or f is provided", {

  raw_signals  <-  NULL
  f_signals    <-  NULL
  testthat::expect_error(stsc(y,
                              raw_signals,
                              f_signals,
                              sample_length,
                              lambda_grid,
                              kappa_grid,
                              burn_in_tvc,
                              gamma_grid,
                              psi_grid,
                              delta,
                              burn_in_dsc,
                              method,
                              risk_aversion,
                              min_weight,
                              max_weight),
            "Assertion failed. One of the following must apply:
 * checkmate::checkMatrix(X): Must be of type 'matrix', not 'NULL'
 * checkmate::checkMatrix(F): Must be of type 'matrix', not 'NULL'",
            fixed = TRUE)
})

test_that("Test whether Code still works with only raw signals / only point forecasts", {

  testthat::expect_no_error(stsc(y,
                                 raw_signals,
                                 NULL,
                                 sample_length,
                                 lambda_grid,
                                 kappa_grid,
                                 burn_in_tvc,
                                 gamma_grid,
                                 psi_grid,
                                 delta,
                                 burn_in_dsc,
                                 method,
                                 risk_aversion,
                                 min_weight,
                                 max_weight))

  testthat::expect_no_error(stsc(y,
                                 NULL,
                                 f_signals,
                                 sample_length,
                                 lambda_grid,
                                 kappa_grid,
                                 burn_in_tvc,
                                 gamma_grid,
                                 psi_grid,
                                 delta,
                                 burn_in_dsc,
                                 method,
                                 risk_aversion,
                                 min_weight,
                                 max_weight))
})

test_that("Test whether Code still works with NA-values", {

  raw_signals[1:20, c(1, 3, 5)] <-  NA
  testthat::expect_no_error(stsc(y,
                                 raw_signals,
                                 f_signals,
                                 sample_length,
                                 lambda_grid,
                                 kappa_grid,
                                 burn_in_tvc,
                                 gamma_grid,
                                 psi_grid,
                                 delta,
                                 burn_in_dsc,
                                 method,
                                 risk_aversion,
                                 min_weight,
                                 max_weight))
})

test_that("Test whether Code still works without dimnames", {

  colnames(raw_signals)  <-  NULL
  colnames(f_signals)    <-  NULL
  testthat::expect_no_error(stsc(y,
                                 raw_signals,
                                 f_signals,
                                 sample_length,
                                 lambda_grid,
                                 kappa_grid,
                                 burn_in_tvc,
                                 gamma_grid,
                                 psi_grid,
                                 delta,
                                 burn_in_dsc,
                                 method,
                                 risk_aversion,
                                 min_weight,
                                 max_weight))
})

### Tests on Methods
test_that("Test whether Code still works with different methods", {

  testthat::expect_no_error(stsc(y,
                                 raw_signals,
                                 f_signals,
                                 sample_length,
                                 lambda_grid,
                                 kappa_grid,
                                 burn_in_tvc,
                                 gamma_grid,
                                 psi_grid,
                                 delta,
                                 burn_in_dsc,
                                 method = 2,
                                 risk_aversion,
                                 min_weight,
                                 max_weight))

  testthat::expect_no_error(stsc(y,
                                 raw_signals,
                                 f_signals,
                                 sample_length,
                                 lambda_grid,
                                 kappa_grid,
                                 burn_in_tvc,
                                 gamma_grid,
                                 psi_grid,
                                 delta,
                                 burn_in_dsc,
                                 method = 3,
                                 risk_aversion,
                                 min_weight,
                                 max_weight))

  testthat::expect_no_error(stsc(y,
                                 raw_signals,
                                 f_signals,
                                 sample_length,
                                 lambda_grid,
                                 kappa_grid,
                                 burn_in_tvc,
                                 gamma_grid,
                                 psi_grid,
                                 delta,
                                 burn_in_dsc,
                                 method = 4,
                                 risk_aversion,
                                 min_weight,
                                 max_weight))

})

test_that("Test whether method is of given set", {

  method  <-  5
  testthat::expect_error(stsc(y,
                              raw_signals,
                              f_signals,
                              sample_length,
                              lambda_grid,
                              kappa_grid,
                              burn_in_tvc,
                              gamma_grid,
                              psi_grid,
                              delta,
                              burn_in_dsc,
                              method,
                              risk_aversion,
                              min_weight,
                              max_weight),
                        "Must be element of set {'1','2','3','4'}",
                        fixed = TRUE)
})

### Tests on Return Parameter
test_that("Test whether relevant parameter are provided", {

  method        <- 4
  risk_aversion <- NULL
  min_weight    <- NULL
  max_weight    <- NULL
  testthat::expect_error(stsc(y,
                              raw_signals,
                              f_signals,
                              sample_length,
                              lambda_grid,
                              kappa_grid,
                              burn_in_tvc,
                              gamma_grid,
                              psi_grid,
                              delta,
                              burn_in_dsc,
                              method,
                              risk_aversion,
                              min_weight,
                              max_weight),
        "Must be of type 'number', not 'NULL'",
            fixed = TRUE)
})


### Output
test_that("Test whether the output has the right format", {

  # Apply TVP-Function
  results  <-  stsc(y,
                    raw_signals,
                    f_signals,
                    sample_length,
                    lambda_grid,
                    kappa_grid,
                    burn_in_tvc,
                    gamma_grid,
                    psi_grid,
                    delta,
                    burn_in_dsc,
                    method,
                    risk_aversion,
                    min_weight,
                    max_weight)

  # List Contains Five Elements
  testthat::expect_equal(length(results), 5)

  # Number of Forecasts
  checkmate::expect_numeric(results[[1]],
                            len = numb_obs,
                            finite = TRUE)

  # Number of Variances
  checkmate::expect_numeric(results[[2]],
                            len = numb_obs,
                            lower = 0,
                            finite = TRUE)
  # Length of Gamma-Vector
  checkmate::expect_numeric(results[[3]],
                            len = numb_obs,
                            lower = 0,
                            finite = TRUE)
  # Length of Psi-Vector
  checkmate::expect_numeric(results[[4]],
                            len = numb_obs,
                            lower = 0,
                            finite = TRUE)

  # Dimension of selected Candidate Forecasts
  checkmate::expect_matrix(results[[5]],
                           mode = "integerish",
                           nrows = numb_obs,
                           ncols = (ncol(raw_signals) + ncol(f_signals)))
})