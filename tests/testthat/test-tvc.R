### Simulate Data
set.seed(123)

# Set Dimensions
numb_obs   <-  500
numb_pred  <-  50
numb_forc  <-  10

# Create Random Target-Variable
target_var  <-  rnorm(n = numb_obs, mean = 0, sd = 1)

# Create Random Simple Signals
raw_signals           <-  replicate(numb_pred, sample(0:10, numb_obs, rep = TRUE), )
raw_names             <-  paste0("X", as.character(seq_len(ncol(raw_signals))))
colnames(raw_signals) <-  raw_names

# Create Random (External) Point Forecasts
f_signals            <-  replicate(10, rnorm(n = numb_obs, mean = 0, sd = 0.5), )
f_names              <-  paste0("F", as.character(seq_len(ncol(f_signals))))
colnames(f_signals)  <-  f_names

# Specify TV-C-Parameter
sample_length  <-  100
lambda_grid    <-  c(0.99, 0.999, 1.000)
kappa_grid     <-  c(0.94)
n_cores        <-  1

### Tests on Y
test_that("Test whether y is Numeric Vector", {

  target_var  <-  as.data.frame(target_var)
  testthat::expect_error(tvc(target_var,
                             raw_signals,
                             f_signals,
                             lambda_grid,
                             kappa_grid,
                             sample_length,
                             n_cores),
              "Must be of type 'numeric', not 'data.frame'.", fixed = TRUE)
})

test_that("Test whether y is not NULL", {

  target_var  <-  NULL
  testthat::expect_error(tvc(target_var,
                             raw_signals,
                             f_signals,
                             lambda_grid,
                             kappa_grid,
                             sample_length,
                             n_cores),
              "Must be of type 'numeric', not 'NULL'.", fixed = TRUE)
})

test_that("Test whether y has only numeric values", {

  target_var[10]  <-  "test"
  testthat::expect_error(tvc(target_var,
                             raw_signals,
                             f_signals,
                             lambda_grid,
                             kappa_grid,
                             sample_length,
                             n_cores),
              "Must be of type 'numeric', not 'character'.", fixed = TRUE)
})

test_that("Test whether y has no NA-Values", {

  target_var[10]  <-  NA
  testthat::expect_error(tvc(target_var,
                             raw_signals,
                             f_signals,
                             lambda_grid,
                             kappa_grid,
                             sample_length,
                             n_cores),
              "Contains missing values", fixed = TRUE)
})

### Tests on X
test_that("Test whether x is matrix", {

  raw_signals  <-  as.data.frame(raw_signals)
  testthat::expect_error(tvc(target_var,
                             raw_signals,
                             f_signals,
                             lambda_grid,
                             kappa_grid,
                             sample_length,
                             n_cores),
              "Must be of type 'matrix' (or 'NULL')", fixed = TRUE)
})

test_that("Test whether x has the same number of observations as y", {

  raw_signals  <-  raw_signals[1:10, ]
  testthat::expect_error(tvc(target_var,
                             raw_signals,
                             f_signals,
                             lambda_grid,
                             kappa_grid,
                             sample_length,
                             n_cores),
                        "Must have exactly", fixed = TRUE)
})

test_that("Test whether exception works when cov_mat cannot be initialised", {

  raw_signals[1:100, 10]  <-  0
  testthat::expect_no_error(tvc(target_var,
                                raw_signals,
                                f_signals,
                                lambda_grid,
                                kappa_grid,
                                sample_length,
                                n_cores))
})

### Tests on f
test_that("Test whether f is matrix", {

  f_signals  <-  as.data.frame(f_signals)
  testthat::expect_error(tvc(target_var,
                             raw_signals,
                             f_signals,
                             lambda_grid,
                             kappa_grid,
                             sample_length,
                             n_cores),
                        "Must be of type 'matrix' (or 'NULL')", fixed = TRUE)
})

test_that("Test whether f has the same number of observations as y", {

  f_signals  <-  f_signals[1:10, ]
  testthat::expect_error(tvc(target_var,
                             raw_signals,
                             f_signals,
                             lambda_grid,
                             kappa_grid,
                             sample_length,
                             n_cores),
                        "Must have exactly", fixed = TRUE)
})

### Tests on x and f
test_that("Test whether either x or f is provided", {

  raw_signals  <-  NULL
  f_signals    <-  NULL
  testthat::expect_error(tvc(target_var,
                             raw_signals,
                             f_signals,
                             lambda_grid,
                             kappa_grid,
                             sample_length,
                             n_cores),
            "Assertion failed. One of the following must apply:
 * checkmate::checkMatrix(X): Must be of type 'matrix', not 'NULL'
 * checkmate::checkMatrix(Ext_F): Must be of type 'matrix', not 'NULL'",
            fixed = TRUE)
})

test_that("Test whether Code still works without dimnames", {

  colnames(raw_signals)  <-  NULL
  colnames(f_signals)    <-  NULL
  testthat::expect_no_error(tvc(target_var,
                                raw_signals,
                                f_signals,
                                lambda_grid,
                                kappa_grid,
                                sample_length,
                                n_cores))
})

### Output
test_that("Test whether the output has the right format", {

  # Apply TVP-Function
  results  <-  tvc(target_var,
                   raw_signals,
                   f_signals,
                   lambda_grid,
                   kappa_grid,
                   sample_length,
                   n_cores)

  # List Contains Two Elements
  testthat::expect_equal(length(results), 2)

  # Number of Models
  numb_mods  <-  length(lambda_grid) * length(kappa_grid) * numb_pred +
                 length(lambda_grid) * length(kappa_grid) * numb_forc

  # Dimension of Forecasts
  checkmate::expect_matrix(results[[1]],
                           mode = "numeric",
                           nrows = numb_obs,
                           ncols = numb_mods)

  # Dimension of Variances
  checkmate::expect_matrix(results[[2]],
                           mode = "numeric",
                           nrows = numb_obs,
                           ncols = numb_mods)

  # Only positive values in Var-Matrix
  checkmate::qassert(results[[2]],
                     c("m+[0,]"))

  # Check Candidate Forecast Names
  checkmate::expect_character(colnames(results[[1]]),
                              any.missing = FALSE,
                              len = numb_mods,
                              unique = TRUE)
})
