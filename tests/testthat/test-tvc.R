### Simulate Data
set.seed(123)

# Set Dimensions
numb_obs   <-  500
numb_pred  <-  50
numb_forc  <-  10

# Create Random Target-Variable
equity_premium  <-  rnorm(n = numb_obs, mean = 0, sd = 1)

# Create Random Text-Variables
raw_preds            <-  replicate(numb_pred, sample(0:10, numb_obs, rep = TRUE), )
raw_names            <-  paste0("X", as.character(seq_len(ncol(raw_preds))))
colnames(raw_preds)  <-  raw_names

# Create Random Point Forecasts
f_preds            <-  replicate(10, rnorm(n = numb_obs, mean = 0, sd = 0.5), )
f_names            <-  paste0("F", as.character(seq_len(ncol(f_preds))))
colnames(f_preds)  <-  f_names

# Specify TV-C-Parameter
sample_length  <-  100
lambda_grid    <-  c(0.9995, 0.9999, 1.0000)
kappa_grid     <-  c(0.94)
n_cores        <-  1

### Tests on Y
test_that("Test whether y is Numeric Vector", {

  equity_premium  <-  as.data.frame(equity_premium)
  testthat::expect_error(tvc(equity_premium,
                             raw_preds,
                             f_preds,
                             lambda_grid,
                             kappa_grid,
                             sample_length,
                             n_cores),
              "Must be of type 'numeric', not 'data.frame'.", fixed = TRUE)
})

test_that("Test whether y is not NULL", {

  equity_premium  <-  NULL
  testthat::expect_error(tvc(equity_premium,
                             raw_preds,
                             f_preds,
                             lambda_grid,
                             kappa_grid,
                             sample_length,
                             n_cores),
              "Must be of type 'numeric', not 'NULL'.", fixed = TRUE)
})

test_that("Test whether y has only numeric values", {

  equity_premium[10]  <-  "test"
  testthat::expect_error(tvc(equity_premium,
                             raw_preds,
                             f_preds,
                             lambda_grid,
                             kappa_grid,
                             sample_length,
                             n_cores),
              "Must be of type 'numeric', not 'character'.", fixed = TRUE)
})

test_that("Test whether y has no NA-Values", {

  equity_premium[10]  <-  NA
  testthat::expect_error(tvc(equity_premium,
                             raw_preds,
                             f_preds,
                             lambda_grid,
                             kappa_grid,
                             sample_length,
                             n_cores),
              "Contains missing values", fixed = TRUE)
})

### Tests on X
test_that("Test whether x is matrix", {

  raw_preds  <-  as.data.frame(raw_preds)
  testthat::expect_error(tvc(equity_premium,
                             raw_preds,
                             f_preds,
                             lambda_grid,
                             kappa_grid,
                             sample_length,
                             n_cores),
              "Must be of type 'matrix' (or 'NULL')", fixed = TRUE)
})

test_that("Test whether x has the same number of observations as y", {

  raw_preds  <-  raw_preds[1:10, ]
  testthat::expect_error(tvc(equity_premium,
                             raw_preds,
                             f_preds,
                             lambda_grid,
                             kappa_grid,
                             sample_length,
                             n_cores),
                        "Must have exactly", fixed = TRUE)
})

### Tests on f
test_that("Test whether f is matrix", {

  f_preds  <-  as.data.frame(f_preds)
  testthat::expect_error(tvc(equity_premium,
                             raw_preds,
                             f_preds,
                             lambda_grid,
                             kappa_grid,
                             sample_length,
                             n_cores),
                        "Must be of type 'matrix' (or 'NULL')", fixed = TRUE)
})

test_that("Test whether f has the same number of observations as y", {

  f_preds  <-  f_preds[1:10, ]
  testthat::expect_error(tvc(equity_premium,
                             raw_preds,
                             f_preds,
                             lambda_grid,
                             kappa_grid,
                             sample_length,
                             n_cores),
                        "Must have exactly", fixed = TRUE)
})

### Tests on x and f
test_that("Test whether either x or f is provided", {

  raw_preds  <-  NULL
  f_preds    <-  NULL
  testthat::expect_error(tvc(equity_premium,
                             raw_preds,
                             f_preds,
                             lambda_grid,
                             kappa_grid,
                             sample_length,
                             n_cores),
            "Assertion failed. One of the following must apply:
 * checkmate::checkMatrix(x): Must be of type 'matrix', not 'NULL'
 * checkmate::checkMatrix(f): Must be of type 'matrix', not 'NULL'",
            fixed = TRUE)
})

### Output
test_that("Test whether the output has the right format", {

  # Apply TVP-Function
  results  <-  tvc(equity_premium,
                   raw_preds,
                   f_preds,
                   lambda_grid,
                   kappa_grid,
                   sample_length,
                   n_cores)

  # List Contains Three Elements
  testthat::expect_equal(length(results), 3)

  # Number of Models
  numb_mods  <-  length(lambda_grid) * length(kappa_grid) * numb_pred +
                 length(kappa_grid) * numb_forc

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

  # Length of Model Names
  checkmate::expect_character(results[[3]],
                              any.missing = FALSE,
                              unique = TRUE,
                              len = numb_mods)
})
