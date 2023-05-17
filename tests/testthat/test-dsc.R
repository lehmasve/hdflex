### Simulate Data
set.seed(123)

# Set Dimensions
numb_obs   <-  500
numb_mods  <-  50

# Create Random Target-Variable
target_var  <-  rnorm(n = numb_obs, mean = 0, sd = 1)

# Create Random Candidate-Forecast-Matrix
forecast_tvc            <-  replicate(numb_mods, rnorm(n = numb_obs, mean = 0, sd = 1), )
f_tvc_names             <-  paste0("X", as.character(seq_len(ncol(forecast_tvc))))
colnames(forecast_tvc)  <-  f_tvc_names

# Create Random Candidate-Variance-Matrix
variance_tvc            <-  replicate(numb_mods, abs(rnorm(n = numb_obs, mean = 0, sd = 1)), )
v_tvc_names             <-  paste0("F", as.character(seq_len(ncol(variance_tvc))))
colnames(variance_tvc)  <-  v_tvc_names

# Set DSC Parameter
nr_mods     <-  numb_mods
gamma_grid  <-  c(0.8, 0.9, 0.95, 0.99, 1)
psi_grid    <-  c(1, 2, 3)
delta       <-  0.99
n_cores     <-  1

### Test for no NULLs
test_that("Test whether every input parameter is specified.", {

  testthat::expect_error(dsc(NULL,
                             psi_grid,
                             target_var,
                             forecast_tvc,
                             variance_tvc,
                             delta,
                             n_cores),
              "not 'NULL'.", fixed = TRUE)

  testthat::expect_error(dsc(gamma_grid,
                             NULL,
                             target_var,
                             forecast_tvc,
                             variance_tvc,
                             delta,
                             n_cores),
            "not 'NULL'.", fixed = TRUE)

  testthat::expect_error(dsc(gamma_grid,
                             psi_grid,
                             NULL,
                             forecast_tvc,
                             variance_tvc,
                             delta,
                             n_cores),
              "not 'NULL'.", fixed = TRUE)

  testthat::expect_error(dsc(gamma_grid,
                             psi_grid,
                             target_var,
                             NULL,
                             variance_tvc,
                             delta,
                             n_cores),
              "not 'NULL'.", fixed = TRUE)

  testthat::expect_error(dsc(gamma_grid,
                             psi_grid,
                             target_var,
                             forecast_tvc,
                             NULL,
                             delta,
                             n_cores),
              "not 'NULL'.", fixed = TRUE)

  testthat::expect_error(dsc(gamma_grid,
                             psi_grid,
                             target_var,
                             forecast_tvc,
                             variance_tvc,
                             NULL,
                             n_cores),
              "not 'NULL'.", fixed = TRUE)

  testthat::expect_error(dsc(gamma_grid,
                             psi_grid,
                             target_var,
                             forecast_tvc,
                             variance_tvc,
                             delta,
                             NULL),
              "not 'NULL'.", fixed = TRUE)
})

### Output
test_that("Test whether the output has the right format", {

  # Apply TVP-Function
  results  <-  dsc(gamma_grid,
                   psi_grid,
                   target_var,
                   forecast_tvc,
                   variance_tvc,
                   delta,
                   n_cores)

  # List Contains Five Elements
  testthat::expect_equal(length(results), 5)

  # Number of Models
  numb_combs  <-  length(gamma_grid) * length(psi_grid)

  # Number of Forecasts
  checkmate::expect_numeric(results[[1]],
                            len = numb_obs,
                            any.missing = FALSE,
                            finite = TRUE)

  # Number of Variances
  checkmate::expect_numeric(results[[2]],
                            len = numb_obs,
                            any.missing = FALSE,
                            lower = 0,
                            finite = TRUE)
  # Length of Gamma-Vector
  checkmate::expect_numeric(results[[3]],
                            len = numb_obs,
                            any.missing = FALSE,
                            lower = 0,
                            finite = TRUE)
  # Length of Psi-Vector
  checkmate::expect_numeric(results[[4]],
                            len = numb_obs,
                            any.missing = FALSE,
                            lower = 0,
                            finite = TRUE)

  # Dimension of selected Candidate Forecasts
  checkmate::expect_matrix(results[[5]],
                           mode = "integerish",
                           nrows = numb_obs,
                           ncols = numb_mods)
})
