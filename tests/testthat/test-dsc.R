context("dsc")

### Simulate Data
set.seed(123)

# Set Dimensions
numb_obs   <-  500
numb_mods  <-  50

# Create Random Target-Variable
equity_premium  <-  rnorm(n = numb_obs, mean = 0, sd = 1)

# Create Forecast-Matrix
forecast_tvc            <-  replicate(numb_mods, rnorm(n = numb_obs, mean = 0, sd = 1), )
f_tvc_names             <-  paste0("X", as.character(seq_len(ncol(forecast_tvc))))
colnames(forecast_tvc)  <-  f_tvc_names

# Create Variance-Matrix
variance_tvc            <-  replicate(numb_mods, abs(rnorm(n = numb_obs, mean = 0, sd = 1)), )
v_tvc_names             <-  paste0("F", as.character(seq_len(ncol(variance_tvc))))
colnames(variance_tvc)  <-  v_tvc_names

# Set DSC Parameter
nr_mods     <-  numb_mods
gamma_grid  <-  c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.995, 0.9999, 1)
psi_grid    <-  c(1, 2, 3)
ew_grid     <-  c(0, 1)
delta       <-  0.9992
n_cores     <-  1

### Test for no NULLs
test_that("Test whether every input parameter is specified.", {

  testthat::expect_error(dsc(NULL,
                             psi_grid,
                             ew_grid,
                             equity_premium,
                             forecast_tvc,
                             variance_tvc,
                             delta,
                             n_cores),
              "not 'NULL'.", fixed = TRUE)

  testthat::expect_error(dsc(gamma_grid,
                             NULL,
                             ew_grid,
                             equity_premium,
                             forecast_tvc,
                             variance_tvc,
                             delta,
                             n_cores),
            "not 'NULL'.", fixed = TRUE)

  testthat::expect_error(dsc(gamma_grid,
                             psi_grid,
                             NULL,
                             equity_premium,
                             forecast_tvc,
                             variance_tvc,
                             delta,
                             n_cores),
              "not 'NULL'.", fixed = TRUE)

  testthat::expect_error(dsc(gamma_grid,
                             psi_grid,
                             ew_grid,
                             NULL,
                             forecast_tvc,
                             variance_tvc,
                             delta,
                             n_cores),
              "not 'NULL'.", fixed = TRUE)

  testthat::expect_error(dsc(gamma_grid,
                             psi_grid,
                             ew_grid,
                             equity_premium,
                             NULL,
                             variance_tvc,
                             delta,
                             n_cores),
              "not 'NULL'.", fixed = TRUE)

  testthat::expect_error(dsc(gamma_grid,
                             psi_grid,
                             ew_grid,
                             equity_premium,
                             forecast_tvc,
                             NULL,
                             delta,
                             n_cores),
              "not 'NULL'.", fixed = TRUE)

  testthat::expect_error(dsc(gamma_grid,
                             psi_grid,
                             ew_grid,
                             equity_premium,
                             forecast_tvc,
                             variance_tvc,
                             NULL,
                             n_cores),
              "not 'NULL'.", fixed = TRUE)

  testthat::expect_error(dsc(gamma_grid,
                             psi_grid,
                             ew_grid,
                             equity_premium,
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
                   ew_grid,
                   equity_premium,
                   forecast_tvc,
                   variance_tvc,
                   delta,
                   n_cores)

  # List Contains Five Elements
  testthat::expect_equal(length(results), 5)

  # Number of Models
  numb_combs  <-  length(gamma_grid) * length(psi_grid) * length(ew_grid)

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
  # Length of Model Names
  checkmate::expect_character(results[[3]],
                              any.missing = FALSE,
                              unique = TRUE,
                              len = numb_combs)
  # Type and Dimension
  checkmate::expect_matrix(results[[4]],
                           mode = "logical",
                           nrows = numb_obs,
                           ncols = numb_combs)

  # Length of List
  checkmate::expect_list(results[[5]],
                         len = numb_combs)
})
