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

# Benchmark
benchmark  <-  dplyr::lag(roll::roll_mean(equity_premium,
                                          width = length(equity_premium),
                                          min_obs = 1), n = 1)

# Specify TV-C-Parameter
sample_length  <-  floor(numb_obs / 10)
lambda_grid    <-  c(0.9995, 0.9999, 1.0000)
kappa_grid     <-  c(0.94)
n_cores        <-  1

# Apply TVP-Function
results  <-  hdflex::tvc(equity_premium,
                         raw_preds,
                         f_preds,
                         lambda_grid,
                         kappa_grid,
                         sample_length,
                         n_cores)

# Assign Results
forecast_tvc      <-  results[[1]]
variance_tvc      <-  results[[2]]
model_names_tvc   <-  results[[3]]

# Define Cut Length
nr_drp              <-  dim(forecast_tvc)[1] - dim(na.omit(forecast_tvc))[1] + 1
sample_period_idx   <-  (nr_drp + sample_length):numb_obs

# Trim Objects
sub_forecast_tvc    <-  forecast_tvc[sample_period_idx, , drop = FALSE]
sub_variance_tvc    <-  variance_tvc[sample_period_idx, , drop = FALSE]
sub_benchmark       <-  benchmark[sample_period_idx]
sub_equity_premium  <-  equity_premium[sample_period_idx]

# Remove Objects
rm(list = c("results", "forecast_tvc", "variance_tvc"))

##### Dynamic Subset Combination #####
# Set DSC Parameter
nr_mods     <-  length(model_names_tvc)
gamma_grid  <-  c(0.9, 0.95, 0.99, 1)
psi_grid    <-  c(1, 2, 3)
delta       <-  0.9992
n_cores     <-  1

# Apply DSC-Function
results  <-  hdflex::dsc(gamma_grid,
                         psi_grid,
                         sub_equity_premium,
                         sub_forecast_tvc,
                         sub_variance_tvc,
                         delta,
                         n_cores)

# Assign Results
sub_forecast_bcomb    <-  results[[1]]
sub_variance_bcomb    <-  results[[2]]
model_names_comb      <-  results[[3]]
sub_chosen_parameter  <-  results[[4]]
sub_models_idx        <-  results[[5]]

# Define Evaluation Period
eval_period_idx      <-  50:length(sub_equity_premium)

# Trim Objects
oos_equity_premium   <-  sub_equity_premium[eval_period_idx]
oos_benchmark        <-  sub_benchmark[eval_period_idx]
oos_forecast_bcomb   <-  sub_forecast_bcomb[eval_period_idx]
oos_variance_bcomb   <-  sub_variance_bcomb[eval_period_idx]
oos_models_idx       <-  lapply(seq_along(model_names_comb), function(i) {
                                    sub_models_idx[[i]][eval_period_idx]})
oos_chosen_parameter <-  sub_chosen_parameter[eval_period_idx, , drop = FALSE]
oos_dates            <-  seq(as.Date("1989-12-01"),
                             by = "day", length.out = length(eval_period_idx))

# Assign Names
names(oos_forecast_bcomb)      <-  oos_dates
names(oos_variance_bcomb)      <-  oos_dates

### Test for no NULLs
test_that("Test whether every input parameter is specified.", {

  testthat::expect_error(eval_dsc(NULL,
                                  oos_benchmark,
                                  oos_forecast_bcomb,
                                  oos_dates,
                                  oos_chosen_parameter,
                                  model_names_tvc,
                                  oos_models_idx),
                        "not 'NULL'.", fixed = TRUE)

  testthat::expect_error(eval_dsc(oos_equity_premium,
                                  NULL,
                                  oos_forecast_bcomb,
                                  oos_dates,
                                  oos_chosen_parameter,
                                  model_names_tvc,
                                  oos_models_idx),
                        "not 'NULL'.", fixed = TRUE)

  testthat::expect_error(eval_dsc(oos_equity_premium,
                                  oos_benchmark,
                                  NULL,
                                  oos_dates,
                                  oos_chosen_parameter,
                                  model_names_tvc,
                                  oos_models_idx),
                        "not 'NULL'.", fixed = TRUE)

  testthat::expect_error(eval_dsc(oos_equity_premium,
                                  oos_benchmark,
                                  oos_forecast_bcomb,
                                  NULL,
                                  oos_chosen_parameter,
                                  model_names_tvc,
                                  oos_models_idx),
                        "not 'NULL'.", fixed = TRUE)

  testthat::expect_error(eval_dsc(oos_equity_premium,
                                  oos_benchmark,
                                  oos_forecast_bcomb,
                                  oos_dates,
                                  NULL,
                                  model_names_tvc,
                                  oos_models_idx),
                        "not 'NULL'.", fixed = TRUE)

  testthat::expect_error(eval_dsc(oos_equity_premium,
                                  oos_benchmark,
                                  oos_forecast_bcomb,
                                  oos_dates,
                                  oos_chosen_parameter,
                                  NULL,
                                  oos_models_idx),
                        "not 'NULL'.", fixed = TRUE)

  testthat::expect_error(eval_dsc(oos_equity_premium,
                                  oos_benchmark,
                                  oos_forecast_bcomb,
                                  oos_dates,
                                  oos_chosen_parameter,
                                  model_names_tvc,
                                  NULL),
                        "not 'NULL'.", fixed = TRUE)
})

### Output
test_that("Test whether the output has the right format", {

  # Apply Statistial-Evaluation-Function
  eval_results  <-  eval_dsc(oos_equity_premium,
                             oos_benchmark,
                             oos_forecast_bcomb,
                             oos_dates,
                             oos_chosen_parameter,
                             model_names_tvc,
                             oos_models_idx)

  # List Contains Four Elements
  testthat::expect_equal(length(eval_results), 4)

  # Clark-West
  checkmate::expect_number(eval_results[[1]],
                           na.ok = FALSE,
                           finite = TRUE)

  # OOS-R2
  checkmate::expect_number(eval_results[[2]],
                           na.ok = FALSE,
                           upper = 1,
                           finite = TRUE)

  # CSED
  checkmate::expect_numeric(eval_results[[3]],
                            len = length(oos_equity_premium),
                            any.missing = FALSE,
                            finite = TRUE)

  # Pockets
  checkmate::expect_matrix(eval_results[[4]],
                           mode = "integerish",
                           nrows = length(oos_equity_premium),
                           ncols = (numb_pred + numb_forc))
})