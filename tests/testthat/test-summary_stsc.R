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
f_signals           <-  replicate(10, rnorm(n = numb_obs, mean = 0, sd = 0.5), )
f_names             <-  paste0("F", as.character(seq_len(ncol(f_signals))))
colnames(f_signals) <-  f_names

# Benchmark
benchmark  <-  dplyr::lag(roll::roll_mean(target_var,
                                          width = length(target_var),
                                          min_obs = 1), n = 1)

# Specify TV-C-Parameter
sample_length  <-  floor(numb_obs / 10)
lambda_grid    <-  c(0.99, 0.999, 1.000)
kappa_grid     <-  c(0.94)
n_cores        <-  1

# Apply TVP-Function
results  <-  hdflex::tvc(target_var,
                         raw_signals,
                         f_signals,
                         lambda_grid,
                         kappa_grid,
                         sample_length,
                         n_cores)

# Assign Results
forecast_tvc      <-  results[[1]]
variance_tvc      <-  results[[2]]

# Define Cut Length
sample_period_idx   <-  (sample_length + 1):numb_obs

# Trim Objects
sub_forecast_tvc    <-  forecast_tvc[sample_period_idx, , drop = FALSE]
sub_variance_tvc    <-  variance_tvc[sample_period_idx, , drop = FALSE]
sub_benchmark       <-  benchmark[sample_period_idx]
sub_target_var      <-  target_var[sample_period_idx]

# Remove Objects
rm(list = c("results", "forecast_tvc", "variance_tvc"))

##### Dynamic Subset Combination #####
# Set DSC Parameter
nr_mods     <-  ncol(sub_forecast_tvc)
gamma_grid  <-  c(0.8, 0.9, 0.95, 0.99, 1)
psi_grid    <-  c(1:25)
delta       <-  0.99
n_cores     <-  1

# Apply DSC-Function
results  <-  hdflex::dsc(gamma_grid,
                         psi_grid,
                         sub_target_var,
                         sub_forecast_tvc,
                         sub_variance_tvc,
                         delta,
                         n_cores)

 # Assign Results
 sub_forecast_stsc    <-  results[[1]]
 sub_variance_stsc    <-  results[[2]]
 sub_chosen_gamma     <-  results[[3]]
 sub_chosen_psi       <-  results[[4]]
 sub_pred_pockets     <-  results[[5]]

 # Define Evaluation Period
 eval_period_idx      <-  51:length(sub_target_var)

 # Trim Objects
 oos_target_var       <-  sub_target_var[eval_period_idx]
 oos_benchmark        <-  sub_benchmark[eval_period_idx]
 oos_forecast_stsc    <-  sub_forecast_stsc[eval_period_idx]
 oos_variance_stsc    <-  sub_variance_stsc[eval_period_idx]
 oos_chosen_gamma     <-  sub_chosen_gamma[eval_period_idx]
 oos_chosen_psi       <-  sub_chosen_psi[eval_period_idx]
 oos_pred_pockets     <-  sub_pred_pockets[eval_period_idx, , drop = FALSE]
 oos_length           <-  length(eval_period_idx)

### Test for no NULLs
test_that("Test whether every input parameter is specified.", {

  testthat::expect_error(summary_stsc(NULL,
                                      oos_benchmark,
                                      oos_forecast_stsc),
                        "not 'NULL'.", fixed = TRUE)

  testthat::expect_error(summary_stsc(oos_target_var,
                                      NULL,
                                      oos_forecast_stsc),
                        "not 'NULL'.", fixed = TRUE)

  testthat::expect_error(summary_stsc(oos_target_var,
                                      oos_benchmark,
                                      NULL),
                        "not 'NULL'.", fixed = TRUE)
})

### Output
test_that("Test whether the output has the right format", {

  # Apply Statistial-Evaluation-Function
  summary_results  <-  summary_stsc(oos_target_var,
                                    oos_benchmark,
                                    oos_forecast_stsc)

  # List Contains Four Elements
  testthat::expect_equal(length(summary_results), 4)

  # Clark-West
  checkmate::expect_number(summary_results[[1]],
                           na.ok = FALSE,
                           finite = TRUE)

  # OOS-R2
  checkmate::expect_number(summary_results[[2]],
                           na.ok = FALSE,
                           upper = 1,
                           finite = TRUE)

  # CSSED
  checkmate::expect_numeric(summary_results[[3]],
                            len = length(oos_target_var),
                            any.missing = FALSE,
                            finite = TRUE)

  # MSE
  checkmate::expect_list(summary_results[[4]],
                         types = "numeric",
                         len = 2)
})