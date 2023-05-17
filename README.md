
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hdflex <a href='https://github.com/lehmasve/hdflex'><img src='man/figures/logo.png' align="right" height="180" /></a>

## About

This package contains the forecasting algorithm developed by [Adämmer,
Lehmann and Schüssler
(2023)](https://www.researchgate.net/publication/367531209_Local_Predictability_in_High_Dimensions).
Please cite the paper when using the package.

## Installation

You can install **hdflex** from
[GitHub](https://github.com/lehmasve/hdflex):

``` r
# install.packages("devtools")
devtools::install_github("https://github.com/lehmasve/hdflex")
```

The package compiles some C++ source code for installation, which is why
you need the appropriate compilers:

- On Windows you need
  [Rtools](https://cran.r-project.org/bin/windows/Rtools/) available
  from CRAN.

- On macOS you need the very least Xcode and a Fortran compiler - for
  more details see [Compiler](https://mac.r-project.org/tools/).

## Main Features

We propose a new forecasting strategy (STSC) that captures
predictability patterns of high-dimensional signals, especially those
that are locally concentrated in time. In each period, an aggregate
density forecast is determined by a subset of candidate density
forecasts, where each candidate is based on a state-space model. Our
approach allows for several stylized features of (macro-)economic and
financial returns, including time-varying predictability patterns and
time-varying volatility. While the approach is capable of exploring vast
model spaces, including datasets that contain mostly irrelevant
predictors, it is parsimoniously parameterized, thereby reducing
estimation error and maintaining transparency. Furthermore, it is very
fast by using empirically reasonable approximations, online updating of
model parameters, and parallel computing.

The package comprises three function:

- `tvc()` can be used to compute density forecasts based on univariate
  time-varying coefficient (TV-C) models in state-space form

- `dsc()` can be used to dynamically generate forecast combinations from
  a subset of candidate density forecasts (dynamic subset combination).

- `summary_stsc()` returns a statistical summary for the forecasting
  results. It provides statistical measures such as
  Clark-West-Statistic, OOS-R2, Mean-Squared-Error and Cumulated Sum of
  Squared-Error-Differences.

## Usage

``` r
#########################################################
######### Forecasting quarterly U.S. inflation ##########
#### Please see Koop & Korobilis (2023) for further  ####
#### details regarding the data & External Forecasts ####
#########################################################

# Load Packages
library("tidyverse")
#> ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
#> ✔ dplyr     1.1.2     ✔ readr     2.1.4
#> ✔ forcats   1.0.0     ✔ stringr   1.5.0
#> ✔ ggplot2   3.4.2     ✔ tibble    3.2.1
#> ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
#> ✔ purrr     1.0.1     
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
#> ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
library("hdflex")

# Set Target-Variables
target_var_names  <- c("GDPCTPI", "PCECTPI", "CPIAUCSL", "CPILFESL")

# Load Data
data  <-  inflation_data

# Loop over Target Variables
results <-  do.call("rbind", lapply(X = seq_along(target_var_names), FUN = function(p) { 

    # Y-Column-Name
    y_target    <-  paste(target_var_names[p], "h_1", sep = "_")
    y_signal    <-  target_var_names[p]
    not_target  <-  setdiff(paste(target_var_names, "h_1", sep = "_"), y_target)

    # Create Forecast-Dataset
    dataset  <-  data                                                                                     %>%
                  select(-any_of(c(y_signal, not_target)))                                                %>%
                  mutate(across(all_of(y_target), .fns = list("lag_1" = ~dplyr::lag(., 1))), .after = 2)  %>%
                  mutate(across(all_of(y_target), .fns = list("lag_2" = ~dplyr::lag(., 2))), .after = 3)  %>%
                  mutate(across(-c(1, 2, 3, 4), dplyr::lag))                                              %>%
                  slice(-c(1:3))                                                                          %>%
                  column_to_rownames("Date")                                                              %>%
                  as.matrix()
    
    # Get Dates & Length
    full_dates   <-  rownames(dataset)
    full_length  <-  length(full_dates)

    # Create Time-Sequence for loop
    T_full      <-  full_length -1
    T_sequence  <-  122:T_full 
    
    ### Benchmark Model ###
    # Create Result Matrices
    preds_ar2  <-  matrix(NA, ncol = 1, nrow = full_length, dimnames = list(full_dates, "AR"))
    se_ar2     <-  matrix(NA, ncol = 1, nrow = full_length, dimnames = list(full_dates, "AR"))
        
    # Loop over t
    for(t in T_sequence) {
            
        ### Pre-Process Data ###
        # Train Data
        x_train     <-  scale(dataset[1:t, -1, drop = FALSE])
        y_train     <-        dataset[1:t,  1, drop = FALSE]
  
        # Predict Data
        x_pred      <-  scale(dataset[1:(t+1), -1, drop = FALSE])[(t+1), , drop = FALSE]
        y_pred      <-        dataset[t+1, 1]

        ### Model 1: AR(2) ###
        # Train Data
        x_train_ar  <-  cbind(int = 1, x_train[, c(1:2), drop = FALSE])

        # Predict Data
        x_pred_ar   <-  cbind(int = 1,  x_pred[, c(1:2), drop = FALSE])

        # Fit Regressions
        model_ar    <-  .lm.fit(x_train_ar, y_train)

        # Predict & Combine
        pred                  <-  model_ar$coefficients %*% x_pred_ar[,]
        preds_ar2[t+1, "AR"]  <-  pred
        se_ar2[t+1, "AR"]     <-  (y_pred - pred) ** 2
    }

    ##### TV-C Models ##### 
    # Set Target Variable
    Y  <-  dataset[,  1, drop = FALSE]
    
    # Set 'Simple' Signals
    X  <-  dataset[, -1, drop = FALSE]

    # Load External Point Forecasts (Koop & Korobilis 2023)
    F  <-  get(paste0("Ext_PF_", target_var_names[p]))

    # Set TV-C-Parameter
    sample_length  <-  4 * 3
    lambda_grid    <-  c(0.90, 0.95, 0.99, 0.999, 1) 
    kappa_grid     <-  0.98
    n_cores        <-  1

    # Apply TV-C-Function
    results  <-  hdflex::tvc(Y, 
                             X,
                             F, 
                             lambda_grid,
                             kappa_grid,
                             sample_length,
                             n_cores)
                             
    # Assign Results
    forecast_tvc      <-  results[[1]]
    variance_tvc      <-  results[[2]]
    model_names_tvc   <-  colnames(forecast_tvc)
    
    # Define Cut Length and Trim Objects
    sample_period_idx  <-  80:full_length   
    sub_forecast_tvc   <-  forecast_tvc[sample_period_idx, , drop = FALSE]
    sub_variance_tvc   <-  variance_tvc[sample_period_idx, , drop = FALSE]
    sub_Y              <-  Y[sample_period_idx, , drop = FALSE]
    sub_dates          <-  full_dates[sample_period_idx]
    sub_length         <-  length(sub_dates)

    ##### Dynamic Subset Combination ##### 
    # Set DSC-Parameter
    nr_mods     <-  ncol(sub_forecast_tvc)
    gamma_grid  <-  c(0.4, 0.5, 0.6, 0.7, 0.8, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.999, 1.00)
    psi_grid    <-  c(1:100)
    delta       <-  0.95
    n_cores     <-  1

    # Apply DSC-Function
    results  <-  hdflex::dsc(gamma_grid,
                             psi_grid,
                             sub_Y,
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
    eval_date_start      <-  "1991-01-01" 
    eval_date_end        <-  "2021-12-31"     
    eval_period_idx      <-  which(sub_dates > eval_date_start & sub_dates <= eval_date_end)

    # Trim Objects
    oos_Y                <-  sub_Y[eval_period_idx, ]
    oos_benchmark        <-  preds_ar2[rownames(preds_ar2) > eval_date_start, "AR"] 
    oos_forecast_stsc    <-  sub_forecast_stsc[eval_period_idx]
    oos_variance_stsc    <-  sub_variance_stsc[eval_period_idx]
    oos_chosen_gamma     <-  sub_chosen_gamma[eval_period_idx]
    oos_chosen_psi       <-  sub_chosen_psi[eval_period_idx]
    oos_pred_pockets     <-  sub_pred_pockets[eval_period_idx, , drop = FALSE]
    oos_length           <-  length(eval_period_idx)
    oos_dates            <-  sub_dates[eval_period_idx]

    # Add Dates
    names(oos_forecast_stsc)   <-  oos_dates
    names(oos_variance_stsc)   <-  oos_dates
    names(oos_chosen_gamma)    <-  oos_dates
    names(oos_chosen_psi)      <-  oos_dates
    rownames(oos_pred_pockets) <-  oos_dates

    ##### Evaluate ##### 
    # Apply Statistial-Evaluation-Function
    summary_results  <-  summary_stsc(oos_Y,
                                      oos_benchmark,
                                      oos_forecast_stsc) 
     # Assign MSE-Results
    mse  <-  summary_results[[4]]

    # Return
    return(c(mse[[2]], mse[[1]]))

}))

# MSE-Results
dimnames(results)  <-  list(target_var_names, c("AR", "STSC"))
round(results / results[, 1], 4)
#>          AR   STSC
#> GDPCTPI   1 0.8818
#> PCECTPI   1 0.6899
#> CPIAUCSL  1 0.9247
#> CPILFESL  1 0.9093
```

### Authors

Philipp Adämmer, Sven Lehmann and Rainer Schüssler

### License

MIT
