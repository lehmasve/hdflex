
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hdflex <a href='https://github.com/lehmasve/hdflex'><img src='man/figures/logo.png' align="right" height="160" /></a>

⁠<!-- badges: start -->⁠ [![CRAN
Version](https://www.r-pkg.org/badges/version/hdflex)](https://CRAN.R-project.org/package=hdflex)
[![DOI:10.2139/ssrn.4342487](http://img.shields.io/badge/DOI-10.2139/ssrn.4342487-163870.svg)](https://dx.doi.org/10.2139/ssrn.4342487)
[![R-CMD-check](https://github.com/lehmasve/hdflex/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lehmasve/hdflex/actions/workflows/R-CMD-check.yaml)
[![Total
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/hdflex?color=orange)](https://CRAN.R-project.org/package=hdflex)
[![codecov](https://codecov.io/gh/lehmasve/hdflex/graph/badge.svg?token=leKtsb0Kub)](https://app.codecov.io/gh/lehmasve/hdflex)
⁠<!-- badges: end -->⁠

## About

This package contains the forecasting algorithm developed by [Adämmer,
Lehmann and Schüssler (2023)](https://dx.doi.org/10.2139/ssrn.4342487).
When using the package, please remember to cite the paper.

The package comprises three functions:

- `stsc()` can be used to directly apply the
  “Signal-Transformed-Subset-Combination” forecasting algorithm
  described in [Adämmer, Lehmann and Schüssler
  (2023)](https://dx.doi.org/10.2139/ssrn.4342487).

- `tvc()` can be used to transform predictive signals into univariate
  density forecasts via time-varying coefficient models, where each
  model generates a conditionally gaussian predictive density for each
  signal at each point in time (first part of the STSC algorithm).

- `dsc()` can be used to dynamically select a subset of candidate
  forecast models for each period based on their past performance in
  terms of density forecast accuracy (second part of the STSC
  algorithm).

## Installation

You can install the released version of **hdflex** from
[CRAN](https://CRAN.R-project.org):

``` r
install.packages("hdflex")
```

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

- On macOS you need Xcode and a Fortran compiler - for more details see
  [Compiler](https://mac.r-project.org/tools/).

## Usage

Example using the `stsc()` function:

``` r
#########################################################
######### Forecasting quarterly U.S. inflation ##########
#### Please see Koop & Korobilis (2023) for further  ####
#### details regarding the data & external forecasts ####
#########################################################

# Load Package
library("hdflex")
library("ggplot2")
library("cowplot")

########## Get Data ##########
# Load Package Data
inflation_data <- inflation_data

# Set Target Variable
y <- inflation_data[,  1]

# Set 'P-Signals'
X <- inflation_data[, 2:442]

# Set 'F-Signals'
Ext_F <- inflation_data[, 443:462]

# Get Dates and Number of Observations
tdates <- rownames(inflation_data)
tlength <- length(tdates)

# First complete observation (no missing values)
first_complete <- which(complete.cases(inflation_data))[1]

########## Rolling AR2-Benchmark ##########
# Set up matrix for predictions
benchmark <- matrix(NA, nrow = tlength,
                    ncol = 1, dimnames = list(tdates, "AR2"))

# Set Window-Size (15 years of quarterly data)
window_size <- 15 * 4

# Time Sequence
t_seq <- seq(window_size, tlength - 1)

# Loop with rolling window
for (t in t_seq) {

  # Split Data for Training Train Data
  x_train <- cbind(int = 1, X[(t - window_size + 1):t, 1:2])
  y_train <- y[(t - window_size + 1):t]

  # Split Data for Prediction
  x_pred <- cbind(int = 1, X[t + 1, 1:2, drop = FALSE])

  # Fit AR-Model
  model_ar <- .lm.fit(x_train, y_train)

  # Predict and store in benchmark matrix
  benchmark[t + 1, ] <- x_pred %*% model_ar$coefficients
}

########## STSC ##########
# Set TV-C-Parameter
init <- 5 * 4
lambda_grid <- c(0.90, 0.95, 1.00)
kappa_grid <- c(0.94, 0.96, 0.98)
bias <- TRUE

# Set DSC-Parameter
gamma_grid <- c(0.40, 0.50, 0.60, 0.70, 0.80, 0.90,
                0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00)
n_tvc <- (ncol(X) + ncol(Ext_F)) * length(lambda_grid) * length(kappa_grid)
psi_grid <- c(1:100, sapply(1:4, function(i) floor(i * n_tvc / 4)))
delta <- 0.95
burn_in <- first_complete + init / 2
burn_in_dsc <- 1
metric <- 5
equal_weight <- TRUE
incl <- NULL
parallel <- FALSE
n_threads <- NULL

# Apply STSC-Function
results <- hdflex::stsc(y,
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
                        NULL)

########## Evaluation ##########
# Define Evaluation Period (OOS-Period)
eval_period <- which(tdates >= "1991-04-01" & tdates <= "2021-12-01")

# Apply Evaluation Summary for STSC
eval_results <- summary(obj = results, eval_period = eval_period)

# Calculate (Mean-)Squared-Errors for AR2-Benchmark
se_ar2 <- (y[eval_period] - benchmark[eval_period, 1])^2
mse_ar2 <- mean(se_ar2)

# Create Cumulative Squared Error Differences (CSSED) Plot
cssed <- cumsum(se_ar2 - eval_results$MSE[[2]])
plot_cssed <- ggplot(
  data.frame(eval_period, cssed),
  aes(x = eval_period, y = cssed)
) +
  geom_line() +
  ylim(-0.0008, 0.0008) +
  ggtitle("Cumulative Squared Error Differences") +
  xlab("Time Index") +
  ylab("CSSED") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    axis.ticks = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5)
  )

# Show Plots
options(repr.plot.width = 15, repr.plot.height = 15)
plots_list <- eval_results$Plots
plots_list <- c(list(plot_cssed), plots_list)
cowplot::plot_grid(plotlist = plots_list, ncol = 2, nrow = 3, align = "hv")

# Relative MSE
print(paste("Relative MSE:", round(eval_results$MSE[[1]] / mse_ar2, 4)))
```

Same example using the `tvc()` and `dsc()` functions:

``` r
#########################################################
######### Forecasting quarterly U.S. inflation ##########
#### Please see Koop & Korobilis (2023) for further  ####
#### details regarding the data & external forecasts ####
#########################################################

# Load Package
library("hdflex")
library("ggplot2")
library("cowplot")

########## Get Data ##########
# Load Package Data
inflation_data <- inflation_data

# Set Target Variable
y <- inflation_data[,  1]

# Set 'P-Signals'
X <- inflation_data[, 2:442]

# Set 'F-Signals'
Ext_F <- inflation_data[, 443:462]

# Get Dates and Number of Observations
tdates <- rownames(inflation_data)
tlength <- length(tdates)

# First complete observation (no missing values)
first_complete <- which(complete.cases(inflation_data))[1]

########## Rolling AR2-Benchmark ##########
# Set up matrix for predictions
benchmark <- matrix(NA, nrow = tlength,
                    ncol = 1, dimnames = list(tdates, "AR2"))

# Set Window-Size (15 years of quarterly data)
window_size <- 15 * 4

# Time Sequence
t_seq <- seq(window_size, tlength - 1)

# Loop with rolling window
for (t in t_seq) {

  # Split Data for Training Train Data
  x_train <- cbind(int = 1, X[(t - window_size + 1):t, 1:2])
  y_train <- y[(t - window_size + 1):t]

  # Split Data for Prediction
  x_pred <- cbind(int = 1, X[t + 1, 1:2, drop = FALSE])

  # Fit AR-Model
  model_ar <- .lm.fit(x_train, y_train)

  # Predict and store in benchmark matrix
  benchmark[t + 1, ] <- x_pred %*% model_ar$coefficients
}

########## STSC ##########
### Part 1: TVC-Function
# Set TV-C-Parameter
init <- 5 * 4
lambda_grid <- c(0.90, 0.95, 1.00)
kappa_grid <- c(0.94, 0.96, 0.98)
bias <- TRUE

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

# First complete forecast period (no missing values)
sub_period <- seq(which(complete.cases(forecast_tvc))[1], tlength)

### Part 2: DSC-Function
# Set DSC-Parameter
gamma_grid <- c(0.40, 0.50, 0.60, 0.70, 0.80, 0.90,
                0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00)
psi_grid <- c(1:100, sapply(1:4, function(i) floor(i * ncol(forecast_tvc) / 4)))
delta <- 0.95
burn_in <- (init / 2) + 1
burn_in_dsc <- 1
metric <- 5
equal_weight <- TRUE
incl <- NULL

# Apply DSC-Function
dsc_results <- hdflex::dsc(y[sub_period],
                           forecast_tvc[sub_period, , drop = FALSE],
                           variance_tvc[sub_period, , drop = FALSE],
                           gamma_grid,
                           psi_grid,
                           delta,
                           burn_in,
                           burn_in_dsc,
                           metric,
                           equal_weight,
                           incl,
                           NULL)

# Assign DSC-Results
pred_stsc <- dsc_results$Forecasts$Point_Forecasts
var_stsc <- dsc_results$Forecasts$Variance_Forecasts

########## Evaluation ##########
# Define Evaluation Period (OOS-Period)
eval_period <- which(tdates[sub_period] >= "1991-04-01" & tdates[sub_period] <= "2021-12-01")

# Get Evaluation Summary for STSC
eval_results <- summary(obj = dsc_results, eval_period = eval_period)

# Calculate (Mean-)Squared-Errors for AR2-Benchmark
oos_y <- y[sub_period][eval_period]
oos_benchmark <- benchmark[sub_period[eval_period], , drop = FALSE]
se_ar2 <- (oos_y - oos_benchmark)^2
mse_ar2 <- mean(se_ar2)

# Create Cumulative Squared Error Differences (CSSED) Plot
cssed <- cumsum(se_ar2 - eval_results$MSE[[2]])
plot_cssed <- ggplot(
  data.frame(eval_period, cssed),
  aes(x = eval_period, y = cssed)
) +
  geom_line() +
  ylim(-0.0008, 0.0008) +
  ggtitle("Cumulative Squared Error Differences") +
  xlab("Time Index") +
  ylab("CSSED") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    axis.ticks = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5)
  )

# Show Plots
options(repr.plot.width = 15, repr.plot.height = 15)
plots_list <- eval_results$Plots
plots_list <- c(list(plot_cssed), plots_list)
cowplot::plot_grid(plotlist = plots_list, ncol = 2, nrow = 3, align = "hv")

# Relative MSE
print(paste("Relative MSE:", round(eval_results$MSE[[1]] / mse_ar2, 4)))
```

### Authors

Philipp Adämmer, Sven Lehmann and Rainer Schüssler

### License

GPL (\>= 2)
