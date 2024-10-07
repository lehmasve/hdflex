#' @name tvc
#' @title Compute density forecasts using univariate time-varying coefficient (TV-C) models
#' @description The `tvc()` function generates density forecasts
#' based on univariate time-varying coefficient models in state-space form.
#' Each forecasting model includes an intercept and one predictive signal,
#' which can either be a 'P-signal' or 'F-signal'.
#' All models are estimated independently and both
#' estimation and forecasting are carried out recursively.
#' @param y A matrix of dimension `T * 1` or numeric vector of length `T`
#' containing the observations of the target variable.
#' @param X A matrix with `T` rows containing the lagged 'P-signals'
#' in each column. Use `NULL` if no (external) 'P-signal' is to be included.
#' @param Ext_F A matrix with `T` rows containing the (external) 'F-signals'
#' in each column. For 'F-Signals', the slope of the TV-C models is fixed to 1.
#' Use `NULL` if no (external) 'F-signal' is to be included.
#' @param init An integer that denotes the number of observations used
#' to initialize the observational variance and the coefficients' variance
#' in the TV-C models.
#' @param lambda_grid A numeric vector which takes values between 0 and 1
#' denoting the discount factor(s) that control the dynamics of the time-varying
#' coefficients. Each signal in combination with each value of
#' lambda provides a separate candidate forecast.
#' Constant coefficients are nested for the case `lambda = 1`.
#' @param kappa_grid A numeric vector which takes values between 0 and 1
#' to accommodate time-varying volatility in the TV-C models.
#' The observational variance is estimated via
#' Exponentially Weighted Moving Average and kappa denotes the underlying
#' decay factor. Constant variance is nested for the case `kappa = 1`.
#' Each signal in combination with each value of
#' kappa provides a separate candidate density forecast.
#' For the values of kappa, we follow the recommendation
#' of RiskMetrics (Reuters, 1996).
#' @param bias A boolean to indicate whether the TV-C-models
#' allow for a bias correction to F-signals.
#' `TRUE` allows for a time-varying intercept, and `FALSE` sets (and fixes)
#' the intercept to 0.
#' @return A list containing:
#' \describe{
#'  \item{Forecasts}{A list containing:
#'   \describe{
#'     \item{Realization: }{
#'     A vector with the actual values of the target variable.
#'     }
#'     \item{Point_Forecasts: }{
#'     A vector with the first moments of the predictive densities.
#'     }
#'     \item{Variance_Forecasts: }{
#'     A vector with the second moments of the predictive densities.
#'     }
#'   }
#' }
#' \item{Model}{A list containing:
#'   \describe{
#'     \item{Lambda_grid}{
#'     The grid of lambda values used in the model.
#'     }
#'     \item{
#'     Kappa_grid}{The grid of kappa values used in the model.
#'     }
#'     \item{Init}{
#'     The init value used in the model.
#'     }
#'     \item{Bias}{
#'     A boolean indicating if bias correct was applied to F-signals.
#'     }
#'   }
#'  }
#' }
#' @seealso \url{https://github.com/lehmasve/hdflex#readme}
#' @author Philipp Adämmer, Sven Lehmann, Rainer Schüssler
#' @references
#' Beckmann, J., Koop, G., Korobilis, D., and Schüssler, R. A. (2020)
#' "Exchange rate predictability and dynamic bayesian learning."
#' \emph{Journal of Applied Econometrics}, 35 (4): 410–421.
#'
#' Dangl, T. and Halling, M. (2012)
#' "Predictive regressions with time-varying coefficients."
#' \emph{Journal of Financial Economics}, 106 (1): 157–181.
#'
#' Del Negro, M., Hasegawa, R. B., and Schorfheide, F. (2016)
#' "Dynamic prediction pools:
#' An investigation of financial frictions and forecasting performance."
#' \emph{Journal of Econometrics}, 192 (2): 391–405.
#'
#' Koop, G. and Korobilis, D. (2012)
#' "Forecasting inflation using dynamic model averaging."
#' \emph{International Economic Review}, 53 (3): 867–886.
#'
#' Koop, G. and Korobilis, D. (2023)
#' "Bayesian dynamic variable selection in high dimensions."
#' \emph{International Economic Review}.
#'
#' Raftery, A. E., Kárn`y, M., and Ettler, P. (2010)
#' "Online prediction under model uncertainty via dynamic model averaging:
#' Application to a cold rolling mill."
#' \emph{Technometrics}, 52 (1): 52–66.
#'
#' West, M. and Harrison, J. (1997)
#' "Bayesian forecasting and dynamic models"
#' \emph{Springer}, 2nd edn.
#' 
#' @export
#' @import checkmate
#' @importFrom stats na.omit
#' @examples
#' \donttest{
#'
#'    #########################################################
#'    ######### Forecasting quarterly U.S. inflation ##########
#'    #### Please see Koop & Korobilis (2023) for further  ####
#'    #### details regarding the data & external forecasts ####
#'    #########################################################
#'
#'    # Load Package
#'    library("hdflex")
#'    library("ggplot2")
#'    library("cowplot")
#'
#'    ########## Get Data ##########
#'    # Load Package Data
#'    inflation_data <- inflation_data
#'
#'    # Set Target Variable
#'    y <- inflation_data[,  1]
#'
#'    # Set 'P-Signals'
#'    X <- inflation_data[, 2:442]
#'
#'    # Set 'F-Signals'
#'    Ext_F <- inflation_data[, 443:462]
#'
#'    # Get Dates and Number of Observations
#'    tdates <- rownames(inflation_data)
#'    tlength <- length(tdates)
#'
#'    # First complete observation (no missing values)
#'    first_complete <- which(complete.cases(inflation_data))[1]
#'
#'    ########## Rolling AR2-Benchmark ##########
#'    # Set up matrix for predictions
#'    benchmark <- matrix(NA, nrow = tlength,
#'                        ncol = 1, dimnames = list(tdates, "AR2"))
#'
#'    # Set Window-Size (15 years of quarterly data)
#'    window_size <- 15 * 4
#'
#'    # Time Sequence
#'    t_seq <- seq(window_size, tlength - 1)
#'
#'    # Loop with rolling window
#'    for (t in t_seq) {
#'
#'      # Split Data for Training Train Data
#'      x_train <- cbind(int = 1, X[(t - window_size + 1):t, 1:2])
#'      y_train <- y[(t - window_size + 1):t]
#'
#'      # Split Data for Prediction
#'      x_pred <- cbind(int = 1, X[t + 1, 1:2, drop = FALSE])
#'
#'      # Fit AR-Model
#'      model_ar <- .lm.fit(x_train, y_train)
#'
#'      # Predict and store in benchmark matrix
#'      benchmark[t + 1, ] <- x_pred %*% model_ar$coefficients
#'    }
#'
#'    ########## STSC ##########
#'    ### Part 1: TVC-Function
#'    # Set TV-C-Parameter
#'    init <- 5 * 4
#'    lambda_grid <- c(0.90, 0.95, 1.00)
#'    kappa_grid <- c(0.94, 0.96, 0.98)
#'    bias <- TRUE
#'
#'    # Apply TVC-Function
#'    tvc_results <- hdflex::tvc(y,
#'                               X,
#'                               Ext_F,
#'                               init,
#'                               lambda_grid,
#'                               kappa_grid,
#'                               bias)
#'
#'    # Assign TVC-Results
#'    forecast_tvc <- tvc_results$Forecasts$Point_Forecasts
#'    variance_tvc <- tvc_results$Forecasts$Variance_Forecasts
#'
#'    # First complete forecast period (no missing values)
#'    sub_period <- seq(which(complete.cases(forecast_tvc))[1], tlength)
#'
#'    ### Part 2: DSC-Function
#'    # Set DSC-Parameter
#'    gamma_grid <- c(0.40, 0.50, 0.60, 0.70, 0.80, 0.90,
#'                    0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00)
#'    psi_grid <- c(1:100, sapply(1:4, function(i) floor(i * ncol(forecast_tvc) / 4)))
#'    delta <- 0.95
#'    burn_in_tvc <- (init / 2) + 1
#'    burn_in_dsc <- 1
#'    metric <- 5
#'    equal_weight <- TRUE
#'    incl <- NULL
#'
#'    # Apply DSC-Function
#'    dsc_results <- hdflex::dsc(y[sub_period],
#'                               forecast_tvc[sub_period, , drop = FALSE],
#'                               variance_tvc[sub_period, , drop = FALSE],
#'                               gamma_grid,
#'                               psi_grid,
#'                               delta,
#'                               burn_in_tvc,
#'                               burn_in_dsc,
#'                               metric,
#'                               equal_weight,
#'                               incl,
#'                               NULL)
#'
#'    # Assign DSC-Results
#'    pred_stsc <- dsc_results$Forecasts$Point_Forecasts
#'    var_stsc <- dsc_results$Forecasts$Variance_Forecasts
#'
#'    ########## Evaluation ##########
#'    # Define Evaluation Period (OOS-Period)
#'    eval_period <- which(tdates[sub_period] >= "1991-04-01" & tdates[sub_period] <= "2021-12-01")
#'
#'    # Get Evaluation Summary for STSC
#'    eval_results <- summary(obj = dsc_results, eval_period = eval_period)
#'
#'    # Calculate (Mean-)Squared-Errors for AR2-Benchmark
#'    oos_y <- y[sub_period][eval_period]
#'    oos_benchmark <- benchmark[sub_period[eval_period], , drop = FALSE]
#'    se_ar2 <- (oos_y - oos_benchmark)^2
#'    mse_ar2 <- mean(se_ar2)
#'
#'    # Create Cumulative Squared Error Differences (CSSED) Plot
#'    cssed <- cumsum(se_ar2 - eval_results$MSE[[2]])
#'    plot_cssed <- ggplot(
#'      data.frame(eval_period, cssed),
#'      aes(x = eval_period, y = cssed)
#'    ) +
#'      geom_line() +
#'      ylim(-0.0008, 0.0008) +
#'      ggtitle("Cumulative Squared Error Differences") +
#'      xlab("Time Index") +
#'      ylab("CSSED") +
#'      geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
#'      theme_minimal(base_size = 15) +
#'      theme(
#'        panel.grid.major = element_blank(),
#'        panel.grid.minor = element_blank(),
#'        panel.border = element_rect(colour = "black", fill = NA),
#'        axis.ticks = element_line(colour = "black"),
#'        plot.title = element_text(hjust = 0.5)
#'      )
#'
#'    # Show Plots
#'    options(repr.plot.width = 15, repr.plot.height = 15)
#'    plots_list <- eval_results$Plots
#'    plots_list <- c(list(plot_cssed), plots_list)
#'    cowplot::plot_grid(plotlist = plots_list,
#'                       ncol = 2,
#'                       nrow = 3,
#'                       align = "hv")
#'
#'    # Relative MSE
#'    print(paste("Relative MSE:", round(eval_results$MSE[[1]] / mse_ar2, 4)))
#'  }

### Time-Varying Coefficient Model
tvc  <- function(y,
                 X,
                 Ext_F,
                 init,
                 lambda_grid,
                 kappa_grid,
                 bias) {

  # Convert y to matrix
  y <- as.matrix(y, ncol = 1)

  ########################################################
  ### Checkmate
  # Check if numeric vector without missing or infinite values
  assertNumeric(y, min.len = 1, any.missing = FALSE, finite = TRUE)

  # Either X or Ext_F must not be null
  assert(checkMatrix(X), checkMatrix(Ext_F), combine = "or")

  # Check if numeric matrices and have same length as y
  assertMatrix(X, mode = "numeric", nrows = nrow(y), null.ok = TRUE)
  assertMatrix(Ext_F, mode = "numeric", nrows = nrow(y), null.ok = TRUE)

  # Check if integer between 2 and the length of y
  assertInt(init, lower = 2, upper = nrow(y))

  # Check if numeric vectors with values between exp(-10) and 1
  assertNumeric(lambda_grid, lower = exp(-10), upper = 1, min.len = 1, any.missing = FALSE, finite = TRUE)
  assertNumeric(kappa_grid, lower = exp(-10), upper = 1, min.len = 1, any.missing = FALSE, finite = TRUE)

  # Check if boolean
  assertLogical(bias, len = 1, any.missing = FALSE)

  # Check if there are any non-consecutive NA values
  non_consec_na <- any(
    apply(
      is.na(cbind(if (exists("X")) X, if (exists("Ext_F")) Ext_F)), 2,
      function(x) sum(abs(diff(x))) > 1 | sum(diff(x)) == 1
    )
  )
  assertFALSE(non_consec_na)

  # Check if any column in X is constant for the first observations
  if (!is.null(X)) {
    if (any(apply(X, 2, function(x) length(unique(na.omit(x)[1:init])) == 1))) {
      print("One or more columns in X are constant for the first 1:init observations.")
    }
  }

  # Check if any column in Ext_F is constant for the first observations
  if (!is.null(Ext_F)) {
    if (any(apply(Ext_F, 2, function(x) length(unique(na.omit(x)[1:init])) == 1))) {
      print("One or more columns in Ext_F are constant for the first 1:init observations.")
    }
  }

  ########################################################
  ### Apply Rcpp-Function
  tvc_results <- tvc_(y,
                      X,
                      Ext_F,
                      init,
                      lambda_grid,
                      kappa_grid,
                      bias)

  ### Assign Results
  forecast_tvc <- tvc_results[[1]]
  variance_tvc <- tvc_results[[2]]

  ### Remove
  rm(list = c("tvc_results"))

  ### Create / Get P-Signal Names
  x_names <- if (!is.null(X)) {
    if (!is.null(colnames(X))) {
      colnames(X)
    } else {
      paste0("X", as.character(seq_len(ncol(X))))
    }
  }

  if (!is.null(X)) {

    # Preallocate TVC-Model Names
    tvc_x_name <- vector("character", length(lambda_grid) * length(kappa_grid) * ncol(X))

    # Create TVC-Model Names
    i <- 1
    for (l in lambda_grid) {
      for (k in kappa_grid) {
        for (j in seq_len(ncol(X))) {
          # Create TVC-Model Name
          tvc_x_name[i] <- paste(x_names[j], l, k, sep = "_")
          i <- i + 1
        }
      }
    }
  }

  ### Create / Get F-Signal Names
  f_names <- if (!is.null(Ext_F)) {
    if (!is.null(colnames(Ext_F))) {
      colnames(Ext_F)
    } else {
      paste0("Ext_F", as.character(seq_len(ncol(Ext_F))))
    }
  }

  if (!is.null(Ext_F)) {

    # Preallocate TVC-Model Names
    tvc_f_name <- vector("character", length(lambda_grid) * length(kappa_grid) * ncol(Ext_F))

    # Create TVC-Model Names
    i <- 1
    for (l in lambda_grid) {
      for (k in kappa_grid) {
        for (j in seq_len(ncol(Ext_F))) {
          # Create TVC-Model Name
          tvc_f_name[i] <- paste(f_names[j], l, k, sep = "_")
          i <- i + 1
        }
      }
    }
  }

  # Combine Signal Names
  model_names_tvc <- c(if (exists("tvc_x_name")) tvc_x_name,
                       if (exists("tvc_f_name")) tvc_f_name)

  # Assign Model Names (-> Column Names)
  colnames(forecast_tvc) <- model_names_tvc
  colnames(variance_tvc) <- model_names_tvc

  # Return Results
  return(
    list(
      Forecasts = list(
        Realization = y,
        Point_Forecasts = forecast_tvc,
        Variance_Forecasts = variance_tvc
      ),
      Model = list(
        Lambda_grid = lambda_grid,
        Kappa_grid = kappa_grid,
        Init = init,
        Bias = bias
      )
    )
  )
}