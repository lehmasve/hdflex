#########################  STSC - S3 Object #########################
#' @name summary.stsc_obj
#' @title Summary for 'stsc' object
#' @description This function plots the evolution of the tuning parameters for an 'stsc' object and returns basic performance metrics.
#' @param object An object of type 'stsc'.
#' @param eval_period (Optional) A vector of indices to specify the evaluation period. Defaults to the entire period after burn-in.
#' @param ... Additional arguments to be consistent with the S3 print() function.
#' @method summary stsc_obj
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom stats complete.cases dnorm pnorm
#' @references
#' Gneiting, T., Raftery, A. E., Westveld, A. H., and Goldman, T. (2005):
#' Calibrated Probabilistic Forecasting Using Ensemble Model Output Statistics and Minimum CRPS Estimation.
#' \emph{Monthly Weather Review}, 133: 1098–1118.
#'
#' Jordan, A., Krueger, F., and Lerch, S. (2019):
#' "Evaluating Probabilistic Forecasts with scoringRules."
#' \emph{Journal of Statistical Software}, 90(12): 1-37.
#'
#' @return A list containing:
#' \describe{
#'   \item{MSE}{A list with the mean squared error (MSE) and squared errors (SE).}
#'   \item{ACRPS}{A list with the average continuous ranked probability score (ACRPS) and CRPS values.}
#'   \item{APLL}{A list with the average predictive log-likelihood (APLL) and predictive log-likelihood (PLL) values.}
#'   \item{Plots}{A list of ggplot objects for visualizing the tuning parameters and selected signals.}
#' }
#' @examples
#' \donttest{
#'
#' # See example for stsc().
#'
#' }
#' @export

summary.stsc_obj <- function(object, eval_period = NULL, ...) {

  ### Data
  # Extract realized values from object
  y <- object$Forecasts$Realization
  point_forecast <- object$Forecasts$Point_Forecasts
  variance_forecast <- object$Forecasts$Variance_Forecasts
  chosen_gamma <- object$Tuning_Parameters$Gamma
  chosen_psi <- object$Tuning_Parameters$Psi
  chosen_lambda <- object$Tuning_Parameters$Lambda
  chosen_kappa <- object$Tuning_Parameters$Kappa
  chosen_signals <- object$Tuning_Parameters$Signals
  gamma_grid <- object$Model$Gamma_grid
  psi_grid <- object$Model$Psi_grid
  init <- object$Model$Init
  burn_in <- object$Model$Burn_in
  burn_in_dsc <- object$Model$Burn_in_dsc

  ### Evaluation
  # Set Evaluation Period
  if (is.null(eval_period)) {
    start <- max(init, burn_in, burn_in_dsc) + 1
    eval_period <- start:nrow(y)
  }

  # Check for NaNs in point_forecast within eval_period
  if (any(is.na(point_forecast[eval_period]))) {
    stop("Invalid 'eval_period': results contain NaNs. Please adjust 'eval_period'.")
  }

  # Cut Objects to evaluation period
  y <- y[eval_period, ]
  point_forecast <- point_forecast[eval_period]
  variance_forecast <- variance_forecast[eval_period]
  chosen_gamma <- chosen_gamma[eval_period]
  chosen_psi <- chosen_psi[eval_period]
  chosen_lambda <- chosen_lambda[eval_period, , drop = FALSE]
  chosen_kappa <- chosen_kappa[eval_period, , drop = FALSE]
  chosen_signals <- chosen_signals[eval_period, , drop = FALSE]

  # Calculate MSE / SE
  se <- (y - point_forecast)^2
  mse <- mean(se)

  # Calculate ACRPS / CRPS
  z <- (y - point_forecast) / sqrt(variance_forecast)
  pdf <- dnorm(z, 0.0, 1.0)
  cdf <- pnorm(z, 0.0, 1.0)
  pi_inv <- 1.0 / sqrt(pi)
  crps <- sqrt(variance_forecast) * (z * (2.0 * cdf - 1.0) + 2.0 * pdf - pi_inv)
  acrps <- mean(crps)

  # (Mean) Predictive Log-Likelihood
  pll <- pnorm(y, point_forecast, sqrt(variance_forecast), log.p = TRUE)
  apll <- mean(pll)

  ### Visualization
  # Function to create a factor plot
  plot_factor <- function(x, y, grid, main, xlab, ylab) {

    # Convert to factor for plotting
    factor_y <- factor(y, levels = grid)
    data <- data.frame(Time = x, Value = factor_y)

    # Ticks / Breaks for y-axis
    max_ticks <- 30
    breaks <- if (length(levels(factor_y)) > max_ticks) {
      pretty(seq_along(levels(factor_y)), n = max_ticks)
    } else {
      seq_along(levels(factor_y))
    }
    breaks <- levels(factor_y)[breaks]

    # Create the ggplot
    ggplot(data, aes(x = Time, y = Value)) +
      geom_point() +
      scale_y_discrete(drop = FALSE, breaks = breaks) +
      labs(title = main, x = xlab, y = ylab) +
      theme_minimal(base_size = 15) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.ticks = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5)
      )
  }

  # Function to create an area plot
  plot_area <- function(x, m, main, xlab, ylab, legend_title) {

    # Normalize matrix (adding up to 100%)
    m_normalized <- sweep(m, 1, rowSums(m), FUN = "/")

    # Convert to data frame for ggplot
    data <- as.data.frame(m_normalized)
    data$Time <- x
    data_long <- reshape2::melt(data,
      id.vars = "Time",
      variable.name = "Variable",
      value.name = "Value"
    )

    # Create the ggplot
    ggplot(data_long, aes(x = Time, y = Value, fill = Variable)) +
      geom_area(position = "fill") +
      labs(title = main, x = xlab, y = ylab, fill = legend_title) +
      scale_fill_grey(start = 0.3, end = 0.9) +
      theme_minimal(base_size = 15) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.ticks = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.box = "horizontal"
      )
  }

  # Function to create a signal plot
  plot_signal <- function(x, signals, main, xlab, ylab) {

    # Dimension restriction
    if (ncol(signals) > 5000) {
      # Select top 5000 signals by non-zero counts
      non_zero_counts <- colSums(signals != 0)
      top_columns <- order(non_zero_counts, decreasing = TRUE)[1:5000]
      signals <- signals[, sort(top_columns), drop = FALSE]
    }

    # Prepare data for ggplot
    mat <- signals %*% diag(seq_len(ncol(signals)))
    mat[mat == 0] <- NA
    data <- as.data.frame(mat)
    data$Time <- x
    data_long <- reshape2::melt(data,
      id.vars = "Time",
      variable.name = "Variable",
      value.name = "Value"
    )

    # Create the ggplot
    ggplot(data_long, aes(x = Time, y = Value)) +
      geom_point(size = 0.5, na.rm = TRUE) +
      labs(title = main, x = xlab, y = ylab) +
      expand_limits(y = c(1, ncol(mat))) +
      theme_minimal(base_size = 15) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.ticks = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5)
      )
  }

  # List to store ggplot objects
  plots_list <- list()

  ### Plot: Tuning Parameter Gamma
  plots_list$Gamma <- plot_factor(eval_period,
                                  chosen_gamma,
                                  gamma_grid,
                                  expression("Discount Factor" ~ gamma),
                                  "Time Index",
                                  expression(gamma))

  ### Plot: Tuning Parameter Psi
  plots_list$Psi <- plot_factor(eval_period,
                                chosen_psi,
                                psi_grid,
                                expression("Subset-Size" ~ psi),
                                "Time Index",
                                expression(psi))

  ### Plot: Evolution of selected signals
  plots_list$Signals <- plot_signal(eval_period,
                                    chosen_signals,
                                    "Selected Signals",
                                    "Time Index",
                                    "Predictive Signal")

  ### Plot: Tuning Parameter Lambda
  plots_list$Lambda <- plot_area(eval_period,
                                 chosen_lambda,
                                 expression("Discount Factor" ~ lambda),
                                 "Time Index",
                                 expression(lambda),
                                 "Values")

  ### Plot: Tuning Parameter Kappa
  plots_list$Kappa <- plot_area(eval_period,
                                chosen_kappa,
                                expression("Discount Factor" ~ kappa),
                                "Time Index",
                                expression(kappa),
                                "Values")

  # Return
  return(
    list(
      MSE = list(mse, se),
      ACRPS = list(acrps, crps),
      APLL = list(apll, pll),
      Plots = plots_list
    )
  )
}

#########################  DSC - S3 Object #########################
#' @name summary.dsc_obj
#' @title Summary for 'dsc' object
#' @description This function plots the evolution of the tuning parameters for a 'dsc' object and returns basic performance metrics.
#' @param object An object of type 'dsc'.
#' @param eval_period (Optional) A vector of indices to specify the evaluation period. Defaults to the entire period after burn-in.
#' @param ... Additional arguments to be consistent with the S3 print() function.
#' @method summary dsc_obj
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom stats complete.cases dnorm na.omit pnorm
#' @references
#' Gneiting, T., Raftery, A. E., Westveld, A. H., and Goldman, T. (2005):
#' Calibrated Probabilistic Forecasting Using Ensemble Model Output Statistics and Minimum CRPS Estimation.
#' \emph{Monthly Weather Review}, 133: 1098–1118.
#'
#' Jordan, A., Krueger, F., and Lerch, S. (2019):
#' "Evaluating Probabilistic Forecasts with scoringRules."
#' \emph{Journal of Statistical Software}, 90(12): 1-37.
#'
#' @return A list containing:
#' \describe{
#'   \item{MSE}{A list with the mean squared error (MSE) and squared errors (SE).}
#'   \item{ACRPS}{A list with the average continuous ranked probability score (ACRPS) and CRPS values.}
#'   \item{APLL}{A list with the average predictive log-likelihood (APLL) and predictive log-likelihood (PLL) values.}
#'   \item{Plots}{A list of ggplot objects for visualizing the tuning parameters and selected CFMs.}
#' }
#' @examples
#' \donttest{
#'
#' # See example for tvc().
#'
#' }
#' @export

summary.dsc_obj <- function(object, eval_period = NULL, ...) {

  ### Data
  # Extract realized values from object
  y <- object$Forecasts$Realization
  point_forecast <- object$Forecasts$Point_Forecasts
  variance_forecast <- object$Forecasts$Variance_Forecasts
  gamma_grid <- object$Model$Gamma_grid
  psi_grid <- object$Model$Psi_grid
  chosen_gamma <- object$Tuning_Parameters$Gamma
  chosen_psi <- object$Tuning_Parameters$Psi
  chosen_cfms <- object$Tuning_Parameters$CFM
  burn_in <- object$Model$Burn_in
  burn_in_dsc <- object$Model$Burn_in_dsc

  ### Evaluation
  # Set Evaluation Period
  if (is.null(eval_period)) {
    start <- max(burn_in, burn_in_dsc) + 1
    eval_period <- start:nrow(y)
  }

  # Check for NaNs in point_forecast within eval_period
  if (any(is.na(point_forecast[eval_period]))) {
    stop("Invalid 'eval_period': results contain NaNs. Please adjust 'eval_period'.")
  }

  # Cut Objects to evaluation period
  y <- y[eval_period, ]
  point_forecast <- point_forecast[eval_period]
  variance_forecast <- variance_forecast[eval_period]
  chosen_gamma <- chosen_gamma[eval_period]
  chosen_psi <- chosen_psi[eval_period]
  chosen_cfms <- chosen_cfms[eval_period, , drop = FALSE]

  # Calculate MSE / SE
  se <- (y - point_forecast)^2
  mse <- mean(se)

  # Calculate ACRPS / CRPS
  z <- (y - point_forecast) / sqrt(variance_forecast)
  pdf <- dnorm(z, 0.0, 1.0)
  cdf <- pnorm(z, 0.0, 1.0)
  pi_inv <- 1.0 / sqrt(pi)
  crps <- sqrt(variance_forecast) * (z * (2.0 * cdf - 1.0) + 2.0 * pdf - pi_inv)
  acrps <- mean(crps)

  # (Mean) Predictive Log-Likelihood
  pll <- pnorm(y, point_forecast, sqrt(variance_forecast), log.p = TRUE)
  apll <- mean(pll)

  ### Visualization
  # Function to create a factor plot
  plot_factor <- function(x, y, grid, main, xlab, ylab) {

    # Convert to factor for plotting
    factor_y <- factor(y, levels = grid)
    data <- data.frame(Time = x, Value = factor_y)

    # Ticks / Breaks for y-axis
    max_ticks <- 30
    breaks <- if (length(levels(factor_y)) > max_ticks) {
      pretty(seq_along(levels(factor_y)), n = max_ticks)
    } else {
      seq_along(levels(factor_y))
    }
    breaks <- levels(factor_y)[breaks]

    # Create the ggplot
    ggplot(data, aes(x = Time, y = Value)) +
      geom_point() +
      scale_y_discrete(drop = FALSE, breaks = breaks) +
      labs(title = main, x = xlab, y = ylab) +
      theme_minimal(base_size = 15) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.ticks = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5)
      )
  }

  # Function to create a cfm plot
  plot_cfm <- function(x, cfms, main, xlab, ylab) {

    # Dimension restriction
    if (ncol(cfms) > 5000) {
      # Select top 5000 CFMs by non-zero counts
      non_zero_counts <- colSums(cfms != 0)
      top_columns <- order(non_zero_counts, decreasing = TRUE)[1:5000]
      cfms <- cfms[, sort(top_columns), drop = FALSE]
    }

    # Prepare data for ggplot
    mat <- cfms %*% diag(seq_len(ncol(cfms)))
    mat[mat == 0] <- NA
    data <- as.data.frame(mat)
    data$Time <- x
    data_long <- reshape2::melt(data,
      id.vars = "Time",
      variable.name = "Variable",
      value.name = "Value"
    )

    # Create the ggplot
    ggplot(data_long, aes(x = Time, y = Value)) +
      geom_point(size = 0.5, na.rm = TRUE) +
      labs(title = main, x = xlab, y = ylab) +
      expand_limits(y = c(1, ncol(mat))) +
      theme_minimal(base_size = 15) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.ticks = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5)
      )
  }

  # List to store ggplot objects
  plots_list <- list()

  ### Plot: Tuning Parameter Gamma
  plots_list$Gamma <- plot_factor(eval_period,
                                  chosen_gamma,
                                  gamma_grid,
                                  expression("Discount Factor" ~ gamma),
                                  "Time Index",
                                  expression(gamma))

  ### Plot: Tuning Parameter Psi
  plots_list$Psi <- plot_factor(eval_period,
                                chosen_psi,
                                psi_grid,
                                expression("Subset-Size" ~ psi),
                                "Time Index",
                                expression(psi))

  ### Plot: Evolution of selected CFMs
  plots_list$CFM <- plot_cfm(eval_period,
                             chosen_cfms,
                             "Selected CFMs",
                             "Time Index",
                             "Predictive Signal")

  # Return
  return(
    list(
      MSE = list(mse, se),
      ACRPS = list(acrps, crps),
      APLL = list(apll, pll),
      Plots = plots_list
    )
  )
}
