#' @title Dynamic Nelson-Siegel Model for Interest Rate Modeling
#' @description A comprehensive toolkit for modeling, fitting, and simulating interest rates
#' using the Dynamic Nelson-Siegel (DNS) model. The package provides functions to calibrate
#' the DNS model from yield curve data, estimate factor dynamics, simulate yield curves,
#' and visualize results.
#'
#' @details
#' The Dynamic Nelson-Siegel model is a popular term structure model used in finance for
#' modeling the yield curve. It represents the yield curve using three latent factors (level,
#' slope, and curvature) that evolve over time according to a mean-reverting stochastic process.
#'
#' This package provides a comprehensive implementation of the DNS model, including:
#'
#' \itemize{
#'   \item Data retrieval and processing of U.S. Treasury yield curve data
#'   \item Cross-sectional fitting of the Nelson-Siegel model to yield curves
#'   \item Time-series estimation of factor dynamics
#'   \item Stochastic simulation of interest rate paths
#'   \item Visualization tools for yield curves and factor analysis
#' }
#'
#' The package is designed to be used by financial analysts, risk managers, and researchers
#' for applications such as:
#'
#' \itemize{
#'   \item Interest rate forecasting and scenario analysis
#'   \item Risk management and stress testing
#'   \item Pricing and valuation of fixed income securities
#'   \item Asset-liability management
#' }
#'
#' @section Key Functions:
#' \describe{
#'   \item{\code{\link{get_treasury_yield_data}}}{Download and process Treasury yield curve data}
#'   \item{\code{\link{get_treasury_yield_combined}}}{Get historical and current Treasury yield data}
#'   \item{\code{\link{fit_ns_model}}}{Fit Nelson-Siegel model to a single yield curve}
#'   \item{\code{\link{fit_dns_model}}}{Fit Dynamic Nelson-Siegel model to time series of yield curves}
#'   \item{\code{\link{estimate_dns_dynamics}}}{Estimate dynamics of the DNS factors}
#'   \item{\code{\link{get_calibrated_dns_params}}}{Get calibrated DNS parameters from market data}
#'   \item{\code{\link{simulate_dns_factors}}}{Simulate DNS factors using YUIMA package}
#'   \item{\code{\link{calculate_dns_yields}}}{Calculate yield curves from DNS factors}
#'   \item{\code{\link{simulate_dns_yield_curves}}}{Simulate DNS yield curves}
#'   \item{\code{\link{plot_yield_curves}}}{Plot yield curves and fitted DNS models}
#'   \item{\code{\link{plot_factor_time_series}}}{Plot DNS factor time series}
#' }
#'
#' @section Model Background:
#' The Nelson-Siegel model expresses the yield curve as:
#' \deqn{y(m) = \beta_1 + \beta_2 \left(\frac{1-e^{-\lambda m}}{\lambda m}\right) + \beta_3 \left(\frac{1-e^{-\lambda m}}{\lambda m} - e^{-\lambda m}\right)}
#' where \eqn{m} is maturity, \eqn{\beta_1}, \eqn{\beta_2}, and \eqn{\beta_3} represent the level, slope, and curvature
#' factors, respectively, and \eqn{\lambda} controls the decay rate of factor loadings.
#'
#' In the Dynamic Nelson-Siegel model, these factors follow a mean-reverting stochastic process:
#' \deqn{d\beta(t) = \kappa(\theta - \beta(t))dt + \sigma dW(t)}
#' where \eqn{\kappa} is the speed of mean reversion, \eqn{\theta} is the long-run mean,
#' \eqn{\sigma} is the volatility, and \eqn{W(t)} is a Wiener process.
#'
#' @examples
#' \dontrun{
#' # Get Treasury yield data
#' yield_data <- get_treasury_yield_combined("2020-01-01")
#'
#' # Fit DNS model
#' dns_result <- fit_dns_model(yield_data, n_days = 100)
#'
#' # Estimate factor dynamics
#' dynamics <- estimate_dns_dynamics(dns_result$parameters)
#'
#' # Simulate future yield curves, starting from the most recent fitted factors
#' latest <- dns_result$parameters[which.max(dns_result$parameters$Date), ]
#' params <- dns_dynamics_to_sim_params(dynamics, lambda = latest$lambda)
#' sim_result <- simulate_dns_yield_curves(
#'   params = params,
#'   initial_factors = as.numeric(latest[c("beta1", "beta2", "beta3")]),
#'   maturities = c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30),
#'   n_paths = 100,
#'   engine = "exact"
#' )
#'
#' # Visualize results
#' plot_yield_curves(yield_data, dns_result)
#' plot_factor_time_series(dns_result)
#' plot_simulated_yields(sim_result)
#' }
#'
#' @name dnsr-package
"_PACKAGE"

# Imports shared across the package. These are declared here rather than in
# each file so there is one place to look; per-file blocks cover what is local
# to a file.
#' @importFrom stats coef cor lm median na.omit optim prcomp predict quantile
#' @importFrom stats resid residuals rnorm sd smooth.spline
#' @importFrom utils head tail
#' @importFrom methods is
NULL

# Columns created inside the function (by reshape2::melt, expand.grid or
# data.frame) and then referenced in ggplot2 aesthetics. R CMD check cannot see
# that they are local, so it reports them as undefined globals.
utils::globalVariables(c(
  # yield data and fitted factors
  "Date", "rmse", "maturity", "yield", "type", "label",
  # melted long-format frames
  "Factor", "Value", "Change", "Loading", "Label",
  "Maturity", "Yield", "Model",
  # simulation frames
  "Time", "Type", "Path", "Mean", "Q10", "Q90", "Forward",
  # PCA frames
  "Component", "Weight", "Score", "Variance", "Cumulative"
))
