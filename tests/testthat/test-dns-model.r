library(testthat)
library(dnsr)

test_that("fit_ns_model handles simple yield curve correctly", {
  # Create a simple yield curve
  maturities <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)
  yields <- c(0.01, 0.012, 0.015, 0.018, 0.02, 0.025, 0.027, 0.03, 0.032, 0.033)
  
  # Fit Nelson-Siegel model
  result <- fit_ns_model(maturities, yields)
  
  # Check result structure
  expect_type(result, "list")
  expect_named(result, c("parameters", "fitted_yields", "rmse", "maturities", 
                        "actual_yields", "residuals", "success"))
  
  # Check parameters
  expect_length(result$parameters, 4)
  expect_named(result$parameters, c("beta1", "beta2", "beta3", "lambda"))
  
  # Check that all parameters are numeric and finite
  expect_true(all(is.finite(result$parameters)))
  
  # Check fitted values
  expect_length(result$fitted_yields, length(maturities))
  expect_true(all(is.finite(result$fitted_yields)))
  
  # Check that RMSE is reasonable
  expect_lt(result$rmse, 0.01)  # RMSE should be less than 1%
  
  # Check success flag
  expect_true(result$success)
})

test_that("fit_ns_model handles missing values correctly", {
  # Create a yield curve with missing values
  maturities <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)
  yields <- c(0.01, NA, 0.015, 0.018, NA, 0.025, 0.027, 0.03, NA, 0.033)
  
  # Fit Nelson-Siegel model
  result <- fit_ns_model(maturities, yields)
  
  # Check result structure
  expect_type(result, "list")
  expect_named(result, c("parameters", "fitted_yields", "rmse", "maturities", 
                        "actual_yields", "residuals", "success"))
  
  # Check that maturities with missing yields are excluded from fitted values
  valid_idx <- !is.na(yields)
  expect_length(result$maturities, sum(valid_idx))
  expect_length(result$fitted_yields, sum(valid_idx))
  
  # Check that filtered maturities match expected
  expect_equal(result$maturities, maturities[valid_idx])
  
  # Check that RMSE is reasonable
  expect_lt(result$rmse, 0.01)  # RMSE should be less than 1%
})

test_that("fit_ns_model handles fixed lambda correctly", {
  # Create a simple yield curve
  maturities <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)
  yields <- c(0.01, 0.012, 0.015, 0.018, 0.02, 0.025, 0.027, 0.03, 0.032, 0.033)
  
  # Fixed lambda value
  fixed_lambda <- 0.1
  
  # Fit Nelson-Siegel model with fixed lambda
  result <- fit_ns_model(maturities, yields, initial_lambda = fixed_lambda, 
                        lambda_constraint = "fixed")
  
  # Check that lambda equals the fixed value
  expect_equal(result$parameters["lambda"], fixed_lambda)
})

test_that("fit_ns_model handles extreme yield curves gracefully", {
  # Create an extreme, but realistic yield curve (inverted)
  maturities <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)
  yields <- c(0.05, 0.049, 0.048, 0.045, 0.043, 0.04, 0.038, 0.035, 0.033, 0.032)
  
  # Fit Nelson-Siegel model
  result <- fit_ns_model(maturities, yields)
  
  # Check result structure
  expect_type(result, "list")
  expect_named(result, c("parameters", "fitted_yields", "rmse", "maturities", 
                        "actual_yields", "residuals", "success"))
  
  # Check parameters - beta2 should be negative for inverted curve
  expect_lt(result$parameters["beta2"], 0)
  
  # Create another extreme case - humped yield curve
  yields_humped <- c(0.01, 0.015, 0.02, 0.025, 0.028, 0.026, 0.024, 0.023, 0.022, 0.021)
  
  # Fit Nelson-Siegel model
  result_humped <- fit_ns_model(maturities, yields_humped)
  
  # Check parameters - beta3 should be non-zero for humped curve
  expect_true(abs(result_humped$parameters["beta3"]) > 0.005)
  
  # Check success flag
  expect_true(result_humped$success)
})

test_that("fit_dns_model handles time series data correctly", {
  # Skip if we don't have an internet connection for downloading data
  skip_if_offline()
  skip_on_cran()
  
  # Create a very simple time series dataset
  dates <- seq(as.Date("2023-01-01"), as.Date("2023-01-10"), by = "day")
  maturities <- c(0.25, 0.5, 1, 2, 5, 10, 30)
  
  # Create a simple yield data frame
  yield_data <- data.frame(Date = dates)
  for (m in maturities) {
    col_name <- as.character(m)
    # Create synthetic yields with some time variation
    yield_data[[col_name]] <- 0.02 + 0.01 * sin(m / 5) + 0.002 * as.numeric(dates - dates[1]) / 10
  }
  
  # Fit DNS model
  dns_result <- fit_dns_model(yield_data, parallel = FALSE, show_progress = FALSE)
  
  # Check result structure
  expect_type(dns_result, "list")
  expect_named(dns_result, c("parameters", "fitted_yields", "maturities", "success_rate"))
  
  # Check parameters data frame
  expect_s3_class(dns_result$parameters, "data.frame")
  expect_true("Date" %in% names(dns_result$parameters))
  expect_true(all(c("beta1", "beta2", "beta3", "lambda", "rmse") %in% names(dns_result$parameters)))
  
  # Check dimensions of fitted yields
  expect_equal(nrow(dns_result$fitted_yields), nrow(dns_result$parameters))
  expect_equal(ncol(dns_result$fitted_yields), length(maturities))
  
  # Check success rate
  expect_gte(dns_result$success_rate, 0)
  expect_lte(dns_result$success_rate, 100)
})

test_that("estimate_dns_dynamics calculates correct parameters", {
  # Create a simulated factor time series
  n_days <- 100
  dates <- seq(as.Date("2023-01-01"), by = "day", length.out = n_days)
  
  # Create factor parameters for simulation
  kappa1 <- 0.02  # Slow mean reversion for level
  kappa2 <- 0.3   # Medium mean reversion for slope
  kappa3 <- 0.5   # Fast mean reversion for curvature
  
  theta1 <- 0.04  # Long-run mean for level
  theta2 <- -0.02 # Long-run mean for slope
  theta3 <- 0.01  # Long-run mean for curvature
  
  sigma1 <- 0.005 # Volatility for level
  sigma2 <- 0.01  # Volatility for slope
  sigma3 <- 0.02  # Volatility for curvature
  
  # Simulate factors
  dt <- 1/252
  beta1 <- numeric(n_days)
  beta2 <- numeric(n_days)
  beta3 <- numeric(n_days)
  
  # Initial values
  beta1[1] <- 0.035
  beta2[1] <- -0.015
  beta3[1] <- 0.005
  
  # Simulate Ornstein-Uhlenbeck processes
  set.seed(123)
  for (i in 2:n_days) {
    beta1[i] <- beta1[i-1] + kappa1 * (theta1 - beta1[i-1]) * dt + sigma1 * sqrt(dt) * rnorm(1)
    beta2[i] <- beta2[i-1] + kappa2 * (theta2 - beta2[i-1]) * dt + sigma2 * sqrt(dt) * rnorm(1)
    beta3[i] <- beta3[i-1] + kappa3 * (theta3 - beta3[i-1]) * dt + sigma3 * sqrt(dt) * rnorm(1)
  }
  
  # Create parameters data frame
  params_df <- data.frame(
    Date = dates,
    beta1 = beta1,
    beta2 = beta2,
    beta3 = beta3
  )
  
  # Calculate factor changes
  params_df$beta1_change <- c(NA, diff(params_df$beta1))
  params_df$beta2_change <- c(NA, diff(params_df$beta2))
  params_df$beta3_change <- c(NA, diff(params_df$beta3))
  
  # Estimate dynamics
  dynamics <- estimate_dns_dynamics(params_df)
  
  # Check result structure
  expect_type(dynamics, "list")
  expect_true(all(c("beta1", "beta2", "beta3", "correlation") %in% names(dynamics)))
  
  # Check estimated parameters (should be close to simulation parameters)
  expect_lt(abs(dynamics$beta1$kappa - kappa1), 0.05)
  expect_lt(abs(dynamics$beta2$kappa - kappa2), 0.1)
  expect_lt(abs(dynamics$beta3$kappa - kappa3), 0.15)
  
  expect_lt(abs(dynamics$beta1$theta - theta1), 0.01)
  expect_lt(abs(dynamics$beta2$theta - theta2), 0.02)
  expect_lt(abs(dynamics$beta3$theta - theta3), 0.02)
  
  expect_lt(abs(dynamics$beta1$sigma - sigma1), 0.003)
  expect_lt(abs(dynamics$beta2$sigma - sigma2), 0.005)
  expect_lt(abs(dynamics$beta3$sigma - sigma3), 0.008)
})
