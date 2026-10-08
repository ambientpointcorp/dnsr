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
  expect_true(all(c("parameters", "fitted_yields", "rmse", "maturities",
                    "actual_yields", "residuals", "success", "valid_idx",
                    "at_bound") %in% names(result)))
  
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
  expect_true(all(c("parameters", "fitted_yields", "rmse", "maturities",
                    "actual_yields", "residuals", "success", "valid_idx",
                    "at_bound") %in% names(result)))
  
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
  expect_equal(unname(result$parameters["lambda"]), fixed_lambda)
})

test_that("fit_ns_model handles extreme yield curves gracefully", {
  # Create an extreme, but realistic yield curve (inverted)
  maturities <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)
  yields <- c(0.05, 0.049, 0.048, 0.045, 0.043, 0.04, 0.038, 0.035, 0.033, 0.032)
  
  # Fit Nelson-Siegel model
  result <- fit_ns_model(maturities, yields)
  expect_false(any(result$at_bound))
  
  # Check result structure
  expect_type(result, "list")
  expect_true(all(c("parameters", "fitted_yields", "rmse", "maturities",
                    "actual_yields", "residuals", "success", "valid_idx",
                    "at_bound") %in% names(result)))
  
  # Check parameters - beta2 should be positive for an inverted curve: the
  # slope loading (1 - exp(-lambda*m)) / (lambda*m) is 1 at the short end
  # and decays toward 0, so short ~ beta1 + beta2 and long ~ beta1, and an
  # inverted curve (short above long) needs beta2 > 0.
  expect_gt(result$parameters["beta2"], 0)
  
  # Create another extreme case - humped yield curve. This hump is sharp enough
  # that the fit wants a faster decay than the default lambda ceiling of 0.5
  # allows, so it comes to rest there. L-BFGS-B reports convergence at a bound,
  # so the fit says so explicitly instead of returning an unqualified success.
  yields_humped <- c(0.01, 0.015, 0.02, 0.025, 0.028, 0.026, 0.024, 0.023, 0.022, 0.021)
  
  # Fit Nelson-Siegel model
  expect_warning(result_humped <- fit_ns_model(maturities, yields_humped),
                 "rests on a parameter bound")
  expect_true(result_humped$at_bound[["lambda"]])
  expect_equal(unname(result_humped$parameters["lambda"]), 0.5)
  
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
  expect_true(all(c("parameters", "fitted_yields", "maturities", "success_rate",
                    "failures") %in% names(dns_result)))
  
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
  # Simulate a noiseless OU path long enough, relative to the half-lives,
  # for the parameters to be identifiable: 10 years against half-lives of
  # 0.2-1.4 years. (Slower reversion, e.g. kappa = 0.02, has a 35-year
  # half-life that no realistic sample can recover.) kappa estimates carry
  # the well-known upward small-sample bias of order 5 / T_years plus
  # sampling noise of order sqrt(2 * kappa / T_years), which the
  # tolerances below reflect; theta and sigma are recovered tightly.
  n_days <- 2520
  dt <- 1 / 252
  kappa <- c(0.5, 1.5, 3.0)
  theta <- c(0.04, -0.02, 0.01)
  sigma <- c(0.005, 0.01, 0.02)

  set.seed(123)
  b <- matrix(NA_real_, n_days, 3)
  b[1, ] <- theta
  for (i in 2:n_days) {
    b[i, ] <- b[i - 1, ] + kappa * (theta - b[i - 1, ]) * dt +
      sigma * sqrt(dt) * rnorm(3)
  }

  params_df <- data.frame(
    Date = seq(as.Date("2014-01-01"), by = "day", length.out = n_days),
    beta1 = b[, 1],
    beta2 = b[, 2],
    beta3 = b[, 3]
  )
  params_df$beta1_change <- c(NA, diff(params_df$beta1))
  params_df$beta2_change <- c(NA, diff(params_df$beta2))
  params_df$beta3_change <- c(NA, diff(params_df$beta3))

  # Estimate dynamics on the daily series as given: this test validates the
  # raw OU regression on a noiseless simulated path, so no subsampling.
  dynamics <- estimate_dns_dynamics(params_df, sampling = "daily")

  # Check result structure
  expect_type(dynamics, "list")
  expect_true(all(c("beta1", "beta2", "beta3", "correlation") %in% names(dynamics)))

  # Check estimated parameters (close to simulation parameters, within the
  # small-sample bias + noise bands above)
  expect_lt(abs(dynamics$beta1$kappa - kappa[1]), 1.5)
  expect_lt(abs(dynamics$beta2$kappa - kappa[2]), 2.5)
  expect_lt(abs(dynamics$beta3$kappa - kappa[3]), 3.0)

  expect_lt(abs(dynamics$beta1$theta - theta[1]), 0.02)
  expect_lt(abs(dynamics$beta2$theta - theta[2]), 0.02)
  expect_lt(abs(dynamics$beta3$theta - theta[3]), 0.02)

  expect_lt(abs(dynamics$beta1$sigma - sigma[1]), 0.003)
  expect_lt(abs(dynamics$beta2$sigma - sigma[2]), 0.005)
  expect_lt(abs(dynamics$beta3$sigma - sigma[3]), 0.008)
})
