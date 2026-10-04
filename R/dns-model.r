#' @title Dynamic Nelson-Siegel Model Implementation
#' @description Functions to calibrate and fit the Dynamic Nelson-Siegel (DNS) model to yield curve data.

#' @importFrom stats optim
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom progress progress_bar
#' @importFrom dplyr arrange filter mutate select

#' @title Fit Nelson-Siegel Model to a Single Yield Curve
#'
#' @description
#' Fits the Nelson-Siegel model to a single yield curve (cross-section) for a specific date.
#'
#' @param maturities Numeric vector. The maturities (in years) for the yield curve points.
#' @param yields Numeric vector. The yield values corresponding to the maturities.
#' @param initial_lambda Numeric. Initial value for the decay parameter lambda.
#' @param lambda_constraint Character. One of "fixed", "constrained", or "free":
#'   - "fixed": lambda is fixed at the initial_lambda value.
#'   - "constrained": lambda is optimized within bounds.
#'   - "free": lambda is optimized without constraints.
#' @param lower_bounds Numeric vector of length 4. Lower bounds for parameters (beta1, beta2, beta3, lambda).
#' @param upper_bounds Numeric vector of length 4. Upper bounds for parameters (beta1, beta2, beta3, lambda).
#' @param optimizer Character. The optimization method to use. Default is "L-BFGS-B".
#'
#' @return A list containing:
#'   - parameters: Named vector of fitted parameters (beta1, beta2, beta3, lambda)
#'   - fitted_yields: Vector of fitted yield values
#'   - rmse: Root mean squared error of the fit
#'   - maturities: The input maturities
#'   - actual_yields: The input yields
#'   - residuals: The fitting residuals (actual - fitted)
#'   - success: Logical indicating if the optimization converged
#'
#' @details
#' The Nelson-Siegel model expresses the yield curve as:
#' \deqn{y(m) = \beta_1 + \beta_2 \left(\frac{1-e^{-\lambda m}}{\lambda m}\right) + \beta_3 \left(\frac{1-e^{-\lambda m}}{\lambda m} - e^{-\lambda m}\right)}
#' where \eqn{m} is maturity, \eqn{\beta_1}, \eqn{\beta_2}, and \eqn{\beta_3} represent the level, slope, and curvature
#' factors, respectively, and \eqn{\lambda} controls the decay rate of factor loadings.
#'
#' @examples
#' \dontrun{
#' maturities <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)
#' yields <- c(0.0050, 0.0052, 0.0055, 0.0065, 0.0078, 0.0105, 0.0128, 0.0152, 0.0215, 0.0225)
#' result <- fit_ns_model(maturities, yields)
#' }
#'
#' @export
fit_ns_model <- function(maturities, yields, initial_lambda = 0.0609,
                         lambda_constraint = "constrained",
                         lower_bounds = c(-0.1, -0.5, -0.5, 0.01),
                         upper_bounds = c(0.2, 0.5, 0.5, 0.5),
                         optimizer = "L-BFGS-B") {

  # Parameter validation. Missing values are permitted here: observations
  # with NA maturities or yields are filtered below, with a graceful
  # success = FALSE result when fewer than 3 valid points remain.
  checkmate::assert_numeric(maturities, lower = 0)
  checkmate::assert_numeric(yields)
  checkmate::assert_numeric(initial_lambda, lower = 0.001, upper = 5, len = 1)
  checkmate::assert_choice(lambda_constraint, c("fixed", "constrained", "free"))
  checkmate::assert_numeric(lower_bounds, len = 4)
  checkmate::assert_numeric(upper_bounds, len = 4)
  checkmate::assert_character(optimizer, len = 1)

  # Make sure inputs are numeric and have same length
  maturities <- as.numeric(maturities)
  yields <- as.numeric(yields)

  # Check for length mismatch
  if (length(maturities) != length(yields)) {
    stop("Length of maturities and yields must be the same")
  }

  # Remove any NA values or zero/negative maturities
  valid_idx <- !is.na(maturities) & !is.na(yields) & maturities > 0
  if (sum(valid_idx) < 3) {
    return(list(
      parameters = c(beta1 = NA, beta2 = NA, beta3 = NA, lambda = NA),
      fitted_yields = rep(NA, length(maturities)),
      rmse = NA,
      maturities = maturities,
      actual_yields = yields,
      residuals = rep(NA, length(maturities)),
      success = FALSE
    ))
  }

  maturities <- maturities[valid_idx]
  yields <- yields[valid_idx]

  # Create the factor loading calculation function
  calculate_loadings <- function(m, lambda) {
    # Avoid numerical issues when m*lambda is very small
    if (m * lambda < 1e-10) {
      loading_2 <- 1.0  # limit of (1-exp(-λm))/(λm) as λm → 0
      loading_3 <- 0.0  # limit of loading_2 - exp(-λm) as λm → 0
    } else {
      # Standard formula
      loading_2 <- (1 - exp(-lambda * m)) / (lambda * m)
      loading_3 <- loading_2 - exp(-lambda * m)
    }
    return(c(loading_2, loading_3))
  }

  # Function to calculate Nelson-Siegel yields
  ns_yields <- function(params, maturities) {
    beta1 <- params[1]
    beta2 <- params[2]
    beta3 <- params[3]
    lambda <- params[4]

    sapply(maturities, function(m) {
      loadings <- calculate_loadings(m, lambda)
      beta1 + beta2 * loadings[1] + beta3 * loadings[2]
    })
  }

  # Error function to minimize (sum of squared errors)
  error_fn <- function(params) {
    # Calculate model yields
    model_yields <- ns_yields(params, maturities)

    # Return sum of squared errors
    sum((yields - model_yields)^2)
  }

  # Initial parameter guess
  initial_params <- c(
    beta1 = mean(yields[maturities >= 5], na.rm = TRUE),  # Level (avg of long-term yields)
    beta2 = yields[1] - yields[length(yields)],           # Slope (short - long)
    beta3 = 0,                                           # Curvature
    lambda = initial_lambda                              # Decay rate
  )

  # Handle NAs in initial guess
  if (is.na(initial_params[1])) initial_params[1] <- mean(yields)
  if (is.na(initial_params[2])) initial_params[2] <- 0

  # Adjust optimization setup based on lambda constraint
  if (lambda_constraint == "fixed") {
    # Only optimize beta parameters
    opt_fn <- function(betas) {
      params <- c(betas, initial_lambda)
      error_fn(params)
    }

    opt_result <- optim(
      par = initial_params[1:3],
      fn = opt_fn,
      method = optimizer,
      lower = lower_bounds[1:3],
      upper = upper_bounds[1:3]
    )

    # Combine optimized betas with fixed lambda
    params <- c(opt_result$par, initial_lambda)
    names(params) <- c("beta1", "beta2", "beta3", "lambda")

  } else if (lambda_constraint == "constrained") {
    # Optimize all parameters with constraints
    opt_result <- optim(
      par = initial_params,
      fn = error_fn,
      method = optimizer,
      lower = lower_bounds,
      upper = upper_bounds
    )

    # Extract optimized parameters
    params <- opt_result$par
    names(params) <- c("beta1", "beta2", "beta3", "lambda")

  } else if (lambda_constraint == "free") {
    # Optimize without bounds
    opt_result <- optim(
      par = initial_params,
      fn = error_fn,
      method = "BFGS"  # Method that doesn't require bounds
    )

    # Extract optimized parameters
    params <- opt_result$par
    names(params) <- c("beta1", "beta2", "beta3", "lambda")
  }

  # Calculate fitted yields with the optimized parameters
  fitted_yields <- ns_yields(params, maturities)

  # Calculate residuals
  residuals <- yields - fitted_yields

  # Calculate RMSE
  rmse <- sqrt(mean(residuals^2))

  # Prepare return values
  result <- list(
    parameters = params,
    fitted_yields = fitted_yields,
    rmse = rmse,
    maturities = maturities,
    actual_yields = yields,
    residuals = residuals,
    success = opt_result$convergence == 0
  )

  return(result)
}

#' @title Fit Dynamic Nelson-Siegel Model to Time Series of Yield Curves
#'
#' @description
#' Fits the DNS model to a time series of yield curves.
#'
#' @param yield_data Data frame. The yield curve data with first column 'Date' and remaining
#'   columns containing yield values for different maturities.
#' @param n_days Integer. Number of most recent days to include in the fitting.
#' @param parallel Logical. If TRUE, use parallel processing to speed up computation.
#' @param lambda_constraint Character. How to handle the lambda parameter (fixed, constrained, free).
#' @param initial_lambda Numeric. Initial value for lambda parameter.
#' @param show_progress Logical. If TRUE, display a progress bar.
#'
#' @return A list containing:
#'   - parameters: Data frame with model parameters for each date.
#'   - fitted_yields: Matrix of fitted yield values.
#'   - maturities: Vector of maturities in years.
#'   - success_rate: Percentage of successful model fits.
#'
#' @examples
#' \dontrun{
#' yield_data <- get_treasury_yield_combined("2020-01-01")
#' dns_result <- fit_dns_model(yield_data, n_days = 30)
#' }
#'
#' @export
fit_dns_model <- function(yield_data, n_days = 60, parallel = FALSE,
                          lambda_constraint = "constrained",
                          initial_lambda = 0.0609,
                          show_progress = TRUE) {

  # Parameter validation
  checkmate::assert_data_frame(yield_data)
  checkmate::assert_names(names(yield_data), must.include = "Date")
  checkmate::assert_numeric(n_days, lower = 1)
  checkmate::assert_logical(parallel, len = 1)
  checkmate::assert_logical(show_progress, len = 1)

  # Extract dates
  dates <- yield_data$Date[(nrow(yield_data) - min(as.integer(n_days), nrow(yield_data)) + 1):nrow(yield_data)]

  # Get maturities from column names
  maturities <- as.numeric(colnames(yield_data)[-1])

  # Initialize storage for parameters
  params_df <- data.frame(
    Date = dates,
    beta1 = numeric(length(dates)),
    beta2 = numeric(length(dates)),
    beta3 = numeric(length(dates)),
    lambda = numeric(length(dates)),
    rmse = numeric(length(dates)),
    stringsAsFactors = FALSE
  )

  # Prepare storage for fitted yields
  fitted_yields_matrix <- matrix(NA, nrow = length(dates), ncol = length(maturities))
  colnames(fitted_yields_matrix) <- colnames(yield_data)[-1]
  rownames(fitted_yields_matrix) <- as.character(dates)

  # Setup parallel processing if requested
  if (parallel) {
    oplan <- future::plan(future::multisession)
    on.exit(future::plan(oplan), add = TRUE)
  } else {
    future::plan(future::sequential)
  }

  # Setup progress bar if requested
  if (show_progress) {
    pb <- progress::progress_bar$new(
      format = "[:bar] :percent Fitted: :current/:total (:eta remaining)",
      total = length(dates)
    )
  }

  # Process each date
  fit_results <- future.apply::future_lapply(seq_along(dates), function(i) {
    date <- dates[i]

    # Update progress bar
    if (show_progress && !parallel) pb$tick()

    # Extract yields for this date
    date_idx <- which(yield_data$Date == date)

    if (length(date_idx) == 0) {
      return(list(success = FALSE, message = "No data found for date"))
    }

    # Extract all columns except the date column for this specific row
    yields <- as.numeric(as.vector(yield_data[date_idx, -1]))

    # Check if we have valid data
    if (length(yields) == 0 || sum(!is.na(yields)) < 3) {
      return(list(success = FALSE, message = "Not enough valid yield data"))
    }

    # Get previous lambda for better initialization if available
    last_lambda <- NULL
    if (i > 1 && exists("last_fit_result") && last_fit_result$success) {
      last_lambda <- last_fit_result$parameters["lambda"]
    }

    # Use the provided initial lambda if no previous fit is available
    lambda_init <- if (!is.null(last_lambda)) last_lambda else initial_lambda

    # Try to fit the model
    result <- tryCatch({
      fit_ns_model(
        maturities = maturities,
        yields = yields,
        initial_lambda = lambda_init,
        lambda_constraint = lambda_constraint
      )
    }, error = function(e) {
      return(list(success = FALSE, message = e$message))
    })

    # Store the result for next iteration if successful
    if (!is.list(result)) {
      return(list(success = FALSE, message = "Fit failed with unknown error"))
    }

    # This result will be used for the next iteration's initialization
    last_fit_result <<- result

    # Create a complete result with date information
    result$date <- date
    result$date_index <- i

    return(result)
  })

  # If parallel and progress bar, update afterward
  if (show_progress && parallel) {
    pb$update(1)
  }

  # Process results
  successful_fits <- 0
  error_messages <- list()

  for (i in seq_along(fit_results)) {
    result <- fit_results[[i]]

    if (!result$success) {
      # Store error message
      error_messages[[length(error_messages) + 1]] <-
        paste("Date", dates[i], ":", result$message)
      next
    }

    # Store parameters
    params_df$beta1[i] <- result$parameters["beta1"]
    params_df$beta2[i] <- result$parameters["beta2"]
    params_df$beta3[i] <- result$parameters["beta3"]
    params_df$lambda[i] <- result$parameters["lambda"]
    params_df$rmse[i] <- result$rmse

    # Store fitted yields
    fitted_yields_matrix[i, ] <- result$fitted_yields

    successful_fits <- successful_fits + 1
  }

  # Report fitting results
  if (length(error_messages) > 0 && length(error_messages) <= 5) {
    for (msg in error_messages) {
      warning(msg)
    }
  } else if (length(error_messages) > 5) {
    warning(sprintf("%d fitting errors occurred. First 5: %s",
                   length(error_messages),
                   paste(error_messages[1:5], collapse = "; ")))
  }

  if (successful_fits == 0) {
    stop("Failed to fit model to any yield curves")
  }

  # Remove rows with NA parameters
  valid_rows <- !is.na(params_df$beta1)
  params_df <- params_df[valid_rows, ]
  fitted_yields_matrix <- fitted_yields_matrix[valid_rows, ]

  # Calculate daily changes in factors
  if (nrow(params_df) > 1) {
    params_df$beta1_change <- c(NA, diff(params_df$beta1))
    params_df$beta2_change <- c(NA, diff(params_df$beta2))
    params_df$beta3_change <- c(NA, diff(params_df$beta3))
  } else {
    params_df$beta1_change <- NA
    params_df$beta2_change <- NA
    params_df$beta3_change <- NA
    warning("Only one valid yield curve - cannot calculate factor changes")
  }

  # Calculate success rate
  success_rate <- successful_fits / length(dates) * 100

  return(list(
    parameters = params_df,
    fitted_yields = fitted_yields_matrix,
    maturities = maturities,
    success_rate = success_rate
  ))
}
