#' @title Dynamic Nelson-Siegel Model Dynamics
#' @description Functions to estimate and analyze the dynamics of DNS model factors.

#' @importFrom stats lm coef sd resid na.omit
#' @importFrom dplyr filter select
#' @importFrom checkmate assert_data_frame assert_subset assert_class
#' @importFrom lubridate today
#' @importFrom utils head

#' @title Estimate Dynamics of the DNS Factors
#'
#' @description
#' Estimates the stochastic dynamics of the factors in the Dynamic Nelson-Siegel model
#' based on historical factor time series.
#'
#' @param params_df Data frame. The parameters data frame from fit_dns_model() containing
#'   columns for Date, beta1, beta2, beta3, and the corresponding change columns.
#' @param min_observations Integer. Minimum number of observations required for estimation.
#' @param dt Numeric. Time step in years. Default is 1/252 (daily data).
#'
#' @return A list containing:
#'   - For each factor (beta1, beta2, beta3):
#'     - kappa: Mean reversion speed
#'     - theta: Long-run mean
#'     - sigma: Volatility
#'     - model: Fitted regression model
#'     - half_life: Half-life of mean reversion (in days)
#'     - summary: Basic statistics about the factor
#'   - correlation: Correlation matrix of factor innovations
#'   - innovations: Matrix of factor innovations (residuals from mean reversion)
#'   - estimation_period: Date range used for estimation
#'
#' @details
#' This function implements the estimation of an Ornstein-Uhlenbeck process for each factor:
#' \deqn{d\beta(t) = \kappa(\theta - \beta(t))dt + \sigma dW(t)}
#' The discretized version used for estimation is:
#' \deqn{\Delta\beta(t) = \kappa(\theta - \beta(t-1))\Delta t + \sigma\epsilon(t)}
#'
#' @examples
#' \dontrun{
#' yield_data <- get_treasury_yield_combined("2020-01-01")
#' dns_result <- fit_dns_model(yield_data, n_days = 100)
#' dynamics <- estimate_dns_dynamics(dns_result$parameters)
#' }
#'
#' @export
estimate_dns_dynamics <- function(params_df, min_observations = 20, dt = 1/252) {
  # Validate inputs
  checkmate::assert_data_frame(params_df)
  checkmate::assert_subset(
    c("Date", "beta1", "beta2", "beta3", "beta1_change", "beta2_change", "beta3_change"),
    names(params_df)
  )
  checkmate::assert_class(params_df$Date, "Date")
  checkmate::assert_numeric(dt, lower = 0, len = 1)
  checkmate::assert_numeric(min_observations, lower = 10, len = 1)

  # Check if we have enough data
  if (nrow(params_df) < as.integer(min_observations)) {
    stop("Not enough observations for reliable dynamics estimation. Need at least ",
         min_observations, " observations, but only have ", nrow(params_df))
  }

  # Extract factor time series (excluding first row which has no changes)
  factors_df <- params_df[2:nrow(params_df), ]
  factors_df <- na.omit(factors_df)  # Remove any rows with NA values

  # Check if we still have enough data after removing NA values
  if (nrow(factors_df) < as.integer(min_observations)) {
    stop("Not enough complete observations after removing NA values")
  }

  # Initialize result storage
  dynamics <- list()
  innovations <- matrix(NA, nrow = nrow(factors_df), ncol = 3)
  colnames(innovations) <- c("beta1", "beta2", "beta3")

  # Extract date range for reporting
  date_range <- range(factors_df$Date)
  dynamics$estimation_period <- date_range

  # For each factor, estimate the mean reversion parameters
  factor_names <- c("beta1", "beta2", "beta3")
  change_names <- c("beta1_change", "beta2_change", "beta3_change")

  for (i in 1:3) {
    factor_name <- factor_names[i]
    change_name <- change_names[i]

    # Prepare regression data
    reg_data <- data.frame(
      factor_level = factors_df[[factor_name]],
      factor_change = factors_df[[change_name]]
    )

    # Calculate basic statistics
    factor_summary <- list(
      mean = mean(reg_data$factor_level),
      sd = sd(reg_data$factor_level),
      min = min(reg_data$factor_level),
      max = max(reg_data$factor_level),
      q25 = quantile(reg_data$factor_level, 0.25),
      q50 = quantile(reg_data$factor_level, 0.50),
      q75 = quantile(reg_data$factor_level, 0.75)
    )

    # Linear regression: Δβ(t) ~ κθΔt - κβ(t-1)Δt
    model <- lm(factor_change ~ factor_level, data = reg_data)

    # FIX: Extract parameters with more robust handling of NA/NULL values
    # Check if coefficients are NA or NULL
    if (is.null(coef(model)) || length(coef(model)) < 2 || any(is.na(coef(model)))) {
      warning(paste("Regression failed for factor", factor_name, "- using fallback values"))
      kappa <- 0.001  # Small positive value as fallback
      theta <- factor_summary$mean
      sigma <- factor_summary$sd * sqrt(dt)
    } else {
      # Normal calculation
      kappa <- -coef(model)[2] / dt

      # FIX: Check if kappa is valid before using it in division
      if (is.na(kappa) || kappa <= 0) {
        # If kappa is NA or not positive, use a fallback approach
        kappa <- max(0.001, abs(kappa))  # Ensure positive value
        theta <- factor_summary$mean
      } else {
        # Normal calculation with valid kappa
        theta <- -coef(model)[1] / coef(model)[2]
      }

      sigma <- sd(residuals(model)) / sqrt(dt)

      # Handle missing or invalid values
      if (is.na(sigma) || sigma <= 0) {
        sigma <- 0.001  # Fallback to small positive value
      }
    }

    # Calculate half-life of mean reversion
    half_life_years <- if (kappa > 0) log(2) / kappa else Inf
    half_life_days <- half_life_years * 252  # Convert to trading days

    # Calculate pseudo-R² (model explanatory power)
    ss_total <- sum((reg_data$factor_change - mean(reg_data$factor_change))^2)
    ss_residual <- sum(residuals(model)^2)
    r_squared <- 1 - ss_residual / ss_total

    # Store results
    dynamics[[factor_name]] <- list(
      kappa = kappa,
      theta = theta,
      sigma = sigma,
      model = model,
      half_life_years = half_life_years,
      half_life_days = half_life_days,
      r_squared = r_squared,
      summary = factor_summary
    )

    # Store innovations (residuals from mean reversion model)
    innovations[, i] <- residuals(model)
  }

  # Estimate correlation between factor innovations
  cor_matrix <- cor(innovations)
  rownames(cor_matrix) <- colnames(cor_matrix) <- factor_names

  # Add correlation to results
  dynamics$correlation <- cor_matrix
  dynamics$innovations <- innovations

  # Analyze persistence of the model
  eigenvalues <- kappa_matrix_eigenvalues(dynamics)
  dynamics$eigenvalues <- eigenvalues
  dynamics$system_stability <- all(Re(eigenvalues) > 0)

  # Return complete dynamics structure
  return(dynamics)
}

#' Calculate the eigenvalues of the mean reversion matrix
#'
#' @param dynamics A dynamics list as returned by estimate_dns_dynamics
#' @return Vector of eigenvalues
#' @keywords internal
kappa_matrix_eigenvalues <- function(dynamics) {
  # Check if the dynamics are in the expected format
  if (!all(c("beta1", "beta2", "beta3") %in% names(dynamics))) {
    stop("Dynamics object does not have the expected structure")
  }

  # Extract kappa values from dynamics
  kappa_matrix <- diag(c(
    dynamics$beta1$kappa,
    dynamics$beta2$kappa,
    dynamics$beta3$kappa
  ))

  # Calculate eigenvalues
  eigenvalues <- eigen(kappa_matrix, only.values = TRUE)$values

  return(eigenvalues)
}

#' @title Convert DNS Dynamics to Simulation Parameters
#'
#' @description
#' Converts the estimated DNS dynamics to parameters suitable for simulation.
#'
#' @param dynamics List. The dynamics object from estimate_dns_dynamics().
#' @param include_correlation Logical. If TRUE, include the correlation structure.
#'
#' @return A list of parameters for DNS simulation.
#'
#' @examples
#' \dontrun{
#' yield_data <- get_treasury_yield_combined("2020-01-01")
#' dns_result <- fit_dns_model(yield_data, n_days = 100)
#' dynamics <- estimate_dns_dynamics(dns_result$parameters)
#' sim_params <- dns_dynamics_to_sim_params(dynamics)
#' }
#'
#' @export
dns_dynamics_to_sim_params <- function(dynamics, include_correlation = TRUE) {
  # Validate input
  if (!all(c("beta1", "beta2", "beta3", "correlation") %in% names(dynamics))) {
    stop("Invalid dynamics object. Must include beta1, beta2, beta3, and correlation components.")
  }

  # Extract parameters
  kappa1 <- dynamics$beta1$kappa
  kappa2 <- dynamics$beta2$kappa
  kappa3 <- dynamics$beta3$kappa

  theta1 <- dynamics$beta1$theta
  theta2 <- dynamics$beta2$theta
  theta3 <- dynamics$beta3$theta

  sigma1 <- dynamics$beta1$sigma
  sigma2 <- dynamics$beta2$sigma
  sigma3 <- dynamics$beta3$sigma

  # Package parameters for simulation
  params <- list(
    kappa = c(kappa1, kappa2, kappa3),
    theta = c(theta1, theta2, theta3),
    sigma = c(sigma1, sigma2, sigma3),
    names = c("Level", "Slope", "Curvature"),
    lambda = 0.0609  # Default lambda if not specified
  )

  # Add correlation if requested
  if (include_correlation) {
    params$correlation <- dynamics$correlation
  }

  # Check for the presence of lambda in dns_result parameters
  if ("lambda" %in% names(dynamics)) {
    params$lambda <- dynamics$lambda
  }

  # Safety check for unreasonable parameters
  if (any(params$kappa < 0)) {
    warning("Negative mean reversion speeds detected. Setting to small positive value.")
    params$kappa[params$kappa < 0] <- 0.001
  }

  if (any(params$sigma <= 0)) {
    warning("Non-positive volatilities detected. Setting to small positive value.")
    params$sigma[params$sigma <= 0] <- 0.001
  }

  return(params)
}

#' @title Get Calibrated DNS Parameters
#'
#' @description
#' Returns calibrated DNS parameters based on the most recent market data,
#' suitable for use in simulation or pricing applications.
#'
#' @param lookback_days Integer. Number of days to use for calibration.
#' @param lambda Numeric. Fixed lambda value to use for calibration.
#' @param start_date Date or character. Start date for the data window.
#' @param end_date Date or character. End date for the data window.
#' @param verbose Logical. If TRUE, print progress information.
#'
#' @return A list containing:
#'   - params: List of DNS parameters for simulation
#'   - initial_factors: Vector of current factor values
#'   - lambda: Fitted lambda value
#'   - calibration_summary: Summary statistics about the calibration
#'
#' @examples
#' \dontrun{
#' dns_params <- get_calibrated_dns_params(lookback_days = 60)
#' }
#'
#' @export
get_calibrated_dns_params <- function(lookback_days = 60, lambda = NULL,
                                     start_date = NULL, end_date = NULL,
                                     verbose = FALSE) {
  # Validate inputs
  checkmate::assert_numeric(lookback_days, lower = 10)
  if (!is.null(lambda)) checkmate::assert_numeric(lambda, lower = 0.001, upper = 1)
  checkmate::assert_logical(verbose, len = 1)

  # Set default end date if not specified
  if (is.null(end_date)) end_date <- Sys.Date()
  end_date <- as.Date(end_date)

  # Set default start date if not specified
  if (is.null(start_date)) {
    # Calculate a start date that will give us enough data
    start_date <- end_date - as.integer(lookback_days) * 1.5
  }
  start_date <- as.Date(start_date)

  # Get yield curve data for the specified period
  if (verbose) cat("Getting yield curve data...\n")
  yield_data <- get_treasury_yield_combined(
    start_date = start_date,
    end_date = end_date,
    verbose = verbose
  )

  # Check if we have enough data
  if (nrow(yield_data) < lookback_days / 2) {
    warning(sprintf("Only %d days of data available. Proceeding with limited data.",
                   nrow(yield_data)))
  }

  # Fit DNS model to the yield curve data
  if (verbose) cat("Fitting DNS model...\n")

  # Set lambda constraint based on input
  lambda_constraint <- if (is.null(lambda)) "constrained" else "fixed"
  initial_lambda <- if (is.null(lambda)) 0.0609 else lambda

  dns_result <- fit_dns_model(
    yield_data = yield_data,
    n_days = min(nrow(yield_data), as.integer(lookback_days)),
    lambda_constraint = lambda_constraint,
    initial_lambda = initial_lambda,
    show_progress = verbose
  )

  # Check if we have valid results
  if (nrow(dns_result$parameters) < 5) {
    warning("Limited DNS fitting results. Model parameters may not be reliable.")
  }

  # Estimate factor dynamics if we have enough data
  if (nrow(dns_result$parameters) >= 20) {
    if (verbose) cat("Estimating factor dynamics...\n")
    dynamics <- tryCatch({
      estimate_dns_dynamics(dns_result$parameters)
    }, error = function(e) {
      warning("Error estimating dynamics: ", e$message,
              ". Using default parameters.")
      NULL
    })
  } else {
    if (verbose) cat("Not enough data for dynamics estimation. Using default values.\n")
    dynamics <- NULL
  }

  # If dynamics estimation failed or we don't have enough data, use default values
  if (is.null(dynamics)) {
    # Use default parameters
    sim_params <- list(
      kappa = c(0.02, 0.3, 0.5),     # Speed of mean reversion
      theta = c(0.04, -0.02, 0.01),  # Long-run mean
      sigma = c(0.005, 0.01, 0.02),  # Volatility
      names = c("Level", "Slope", "Curvature")
    )
  } else {
    # Convert dynamics to simulation parameters
    sim_params <- dns_dynamics_to_sim_params(dynamics)
  }

  # Get latest factors (for initial conditions)
  latest_factors <- as.numeric(dns_result$parameters[1, c("beta1", "beta2", "beta3")])

  # Get fitted lambda
  lambda_value <- dns_result$parameters$lambda[1]

  # Create summary information
  calibration_summary <- list(
    calibration_date = Sys.Date(),
    data_date_range = range(yield_data$Date),
    n_observations = nrow(dns_result$parameters),
    maturities = dns_result$maturities,
    rmse_bp = median(dns_result$parameters$rmse) * 10000  # In basis points
  )

  # Return the calibrated parameters
  return(list(
    params = sim_params,
    initial_factors = latest_factors,
    lambda = lambda_value,
    calibration_summary = calibration_summary
  ))
}
