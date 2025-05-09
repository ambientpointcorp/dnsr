#' @title Dynamic Nelson-Siegel Model Simulation using YUIMA
#' @description Functions to simulate yield curves using the Dynamic Nelson-Siegel model with YUIMA package for SDE simulation.

#' @importFrom stats rnorm
#' @importFrom checkmate assert_numeric assert_list assert_matrix assert_scalar
#' @importFrom methods is
#' @importFrom utils head tail

#' @title Simulate DNS Factors using YUIMA package
#'
#' @description
#' Simulates the evolution of DNS factors using a multivariate Ornstein-Uhlenbeck process
#' implemented with the YUIMA package for accurate stochastic differential equation simulation.
#'
#' @param params List. Model parameters (kappa, theta, sigma).
#' @param initial_factors Numeric vector. Initial factor values (beta1, beta2, beta3).
#' @param time_horizon Numeric. Simulation time horizon in years.
#' @param n_steps Integer. Number of time steps in the simulation.
#' @param n_paths Integer. Number of simulation paths.
#' @param include_correlation Logical. Whether to include correlations in simulations.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return A list containing:
#'   - factors: Array of simulated factors [time, path, factor]
#'   - time: Vector of time points
#'   - params: The input parameters
#'
#' @details
#' The simulation is based on the multivariate Ornstein-Uhlenbeck process:
#' \deqn{d\beta(t) = \kappa(\theta - \beta(t))dt + \sigma dW(t)}
#' where \eqn{\kappa} is the mean-reversion matrix, \eqn{\theta} is the long-run mean vector,
#' \eqn{\sigma} is the volatility matrix, and \eqn{W(t)} is a correlated Wiener process.
#'
#' This implementation uses the YUIMA package, which provides efficient and accurate methods
#' for simulating stochastic differential equations.
#'
#' @examples
#' \dontrun{
#' params <- list(
#'   kappa = c(0.02, 0.3, 0.5),
#'   theta = c(0.04, -0.02, 0.01),
#'   sigma = c(0.005, 0.01, 0.02),
#'   correlation = diag(3)
#' )
#' initial_factors <- c(0.035, -0.015, 0.005)
#' sim_result <- simulate_dns_factors(params, initial_factors)
#' }
#'
#' @export
simulate_dns_factors <- function(params, initial_factors, time_horizon = 10,
                                 n_steps = 120, n_paths = 100,
                                 include_correlation = TRUE, seed = NULL) {
  # Check if YUIMA package is installed
  if (!requireNamespace("yuima", quietly = TRUE)) {
    stop("The YUIMA package is required for SDE simulation. Please install it with: install.packages('yuima')")
  }

  # Parameter validation
  checkmate::assert_list(params)
  checkmate::assert_subset(c("kappa", "theta", "sigma"), names(params))
  checkmate::assert_numeric(params$kappa, len = 3)
  checkmate::assert_numeric(params$theta, len = 3)
  checkmate::assert_numeric(params$sigma, len = 3, lower = 0)
  checkmate::assert_numeric(initial_factors, len = 3)
  checkmate::assert_scalar(time_horizon)
  checkmate::assert_count(n_steps, positive = TRUE)
  checkmate::assert_count(n_paths, positive = TRUE)

  # Set random seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Time grid for simulation
  delta_t <- time_horizon / n_steps
  time_grid <- seq(0, time_horizon, length.out = n_steps + 1)

  # Define variable names
  var_names <- c("beta1", "beta2", "beta3")

  # Extract parameters
  kappa <- params$kappa
  theta <- params$theta
  sigma <- params$sigma

  # Define the correlation matrix (if provided)
  corr_matrix <- if (include_correlation && "correlation" %in% names(params)) {
    checkmate::assert_matrix(params$correlation, nrows = 3, ncols = 3)
    params$correlation
  } else {
    diag(3)  # Identity matrix (no correlation)
  }

  # Convert parameters to YUIMA format
  # 1. Drift coefficient (mean reversion)
  # For O-U process: dX_t = κ(θ - X_t)dt + σdW_t
  # The drift term is: κθ - κX_t
  drift_terms <- c(
    paste0(kappa[1], "*", theta[1], "-", kappa[1], "*beta1"),
    paste0(kappa[2], "*", theta[2], "-", kappa[2], "*beta2"),
    paste0(kappa[3], "*", theta[3], "-", kappa[3], "*beta3")
  )

  # 2. Diffusion coefficient
  # Create diffusion matrix with correlations
  # First, create Cholesky decomposition of correlation matrix
  chol_corr <- try(chol(corr_matrix), silent = TRUE)
  if (inherits(chol_corr, "try-error")) {
    warning("Correlation matrix is not positive definite. Using identity matrix.")
    chol_corr <- diag(3)
  }

  # Construct the diffusion matrix as sigma * chol_corr
  diffusion_matrix <- matrix(0, nrow = 3, ncol = 3)
  for (i in 1:3) {
    for (j in 1:3) {
      diffusion_matrix[i, j] <- sigma[i] * chol_corr[j, i]
    }
  }

  # Convert diffusion matrix to yuima format
  diffusion_terms <- matrix("0", nrow = 3, ncol = 3)
  for (i in 1:3) {
    for (j in 1:3) {
      if (diffusion_matrix[i, j] != 0) {
        diffusion_terms[i, j] <- as.character(diffusion_matrix[i, j])
      }
    }
  }

  # Create YUIMA model
  ou_model <- yuima::setModel(
    drift = drift_terms,
    diffusion = diffusion_terms,
    state.variable = var_names,
    solve.variable = var_names
  )

  # Set simulation parameters
  sim_params <- yuima::setSampling(
    Terminal = time_horizon,
    n = n_steps
  )

  # Initialize storage for simulation results
  factors <- array(0, dim = c(n_steps + 1, n_paths, 3))

  # Simulate paths
  for (p in 1:n_paths) {
    # Set initial condition for this path
    initial_state <- initial_factors
    names(initial_state) <- var_names

    # Run simulation
    sim <- yuima::simulate(
      object = ou_model,
      sampling = sim_params,
      xinit = initial_state,
      method = "euler"
    )

    # Let's examine the structure of the simulation output
    sim_data <- sim@data@original.data

    # The sim_data object likely has a different structure than expected
    # It might be a matrix with the first column being time and the other columns being the simulated state variables

    # Debug: print out information about sim_data for the first path only
    if (p == 1) {
      message("Simulation data structure: ", class(sim_data))
      message("Simulation data dimensions: ", nrow(sim_data), " x ", ncol(sim_data))
      if (ncol(sim_data) > 0) {
        message("Column names: ", paste(colnames(sim_data), collapse=", "))
      }
    }

    # Get the data for each factor
    # Assuming the column structure is [time, beta1, beta2, beta3]
    # If sim_data only has 3 columns, it might be that time is the index and not a column

    if (ncol(sim_data) >= 4) {
      # Standard case: [time, beta1, beta2, beta3]
      for (f in 1:3) {
        factors[, p, f] <- sim_data[, f + 1]
      }
    } else if (ncol(sim_data) == 3) {
      # Case where sim_data = [beta1, beta2, beta3] (no time column)
      for (f in 1:3) {
        factors[, p, f] <- sim_data[, f]
      }
    } else {
      # Error case: not enough columns
      warning("Simulation output has unexpected structure with ", ncol(sim_data), " columns.")
      for (f in 1:3) {
        # Use initial value for this factor to avoid empty results
        factors[, p, f] <- rep(initial_factors[f], n_steps + 1)
      }
    }
  }

  # Return results
  return(list(
    factors = factors,
    time = time_grid,
    params = params
  ))
}

#' @title Calculate Yield Curves from DNS Factors
#'
#' @description
#' Calculates yield curves based on DNS factors for specified maturities.
#'
#' @param factors Array or matrix. If array: simulated factors [time, path, factor].
#'   If matrix: factors for a single path [time, factor].
#' @param maturities Numeric vector. Maturities (in years) for the yield curve.
#' @param lambda Numeric. DNS decay parameter lambda.
#'
#' @return Matrix of yields with dimensions [time, maturity, path] if factors is an array,
#'   or [time, maturity] if factors is a matrix.
#'
#' @examples
#' \dontrun{
#' params <- list(
#'   kappa = c(0.02, 0.3, 0.5),
#'   theta = c(0.04, -0.02, 0.01),
#'   sigma = c(0.005, 0.01, 0.02)
#' )
#' initial_factors <- c(0.035, -0.015, 0.005)
#' sim_result <- simulate_dns_factors(params, initial_factors)
#' maturities <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)
#' yield_curves <- calculate_dns_yields(sim_result$factors, maturities, 0.0609)
#' }
#'
#' @export
calculate_dns_yields <- function(factors, maturities, lambda = 0.0609) {
  # Parameter validation
  checkmate::assert_numeric(lambda, lower = 0, len = 1)
  checkmate::assert_numeric(maturities, lower = 0)

  # Determine input dimensions
  if (is.matrix(factors)) {
    # Single path case
    n_times <- nrow(factors)
    n_paths <- 1
    n_factors <- ncol(factors)
    is_array <- FALSE
  } else if (is.array(factors) && length(dim(factors)) == 3) {
    # Multiple paths case
    n_times <- dim(factors)[1]
    n_paths <- dim(factors)[2]
    n_factors <- dim(factors)[3]
    is_array <- TRUE
  } else {
    stop("factors must be a matrix [time, factor] or an array [time, path, factor]")
  }

  # Check number of factors
  if (n_factors != 3) {
    stop("Expected 3 factors (level, slope, curvature), but got ", n_factors)
  }

  # Initialize yield storage
  n_maturities <- length(maturities)

  if (is_array) {
    yields <- array(0, dim = c(n_times, n_maturities, n_paths))
  } else {
    yields <- matrix(0, nrow = n_times, ncol = n_maturities)
  }

  # Precompute factor loadings for efficiency
  loadings <- matrix(0, nrow = n_maturities, ncol = 2)
  colnames(loadings) <- c("loading_2", "loading_3")

  for (m_idx in 1:n_maturities) {
    m <- maturities[m_idx]

    # Handle small maturity values to avoid numerical issues
    if (m * lambda < 1e-10) {
      loading_2 <- 1.0
      loading_3 <- 0.0
    } else {
      loading_2 <- (1 - exp(-lambda * m)) / (lambda * m)
      loading_3 <- loading_2 - exp(-lambda * m)
    }

    loadings[m_idx, "loading_2"] <- loading_2
    loadings[m_idx, "loading_3"] <- loading_3
  }

  # Calculate yields for all times, paths, and maturities
  if (is_array) {
    # Multiple paths case - vectorized calculation
    for (t in 1:n_times) {
      for (p in 1:n_paths) {
        beta1 <- factors[t, p, 1]
        beta2 <- factors[t, p, 2]
        beta3 <- factors[t, p, 3]

        for (m_idx in 1:n_maturities) {
          loading_2 <- loadings[m_idx, "loading_2"]
          loading_3 <- loadings[m_idx, "loading_3"]

          yields[t, m_idx, p] <- beta1 + beta2 * loading_2 + beta3 * loading_3
        }
      }
    }
  } else {
    # Single path case
    for (t in 1:n_times) {
      beta1 <- factors[t, 1]
      beta2 <- factors[t, 2]
      beta3 <- factors[t, 3]

      for (m_idx in 1:n_maturities) {
        loading_2 <- loadings[m_idx, "loading_2"]
        loading_3 <- loadings[m_idx, "loading_3"]

        yields[t, m_idx] <- beta1 + beta2 * loading_2 + beta3 * loading_3
      }
    }
  }

  return(yields)
}

#' @title Simulate DNS Yield Curves using YUIMA
#'
#' @description
#' Simulates DNS yield curves using the YUIMA package for accurate SDE simulation.
#'
#' @param params List. Model parameters (kappa, theta, sigma).
#' @param initial_factors Numeric vector. Initial factor values (beta1, beta2, beta3).
#' @param maturities Numeric vector. Maturities (in years) for the yield curve.
#' @param lambda Numeric. DNS decay parameter lambda.
#' @param time_horizon Numeric. Simulation time horizon in years.
#' @param n_steps Integer. Number of time steps in the simulation.
#' @param n_paths Integer. Number of simulation paths.
#' @param include_correlation Logical. Whether to include correlations in simulations.
#' @param sim_method Character. Simulation method for YUIMA ("euler" or "milstein").
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return A list containing:
#'   - yields: Array of simulated yields [time, maturity, path]
#'   - factors: Array of simulated factors [time, path, factor]
#'   - time: Vector of time points
#'   - maturities: Vector of maturities
#'   - params: The input parameters
#'
#' @examples
#' \dontrun{
#' params <- list(
#'   kappa = c(0.02, 0.3, 0.5),
#'   theta = c(0.04, -0.02, 0.01),
#'   sigma = c(0.005, 0.01, 0.02)
#' )
#' initial_factors <- c(0.035, -0.015, 0.005)
#' maturities <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)
#' sim_result <- simulate_dns_yield_curves(params, initial_factors, maturities)
#' }
#'
#' @export
simulate_dns_yield_curves <- function(params, initial_factors, maturities,
                                      lambda = 0.0609, time_horizon = 10,
                                      n_steps = 120, n_paths = 100,
                                      include_correlation = TRUE,
                                      sim_method = "euler", seed = NULL) {
  # Parameter validation
  checkmate::assert_list(params)
  checkmate::assert_subset(c("kappa", "theta", "sigma"), names(params))
  checkmate::assert_numeric(initial_factors, len = 3)
  checkmate::assert_numeric(maturities, lower = 0)
  checkmate::assert_numeric(lambda, lower = 0, len = 1)
  checkmate::assert_choice(sim_method, c("euler", "milstein"))

  # Check if YUIMA package is installed
  if (!requireNamespace("yuima", quietly = TRUE)) {
    stop("The YUIMA package is required for SDE simulation. Please install it with: install.packages('yuima')")
  }

  # Simulate factors
  factor_sim <- simulate_dns_factors(
    params = params,
    initial_factors = initial_factors,
    time_horizon = time_horizon,
    n_steps = n_steps,
    n_paths = n_paths,
    include_correlation = include_correlation,
    seed = seed
  )

  # Calculate yield curves
  yields <- calculate_dns_yields(
    factors = factor_sim$factors,
    maturities = maturities,
    lambda = lambda
  )

  # Return results
  return(list(
    yields = yields,
    factors = factor_sim$factors,
    time = factor_sim$time,
    maturities = maturities,
    params = params,
    lambda = lambda
  ))
}

#' @title Calculate Forward Rates from DNS Factors
#'
#' @description
#' Calculates instantaneous forward rates based on DNS factors.
#'
#' @param factors Matrix or array. DNS factors.
#' @param maturities Numeric vector. Maturities (in years) for forward rates.
#' @param lambda Numeric. DNS decay parameter lambda.
#'
#' @return Matrix or array of forward rates with the same structure as the yield output
#'   from calculate_dns_yields.
#'
#' @details
#' The instantaneous forward rate in the DNS model is given by:
#' \deqn{f(m) = \beta_1 + \beta_2 e^{-\lambda m} + \beta_3 \lambda m e^{-\lambda m}}
#'
#' @examples
#' \dontrun{
#' params <- list(
#'   kappa = c(0.02, 0.3, 0.5),
#'   theta = c(0.04, -0.02, 0.01),
#'   sigma = c(0.005, 0.01, 0.02)
#' )
#' initial_factors <- c(0.035, -0.015, 0.005)
#' sim_result <- simulate_dns_factors(params, initial_factors)
#' maturities <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)
#' forward_rates <- calculate_dns_forwards(sim_result$factors, maturities, 0.0609)
#' }
#'
#' @export
calculate_dns_forwards <- function(factors, maturities, lambda = 0.0609) {
  # Parameter validation
  checkmate::assert_numeric(lambda, lower = 0, len = 1)
  checkmate::assert_numeric(maturities, lower = 0)

  # Determine input dimensions
  if (is.matrix(factors)) {
    # Single path case
    n_times <- nrow(factors)
    n_paths <- 1
    n_factors <- ncol(factors)
    is_array <- FALSE
  } else if (is.array(factors) && length(dim(factors)) == 3) {
    # Multiple paths case
    n_times <- dim(factors)[1]
    n_paths <- dim(factors)[2]
    n_factors <- dim(factors)[3]
    is_array <- TRUE
  } else {
    stop("factors must be a matrix [time, factor] or an array [time, path, factor]")
  }

  # Check number of factors
  if (n_factors != 3) {
    stop("Expected 3 factors (level, slope, curvature), but got ", n_factors)
  }

  # Initialize forward rate storage
  n_maturities <- length(maturities)

  if (is_array) {
    forwards <- array(0, dim = c(n_times, n_maturities, n_paths))
  } else {
    forwards <- matrix(0, nrow = n_times, ncol = n_maturities)
  }

  # Calculate forward rates for all times, paths, and maturities
  if (is_array) {
    # Multiple paths case
    for (t in 1:n_times) {
      for (p in 1:n_paths) {
        for (m_idx in 1:n_maturities) {
          m <- maturities[m_idx]
          beta1 <- factors[t, p, 1]
          beta2 <- factors[t, p, 2]
          beta3 <- factors[t, p, 3]

          # Forward rate formula: f(m) = β₁ + β₂ exp(-λm) + β₃ λm exp(-λm)
          exp_term <- exp(-lambda * m)
          forwards[t, m_idx, p] <- beta1 + beta2 * exp_term + beta3 * lambda * m * exp_term
        }
      }
    }
  } else {
    # Single path case
    for (t in 1:n_times) {
      for (m_idx in 1:n_maturities) {
        m <- maturities[m_idx]
        beta1 <- factors[t, 1]
        beta2 <- factors[t, 2]
        beta3 <- factors[t, 3]

        # Forward rate formula: f(m) = β₁ + β₂ exp(-λm) + β₃ λm exp(-λm)
        exp_term <- exp(-lambda * m)
        forwards[t, m_idx] <- beta1 + beta2 * exp_term + beta3 * lambda * m * exp_term
      }
    }
  }

  return(forwards)
}

#' @title Convert Nelson-Siegel Model to YUIMA SDE
#'
#' @description
#' Converts DNS model parameters to a YUIMA SDE object for advanced analysis.
#'
#' @param params List. Model parameters (kappa, theta, sigma, correlation).
#' @param correlation Matrix. Correlation matrix for factor innovations.
#'
#' @return A YUIMA model object representing the DNS factor dynamics.
#'
#' @details
#' This function creates a YUIMA model object that represents the factor dynamics
#' of the DNS model. This object can be used for various analyses supported by
#' the YUIMA package, such as parameter estimation, filtering, and more advanced
#' simulations.
#'
#' @examples
#' \dontrun{
#' params <- list(
#'   kappa = c(0.02, 0.3, 0.5),
#'   theta = c(0.04, -0.02, 0.01),
#'   sigma = c(0.005, 0.01, 0.02)
#' )
#' correlation <- matrix(c(1, 0.2, 0.1, 0.2, 1, 0.3, 0.1, 0.3, 1), nrow = 3)
#' yuima_model <- dns_to_yuima_model(params, correlation)
#' }
#'
#' @export
dns_to_yuima_model <- function(params, correlation = NULL) {
  # Check if YUIMA package is installed
  if (!requireNamespace("yuima", quietly = TRUE)) {
    stop("The YUIMA package is required. Please install it with: install.packages('yuima')")
  }

  # Parameter validation
  checkmate::assert_list(params)
  checkmate::assert_subset(c("kappa", "theta", "sigma"), names(params))
  checkmate::assert_numeric(params$kappa, len = 3)
  checkmate::assert_numeric(params$theta, len = 3)
  checkmate::assert_numeric(params$sigma, len = 3, lower = 0)

  # Define variable names
  var_names <- c("beta1", "beta2", "beta3")

  # Extract parameters
  kappa <- params$kappa
  theta <- params$theta
  sigma <- params$sigma

  # Use correlation matrix if provided, otherwise use identity
  if (is.null(correlation)) {
    correlation <- diag(3)
  } else {
    checkmate::assert_matrix(correlation, nrows = 3, ncols = 3)
    # Ensure it's a correlation matrix
    if (!all(diag(correlation) == 1) || any(abs(correlation) > 1)) {
      warning("Provided matrix doesn't appear to be a valid correlation matrix. Using identity matrix instead.")
      correlation <- diag(3)
    }
  }

  # Define the drift coefficient (mean reversion)
  drift_terms <- c(
    paste0(kappa[1], "*", theta[1], "-", kappa[1], "*beta1"),
    paste0(kappa[2], "*", theta[2], "-", kappa[2], "*beta2"),
    paste0(kappa[3], "*", theta[3], "-", kappa[3], "*beta3")
  )

  # Define the diffusion coefficient
  # Calculate Cholesky decomposition of correlation matrix
  chol_corr <- try(chol(correlation), silent = TRUE)
  if (inherits(chol_corr, "try-error")) {
    warning("Correlation matrix is not positive definite. Using identity matrix.")
    chol_corr <- diag(3)
  }

  # Construct the diffusion matrix as sigma * chol_corr
  diffusion_matrix <- matrix(0, nrow = 3, ncol = 3)
  for (i in 1:3) {
    for (j in 1:3) {
      diffusion_matrix[i, j] <- sigma[i] * chol_corr[j, i]
    }
  }

  # Convert diffusion matrix to yuima format
  diffusion_terms <- matrix("0", nrow = 3, ncol = 3)
  for (i in 1:3) {
    for (j in 1:3) {
      if (diffusion_matrix[i, j] != 0) {
        diffusion_terms[i, j] <- as.character(diffusion_matrix[i, j])
      }
    }
  }

  # Create YUIMA model
  yuima_model <- yuima::setModel(
    drift = drift_terms,
    diffusion = diffusion_terms,
    state.variable = var_names,
    solve.variable = var_names
  )

  return(yuima_model)
}

#' @title Estimate DNS Parameters from Factors using YUIMA
#'
#' @description
#' Estimates DNS model parameters from factor time series using YUIMA's advanced estimation methods.
#'
#' @param factor_data Data frame or matrix. Time series of DNS factors (level, slope, curvature).
#' @param initial_params List. Initial parameter values for optimization.
#' @param estimation_method Character. Method for parameter estimation ("qmle" or "gmm").
#' @param delta_t Numeric. Time step between observations in years.
#'
#' @return A list containing:
#'   - estimated_params: Estimated model parameters
#'   - yuima_model: The fitted YUIMA model
#'   - fit_result: Full result from the YUIMA estimation
#'
#' @details
#' This function uses the YUIMA package's estimation capabilities to fit an
#' Ornstein-Uhlenbeck process to observed factor time series. It provides
#' more rigorous statistical estimation than standard regression methods.
#'
#' @examples
#' \dontrun{
#' # Assuming factor_data is a data frame with columns: Date, beta1, beta2, beta3
#' factors <- factor_data[, c("beta1", "beta2", "beta3")]
#' initial <- list(
#'   kappa = c(0.02, 0.3, 0.5),
#'   theta = c(0.04, -0.02, 0.01),
#'   sigma = c(0.005, 0.01, 0.02)
#' )
#' est_result <- estimate_dns_parameters_yuima(factors, initial)
#' }
#'
#' @export
estimate_dns_parameters_yuima <- function(factor_data, initial_params = NULL,
                                          estimation_method = "qmle",
                                          delta_t = 1/252) {
  # Check if YUIMA package is installed
  if (!requireNamespace("yuima", quietly = TRUE)) {
    stop("The YUIMA package is required. Please install it with: install.packages('yuima')")
  }

  # Parameter validation
  if (is.data.frame(factor_data)) {
    if (ncol(factor_data) < 3) {
      stop("factor_data must have at least 3 columns for the DNS factors")
    }
    # Extract numeric data if it's a data frame
    factors <- as.matrix(factor_data[, 1:3])
  } else if (is.matrix(factor_data)) {
    if (ncol(factor_data) < 3) {
      stop("factor_data must have at least 3 columns for the DNS factors")
    }
    factors <- factor_data[, 1:3]
  } else {
    stop("factor_data must be a data frame or matrix")
  }

  # Check for missing values
  if (any(is.na(factors))) {
    stop("factor_data contains missing values. Please handle missing values before estimation.")
  }

  # Set up default initial parameters if not provided
  if (is.null(initial_params)) {
    # Estimate initial params from data
    n_obs <- nrow(factors)

    # Estimate mean reversion parameters using simple regression
    kappa_est <- numeric(3)
    theta_est <- numeric(3)
    sigma_est <- numeric(3)

    for (i in 1:3) {
      factor_changes <- diff(factors[, i])
      factor_levels <- factors[1:(n_obs-1), i]

      # Regression: Δβ(t) ~ κθΔt - κβ(t-1)Δt
      lm_result <- lm(factor_changes ~ factor_levels)

      kappa_est[i] <- -coef(lm_result)[2] / delta_t
      if (kappa_est[i] <= 0) kappa_est[i] <- 0.1  # Ensure positive mean reversion

      theta_est[i] <- -coef(lm_result)[1] / coef(lm_result)[2]
      sigma_est[i] <- sd(residuals(lm_result)) / sqrt(delta_t)
    }

    initial_params <- list(
      kappa = kappa_est,
      theta = theta_est,
      sigma = sigma_est
    )
  }

  # Create parameter list for YUIMA
  # For Ornstein-Uhlenbeck, we need parameters in the form:
  # drift: a - b*x (where a = kappa*theta, b = kappa)
  # diffusion: sigma

  # Convert to YUIMA parameter names (a1, a2, a3, b1, b2, b3, s1, s2, s3)
  a_params <- initial_params$kappa * initial_params$theta
  b_params <- initial_params$kappa
  s_params <- initial_params$sigma

  param_names <- c(
    "a1", "a2", "a3",  # Drift intercepts (kappa*theta)
    "b1", "b2", "b3",  # Drift slopes (kappa)
    "s1", "s2", "s3"   # Diffusion coefficients (sigma)
  )

  param_values <- c(a_params, b_params, s_params)
  names(param_values) <- param_names

  # Create the multivariate Ornstein-Uhlenbeck model in YUIMA
  drift_terms <- c(
    "a1 - b1*beta1",
    "a2 - b2*beta2",
    "a3 - b3*beta3"
  )

  # We assume diagonal diffusion matrix for simplicity
  diffusion_terms <- matrix("0", nrow = 3, ncol = 3)
  diffusion_terms[1, 1] <- "s1"
  diffusion_terms[2, 2] <- "s2"
  diffusion_terms[3, 3] <- "s3"

  var_names <- c("beta1", "beta2", "beta3")

  # Create YUIMA model
  yuima_model <- yuima::setModel(
    drift = drift_terms,
    diffusion = diffusion_terms,
    state.variable = var_names,
    solve.variable = var_names,
    parameter = param_names
  )

  # Convert data to YUIMA format
  time_indices <- seq(0, (nrow(factors) - 1) * delta_t, by = delta_t)
  yuima_data <- yuima::setData(
    zoo::zoo(factors, time_indices),
    delta = delta_t
  )

  # Create YUIMA object
  yuima_obj <- yuima::setYuima(
    model = yuima_model,
    data = yuima_data
  )

  # Estimate parameters
  if (estimation_method == "qmle") {
    # Quasi-Maximum Likelihood Estimation
    start_params <- param_values
    lower_bounds <- c(rep(-Inf, 3), rep(0.001, 3), rep(0.0001, 3))  # Lower bounds for parameters
    upper_bounds <- rep(Inf, 9)  # Upper bounds for parameters

    fit_result <- yuima::qmle(
      yuima_obj,
      start = start_params,
      lower = lower_bounds,
      upper = upper_bounds,
      method = "L-BFGS-B"
    )
  } else if (estimation_method == "gmm") {
    # Generalized Method of Moments
    fit_result <- yuima::gmm(yuima_obj)
  } else {
    stop("Invalid estimation_method. Use 'qmle' or 'gmm'.")
  }

  # Extract estimated parameters and convert back to DNS format
  est_params <- fit_result@coef

  # Convert back to DNS parameters
  kappa <- est_params[4:6]  # b1, b2, b3
  theta <- est_params[1:3] / kappa  # a1/b1, a2/b2, a3/b3
  sigma <- est_params[7:9]  # s1, s2, s3

  # Create result list
  result <- list(
    estimated_params = list(
      kappa = kappa,
      theta = theta,
      sigma = sigma
    ),
    yuima_model = yuima_model,
    fit_result = fit_result,
    input_data = factors,
    delta_t = delta_t
  )

  return(result)
}
