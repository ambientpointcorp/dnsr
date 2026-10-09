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
#' @param sim_method Character. Discretisation scheme passed to YUIMA, either
#'   \code{"euler"} (the default) or \code{"milstein"}. Ignored by
#'   \code{engine = "exact"}, which does not discretise.
#' @param engine Character. \code{"yuima"} (the default) simulates through the
#'   YUIMA package, one call per path. \code{"exact"} samples the exact
#'   Gaussian transition law of the factor system directly: it needs no extra
#'   package, carries no discretisation error at any step size, and is roughly
#'   two orders of magnitude faster. The default is expected to change to
#'   \code{"exact"} in a future release.
#'
#' @return A list containing:
#'   - factors: Array of simulated factors (time, path, factor)
#'   - time: Vector of time points
#'   - params: The input parameters
#'   - sim_method: The scheme used
#'   - engine: The engine used
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
#' The diffusion coefficient of this model does not depend on the state, so the
#' Milstein correction term vanishes and \code{sim_method = "milstein"} coincides
#' with the Euler-Maruyama scheme. The argument is honoured for completeness;
#' accuracy is governed by \code{n_steps}.
#'
#' That accuracy ceiling is what \code{engine = "exact"} removes. The factor
#' system is a linear SDE with state-independent diffusion, so its transition
#' density over any step is Gaussian with a mean and covariance available in
#' closed form. Sampling that law directly is exact in distribution at every
#' step size, and replaces one YUIMA call per path with vectorised draws over
#' all paths at once. The two engines agree in distribution; they do not
#' produce the same paths from the same seed.
#'
#' @examples
#' params <- list(
#'   kappa = c(0.02, 0.3, 0.5),
#'   theta = c(0.04, -0.02, 0.01),
#'   sigma = c(0.005, 0.01, 0.02),
#'   correlation = diag(3)
#' )
#' initial_factors <- c(0.035, -0.015, 0.005)
#'
#' # The exact engine needs no extra package
#' sim_result <- simulate_dns_factors(params, initial_factors, n_paths = 50,
#'                                    engine = "exact", seed = 1)
#' sim_result
#'
#' \dontrun{
#' # The YUIMA engine, for comparison
#' simulate_dns_factors(params, initial_factors, n_paths = 50, engine = "yuima")
#' }
#'
#' @export
simulate_dns_factors <- function(params, initial_factors, time_horizon = 10,
                                 n_steps = 120, n_paths = 100,
                                 include_correlation = TRUE, seed = NULL,
                                 sim_method = "euler",
                                 engine = c("yuima", "exact")) {
  engine <- match.arg(engine)

  # YUIMA is only needed by the engine that uses it
  if (engine == "yuima" && !requireNamespace("yuima", quietly = TRUE)) {
    stop("The YUIMA package is required for engine = \"yuima\". Install it ",
         "with install.packages('yuima'), or pass engine = \"exact\", which ",
         "needs no extra package and samples the exact transition law.")
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
  checkmate::assert_choice(sim_method, c("euler", "milstein"))

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

  # The exact engine needs nothing beyond the parameters resolved above.
  if (engine == "exact") {
    factors <- simulate_ou_exact(
      kappa = kappa, theta = theta, sigma = sigma, corr_matrix = corr_matrix,
      x0 = initial_factors, dt = delta_t, n_steps = n_steps, n_paths = n_paths
    )

    result <- list(
      factors = factors,
      time = time_grid,
      params = params,
      sim_method = sim_method,
      engine = engine
    )
    class(result) <- c("dns_simulation", "list")
    return(result)
  }

  # Convert parameters to YUIMA format
  # 1. Drift coefficient (mean reversion)
  # For O-U process: dX_t = kappa*(theta - X_t)dt + sigma*dW_t
  # The drift term is: kappa*theta - kappa*X_t
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
      method = sim_method
    )

    sim_data <- sim@data@original.data

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
  result <- list(
    factors = factors,
    time = time_grid,
    params = params,
    sim_method = sim_method,
    engine = engine
  )
  class(result) <- c("dns_simulation", "list")

  return(result)
}

#' @title Calculate Yield Curves from DNS Factors
#'
#' @description
#' Calculates yield curves based on DNS factors for specified maturities.
#'
#' @param factors Array or matrix. If an array, simulated factors indexed by
#'   (time, path, factor); if a matrix, factors for a single path indexed by
#'   (time, factor).
#' @param maturities Numeric vector. Maturities (in years) for the yield curve.
#' @param lambda Numeric. DNS decay parameter lambda, per unit of maturity.
#'   Maturities here are in years, and the default places the curvature hump at
#'   3.5 years, near the middle of a Treasury curve. It is written as
#'   \code{ns_lambda(3.5)} so the number states its own meaning; see
#'   \code{\link{ns_lambda}}. Diebold and Li's 0.0609 is a month-denominated
#'   constant and in years would place that hump at 29.4 years.
#'
#' @return Array of yields indexed by (time, maturity, path) if \code{factors} is
#'   an array, or a matrix indexed by (time, maturity) if it is a matrix.
#'
#' @examples
#' maturities <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)
#'
#' # A single curve: one row of (level, slope, curvature)
#' factors <- matrix(c(0.035, -0.015, 0.005), nrow = 1)
#' round(calculate_dns_yields(factors, maturities, ns_lambda(3.5))[1, ], 5)
#'
#' # Three paths over four dates, as returned by simulate_dns_factors()
#' paths <- array(rep(c(0.035, -0.015, 0.005), each = 12), dim = c(4, 3, 3))
#' dim(calculate_dns_yields(paths, maturities, ns_lambda(3.5)))
#'
#' @export
calculate_dns_yields <- function(factors, maturities,
                                 lambda = ns_lambda(3.5)) {
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
  loadings <- ns_loadings(maturities, lambda)

  # Calculate yields for all times, paths, and maturities. The loop runs over
  # maturities only; the time and path dimensions are handled by vectorised
  # arithmetic. The expression and its left-to-right operation order are
  # unchanged from the scalar form, so results are identical to the last bit.
  if (is_array) {
    # matrix() keeps each factor two-dimensional when n_times or n_paths is 1
    beta1 <- matrix(factors[, , 1], nrow = n_times, ncol = n_paths)
    beta2 <- matrix(factors[, , 2], nrow = n_times, ncol = n_paths)
    beta3 <- matrix(factors[, , 3], nrow = n_times, ncol = n_paths)

    for (m_idx in seq_len(n_maturities)) {
      loading_2 <- loadings[m_idx, 1]
      loading_3 <- loadings[m_idx, 2]

      yields[, m_idx, ] <- beta1 + beta2 * loading_2 + beta3 * loading_3
    }
  } else {
    # Single path case
    beta1 <- factors[, 1]
    beta2 <- factors[, 2]
    beta3 <- factors[, 3]

    for (m_idx in seq_len(n_maturities)) {
      loading_2 <- loadings[m_idx, 1]
      loading_3 <- loadings[m_idx, 2]

      yields[, m_idx] <- beta1 + beta2 * loading_2 + beta3 * loading_3
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
#' @param lambda Numeric. DNS decay parameter lambda, per unit of maturity.
#'   Maturities here are in years, and the default places the curvature hump at
#'   3.5 years, near the middle of a Treasury curve. It is written as
#'   \code{ns_lambda(3.5)} so the number states its own meaning; see
#'   \code{\link{ns_lambda}}. Diebold and Li's 0.0609 is a month-denominated
#'   constant and in years would place that hump at 29.4 years.
#' @param time_horizon Numeric. Simulation time horizon in years.
#' @param n_steps Integer. Number of time steps in the simulation.
#' @param n_paths Integer. Number of simulation paths.
#' @param include_correlation Logical. Whether to include correlations in simulations.
#' @param sim_method Character. Simulation method for YUIMA ("euler" or "milstein").
#' @param seed Integer. Random seed for reproducibility.
#' @param engine Character. \code{"yuima"} (the default) or \code{"exact"};
#'   see \code{\link{simulate_dns_factors}}.
#'
#' @return A list containing:
#'   - yields: Array of simulated yields, indexed by (time, maturity, path)
#'   - factors: Array of simulated factors, indexed by (time, path, factor)
#'   - time: Vector of time points
#'   - maturities: Vector of maturities
#'   - params: The input parameters
#'   - sim_method, engine: The scheme and engine used
#'
#' @examples
#' params <- list(
#'   kappa = c(0.02, 0.3, 0.5),
#'   theta = c(0.04, -0.02, 0.01),
#'   sigma = c(0.005, 0.01, 0.02)
#' )
#' initial_factors <- c(0.035, -0.015, 0.005)
#' maturities <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)
#'
#' sim_result <- simulate_dns_yield_curves(params, initial_factors, maturities,
#'                                         n_paths = 50, engine = "exact",
#'                                         seed = 1)
#' sim_result
#'
#' @export
simulate_dns_yield_curves <- function(params, initial_factors, maturities,
                                      lambda = ns_lambda(3.5), time_horizon = 10,
                                      n_steps = 120, n_paths = 100,
                                      include_correlation = TRUE,
                                      sim_method = "euler", seed = NULL,
                                      engine = c("yuima", "exact")) {
  # Parameter validation
  engine <- match.arg(engine)
  checkmate::assert_list(params)
  checkmate::assert_subset(c("kappa", "theta", "sigma"), names(params))
  checkmate::assert_numeric(initial_factors, len = 3)
  checkmate::assert_numeric(maturities, lower = 0)
  checkmate::assert_numeric(lambda, lower = 0, len = 1)
  checkmate::assert_choice(sim_method, c("euler", "milstein"))

  # Simulate factors (simulate_dns_factors() checks for YUIMA if it needs it)
  factor_sim <- simulate_dns_factors(
    params = params,
    initial_factors = initial_factors,
    time_horizon = time_horizon,
    n_steps = n_steps,
    n_paths = n_paths,
    include_correlation = include_correlation,
    seed = seed,
    sim_method = sim_method,
    engine = engine
  )

  # Calculate yield curves
  yields <- calculate_dns_yields(
    factors = factor_sim$factors,
    maturities = maturities,
    lambda = lambda
  )

  # Return results
  result <- list(
    yields = yields,
    factors = factor_sim$factors,
    time = factor_sim$time,
    maturities = maturities,
    params = params,
    lambda = lambda,
    sim_method = sim_method,
    engine = engine
  )
  class(result) <- c("dns_simulation", "list")

  return(result)
}

#' @title Calculate Forward Rates from DNS Factors
#'
#' @description
#' Calculates instantaneous forward rates based on DNS factors.
#'
#' @param factors Matrix or array. DNS factors.
#' @param maturities Numeric vector. Maturities (in years) for forward rates.
#' @param lambda Numeric. DNS decay parameter lambda, per unit of maturity.
#'   Maturities here are in years, and the default places the curvature hump at
#'   3.5 years, near the middle of a Treasury curve. It is written as
#'   \code{ns_lambda(3.5)} so the number states its own meaning; see
#'   \code{\link{ns_lambda}}. Diebold and Li's 0.0609 is a month-denominated
#'   constant and in years would place that hump at 29.4 years.
#'
#' @return Matrix or array of forward rates with the same structure as the yield output
#'   from calculate_dns_yields.
#'
#' @details
#' The instantaneous forward rate in the DNS model is given by:
#' \deqn{f(m) = \beta_1 + \beta_2 e^{-\lambda m} + \beta_3 \lambda m e^{-\lambda m}}
#'
#' @examples
#' maturities <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)
#' factors <- matrix(c(0.035, -0.015, 0.005), nrow = 1)
#'
#' # Forwards sit above spots where the curve is rising
#' rbind(
#'   spot = calculate_dns_yields(factors, maturities, ns_lambda(3.5))[1, ],
#'   forward = calculate_dns_forwards(factors, maturities, ns_lambda(3.5))[1, ]
#' )
#'
#' @export
calculate_dns_forwards <- function(factors, maturities,
                                   lambda = ns_lambda(3.5)) {
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

  # Calculate forward rates for all times, paths, and maturities. As in
  # calculate_dns_yields(), only the maturity dimension is looped; the
  # arithmetic expression is unchanged, including its grouping, so results are
  # identical to the last bit.
  #
  # Forward rate formula: f(m) = b1 + b2 exp(-lambda m) + b3 lambda m exp(-lambda m)
  if (is_array) {
    beta1 <- matrix(factors[, , 1], nrow = n_times, ncol = n_paths)
    beta2 <- matrix(factors[, , 2], nrow = n_times, ncol = n_paths)
    beta3 <- matrix(factors[, , 3], nrow = n_times, ncol = n_paths)

    for (m_idx in seq_len(n_maturities)) {
      m <- maturities[m_idx]
      exp_term <- exp(-lambda * m)
      forwards[, m_idx, ] <- beta1 + beta2 * exp_term + beta3 * lambda * m * exp_term
    }
  } else {
    # Single path case
    beta1 <- factors[, 1]
    beta2 <- factors[, 2]
    beta3 <- factors[, 3]

    for (m_idx in seq_len(n_maturities)) {
      m <- maturities[m_idx]
      exp_term <- exp(-lambda * m)
      forwards[, m_idx] <- beta1 + beta2 * exp_term + beta3 * lambda * m * exp_term
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
    if (!isTRUE(all.equal(diag(correlation), rep(1, 3))) ||
        any(abs(correlation) > 1 + 1e-8)) {
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
#' @param factor_data Data frame or matrix. Time series of DNS factors (level,
#'   slope, curvature). Columns named \code{beta1}, \code{beta2} and
#'   \code{beta3} are used when present, so a \code{parameters} frame from
#'   \code{\link{fit_dns_model}} can be passed directly; otherwise the first
#'   three numeric columns are taken, after any \code{Date} column is dropped.
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
  checkmate::assert_choice(estimation_method, c("qmle", "gmm"))
  checkmate::assert_numeric(delta_t, lower = 0, len = 1, any.missing = FALSE)

  if (!is.data.frame(factor_data) && !is.matrix(factor_data)) {
    stop("factor_data must be a data frame or matrix")
  }

  # Select the three factor columns by name when they are present: a
  # `parameters` frame from fit_dns_model() leads with Date, and taking columns
  # 1:3 positionally would feed dates into the regression.
  beta_names <- c("beta1", "beta2", "beta3")
  cols <- colnames(factor_data)

  if (!is.null(cols) && all(beta_names %in% cols)) {
    factors <- as.matrix(as.data.frame(factor_data)[, beta_names, drop = FALSE])
  } else {
    candidate <- as.data.frame(factor_data, stringsAsFactors = FALSE)
    candidate <- candidate[, !names(candidate) %in% "Date", drop = FALSE]
    numeric_cols <- names(candidate)[vapply(candidate, is.numeric, logical(1))]
    if (length(numeric_cols) < 3) {
      stop("factor_data must have at least 3 numeric columns for the DNS ",
           "factors, or columns named beta1, beta2 and beta3")
    }
    factors <- as.matrix(candidate[, numeric_cols[1:3], drop = FALSE])
  }
  colnames(factors) <- beta_names

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

      # Regression: d_beta(t) ~ kappa*theta*dt - kappa*beta(t-1)*dt
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

  # Create YUIMA model (parameter names are inferred from the drift and
  # diffusion expressions; setModel() has no `parameter` argument)
  yuima_model <- yuima::setModel(
    drift = drift_terms,
    diffusion = diffusion_terms,
    state.variable = var_names,
    solve.variable = var_names
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
    # Quasi-Maximum Likelihood Estimation. yuima::qmle requires `start`
    # (and any bounds) as named lists keyed by parameter name.
    start_params <- as.list(param_values)
    lower_bounds <- stats::setNames(
      as.list(c(rep(-Inf, 3), rep(0.001, 3), rep(0.0001, 3))), param_names)
    upper_bounds <- stats::setNames(as.list(rep(Inf, 9)), param_names)

    fit_result <- yuima::qmle(
      yuima_obj,
      start = start_params,
      lower = lower_bounds,
      upper = upper_bounds,
      method = "L-BFGS-B"
    )
  } else {
    # Generalized Method of Moments
    fit_result <- yuima::gmm(yuima_obj)
  }

  # Extract estimated parameters and convert back to DNS format.
  # Index by name: yuima's coefficient ordering follows the model's
  # internal parameter order, not the order the names were declared in.
  est_params <- fit_result@coef

  a_hat <- est_params[c("a1", "a2", "a3")]
  kappa <- stats::setNames(est_params[c("b1", "b2", "b3")],
                           c("kappa1", "kappa2", "kappa3"))
  theta <- stats::setNames(a_hat / as.numeric(kappa),
                           c("theta1", "theta2", "theta3"))
  sigma <- stats::setNames(est_params[c("s1", "s2", "s3")],
                           c("sigma1", "sigma2", "sigma3"))

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


# Exact simulation of the DNS factor system.
#
# Each factor follows dX_i = kappa_i (theta_i - X_i) dt + sigma_i dW_i, with
# corr(dW_i, dW_j) = rho_ij. The system is linear and its diffusion does not
# depend on the state, so the transition law over a step of length dt is
# Gaussian and known in closed form:
#
#   E[X_i(t+dt) | X(t)] = theta_i + (X_i(t) - theta_i) exp(-kappa_i dt)
#   Cov_ij              = rho_ij sigma_i sigma_j *
#                         (1 - exp(-(kappa_i + kappa_j) dt)) / (kappa_i + kappa_j)
#
# The covariance couples the two reversion speeds. Using each factor's own
# step variance together with a Cholesky factor of the correlation matrix is
# only correct when every kappa is equal, which is why the full matrix is
# built here. Sampling this law carries no discretisation error at any step
# size, unlike an Euler scheme.
simulate_ou_exact <- function(kappa, theta, sigma, corr_matrix, x0,
                              dt, n_steps, n_paths) {
  decay <- exp(-kappa * dt)

  kappa_sum <- outer(kappa, kappa, "+")
  # (1 - exp(-s dt)) / s, taking its limit dt as s -> 0
  step_scale <- ifelse(
    abs(kappa_sum) > 1e-10,
    (1 - exp(-kappa_sum * dt)) / kappa_sum,
    dt
  )
  cov_step <- corr_matrix * outer(sigma, sigma) * step_scale

  chol_step <- tryCatch(chol(cov_step), error = function(e) NULL)
  if (is.null(chol_step)) {
    warning("Step covariance matrix is not positive definite. ",
            "Simulating the factors without correlation.")
    chol_step <- diag(sqrt(pmax(diag(cov_step), 0)), 3)
  }

  out <- array(0, dim = c(n_steps + 1, n_paths, 3))
  theta_rep <- rep(theta, each = n_paths)
  decay_rep <- rep(decay, each = n_paths)

  state <- matrix(x0, nrow = n_paths, ncol = 3, byrow = TRUE)
  out[1, , ] <- state

  for (step in seq_len(n_steps)) {
    shock <- matrix(rnorm(n_paths * 3), nrow = n_paths, ncol = 3) %*% chol_step
    state <- theta_rep + (state - theta_rep) * decay_rep + shock
    out[step + 1, , ] <- state
  }

  out
}
