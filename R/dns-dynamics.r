#' @title Dynamic Nelson-Siegel Model Dynamics
#' @description Functions to estimate and analyze the dynamics of DNS model factors.

#' @importFrom stats lm coef sd resid na.omit vcov df.residual qnorm
#' @importFrom checkmate assert_data_frame assert_subset assert_class

#' @title Estimate Dynamics of the DNS Factors
#'
#' @description
#' Estimates the stochastic dynamics of the factors in the Dynamic Nelson-Siegel model
#' based on historical factor time series.
#'
#' @param params_df Data frame. The parameters data frame from fit_dns_model() containing
#'   columns for Date, beta1, beta2, beta3. With \code{sampling = "daily"} the
#'   corresponding change columns are also required; with \code{sampling = "monthly"}
#'   they are recomputed internally on the month-end grid and need not be present.
#' @param min_observations Integer. Minimum number of observations required for
#'   estimation, counted on the series the regression actually runs on (month-ends
#'   when \code{sampling = "monthly"}).
#' @param dt Numeric. Time step in years of the series the regression runs on.
#'   Defaults to the sampling mode's spacing: 1/12 for \code{"monthly"}, 1/252 for
#'   \code{"daily"}. With \code{sampling = "monthly"} the step is month-end to
#'   month-end by construction, so a conflicting value is ignored with a warning.
#' @param sampling Character. \code{"monthly"} (default) subsamples the factor
#'   series to its last observation per month before estimating; \code{"daily"}
#'   estimates on the series as given, preserving the pre-0.2.0 behavior.
#' @param structure Character. \code{"diagonal"} (the default) estimates three
#'   independent Ornstein-Uhlenbeck processes with correlated innovations, as
#'   earlier versions did. \code{"full"} regresses each factor change on the
#'   levels of all three factors, giving the full mean-reversion matrix of the
#'   Diebold-Li VAR(1), in which a shock to one factor can pull the others.
#'
#' @return A list containing:
#'   - For each factor (beta1, beta2, beta3):
#'     - kappa: Mean reversion speed
#'     - theta: Long-run mean
#'     - sigma: Volatility
#'     - model: Fitted regression model
#'     - half_life: Half-life of mean reversion (in days)
#'     - summary: Basic statistics about the factor
#'     - se: Standard errors of kappa, theta and sigma
#'     - t_value: t statistics for kappa and theta against zero
#'     - half_life_ci: 95 percent confidence interval for the half-life, in years
#'   - correlation: Correlation matrix of factor innovations
#'   - innovations: Matrix of factor innovations (residuals from mean reversion)
#'   - innovation_dates: Date vector aligned with the rows of innovations
#'   - estimation_period: Date range used for estimation
#'   - sampling, dt: The sampling mode and time step actually used
#'   - structure: The mean-reversion structure estimated
#'   - kappa_matrix: The 3 by 3 mean-reversion matrix. Diagonal under
#'     \code{structure = "diagonal"}; the estimated VAR(1) matrix under
#'     \code{"full"}
#'   - theta_vector: Long-run means as a vector, named by factor
#'   - eigenvalues, system_stability: Eigenvalues of \code{kappa_matrix} and
#'     whether every real part is positive
#'
#' @details
#' This function implements the estimation of an Ornstein-Uhlenbeck process for each factor:
#' \deqn{d\beta(t) = \kappa(\theta - \beta(t))dt + \sigma dW(t)}
#' The discretized version used for estimation is:
#' \deqn{\Delta\beta(t) = \kappa(\theta - \beta(t-1))\Delta t + \sigma\epsilon(t)}
#'
#' The default month-end sampling follows Diebold-Li practice. Daily
#' cross-sectional refits carry measurement noise that is serially
#' uncorrelated; estimating the OU regression on daily changes reads that
#' noise as mean reversion and inflates kappa by an order of magnitude
#' (errors-in-variables), which in turn makes fixed-step simulation of the
#' resulting dynamics numerically unstable. Month-end sampling averages the
#' noise down and recovers reversion speeds on the economically meaningful
#' time scale. Subsampling is idempotent: passing an already month-end
#' series through \code{sampling = "monthly"} reproduces it unchanged.
#'
#' Standard errors come from the fitted regression. For kappa the mapping from
#' the slope is linear, so \eqn{se(\kappa) = se(b) / \Delta t}. For theta,
#' which is a ratio of coefficients, they come from the delta method. The one
#' for sigma is the asymptotic \eqn{\sigma / \sqrt{2(n-k)}}. Note that kappa
#' carries a well-known upward small-sample bias of order \eqn{5 / T} in a
#' sample of \eqn{T} years, which these standard errors do not capture: they
#' describe the sampling noise, not the bias.
#'
#' With \code{structure = "full"} the system estimated is
#' \deqn{\Delta\beta(t) = K(\theta - \beta(t-1))\Delta t + \Sigma\epsilon(t),}
#' so each change is regressed on all three lagged levels. Writing the fitted
#' intercepts as \eqn{a} and slopes as \eqn{B}, the mean-reversion matrix is
#' \eqn{K = -B / \Delta t} and the long-run mean is
#' \eqn{\theta = K^{-1} a / \Delta t}. The per-factor \code{kappa} entries
#' then hold the diagonal of \eqn{K}; the cross terms are only in
#' \code{kappa_matrix}, and \code{\link{simulate_dns_factors}} uses the
#' diagonal.
#'
#' @examples
#' data(treasury_yields)
#'
#' # The default monthly sampling estimates on month-end observations, so
#' # fitting only those is enough and is an order of magnitude quicker.
#' # warm_start = FALSE makes each cross-section independent of which other
#' # dates are in the panel, so these fits are identical to the ones fitting
#' # every day would produce on the same month-ends.
#' month_ends <- treasury_yields[
#'   !duplicated(format(treasury_yields$Date, "%Y-%m"), fromLast = TRUE),
#' ]
#' dns_result <- fit_dns_model(month_ends, n_days = nrow(month_ends),
#'                             show_progress = FALSE, warm_start = FALSE)
#'
#' dynamics <- estimate_dns_dynamics(dns_result$parameters)
#' dynamics
#'
#' # Standard errors show how weakly kappa is identified in a short sample
#' dynamics$beta1$se
#' dynamics$beta1$half_life_ci
#'
#' # The full VAR(1), in which a shock to one factor pulls the others
#' full <- estimate_dns_dynamics(dns_result$parameters, structure = "full")
#' round(full$kappa_matrix, 3)
#'
#' @export
estimate_dns_dynamics <- function(params_df, min_observations = 20, dt = NULL,
                                  sampling = c("monthly", "daily"),
                                  structure = c("diagonal", "full")) {
  # Validate inputs
  sampling <- match.arg(sampling)
  structure <- match.arg(structure)
  checkmate::assert_data_frame(params_df)
  checkmate::assert_subset(c("Date", "beta1", "beta2", "beta3"), names(params_df))
  if (sampling == "daily") {
    checkmate::assert_subset(
      c("beta1_change", "beta2_change", "beta3_change"), names(params_df)
    )
  }
  checkmate::assert_class(params_df$Date, "Date")
  checkmate::assert_numeric(min_observations, lower = 10, len = 1)

  if (sampling == "monthly") {
    if (!is.null(dt) && !isTRUE(all.equal(dt, 1 / 12))) {
      warning("sampling = \"monthly\" estimates on month-end observations; ",
              "ignoring dt = ", format(dt), " and using dt = 1/12.")
    }
    dt <- 1 / 12
    params_df <- params_df[order(params_df$Date), ]
    month_key <- format(params_df$Date, "%Y-%m")
    params_df <- params_df[!duplicated(month_key, fromLast = TRUE), ]
    for (f in c("beta1", "beta2", "beta3")) {
      params_df[[paste0(f, "_change")]] <- c(NA, diff(params_df[[f]]))
    }
  } else if (is.null(dt)) {
    dt <- 1 / 252
  }
  checkmate::assert_numeric(dt, lower = 0, len = 1)

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

  # The OU regression is on the LAGGED level: the model is
  # d_beta(t) = kappa (theta - beta(t-1)) dt + sigma eps(t). Each row of
  # factors_df holds beta(t) and beta(t) - beta(t-1), so beta(t-1) is their
  # difference. Deriving it this way is exact in both sampling modes and makes
  # no assumption that the retained rows are contiguous.
  #
  # Regressing on beta(t) instead would put it on both sides of the equation
  # and return the slope with the opposite sign, since
  # cov(beta(t) - beta(t-1), beta(t)) = var(beta) - gamma(1) > 0.
  lagged <- vapply(
    factor_names,
    function(f) factors_df[[f]] - factors_df[[paste0(f, "_change")]],
    numeric(nrow(factors_df))
  )
  lagged <- matrix(lagged, nrow = nrow(factors_df), ncol = 3,
                   dimnames = list(NULL, factor_names))

  # Under structure = "full" each change is regressed on all three lagged
  # levels, so the regressions are built from a shared design.
  full_design <- NULL
  if (structure == "full") {
    full_design <- data.frame(
      beta1_level = lagged[, "beta1"],
      beta2_level = lagged[, "beta2"],
      beta3_level = lagged[, "beta3"]
    )
    slope_matrix <- matrix(NA_real_, 3, 3,
                           dimnames = list(factor_names, factor_names))
    intercepts <- stats::setNames(rep(NA_real_, 3), factor_names)
  }

  for (i in 1:3) {
    factor_name <- factor_names[i]
    change_name <- change_names[i]

    # Prepare regression data. factor_level is the lagged level; the summary
    # statistics below describe the factor's own distribution and so use the
    # contemporaneous series.
    reg_data <- data.frame(
      factor_level = lagged[, factor_name],
      factor_change = factors_df[[change_name]]
    )
    observed_level <- factors_df[[factor_name]]

    # Calculate basic statistics
    factor_summary <- list(
      mean = mean(observed_level),
      sd = sd(observed_level),
      min = min(observed_level),
      max = max(observed_level),
      q25 = quantile(observed_level, 0.25),
      q50 = quantile(observed_level, 0.50),
      q75 = quantile(observed_level, 0.75)
    )

    # Linear regression: d_beta(t) ~ kappa*theta*dt - kappa*beta(t-1)*dt
    if (structure == "full") {
      model <- lm(factors_df[[change_name]] ~ beta1_level + beta2_level +
                    beta3_level, data = full_design)
      coefs <- coef(model)
      intercepts[i] <- coefs[1]
      slope_matrix[i, ] <- coefs[2:4]
    } else {
      model <- lm(factor_change ~ factor_level, data = reg_data)
    }

    # FIX: Extract parameters with more robust handling of NA/NULL values
    # Check if coefficients are NA or NULL
    if (is.null(coef(model)) || length(coef(model)) < 2 || any(is.na(coef(model)))) {
      warning(paste("Regression failed for factor", factor_name, "- using fallback values"))
      kappa <- 0.001  # Small positive value as fallback
      theta <- factor_summary$mean
      sigma <- factor_summary$sd * sqrt(dt)
    } else {
      # kappa from the own slope: b = -kappa * dt, so kappa = -b / dt.
      own_slope <- coef(model)[[if (structure == "full") 1L + i else 2L]]
      kappa <- -own_slope / dt

      if (is.na(kappa) || kappa <= 0) {
        # No mean reversion detected in this sample. Taking the magnitude keeps
        # the object usable by the simulator, but the sign is information and
        # is reported rather than swallowed: it means the factor drifted
        # without reverting over the estimation window, so theta is not
        # identified from the regression and the sample mean stands in.
        #
        # Only under the diagonal structure, though. There a non-positive
        # kappa really does mean the factor does not revert. In the full
        # system it means only that this factor's own coefficient is negative,
        # which a VAR is free to do while remaining stable overall: what
        # decides that is the eigenvalues of the whole matrix, checked once it
        # has been assembled. Both the kappa and the theta set here are
        # replaced by the matrix solution below.
        if (structure == "diagonal") {
          warning(sprintf(
            paste0("No mean reversion detected for %s (estimated kappa = %s). ",
                   "Using |kappa| and the sample mean for theta; treat the ",
                   "dynamics of this factor as unidentified over %s to %s."),
            factor_name, format(signif(kappa, 3)),
            format(date_range[1]), format(date_range[2])
          ))
        }
        kappa <- max(0.001, abs(kappa))
        theta <- factor_summary$mean
      } else {
        # theta from the intercept: a = kappa * theta * dt, so theta = -a / b.
        theta <- -coef(model)[[1]] / own_slope
      }

      sigma <- sd(residuals(model)) / sqrt(dt)

      # Handle missing or invalid values
      if (is.na(sigma) || sigma <= 0) {
        sigma <- 0.001  # Fallback to small positive value
      }
    }

    # Standard errors. kappa is a linear rescaling of the own slope, so its
    # standard error rescales with it; theta is a ratio of coefficients and
    # needs the delta method; sigma uses the asymptotic sd of a sd estimate.
    own_slope_idx <- if (structure == "full") 1L + i else 2L
    se_kappa <- NA_real_
    se_theta <- NA_real_
    se_sigma <- NA_real_
    vc <- tryCatch(stats::vcov(model), error = function(e) NULL)

    if (!is.null(vc) && all(dim(vc) >= own_slope_idx) &&
        !any(is.na(vc[c(1L, own_slope_idx), c(1L, own_slope_idx)]))) {
      coefs_all <- coef(model)
      a_hat <- coefs_all[[1]]
      b_hat <- coefs_all[[own_slope_idx]]
      var_a <- vc[1L, 1L]
      var_b <- vc[own_slope_idx, own_slope_idx]
      cov_ab <- vc[1L, own_slope_idx]

      se_kappa <- sqrt(var_b) / dt

      # theta = -a / b, so d(theta)/da = -1/b and d(theta)/db = a / b^2.
      if (is.finite(b_hat) && abs(b_hat) > 0) {
        var_theta <- (var_a + (a_hat / b_hat)^2 * var_b -
                        2 * (a_hat / b_hat) * cov_ab) / b_hat^2
        se_theta <- if (is.finite(var_theta) && var_theta >= 0) {
          sqrt(var_theta)
        } else {
          NA_real_
        }
      }

      df_resid <- stats::df.residual(model)
      if (!is.null(df_resid) && df_resid > 0) {
        se_sigma <- sigma / sqrt(2 * df_resid)
      }
    }

    # Calculate half-life of mean reversion
    half_life_years <- if (kappa > 0) log(2) / kappa else Inf
    half_life_days <- half_life_years * 252  # Convert to trading days

    # Half-life is log(2)/kappa, so its standard error follows from kappa's.
    half_life_ci <- c(lower = NA_real_, upper = NA_real_)
    if (is.finite(kappa) && kappa > 0 && is.finite(se_kappa)) {
      kappa_ci <- kappa + c(1, -1) * stats::qnorm(0.975) * se_kappa
      half_life_ci <- c(
        lower = if (kappa_ci[1] > 0) log(2) / kappa_ci[1] else 0,
        upper = if (kappa_ci[2] > 0) log(2) / kappa_ci[2] else Inf
      )
    }

    # Calculate pseudo-R^2 (model explanatory power)
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
      half_life_ci = half_life_ci,
      r_squared = r_squared,
      se = c(kappa = se_kappa, theta = se_theta, sigma = se_sigma),
      t_value = c(kappa = unname(kappa / se_kappa),
                  theta = unname(theta / se_theta)),
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
  dynamics$innovation_dates <- factors_df$Date
  dynamics$sampling <- sampling
  dynamics$dt <- dt
  dynamics$structure <- structure

  # Mean-reversion matrix. Delta_beta = K theta dt - K beta dt, so the fitted
  # slopes B give K = -B / dt and the intercepts a give theta = K^-1 a / dt.
  if (structure == "full") {
    kappa_matrix <- -slope_matrix / dt
    theta_vector <- tryCatch(
      stats::setNames(as.numeric(solve(kappa_matrix, intercepts / dt)),
                      factor_names),
      error = function(e) {
        warning("Mean-reversion matrix is singular; long-run means cannot be ",
                "solved for. Falling back to the sample means.")
        stats::setNames(
          vapply(factor_names, function(f) mean(factors_df[[f]]), numeric(1)),
          factor_names
        )
      }
    )

    # The per-factor entries carry the own-reversion speed and the full
    # system's long-run mean, so this object still works wherever a diagonal
    # one did. Cross terms live only in kappa_matrix.
    for (i in seq_along(factor_names)) {
      dynamics[[factor_names[i]]]$kappa <- kappa_matrix[i, i]
      dynamics[[factor_names[i]]]$theta <- theta_vector[[i]]
      dynamics[[factor_names[i]]]$half_life_years <-
        if (kappa_matrix[i, i] > 0) log(2) / kappa_matrix[i, i] else Inf
      dynamics[[factor_names[i]]]$half_life_days <-
        dynamics[[factor_names[i]]]$half_life_years * 252
    }
  } else {
    kappa_matrix <- diag(vapply(factor_names,
                                function(f) dynamics[[f]]$kappa, numeric(1)))
    dimnames(kappa_matrix) <- list(factor_names, factor_names)
    theta_vector <- stats::setNames(
      vapply(factor_names, function(f) dynamics[[f]]$theta, numeric(1)),
      factor_names
    )
  }

  dynamics$kappa_matrix <- kappa_matrix
  dynamics$theta_vector <- theta_vector

  # Analyze persistence of the model. The mean-reverting system is stable when
  # every eigenvalue of the mean-reversion matrix has a positive real part; for
  # a full matrix that is not the same question as whether each diagonal entry
  # is positive.
  eigenvalues <- kappa_matrix_eigenvalues(dynamics)
  dynamics$eigenvalues <- eigenvalues
  dynamics$system_stability <- all(Re(eigenvalues) > 0)

  if (!dynamics$system_stability) {
    warning(sprintf(
      paste0("The estimated factor system is not mean reverting: the ",
             "mean-reversion matrix has eigenvalues with non-positive real ",
             "parts (%s). Simulating from it will not be stationary; treat ",
             "the dynamics over %s to %s as unidentified."),
      paste(format(signif(Re(eigenvalues), 3)), collapse = ", "),
      format(date_range[1]), format(date_range[2])
    ))
  }

  class(dynamics) <- c("dns_dynamics", "list")

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

  # Use the estimated matrix when there is one. Under structure = "full" the
  # cross terms change the eigenvalues, and so the stability verdict, which a
  # diagonal reconstruction would miss.
  kappa_matrix <- dynamics$kappa_matrix
  if (is.null(kappa_matrix)) {
    kappa_matrix <- diag(c(
      dynamics$beta1$kappa,
      dynamics$beta2$kappa,
      dynamics$beta3$kappa
    ))
  }

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
#' @param lambda Numeric. The Nelson-Siegel decay parameter to carry into the
#'   simulation parameters. If NULL (the default) it is taken from
#'   \code{dynamics$lambda} when present, and otherwise falls back to 0.0609.
#'   \code{estimate_dns_dynamics()} does not estimate lambda, so pass the value
#'   fitted by \code{\link{fit_dns_model}} to avoid simulating off the default.
#'
#' @return A list of parameters for DNS simulation.
#'
#' @details
#' \code{kappa}, \code{theta} and \code{sigma} are the per-factor values, so a
#' dynamics object estimated with \code{structure = "full"} contributes only
#' the diagonal of its mean-reversion matrix here: \code{simulate_dns_factors()}
#' reverts each factor to its own long-run mean. The full matrix is carried
#' through as \code{kappa_matrix} for inspection.
#'
#' @examples
#' data(treasury_yields)
#'
#' # Month-ends only: the dynamics are estimated on them anyway
#' month_ends <- treasury_yields[
#'   !duplicated(format(treasury_yields$Date, "%Y-%m"), fromLast = TRUE),
#' ]
#' dns_result <- fit_dns_model(month_ends, n_days = nrow(month_ends),
#'                             show_progress = FALSE, warm_start = FALSE)
#' dynamics <- estimate_dns_dynamics(dns_result$parameters)
#'
#' # Pass the fitted lambda through, rather than leaving the 0.0609 default
#' latest <- dns_result$parameters[which.max(dns_result$parameters$Date), ]
#' sim_params <- dns_dynamics_to_sim_params(dynamics, lambda = latest$lambda)
#' str(sim_params)
#'
#' @export
dns_dynamics_to_sim_params <- function(dynamics, include_correlation = TRUE,
                                       lambda = NULL) {
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

  # Resolve lambda: an explicit argument wins, then a lambda carried on the
  # dynamics object, then the conventional Diebold-Li value.
  if (is.null(lambda)) {
    lambda <- if ("lambda" %in% names(dynamics)) dynamics$lambda else 0.0609
  }
  checkmate::assert_numeric(lambda, lower = 0, len = 1, any.missing = FALSE)

  # Package parameters for simulation
  params <- list(
    kappa = c(kappa1, kappa2, kappa3),
    theta = c(theta1, theta2, theta3),
    sigma = c(sigma1, sigma2, sigma3),
    names = c("Level", "Slope", "Curvature"),
    lambda = lambda
  )

  # Add correlation if requested
  if (include_correlation) {
    params$correlation <- dynamics$correlation
  }

  # Carry the mean-reversion matrix through for inspection. The simulator
  # reverts each factor to its own mean, so only the diagonal is used; with
  # structure = "full" the cross terms are dropped at this point.
  if (!is.null(dynamics$kappa_matrix)) {
    params$kappa_matrix <- dynamics$kappa_matrix
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
#' @param lambda Numeric. Fixed lambda value to use for calibration. See
#'   \code{\link{ns_lambda}} on its units; the 0.0609 fallback is a
#'   month-denominated constant and maturities here are in years.
#' @param start_date Date or character. Start date for the data window.
#' @param end_date Date or character. End date for the data window.
#' @param verbose Logical. If TRUE, print progress information.
#' @param api_key Character. FRED API key, passed through to the data fetch. If
#'   NULL the \code{fred_api_key} environment variable is used.
#'
#' @return A list containing:
#'   - params: List of DNS parameters for simulation, including the lambda
#'     actually fitted to the most recent curve
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
                                     verbose = FALSE, api_key = NULL) {
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
    api_key = api_key,
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

  # Estimate factor dynamics. estimate_dns_dynamics() owns the sampling
  # policy: its default monthly mode subsamples the factor series to
  # month-ends and estimates at dt = 1/12 (see its documentation for why
  # daily changes cannot be used directly). Windows with fewer than 20
  # month-ends fall back to the documented default dynamics, which short
  # samples cannot improve on.
  pars <- dns_result$parameters
  n_months <- length(unique(format(pars$Date, "%Y-%m")))

  if (n_months >= 20) {
    if (verbose) cat("Estimating factor dynamics (monthly sampling)...\n")
    dynamics <- tryCatch({
      estimate_dns_dynamics(pars, sampling = "monthly")
    }, error = function(e) {
      warning("Error estimating dynamics: ", e$message,
              ". Using default parameters.")
      NULL
    })
  } else {
    if (verbose) cat("Fewer than 20 month-end observations; dynamics cannot be",
                     "identified from this window. Using default values.\n")
    dynamics <- NULL
  }

  # The lambda fitted to the most recent curve, which the simulation
  # parameters below are built around.
  latest_row <- which.max(dns_result$parameters$Date)
  lambda_fitted <- dns_result$parameters$lambda[latest_row]

  # If dynamics estimation failed or we don't have enough data, use default values
  if (is.null(dynamics)) {
    # Use default parameters
    sim_params <- list(
      kappa = c(0.02, 0.3, 0.5),     # Speed of mean reversion
      theta = c(0.04, -0.02, 0.01),  # Long-run mean
      sigma = c(0.005, 0.01, 0.02),  # Volatility
      names = c("Level", "Slope", "Curvature"),
      lambda = lambda_fitted
    )
  } else {
    # Convert dynamics to simulation parameters, carrying the fitted lambda
    # rather than letting it fall back to the 0.0609 default.
    sim_params <- dns_dynamics_to_sim_params(dynamics, lambda = lambda_fitted)
  }

  # Get latest factors (for initial conditions). Selected by date rather than
  # by row position so the result does not depend on the frame's row order.
  latest_factors <- as.numeric(
    dns_result$parameters[latest_row, c("beta1", "beta2", "beta3")]
  )
  lambda_value <- lambda_fitted

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
