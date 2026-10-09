# One-step state-space estimation.
#
# The parameters divide into those a sample of realistic length determines
# well (lambda, the shock volatilities, the measurement error) and those it
# does not (the mean reversion speeds and long-run means). Assertions are
# tight on the first group and loose on the second; see
# ?estimate_dns_state_space for the measured errors.

sim_state_space <- function(n = 240, seed = 1,
                            kappa = c(0.20, 0.60, 1.20),
                            theta = c(0.040, -0.015, -0.010),
                            sigma = c(0.012, 0.018, 0.030),
                            rho12 = -0.35, rho13 = 0.15, rho23 = -0.25,
                            lambda = 0.55, meas_sd = 8e-4,
                            maturities = c(0.25, 1, 2, 5, 10, 30),
                            dt = 1 / 12) {
  set.seed(seed)
  corr <- matrix(c(1, rho12, rho13, rho12, 1, rho23, rho13, rho23, 1), 3, 3)
  A <- diag(exp(-kappa * dt), 3)
  ks <- outer(kappa, kappa, "+")
  Q <- corr * outer(sigma, sigma) * (1 - exp(-ks * dt)) / ks
  P0 <- corr * outer(sigma, sigma) / ks

  f <- matrix(NA_real_, n, 3)
  f[1, ] <- theta + as.numeric(t(chol(P0)) %*% stats::rnorm(3))
  cq <- t(chol(Q))
  for (t in 2:n) {
    f[t, ] <- theta + as.numeric(A %*% (f[t - 1, ] - theta)) +
      as.numeric(cq %*% stats::rnorm(3))
  }

  L <- cbind(1, (1 - exp(-lambda * maturities)) / (lambda * maturities),
             (1 - exp(-lambda * maturities)) / (lambda * maturities) -
               exp(-lambda * maturities))
  Y <- f %*% t(L) +
    matrix(stats::rnorm(n * length(maturities), sd = meas_sd),
           n, length(maturities))

  panel <- data.frame(Date = seq(as.Date("1990-01-31"), by = "month",
                                 length.out = n))
  for (k in seq_along(maturities)) {
    panel[[as.character(maturities[k])]] <- Y[, k]
  }
  list(panel = panel, factors = f, truth = list(
    kappa = kappa, theta = theta, sigma = sigma, corr = corr,
    lambda = lambda, meas_sd = meas_sd), dt = dt)
}


# One fit, computed on first use and reused. Each call to the estimator takes
# several seconds: BFGS with numerical gradients over a dozen parameters, with
# every trial point running the filter across the whole panel.
shared <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      sim <- sim_state_space(n = 84, seed = 13)
      cached <<- list(sim = sim, fit = estimate_dns_state_space(
        sim$panel, sampling = "all", dt = sim$dt, lambda = 0.55,
        measurement = "common"))
    }
    cached
  }
})


# ---- the pieces the filter is built on --------------------------------------

test_that("the correlation transform reaches only valid matrices, and round trips", {
  corr_from <- getFromNamespace("corr_from_unconstrained", "dnsr")
  to_unc <- getFromNamespace("unconstrained_from_corr", "dnsr")

  set.seed(4)
  for (i in 1:100) {
    R <- corr_from(stats::rnorm(3, sd = 1.2))
    expect_equal(diag(R), rep(1, 3))
    expect_equal(R, t(R))
    expect_gt(min(eigen(R, only.values = TRUE)$values), 0)
  }

  # Exact in the direction it is used: from a sample correlation matrix to
  # coordinates the optimiser starts from.
  for (i in 1:100) {
    R <- stats::cor(matrix(stats::rnorm(60), 20, 3))
    expect_equal(corr_from(to_unc(R)), R)
  }
})

test_that("the OU transition matches the simulator's own transition law", {
  ou <- getFromNamespace("ou_transition", "dnsr")
  kappa <- c(0.2, 0.6, 1.2); sigma <- c(0.012, 0.018, 0.03)
  corr <- matrix(c(1, -0.3, 0.1, -0.3, 1, -0.2, 0.1, -0.2, 1), 3, 3)
  dt <- 1 / 12
  tr <- ou(kappa, sigma, corr, dt)

  expect_equal(diag(tr$A), exp(-kappa * dt))
  ks <- outer(kappa, kappa, "+")
  expect_equal(tr$Q, corr * outer(sigma, sigma) * (1 - exp(-ks * dt)) / ks)

  # P0 is the stationary covariance, so A P0 A' + Q = P0. If that identity
  # fails the filter starts from the wrong prior.
  expect_equal(tr$A %*% tr$P0 %*% t(tr$A) + tr$Q, tr$P0)
  expect_gt(min(eigen(tr$Q, only.values = TRUE)$values), 0)
})

test_that("the filter's likelihood agrees with a direct Gaussian evaluation", {
  # Three dates and three maturities, small enough to write the joint
  # density out directly and compare.
  kf <- getFromNamespace("dns_kalman", "dnsr")
  set.seed(11)
  mats <- c(1, 5, 10)
  lam <- 0.5
  L <- getFromNamespace("ns_design", "dnsr")(mats, lam)
  kappa <- c(0.5, 0.8, 1.1); sigma <- c(0.01, 0.02, 0.03)
  theta <- c(0.03, -0.01, 0); corr <- diag(3); dt <- 1 / 12
  Y <- matrix(stats::rnorm(9, 0.03, 0.01), 3, 3)
  mv <- rep(1e-6, 3)

  got <- kf(Y, L, kappa, theta, sigma, corr, mv, dt)$logLik

  # Direct: the stacked observations are jointly Gaussian with mean L theta
  # and covariance built from the stationary and step covariances.
  ou <- getFromNamespace("ou_transition", "dnsr")
  tr <- ou(kappa, sigma, corr, dt)
  n <- nrow(Y)
  big_mean <- rep(as.numeric(L %*% theta), n)
  S <- matrix(0, 3 * n, 3 * n)
  for (i in 1:n) for (j in 1:n) {
    # For a stationary AR(1), Cov(f_i, f_j) = P0 (A^(j-i))' for j >= i.
    cov_f <- if (i <= j) {
      tr$P0 %*% diag(exp(-kappa * dt * (j - i)), 3)
    } else {
      diag(exp(-kappa * dt * (i - j)), 3) %*% tr$P0
    }
    block <- L %*% cov_f %*% t(L)
    if (i == j) block <- block + diag(mv, 3)
    S[((i - 1) * 3 + 1):(i * 3), ((j - 1) * 3 + 1):(j * 3)] <- block
  }
  v <- as.numeric(t(Y)) - big_mean
  ch <- chol(S)
  direct <- -0.5 * (length(v) * log(2 * pi) + 2 * sum(log(diag(ch))) +
                      sum(backsolve(ch, v, transpose = TRUE)^2))

  expect_equal(got, direct, tolerance = 1e-8)
})


# ---- the estimator ----------------------------------------------------------

test_that("the one-step fit recovers the parameters a sample can identify", {
  skip_on_cran()
  sim <- sim_state_space(n = 300, seed = 3)
  fit <- estimate_dns_state_space(sim$panel, sampling = "all", dt = sim$dt,
                                  lambda = NULL, measurement = "common")

  expect_s3_class(fit, "dns_state_space")
  expect_equal(fit$convergence, 0L)

  # Well determined: lambda, the shock volatilities, the measurement error.
  expect_equal(fit$lambda, sim$truth$lambda, tolerance = 0.05)
  got_sigma <- c(fit$beta1$sigma, fit$beta2$sigma, fit$beta3$sigma)
  expect_equal(got_sigma, sim$truth$sigma, tolerance = 0.15)
  expect_equal(stats::median(fit$measurement_sd), sim$truth$meas_sd,
               tolerance = 0.1)

  # The filtered states.
  for (k in 1:3) {
    expect_gt(stats::cor(fit$factors[[k + 1]], sim$factors[, k]), 0.97)
  }

  # Mean reversion speeds and long-run means are weakly determined at this
  # sample length: kappa is biased upward, and theta has a standard error of
  # the same order as its own value. These are bounds, not recovery claims.
  # See ?estimate_dns_state_space.
  got_kappa <- c(fit$beta1$kappa, fit$beta2$kappa, fit$beta3$kappa)
  expect_true(all(got_kappa > 0))
  expect_true(all(got_kappa < 5))
  expect_true(fit$system_stability)
})

test_that("a fixed lambda is honoured and estimation is reported", {
  # Estimating lambda is the most expensive path: it moves the loadings, so
  # every trial point rebuilds the design matrix and reruns the filter.
  skip_on_cran()
  sim <- sim_state_space(n = 72, seed = 5)
  fixed <- estimate_dns_state_space(sim$panel, sampling = "all", dt = sim$dt,
                                    lambda = 0.4, measurement = "common")
  expect_equal(fixed$lambda, 0.4)
  expect_false(fixed$lambda_estimated)

  free <- estimate_dns_state_space(sim$panel, sampling = "all", dt = sim$dt,
                                   lambda = NULL, measurement = "common")
  expect_true(free$lambda_estimated)
  # Estimating lambda cannot give a lower maximised likelihood than holding
  # it at a fixed value.
  expect_gte(free$logLik, fixed$logLik)
  expect_equal(free$n_params, fixed$n_params + 1L)
})

test_that("missing yields shrink the measurement equation, not the panel", {
  sim <- sim_state_space(n = 84, seed = 7)
  holed <- sim$panel
  # One maturity removed on a third of the dates, and one date removed
  # entirely.
  holed[seq(2, 60, by = 3), "30"] <- NA
  holed[seq(5, 60, by = 9), "0.25"] <- NA
  holed[60, setdiff(names(holed), "Date")] <- NA

  fit <- estimate_dns_state_space(holed, sampling = "all", dt = sim$dt,
                                  lambda = 0.55, measurement = "common")

  # Dates with at least three observations are retained; the blank one is
  # dropped.
  expect_equal(fit$n_dates, nrow(holed) - 1L)
  expect_false(as.Date(holed$Date[60]) %in% fit$factors$Date)
  expect_false(anyNA(fit$factors$beta1))
  expect_equal(fit$convergence, 0L)

  # Residuals are NA exactly where the yields were.
  expect_true(all(is.na(fit$residuals[is.na(as.matrix(
    holed[-60, setdiff(names(holed), "Date")]))])))
})

test_that("the fit is accepted by the simulation parameter converter", {
  fit <- shared()$fit

  # The continuous-time parameterisation needs no conversion step here.
  params <- dns_dynamics_to_sim_params(fit)
  expect_named(params, c("kappa", "theta", "sigma", "names", "lambda",
                         "correlation", "kappa_matrix"))
  expect_equal(params$lambda, 0.55)
  expect_equal(params$kappa,
               c(fit$beta1$kappa, fit$beta2$kappa, fit$beta3$kappa))

  # And it simulates forward.
  simd <- simulate_dns_factors(params, initial_factors = params$theta,
                               time_horizon = 1, n_steps = 12, n_paths = 5,
                               engine = "exact")
  expect_equal(dim(simd$factors), c(13L, 5L, 3L))
})

test_that("the smoothed factors use the whole sample and the filtered ones do not", {
  sh <- shared(); sim <- sh$sim; fit <- sh$fit

  expect_equal(nrow(fit$factors), nrow(fit$filtered))
  # They agree at the last date, where there is no future to smooth with,
  # and differ before it.
  n <- nrow(fit$factors)
  expect_equal(as.numeric(fit$factors[n, -1]),
               as.numeric(fit$filtered[n, -1]), tolerance = 1e-6)
  expect_false(isTRUE(all.equal(as.numeric(fit$factors[1, -1]),
                                as.numeric(fit$filtered[1, -1]))))

  # Smoothing uses more information, so it tracks the simulated factors
  # more closely.
  rmse <- function(a, b) sqrt(mean((a - b)^2))
  expect_lt(rmse(fit$factors$beta1, sim$factors[, 1]),
            rmse(fit$filtered$beta1, sim$factors[, 1]))
})

test_that("per-maturity measurement error is estimated separately", {
  skip_on_cran()
  sim <- sim_state_space(n = 72, seed = 17)
  per <- estimate_dns_state_space(sim$panel, sampling = "all", dt = sim$dt,
                                  lambda = 0.55, measurement = "per_maturity")
  common <- estimate_dns_state_space(sim$panel, sampling = "all", dt = sim$dt,
                                     lambda = 0.55, measurement = "common")

  expect_length(per$measurement_sd, ncol(sim$panel) - 1L)
  expect_gt(stats::sd(per$measurement_sd), 0)
  expect_equal(stats::sd(common$measurement_sd), 0)
  expect_gt(per$n_params, common$n_params)
  # More parameters on the same data cannot lower the maximised likelihood.
  expect_gte(per$logLik, common$logLik - 1e-6)
})

test_that("print and summary report without error", {
  fit <- shared()$fit

  out <- capture.output(print(fit))
  expect_true(any(grepl("one-step state-space", out)))
  expect_true(any(grepl("level \\(b1\\)", out)))
  expect_true(any(grepl("logLik", out)))

  s <- summary(fit)
  expect_s3_class(s, "summary.dns_state_space")
  sout <- capture.output(print(s))
  expect_true(any(grepl("Shock correlation", sout)))
  expect_true(any(grepl("Measurement error by maturity", sout)))
  expect_true(any(grepl("two-step estimates", sout)))
})

test_that("the two-step starting point is reported alongside the result", {
  fit <- shared()$fit
  expect_false(is.null(fit$two_step))
  expect_true(all(c("beta1", "beta2", "beta3", "correlation") %in%
                    names(fit$two_step)))
})

test_that("it rejects what it cannot fit", {
  expect_error(estimate_dns_state_space(data.frame(a = 1:5)), "Date")
  bad <- data.frame(Date = seq(as.Date("2020-01-31"), by = "month",
                               length.out = 10))
  bad$notamaturity <- 0.03
  bad$`5` <- 0.04
  expect_error(estimate_dns_state_space(bad), "maturity in years")

  short <- data.frame(Date = seq(as.Date("2020-01-31"), by = "month",
                                 length.out = 4))
  short$`1` <- 0.03; short$`5` <- 0.04; short$`10` <- 0.045
  expect_error(estimate_dns_state_space(short, sampling = "all"),
               "at least 5 are needed")

  # The monthly default can thin a daily panel below that minimum, so the
  # message says how many dates survived and how to change it.
  daily <- treasury_yields[1:60, , drop = FALSE]
  expect_error(estimate_dns_state_space(daily),
               "sampling = \"monthly\" kept 3 of the 60 dates")

  sim <- sim_state_space(n = 60, seed = 29)
  expect_error(estimate_dns_state_space(sim$panel, lambda = -1), "lambda")
  expect_error(estimate_dns_state_space(sim$panel, maxit = 0), "maxit")
})

test_that("monthly sampling thins a daily panel to month-ends", {
  fit <- estimate_dns_state_space(treasury_yields, n_days = 200,
                                  lambda = ns_lambda(3.5),
                                  measurement = "common")
  expect_equal(fit$sampling, "monthly")
  # Roughly one date per month across the 200 trading days used.
  expect_lt(fit$n_dates, 15L)
  expect_gt(fit$n_dates, 6L)
  expect_equal(format(fit$factors$Date, "%Y-%m"),
               unique(format(fit$factors$Date, "%Y-%m")))
})
