# The 0.3.0 performance work must not move any number. P-01 and P-03 are pure
# rewrites and are pinned against reference implementations here; P-02 adds an
# opt-in engine, so what is pinned there is the distribution it samples and the
# fact that the default is untouched.


# ---- P-01 -------------------------------------------------------------------

# The scalar triple loop that calculate_dns_yields() used before 0.3.0.
yields_scalar <- function(factors, maturities, lambda) {
  n_times <- dim(factors)[1]
  n_paths <- dim(factors)[2]
  n_mat <- length(maturities)
  out <- array(0, dim = c(n_times, n_mat, n_paths))

  loadings <- matrix(0, nrow = n_mat, ncol = 2)
  for (m_idx in seq_len(n_mat)) {
    m <- maturities[m_idx]
    if (m * lambda < 1e-10) {
      loadings[m_idx, ] <- c(1, 0)
    } else {
      l2 <- (1 - exp(-lambda * m)) / (lambda * m)
      loadings[m_idx, ] <- c(l2, l2 - exp(-lambda * m))
    }
  }

  for (t in seq_len(n_times)) {
    for (p in seq_len(n_paths)) {
      beta1 <- factors[t, p, 1]
      beta2 <- factors[t, p, 2]
      beta3 <- factors[t, p, 3]
      for (m_idx in seq_len(n_mat)) {
        out[t, m_idx, p] <- beta1 + beta2 * loadings[m_idx, 1] +
          beta3 * loadings[m_idx, 2]
      }
    }
  }
  out
}

# The scalar triple loop that calculate_dns_forwards() used before 0.3.0.
forwards_scalar <- function(factors, maturities, lambda) {
  n_times <- dim(factors)[1]
  n_paths <- dim(factors)[2]
  n_mat <- length(maturities)
  out <- array(0, dim = c(n_times, n_mat, n_paths))

  for (t in seq_len(n_times)) {
    for (p in seq_len(n_paths)) {
      for (m_idx in seq_len(n_mat)) {
        m <- maturities[m_idx]
        beta1 <- factors[t, p, 1]
        beta2 <- factors[t, p, 2]
        beta3 <- factors[t, p, 3]
        exp_term <- exp(-lambda * m)
        out[t, m_idx, p] <- beta1 + beta2 * exp_term +
          beta3 * lambda * m * exp_term
      }
    }
  }
  out
}

test_that("vectorised yields and forwards are identical to the scalar form", {
  set.seed(99)
  # Includes m = 0, which takes the limit branch of the loading formula.
  maturities <- c(0, 0.0833, 0.25, 1, 2, 7, 10, 30)
  factors <- array(stats::rnorm(13 * 5 * 3, 0.02, 0.01), dim = c(13, 5, 3))

  for (lambda in c(0, 0.0609, 0.3, 1.5)) {
    expect_identical(calculate_dns_yields(factors, maturities, lambda),
                     yields_scalar(factors, maturities, lambda))
    expect_identical(calculate_dns_forwards(factors, maturities, lambda),
                     forwards_scalar(factors, maturities, lambda))
  }

  # Single-path matrix input takes the other branch of both functions.
  single <- matrix(stats::rnorm(13 * 3, 0.02, 0.01), nrow = 13, ncol = 3)
  as_array <- array(single, dim = c(13, 1, 3))
  expect_equal(calculate_dns_yields(single, maturities, 0.0609),
               yields_scalar(as_array, maturities, 0.0609)[, , 1])
  expect_equal(calculate_dns_forwards(single, maturities, 0.0609),
               forwards_scalar(as_array, maturities, 0.0609)[, , 1])
})

test_that("degenerate shapes keep their dimensions", {
  maturities <- c(1, 5, 10)

  one_path <- array(stats::rnorm(7 * 1 * 3, 0.02, 0.01), dim = c(7, 1, 3))
  expect_equal(dim(calculate_dns_yields(one_path, maturities)), c(7L, 3L, 1L))

  one_time <- array(stats::rnorm(1 * 4 * 3, 0.02, 0.01), dim = c(1, 4, 3))
  expect_equal(dim(calculate_dns_yields(one_time, maturities)), c(1L, 3L, 4L))

  one_mat <- array(stats::rnorm(7 * 4 * 3, 0.02, 0.01), dim = c(7, 4, 3))
  expect_equal(dim(calculate_dns_yields(one_mat, 10)), c(7L, 1L, 4L))
})


# ---- P-02 -------------------------------------------------------------------

exact_params <- function() {
  list(
    kappa = c(0.3, 1.2, 2.5),
    theta = c(0.04, -0.02, 0.01),
    sigma = c(0.005, 0.01, 0.02),
    correlation = matrix(c(1, 0.4, 0.2, 0.4, 1, -0.3, 0.2, -0.3, 1), 3, 3)
  )
}

test_that("the exact engine samples the closed-form transition law", {
  params <- exact_params()
  x0 <- c(0.055, -0.040, 0.030)
  dt <- 0.25
  n_paths <- 2e5

  set.seed(11)
  sim <- simulate_dns_factors(params, x0, time_horizon = dt, n_steps = 1,
                              n_paths = n_paths, engine = "exact")
  drawn <- sim$factors[2, , ]

  mean_theory <- params$theta + (x0 - params$theta) * exp(-params$kappa * dt)

  # The step covariance couples the two reversion speeds. Each factor's own
  # variance combined with a Cholesky factor of the correlation matrix would
  # only be right if every kappa were equal.
  kappa_sum <- outer(params$kappa, params$kappa, "+")
  cov_theory <- params$correlation * outer(params$sigma, params$sigma) *
    (1 - exp(-kappa_sum * dt)) / kappa_sum

  # Compared as z-scores against the Monte Carlo standard error, which is what
  # the sampling noise here actually is. (expect_equal's tolerance is
  # relative, so it tightens on factors whose mean sits near zero.)
  std_error <- sqrt(diag(cov_theory) / n_paths)
  expect_lt(max(abs(colMeans(drawn) - mean_theory) / std_error), 5)

  expect_lt(max(abs(stats::cov(drawn) - cov_theory) / abs(cov_theory)), 0.05)
})

test_that("the exact engine reaches the stationary distribution", {
  params <- exact_params()

  set.seed(12)
  sim <- simulate_dns_factors(params, params$theta, time_horizon = 60,
                              n_steps = 720, n_paths = 3000, engine = "exact")
  terminal <- sim$factors[dim(sim$factors)[1], , ]

  # Stationary sd of an OU process is sigma / sqrt(2 kappa).
  sd_stationary <- params$sigma / sqrt(2 * params$kappa)
  expect_lt(max(abs(apply(terminal, 2, stats::sd) / sd_stationary - 1)), 0.08)

  # Paths are independent, so the mean's standard error is sd / sqrt(n_paths);
  # `terminal` is n_paths by 3.
  std_error <- sd_stationary / sqrt(nrow(terminal))
  expect_lt(max(abs(colMeans(terminal) - params$theta) / std_error), 5)
})

test_that("the exact engine agrees with yuima in distribution", {
  skip_if_not_installed("yuima")
  skip_on_cran()

  params <- exact_params()
  x0 <- c(0.055, -0.040, 0.030)

  set.seed(21)
  exact <- simulate_dns_factors(params, x0, time_horizon = 2, n_steps = 400,
                                n_paths = 300, engine = "exact")
  set.seed(21)
  viayuima <- suppressWarnings(
    simulate_dns_factors(params, x0, time_horizon = 2, n_steps = 400,
                         n_paths = 300, engine = "yuima")
  )

  last <- dim(exact$factors)[1]
  for (f in 1:3) {
    p_value <- suppressWarnings(stats::ks.test(
      exact$factors[last, , f], viayuima$factors[last, , f]
    )$p.value)
    expect_gt(p_value, 0.01)
  }

  # Same law, not the same paths: the engines consume randomness differently.
  expect_false(identical(exact$factors, viayuima$factors))
})

test_that("the exact engine needs no extra package and handles edge shapes", {
  params <- exact_params()
  x0 <- c(0.055, -0.040, 0.030)

  # Reaches no yuima code path at all.
  single <- simulate_dns_factors(params, x0, time_horizon = 1, n_steps = 5,
                                 n_paths = 1, engine = "exact")
  expect_equal(dim(single$factors), c(6L, 1L, 3L))
  expect_equal(as.numeric(single$factors[1, 1, ]), x0)
  expect_true(all(is.finite(single$factors)))
  expect_equal(single$engine, "exact")

  # kappa = 0 is the no-reversion limit: a driftless random walk whose
  # step variance is sigma^2 dt.
  flat <- params
  flat$kappa <- c(0, 0, 0)
  flat$correlation <- diag(3)
  set.seed(31)
  walk <- simulate_dns_factors(flat, x0, time_horizon = 1, n_steps = 10,
                               n_paths = 2000, engine = "exact")
  terminal <- walk$factors[11, , 1]
  expect_equal(stats::sd(terminal), flat$sigma[1], tolerance = 0.06)
  expect_equal(mean(terminal), x0[1], tolerance = 0.01)
})

test_that("the default engine is unchanged", {
  expect_equal(eval(formals(simulate_dns_factors)$engine)[1], "yuima")
  expect_equal(eval(formals(simulate_dns_yield_curves)$engine)[1], "yuima")
})

test_that("engine reaches simulate_dns_yield_curves", {
  params <- exact_params()
  sim <- simulate_dns_yield_curves(
    params, c(0.055, -0.040, 0.030), maturities = c(1, 5, 10),
    time_horizon = 2, n_steps = 10, n_paths = 5, engine = "exact"
  )
  expect_equal(sim$engine, "exact")
  expect_equal(dim(sim$yields), c(11L, 3L, 5L))

  # Yields are still the loading combination of the simulated factors.
  expect_equal(sim$yields,
               calculate_dns_yields(sim$factors, c(1, 5, 10), sim$lambda))
})


# ---- P-03 -------------------------------------------------------------------

test_that("plot frames keep their row order after the rbind rewrite", {
  set.seed(5)
  n_times <- 21
  maturities <- c(0.25, 1, 2, 5, 10, 30)
  factors <- array(stats::rnorm(n_times * 6 * 3, 0.02, 0.004),
                   dim = c(n_times, 6, 3))
  sim <- list(
    factors = factors,
    time = seq(0, 10, length.out = n_times),
    params = list(names = c("Level", "Slope", "Curvature")),
    maturities = maturities,
    lambda = 0.0609,
    yields = calculate_dns_yields(factors, maturities, 0.0609)
  )

  # Trajectories: time varies fastest within each path.
  traj <- plot_simulated_factors(sim, path_subset = 1:4,
                                 plot_trajectories = TRUE)$level$data
  expect_equal(nrow(traj), n_times * 4)
  expect_equal(traj$Time[seq_len(n_times)], sim$time)
  expect_equal(unique(traj$Path), paste("Path", 1:4))
  expect_equal(traj$Value, as.vector(factors[, 1:4, 1]))

  # Yield curves: maturity inner, path middle, time outer.
  curves <- plot_simulated_yields(sim, time_points = c(0, 1, 5, 10),
                                  path_subset = 1:3)$yield_curves$data
  expect_equal(nrow(curves), 4 * 3 * length(maturities))
  expect_equal(curves$Maturity[seq_along(maturities)], maturities)
  expect_equal(as.character(curves$Time[1]), "Year 0")
  expect_equal(curves$Path[1], "Path 1")

  # Spot-check a value against the source array.
  t_idx <- which.min(abs(sim$time - 5))
  row <- which(as.character(curves$Time) == "Year 5" &
                 curves$Path == "Path 2" & curves$Maturity == 10)
  expect_equal(curves$Yield[row], sim$yields[t_idx, 5, 2] * 100)

  # Term structure: time inner, maturity outer.
  term <- plot_simulated_yields(sim, time_points = 5)$term_structure$data
  expect_equal(nrow(term), n_times * length(maturities))
  expect_equal(term$Time[seq_len(n_times)], sim$time)
  expect_equal(term$Mean[1], mean(sim$yields[1, 1, ]) * 100)
})
