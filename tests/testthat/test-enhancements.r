# Tests for the 0.4.0 enhancements, and for the lagged-level regression fix in
# estimate_dns_dynamics().

ns_curve <- function(maturities, betas, lambda) {
  l2 <- (1 - exp(-lambda * maturities)) / (lambda * maturities)
  l3 <- l2 - exp(-lambda * maturities)
  betas[1] + betas[2] * l2 + betas[3] * l3
}

# An OU system simulated at a stated step, returned in the shape
# estimate_dns_dynamics() expects.
ou_panel <- function(n, kappa, theta, sigma, dt, seed = 1, by = "month") {
  set.seed(seed)
  k_is_matrix <- is.matrix(kappa)
  b <- matrix(NA_real_, n, 3)
  b[1, ] <- theta
  for (i in 2:n) {
    pull <- if (k_is_matrix) {
      as.numeric(kappa %*% (theta - b[i - 1, ]))
    } else {
      kappa * (theta - b[i - 1, ])
    }
    b[i, ] <- b[i - 1, ] + pull * dt + sigma * sqrt(dt) * stats::rnorm(3)
  }
  out <- data.frame(
    Date = seq(as.Date("1960-01-28"), by = by, length.out = n),
    beta1 = b[, 1], beta2 = b[, 2], beta3 = b[, 3]
  )
  for (f in c("beta1", "beta2", "beta3")) {
    out[[paste0(f, "_change")]] <- c(NA, diff(out[[f]]))
  }
  out
}


# ---- the OU regression uses the lagged level ---------------------------------

test_that("the OU regression is on the lagged level, giving a negative slope", {
  panel <- ou_panel(600, kappa = c(0.8, 1.5, 2.5), theta = c(0.04, -0.02, 0.01),
                    sigma = c(0.01, 0.014, 0.03), dt = 1 / 12)

  dynamics <- estimate_dns_dynamics(panel, sampling = "monthly")

  for (f in c("beta1", "beta2", "beta3")) {
    # b = -kappa * dt, so a mean-reverting series must give a negative slope.
    # Regressing on the contemporaneous level puts beta(t) on both sides and
    # returns +(1 - phi) instead, which made kappa come out negative and get
    # rescued by abs().
    expect_lt(stats::coef(dynamics[[f]]$model)[[2]], 0)
    expect_gt(dynamics[[f]]$kappa, 0)
  }

  # theta now comes from -a/b rather than always falling back to the sample
  # mean, which is what happened while the slope's sign was inverted.
  for (f in c("beta1", "beta2", "beta3")) {
    model <- dynamics[[f]]$model
    expect_equal(
      dynamics[[f]]$theta,
      -stats::coef(model)[[1]] / stats::coef(model)[[2]]
    )
  }

  # The lagged level is the contemporaneous level minus the change.
  reg_levels <- stats::model.frame(dynamics$beta1$model)$factor_level
  expect_false(anyNA(reg_levels))
})

test_that("a series with no mean reversion is reported, not silently flipped", {
  # A pure random walk: no reversion to find.
  set.seed(7)
  n <- 400
  walk <- cumsum(stats::rnorm(n, 0, 0.01)) + 0.04
  panel <- data.frame(Date = seq(as.Date("1960-01-28"), by = "month",
                                 length.out = n),
                      beta1 = walk, beta2 = walk, beta3 = walk)
  for (f in c("beta1", "beta2", "beta3")) {
    panel[[paste0(f, "_change")]] <- c(NA, diff(panel[[f]]))
  }

  # Whether this particular sample shows a positive or negative slope is luck;
  # what matters is that a non-reverting estimate warns rather than being
  # turned positive in silence.
  warnings_seen <- character(0)
  dynamics <- withCallingHandlers(
    estimate_dns_dynamics(panel, sampling = "monthly"),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  if (dynamics$beta1$kappa <= 0.002) {
    expect_true(any(grepl("No mean reversion detected", warnings_seen)))
  }
  expect_true(all(vapply(c("beta1", "beta2", "beta3"),
                         function(f) dynamics[[f]]$kappa > 0, logical(1))))
})


# ---- E-02 --------------------------------------------------------------------

test_that("the profiled fit is the exact least-squares solution at fixed lambda", {
  maturities <- c(0.25, 1, 2, 5, 10, 30)
  yields <- ns_curve(maturities, c(0.045, -0.02, 0.01), 0.3) +
    c(1, -1, 1, -1, 1, -1) * 2e-4

  fit <- fit_ns_model(maturities, yields, lambda_constraint = "fixed",
                      initial_lambda = 0.0609, method = "profile")

  l2 <- (1 - exp(-0.0609 * maturities)) / (0.0609 * maturities)
  l3 <- l2 - exp(-0.0609 * maturities)
  reference <- stats::coef(stats::lm(yields ~ l2 + l3))

  expect_equal(unname(fit$parameters[1:3]), unname(reference))
  expect_equal(fit$parameters[["lambda"]], 0.0609)
})

test_that("the profiled fit recovers a synthetic curve far more closely", {
  maturities <- c(0.25, 1, 2, 5, 10, 30)
  truth <- c(0.045, -0.02, 0.01)
  yields <- ns_curve(maturities, truth, 0.3)

  profiled <- fit_ns_model(maturities, yields, method = "profile")
  numeric_fit <- fit_ns_model(maturities, yields, method = "optim")

  # The curve is exactly Nelson-Siegel, so an exact solver should essentially
  # reproduce it.
  expect_equal(profiled$parameters[["lambda"]], 0.3, tolerance = 1e-3)
  expect_equal(unname(profiled$parameters[1:3]), truth, tolerance = 1e-4)
  expect_lt(profiled$rmse, 1e-8)
  expect_lt(profiled$rmse, numeric_fit$rmse)
})

test_that("the profiled fit reports a lambda resting on its bound", {
  maturities <- c(0.25, 1, 2, 5, 10, 30)
  # Generated with a decay well above the default ceiling of 0.5.
  yields <- ns_curve(maturities, c(0.045, -0.02, 0.03), 1.2)

  expect_warning(
    clamped <- fit_ns_model(maturities, yields, method = "profile"),
    "rests on a parameter bound"
  )
  expect_true(clamped$at_bound[["lambda"]])
  expect_equal(clamped$parameters[["lambda"]], 0.5, tolerance = 1e-5)

  # Given room, it finds the interior optimum and reports no bound.
  free <- fit_ns_model(maturities, yields, method = "profile",
                       lower_bounds = c(-1, -1, -1, 0.01),
                       upper_bounds = c(1, 1, 1, 5))
  expect_false(any(free$at_bound))
  expect_equal(free$parameters[["lambda"]], 1.2, tolerance = 1e-3)
})

test_that("fit_dns_model reports bound hits once for the whole panel", {
  maturities <- c(0.25, 1, 2, 5, 10, 30)
  dates <- seq(as.Date("2023-01-02"), by = "day", length.out = 12)
  panel <- data.frame(Date = dates)
  for (k in seq_along(maturities)) {
    panel[[as.character(maturities[k])]] <-
      ns_curve(maturities[k], c(0.045, -0.02, 0.03), 1.2) +
      1e-4 * seq_along(dates)
  }

  # One warning naming the count, not one per date.
  expect_warning(
    fit <- fit_dns_model(panel, n_days = 12, show_progress = FALSE,
                         method = "profile"),
    "rest on a bound in some fits"
  )
  expect_true(is.matrix(fit$at_bound))
  expect_equal(nrow(fit$at_bound), nrow(fit$parameters))
  expect_true(all(fit$at_bound[, "lambda"]))

  # parameters stays numeric apart from Date, so coercion still works.
  numeric_cols <- setdiff(names(fit$parameters), "Date")
  expect_true(all(vapply(fit$parameters[numeric_cols], is.numeric, logical(1))))
})


# ---- E-03 --------------------------------------------------------------------

test_that("the full structure recovers a known mean-reversion matrix", {
  kappa <- matrix(c(0.80, -0.30, 0.10,
                    -0.20, 1.20, 0.25,
                    0.05, -0.15, 2.00), 3, 3, byrow = TRUE)
  theta <- c(0.040, -0.015, 0.008)
  panel <- ou_panel(1200, kappa = kappa, theta = theta,
                    sigma = c(0.010, 0.014, 0.030), dt = 1 / 12, seed = 4)

  dynamics <- estimate_dns_dynamics(panel, structure = "full",
                                    sampling = "monthly")

  expect_equal(dynamics$structure, "full")
  expect_equal(dim(dynamics$kappa_matrix), c(3L, 3L))
  expect_lt(max(abs(dynamics$kappa_matrix - kappa)), 0.6)
  # Compared on an absolute scale: theta3 is 0.008, so a relative tolerance
  # would demand far more precision of it than of theta1.
  expect_lt(max(abs(unname(dynamics$theta_vector) - theta)), 0.004)

  # The true system is stable, so every eigenvalue has a positive real part.
  expect_true(all(Re(dynamics$eigenvalues) > 0))
  expect_true(dynamics$system_stability)

  # Per-factor entries carry the diagonal, so the object still feeds the
  # simulator.
  expect_equal(
    c(dynamics$beta1$kappa, dynamics$beta2$kappa, dynamics$beta3$kappa),
    diag(dynamics$kappa_matrix), ignore_attr = TRUE
  )
  sim_params <- dns_dynamics_to_sim_params(dynamics)
  expect_equal(sim_params$kappa, diag(dynamics$kappa_matrix),
               ignore_attr = TRUE)
  expect_equal(sim_params$kappa_matrix, dynamics$kappa_matrix)
})

test_that("the diagonal structure stays diagonal and is the default", {
  panel <- ou_panel(400, kappa = c(0.8, 1.5, 2.5), theta = c(0.04, -0.02, 0.01),
                    sigma = c(0.01, 0.014, 0.03), dt = 1 / 12)

  dynamics <- estimate_dns_dynamics(panel)
  expect_equal(dynamics$structure, "diagonal")
  expect_true(all(dynamics$kappa_matrix[upper.tri(dynamics$kappa_matrix)] == 0))
  expect_true(all(dynamics$kappa_matrix[lower.tri(dynamics$kappa_matrix)] == 0))

  # Eigenvalues of a diagonal matrix are its diagonal.
  expect_setequal(Re(dynamics$eigenvalues), diag(dynamics$kappa_matrix))

  # Each factor's own regression has one regressor.
  expect_length(stats::coef(dynamics$beta1$model), 2)
})


# ---- E-04 --------------------------------------------------------------------

test_that("standard errors follow from the fitted regression", {
  panel <- ou_panel(500, kappa = c(0.8, 1.5, 2.5), theta = c(0.04, -0.02, 0.01),
                    sigma = c(0.01, 0.014, 0.03), dt = 1 / 12)
  dynamics <- estimate_dns_dynamics(panel, sampling = "monthly")
  dt <- dynamics$dt

  for (f in c("beta1", "beta2", "beta3")) {
    fac <- dynamics[[f]]
    vc <- stats::vcov(fac$model)

    # kappa = -b/dt, a linear rescaling, so its standard error rescales too.
    expect_equal(fac$se[["kappa"]], sqrt(vc[2, 2]) / dt)
    expect_true(is.finite(fac$se[["theta"]]))
    expect_equal(fac$se[["sigma"]],
                 fac$sigma / sqrt(stats::df.residual(fac$model)* 2))
    expect_equal(fac$t_value[["kappa"]], fac$kappa / fac$se[["kappa"]])

    # The half-life interval brackets the point estimate.
    expect_lte(fac$half_life_ci[["lower"]], fac$half_life_years)
    expect_gte(fac$half_life_ci[["upper"]], fac$half_life_years)
  }
})


# ---- E-05 --------------------------------------------------------------------

test_that("classes are additive and the lists still behave as lists", {
  data(treasury_yields)
  # The recent curve wants a faster decay than the default lambda ceiling, so
  # some of these fits sit on it; that is reported and is not what this test is
  # about.
  fit <- suppressWarnings(
    fit_dns_model(treasury_yields, n_days = 40, show_progress = FALSE)
  )

  expect_s3_class(fit, "dns_fit")
  expect_true(is.list(fit))
  expect_true(inherits(fit, "list"))

  # Every element is reachable under the name it always had.
  for (nm in c("parameters", "fitted_yields", "maturities", "success_rate")) {
    expect_false(is.null(fit[[nm]]))
  }
  expect_silent(checkmate::assert_list(fit))

  # Functions that validate on names still accept it.
  expect_s3_class(plot_factor_time_series(fit)$factors, "ggplot")

  # print and summary return their input invisibly.
  expect_output(print(fit), "Dynamic Nelson-Siegel fit")
  expect_s3_class(summary(fit), "summary.dns_fit")
  expect_output(print(summary(fit)), "Fit error")

  curve <- fit_ns_model(c(0.25, 1, 2, 5, 10, 30),
                        ns_curve(c(0.25, 1, 2, 5, 10, 30),
                                 c(0.045, -0.02, 0.01), 0.3))
  expect_s3_class(curve, "ns_fit")
  expect_output(print(curve), "Nelson-Siegel curve fit")

  params <- list(kappa = c(0.5, 1, 2), theta = c(0.04, -0.02, 0.01),
                 sigma = c(0.005, 0.01, 0.02))
  sim <- simulate_dns_factors(params, c(0.05, -0.01, 0.01), n_steps = 10,
                              n_paths = 5, engine = "exact", seed = 1)
  expect_s3_class(sim, "dns_simulation")
  expect_output(print(sim), "DNS simulation")
})


# ---- E-06 --------------------------------------------------------------------

test_that("discount factors are the discounted DNS curve", {
  maturities <- c(0.25, 1, 2, 5, 10, 30)
  factors <- matrix(c(0.045, -0.02, 0.01), nrow = 1)

  discount <- dns_discount_factors(factors, maturities, 0.3)
  yields <- calculate_dns_yields(factors, maturities, 0.3)

  expect_equal(as.numeric(discount), as.numeric(exp(-yields * maturities)))
  expect_equal(dns_discount_factors(factors, 0, 0.3)[1, 1], 1)
  expect_true(all(diff(as.numeric(discount)) < 0))

  # Zero prices are the discount factors times face.
  expect_equal(dns_zero_prices(factors, maturities, 0.3),
               100 * discount)
  expect_equal(dns_zero_prices(factors, maturities, 0.3, face = 1), discount)

  # Shapes follow calculate_dns_yields.
  paths <- array(rep(c(0.045, -0.02, 0.01), each = 12), dim = c(4, 3, 3))
  expect_equal(dim(dns_discount_factors(paths, maturities, 0.3)),
               c(4L, length(maturities), 3L))
})

test_that("par yields price their bond at par", {
  maturities <- c(0.5, 1, 2, 5, 10, 30)
  factors <- matrix(c(0.045, -0.02, 0.01), nrow = 1)
  lambda <- 0.3

  par <- dns_par_yields(factors, maturities, lambda)

  for (k in seq_along(maturities)) {
    horizon <- maturities[k]
    coupon_dates <- seq_len(horizon * 2) / 2
    discount <- as.numeric(dns_discount_factors(factors, coupon_dates, lambda))

    # Coupons plus redemption, per unit of face.
    price <- sum(par[1, k] / 2 * discount) + utils::tail(discount, 1)
    expect_equal(price, 1, tolerance = 1e-12)
  }
})

test_that("par yields are undefined below one coupon period", {
  factors <- matrix(c(0.045, -0.02, 0.01), nrow = 1)

  expect_warning(
    par <- dns_par_yields(factors, c(0.0833, 0.25, 1), 0.3),
    "not a whole number"
  )
  expect_true(all(is.na(par[1, 1:2])))
  expect_false(is.na(par[1, 3]))

  # A quarterly bond makes a 3-month maturity well defined.
  expect_silent(dns_par_yields(factors, c(0.25, 1), 0.3, frequency = 4))

  # Every maturity undefined: all NA, still the right shape.
  expect_warning(none <- dns_par_yields(factors, c(0.0833, 0.25), 0.3))
  expect_true(all(is.na(none)))
  expect_equal(dim(none), c(1L, 2L))
})

test_that("par and zero yields differ by a convention term and a shape term", {
  maturities <- c(1, 2, 5, 10, 30)

  # Upward sloping zero curve: beta2 negative lifts the long end.
  rising <- matrix(c(0.05, -0.03, 0.00), nrow = 1)
  zero_r <- calculate_dns_yields(rising, maturities, 0.3)[1, ]
  par_r <- dns_par_yields(rising, maturities, 0.3)[1, ]

  # Far out, the shape term dominates: coupons are discounted at the lower
  # rates earlier in the curve, so the par yield falls below the zero yield.
  expect_lt(par_r[5], zero_r[5])
  expect_lt(par_r[4], zero_r[4])

  # At one year the shape term is tiny and the compounding convention
  # dominates: a continuously compounded zero against a semiannual par rate
  # differs by about y^2 / (2f).
  expect_gt(par_r[1], zero_r[1])
  expect_equal(par_r[1] - zero_r[1], zero_r[1]^2 / 4, tolerance = 0.1)

  # Inverted zero curve: both terms push the par yield up, at every maturity.
  inverted <- matrix(c(0.03, 0.03, 0.00), nrow = 1)
  zero_i <- calculate_dns_yields(inverted, maturities, 0.3)[1, ]
  par_i <- dns_par_yields(inverted, maturities, 0.3)[1, ]
  expect_true(all(par_i > zero_i))
})


# ---- E-07 --------------------------------------------------------------------

test_that("the bundled Treasury panel is well formed", {
  data(treasury_yields)

  expect_s3_class(treasury_yields, "data.frame")
  expect_equal(names(treasury_yields)[1], "Date")
  expect_s3_class(treasury_yields$Date, "Date")
  expect_false(anyNA(treasury_yields))
  expect_false(is.unsorted(treasury_yields$Date))
  expect_gt(nrow(treasury_yields), 1400)

  # Maturity columns are named by maturity in years, which is what the fitting
  # functions read.
  maturities <- suppressWarnings(
    as.numeric(setdiff(names(treasury_yields), "Date"))
  )
  expect_false(anyNA(maturities))
  expect_false(is.unsorted(maturities))

  # Decimal form, not percentage points.
  yields <- unlist(treasury_yields[, -1])
  expect_lt(max(yields), 0.12)
  expect_gte(min(yields), 0)

  # Enough month-ends for the default monthly dynamics estimation.
  expect_gte(length(unique(format(treasury_yields$Date, "%Y-%m"))), 20)
})


# ---- E-08 --------------------------------------------------------------------

test_that("a cached response is served without a network call or an API key", {
  cache_dir <- file.path(tempdir(), "dnsr-cache-test")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  reference <- data.frame(date = as.Date("2024-01-01") + 0:2,
                          value = c(3.9, 4.0, 4.1))
  saveRDS(reference,
          file.path(cache_dir, "fred_DGS10_2024-01-01_2024-01-03.rds"))

  old_key <- Sys.getenv("fred_api_key")
  Sys.unsetenv("fred_api_key")
  on.exit(if (nzchar(old_key)) Sys.setenv(fred_api_key = old_key), add = TRUE)

  # A cache hit sends nothing, so it must not demand a key.
  got <- fetch_fred_series("DGS10", start_date = "2024-01-01",
                           end_date = "2024-01-03", cache_dir = cache_dir)
  expect_equal(got, reference)

  # A request with nothing cached still needs one.
  expect_error(
    fetch_fred_series("DGS2", start_date = "2024-01-01",
                      end_date = "2024-01-03", cache_dir = cache_dir),
    "API key not found"
  )
})

test_that("retry and cache arguments are exposed end to end", {
  expect_true(all(c("retries", "cache_dir") %in%
                    names(formals(fetch_fred_series))))
  expect_true(all(c("retries", "cache_dir") %in%
                    names(formals(get_treasury_yield_data))))

  # Appended, so existing positional calls still land where they did.
  expect_equal(names(formals(fetch_fred_series))[1:5],
               c("series_id", "api_key", "start_date", "end_date", "verbose"))
})


# ---- Nelson-Siegel-Svensson -------------------------------------------------

test_that("the Svensson fitter recovers curves generated by its own model", {
  fit_svensson <- getFromNamespace("fit_svensson", "dnsr")
  nss_design <- getFromNamespace("nss_design", "dnsr")
  maturities <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)

  cases <- list(
    list(betas = c(0.045, -0.020, 0.010, 0.030), lambdas = c(0.60, 4.0)),
    list(betas = c(0.030, 0.020, -0.040, 0.050), lambdas = c(0.25, 1.5)),
    list(betas = c(0.050, -0.030, 0.020, -0.020), lambdas = c(1.00, 6.0)),
    list(betas = c(0.035, 0.010, -0.020, 0.030), lambdas = c(0.80, 10.0))
  )

  for (case in cases) {
    yields <- as.numeric(
      nss_design(maturities, case$lambdas[1], case$lambdas[2]) %*% case$betas
    )
    fit <- fit_svensson(maturities, yields)

    expect_equal(fit$parameters[["lambda1"]], case$lambdas[1],
                 tolerance = 1e-4)
    expect_equal(fit$parameters[["lambda2"]], case$lambdas[2],
                 tolerance = 1e-4)
    expect_lt(max(abs(fit$parameters[1:4] - case$betas)), 1e-6)
    expect_lt(fit$rmse, 1e-9)
  }
})

test_that("the Svensson fit stays out of its degenerate corners", {
  fit_svensson <- getFromNamespace("fit_svensson", "dnsr")
  data(treasury_yields)
  maturities <- as.numeric(setdiff(names(treasury_yields), "Date"))

  # A spread of dates across the panel's regimes. Unbounded, the search found
  # lambda1 = 1e-4 with betas of +/-16,800 on one of these and lambda1 = 40 on
  # others -- degenerate parameterisations with flattering residuals. Every
  # tenth date was swept while the bounds were being worked out; the stride is
  # wider here so the suite stays quick under coverage instrumentation.
  for (i in seq(1, nrow(treasury_yields), by = 75)) {
    fit <- fit_svensson(maturities, as.numeric(treasury_yields[i, -1]))

    expect_gte(fit$parameters[["lambda1"]], 0.02)
    expect_lte(fit$parameters[["lambda1"]], 2)
    expect_gt(fit$parameters[["lambda2"]], fit$parameters[["lambda1"]])

    # Yields are decimals, so a factor above 1 in magnitude means the betas are
    # cancelling rather than describing the curve.
    expect_lt(max(abs(fit$parameters[1:4])), 1)
    expect_lt(fit$rmse, 0.002)
  }
})

test_that("Svensson nests Nelson-Siegel, so it fits at least as well", {
  data(treasury_yields)

  comparison <- suppressWarnings(compare_term_structure_models(
    treasury_yields, models = c("NS", "NSS")
  ))
  expect_equal(nrow(comparison$metrics), 2)
  expect_lte(comparison$metrics$RMSE[comparison$metrics$Model ==
                                       "Nelson-Siegel-Svensson"],
             comparison$metrics$RMSE[comparison$metrics$Model ==
                                       "Nelson-Siegel"])
})

test_that("Svensson needs enough maturities to be identified", {
  data(treasury_yields)
  narrow <- treasury_yields[, c("Date", "1", "10", "30")]

  expect_warning(
    comparison <- compare_term_structure_models(narrow, models = c("NS", "NSS")),
    "at least 5 observed maturities"
  )
  expect_equal(nrow(comparison$metrics), 1)
})


# ---- Smith-Wilson -----------------------------------------------------------

test_that("Smith-Wilson interpolates the observed curve exactly", {
  skip_if_not_installed("SmithWilsonYieldCurve")
  data(treasury_yields)

  comparison <- suppressWarnings(compare_term_structure_models(
    treasury_yields, models = c("NS", "SmithWilson")
  ))
  smith_wilson <- comparison$metrics[comparison$metrics$Model == "Smith-Wilson", ]

  expect_equal(nrow(smith_wilson), 1)
  # One kernel term per instrument, so it reproduces every input. The error is
  # zero by construction and is not comparable with a parametric fit's.
  expect_lt(smith_wilson$RMSE, 1e-6)
  expect_lt(smith_wilson$MaxError, 1e-6)
})

test_that("Smith-Wilson is skipped cleanly when its package is absent", {
  data(treasury_yields)

  if (requireNamespace("SmithWilsonYieldCurve", quietly = TRUE)) {
    skip("SmithWilsonYieldCurve is installed, so the guard cannot be exercised")
  }
  expect_warning(
    comparison <- compare_term_structure_models(treasury_yields,
                                                models = c("NS", "SmithWilson")),
    "SmithWilsonYieldCurve"
  )
  expect_equal(nrow(comparison$metrics), 1)
})

test_that("SmithWilson is opt-in and its arguments are appended", {
  args <- names(formals(compare_term_structure_models))
  expect_equal(args[1:4],
               c("yield_data", "date_to_compare", "models", "spline_df"))
  expect_true(all(c("ufr", "alpha") %in% args))

  # Not in the default set, so no optional dependency is pulled in by default.
  expect_false("SmithWilson" %in% eval(formals(
    compare_term_structure_models)$models))
})


# ---- which diagnostic belongs to which structure ---------------------------

test_that("a negative own-diagonal in a stable full system is not a warning", {
  # A full mean-reversion matrix is free to have a negative entry on its own
  # diagonal while the system reverts overall: what decides that is the
  # eigenvalues of the whole matrix. Reporting "no mean reversion" per factor
  # here was wrong, and fired on the bundled panel.
  data(treasury_yields)
  month_ends <- treasury_yields[
    !duplicated(format(treasury_yields$Date, "%Y-%m"), fromLast = TRUE),
  ]
  fit <- fit_dns_model(month_ends, n_days = nrow(month_ends),
                       show_progress = FALSE, warm_start = FALSE,
                       lambda_constraint = "fixed", initial_lambda = 0.0609)

  seen <- character(0)
  dynamics <- withCallingHandlers(
    estimate_dns_dynamics(fit$parameters, structure = "full"),
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_true(dynamics$system_stability)
  expect_true(all(Re(dynamics$eigenvalues) > 0))
  expect_false(any(grepl("No mean reversion detected", seen)))
  expect_false(any(grepl("not mean reverting", seen)))

  # The point of the test: at least one own-diagonal entry is negative, so the
  # old per-factor check would have fired.
  expect_true(any(diag(dynamics$kappa_matrix) < 0))
})

test_that("a genuinely explosive system is reported by each structure", {
  # phi > 1, so the change loads positively on the lagged level and kappa
  # really is negative.
  set.seed(3)
  n <- 240
  phi <- 1.04
  betas <- matrix(NA_real_, n, 3)
  betas[1, ] <- c(0.04, -0.01, 0.01)
  for (i in 2:n) {
    betas[i, ] <- phi * betas[i - 1, ] + stats::rnorm(3, 0, 0.002)
  }
  panel <- data.frame(
    Date = seq(as.Date("1990-01-31"), by = "month", length.out = n),
    beta1 = betas[, 1], beta2 = betas[, 2], beta3 = betas[, 3]
  )

  collect <- function(structure) {
    seen <- character(0)
    dynamics <- withCallingHandlers(
      estimate_dns_dynamics(panel, structure = structure),
      warning = function(w) {
        seen <<- c(seen, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    list(dynamics = dynamics, warnings = seen)
  }

  # Diagonal: each factor fails to revert on its own terms, and says so.
  diagonal <- collect("diagonal")
  expect_true(all(grepl("No mean reversion detected", diagonal$warnings)))
  expect_length(diagonal$warnings, 3)

  # Full: the system-level verdict, from the eigenvalues.
  full <- collect("full")
  expect_false(full$dynamics$system_stability)
  expect_true(any(Re(full$dynamics$eigenvalues) <= 0))
  expect_true(any(grepl("not mean reverting", full$warnings)))
})


# ---- ns_lambda ---------------------------------------------------------------

test_that("ns_lambda puts the curvature peak where it says", {
  for (peak in c(0.5, 1, 2.5, 10, 30)) {
    lambda <- ns_lambda(peak)
    grid <- seq(peak / 20, peak * 20, length.out = 20001)
    loading <- (1 - exp(-lambda * grid)) / (lambda * grid) - exp(-lambda * grid)
    expect_equal(grid[which.max(loading)], peak, tolerance = 1e-3)
  }
})

test_that("ns_lambda and ns_curvature_peak are inverses, and vectorised", {
  peaks <- c(0.5, 1, 2.5, 4, 10, 30)
  expect_equal(ns_curvature_peak(ns_lambda(peaks)), peaks)
  expect_equal(ns_lambda(ns_curvature_peak(peaks)), peaks)
  expect_length(ns_lambda(peaks), length(peaks))

  # Lambda is per unit of maturity, which is the whole point: the same curve
  # shape in months and in years differs by the unit conversion.
  expect_equal(ns_lambda(2.5) / ns_lambda(30), 12, tolerance = 1e-8)
})

test_that("ns_curvature_peak explains the package default", {
  # 0.0609 is the Diebold-Li value for maturities in months, where it peaks
  # just short of 30. The package measures maturities in years.
  expect_equal(ns_curvature_peak(0.0609), 29.4, tolerance = 0.01)
  expect_equal(ns_lambda(30), 0.0598, tolerance = 1e-3)

  # The default of fit_ns_model() is that same number.
  expect_equal(eval(formals(fit_ns_model)$initial_lambda), 0.0609)
})

test_that("ns_lambda rejects nonsense", {
  expect_error(ns_lambda(-1), "peak_maturity")
  expect_error(ns_lambda(NA_real_), "peak_maturity")
  expect_error(ns_curvature_peak(numeric(0)), "lambda")
})
