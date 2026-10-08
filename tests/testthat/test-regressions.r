# One test per defect fixed in 0.3.0. Each reproduces the original failure, so
# a regression shows up as a named failure rather than a subtle shift in
# calibration output.

# A Nelson-Siegel panel with a known generating lambda and a rising level, so
# the sign of the factor changes is unambiguous.
ns_panel <- function(n = 20, lambda = 0.3, drift = 4e-4,
                     maturities = c(0.25, 1, 2, 5, 10, 30),
                     start = as.Date("2023-01-02")) {
  dates <- seq(start, by = "day", length.out = n)
  out <- data.frame(Date = dates)
  loading <- (1 - exp(-lambda * maturities)) / (lambda * maturities)
  for (k in seq_along(maturities)) {
    out[[as.character(maturities[k])]] <-
      0.02 + drift * seq_len(n) + 0.01 * loading[k]
  }
  out
}

# As returned by get_treasury_yield_*(): most recent date first.
as_fred_order <- function(x) x[order(x$Date, decreasing = TRUE), , drop = FALSE]


# ---- C-01 -------------------------------------------------------------------

test_that("gaps are filled from the preceding date, not the following one", {
  locf <- getFromNamespace("locf_by_date", "dnsr")

  dates <- as.Date("2023-01-02") + 0:3
  x <- c(0.010, NA, NA, 0.040)

  # The only correct answer carries 0.010 forward. Taking 0.040 would mean
  # reading a value from a date that had not happened yet.
  expect_equal(locf(x, dates), c(0.010, 0.010, 0.010, 0.040))

  # A leading NA has nothing earlier to carry and must stay NA.
  expect_equal(locf(c(NA, 0.02, NA), dates[1:3]), c(NA, 0.02, 0.02))

  # max_gap is measured in calendar days.
  expect_equal(locf(x, dates, max_gap = 1), c(0.010, 0.010, NA, 0.040))
  expect_equal(locf(x, dates, max_gap = 0), x)

  # All-missing and empty inputs pass through.
  expect_equal(locf(rep(NA_real_, 3), dates[1:3]), rep(NA_real_, 3))
  expect_equal(locf(numeric(0), dates[0]), numeric(0))
})


# ---- C-02 -------------------------------------------------------------------

test_that("n_days selects the most recent dates whatever the row order", {
  panel <- ns_panel(20)
  recent <- utils::tail(sort(panel$Date), 5)

  for (input in list(panel, as_fred_order(panel))) {
    fit <- fit_dns_model(input, n_days = 5, show_progress = FALSE)
    expect_setequal(fit$parameters$Date, recent)
  }
})

test_that("date_order controls the row order of the result", {
  panel <- ns_panel(20)
  fred <- as_fred_order(panel)

  # "input" keeps the caller's direction, so parameters[1, ] stays the most
  # recent row for a most-recent-first frame. README and
  # get_calibrated_dns_params() both read it that way.
  fit <- fit_dns_model(fred, n_days = 5, show_progress = FALSE)
  expect_equal(fit$parameters$Date[1], max(fit$parameters$Date))

  asc <- fit_dns_model(fred, n_days = 5, show_progress = FALSE,
                       date_order = "ascending")
  expect_false(is.unsorted(asc$parameters$Date))

  desc <- fit_dns_model(panel, n_days = 5, show_progress = FALSE,
                        date_order = "descending")
  expect_false(is.unsorted(rev(desc$parameters$Date)))
})


# ---- C-03 -------------------------------------------------------------------

test_that("unfittable dates are dropped, not left as zero-factor rows", {
  panel <- ns_panel(20)
  panel[15, -1] <- NA
  bad_date <- panel$Date[15]

  fit <- suppressWarnings(
    fit_dns_model(panel, n_days = 10, show_progress = FALSE)
  )

  expect_false(bad_date %in% fit$parameters$Date)
  expect_equal(nrow(fit$parameters), 9)
  expect_equal(fit$success_rate, 90)

  # Nothing zeroed: a pre-sized numeric frame would leave beta = 0 here and
  # every is.na() filter downstream would pass it through.
  expect_false(any(fit$parameters$beta1 == 0))
  expect_equal(nrow(fit$failures), 1)
  expect_equal(fit$failures$Date, bad_date)

  kept <- suppressWarnings(
    fit_dns_model(panel, n_days = 10, show_progress = FALSE,
                  drop_failed = FALSE)
  )
  expect_equal(nrow(kept$parameters), 10)
  expect_true(is.na(kept$parameters$beta1[kept$parameters$Date == bad_date]))
})


# ---- C-04 -------------------------------------------------------------------

test_that("a single missing maturity does not abort the fit", {
  panel <- ns_panel(20)
  panel[["30"]][15] <- NA
  gap_date <- panel$Date[15]

  fit <- fit_dns_model(panel, n_days = 10, show_progress = FALSE)

  # Full-width row, with the unobserved maturity marked NA rather than the
  # fitted vector being written short.
  expect_equal(ncol(fit$fitted_yields), length(fit$maturities))
  row <- fit$fitted_yields[fit$parameters$Date == gap_date, ]
  expect_true(is.na(row[["30"]]))
  expect_false(anyNA(row[c("0.25", "1", "2", "5", "10")]))

  # Other dates are untouched.
  other <- fit$fitted_yields[fit$parameters$Date == panel$Date[14], ]
  expect_false(anyNA(other))
})


# ---- C-05 -------------------------------------------------------------------

test_that("factor changes are calendar-order differences", {
  panel <- ns_panel(20)

  from_asc <- fit_dns_model(panel, n_days = 10, show_progress = FALSE,
                            date_order = "ascending")
  from_desc <- fit_dns_model(as_fred_order(panel), n_days = 10,
                             show_progress = FALSE, date_order = "ascending")

  # Same panel, either input order, same changes.
  expect_equal(from_asc$parameters$beta1_change,
               from_desc$parameters$beta1_change)

  # change[t] is beta(t) - beta(t-1), with no leading value.
  expect_true(is.na(from_asc$parameters$beta1_change[1]))
  expect_equal(from_asc$parameters$beta1_change[-1],
               diff(from_asc$parameters$beta1))

  # Written back in the output's own row order.
  rev_out <- fit_dns_model(panel, n_days = 10, show_progress = FALSE,
                           date_order = "descending")
  expect_true(is.na(rev_out$parameters$beta1_change[nrow(rev_out$parameters)]))
  matched <- match(rev_out$parameters$Date, from_asc$parameters$Date)
  expect_equal(rev_out$parameters$beta1_change,
               from_asc$parameters$beta1_change[matched])
})


# ---- C-06 -------------------------------------------------------------------

test_that("fitting leaves the caller's future plan alone", {
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  future::plan(future::multisession, workers = 2)
  expected <- class(future::plan())

  panel <- ns_panel(12)
  invisible(fit_dns_model(panel, n_days = 10, parallel = FALSE,
                          show_progress = FALSE))
  expect_equal(class(future::plan()), expected)

  invisible(fit_dns_model(panel, n_days = 10, parallel = TRUE,
                          show_progress = FALSE))
  expect_equal(class(future::plan()), expected)
})


# ---- C-07 -------------------------------------------------------------------

test_that("fitting writes nothing into the global environment", {
  if (exists("last_fit_result", envir = globalenv(), inherits = FALSE)) {
    rm("last_fit_result", envir = globalenv())
  }
  invisible(fit_dns_model(ns_panel(12), n_days = 10, show_progress = FALSE))
  expect_false(exists("last_fit_result", envir = globalenv(), inherits = FALSE))
})

test_that("parallel and sequential fits agree", {
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  panel <- ns_panel(14)
  seq_fit <- fit_dns_model(panel, n_days = 10, parallel = FALSE,
                           show_progress = FALSE, warm_start = FALSE)
  par_fit <- fit_dns_model(panel, n_days = 10, parallel = TRUE,
                           show_progress = FALSE)

  expect_equal(seq_fit$parameters$Date, par_fit$parameters$Date)
  expect_equal(seq_fit$parameters$beta1, par_fit$parameters$beta1)
})


# ---- H-01 -------------------------------------------------------------------

test_that("model comparison errors are genuinely in basis points", {
  panel <- ns_panel(20)
  curve <- as.numeric(panel[nrow(panel), -1])
  maturities <- as.numeric(setdiff(names(panel), "Date"))
  reference <- fit_ns_model(maturities, curve)$rmse

  bp <- suppressWarnings(
    compare_term_structure_models(panel, models = "NS")
  )
  expect_equal(bp$metrics$RMSE, reference * 10000)

  dec <- suppressWarnings(
    compare_term_structure_models(panel, models = "NS", units = "decimal")
  )
  expect_equal(dec$metrics$RMSE, reference)
})

test_that("model comparison survives an unobserved maturity", {
  panel <- ns_panel(20)
  panel[["30"]][nrow(panel)] <- NA

  expect_warning(
    cmp <- compare_term_structure_models(panel, models = c("NS", "Spline")),
    "no observation"
  )
  expect_equal(nrow(cmp$metrics), 2)
  expect_true(all(is.finite(cmp$metrics$RMSE)))
})


# ---- H-02 -------------------------------------------------------------------

test_that("sim_method reaches the simulator", {
  expect_true("sim_method" %in% names(formals(simulate_dns_factors)))

  skip_if_not_installed("yuima")
  params <- list(kappa = c(0.02, 0.3, 0.5), theta = c(0.04, -0.02, 0.01),
                 sigma = c(0.005, 0.01, 0.02))

  sim <- suppressWarnings(simulate_dns_yield_curves(
    params, c(0.035, -0.015, 0.005), maturities = c(1, 5, 10),
    n_steps = 10, n_paths = 2, sim_method = "milstein", seed = 1
  ))
  expect_equal(sim$sim_method, "milstein")
})

test_that("simulation does not print debug diagnostics", {
  skip_if_not_installed("yuima")
  params <- list(kappa = c(0.02, 0.3, 0.5), theta = c(0.04, -0.02, 0.01),
                 sigma = c(0.005, 0.01, 0.02))

  msgs <- capture.output(
    invisible(suppressWarnings(simulate_dns_factors(
      params, c(0.035, -0.015, 0.005), n_steps = 5, n_paths = 2, seed = 1
    ))),
    type = "message"
  )
  expect_false(any(grepl("Simulation data", msgs)))
})


# ---- H-03 -------------------------------------------------------------------

test_that("yield curve facets are labelled with dates, in calendar order", {
  panel <- ns_panel(20)
  fit <- fit_dns_model(panel, n_days = 10, show_progress = FALSE)

  p <- plot_yield_curves(panel, fit, max_plots = 3)
  labels <- levels(p$data$date)

  # Iterating a Date vector drops the class, so these used to come out as day
  # numbers ("19369").
  expect_true(all(grepl("^\\d{4}-\\d{2}-\\d{2}$", labels)))
  expect_equal(labels, sort(labels))
  expect_true(all(as.Date(labels) %in% panel$Date))
})


# ---- H-04 -------------------------------------------------------------------

test_that("plot_simulated_factors handles a single path", {
  skip_if_not_installed("yuima")
  params <- list(kappa = c(0.02, 0.3, 0.5), theta = c(0.04, -0.02, 0.01),
                 sigma = c(0.005, 0.01, 0.02))
  sim <- suppressWarnings(simulate_dns_factors(
    params, c(0.035, -0.015, 0.005), n_steps = 10, n_paths = 1, seed = 1
  ))
  expect_no_error(plot_simulated_factors(sim))
})

test_that("plot_yield_heatmap reverse_time reverses the date axis", {
  panel <- ns_panel(20)

  normal <- ggplot2::ggplot_build(
    plot_yield_heatmap(panel, reverse_time = FALSE)
  )$layout$panel_params[[1]]$y.range
  reversed <- ggplot2::ggplot_build(
    plot_yield_heatmap(panel, reverse_time = TRUE)
  )$layout$panel_params[[1]]$y.range

  expect_lt(normal[1], normal[2])
  expect_gt(reversed[1], reversed[2])
})

test_that("plot_yield_pca returns var_explained with every plot flag off", {
  panel <- ns_panel(20)
  res <- plot_yield_pca(panel, plot_weights = FALSE, plot_scores = FALSE,
                        plot_variance = FALSE)

  expect_length(res$plots, 0)
  expect_true(length(res$var_explained) > 0)
  expect_equal(sum(res$var_explained), 1)
})


# ---- H-05 -------------------------------------------------------------------

test_that("simulated yield facets follow time, not the alphabet", {
  skip_if_not_installed("yuima")
  params <- list(kappa = c(0.02, 0.3, 0.5), theta = c(0.04, -0.02, 0.01),
                 sigma = c(0.005, 0.01, 0.02))
  sim <- suppressWarnings(simulate_dns_yield_curves(
    params, c(0.035, -0.015, 0.005), maturities = c(1, 5, 10),
    time_horizon = 10, n_steps = 20, n_paths = 3, seed = 1
  ))

  plots <- plot_simulated_yields(sim, time_points = c(0, 1, 5, 10),
                                 show_forwards = TRUE)

  # Character labels sort "Year 10" between "Year 1" and "Year 5".
  expect_equal(levels(plots$yield_curves$data$Time),
               c("Year 0", "Year 1", "Year 5", "Year 10"))
  expect_equal(levels(plots$forwards$data$Time),
               c("Year 0", "Year 1", "Year 5", "Year 10"))
})


# ---- H-06 -------------------------------------------------------------------

test_that("a fitted lambda can reach the simulation parameters", {
  dynamics <- list(
    beta1 = list(kappa = 1, theta = 0.04, sigma = 0.01),
    beta2 = list(kappa = 1, theta = 0.00, sigma = 0.01),
    beta3 = list(kappa = 1, theta = 0.00, sigma = 0.01),
    correlation = diag(3)
  )

  expect_equal(dns_dynamics_to_sim_params(dynamics, lambda = 0.09)$lambda, 0.09)

  # Unchanged default, and a lambda carried on the dynamics object is honoured.
  expect_equal(dns_dynamics_to_sim_params(dynamics)$lambda, 0.0609)
  dynamics$lambda <- 0.07
  expect_equal(dns_dynamics_to_sim_params(dynamics)$lambda, 0.07)
})


# ---- H-07 -------------------------------------------------------------------

test_that("default bounds cover the data scale instead of clamping", {
  maturities <- c(0.25, 1, 2, 5, 10, 30)
  loading <- (1 - exp(-0.3 * maturities)) / (0.3 * maturities)
  decimal <- 0.028 + 0.01 * loading

  # Percentage-point yields used to pin all four parameters at their bounds
  # (RMSE 2.7 percentage points) and still report success = TRUE.
  pct <- fit_ns_model(maturities, decimal * 100)
  expect_false(any(pct$at_bound))
  expect_true(pct$success)
  expect_lt(pct$rmse, 0.01)
  expect_equal(pct$parameters[["lambda"]], 0.3, tolerance = 0.1)
})

test_that("a boundary solution is reported rather than passed off as success", {
  maturities <- c(0.25, 1, 2, 5, 10, 30)
  loading <- (1 - exp(-0.3 * maturities)) / (0.3 * maturities)
  decimal <- 0.028 + 0.01 * loading

  expect_warning(
    clamped <- fit_ns_model(maturities, decimal,
                            lower_bounds = c(0.035, -0.5, -0.5, 0.01),
                            upper_bounds = c(0.200, 0.5, 0.5, 0.50)),
    "rests on a parameter bound"
  )
  expect_true(clamped$at_bound[["beta1"]])
  expect_equal(clamped$parameters[["beta1"]], 0.035)
})


# ---- H-08 -------------------------------------------------------------------

test_that("estimate_dns_parameters_yuima selects factor columns by name", {
  skip_if_not_installed("yuima")

  panel <- ns_panel(60)
  fit <- fit_dns_model(panel, n_days = 60, show_progress = FALSE)

  # Passing the parameters frame straight in used to feed the Date column into
  # diff() as factor 1.
  expect_no_error(
    suppressWarnings(estimate_dns_parameters_yuima(fit$parameters,
                                                   delta_t = 1 / 252))
  )
})

test_that("get_calibrated_dns_params takes an api_key", {
  args <- names(formals(get_calibrated_dns_params))
  expect_true("api_key" %in% args)

  # Appended, so existing positional calls still land on the same arguments.
  expect_equal(args[1:5],
               c("lookback_days", "lambda", "start_date", "end_date", "verbose"))
})
