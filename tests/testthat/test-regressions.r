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

  # A tweaked sequential plan, which costs nothing to install and is still
  # distinguishable from the plain sequential one the bug used to force. That
  # is precisely what has to be detected, and it needs no workers: starting
  # those was most of this test's runtime.
  future::plan(future::sequential, split = TRUE)
  expected <- class(future::plan())
  expect_true("tweaked" %in% expected)

  panel <- ns_panel(14)
  invisible(fit_dns_model(panel, n_days = 10, parallel = FALSE,
                          show_progress = FALSE))
  expect_equal(class(future::plan()), expected)
})

test_that("parallel and sequential fits agree, and restore the plan", {
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  future::plan(future::multisession, workers = 2)
  expected <- class(future::plan())

  panel <- ns_panel(14)
  sequential_fit <- fit_dns_model(panel, n_days = 10, parallel = FALSE,
                                  show_progress = FALSE, warm_start = FALSE)
  parallel_fit <- fit_dns_model(panel, n_days = 10, parallel = TRUE,
                                show_progress = FALSE)

  expect_equal(class(future::plan()), expected)

  # C-07: the two paths agree, which they could not while the warm start was
  # carried through a global variable the workers never saw.
  expect_equal(sequential_fit$parameters$Date, parallel_fit$parameters$Date)
  expect_equal(sequential_fit$parameters$beta1, parallel_fit$parameters$beta1)
})


# ---- C-07 -------------------------------------------------------------------

test_that("fitting writes nothing into the global environment", {
  if (exists("last_fit_result", envir = globalenv(), inherits = FALSE)) {
    rm("last_fit_result", envir = globalenv())
  }
  invisible(fit_dns_model(ns_panel(12), n_days = 10, show_progress = FALSE))
  expect_false(exists("last_fit_result", envir = globalenv(), inherits = FALSE))
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

  # Which calendar date ends up at the top of the panel. Asserted through the
  # built y positions rather than the panel range: the range is an internal
  # that ggplot2 4.x normalises, and an earlier fix here reversed the axis on
  # 3.x while silently doing nothing on 4.x.
  top_date <- function(reverse) {
    built <- ggplot2::ggplot_build(
      plot_yield_heatmap(panel, reverse_time = reverse)
    )
    highest <- max(built$data[[1]]$y)
    as.Date(if (reverse) -highest else highest, origin = "1970-01-01")
  }

  # A date axis increases upward, so the default puts the most recent date at
  # the top and reverse_time flips that.
  expect_equal(top_date(FALSE), max(panel$Date))
  expect_equal(top_date(TRUE), min(panel$Date))
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

  # QMLE's cost grows with the series length and this asserts only that the
  # factor columns are picked out by name, so a short panel is enough.
  panel <- ns_panel(30)
  fit <- fit_dns_model(panel, n_days = 30, show_progress = FALSE)

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


# ---- plot rendering ---------------------------------------------------------

test_that("every plot actually draws, not just builds", {
  # ggplot_build() stops short of text measurement, so it cannot catch a label
  # the device is unable to render. That gap is how literal Greek in these
  # labels survived: on Windows the default device converts text through a
  # single-byte charset and failed with "conversion failure on 'lambda =
  # 0.0609' in mbcsToSbcs", which only surfaced once the examples were made
  # runnable and R CMD check started drawing them. Printing to a null device
  # exercises the whole pipeline.
  data(treasury_yields)
  panel <- treasury_yields[1:30, ]
  fit <- suppressWarnings(
    fit_dns_model(panel, n_days = 20, show_progress = FALSE)
  )
  # One decomposition, both of its plots.
  pca <- plot_yield_pca(panel, scale = FALSE)

  draws <- list(
    loadings = plot_dns_loadings(),
    loadings_plain = plot_dns_loadings(show_labels = FALSE),
    curves = plot_yield_curves(panel, fit, max_plots = 3),
    heatmap = plot_yield_heatmap(panel),
    heatmap_reversed = plot_yield_heatmap(panel, reverse_time = TRUE),
    pca_weights = pca$plots$weights,
    pca_variance = pca$plots$variance,
    comparison = suppressWarnings(
      compare_term_structure_models(panel, models = c("NS", "NSS", "Spline"))
    )$plot
  )

  factor_plots <- plot_factor_time_series(fit, show_changes = TRUE)
  draws$factors <- factor_plots$factors
  draws$changes <- factor_plots$changes
  draws$rmse <- factor_plots$rmse

  params <- list(kappa = c(0.02, 0.3, 0.5), theta = c(0.04, -0.02, 0.01),
                 sigma = c(0.005, 0.01, 0.02))
  simulation <- simulate_dns_yield_curves(
    params, c(0.05, -0.01, 0.01), maturities = c(1, 5, 10, 30),
    n_paths = 20, engine = "exact", seed = 1
  )
  draws$simulated_factors <- plot_simulated_factors(simulation)$level
  draws$simulated_yields <- plot_simulated_yields(
    simulation, time_points = c(0, 1, 5, 10)
  )$yield_curves

  grDevices::pdf(NULL)
  on.exit({
    while (grDevices::dev.cur() > 1) grDevices::dev.off()
  }, add = TRUE)

  for (name in names(draws)) {
    # expect_no_error() takes no info argument, so the plot's name is carried
    # in the label instead -- otherwise a failure here says only that some
    # plot among a dozen did not draw.
    outcome <- tryCatch({
      # The changes plot drops the leading NA of each change column, which
      # ggplot2 reports; that is expected and is not what this test is about.
      suppressWarnings(print(draws[[name]]))
      TRUE
    }, error = function(e) conditionMessage(e))

    expect_true(
      isTRUE(outcome),
      label = if (isTRUE(outcome)) {
        paste0(name, " draws")
      } else {
        paste0(name, " failed to draw: ", outcome)
      }
    )
  }
})


# ---- panel independence -----------------------------------------------------

test_that("without warm start, a date's fit is independent of the panel", {
  data(treasury_yields)
  # A short window on purpose: this asserts an invariant about how the fits are
  # initialised, which does not need a long panel, and each extra date here
  # costs four cross-sectional fits.
  window <- utils::tail(treasury_yields, 60)
  month_ends <- window[
    !duplicated(format(window$Date, "%Y-%m"), fromLast = TRUE),
  ]

  fit_subset <- suppressWarnings(fit_dns_model(
    month_ends, n_days = nrow(month_ends), show_progress = FALSE,
    warm_start = FALSE
  ))
  fit_full <- suppressWarnings(fit_dns_model(
    window, n_days = nrow(window), show_progress = FALSE, warm_start = FALSE
  ))

  shared <- match(fit_subset$parameters$Date, fit_full$parameters$Date)
  expect_false(anyNA(shared))

  # Identical, which is what lets the examples fit only the month-ends the
  # monthly dynamics estimation uses.
  for (column in c("beta1", "beta2", "beta3", "lambda", "rmse")) {
    expect_equal(fit_subset$parameters[[column]],
                 fit_full$parameters[[column]][shared],
                 info = column)
  }

  # Warm starting changes where each fit begins -- the previous date's lambda,
  # and the month-end chain is not the daily one -- but no longer changes where
  # it ends. While the optimiser's relative convergence test was degenerate it
  # stopped close enough to its start value for the chain to leak into the
  # answer, and these two disagreed. With the objective scaled so that test
  # means what it says, both searches reach the same optimum and warm starting
  # buys speed only. That is the property worth pinning.
  warm_subset <- suppressWarnings(fit_dns_model(
    month_ends, n_days = nrow(month_ends), show_progress = FALSE
  ))
  warm_full <- suppressWarnings(fit_dns_model(
    window, n_days = nrow(window), show_progress = FALSE
  ))
  warm_shared <- match(warm_subset$parameters$Date, warm_full$parameters$Date)
  expect_equal(warm_subset$parameters$lambda,
               warm_full$parameters$lambda[warm_shared])
})


# ---- bounds reach the panel fitter ------------------------------------------

test_that("fit_dns_model can widen the bounds its at_bound report complains about", {
  data(treasury_yields)
  month_ends <- treasury_yields[
    !duplicated(format(treasury_yields$Date, "%Y-%m"), fromLast = TRUE),
  ]

  # The subject here is whether bounds reach the cross-sectional fitter, so the
  # starting lambda is pinned rather than left to the default: the two
  # questions are separable, and a default that sits in the degenerate
  # small-lambda basin would otherwise decide this test's outcome.
  start <- ns_lambda(3.5)

  # Default ceiling on lambda is 0.5, and real curves reach it.
  expect_warning(
    clamped <- fit_dns_model(month_ends, n_days = nrow(month_ends),
                             show_progress = FALSE, warm_start = FALSE,
                             initial_lambda = start),
    "rest on a bound in some fits"
  )
  expect_true(any(clamped$at_bound[, "lambda"]))
  expect_true(all(clamped$parameters$lambda <= 0.5 + 1e-8))

  # Widening them clears it, which previously could not be done from here at
  # all: the report named a problem the function gave no way to act on.
  freed <- fit_dns_model(month_ends, n_days = nrow(month_ends),
                         show_progress = FALSE, warm_start = FALSE,
                         initial_lambda = start,
                         lower_bounds = c(-1, -1, -1, 0.01),
                         upper_bounds = c(1, 1, 1, 5))
  expect_false(any(freed$at_bound))
  expect_gt(max(freed$parameters$lambda), 0.5)

  # And the freed fit is better, which is why the clamp mattered.
  expect_lt(stats::median(freed$parameters$rmse),
            stats::median(clamped$parameters$rmse))
})

test_that("omitting the bounds leaves fit_ns_model's own defaults in force", {
  # fit_ns_model() widens its beta box to the data scale only when both bounds
  # are missing(), so forwarding its defaults explicitly would switch that off.
  # Percentage-point yields are the case that needs it.
  maturities <- c(0.25, 1, 2, 5, 10, 30)
  loading <- (1 - exp(-0.3 * maturities)) / (0.3 * maturities)
  dates <- seq(as.Date("2023-01-02"), by = "day", length.out = 6)

  panel <- data.frame(Date = dates)
  for (k in seq_along(maturities)) {
    panel[[as.character(maturities[k])]] <-
      (0.028 + 0.01 * loading[k]) * 100 + 1e-3 * seq_along(dates)
  }

  fit <- fit_dns_model(panel, n_days = nrow(panel), show_progress = FALSE)

  # Beta1 is around 2.8 on this scale, far outside the default box of
  # [-0.1, 0.2]; the widening is what keeps it off a bound.
  expect_false(any(fit$at_bound[, c("beta1", "beta2", "beta3")]))
  expect_gt(max(fit$parameters$beta1), 1)
})
