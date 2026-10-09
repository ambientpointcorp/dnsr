#' @title Dynamic Nelson-Siegel Model Implementation
#' @description Functions to calibrate and fit the Dynamic Nelson-Siegel (DNS) model to yield curve data.

#' @importFrom stats optim
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
#' @importFrom progress progress_bar

#' @title Fit Nelson-Siegel Model to a Single Yield Curve
#'
#' @description
#' Fits the Nelson-Siegel model to a single yield curve (cross-section) for a specific date.
#'
#' @param maturities Numeric vector. The maturities (in years) for the yield curve points.
#' @param yields Numeric vector. The yield values corresponding to the maturities.
#' @param initial_lambda Numeric. Initial value for the decay parameter lambda,
#'   per unit of maturity, and the value used outright when
#'   \code{lambda_constraint = "fixed"}. Maturities here are in years, and the
#'   default places the curvature hump at 3.5 years, near the middle of a
#'   Treasury curve; it is written as \code{ns_lambda(3.5)} so the number
#'   states its own meaning. Diebold and Li's 0.0609 is a month-denominated
#'   constant: in years it places the hump at 29.4 years, off the end of most
#'   curves, where the slope loading is nearly parallel to the level loading
#'   and the two factors are not separately identified. See
#'   \code{\link{ns_lambda}}.
#' @param lambda_constraint Character. One of "fixed", "constrained", or "free":
#'   - "fixed": lambda is fixed at the initial_lambda value.
#'   - "constrained": lambda is optimized within bounds.
#'   - "free": lambda is optimized without constraints.
#'
#'   With \code{method = "optim"} the search over lambda is weak: L-BFGS-B's
#'   numerical gradient is poorly scaled when the betas are around 0.05 and
#'   lambda around 0.5, and in practice lambda moves very little from
#'   \code{initial_lambda} even when unconstrained. \code{method = "profile"}
#'   searches lambda properly and reaches the same optimum from any starting
#'   value, so prefer it when lambda is meant to be estimated rather than
#'   supplied.
#' @param lower_bounds Numeric vector of length 4. Lower bounds for parameters (beta1, beta2, beta3, lambda).
#' @param upper_bounds Numeric vector of length 4. Upper bounds for parameters
#'   (beta1, beta2, beta3, lambda). The lambda ceiling is 2, which places the
#'   curvature hump at 0.9 years; a ceiling of 0.5 bound the fit on roughly a
#'   third of real month-end curves, and a clamped lambda is reported through
#'   \code{at_bound} rather than silently. Raising it further buys in-sample
#'   fit at the cost of comparability across dates, since factors estimated at
#'   different lambda are not on the same scale.
#' @param optimizer Character. The optimization method to use. Default is
#'   "L-BFGS-B". Used by \code{method = "optim"} only.
#' @param control List. Passed to \code{\link[stats]{optim}}'s \code{control}
#'   argument by \code{method = "optim"}; ignored by \code{method = "profile"},
#'   which uses \code{\link[stats]{optimize}}. The objective is minimised in
#'   squared basis points so that L-BFGS-B's relative convergence test is
#'   genuinely relative, which means \code{factr} keeps its documented meaning
#'   here and rarely needs changing.
#' @param method Character. How to solve the fit. \code{"optim"} (the default)
#'   searches all four parameters numerically, as in earlier versions.
#'   \code{"profile"} exploits the fact that the curve is linear in the betas
#'   once lambda is fixed: the betas are solved by least squares in closed form
#'   and only lambda is searched, over one dimension. The profiled fit does not
#'   depend on beta start values and does not need beta bounds, which removes
#'   the boundary solutions \code{"optim"} can return; it is the recommended
#'   setting and is expected to become the default.
#'
#' @return A list containing:
#'   - parameters: Named vector of fitted parameters (beta1, beta2, beta3, lambda)
#'   - fitted_yields: Vector of fitted yield values, one per retained maturity
#'   - rmse: Root mean squared error of the fit
#'   - maturities: The retained maturities (input maturities, minus any dropped
#'     for being missing or non-positive)
#'   - actual_yields: The retained input yields
#'   - residuals: The fitting residuals (actual - fitted)
#'   - success: Logical indicating if the optimization converged
#'   - valid_idx: Logical vector over the *input* maturities marking which
#'     points were retained, so callers can place the fitted values back into a
#'     full-width row
#'   - at_bound: Named logical vector marking parameters that came to rest on
#'     one of their bounds
#'
#' @details
#' The Nelson-Siegel model expresses the yield curve as:
#' \deqn{y(m) = \beta_1 + \beta_2 \left(\frac{1-e^{-\lambda m}}{\lambda m}\right) + \beta_3 \left(\frac{1-e^{-\lambda m}}{\lambda m} - e^{-\lambda m}\right)}
#' where \eqn{m} is maturity, \eqn{\beta_1}, \eqn{\beta_2}, and \eqn{\beta_3} represent the level, slope, and curvature
#' factors, respectively, and \eqn{\lambda} controls the decay rate of factor loadings.
#'
#' With \code{method = "profile"} the sum of squared errors is concentrated:
#' for any lambda the optimal betas are the least-squares solution of
#' \eqn{y = X(\lambda)\beta}, so the objective reduces to a function of lambda
#' alone, minimised over \code{[lower_bounds[4], upper_bounds[4]]} by
#' \code{\link[stats]{optimize}}. Beta bounds are not consulted, since the
#' linear solve is exact.
#'
#' The default bounds are written for yields in decimal form. When they are
#' left at their defaults and the data is evidently on another scale (yields in
#' percentage points, say), the beta box is widened to cover the data rather
#' than clamping the level factor: L-BFGS-B reports convergence at a bound, so
#' a clamped fit would otherwise be returned with \code{success = TRUE}. Any
#' parameter that does come to rest on a bound is reported in \code{at_bound}
#' and raises a warning.
#'
#' @examples
#' maturities <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)
#' yields <- c(0.0050, 0.0052, 0.0055, 0.0065, 0.0078, 0.0105, 0.0128, 0.0152, 0.0215, 0.0225)
#' fit_ns_model(maturities, yields)
#'
#' # The profiled fit solves the betas in closed form at each lambda, so it
#' # needs no start values and cannot return a clamped beta.
#' fit_ns_model(maturities, yields, method = "profile")
#'
#' # A real curve from the bundled panel
#' data(treasury_yields)
#' latest <- treasury_yields[nrow(treasury_yields), ]
#' fit_ns_model(as.numeric(setdiff(names(treasury_yields), "Date")),
#'              as.numeric(latest[-1]), method = "profile")
#'
#' @export
fit_ns_model <- function(maturities, yields,
                         initial_lambda = ns_lambda(3.5),
                         lambda_constraint = "constrained",
                         lower_bounds = c(-0.1, -0.5, -0.5, 0.01),
                         upper_bounds = c(0.2, 0.5, 0.5, 2),
                         optimizer = "L-BFGS-B",
                         method = c("optim", "profile"),
                         control = list()) {
  method <- match.arg(method)

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
  checkmate::assert_list(control, names = "unique")

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
      success = FALSE,
      valid_idx = valid_idx,
      at_bound = c(beta1 = NA, beta2 = NA, beta3 = NA, lambda = NA),
      message = "Fewer than 3 valid yield observations"
    ))
  }

  maturities <- maturities[valid_idx]
  yields <- yields[valid_idx]

  # Function to calculate Nelson-Siegel yields
  ns_yields <- function(params, maturities) {
    beta1 <- params[1]
    beta2 <- params[2]
    beta3 <- params[3]
    lambda <- params[4]

    loadings <- ns_loadings(maturities, lambda)
    beta1 + beta2 * loadings[, 1] + beta3 * loadings[, 2]
  }

  # Error function to minimize (sum of squared errors), scaled for the
  # optimiser only.
  #
  # L-BFGS-B stops when the achieved reduction in f, divided by
  # max(|f_old|, |f_new|, 1), falls below factr * .Machine$double.eps. A sum of
  # squared errors in decimal yields is around 1e-5, so that denominator is the
  # constant 1 and the "relative" test degenerates into an absolute one of
  # roughly 2.2e-9 -- larger than most of the progress still available. The
  # search then reports CONVERGENCE after a handful of evaluations without
  # having moved off its start value.
  #
  # Working in squared basis points puts the objective above 1, where the test
  # measures what it is documented to measure. A positive constant does not
  # move the minimiser, so the parameters, rmse and residuals below are those
  # of the unscaled problem; only the stopping decision changes.
  SSE_SCALE <- 1e8

  error_fn <- function(params) {
    # Calculate model yields
    model_yields <- ns_yields(params, maturities)

    # Return sum of squared errors
    sum((yields - model_yields)^2) * SSE_SCALE
  }

  # Initial parameter guess. The slope start value is short minus long, taken
  # by maturity rather than by position, so an unsorted maturity vector does
  # not invert it.
  short_i <- which.min(maturities)
  long_i <- which.max(maturities)
  initial_params <- c(
    beta1 = mean(yields[maturities >= 5], na.rm = TRUE),  # Level (avg of long-term yields)
    beta2 = yields[short_i] - yields[long_i],             # Slope (short - long)
    beta3 = 0,                                           # Curvature
    lambda = initial_lambda                              # Decay rate
  )

  # Handle NAs in initial guess
  if (is.na(initial_params[1])) initial_params[1] <- mean(yields)
  if (is.na(initial_params[2])) initial_params[2] <- 0

  # The default beta box assumes decimal yields. If it was left at the default
  # and cannot contain this data, widen it: L-BFGS-B returns convergence == 0
  # at a bound, so a box that is too tight yields a clamped level factor
  # reported as a successful fit.
  if (missing(lower_bounds) && missing(upper_bounds)) {
    yield_scale <- max(abs(yields), na.rm = TRUE)
    if (is.finite(yield_scale) && yield_scale > 0) {
      lower_bounds[1:3] <- pmin(lower_bounds[1:3], -2 * yield_scale)
      upper_bounds[1:3] <- pmax(upper_bounds[1:3],  2 * yield_scale)
    }
  }

  # Keep the start values inside the box; L-BFGS-B errors outright otherwise.
  if (lambda_constraint != "free") {
    initial_params <- pmin(pmax(initial_params, lower_bounds), upper_bounds)
  }

  if (method == "profile") {
    # Concentrated least squares: solve the betas exactly at each lambda and
    # search only over lambda.
    profile_sse <- function(lambda) {
      design <- ns_design(maturities, lambda)
      sum(qr.resid(qr(design), yields)^2)
    }

    interval <- c(lower_bounds[4], upper_bounds[4])
    if (lambda_constraint == "fixed") {
      lambda_hat <- initial_lambda
    } else {
      # "free" has no box to search, so use the range initial_lambda itself is
      # validated against.
      interval <- if (lambda_constraint == "constrained") {
        c(lower_bounds[4], upper_bounds[4])
      } else {
        c(0.001, 5)
      }
      # Grid first, then refine. The profiled objective is not unimodal in
      # lambda: it has a local basin near zero, where the level and slope
      # loadings both tend to 1, separated from the basin a real curve wants by
      # a ridge. optimize() is a local bracketing search, so on its own it
      # returns whichever basin its first three probes happen to straddle --
      # which is how it came to miss the optimum of an exactly representable
      # curve once the search interval was widened. A log-spaced scan locates
      # the right basin, and optimize() then resolves it between the
      # neighbouring grid points, which bracket the minimum.
      #
      # The scan is log-spaced because lambda acts as a reciprocal maturity: a
      # step from 0.02 to 0.04 moves the curvature hump by 45 years, the same
      # step from 1 to 1.02 moves it by half a month.
      grid <- exp(seq(log(interval[1]), log(interval[2]),
                      length.out = PROFILE_GRID_N))
      grid_sse <- vapply(grid, profile_sse, numeric(1))
      best <- which.min(grid_sse)

      # A tolerance well inside the edge test below, so a lambda that wants to
      # sit on a bound is resolved closely enough to be recognised as such.
      lambda_hat <- stats::optimize(
        profile_sse,
        interval = c(grid[max(1L, best - 1L)],
                     grid[min(length(grid), best + 1L)]),
        tol = 1e-9
      )$minimum

      # The scan is evaluated at the endpoints, so it can see a minimum that
      # the refinement, confined to the interior of a bracket, cannot reach.
      if (grid_sse[best] < profile_sse(lambda_hat)) lambda_hat <- grid[best]
    }

    design <- ns_design(maturities, lambda_hat)
    qr_design <- qr(design)
    if (qr_design$rank < 3) {
      warning("Nelson-Siegel design matrix is rank deficient at lambda = ",
              format(lambda_hat), "; the factors are not separately identified.")
    }
    betas <- qr.coef(qr_design, yields)
    betas[is.na(betas)] <- 0

    params <- c(betas, lambda_hat)
    names(params) <- c("beta1", "beta2", "beta3", "lambda")

    # Shaped like an optim() result so the shared code below needs no branch.
    opt_result <- list(par = params, convergence = 0L,
                       value = sum(qr.resid(qr_design, yields)^2))

    # Only lambda could have come to rest on a bound here; the betas are
    # unconstrained by construction. optimize() stops within its tolerance of
    # an endpoint rather than exactly on it, so the test is against a window
    # scaled to the search interval, not an equality.
    lambda_at_bound <- FALSE
    if (lambda_constraint == "constrained") {
      edge_tol <- 1e-6 * diff(interval)
      lambda_at_bound <- (lambda_hat - interval[1]) <= edge_tol ||
        (interval[2] - lambda_hat) <= edge_tol
    }

  } else if (lambda_constraint == "fixed") {
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
      upper = upper_bounds[1:3],
      control = control
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
      upper = upper_bounds,
      control = control
    )

    # Extract optimized parameters
    params <- opt_result$par
    names(params) <- c("beta1", "beta2", "beta3", "lambda")

  } else if (lambda_constraint == "free") {
    # Optimize without bounds. BFGS stays the default here because L-BFGS-B,
    # the function-wide default, is a bounded method; an explicitly supplied
    # optimizer is honoured.
    opt_result <- optim(
      par = initial_params,
      fn = error_fn,
      method = if (missing(optimizer)) "BFGS" else optimizer,
      control = control
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

  # Flag parameters resting on a bound. L-BFGS-B reports convergence == 0 in
  # that case, so `success` alone cannot distinguish a genuine optimum from a
  # clamped one.
  at_bound <- c(beta1 = FALSE, beta2 = FALSE, beta3 = FALSE, lambda = FALSE)
  if (method == "profile") {
    at_bound[["lambda"]] <- lambda_at_bound
    if (lambda_at_bound) {
      warning("Nelson-Siegel fit rests on a parameter bound (lambda); the ",
              "factor decomposition is not identified. Widen ",
              "lower_bounds/upper_bounds or check the yield scale.")
    }
  } else if (lambda_constraint != "free") {
    n_free <- if (lambda_constraint == "fixed") 3L else 4L
    tol <- 1e-8 * pmax(1, abs(upper_bounds - lower_bounds))
    idx <- seq_len(n_free)
    at_bound[idx] <- abs(params[idx] - lower_bounds[idx]) <= tol[idx] |
                     abs(params[idx] - upper_bounds[idx]) <= tol[idx]
    if (any(at_bound)) {
      warning("Nelson-Siegel fit rests on a parameter bound (",
              paste(names(at_bound)[at_bound], collapse = ", "),
              "); the factor decomposition is not identified. Widen ",
              "lower_bounds/upper_bounds or check the yield scale.")
    }
  }

  # Prepare return values
  result <- list(
    parameters = params,
    fitted_yields = fitted_yields,
    rmse = rmse,
    maturities = maturities,
    actual_yields = yields,
    residuals = residuals,
    success = opt_result$convergence == 0,
    valid_idx = valid_idx,
    at_bound = at_bound
  )

  # A class for display only: every element above keeps its name, and the
  # object is still a plain list to is.list(), `$` and checkmate.
  class(result) <- c("ns_fit", "list")

  return(result)
}

#' @title Fit Dynamic Nelson-Siegel Model to Time Series of Yield Curves
#'
#' @description
#' Fits the DNS model to a time series of yield curves.
#'
#' @param yield_data Data frame. The yield curve data with a 'Date' column and
#'   remaining columns containing yield values for different maturities, named
#'   by maturity in years.
#' @param n_days Integer. Number of most recent dates to include in the fitting.
#'   The window is selected by date, so it does not depend on the row order of
#'   \code{yield_data}.
#' @param parallel Logical. If TRUE, use parallel processing to speed up computation.
#' @param lambda_constraint Character. How to handle the lambda parameter (fixed, constrained, free).
#' @param initial_lambda Numeric. Initial value for lambda, passed to
#'   \code{\link{fit_ns_model}}. The default places the curvature hump at 3.5
#'   years; see \code{\link{ns_lambda}} on its units.
#' @param show_progress Logical. If TRUE, display a progress bar.
#' @param date_order Character. Row order of the returned \code{parameters} and
#'   \code{fitted_yields}: \code{"input"} (the default) keeps the order the
#'   selected rows had in \code{yield_data}, so a most-recent-first frame yields
#'   a most-recent-first result; \code{"ascending"} and \code{"descending"} sort
#'   explicitly.
#' @param drop_failed Logical. If TRUE (the default) dates whose cross-section
#'   could not be fitted are removed from the result. If FALSE they are kept
#'   with \code{NA} parameters.
#' @param warm_start Logical. If TRUE (the default) each cross-section is
#'   initialised from the previous date's fitted lambda, walking forward in
#'   calendar time. Ignored when \code{parallel = TRUE}, where the fits are not
#'   ordered, and by \code{method = "profile"}, which needs no start values.
#'
#'   Note that warm starting makes a date's fit depend on which \emph{other}
#'   dates are in the panel: fitting a month-end series directly and fitting
#'   the daily series that contains it give slightly different parameters on
#'   the shared dates, because the chain of start values differs. Set
#'   \code{warm_start = FALSE} for cross-sections that are independent of the
#'   panel they arrived in.
#' @param method Character. Passed to \code{\link{fit_ns_model}}:
#'   \code{"optim"} (the default) or \code{"profile"}.
#' @param lower_bounds,upper_bounds Numeric vectors of length 4, passed to
#'   \code{\link{fit_ns_model}} as the bounds on
#'   \code{(beta1, beta2, beta3, lambda)}. \code{NULL}, the default, leaves
#'   that function's own defaults in force, including its widening of the beta
#'   box to the data scale -- which it only applies when both bounds are left
#'   alone.
#'
#'   Supply them when the \code{at_bound} report says a fit is clamped. The
#'   default ceiling on lambda is 0.5, and real Treasury curves reach it: on
#'   the panel in \code{\link{treasury_yields}} a freely estimated lambda
#'   rests there on roughly a quarter of dates. Without these arguments there
#'   was no way to act on that from here.
#'
#' @return A list containing:
#'   - parameters: Data frame with model parameters for each date, plus
#'     \code{beta1_change}, \code{beta2_change} and \code{beta3_change}.
#'   - fitted_yields: Matrix of fitted yield values, one column per maturity
#'     column of \code{yield_data}. Maturities that were missing on a given
#'     date are \code{NA}.
#'   - maturities: Vector of maturities in years.
#'   - success_rate: Percentage of successful model fits.
#'   - failures: Data frame of dates that could not be fitted, with the reason.
#'   - at_bound: Logical matrix, one row per retained date and one column per
#'     parameter, marking fits that came to rest on a bound. Those dates'
#'     factor decompositions are not identified.
#'
#' @details
#' The factor change columns are first differences in calendar order
#' (\eqn{\beta(t) - \beta(t-1)}) regardless of the row order of the output, so
#' they can be fed directly to \code{\link{estimate_dns_dynamics}}.
#'
#' @examples
#' data(treasury_yields)
#'
#' # The 60 most recent curves. The window is chosen by date, so passing the
#' # panel in either row order selects the same dates.
#' dns_result <- fit_dns_model(treasury_yields, n_days = 60, show_progress = FALSE)
#' dns_result
#' summary(dns_result)
#'
#' head(dns_result$parameters)
#'
#' @export
fit_dns_model <- function(yield_data, n_days = 60, parallel = FALSE,
                          lambda_constraint = "constrained",
                          initial_lambda = ns_lambda(3.5),
                          show_progress = TRUE,
                          date_order = c("input", "ascending", "descending"),
                          drop_failed = TRUE,
                          warm_start = TRUE,
                          method = c("optim", "profile"),
                          lower_bounds = NULL, upper_bounds = NULL) {

  # Parameter validation
  date_order <- match.arg(date_order)
  method <- match.arg(method)
  if (!is.null(lower_bounds)) checkmate::assert_numeric(lower_bounds, len = 4)
  if (!is.null(upper_bounds)) checkmate::assert_numeric(upper_bounds, len = 4)
  checkmate::assert_data_frame(yield_data, min.rows = 1)
  checkmate::assert_names(names(yield_data), must.include = "Date")
  checkmate::assert_numeric(n_days, lower = 1, len = 1)
  checkmate::assert_logical(parallel, len = 1)
  checkmate::assert_logical(show_progress, len = 1)
  checkmate::assert_logical(drop_failed, len = 1)
  checkmate::assert_logical(warm_start, len = 1)

  # Locate the maturity columns by name rather than assuming Date is first.
  mat_cols <- setdiff(names(yield_data), "Date")
  if (length(mat_cols) == 0) {
    stop("yield_data has no maturity columns besides 'Date'")
  }
  maturities <- suppressWarnings(as.numeric(mat_cols))
  if (all(is.na(maturities))) {
    stop("No maturity column names of yield_data could be read as numbers. ",
         "Columns other than 'Date' must be named by maturity in years, ",
         "e.g. \"0.25\", \"2\", \"10\".")
  }

  dates_all <- as.Date(yield_data$Date)

  # Select the n_days most recent dates BY DATE. Taking the last n rows is only
  # "most recent" for an ascending frame, and this package's own data functions
  # return most-recent-first.
  n_keep <- min(as.integer(n_days), nrow(yield_data))
  by_recency <- order(dates_all, decreasing = TRUE)
  keep_rows <- sort(by_recency[seq_len(n_keep)])

  # Emit rows in the requested order.
  panel <- yield_data[keep_rows, , drop = FALSE]
  panel_dates <- dates_all[keep_rows]
  out_order <- switch(
    date_order,
    input      = seq_along(panel_dates),
    ascending  = order(panel_dates),
    descending = order(panel_dates, decreasing = TRUE)
  )
  panel <- panel[out_order, , drop = FALSE]
  dates <- panel_dates[out_order]

  # Yields as a plain numeric matrix, rows aligned with `dates`.
  yield_matrix <- as.matrix(
    vapply(mat_cols, function(cn) as.numeric(panel[[cn]]),
           numeric(nrow(panel)))
  )
  dim(yield_matrix) <- c(nrow(panel), length(mat_cols))
  colnames(yield_matrix) <- mat_cols

  # Setup parallel processing if requested. The plan is always restored: a
  # fitting call must not leave the caller's future plan changed.
  if (parallel) {
    oplan <- future::plan(future::multisession)
    on.exit(future::plan(oplan), add = TRUE)
  }

  fit_one <- function(i, lambda_init) {
    yields <- yield_matrix[i, ]

    if (sum(!is.na(yields)) < 3) {
      return(list(success = FALSE, message = "Not enough valid yield data"))
    }

    # Bounds are forwarded only when the caller supplied them. fit_ns_model()
    # decides whether to widen its beta box by asking whether both bounds are
    # missing(), so passing its own defaults through explicitly would silently
    # switch that off.
    args <- list(
      maturities = maturities,
      yields = yields,
      initial_lambda = lambda_init,
      lambda_constraint = lambda_constraint,
      method = method
    )
    if (!is.null(lower_bounds)) args$lower_bounds <- lower_bounds
    if (!is.null(upper_bounds)) args$upper_bounds <- upper_bounds

    tryCatch(
      # The per-date bound warning is muffled here and reported once for the
      # whole panel below: a run over several hundred dates can hit a bound on
      # many of them, and one warning each would bury everything else. The
      # information is not lost -- it is in each fit's at_bound element.
      withCallingHandlers(
        do.call(fit_ns_model, args),
        warning = function(w) {
          if (grepl("rests on a parameter bound", conditionMessage(w))) {
            invokeRestart("muffleWarning")
          }
        }
      ),
      error = function(e) list(success = FALSE, message = conditionMessage(e))
    )
  }

  pb <- NULL
  if (show_progress) {
    pb <- progress::progress_bar$new(
      format = "[:bar] :percent Fitted: :current/:total (:eta remaining)",
      total = length(dates)
    )
  }

  if (parallel) {
    # No warm start here: the fits are independent and unordered.
    fit_results <- future.apply::future_lapply(
      seq_along(dates),
      function(i) fit_one(i, initial_lambda),
      future.seed = TRUE
    )
    if (!is.null(pb)) pb$update(1)
  } else {
    # Sequential: walk forward in calendar time so the warm start carries the
    # previous date's lambda, in a local variable rather than a global one.
    calendar_order <- order(dates)
    fit_results <- vector("list", length(dates))
    lambda_init <- initial_lambda
    for (i in calendar_order) {
      res <- fit_one(i, lambda_init)
      if (warm_start && isTRUE(res$success) &&
          !is.na(res$parameters[["lambda"]])) {
        lambda_init <- res$parameters[["lambda"]]
      }
      fit_results[[i]] <- res
      if (!is.null(pb)) pb$tick()
    }
  }

  # ---- assemble -------------------------------------------------------------
  n <- length(dates)
  params_df <- data.frame(
    Date = dates,
    beta1 = rep(NA_real_, n),
    beta2 = rep(NA_real_, n),
    beta3 = rep(NA_real_, n),
    lambda = rep(NA_real_, n),
    rmse = rep(NA_real_, n),
    stringsAsFactors = FALSE
  )

  fitted_yields_matrix <- matrix(
    NA_real_, nrow = n, ncol = length(mat_cols),
    dimnames = list(as.character(dates), mat_cols)
  )

  successful_fits <- 0L
  failure_dates <- as.Date(character(0))
  failure_reasons <- character(0)
  at_bound_matrix <- matrix(
    FALSE, nrow = n, ncol = 4L,
    dimnames = list(NULL, c("beta1", "beta2", "beta3", "lambda"))
  )

  for (i in seq_along(fit_results)) {
    result <- fit_results[[i]]

    if (!isTRUE(result$success)) {
      failure_dates <- c(failure_dates, dates[i])
      msg <- result$message
      failure_reasons <- c(
        failure_reasons,
        if (is.null(msg)) "Fit did not converge" else as.character(msg)[1]
      )
      next
    }

    params_df$beta1[i] <- result$parameters[["beta1"]]
    params_df$beta2[i] <- result$parameters[["beta2"]]
    params_df$beta3[i] <- result$parameters[["beta3"]]
    params_df$lambda[i] <- result$parameters[["lambda"]]
    params_df$rmse[i] <- result$rmse

    # Place fitted values by maturity: fit_ns_model drops missing or
    # non-positive maturities, so the vector is generally shorter than the row.
    fitted_yields_matrix[i, result$valid_idx] <- result$fitted_yields

    if (!is.null(result$at_bound)) {
      at_bound_matrix[i, ] <- as.logical(result$at_bound)
    }

    successful_fits <- successful_fits + 1L
  }

  failures <- data.frame(
    Date = failure_dates, reason = failure_reasons, stringsAsFactors = FALSE
  )

  # Report fitting results
  if (nrow(failures) > 0 && nrow(failures) <= 5) {
    for (k in seq_len(nrow(failures))) {
      warning(paste("Date", failures$Date[k], ":", failures$reason[k]))
    }
  } else if (nrow(failures) > 5) {
    warning(sprintf("%d fitting errors occurred. First 5: %s",
                    nrow(failures),
                    paste(paste("Date", failures$Date[1:5], ":",
                                failures$reason[1:5]), collapse = "; ")))
  }

  if (successful_fits == 0) {
    stop("Failed to fit model to any yield curves")
  }

  # One report for the whole panel, naming each parameter that came to rest on
  # a bound and on how many dates.
  bound_counts <- colSums(at_bound_matrix, na.rm = TRUE)
  if (any(bound_counts > 0)) {
    warning(sprintf(
      paste0("Parameters rest on a bound in some fits: %s (of %d dates). ",
             "The factor decomposition is not identified on those dates; ",
             "see the at_bound element of the result, pass wider ",
             "lower_bounds/upper_bounds, or check the yield scale."),
      paste(sprintf("%s on %d", names(bound_counts)[bound_counts > 0],
                    bound_counts[bound_counts > 0]), collapse = ", "),
      length(dates)
    ))
  }

  # Success rate over every date attempted, before any rows are dropped.
  success_rate <- successful_fits / length(dates) * 100

  # Drop unfitted dates. These rows carry NA, not zero: a pre-sized numeric
  # frame would leave them as 0 and they would survive every NA filter.
  # Kept beside the parameters rather than inside them: `parameters` is a
  # numeric frame apart from its Date column, and callers coerce it.
  rownames(at_bound_matrix) <- as.character(dates)

  if (drop_failed) {
    valid_rows <- !is.na(params_df$beta1)
    params_df <- params_df[valid_rows, , drop = FALSE]
    fitted_yields_matrix <- fitted_yields_matrix[valid_rows, , drop = FALSE]
    at_bound_matrix <- at_bound_matrix[valid_rows, , drop = FALSE]
    rownames(params_df) <- NULL
  }

  # Factor changes as first differences in CALENDAR order, written back in the
  # row order of the output.
  if (nrow(params_df) > 1) {
    asc <- order(params_df$Date)
    for (f in c("beta1", "beta2", "beta3")) {
      change <- rep(NA_real_, nrow(params_df))
      change[asc] <- c(NA_real_, diff(params_df[[f]][asc]))
      params_df[[paste0(f, "_change")]] <- change
    }
  } else {
    params_df$beta1_change <- NA_real_
    params_df$beta2_change <- NA_real_
    params_df$beta3_change <- NA_real_
    warning("Only one valid yield curve - cannot calculate factor changes")
  }

  result <- list(
    parameters = params_df,
    fitted_yields = fitted_yields_matrix,
    maturities = maturities,
    success_rate = success_rate,
    failures = failures,
    at_bound = at_bound_matrix
  )
  class(result) <- c("dns_fit", "list")

  return(result)
}
