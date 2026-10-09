#' @importFrom stats optim na.omit
NULL

# ---- parameter transforms ---------------------------------------------------
#
# The optimiser works in unconstrained coordinates. These maps impose the
# constraints the state space requires at every trial point: kappa and sigma
# positive, the correlation matrix positive definite, the measurement
# variances positive. The alternative, rejecting illegal trial points with a
# penalty, introduces a discontinuity at the boundary.

# Three unconstrained reals to a 3-by-3 correlation matrix, through the
# canonical partial correlations of its Cholesky factor. Every real triple
# maps to a positive definite matrix with a unit diagonal, and every such
# matrix is reachable.
corr_from_unconstrained <- function(y) {
  z <- tanh(y)
  n <- 3L
  L <- matrix(0, n, n)
  L[1, 1] <- 1
  k <- 1L
  for (i in 2:n) {
    L[i, 1] <- z[k]; k <- k + 1L
    acc <- sqrt(max(0, 1 - L[i, 1]^2))
    if (i > 2) {
      for (j in 2:(i - 1)) {
        L[i, j] <- z[k] * acc; k <- k + 1L
        acc <- acc * sqrt(max(0, 1 - z[k - 1L]^2))
      }
    }
    L[i, i] <- acc
  }
  R <- L %*% t(L)
  # Guard the diagonal against accumulated rounding; the filter divides by it.
  diag(R) <- 1
  R
}

# The inverse, used to start the optimiser from the two-step correlation.
unconstrained_from_corr <- function(R) {
  L <- t(chol(R))
  y <- numeric(3)
  y[1] <- L[2, 1]
  y[2] <- L[3, 1]
  acc3 <- sqrt(max(0, 1 - L[3, 1]^2))
  y[3] <- if (acc3 > 0) L[3, 2] / acc3 else 0
  atanh(pmin(pmax(y, -0.999), 0.999))
}


# ---- the state space --------------------------------------------------------
#
# Measurement:  y_t = Lambda(lambda) f_t + e_t,     e_t ~ N(0, H), H diagonal
# Transition:   f_t = theta + A (f_{t-1} - theta) + n_t,  n_t ~ N(0, Q)
#
# A and Q are not free parameters. They are the exact transition law of the
# same Ornstein-Uhlenbeck system the rest of the package simulates:
#
#   A_ii = exp(-kappa_i dt)
#   Q_ij = rho_ij sigma_i sigma_j (1 - exp(-(kappa_i + kappa_j) dt)) /
#                                 (kappa_i + kappa_j)
#
# Deriving them rather than estimating them freely has three consequences.
# Q is a covariance matrix at every trial point. kappa > 0 makes the system
# stationary, so the filter has a proper initial covariance. And the result is
# in the same continuous-time parameters dns_dynamics_to_sim_params() and
# simulate_dns_factors() take, so a fit can be simulated forward without a
# conversion step. The cost is that cross-factor terms cannot be represented.
ou_transition <- function(kappa, sigma, corr, dt) {
  kappa_sum <- outer(kappa, kappa, "+")
  scale <- ifelse(abs(kappa_sum) > 1e-10,
                  (1 - exp(-kappa_sum * dt)) / kappa_sum,
                  dt)
  list(A = diag(exp(-kappa * dt), 3),
       Q = corr * outer(sigma, sigma) * scale,
       # Stationary covariance of the same system, the limit of Q as dt
       # grows. Used as the prior for f_0, which assumes the process starts
       # in its stationary distribution.
       P0 = corr * outer(sigma, sigma) / kappa_sum)
}

# Kalman filter log-likelihood, with an optional smoothing pass.
#
# Missing yields shrink the measurement equation for that date rather than
# dropping the date, so a curve missing one maturity still contributes the
# others and leaves no gap in the series. Dates with no observations
# contribute no measurement update and propagate the state.
dns_kalman <- function(Y, loadings, kappa, theta, sigma, corr, meas_var, dt,
                       smooth = FALSE) {
  tr <- ou_transition(kappa, sigma, corr, dt)
  A <- tr$A; Q <- tr$Q

  n_t <- nrow(Y)
  x <- theta
  P <- tr$P0
  ll <- 0

  # Kept for the smoother.
  xp <- matrix(NA_real_, n_t, 3); Pp <- vector("list", n_t)
  xf <- matrix(NA_real_, n_t, 3); Pf <- vector("list", n_t)

  for (t in seq_len(n_t)) {
    # Predict
    x_pred <- as.numeric(theta + A %*% (x - theta))
    P_pred <- A %*% P %*% t(A) + Q
    P_pred <- (P_pred + t(P_pred)) / 2
    xp[t, ] <- x_pred; Pp[[t]] <- P_pred

    obs <- which(!is.na(Y[t, ]))
    if (length(obs) == 0L) {
      x <- x_pred; P <- P_pred
      xf[t, ] <- x; Pf[[t]] <- P
      next
    }

    L_t <- loadings[obs, , drop = FALSE]
    v <- Y[t, obs] - as.numeric(L_t %*% x_pred)
    F_t <- L_t %*% P_pred %*% t(L_t) + diag(meas_var[obs], length(obs))
    F_t <- (F_t + t(F_t)) / 2

    # The Cholesky factor gives the log determinant and the solve from one
    # decomposition, and errors if F_t is not positive definite.
    ch <- tryCatch(chol(F_t), error = function(e) NULL)
    if (is.null(ch)) return(list(logLik = -Inf))

    half <- backsolve(ch, v, transpose = TRUE)
    ll <- ll - 0.5 * (length(obs) * log(2 * pi) +
                      2 * sum(log(diag(ch))) + sum(half^2))

    # Update
    PL <- P_pred %*% t(L_t)
    K <- t(backsolve(ch, backsolve(ch, t(PL), transpose = TRUE)))
    x <- as.numeric(x_pred + K %*% v)
    P <- P_pred - K %*% L_t %*% P_pred
    P <- (P + t(P)) / 2
    xf[t, ] <- x; Pf[[t]] <- P
  }

  out <- list(logLik = ll, filtered = xf)
  if (!smooth) return(out)

  # Rauch-Tung-Striebel backward pass. The smoothed factors use the whole
  # sample at every date; the two-step betas use only that date's
  # cross-section.
  xs <- xf; Ps <- Pf
  for (t in (n_t - 1):1) {
    Pn <- Pp[[t + 1]]
    J <- tryCatch(Pf[[t]] %*% t(A) %*% solve(Pn), error = function(e) NULL)
    if (is.null(J)) break
    xs[t, ] <- xf[t, ] + as.numeric(J %*% (xs[t + 1, ] - xp[t + 1, ]))
    Ps[[t]] <- Pf[[t]] + J %*% (Ps[[t + 1]] - Pn) %*% t(J)
  }
  out$smoothed <- xs
  out$smoothed_var <- Ps
  out
}


# Two-step estimates of the factor system, at a given dt.
#
# Computed here rather than through estimate_dns_dynamics() for two reasons:
# that function re-subsamples the panel and selects its own dt, which would
# put the starting values on a different time scale from the filter's, and it
# requires 20 observations. Reported on the result as $two_step.
#
# The mapping back to continuous time agrees with ou_transition():
# Phi_ii = exp(-kappa_i dt), and the AR(1) residual variance is
# sigma_i^2 (1 - exp(-2 kappa_i dt)) / (2 kappa_i).
two_step_estimates <- function(betas, dt) {
  n <- nrow(betas)
  if (n < 4L) return(NULL)
  lagged <- betas[-n, , drop = FALSE]
  current <- betas[-1L, , drop = FALSE]

  kappa <- theta <- sigma <- numeric(3)
  resid <- matrix(NA_real_, n - 1L, 3)
  for (i in 1:3) {
    fit <- stats::lm(current[, i] ~ lagged[, i])
    phi <- stats::coef(fit)[[2]]
    intercept <- stats::coef(fit)[[1]]
    # A sample can show no mean reversion, or explosive dynamics. Clamped
    # into the stationary interior; these are starting values only.
    phi <- min(max(phi, 1e-4), 0.9999)
    kappa[i] <- -log(phi) / dt
    theta[i] <- intercept / (1 - phi)
    resid[, i] <- stats::resid(fit)
    step_var <- stats::var(resid[, i])
    scale <- (1 - exp(-2 * kappa[i] * dt)) / (2 * kappa[i])
    sigma[i] <- sqrt(max(step_var, .Machine$double.eps) / max(scale, 1e-12))
  }

  corr <- stats::cor(resid)
  corr[!is.finite(corr)] <- 0
  diag(corr) <- 1
  factor_names <- c("beta1", "beta2", "beta3")
  list(
    beta1 = list(kappa = kappa[1], theta = theta[1], sigma = sigma[1],
                 half_life = log(2) / kappa[1]),
    beta2 = list(kappa = kappa[2], theta = theta[2], sigma = sigma[2],
                 half_life = log(2) / kappa[2]),
    beta3 = list(kappa = kappa[3], theta = theta[3], sigma = sigma[3],
                 half_life = log(2) / kappa[3]),
    correlation = `dimnames<-`(corr, list(factor_names, factor_names)),
    dt = dt
  )
}


# ---- the estimator ----------------------------------------------------------

#' @title Estimate the Dynamic Nelson-Siegel Model in One Step
#'
#' @description
#' Fits the DNS model to a panel of yield curves by maximum likelihood, with
#' the factors treated as an unobserved state and filtered out of the yields,
#' rather than fitted curve by curve and then modelled. This is the one-step
#' state-space approach of Diebold, Rudebusch and Aruoba (2006), as against
#' the two-step approach of Diebold and Li (2006) that
#' \code{\link{fit_dns_model}} and \code{\link{estimate_dns_dynamics}}
#' implement between them.
#'
#' @details
#' The model is
#' \deqn{y_t = \Lambda(\lambda) f_t + e_t, \quad e_t \sim N(0, H)}
#' \deqn{f_t = \theta + A (f_{t-1} - \theta) + \eta_t, \quad \eta_t \sim N(0, Q)}
#' with \eqn{\Lambda} the Nelson-Siegel loadings and \eqn{H} diagonal.
#' \eqn{A} and \eqn{Q} are the transition law of the Ornstein-Uhlenbeck system
#' the rest of the package simulates:
#' \deqn{A_{ii} = e^{-\kappa_i \Delta t}}
#' \deqn{Q_{ij} = \rho_{ij} \sigma_i \sigma_j (1 - e^{-(\kappa_i+\kappa_j)\Delta t}) / (\kappa_i + \kappa_j)}
#'
#' \eqn{A} and \eqn{Q} are therefore functions of \eqn{\kappa}, \eqn{\sigma}
#' and \eqn{\rho} rather than free parameters. This has three consequences.
#' \eqn{Q} is a covariance matrix at every trial point. \eqn{\kappa > 0} makes
#' the system stationary, so the filter has a proper initial covariance. And
#' the result is in the parameters \code{\link{dns_dynamics_to_sim_params}} and
#' \code{\link{simulate_dns_factors}} accept, so a fit can be simulated
#' forward without conversion. The cost is that cross-factor terms cannot be
#' represented.
#'
#' \strong{One step against two.} The two-step estimator uses the fitted betas
#' as data. They are estimates carrying sampling error, and using them as
#' observed in a time-series regression biases mean reversion upward, the
#' errors-in-variables direction. Filtering the factors propagates that
#' uncertainty into the likelihood and estimates \eqn{\lambda} jointly with
#' the dynamics. The costs are run time, needing the whole panel at once, and
#' a harder numerical problem: a dozen or more parameters in one optimisation
#' rather than a sequence of independent three-parameter fits.
#'
#' \strong{Limitations.} The factor system is diagonal: each factor reverts to
#' its own mean at its own speed, with contemporaneously correlated shocks. A
#' full mean-reversion matrix is not supported. Supporting it would require a
#' matrix exponential for \eqn{A} and a Sylvester solve for \eqn{Q}, and
#' \code{\link{dns_dynamics_to_sim_params}} discards cross terms, so the
#' result could not be simulated by the rest of the package. Use
#' \code{estimate_dns_dynamics(structure = "full")} on two-step betas for the
#' cross terms.
#'
#' Missing yields shrink the measurement equation for that date rather than
#' dropping the date, so a curve missing one maturity still contributes its
#' other maturities and leaves no gap in the series.
#'
#' \strong{Precision of the estimates.} The parameters divide into two groups.
#' Simulating from a known DNS state space with monthly data and ten
#' maturities, then re-estimating, gives these percentage errors:
#'
#' \itemize{
#'   \item \eqn{\lambda}: within 1 percent at 20, 40 and 100 years of data.
#'   \item \eqn{\sigma}: within 5 percent at all three lengths.
#'   \item measurement standard deviation: within 1 percent at all three.
#'   \item \eqn{\kappa} for the level factor: 167 percent too high at 20
#'     years, 77 percent at 40, 58 percent at 100.
#'   \item \eqn{\theta} for the slope factor: 54 percent low at 20 years,
#'     66 percent high at 40, 2 percent low at 100.
#' }
#'
#' \eqn{\lambda}, the shock volatilities and the measurement error are
#' recovered precisely at realistic sample lengths, and the smoothed factors
#' correlate above 0.99 with the simulated ones. The mean reversion speeds and
#' long-run means are not. \eqn{\kappa} carries an upward bias that decays as
#' the sample lengthens, and \eqn{\theta} has a standard error of the same
#' order as its own value, because a factor with a multi-year half-life
#' supplies few independent observations about its long-run mean however
#' finely it is sampled.
#'
#' The errors shrink toward zero as the sample grows, so this reflects the
#' difficulty of estimating mean reversion rather than an error in the
#' implementation, and it applies to the two-step route as well. Treat
#' \eqn{\kappa} and \eqn{\theta} as approximate.
#'
#' @param yield_data Data frame. Yield curve panel with a \code{Date} column
#'   and one column per maturity, named by maturity in years -- the same shape
#'   \code{\link{fit_dns_model}} takes.
#' @param n_days Integer or NULL. Use only the most recent \code{n_days}
#'   dates. NULL (the default) uses all of them.
#' @param sampling Character. \code{"monthly"} (the default) keeps the last
#'   observation in each calendar month, the frequency this model is usually
#'   estimated at; \code{"all"} uses every date.
#' @param lambda Numeric or NULL. NULL (the default) estimates lambda jointly
#'   with the dynamics. A number holds it fixed, which is faster.
#' @param measurement Character. \code{"per_maturity"} (the default) gives each
#'   maturity its own measurement error variance, as in Diebold, Rudebusch and
#'   Aruoba; \code{"common"} estimates a single variance for all maturities,
#'   which is better determined on short panels.
#' @param dt Numeric or NULL. Time step between observations, in years. NULL
#'   (the default) infers it from the median spacing of the dates used.
#' @param start List or NULL. Optional starting values, with any of the
#'   elements \code{kappa}, \code{theta}, \code{sigma}, \code{correlation},
#'   \code{lambda}, \code{measurement_sd}. NULL (the default) starts from the
#'   two-step estimates, which are also reported on the result.
#' @param maxit Integer. Maximum optimiser iterations. Default 500.
#' @param verbose Logical. If TRUE, report the log-likelihood as it improves.
#'
#' @return A list with class \code{"dns_state_space"}. Its
#'   \code{beta1}/\code{beta2}/\code{beta3}, \code{correlation} and
#'   \code{lambda} elements are shaped so that
#'   \code{\link{dns_dynamics_to_sim_params}} accepts it directly:
#'   \itemize{
#'     \item beta1, beta2, beta3: each a list of \code{kappa}, \code{theta}
#'       and \code{sigma}, plus the implied \code{half_life}
#'     \item correlation: the 3-by-3 shock correlation matrix
#'     \item lambda: the decay parameter, estimated or as supplied
#'     \item factors: data frame of smoothed factors by date, directly
#'       comparable with the \code{parameters} frame \code{fit_dns_model()}
#'       returns
#'     \item filtered: the filtered factors, using only data up to each date
#'     \item measurement_sd: fitted measurement error per maturity, in the
#'       units of the input yields
#'     \item logLik, n_params, aic, bic: fit and its cost in parameters
#'     \item fitted, residuals: fitted yields from the smoothed factors, and
#'       observed minus fitted
#'     \item rmse: root mean squared measurement residual
#'     \item kappa_matrix, theta_vector, eigenvalues, system_stability: as
#'       \code{\link{estimate_dns_dynamics}} reports them
#'     \item two_step: the two-step estimates used as starting values
#'     \item convergence, message: from \code{\link[stats]{optim}}
#'     \item dt, sampling, structure, n_dates, maturities
#'   }
#'
#' @references
#' Diebold, F.X., Rudebusch, G.D. and Aruoba, S.B. (2006). The macroeconomy and
#' the yield curve: a dynamic latent factor approach. \emph{Journal of
#' Econometrics} 131, 309-338.
#'
#' Diebold, F.X. and Li, C. (2006). Forecasting the term structure of
#' government bond yields. \emph{Journal of Econometrics} 130, 337-364.
#'
#' @seealso \code{\link{fit_dns_model}} and
#'   \code{\link{estimate_dns_dynamics}} for the two-step route,
#'   \code{\link{dns_dynamics_to_sim_params}} to simulate a fit forward.
#'
#' @examples
#' # Maximising the likelihood takes a few seconds, so the example is not run
#' # during checks. The run time is dominated by the number of optimiser
#' # iterations rather than by the size of the panel.
#' \donttest{
#' # Two years of the bundled panel, thinned to month-ends. lambda is held
#' # fixed here; pass lambda = NULL to estimate it jointly with the dynamics.
#' # n_days counts rows of the input, which is daily, so 500 trading days
#' # leaves 24 month-ends to filter.
#' fit <- estimate_dns_state_space(treasury_yields, n_days = 500,
#'                                 lambda = ns_lambda(3.5),
#'                                 measurement = "common")
#' fit
#'
#' # The shock correlation, the measurement error by maturity, and the
#' # two-step estimates alongside the one-step ones
#' summary(fit)
#'
#' # The smoothed factors have the same shape as the two-step betas
#' head(fit$factors)
#'
#' # The fit is accepted by the simulation parameter converter
#' params <- dns_dynamics_to_sim_params(fit)
#' str(params)
#' }
#'
#' @export
estimate_dns_state_space <- function(yield_data, n_days = NULL,
                                     sampling = c("monthly", "all"),
                                     lambda = NULL,
                                     measurement = c("per_maturity", "common"),
                                     dt = NULL, start = NULL,
                                     maxit = 500L, verbose = FALSE) {
  sampling <- match.arg(sampling)
  measurement <- match.arg(measurement)

  checkmate::assert_data_frame(yield_data, min.rows = 3)
  if (!"Date" %in% names(yield_data)) {
    stop("yield_data must have a 'Date' column")
  }
  if (!is.null(lambda)) {
    checkmate::assert_numeric(lambda, lower = 1e-4, upper = 50, len = 1,
                              any.missing = FALSE)
  }
  checkmate::assert_count(maxit, positive = TRUE)
  checkmate::assert_logical(verbose, len = 1)
  if (!is.null(start)) checkmate::assert_list(start, names = "unique")

  panel <- yield_data
  panel$Date <- as.Date(panel$Date)
  panel <- panel[!is.na(panel$Date), , drop = FALSE]
  # The filter runs forward in time regardless of the input row order.
  panel <- panel[order(panel$Date), , drop = FALSE]

  if (!is.null(n_days)) {
    checkmate::assert_count(n_days, positive = TRUE)
    keep <- utils::tail(seq_len(nrow(panel)), n_days)
    panel <- panel[keep, , drop = FALSE]
  }

  n_supplied <- nrow(panel)
  if (sampling == "monthly") {
    month <- format(panel$Date, "%Y-%m")
    panel <- panel[!duplicated(month, fromLast = TRUE), , drop = FALSE]
  }

  maturity_names <- setdiff(names(panel), "Date")
  maturities <- suppressWarnings(as.numeric(maturity_names))
  if (anyNA(maturities) || any(maturities <= 0)) {
    stop("Every column other than 'Date' must be named by a positive ",
         "maturity in years; got ",
         paste(maturity_names[is.na(maturities) | maturities <= 0],
               collapse = ", "))
  }

  Y <- as.matrix(panel[, maturity_names, drop = FALSE])
  storage.mode(Y) <- "double"
  usable <- rowSums(!is.na(Y)) >= 3L
  if (sum(usable) < 5L) {
    stop("Only ", sum(usable), " date(s) have 3 or more observed yields, ",
         "and at least 5 are needed to identify a state-space model. ",
         if (sampling == "monthly")
           paste0("sampling = \"monthly\" kept ", nrow(panel),
                  " of the ", n_supplied, " dates supplied, one per calendar",
                  " month; pass sampling = \"all\" or a larger n_days.")
         else "")
  }
  Y <- Y[usable, , drop = FALSE]
  dates <- panel$Date[usable]

  if (is.null(dt)) {
    gaps <- as.numeric(diff(dates))
    dt <- stats::median(gaps[gaps > 0]) / 365.25
  }
  checkmate::assert_numeric(dt, lower = 1e-8, len = 1, any.missing = FALSE)

  # ---- starting values, from the two-step estimates ----
  #
  # The two-step fit is cheap and lands near the optimum, and is reported on
  # the result so the one-step estimates can be compared against it.
  lambda_start <- if (!is.null(lambda)) lambda else
    if (!is.null(start$lambda)) start$lambda else ns_lambda(3.5)

  two_step <- NULL
  ts_fit <- try(suppressWarnings(fit_dns_model(
    panel[usable, , drop = FALSE], n_days = nrow(Y), show_progress = FALSE,
    lambda_constraint = "fixed", initial_lambda = lambda_start
  )), silent = TRUE)
  if (!inherits(ts_fit, "try-error") && nrow(ts_fit$parameters) >= 4L) {
    # Sorted by date before lagging. fit_dns_model() returns rows in the
    # order it received them, and an AR(1) fitted to a reversed series
    # estimates the time-reverse of the process.
    ts_params <- ts_fit$parameters[order(ts_fit$parameters$Date), ,
                                   drop = FALSE]
    betas <- as.matrix(ts_params[, c("beta1", "beta2", "beta3")])
    two_step <- two_step_estimates(betas, dt)
  }

  pick <- function(user, from_two_step, fallback) {
    if (!is.null(user)) return(user)
    if (!is.null(from_two_step) && all(is.finite(from_two_step))) {
      return(from_two_step)
    }
    fallback
  }
  ts_get <- function(field) {
    if (is.null(two_step)) return(NULL)
    vapply(c("beta1", "beta2", "beta3"),
           function(b) two_step[[b]][[field]], numeric(1))
  }

  kappa0 <- pmax(pick(start$kappa, ts_get("kappa"), c(0.3, 0.5, 1.0)), 1e-3)
  theta0 <- pick(start$theta, ts_get("theta"),
                 c(mean(Y, na.rm = TRUE), 0, 0))
  sigma0 <- pmax(pick(start$sigma, ts_get("sigma"), rep(0.01, 3)), 1e-6)
  corr0  <- pick(start$correlation,
                 if (!is.null(two_step)) two_step$correlation else NULL,
                 diag(3))
  # A sample correlation can sit on the boundary; moved inside so the
  # transform to unconstrained coordinates is finite.
  corr0 <- as.matrix(corr0)
  if (min(eigen(corr0, only.values = TRUE)$values) < 1e-6) {
    corr0 <- 0.95 * corr0 + 0.05 * diag(3)
  }

  n_mat <- length(maturities)
  meas0 <- if (!is.null(start$measurement_sd)) {
    rep_len(start$measurement_sd, if (measurement == "common") 1L else n_mat)
  } else if (!inherits(ts_fit, "try-error")) {
    r <- stats::sd(as.numeric(ts_fit$fitted_yields) - as.numeric(Y),
                   na.rm = TRUE)
    rep(max(r, 1e-6), if (measurement == "common") 1L else n_mat)
  } else {
    rep(5e-4, if (measurement == "common") 1L else n_mat)
  }

  # ---- pack, unpack, objective ----
  estimate_lambda <- is.null(lambda)
  n_meas <- length(meas0)

  unpack <- function(p) {
    i <- 0L
    # The packed vector picks up names from the two-step starting values.
    # Dropped here, so the returned scalars carry none.
    nxt <- function(k) { out <- unname(p[i + seq_len(k)]); i <<- i + k; out }
    theta <- nxt(3L)
    kappa <- exp(nxt(3L))
    sigma <- exp(nxt(3L))
    corr  <- corr_from_unconstrained(nxt(3L))
    msd   <- exp(nxt(n_meas))
    lam   <- if (estimate_lambda) exp(nxt(1L)) else lambda
    list(theta = theta, kappa = kappa, sigma = sigma, corr = corr,
         meas_sd = msd, lambda = lam)
  }

  p0 <- unname(c(theta0, log(kappa0), log(sigma0),
                 unconstrained_from_corr(corr0), log(meas0)))
  if (estimate_lambda) p0 <- c(p0, log(unname(lambda_start)))

  best <- -Inf
  neg_ll <- function(p) {
    if (any(!is.finite(p))) return(1e10)
    th <- unpack(p)
    # Outside this range the loadings are numerically flat and the
    # likelihood is uninformative about lambda.
    if (th$lambda < 1e-4 || th$lambda > 50) return(1e10)
    loadings <- ns_design(maturities, th$lambda)
    mv <- rep_len(th$meas_sd, n_mat)^2
    out <- dns_kalman(Y, loadings, th$kappa, th$theta, th$sigma, th$corr,
                      mv, dt)
    if (!is.finite(out$logLik)) return(1e10)
    if (verbose && out$logLik > best) {
      best <<- out$logLik
      message(sprintf("logLik %.3f", out$logLik))
    }
    -out$logLik
  }

  opt <- stats::optim(p0, neg_ll, method = "BFGS",
                      control = list(maxit = maxit))

  th <- unpack(opt$par)
  loadings <- ns_design(maturities, th$lambda)
  mv <- rep_len(th$meas_sd, n_mat)^2
  filt <- dns_kalman(Y, loadings, th$kappa, th$theta, th$sigma, th$corr,
                     mv, dt, smooth = TRUE)

  states <- if (!is.null(filt$smoothed)) filt$smoothed else filt$filtered
  fitted <- states %*% t(loadings)
  dimnames(fitted) <- dimnames(Y)
  residuals <- Y - fitted

  n_params <- length(opt$par)
  n_obs <- sum(!is.na(Y))
  half_life <- log(2) / th$kappa

  factor_frame <- data.frame(Date = dates, beta1 = states[, 1],
                             beta2 = states[, 2], beta3 = states[, 3])
  filtered_frame <- data.frame(Date = dates, beta1 = filt$filtered[, 1],
                               beta2 = filt$filtered[, 2],
                               beta3 = filt$filtered[, 3])

  factor_names <- c("beta1", "beta2", "beta3")
  result <- list(
    beta1 = list(kappa = th$kappa[1], theta = th$theta[1],
                 sigma = th$sigma[1], half_life = half_life[1]),
    beta2 = list(kappa = th$kappa[2], theta = th$theta[2],
                 sigma = th$sigma[2], half_life = half_life[2]),
    beta3 = list(kappa = th$kappa[3], theta = th$theta[3],
                 sigma = th$sigma[3], half_life = half_life[3]),
    correlation = `dimnames<-`(th$corr, list(factor_names, factor_names)),
    lambda = th$lambda,
    lambda_estimated = estimate_lambda,
    factors = factor_frame,
    filtered = filtered_frame,
    measurement_sd = stats::setNames(rep_len(th$meas_sd, n_mat),
                                     maturity_names),
    logLik = filt$logLik,
    n_params = n_params,
    n_obs = n_obs,
    aic = -2 * filt$logLik + 2 * n_params,
    bic = -2 * filt$logLik + log(n_obs) * n_params,
    fitted = fitted,
    residuals = residuals,
    rmse = sqrt(mean(residuals^2, na.rm = TRUE)),
    kappa_matrix = `dimnames<-`(diag(th$kappa, 3),
                                list(factor_names, factor_names)),
    theta_vector = stats::setNames(th$theta, factor_names),
    eigenvalues = th$kappa,
    system_stability = all(th$kappa > 0),
    two_step = two_step,
    estimation_period = range(dates),
    n_dates = nrow(Y),
    maturities = maturities,
    dt = dt,
    sampling = sampling,
    structure = "diagonal",
    convergence = opt$convergence,
    message = if (is.null(opt$message)) NA_character_ else opt$message
  )

  if (opt$convergence != 0) {
    warning("The state-space likelihood did not converge (optim code ",
            opt$convergence, "). Treat the estimates as provisional: try a ",
            "larger maxit, measurement = \"common\", or holding lambda ",
            "fixed.")
  }

  class(result) <- c("dns_state_space", "list")
  result
}
