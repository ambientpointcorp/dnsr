# Nelson-Siegel factor loadings.
#
# The curve is y(m) = beta1 + beta2 * L2(m) + beta3 * L3(m), with
#
#   L2(m) = (1 - exp(-lambda m)) / (lambda m)
#   L3(m) = L2(m) - exp(-lambda m)
#
# Both expressions are 0/0 at lambda * m = 0 and are replaced there by their
# limits, L2 -> 1 and L3 -> 0.
#
# This is the single definition used by the fitter, the yield and forward
# calculators, the pricing helpers and the loadings plot. The arithmetic is
# written exactly as it was when each of those carried its own copy, so
# results are unchanged to the last bit.
ns_loadings <- function(maturities, lambda) {
  n <- length(maturities)
  out <- matrix(0, nrow = n, ncol = 2,
                dimnames = list(NULL, c("slope", "curvature")))

  for (i in seq_len(n)) {
    m <- maturities[i]
    if (m * lambda < 1e-10) {
      out[i, 1] <- 1.0  # limit of (1-exp(-lambda*m))/(lambda*m) as lambda*m -> 0
      out[i, 2] <- 0.0  # limit of L2 - exp(-lambda*m) as lambda*m -> 0
    } else {
      loading_2 <- (1 - exp(-lambda * m)) / (lambda * m)
      out[i, 1] <- loading_2
      out[i, 2] <- loading_2 - exp(-lambda * m)
    }
  }

  out
}

# Design matrix of the Nelson-Siegel curve at a given lambda: the three
# columns multiply beta1, beta2 and beta3 respectively. With lambda held
# fixed the model is linear in the betas, which is what the profiled fit in
# fit_ns_model(method = "profile") exploits.
ns_design <- function(maturities, lambda) {
  loadings <- ns_loadings(maturities, lambda)
  cbind(beta1 = rep(1, length(maturities)),
        beta2 = loadings[, 1],
        beta3 = loadings[, 2])
}


# Design matrix of the Nelson-Siegel-Svensson curve, which adds a second
# curvature term with its own decay:
#
#   y(m) = b1 + b2 L2(l1,m) + b3 L3(l1,m) + b4 L3(l2,m)
#
# Like the three-factor model it is linear in the betas once the decays are
# fixed, so the same profiled solve applies with a two-dimensional search.
nss_design <- function(maturities, lambda1, lambda2) {
  first <- ns_loadings(maturities, lambda1)
  second <- ns_loadings(maturities, lambda2)
  cbind(beta1 = rep(1, length(maturities)),
        beta2 = first[, 1],
        beta3 = first[, 2],
        beta4 = second[, 2])
}


# Fit the Svensson curve by concentrated least squares: solve the four betas
# exactly at each (lambda1, lambda2) and search only over the two decays.
#
# The search is parameterised as (lambda1, ratio) with lambda2 = lambda1 *
# ratio, so that both constraints the model needs are plain box bounds:
#
#   lambda1 in [0.02, 2]   ratio in [1.2, 30]
#
# Both are needed, because the profiled residual has degenerate corners that
# fit beautifully and mean nothing:
#
#   - ratio -> 1: the third and fourth loadings coincide;
#   - lambda1 -> 0: the slope loading flattens to a constant 1 and goes
#     collinear with the intercept;
#   - lambda1 large: every loading but the first decays to zero inside the
#     shortest maturity.
#
# In each corner the design goes ill-conditioned and the betas diverge to
# cancel one another. Left unbounded, a search on this panel returned
# lambda1 = 1e-4 with betas of +/-16,800 on one date and lambda1 = 40 on
# others, both with flattering residuals.
#
# The surface also has local minima, so the global step is a coarse grid --
# each evaluation is one QR of a small matrix, so a couple of hundred are cheap
# and, unlike a multi-start search, give the same answer every time -- and only
# the best cell is refined. The refinement runs Nelder-Mead on a logistic
# transform of the box rather than L-BFGS-B on the box itself: the latter's
# default finite-difference step is coarse against this surface's curvature and
# left the decays four orders of magnitude less accurate on curves generated
# exactly by the model.
fit_svensson <- function(maturities, yields,
                         lambda1_range = c(0.02, 2),
                         ratio_range = c(1.2, 30)) {
  solve_at <- function(lambda1, ratio) {
    design <- nss_design(maturities, lambda1, lambda1 * ratio)
    if (!all(is.finite(design)) || kappa(design) > 1e6) return(NULL)
    fit <- qr(design)
    if (fit$rank < ncol(design)) return(NULL)
    fit
  }

  sse <- function(par) {
    fit <- solve_at(par[1], par[2])
    if (is.null(fit)) return(.Machine$double.xmax)
    sum(qr.resid(fit, yields)^2)
  }

  # Global step: a log-spaced grid over the box.
  grid <- expand.grid(
    lambda1 = exp(seq(log(lambda1_range[1]), log(lambda1_range[2]),
                      length.out = 16)),
    ratio = exp(seq(log(ratio_range[1]), log(ratio_range[2]),
                    length.out = 12)),
    KEEP.OUT.ATTRS = FALSE
  )
  values <- vapply(seq_len(nrow(grid)),
                   function(i) sse(c(grid$lambda1[i], grid$ratio[i])),
                   numeric(1))
  if (all(values >= .Machine$double.xmax)) {
    stop("No well-conditioned Svensson design over the search range")
  }
  start <- as.numeric(grid[which.min(values), ])

  # Local step: Nelder-Mead on a logistic map of the same box, so the search is
  # unconstrained and smooth while the decays stay inside it.
  to_box <- function(u) {
    c(lambda1_range[1] + diff(lambda1_range) * stats::plogis(u[1]),
      ratio_range[1] + diff(ratio_range) * stats::plogis(u[2]))
  }
  from_box <- function(par) {
    c(stats::qlogis((par[1] - lambda1_range[1]) / diff(lambda1_range)),
      stats::qlogis((par[2] - ratio_range[1]) / diff(ratio_range)))
  }
  clamp <- function(x, lo, hi) pmin(pmax(x, lo + 1e-9), hi - 1e-9)

  opt <- tryCatch(
    stats::optim(from_box(clamp(start, c(lambda1_range[1], ratio_range[1]),
                                c(lambda1_range[2], ratio_range[2]))),
                 function(u) sse(to_box(u)), method = "Nelder-Mead",
                 control = list(reltol = 1e-13, maxit = 3000)),
    error = function(e) NULL
  )
  par <- if (!is.null(opt) && opt$value < min(values)) to_box(opt$par) else start

  fit <- solve_at(par[1], par[2])
  betas <- qr.coef(fit, yields)
  betas[is.na(betas)] <- 0
  design <- nss_design(maturities, par[1], par[1] * par[2])
  fitted <- as.numeric(design %*% betas)

  list(
    parameters = c(betas, lambda1 = par[1], lambda2 = par[1] * par[2]),
    fitted_yields = fitted,
    residuals = yields - fitted,
    rmse = sqrt(mean((yields - fitted)^2)),
    success = TRUE
  )
}


#' @title Choose the Nelson-Siegel Decay from the Curvature Peak
#'
#' @description
#' Converts between the decay parameter lambda and the maturity at which the
#' curvature loading peaks. \code{ns_lambda()} goes from a maturity to a
#' lambda; \code{ns_curvature_peak()} goes back.
#'
#' @param peak_maturity Numeric. The maturity at which the curvature loading
#'   should peak, in whatever units the maturities are measured in.
#' @param lambda Numeric. A decay parameter, per unit of maturity.
#'
#' @return A numeric vector of the same length as the input:
#'   \code{ns_lambda()} returns decays, \code{ns_curvature_peak()} maturities.
#'
#' @details
#' The curvature loading
#' \deqn{L_3(m) = \frac{1-e^{-\lambda m}}{\lambda m} - e^{-\lambda m}}
#' is maximised where \eqn{\lambda m = 1.793282}, the root of
#' \eqn{\partial L_3 / \partial m = 0}. So one defensible way to pick lambda is
#' to decide where the hump belongs and divide.
#'
#' The point of doing it this way is units. Lambda is per unit of maturity, so
#' its numerical value depends on whether maturities are years or months, and
#' nothing about a lambda on its own says which. Stating a maturity instead
#' settles it: \code{ns_lambda(30)} on a curve measured in months and
#' \code{ns_lambda(2.5)} on the same curve measured in years describe the same
#' shape.
#'
#' That also explains why this package does not use Diebold and Li's 0.0609 as
#' its default. That is the value they published for maturities in
#' \strong{months}, where it places the hump at 29.4 months. The functions here
#' take maturities in \strong{years}, so the same number places it at 29.4
#' \emph{years} -- off the end of most curves, where the slope loading is
#' nearly parallel to the level loading and the two factors are not separately
#' identified. \code{ns_curvature_peak()} reports that directly, and the
#' package default is written as \code{ns_lambda(3.5)} instead.
#'
#' Two cautions. The peak is flat -- at half or double the nominal peak
#' maturity the loading is still about 85 percent of its maximum -- so lambda
#' is not tightly determined by this choice, and the apparent precision of the
#' result is misleading. And lambda sets the slope loading's decay at the same
#' time, which is often the more consequential effect: a hump at 2.5 years
#' leaves the slope loading at 0.53 by two years, while a hump at 10 years
#' leaves it at 0.84. Where lambda can be estimated rather than chosen, for
#' instance with \code{fit_ns_model(method = "profile")}, estimating it is
#' usually better. Fixing it is still the right move for a panel, where the
#' factors have to be comparable across dates.
#'
#' @examples
#' # A hump at two and a half years, for maturities measured in years
#' ns_lambda(2.5)
#'
#' # Diebold and Li aimed for a hump at 30 months, measuring maturities in
#' # months. That criterion gives 0.0598; the value they published, 0.0609,
#' # peaks a little short of it.
#' ns_lambda(30)
#' ns_curvature_peak(0.0609)
#'
#' # What the package default implies when maturities are in years
#' ns_curvature_peak(0.0609)
#'
#' # They are inverses
#' ns_curvature_peak(ns_lambda(4))
#'
#' @seealso \code{\link{plot_dns_loadings}} to see the loadings a lambda
#'   produces, and \code{\link{fit_ns_model}} to estimate one instead.
#'
#' @export
ns_lambda <- function(peak_maturity) {
  checkmate::assert_numeric(peak_maturity, lower = 0, any.missing = FALSE,
                            min.len = 1)
  NS_CURVATURE_PEAK_PRODUCT / peak_maturity
}

#' @rdname ns_lambda
#' @export
ns_curvature_peak <- function(lambda) {
  checkmate::assert_numeric(lambda, lower = 0, any.missing = FALSE,
                            min.len = 1)
  NS_CURVATURE_PEAK_PRODUCT / lambda
}

# The product lambda * m at which the curvature loading is maximised: the root
# of d/dm [ (1-exp(-lambda m))/(lambda m) - exp(-lambda m) ] = 0, which has no
# closed form. Computed once to 6 figures rather than solved on every call.
NS_CURVATURE_PEAK_PRODUCT <- 1.793282

# Points in the log-spaced lambda scan that method = "profile" runs before it
# refines. Enough that no basin of the profiled objective is stepped over on a
# plausible interval, and cheap: each point costs one QR of an n-by-3 matrix.
PROFILE_GRID_N <- 96L
