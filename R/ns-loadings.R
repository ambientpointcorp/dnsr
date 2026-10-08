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
