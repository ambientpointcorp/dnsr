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
