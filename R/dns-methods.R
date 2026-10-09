#' @title Print and Summary Methods for DNS Objects
#' @description
#' Compact displays for the objects the fitting, estimation and simulation
#' functions return.
#'
#' @details
#' These are classes layered onto the lists the package already returned: every
#' element is still present under the same name, and \code{is.list()} and
#' \code{$} access behave as before. Only the default display changes.
#'
#' @name dns-methods
#' @keywords internal
NULL

# Formats a numeric vector for one-line display.
fmt_num <- function(x, digits = 4) {
  paste(formatC(x, format = "f", digits = digits), collapse = "  ")
}


#' @rdname dns-methods
#' @param x An object returned by the corresponding dnsr function.
#' @param ... Ignored, present for consistency with the generic.
#' @return \code{x}, invisibly. Called for the display it produces.
#' @export
print.ns_fit <- function(x, ...) {
  cat("Nelson-Siegel curve fit\n")
  if (!isTRUE(x$success)) {
    cat("  fit failed",
        if (!is.null(x$message)) paste0(": ", x$message) else "", "\n", sep = "")
    return(invisible(x))
  }

  cat(sprintf("  beta1 (level)     %s\n", fmt_num(x$parameters[["beta1"]], 5)))
  cat(sprintf("  beta2 (slope)     %s\n", fmt_num(x$parameters[["beta2"]], 5)))
  cat(sprintf("  beta3 (curvature) %s\n", fmt_num(x$parameters[["beta3"]], 5)))
  cat(sprintf("  lambda            %s\n", fmt_num(x$parameters[["lambda"]], 5)))
  cat(sprintf("  rmse              %.2f bp over %d maturities\n",
              x$rmse * 10000, length(x$maturities)))

  if (!is.null(x$at_bound) && any(x$at_bound, na.rm = TRUE)) {
    cat(sprintf("  at bound          %s (decomposition not identified)\n",
                paste(names(x$at_bound)[which(x$at_bound)], collapse = ", ")))
  }
  invisible(x)
}


#' @rdname dns-methods
#' @export
print.dns_fit <- function(x, ...) {
  params <- x$parameters
  cat("Dynamic Nelson-Siegel fit\n")
  cat(sprintf("  %d curves, %s to %s\n", nrow(params),
              format(min(params$Date)), format(max(params$Date))))
  cat(sprintf("  maturities        %s\n",
              paste(format(x$maturities), collapse = " ")))
  cat(sprintf("  rmse              %.2f bp median, %.2f bp worst\n",
              stats::median(params$rmse, na.rm = TRUE) * 10000,
              max(params$rmse, na.rm = TRUE) * 10000))
  cat(sprintf("  lambda            %.4f median\n",
              stats::median(params$lambda, na.rm = TRUE)))
  cat(sprintf("  success rate      %.1f%%\n", x$success_rate))

  if (!is.null(x$failures) && nrow(x$failures) > 0) {
    cat(sprintf("  unfitted dates    %d (see $failures)\n", nrow(x$failures)))
  }
  if (!is.null(x$at_bound) && any(x$at_bound)) {
    hits <- colSums(x$at_bound)
    cat(sprintf("  at bound          %s (see $at_bound)\n",
                paste(sprintf("%s on %d", names(hits)[hits > 0], hits[hits > 0]),
                      collapse = ", ")))
  }
  invisible(x)
}


#' @rdname dns-methods
#' @export
summary.dns_fit <- function(object, ...) {
  params <- object$parameters
  factor_summary <- t(vapply(
    c("beta1", "beta2", "beta3", "lambda"),
    function(f) c(min = min(params[[f]], na.rm = TRUE),
                  mean = mean(params[[f]], na.rm = TRUE),
                  max = max(params[[f]], na.rm = TRUE),
                  sd = stats::sd(params[[f]], na.rm = TRUE)),
    numeric(4)
  ))

  out <- list(
    n_curves = nrow(params),
    period = range(params$Date),
    maturities = object$maturities,
    success_rate = object$success_rate,
    factors = factor_summary,
    rmse_bp = stats::quantile(params$rmse * 10000, c(0, 0.5, 0.95, 1),
                              na.rm = TRUE),
    n_failures = if (is.null(object$failures)) 0L else nrow(object$failures)
  )
  class(out) <- "summary.dns_fit"
  out
}


#' @rdname dns-methods
#' @export
print.summary.dns_fit <- function(x, ...) {
  cat("Dynamic Nelson-Siegel fit\n")
  cat(sprintf("  %d curves, %s to %s, %.1f%% fitted\n", x$n_curves,
              format(x$period[1]), format(x$period[2]), x$success_rate))
  cat("\nFactors:\n")
  print(round(x$factors, 5))
  cat("\nFit error (basis points):\n")
  print(round(x$rmse_bp, 2))
  if (x$n_failures > 0) {
    cat(sprintf("\n%d dates could not be fitted.\n", x$n_failures))
  }
  invisible(x)
}


#' @rdname dns-methods
#' @export
print.dns_dynamics <- function(x, ...) {
  cat("DNS factor dynamics\n")
  cat(sprintf("  estimated on     %s to %s\n",
              format(x$estimation_period[1]), format(x$estimation_period[2])))
  cat(sprintf("  sampling         %s (dt = 1/%g), %d observations\n",
              x$sampling, 1 / x$dt, nrow(x$innovations)))
  cat(sprintf("  structure        %s\n",
              if (is.null(x$structure)) "diagonal" else x$structure))

  rows <- lapply(c("beta1", "beta2", "beta3"), function(f) {
    fac <- x[[f]]
    se <- if (is.null(fac$se)) c(kappa = NA_real_) else fac$se
    c(kappa = fac$kappa, se_kappa = se[["kappa"]], theta = fac$theta,
      sigma = fac$sigma, half_life_y = fac$half_life_years)
  })
  table <- do.call(rbind, rows)
  rownames(table) <- c("level (b1)", "slope (b2)", "curvature (b3)")
  cat("\n")
  print(round(table, 4))

  cat(sprintf("\n  stable           %s (eigenvalues %s)\n",
              x$system_stability,
              paste(formatC(Re(x$eigenvalues), format = "f", digits = 3),
                    collapse = " ")))
  if (!is.null(x$structure) && identical(x$structure, "full")) {
    cat("  cross terms are in $kappa_matrix\n")
  }
  invisible(x)
}


#' @rdname dns-methods
#' @export
print.dns_simulation <- function(x, ...) {
  dims <- dim(x$factors)
  cat("DNS simulation\n")
  cat(sprintf("  %d paths, %d steps over %g years\n",
              dims[2], dims[1] - 1L, max(x$time)))
  cat(sprintf("  engine           %s%s\n",
              if (is.null(x$engine)) "yuima" else x$engine,
              if (identical(x$engine, "exact")) "" else
                paste0(" (", x$sim_method, ")")))

  terminal <- matrix(x$factors[dims[1], , ], nrow = dims[2], ncol = 3)
  cat(sprintf("  terminal factors %s\n", fmt_num(colMeans(terminal), 4)))
  cat(sprintf("    (sd)           %s\n",
              fmt_num(apply(terminal, 2, stats::sd), 4)))

  if (!is.null(x$yields)) {
    cat(sprintf("  maturities       %s\n",
                paste(format(x$maturities), collapse = " ")))
    last <- x$yields[dim(x$yields)[1], , , drop = FALSE]
    cat(sprintf("  terminal yields  %s (%% mean)\n",
                fmt_num(apply(last, 2, mean) * 100, 2)))
  }
  invisible(x)
}


#' @title Print a one-step DNS state-space fit
#' @param x A \code{dns_state_space} object.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.dns_state_space <- function(x, ...) {
  cat("DNS one-step state-space fit (Kalman filter, maximum likelihood)\n")
  cat(sprintf("  estimated on     %s to %s\n",
              format(x$estimation_period[1]), format(x$estimation_period[2])))
  cat(sprintf("  sampling         %s (dt = 1/%.1f), %d dates, %d maturities\n",
              x$sampling, 1 / x$dt, x$n_dates, length(x$maturities)))
  cat(sprintf("  lambda           %.4f%s (hump at %.2f years)\n",
              x$lambda,
              if (isTRUE(x$lambda_estimated)) " estimated" else " fixed",
              ns_curvature_peak(x$lambda)))

  rows <- lapply(c("beta1", "beta2", "beta3"), function(f) {
    fac <- x[[f]]
    c(kappa = fac$kappa, theta = fac$theta, sigma = fac$sigma,
      half_life_y = fac$half_life)
  })
  table <- do.call(rbind, rows)
  rownames(table) <- c("level (b1)", "slope (b2)", "curvature (b3)")
  cat("\n")
  print(round(table, 4))

  cat(sprintf("\n  measurement sd   %.2f bp median across maturities\n",
              stats::median(x$measurement_sd) * 10000))
  cat(sprintf("  fit              %.2f bp rmse\n", x$rmse * 10000))
  cat(sprintf("  logLik           %.2f on %d parameters (AIC %.1f, BIC %.1f)\n",
              x$logLik, x$n_params, x$aic, x$bic))
  cat(sprintf("  stable           %s\n", x$system_stability))
  if (!identical(x$convergence, 0L) && !identical(x$convergence, 0)) {
    cat(sprintf("  NOT CONVERGED    optim code %s; estimates are provisional\n",
                x$convergence))
  }
  if (!is.null(x$two_step)) {
    cat("  $two_step holds the two-step estimates used as starting values\n")
  }
  invisible(x)
}

#' @title Summarise a one-step DNS state-space fit
#'
#' @description
#' Reports the shock correlation matrix, the measurement error and residual
#' error by maturity, and the two-step estimates alongside the one-step ones.
#'
#' @param object A \code{dns_state_space} object.
#' @param ... Ignored.
#' @return A list with class \code{"summary.dns_state_space"}, invisibly
#'   printed.
#' @export
summary.dns_state_space <- function(object, ...) {
  comparison <- NULL
  if (!is.null(object$two_step)) {
    fields <- c("kappa", "theta", "sigma")
    comparison <- do.call(rbind, lapply(c("beta1", "beta2", "beta3"),
      function(f) {
        one <- vapply(fields, function(k) object[[f]][[k]], numeric(1))
        two <- vapply(fields, function(k) object$two_step[[f]][[k]], numeric(1))
        stats::setNames(c(one, two), c(paste0("one_", fields),
                                       paste0("two_", fields)))
      }))
    rownames(comparison) <- c("level (b1)", "slope (b2)", "curvature (b3)")
  }

  out <- list(
    fit = object,
    correlation = object$correlation,
    measurement_bp = object$measurement_sd * 10000,
    residual_bp_by_maturity = apply(object$residuals, 2,
                                    function(r) sqrt(mean(r^2, na.rm = TRUE))) * 10000,
    comparison = comparison
  )
  class(out) <- c("summary.dns_state_space", "list")
  out
}

#' @title Print a one-step DNS state-space summary
#' @param x A \code{summary.dns_state_space} object.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.summary.dns_state_space <- function(x, ...) {
  print(x$fit)

  cat("\nShock correlation\n")
  print(round(x$correlation, 3))

  cat("\nMeasurement error by maturity (bp)\n")
  print(round(x$measurement_bp, 2))

  cat("\nResidual rmse by maturity (bp)\n")
  print(round(x$residual_bp_by_maturity, 2))

  if (!is.null(x$comparison)) {
    cat("\nOne step against the two-step estimates it started from\n")
    print(round(x$comparison, 4))
    cat("\nThe two-step estimates treat the fitted betas as observed when\n")
    cat("they carry sampling error, which biases kappa upward. See\n")
    cat("?estimate_dns_state_space on the precision of kappa and theta.\n")
  }
  invisible(x)
}
