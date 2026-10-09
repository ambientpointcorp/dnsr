#' @title Discount Factors, Zero Prices and Par Yields from DNS Factors
#' @description
#' Converts DNS factors into the quantities a cashflow valuation actually
#' consumes: discount factors, zero-coupon bond prices and par coupon yields.

#' @title Discount Factors from DNS Factors
#'
#' @description
#' Converts DNS factors into discount factors.
#'
#' @param factors Array or matrix. If an array, simulated factors indexed by
#'   (time, path, factor); if a matrix, factors for a single path indexed by
#'   (time, factor). The same shapes \code{\link{calculate_dns_yields}} accepts.
#' @param maturities Numeric vector. Maturities in years.
#' @param lambda Numeric. DNS decay parameter lambda, per unit of maturity.
#'   Maturities here are in years, and the default places the curvature hump at
#'   3.5 years, near the middle of a Treasury curve. It is written as
#'   \code{ns_lambda(3.5)} so the number states its own meaning; see
#'   \code{\link{ns_lambda}}. Diebold and Li's 0.0609 is a month-denominated
#'   constant and in years would place that hump at 29.4 years.
#'
#' @return An array of discount factors indexed by (time, maturity, path), or a
#'   matrix indexed by (time, maturity) when \code{factors} is a matrix --
#'   the same shape \code{\link{calculate_dns_yields}} returns.
#'
#' @details
#' The DNS curve gives continuously compounded zero rates, so the discount
#' factor at maturity \eqn{m} is
#' \deqn{P(m) = \exp(-y(m) \, m).}
#' A zero maturity discounts by 1.
#'
#' @examples
#' data(treasury_yields)
#' maturities <- as.numeric(setdiff(names(treasury_yields), "Date"))
#'
#' # Fit the most recent curve, then discount off it
#' curve <- fit_ns_model(maturities, as.numeric(treasury_yields[nrow(treasury_yields), -1]))
#' factors <- matrix(curve$parameters[1:3], nrow = 1)
#'
#' discount <- dns_discount_factors(factors, maturities, curve$parameters[["lambda"]])
#' round(discount[1, ], 4)
#'
#' @seealso \code{\link{dns_zero_prices}}, \code{\link{dns_par_yields}}
#'
#' @export
dns_discount_factors <- function(factors, maturities,
                                 lambda = ns_lambda(3.5)) {
  checkmate::assert_numeric(maturities, lower = 0, any.missing = FALSE,
                            min.len = 1)
  checkmate::assert_numeric(lambda, lower = 0, len = 1, any.missing = FALSE)

  yields <- calculate_dns_yields(factors, maturities, lambda)

  # Maturity varies along the second dimension, so recycling a vector that
  # repeats each maturity n_times over lines up with the column-major layout
  # of both the matrix and the array case.
  n_times <- dim(yields)[1]
  exp(-yields * rep(maturities, each = n_times))
}


#' @title Zero-Coupon Bond Prices from DNS Factors
#'
#' @description
#' Prices zero-coupon bonds off a DNS curve.
#'
#' @param factors Array or matrix. DNS factors, as in
#'   \code{\link{dns_discount_factors}}.
#' @param maturities Numeric vector. Maturities in years.
#' @param lambda Numeric. DNS decay parameter lambda, per unit of maturity.
#'   Maturities here are in years, and the default places the curvature hump at
#'   3.5 years, near the middle of a Treasury curve. It is written as
#'   \code{ns_lambda(3.5)} so the number states its own meaning; see
#'   \code{\link{ns_lambda}}. Diebold and Li's 0.0609 is a month-denominated
#'   constant and in years would place that hump at 29.4 years.
#' @param face Numeric. Face value of the bond. Defaults to 100.
#'
#' @return Prices in the same shape as \code{\link{dns_discount_factors}}.
#'
#' @details
#' A zero-coupon bond is worth its face value discounted to today, so this is
#' \code{face} times the discount factor. It is separated out because quoting
#' per 100 of face is the market convention.
#'
#' @examples
#' data(treasury_yields)
#' maturities <- as.numeric(setdiff(names(treasury_yields), "Date"))
#' curve <- fit_ns_model(maturities, as.numeric(treasury_yields[nrow(treasury_yields), -1]))
#' factors <- matrix(curve$parameters[1:3], nrow = 1)
#'
#' round(dns_zero_prices(factors, maturities, curve$parameters[["lambda"]])[1, ], 3)
#'
#' @seealso \code{\link{dns_discount_factors}}, \code{\link{dns_par_yields}}
#'
#' @export
dns_zero_prices <- function(factors, maturities, lambda = ns_lambda(3.5),
                            face = 100) {
  checkmate::assert_numeric(face, len = 1, any.missing = FALSE)
  face * dns_discount_factors(factors, maturities, lambda)
}


#' @title Par Coupon Yields from DNS Factors
#'
#' @description
#' Computes the par coupon yield at each maturity: the coupon rate that prices
#' a bullet bond at par off the DNS curve.
#'
#' @param factors Array or matrix. DNS factors, as in
#'   \code{\link{dns_discount_factors}}.
#' @param maturities Numeric vector. Maturities in years. Each must be a whole
#'   number of coupon periods; see Details.
#' @param lambda Numeric. DNS decay parameter lambda, per unit of maturity.
#'   Maturities here are in years, and the default places the curvature hump at
#'   3.5 years, near the middle of a Treasury curve. It is written as
#'   \code{ns_lambda(3.5)} so the number states its own meaning; see
#'   \code{\link{ns_lambda}}. Diebold and Li's 0.0609 is a month-denominated
#'   constant and in years would place that hump at 29.4 years.
#' @param frequency Integer. Coupon payments per year. Defaults to 2, the
#'   U.S. Treasury convention.
#'
#' @return Par yields in the same shape as
#'   \code{\link{dns_discount_factors}}, as annualised decimal rates.
#'
#' @details
#' For a bond maturing at \eqn{T} with \eqn{n = T f} coupon dates
#' \eqn{t_i = i / f}, the coupon rate that makes its price equal its face value
#' is
#' \deqn{c = f \, \frac{1 - P(T)}{\sum_{i=1}^{n} P(t_i)}.}
#' The denominator is the annuity factor: the present value of one unit paid at
#' every coupon date.
#'
#' Par yields are a coupon-bond quantity, so they are only defined on
#' maturities that are a whole number of coupon periods. At the default
#' semiannual frequency a 6-month maturity qualifies but a 1-month or 3-month
#' one does not; those return \code{NA} with a warning rather than being forced
#' onto a convention they do not have. Bill maturities belong in
#' \code{\link{calculate_dns_yields}} or \code{\link{dns_discount_factors}}.
#'
#' The gap between the par yield and the zero yield at the same maturity has
#' two parts, which can pull in opposite directions:
#'
#' \itemize{
#'   \item A compounding convention term. The DNS curve is continuously
#'     compounded while the par yield is a rate payable \code{frequency} times
#'     a year, which adds roughly \eqn{y^2 / (2f)} and is always positive.
#'   \item A curve shape term. Coupons are discounted at the rates earlier in
#'     the curve, not at the terminal rate alone, so an upward-sloping zero
#'     curve pulls the par yield down and an inverted one pushes it up.
#' }
#'
#' The convention term dominates at short maturities, where the shape term is
#' small; the shape term dominates further out. On a steeply upward-sloping
#' curve the par yield is therefore slightly above the zero yield at one year
#' and clearly below it at thirty.
#'
#' @examples
#' data(treasury_yields)
#' maturities <- as.numeric(setdiff(names(treasury_yields), "Date"))
#' curve <- fit_ns_model(maturities, as.numeric(treasury_yields[nrow(treasury_yields), -1]))
#' factors <- matrix(curve$parameters[1:3], nrow = 1)
#' lambda <- curve$parameters[["lambda"]]
#'
#' coupon_maturities <- c(0.5, 1, 2, 5, 10, 30)
#' zero <- calculate_dns_yields(factors, coupon_maturities, lambda)
#' par <- dns_par_yields(factors, coupon_maturities, lambda)
#'
#' data.frame(
#'   maturity = coupon_maturities,
#'   zero_pct = round(zero[1, ] * 100, 3),
#'   par_pct = round(par[1, ] * 100, 3)
#' )
#'
#' @seealso \code{\link{dns_discount_factors}}, \code{\link{dns_zero_prices}}
#'
#' @export
dns_par_yields <- function(factors, maturities, lambda = ns_lambda(3.5),
                           frequency = 2) {
  checkmate::assert_numeric(maturities, lower = 0, any.missing = FALSE,
                            min.len = 1)
  checkmate::assert_numeric(lambda, lower = 0, len = 1, any.missing = FALSE)
  checkmate::assert_count(frequency, positive = TRUE)

  period <- 1 / frequency
  n_coupons <- maturities * frequency
  whole <- abs(n_coupons - round(n_coupons)) < 1e-8 & n_coupons >= 1
  n_coupons <- round(n_coupons)

  if (any(!whole)) {
    warning(sprintf(
      "Par yields are undefined for %s: not a whole number of %g-month coupon periods. Returning NA.",
      paste(maturities[!whole], collapse = ", "), 12 / frequency
    ))
  }

  # Discount factors on every coupon date needed by any requested maturity,
  # computed once.
  if (!any(whole)) {
    out <- calculate_dns_yields(factors, maturities, lambda)
    out[] <- NA_real_
    return(out)
  }

  coupon_grid <- seq_len(max(n_coupons[whole])) * period
  discount <- dns_discount_factors(factors, coupon_grid, lambda)

  is_array <- length(dim(discount)) == 3L
  n_times <- dim(discount)[1]

  # Cumulative annuity factor along the maturity dimension: column i is the
  # present value of one unit paid at each of the first i coupon dates.
  annuity <- if (is_array) {
    aperm(apply(discount, c(1, 3), cumsum), c(2, 1, 3))
  } else {
    t(apply(discount, 1, cumsum))
  }

  # Assemble the output on the requested maturity grid.
  template <- calculate_dns_yields(factors, maturities, lambda)
  template[] <- NA_real_

  for (k in which(whole)) {
    idx <- n_coupons[k]
    if (is_array) {
      template[, k, ] <- frequency * (1 - discount[, idx, ]) / annuity[, idx, ]
    } else {
      template[, k] <- frequency * (1 - discount[, idx]) / annuity[, idx]
    }
  }

  template
}
