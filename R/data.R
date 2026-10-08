#' @title U.S. Treasury Yield Curves, 2019-2024
#'
#' @description
#' Daily U.S. Treasury constant maturity yields across eleven maturities, in
#' the shape the rest of the package expects: a \code{Date} column followed by
#' one column per maturity, named by maturity in years.
#'
#' Bundled so that the examples, vignettes and tests have a real yield panel to
#' work on without a FRED API key.
#'
#' @format A data frame with 1,501 rows and 12 columns:
#' \describe{
#'   \item{Date}{Observation date, ascending.}
#'   \item{0.0833, 0.25, 0.5}{1-month, 3-month and 6-month yields.}
#'   \item{1, 2, 3, 5, 7, 10, 20, 30}{1- to 30-year yields.}
#' }
#' Yields are in decimal form, so 0.0425 is 4.25 percent. There are no missing
#' values: market holidays, which FRED reports as rows with no observation at
#' any maturity, are excluded rather than carried as \code{NA}.
#'
#' @source
#' Federal Reserve Economic Data (FRED), series \code{DGS1MO}, \code{DGS3MO},
#' \code{DGS6MO}, \code{DGS1}, \code{DGS2}, \code{DGS3}, \code{DGS5},
#' \code{DGS7}, \code{DGS10}, \code{DGS20} and \code{DGS30}. These are the
#' Treasury constant maturity series of the Federal Reserve's H.15 release,
#' retrieved 2026-10-08 and covering 2 January 2019 to 31 December 2024. As
#' U.S. government works they are in the public domain.
#'
#' The script that builds this dataset is \code{data-raw/treasury_yields.R} in
#' the package sources; it calls \code{\link{get_treasury_yield_data}}.
#'
#' @details
#' Six years of daily observations spans 72 month-ends, which is enough for
#' \code{\link{estimate_dns_dynamics}} to run in its default monthly sampling
#' mode. The window covers the 2020 collapse to near-zero short rates, the
#' 2022-2023 inversion, and the start of the subsequent normalisation, so it
#' contains level, slope and curvature variation rather than a single regime.
#'
#' @examples
#' data(treasury_yields)
#'
#' # The most recent curve
#' tail(treasury_yields, 1)
#'
#' # Dates on which the 3-month yield exceeded the 10-year
#' inverted <- treasury_yields$`0.25` > treasury_yields$`10`
#' range(treasury_yields$Date[inverted])
"treasury_yields"
