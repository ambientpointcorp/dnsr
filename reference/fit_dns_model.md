# Fit Dynamic Nelson-Siegel Model to Time Series of Yield Curves

Fits the DNS model to a time series of yield curves.

## Usage

``` r
fit_dns_model(
  yield_data,
  n_days = 60,
  parallel = FALSE,
  lambda_constraint = "constrained",
  initial_lambda = 0.0609,
  show_progress = TRUE,
  date_order = c("input", "ascending", "descending"),
  drop_failed = TRUE,
  warm_start = TRUE,
  method = c("optim", "profile"),
  lower_bounds = NULL,
  upper_bounds = NULL
)
```

## Arguments

- yield_data:

  Data frame. The yield curve data with a 'Date' column and remaining
  columns containing yield values for different maturities, named by
  maturity in years.

- n_days:

  Integer. Number of most recent dates to include in the fitting. The
  window is selected by date, so it does not depend on the row order of
  `yield_data`.

- parallel:

  Logical. If TRUE, use parallel processing to speed up computation.

- lambda_constraint:

  Character. How to handle the lambda parameter (fixed, constrained,
  free).

- initial_lambda:

  Numeric. Initial value for lambda, passed to
  [`fit_ns_model`](https://ambientpointcorp.github.io/dnsr/reference/fit_ns_model.md).
  See
  [`ns_lambda`](https://ambientpointcorp.github.io/dnsr/reference/ns_lambda.md)
  on its units: the 0.0609 default is a month-denominated constant and
  maturities here are in years.

- show_progress:

  Logical. If TRUE, display a progress bar.

- date_order:

  Character. Row order of the returned `parameters` and `fitted_yields`:
  `"input"` (the default) keeps the order the selected rows had in
  `yield_data`, so a most-recent-first frame yields a most-recent-first
  result; `"ascending"` and `"descending"` sort explicitly.

- drop_failed:

  Logical. If TRUE (the default) dates whose cross-section could not be
  fitted are removed from the result. If FALSE they are kept with `NA`
  parameters.

- warm_start:

  Logical. If TRUE (the default) each cross-section is initialised from
  the previous date's fitted lambda, walking forward in calendar time.
  Ignored when `parallel = TRUE`, where the fits are not ordered, and by
  `method = "profile"`, which needs no start values.

  Note that warm starting makes a date's fit depend on which *other*
  dates are in the panel: fitting a month-end series directly and
  fitting the daily series that contains it give slightly different
  parameters on the shared dates, because the chain of start values
  differs. Set `warm_start = FALSE` for cross-sections that are
  independent of the panel they arrived in.

- method:

  Character. Passed to
  [`fit_ns_model`](https://ambientpointcorp.github.io/dnsr/reference/fit_ns_model.md):
  `"optim"` (the default) or `"profile"`.

- lower_bounds, upper_bounds:

  Numeric vectors of length 4, passed to
  [`fit_ns_model`](https://ambientpointcorp.github.io/dnsr/reference/fit_ns_model.md)
  as the bounds on `(beta1, beta2, beta3, lambda)`. `NULL`, the default,
  leaves that function's own defaults in force, including its widening
  of the beta box to the data scale – which it only applies when both
  bounds are left alone.

  Supply them when the `at_bound` report says a fit is clamped. The
  default ceiling on lambda is 0.5, and real Treasury curves reach it:
  on the panel in
  [`treasury_yields`](https://ambientpointcorp.github.io/dnsr/reference/treasury_yields.md)
  a freely estimated lambda rests there on roughly a quarter of dates.
  Without these arguments there was no way to act on that from here.

## Value

A list containing:

- parameters: Data frame with model parameters for each date, plus
  `beta1_change`, `beta2_change` and `beta3_change`.

- fitted_yields: Matrix of fitted yield values, one column per maturity
  column of `yield_data`. Maturities that were missing on a given date
  are `NA`.

- maturities: Vector of maturities in years.

- success_rate: Percentage of successful model fits.

- failures: Data frame of dates that could not be fitted, with the
  reason.

- at_bound: Logical matrix, one row per retained date and one column per
  parameter, marking fits that came to rest on a bound. Those dates'
  factor decompositions are not identified.

## Details

The factor change columns are first differences in calendar order
(\\\beta(t) - \beta(t-1)\\) regardless of the row order of the output,
so they can be fed directly to
[`estimate_dns_dynamics`](https://ambientpointcorp.github.io/dnsr/reference/estimate_dns_dynamics.md).

## Examples

``` r
data(treasury_yields)

# The 60 most recent curves. The window is chosen by date, so passing the
# panel in either row order selects the same dates.
dns_result <- fit_dns_model(treasury_yields, n_days = 60, show_progress = FALSE)
#> Warning: Parameters rest on a bound in some fits: lambda on 38 (of 60 dates). The factor decomposition is not identified on those dates; see the at_bound element of the result, pass wider lower_bounds/upper_bounds, or check the yield scale.
dns_result
#> Dynamic Nelson-Siegel fit
#>   60 curves, 2024-10-03 to 2024-12-31
#>   maturities         0.0833  0.2500  0.5000  1.0000  2.0000  3.0000  5.0000  7.0000 10.0000 20.0000 30.0000
#>   rmse              5.10 bp median, 8.04 bp worst
#>   lambda            0.5000 median
#>   success rate      100.0%
#>   at bound          lambda on 38 (see $at_bound)
summary(dns_result)
#> Dynamic Nelson-Siegel fit
#>   60 curves, 2024-10-03 to 2024-12-31, 100.0% fitted
#> 
#> Factors:
#>             min     mean      max      sd
#> beta1   0.04489  0.04717  0.05004 0.00129
#> beta2  -0.00638 -0.00041  0.00411 0.00274
#> beta3  -0.03911 -0.02026 -0.00849 0.00742
#> lambda  0.49157  0.49818  0.50000 0.00330
#> 
#> Fit error (basis points):
#>   0%  50%  95% 100% 
#> 3.54 5.10 6.89 8.04 

head(dns_result$parameters)
#>         Date      beta1       beta2       beta3    lambda         rmse
#> 1 2024-10-03 0.04521919 0.004054147 -0.03910520 0.5000000 0.0008035182
#> 2 2024-10-04 0.04535422 0.004113613 -0.03217322 0.4999392 0.0007049269
#> 3 2024-10-07 0.04564922 0.003862919 -0.03097348 0.5000000 0.0006886792
#> 4 2024-10-08 0.04589329 0.003333767 -0.03120564 0.5000000 0.0006464879
#> 5 2024-10-09 0.04603303 0.003073846 -0.03007115 0.5000000 0.0005936956
#> 6 2024-10-10 0.04663458 0.002660000 -0.03174549 0.4999100 0.0006832128
#>   beta1_change  beta2_change  beta3_change
#> 1           NA            NA            NA
#> 2 0.0001350378  5.946613e-05  0.0069319757
#> 3 0.0002949941 -2.506935e-04  0.0011997457
#> 4 0.0002440761 -5.291528e-04 -0.0002321662
#> 5 0.0001397336 -2.599206e-04  0.0011344893
#> 6 0.0006015521 -4.138464e-04 -0.0016743405
```
