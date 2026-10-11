# Fit Dynamic Nelson-Siegel Model to Time Series of Yield Curves

Fits the DNS model to a time series of yield curves.

## Usage

``` r
fit_dns_model(
  yield_data,
  n_days = 60,
  parallel = FALSE,
  lambda_constraint = "constrained",
  initial_lambda = ns_lambda(3.5),
  show_progress = TRUE,
  date_order = c("input", "ascending", "descending"),
  drop_failed = TRUE,
  warm_start = TRUE,
  method = c("profile", "optim"),
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
  The default places the curvature hump at 3.5 years; see
  [`ns_lambda`](https://ambientpointcorp.github.io/dnsr/reference/ns_lambda.md)
  on its units.

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
  `"profile"` (the default) or `"optim"`.

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
dns_result
#> Dynamic Nelson-Siegel fit
#>   60 curves, 2024-10-03 to 2024-12-31
#>   maturities         0.0833  0.2500  0.5000  1.0000  2.0000  3.0000  5.0000  7.0000 10.0000 20.0000 30.0000
#>   rmse              4.66 bp median, 6.00 bp worst
#>   lambda            0.6169 median
#>   success rate      100.0%
summary(dns_result)
#> Dynamic Nelson-Siegel fit
#>   60 curves, 2024-10-03 to 2024-12-31, 100.0% fitted
#> 
#> Factors:
#>             min     mean      max      sd
#> beta1   0.04370  0.04662  0.04946 0.00138
#> beta2  -0.00541  0.00051  0.00663 0.00317
#> beta3  -0.03712 -0.01996 -0.00849 0.00676
#> lambda  0.43897  0.59613  0.71243 0.07070
#> 
#> Fit error (basis points):
#>   0%  50%  95% 100% 
#> 3.43 4.66 5.66 6.00 

head(dns_result$parameters)
#>         Date      beta1       beta2       beta3    lambda         rmse
#> 1 2024-10-03 0.04369712 0.006632190 -0.03712005 0.6637689 0.0005648560
#> 2 2024-10-04 0.04426872 0.005914538 -0.03051946 0.6350888 0.0005868570
#> 3 2024-10-07 0.04469211 0.005432235 -0.02950361 0.6196247 0.0006002489
#> 4 2024-10-08 0.04490083 0.004971963 -0.02983858 0.6250754 0.0005367854
#> 5 2024-10-09 0.04513430 0.004548608 -0.02885400 0.6155153 0.0004960210
#> 6 2024-10-10 0.04546945 0.004615491 -0.03037229 0.6503208 0.0005303682
#>   beta1_change  beta2_change  beta3_change
#> 1           NA            NA            NA
#> 2 0.0005716005 -7.176527e-04  0.0066005937
#> 3 0.0004233857 -4.823028e-04  0.0010158447
#> 4 0.0002087230 -4.602718e-04 -0.0003349673
#> 5 0.0002334720 -4.233556e-04  0.0009845737
#> 6 0.0003351439  6.688371e-05 -0.0015182899
```
