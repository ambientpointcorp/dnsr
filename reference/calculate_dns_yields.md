# Calculate Yield Curves from DNS Factors

Calculates yield curves based on DNS factors for specified maturities.

## Usage

``` r
calculate_dns_yields(factors, maturities, lambda = ns_lambda(3.5))
```

## Arguments

- factors:

  Array or matrix. If an array, simulated factors indexed by (time,
  path, factor); if a matrix, factors for a single path indexed by
  (time, factor).

- maturities:

  Numeric vector. Maturities (in years) for the yield curve.

- lambda:

  Numeric. DNS decay parameter lambda, per unit of maturity. Maturities
  here are in years, and the default places the curvature hump at 3.5
  years, near the middle of a Treasury curve. It is written as
  `ns_lambda(3.5)` so the number states its own meaning; see
  [`ns_lambda`](https://ambientpointcorp.github.io/dnsr/reference/ns_lambda.md).
  Diebold and Li's 0.0609 is a month-denominated constant and in years
  would place that hump at 29.4 years.

## Value

Array of yields indexed by (time, maturity, path) if `factors` is an
array, or a matrix indexed by (time, maturity) if it is a matrix.

## Examples

``` r
maturities <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)

# A single curve: one row of (level, slope, curvature)
factors <- matrix(c(0.035, -0.015, 0.005), nrow = 1)
round(calculate_dns_yields(factors, maturities, ns_lambda(3.5))[1, ], 5)
#>  [1] 0.02122 0.02231 0.02418 0.02695 0.02882 0.03101 0.03215 0.03303 0.03402
#> [10] 0.03435

# Three paths over four dates, as returned by simulate_dns_factors()
paths <- array(rep(c(0.035, -0.015, 0.005), each = 12), dim = c(4, 3, 3))
dim(calculate_dns_yields(paths, maturities, ns_lambda(3.5)))
#> [1]  4 10  3
```
