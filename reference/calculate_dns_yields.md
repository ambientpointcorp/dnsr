# Calculate Yield Curves from DNS Factors

Calculates yield curves based on DNS factors for specified maturities.

## Usage

``` r
calculate_dns_yields(factors, maturities, lambda = 0.0609)
```

## Arguments

- factors:

  Array or matrix. If an array, simulated factors indexed by (time,
  path, factor); if a matrix, factors for a single path indexed by
  (time, factor).

- maturities:

  Numeric vector. Maturities (in years) for the yield curve.

- lambda:

  Numeric. DNS decay parameter lambda.

## Value

Array of yields indexed by (time, maturity, path) if `factors` is an
array, or a matrix indexed by (time, maturity) if it is a matrix.

## Examples

``` r
maturities <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)

# A single curve: one row of (level, slope, curvature)
factors <- matrix(c(0.035, -0.015, 0.005), nrow = 1)
round(calculate_dns_yields(factors, maturities, 0.0609)[1, ], 5)
#>  [1] 0.02015 0.02030 0.02059 0.02116 0.02170 0.02269 0.02359 0.02479 0.02774
#> [10] 0.02960

# Three paths over four dates, as returned by simulate_dns_factors()
paths <- array(rep(c(0.035, -0.015, 0.005), each = 12), dim = c(4, 3, 3))
dim(calculate_dns_yields(paths, maturities, 0.0609))
#> [1]  4 10  3
```
