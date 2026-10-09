# Calculate Forward Rates from DNS Factors

Calculates instantaneous forward rates based on DNS factors.

## Usage

``` r
calculate_dns_forwards(factors, maturities, lambda = ns_lambda(3.5))
```

## Arguments

- factors:

  Matrix or array. DNS factors.

- maturities:

  Numeric vector. Maturities (in years) for forward rates.

- lambda:

  Numeric. DNS decay parameter lambda, per unit of maturity. Maturities
  here are in years, and the default places the curvature hump at 3.5
  years, near the middle of a Treasury curve. It is written as
  `ns_lambda(3.5)` so the number states its own meaning; see
  [`ns_lambda`](https://ambientpointcorp.github.io/dnsr/reference/ns_lambda.md).
  Diebold and Li's 0.0609 is a month-denominated constant and in years
  would place that hump at 29.4 years.

## Value

Matrix or array of forward rates with the same structure as the yield
output from calculate_dns_yields.

## Details

The instantaneous forward rate in the DNS model is given by: \$\$f(m) =
\beta_1 + \beta_2 e^{-\lambda m} + \beta_3 \lambda m e^{-\lambda m}\$\$

## Examples

``` r
maturities <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)
factors <- matrix(c(0.035, -0.015, 0.005), nrow = 1)

# Forwards sit above spots where the curve is rising
rbind(
  spot = calculate_dns_yields(factors, maturities, ns_lambda(3.5))[1, ],
  forward = calculate_dns_forwards(factors, maturities, ns_lambda(3.5))[1, ]
)
#>               [,1]       [,2]       [,3]       [,4]       [,5]       [,6]
#> spot    0.02121510 0.02230819 0.02417968 0.02694920 0.02881798 0.03101193
#> forward 0.02236687 0.02438143 0.02754859 0.03145546 0.03342735 0.03483095
#>               [,7]       [,8]       [,9]      [,10]
#> spot    0.03215056 0.03303012 0.03402399 0.03434942
#> forward 0.03508122 0.03506322 0.03500128 0.03500001
```
