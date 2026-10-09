# Calculate Forward Rates from DNS Factors

Calculates instantaneous forward rates based on DNS factors.

## Usage

``` r
calculate_dns_forwards(factors, maturities, lambda = 0.0609)
```

## Arguments

- factors:

  Matrix or array. DNS factors.

- maturities:

  Numeric vector. Maturities (in years) for forward rates.

- lambda:

  Numeric. DNS decay parameter lambda, per unit of maturity. Maturities
  here are in years, so the 0.0609 default is the Diebold-Li value for
  maturities in months and places the curvature hump at 29.4 years
  rather than months; see
  [`ns_lambda`](https://ambientpointcorp.github.io/dnsr/reference/ns_lambda.md).

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
  spot = calculate_dns_yields(factors, maturities, 0.0609)[1, ],
  forward = calculate_dns_forwards(factors, maturities, 0.0609)[1, ]
)
#>               [,1]       [,2]       [,3]       [,4]       [,5]       [,6]
#> spot    0.02015129 0.02030067 0.02059382 0.02115838 0.02169523 0.02269159
#> forward 0.02030162 0.02059755 0.02117275 0.02225928 0.02326569 0.02506045
#>               [,7]       [,8]       [,9]      [,10]
#> spot    0.02359372 0.02479111 0.02773946 0.02960272
#> forward 0.02659789 0.02849774 0.03236423 0.03405635
```
