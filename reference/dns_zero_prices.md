# Zero-Coupon Bond Prices from DNS Factors

Prices zero-coupon bonds off a DNS curve.

## Usage

``` r
dns_zero_prices(factors, maturities, lambda = ns_lambda(3.5), face = 100)
```

## Arguments

- factors:

  Array or matrix. DNS factors, as in
  [`dns_discount_factors`](https://ambientpointcorp.github.io/dnsr/reference/dns_discount_factors.md).

- maturities:

  Numeric vector. Maturities in years.

- lambda:

  Numeric. DNS decay parameter lambda, per unit of maturity. Maturities
  here are in years, and the default places the curvature hump at 3.5
  years, near the middle of a Treasury curve. It is written as
  `ns_lambda(3.5)` so the number states its own meaning; see
  [`ns_lambda`](https://ambientpointcorp.github.io/dnsr/reference/ns_lambda.md).
  Diebold and Li's 0.0609 is a month-denominated constant and in years
  would place that hump at 29.4 years.

- face:

  Numeric. Face value of the bond. Defaults to 100.

## Value

Prices in the same shape as
[`dns_discount_factors`](https://ambientpointcorp.github.io/dnsr/reference/dns_discount_factors.md).

## Details

A zero-coupon bond is worth its face value discounted to today, so this
is `face` times the discount factor. It is separated out because quoting
per 100 of face is the market convention.

## See also

[`dns_discount_factors`](https://ambientpointcorp.github.io/dnsr/reference/dns_discount_factors.md),
[`dns_par_yields`](https://ambientpointcorp.github.io/dnsr/reference/dns_par_yields.md)

## Examples

``` r
data(treasury_yields)
maturities <- as.numeric(setdiff(names(treasury_yields), "Date"))
curve <- fit_ns_model(maturities, as.numeric(treasury_yields[nrow(treasury_yields), -1]))
factors <- matrix(curve$parameters[1:3], nrow = 1)

round(dns_zero_prices(factors, maturities, curve$parameters[["lambda"]])[1, ], 3)
#>  [1] 99.635 98.921 97.882 95.872 91.951 88.043 80.307 72.955 62.997 38.497
#> [11] 23.519
```
