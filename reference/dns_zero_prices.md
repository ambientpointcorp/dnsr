# Zero-Coupon Bond Prices from DNS Factors

Prices zero-coupon bonds off a DNS curve.

## Usage

``` r
dns_zero_prices(factors, maturities, lambda = 0.0609, face = 100)
```

## Arguments

- factors:

  Array or matrix. DNS factors, as in
  [`dns_discount_factors`](https://ambientpointcorp.github.io/dnsr/reference/dns_discount_factors.md).

- maturities:

  Numeric vector. Maturities in years.

- lambda:

  Numeric. DNS decay parameter lambda, per unit of maturity. Maturities
  here are in years, so the 0.0609 default is the Diebold-Li value for
  maturities in months and places the curvature hump at 29.4 years
  rather than months; see
  [`ns_lambda`](https://ambientpointcorp.github.io/dnsr/reference/ns_lambda.md).

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
#>  [1] 99.647 98.942 97.891 95.808 91.721 87.748 80.161 73.072 63.390 38.755
#> [11] 23.376
```
