# Discount Factors, Zero Prices and Par Yields from DNS Factors

Converts DNS factors into the quantities a cashflow valuation actually
consumes: discount factors, zero-coupon bond prices and par coupon
yields.

Converts DNS factors into discount factors.

## Usage

``` r
dns_discount_factors(factors, maturities, lambda = 0.0609)
```

## Arguments

- factors:

  Array or matrix. If an array, simulated factors indexed by (time,
  path, factor); if a matrix, factors for a single path indexed by
  (time, factor). The same shapes
  [`calculate_dns_yields`](https://ambientpointcorp.github.io/dnsr/reference/calculate_dns_yields.md)
  accepts.

- maturities:

  Numeric vector. Maturities in years.

- lambda:

  Numeric. DNS decay parameter lambda, per unit of maturity. Maturities
  here are in years, so the 0.0609 default is the Diebold-Li value for
  maturities in months and places the curvature hump at 29.4 years
  rather than months; see
  [`ns_lambda`](https://ambientpointcorp.github.io/dnsr/reference/ns_lambda.md).

## Value

An array of discount factors indexed by (time, maturity, path), or a
matrix indexed by (time, maturity) when `factors` is a matrix – the same
shape
[`calculate_dns_yields`](https://ambientpointcorp.github.io/dnsr/reference/calculate_dns_yields.md)
returns.

## Details

The DNS curve gives continuously compounded zero rates, so the discount
factor at maturity \\m\\ is \$\$P(m) = \exp(-y(m) \\ m).\$\$ A zero
maturity discounts by 1.

## See also

[`dns_zero_prices`](https://ambientpointcorp.github.io/dnsr/reference/dns_zero_prices.md),
[`dns_par_yields`](https://ambientpointcorp.github.io/dnsr/reference/dns_par_yields.md)

## Examples

``` r
data(treasury_yields)
maturities <- as.numeric(setdiff(names(treasury_yields), "Date"))

# Fit the most recent curve, then discount off it
curve <- fit_ns_model(maturities, as.numeric(treasury_yields[nrow(treasury_yields), -1]))
factors <- matrix(curve$parameters[1:3], nrow = 1)

discount <- dns_discount_factors(factors, maturities, curve$parameters[["lambda"]])
round(discount[1, ], 4)
#>  [1] 0.9965 0.9894 0.9789 0.9581 0.9172 0.8775 0.8016 0.7307 0.6339 0.3876
#> [11] 0.2338
```
