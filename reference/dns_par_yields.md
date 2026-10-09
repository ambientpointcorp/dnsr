# Par Coupon Yields from DNS Factors

Computes the par coupon yield at each maturity: the coupon rate that
prices a bullet bond at par off the DNS curve.

## Usage

``` r
dns_par_yields(factors, maturities, lambda = 0.0609, frequency = 2)
```

## Arguments

- factors:

  Array or matrix. DNS factors, as in
  [`dns_discount_factors`](https://ambientpointcorp.github.io/dnsr/reference/dns_discount_factors.md).

- maturities:

  Numeric vector. Maturities in years. Each must be a whole number of
  coupon periods; see Details.

- lambda:

  Numeric. DNS decay parameter lambda, per unit of maturity. Maturities
  here are in years, so the 0.0609 default is the Diebold-Li value for
  maturities in months and places the curvature hump at 29.4 years
  rather than months; see
  [`ns_lambda`](https://ambientpointcorp.github.io/dnsr/reference/ns_lambda.md).

- frequency:

  Integer. Coupon payments per year. Defaults to 2, the U.S. Treasury
  convention.

## Value

Par yields in the same shape as
[`dns_discount_factors`](https://ambientpointcorp.github.io/dnsr/reference/dns_discount_factors.md),
as annualised decimal rates.

## Details

For a bond maturing at \\T\\ with \\n = T f\\ coupon dates \\t_i = i /
f\\, the coupon rate that makes its price equal its face value is \$\$c
= f \\ \frac{1 - P(T)}{\sum\_{i=1}^{n} P(t_i)}.\$\$ The denominator is
the annuity factor: the present value of one unit paid at every coupon
date.

Par yields are a coupon-bond quantity, so they are only defined on
maturities that are a whole number of coupon periods. At the default
semiannual frequency a 6-month maturity qualifies but a 1-month or
3-month one does not; those return `NA` with a warning rather than being
forced onto a convention they do not have. Bill maturities belong in
[`calculate_dns_yields`](https://ambientpointcorp.github.io/dnsr/reference/calculate_dns_yields.md)
or
[`dns_discount_factors`](https://ambientpointcorp.github.io/dnsr/reference/dns_discount_factors.md).

The gap between the par yield and the zero yield at the same maturity
has two parts, which can pull in opposite directions:

- A compounding convention term. The DNS curve is continuously
  compounded while the par yield is a rate payable `frequency` times a
  year, which adds roughly \\y^2 / (2f)\\ and is always positive.

- A curve shape term. Coupons are discounted at the rates earlier in the
  curve, not at the terminal rate alone, so an upward-sloping zero curve
  pulls the par yield down and an inverted one pushes it up.

The convention term dominates at short maturities, where the shape term
is small; the shape term dominates further out. On a steeply
upward-sloping curve the par yield is therefore slightly above the zero
yield at one year and clearly below it at thirty.

## See also

[`dns_discount_factors`](https://ambientpointcorp.github.io/dnsr/reference/dns_discount_factors.md),
[`dns_zero_prices`](https://ambientpointcorp.github.io/dnsr/reference/dns_zero_prices.md)

## Examples

``` r
data(treasury_yields)
maturities <- as.numeric(setdiff(names(treasury_yields), "Date"))
curve <- fit_ns_model(maturities, as.numeric(treasury_yields[nrow(treasury_yields), -1]))
factors <- matrix(curve$parameters[1:3], nrow = 1)
lambda <- curve$parameters[["lambda"]]

coupon_maturities <- c(0.5, 1, 2, 5, 10, 30)
zero <- calculate_dns_yields(factors, coupon_maturities, lambda)
par <- dns_par_yields(factors, coupon_maturities, lambda)

data.frame(
  maturity = coupon_maturities,
  zero_pct = round(zero[1, ] * 100, 3),
  par_pct = round(par[1, ] * 100, 3)
)
#>   maturity zero_pct par_pct
#> 1      0.5    4.263   4.309
#> 2      1.0    4.283   4.329
#> 3      2.0    4.321   4.367
#> 4      5.0    4.423   4.466
#> 5     10.0    4.559   4.590
#> 6     30.0    4.845   4.808
```
