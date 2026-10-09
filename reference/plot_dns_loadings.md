# Plot DNS Factor Loadings

Creates a plot showing the factor loadings of the DNS model.

## Usage

``` r
plot_dns_loadings(
  maturities = seq(0, 30, length.out = 100),
  lambda = ns_lambda(3.5),
  show_labels = TRUE
)
```

## Arguments

- maturities:

  Numeric vector. Maturities (in years) for which to calculate loadings.

- lambda:

  Numeric. DNS decay parameter lambda, per unit of maturity. Maturities
  here are in years, and the default places the curvature hump at 3.5
  years, near the middle of a Treasury curve. It is written as
  `ns_lambda(3.5)` so the number states its own meaning; see
  [`ns_lambda`](https://ambientpointcorp.github.io/dnsr/reference/ns_lambda.md).
  Diebold and Li's 0.0609 is a month-denominated constant and in years
  would place that hump at 29.4 years.

- show_labels:

  Logical. If TRUE, add labels to the plot.

## Value

A ggplot object with the factor loadings plot.

## Examples

``` r
plot_dns_loadings()


# A faster decay concentrates the curvature loading at shorter maturities
plot_dns_loadings(lambda = 0.4)

```
