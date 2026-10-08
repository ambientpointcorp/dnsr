# Plot DNS Factor Loadings

Creates a plot showing the factor loadings of the DNS model.

## Usage

``` r
plot_dns_loadings(
  maturities = seq(0, 30, length.out = 100),
  lambda = 0.0609,
  show_labels = TRUE
)
```

## Arguments

- maturities:

  Numeric vector. Maturities (in years) for which to calculate loadings.

- lambda:

  Numeric. DNS decay parameter lambda.

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
