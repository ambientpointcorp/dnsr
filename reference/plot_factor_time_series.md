# Plot DNS Factor Time Series

Creates plots of DNS factor time series and fit quality.

## Usage

``` r
plot_factor_time_series(dns_result, show_changes = FALSE)
```

## Arguments

- dns_result:

  List. The DNS fitting result from fit_dns_model().

- show_changes:

  Logical. If TRUE, include plots of factor changes.

## Value

A list of ggplot objects:

- factors: Plot of factor time series

- rmse: Plot of RMSE time series

- changes: Plot of factor changes (if show_changes = TRUE)

## Examples

``` r
data(treasury_yields)
dns_result <- fit_dns_model(treasury_yields, n_days = 250, show_progress = FALSE)
#> Warning: Parameters rest on a bound in some fits: lambda on 107 (of 250 dates). The factor decomposition is not identified on those dates; see the at_bound element of the result, widen lower_bounds/upper_bounds, or check the yield scale.
plots <- plot_factor_time_series(dns_result)
plots$factors

plots$rmse

```
