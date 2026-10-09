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
plots <- plot_factor_time_series(dns_result)
plots$factors

plots$rmse

```
