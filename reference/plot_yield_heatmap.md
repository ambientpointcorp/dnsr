# Create Heat Map of Yield Curve Evolution

Creates a heat map visualization of yield curve evolution over time.

## Usage

``` r
plot_yield_heatmap(
  yield_data,
  reverse_time = FALSE,
  plot_type = "yields",
  smooth_data = FALSE,
  color_palette = NULL
)
```

## Arguments

- yield_data:

  Data frame. Historical yield curve data with 'Date' column.

- reverse_time:

  Logical. Flips the date axis. A date axis increases upward, so by
  default the oldest date is at the bottom and the most recent at the
  top; `TRUE` puts the most recent at the bottom. (Before 0.4.0 this
  argument did nothing, and its documentation described the default
  orientation rather than the effect of setting it.)

- plot_type:

  Character. One of "yields", "changes", or "curvature".

- smooth_data:

  Logical. If TRUE, apply smoothing to the heat map data.

- color_palette:

  Character vector. Colors for the heat map.

## Value

A ggplot object with the heat map.

## Examples

``` r
data(treasury_yields)
plot_yield_heatmap(treasury_yields)


# Most recent dates at the top
plot_yield_heatmap(treasury_yields, reverse_time = TRUE)

```
