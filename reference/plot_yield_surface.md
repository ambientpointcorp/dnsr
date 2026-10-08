# Plot 3D Yield Surface

Creates a 3D visualization of the yield curve surface.

## Usage

``` r
plot_yield_surface(
  yield_data,
  interactive = TRUE,
  sample_dates = 50,
  color_scale = "Viridis"
)
```

## Arguments

- yield_data:

  Data frame. Historical yield curve data with 'Date' column.

- interactive:

  Logical. If TRUE, create an interactive plot with plotly.

- sample_dates:

  Integer. Number of dates to sample for plotting.

- color_scale:

  Character. Color scale for the surface.

## Value

A 3D plot object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Needs plotly (interactive) or lattice (static)
data(treasury_yields)
plot_yield_surface(treasury_yields)
} # }
```
