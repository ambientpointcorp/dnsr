# Plot Simulated DNS Factors

Creates plots of simulated DNS factors from a DNS simulation.

## Usage

``` r
plot_simulated_factors(
  sim_result,
  path_subset = NULL,
  include_quantiles = TRUE,
  plot_trajectories = FALSE
)
```

## Arguments

- sim_result:

  List. The simulation result from simulate_dns_factors().

- path_subset:

  Integer vector. Indices of paths to plot. If NULL, plots summary
  statistics.

- include_quantiles:

  Logical. If TRUE, include quantile paths.

- plot_trajectories:

  Logical. If TRUE, plot individual trajectories instead of summary
  statistics.

## Value

A list of ggplot objects, one for each factor.

## Examples

``` r
params <- list(
  kappa = c(0.02, 0.3, 0.5),
  theta = c(0.04, -0.02, 0.01),
  sigma = c(0.005, 0.01, 0.02)
)
sim_result <- simulate_dns_factors(params, c(0.035, -0.015, 0.005),
                                   n_paths = 50, engine = "exact", seed = 1)
plots <- plot_simulated_factors(sim_result)
plots$level


plot_simulated_factors(sim_result, plot_trajectories = TRUE)$level

```
