# Simulate DNS Yield Curves using YUIMA

Simulates DNS yield curves using the YUIMA package for accurate SDE
simulation.

## Usage

``` r
simulate_dns_yield_curves(
  params,
  initial_factors,
  maturities,
  lambda = 0.0609,
  time_horizon = 10,
  n_steps = 120,
  n_paths = 100,
  include_correlation = TRUE,
  sim_method = "euler",
  seed = NULL,
  engine = c("yuima", "exact")
)
```

## Arguments

- params:

  List. Model parameters (kappa, theta, sigma).

- initial_factors:

  Numeric vector. Initial factor values (beta1, beta2, beta3).

- maturities:

  Numeric vector. Maturities (in years) for the yield curve.

- lambda:

  Numeric. DNS decay parameter lambda.

- time_horizon:

  Numeric. Simulation time horizon in years.

- n_steps:

  Integer. Number of time steps in the simulation.

- n_paths:

  Integer. Number of simulation paths.

- include_correlation:

  Logical. Whether to include correlations in simulations.

- sim_method:

  Character. Simulation method for YUIMA ("euler" or "milstein").

- seed:

  Integer. Random seed for reproducibility.

- engine:

  Character. `"yuima"` (the default) or `"exact"`; see
  [`simulate_dns_factors`](https://ambientpointcorp.github.io/dnsr/reference/simulate_dns_factors.md).

## Value

A list containing:

- yields: Array of simulated yields, indexed by (time, maturity, path)

- factors: Array of simulated factors, indexed by (time, path, factor)

- time: Vector of time points

- maturities: Vector of maturities

- params: The input parameters

- sim_method, engine: The scheme and engine used

## Examples

``` r
params <- list(
  kappa = c(0.02, 0.3, 0.5),
  theta = c(0.04, -0.02, 0.01),
  sigma = c(0.005, 0.01, 0.02)
)
initial_factors <- c(0.035, -0.015, 0.005)
maturities <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)

sim_result <- simulate_dns_yield_curves(params, initial_factors, maturities,
                                        n_paths = 50, engine = "exact",
                                        seed = 1)
sim_result
#> DNS simulation
#>   50 paths, 120 steps over 10 years
#>   engine           exact
#>   terminal factors 0.0353  -0.0181  0.0107
#>     (sd)           0.0152  0.0143  0.0193
#>   maturities        0.25  0.50  1.00  2.00  3.00  5.00  7.00 10.00 20.00 30.00
#>   terminal yields  1.73  1.76  1.80  1.88  1.96  2.10  2.22  2.39  2.78  3.01 (% mean)
```
