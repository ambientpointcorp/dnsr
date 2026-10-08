# Dynamic Nelson-Siegel Model Simulation using YUIMA

Functions to simulate yield curves using the Dynamic Nelson-Siegel model
with YUIMA package for SDE simulation.

Simulates the evolution of DNS factors using a multivariate
Ornstein-Uhlenbeck process implemented with the YUIMA package for
accurate stochastic differential equation simulation.

## Usage

``` r
simulate_dns_factors(
  params,
  initial_factors,
  time_horizon = 10,
  n_steps = 120,
  n_paths = 100,
  include_correlation = TRUE,
  seed = NULL,
  sim_method = "euler",
  engine = c("yuima", "exact")
)
```

## Arguments

- params:

  List. Model parameters (kappa, theta, sigma).

- initial_factors:

  Numeric vector. Initial factor values (beta1, beta2, beta3).

- time_horizon:

  Numeric. Simulation time horizon in years.

- n_steps:

  Integer. Number of time steps in the simulation.

- n_paths:

  Integer. Number of simulation paths.

- include_correlation:

  Logical. Whether to include correlations in simulations.

- seed:

  Integer. Random seed for reproducibility.

- sim_method:

  Character. Discretisation scheme passed to YUIMA, either `"euler"`
  (the default) or `"milstein"`. Ignored by `engine = "exact"`, which
  does not discretise.

- engine:

  Character. `"yuima"` (the default) simulates through the YUIMA
  package, one call per path. `"exact"` samples the exact Gaussian
  transition law of the factor system directly: it needs no extra
  package, carries no discretisation error at any step size, and is
  roughly two orders of magnitude faster. The default is expected to
  change to `"exact"` in a future release.

## Value

A list containing:

- factors: Array of simulated factors (time, path, factor)

- time: Vector of time points

- params: The input parameters

- sim_method: The scheme used

- engine: The engine used

## Details

The simulation is based on the multivariate Ornstein-Uhlenbeck process:
\$\$d\beta(t) = \kappa(\theta - \beta(t))dt + \sigma dW(t)\$\$ where
\\\kappa\\ is the mean-reversion matrix, \\\theta\\ is the long-run mean
vector, \\\sigma\\ is the volatility matrix, and \\W(t)\\ is a
correlated Wiener process.

This implementation uses the YUIMA package, which provides efficient and
accurate methods for simulating stochastic differential equations.

The diffusion coefficient of this model does not depend on the state, so
the Milstein correction term vanishes and `sim_method = "milstein"`
coincides with the Euler-Maruyama scheme. The argument is honoured for
completeness; accuracy is governed by `n_steps`.

That accuracy ceiling is what `engine = "exact"` removes. The factor
system is a linear SDE with state-independent diffusion, so its
transition density over any step is Gaussian with a mean and covariance
available in closed form. Sampling that law directly is exact in
distribution at every step size, and replaces one YUIMA call per path
with vectorised draws over all paths at once. The two engines agree in
distribution; they do not produce the same paths from the same seed.

## Examples

``` r
params <- list(
  kappa = c(0.02, 0.3, 0.5),
  theta = c(0.04, -0.02, 0.01),
  sigma = c(0.005, 0.01, 0.02),
  correlation = diag(3)
)
initial_factors <- c(0.035, -0.015, 0.005)

# The exact engine needs no extra package
sim_result <- simulate_dns_factors(params, initial_factors, n_paths = 50,
                                   engine = "exact", seed = 1)
sim_result
#> DNS simulation
#>   50 paths, 120 steps over 10 years
#>   engine           exact
#>   terminal factors 0.0353  -0.0181  0.0107
#>     (sd)           0.0152  0.0143  0.0193

if (FALSE) { # \dontrun{
# The YUIMA engine, for comparison
simulate_dns_factors(params, initial_factors, n_paths = 50, engine = "yuima")
} # }
```
