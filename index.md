# dnsr: Dynamic Nelson-Siegel Model in R

Reference documentation: <https://ambientpointcorp.github.io/dnsr>

## Overview

The `dnsr` package provides a comprehensive implementation of the
Dynamic Nelson-Siegel (DNS) model for yield curve modeling and analysis.
It offers a complete workflow from data acquisition to model fitting,
factor dynamics estimation, advanced stochastic simulation, and
visualization.

## Features

- **Data Acquisition**: Download and process U.S. Treasury yield curve
  data, with retries and an optional on-disk cache
- **Model Calibration**: Robust fitting of the DNS model to yield curve
  time series, numerically or by concentrated least squares
- **Factor Dynamics**: Estimate the stochastic processes governing yield
  curve evolution, as independent factors or as a full VAR(1), with
  standard errors
- **Advanced Simulation**: Generate future yield curve scenarios, either
  through YUIMA or by sampling the exact transition law
- **Pricing**: Discount factors, zero-coupon prices and par yields off a
  fitted curve
- **Visualization**: Create informative plots of yield curves and model
  results
- **Parallel Processing**: Support for faster model fitting with
  multiple cores
- **Bundled data**: A six-year Treasury panel, so everything above can
  be tried without an API key

## Installation

`dnsr` is not on CRAN; install it from GitHub:

``` r

# install.packages("remotes")
remotes::install_github("ambientpointcorp/dnsr")

# For the SDE simulation paths, also install YUIMA (a Suggests dependency)
install.packages("yuima")
```

## Setup

Before using the package, you need a free FRED API key to download
Treasury data:

1.  **Get a FRED API key**: Visit <https://fred.stlouisfed.org/> and
    create a free account
2.  **Request an API key**: Go to
    <https://fred.stlouisfed.org/docs/api/api_key.html>
3.  **Set your key**:

``` r

# Option 1: Set for current session
Sys.setenv(fred_api_key = "your_api_key_here")

# Option 2: Add to .Renviron for permanent setup
usethis::edit_r_environ()
# Add line: fred_api_key=your_api_key_here
```

## Quick Start

The package ships a six-year Treasury panel, so this runs as written:

``` r

library(dnsr)
data(treasury_yields)

# Fit the 750 most recent curves -- about three years, which is 36 month-ends,
# enough for the dynamics estimation below. "profile" solves the betas in
# closed form at each lambda, so it needs no start values.
dns_result <- fit_dns_model(treasury_yields, n_days = 750, method = "profile")
dns_result

# Factor dynamics, month-end sampled by default
dynamics <- estimate_dns_dynamics(dns_result$parameters)
dynamics

# Simulate forward from the most recent curve
latest <- dns_result$parameters[which.max(dns_result$parameters$Date), ]
params <- dns_dynamics_to_sim_params(dynamics, lambda = latest$lambda)

sim <- simulate_dns_yield_curves(
  params = params,
  initial_factors = as.numeric(latest[c("beta1", "beta2", "beta3")]),
  maturities = c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30),
  n_paths = 1000,
  engine = "exact"   # exact transition law; no YUIMA needed
)

# Discount factors and par yields off the fitted curve
factors <- matrix(as.numeric(latest[c("beta1", "beta2", "beta3")]), nrow = 1)
dns_discount_factors(factors, c(1, 5, 10, 30), latest$lambda)
dns_par_yields(factors, c(1, 5, 10, 30), latest$lambda)
```

## Working with live data

``` r

# Set FRED API key (do this once per session)
Sys.setenv(fred_api_key = "your_api_key_here")

library(dnsr)
yield_data <- get_treasury_yield_combined(start_date = "2020-01-01")

# Fit the DNS model
dns_result <- fit_dns_model(yield_data, n_days = 100)

# Visualize the results
plot_yield_curves(yield_data, dns_result)
factor_plots <- plot_factor_time_series(dns_result)
print(factor_plots$factors)

# Estimate factor dynamics
# (month-end sampling by default since 0.2.0; use sampling = "daily" for the
#  0.1.0 behavior)
dynamics <- estimate_dns_dynamics(dns_result$parameters)

# Simulate future yield curves with advanced SDE methods (requires YUIMA)
# get_treasury_yield_*() returns rows most-recent-first and fit_dns_model()
# keeps that order by default, so row 1 is the latest curve. Pass
# date_order = "ascending" if you would rather read it the other way round.
params <- dns_dynamics_to_sim_params(dynamics, lambda = dns_result$parameters$lambda[1])
initial_factors <- dns_result$parameters[1, c("beta1", "beta2", "beta3")]
maturities <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)

sim_result <- simulate_dns_yield_curves(
  params = params, 
  initial_factors = as.numeric(initial_factors),
  maturities = maturities,
  sim_method = "euler"  # Using YUIMA's Euler-Maruyama scheme
)

# Plot simulation results
yield_plots <- plot_simulated_yields(sim_result, time_points = c(0, 1, 5, 10))
print(yield_plots$yield_curves)
```

## Advanced Stochastic Simulation

The package leverages the YUIMA framework for accurate simulation of
stochastic differential equations. This provides:

- **Mathematical Rigor**: Proper implementation of numerical SDE schemes
  (Euler-Maruyama, Milstein). The DNS diffusion does not depend on the
  state, so the two schemes coincide here; accuracy is set by `n_steps`.
- **Correlation Handling**: Accurate simulation of correlated factors
- **Statistical Estimation**: Advanced methods for parameter estimation

The factor system is a linear SDE, so its transition law over any step
is a known Gaussian. `engine = "exact"` samples that law directly — no
discretisation error at any step size, no YUIMA dependency, and about
two orders of magnitude faster than one
[`yuima::simulate()`](https://rdrr.io/pkg/yuima/man/simulate.html) call
per path:

``` r

params <- list(
  kappa = c(0.02, 0.3, 0.5),
  theta = c(0.04, -0.02, 0.01),
  sigma = c(0.005, 0.01, 0.02),
  correlation = diag(3)
)

sim_result <- simulate_dns_yield_curves(
  params = params,
  initial_factors = c(0.045, -0.005, -0.015),
  maturities = c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30),
  n_paths = 1000,
  engine = "exact"
)
```

The default is still `engine = "yuima"`. The two agree in distribution,
but they consume randomness differently, so a given seed does not
reproduce the same paths across engines.

YUIMA also provides the quasi-maximum-likelihood estimator, as an
alternative to the regression in
[`estimate_dns_dynamics()`](https://ambientpointcorp.github.io/dnsr/reference/estimate_dns_dynamics.md):

``` r

# A YUIMA model object, for analysis the package does not wrap
yuima_model <- dns_to_yuima_model(params)

# QMLE on the fitted factor series. A parameters frame can be passed straight
# in; its beta1/beta2/beta3 columns are picked out by name.
est_result <- estimate_dns_parameters_yuima(dns_result$parameters,
                                            delta_t = 1 / 252)
est_result$estimated_params
```

## The Nelson-Siegel Model

The Nelson-Siegel model expresses the yield curve as:

``` math
y(m) = \beta_1 + \beta_2 \left(\frac{1-e^{-\lambda m}}{\lambda m}\right) + \beta_3 \left(\frac{1-e^{-\lambda m}}{\lambda m} - e^{-\lambda m}\right)
```

where: - $`y(m)`$ is the yield at maturity $`m`$ - $`\beta_1`$
represents the long-term level of interest rates - $`\beta_2`$
represents the yield curve slope - $`\beta_3`$ represents the
curvature - $`\lambda`$ controls the decay rate of factor loadings

`λ` is measured per unit of maturity, and every function here takes
maturities in years. The conventional 0.0609 comes from Diebold and Li,
who worked in months, where it places the curvature hump just short of
30 months. In years the same number places it at 29.4, beyond the end of
most curves, where the slope loading is close to the level loading and
the two factors are not separately identified. Since 0.6.0 the default
is `ns_lambda(3.5)`, a hump at three and a half years.
[`ns_lambda()`](https://ambientpointcorp.github.io/dnsr/reference/ns_lambda.md)
and
[`ns_curvature_peak()`](https://ambientpointcorp.github.io/dnsr/reference/ns_lambda.md)
convert between a decay and its hump maturity:

``` r

ns_curvature_peak(0.0609)   # 29.4 -- in years, the hump is off the curve
ns_lambda(2.5)              # 0.717 -- a hump at two and a half years
ns_lambda(3.5)              # 0.512 -- the package default
```

In the dynamic extension, the factors $`\beta_1`$, $`\beta_2`$, and
$`\beta_3`$ are modeled as mean-reverting Ornstein-Uhlenbeck processes:

``` math
d\beta_i(t) = \kappa_i(\theta_i - \beta_i(t))dt + \sigma_i dW_i(t)
```

Two estimators are available.
[`fit_dns_model()`](https://ambientpointcorp.github.io/dnsr/reference/fit_dns_model.md)
followed by
[`estimate_dns_dynamics()`](https://ambientpointcorp.github.io/dnsr/reference/estimate_dns_dynamics.md)
is the two-step approach of Diebold and Li: fit each cross-section, then
model the resulting factors.
[`estimate_dns_state_space()`](https://ambientpointcorp.github.io/dnsr/reference/estimate_dns_state_space.md)
is the one-step state-space approach of Diebold, Rudebusch and Aruoba:
the factors are an unobserved state, filtered out of the yields by
maximum likelihood, with lambda estimated jointly with the dynamics.

``` r

one_step <- estimate_dns_state_space(treasury_yields, lambda = ns_lambda(3.5))
one_step
params <- dns_dynamics_to_sim_params(one_step)   # simulates without conversion
```

## Visualizations

Each plotting function returns `ggplot` objects (or a named list of
them), so output can be themed, combined, or saved by the caller:

| Function | Shows |
|:---|:---|
| [`plot_yield_curves()`](https://ambientpointcorp.github.io/dnsr/reference/plot_yield_curves.md) | Fitted curves against observed yields, with RMSE |
| [`plot_dns_loadings()`](https://ambientpointcorp.github.io/dnsr/reference/plot_dns_loadings.md) | Level, slope, and curvature loadings against maturity |
| [`plot_factor_time_series()`](https://ambientpointcorp.github.io/dnsr/reference/plot_factor_time_series.md) | Estimated factors over time, optionally as changes |
| [`plot_yield_surface()`](https://ambientpointcorp.github.io/dnsr/reference/plot_yield_surface.md) | Yield surface across maturity and time |
| [`plot_yield_heatmap()`](https://ambientpointcorp.github.io/dnsr/reference/plot_yield_heatmap.md) | Yield evolution as a maturity-by-time heatmap |
| [`plot_yield_pca()`](https://ambientpointcorp.github.io/dnsr/reference/plot_yield_pca.md) | Principal components of the yield panel |
| [`plot_simulated_factors()`](https://ambientpointcorp.github.io/dnsr/reference/plot_simulated_factors.md) | Simulated factor paths |
| [`plot_simulated_yields()`](https://ambientpointcorp.github.io/dnsr/reference/plot_simulated_yields.md) | Simulated curves at selected horizons |

## Data

`treasury_yields` holds 1,501 daily U.S. Treasury curves across eleven
maturities, 2 January 2019 to 31 December 2024, in decimal form. The
source is FRED’s H.15 constant maturity series (`DGS1MO` through
`DGS30`), which are U.S. government works in the public domain.
`data-raw/treasury_yields.R` rebuilds it.

## Fitting, Dynamics and Pricing

| Function | Does |
|:---|:---|
| [`fit_ns_model()`](https://ambientpointcorp.github.io/dnsr/reference/fit_ns_model.md) | Fit one cross-section; `method = "profile"` solves the betas in closed form |
| [`fit_dns_model()`](https://ambientpointcorp.github.io/dnsr/reference/fit_dns_model.md) | Fit a panel, selecting the window by date |
| [`ns_lambda()`](https://ambientpointcorp.github.io/dnsr/reference/ns_lambda.md), [`ns_curvature_peak()`](https://ambientpointcorp.github.io/dnsr/reference/ns_lambda.md) | Convert between the decay and the maturity its curvature hump sits at |
| [`estimate_dns_dynamics()`](https://ambientpointcorp.github.io/dnsr/reference/estimate_dns_dynamics.md) | Factor dynamics, diagonal or as a full VAR(1), with standard errors |
| [`dns_dynamics_to_sim_params()`](https://ambientpointcorp.github.io/dnsr/reference/dns_dynamics_to_sim_params.md) | Package dynamics for the simulator |
| [`get_calibrated_dns_params()`](https://ambientpointcorp.github.io/dnsr/reference/get_calibrated_dns_params.md) | Download, fit and return simulation-ready parameters in one call |
| [`simulate_dns_factors()`](https://ambientpointcorp.github.io/dnsr/reference/simulate_dns_factors.md), [`simulate_dns_yield_curves()`](https://ambientpointcorp.github.io/dnsr/reference/simulate_dns_yield_curves.md) | Project forward, via YUIMA or the exact transition law |
| [`calculate_dns_yields()`](https://ambientpointcorp.github.io/dnsr/reference/calculate_dns_yields.md), [`calculate_dns_forwards()`](https://ambientpointcorp.github.io/dnsr/reference/calculate_dns_forwards.md) | Curves and instantaneous forwards from factors |
| [`dns_discount_factors()`](https://ambientpointcorp.github.io/dnsr/reference/dns_discount_factors.md), [`dns_zero_prices()`](https://ambientpointcorp.github.io/dnsr/reference/dns_zero_prices.md), [`dns_par_yields()`](https://ambientpointcorp.github.io/dnsr/reference/dns_par_yields.md) | Discount off a fitted curve |
| [`compare_term_structure_models()`](https://ambientpointcorp.github.io/dnsr/reference/compare_term_structure_models.md) | Nelson-Siegel against Svensson, a spline and Smith-Wilson on one date |
| [`dns_to_yuima_model()`](https://ambientpointcorp.github.io/dnsr/reference/dns_to_yuima_model.md), [`estimate_dns_parameters_yuima()`](https://ambientpointcorp.github.io/dnsr/reference/estimate_dns_parameters_yuima.md) | A YUIMA model object, and QMLE on the factors |
| [`fetch_fred_series()`](https://ambientpointcorp.github.io/dnsr/reference/fetch_fred_series.md), [`get_treasury_yield_data()`](https://ambientpointcorp.github.io/dnsr/reference/get_treasury_yield_data.md), [`get_treasury_yield_combined()`](https://ambientpointcorp.github.io/dnsr/reference/get_treasury_yield_combined.md) | Treasury data from FRED, with retries and caching |

`treasury_yields` is the bundled panel: 1,501 daily curves across eleven
maturities, 2019 to 2024.

## Documentation

For detailed documentation and examples, see the package vignettes:

``` r

# Basic introduction to the DNS model
vignette("dns-model-introduction", package = "dnsr")

# Advanced simulation with YUIMA integration
vignette("dnsr-yuima-integration", package = "dnsr")

# Obtaining and configuring a FRED API key
vignette("fred-api-setup", package = "dnsr")
```

## Applications

The `dnsr` package is useful for:

- **Interest Rate Risk Management**: Model yield curve shifts for risk
  analysis
- **Economic Scenario Generation**: Create consistent interest rate
  scenarios
- **Asset-Liability Management**: Simulate future yield curves for ALM
  studies
- **Fixed Income Research**: Analyze yield curve dynamics and factor
  behavior
- **Stress Testing**: Generate extreme but plausible scenarios for risk
  management

## Citation

If you use `dnsr` in your research, please cite:

    @Manual{,
      title = {dnsr: Dynamic Nelson-Siegel Model in R},
      author = {Ariel R. H.},
      year = {2026},
      note = {R package version 0.4.0},
      url = {https://github.com/ambientpointcorp/dnsr},
    }

## Related

- [`oasr`](https://github.com/ambientpointcorp/oasr) — option-adjusted
  spread analysis; consumes `dnsr` rate paths for cashflow projection
  and OAS.
- *Monte Carlo Fair Value Analysis for Loan Pools* — develops the DNS
  theory this package implements, from the Nelson-Siegel decomposition
  through SDE simulation. Its worked examples run on `dnsr`; the chapter
  code is at
  [mc-fair-value-book](https://github.com/ambientpointcorp/mc-fair-value-book).

## License

MIT © Ambient Point Corp. `LICENSE` carries the copyright stub R
expects; `LICENSE.note` has the full text.
