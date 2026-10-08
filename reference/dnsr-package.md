# Dynamic Nelson-Siegel Model for Interest Rate Modeling

A comprehensive toolkit for modeling, fitting, and simulating interest
rates using the Dynamic Nelson-Siegel (DNS) model. The package provides
functions to calibrate the DNS model from yield curve data, estimate
factor dynamics, simulate yield curves, and visualize results.

## Details

The Dynamic Nelson-Siegel model is a popular term structure model used
in finance for modeling the yield curve. It represents the yield curve
using three latent factors (level, slope, and curvature) that evolve
over time according to a mean-reverting stochastic process.

This package provides a comprehensive implementation of the DNS model,
including:

- Data retrieval and processing of U.S. Treasury yield curve data

- Cross-sectional fitting of the Nelson-Siegel model to yield curves

- Time-series estimation of factor dynamics

- Stochastic simulation of interest rate paths

- Visualization tools for yield curves and factor analysis

The package is designed to be used by financial analysts, risk managers,
and researchers for applications such as:

- Interest rate forecasting and scenario analysis

- Risk management and stress testing

- Pricing and valuation of fixed income securities

- Asset-liability management

## Key Functions

- [`get_treasury_yield_data`](https://ambientpointcorp.github.io/dnsr/reference/get_treasury_yield_data.md):

  Download and process Treasury yield curve data

- [`get_treasury_yield_combined`](https://ambientpointcorp.github.io/dnsr/reference/get_treasury_yield_combined.md):

  Get historical and current Treasury yield data

- [`fit_ns_model`](https://ambientpointcorp.github.io/dnsr/reference/fit_ns_model.md):

  Fit Nelson-Siegel model to a single yield curve

- [`fit_dns_model`](https://ambientpointcorp.github.io/dnsr/reference/fit_dns_model.md):

  Fit Dynamic Nelson-Siegel model to time series of yield curves

- [`estimate_dns_dynamics`](https://ambientpointcorp.github.io/dnsr/reference/estimate_dns_dynamics.md):

  Estimate dynamics of the DNS factors

- [`get_calibrated_dns_params`](https://ambientpointcorp.github.io/dnsr/reference/get_calibrated_dns_params.md):

  Get calibrated DNS parameters from market data

- [`simulate_dns_factors`](https://ambientpointcorp.github.io/dnsr/reference/simulate_dns_factors.md):

  Simulate DNS factors using YUIMA package

- [`calculate_dns_yields`](https://ambientpointcorp.github.io/dnsr/reference/calculate_dns_yields.md):

  Calculate yield curves from DNS factors

- [`simulate_dns_yield_curves`](https://ambientpointcorp.github.io/dnsr/reference/simulate_dns_yield_curves.md):

  Simulate DNS yield curves

- [`plot_yield_curves`](https://ambientpointcorp.github.io/dnsr/reference/plot_yield_curves.md):

  Plot yield curves and fitted DNS models

- [`plot_factor_time_series`](https://ambientpointcorp.github.io/dnsr/reference/plot_factor_time_series.md):

  Plot DNS factor time series

## Model Background

The Nelson-Siegel model expresses the yield curve as: \$\$y(m) =
\beta_1 + \beta_2 \left(\frac{1-e^{-\lambda m}}{\lambda m}\right) +
\beta_3 \left(\frac{1-e^{-\lambda m}}{\lambda m} - e^{-\lambda
m}\right)\$\$ where \\m\\ is maturity, \\\beta_1\\, \\\beta_2\\, and
\\\beta_3\\ represent the level, slope, and curvature factors,
respectively, and \\\lambda\\ controls the decay rate of factor
loadings.

In the Dynamic Nelson-Siegel model, these factors follow a
mean-reverting stochastic process: \$\$d\beta(t) = \kappa(\theta -
\beta(t))dt + \sigma dW(t)\$\$ where \\\kappa\\ is the speed of mean
reversion, \\\theta\\ is the long-run mean, \\\sigma\\ is the
volatility, and \\W(t)\\ is a Wiener process.

## See also

Useful links:

- <https://github.com/ambientpointcorp/dnsr>

- Report bugs at <https://github.com/ambientpointcorp/dnsr/issues>

## Author

**Maintainer**: Ariel R. H. <ariel@ambientpointcorp.com>

## Examples

``` r
if (FALSE) { # \dontrun{
# Get Treasury yield data
yield_data <- get_treasury_yield_combined("2020-01-01")

# Fit DNS model
dns_result <- fit_dns_model(yield_data, n_days = 100)

# Estimate factor dynamics
dynamics <- estimate_dns_dynamics(dns_result$parameters)

# Simulate future yield curves, starting from the most recent fitted factors
latest <- dns_result$parameters[which.max(dns_result$parameters$Date), ]
params <- dns_dynamics_to_sim_params(dynamics, lambda = latest$lambda)
sim_result <- simulate_dns_yield_curves(
  params = params,
  initial_factors = as.numeric(latest[c("beta1", "beta2", "beta3")]),
  maturities = c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30),
  n_paths = 100,
  engine = "exact"
)

# Visualize results
plot_yield_curves(yield_data, dns_result)
plot_factor_time_series(dns_result)
plot_simulated_yields(sim_result)
} # }
```
