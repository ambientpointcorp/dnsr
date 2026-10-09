# Package index

## Overview

- [`dnsr`](https://ambientpointcorp.github.io/dnsr/reference/dnsr-package.md)
  [`dnsr-package`](https://ambientpointcorp.github.io/dnsr/reference/dnsr-package.md)
  : Dynamic Nelson-Siegel Model for Interest Rate Modeling

## Treasury data

Fetch and assemble yield curve panels from FRED.

- [`fetch_fred_series()`](https://ambientpointcorp.github.io/dnsr/reference/fetch_fred_series.md)
  : Treasury Yield Curve Data Utilities
- [`get_treasury_yield_data()`](https://ambientpointcorp.github.io/dnsr/reference/get_treasury_yield_data.md)
  : Get Treasury Yield Curve Data from FRED
- [`get_treasury_yield_combined()`](https://ambientpointcorp.github.io/dnsr/reference/get_treasury_yield_combined.md)
  : Get Historical Treasury Yield Data from FRED

## Cross-sectional fitting

Fit the Nelson-Siegel curve to one date, or to a panel of dates.

- [`ns_lambda()`](https://ambientpointcorp.github.io/dnsr/reference/ns_lambda.md)
  [`ns_curvature_peak()`](https://ambientpointcorp.github.io/dnsr/reference/ns_lambda.md)
  : Choose the Nelson-Siegel Decay from the Curvature Peak
- [`fit_ns_model()`](https://ambientpointcorp.github.io/dnsr/reference/fit_ns_model.md)
  : Dynamic Nelson-Siegel Model Implementation
- [`fit_dns_model()`](https://ambientpointcorp.github.io/dnsr/reference/fit_dns_model.md)
  : Fit Dynamic Nelson-Siegel Model to Time Series of Yield Curves
- [`compare_term_structure_models()`](https://ambientpointcorp.github.io/dnsr/reference/compare_term_structure_models.md)
  : Compare Nelson-Siegel Model to Other Term Structure Models

## Factor dynamics

Estimate the Ornstein-Uhlenbeck dynamics of the fitted factors.

- [`estimate_dns_dynamics()`](https://ambientpointcorp.github.io/dnsr/reference/estimate_dns_dynamics.md)
  : Dynamic Nelson-Siegel Model Dynamics
- [`dns_dynamics_to_sim_params()`](https://ambientpointcorp.github.io/dnsr/reference/dns_dynamics_to_sim_params.md)
  : Convert DNS Dynamics to Simulation Parameters
- [`get_calibrated_dns_params()`](https://ambientpointcorp.github.io/dnsr/reference/get_calibrated_dns_params.md)
  : Get Calibrated DNS Parameters
- [`estimate_dns_parameters_yuima()`](https://ambientpointcorp.github.io/dnsr/reference/estimate_dns_parameters_yuima.md)
  : Estimate DNS Parameters from Factors using YUIMA

## Simulation

Project factors forward and convert them to curves.

- [`simulate_dns_factors()`](https://ambientpointcorp.github.io/dnsr/reference/simulate_dns_factors.md)
  : Dynamic Nelson-Siegel Model Simulation using YUIMA
- [`simulate_dns_yield_curves()`](https://ambientpointcorp.github.io/dnsr/reference/simulate_dns_yield_curves.md)
  : Simulate DNS Yield Curves using YUIMA
- [`calculate_dns_yields()`](https://ambientpointcorp.github.io/dnsr/reference/calculate_dns_yields.md)
  : Calculate Yield Curves from DNS Factors
- [`calculate_dns_forwards()`](https://ambientpointcorp.github.io/dnsr/reference/calculate_dns_forwards.md)
  : Calculate Forward Rates from DNS Factors
- [`dns_to_yuima_model()`](https://ambientpointcorp.github.io/dnsr/reference/dns_to_yuima_model.md)
  : Convert Nelson-Siegel Model to YUIMA SDE

## Pricing

Discount factors, zero prices and par yields off a fitted curve.

- [`dns_discount_factors()`](https://ambientpointcorp.github.io/dnsr/reference/dns_discount_factors.md)
  : Discount Factors, Zero Prices and Par Yields from DNS Factors
- [`dns_zero_prices()`](https://ambientpointcorp.github.io/dnsr/reference/dns_zero_prices.md)
  : Zero-Coupon Bond Prices from DNS Factors
- [`dns_par_yields()`](https://ambientpointcorp.github.io/dnsr/reference/dns_par_yields.md)
  : Par Coupon Yields from DNS Factors

## Data

- [`treasury_yields`](https://ambientpointcorp.github.io/dnsr/reference/treasury_yields.md)
  : U.S. Treasury Yield Curves, 2019-2024

## Visualization

- [`plot_yield_curves()`](https://ambientpointcorp.github.io/dnsr/reference/plot_yield_curves.md)
  : Dynamic Nelson-Siegel Model Visualization Functions
- [`plot_dns_loadings()`](https://ambientpointcorp.github.io/dnsr/reference/plot_dns_loadings.md)
  : Plot DNS Factor Loadings
- [`plot_factor_time_series()`](https://ambientpointcorp.github.io/dnsr/reference/plot_factor_time_series.md)
  : Plot DNS Factor Time Series
- [`plot_yield_heatmap()`](https://ambientpointcorp.github.io/dnsr/reference/plot_yield_heatmap.md)
  : Create Heat Map of Yield Curve Evolution
- [`plot_yield_surface()`](https://ambientpointcorp.github.io/dnsr/reference/plot_yield_surface.md)
  : Plot 3D Yield Surface
- [`plot_yield_pca()`](https://ambientpointcorp.github.io/dnsr/reference/plot_yield_pca.md)
  : Plot Principal Component Analysis of Yield Curve Movements
- [`plot_simulated_factors()`](https://ambientpointcorp.github.io/dnsr/reference/plot_simulated_factors.md)
  : Plot Simulated DNS Factors
- [`plot_simulated_yields()`](https://ambientpointcorp.github.io/dnsr/reference/plot_simulated_yields.md)
  : Plot Simulated Yield Curves
