# Convert Nelson-Siegel Model to YUIMA SDE

Converts DNS model parameters to a YUIMA SDE object for advanced
analysis.

## Usage

``` r
dns_to_yuima_model(params, correlation = NULL)
```

## Arguments

- params:

  List. Model parameters (kappa, theta, sigma, correlation).

- correlation:

  Matrix. Correlation matrix for factor innovations.

## Value

A YUIMA model object representing the DNS factor dynamics.

## Details

This function creates a YUIMA model object that represents the factor
dynamics of the DNS model. This object can be used for various analyses
supported by the YUIMA package, such as parameter estimation, filtering,
and more advanced simulations.

## Examples

``` r
if (FALSE) { # \dontrun{
params <- list(
  kappa = c(0.02, 0.3, 0.5),
  theta = c(0.04, -0.02, 0.01),
  sigma = c(0.005, 0.01, 0.02)
)
correlation <- matrix(c(1, 0.2, 0.1, 0.2, 1, 0.3, 0.1, 0.3, 1), nrow = 3)
yuima_model <- dns_to_yuima_model(params, correlation)
} # }
```
