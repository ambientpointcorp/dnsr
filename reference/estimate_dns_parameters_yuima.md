# Estimate DNS Parameters from Factors using YUIMA

Estimates DNS model parameters from factor time series using YUIMA's
advanced estimation methods.

## Usage

``` r
estimate_dns_parameters_yuima(
  factor_data,
  initial_params = NULL,
  estimation_method = "qmle",
  delta_t = 1/252
)
```

## Arguments

- factor_data:

  Data frame or matrix. Time series of DNS factors (level, slope,
  curvature). Columns named `beta1`, `beta2` and `beta3` are used when
  present, so a `parameters` frame from
  [`fit_dns_model`](https://ambientpointcorp.github.io/dnsr/reference/fit_dns_model.md)
  can be passed directly; otherwise the first three numeric columns are
  taken, after any `Date` column is dropped.

- initial_params:

  List. Initial parameter values for optimization.

- estimation_method:

  Character. Method for parameter estimation ("qmle" or "gmm").

- delta_t:

  Numeric. Time step between observations in years.

## Value

A list containing:

- estimated_params: Estimated model parameters

- yuima_model: The fitted YUIMA model

- fit_result: Full result from the YUIMA estimation

## Details

This function uses the YUIMA package's estimation capabilities to fit an
Ornstein-Uhlenbeck process to observed factor time series. It provides
more rigorous statistical estimation than standard regression methods.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming factor_data is a data frame with columns: Date, beta1, beta2, beta3
factors <- factor_data[, c("beta1", "beta2", "beta3")]
initial <- list(
  kappa = c(0.02, 0.3, 0.5),
  theta = c(0.04, -0.02, 0.01),
  sigma = c(0.005, 0.01, 0.02)
)
est_result <- estimate_dns_parameters_yuima(factors, initial)
} # }
```
