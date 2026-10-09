# Dynamic Nelson-Siegel Model Implementation

Functions to calibrate and fit the Dynamic Nelson-Siegel (DNS) model to
yield curve data.

Fits the Nelson-Siegel model to a single yield curve (cross-section) for
a specific date.

## Usage

``` r
fit_ns_model(
  maturities,
  yields,
  initial_lambda = 0.0609,
  lambda_constraint = "constrained",
  lower_bounds = c(-0.1, -0.5, -0.5, 0.01),
  upper_bounds = c(0.2, 0.5, 0.5, 0.5),
  optimizer = "L-BFGS-B",
  method = c("optim", "profile")
)
```

## Arguments

- maturities:

  Numeric vector. The maturities (in years) for the yield curve points.

- yields:

  Numeric vector. The yield values corresponding to the maturities.

- initial_lambda:

  Numeric. Initial value for the decay parameter lambda, per unit of
  maturity, and the value used outright when
  `lambda_constraint = "fixed"`. Maturities here are in years, so the
  0.0609 default is the Diebold-Li value for maturities in months and
  places the curvature hump at 29.4 years rather than months; see
  [`ns_lambda`](https://ambientpointcorp.github.io/dnsr/reference/ns_lambda.md).

- lambda_constraint:

  Character. One of "fixed", "constrained", or "free":

  - "fixed": lambda is fixed at the initial_lambda value.

  - "constrained": lambda is optimized within bounds.

  - "free": lambda is optimized without constraints.

  With `method = "optim"` the search over lambda is weak: L-BFGS-B's
  numerical gradient is poorly scaled when the betas are around 0.05 and
  lambda around 0.5, and in practice lambda moves very little from
  `initial_lambda` even when unconstrained. `method = "profile"`
  searches lambda properly and reaches the same optimum from any
  starting value, so prefer it when lambda is meant to be estimated
  rather than supplied.

- lower_bounds:

  Numeric vector of length 4. Lower bounds for parameters (beta1, beta2,
  beta3, lambda).

- upper_bounds:

  Numeric vector of length 4. Upper bounds for parameters (beta1, beta2,
  beta3, lambda).

- optimizer:

  Character. The optimization method to use. Default is "L-BFGS-B". Used
  by `method = "optim"` only.

- method:

  Character. How to solve the fit. `"optim"` (the default) searches all
  four parameters numerically, as in earlier versions. `"profile"`
  exploits the fact that the curve is linear in the betas once lambda is
  fixed: the betas are solved by least squares in closed form and only
  lambda is searched, over one dimension. The profiled fit does not
  depend on beta start values and does not need beta bounds, which
  removes the boundary solutions `"optim"` can return; it is the
  recommended setting and is expected to become the default.

## Value

A list containing:

- parameters: Named vector of fitted parameters (beta1, beta2, beta3,
  lambda)

- fitted_yields: Vector of fitted yield values, one per retained
  maturity

- rmse: Root mean squared error of the fit

- maturities: The retained maturities (input maturities, minus any
  dropped for being missing or non-positive)

- actual_yields: The retained input yields

- residuals: The fitting residuals (actual - fitted)

- success: Logical indicating if the optimization converged

- valid_idx: Logical vector over the *input* maturities marking which
  points were retained, so callers can place the fitted values back into
  a full-width row

- at_bound: Named logical vector marking parameters that came to rest on
  one of their bounds

## Details

The Nelson-Siegel model expresses the yield curve as: \$\$y(m) =
\beta_1 + \beta_2 \left(\frac{1-e^{-\lambda m}}{\lambda m}\right) +
\beta_3 \left(\frac{1-e^{-\lambda m}}{\lambda m} - e^{-\lambda
m}\right)\$\$ where \\m\\ is maturity, \\\beta_1\\, \\\beta_2\\, and
\\\beta_3\\ represent the level, slope, and curvature factors,
respectively, and \\\lambda\\ controls the decay rate of factor
loadings.

With `method = "profile"` the sum of squared errors is concentrated: for
any lambda the optimal betas are the least-squares solution of \\y =
X(\lambda)\beta\\, so the objective reduces to a function of lambda
alone, minimised over `[lower_bounds[4], upper_bounds[4]]` by
[`optimize`](https://rdrr.io/r/stats/optimize.html). Beta bounds are not
consulted, since the linear solve is exact.

The default bounds are written for yields in decimal form. When they are
left at their defaults and the data is evidently on another scale
(yields in percentage points, say), the beta box is widened to cover the
data rather than clamping the level factor: L-BFGS-B reports convergence
at a bound, so a clamped fit would otherwise be returned with
`success = TRUE`. Any parameter that does come to rest on a bound is
reported in `at_bound` and raises a warning.

## Examples

``` r
maturities <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)
yields <- c(0.0050, 0.0052, 0.0055, 0.0065, 0.0078, 0.0105, 0.0128, 0.0152, 0.0215, 0.0225)
fit_ns_model(maturities, yields)
#> Nelson-Siegel curve fit
#>   beta1 (level)     0.02405
#>   beta2 (slope)     -0.02000
#>   beta3 (curvature) 0.02439
#>   lambda            0.07177
#>   rmse              4.66 bp over 10 maturities

# The profiled fit solves the betas in closed form at each lambda, so it
# needs no start values and cannot return a clamped beta.
fit_ns_model(maturities, yields, method = "profile")
#> Nelson-Siegel curve fit
#>   beta1 (level)     0.02707
#>   beta2 (slope)     -0.02216
#>   beta3 (curvature) -0.02021
#>   lambda            0.33811
#>   rmse              2.89 bp over 10 maturities

# A real curve from the bundled panel
data(treasury_yields)
latest <- treasury_yields[nrow(treasury_yields), ]
fit_ns_model(as.numeric(setdiff(names(treasury_yields), "Date")),
             as.numeric(latest[-1]), method = "profile")
#> Warning: Nelson-Siegel fit rests on a parameter bound (lambda); the factor decomposition is not identified. Widen lower_bounds/upper_bounds or check the yield scale.
#> Nelson-Siegel curve fit
#>   beta1 (level)     0.04988
#>   beta2 (slope)     -0.00617
#>   beta3 (curvature) -0.01404
#>   lambda            0.50000
#>   rmse              4.95 bp over 11 maturities
#>   at bound          lambda (decomposition not identified)
```
