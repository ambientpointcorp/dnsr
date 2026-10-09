# Dynamic Nelson-Siegel Model Dynamics

Functions to estimate and analyze the dynamics of DNS model factors.

Estimates the stochastic dynamics of the factors in the Dynamic
Nelson-Siegel model based on historical factor time series.

## Usage

``` r
estimate_dns_dynamics(
  params_df,
  min_observations = 20,
  dt = NULL,
  sampling = c("monthly", "daily"),
  structure = c("diagonal", "full")
)
```

## Arguments

- params_df:

  Data frame. The parameters data frame from fit_dns_model() containing
  columns for Date, beta1, beta2, beta3. With `sampling = "daily"` the
  corresponding change columns are also required; with
  `sampling = "monthly"` they are recomputed internally on the month-end
  grid and need not be present.

- min_observations:

  Integer. Minimum number of observations required for estimation,
  counted on the series the regression actually runs on (month-ends when
  `sampling = "monthly"`).

- dt:

  Numeric. Time step in years of the series the regression runs on.
  Defaults to the sampling mode's spacing: 1/12 for `"monthly"`, 1/252
  for `"daily"`. With `sampling = "monthly"` the step is month-end to
  month-end by construction, so a conflicting value is ignored with a
  warning.

- sampling:

  Character. `"monthly"` (default) subsamples the factor series to its
  last observation per month before estimating; `"daily"` estimates on
  the series as given, preserving the pre-0.2.0 behavior.

- structure:

  Character. `"diagonal"` (the default) estimates three independent
  Ornstein-Uhlenbeck processes with correlated innovations, as earlier
  versions did. `"full"` regresses each factor change on the levels of
  all three factors, giving the full mean-reversion matrix of the
  Diebold-Li VAR(1), in which a shock to one factor can pull the others.

## Value

A list containing:

- For each factor (beta1, beta2, beta3):

  - kappa: Mean reversion speed

  - theta: Long-run mean

  - sigma: Volatility

  - model: Fitted regression model

  - half_life: Half-life of mean reversion (in days)

  - summary: Basic statistics about the factor

  - se: Standard errors of kappa, theta and sigma

  - t_value: t statistics for kappa and theta against zero

  - half_life_ci: 95 percent confidence interval for the half-life, in
    years

- correlation: Correlation matrix of factor innovations

- innovations: Matrix of factor innovations (residuals from mean
  reversion)

- innovation_dates: Date vector aligned with the rows of innovations

- estimation_period: Date range used for estimation

- sampling, dt: The sampling mode and time step actually used

- structure: The mean-reversion structure estimated

- kappa_matrix: The 3 by 3 mean-reversion matrix. Diagonal under
  `structure = "diagonal"`; the estimated VAR(1) matrix under `"full"`

- theta_vector: Long-run means as a vector, named by factor

- eigenvalues, system_stability: Eigenvalues of `kappa_matrix` and
  whether every real part is positive

## Details

This function implements the estimation of an Ornstein-Uhlenbeck process
for each factor: \$\$d\beta(t) = \kappa(\theta - \beta(t))dt + \sigma
dW(t)\$\$ The discretized version used for estimation is:
\$\$\Delta\beta(t) = \kappa(\theta - \beta(t-1))\Delta t +
\sigma\epsilon(t)\$\$

The default month-end sampling follows Diebold-Li practice. Daily
cross-sectional refits carry measurement noise that is serially
uncorrelated; estimating the OU regression on daily changes reads that
noise as mean reversion and inflates kappa by an order of magnitude
(errors-in-variables), which in turn makes fixed-step simulation of the
resulting dynamics numerically unstable. Month-end sampling averages the
noise down and recovers reversion speeds on the economically meaningful
time scale. Subsampling is idempotent: passing an already month-end
series through `sampling = "monthly"` reproduces it unchanged.

Standard errors come from the fitted regression. For kappa the mapping
from the slope is linear, so \\se(\kappa) = se(b) / \Delta t\\. For
theta, which is a ratio of coefficients, they come from the delta
method. The one for sigma is the asymptotic \\\sigma / \sqrt{2(n-k)}\\.
Note that kappa carries a well-known upward small-sample bias of order
\\5 / T\\ in a sample of \\T\\ years, which these standard errors do not
capture: they describe the sampling noise, not the bias.

With `structure = "full"` the system estimated is \$\$\Delta\beta(t) =
K(\theta - \beta(t-1))\Delta t + \Sigma\epsilon(t),\$\$ so each change
is regressed on all three lagged levels. Writing the fitted intercepts
as \\a\\ and slopes as \\B\\, the mean-reversion matrix is \\K = -B /
\Delta t\\ and the long-run mean is \\\theta = K^{-1} a / \Delta t\\.
The per-factor `kappa` entries then hold the diagonal of \\K\\; the
cross terms are only in `kappa_matrix`, and
[`simulate_dns_factors`](https://ambientpointcorp.github.io/dnsr/reference/simulate_dns_factors.md)
uses the diagonal.

## Examples

``` r
data(treasury_yields)

# The default monthly sampling estimates on month-end observations, so
# fitting only those is enough and is an order of magnitude quicker.
# warm_start = FALSE makes each cross-section independent of which other
# dates are in the panel, so these fits are identical to the ones fitting
# every day would produce on the same month-ends.
month_ends <- treasury_yields[
  !duplicated(format(treasury_yields$Date, "%Y-%m"), fromLast = TRUE),
]
dns_result <- fit_dns_model(month_ends, n_days = nrow(month_ends),
                            show_progress = FALSE, warm_start = FALSE)
#> Warning: Parameters rest on a bound in some fits: lambda on 13 (of 72 dates). The factor decomposition is not identified on those dates; see the at_bound element of the result, widen lower_bounds/upper_bounds, or check the yield scale.

dynamics <- estimate_dns_dynamics(dns_result$parameters)
dynamics
#> DNS factor dynamics
#>   estimated on     2019-02-28 to 2024-12-31
#>   sampling         monthly (dt = 1/12), 71 observations
#>   structure        diagonal
#> 
#>                 kappa se_kappa   theta  sigma half_life_y
#> level (b1)     2.7753   0.9339  0.0346 0.0304      0.2498
#> slope (b2)     2.6320   0.8909 -0.0083 0.0330      0.2634
#> curvature (b3) 2.1218   0.8195 -0.0043 0.0579      0.3267
#> 
#>   stable           TRUE (eigenvalues 2.775 2.632 2.122)

# Standard errors show how weakly kappa is identified in a short sample
dynamics$beta1$se
#>       kappa       theta       sigma 
#> 0.933902253 0.004531378 0.002584906 
dynamics$beta1$half_life_ci
#>     lower     upper 
#> 0.1504970 0.7335731 

# The full VAR(1), in which a shock to one factor pulls the others
full <- estimate_dns_dynamics(dns_result$parameters, structure = "full")
round(full$kappa_matrix, 3)
#>        beta1  beta2  beta3
#> beta1  3.746 -3.289 -0.031
#> beta2 -5.999  3.251 -1.160
#> beta3  4.821  4.394  4.768
```
