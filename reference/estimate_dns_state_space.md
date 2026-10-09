# Estimate the Dynamic Nelson-Siegel Model in One Step

Fits the DNS model to a panel of yield curves by maximum likelihood,
with the factors treated as an unobserved state and filtered out of the
yields, rather than fitted curve by curve and then modelled. This is the
one-step state-space approach of Diebold, Rudebusch and Aruoba (2006),
as against the two-step approach of Diebold and Li (2006) that
[`fit_dns_model`](https://ambientpointcorp.github.io/dnsr/reference/fit_dns_model.md)
and
[`estimate_dns_dynamics`](https://ambientpointcorp.github.io/dnsr/reference/estimate_dns_dynamics.md)
implement between them.

## Usage

``` r
estimate_dns_state_space(
  yield_data,
  n_days = NULL,
  sampling = c("monthly", "all"),
  lambda = NULL,
  measurement = c("per_maturity", "common"),
  dt = NULL,
  start = NULL,
  maxit = 500L,
  verbose = FALSE
)
```

## Arguments

- yield_data:

  Data frame. Yield curve panel with a `Date` column and one column per
  maturity, named by maturity in years – the same shape
  [`fit_dns_model`](https://ambientpointcorp.github.io/dnsr/reference/fit_dns_model.md)
  takes.

- n_days:

  Integer or NULL. Use only the most recent `n_days` dates. NULL (the
  default) uses all of them.

- sampling:

  Character. `"monthly"` (the default) keeps the last observation in
  each calendar month, the frequency this model is usually estimated at;
  `"all"` uses every date.

- lambda:

  Numeric or NULL. NULL (the default) estimates lambda jointly with the
  dynamics. A number holds it fixed, which is faster.

- measurement:

  Character. `"per_maturity"` (the default) gives each maturity its own
  measurement error variance, as in Diebold, Rudebusch and Aruoba;
  `"common"` estimates a single variance for all maturities, which is
  better determined on short panels.

- dt:

  Numeric or NULL. Time step between observations, in years. NULL (the
  default) infers it from the median spacing of the dates used.

- start:

  List or NULL. Optional starting values, with any of the elements
  `kappa`, `theta`, `sigma`, `correlation`, `lambda`, `measurement_sd`.
  NULL (the default) starts from the two-step estimates, which are also
  reported on the result.

- maxit:

  Integer. Maximum optimiser iterations. Default 500.

- verbose:

  Logical. If TRUE, report the log-likelihood as it improves.

## Value

A list with class `"dns_state_space"`. Its `beta1`/`beta2`/`beta3`,
`correlation` and `lambda` elements are shaped so that
[`dns_dynamics_to_sim_params`](https://ambientpointcorp.github.io/dnsr/reference/dns_dynamics_to_sim_params.md)
accepts it directly:

- beta1, beta2, beta3: each a list of `kappa`, `theta` and `sigma`, plus
  the implied `half_life`

- correlation: the 3-by-3 shock correlation matrix

- lambda: the decay parameter, estimated or as supplied

- factors: data frame of smoothed factors by date, directly comparable
  with the `parameters` frame
  [`fit_dns_model()`](https://ambientpointcorp.github.io/dnsr/reference/fit_dns_model.md)
  returns

- filtered: the filtered factors, using only data up to each date

- measurement_sd: fitted measurement error per maturity, in the units of
  the input yields

- logLik, n_params, aic, bic: fit and its cost in parameters

- fitted, residuals: fitted yields from the smoothed factors, and
  observed minus fitted

- rmse: root mean squared measurement residual

- kappa_matrix, theta_vector, eigenvalues, system_stability: as
  [`estimate_dns_dynamics`](https://ambientpointcorp.github.io/dnsr/reference/estimate_dns_dynamics.md)
  reports them

- two_step: the two-step estimates used as starting values

- convergence, message: from
  [`optim`](https://rdrr.io/r/stats/optim.html)

- dt, sampling, structure, n_dates, maturities

## Details

The model is \$\$y_t = \Lambda(\lambda) f_t + e_t, \quad e_t \sim N(0,
H)\$\$ \$\$f_t = \theta + A (f\_{t-1} - \theta) + \eta_t, \quad \eta_t
\sim N(0, Q)\$\$ with \\\Lambda\\ the Nelson-Siegel loadings and \\H\\
diagonal. \\A\\ and \\Q\\ are the transition law of the
Ornstein-Uhlenbeck system the rest of the package simulates: \$\$A\_{ii}
= e^{-\kappa_i \Delta t}\$\$ \$\$Q\_{ij} = \rho\_{ij} \sigma_i \sigma_j
(1 - e^{-(\kappa_i+\kappa_j)\Delta t}) / (\kappa_i + \kappa_j)\$\$

\\A\\ and \\Q\\ are therefore functions of \\\kappa\\, \\\sigma\\ and
\\\rho\\ rather than free parameters. This has three consequences. \\Q\\
is a covariance matrix at every trial point. \\\kappa \> 0\\ makes the
system stationary, so the filter has a proper initial covariance. And
the result is in the parameters
[`dns_dynamics_to_sim_params`](https://ambientpointcorp.github.io/dnsr/reference/dns_dynamics_to_sim_params.md)
and
[`simulate_dns_factors`](https://ambientpointcorp.github.io/dnsr/reference/simulate_dns_factors.md)
accept, so a fit can be simulated forward without conversion. The cost
is that cross-factor terms cannot be represented.

**One step against two.** The two-step estimator uses the fitted betas
as data. They are estimates carrying sampling error, and using them as
observed in a time-series regression biases mean reversion upward, the
errors-in-variables direction. Filtering the factors propagates that
uncertainty into the likelihood and estimates \\\lambda\\ jointly with
the dynamics. The costs are run time, needing the whole panel at once,
and a harder numerical problem: a dozen or more parameters in one
optimisation rather than a sequence of independent three-parameter fits.

**Limitations.** The factor system is diagonal: each factor reverts to
its own mean at its own speed, with contemporaneously correlated shocks.
A full mean-reversion matrix is not supported. Supporting it would
require a matrix exponential for \\A\\ and a Sylvester solve for \\Q\\,
and
[`dns_dynamics_to_sim_params`](https://ambientpointcorp.github.io/dnsr/reference/dns_dynamics_to_sim_params.md)
discards cross terms, so the result could not be simulated by the rest
of the package. Use `estimate_dns_dynamics(structure = "full")` on
two-step betas for the cross terms.

Missing yields shrink the measurement equation for that date rather than
dropping the date, so a curve missing one maturity still contributes its
other maturities and leaves no gap in the series.

**Precision of the estimates.** The parameters divide into two groups.
Simulating from a known DNS state space with monthly data and ten
maturities, then re-estimating, gives these percentage errors:

- \\\lambda\\: within 1 percent at 20, 40 and 100 years of data.

- \\\sigma\\: within 5 percent at all three lengths.

- measurement standard deviation: within 1 percent at all three.

- \\\kappa\\ for the level factor: 167 percent too high at 20 years, 77
  percent at 40, 58 percent at 100.

- \\\theta\\ for the slope factor: 54 percent low at 20 years, 66
  percent high at 40, 2 percent low at 100.

\\\lambda\\, the shock volatilities and the measurement error are
recovered precisely at realistic sample lengths, and the smoothed
factors correlate above 0.99 with the simulated ones. The mean reversion
speeds and long-run means are not. \\\kappa\\ carries an upward bias
that decays as the sample lengthens, and \\\theta\\ has a standard error
of the same order as its own value, because a factor with a multi-year
half-life supplies few independent observations about its long-run mean
however finely it is sampled.

The errors shrink toward zero as the sample grows, so this reflects the
difficulty of estimating mean reversion rather than an error in the
implementation, and it applies to the two-step route as well. Treat
\\\kappa\\ and \\\theta\\ as approximate.

## References

Diebold, F.X., Rudebusch, G.D. and Aruoba, S.B. (2006). The macroeconomy
and the yield curve: a dynamic latent factor approach. *Journal of
Econometrics* 131, 309-338.

Diebold, F.X. and Li, C. (2006). Forecasting the term structure of
government bond yields. *Journal of Econometrics* 130, 337-364.

## See also

[`fit_dns_model`](https://ambientpointcorp.github.io/dnsr/reference/fit_dns_model.md)
and
[`estimate_dns_dynamics`](https://ambientpointcorp.github.io/dnsr/reference/estimate_dns_dynamics.md)
for the two-step route,
[`dns_dynamics_to_sim_params`](https://ambientpointcorp.github.io/dnsr/reference/dns_dynamics_to_sim_params.md)
to simulate a fit forward.

## Examples

``` r
# Maximising the likelihood takes a few seconds, so the example is not run
# during checks. The run time is dominated by the number of optimiser
# iterations rather than by the size of the panel.
# \donttest{
# Two years of the bundled panel, thinned to month-ends. lambda is held
# fixed here; pass lambda = NULL to estimate it jointly with the dynamics.
# n_days counts rows of the input, which is daily, so 500 trading days
# leaves 24 month-ends to filter.
fit <- estimate_dns_state_space(treasury_yields, n_days = 500,
                                lambda = ns_lambda(3.5),
                                measurement = "common")
fit
#> DNS one-step state-space fit (Kalman filter, maximum likelihood)
#>   estimated on     2023-01-31 to 2024-12-31
#>   sampling         monthly (dt = 1/11.8), 24 dates, 11 maturities
#>   lambda           0.5124 fixed (hump at 3.50 years)
#> 
#>                 kappa   theta  sigma  half_life_y
#> level (b1)     0.0000  0.0436 0.0000 1.819697e+09
#> slope (b2)     0.7682  0.0064 0.0066 9.023000e-01
#> curvature (b3) 7.9965 -0.0200 0.0535 8.670000e-02
#> 
#>   measurement sd   20.84 bp median across maturities
#>   fit              19.39 bp rmse
#>   logLik           1209.63 on 13 parameters (AIC -2393.3, BIC -2346.8)
#>   stable           TRUE
#>   $two_step holds the two-step estimates used as starting values

# The shock correlation, the measurement error by maturity, and the
# two-step estimates alongside the one-step ones
summary(fit)
#> DNS one-step state-space fit (Kalman filter, maximum likelihood)
#>   estimated on     2023-01-31 to 2024-12-31
#>   sampling         monthly (dt = 1/11.8), 24 dates, 11 maturities
#>   lambda           0.5124 fixed (hump at 3.50 years)
#> 
#>                 kappa   theta  sigma  half_life_y
#> level (b1)     0.0000  0.0436 0.0000 1.819697e+09
#> slope (b2)     0.7682  0.0064 0.0066 9.023000e-01
#> curvature (b3) 7.9965 -0.0200 0.0535 8.670000e-02
#> 
#>   measurement sd   20.84 bp median across maturities
#>   fit              19.39 bp rmse
#>   logLik           1209.63 on 13 parameters (AIC -2393.3, BIC -2346.8)
#>   stable           TRUE
#>   $two_step holds the two-step estimates used as starting values
#> 
#> Shock correlation
#>        beta1  beta2  beta3
#> beta1  1.000  0.372 -1.000
#> beta2  0.372  1.000 -0.374
#> beta3 -1.000 -0.374  1.000
#> 
#> Measurement error by maturity (bp)
#> 0.0833   0.25    0.5      1      2      3      5      7     10     20     30 
#>  20.84  20.84  20.84  20.84  20.84  20.84  20.84  20.84  20.84  20.84  20.84 
#> 
#> Residual rmse by maturity (bp)
#> 0.0833   0.25    0.5      1      2      3      5      7     10     20     30 
#>  27.10   7.16  14.05  21.47  14.25   8.65   6.50   9.94  20.19  29.87  31.21 
#> 
#> One step against the two-step estimates it started from
#>                one_kappa one_theta one_sigma two_kappa two_theta two_sigma
#> level (b1)        0.0000    0.0436    0.0000    2.6700    0.0463    0.0096
#> slope (b2)        0.7682    0.0064    0.0066    0.0012  -11.5237    0.0114
#> curvature (b3)    7.9965   -0.0200    0.0535    7.1414   -0.0210    0.0413
#> 
#> The two-step estimates treat the fitted betas as observed when
#> they carry sampling error, which biases kappa upward. See
#> ?estimate_dns_state_space on the precision of kappa and theta.

# The smoothed factors have the same shape as the two-step betas
head(fit$factors)
#>         Date      beta1       beta2       beta3
#> 1 2023-01-31 0.04363048 0.007135860 -0.03178261
#> 2 2023-02-28 0.04363048 0.007546616 -0.01451132
#> 3 2023-03-31 0.04363048 0.008568861 -0.03670761
#> 4 2023-04-28 0.04363048 0.009862089 -0.04029340
#> 5 2023-05-31 0.04363048 0.013336431 -0.03555936
#> 6 2023-06-30 0.04363048 0.013779761 -0.02407746

# The fit is accepted by the simulation parameter converter
params <- dns_dynamics_to_sim_params(fit)
str(params)
#> List of 7
#>  $ kappa       : num [1:3] 3.81e-10 7.68e-01 8.00
#>  $ theta       : num [1:3] 0.0436 0.0064 -0.02
#>  $ sigma       : num [1:3] 2.28e-10 6.59e-03 5.35e-02
#>  $ names       : chr [1:3] "Level" "Slope" "Curvature"
#>  $ lambda      : num 0.512
#>  $ correlation : num [1:3, 1:3] 1 0.372 -1 0.372 1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:3] "beta1" "beta2" "beta3"
#>   .. ..$ : chr [1:3] "beta1" "beta2" "beta3"
#>  $ kappa_matrix: num [1:3, 1:3] 3.81e-10 0.00 0.00 0.00 7.68e-01 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:3] "beta1" "beta2" "beta3"
#>   .. ..$ : chr [1:3] "beta1" "beta2" "beta3"
# }
```
