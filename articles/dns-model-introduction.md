# Introduction to the Dynamic Nelson-Siegel Model

## Overview

The `dnsr` package implements the Dynamic Nelson-Siegel (DNS) model for
yield curve modelling. This vignette covers the workflow: fitting the
curve, estimating how its factors move, simulating forward, and pricing
off the result.

Everything here runs against `treasury_yields`, a six-year panel of real
U.S. Treasury curves that ships with the package, so none of it needs a
FRED API key. For live data see
[`vignette("fred-api-setup", package = "dnsr")`](https://ambientpointcorp.github.io/dnsr/articles/fred-api-setup.md),
and for the YUIMA simulation engine see
[`vignette("dnsr-yuima-integration", package = "dnsr")`](https://ambientpointcorp.github.io/dnsr/articles/dnsr-yuima-integration.md).

``` r

library(dnsr)
data(treasury_yields)

dim(treasury_yields)
#> [1] 1501   12
range(treasury_yields$Date)
#> [1] "2019-01-02" "2024-12-31"
```

The frame is the shape every function in the package expects: a `Date`
column followed by one column per maturity, named by maturity in years,
with yields in decimal form.

``` r

treasury_yields[1:3, 1:6]
#>         Date 0.0833   0.25    0.5      1      2
#> 1 2019-01-02 0.0240 0.0242 0.0251 0.0260 0.0250
#> 2 2019-01-03 0.0242 0.0241 0.0247 0.0250 0.0239
#> 3 2019-01-04 0.0240 0.0242 0.0251 0.0257 0.0250
```

## The Nelson-Siegel Model

The Nelson-Siegel model is a parsimonious parametric model for the yield
curve. It expresses the yield at each maturity through three factors:

``` math
y(m) = \beta_1 + \beta_2 \left(\frac{1-e^{-\lambda m}}{\lambda m}\right) + \beta_3 \left(\frac{1-e^{-\lambda m}}{\lambda m} - e^{-\lambda m}\right)
```

where

- $`y(m)`$ is the yield at maturity $`m`$,
- $`\beta_1`$ is the long-term level of interest rates,
- $`\beta_2`$ is the slope, short-term minus long-term,
- $`\beta_3`$ is the curvature, the size of the hump,
- $`\lambda`$ controls how fast the factor loadings decay.

The loadings are what give each factor its interpretation: the level
loads equally on every maturity, the slope loads most at the short end
and decays to zero, and the curvature peaks somewhere in the middle.

### The units of lambda

$`\lambda`$ appears only as the product $`\lambda m`$, so it is measured
per unit of maturity. Its numerical value therefore depends on whether
maturities are years or months, and a $`\lambda`$ quoted on its own does
not say which. Every function in this package takes maturities in years.

Diebold and Li chose $`\lambda = 0.0609`$ working in *months*, where it
places the curvature hump just short of 30 months. Against maturities in
years the same number places it elsewhere:

``` r

ns_curvature_peak(0.0609)   # where the hump sits, in the maturity's own units
#> [1] 29.44634
```

Twenty-nine years, beyond the end of most curves, where the curvature
loading is close to flat.
[`ns_lambda()`](https://ambientpointcorp.github.io/dnsr/reference/ns_lambda.md)
maps in the other direction, from a hump maturity to the decay that
places it there:

``` r

ns_lambda(2.5)   # a hump at two and a half years
#> [1] 0.7173128
ns_lambda(30)    # the Diebold-Li criterion, in the months they worked in
#> [1] 0.05977607
```

`0.0609` is therefore correct for month-denominated maturities and
incorrect for maturities in years. The package default is
`ns_lambda(3.5)`, a hump at three and a half years, near the middle of a
Treasury curve. The two loadings:

``` r

plot_dns_loadings(lambda = 0.0609)
```

![](dns-model-introduction_files/figure-html/loadings-1.png)

``` r

plot_dns_loadings(lambda = ns_lambda(2.5))
```

![](dns-model-introduction_files/figure-html/loadings_fast-1.png)

The first has no hump within thirty years. The second has the shape the
three-factor description implies.

The consequence is in the fit, not only the picture. At
$`\lambda = 0.0609`$ the slope loading falls from 0.992 at three months
to 0.749 at ten years, which is close to the level loading of 1 at every
maturity. The two factors are then nearly collinear and not separately
identified, and $`\beta_2`$ no longer behaves like a slope. Fitting the
same panel both ways and comparing $`\beta_2`$ with the quoted spread:

| $`\lambda`$ | hump at | cor($`\beta_2`$, 3m-10y) | same sign | median fit |
|------------:|--------:|-------------------------:|----------:|-----------:|
|      0.0609 |  29.4 y |                   -0.928 |        5% |    23.0 bp |
|       0.448 |   4.0 y |                   +0.941 |       90% |     5.3 bp |

At the month-denominated value the slope factor is negatively correlated
with the spread and carries the same sign on 5 percent of dates. `dnsr`
therefore does not use 0.0609 as its default.

Two qualifications on
[`ns_lambda()`](https://ambientpointcorp.github.io/dnsr/reference/ns_lambda.md).
The curvature peak is flat: at half or double the nominal peak maturity
the loading is still about 85 percent of its maximum, so $`\lambda`$ is
not tightly determined by this choice and the precision of the returned
number overstates what it fixes. And $`\lambda`$ also sets the decay of
the slope loading, which is often the larger effect: a hump at 2.5 years
leaves the slope loading at 0.53 by two years, a hump at 10 years leaves
it at 0.84. Estimate $`\lambda`$ where that is possible.

## Dynamic Extension

The DNS model lets the three factors vary over time, each as an
Ornstein-Uhlenbeck process:

``` math
d\beta_i(t) = \kappa_i(\theta_i - \beta_i(t))dt + \sigma_i dW_i(t)
```

with $`\kappa_i`$ the speed of mean reversion, $`\theta_i`$ the long-run
mean, $`\sigma_i`$ the volatility and $`W_i(t)`$ a Wiener process.
Fitting the curve is therefore two steps: a cross-sectional fit on each
date, then a time-series estimate of how the resulting factors move.

## Fitting a Single Curve

[`fit_ns_model()`](https://ambientpointcorp.github.io/dnsr/reference/fit_ns_model.md)
fits one cross-section. The default `"profile"` method solves the three
betas in closed form at each $`\lambda`$ and searches only over
$`\lambda`$, so it needs no starting values for them:

``` r

maturities <- as.numeric(setdiff(names(treasury_yields), "Date"))
latest <- as.numeric(treasury_yields[nrow(treasury_yields), -1])

curve <- fit_ns_model(maturities, latest)
curve
#> Nelson-Siegel curve fit
#>   beta1 (level)     0.04927
#>   beta2 (slope)     -0.00515
#>   beta3 (curvature) -0.01551
#>   lambda            0.66877
#>   rmse              4.43 bp over 11 maturities
```

The fit reports whether any parameter came to rest on one of its bounds.
A clamped fit can have a low RMSE while the factor decomposition is not
identified, and a bounded optimiser reports convergence at a bound, so
`success` alone does not distinguish the two:

``` r

curve$at_bound
#>  beta1  beta2  beta3 lambda 
#>  FALSE  FALSE  FALSE  FALSE
```

Nothing is clamped here. Imposing the ceiling of 0.5 used before version
0.6.0 clamps this curve, which fits a faster decay than that allows:

``` r

pinched <- fit_ns_model(maturities, latest,
                        upper_bounds = c(0.2, 0.5, 0.5, 0.5))
#> Warning in fit_ns_model(maturities, latest, upper_bounds = c(0.2, 0.5, 0.5, :
#> Nelson-Siegel fit rests on a parameter bound (lambda); the factor decomposition
#> is not identified. Widen lower_bounds/upper_bounds or check the yield scale.
c(lambda = pinched$parameters[["lambda"]], rmse_bp = pinched$rmse * 10000)
#>   lambda  rmse_bp 
#> 0.500000 4.948436
pinched$at_bound
#>  beta1  beta2  beta3 lambda 
#>  FALSE  FALSE  FALSE   TRUE
```

The clamped fit returns $`\lambda`$ exactly on 0.5 with a higher error,
and flags it through `at_bound`.

## Fitting the Panel

[`fit_dns_model()`](https://ambientpointcorp.github.io/dnsr/reference/fit_dns_model.md)
fits every date. The window is selected by date, so it does not matter
which row order the panel arrives in.

A panel raises a choice that a single curve does not. If $`\lambda`$ is
estimated freely it varies from date to date, on this data from 0.17 to
2.00, the ceiling. $`\beta_2`$ then multiplies a differently shaped
loading at each date, so the betas are not a comparable time series,
which is what the next step requires.

Standard DNS practice, following Diebold and Li, therefore fixes
$`\lambda`$ at a single value for the whole panel. The cost is
cross-sectional fit, about 6.4 basis points median here against 4.5 with
$`\lambda`$ free; the gain is factors that are comparable across dates
and can be modelled over time.

`ns_lambda(4)` places the curvature hump at four years, close to the
panel optimum, and keeps the loadings distinct:

``` r

recent <- fit_dns_model(treasury_yields, n_days = 250, show_progress = FALSE,
                        lambda_constraint = "fixed",
                        initial_lambda = ns_lambda(4))
recent
#> Dynamic Nelson-Siegel fit
#>   250 curves, 2024-01-02 to 2024-12-31
#>   maturities         0.0833  0.2500  0.5000  1.0000  2.0000  3.0000  5.0000  7.0000 10.0000 20.0000 30.0000
#>   rmse              5.28 bp median, 12.08 bp worst
#>   lambda            0.4483 median
#>   success rate      100.0%
```

[`summary()`](https://rdrr.io/r/base/summary.html) gives the
distribution of each factor and of the fit error:

``` r

summary(recent)
#> Dynamic Nelson-Siegel fit
#>   250 curves, 2024-01-02 to 2024-12-31, 100.0% fitted
#> 
#> Factors:
#>             min     mean      max      sd
#> beta1   0.04314  0.04649  0.05018 0.00135
#> beta2  -0.00663  0.00680  0.01352 0.00484
#> beta3  -0.04877 -0.02768 -0.00824 0.01029
#> lambda  0.44832  0.44832  0.44832 0.00000
#> 
#> Fit error (basis points):
#>    0%   50%   95%  100% 
#>  3.43  5.28 10.24 12.08
```

Fitted curves against the observations:

``` r

plot_yield_curves(treasury_yields, recent, max_plots = 4)
```

![](dns-model-introduction_files/figure-html/plot_fits-1.png)

And the factors over time:

``` r

factor_plots <- plot_factor_time_series(recent)
factor_plots$factors
```

![](dns-model-introduction_files/figure-html/plot_factors-1.png)

``` r

factor_plots$rmse
```

![](dns-model-introduction_files/figure-html/plot_rmse-1.png)

### Checking the factor interpretation

The three factors are described as level, slope and curvature. Each
claim is testable. The slope factor should track the quoted spread:

``` r

aligned <- treasury_yields[match(recent$parameters$Date, treasury_yields$Date), ]
observed_spread <- aligned[["0.25"]] - aligned[["10"]]

c(
  correlation = cor(recent$parameters$beta2, observed_spread),
  same_sign = mean(sign(recent$parameters$beta2) == sign(observed_spread))
)
#> correlation   same_sign 
#>   0.9410156   0.9040000
```

The correlation is high and the sign follows from the loadings: the
slope loading is 1 at the short end and decays to zero, so
$`y(\text{short}) \approx \beta_1 + \beta_2`$ and
$`y(\text{long}) \approx \beta_1`$, making $`\beta_2`$ the
short-minus-long spread. A positive $`\beta_2`$ corresponds to an
inverted curve.

The level factor tracks the policy cycle:

``` r

range(recent$parameters$beta1)
#> [1] 0.04314019 0.05018226
```

### Residuals by maturity

The median fit error summarises the residuals; their distribution across
maturities is more informative. Three factors constrain the shapes the
curve can take:

``` r

fitted_rows <- recent$fitted_yields
residual <- as.matrix(aligned[, -1]) - fitted_rows

data.frame(
  maturity = as.numeric(colnames(fitted_rows)),
  mean_bp = round(colMeans(residual, na.rm = TRUE) * 10000, 1),
  rmse_bp = round(sqrt(colMeans(residual^2, na.rm = TRUE)) * 10000, 1)
)
#>        maturity mean_bp rmse_bp
#> 0.0833   0.0833     0.1     7.9
#> 0.25     0.2500     3.3     5.0
#> 0.5      0.5000     1.1     6.1
#> 1        1.0000    -4.7     8.5
#> 2        2.0000    -2.7     5.5
#> 3        3.0000    -0.5     1.9
#> 5        5.0000     1.5     4.3
#> 7        7.0000     3.1     5.0
#> 10      10.0000    -1.2     3.2
#> 20      20.0000     8.6     8.7
#> 30      30.0000    -8.7     9.4
```

The errors are not uniform across maturities. The middle of the curve,
where the curvature loading is largest, fits to a few basis points. Both
ends are three to five times worse, and the long end shows a systematic
pattern: twenty years is fitted low and thirty high, which is what a
single exponential decay produces on a curve that flattens faster than
it can represent. The one-month point has a different cause: it is not
well described by the same curve as the rest, and no three-factor shape
fits both the bill market and the thirty-year point.

Which maturities matter depends on the application. For discounting long
cashflows, the median RMSE is not the relevant summary.

## Estimating Factor Dynamics

[`estimate_dns_dynamics()`](https://ambientpointcorp.github.io/dnsr/reference/estimate_dns_dynamics.md)
regresses each factor’s change on its lagged level. By default it
samples the factor series at month-ends and estimates at
$`\Delta t = 1/12`$. Daily cross-sectional refits carry serially
uncorrelated measurement noise, and an OU regression on daily changes
attributes that noise to mean reversion, inflating $`\kappa`$ by about
an order of magnitude.

Since the estimation uses month-ends, fitting only those is sufficient
and faster. `warm_start = FALSE` makes each cross-section independent of
the other dates in the panel, so these fits match what fitting every day
would produce:

``` r

month_ends <- treasury_yields[
  !duplicated(format(treasury_yields$Date, "%Y-%m"), fromLast = TRUE),
]
monthly_fit <- fit_dns_model(month_ends, n_days = nrow(month_ends),
                             show_progress = FALSE, warm_start = FALSE,
                             lambda_constraint = "fixed",
                             initial_lambda = ns_lambda(4))

dynamics <- estimate_dns_dynamics(monthly_fit$parameters)
dynamics
#> DNS factor dynamics
#>   estimated on     2019-02-28 to 2024-12-31
#>   sampling         monthly (dt = 1/12), 71 observations
#>   structure        diagonal
#> 
#>                 kappa se_kappa   theta  sigma half_life_y
#> level (b1)     0.0604   0.3378  0.0862 0.0088     11.4681
#> slope (b2)     0.3534   0.3484 -0.0058 0.0114      1.9615
#> curvature (b3) 1.5240   0.7039 -0.0150 0.0306      0.4548
#> 
#>   stable           TRUE (eigenvalues 1.524 0.353 0.060)
```

Mean reversion is weakly determined in a sample this short. The standard
errors and the half-life interval show the range involved:

``` r

dynamics$beta1$se
#>        kappa        theta        sigma 
#> 0.3378346439 0.3155525922 0.0007481728
dynamics$beta1$t_value
#>     kappa     theta 
#> 0.1789075 0.2731285
dynamics$beta1$half_life_ci
#>     lower     upper 
#> 0.9592606       Inf
```

$`\kappa`$ also carries an upward small-sample bias of order $`5/T`$ in
a sample of $`T`$ years. The standard errors above measure sampling
noise, not that bias.

Factor innovations are correlated, and the simulation uses the
correlation:

``` r

round(dynamics$correlation, 3)
#>        beta1  beta2  beta3
#> beta1  1.000 -0.548  0.132
#> beta2 -0.548  1.000 -0.302
#> beta3  0.132 -0.302  1.000
```

### A full mean-reversion matrix

By default the three factors revert independently. `structure = "full"`
estimates the Diebold-Li VAR(1) instead, in which a shock to one factor
affects the others:

``` r

full <- estimate_dns_dynamics(monthly_fit$parameters, structure = "full")
round(full$kappa_matrix, 3)
#>        beta1  beta2  beta3
#> beta1  0.436 -0.494 -0.391
#> beta2 -1.190  0.866 -0.843
#> beta3  0.921  0.635  1.675
```

The eigenvalues of that matrix determine whether the system is stable. A
diagonal reconstruction gives a different answer when the cross terms
are non-negligible:

``` r

round(Re(full$eigenvalues), 3)
#> [1] 1.176 1.176 0.626
full$system_stability
#> [1] TRUE
```

## Diagnostics

The package reports fit diagnostics across several objects:

``` r

# 1. Did every cross-section fit at all?
recent$success_rate
#> [1] 100
nrow(recent$failures)
#> [1] 0

# 2. Did any fit come to rest on a bound? Those factor decompositions are not
#    identified, however good their RMSE looks.
colSums(recent$at_bound)
#>  beta1  beta2  beta3 lambda 
#>      0      0      0      0

# 3. How bad is the worst curve, not just the median?
round(quantile(recent$parameters$rmse * 10000, c(0.5, 0.95, 1)), 2)
#>   50%   95%  100% 
#>  5.28 10.24 12.08

# 4. Is the estimated factor system actually mean reverting?
dynamics$system_stability
#> [1] TRUE
round(Re(dynamics$eigenvalues), 3)
#> [1] 1.524 0.353 0.060

# 5. How well identified is the reversion speed?
dynamics$beta1$half_life_ci
#>     lower     upper 
#> 0.9592606       Inf
```

What to look for: a high success rate with no bound hits, an acceptable
worst-case RMSE, a stable system, and a finite half-life interval. If
the half-life interval runs to infinity, as it does on this six-year
sample, $`\kappa`$ is not determined by the data, and a simulation
horizon longer than the sample extrapolates beyond what was estimated.

## Simulating Future Yield Curves

The factor system is a linear SDE whose diffusion does not depend on the
state, so its transition law over any step is a known Gaussian.
`engine = "exact"` samples that law directly, which carries no
discretisation error at any step size and needs no additional package.

``` r

latest_row <- monthly_fit$parameters[which.max(monthly_fit$parameters$Date), ]
initial_factors <- as.numeric(latest_row[c("beta1", "beta2", "beta3")])

sim_params <- dns_dynamics_to_sim_params(dynamics, lambda = latest_row$lambda)

simulation <- simulate_dns_yield_curves(
  params = sim_params,
  initial_factors = initial_factors,
  maturities = maturities,
  time_horizon = 5,
  n_steps = 60,
  n_paths = 500,
  engine = "exact",
  seed = 1
)
simulation
#> DNS simulation
#>   500 paths, 60 steps over 5 years
#>   engine           exact
#>   terminal factors 0.0597  -0.0065  -0.0145
#>     (sd)           0.0177  0.0135  0.0173
#>   maturities        0.0833  0.2500  0.5000  1.0000  2.0000  3.0000  5.0000  7.0000 10.0000 20.0000 30.0000
#>   terminal yields  5.30  5.27  5.24  5.19  5.18  5.21  5.33  5.44  5.57  5.77  5.83 (% mean)
```

Simulated curves at a few horizons. The facets follow the calendar
rather than sorting lexically, so “Year 5” comes after “Year 3”:

``` r

yield_plots <- plot_simulated_yields(
  simulation,
  time_points = c(0, 1, 3, 5),
  path_subset = 1:5
)
yield_plots$yield_curves
```

![](dns-model-introduction_files/figure-html/plot_simulated_yields-1.png)

The factor paths, as a fan of quantiles across all 500 of them:

``` r

plot_simulated_factors(simulation)$level
```

![](dns-model-introduction_files/figure-html/plot_simulated_factors-1.png)

The YUIMA engine remains available and is the default, for the
estimation functions it provides. See
[`vignette("dnsr-yuima-integration", package = "dnsr")`](https://ambientpointcorp.github.io/dnsr/articles/dnsr-yuima-integration.md).

## Pricing off a Fitted Curve

The DNS curve gives continuously compounded zero rates, so the discount
factor at maturity $`m`$ is $`P(m) = \exp(-y(m)\,m)`$:

``` r

factors <- matrix(curve$parameters[1:3], nrow = 1)
lambda <- curve$parameters[["lambda"]]

coupon_maturities <- c(0.5, 1, 2, 5, 10, 30)
data.frame(
  maturity = coupon_maturities,
  zero_pct = round(calculate_dns_yields(factors, coupon_maturities, lambda)[1, ] * 100, 3),
  discount = round(dns_discount_factors(factors, coupon_maturities, lambda)[1, ], 4),
  price_100 = round(dns_zero_prices(factors, coupon_maturities, lambda)[1, ], 2),
  par_pct = round(dns_par_yields(factors, coupon_maturities, lambda)[1, ] * 100, 3)
)
#>   maturity zero_pct discount price_100 par_pct
#> 1      0.5    4.282   0.9788     97.88   4.328
#> 2      1.0    4.216   0.9587     95.87   4.261
#> 3      2.0    4.196   0.9195     91.95   4.240
#> 4      5.0    4.386   0.8031     80.31   4.423
#> 5     10.0    4.621   0.6300     63.00   4.639
#> 6     30.0    4.824   0.2352     23.52   4.812
```

Par yields are the coupon rate that prices a bullet bond at par. They
differ from zero yields for two reasons with opposite signs: a
compounding term, since the curve is continuously compounded while the
par rate pays semi-annually, and a shape term, since coupons are
discounted at the earlier rates rather than the terminal one. The
compounding term dominates at the short end, the shape term further out.

Forward rates come off the same factors:

``` r

rbind(
  spot = round(calculate_dns_yields(factors, coupon_maturities, lambda)[1, ] * 100, 3),
  forward = round(calculate_dns_forwards(factors, coupon_maturities, lambda)[1, ] * 100, 3)
)
#>          [,1]  [,2]  [,3]  [,4]  [,5]  [,6]
#> spot    4.282 4.216 4.196 4.386 4.621 4.824
#> forward 4.188 4.132 4.248 4.726 4.914 4.927
```

## Looking at the Panel

A heatmap of the whole panel shows the 2020 collapse to near-zero short
rates and the 2022-2023 inversion:

``` r

plot_yield_heatmap(treasury_yields)
```

![](dns-model-introduction_files/figure-html/heatmap-1.png)

Principal components of the daily changes. Decomposing the covariance
matrix rather than standardising each maturity first is the usual
convention for yield curves, and produces the level, slope and curvature
shapes:

``` r

pca <- plot_yield_pca(treasury_yields, scale = FALSE)
round(pca$var_explained[1:3], 4)
#> [1] 0.7160 0.1087 0.0991
pca$plots$weights
```

![](dns-model-introduction_files/figure-html/pca-1.png)

The first three components are the empirical counterpart of the model’s
three factors.

## Comparing Term Structure Models

[`compare_term_structure_models()`](https://ambientpointcorp.github.io/dnsr/reference/compare_term_structure_models.md)
fits several curves to one cross-section and reports the error of each,
in basis points:

``` r

comparison <- compare_term_structure_models(
  treasury_yields,
  models = c("NS", "NSS", "Spline")
)
comparison$metrics
#>                    Model     RMSE      MAE  MaxError
#> 1          Nelson-Siegel 4.431432 3.887569  8.700489
#> 2 Nelson-Siegel-Svensson 3.542643 2.502529  8.621915
#> 3           Cubic Spline 5.561295 3.995767 11.985578
```

``` r

comparison$plot
```

![](dns-model-introduction_files/figure-html/comparison_plot-1.png)

Nelson-Siegel-Svensson adds a second curvature term with its own decay.
It nests Nelson-Siegel and so cannot fit worse in sample, at the cost of
two further parameters. A smoothing spline usually has a lower in-sample
error than either and imposes no shape, so its parameters do not carry
an interpretation.

## One-Step Estimation

Everything above estimates in two steps: each cross-section is fitted,
and the resulting factors are then treated as data for a time-series
model. The factors are estimates carrying sampling error, and using them
as observed biases $`\kappa`$ upward.

[`estimate_dns_state_space()`](https://ambientpointcorp.github.io/dnsr/reference/estimate_dns_state_space.md)
estimates the same model in one step, as the state-space system of
Diebold, Rudebusch and Aruoba, by maximum likelihood with a Kalman
filter. The factors become an unobserved state, their measurement
uncertainty enters the likelihood, $`\lambda`$ is estimated jointly with
the dynamics, and missing maturities are handled by shrinking the
measurement equation for that date.

``` r

one_step <- estimate_dns_state_space(treasury_yields,
                                     lambda = ns_lambda(3.5),
                                     measurement = "common")
one_step
#> DNS one-step state-space fit (Kalman filter, maximum likelihood)
#>   estimated on     2019-01-31 to 2024-12-31
#>   sampling         monthly (dt = 1/11.8), 72 dates, 11 maturities
#>   lambda           0.5124 fixed (hump at 3.50 years)
#> 
#>                 kappa   theta  sigma half_life_y
#> level (b1)     0.1214  0.0365 0.0078      5.7109
#> slope (b2)     0.4663 -0.0034 0.0109      1.4865
#> curvature (b3) 0.4533 -0.0160 0.0259      1.5291
#> 
#>   measurement sd   13.37 bp median across maturities
#>   fit              11.92 bp rmse
#>   logLik           3863.72 on 13 parameters (AIC -7701.4, BIC -7640.7)
#>   stable           TRUE
#>   $two_step holds the two-step estimates used as starting values
```

The smoothed factors have the same shape as the two-step `parameters`
frame, and [`summary()`](https://rdrr.io/r/base/summary.html) reports
both sets side by side:

``` r

head(one_step$factors)
#>         Date      beta1        beta2       beta3
#> 1 2019-01-31 0.03007925 -0.005008562 -0.01156606
#> 2 2019-02-28 0.03077127 -0.005530049 -0.01055434
#> 3 2019-03-29 0.02872613 -0.003504591 -0.01647100
#> 4 2019-04-30 0.02929922 -0.003959725 -0.01651064
#> 5 2019-05-31 0.02671109 -0.001949996 -0.02249331
#> 6 2019-06-28 0.02555954 -0.002923954 -0.02310094
```

The costs are run time, needing the whole panel at once, and a harder
numerical problem. There is also a limit on what either estimator
determines: $`\lambda`$, the shock volatilities and the measurement
error are recovered precisely, but $`\kappa`$ and $`\theta`$ are not
well determined at realistic sample lengths.
[`?estimate_dns_state_space`](https://ambientpointcorp.github.io/dnsr/reference/estimate_dns_state_space.md)
gives the measured errors from simulated data.

The factor system here is diagonal. For a full mean-reversion matrix,
use `estimate_dns_dynamics(structure = "full")` on two-step betas.

## What This Package Does Not Do

**It carries no arbitrage-free restriction.** Nelson-Siegel is a
flexible curve rather than a term structure model derived from a
no-arbitrage condition, so a fitted or simulated curve can imply
arbitrage. The arbitrage-free variant, AFNS, adds a yield-adjustment
term that ties the dynamics and the cross-section together. The
distinction matters for pricing derivatives off these curves, less so
for discounting cashflows.

## Where to Go Next

``` r

# Live Treasury data, caching and the API key
vignette("fred-api-setup", package = "dnsr")

# The YUIMA engine, and quasi-maximum-likelihood estimation
vignette("dnsr-yuima-integration", package = "dnsr")
```

[`get_calibrated_dns_params()`](https://ambientpointcorp.github.io/dnsr/reference/get_calibrated_dns_params.md)
wraps the two-step path above into one call, from downloading data to
returning simulation-ready parameters, for when the intermediate objects
are not needed.
