# Convert DNS Dynamics to Simulation Parameters

Converts the estimated DNS dynamics to parameters suitable for
simulation.

## Usage

``` r
dns_dynamics_to_sim_params(dynamics, include_correlation = TRUE, lambda = NULL)
```

## Arguments

- dynamics:

  List. The dynamics object from estimate_dns_dynamics().

- include_correlation:

  Logical. If TRUE, include the correlation structure.

- lambda:

  Numeric. The Nelson-Siegel decay parameter to carry into the
  simulation parameters. If NULL (the default) it is taken from
  `dynamics$lambda` when present, and otherwise falls back to 0.0609.
  [`estimate_dns_dynamics()`](https://ambientpointcorp.github.io/dnsr/reference/estimate_dns_dynamics.md)
  does not estimate lambda, so pass the value fitted by
  [`fit_dns_model`](https://ambientpointcorp.github.io/dnsr/reference/fit_dns_model.md)
  to avoid simulating off the default.

## Value

A list of parameters for DNS simulation.

## Details

`kappa`, `theta` and `sigma` are the per-factor values, so a dynamics
object estimated with `structure = "full"` contributes only the diagonal
of its mean-reversion matrix here:
[`simulate_dns_factors()`](https://ambientpointcorp.github.io/dnsr/reference/simulate_dns_factors.md)
reverts each factor to its own long-run mean. The full matrix is carried
through as `kappa_matrix` for inspection.

## Examples

``` r
data(treasury_yields)

# Month-ends only: the dynamics are estimated on them anyway
month_ends <- treasury_yields[
  !duplicated(format(treasury_yields$Date, "%Y-%m"), fromLast = TRUE),
]
dns_result <- fit_dns_model(month_ends, n_days = nrow(month_ends),
                            show_progress = FALSE, warm_start = FALSE)
#> Warning: Parameters rest on a bound in some fits: lambda on 13 (of 72 dates). The factor decomposition is not identified on those dates; see the at_bound element of the result, widen lower_bounds/upper_bounds, or check the yield scale.
dynamics <- estimate_dns_dynamics(dns_result$parameters)

# Pass the fitted lambda through, rather than leaving the 0.0609 default
latest <- dns_result$parameters[which.max(dns_result$parameters$Date), ]
sim_params <- dns_dynamics_to_sim_params(dynamics, lambda = latest$lambda)
str(sim_params)
#> List of 7
#>  $ kappa       : num [1:3] 2.78 2.63 2.12
#>  $ theta       : num [1:3] 0.03456 -0.00831 -0.00431
#>  $ sigma       : num [1:3] 0.0304 0.033 0.0579
#>  $ names       : chr [1:3] "Level" "Slope" "Curvature"
#>  $ lambda      : num 0.0614
#>  $ correlation : num [1:3, 1:3] 1 -0.791 -0.549 -0.791 1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:3] "beta1" "beta2" "beta3"
#>   .. ..$ : chr [1:3] "beta1" "beta2" "beta3"
#>  $ kappa_matrix: num [1:3, 1:3] 2.78 0 0 0 2.63 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:3] "beta1" "beta2" "beta3"
#>   .. ..$ : chr [1:3] "beta1" "beta2" "beta3"
```
