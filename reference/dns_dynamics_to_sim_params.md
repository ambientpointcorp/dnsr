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
  `dynamics$lambda` when present, and otherwise falls back to
  `ns_lambda(3.5)`.
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
#> Warning: Parameters rest on a bound in some fits: lambda on 13 (of 72 dates). The factor decomposition is not identified on those dates; see the at_bound element of the result, pass wider lower_bounds/upper_bounds, or check the yield scale.
dynamics <- estimate_dns_dynamics(dns_result$parameters)

# Pass the fitted lambda through, rather than leaving the default
latest <- dns_result$parameters[which.max(dns_result$parameters$Date), ]
sim_params <- dns_dynamics_to_sim_params(dynamics, lambda = latest$lambda)
str(sim_params)
#> List of 7
#>  $ kappa       : num [1:3] 0.351 0.573 4.514
#>  $ theta       : num [1:3] 0.039 -0.00605 -0.01289
#>  $ sigma       : num [1:3] 0.0111 0.0146 0.0578
#>  $ names       : chr [1:3] "Level" "Slope" "Curvature"
#>  $ lambda      : num 0.669
#>  $ correlation : num [1:3, 1:3] 1 -0.594 -0.184 -0.594 1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:3] "beta1" "beta2" "beta3"
#>   .. ..$ : chr [1:3] "beta1" "beta2" "beta3"
#>  $ kappa_matrix: num [1:3, 1:3] 0.351 0 0 0 0.573 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:3] "beta1" "beta2" "beta3"
#>   .. ..$ : chr [1:3] "beta1" "beta2" "beta3"
```
