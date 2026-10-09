# Choose the Nelson-Siegel Decay from the Curvature Peak

Converts between the decay parameter lambda and the maturity at which
the curvature loading peaks. `ns_lambda()` goes from a maturity to a
lambda; `ns_curvature_peak()` goes back.

## Usage

``` r
ns_lambda(peak_maturity)

ns_curvature_peak(lambda)
```

## Arguments

- peak_maturity:

  Numeric. The maturity at which the curvature loading should peak, in
  whatever units the maturities are measured in.

- lambda:

  Numeric. A decay parameter, per unit of maturity.

## Value

A numeric vector of the same length as the input: `ns_lambda()` returns
decays, `ns_curvature_peak()` maturities.

## Details

The curvature loading \$\$L_3(m) = \frac{1-e^{-\lambda m}}{\lambda m} -
e^{-\lambda m}\$\$ is maximised where \\\lambda m = 1.793282\\, the root
of \\\partial L_3 / \partial m = 0\\. So one defensible way to pick
lambda is to decide where the hump belongs and divide.

The point of doing it this way is units. Lambda is per unit of maturity,
so its numerical value depends on whether maturities are years or
months, and nothing about a lambda on its own says which. Stating a
maturity instead settles it: `ns_lambda(30)` on a curve measured in
months and `ns_lambda(2.5)` on the same curve measured in years describe
the same shape.

That also explains the package default of 0.0609, which is the value
Diebold and Li published for maturities in **months**, where it places
the hump at 29.4 months. The functions in this package take maturities
in **years**, so the same number places it at 29.4 *years*;
`ns_curvature_peak()` reports that directly.

Two cautions. The peak is flat – at half or double the nominal peak
maturity the loading is still about 85 percent of its maximum – so
lambda is not tightly determined by this choice, and the apparent
precision of the result is misleading. And lambda sets the slope
loading's decay at the same time, which is often the more consequential
effect: a hump at 2.5 years leaves the slope loading at 0.53 by two
years, while a hump at 10 years leaves it at 0.84. Where lambda can be
estimated rather than chosen, for instance with
`fit_ns_model(method = "profile")`, estimating it is usually better.
Fixing it is still the right move for a panel, where the factors have to
be comparable across dates.

## See also

[`plot_dns_loadings`](https://ambientpointcorp.github.io/dnsr/reference/plot_dns_loadings.md)
to see the loadings a lambda produces, and
[`fit_ns_model`](https://ambientpointcorp.github.io/dnsr/reference/fit_ns_model.md)
to estimate one instead.

## Examples

``` r
# A hump at two and a half years, for maturities measured in years
ns_lambda(2.5)
#> [1] 0.7173128

# Diebold and Li aimed for a hump at 30 months, measuring maturities in
# months. That criterion gives 0.0598; the value they published, 0.0609,
# peaks a little short of it.
ns_lambda(30)
#> [1] 0.05977607
ns_curvature_peak(0.0609)
#> [1] 29.44634

# What the package default implies when maturities are in years
ns_curvature_peak(0.0609)
#> [1] 29.44634

# They are inverses
ns_curvature_peak(ns_lambda(4))
#> [1] 4
```
