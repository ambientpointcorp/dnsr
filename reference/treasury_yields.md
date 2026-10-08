# U.S. Treasury Yield Curves, 2019-2024

Daily U.S. Treasury constant maturity yields across eleven maturities,
in the shape the rest of the package expects: a `Date` column followed
by one column per maturity, named by maturity in years.

Bundled so that the examples, vignettes and tests have a real yield
panel to work on without a FRED API key.

## Usage

``` r
treasury_yields
```

## Format

A data frame with 1,501 rows and 12 columns:

- Date:

  Observation date, ascending.

- 0.0833, 0.25, 0.5:

  1-month, 3-month and 6-month yields.

- 1, 2, 3, 5, 7, 10, 20, 30:

  1- to 30-year yields.

Yields are in decimal form, so 0.0425 is 4.25 percent. There are no
missing values: market holidays, which FRED reports as rows with no
observation at any maturity, are excluded rather than carried as `NA`.

## Source

Federal Reserve Economic Data (FRED), series `DGS1MO`, `DGS3MO`,
`DGS6MO`, `DGS1`, `DGS2`, `DGS3`, `DGS5`, `DGS7`, `DGS10`, `DGS20` and
`DGS30`. These are the Treasury constant maturity series of the Federal
Reserve's H.15 release, retrieved 2026-10-08 and covering 2 January 2019
to 31 December 2024. As U.S. government works they are in the public
domain.

The script that builds this dataset is `data-raw/treasury_yields.R` in
the package sources; it calls
[`get_treasury_yield_data`](https://ambientpointcorp.github.io/dnsr/reference/get_treasury_yield_data.md).

## Details

Six years of daily observations spans 72 month-ends, which is enough for
[`estimate_dns_dynamics`](https://ambientpointcorp.github.io/dnsr/reference/estimate_dns_dynamics.md)
to run in its default monthly sampling mode. The window covers the 2020
collapse to near-zero short rates, the 2022-2023 inversion, and the
start of the subsequent normalisation, so it contains level, slope and
curvature variation rather than a single regime.

## Examples

``` r
data(treasury_yields)

# The most recent curve
tail(treasury_yields, 1)
#>            Date 0.0833   0.25    0.5      1      2      3      5      7     10
#> 1501 2024-12-31  0.044 0.0437 0.0424 0.0416 0.0425 0.0427 0.0438 0.0448 0.0458
#>          20     30
#> 1501 0.0486 0.0478

# Dates on which the 3-month yield exceeded the 10-year
inverted <- treasury_yields$`0.25` > treasury_yields$`10`
range(treasury_yields$Date[inverted])
#> [1] "2019-03-22" "2024-12-12"
```
