# Get Treasury Yield Curve Data from FRED

Downloads U.S. Treasury yield curve data from FRED API for multiple
maturities and processes it into a consistent format suitable for DNS
model calibration.

## Usage

``` r
get_treasury_yield_data(
  api_key = NULL,
  maturities = c("1MO", "3MO", "6MO", "1", "2", "3", "5", "7", "10", "20", "30"),
  start_date = NULL,
  end_date = NULL,
  convert_to_decimal = TRUE,
  verbose = FALSE,
  fill_missing = TRUE,
  max_gap = Inf,
  date_order = c("descending", "ascending"),
  retries = 3L,
  cache_dir = NULL
)
```

## Arguments

- api_key:

  Character. FRED API key. If NULL, will look for fred_api_key
  environment variable.

- maturities:

  Character vector. Maturities to download. Options are "1MO", "2MO",
  "3MO", "6MO", "1", "2", "3", "5", "7", "10", "20", "30".

- start_date:

  Date or character in "YYYY-MM-DD" format. The earliest date to
  include.

- end_date:

  Date or character in "YYYY-MM-DD" format. The latest date to include.

- convert_to_decimal:

  Logical. If TRUE, convert yields from percentages to decimal format.

- verbose:

  Logical. If TRUE, print processing information.

- fill_missing:

  Logical, or one of `"forward"` / `"none"`. If TRUE (the default) gaps
  are filled by carrying the last observation forward in calendar time.

- max_gap:

  Numeric. Maximum number of calendar days a value may be carried across
  when `fill_missing` is in effect. Gaps wider than this are left as
  `NA`. Defaults to `Inf`, which fills every gap.

- date_order:

  Character. Row order of the returned frame: `"descending"` (the
  default, most recent first) or `"ascending"`.

- retries:

  Integer. HTTP attempts per series, passed to
  [`fetch_fred_series`](https://ambientpointcorp.github.io/dnsr/reference/fetch_fred_series.md).

- cache_dir:

  Character. Directory in which to cache responses, passed to
  [`fetch_fred_series`](https://ambientpointcorp.github.io/dnsr/reference/fetch_fred_series.md).
  `NULL` disables caching.

## Value

A data.frame with processed yield curve data where:

- The first column is 'Date' containing the observation dates.

- Subsequent columns contain yield values for different maturities.

- Column names represent maturities in years.

## Details

This function fetches data from FRED for the following series:

- DGS1MO (1-Month)

- DGS3MO (3-Month)

- DGS6MO (6-Month)

- DGS1 (1-Year)

- DGS2 (2-Year)

- DGS3 (3-Year)

- DGS5 (5-Year)

- DGS7 (7-Year)

- DGS10 (10-Year)

- DGS20 (20-Year)

- DGS30 (30-Year)

Requires a FRED API key which can be obtained free at
https://fred.stlouisfed.org/docs/api/api_key.html

Gap filling runs on an ascending date index, so a missing yield takes
the value of the most recent *earlier* observation. This matters because
the frame is returned newest-first: filling in row order would carry
values backwards in time and introduce look-ahead bias into anything
calibrated on the result. Long series have genuine multi-year holes
(DGS20 before 1993, DGS30 between 2002 and 2006); set `max_gap` to stop
those being bridged.

## Examples

``` r
if (FALSE) { # \dontrun{
# Set your API key first
Sys.setenv(fred_api_key = "your_api_key_here")
yield_data <- get_treasury_yield_data(start_date = "2020-01-01")
head(yield_data)
} # }
```
