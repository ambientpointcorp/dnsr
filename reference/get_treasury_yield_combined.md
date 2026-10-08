# Get Historical Treasury Yield Data from FRED

Downloads U.S. Treasury yield curve data from FRED API for a specified
date range. This is a simplified wrapper around
get_treasury_yield_data() for backward compatibility.

## Usage

``` r
get_treasury_yield_combined(
  api_key = NULL,
  start_date = NULL,
  end_date = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- api_key:

  Character. FRED API key. If NULL, will look for fred_api_key
  environment variable.

- start_date:

  Date or character in "YYYY-MM-DD" format. The earliest date to
  include. If NULL, defaults to "1990-01-01".

- end_date:

  Date or character in "YYYY-MM-DD" format. The latest date to include.
  If NULL, defaults to current date.

- verbose:

  Logical. If TRUE, print processing information.

- ...:

  Further arguments passed to
  [`get_treasury_yield_data`](https://ambientpointcorp.github.io/dnsr/reference/get_treasury_yield_data.md),
  for example `max_gap`, `date_order` or `fill_missing`.

## Value

A data.frame with combined yield curve data, most recent date first.

## Details

This function is equivalent to calling get_treasury_yield_data() with
default maturities. It provides a simpler interface for users who just
want all available maturities.

## Examples

``` r
if (FALSE) { # \dontrun{
# Set your API key first
Sys.setenv(fred_api_key = "your_api_key_here")
yield_data <- get_treasury_yield_combined(start_date = "2018-01-01")
head(yield_data)
} # }
```
