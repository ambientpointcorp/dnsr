# Get Calibrated DNS Parameters

Returns calibrated DNS parameters based on the most recent market data,
suitable for use in simulation or pricing applications.

## Usage

``` r
get_calibrated_dns_params(
  lookback_days = 60,
  lambda = NULL,
  start_date = NULL,
  end_date = NULL,
  verbose = FALSE,
  api_key = NULL
)
```

## Arguments

- lookback_days:

  Integer. Number of days to use for calibration.

- lambda:

  Numeric. Fixed lambda value to use for calibration. See
  [`ns_lambda`](https://ambientpointcorp.github.io/dnsr/reference/ns_lambda.md)
  on its units; the 0.0609 fallback is a month-denominated constant and
  maturities here are in years.

- start_date:

  Date or character. Start date for the data window.

- end_date:

  Date or character. End date for the data window.

- verbose:

  Logical. If TRUE, print progress information.

- api_key:

  Character. FRED API key, passed through to the data fetch. If NULL the
  `fred_api_key` environment variable is used.

## Value

A list containing:

- params: List of DNS parameters for simulation, including the lambda
  actually fitted to the most recent curve

- initial_factors: Vector of current factor values

- lambda: Fitted lambda value

- calibration_summary: Summary statistics about the calibration

## Examples

``` r
if (FALSE) { # \dontrun{
dns_params <- get_calibrated_dns_params(lookback_days = 60)
} # }
```
