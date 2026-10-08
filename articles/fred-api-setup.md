# FRED API Setup

## FRED API Setup

Before using the `dnsr` package to fetch Treasury yield data, you need
to set up access to the Federal Reserve Economic Data (FRED) API.

### Getting a FRED API Key

1.  **Create a FRED account** (free):
    - Visit <https://fred.stlouisfed.org/>
    - Click “Sign In” and create a free account
2.  **Request an API key**:
    - Once logged in, go to
      <https://fred.stlouisfed.org/docs/api/api_key.html>
    - Click “Request API Key”
    - Fill out a brief form explaining your use case
    - Your API key will be generated instantly
3.  **Set up your API key in R**:

``` r

# Option 1: Set for current R session only
Sys.setenv(fred_api_key = "your_api_key_here")

# Option 2: Add to .Renviron for permanent setup
usethis::edit_r_environ()
# In the file that opens, add this line:
# fred_api_key=your_api_key_here
# Save and restart R
```

### Verifying Your Setup

After setting your API key, verify it works:

``` r

# Try fetching a small amount of data
test_data <- dnsr::fetch_fred_series(
  series_id = "DGS10",
  start_date = "2023-01-01",
  end_date = "2023-01-31"
)

# If successful, you should see data:
head(test_data)
```

### API Rate Limits

FRED API has the following limits: - **120 requests per minute** per IP
address - **No daily limit** for registered users

The `dnsr` package typically makes 10-12 API calls when fetching a full
yield curve (one per maturity), which is well within these limits.

### Troubleshooting

#### “FRED API key not found” error

If you see this error:

    Error: FRED API key not found. Set fred_api_key environment variable or pass api_key parameter.

**Solution**: Make sure you’ve set the environment variable. Check with:

``` r

Sys.getenv("fred_api_key")
```

If it returns an empty string, your key isn’t set. Use
[`Sys.setenv()`](https://rdrr.io/r/base/Sys.setenv.html) to set it.

#### “FRED API request failed” error

If you see this error, possible causes: 1. **Invalid API key**:
Double-check your key is correct 2. **Network issues**: Verify your
internet connection 3. **FRED service down**: Check
<https://fred.stlouisfed.org/> 4. **Invalid series ID**: Make sure the
series exists

#### Rate limit exceeded

If you hit rate limits:

    Error: FRED API request failed with status 429: Too Many Requests

**Solution**: Wait one minute before retrying, or implement caching:

``` r

# Cache data to avoid repeated API calls
cache_file <- "treasury_data_cache.rds"

if (file.exists(cache_file)) {
  # Load from cache if recent (< 24 hours old)
  cache <- readRDS(cache_file)
  if (difftime(Sys.time(), cache$timestamp, units = "hours") < 24) {
    yield_data <- cache$data
  }
} else {
  # Fetch fresh data
  yield_data <- get_treasury_yield_combined(start_date = "2020-01-01")
  saveRDS(list(data = yield_data, timestamp = Sys.time()), cache_file)
}
```

### Security Best Practices

⚠️ **Never commit your API key to version control!**

- Don’t hardcode your API key in scripts
- Use `.Renviron` for permanent storage
- Add `.Renviron` to your `.gitignore` file
- When sharing code, provide instructions for others to set up their own
  keys

------------------------------------------------------------------------

## Loading the Package

Now that FRED API access is set up, load the package:

``` r

library(dnsr)
```

------------------------------------------------------------------------

\[Rest of vignette content continues here…\]

------------------------------------------------------------------------

## Example: Fetching Recent Treasury Data

With your API key configured, fetching Treasury data is simple:

``` r

# Get the last 2 years of data for all maturities
yield_data <- get_treasury_yield_combined(
  start_date = Sys.Date() - 365 * 2,
  verbose = TRUE
)

# View the structure
head(yield_data)
str(yield_data)
```

### Selecting Specific Maturities

You can also fetch only specific maturities:

``` r

# Get only 2, 5, 10, and 30-year yields
yield_data <- get_treasury_yield_data(
  maturities = c("2", "5", "10", "30"),
  start_date = "2020-01-01",
  end_date = "2023-12-31",
  verbose = TRUE
)
```

Available maturity codes: - Short-term: `"1MO"`, `"2MO"`, `"3MO"`,
`"6MO"` - Medium-term: `"1"`, `"2"`, `"3"`, `"5"`, `"7"` - Long-term:
`"10"`, `"20"`, `"30"`

------------------------------------------------------------------------

## Data Format

The fetched data is returned as a data frame with: - **First column**:
`Date` (most recent dates first) - **Remaining columns**: Yield values
for each maturity (in decimal form) - **Column names**: Maturity in
years (e.g., `"0.25"` for 3-month, `"10"` for 10-year)

``` r

# Example data structure:
#   Date         0.25    0.5     1       2       5       10      30
#   2023-12-31  0.0521  0.0534  0.0489  0.0452  0.0425  0.0398  0.0412
#   2023-12-29  0.0523  0.0536  0.0491  0.0454  0.0427  0.0400  0.0414
#   ...
```

This format is designed for direct use with the DNS model fitting
functions in the package.

------------------------------------------------------------------------
