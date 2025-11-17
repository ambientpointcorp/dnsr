#' @title Treasury Yield Curve Data Utilities
#' @description Functions for downloading and processing U.S. Treasury yield curve data from FRED API.

#' @importFrom httr GET content
#' @importFrom dplyr %>% mutate across where lag filter select all_of bind_rows arrange
#' @importFrom lubridate mdy date
#' @importFrom tidyr pivot_wider replace_na
#' @importFrom stats na.omit

#' @title Download Treasury Yield Data from FRED
#'
#' @description
#' Downloads U.S. Treasury yield data for specific maturities from the Federal Reserve Economic Data (FRED) API.
#'
#' @param series_id Character. FRED series ID for the treasury yield (e.g., "DGS10" for 10-year).
#' @param api_key Character. FRED API key. If NULL, will look for fred_api_key environment variable.
#' @param start_date Date or character in "YYYY-MM-DD" format. The earliest date to include.
#' @param end_date Date or character in "YYYY-MM-DD" format. The latest date to include.
#' @param verbose Logical. If TRUE, print processing information.
#'
#' @return A data.frame with date and yield value columns.
#'
#' @details
#' Requires a FRED API key which can be obtained free at https://fred.stlouisfed.org/docs/api/api_key.html
#' Set the environment variable fred_api_key or pass it directly to the function.
#'
#' @examples
#' \dontrun{
#' # Set your API key first
#' Sys.setenv(fred_api_key = "your_api_key_here")
#' yield_10y <- fetch_fred_series("DGS10", start_date = "2020-01-01")
#' head(yield_10y)
#' }
#'
#' @export
fetch_fred_series <- function(series_id, api_key = NULL, start_date = NULL,
                              end_date = NULL, verbose = FALSE) {
  # Parameter validation
  checkmate::assert_string(series_id, min.chars = 1)
  checkmate::assert_logical(verbose, len = 1)

  # Get API key
  if (is.null(api_key)) {
    api_key <- Sys.getenv("fred_api_key")
    if (api_key == "") {
      stop("FRED API key not found. Set fred_api_key environment variable or pass api_key parameter.")
    }
  }

  # Create informative logging function
  log_info <- function(...) {
    if (verbose) cat(..., "\n")
  }

  # Build API URL
  base_url <- "https://api.stlouisfed.org/fred/series/observations"

  params <- list(
    series_id = series_id,
    api_key = api_key,
    file_type = "json"
  )

  if (!is.null(start_date)) {
    params$observation_start <- as.character(as.Date(start_date))
  }

  if (!is.null(end_date)) {
    params$observation_end <- as.character(as.Date(end_date))
  }

  log_info(sprintf("Fetching FRED series: %s", series_id))

  # Make API request
  response <- tryCatch({
    httr::GET(base_url, query = params)
  }, error = function(e) {
    stop(sprintf("Failed to connect to FRED API: %s", e$message))
  })

  # Check response status
  if (httr::http_error(response)) {
    stop(sprintf("FRED API request failed with status %d: %s",
                 httr::status_code(response),
                 httr::http_status(response)$message))
  }

  # Parse JSON response
  content <- httr::content(response, as = "parsed")

  if (is.null(content$observations) || length(content$observations) == 0) {
    warning(sprintf("No data returned for series %s", series_id))
    return(data.frame(date = as.Date(character(0)), value = numeric(0)))
  }

  # Extract data
  dates <- sapply(content$observations, function(x) x$date)
  values <- sapply(content$observations, function(x) {
    val <- x$value
    if (val == ".") return(NA_real_)
    return(as.numeric(val))
  })

  result <- data.frame(
    date = as.Date(dates),
    value = values,
    stringsAsFactors = FALSE
  )

  log_info(sprintf("Retrieved %d observations for %s", nrow(result), series_id))

  return(result)
}


#' @title Get Treasury Yield Curve Data from FRED
#'
#' @description
#' Downloads U.S. Treasury yield curve data from FRED API for multiple maturities
#' and processes it into a consistent format suitable for DNS model calibration.
#'
#' @param api_key Character. FRED API key. If NULL, will look for fred_api_key environment variable.
#' @param maturities Character vector. Maturities to download. Options include:
#'   "1MO", "3MO", "6MO", "1", "2", "3", "5", "7", "10", "20", "30"
#' @param start_date Date or character in "YYYY-MM-DD" format. The earliest date to include.
#' @param end_date Date or character in "YYYY-MM-DD" format. The latest date to include.
#' @param convert_to_decimal Logical. If TRUE, convert yields from percentages to decimal format.
#' @param verbose Logical. If TRUE, print processing information.
#' @param fill_missing Logical. If TRUE, attempt to fill in missing values using forward fill.
#'
#' @return A data.frame with processed yield curve data where:
#'   - The first column is 'Date' containing the observation dates.
#'   - Subsequent columns contain yield values for different maturities.
#'   - Column names represent maturities in years.
#'
#' @details
#' This function fetches data from FRED for the following series:
#' - DGS1MO (1-Month)
#' - DGS3MO (3-Month)
#' - DGS6MO (6-Month)
#' - DGS1 (1-Year)
#' - DGS2 (2-Year)
#' - DGS3 (3-Year)
#' - DGS5 (5-Year)
#' - DGS7 (7-Year)
#' - DGS10 (10-Year)
#' - DGS20 (20-Year)
#' - DGS30 (30-Year)
#'
#' Requires a FRED API key which can be obtained free at https://fred.stlouisfed.org/docs/api/api_key.html
#'
#' @examples
#' \dontrun{
#' # Set your API key first
#' Sys.setenv(fred_api_key = "your_api_key_here")
#' yield_data <- get_treasury_yield_data(start_date = "2020-01-01")
#' head(yield_data)
#' }
#'
#' @export
get_treasury_yield_data <- function(api_key = NULL,
                                    maturities = c("1MO", "3MO", "6MO", "1", "2", "3",
                                                   "5", "7", "10", "20", "30"),
                                    start_date = NULL,
                                    end_date = NULL,
                                    convert_to_decimal = TRUE,
                                    verbose = FALSE,
                                    fill_missing = TRUE) {
  # Parameter validation
  checkmate::assert_logical(convert_to_decimal, len = 1)
  checkmate::assert_logical(verbose, len = 1)
  checkmate::assert_logical(fill_missing, len = 1)
  checkmate::assert_character(maturities, min.len = 1)

  # Create informative logging function
  log_info <- function(...) {
    if (verbose) cat(..., "\n")
  }

  # Map maturity codes to FRED series IDs
  maturity_map <- c(
    "1MO" = "DGS1MO",
    "2MO" = "DGS2MO",
    "3MO" = "DGS3MO",
    "6MO" = "DGS6MO",
    "1" = "DGS1",
    "2" = "DGS2",
    "3" = "DGS3",
    "5" = "DGS5",
    "7" = "DGS7",
    "10" = "DGS10",
    "20" = "DGS20",
    "30" = "DGS30"
  )

  # Map maturity codes to column names (years)
  maturity_years <- c(
    "1MO" = "0.0833",
    "2MO" = "0.1667",
    "3MO" = "0.25",
    "6MO" = "0.5",
    "1" = "1",
    "2" = "2",
    "3" = "3",
    "5" = "5",
    "7" = "7",
    "10" = "10",
    "20" = "20",
    "30" = "30"
  )

  # Validate maturities
  invalid_maturities <- setdiff(maturities, names(maturity_map))
  if (length(invalid_maturities) > 0) {
    stop(sprintf("Invalid maturities: %s. Valid options are: %s",
                 paste(invalid_maturities, collapse = ", "),
                 paste(names(maturity_map), collapse = ", ")))
  }

  log_info("Fetching Treasury yield data from FRED...")
  log_info(sprintf("Maturities: %s", paste(maturities, collapse = ", ")))

  # Fetch data for each maturity
  all_data <- list()

  for (mat in maturities) {
    series_id <- maturity_map[[mat]]
    col_name <- maturity_years[[mat]]

    tryCatch({
      df <- fetch_fred_series(
        series_id = series_id,
        api_key = api_key,
        start_date = start_date,
        end_date = end_date,
        verbose = verbose
      )

      if (nrow(df) > 0) {
        names(df)[names(df) == "value"] <- col_name
        all_data[[mat]] <- df
      }
    }, error = function(e) {
      warning(sprintf("Failed to fetch %s: %s", series_id, e$message))
    })
  }

  if (length(all_data) == 0) {
    stop("Failed to fetch any data from FRED")
  }

  # Merge all series by date
  log_info("Merging data across maturities...")

  yield_data <- all_data[[1]]
  if (length(all_data) > 1) {
    for (i in 2:length(all_data)) {
      yield_data <- merge(yield_data, all_data[[i]], by = "date", all = TRUE)
    }
  }

  # Rename date column to Date for consistency
  names(yield_data)[names(yield_data) == "date"] <- "Date"

  # Convert to decimal if requested (FRED data is already in percentage form)
  if (convert_to_decimal) {
    log_info("Converting yields from percentage to decimal...")
    for (col in setdiff(names(yield_data), "Date")) {
      yield_data[[col]] <- yield_data[[col]] / 100
    }
  }

  # Sort data by date (most recent first)
  yield_data <- yield_data[order(yield_data$Date, decreasing = TRUE), ]

  # Fill missing values if requested
  if (fill_missing) {
    log_info("Filling missing values using forward fill...")
    yield_data <- yield_data %>%
      arrange(desc(Date)) %>%
      mutate(across(where(is.numeric), ~{
        # Forward fill (carry last observation forward)
        for(i in seq_along(.)) {
          if(is.na(.)[i] && i > 1) {
            .[i] <- .[i-1]
          }
        }
        .
      }))
  }

  # Final log output
  if (verbose) {
    log_info("Final yield data structure:")
    str(yield_data)
    log_info("First few rows of final data:")
    print(head(yield_data, 3))
    log_info(sprintf("Data range: %s to %s",
                     min(yield_data$Date, na.rm = TRUE),
                     max(yield_data$Date, na.rm = TRUE)))
    log_info(sprintf("Total observations: %d", nrow(yield_data)))
  }

  return(yield_data)
}


#' @title Get Historical Treasury Yield Data from FRED
#'
#' @description
#' Downloads U.S. Treasury yield curve data from FRED API for a specified date range.
#' This is a simplified wrapper around get_treasury_yield_data() for backward compatibility.
#'
#' @param api_key Character. FRED API key. If NULL, will look for fred_api_key environment variable.
#' @param start_date Date or character in "YYYY-MM-DD" format. The earliest date to include.
#'   If NULL, defaults to "1990-01-01".
#' @param end_date Date or character in "YYYY-MM-DD" format. The latest date to include.
#'   If NULL, defaults to current date.
#' @param verbose Logical. If TRUE, print processing information.
#'
#' @return A data.frame with combined yield curve data.
#'
#' @details
#' This function is equivalent to calling get_treasury_yield_data() with default maturities.
#' It provides a simpler interface for users who just want all available maturities.
#'
#' @examples
#' \dontrun{
#' # Set your API key first
#' Sys.setenv(fred_api_key = "your_api_key_here")
#' yield_data <- get_treasury_yield_combined(start_date = "2018-01-01")
#' head(yield_data)
#' }
#'
#' @export
get_treasury_yield_combined <- function(api_key = NULL,
                                       start_date = NULL,
                                       end_date = NULL,
                                       verbose = FALSE) {
  # Set defaults
  if (is.null(start_date)) {
    start_date <- "1990-01-01"
  }

  if (is.null(end_date)) {
    end_date <- Sys.Date()
  }

  if (verbose) {
    cat(sprintf("Fetching Treasury yield data from FRED: %s to %s\n",
                start_date, end_date))
  }

  # Use the main function with all maturities
  yield_data <- get_treasury_yield_data(
    api_key = api_key,
    maturities = c("1MO", "3MO", "6MO", "1", "2", "3", "5", "7", "10", "20", "30"),
    start_date = start_date,
    end_date = end_date,
    convert_to_decimal = TRUE,
    verbose = verbose,
    fill_missing = TRUE
  )

  return(yield_data)
}
