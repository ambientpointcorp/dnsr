#' @title Treasury Yield Curve Data Utilities
#' @description Functions for downloading and processing U.S. Treasury yield curve data from FRED API.

#' @importFrom httr GET RETRY content
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
#' @param retries Integer. Number of attempts for the HTTP request. Transient
#'   failures and 5xx responses are retried with exponential backoff; a 4xx
#'   response (a bad key or series id) is not, since retrying cannot help.
#' @param cache_dir Character. Directory in which to cache responses. If given,
#'   a request already on disk is read from there instead of being re-fetched.
#'   \code{NULL} (the default) disables caching.
#'
#' @return A data.frame with date and yield value columns.
#'
#' @details
#' Requires a FRED API key which can be obtained free at https://fred.stlouisfed.org/docs/api/api_key.html
#' Set the environment variable fred_api_key or pass it directly to the function.
#'
#' Cache files are keyed on the series id and the requested date range only.
#' The API key is never part of the key and is never written to disk. A cache
#' hit is served before the key is looked up, so a cached corpus can be
#' replayed with no key configured at all.
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
                              end_date = NULL, verbose = FALSE,
                              retries = 3L, cache_dir = NULL) {
  # Parameter validation
  checkmate::assert_string(series_id, min.chars = 1)
  checkmate::assert_logical(verbose, len = 1)
  checkmate::assert_count(retries, positive = TRUE)
  if (!is.null(cache_dir)) checkmate::assert_string(cache_dir, min.chars = 1)

  # Create informative logging function
  log_info <- function(...) {
    if (verbose) cat(..., "\n")
  }

  # Serve from cache before anything else. The key covers the series and the
  # date range; the API key is deliberately excluded from it and is never
  # written to disk. A cache hit sends no request, so it needs no key -- which
  # is what lets a cached corpus be replayed offline.
  cache_file <- NULL
  if (!is.null(cache_dir)) {
    key <- paste(
      series_id,
      if (is.null(start_date)) "min" else as.character(as.Date(start_date)),
      if (is.null(end_date)) "max" else as.character(as.Date(end_date)),
      sep = "_"
    )
    key <- gsub("[^A-Za-z0-9._-]", "-", key)
    cache_file <- file.path(cache_dir, paste0("fred_", key, ".rds"))

    if (file.exists(cache_file)) {
      log_info(sprintf("Reading %s from cache", series_id))
      return(readRDS(cache_file))
    }
  }

  # Get API key
  if (is.null(api_key)) {
    api_key <- Sys.getenv("fred_api_key")
    if (api_key == "") {
      stop("FRED API key not found. Set fred_api_key environment variable or pass api_key parameter.")
    }
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

  # Make API request. RETRY backs off on transient failures; a 4xx is returned
  # as-is for the status check below, since retrying a rejected key or an
  # unknown series id cannot succeed.
  response <- tryCatch({
    httr::RETRY(
      "GET", base_url, query = params,
      times = as.integer(retries),
      terminate_on = c(400L, 401L, 403L, 404L, 422L),
      quiet = !verbose
    )
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

  if (!is.null(cache_file)) {
    dir.create(dirname(cache_file), showWarnings = FALSE, recursive = TRUE)
    saveRDS(result, cache_file)
    log_info(sprintf("Cached %s", series_id))
  }

  return(result)
}


#' @title Get Treasury Yield Curve Data from FRED
#'
#' @description
#' Downloads U.S. Treasury yield curve data from FRED API for multiple maturities
#' and processes it into a consistent format suitable for DNS model calibration.
#'
#' @param api_key Character. FRED API key. If NULL, will look for fred_api_key environment variable.
#' @param maturities Character vector. Maturities to download. Options are
#'   "1MO", "2MO", "3MO", "6MO", "1", "2", "3", "5", "7", "10", "20", "30".
#' @param start_date Date or character in "YYYY-MM-DD" format. The earliest date to include.
#' @param end_date Date or character in "YYYY-MM-DD" format. The latest date to include.
#' @param convert_to_decimal Logical. If TRUE, convert yields from percentages to decimal format.
#' @param verbose Logical. If TRUE, print processing information.
#' @param fill_missing Logical, or one of \code{"forward"} / \code{"none"}. If TRUE
#'   (the default) gaps are filled by carrying the last observation forward in
#'   calendar time.
#' @param max_gap Numeric. Maximum number of calendar days a value may be
#'   carried across when \code{fill_missing} is in effect. Gaps wider than this
#'   are left as \code{NA}. Defaults to \code{Inf}, which fills every gap.
#' @param date_order Character. Row order of the returned frame:
#'   \code{"descending"} (the default, most recent first) or \code{"ascending"}.
#' @param retries Integer. HTTP attempts per series, passed to
#'   \code{\link{fetch_fred_series}}.
#' @param cache_dir Character. Directory in which to cache responses, passed to
#'   \code{\link{fetch_fred_series}}. \code{NULL} disables caching.
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
#' Gap filling runs on an ascending date index, so a missing yield takes the
#' value of the most recent *earlier* observation. This matters because the
#' frame is returned newest-first: filling in row order would carry values
#' backwards in time and introduce look-ahead bias into anything calibrated on
#' the result. Long series have genuine multi-year holes (DGS20 before 1993,
#' DGS30 between 2002 and 2006); set \code{max_gap} to stop those being bridged.
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
                                    fill_missing = TRUE,
                                    max_gap = Inf,
                                    date_order = c("descending", "ascending"),
                                    retries = 3L,
                                    cache_dir = NULL) {
  # Parameter validation
  date_order <- match.arg(date_order)
  checkmate::assert_logical(convert_to_decimal, len = 1)
  checkmate::assert_logical(verbose, len = 1)
  checkmate::assert_character(maturities, min.len = 1)
  checkmate::assert_numeric(max_gap, lower = 0, len = 1)

  # fill_missing accepts the historical logical as well as a named method.
  if (is.character(fill_missing)) {
    fill_missing <- match.arg(fill_missing, c("forward", "none"))
    do_fill <- identical(fill_missing, "forward")
  } else {
    checkmate::assert_logical(fill_missing, len = 1, any.missing = FALSE)
    do_fill <- isTRUE(fill_missing)
  }

  # Resolve the API key once up front: a missing key is a configuration
  # error, not a per-series fetch failure, and must not be swallowed by the
  # per-maturity tryCatch below.
  if (is.null(api_key)) {
    api_key <- Sys.getenv("fred_api_key")
    if (api_key == "") {
      stop("FRED API key not found. Set fred_api_key environment variable ",
           "or pass api_key parameter.")
    }
  }

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
        verbose = verbose,
        retries = retries,
        cache_dir = cache_dir
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

  # Fill gaps on an ASCENDING date index. The carry has to run forward in
  # calendar time, which is the opposite direction from the newest-first row
  # order this function returns by default; filling in row order would take
  # each missing yield from a later date.
  if (do_fill) {
    log_info("Filling missing values (last observation carried forward)...")
    yield_data <- yield_data[order(yield_data$Date), , drop = FALSE]
    for (col in setdiff(names(yield_data), "Date")) {
      if (is.numeric(yield_data[[col]])) {
        yield_data[[col]] <- locf_by_date(
          yield_data[[col]], yield_data$Date, max_gap = max_gap
        )
      }
    }
  }

  # Order the output as requested; the default reproduces the historical
  # most-recent-first convention.
  yield_data <- yield_data[
    order(yield_data$Date, decreasing = identical(date_order, "descending")), ,
    drop = FALSE
  ]
  rownames(yield_data) <- NULL

  # Final log output
  if (verbose) {
    log_info(sprintf("Final yield data: %d columns (%s)",
                     ncol(yield_data), paste(names(yield_data), collapse = ", ")))
    log_info("First few rows of final data:")
    print(utils::head(yield_data, 3))
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
#' @param ... Further arguments passed to \code{\link{get_treasury_yield_data}},
#'   for example \code{max_gap}, \code{date_order} or \code{fill_missing}.
#'
#' @return A data.frame with combined yield curve data, most recent date first.
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
                                       verbose = FALSE,
                                       ...) {
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
    ...
  )

  return(yield_data)
}


# Carry the last observation forward along an ascending date index.
#
# `x` and `dates` must be ordered by date, ascending. `max_gap` caps how many
# days a value may be carried across, so a multi-year hole in a series (DGS20
# before 1993, DGS30 between 2002 and 2006) is not silently bridged.
locf_by_date <- function(x, dates, max_gap = Inf) {
  n <- length(x)
  if (n == 0L) return(x)

  observed <- which(!is.na(x))
  if (length(observed) == 0L) return(x)

  # For each position, the index of the most recent observed value at or
  # before it. Positions before the first observation get NA and stay NA:
  # there is nothing earlier to carry.
  rank_of_last <- cumsum(!is.na(x))
  rank_of_last[rank_of_last == 0L] <- NA_integer_
  source_idx <- observed[rank_of_last]

  gap_days <- as.numeric(dates - dates[source_idx])
  fillable <- is.na(x) & !is.na(source_idx) & gap_days <= max_gap
  x[fillable] <- x[source_idx[fillable]]
  x
}
