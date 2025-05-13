#' @title Treasury Yield Curve Data Utilities
#' @description Functions for downloading and processing U.S. Treasury yield curve data.

#' @importFrom httr GET content
#' @importFrom readr read_csv
#' @importFrom dplyr %>% mutate across where lag filter select all_of bind_rows
#' @importFrom lubridate mdy date
#' @importFrom tidyr replace_na
#' @importFrom stats na.omit

#' @title Download and Process Treasury Yield Curve Data
#'
#' @description
#' Downloads U.S. Treasury yield curve data from the Treasury website and processes it into
#' a consistent format suitable for DNS model calibration.
#'
#' @param url Character. URL for the Treasury yield curve data CSV file.
#' @param convert_to_decimal Logical. If TRUE, convert yields from percentages to decimal format.
#' @param verbose Logical. If TRUE, print processing information.
#' @param fill_missing Logical. If TRUE, attempt to fill in missing values.
#'
#' @return A data.frame with processed yield curve data where:
#'   - The first column is 'Date' containing the observation dates.
#'   - Subsequent columns contain yield values for different maturities.
#'   - Column names represent maturities in years.
#'
#' @details
#' This function handles the inconsistencies in Treasury data formatting, cleans column names
#' to represent maturities in years, and processes the data for DNS model fitting.
#'
#' @examples
#' \dontrun{
#' url <- "https://home.treasury.gov/system/files/276/yield-curve-rates-1990-2024.csv"
#' yield_data <- get_treasury_yield_data(url)
#' head(yield_data)
#' }
#'
#' @export
get_treasury_yield_data <- function(url, convert_to_decimal = TRUE, verbose = FALSE, fill_missing = TRUE) {
  # Parameter validation
  checkmate::assert_string(url, min.chars = 10)
  checkmate::assert_logical(convert_to_decimal, len = 1)
  checkmate::assert_logical(verbose, len = 1)
  checkmate::assert_logical(fill_missing, len = 1)

  # Create informative logging function
  log_info <- function(...) {
    if (verbose) cat(..., "\n")
  }

  # Download and read the file with error handling
  log_info("Downloading data from:", url)

  temp_file <- tempfile(fileext = '.csv')

  # Download the file with comprehensive error handling
  tryCatch({
    download.file(url, temp_file, quiet = !verbose, mode = "wb")

    # Check if file exists and has content
    if (!file.exists(temp_file) || file.size(temp_file) == 0) {
      stop("Failed to download data or empty file received")
    }

    # Try to read the CSV
    yield_data <- tryCatch({
      read_csv(temp_file, show_col_types = FALSE)
    }, error = function(e) {
      stop("Failed to read CSV data: ", e$message)
    })
  },
  warning = function(w) {
    message("Warning in download: ", w$message)
  },
  error = function(e) {
    if (file.exists(temp_file)) file.remove(temp_file)
    stop("Failed to download or process data: ", e$message)
  },
  finally = {
    # Make sure we clean up the temporary file
    if (file.exists(temp_file)) file.remove(temp_file)
  })

  # Show structure information in verbose mode
  if (verbose) {
    log_info("Structure of loaded data:")
    str(yield_data)
    log_info("First few rows of data:")
    print(head(yield_data, 3))
    log_info("Number of columns:", ncol(yield_data))
    log_info("Column names:", paste(names(yield_data), collapse = ", "))
  }

  # Ensure Date column exists
  if (!"Date" %in% names(yield_data)) {
    # Try to find a date column with different name
    date_col_idx <- grep("date", names(yield_data), ignore.case = TRUE)[1]
    if (!is.na(date_col_idx)) {
      names(yield_data)[date_col_idx] <- "Date"
    } else {
      stop("Cannot find a date column in the data")
    }
  }

  # Identify maturity columns (all columns except Date)
  maturity_cols <- setdiff(names(yield_data), "Date")

  # Create clean maturity labels based on column names
  clean_maturities <- sapply(maturity_cols, function(col) {
    # Try to extract maturity in years
    if (grepl("mo", col, ignore.case = TRUE)) {
      # Column has months
      months <- as.numeric(gsub("[^0-9]", "", col))
      if (!is.na(months)) {
        return(round(months / 12, 4))
      }
    } else if (grepl("yr", col, ignore.case = TRUE) || grepl("year", col, ignore.case = TRUE)) {
      # Column has years
      years <- as.numeric(gsub("[^0-9]", "", col))
      if (!is.na(years)) {
        return(years)
      }
    }

    # Try to extract numeric part as a fallback
    num <- as.numeric(gsub("[^0-9.]", "", col))
    if (!is.na(num)) {
      if (num < 12 && !grepl("yr", col, ignore.case = TRUE)) {
        # Likely months if small number and not explicitly years
        return(round(num / 12, 4))
      } else {
        return(num)
      }
    }

    # If all else fails, return NA
    return(NA)
  })

  # Check if we have any valid maturities
  if (all(is.na(clean_maturities))) {
    stop("Failed to extract any valid maturity information from column names")
  }

  # Create new column names
  new_names <- c("Date", as.character(clean_maturities))

  # Handle duplicates in the new column names
  if (length(new_names) != length(unique(new_names))) {
    log_info("Warning: Duplicate column names detected after cleaning")

    # Find and resolve duplicates
    for (i in 2:length(new_names)) {
      if (new_names[i] %in% new_names[1:(i-1)]) {
        # If duplicate, slightly modify the value
        new_names[i] <- sprintf("%.4f", as.numeric(new_names[i]) + 0.0001)
      }
    }
  }

  # Apply the new column names
  names(yield_data) <- new_names

  # Convert yield values from percentages to decimals if requested
  if (convert_to_decimal) {
    for (col in names(yield_data)[-1]) {
      if (is.numeric(yield_data[[col]])) {
        yield_data[[col]] <- yield_data[[col]] / 100
      }
    }
  }

  # Handle date formatting
  if (is.character(yield_data$Date)) {
    # Try to parse date with lubridate
    tryCatch({
      # First try mdy format (common for Treasury data)
      yield_data$Date <- mdy(yield_data$Date)
    }, error = function(e) {
      # If that fails, try other formats
      tryCatch({
        yield_data$Date <- as.Date(yield_data$Date)
      }, error = function(e2) {
        warning("Could not parse dates in expected formats. Date column may not be properly formatted.")
      })
    })
  }

  # Sort data by date (most recent first)
  yield_data <- yield_data[order(yield_data$Date, decreasing = TRUE), ]

  # Fill missing values if requested
  if (fill_missing) {
    # Replace missing values with previous value
    yield_data <- yield_data %>%
      arrange(desc(Date)) %>%
      mutate(across(where(is.numeric), ~ifelse(is.na(.), lag(.), .)))
  }

  # Convert any character columns to numeric (except Date)
  for (col in setdiff(names(yield_data), "Date")) {
    if (is.character(yield_data[[col]])) {
      yield_data[[col]] <- as.numeric(yield_data[[col]])
    }
  }

  # Final log output
  if (verbose) {
    log_info("Final yield data structure:")
    str(yield_data)
    log_info("First few rows of final data:")
    print(head(yield_data, 3))
    log_info("Data range:", min(yield_data$Date, na.rm = TRUE), "to",
             max(yield_data$Date, na.rm = TRUE))
  }

  return(yield_data)
}

#' @title Get Historical and Current Treasury Yield Data
#'
#' @description
#' Downloads and combines historical and current U.S. Treasury yield curve data.
#'
#' @param start_date Date or character in "YYYY-MM-DD" format. The earliest date to include.
#' @param end_date Date or character in "YYYY-MM-DD" format. The latest date to include.
#' @param verbose Logical. If TRUE, print processing information.
#'
#' @return A data.frame with combined yield curve data.
#'
#' @examples
#' \dontrun{
#' yield_data <- get_treasury_yield_combined(start_date = "2018-01-01")
#' head(yield_data)
#' }
#'
#' @export
get_treasury_yield_combined <- function(start_date = NULL, end_date = NULL, verbose = FALSE) {
  # URLs for Treasury data sources
  yc_url_1990_2024 <- "https://home.treasury.gov/system/files/276/yield-curve-rates-1990-2024.csv"
  yc_url_2025 <- "https://home.treasury.gov/resource-center/data-chart-center/interest-rates/daily-treasury-rates.csv/2025/all?type=daily_treasury_yield_curve&field_tdr_date_value=2025&page&_format=csv"

  # Get historical data
  if (verbose) cat("Downloading historical Treasury yield curve data...\n")
  yield_curve_daily_hist <- get_treasury_yield_data(yc_url_1990_2024, verbose = verbose)

  # Get current year data
  if (verbose) cat("Downloading current year Treasury yield curve data...\n")
  yield_curve_daily_current <- try(
    get_treasury_yield_data(yc_url_2025, verbose = verbose),
    silent = TRUE
  )

  # Combine data
  if (verbose) cat("Combining data sets...\n")
  if (inherits(yield_curve_daily_current, "try-error")) {
    if (verbose) cat("Warning: Could not download current year data. Using only historical data.\n")
    yield_data <- yield_curve_daily_hist
  } else {
    # Find common columns to ensure compatibility
    common_cols <- intersect(
      names(yield_curve_daily_hist),
      names(yield_curve_daily_current)
    )

    # Ensure 'Date' is included
    if (!"Date" %in% common_cols) {
      stop("Date column missing from one or both data sets")
    }

    # Combine data sets
    yield_data <- bind_rows(
      yield_curve_daily_current[, common_cols],
      yield_curve_daily_hist[, common_cols]
    )
  }

  # Filter by date range if specified
  if (!is.null(start_date)) {
    start_date <- as.Date(start_date)
    yield_data <- yield_data[yield_data$Date >= start_date, ]
  }

  if (!is.null(end_date)) {
    end_date <- as.Date(end_date)
    yield_data <- yield_data[yield_data$Date <= end_date, ]
  }

  # Fill missing values (can occur around weekends or holidays)
  yield_data <- yield_data %>%
    arrange(Date) %>%
    mutate(across(where(is.numeric), ~ifelse(is.na(.), lag(.), .)))

  # Clean up unused columns (uncommon maturities)
  col_statistics <- sapply(yield_data[,-1], function(col) {
    sum(!is.na(col)) / length(col)
  })

  # Remove columns with more than 50% missing values
  cols_to_keep <- c("Date", names(col_statistics)[col_statistics > 0.5])
  yield_data <- yield_data[, cols_to_keep]

  # Final data summary
  if (verbose) {
    cat("\nYield Data Summary:\n")
    cat("Number of days:", nrow(yield_data), "\n")
    cat("Date range:", min(yield_data$Date), "to", max(yield_data$Date), "\n")
    cat("Maturities:", paste(names(yield_data)[-1], collapse = ", "), "\n\n")
  }

  return(yield_data)
}
