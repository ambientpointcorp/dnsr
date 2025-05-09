library(testthat)
library(dnsr)
library(mockery)

test_that("get_treasury_yield_data handles well-formatted data", {
  # Create mock data
  mock_csv <- "Date,1 Mo,2 Mo,3 Mo,6 Mo,1 Yr,2 Yr,3 Yr,5 Yr,7 Yr,10 Yr,20 Yr,30 Yr
  05/04/2022,0.75,0.99,1.32,1.93,2.32,2.71,2.91,3.07,3.07,3.00,3.18,3.04
  05/03/2022,0.79,0.97,1.36,1.96,2.28,2.83,2.99,3.10,3.08,2.96,3.16,3.00"
  
  # Create temporary file
  temp_file <- tempfile(fileext = ".csv")
  writeLines(mock_csv, temp_file)
  
  # Mock download.file function
  with_mock(
    `download.file` = function(url, destfile, ...) {
      file.copy(temp_file, destfile)
      return(0)
    },
    {
      # Test function
      result <- get_treasury_yield_data("mock_url", verbose = FALSE)
      
      # Check result structure
      expect_s3_class(result, "data.frame")
      expect_true("Date" %in% names(result))
      expect_equal(nrow(result), 2)
      expect_equal(ncol(result), 13)  # Date + 12 maturity columns
      
      # Check maturity conversion
      expected_columns <- c("Date", "0.0833", "0.1667", "0.25", "0.5", "1", 
                           "2", "3", "5", "7", "10", "20", "30")
      expect_equal(names(result), expected_columns)
      
      # Check decimal conversion
      expect_true(all(result[1, -1] < 1))  # All yield values converted to decimals
    }
  )
  
  # Clean up
  unlink(temp_file)
})

test_that("get_treasury_yield_data handles missing values", {
  # Create mock data with missing values
  mock_csv <- "Date,1 Mo,2 Mo,3 Mo,6 Mo,1 Yr,2 Yr,3 Yr,5 Yr,7 Yr,10 Yr,20 Yr,30 Yr
  05/04/2022,0.75,,1.32,1.93,2.32,2.71,2.91,3.07,3.07,3.00,3.18,3.04
  05/03/2022,0.79,0.97,1.36,1.96,,2.83,2.99,3.10,3.08,2.96,3.16,3.00"
  
  # Create temporary file
  temp_file <- tempfile(fileext = ".csv")
  writeLines(mock_csv, temp_file)
  
  # Mock download.file function
  with_mock(
    `download.file` = function(url, destfile, ...) {
      file.copy(temp_file, destfile)
      return(0)
    },
    {
      # Test function with fill_missing = TRUE
      result <- get_treasury_yield_data("mock_url", verbose = FALSE, fill_missing = TRUE)
      
      # Check that gaps were filled
      expect_false(any(is.na(result[1, -1])))  # No NAs in first row yields
      expect_false(any(is.na(result[2, -1])))  # No NAs in second row yields
      
      # Test function with fill_missing = FALSE
      result_no_fill <- get_treasury_yield_data("mock_url", verbose = FALSE, fill_missing = FALSE)
      
      # Check that gaps were not filled
      expect_true(is.na(result_no_fill[1, "0.1667"]))  # 2 Mo in first row should be NA
      expect_true(is.na(result_no_fill[2, "1"]))       # 1 Yr in second row should be NA
    }
  )
  
  # Clean up
  unlink(temp_file)
})

test_that("get_treasury_yield_data handles download errors gracefully", {
  # Mock failed download
  with_mock(
    `download.file` = function(url, destfile, ...) {
      stop("Simulated download error")
    },
    {
      # Test function with error
      expect_error(get_treasury_yield_data("mock_url", verbose = FALSE), 
                   "Failed to download or process data")
    }
  )
})

test_that("get_treasury_yield_combined combines datasets correctly", {
  # Create mock historical data
  mock_hist <- "Date,1 Mo,2 Mo,3 Mo,6 Mo,1 Yr,2 Yr,3 Yr,5 Yr,7 Yr,10 Yr,20 Yr,30 Yr
  12/31/2023,5.55,5.53,5.41,5.29,4.79,4.25,3.97,3.84,3.87,3.88,4.27,4.03
  12/30/2023,5.54,5.52,5.40,5.28,4.78,4.24,3.96,3.83,3.86,3.87,4.26,4.02"
  
  # Create mock current data
  mock_current <- "Date,1 Mo,2 Mo,3 Mo,6 Mo,1 Yr,2 Yr,3 Yr,5 Yr,7 Yr,10 Yr,20 Yr,30 Yr
  01/02/2024,5.56,5.54,5.42,5.30,4.80,4.26,3.98,3.85,3.88,3.89,4.28,4.04
  01/01/2024,5.57,5.55,5.43,5.31,4.81,4.27,3.99,3.86,3.89,3.90,4.29,4.05"
  
  # Create temporary files
  temp_hist <- tempfile(fileext = ".csv")
  temp_current <- tempfile(fileext = ".csv")
  writeLines(mock_hist, temp_hist)
  writeLines(mock_current, temp_current)
  
  # Create mock functions
  with_mock(
    `get_treasury_yield_data` = function(url, ...) {
      if (grepl("1990-2024", url)) {
        # Read the historical data
        read.csv(temp_hist, stringsAsFactors = FALSE)
      } else {
        # Read the current year data
        read.csv(temp_current, stringsAsFactors = FALSE)
      }
    },
    {
      # Test function
      result <- get_treasury_yield_combined(verbose = FALSE)
      
      # Check result structure
      expect_s3_class(result, "data.frame")
      expect_equal(nrow(result), 4)  # 2 rows from each file
      
      # Test with date range
      result_filtered <- get_treasury_yield_combined(
        start_date = "2023-12-31", 
        end_date = "2024-01-01",
        verbose = FALSE
      )
      
      expect_equal(nrow(result_filtered), 2)  # Only rows within date range
    }
  )
  
  # Clean up
  unlink(temp_hist)
  unlink(temp_current)
})
