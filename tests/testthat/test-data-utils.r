library(testthat)
library(dnsr)
library(mockery)

test_that("fetch_fred_series works with valid API key", {
  skip_if_not(nzchar(Sys.getenv("fred_api_key")), "No FRED API key set")
  skip_on_cran()

  result <- fetch_fred_series(
    series_id = "DGS10",
    start_date = "2023-01-01",
    end_date = "2023-01-31"
  )

  expect_s3_class(result, "data.frame")
  expect_true("date" %in% names(result))
  expect_true("value" %in% names(result))
  expect_gt(nrow(result), 0)
})

test_that("get_treasury_yield_data fetches multiple maturities", {
  skip_if_not(nzchar(Sys.getenv("fred_api_key")), "No FRED API key set")
  skip_on_cran()

  result <- get_treasury_yield_data(
    maturities = c("2", "10"),
    start_date = "2023-01-01",
    end_date = "2023-01-31"
  )

  expect_s3_class(result, "data.frame")
  expect_true("Date" %in% names(result))
  expect_true("2" %in% names(result))
  expect_true("10" %in% names(result))
})

test_that("get_treasury_yield_data handles missing API key", {
  # Temporarily unset API key
  old_key <- Sys.getenv("fred_api_key")
  Sys.unsetenv("fred_api_key")

  expect_error(
    get_treasury_yield_data(maturities = "10"),
    "FRED API key not found"
  )

  # Restore API key
  if (nzchar(old_key)) Sys.setenv(fred_api_key = old_key)
})

test_that("get_treasury_yield_data validates maturity codes", {
  skip_if_not(nzchar(Sys.getenv("fred_api_key")), "No FRED API key set")

  expect_error(
    get_treasury_yield_data(maturities = c("99", "invalid")),
    "Invalid maturities"
  )
})

test_that("get_treasury_yield_combined uses defaults correctly", {
  skip_if_not(nzchar(Sys.getenv("fred_api_key")), "No FRED API key set")
  skip_on_cran()

  result <- get_treasury_yield_combined(
    start_date = "2023-01-01",
    end_date = "2023-01-31"
  )

  expect_s3_class(result, "data.frame")
  expect_true("Date" %in% names(result))
  # Should have multiple maturity columns
  expect_gt(ncol(result), 5)
})
