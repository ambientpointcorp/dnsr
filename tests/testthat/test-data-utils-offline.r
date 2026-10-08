# Offline coverage of the FRED data path, which the live tests skip and which
# was the least-covered file in the package.
#
# Two routes, neither needing a key or a network:
#   - the on-disk cache, which fetch_fred_series() consults before it looks for
#     a key, so a prepared cache exercises the whole assembly in
#     get_treasury_yield_data() for real;
#   - mockery, for the request and parse layer the cache skips over.

# A FRED observations payload as httr::content(as = "parsed") returns it.
fred_payload <- function(dates, values) {
  list(observations = Map(
    function(d, v) list(date = d, value = v),
    as.character(dates), as.character(values)
  ))
}

# A cache directory holding one series per file, keyed the way
# fetch_fred_series() keys them.
seed_cache <- function(series, start, end) {
  dir <- file.path(tempdir(), paste0("dnsr-cache-", basename(tempfile(""))))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  for (id in names(series)) {
    key <- gsub("[^A-Za-z0-9._-]", "-", paste(id, start, end, sep = "_"))
    saveRDS(series[[id]], file.path(dir, paste0("fred_", key, ".rds")))
  }
  dir
}


# ---- the request and parse layer --------------------------------------------

test_that("fetch_fred_series parses observations and maps '.' to NA", {
  skip_if_not_installed("mockery")

  dates <- as.Date("2024-01-02") + 0:3
  payload <- fred_payload(dates, c("3.95", ".", "4.01", "4.05"))

  mockery::stub(fetch_fred_series, "httr::RETRY", structure(list(), class = "response"))
  mockery::stub(fetch_fred_series, "httr::http_error", FALSE)
  mockery::stub(fetch_fred_series, "httr::content", payload)

  result <- fetch_fred_series("DGS10", api_key = "test-key",
                              start_date = "2024-01-02", end_date = "2024-01-05")

  expect_s3_class(result, "data.frame")
  expect_named(result, c("date", "value"))
  expect_s3_class(result$date, "Date")
  expect_equal(result$date, dates)
  # A FRED "." means no observation that day, not zero.
  expect_equal(result$value, c(3.95, NA, 4.01, 4.05))
})

test_that("fetch_fred_series reports an HTTP failure rather than returning it", {
  skip_if_not_installed("mockery")

  mockery::stub(fetch_fred_series, "httr::RETRY", structure(list(), class = "response"))
  mockery::stub(fetch_fred_series, "httr::http_error", TRUE)
  mockery::stub(fetch_fred_series, "httr::status_code", 500L)
  mockery::stub(fetch_fred_series, "httr::http_status",
                list(message = "Server error: (500) Internal Server Error"))

  expect_error(
    fetch_fred_series("DGS10", api_key = "test-key"),
    "FRED API request failed with status 500"
  )
})

test_that("fetch_fred_series warns and returns an empty frame for no data", {
  skip_if_not_installed("mockery")

  mockery::stub(fetch_fred_series, "httr::RETRY", structure(list(), class = "response"))
  mockery::stub(fetch_fred_series, "httr::http_error", FALSE)
  mockery::stub(fetch_fred_series, "httr::content", list(observations = list()))

  expect_warning(
    empty <- fetch_fred_series("DGS10", api_key = "test-key"),
    "No data returned"
  )
  expect_s3_class(empty, "data.frame")
  expect_equal(nrow(empty), 0)
  expect_s3_class(empty$date, "Date")
})

test_that("fetch_fred_series writes and then reads its cache", {
  skip_if_not_installed("mockery")

  cache_dir <- file.path(tempdir(), "dnsr-cache-roundtrip")
  unlink(cache_dir, recursive = TRUE)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  dates <- as.Date("2024-01-02") + 0:2
  mockery::stub(fetch_fred_series, "httr::RETRY", structure(list(), class = "response"))
  mockery::stub(fetch_fred_series, "httr::http_error", FALSE)
  mockery::stub(fetch_fred_series, "httr::content",
                fred_payload(dates, c("3.95", "4.01", "4.05")))

  first <- fetch_fred_series("DGS10", api_key = "test-key",
                             start_date = "2024-01-02", end_date = "2024-01-04",
                             cache_dir = cache_dir)
  expect_equal(nrow(first), 3)
  expect_length(list.files(cache_dir), 1)

  # The key covers the series and date range only, so nothing secret is on disk.
  written <- readLines(list.files(cache_dir, full.names = TRUE)[1],
                       warn = FALSE, n = -1L)
  expect_false(any(grepl("test-key", written, fixed = TRUE)))

  # Second call is served from disk: no key needed, because no request is made.
  second <- fetch_fred_series("DGS10", start_date = "2024-01-02",
                              end_date = "2024-01-04", cache_dir = cache_dir)
  expect_equal(second, first)
})

test_that("fetch_fred_series validates its arguments", {
  expect_error(fetch_fred_series(""), "series_id")
  expect_error(fetch_fred_series("DGS10", api_key = "k", retries = 0), "retries")

  old_key <- Sys.getenv("fred_api_key")
  Sys.unsetenv("fred_api_key")
  on.exit(if (nzchar(old_key)) Sys.setenv(fred_api_key = old_key), add = TRUE)
  expect_error(fetch_fred_series("DGS10"), "API key not found")
})


# ---- the assembly, through the cache ----------------------------------------

test_that("get_treasury_yield_data assembles, converts and orders the panel", {
  dates <- as.Date("2024-01-02") + 0:3
  cache_dir <- seed_cache(
    list(
      DGS2 = data.frame(date = dates, "2" = c(4.30, 4.32, 4.28, 4.31),
                        check.names = FALSE),
      DGS10 = data.frame(date = dates, "10" = c(3.95, 4.01, 4.05, 4.02),
                         check.names = FALSE)
    ),
    "2024-01-02", "2024-01-05"
  )
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  panel <- get_treasury_yield_data(
    api_key = "unused-because-cached", maturities = c("2", "10"),
    start_date = "2024-01-02", end_date = "2024-01-05", cache_dir = cache_dir
  )

  expect_named(panel, c("Date", "2", "10"))
  expect_s3_class(panel$Date, "Date")
  expect_equal(nrow(panel), 4)

  # Percentages become decimals.
  expect_equal(panel[["10"]][panel$Date == dates[1]], 0.0395)

  # Most recent first, which is the documented default.
  expect_equal(panel$Date, sort(dates, decreasing = TRUE))

  ascending <- get_treasury_yield_data(
    api_key = "unused", maturities = c("2", "10"),
    start_date = "2024-01-02", end_date = "2024-01-05",
    cache_dir = cache_dir, date_order = "ascending"
  )
  expect_equal(ascending$Date, sort(dates))

  # convert_to_decimal = FALSE leaves percentage points alone.
  raw <- get_treasury_yield_data(
    api_key = "unused", maturities = c("10"),
    start_date = "2024-01-02", end_date = "2024-01-05",
    cache_dir = cache_dir, convert_to_decimal = FALSE
  )
  expect_equal(max(raw[["10"]]), 4.05)
})

test_that("get_treasury_yield_data fills gaps forward and honours max_gap", {
  dates <- as.Date("2024-01-02") + 0:4
  cache_dir <- seed_cache(
    list(DGS10 = data.frame(date = dates,
                            "10" = c(3.95, NA, NA, NA, 4.05),
                            check.names = FALSE)),
    "2024-01-02", "2024-01-06"
  )
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  filled <- get_treasury_yield_data(
    api_key = "unused", maturities = "10", start_date = "2024-01-02",
    end_date = "2024-01-06", cache_dir = cache_dir, date_order = "ascending"
  )
  # Carried from the last earlier observation, never from a later one.
  expect_equal(filled[["10"]], c(0.0395, 0.0395, 0.0395, 0.0395, 0.0405))

  capped <- get_treasury_yield_data(
    api_key = "unused", maturities = "10", start_date = "2024-01-02",
    end_date = "2024-01-06", cache_dir = cache_dir, date_order = "ascending",
    max_gap = 2
  )
  expect_equal(capped[["10"]], c(0.0395, 0.0395, 0.0395, NA, 0.0405))

  unfilled <- get_treasury_yield_data(
    api_key = "unused", maturities = "10", start_date = "2024-01-02",
    end_date = "2024-01-06", cache_dir = cache_dir, date_order = "ascending",
    fill_missing = FALSE
  )
  expect_equal(sum(is.na(unfilled[["10"]])), 3)

  # The character form of the argument is accepted too.
  named <- get_treasury_yield_data(
    api_key = "unused", maturities = "10", start_date = "2024-01-02",
    end_date = "2024-01-06", cache_dir = cache_dir, date_order = "ascending",
    fill_missing = "none"
  )
  expect_equal(named[["10"]], unfilled[["10"]])
})

test_that("get_treasury_yield_data rejects unknown maturities before fetching", {
  expect_error(
    get_treasury_yield_data(api_key = "k", maturities = c("99", "invalid")),
    "Invalid maturities"
  )
  # 2MO is documented and must pass validation. It has no cache entry here, so
  # the call still fails -- but on the fetch, not on the maturity name.
  failure <- tryCatch(
    suppressWarnings(get_treasury_yield_data(
      api_key = "k", maturities = "2MO", cache_dir = tempfile()
    )),
    error = function(e) conditionMessage(e)
  )
  expect_false(grepl("Invalid maturities", failure))
})

test_that("get_treasury_yield_data reports when every series fails", {
  skip_if_not_installed("mockery")

  mockery::stub(get_treasury_yield_data, "fetch_fred_series",
                function(...) stop("network down"))

  warnings_seen <- character(0)
  outcome <- withCallingHandlers(
    tryCatch(
      get_treasury_yield_data(api_key = "k", maturities = c("2", "10")),
      error = function(e) conditionMessage(e)
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_match(outcome, "Failed to fetch any data")
  # One warning per series, naming each, before the overall failure.
  expect_length(warnings_seen, 2)
  expect_true(all(grepl("network down", warnings_seen)))
})

test_that("get_treasury_yield_combined passes its arguments through", {
  dates <- as.Date("2024-01-02") + 0:2
  series <- stats::setNames(
    lapply(c("1MO", "3MO", "6MO", "1", "2", "3", "5", "7", "10", "20", "30"),
           function(m) NULL),
    c("DGS1MO", "DGS3MO", "DGS6MO", "DGS1", "DGS2", "DGS3", "DGS5", "DGS7",
      "DGS10", "DGS20", "DGS30")
  )
  columns <- c("0.0833", "0.25", "0.5", "1", "2", "3", "5", "7", "10", "20", "30")
  for (i in seq_along(series)) {
    frame <- data.frame(date = dates, v = c(4.0, 4.1, 4.2))
    names(frame)[2] <- columns[i]
    series[[i]] <- frame
  }
  cache_dir <- seed_cache(series, "2024-01-02", "2024-01-04")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  panel <- get_treasury_yield_combined(
    api_key = "unused", start_date = "2024-01-02", end_date = "2024-01-04",
    cache_dir = cache_dir
  )

  expect_equal(ncol(panel), 12)
  expect_equal(names(panel)[1], "Date")
  expect_equal(panel$Date, sort(dates, decreasing = TRUE))
  expect_true(all(unlist(panel[, -1]) < 1))
})
