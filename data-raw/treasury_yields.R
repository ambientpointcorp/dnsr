# Builds data/treasury_yields.rda.
#
# Source: Federal Reserve Economic Data (FRED), series DGS1MO, DGS3MO, DGS6MO,
# DGS1, DGS2, DGS3, DGS5, DGS7, DGS10, DGS20, DGS30 -- the Treasury constant
# maturity yields published in the Federal Reserve's H.15 release. These are
# U.S. government works and are in the public domain.
#
# Re-run with a FRED API key in the fred_api_key environment variable:
#   Rscript data-raw/treasury_yields.R

library(dnsr)

treasury_yields <- get_treasury_yield_data(
  maturities = c("1MO", "3MO", "6MO", "1", "2", "3", "5", "7", "10", "20", "30"),
  start_date = "2019-01-01",
  end_date = "2024-12-31",
  convert_to_decimal = TRUE,
  fill_missing = FALSE,
  date_order = "ascending"
)

# FRED emits a row for every weekday, with "." for market holidays. Those rows
# carry no observation at any maturity, so they are dropped rather than shipped
# as 65 rows of NA.
yield_cols <- setdiff(names(treasury_yields), "Date")
all_missing <- rowSums(is.na(treasury_yields[, yield_cols])) == length(yield_cols)
treasury_yields <- treasury_yields[!all_missing, , drop = FALSE]
rownames(treasury_yields) <- NULL

stopifnot(
  nrow(treasury_yields) > 1400,
  !anyNA(treasury_yields),
  !is.unsorted(treasury_yields$Date)
)

save(treasury_yields, file = "data/treasury_yields.rda",
     compress = "xz", version = 2)
