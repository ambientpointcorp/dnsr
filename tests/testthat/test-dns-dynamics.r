# Tests for estimate_dns_dynamics() sampling policy

# Simulate an OU path at daily spacing, optionally with iid measurement
# noise on the observed level (the situation month-end sampling corrects).
.sim_ou_params_df <- function(n_days, kappa = c(1.5, 2.0, 2.5),
                              theta = c(5, -2, 1),
                              sigma = c(1.0, 1.5, 2.0),
                              noise_sd = 0, seed = 42) {
  set.seed(seed)
  dt <- 1 / 252
  dates <- seq(as.Date("2018-01-02"), by = "day", length.out = ceiling(n_days * 1.45))
  dates <- dates[!weekdays(dates) %in% c("Saturday", "Sunday")][seq_len(n_days)]

  betas <- matrix(NA_real_, nrow = n_days, ncol = 3)
  betas[1, ] <- theta
  for (i in 2:n_days) {
    betas[i, ] <- betas[i - 1, ] + kappa * (theta - betas[i - 1, ]) * dt +
      sigma * sqrt(dt) * rnorm(3)
  }
  obs <- betas + matrix(rnorm(3 * n_days, 0, noise_sd), ncol = 3)

  df <- data.frame(Date = dates, beta1 = obs[, 1], beta2 = obs[, 2], beta3 = obs[, 3])
  for (f in c("beta1", "beta2", "beta3")) {
    df[[paste0(f, "_change")]] <- c(NA, diff(df[[f]]))
  }
  df
}

test_that("monthly sampling deflates the kappa inflation from measurement noise", {
  df <- .sim_ou_params_df(n_days = 1512, noise_sd = 0.15)  # ~6 years

  daily   <- estimate_dns_dynamics(df, sampling = "daily")
  monthly <- estimate_dns_dynamics(df, sampling = "monthly")

  for (f in c("beta1", "beta2", "beta3")) {
    # Daily changes read the iid noise as strong mean reversion
    expect_gt(daily[[f]]$kappa, 8)
    # Month-end sampling recovers the economically meaningful scale
    expect_lt(monthly[[f]]$kappa, daily[[f]]$kappa / 2)
    expect_lt(monthly[[f]]$kappa, 8)
  }
})

test_that("monthly subsampling is idempotent", {
  df <- .sim_ou_params_df(n_days = 1008, noise_sd = 0.1)

  # Pre-subsample exactly the way 0.1.x callers did
  month_key <- format(df$Date, "%Y-%m")
  monthly_df <- df[!duplicated(month_key, fromLast = TRUE), ]
  for (f in c("beta1", "beta2", "beta3")) {
    monthly_df[[paste0(f, "_change")]] <- c(NA, diff(monthly_df[[f]]))
  }

  from_daily   <- estimate_dns_dynamics(df, sampling = "monthly")
  from_monthly <- estimate_dns_dynamics(monthly_df, min_observations = 20, dt = 1 / 12,
                                        sampling = "monthly")

  for (f in c("beta1", "beta2", "beta3")) {
    expect_equal(from_daily[[f]]$kappa, from_monthly[[f]]$kappa)
    expect_equal(from_daily[[f]]$theta, from_monthly[[f]]$theta)
    expect_equal(from_daily[[f]]$sigma, from_monthly[[f]]$sigma)
  }
  expect_equal(from_daily$innovations, from_monthly$innovations)
  expect_equal(from_daily$innovation_dates, from_monthly$innovation_dates)
})

test_that("innovation_dates align with innovations and are month-ends of the input", {
  df <- .sim_ou_params_df(n_days = 1008)
  d <- estimate_dns_dynamics(df)  # default monthly

  expect_equal(length(d$innovation_dates), nrow(d$innovations))
  expect_true(all(d$innovation_dates %in% df$Date))
  # One observation per month, in order
  expect_false(any(duplicated(format(d$innovation_dates, "%Y-%m"))))
  expect_true(!is.unsorted(d$innovation_dates))
  expect_equal(d$sampling, "monthly")
  expect_equal(d$dt, 1 / 12)
})

test_that("monthly mode ignores a conflicting dt with a warning", {
  df <- .sim_ou_params_df(n_days = 1008)

  expect_warning(d <- estimate_dns_dynamics(df, dt = 1 / 252, sampling = "monthly"),
                 "month-end")
  expect_equal(d$dt, 1 / 12)
  ref <- estimate_dns_dynamics(df, sampling = "monthly")
  expect_equal(d$beta1$kappa, ref$beta1$kappa)

  # dt = 1/12 is consistent and accepted silently
  expect_silent(estimate_dns_dynamics(df, dt = 1 / 12, sampling = "monthly"))
})

test_that("change columns are only required in daily mode", {
  df <- .sim_ou_params_df(n_days = 1008)
  bare <- df[, c("Date", "beta1", "beta2", "beta3")]

  expect_silent(estimate_dns_dynamics(bare, sampling = "monthly"))
  expect_error(estimate_dns_dynamics(bare, sampling = "daily"),
               "beta1_change")
})

test_that("min_observations counts the estimation series, not the daily input", {
  # ~10 months of daily data: plenty of daily rows, too few month-ends
  df <- .sim_ou_params_df(n_days = 210)

  expect_error(estimate_dns_dynamics(df, sampling = "monthly"),
               "Not enough observations")
  expect_silent(estimate_dns_dynamics(df, sampling = "daily"))
})

test_that("daily mode reproduces the legacy default dt", {
  df <- .sim_ou_params_df(n_days = 504)

  implicit <- estimate_dns_dynamics(df, sampling = "daily")
  explicit <- estimate_dns_dynamics(df, dt = 1 / 252, sampling = "daily")
  expect_equal(implicit$beta1$kappa, explicit$beta1$kappa)
  expect_equal(implicit$dt, 1 / 252)
})
