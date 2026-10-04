# dnsr 0.2.0

## `estimate_dns_dynamics()` owns the sampling policy

* New `sampling` argument: `"monthly"` (the default) subsamples the factor
  series to its last observation per month, recomputes the change columns on
  the month-end grid, and estimates at `dt = 1/12`, per Diebold-Li practice.
  Daily cross-sectional refits carry serially uncorrelated measurement noise
  that the OU change-on-level regression reads as mean reversion, inflating
  kappa by an order of magnitude; month-end sampling averages the noise down
  and recovers reversion speeds on the economically meaningful time scale.
  Every caller now gets this policy from one place instead of replicating it.
* `sampling = "daily"` estimates on the series as given and reproduces the
  0.1.0 behavior exactly (including the `dt = 1/252` default).
* `dt` now defaults to the sampling mode's spacing. In monthly mode the step
  is month-end to month-end by construction, so a conflicting explicit `dt`
  is ignored with a warning.
* Monthly subsampling is idempotent: a series that is already month-end
  passes through unchanged, so callers that pre-subsample (as
  `get_calibrated_dns_params()` did in 0.1.x) get identical estimates.
* In monthly mode the `beta*_change` input columns are no longer required;
  they are derived internally.
* The returned dynamics gain `$innovation_dates` (the date axis of the
  `$innovations` rows, so consumers can align innovations across models
  estimated on different calendars), plus `$sampling` and `$dt`.
* `get_calibrated_dns_params()` now delegates the subsampling to
  `estimate_dns_dynamics()`; its results are unchanged.

# dnsr 0.1.0

* Initial release: Treasury data fetch, DNS cross-sectional fitting, OU
  factor dynamics, simulation, and visualization.
