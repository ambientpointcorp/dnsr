# Changelog

## dnsr 0.4.2

The vignettes, which had been left at 0.1.0-era content while everything
else moved, plus one defect they surfaced.

### Fix

- `estimate_dns_dynamics(structure = "full")` warned “No mean reversion
  detected” whenever a factor’s own diagonal entry came out
  non-positive. That is the right diagnostic for the diagonal structure,
  where a non-positive kappa does mean the factor never reverts, and the
  wrong one for a full matrix, where a negative own coefficient is
  perfectly compatible with a stable system: what settles that is the
  eigenvalues of the whole matrix. On the bundled panel the full VAR has
  `kappa[1,1] = -2.6` and eigenvalues of 1.40, 1.40 and 0.49 – stable,
  and warned about anyway. The per-factor warning is now raised only
  under `structure = "diagonal"`, and a system whose eigenvalues are not
  all positive raises one warning naming them.

### Vignettes

- `dns-model-introduction` ran almost nothing: 14 of its 16 chunks were
  `eval = FALSE`, and the two that ran fabricated yields with
  [`rnorm()`](https://rdrr.io/r/stats/Normal.html). It now runs end to
  end against the bundled `treasury_yields` panel, so the output is real
  and `R CMD check` verifies it. It covers the profiled fit,
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html), the standard
  errors and half-life interval, the full mean-reversion matrix, the
  exact simulation engine, the pricing helpers, and the model
  comparison.

  It also now explains the choice a panel forces and a single curve does
  not: with lambda free it ranges from 0.04 to 0.50 across these dates,
  so beta2 on one date does not mean what it means on the next and the
  factors stop being a comparable time series. Fixing lambda, as Diebold
  and Li do, costs about 14 basis points of median fit against 7 and
  buys factors that can be modelled over time.

- `dnsr-yuima-integration` claimed `"milstein"` gives “higher accuracy”.
  It does not: the Milstein correction is proportional to the derivative
  of the diffusion with respect to the state, and this diffusion is
  constant, so the term vanishes and the two schemes coincide. Removed,
  and the two engines now have a section explaining what each is
  actually for.

  Its 41-line hand-rolled Euler fallback for when YUIMA is absent is
  gone; `engine = "exact"` does that properly, with no discretisation
  error.

  The QMLE section printed a true-versus-estimated table with no
  comment, and the table is alarming on its face: kappa for the level
  factor comes back 120 times too large. The vignette now explains why,
  which makes it the most useful passage in it. Sigma is recovered to
  within a few percent, because volatility is identified by
  high-frequency variation. Kappa is not, and the error tracks the ratio
  of each factor’s half-life to the sample: 35 years against 2 gives
  120x, 2.3 years gives 7x, 1.4 years gives 1.3x. No estimator recovers
  a 35-year half-life from two years of data, and theta inherits the
  problem, since the level a process reverts to is identified only if
  the reversion is.

- `fred-api-setup` told readers to hand-roll a caching layer and to wait
  out rate limits. Both are built in now, so it documents `cache_dir`
  and `retries` instead – including that cache keys exclude the API key,
  that a hit is served before the key is looked up, and that the cache
  does not expire. It gains a section on missing observations covering
  `max_gap`, `fill_missing` and why the fill direction matters, a note
  on `date_order`, and a pointer to the bundled data for readers who do
  not need a key yet.

- Fixed a broken cross-reference: the introduction pointed at
  `vignette("yuima-integration")`, which does not exist.

## dnsr 0.4.1

Follow-ups to 0.4.0: no change to what any function computes, except one
documentation correction.

- [`fit_dns_model()`](https://ambientpointcorp.github.io/dnsr/reference/fit_dns_model.md)’s
  warm start makes a date’s fit depend on which *other* dates are in the
  panel, because each cross-section is initialised from the previous
  date’s lambda and the chain differs. Fitting a month-end series
  directly and fitting the daily series containing it give lambdas
  differing by up to 0.096 on the shared dates. That is inherent to warm
  starting rather than a defect, but it was undocumented, and
  `warm_start = FALSE` makes the fits panel-independent – verified
  identical, not merely close. Now documented, and pinned by a test.
- The
  [`estimate_dns_dynamics()`](https://ambientpointcorp.github.io/dnsr/reference/estimate_dns_dynamics.md)
  and
  [`dns_dynamics_to_sim_params()`](https://ambientpointcorp.github.io/dnsr/reference/dns_dynamics_to_sim_params.md)
  examples fitted 750 cross-sections to reach the 20 month-ends the
  default monthly sampling needs, taking about 9 seconds each and
  tripping the CRAN check note for examples over 5 seconds. They now fit
  only the month-ends the estimation actually uses, with
  `warm_start = FALSE` so the result matches what fitting every day
  would give. About nine times quicker.
- Offline coverage for the FRED data path, which the live tests skip and
  which was the least-covered file at 22.7 percent. Two routes, neither
  needing a key or a network: a prepared cache, which
  [`fetch_fred_series()`](https://ambientpointcorp.github.io/dnsr/reference/fetch_fred_series.md)
  consults before looking for a key, so the whole assembly in
  [`get_treasury_yield_data()`](https://ambientpointcorp.github.io/dnsr/reference/get_treasury_yield_data.md)
  runs for real; and `mockery` for the request and parse layer the cache
  skips. Covers the “.” to `NA` mapping, HTTP failure reporting, the
  empty-payload warning, the cache round trip including that no key
  reaches disk, decimal conversion, row ordering, gap filling with
  `max_gap`, and per-series fetch failures. `R/data-utils.r` goes from
  22.7 to 92.8 percent and the package from 76.2 to 82.4.
- Trimmed the test suite from 73 to 40 seconds without weakening an
  assertion: the exact simulation engine is exact at any step size, so
  the stationary test uses 120 steps rather than 720; the transition-law
  test needs fewer paths than it had, since its assertions are z-scores
  against the Monte Carlo standard error; the plan-restoration test uses
  a tweaked sequential plan, which is distinguishable from the plain one
  the bug installed and needs no workers; and the Svensson sweep uses a
  wider stride over the same range. The coverage job was taking 17
  minutes.
- New `pkgdown` workflow. `_pkgdown.yml` arrived in 0.4.0 but nothing
  rendered it: `docs/` was a gitignored local build, never committed,
  ten months stale, documenting three functions that no longer existed
  and missing six that did. The site is now built on a runner that has
  pandoc and published to `gh-pages`. Standing the workflow up
  immediately showed the index was incomplete – it omitted the
  `dnsr-package` topic, which pkgdown treats as an error – a gap that
  had gone unnoticed because nothing had ever rendered the config.

## dnsr 0.4.0

Enhancements, all reachable only through new arguments and new exports,
plus one correctness fix in the factor dynamics. Nothing that existed
before changes its name, arguments or argument order.

### Correctness

- [`estimate_dns_dynamics()`](https://ambientpointcorp.github.io/dnsr/reference/estimate_dns_dynamics.md)
  regressed the factor change on the **contemporaneous** level rather
  than the lagged level its own documented model specifies. Since
  beta(t) appears on both sides of beta(t) - beta(t-1), the ordinary
  least squares slope came back with the opposite sign: +(1 - phi)
  instead of -(1 - phi). The consequences were systematic rather than
  occasional:
  - kappa was computed negative for every factor on every call, and the
    `kappa <= 0` guard then rescued it with
    [`abs()`](https://rdrr.io/r/base/MathFun.html), which recovers
    roughly the right magnitude and so hid the problem.
  - because that guard always fired, theta **always** fell back to the
    sample mean, and the documented `-a/b` estimator never once
    executed.
  - `system_stability` was always TRUE, since
    [`abs()`](https://rdrr.io/r/base/MathFun.html) forced kappa
    positive.

  The regression is now on beta(t-1), derived as the level minus the
  change so it is exact under either sampling mode. kappa comes out
  positive on its own, theta comes from the intercept and slope, and a
  genuinely non-reverting factor now raises a warning naming the factor
  instead of being silently turned around. On the bundled Treasury panel
  this moves kappa by about 0.1 percent and leaves theta unchanged,
  because the two slopes are mirror images; it is decisive for the full
  mean-reversion matrix below, where there is no scalar sign for
  [`abs()`](https://rdrr.io/r/base/MathFun.html) to flip.

### New capabilities

- [`fit_ns_model()`](https://ambientpointcorp.github.io/dnsr/reference/fit_ns_model.md)
  and
  [`fit_dns_model()`](https://ambientpointcorp.github.io/dnsr/reference/fit_dns_model.md)
  gain `method`. `method = "profile"` concentrates the sum of squared
  errors: with lambda fixed the curve is linear in the betas, so they
  are solved exactly by least squares and only lambda is searched, over
  one dimension. It needs no start values, consults no beta bounds, and
  cannot return a clamped beta. On a synthetic Nelson-Siegel curve it
  recovers lambda to five decimals where the numerical search reaches
  three. `"optim"` remains the default.

- [`estimate_dns_dynamics()`](https://ambientpointcorp.github.io/dnsr/reference/estimate_dns_dynamics.md)
  gains `structure = "full"`, which regresses each factor change on all
  three lagged levels and returns the full mean-reversion matrix of the
  Diebold-Li VAR(1) in `kappa_matrix`, with long-run means solved from
  it in `theta_vector`.
  [`kappa_matrix_eigenvalues()`](https://ambientpointcorp.github.io/dnsr/reference/kappa_matrix_eigenvalues.md)
  now uses that matrix, so the stability verdict accounts for the cross
  terms. `"diagonal"` remains the default and still produces a diagonal
  `kappa_matrix`.

- Each factor’s dynamics now carry `se` (standard errors of kappa, theta
  and sigma), `t_value` and `half_life_ci`. kappa’s standard error
  rescales from the slope’s; theta’s comes from the delta method, since
  it is a ratio of coefficients. These describe sampling noise, not the
  upward small-sample bias in kappa, which is documented separately.

- New
  [`dns_discount_factors()`](https://ambientpointcorp.github.io/dnsr/reference/dns_discount_factors.md),
  [`dns_zero_prices()`](https://ambientpointcorp.github.io/dnsr/reference/dns_zero_prices.md)
  and
  [`dns_par_yields()`](https://ambientpointcorp.github.io/dnsr/reference/dns_par_yields.md)
  convert factors into the quantities a cashflow valuation consumes. Par
  yields are computed on the annuity factor and are defined only on
  whole numbers of coupon periods; shorter maturities return `NA` with a
  warning rather than being forced onto a convention they do not have.

- New `treasury_yields` dataset: 1,501 daily curves across eleven
  maturities, 2019 to 2024, from FRED’s H.15 Treasury constant maturity
  series. Six years spans 72 month-ends, enough for the default monthly
  dynamics estimation. Most examples now run against it instead of
  sitting inside `\dontrun{}`, and the test suite can cover the data
  path offline. Built by `data-raw/treasury_yields.R`.

- [`fetch_fred_series()`](https://ambientpointcorp.github.io/dnsr/reference/fetch_fred_series.md)
  and
  [`get_treasury_yield_data()`](https://ambientpointcorp.github.io/dnsr/reference/get_treasury_yield_data.md)
  gain `retries` and `cache_dir`. Requests are retried with exponential
  backoff on transient failures but not on a 4xx, where retrying cannot
  help. Cache keys cover the series and date range only: the API key is
  never part of the key and never written to disk, and a cache hit is
  served before the key is looked up, so a cached corpus can be replayed
  with no key configured.

- [`fit_ns_model()`](https://ambientpointcorp.github.io/dnsr/reference/fit_ns_model.md),
  [`fit_dns_model()`](https://ambientpointcorp.github.io/dnsr/reference/fit_dns_model.md),
  [`estimate_dns_dynamics()`](https://ambientpointcorp.github.io/dnsr/reference/estimate_dns_dynamics.md)
  and the two simulators now return classed objects with
  [`print()`](https://rdrr.io/r/base/print.html) methods, and
  [`fit_dns_model()`](https://ambientpointcorp.github.io/dnsr/reference/fit_dns_model.md)
  also has [`summary()`](https://rdrr.io/r/base/summary.html). These are
  classes layered onto the lists already returned: every element keeps
  its name, and [`is.list()`](https://rdrr.io/r/base/list.html), `$` and
  [`checkmate::assert_list()`](https://mllg.github.io/checkmate/reference/checkList.html)
  behave as before. Only the default display changes.

- The Nelson-Siegel-Svensson comparison in
  [`compare_term_structure_models()`](https://ambientpointcorp.github.io/dnsr/reference/compare_term_structure_models.md)
  could never have produced a column. It called
  `SmithWilson::Nelson.Siegel.Svensson()`, and `SmithWilson` is not a
  package: the one that exists is `SmithWilsonYieldCurve`. So
  [`requireNamespace("SmithWilson")`](https://rdrr.io/r/base/ns-load.html)
  returned FALSE on every machine and the branch always skipped with a
  warning. Correcting the name would not have helped either, because
  `SmithWilsonYieldCurve` exports only `fFitSmithWilsonYieldCurve()` and
  `fFitSmithWilsonYieldCurveToInstruments()` – no Svensson fitter.

  Svensson is now fitted directly. The curve is linear in its four betas
  once the two decays are fixed, so it uses the same concentrated least
  squares as `method = "profile"`, with the decays searched over a
  coarse grid and then refined. The comparison is therefore available to
  everyone with no optional dependency, both decays are estimated rather
  than fixed at 0.5 and 3, and on the bundled panel it fits better than
  Nelson-Siegel on every date, as a nested model must.

  The search is bounded because the profiled residual has three
  degenerate corners that fit beautifully and mean nothing: the two
  decays coinciding, the first approaching zero (the slope loading
  flattens to a constant and goes collinear with the intercept), and the
  first growing large (every loading but the level decays to nothing
  inside the shortest maturity). Unbounded, the search returned lambda1
  = 1e-4 with betas of plus and minus 16,800 on one real date. Bounded,
  no date in the panel produces a factor above 1 in magnitude.

- New `"SmithWilson"` option in
  [`compare_term_structure_models()`](https://ambientpointcorp.github.io/dnsr/reference/compare_term_structure_models.md),
  which is what `SmithWilsonYieldCurve` actually provides, with `ufr`
  and `alpha` arguments. It is opt-in rather than in the default
  `models`, and documented as an interpolation rather than a fit: it
  carries one kernel term per instrument, so it reproduces every
  observed yield to within 4e-12 basis points and its error metrics are
  zero by construction, not comparable with a three- or four-factor
  curve’s. The useful part is its extrapolation beyond the last
  maturity, which converges to `ufr`.

- Dropped `gganimate` and `gifski` from Suggests; nothing has referenced
  them since `animate_yield_curves()` was removed.

- [`plot_yield_pca()`](https://ambientpointcorp.github.io/dnsr/reference/plot_yield_pca.md)
  gains `scale`. `scale = FALSE` decomposes the covariance matrix rather
  than standardising each maturity first, which is the usual convention
  for yield curves and is what produces the familiar level, slope and
  curvature shapes. `TRUE` remains the default.

### Other fixes

- [`fit_dns_model()`](https://ambientpointcorp.github.io/dnsr/reference/fit_dns_model.md)
  reports parameters resting on a bound once for the whole panel, naming
  each parameter and the number of dates, instead of emitting one
  warning per date; a run over several hundred dates could previously
  produce hundreds. The per-date detail is in the new `at_bound`
  element, a logical matrix kept beside `parameters` rather than inside
  it, so that frame stays numeric apart from its `Date` column.
- The Greek letters in plot labels made several plots unrenderable on
  Windows. The default device there converts text through a single-byte
  charset and fails outright on U+03BB:
  `conversion failure on 'lambda = 0.0609' in 'mbcsToSbcs'`.
  [`plot_dns_loadings()`](https://ambientpointcorp.github.io/dnsr/reference/plot_dns_loadings.md)
  and
  [`plot_factor_time_series()`](https://ambientpointcorp.github.io/dnsr/reference/plot_factor_time_series.md)
  now use plotmath, which goes through R’s symbol font and renders
  anywhere, and the package source is fully ASCII. The bug predates this
  release; it surfaced only once the examples became runnable and
  `R CMD check` started drawing them. A test now prints every plot to a
  null device, since
  [`ggplot_build()`](https://ggplot2.tidyverse.org/reference/ggplot_build.html)
  stops short of text measurement and cannot catch this class of fault.
- `geom_line(size = )` replaced with `linewidth =`, which ggplot2 has
  wanted since 3.4.
- [`plot_yield_heatmap()`](https://ambientpointcorp.github.io/dnsr/reference/plot_yield_heatmap.md)
  decided whether to rescale yields to percent one column at a time, so
  a panel whose short end was below 1 and long end above it could be
  drawn in two different units. The decision is now made once from the
  whole frame. `plot_type = "curvature"` with fewer than three
  maturities counted backwards through its loop and now reports the
  requirement.
- [`dns_to_yuima_model()`](https://ambientpointcorp.github.io/dnsr/reference/dns_to_yuima_model.md)
  tested the correlation diagonal for exact equality with 1, rejecting
  matrices that are correct to within rounding.
- The package-level example simulated from `dynamics$initial_factors`,
  which
  [`estimate_dns_dynamics()`](https://ambientpointcorp.github.io/dnsr/reference/estimate_dns_dynamics.md)
  has never returned, so it errored. It now starts from the most recent
  fitted factors.
- The internal Nelson-Siegel loading formula had four copies; it is now
  defined once. The arithmetic is unchanged, so yields and forwards
  remain identical to the last bit.
- Dropped an unused local in
  [`plot_simulated_factors()`](https://ambientpointcorp.github.io/dnsr/reference/plot_simulated_factors.md).

## dnsr 0.3.0

Correctness, performance and packaging release. Every exported function
keeps its name, arguments and argument order; new behaviour is reachable
only through new trailing arguments whose defaults reproduce 0.2.0. Four
returned values change, each because the old value was wrong — see
*Values that change* below.

### Critical fixes

- [`get_treasury_yield_data()`](https://ambientpointcorp.github.io/dnsr/reference/get_treasury_yield_data.md)
  filled missing yields **backwards in calendar time**. The frame is
  sorted most-recent-first and the fill ran in row order, so each gap
  took the value of the next *later* date — look-ahead bias, on by
  default, in every panel with holes (DGS20 before 1993, DGS30 between
  2002 and 2006). Gap filling now runs on an ascending date index. New
  `max_gap` caps how many calendar days a value may be carried, and new
  `date_order` chooses the row order of the result.
- [`fit_dns_model()`](https://ambientpointcorp.github.io/dnsr/reference/fit_dns_model.md)
  selected the **oldest** `n_days` rows, not the most recent. It took
  the last rows of the frame, which is only “most recent” for ascending
  input, while this package’s own data functions return descending. The
  window is now chosen by date and does not depend on row order. New
  `date_order` controls the output order and defaults to `"input"`, so
  `parameters[1, ]` remains the most recent row for a most-recent-first
  frame.
- Cross-sections that could not be fitted **survived as zero-factor
  rows**. `params_df` was pre-allocated with
  [`numeric()`](https://rdrr.io/r/base/numeric.html), so failed dates
  kept `beta1 = beta2 = beta3 = rmse = 0` and the `!is.na(beta1)` filter
  meant to drop them never matched. Failed dates are now `NA` and are
  removed; new `drop_failed = FALSE` keeps them as `NA`, and the result
  gains a `failures` data frame naming each date and reason.
- A **single missing maturity on a single date aborted the whole fit**
  with
  `number of items to replace is not a multiple of replacement length`.
  [`fit_ns_model()`](https://ambientpointcorp.github.io/dnsr/reference/fit_ns_model.md)
  returns fitted values only for the maturities it retained; these are
  now placed back into a full-width row by index, leaving absent
  maturities as `NA`.
  [`fit_ns_model()`](https://ambientpointcorp.github.io/dnsr/reference/fit_ns_model.md)
  gains `valid_idx` for this purpose.
- `beta*_change` columns were differenced **in row order rather than
  calendar order**, so for descending input every change had the wrong
  sign and the change-on-level slope — and therefore `kappa` — flipped.
  Changes are now first differences in calendar order, written back in
  the output’s row order.
- `fit_dns_model(parallel = FALSE)` **permanently reset the caller’s
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)**
  to sequential, tearing down a multisession pool set up for the rest of
  the session. No plan is now changed as a side effect, and the
  `parallel = TRUE` path restores the previous plan on exit.
- The lambda warm start was carried in `last_fit_result <<- result`,
  which **wrote into the global environment** and could not work under
  `future_lapply()` in any case. Warm starting is now explicit, local,
  and walks forward in calendar time; new `warm_start = FALSE` disables
  it.

### Other fixes

- [`compare_term_structure_models()`](https://ambientpointcorp.github.io/dnsr/reference/compare_term_structure_models.md)
  scaled its error metrics by 100 while labelling them basis points —
  every reported RMSE, MAE and MaxError was 100x too small. Decimal
  yields now scale by 10,000, matching
  [`fit_dns_model()`](https://ambientpointcorp.github.io/dnsr/reference/fit_dns_model.md).
  New `units = "decimal"` returns them in the units of the input. The
  function also no longer breaks on a maturity with no observation on
  the chosen date.
- [`simulate_dns_yield_curves()`](https://ambientpointcorp.github.io/dnsr/reference/simulate_dns_yield_curves.md)
  validated `sim_method` and then discarded it, so `"milstein"` silently
  ran Euler.
  [`simulate_dns_factors()`](https://ambientpointcorp.github.io/dnsr/reference/simulate_dns_factors.md)
  gains `sim_method` and the value is forwarded and recorded on the
  result. (For this model the diffusion is state-independent, so the two
  schemes coincide; the documentation now says so.)
- [`plot_yield_curves()`](https://ambientpointcorp.github.io/dnsr/reference/plot_yield_curves.md)
  labelled facets with integer date codes (“19369”):
  [`intersect()`](https://rdrr.io/r/base/sets.html) and
  `for (date in dates)` both drop the `Date` class. Facets are now dated
  and ordered by calendar. RMSE annotations moved to their own layer
  instead of being smuggled through the colour scale as a fake series.
- [`plot_simulated_factors()`](https://ambientpointcorp.github.io/dnsr/reference/plot_simulated_factors.md)
  failed on a single-path simulation
  (`dim(X) must have a positive length`), as did the term-structure
  branch of
  [`plot_simulated_yields()`](https://ambientpointcorp.github.io/dnsr/reference/plot_simulated_yields.md).
  Both now keep the path dimension.
- `plot_yield_heatmap(reverse_time = TRUE)` was a no-op, and its
  documentation described the *default* orientation rather than the
  effect of setting it: a date axis increases upward, so the most recent
  date was already at the top. The argument now flips the axis, putting
  the most recent date at the bottom, and says so. It is implemented on
  the data rather than through the scale or the coordinate system,
  because neither is portable –
  [`scale_y_reverse()`](https://ggplot2.tidyverse.org/reference/scale_continuous.html)
  negates the variable and unary minus is undefined for `Date`, while
  `coord_cartesian(ylim = rev(range))` reverses the panel on ggplot2 3.x
  but is normalised away on 4.x.
- [`plot_yield_pca()`](https://ambientpointcorp.github.io/dnsr/reference/plot_yield_pca.md)
  errored with `object 'var_explained' not found` when all three
  `plot_*` flags were `FALSE`.
- [`plot_simulated_yields()`](https://ambientpointcorp.github.io/dnsr/reference/plot_simulated_yields.md)
  facets sorted lexically, putting “Year 10” between “Year 1” and “Year
  5” for the default `time_points`.
- [`dns_dynamics_to_sim_params()`](https://ambientpointcorp.github.io/dnsr/reference/dns_dynamics_to_sim_params.md)
  always returned `lambda = 0.0609`; the branch meant to override it
  could never fire, because
  [`estimate_dns_dynamics()`](https://ambientpointcorp.github.io/dnsr/reference/estimate_dns_dynamics.md)
  does not return a lambda. New `lambda` argument, and
  [`get_calibrated_dns_params()`](https://ambientpointcorp.github.io/dnsr/reference/get_calibrated_dns_params.md)
  now passes the lambda it actually fitted.
- [`fit_ns_model()`](https://ambientpointcorp.github.io/dnsr/reference/fit_ns_model.md)’s
  default bounds assumed decimal yields and clamped silently. L-BFGS-B
  reports `convergence == 0` at a bound, so a clamped fit was returned
  with `success = TRUE`: percentage-point yields pinned all four
  parameters at their bounds with an RMSE of 2.7 percentage points.
  Default beta bounds now widen to cover the data scale, start values
  are kept inside the box, and any parameter resting on a bound is
  reported in the new `at_bound` field and raises a warning. An
  explicitly supplied `optimizer` is now honoured by
  `lambda_constraint = "free"`, which previously always used BFGS. The
  slope start value is taken by maturity rather than by position.
- [`estimate_dns_parameters_yuima()`](https://ambientpointcorp.github.io/dnsr/reference/estimate_dns_parameters_yuima.md)
  took factor columns positionally, so passing the `parameters` frame
  the documentation describes fed `Date` into
  [`diff()`](https://rdrr.io/r/base/diff.html). Columns named
  `beta1`/`beta2`/`beta3` are now used when present.
- [`get_calibrated_dns_params()`](https://ambientpointcorp.github.io/dnsr/reference/get_calibrated_dns_params.md)
  gains `api_key`; the calibration path was previously reachable only
  through the environment variable. Its “latest” factors and lambda are
  now selected by date rather than by row position.
- `tests/testthat.R` was missing, so `R CMD check` skipped the entire
  test suite. Added, along with `mockery` in Suggests and a regression
  test for each fix above.

### Performance

- [`calculate_dns_yields()`](https://ambientpointcorp.github.io/dnsr/reference/calculate_dns_yields.md)
  and
  [`calculate_dns_forwards()`](https://ambientpointcorp.github.io/dnsr/reference/calculate_dns_forwards.md)
  were triple loops over scalars. Only the maturity dimension is looped
  now; time and path are handled by vectorised arithmetic, with the
  expression and its left-to-right grouping unchanged, so results are
  identical to the last bit. On 601 steps x 500 paths x 11 maturities:
  yields 24.0s to 0.079s, forwards 2.33s to 0.089s. (The gap between the
  two was `loadings[m_idx, "loading_2"]` – indexing a matrix by column
  name inside the innermost loop.)
- New `engine` argument on
  [`simulate_dns_factors()`](https://ambientpointcorp.github.io/dnsr/reference/simulate_dns_factors.md)
  and
  [`simulate_dns_yield_curves()`](https://ambientpointcorp.github.io/dnsr/reference/simulate_dns_yield_curves.md).
  `engine = "exact"` samples the factor system’s exact Gaussian
  transition law instead of calling
  [`yuima::simulate()`](https://rdrr.io/pkg/yuima/man/simulate.html)
  once per path: no discretisation error at any step size, no YUIMA
  dependency, and about 120x faster (100 paths x 120 steps: 1.78s to
  0.015s). The step covariance is built as
  `rho_ij sigma_i sigma_j (1 - exp(-(kappa_i + kappa_j) dt)) / (kappa_i + kappa_j)`,
  which couples the two reversion speeds – using each factor’s own
  variance with a Cholesky factor of the correlation matrix is only
  correct when every kappa is equal. The default remains `"yuima"`; the
  two engines agree in distribution but do not produce the same paths
  from the same seed. The default is expected to change to `"exact"` in
  a future release.
- [`fit_dns_model()`](https://ambientpointcorp.github.io/dnsr/reference/fit_dns_model.md)
  no longer searches for each date’s row with `which(Date == date)`,
  which made the fit quadratic in the number of dates.
- The
  [`plot_simulated_factors()`](https://ambientpointcorp.github.io/dnsr/reference/plot_simulated_factors.md)
  and
  [`plot_simulated_yields()`](https://ambientpointcorp.github.io/dnsr/reference/plot_simulated_yields.md)
  frames are built in one pass instead of by
  [`rbind()`](https://rdrr.io/r/base/cbind.html) inside nested loops.
  Row order and values are unchanged.

### Package health

`R CMD check` went from 4 warnings and 4 notes to zero of either, with
the test suite running for the first time. Remaining output on a given
machine is limited to uninstalled `Suggests` packages and, where pandoc
is absent, unbuilt vignettes.

- `R/dnsr-package.R` was at the repository root as `dnsr-package.r`, so
  its roxygen block never ran and the file was never sourced: the
  package-level help page and its `@import` directives were both dead.
  Moved into `R/`, and `?dnsr-package` now exists. `dnsr-structure.txt`
  and `dnsr-package-summary.Rmd` are no longer shipped in the tarball.
- Non-ASCII characters removed from `R/`: the Greek letters and
  subscripts in plot labels are now `\uxxxx` escapes, which renders
  identically, and the mathematical notation in comments is plain ASCII.
- `LICENSE` is now the two-field DCF stub that `MIT + file LICENSE`
  requires; the full text moved to `LICENSE.note`, which still ships.
- Added the missing `importFrom(stats, ...)` entries for `cor`, `median`
  and `residuals`, and a
  [`utils::globalVariables()`](https://rdrr.io/r/utils/globalVariables.html)
  registration for the columns created inside functions and then
  referenced in ggplot2 aesthetics. The “no visible binding” and “no
  visible global function definition” notes – 40 of them – are gone.
- Dropped `readr`, `data.table`, `Matrix`, `tidyr` and `lubridate` from
  Imports: nothing in the package referenced them. The dplyr verb
  imports were also unused; only the pipe is kept, for
  [`plot_yield_surface()`](https://ambientpointcorp.github.io/dnsr/reference/plot_yield_surface.md).
- `[time, path, factor]` in `@return` and `@param` text was being parsed
  as an Rd link, producing cross-reference warnings. Rewritten as prose.
- The YUIMA vignette could not be built on a machine without YUIMA:
  `has_yuima` was defined inside a chunk gated on
  `eval = requireNamespace( "yuima")`, so when that was FALSE the chunk
  never ran and the fallback chunks testing `eval = !has_yuima` died
  with `object 'has_yuima' not found`. It is now set in the always-run
  setup chunk.
- `DESCRIPTION` gains `URL`, `BugReports`, `Language` and
  `Depends: R (>= 4.1)`.
- Added `_pkgdown.yml`, giving the reference index an organised
  structure rather than one flat alphabetical list, and GitHub Actions
  workflows for `R CMD check` (Linux/macOS/Windows, release and oldrel)
  and test coverage. There was no CI before, which is how a missing
  `tests/testthat.R` went unnoticed.

### Values that change

- **Fill direction.** Gaps are filled from the last earlier observation
  rather than the next later one. Opt out with `fill_missing = FALSE`.
- **Change sign.** `beta*_change` is always `beta(t) - beta(t-1)`. This
  reaches only callers using
  `estimate_dns_dynamics(sampling = "daily")`; the monthly default
  already recomputed the changes on a sorted grid.
- **Error units.** `compare_term_structure_models()$metrics` is now
  genuinely in basis points, so the values are 100x their 0.2.0
  counterparts. Pass `units = "decimal"` for the input’s own units.
- **Simulation lambda.** `get_calibrated_dns_params()$params$lambda` is
  the fitted lambda rather than a hard-coded 0.0609.
  [`dns_dynamics_to_sim_params()`](https://ambientpointcorp.github.io/dnsr/reference/dns_dynamics_to_sim_params.md)
  called directly, with no `lambda` argument, still returns 0.0609.

Additive changes to returned objects:
[`fit_ns_model()`](https://ambientpointcorp.github.io/dnsr/reference/fit_ns_model.md)
gains `valid_idx` and `at_bound`;
[`fit_dns_model()`](https://ambientpointcorp.github.io/dnsr/reference/fit_dns_model.md)
gains `failures`;
[`simulate_dns_factors()`](https://ambientpointcorp.github.io/dnsr/reference/simulate_dns_factors.md)
and
[`simulate_dns_yield_curves()`](https://ambientpointcorp.github.io/dnsr/reference/simulate_dns_yield_curves.md)
gain `sim_method` and `engine`. Code reading these lists by name is
unaffected; code asserting an exact set of element names is not.

## dnsr 0.2.0

### `estimate_dns_dynamics()` owns the sampling policy

- New `sampling` argument: `"monthly"` (the default) subsamples the
  factor series to its last observation per month, recomputes the change
  columns on the month-end grid, and estimates at `dt = 1/12`, per
  Diebold-Li practice. Daily cross-sectional refits carry serially
  uncorrelated measurement noise that the OU change-on-level regression
  reads as mean reversion, inflating kappa by an order of magnitude;
  month-end sampling averages the noise down and recovers reversion
  speeds on the economically meaningful time scale. Every caller now
  gets this policy from one place instead of replicating it.
- `sampling = "daily"` estimates on the series as given and reproduces
  the 0.1.0 behavior exactly (including the `dt = 1/252` default).
- `dt` now defaults to the sampling mode’s spacing. In monthly mode the
  step is month-end to month-end by construction, so a conflicting
  explicit `dt` is ignored with a warning.
- Monthly subsampling is idempotent: a series that is already month-end
  passes through unchanged, so callers that pre-subsample (as
  [`get_calibrated_dns_params()`](https://ambientpointcorp.github.io/dnsr/reference/get_calibrated_dns_params.md)
  did in 0.1.x) get identical estimates.
- In monthly mode the `beta*_change` input columns are no longer
  required; they are derived internally.
- The returned dynamics gain `$innovation_dates` (the date axis of the
  `$innovations` rows, so consumers can align innovations across models
  estimated on different calendars), plus `$sampling` and `$dt`.
- [`get_calibrated_dns_params()`](https://ambientpointcorp.github.io/dnsr/reference/get_calibrated_dns_params.md)
  now delegates the subsampling to
  [`estimate_dns_dynamics()`](https://ambientpointcorp.github.io/dnsr/reference/estimate_dns_dynamics.md);
  its results are unchanged.

## dnsr 0.1.0

- Initial release: Treasury data fetch, DNS cross-sectional fitting, OU
  factor dynamics, simulation, and visualization.
