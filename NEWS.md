# News

_This NEWS file is best viewed on our [website](https://grantmcdermott.com/dbreg/NEWS.html)._

## dbreg 0.0.3.9903

- Use bounded FE moments for single-intercept-FE GLMs on DuckDB, reusing the diagonal projection and coefficient-state IRLS path. Preserve native Databricks, varying-slope, forced-SQL and capacity fallbacks.

## dbreg 0.0.3.9901

Local research release.

- Raw compressed predictions use saved column identities, and compound fixed-effect keys preserve native component precision.

- Keep fitted factor dictionaries distinct from display labels when rebuilding OLS and IV predictions.
- Drop compressed regressors absorbed by fixed effects before counting parameters for covariance corrections.
- Preserve raw values when a source supplies regressor, fixed-effect, weight or offset roles, including nearby numeric FE levels.
- Restore compressed composite-FE prediction and support nonstandard FE names.
- Include fitted levels and their joint covariance in constrained spline plots with fixed effects.

- Retained clustering accepts nonstandard column names, and source names
  containing `_x_` or `:` remain literal in sample and fixed-effect metadata.
- GLM fitting and R/SQL prediction handle nonstandard source names consistently,
  including retained samples and later clustering on original columns.
- Learned fixed-effect predictions preserve distinct columns with the same
  displayed coefficient name, including dropped-column estimability checks.
- GLMs report nonfinite initial information as an error instead of silently
  dropping an overflowed predictor as collinear.
- Compressed linear fits handle an ordinary categorical regressor with one
  observed level without changing the compressed grouping keys.
- IV fits keep structural and instrument coordinates separate. R prediction
  reconstructs structural terms from raw columns when public labels collide.
- SQL Server and MySQL `i()` coding reads native column metadata through scoped
  adapters; other backend metadata paths retain their existing behavior.
- `i()` level discovery excludes missing keys and orders integer keys by their
  native values, including when the source column is named `value`.
- HC1 score queries project residual and regressor expressions once, reducing
  repeated SQL while preserving the covariance multiplication order.
- Coefficient labels come from fitted source and basis identities. Literal
  `_x_` text in names or factor levels stays unchanged; products use `:`.
- A categorical level named `_constant` remains distinct from the compiler's
  synthetic constant column.
- Generated SQL design columns avoid outcome, regressor, weight, offset and
  cluster names, including across IV blocks. Distinct interaction and `i()`
  columns retain separate coordinates when their initial aliases collide.
- SQL builders emit backend-specific text casts and column exclusions directly;
  complete statements no longer pass through a global token rewriter.
- Spark AP and fitted-table creation validate physical-table capabilities before
  acquiring resources. Read-only fitting remains available without native mode.
- Fixed-effect coefficient summaries and covariance corrections preserve distinct
  fitted columns with the same public label, including a factor dummy and a raw
  regressor whose names coincide.
- `dbglm(null_model = TRUE)` optionally stores global null likelihood diagnostics
  using the final estimation sample, weights and offset. The default performs no
  null-model queries; printed and extracted diagnostics use saved scalar results.

- FE binscatter keeps explicit level coefficients and covariance under automatic
  strategy selection. Generated bin and spline names avoid source columns,
  including after restoring public coefficient names.
- GLM prediction reconstructs the fitted formula even when newdata contains
  columns matching internal SQL aliases. OLS prediction intervals honor the
  covariance's testing degrees of freedom while retaining residual degrees of
  freedom for prediction-noise variance.
- Native validation uses one checked suite inventory shared by tests, its runner
  and CI; `--list-suites` validates it without connecting to a database.

- Direct SQLite and PostgreSQL connections are no longer supported; their
  dialect branches, optional dependency and backend tests have been removed.

- `dbglm()` supports automatic or explicit exact cell compression before IRLS,
  alongside original-row moments and SQL AP strategies. Compression preserves
  fractional-weight likelihoods, original observation counts, IID/HC1/cluster
  covariance and requested prediction/inference retention across all supported
  families. Retained original-row scores reuse the compact final projector.
- All GLM families normalize very small or large weights internally, restoring
  likelihoods and IID inference to their original units. This prevents false
  nonconvergence with uniformly large binomial or Poisson weights. A shared
  raw-key preflight avoids recursive pruning when no FE group needs removal.

- DuckDB linear fixed-effect fits compute within moments from compact group
  information using FWL, without writing observation-sized demeaning tables.
  Robust covariance reads projected rows through an inline join to compact FE
  coefficients. Near-absorbed columns use a read-only precision check. Requested
  inference retention still saves the original sample for later covariance changes.
- Unweighted DuckDB moments and compressed cells use stable covariance aggregates
  in a single source pass. Demeaned frames reuse known counts and original scales;
  tiny compressed designs use bounded ordinary model matrices.

- Fixed-effect systems now keep cross blocks sparse and use component-aware
  Cholesky for sparse two-way Schur systems, reusing the symbolic factorization
  across IRLS updates. Small systems retain dense factorization. Component null
  projections avoid observation-by-component dense allocations in compact fits;
  numerical zeros rebuild component pins, and direct row scaling handles very
  small information weights without reciprocal overflow. The iterative fallback
  deflates connected-component null directions and uses bounded Schur SGS or
  block SGS preconditioning. It centers source columns before solving while
  preserving the original moment reductions used by GLM inference.

- Linear moments now center variables before forming products, normalize weights,
  and scale rank/solve calculations. Robust score passes retain centered
  coordinates through covariance calculation. Coefficients and covariance return
  in the original units; SQL integer products, interactions and powers are promoted
  before multiplication in both linear and GLM designs.
- Native clustered tables and intervals now use the same testing degrees of
  freedom as post-fit inference. Saturated fits and one-cluster covariance are
  rejected explicitly. DuckDB/Spark sample filtering treats numeric NaNs as missing.
- Compressed factor designs stay sparse during centering; an allocation-bound
  regression test protects against dense cells-by-FE matrices.
- Iterative FE projection checks group means relative to projected weighted RMS as well
  as in absolute units. FE-only combinations decay below the original-value
  roundoff bound before rank detection, preventing spurious identified columns.
- Cluster-score cross-products stay in SQL for ordinary and compressed OLS/IV.
  Unconstrained binscatter builds its basis in SQL and collects only bin geometry;
  sampled-table cleanup is installed immediately.
- Unconstrained binscatter uses automatic strategy selection without fixed
  effects, and compression for FE level plots to retain explicit coefficients
  and joint covariance.
  Plot intervals inherit the fitted testing degrees of freedom, and simultaneous bands use the
  covariance of the plotted contrasts with column-wise standardization. Bands
  also work with pointwise intervals disabled.
- Automatic no-FE fitting skips compression probes for known continuous R input.
  Exact prediction keys are cast before collection, avoiding a second full fetch.
  Added numerical-invariance and DBI-transfer regression tests.

- DuckDB/Spark SQL weights now reject infinities during input validation for
  OLS/IV, before they reach numerical factorization.

- Compression now applies the same declared full-FE degrees-of-freedom policy
  as within estimation on disconnected/nested graphs; numerical dummy rank no
  longer silently changes covariance corrections across strategies.

- Independent references are protected against compatibility-export masking.
- FE projection, backend SQL and resource helpers are shared across estimators.
  Compact marginal/pair/global transport avoids mostly-NULL transfers and needs no
  optional DuckDB JSON extension. Prediction state and recovery transfer have
  separate bounds. Automatic two-way sparse fallback on DuckDB is memory planned.
- Projection outcomes and covariance repair expose diagnostics. Cleanup forgets
  successful drops and keeps failed drops retryable. Retained samples preserve immutable row identity.
- Local Spark 4.2.0 has a small real-SQL syntax/algorithm tier, separate from
  native Databricks ODBC validation and performance benchmarks.

- R prediction explicitly sums potentially overflowing rows, avoiding BLAS- and
  batch-dependent `Inf - Inf` results. Ordinary predictions retain the matrix
  product; real infinite Poisson means and its lower bound remain supported.

- Added supported additive `i()` indicators/numeric group slopes with fitted
  references and literal `ref`/`keep` controls. Bounded typed dictionaries drive
  R/SQL prediction, preserving large integer and adjacent-double categories.
  Malformed/unseen keys and ambiguous public names fail explicitly. Linear models retain native OLS classes, goodness-of-fit statistics and
  moments/compress/demean strategies, including intercept-only, factor and
  varying-slope formulas.
- Added a support/difference table and executable R/SQL lifecycle examples.

- Shared SQL Boolean literals now use `TRUE`/`FALSE`, avoiding native ODBC
  numeric quoting in subset/cluster predicates. Native estimators share string/list SSC policies with retained inference.

- Removed the fixest-shaped estimator facade and colliding helper exports.
  `dbreg()`/`dbglm()` keep `retain="none"` defaults. Delayed S3 methods support
  fixest coefficient summaries, fixed effects and coefficient plots, generics
  augmentation, and insight/modelsummary covariance requests. No fixest global
  helper is changed. Use `case.names()` for sample IDs.
- `tidy()` accepts supported covariance requests or a supplied matrix, including
  the matrix passed by modelsummary, with matching standard errors and intervals.
- `dbglm()` defaults to Gaussian identity and accepts family strings, functions
  and family objects. Formula weights/offsets/subsets and lazy-table `data=` work
  through shared input handling; IV also accepts offsets and subsets.
- Caller-owned GC cleanup queues retired relations for synchronous cleanup;
  private retained DuckDB sessions use one thread. `close()` releases owned state.

- `retain="inference"` saves fitted rows/scores without recovering learned FE
  dictionaries. It supports saved-sample operations and post-fit inference.

- `dbglm(family=gaussian())` reuses the shared FE/varying-slope projectors with
  one certified WLS update. Gaussian dispersion, t inference and likelihood
  follow `fixest::feglm()`. OLS models continue to use `dbreg()`; binomial fits specify their family.
- Post-fit OLS IID covariance with `ssc=list(K.adj=FALSE)` now uses the reference
  `N-1` denominator. The default adjusted covariance is unchanged.

- `dbreg(offset="column")` supports fixed-coefficient numeric offsets through
  the existing strategies on DuckDB/native Databricks. Saved-sample and newdata
  prediction add offsets once; compact offset models require explicit newdata.

- Explicit `retain="all"` saves an immutable database sample and projected scores
  for default/original-sample prediction, fitted values, residuals and `stats::case.names()`.
  SQL inputs return keyed relations; R inputs preserve original positions.
- `vcov()`, `summary()` and `confint()` can recompute IID, HC1, one-/two-way OLS
  and one-way GLM covariance on the retained sample, with supported SSC controls.
  Weighted OLS likelihood matches fixest's unweighted residual-RSS convention.
  Compact native defaults remain unchanged; full retention uses row-sized SQL
  storage and requires its owning session for sample-dependent operations.

- `predict()` can return keyed database output on DuckDB/native Databricks using
  `row_id` and `conn`/`table`/`path` or `newdata`. OLS and all four GLM families
  use saved factor coding and learned effects. `head()` previews output and
  `as.data.frame()` explicitly collects it. Output leases survive model close;
  `close(output)` releases owned tables and preserves caller-owned exports.
- SQL factor literals preserve 17-digit numeric values and native text containing
  quotes, backslashes or UTF-8. Prediction uses original typed FE keys and keeps
  exact large output keys when a driver would otherwise round them.

- Oversized one- or two-FE prediction dictionaries can stay in SQL, with bounded
  key matching and estimability checks. `close(model)` releases owned state;
  caller connections/exports survive. These live dictionaries cannot be restored
  by serializing their handles; small portable dictionaries remain supported.
  Two-FE recovery checks connected components and regressor relations before
  predicting new combinations. FE-only GLMs support zero-column designs.

- `dbglm()` accepts nonnegative fractional observation weights for all supported
  families. Fitting and inference use the original weights; `nobs()` counts
  retained observations. For fractional binomial weights, `logLik()` matches
  fixest/R's rounded-weight reporting convention. Use
  `logLik(model, type = "weighted")` for the unrounded objective.
- Estimator help pages and the introductory vignette describe supported
  prediction, formula, inference and fitted-state requests.
- `dbreg()` rejects unused arguments.
  Coefficient, covariance, confidence-interval and prediction methods reject
  unsupported requests rather than silently ignoring them. Existing `query_only`
  and `fes` aliases remain supported.
- Saved linear-model formulas and compressed prediction formulas no longer
  capture the caller's data or fit-time connection environment. Existing fitted
  coding and predictions are preserved when serializing compact objects.
- OLS now supports `row_id` and `sample_table` for retained-key exports on
  DuckDB and native Databricks. Requested exports fit from an owned snapshot,
  preserve typed composite/integer64 keys and release temporary fit resources.
  OLS and GLM share key validation and export helpers. Databricks exports retain
  nonstandard key names through scoped column mapping and recognize existing
  session tables before fitting.
- `dbglm(retain = "prediction")` saves bounded portable learned-effect state
  for link/response prediction without outcomes or fitting weights. Prediction
  includes varying slopes and checks the identified joint regressor/FE span;
  unknown keys and non-estimable combinations return `NA`. Native retention
  remains off by default. Database-backed large state and OLS support are pending.
- Absorbed GLMs also certify group equations relative to group information,
  preventing premature convergence for extremely small group-specific weights.
- `dbreg(retain = "prediction")` supports bounded portable learned-effect
  predictions for demean/compress fits, including typed composite keys and
  estimability checks. Newdata need not include outcomes or fitting weights.
  Existing stored-cell predictions for compression remain available. Mundlak
  retention and intervals using learned effects remain unsupported.
- Compression preserves exact large integer FE labels and adjacent double
  levels. Cluster scores match typed cell tuples, including embedded delimiters,
  and preserve precise cluster IDs. Only exceptional lossy integer results need
  another SQL fetch; caller connection conversion settings are unchanged.

### Additional linear-model improvements

- OLS supports restricted square/cube terms, explicit factor references and
  composite fixed effects. Predictions retain factor coding and missing-row
  positions, including integer64 inputs; weighted multiway demeaning is fixed.
- Covariance extraction matches visible coefficients. Clustered compression
  shares its expanded design with the score query, and generated SQL aliases
  avoid existing column names.
- OLS/IV gain optional sample audits, per-fit projection controls, and executable
  SQL inspection plans. These additions build on PR80 without changing its head.

### Upstream development changes

New features

- Alternating projections on Spark/Databricks now use sequential SQL window
  functions (both unweighted and weighted) and batch three cycles per temporary
  table write, controlled by the new `dbreg.ap_method` and
  `dbreg.ap_batch_size` options.
- New `?dbreg-options` help page documenting all global options.
- Added `dbivreg()` for instrumental-variables (2SLS) regression on database
  backends. Uses `fixest`-style IV formulas, e.g. `y ~ x | d ~ z` or
  `y ~ x | fe1 + fe2 | d ~ z`, and supports multiple endogenous regressors,
  factors/interactions, weights, and iid, HC1, or clustered SEs (incl.
  `ssc`). Estimation is exact via the `"moments"`, `"demean"` (any number of
  FE), and `"compress"` (no FE) strategies. (#66)
- Added support for weighted regression via the new `weights` argument. Accepts
  a column name string and is compatible with all strategies and SE types (iid,
  HC1, clustered). (#12)
- The `"demean"` strategy now supports three or more fixed effects via
  alternating projections (AP). Previously this was limited to at most two FE.
  The `"auto"` heuristic will select `"demean"` (AP) for 3+ FE when compression
  is not viable. (#12)
- For unbalanced two-way panels, `strategy = "auto"` now routes to `"demean"`
  via alternating projections instead of erroring. Users no longer need to
  explicitly choose between `"compress"` and `"mundlak"` for these cases.
- Formulas now support `fixest`-style `i()` terms, e.g.
  `dbreg(y ~ i(period, treat, ref = 5) | id + period, ...)` for event studies,
  in `dbreg()` and `dbivreg()` and across all strategies. Supports `i(f)`,
  `i(f, x)` with numeric or logical `x`, and the `ref` and `keep` arguments.
  Coefficients are named as in `fixest` (e.g. `period::1:treat`). (#82)

Bug fixes

- Alternating projections no longer fail with an SQL error when a fixed-effect
  or cluster variable also appears as a regressor. Invalid `dbreg.ap_max_iter`
  and `dbreg.ap_tol` options now produce a clear error.
- Fixed a typo in `<dbreg.object>$coeftable` column names that was obscuring
  some downstream functions, e.g. stars in `modelsummary::modelsummary`. (#62) 
- Fixed `dbreg_is_balanced_panel()` to correctly detect unbalanced panels with
  missing cells (previously only checked for unequal cell counts).
- The CR1 small-sample correction for clustered SEs now uses
  `(G/(G-1)) * ((N-1)/(N-K))`, matching `fixest`, instead of
  `(G/(G-1)) * (N/(N-K))`. Similarly, `ssc = "nested"` now excludes `G-1` (not
  `G`) parameters for each fixed effect nested within the cluster variable. In
  practice, clustered SEs shift by a relative factor of roughly `1/(2N)` under
  the default `ssc = "full"` (and are essentially unchanged under `"nested"`),
  which is negligible for large datasets. But it does mean that `dbreg` now
  matches `fixest` to numerical precision in small samples too. (#66)
- With `drop_missings = TRUE` (the default), rows with a missing cluster
  variable are now dropped, as in `fixest`. Previously they were kept and
  treated as their own cluster. (#66)
- Clustered SEs for the `"compress"` strategy no longer stall when there are
  many grouping columns. (#66)
- Transformed formula terms, e.g. `log(y) ~ x` or `y ~ log(x) + I(x^2)`, now
  throw an informative error. Previously these were silently reduced to their
  underlying variables (e.g., regressing on `x` instead of `log(x)`). Please
  create such variables in the data (or a database view) first. (#66)
- Requests to omit the intercept (`y ~ 0 + x` or `y ~ x - 1`) are now honoured
  for models without fixed effects, in both `dbreg()` and `dbivreg()`.
  Previously an intercept was always added. As in `lm()`, the first factor
  keeps all of its levels and R2 is uncentred. (#68, #73)
- `predict()` now rebuilds the design matrix from the model formula, fixing
  predictions for models with interactions (e.g., `"mundlak"`), factors, or no
  intercept. (#67, #73)
- Variable names that aren't valid SQL identifiers (e.g., `Sepal.Length`,
  names with spaces, or reserved words) and factor levels with spaces or other
  special characters now work in `dbreg()`, `dbivreg()` and `dbbinsreg()`.
  Results use the original names. (#58)
- Plain factor regressors (without interactions) now work for the
  `"moments"`, `"demean"` and `"mundlak"` strategies, and factor levels are
  read from all rows rather than the first 100,000. (#76)
- `predict()` now errors for `"demean"` models with more than two fixed
  effects. Previously it silently ignored the third and later fixed effects,
  returning incorrect predictions. (#66)
- `sql_model_matrix()` (and thus interaction/factor terms) now generates
  backend-portable schema and factor-level queries, fixing failures on SQL
  Server and other non-DuckDB backends. (#66)
- Sufficient statistics and HC1/clustered scores are now computed in double
  precision on all backends. Previously, they were cast to `FLOAT`, which is
  single precision on DuckDB, MySQL, Spark and Athena, limiting the
  `"demean"` and `"mundlak"` strategies (and robust SEs) to ~1e-8 accuracy. (#78)
- The `"mundlak"` strategy now includes group means of every regressor column,
  incl. factor dummies and interactions. Previously only raw numeric variables
  were averaged, so models with factors or interactions didn't match the
  equivalent fixed effects estimates. Time-invariant regressors now keep their
  coefficient (dropping their collinear group mean) instead of erroring. (#77)
- Clustered SEs for the `"compress"` strategy no longer error for models with
  interaction terms.

Internal

- Update maintainer email.
- Added a `CLAUDE.md` context file for sharper AI-assisted development, as well
  as a `Makefile` for automating dev workflows. (#64)
- Refactored and modularized the codebase to make it more efficient and safer.
  (#65)
- Much faster cluster-robust SEs for the `"compress"` strategy with many
  clusters. (#72)
- Internal arguments now flow through a single `inputs` environment (filled in
  place via `env2env()`), simplifying function signatures. (#74)
- Estimation internals now share a common set of building blocks (estimation
  frames, a cross-moment matrix, estimator fits, and a SQL score pass for
  robust/clustered meat matrices), which `dbreg()` and `dbivreg()` both use.
  This removes a lot of duplicated code and should make it easier to add new
  estimators. (#66)

## dbreg 0.0.3

_**dbreg** v0.0.3 is a big release, containing many new features and some
important breaking changes. We have closed out most of items on our
[TODO](https://github.com/grantmcdermott/dbreg/issues/5) list and this will
probably be our final "beta" release before submitting to CRAN._

Website

- We finally have a package website:
[grantmcdermott.com/dbreg](https://grantmcdermott.com/dbreg/index.html) 🎉

Breaking changes

- The behaviour of the `"mundlak"` strategy has been changed. (#24)
  - The old `"mundlak"` strategy is remapped to the new `"demean"` (alias
    `"within"`) strategy, to better reflect the fact that this strategy invokes
    a (double) demeaning transformation. Users who want the old behaviour should
    thus use `"demean"` instead of `"mundlak"`.
  - Simultaneously, we now provide a revised `"mundlak"` strategy that
    implements a "true" Mundlak/CRE estimator; see New features below.
- For estimations with two fixed effects on _unbalanced_ panels,
  `strategy="auto"` now errors when the compression limits are exceeded. It does
  this to avoid silently selecting a different estimand (i.e., Mundlak/CRE
  instead of TWFE). For these ambiguous cases, users will now be prompted to
  explicitly choose `"compress"` (with higher limits) or `"mundlak"` (different
  model and thus potentially different coefficients). (#24)
- The default `verbose` behaviour is changed to `FALSE`. Users can revert to
  the old behaviour for a single call (i.e., `dbreg(..., verbose = TRUE)`), or
  set it globally (i.e., `options(dbreg.verbose = TRUE)`). (#33)
- Technically not a breaking change, since we currently support backwards
  compatibility, but several minor arguments have been renamed/superseded. (#34)
  - `query_only` -> `sql_only` (in `dbreg`)
  - `fes` -> `fe` (in `print.dbreg`, `coef.dbreg`, `confint.dbreg`, etc.)

Major new features

- We have added new/revised acceleration strategies. (#24)
  - The `"demean"` (alias `"within"`) strategy implements a (double) demeaning
    transformation and is particularly suited to (balanced) panels with one or
    two fixed effects. Note the underlying query is the same as the _old_
    `"mundlak"` strategy, which was somewhat erroneously named. Speaking of
    which...
  - The revised `"mundlak"` strategy now implements a "true" Mundlak/CRE
    (correlated random effects) estimator by regressing Y on X plus group
    means of X. Unlike the `"demean"` strategy (above), this revised `"mundlak"`
    model obtains consistent coefficients regardless of panel structure (incl.
    unbalanced panels) and supports any number of fixed effects. However, users
    should note that Mundlak/CRE is a different model from "vanilla" fixed
    effects---albeit asymptotically equivalent under certain assumptions---and
    may obtain different coefficients as a result.
  - Please consult the expanded `Acceleration Strategies` section in the
    `?dbreg` helpfile for technical details.
- Add support for clustered standard errors. Follows the `fixest` API:

  ```r
  dbreg(..., vcov = ~cluster_var)
  ```

  Please note that these clustered SEs are computed _analytically_ (not
  bootstrapped) and should thus add minimal overhead to your regressions. See
  the updated README for examples. (#29)
- Added support for binscatter regressions on database backends via the
  `dbbinsreg()` function. This function attempts to mimic the main API of the
  **binsreg** package (modulo some simplifications), but should be much faster
  on big datasets. (#32)
- Added support for interaction terms in formulas (e.g., `y ~ x1 * x2 | fe`).
  Interactions are expanded to SQL via the new (user-facing)
  `sql_model_matrix()` function and work with all acceleration strategies.
  Collinear terms (e.g., treatment main effects absorbed by unit FEs in DiD
  models) are automatically detected and dropped. (#48)

  - Note that interactions support also enables basic event-study estimations,
    although the reference period needs to be specified manually, e.g.:

    ```r
    data("base_did", package = "fixest")
    did = transform(base_did, period = factor(period))
    did$period = relevel(did$period, ref = "5")
    dbreg(y ~ x1 + treat:period | id + period, data = did)
    ```

Other new features

- The `"auto"` strategy logic now considers a `compress_nmax` threshold, which
  governs the maximum allowable size of the compressed data object (default
  threshold = 1 million rows). This additional guardrail is intended to avoid
  cases where the `"compress"` strategy satisfies the `compress_ratio`
  threshold, but could still return a prohibitively large dataset. The most
  common example would be querying a massive dataset on a remote database, where
  network latency makes data I/O transfer expensive, even though we've achieved
  good compression relative to the original data size. (#10)
  - Aside: Improved documentation and messaging (when `verbose = TRUE`) should
    also help users understand the `"auto"` strategy decision tree.
- Enabled weights for double demean (within) specification. (#13)
- Estimations now report some goodness-of-fit statistics like R2 and RMSE,
  powered by the (user-facing) `gof()` function. (#21) 
- Added support for various `*.dbreg` methods (#21, #30):
  - From **stats**: `coef()`, `confint()`, `predict()`, and `vcov()`.
  - From **broom**/**generics**: `tidy()` and `glance()`. These also enable
  post-processing operations like exporting results to coefficient tables via
  `modelsummary::msummary()`. Thanks to @HariharanJayashankar for the request in
  #20.
- Better documentation. (#28, #36, #43, and various other PRs)

Bug fixes

- Added QR decomposition fallback for regression calculations, for cases where
  the default Cholesky solver fails. (#7)
- Improved integration for running regressions on AWS Athena datasets via the
  **noctua** package/driver. (#8)
- Automatically drop incomplete cases (i.e., missing values) prior to any
  aggregation steps, avoiding mismatched matrices during estimation. (#19)
- User-specified `compress_ratio` values should now bind in all cases.
  Previously, these could sometimes be silently ignored due to internal
  overrides. Also, clarify in the argument documentation that the default
  (automatic) `compress_ratio` threshold can vary based on heuristics related to
  model structure. (#25)
- Correctly estimate HC1 standard errors for the `"moments"`, `"demean"`, and
  `"mundlak"` strategies. While the old analytic HC1 approach worked (and still
  does) for the `"compress"` case, it led to misleading SEs for these other
  strategies. The fix does impose some additional computational overhead, since
  it requires a second pass over the data to calculate the individual errors
  and "meat" of the sandwich matrix. But testing suggests that this leads to a
  <2x increase in total estimation time, which seems a reasonable tradeoff for
  heteroskedastic-robust SEs. (#27)
- Fixed a bug for in-memory `dbreg(..., data = <data>)` cases, where factor
  variables in X could cause the `"compress"` strategy to fail. (#41) 

Internals

- Added unit testing framework using **tinytest**. (#16)
- Added GitHub Actions CI. (#18)

## dbreg 0.0.2

**IMPORTANT BREAKING CHANGE:**

The package has been renamed to **dbreg** to better reflect the fact that it
supports multiple database backends. (#4)

Other breaking changes

- The default `vcov` is now "iid". (#2 @grantmcdermott)

New features 

- The new `dbreg(..., strategy = <strategy>)` argument allows users to choose
  between different acceleration strategies for efficient computation of the
  regression coefficients and standard errors. This includes `"compress"` (the
  old default), as well as `"mundlak"` or `"moments"`. The latter two strategies 
  are newly introduced in **dbreg** v0.0.2 and offer important advantages in the 
  case of true panel data. If an explicit strategy is not provided by the user, 
  then `dbreg()` will invoke some internal heuristics to determine the optimal
  strategy based on the size and structure of the data. (#2 @jamesbrandecon
  and @grantmcdermott)

Project

- @jamesbrandecon has joined the project as a core contributor.

## duckreg 0.0.1

* Initial GitHub release.
