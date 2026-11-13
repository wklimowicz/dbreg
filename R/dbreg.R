#' Run a regression on a database backend
#'
#' @md
#' @description
#' Leverages the power of databases to run regressions on very large datasets,
#' which may not fit into R's memory. Various acceleration strategies allow for
#' highly efficient computation, while robust standard errors are computed from
#' sufficient statistics.
#'
#' @section Varying slopes:
#' OLS `group[trend]` and `group[trend1, trend2]` absorb a group intercept plus
#' numeric group-specific slopes; `(group1^group2)[trend]` uses native tuples.
#' Repeated group declarations merge their trend columns. Both `"demean"` and
#' `"compress"` are supported. Auto includes trend values in its compression
#' estimate and uses demeaning after a compressed capacity decline. Mundlak
#' and SQL-only within projections are unavailable. IV varying slopes
#' remain unsupported. Trend columns participate in finite/missing filtering
#' and sample audits, independently of their RHS or native-key roles.
#'
#' Within fits prefer bounded signed coordinate moments, with SQL coefficient-CG
#' fallback after a numerical/capacity decline or with `dbreg.ap_solver="iterative"`.
#' Zero/redundant axes preserve the actual projected span. Resource bounds and
#' independent count/finite/absolute/RMS coordinate certificates apply to both
#' engines. Expected moment declines never use additive-only AP for these fits.
#' An explicit `fe_max_iter` (argument or option) bounds CG work, otherwise
#' `dbreg.cg_max_iter` does; either way at most eight refinement steps. Under the
#' declared `K.exact=FALSE` convention, each trend contributes one coefficient
#' per declared group level, even for a locally redundant trend. `K.fixef`
#' changes FE intercept counts; varying coefficients remain in K. Within predictions reconstruct the weighted projection of the supplied cohort,
#' including trends, without recovering training effects.
#'
#' Compressed fits group by the actual design, FE and trend values. Each cell
#' carries one representative native raw tuple and moments in source units.
#' `data_only` returns these cells; `sql_only` returns the runnable cell query.
#' Row, scalar, key byte, allocation and score upload budgets apply. A large
#' score upload can use bounded cluster-by-cell statistics instead. Trend
#' nuisance coefficients use per-group centered coordinates. Compressed
#' predictions rebuild those fitted coordinates and return `NA` for unknown
#' groups or rows outside the identified training span.
#'
#' @param fml A \code{\link[stats]{formula}} representing the relation to be
#' estimated. Fixed effects should be included after a pipe, e.g
#' `fml = y ~ x1 + x2 | fe1 + f2`. Interaction terms are supported using
#' standard R syntax (`:` for interactions, `*` for main effects plus
#' interaction). `fixest`-style `i()` terms are also supported, e.g.
#' `i(period, treat, ref = 5)` for event studies; see
#' [sql_model_matrix()] for exact typed selectors and the `dbreg.i_max_columns`
#' resource limit. Explicit `factor(column)`, `I(column^2)` and `I(column^3)`
#' terms and their interactions are supported; other transformations and
#' literal regressors are rejected. Treatment coding is fixed at estimation and retained for R
#' prediction; changing global contrast options does not recode a fitted model.
#' An intercept-only model (`y ~ 1`) is supported without fixed effects.
#'   Fixed effects accept native tuples such as `school^year`, combined with
#'   additive families (`school^year + teacher`). Components must be source
#'   column names. DuckDB/native Databricks fits freeze tuple ranks in an owned
#'   snapshot and filter missing components before ranking. Tuple identity is
#'   independent of delimiters and numeric casts. Compressed prediction stores
#'   a bounded training-tuple dictionary under the projection row/cell/memory
#'   limits; individual key strings must stay below 65535 bytes. Compressed raw newdata
#'   must provide all components, with integer64 or exact decimal strings for
#'   BIGINT values beyond double precision. Unknown training tuples produce NA.
#'   In compressed output, source columns that conflict with count/moment
#'   fields get separate `_source` names; the `dbreg_source_columns` attribute
#'   records their actual columns. Counts and raw key values remain distinct.
#' @param conn Database connection, e.g. created with
#' \code{\link[DBI]{dbConnect}}. Can be either persistent (disk-backed) or
#' ephemeral (in-memory). If no connection is provided, then an ephemeral
#' \code{\link[duckdb]{duckdb}} connection will be created automatically and
#' closed before the function exits. Note that a persistent (disk-backed)
#' database connection is required for larger-than-RAM datasets in order to take
#' advantage of out-of-core functionality like streaming (where supported).
#' Default connections for integer64 data use exact BIGINT reads. Native BIGINT
#' registration is verified before fitting; failed setup releases its resources.
#' Caller connection settings are preserved. On DuckDB, large numeric BIGINT
#' category/cell keys are re-read exactly using bit64; numeric regression
#' arithmetic still uses double precision.
#' @param table,data,path Arguments for specifying the data
#' table (object) to be queried. In order of precedence:
#' - `table`: Character string giving the name of the data table in an
#' existing (open) database connection.
#' - `data`: R dataframe registered with DuckDB, or a lazy table with a
#' connection that can be inferred. Ignored if `table` is provided.
#' - `path`: Character string giving a path to the data file(s) on disk, which
#' will be read into `conn`. Internally, this string is passed to the `FROM`
#' query statement, so could (should) include file globbing for
#' Hive-partitioned datasets, e.g. `"mydata/**/.*parquet"`. For more precision,
#' however, it is recommended to pass the desired database reader function as
#' part of this string, e.g. `"read_parquet('mydata/**/*.parquet')"` for DuckDB;
#' note the use of single quotes.
#' Ignored if either `table` or `data` is provided.
#'
#' Lazy tables supplied as `data` are treated as `table`, including connection inference.
#' Unused lower-priority sources are ignored before validation and heuristics.
#' @param weights Column name or one-sided formula such as `~w` naming weights,
#' or NULL (default) for unweighted regression. Weights must be non-negative;
#' rows with zero weight are dropped. Weighted regressions support `"iid"`,
#' `"hc1"`, and clustered SEs, and are compatible with all strategies.
#' @param factor_refs Optional named reference levels for ordinary categorical
#'   regressors. Naming a numeric-coded column treats it as categorical. Each
#'   entry must identify one observed level; i() terms use their own selectors.
#'   Exact fitted coding is retained for prediction across strategies.
#' @param vcov Character string or formula denoting the desired type of variance-
#' covariance correction / standard errors. Options are `"iid"` (default),
#' `"hc1"` (heteroskedasticity-consistent), or a one-sided formula like
#' `~cluster_var` for cluster-robust standard errors. Note that `"hc1"` and
#' clustered SEs require a second pass over the data unless
#' `strategy = "compress"` to construct the residuals.
#' @param strategy Character string indicating the preferred acceleration
#'   strategy. The default `"auto"` will pick an optimal strategy based on
#'   internal heuristics. Users can also override with one of the following
#'   explicit strategies: `"compress"`, `"demean"` (alias: `"within"`),
#'   `"mundlak"`, or `"moments"`. See the Acceleration Strategies section below
#'   for details.
#' @param compress_ratio,compress_nmax Numeric(s). Parameters that help to
#'   determine the acceleration `strategy` under the default `"auto"` option.
#'
#'   - `compress_ratio` defines the compression ratio threshold, i.e. numeric
#'     in the range `[0,1]` defining the minimum acceptable compressed versus
#'     the original data size. Default value of `NULL` means that the threshold
#'     will be automatically determined based on some internal heuristic
#'     (e.g., 0.01 for models without fixed effects).
#'   - `compress_nmax` defines the maximum allowable size (in rows) of the
#'     compressed dataset that can be serialized into R. Pays heed to the idea
#'     that big data serialization can be costly (esp. for remote databases),
#'     even if we have achieved good compression on top of the original dataset.
#'     Default value is 1e6 (i.e., a million rows).
#'
#'   See the Acceleration Strategies section below for further details.
#' @param cluster Optional. Provides an alternative way to specify
#' cluster-robust standard errors (i.e., instead of `vcov = ~cluster_var`).
#' Either a one-sided formula (e.g., `~firm`) or character string giving the
#' variable name. Only single-variable clustering is currently supported.
#' @param ssc Small-sample policy: `"full"` (default), `"nested"`, `"none"`,
#'   or a named list with `K.adj`, `K.fixef`, `K.exact`, `G.adj`, `G.df`, and
#'   `t.df`, as in `fixest::ssc()`. `"full"` counts all fixed-effect parameters;
#'   `"nested"` excludes effects nested in the cluster; `"none"` excludes
#'   fixed-effect parameters. Mundlak group means remain ordinary regressors.
#'   List defaults are `K.adj=TRUE`, `K.fixef="nonnested"`, `K.exact=FALSE`,
#'   `G.adj=TRUE`, `G.df="min"`, `t.df="min"`. `K.exact=TRUE` is unsupported.
#'   A finite numeric `t.df` of at least one overrides the test degrees of
#'   freedom; otherwise clustered tests use G-1 (`"min"`) or N-K
#'   (`"conventional"`). Coefficient tables and intervals share this test df.
#'   Residual df always counts the full declared model. Inference requires
#'   positive residual df and, for clustering, at least two clusters. Under
#'   `K.exact=FALSE`, clustering directly on an FE uses the declared matching-FE
#'   nesting convention; external cluster columns use functional nesting checks.
#' @param sql_only Logical indicating whether to return only the underlying SQL
#'   query (default `FALSE`). For iterative AP, returns a `dbreg_query_plan`
#'   of ordinary-projection SQL templates and certificates without fitting.
#'   See [dbreg-options].
#' @param data_only Logical indicating whether only the compressed dataset
#'   should be returned (i.e., no regression is run). Default is `FALSE`.
#'   Compressed cells use source coordinates. Moment rows use centered, weight-
#'   normalized coordinates and include `dbreg_moment_center_*` and
#'   `dbreg_moment_weight_scale` metadata. Raw compressed sums returned here
#'   remain subject to the numeric range of R doubles.
#' @param drop_missings Logical indicating whether incomplete cases (i.e., rows
#'   where any of the dependent, independent or FE variables are
#'   missing) should be dropped. The default is `TRUE`, according with standard
#'   regression software. It is *strongly* recommended not to change this value
#'   unless you are absolutely sure that your data have no missings and you wish
#'   to skip some internal checks. (Even then, it probably isn't worth it.)
#' @param row_id Optional character vector of unique, nonmissing source key
#'   columns. Keys are validated before subset/missing/weight filtering, using
#'   native SQL tuples; no key values are collected for validation.
#' @param sample_table Optional new temporary table receiving the fitted source
#'   keys. Requires `row_id` and a caller-owned connection (including inference
#'   from a lazy table). DuckDB/native Databricks exports use an owned snapshot
#'   for consistent fitting and keys. Existing relations are never replaced.
#'   The caller owns the export; private intermediates are cleaned on errors
#'   and success. Unacknowledged cleanup remains queued for synchronous retry.
#'   Cannot be combined with `sql_only` or `data_only`.
#' @param diagnostics Logical. Run an additional source aggregate for input,
#'   subset-selected, missing-by-role, nonpositive-weight and invalid-offset
#'   counts? Default is `FALSE`. Counts may overlap. Fitted objects return
#'   `sample_diagnostics` and `nobs_input`; inspection modes skip the audit.
#' @param offset Optional numeric column name or one-sided formula such as `~o`,
#'   added to the predictor with coefficient one. Estimation uses outcome minus
#'   offset. Rows with nonfinite outcome, offset or adjusted difference are
#'   excluded. Supports DuckDB/native Databricks. Prediction requires newdata
#'   containing the offset; intervals condition on its known value.
#' @param subset Optional logical column name or restricted one-sided predicate
#'   formula, evaluated before missing-value/weight filtering. Supports column
#'   comparisons, arithmetic, boolean operators, `is.na()` and literal `%in%`.
#'   Supports DuckDB and native Databricks; positional vectors are rejected.
#'   Weights are validated on the source before subsetting.
#' @param verbose Logical. Print auto strategy and progress messages to the
#'   console? Defaults to `FALSE`. This can be overridden for a single call
#'   by supplying `verbose = TRUE`, or set globally via
#'   `options(dbreg.verbose = TRUE)`.
#' @param factor_refs Optional named reference levels for ordinary categorical
#'   regressors. Naming a numeric-coded column treats it as categorical. Each
#'   entry must identify one observed level; `i()` terms use their own selectors.
#'   Exact fitted coding is retained for prediction across strategies.
#' @param fixef_rm Character string controlling singleton fixed-effect removal:
#'   `"none"` (default) keeps all observations; `"singletons"` iteratively
#'   prunes observations belonging to singleton fixed-effect groups.
#' @param fe_tol Positive per-fit tolerance for fixed-effect projection. Defaults
#'   to `getOption("dbreg.ap_tol", 1e-10)` without changing that option.
#' @param fe_max_iter Positive maximum AP cycles or bounded solve/refinement
#'   steps for this fit. Defaults to `getOption("dbreg.ap_max_iter", 100L)`
#'   without changing it. Bounded refinement uses at most eight steps; inner
#'   PCG has a separate cap. See [dbreg-options].
#' @param ... Additional arguments. Currently ignored, except to handle
#'   superseded arguments for backwards compatibility.
#'
#' @section Bounded fixed-effect projection:
#' DuckDB/native Databricks can use bounded FE moment projection before AP;
#' see [dbreg-options] for capacity limits, explicit solver selection and
#' projection diagnostics. Source keys stay in SQL and private tables are
#' released on success, decline and errors.
#'
#' @section Sample audits:
#' With `diagnostics = TRUE`, `sample_diagnostics` reports original source rows
#' (`n_input`), subset-selected rows (`n_selected`), final original observations
#' (`n_fitted`) and their difference (`n_dropped`). `detail` gives missing counts
#' by public source role before filtering; DuckDB/Spark numeric NaNs are missing.
#' Nonpositive-weight and invalid-offset counts also refer to the original
#' source. An invalid offset has a nonfinite operand or adjusted difference.
#' These counts overlap and cannot be added to reconstruct the fitted count.
#' `query_string` reproduces the source aggregate, with columns mapped by
#' `count_aliases`. Only metadata and counts are collected. Without an audit,
#' `nobs_input` is `NA`; `sql_only` and `data_only` skip auditing.
#'
#' @return A list of class "dbreg" containing various slots, including a table
#' of coefficients (which the associated `print` method will display).
#'
#' @section Acceleration Strategies:
#'
#' `dbreg` offers four primary acceleration strategies for estimating regression
#' results from simplified data representations. Below we use the shorthand
#' Y (outcome), X (explanatory variables), and FE (fixed effects) for exposition
#' purposes:
#'
#' 1. `"compress"`: compresses the data via a `GROUP BY` operation (using X and
#'    the FE as groups), before running weighted least squares on this much
#'    smaller dataset:
#'    \deqn{\hat{\beta} = (X_c' W X_c)^{-1} X_c' W Y_c}
#'    where \eqn{W = \text{diag}(n_g)} are the group frequencies. This procedure
#'    follows Wong et al. (2021).
#' 2. `"moments"`: computes sufficient statistics (\eqn{X'X, X'y}) directly via
#'    SQL aggregation, returning a single-row result. This solves the standard
#'    OLS normal equations \eqn{\hat{\beta} = (X'X)^{-1}X'y}. Limited to cases
#'    without FE.
#' 3. `"demean"` (alias `"within"`): subtracts group-level means from both Y and
#'    X before computing sufficient statistics (per the `"moments"` strategy).
#'    For example, given unit \eqn{i} and time \eqn{t} FE, we apply double
#'    demeaning:
#'    \deqn{\ddot{Y}_{it} = \beta \ddot{X}_{it} + \varepsilon_{it}}
#'    where \eqn{\ddot{X} = X - \bar{X}_i - \bar{X}_t + \bar{X}}. This
#'    (single-pass) within transformation is algebraically equivalent to the
#'    fixed effects projection---i.e., Frisch-Waugh-Lovell partialling out---in
#'    the presence of a single FE. It is also identical for the two-way FE
#'    (TWFE) case if your panel is balanced. For unbalanced two-way panels,
#'    however, the double demeaning strategy is not algebraically equivalent to
#'    the fixed effects projection and therefore does not recover the exact TWFE
#'    coefficients. In such cases, and also for weighted two-way FE, `dbreg`
#'    uses alternating projections (AP) to recover the exact FE coefficients,
#'    at the cost of extra passes over the data. AP also generalizes the
#'    `"demean"` strategy to three or more FE. See [dbreg-options] for AP
#'    tuning options.
#' 4. `"mundlak"`: a generalized Mundlak (1978), or correlated random effects
#'    (CRE) estimator that regresses Y on X plus group means of X:
#'    \deqn{Y_{it} = \alpha + \beta X_{it} + \gamma \bar{X}_i + \varepsilon_{it} \quad \text{(one-way)}}
#'    \deqn{Y_{it} = \alpha + \beta X_{it} + \gamma \bar{X}_{i} + \delta \bar{X}_{t} + \varepsilon_{it} \quad \text{(two-way, etc.)}}
#'    Unlike `"demean"`, Y is not transformed, so predictions are on the
#'    original scale. Supports any number of FE and works correctly for any
#'    panel structure (balanced or unbalanced). However, note that CRE is a
#'    *different model* from FE: while coefficients are asymptotically
#'    equivalent under certain assumptions, they will generally differ in
#'    finite samples.
#'
#' The relative efficiency of each of these strategies depends on the size and
#' structure of the data, as well as the number of unique regressors and FE. For
#' (quote unquote) "standard" cases, the `"compress"` strategy can yield
#' remarkable performance gains and should justifiably be viewed as a good
#' default. However, the compression approach tends to be less efficient for
#' true panels (repeated cross-sections over time), where N >> T. In such
#' cases, it can be more efficient to use a demeaning strategy that first
#' controls for (e.g. subtracts) group means, before computing sufficient
#' statistics on the aggregated data. The reason for this is that time and unit
#' FE are typically high dimensional, but covariate averages are not; see
#' Arkhangelsky & Imbens (2024).
#' 
#' However, the demeaning approaches invite tradeoffs of their own. For example,
#' the single-pass double demeaning transformation only obtains exact TWFE
#' results for balanced panels with two FE. For unbalanced panels, weighted
#' regressions, or three or more FE, `dbreg` uses alternating projections
#' (iterative demeaning) which is exact but requires multiple passes and may be
#' slower to converge on very large datasets. In such cases, `"mundlak"` (CRE)
#' may be preferable as it is a single-pass estimator that obtains consistent
#' coefficients regardless of panel structure and FE count, but at the "cost" of
#' recovering a different estimand. (It is a different model to TWFE, after
#' all.) See Wooldridge (2025) for an extended discussion of these issues.
#' 
#' Users should weigh these tradeoffs when choosing their acceleration strategy.
#' Summarising, we can provide a few guiding principles. `"compress"` is a good
#' default that guarantees the "exact" FE estimates and is usually very 
#' efficient (barring data I/O costs and high FE dimensionality). `"mundlak"` is
#' another efficient alternative provided that the CRE estimand is acceptable
#' (don't be alarmed if your coefficients are not identical). Finally, the
#' `"demean"` and `"moments"` strategies are great for particular use cases
#' (i.e., balanced panels and cases without FE, respectively).
#' 
#' If this all sounds like too much to think about, don't fret. The good news
#' is that `dbreg` can do a lot (all?) of the deciding for you. Specifically, it
#' will invoke an `"auto"` heuristic behind the scenes if a user does not
#' provide an explicit acceleration strategy. Working through the heuristic
#' logic does impose some additional overhead, but this should be negligible in
#' most cases (certainly compared to the overall time savings). The `"auto"`
#' heuristic is as follows:
#'
#' - IF no FE AND (any continuous regressor OR poor compression ratio OR too big
#'   compressed data) THEN `"moments"`.
#' - ELSE IF 1 FE AND (poor compression ratio OR too big compressed data) THEN
#'   `"demean"`.
#' - ELSE IF 2 FE AND (poor compression ratio OR too big compressed data):
#'   - IF balanced panel THEN `"demean"`.
#'   - ELSE `"demean"` via alternating projections.
#' - ELSE IF 3+ FE AND (poor compression ratio OR too big compressed data)
#'   THEN `"demean"` via alternating projections.
#' - ELSE THEN `"compress"`.
#' 
#' _Tip: set `dbreg(..., verbose = TRUE)` to print information about the auto
#' strategy decision criteria._
#'
#' @references
#' Arkhangelsky, D. & Imbens, G. (2024)
#' \cite{Fixed Effects and the Generalized Mundlak Estimator}.
#' The Review of Economic Studies, 91(5), pp. 2545–2571.
#' Available: https://doi.org/10.1093/restud/rdad089
#'
#' Mundlak, Y. (1978)
#' \cite{On the Pooling of Time Series and Cross Section Data}.
#' Econometrica, 46(1), pp. 69–85.
#' Available: https://doi.org/10.2307/1913646
#'
#' Wong, J., Forsell, E., Lewis, R., Mao, T., & Wardrop, M. (2021).
#' \cite{You Only Compress Once: Optimal Data Compression for Estimating Linear Models.}
#' arXiv preprint arXiv:2102.11297.
#' Available: https://doi.org/10.48550/arXiv.2102.11297
#'
#' Wooldridge, J.M. (2025)
#' \cite{Two-way fixed effects, the two-way mundlak regression, and difference-in-differences estimators}.
#' Empirical Economics, 69, pp. 2545–2587.
#' Available: https://doi.org/10.1007/s00181-025-02807-z
#'
#' @seealso \code{\link[DBI]{dbConnect}} for creating database connections,
#' \code{\link[duckdb]{duckdb}} for DuckDB-specific connections
#'
#' @importFrom DBI dbConnect dbDisconnect dbExecute dbGetInfo dbGetQuery dbIsValid dbRemoveTable dbWriteTable
#' @importFrom duckdb duckdb duckdb_register duckdb_unregister
#' @importFrom Formula Formula
#' @importFrom Matrix chol2inv crossprod Diagonal sparse.model.matrix
#' @importFrom stats aggregate as.formula formula pt reformulate setNames terms
#' @importFrom glue glue glue_sql
#'
#' @examples
#' ## In-memory data ----
#' 
#' # We can pass in-memory R data frames to an ephemeral DuckDB connection via
#' # the `data` argument. This is convenient for small(er) datasets and demos.
#'
#' # Default "compress" strategy reduces the data to 4 rows before running OLS
#' dbreg(weight ~ Diet, data = ChickWeight)
#' 
#' # Compare with lm
#' summary(lm(weight ~ Diet, data = ChickWeight))$coefficients
#' 
#' # Add "fixed effects" after a `|` 
#' dbreg(weight ~ Time | Diet, data = ChickWeight)
#' 
#' # "robust" SEs can also be computed using a sufficient statistics approach
#' dbreg(weight ~ Time | Diet, data = ChickWeight, vcov = "hc1")
#' dbreg(weight ~ Time | Diet, data = ChickWeight, vcov = ~Chick)
#' 
#' # Different acceleration strategies + specifications
#' dbreg(weight ~ Time | Diet, data = ChickWeight, strategy = "demean")
#' dbreg(weight ~ Time | Diet, data = ChickWeight, strategy = "mundlak")
#' dbreg(weight ~ Time | Diet + Chick, data = ChickWeight, strategy = "mundlak") # two-way Mundlak
#' dbreg(weight ~ Time, data = ChickWeight, strategy = "moments") # no FEs
#' # etc.
#' 
#' # Interactions: does the effect of Time vary by Diet?
#' # (Diet main effects are collinear with Chick FE, so these drop out)
#' dbreg(weight ~ Time * Diet | Chick, data = ChickWeight)
#' 
#' # Event study with fixest-style i() terms (reference period 5)
#' if (requireNamespace("fixest", quietly = TRUE)) {
#'   data("base_did", package = "fixest")
#'   dbreg(y ~ x1 + i(period, treat, ref = 5) | id + period, data = base_did)
#' }
#' 
#' #
#' ## DBI connection ----
#' 
#' # For persistent databases or more control, use the `conn` + `table` args.
#' # Again, we use DuckDB below but any other DBI-supported backend should work
#' # too (e.g., odbc, bigrquery, noctua (AWS Athena),  etc.) See:
#' # https://r-dbi.org/backends/
#' 
#' library(DBI)
#' con = dbConnect(duckdb::duckdb())
#' dbWriteTable(con, "cw", as.data.frame(ChickWeight))
#' 
#' dbreg(weight ~ Time | Diet, conn = con, table = "cw")
#' 
#' # Tip: Rather than creating or writing (temp) tables, use CREATE VIEW to
#' # define subsets or computed columns without materializing data. This is more
#' # efficient and especially useful for filtering or adding variables.
#' dbExecute(
#'   con,
#'   "
#'   CREATE VIEW cw1 AS
#'   SELECT *
#'   FROM cw
#'   WHERE Diet = 1
#'   "
#' )
#' dbreg(weight ~ Time | Chick, conn = con, table = "cw1")
#' 
#' #
#' ## Path to file ----
#' #
#' # For file-based data (e.g., parquet), use the path argument.
#' 
#' tmp = tempfile(fileext = ".parquet")
#' dbExecute(con, sprintf("COPY cw TO '%s' (FORMAT PARQUET)", tmp))
#' 
#' dbreg(weight ~ Time | Diet, path = tmp)
#' 
#' # Cleanup
#' dbDisconnect(con)
#' unlink(tmp)
#' 
#' #
#' ## Big dataset ----
#' 
#' # For a more compelling and appropriate dbreg use-case, i.e. regression on a
#' # big (~180 million row) dataset of Hive-partioned parquet files, see the
#' # package website:
#' # https://grantmcdermott.com/dbreg/
#' @export
dbreg = function(
  fml,
  conn = NULL,
  table = NULL,
  data = NULL,
  path = NULL,
  weights = NULL,
  vcov = c("iid", "hc1"),
  strategy = c("auto", "compress", "moments", "demean", "within", "mundlak"),
  compress_ratio = NULL,
  compress_nmax = 1e6,
  cluster = NULL,
  ssc = "full",
  sql_only = FALSE,
  data_only = FALSE,
  drop_missings = TRUE,
  verbose = getOption("dbreg.verbose", FALSE),
  factor_refs = NULL,
  subset = NULL,
  offset = NULL,
  diagnostics = FALSE,
  fixef_rm = c("none", "singletons"),
  row_id = NULL,
  sample_table = NULL,
  fe_tol = getOption("dbreg.ap_tol", 1e-10),
  fe_max_iter = getOption("dbreg.ap_max_iter", 100L),
  ...
) {

  fixef_rm = match.arg(fixef_rm)
  # superseded args handled through ...
  dots = list(...)
  if (length(dots)) {
    if (!is.null(dots[["query_only"]]) && !identical(sql_only, dots[["query_only"]])) {
      sql_only = dots[["query_only"]]
      warning("The `query_only` argument has been superseded by `sql_only` and will be deprecated in a future `dbreg` release.\n")
    }
  }
  if (fixef_rm == "singletons" && (isTRUE(sql_only) || !isTRUE(drop_missings))) {
    stop("Singleton removal requires execution and missing-row filtering.", call. = FALSE)
  }
  source = normalize_data_source(table, data, path)
  table = source[["table"]]
  data = source[["data"]]
  path = source[["path"]]
  weights = formula_column(weights, "weights")
  validate_diagnostics(diagnostics)
  offset = linear_offset_column(offset)
  subset_variables(subset)
  verbose = isTRUE(verbose)
  ssc = normalize_ssc(ssc)
  vcov_parsed = parse_vcov_args(vcov, cluster)
  vcov_type_req = vcov_parsed[["vcov_type"]]
  cluster_var = vcov_parsed[["cluster_var"]]
  strategy = match.arg(strategy)
  if (strategy == "within") strategy = "demean"  # alias


  # Collect the (normalised) arguments into the `inputs` environment, which
  # the rest of the pipeline reads from and writes to (reference semantics)
  inputs = new.env(parent = emptyenv())
  env2env(
    environment(),
    inputs,
    c(
      "fml", "conn", "table", "data", "path", "weights", "vcov_type_req",
      "cluster_var", "ssc", "strategy", "sql_only", "data_only", "compress_ratio",
      "compress_nmax", "drop_missings", "verbose"
    )
  )
  inputs[["factor_refs"]] = factor_refs
  validate_sample_args(row_id, sample_table, conn, table, sql_only, data_only)
  inputs[["row_id"]] = row_id
  inputs[["sample_table"]] = sample_table
  inputs[["resources"]] = fit_resources()
  inputs[["subset"]] = subset
  inputs[["offset"]] = offset
  inputs[["diagnostics"]] = diagnostics
  inputs[["fixef_rm"]] = fixef_rm
  inputs[["fe_tol"]] = fe_tol
  inputs[["fe_max_iter"]] = fe_max_iter
  inputs[["fe_max_iter_explicit"]] = !missing(fe_max_iter) || !is.null(getOption("dbreg.ap_max_iter"))
  validate_fe_controls(inputs)
  # Release any connection/registration/temp tables we create, including on errors during
  # input processing (process_dbreg_inputs() records them in `inputs` first)
  on.exit({
    if (length(inputs[["temp_tables"]]) && !is.null(inputs[["conn"]])) {
      for (tbl in rev(inputs[["temp_tables"]])) {
        try(DBI::dbExecute(inputs[["conn"]], paste("DROP TABLE IF EXISTS", tbl)), silent = TRUE)
      }
    }
    release_fit_resources(inputs[["resources"]])
  }, add = TRUE)

  # Process and validate inputs (writes parsed/derived fields into `inputs`)
  process_dbreg_inputs(inputs)

  # Choose strategy (mutates inputs[["is_balanced"]], inputs[["compression_ratio_est"]])
  chosen_strategy = choose_strategy(inputs)

  # Auto may select small varying cells; a documented capacity decline can
  # still use the certified within engine, after releasing the failed snapshot.
  if (chosen_strategy == "compress" && strategy == "auto" && !sql_only && length(unlist(inputs[["fe_slopes"]]))) {
    attempt = projection_attempt("compressed_coordinates", execute_compress_strategy(inputs))
    result = record_projection_attempt(inputs, attempt)
    if (is.null(result)) {
      chosen_strategy = "demean"
      result = execute_demean_strategy(inputs)
    }
  } else result = switch(
    chosen_strategy,
    # sufficient statistics with no fixed effects
    "moments" = execute_moments_strategy(inputs),
    # fixed effects via demeaning (1 FE: analytic; 2+ FE: alternating projections)
    "demean" = execute_demean_strategy(inputs),
    # true Mundlak/CRE: Y ~ X + group means of X
    "mundlak" = execute_mundlak_strategy(inputs),
    # group by regressors (+ fixed effects) -> frequency-weighted rows -> WLS
    "compress" = execute_compress_strategy(inputs),
    stop("Unknown strategy: ", chosen_strategy)
  )
  # Finalize result
  finalize_dbreg_result(result, inputs, chosen_strategy)
}

#' Process and validate dbreg inputs
#'
#' Reads the user arguments from the `inputs` environment (see [dbreg()]) and
#' writes the parsed and derived fields back into it: the database connection,
#' the formula components, the filtered `from_statement`, etc.
#'
#' @param inputs Environment created by [dbreg()].
#' @keywords internal
process_dbreg_inputs = function(inputs) {
  env2env(
    inputs,
    environment(),
    c(
      "fml", "conn", "table", "data", "path", "weights", "cluster_var",
      "compress_ratio", "drop_missings"
    )
  )

  db_setup = setup_db_connection(conn, table, data, path, caller = "dbreg", resources = inputs[["resources"]])
  conn = db_setup[["conn"]]
  own_conn = db_setup[["own_conn"]]
  registered_table = db_setup[["registered_table"]]
  # Record the connection straight away so the caller can release it, even if
  # validation below fails
  env2env(environment(), inputs, c("conn", "own_conn", "registered_table"))
  from_statement = db_setup[["from_statement"]]

  # Parse (and validate) the formula with the user's names
  fml_parsed = parse_regression_formula(fml, allow_empty = TRUE)
  any_continuous = any_continuous_vars(data, fml_parsed[["xvars"]])
  weights_display = weights
  factor_refs = validate_factor_refs(inputs[["factor_refs"]], fml_parsed[["xvars"]])

  # Give model variables SQL-safe names (see safe_source()) and re-parse;
  # results are mapped back to display names at the end
  from_statement = prepare_linear_sample(inputs, from_statement, c(model_vars(fml_parsed[["fml"]]), weights, cluster_var, inputs[["offset"]], subset_variables(inputs[["subset"]])))
  src = safe_source(conn, from_statement,
    c(model_vars(fml_parsed[["fml"]]), weights, cluster_var, inputs[["offset"]],
      subset_variables(inputs[["subset"]]), inputs[["row_id"]]),
    reserved = compressed_moment_fields())
  from_statement = src[["from_statement"]]
  name_map = src[["map"]]
  inputs[["row_id_safe"]] = src[["to_safe"]](inputs[["row_id"]])
  source_map = name_map
  names(factor_refs) = src[["to_safe"]](names(factor_refs))
  inputs[["factor_refs"]] = factor_refs
  fml_parsed = parse_regression_formula(src[["to_safe"]](fml), allow_empty = TRUE)
  weights = src[["to_safe"]](weights)
  cluster_var = src[["to_safe"]](cluster_var)
  fml = fml_parsed[["fml"]]
  if (length(unlist(fml_parsed[["fe_slopes"]]))) inputs[["varying_fml"]] = fml
  yvar = fml_parsed[["yvar"]]
  xvars = fml_parsed[["xvars"]]
  term_labels = fml_parsed[["term_labels"]]
  has_interactions = fml_parsed[["has_interactions"]] || length(factor_refs) > 0L
  fe = fml_parsed[["fe"]]
  inputs[["fe_slopes"]] = fml_parsed[["fe_slopes"]]
  trend_vars = unique(unlist(inputs[["fe_slopes"]], use.names = FALSE))
  if (length(trend_vars)) {
    schema = DBI::dbGetQuery(conn, sql_limit(paste("SELECT", paste(as.character(DBI::dbQuoteIdentifier(conn, trend_vars)), collapse = ","), from_statement), 0, detect_backend(conn)[["name"]]))
    if (any(!vapply(schema, is.numeric, logical(1)))) stop("Varying-slope columns in fixed effects must be numeric.", call. = FALSE)
  }
  # With FE the intercept is absorbed, so `- 1` / `0 +` only matters without FE
  intercept = fml_parsed[["has_intercept"]] || length(fe) > 0

  # In-memory data are checked by (display) column name; otherwise the database
  # is probed using the (safe) source names
  validate_weights(if (is.null(data)) weights else weights_display, data, conn, from_statement)
  audit_from_statement = from_statement
  from_statement = apply_source_subset(conn, from_statement, inputs[["subset"]], src[["to_safe"]])
  adjusted = prepare_linear_offset(inputs, from_statement, yvar, src)
  from_statement = adjusted[["from_statement"]]
  yvar = adjusted[["yvar"]]
  inputs[["sample_diagnostics"]] = sample_diagnostics(
    inputs, audit_from_statement, list(outcome = fml_parsed[["yvar"]], regressors = xvars,
      fe = fe, trends = trend_vars, weights = weights, cluster = cluster_var, offset = src[["to_safe"]](inputs[["offset"]])), src
  )

  compress_ratio = resolve_compress_ratio(compress_ratio, fe)

  from_statement = filter_estimation_sample(
    from_statement,
    vars = c(yvar, xvars, fe, trend_vars, cluster_var),
    weights = weights,
    drop_missings = drop_missings,
    conn = conn
  )
  inputs[["fe_components"]] = fml_parsed[["fe_components"]]
  composite = prepare_composite_fe(inputs, from_statement, fml_parsed[["fe_components"]])
  from_statement = composite[["from"]]
  fe = composite[["fe"]]
  inputs[["composite_fes"]] = composite[["specs"]]
  inputs[["composite_types"]] = composite[["types"]]
  if (length(composite[["specs"]]) || length(trend_vars)) {
    fml = replace_fe_formula(fml, fe)
    name_map = c(name_map, setNames(fe_component_labels(lapply(composite[["specs"]], restore_names, map = source_map)), names(composite[["specs"]])))
  }
  if (identical(inputs[["fixef_rm"]], "singletons") && length(fe)) {
    pruned = linear_prune_singletons(inputs, from_statement, fe, yvar)
    from_statement = paste("FROM", pruned[["table"]])
    inputs[["pruning"]] = pruned
  }
  is_balanced = NULL

  env2env(
    environment(),
    inputs,
    c(
      "fml", "yvar", "xvars", "term_labels", "has_interactions", "intercept", "fe",
      "from_statement", "compress_ratio", "any_continuous", "is_balanced",
      "weights", "cluster_var", "name_map", "source_map"
    )
  )
  invisible(inputs)
}

#' Heuristic: does any regressor look continuous?
#'
#' Only checkable for in-memory `data`; returns `FALSE` otherwise. Used by
#' the `"auto"` strategy to steer away from `"compress"`.
#' @keywords internal
any_continuous_vars = function(data, vars) {
  if (is.null(data)) {
    return(FALSE)
  }
  is_continuous = function(v) {
    xv = data[[v]]
    if (is.factor(xv) || is.character(xv) || is.logical(xv) || is.integer(xv)) {
      return(FALSE)
    }
    if (is.numeric(xv)) {
      return(length(unique(xv)) > min(50, 0.2 * length(xv)))
    }
    FALSE
  }
  any(vapply(vars, is_continuous, logical(1)))
}

#' Validate or default the `compress_ratio` threshold
#' @keywords internal
resolve_compress_ratio = function(compress_ratio, fe) {
  if (is.null(compress_ratio)) {
    # Stricter compress_ratio logic for 1 and 2 FE cases
    return(if (length(fe) %in% 1:2) 0.6 else 0.01)
  }
  if (!(is.numeric(compress_ratio) && compress_ratio >= 0 && compress_ratio <= 1)) {
    stop("Argument `compress_ratio` ratio must be a numeric in the range [0, 1]\n.")
  }
  compress_ratio
}

#' Validate a weights argument
#'
#' Checks that `weights` names a single column and, where cheap, that its
#' values are finite (in-memory data, DuckDB and Spark) and non-negative.
#' @keywords internal
validate_weights = function(weights, data, conn, from_statement) {
  if (is.null(weights)) {
    return(invisible(NULL))
  }
  if (!is.character(weights) || length(weights) != 1) {
    stop("`weights` must be a single character string (column name) or NULL.")
  }
  if (!is.null(data)) {
    if (!weights %in% names(data)) {
      stop("Weight column '", weights, "' not found in data.")
    }
    wv = data[[weights]]
    if (any(!is.finite(wv) & !is.na(wv))) {
      stop("Weights must be finite.")
    }
    if (any(wv < 0, na.rm = TRUE)) {
      stop("Weights must be non-negative.")
    }
  } else {
    backend = detect_backend(conn)[["name"]]
    invalid = glue("{weights} < 0")
    if (backend %in% c("duckdb", "spark")) {
      invalid = glue("({invalid}) OR ABS({weights}) = CAST('Infinity' AS DOUBLE)")
    }
    sql = sql_limit(
      glue("SELECT CASE WHEN {weights} < 0 THEN 1 ELSE 2 END AS invalid FROM (SELECT * {from_statement}) t WHERE {invalid}"),
      1,
      backend
    )
    bad = tryCatch(dbGetQuery(conn, sql), error = function(e) {
      stop("Could not validate weight column '", weights, "': ", conditionMessage(e), call. = FALSE)
    })
    if (nrow(bad)) {
      if (bad[["invalid"]][1] == 1) stop("Weights must be non-negative.")
      stop("Weights must be finite.")
    }
  }
  invisible(NULL)
}

#' Restrict the estimation sample
#'
#' Appends a WHERE clause dropping rows with missing `vars` (if
#' `drop_missings`) and non-positive weights. Native DuckDB/Spark numeric NaNs
#' are treated as missing, independently of DBI upload conversion.
#' @keywords internal
filter_estimation_sample = function(from_statement, vars, weights = NULL, drop_missings = TRUE, conn = NULL) {
  vars = unique(vars[!is.na(vars) & nzchar(vars)])
  where_clauses = c(
    if (isTRUE(drop_missings)) paste0(vars, " IS NOT NULL"),
    if (!is.null(weights)) paste0(weights, " > 0")
  )
  backend = if (!is.null(conn)) detect_backend(conn)[["name"]] else NULL
  missing = unique(c(if (isTRUE(drop_missings)) vars, weights))
  if (identical(backend, "duckdb") && length(missing)) {
    # Constant typeof() branches preserve text categories such as "NaN".
    where_clauses = c(where_clauses, paste0(
      "(CASE WHEN typeof(", missing,
      ") IN ('DOUBLE','FLOAT') THEN NOT isnan(TRY_CAST(", missing, " AS DOUBLE)) ELSE TRUE END)"
    ))
  } else if (identical(backend, "spark") && length(missing)) {
    schema = DBI::dbGetQuery(conn, sql_limit(paste("SELECT *", from_statement), 0, backend))
    numeric = missing[vapply(missing, function(name) is.double(schema[[name]]) && !inherits(schema[[name]], "integer64"), logical(1))]
    if (length(numeric)) where_clauses = c(where_clauses, paste0("NOT isnan(CAST(", numeric, " AS DOUBLE))"))
  }
  if (!length(where_clauses)) {
    return(from_statement)
  }
  # Wrap in subquery if from_statement contains clauses that must come after WHERE
  if (grepl("WHERE|LIMIT|ORDER\\s+BY|GROUP\\s+BY|HAVING", from_statement, ignore.case = TRUE)) {
    from_statement = glue("FROM (SELECT * {from_statement}) AS subq")
  }
  glue("
  {from_statement}
  WHERE {paste(where_clauses, collapse = ' AND ')}
  ")
}

# sql_weight_expr: returns a weight expression or NULL if no weights
sql_weight_expr = function(weights) {
  if (is.null(weights)) {
    return(NULL)
  }
  glue("1.0 * {weights}")
}

# sql_weighted_sum: SUM(expr) or SUM(w * expr) with alias
sql_weighted_sum = function(expr, weights_expr, alias) {
  if (is.null(weights_expr)) {
    return(glue("SUM({expr}) AS {alias}"))
  }
  glue("SUM(({weights_expr}) * ({expr})) AS {alias}")
}

# sql_weighted_mean: AVG(expr) or weighted mean with alias
sql_weighted_mean = function(expr, weights_expr, alias) {
  if (is.null(weights_expr)) {
    return(glue("AVG({expr}) AS {alias}"))
  }
  glue("SUM(({weights_expr}) * ({expr})) / SUM({weights_expr}) AS {alias}")
}

# build_weighted_moment_terms: aggregate terms for weighted sufficient stats.
# Aliases use variable indices (sum_w{i}, sum_w{i}_y, sum_w{i}_{j} for i <= j)
# and are read back into a matrix by cross_moments(). Set intercept = FALSE to
# skip the linear sums (e.g., for demeaned data).
build_weighted_moment_terms = function(
  y_sql,
  x_sql = character(0),
  weights_expr = NULL,
  intercept = TRUE,
  prefix_terms = NULL,
  include_w_sq = FALSE,
  w_sq_alias_base = "sum_w2",
  conn
) {
  # Cast operands, not their products: integer/DECIMAL multiplication can
  # overflow before an outer SUM or cast is evaluated.
  y_sql = sql_cast_double(y_sql, conn)
  x_sql = sql_cast_double(x_sql, conn)
  if (!is.null(weights_expr)) weights_expr = sql_cast_double(weights_expr, conn)
  weights_sq_expr = if (is.null(weights_expr)) NULL else glue("({weights_expr}) * ({weights_expr})")
  moment_terms = c(
    prefix_terms,
    if (is.null(weights_expr)) "COUNT(*) AS sum_w" else glue("SUM({weights_expr}) AS sum_w"),
    sql_weighted_sum(y_sql, weights_expr, "sum_wy"),
    sql_weighted_sum(glue("({y_sql}) * ({y_sql})"), weights_expr, "sum_wy_sq")
  )

  if (isTRUE(include_w_sq)) {
    moment_terms = c(
      moment_terms,
      if (is.null(weights_sq_expr)) {
        glue("COUNT(*) AS {w_sq_alias_base}")
      } else {
        glue("SUM({weights_sq_expr}) AS {w_sq_alias_base}")
      },
      sql_weighted_sum(y_sql, weights_sq_expr, paste0(w_sq_alias_base, "y")),
      sql_weighted_sum(glue("({y_sql}) * ({y_sql})"), weights_sq_expr, paste0(w_sq_alias_base, "y_sq"))
    )
  }

  for (i in seq_along(x_sql)) {
    moment_terms = c(
      moment_terms,
      if (intercept) sql_weighted_sum(x_sql[i], weights_expr, sprintf("sum_w%d", i)),
      sql_weighted_sum(glue("({x_sql[i]}) * ({y_sql})"), weights_expr, sprintf("sum_w%d_y", i))
    )
    for (j in seq_len(i)) {
      moment_terms = c(
        moment_terms,
        sql_weighted_sum(glue("({x_sql[j]}) * ({x_sql[i]})"), weights_expr, sprintf("sum_w%d_%d", j, i))
      )
    }
  }

  moment_terms
}

#' Check if a two-way panel is balanced
#' @keywords internal
dbreg_is_balanced_panel = function(inputs) {
  env2env(inputs, environment(), c("conn", "from_statement", "fe"))
  if (length(fe) != 2) {
    return(NA)
  }
  fe_expr = paste(fe, collapse = ", ")
  balance_sql = glue(
    "SELECT COUNT(DISTINCT cnt) AS n_distinct_counts, COUNT(*) AS n_cells, ",
    "(COUNT(DISTINCT {fe[1]}) * COUNT(DISTINCT {fe[2]})) AS n_expected ",
    "FROM (SELECT COUNT(*) AS cnt, {fe[1]}, {fe[2]} {from_statement} GROUP BY {fe_expr}) t"
  )
  res = tryCatch(dbGetQuery(conn, balance_sql), error = function(e) NULL)
  if (is.null(res)) {
    return(NA)
  }
  res[["n_distinct_counts"]] == 1 && res[["n_cells"]] == res[["n_expected"]]
}

# AP implementation: "window" (SQL window sweeps) or "join" (group-mean
# table joined back per FE). Windows are the default on Spark/Databricks,
# where each table write is costly; override with options(dbreg.ap_method).
ap_method = function(backend) {
  method = getOption("dbreg.ap_method", if (backend == "spark") "window" else "join")
  if (!identical(method, "window") && !identical(method, "join")) {
    stop('[dbreg] dbreg.ap_method must be "window" or "join".', call. = FALSE)
  }
  method
}

# SQL for `cycles` full AP cycles using window functions. Each sweep subtracts
# the weighted mean of every `_tilde` column within one FE, using the preceding
# sweep's result. One subquery per sweep, because SQL cannot nest window aggregates.
sql_ap_window_sweeps = function(table, keep_cols, tilde_names, fe, cycles) {
  sql = paste0("SELECT * FROM ", table)
  for (fe_k in rep(fe, cycles)) {
    projected = sprintf(
      "%1$s - SUM(__w * %1$s) OVER (PARTITION BY %2$s) / SUM(__w) OVER (PARTITION BY %2$s) AS %1$s",
      tilde_names, fe_k
    )
    sql = paste0("SELECT ", paste(c(keep_cols, projected), collapse = ", "), " FROM (", sql, ") ap")
  }
  sql
}

#' Alternating projections (AP) for exact multi-way FE demeaning
#' @keywords internal
dbreg_alternating_projections = function(
  inputs,
  xvars_sql,
  xvar_names,
  max_iter = if (is.null(inputs[["fe_max_iter"]])) getOption("dbreg.ap_max_iter", 100L) else inputs[["fe_max_iter"]],
  tol = if (is.null(inputs[["fe_tol"]])) getOption("dbreg.ap_tol", 1e-10) else inputs[["fe_tol"]]
) {
  env2env(inputs, environment(), c(
    "conn", "from_statement", "fe", "yvar", "weights", "cluster_var",
    "verbose"
  ))
  positive = function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0
  if (!positive(max_iter) || max_iter %% 1 != 0) {
    stop("[dbreg] dbreg.ap_max_iter must be a positive integer.", call. = FALSE)
  }
  if (!positive(tol)) {
    stop("[dbreg] dbreg.ap_tol must be a positive number.", call. = FALSE)
  }
  backend_info = detect_backend(conn)
  backend = backend_info[["name"]]
  if (length(unlist(inputs[["fe_slopes"]]))) return(dbreg_varying_projections(inputs, xvars_sql, xvar_names, max_iter, tol))
  solver = getOption("dbreg.ap_solver", "auto")
  if (!is.character(solver) || length(solver) != 1L || is.na(solver) || !solver %in% c("auto", "moments", "iterative")) {
    stop("dbreg.ap_solver must be auto, moments or iterative.", call. = FALSE)
  }
  inputs[["projection_limits"]] = projection_limits()
  if (backend == "duckdb") {
    method = getOption("dbreg.duckdb_ap_method", "join")
    if (!is.character(method) || length(method) != 1L || is.na(method) || !method %in% c("join", "cg")) {
      stop("dbreg.duckdb_ap_method must be join or cg.", call. = FALSE)
    }
    if (method == "cg" && length(fe) == 2L && solver != "moments") {
      attempt = projection_attempt("bounded_pcg", dbreg_moment_projections(inputs, xvars_sql, xvar_names, max_iter, tol,
        materialize = FALSE, two_family_pcg = TRUE))
      projected = record_projection_attempt(inputs, attempt)
      if (!is.null(projected)) return(projected)
    }
  }
  explicit_method = !is.null(getOption("dbreg.ap_method")) ||
    (backend == "duckdb" && !is.null(getOption("dbreg.duckdb_ap_method")))
  if (solver == "moments" || (solver == "auto" && !explicit_method)) {
    attempt = projection_attempt("bounded_moments", dbreg_moment_projections(inputs, xvars_sql, xvar_names, max_iter, tol, materialize = FALSE))
    projected = record_projection_attempt(inputs, attempt)
    if (!is.null(projected)) return(projected)
    if (solver == "auto" && !explicit_method && backend == "duckdb" && length(fe) >= 3L) {
      attempt = projection_attempt("sql_cg", dbreg_cg_projections(inputs, xvars_sql, xvar_names, tol))
      projected = record_projection_attempt(inputs, attempt)
      if (!is.null(projected)) return(projected)
    }
    if (solver == "auto" && !explicit_method && backend == "duckdb" && length(fe) == 2L) {
      attempt = projection_attempt("bounded_pcg", dbreg_moment_projections(inputs, xvars_sql, xvar_names, max_iter, tol,
        materialize = FALSE, two_family_pcg = TRUE))
      projected = record_projection_attempt(inputs, attempt)
      if (!is.null(projected)) return(projected)
    }
  }
  use_windows = ap_method(backend) == "window"
  batch_size = 1L
  if (use_windows) {
    batch_size = getOption("dbreg.ap_batch_size", 3L)
    if (!positive(batch_size) || batch_size %% 1 != 0) {
      stop("[dbreg] dbreg.ap_batch_size must be a positive integer.", call. = FALSE)
    }
  }
  weights_expr = sql_weight_expr(weights)
  if (is.null(weights_expr)) {
    weights_expr = "1.0"
  }

  # FE and cluster identifiers must survive every AP table, including when
  # they are also regressors (which are otherwise replaced by `_tilde` columns).
  id_cols = unique(c(fe, cluster_var))
  id_cols = id_cols[!is.na(id_cols) & id_cols != ""]

  seed = paste0(
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    "_",
    sprintf("%06d", sample.int(1e6, 1))
  )
  base_table = temp_table_name(paste0("dbreg_ap_base_", seed), backend)
  cur_table = temp_table_name(paste0("dbreg_ap_cur_", seed), backend)
  alt_table = temp_table_name(paste0("dbreg_ap_alt_", seed), backend)
  mean_table = temp_table_name(paste0("dbreg_ap_m_", seed), backend)

  success = FALSE
  on.exit({
    if (!success) {
      drop_table_if_exists(conn, cur_table, backend)
      drop_table_if_exists(conn, alt_table, backend)
      drop_table_if_exists(conn, mean_table, backend)
      drop_table_if_exists(conn, base_table, backend)
    }
  }, add = TRUE)

  base_select = c(
    setdiff(id_cols, c(yvar, xvar_names)), # overlaps are selected below
    sprintf("%s AS %s", yvar, yvar),
    sprintf("%s AS %s", xvars_sql, xvar_names),
    sprintf("%s AS __w", weights_expr)
  )
  base_sql = paste0("SELECT ", paste(base_select, collapse = ", "), " ", from_statement)
  create_projection_table_as(conn, base_table, base_sql, backend)

  vars_all = c(yvar, xvar_names)
  original = dbGetQuery(conn, sql_ap_original_scale(base_table, vars_all, conn))
  original_scale = vapply(original, as.numeric, numeric(1))
  roundoff = 128 * .Machine[["double.eps"]] * original_scale
  tilde_names = paste0(vars_all, "_tilde")

  # Iterate in double precision (e.g. averaging DECIMAL inputs keeps their
  # fixed scale on some backends)
  init_sql = sql_ap_initialization(base_table, vars_all, id_cols, conn)
  create_projection_table_as(conn, cur_table, init_sql, backend)

  # Windows run `batch_size` cycles per table write (fewer writes, but up to
  # batch_size - 1 extra cycles past convergence); joins check every cycle.
  iter = 0L
  accelerate = ap_accelerate()
  plain_next = FALSE
  previous_abs = Inf
  accelerated_steps = 0L
  while (iter < max_iter) {
    accelerated = accelerate && !plain_next && max_iter - iter >= 2L
    if (accelerated) {
      # Irons-Tuck: two cycles keeping X and F(X), then one extrapolation.
      iter = iter + 2L
      accelerated_steps = accelerated_steps + 1L
      if (use_windows) {
        create_projection_table_as(conn, alt_table, sql_ap_window_irons_tuck(cur_table, c(id_cols, "__w"), tilde_names, fe), backend)
        drop_table_if_exists(conn, cur_table, backend)
        tmp = cur_table
        cur_table = alt_table
        alt_table = tmp
      } else {
        carried = character()
        for (snapshot in c("ita", "itb")) {
          for (k in seq_along(fe)) {
            create_projection_table_as(conn, mean_table, sql_ap_group_means(cur_table, vars_all, fe[[k]]), backend)
            update_sql = sql_ap_join_projection(cur_table, mean_table, c(id_cols, carried), vars_all, fe[[k]],
              snapshot = if (k == 1L) snapshot else NULL)
            create_projection_table_as(conn, alt_table, update_sql, backend)
            drop_table_if_exists(conn, cur_table, backend)
            drop_table_if_exists(conn, mean_table, backend)
            tmp = cur_table
            cur_table = alt_table
            alt_table = tmp
            if (k == 1L) carried = c(carried, paste0(vars_all, "_", snapshot))
          }
        }
      }
      alpha = ap_irons_tuck_alpha(dbGetQuery(conn, sql_ap_irons_tuck_moments(cur_table, tilde_names)), length(tilde_names))
      create_projection_table_as(conn, alt_table, sql_ap_irons_tuck_step(cur_table, c(id_cols, "__w"), tilde_names, alpha), backend)
      drop_table_if_exists(conn, cur_table, backend)
      tmp = cur_table
      cur_table = alt_table
      alt_table = tmp
    } else {
      cycles = min(batch_size, max_iter - iter)
      iter = iter + cycles
      if (use_windows) {
        window_sql = sql_ap_window_sweeps(cur_table, c(id_cols, "__w"), tilde_names, fe, cycles)
        create_projection_table_as(conn, alt_table, window_sql, backend)
        drop_table_if_exists(conn, cur_table, backend)
        tmp = cur_table
        cur_table = alt_table
        alt_table = tmp
      } else {
        for (fe_k in fe) {
          mean_sql = sql_ap_group_means(cur_table, vars_all, fe_k)
          create_projection_table_as(conn, mean_table, mean_sql, backend)

          update_sql = sql_ap_join_projection(cur_table, mean_table, id_cols, vars_all, fe_k)
          create_projection_table_as(conn, alt_table, update_sql, backend)

          drop_table_if_exists(conn, cur_table, backend)
          drop_table_if_exists(conn, mean_table, backend)
          tmp = cur_table
          cur_table = alt_table
          alt_table = tmp
        }
      }
    }

    res = dbGetQuery(conn, sql_ap_convergence(cur_table, vars_all, fe))
    column_abs = vapply(res, as.numeric, numeric(1))
    max_abs = if (length(column_abs) && all(is.finite(column_abs))) max(column_abs) else Inf

    # Follow an extrapolation that did not reduce the residual with plain AP.
    plain_next = accelerated && !(max_abs < previous_abs)
    previous_abs = max_abs

    if (isTRUE(verbose)) {
      message("[AP] iter ", iter, ": max abs mean = ", sprintf("%.4e", max_abs))
    }

    # Group-mean residuals in units of max(1, column half-range), as in the
    # bounded certificates: absolute for unit-scale columns, relative for dollars.
    max_level = if (is.finite(max_abs)) max(column_abs / pmax(1, original_scale)) else Inf
    max_relative = Inf
    if (is.finite(max_level) && max_level < tol) {
      rms = vapply(dbGetQuery(conn, sql_ap_rms(cur_table, vars_all, conn)), as.numeric, numeric(1))
      # Below original-scale roundoff a column is absorbed. Otherwise require
      # its group equations to converge relative to its remaining variation.
      if (all(is.finite(rms)) && all(is.finite(roundoff))) {
        identified = rms > roundoff
        max_relative = if (any(identified)) max(column_abs[identified] / rms[identified]) else 0
      }
    }
    if (is.finite(max_level) && max_level < tol && is.finite(max_relative) && max_relative <= tol) {
      success = TRUE
      return(list(table = cur_table, base_table = base_table,
        original_scale = setNames(original_scale, vars_all),
        iterations = iter, accelerated_steps = accelerated_steps,
        max_abs = max_abs, max_relative = max_relative,
        method = if (use_windows) "window" else "row_tables"))
    }
  }

  stop(
    "[dbreg] Alternating projections did not converge within ",
    max_iter,
    " iterations (max abs mean = ",
    sprintf("%.4e", max_abs),
    "; tolerance ", format(tol), " in units of max(1, column range) and relative to projected RMS).\n\n",
    "Options:\n",
    "  - Increase tolerance: options(dbreg.ap_tol = 1e-7)\n",
    "  - Increase iterations: options(dbreg.ap_max_iter = 500)\n",
    "  - Use strategy = 'mundlak' (single-pass CRE estimator, no iteration needed)\n",
    "See ?`dbreg-options` for details.",
    call. = FALSE
  )
}

#' Choose regression strategy based on inputs and auto logic
#' @keywords internal
choose_strategy = function(inputs) {
  env2env(inputs, environment(), c(
    "strategy", "fe", "verbose", "any_continuous", "compress_ratio",
    "compress_nmax", "conn", "from_statement", "xvars", "weights", "caller"
  ))
  tag = if (is.null(caller)) "dbreg" else caller

  varying = length(unlist(inputs[["fe_slopes"]])) > 0L
  if (varying && strategy == "mundlak") stop("Varying slopes do not support the Mundlak strategy.", call. = FALSE)

  # Compression ratio estimator
  estimate_compression = function() {
    key_cols = unique(c(xvars, fe, unlist(inputs[["fe_slopes"]], use.names = FALSE)))
    # Total rows (safe: COUNT(*) is supported pretty much everywhere)
    total_sql = glue(
      "SELECT CAST(COUNT(*) AS BIGINT) AS n FROM (SELECT * {from_statement}) t"
    )
    total_n = as.numeric(dbGetQuery(conn, total_sql)[["n"]])
    if (!length(key_cols)) return(structure(1/max(total_n, 1), comp_size = 1))

    # Helper to count distinct tuples (works for single or multi-column)
    backend = detect_backend(conn)[["name"]]
    count_distinct_tuples = function(cols, limit = NULL) {
      cols_expr = paste(cols, collapse = ", ")
      # Use subquery counting distinct tuples (portable and works in DuckDB/SQL Server/etc)
      inner_sql = glue("SELECT DISTINCT {cols_expr} {from_statement}")
      if (!is.null(limit) && is.finite(limit) && limit >= 1) {
        inner_sql = sql_limit(inner_sql, limit, backend)
      }
      sql = glue("SELECT CAST(COUNT(*) AS BIGINT) AS g FROM ({inner_sql}) t")
      dbGetQuery(conn, sql)[["g"]]
    }

    if (verbose && length(fe)) {
      # count unique FE groups (may be single or multi-column)
      n_groups_fe = tryCatch(count_distinct_tuples(fe), error = function(e) {
        NA_integer_
      })
    } else {
      n_groups_fe = NA_integer_
    }

    # count unique keys over regressors + FEs (may be multi-column);
    # when not verbose, cap at the compression threshold + 1 so large
    # continuous tables do not hash every distinct row.
    key_limit = if (!verbose && is.finite(total_n)) {
      cap = floor(min(total_n * compress_ratio, compress_nmax))
      if (is.finite(cap) && cap >= 0) cap + 1 else NULL
    } else NULL
    n_groups_total = tryCatch(
      count_distinct_tuples(key_cols, limit = key_limit),
      error = function(e) NA_integer_
    )

    if (verbose) {
      data_msg = paste0(
        "        - ", 
        "data has ",
        format(total_n, big.mark = ","), " rows"
      )
      if (length(fe) && !is.na(n_groups_fe)) {
        data_msg = paste0(
          data_msg,
          " with ",
          length(fe), " FE ",
          "(", format(n_groups_fe, big.mark = ","), " unique groups)"
        )
      } else if (length(fe) == 0) {
        data_msg = paste0(data_msg, " with 0 FE")
      }
      message(data_msg)
    }

    comp_rat = n_groups_total / max(total_n, 1)
    attr(comp_rat, "comp_size") = n_groups_total

    return(comp_rat)
  }

  chosen_strategy = strategy
  est_cr = NA_real_

  # Auto logic
  if (strategy == "auto") {
    if (verbose) {
      message(sprintf("[%s] Auto strategy:", tag))
    }
    est_cr = tryCatch(estimate_compression(), error = function(e) {
      NA_real_
    })
    comp_size = attr(est_cr, "comp_size")
    fail_compress_ratio = !is.na(est_cr) && est_cr > compress_ratio
    fail_compress_nmax = !is.na(est_cr) && comp_size > compress_nmax

    if (verbose) {
      compress_ratio_msg_sign = if (fail_compress_ratio) " exceeds " else " satisfies "
      message(paste0(
        "        - ",
        "compression ratio (", sprintf("%.2f", est_cr), ")",
        compress_ratio_msg_sign,
        "threshold (", compress_ratio, ")"
      ))
      # only print compress nmax message if it fails (edge case)
      if (fail_compress_nmax) {
        compress_nmax_msg_sign = if (fail_compress_nmax) " exceeds " else " satisfies "
        message(paste0(
          "        - ",
          "compressed data size (", prettyNum(comp_size, big.mark = ","), " rows)",
          compress_nmax_msg_sign,
          "threshold (", prettyNum(compress_nmax, big.mark = ","), " rows)"
        ))
      }
    }

    if (length(fe) == 0) {
      if (verbose) {
        if (any_continuous) {
          message("        - continuous variables detected")
        }
      }
      if (any_continuous || (fail_compress_ratio || fail_compress_nmax)) {
        chosen_strategy = "moments"
      } else {
        chosen_strategy = "compress"
      }
    } else if (length(fe) >= 1) {
      if (fail_compress_ratio || fail_compress_nmax) {
        chosen_strategy = "demean"
        if (length(fe) == 2 && !varying && is.null(weights)) {
          is_balanced = dbreg_is_balanced_panel(inputs)
          if (verbose) {
            if (isTRUE(is_balanced)) {
              message("        - panel is balanced")
            } else {
              message("        - panel is unbalanced (using alternating projections)")
            }
          }
        } else if (length(fe) > 2 && verbose) {
          message("        - more than 2 FEs, using alternating projections")
        }
      } else {
        chosen_strategy = "compress"
      }
    }
    if (verbose) {
      message("        - decision: ", chosen_strategy)
    }
    } else {
    chosen_strategy = strategy
    if (verbose) {
      message(sprintf("[%s] Using strategy: ", tag), chosen_strategy)
    }
  }

  # Guard unsupported combos
  if (chosen_strategy == "moments" && length(fe) > 0) {
    warning(
      "[dbreg] FE present; moments (no-FE) not applicable. Using compress."
    )
    chosen_strategy = "compress"
  }
  if (chosen_strategy == "demean" && length(fe) == 0) {
    warning("[dbreg] demean requires at least 1 FE. Using moments.")
    chosen_strategy = "moments"
  }

  inputs[["compression_ratio_est"]] = est_cr
  if (exists("is_balanced", inherits = FALSE)) {
    inputs[["is_balanced"]] = is_balanced
  }

  chosen_strategy
}

#' Finalize dbreg result object
#' @keywords internal
finalize_dbreg_result = function(result, inputs, chosen_strategy, class = "dbreg") {
  if (inputs[["sql_only"]]) {
    if (inherits(result, "dbreg_query_plan")) print(result) else cat(result)
    return(invisible(result))
  }
  result = restore_result_names(result, inputs)
  if (inputs[["data_only"]]) {
    return(result)
  }
  result[["fe_components"]] = lapply(inputs[["fe_components"]], function(x) unname(restore_names(x, inputs[["source_map"]])))
  result[["fe_component_types"]] = lapply(inputs[["fe_components"]], function(columns) {
    at = which(vapply(inputs[["composite_fes"]], identical, logical(1), columns))
    if (length(at)) inputs[["composite_types"]][[at[1L]]] else NULL
  })
  result[["varying_slopes"]] = lapply(inputs[["fe_slopes"]], function(x) unname(restore_names(x, inputs[["source_map"]])))
  result[["composite_dictionaries"]] = inputs[["composite_dictionaries"]]
  if (!is.null(inputs[["varying_fml"]])) result[["fml"]] = Formula(rename_expr_vars(stats::formula(inputs[["varying_fml"]]), inputs[["source_map"]]))
  result = finalize_linear_layout(result, inputs)
  result = finalize_sample_diagnostics(result, inputs)
  if (!is.null(inputs[["offset"]])) {
    result[["offset"]] = inputs[["offset"]]
    result[["yvar"]] = inputs[["offset_outcome"]]
    public_formula = stats::formula(result[["fml"]])
    public_formula[[2L]] = as.name(result[["yvar"]])
    result[["fml"]] = Formula(public_formula)
  }
  result[["fixef_rm"]] = if (is.null(inputs[["fixef_rm"]])) "none" else inputs[["fixef_rm"]]
  result[["pruning_rounds"]] = if (is.null(inputs[["pruning"]])) 0L else inputs[["pruning"]][["rounds"]]
  result[["nobs_pruned"]] = if (is.null(inputs[["pruning"]])) 0 else inputs[["pruning"]][["counts"]][["singleton"]]
  result[["strategy"]] = chosen_strategy
  result[["ap_tolerance"]] = inputs[["fe_tol"]]
  result[["ap_max_iter"]] = inputs[["fe_max_iter"]]
  result[["projection_diagnostics"]] = inputs[["projection_diagnostics"]]
  diagnostics = inputs[["ap_diagnostics"]]
  for (field in names(diagnostics)) result[[paste0("ap_", field)]] = diagnostics[[field]]
  # i() columns (display names), for predict()
  result[["i_columns"]] = inputs[["i_columns"]]
  class(result) = c(class, class(result))
  finalize_linear_sample(result, inputs)
}
