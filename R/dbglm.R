#' Controls for database GLMs
#'
#' @param glm_tol Relative deviance tolerance.
#' @param score_tol Maximum absolute score scaled by the information diagonal.
#' @param glm_max_iter Maximum number of IRLS iterations.
#' @param step_max_halvings Maximum number of step halvings per iteration.
#' @param collin_tol Tolerance for structural rank selection on scaled moments.
#' @param fe_tol Relative group-equation tolerance for fixed-effect projection.
#'   Inner solves tighten this further when needed for outer score convergence.
#' @param fe_max_iter Maximum projection iterations per WLS solve.
#' @param projection_method `"auto"` permits bounded FE-summary projection for
#'   large two-way categorical effects; `"sql"` uses the existing SQL projector.
#' @return A validated list of numerical controls.
#' @export
dbglm_control = function(glm_tol = 1e-8, score_tol = 1e-8, glm_max_iter = 100L,
                          step_max_halvings = 25L, collin_tol = 1e-10,
                          fe_tol = 1e-10, fe_max_iter = 10000L,
                          projection_method = c("auto", "sql")) {
  controls = as.list(environment())
  projection_method = match.arg(projection_method)
  controls[["projection_method"]] = NULL
  for (name in names(controls)) {
    value = controls[[name]]
    if (!is.numeric(value) || length(value) != 1 || !is.finite(value) || value <= 0) {
      stop("`", name, "` must be one finite, positive number.", call. = FALSE)
    }
  }
  for (name in c("glm_max_iter", "step_max_halvings", "fe_max_iter")) {
    value = controls[[name]]
    if (value != floor(value) || value > .Machine[["integer.max"]]) {
      stop("`", name, "` must be a positive integer.", call. = FALSE)
    }
    controls[[name]] = as.integer(value)
  }
  controls[["projection_method"]] = projection_method
  controls
}

#' Fit a binomial or Poisson GLM on a database relation
#'
#' Fits binomial logit, probit, complementary-log-log and log-link Poisson models
#' using SQL-backed IRLS. Only bounded sufficient statistics are transferred to R.
#' Supports additive and joint
#' categorical fixed effects and group-specific linear/quadratic trends on native
#' DuckDB and Databricks connections, with recursive fixed-effect sample exclusions.
#'
#' @param fml Regression formula. Supports the existing SQL design syntax:
#'   column names, categorical columns, interactions, and `I(x^2)`/`I(x^3)`.
#'   Intercept-only models and numeric no-intercept models are supported. Use
#'   `factor(column)` for categorical coercion. After `|`, column names are
#'   categorical effects, `a^b` is a joint effect, and `group[t]` or
#'   `group[t,t2]` adds an intercept and group-specific continuous slopes.
#' @param conn,table,data,path Connection and input relation, with the same
#'   precedence as [dbreg()]. A default connection is in memory with spilling
#'   disabled. Caller connections retain their settings and remain open.
#' @param family A built-in `binomial("cloglog")` (default), `binomial("logit")`,
#'   `binomial("probit")`, or `poisson("log")`. Poisson responses must be
#'   nonnegative integer counts; binomial responses must be zero or one.
#' @param weights Optional column of nonnegative observation weights, including
#'   fractional weights. Zero weights are omitted. Weights multiply each row's
#'   log likelihood and score; observation counts are counts of retained rows,
#'   not sums of weights. See the likelihood conventions below.
#' @param offset Optional numeric column entering the predictor with coefficient one.
#' @param factor_refs Optional named categorical reference levels, as in [dbreg()].
#' @param vcov `"iid"`, `"hc1"`, or a one-sided formula naming one cluster column.
#' @param ssc Optional named list overriding fixest-style defaults. `K.adj=TRUE`
#'   by default; use `list(K.adj=FALSE)` for base `glm` covariance comparisons.
#'   `G.adj` controls the cluster correction. Binomial and Poisson inference is normal.
#' @param start Optional starting slope vector, including the intercept if present.
#'   Named vectors must identify the retained design coefficients.
#' @param control A list from [dbglm_control()].
#' @param verbose Print iteration diagnostics?
#' @param fixef_rm Recursive FE exclusion policy, using fixest conventions.
#'   If perfect-prediction groups remain under the selected policy, fitting errors.
#' @param row_id Optional character vector of unique, nonmissing source-key
#'   columns. Otherwise opaque snapshot keys are assigned once.
#' @param sample_table Optional new temporary table on a caller-supplied connection
#'   receiving retained keys after a successful fit. It remains caller-owned.
#' @param separation Default `"error"` withholds fits at numerical boundaries.
#'   Explicit `"prediction"` allows a no-FE boundary probability approximation,
#'   with a warning, `finite_mle=FALSE`, `converged=FALSE`, and unavailable SEs,
#'   tests and AIC/BIC. Coefficients are finite approximations, not a finite MLE.
#' @param retain Default `"none"` keeps the compact fit. `"prediction"` saves
#'   bounded portable FE dictionaries and learned coefficients for new-row
#'   prediction. It adds SQL aggregation and dictionary transfer, without keeping
#'   input rows or a connection. Oversized state errors instead of collecting
#'   the panel; `options(dbreg.prediction_max_entries = 2e6)` bounds matrix and
#'   transfer entries alongside the existing FE moment limits.
#' @return A compact object of class `c("dbglm", "dbreg")`, containing slopes,
#'   covariance, log likelihood, sample counts and convergence diagnostics. Input
#'   rows, database connections and temporary relations are not retained.
#' @details Missing/nonfinite model values and zero weights are omitted. Model
#'   specification is fixed before IRLS, and rank decisions are held constant.
#'   Fits that fail convergence or reach the numerical mean boundary
#'   error rather than returning ordinary successful inference. These checks
#'   do not constitute a complete general separation detector. Absorbed fits
#'   support new-row prediction when fitted with `retain = "prediction"`.
#'   Unknown FE keys and combinations outside the identified training span return
#'   `NA`. The saved state survives serialization and database disconnection.
#'   Grouped binomial responses and categorical
#'   no-intercept models remain unsupported.
#'   For binomial models, [logLik()] follows the `fixest::feglm()` / R binomial
#'   AIC convention: original weights are rounded to integers (ties to even)
#'   when reporting the likelihood. Fitting, deviance and covariance always use
#'   the unrounded weights. Use `logLik(model, type = "weighted")` for the
#'   weighted Bernoulli objective actually optimized. For Poisson models and
#'   integer-weight binomial models the two likelihoods coincide.
#'   If all retained weights are below one, fitting rescales them by a power
#'   of two to avoid premature convergence and loss of precision for very small
#'   weights. The divisor is recorded in `weight_scale`. Likelihood, deviance
#'   and IID covariance are returned on the original weight scale; sandwich
#'   covariances are invariant to this rescaling. Score diagnostics use the
#'   internal scale.
#'   With absorbed effects, expanded regressor columns are centered internally
#'   before forming projection moments. Their slopes, units and covariance are
#'   unchanged; the retained column centers are recorded in `design_center`.
#'   Large two-way categorical FE problems use the shared bounded sparse
#'   projector when possible. It rebuilds weighted sufficient statistics at each
#'   IRLS step and certifies normalized group equations on the transformed rows.
#'   Transfer limits use `dbreg.duckdb_ap_max_edges` and
#'   `dbreg.duckdb_ap_max_entries`. Exceeding limits or numerical certification
#'   falls back to the existing SQL projection. Use
#'   `dbglm_control(projection_method = "sql")` to select that path directly.
#'   Databricks requires a session-preserving native-query DBI connection and
#'   physical temporary-table support. Bounded categorical FE problems use
#'   parameter-space IRLS, grouped moment queries and compact global summaries.
#'   Small dense FE dictionaries use coefficient arrays; large or sparse
#'   dictionaries retain coefficient joins. Working states are inline SQL
#'   expressions on the moments path, avoiding row-sized writes during IRLS.
#'   Two-way categorical effects can also use the shared bounded sparse solver;
#'   the SQL fallback reuses batched window AP with the same normalized tolerance.
#'   `dbreg.spark_ap_chunk_cycles` controls its batch size (default three).
#'   Databricks uses double-precision SUM and a bounded count histogram (at most
#'   one million distinct counts) for the Poisson likelihood constant.
#'   `options(dbreg.glm_databricks_materialize_intermediates = TRUE)` permits
#'   diagnostic comparisons with fully materialized intermediates.
#'   FE lookups use compact JSON arrays with 17-digit double serialization.
#'   `options(dbreg.glm_databricks_lookup = "values")` selects SQL VALUES for
#'   diagnostic comparisons; both representations retain the transfer caps.
#' @examples
#' d = data.frame(y = c(0, 1, 0, 0, 1, 1), x = c(-2, -1, 0, 1, 2, 3))
#' m = dbglm(y ~ x, data = d, family = binomial("logit"))
#' coef(m)
#' predict(m, newdata = d, type = "response")
#' @export
dbglm = function(fml, conn = NULL, table = NULL, data = NULL, path = NULL,
                  family = stats::binomial("cloglog"), weights = NULL, offset = NULL,
                  factor_refs = NULL, vcov = "iid", ssc = NULL, start = NULL,
                  control = dbglm_control(), verbose = FALSE,
                  fixef_rm = c("perfect_fit", "singletons", "infinite_coef", "none"),
                  row_id = NULL, sample_table = NULL, separation = c("error", "prediction"),
                  retain = c("none", "prediction")) {
  family = dbglm_family(family)
  if (is.character(fml) && length(fml) == 1 && !is.na(fml)) fml = stats::as.formula(fml, env = baseenv())
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  if (!is.list(control) || is.null(names(control)) || anyDuplicated(names(control))) {
    stop("`control` must be a named list from dbglm_control().", call. = FALSE)
  }
  control = do.call(dbglm_control, control)
  ssc = dbglm_ssc(ssc)
  covariance = parse_vcov_args(vcov)
  fixef_rm = match.arg(fixef_rm)
  separation = match.arg(separation)
  retain = match.arg(retain)
  if (retain == "prediction") prediction_state_limit()
  if (family[["family"]] == "poisson" && separation == "prediction") {
    stop("Boundary prediction mode supports binomial families only.", call. = FALSE)
  }
  normalized = dbglm_formula(fml)
  # Replace only display names; fitted SQL/R design continues to use aliases.
  display = function(names) {
    for (a in names(normalized[["factors"]])) names = gsub(a, paste0("factor(", normalized[["factors"]][[a]], ")"), names, fixed = TRUE)
    names
  }
  parsed = parse_regression_formula(normalized[["fml"]], allow_empty = TRUE)
  unique_fe = !duplicated(parsed[["fe_labels"]])
  parsed[["fe"]] = parsed[["fe"]][unique_fe]
  parsed[["fe_labels"]] = parsed[["fe_labels"]][unique_fe]
  parsed[["fe_components"]] = parsed[["fe_components"]][unique_fe]
  fe = paste0("__dbglm_fe", seq_along(parsed[["fe"]]))
  if (!length(parsed[["fe"]])) fe = character()
  if (length(fe) && separation == "prediction") stop("Boundary prediction mode currently requires a no-FE model.", call. = FALSE)
  trend_sources = unique(unlist(normalized[["slopes"]]))
  trend_aliases = if (length(trend_sources)) stats::setNames(paste0("__dbglm_t", seq_along(trend_sources)), trend_sources) else character()
  attr(fe, "slopes") = lapply(parsed[["fe_labels"]], function(label) unname(trend_aliases[normalized[["slopes"]][[label]]]))
  validate_row_id_arg(row_id)
  if (!is.null(sample_table) && (is.null(conn) || !is.character(sample_table) || length(sample_table) != 1 || is.na(sample_table))) {
    stop("`sample_table` requires a caller connection and one new table name.", call. = FALSE)
  }
  lhs = stats::formula(parsed[["fml"]], lhs = 1, rhs = 0)[[2]]
  if (!is.symbol(lhs)) stop("The GLM outcome must be a column name.", call. = FALSE)
  for (name in c("weights", "offset")) {
    value = get(name)
    if (!is.null(value) && (!is.character(value) || length(value) != 1 || is.na(value) || !nzchar(value))) {
      stop("`", name, "` must name one column.", call. = FALSE)
    }
  }
  environment(parsed[["fml"]]) = baseenv()
  has_intercept = attr(stats::terms(stats::formula(parsed[["fml"]], lhs = 0, rhs = 1)), "intercept") == 1L
  if (!is.null(factor_refs)) {
    factor_refs = as.list(factor_refs)
    if (is.null(names(factor_refs)) || anyDuplicated(names(factor_refs)) ||
      any(!names(factor_refs) %in% c(parsed[["xvars"]], unlist(normalized[["factors"]]))) || any(lengths(factor_refs) != 1L) ||
      any(vapply(factor_refs, anyNA, logical(1)))) {
      stop("`factor_refs` must name regressor columns with one nonmissing reference each.", call. = FALSE)
    }
  }

  resources = dbglm_resources()
  on.exit(dbglm_release_resources(resources), add = TRUE)
  setup = setup_db_connection(conn, table, data, path, caller = "dbglm", resources = resources,
    duckdb_config = list(temp_directory = ""))
  conn = setup[["conn"]]
  backend_info = detect_backend(conn)
  if (backend_info[["name"]] != "duckdb" && !isTRUE(backend_info[["is_databricks"]])) {
    stop("dbglm() supports native DuckDB or Databricks connections.", call. = FALSE)
  }
  if (isTRUE(backend_info[["is_databricks"]])) validate_databricks_ap_backend(backend_info)
  if (!is.null(sample_table) && sample_table_exists(conn, sample_table)) stop("Sample export table already exists.", call. = FALSE)
  quote_id = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  tables = dbglm_table_manager(conn, resources, backend_info)
  create = tables[["create"]]
  drop = tables[["drop"]]
  variables = unique(c(parsed[["yvar"]], parsed[["xvars"]], unlist(parsed[["fe_components"]]), weights, offset, covariance[["cluster_var"]], unlist(normalized[["slopes"]])))
  real_vars = unique(c(setdiff(variables, names(normalized[["factors"]])), unlist(normalized[["factors"]]), row_id))
  if (any(startsWith(real_vars, "__dbglm_")) || any(startsWith(real_vars, "dbglm_factor_"))) stop("Source names use a reserved GLM prefix.", call. = FALSE)
  snapshot = create(paste("SELECT", paste(c(quote_id(real_vars),
    vapply(names(normalized[["factors"]]), function(a) paste(quote_id(normalized[["factors"]][[a]]), "AS", a), character(1)),
    if (isTRUE(backend_info[["is_databricks"]])) "monotonically_increasing_id() AS __dbglm_key" else "ROW_NUMBER() OVER () AS __dbglm_key"), collapse = ", "), setup[["from_statement"]]), materialize = !is.null(sample_table) && is.null(row_id))
  validate_source_row_id(conn, paste("FROM", snapshot), row_id)
  schema = dbglm_get_query(conn, paste("SELECT * FROM", snapshot, "LIMIT 0"))
  for (name in c(parsed[["yvar"]], weights, offset, unlist(normalized[["slopes"]]))) {
    if (!is.numeric(schema[[name]]) && !(name == parsed[["yvar"]] && is.logical(schema[[name]]))) {
      stop("Outcome, weights and offset must be numeric (the outcome may be logical).", call. = FALSE)
    }
  }
  eligible = paste0(quote_id(variables), " IS NOT NULL")
  numeric_vars = variables[vapply(schema[variables], is.numeric, logical(1))]
  eligible = c(eligible, paste0("isfinite(", quote_id(numeric_vars), ")"))
  if (!length(numeric_vars)) eligible = paste0(quote_id(variables), " IS NOT NULL")
  if (!is.null(weights)) eligible = c(eligible, paste(quote_id(weights), "> 0"))
  elig_cond = paste(eligible, collapse = " AND ")
  outcome = quote_id(parsed[["yvar"]])
  invalid_y = if (family[["family"]] == "poisson") paste0("(", outcome, "<0 OR ", outcome, "!=FLOOR(", outcome, "))") else
    paste(outcome, "NOT IN (0,1)")
  invalid_message = if (family[["family"]] == "poisson") "Poisson outcomes must be nonnegative integers." else "Binomial outcomes must be zero or one."
  w_bad_sql = if (!is.null(weights)) {
    w = quote_id(weights)
    paste0("SUM(CASE WHEN ", w, " < 0 THEN 1 ELSE 0 END) AS bad_w, ",
      "SUM(CASE WHEN ", elig_cond, " AND ", w, " != FLOOR(", w, ") THEN 1 ELSE 0 END) AS fractional_w, ")
  } else "0 AS bad_w, 0 AS fractional_w, "
  snap_info = dbglm_get_query(conn, paste0(
    "SELECT COUNT(*) AS n_input, ", w_bad_sql,
    "MAX(CASE WHEN ", elig_cond, " THEN ", if (is.null(weights)) "1.0" else quote_id(weights), " END) AS weight_max, ",
    "SUM(CASE WHEN ", elig_cond, " THEN 1 ELSE 0 END) AS n_elig, ",
    "SUM(CASE WHEN (", elig_cond, ") AND (", invalid_y, ") THEN 1 ELSE 0 END) AS bad_y, ",
    "MIN(CASE WHEN ", elig_cond, " THEN ", outcome, " END) AS ymin, ",
    "MAX(CASE WHEN ", elig_cond, " THEN ", outcome, " END) AS ymax FROM ", snapshot
  ))
  count_input = as.numeric(snap_info[["n_input"]])
  if (!is.na(snap_info[["bad_w"]]) && snap_info[["bad_w"]] > 0) stop("Original weights must be nonnegative.", call. = FALSE)
  n_eligible = if (is.na(snap_info[["n_elig"]])) 0 else as.numeric(snap_info[["n_elig"]])
  if (!is.na(snap_info[["bad_y"]]) && snap_info[["bad_y"]] > 0) stop(invalid_message, call. = FALSE)
  retained = create(paste("SELECT * FROM", snapshot, "WHERE", elig_cond), materialize = FALSE)
  drop(snapshot)
  prune_mode = fixef_rm
  if (isTRUE(detect_backend(conn)[["is_databricks"]]) && length(fe) && fixef_rm != "none") {
    groups = lapply(parsed[["fe_components"]], function(g) quote_id(unique(g)))
    if (!dbglm_databricks_pruning_needed(conn, retained, groups, outcome, fixef_rm, family[["family"]])) prune_mode = "none"
  }
  # Rank small category dictionaries, not the full N-row relation once per FE.
  if (length(fe)) {
    next_retained = create(sql_rank_fixed_effects(conn, retained, parsed[["fe_components"]], fe), materialize = FALSE)
    drop(retained)
    retained = next_retained
  }
  pruning = dbglm_prune(conn, retained, fe, quote_id(parsed[["yvar"]]), prune_mode, create, drop, family[["family"]])
  retained = pruning[["table"]]
  from = paste("FROM", retained)
  if (pruning[["rounds"]] == 0L) {
    sample = list(n = n_eligible, ymin = snap_info[["ymin"]], ymax = snap_info[["ymax"]],
      weight_max = snap_info[["weight_max"]], bad = 0)
  } else {
    sample = dbglm_get_query(conn, paste("SELECT COUNT(*) AS n, MIN(", quote_id(parsed[["yvar"]]),
      ") AS ymin, MAX(", quote_id(parsed[["yvar"]]), ") AS ymax, SUM(CASE WHEN",
      invalid_y, "THEN 1 ELSE 0 END) AS bad, MAX(",
      if (is.null(weights)) "1.0" else quote_id(weights), ") AS weight_max", from))
  }
  n = as.numeric(sample[["n"]])
  if (n == 0) stop("No complete positive-weight observations remain.", call. = FALSE)
  # A binary power preserves representable input weights exactly, including
  # the half-integer ties used only for binomial likelihood reporting. Scaling
  # before the first WLS solve also protects SQL aggregate/FE certification.
  weight_max = as.numeric(sample[["weight_max"]])
  weight_scale = if (weight_max < 1) 2^floor(log2(weight_max)) else 1
  if (sample[["bad"]] > 0) stop(invalid_message, call. = FALSE)
  if (has_intercept && (family[["family"]] == "poisson" && sample[["ymax"]] == 0 ||
    family[["family"]] == "binomial" && sample[["ymin"]] == sample[["ymax"]])) {
    stop("Constant outcome with an intercept has no finite ", family[["family"]], " solution.", call. = FALSE)
  }
  for (a in names(normalized[["factors"]])) {
    source = normalized[["factors"]][[a]]
    ref = if (source %in% names(factor_refs)) factor_refs[[source]] else dbglm_get_query(conn,
      paste("SELECT", a, from, "GROUP BY", a, "ORDER BY", a, "LIMIT 1"))[[1]][1]
    factor_refs[[a]] = ref
  }
  if (length(fe) && fixef_rm != "perfect_fit" && fixef_rm != "infinite_coef") {
    pure_condition = if (family[["family"]] == "poisson") paste0("MAX(", outcome, ")=0") else
      paste0("MIN(", outcome, ")=MAX(", outcome, ")")
    pure_queries = vapply(fe, function(f) {
      paste0("SELECT 1 AS p ", from, " GROUP BY ", f, " HAVING ", pure_condition)
    }, character(1))
    pure = dbglm_get_query(conn, paste0("SELECT COUNT(*) AS n FROM (", paste(pure_queries, collapse = " UNION ALL "), ") t"))[["n"]]
    if (pure) stop("GLM fixed-effect separation remains under fixef_rm='", fixef_rm, "'; no finite solution.", call. = FALSE)
  }
  spec = build_model_spec(parsed[["fml"]], conn, from, factor_refs = factor_refs, reserved_names = variables)
  if (!has_intercept && length(spec[["factor_levels"]])) {
    stop("Categorical no-intercept designs are not implemented yet.", call. = FALSE)
  }
  has_intercept = has_intercept && !length(fe)
  public = c(if (has_intercept) "(Intercept)", spec[["columns"]][["name"]])
  x = c(if (has_intercept) "__dbglm_intercept", spec[["columns"]][["alias"]])
  if (!length(x) && !length(fe)) stop("No GLM design columns were specified.", call. = FALSE)
  prediction_keys = if (retain == "prediction" && length(fe)) prediction_key_fields(conn, parsed[["fe_components"]]) else
    list(fields = character(), columns = list())
  trend_fields = vapply(trend_sources, function(v) paste0("CAST(", quote_id(v), " AS DOUBLE) AS ", trend_aliases[[v]]), character(1))
  fields = c("__dbglm_key", fe, prediction_keys[["fields"]], trend_fields, sprintf("CAST((%s) AS DOUBLE) AS %s", quote_id(parsed[["yvar"]]), "__dbglm_y"),
    sprintf("CAST((%s) AS DOUBLE)/(%s) AS __dbglm_a", if (is.null(weights)) "1.0" else quote_id(weights), dbglm_number(weight_scale)),
    sprintf("CAST((%s) AS DOUBLE) AS __dbglm_o", if (is.null(offset)) "0.0" else quote_id(offset)),
    if (has_intercept) "CAST(1.0 AS DOUBLE) AS __dbglm_intercept",
    sprintf("CAST((%s) AS DOUBLE) AS %s", spec[["columns"]][["sql"]], spec[["columns"]][["alias"]]),
    if (!is.null(covariance[["cluster_var"]])) paste(quote_id(covariance[["cluster_var"]]), "AS __dbglm_cluster"))
  base = create(paste("SELECT", paste(fields, collapse = ", "), from))
  if (is.null(sample_table)) drop(retained)
  base_terms = c(
    if (length(x)) paste0("SUM(CASE WHEN ", paste(paste0("NOT isfinite(", x, ")"), collapse = " OR "), " THEN 1 ELSE 0 END) AS nonfinite"),
    if (length(fe)) sprintf("COUNT(DISTINCT %s) AS fe_size_%d", fe, seq_along(fe)),
    if (length(fe) && length(x)) sprintf("AVG(%s) AS center_%d", x, seq_along(x)),
    dbglm_moment_terms(x, "__dbglm_a"),
    "SUM(__dbglm_a*__dbglm_y)/SUM(__dbglm_a) AS ybar"
  )
  row = dbglm_get_query(conn, paste("SELECT", paste(base_terms, collapse = ", "), "FROM", base))
  if (length(x) && !is.na(row[["nonfinite"]]) && row[["nonfinite"]] > 0) {
    stop("GLM design contains nonfinite transformed values.", call. = FALSE)
  }
  precomputed_fe_sizes = if (length(fe)) vapply(paste0("fe_size_", seq_along(fe)), function(nm) as.numeric(row[[nm]]), numeric(1)) else integer()
  resources[["fe_sizes"]] = precomputed_fe_sizes
  design_center = stats::setNames(rep(0, length(x)), public)
  if (length(fe) && length(x)) {
    # Every FE space contains the constant. Center expanded design columns
    # before forming Schur differences, preserving slope units and inference.
    # Recompute moments on centered rows instead of subtracting nearly equal
    # raw second moments in R. No observation-sized result leaves the database.
    design_center[] = vapply(seq_along(x), function(i) as.numeric(row[[paste0("center_", i)]]), numeric(1))
    if (any(!is.finite(design_center))) stop("GLM design centering is nonfinite.", call. = FALSE)
    centered_sql = paste0("SELECT * EXCLUDE(", paste(x, collapse = ","), "),",
      paste(sprintf("(%s-(%s)) AS %s", x, dbglm_number(design_center), x), collapse = ","), " FROM ", base)
    base = if (isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates", FALSE))) create(centered_sql) else
      paste0("(", dbglm_sql(conn, centered_sql), ")")
    centered_row = dbglm_get_query(conn, paste("SELECT", paste(dbglm_moment_terms(x, "__dbglm_a"), collapse = ","), "FROM", base))
    for (name in names(centered_row)) row[[name]] = centered_row[[name]]
  }
  if (length(fe)) {
    original = if (!isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates", FALSE))) paste0("(SELECT *, __dbglm_a AS w FROM ", base, ")") else
      create(paste("SELECT *, __dbglm_a AS w FROM", base), materialize = FALSE)
    raw_H = if (length(x)) dbglm_read_information(row, length(x)) else matrix(0, 0, 0)
    raw_scales = if (length(x)) sqrt(pmax(diag(raw_H), .Machine[["double.xmin"]])) else numeric()
    projection = solve_partial_wls(conn, original, x, fe, control, create, drop, response = FALSE, cache = resources, materialize = FALSE, scales = raw_scales, need_table = FALSE)
    if (length(x)) {
      H = if (!is.null(projection[["cross_reduction"]])) raw_H - projection[["cross_reduction"]] else
        dbglm_read_information(dbglm_get_query(conn, paste("SELECT", paste(dbglm_moment_terms(x, "w"), collapse = ","), "FROM", projection[["table"]])), length(x))
      original_H = raw_H
      absorbed = diag(H) <= control[["collin_tol"]] * pmax(diag(original_H), .Machine[["double.xmin"]])
      candidates = which(!absorbed)
      keep = if (length(candidates)) candidates[dbglm_basis(H[candidates, candidates, drop = FALSE], control[["collin_tol"]])] else integer()
    } else keep = integer()
    drop(projection[["table"]])
    for (st in projection[["sub_tables"]]) drop(st)
    drop(original)
  } else keep = if (length(x)) dbglm_basis(dbglm_read_information(row, length(x)), control[["collin_tol"]]) else integer()
  all_names = public
  all_x = x
  public = public[keep]
  x = x[keep]
  beta = stats::setNames(rep(0, length(x)), public)
  if (is.null(start)) {
    if ("(Intercept)" %in% public) beta["(Intercept)"] = family[["linkfun"]](row[["ybar"]])
  } else {
    if (!is.numeric(start) || length(start) != length(beta) || any(!is.finite(start))) {
      stop("`start` must contain one finite value per retained coefficient.", call. = FALSE)
    }
    if (!is.null(names(start))) {
      if (anyDuplicated(names(start)) || !setequal(names(start), display(public))) {
        stop("Named `start` values must match retained coefficients.", call. = FALSE)
      }
      start = start[display(public)]
    }
    beta[] = start
  }
  counts = dbglm_fe_counts(conn, base, fe, length(x), covariance[["vcov_type"]], ssc, parsed[["fe_labels"]], covariance[["cluster_var"]], sizes = precomputed_fe_sizes)
  if (length(fe)) {
    fit = dbglm_fe_fit(conn, base, x, fe, beta, family, control, create, drop, verbose, cache = resources, ybar = row[["ybar"]], vcov_type = covariance[["vcov_type"]])
    covariance_result = dbglm_fe_covariance(conn, x, fit, covariance[["vcov_type"]], ssc, n, counts)
  } else {
    fit = dbglm_irls(conn, base, x, beta, family, control, verbose, separation)
    covariance_result = if (isFALSE(fit[["finite_mle"]])) list(matrix = matrix(NA_real_, length(x), length(x)), clusters = NULL) else
      dbglm_covariance(conn, base, x, fit, family, covariance[["vcov_type"]], ssc, n)
  }
  public = display(public)
  all_names = display(all_names)
  V = covariance_result[["matrix"]]
  if (covariance[["vcov_type"]] == "iid") V = V / weight_scale
  if (!isFALSE(fit[["finite_mle"]]) && (any(!is.finite(V)) || any(diag(V) < 0))) {
    stop("GLM covariance is invalid; check identification and weight scale.", call. = FALSE)
  }
  dimnames(V) = list(public, public)
  beta = stats::setNames(fit[["beta"]], public)
  se = sqrt(diag(V))
  z = beta / se
  coeftable = cbind(estimate = beta, std.error = se, statistic = z,
    p.value = 2 * stats::pnorm(-abs(z)))
  full = stats::setNames(rep(NA_real_, length(all_names)), all_names)
  full[public] = beta
  saved_formula = Formula::Formula(stats::as.formula(paste(deparse(Formula::Formula(fml), width.cutoff = 500L), collapse = ""), env = baseenv()))
  # A do.call() call can contain literal data frames, connections and formulas
  # whose environments capture rows. Record resolved compact metadata instead.
  # Build arguments as a list: removing an absent optional argument with
  # [[<- NULL is safe for lists, but not for language objects.
  call = as.list(call("dbglm", fml = saved_formula))
  call[["family"]] = if (family[["family"]] == "poisson") quote(stats::poisson("log")) else
    substitute(stats::binomial(link), list(link = family[["link"]]))
  if (!setup[["own_conn"]]) call[["conn"]] = as.name(".connection")
  if (!is.null(data)) call[["data"]] = as.name(".data")
  if (!is.null(table)) call[["table"]] = if (is.character(table)) as.character(table) else as.name(".relation")
  if (!is.null(path)) call[["path"]] = as.character(path)
  call[["weights"]] = weights
  call[["offset"]] = offset
  call[["vcov"]] = if (covariance[["vcov_type"]] == "cluster") stats::reformulate(covariance[["cluster_var"]], env = baseenv()) else covariance[["vcov_type"]]
  call[["ssc"]] = ssc
  call[["control"]] = control
  call[["fixef_rm"]] = fixef_rm
  call[["separation"]] = separation
  call[["row_id"]] = if (length(row_id)) as.character(row_id) else NULL
  call[["sample_table"]] = sample_table
  call[["retain"]] = retain
  call = as.call(call)
  saturated_loglik = if (family[["family"]] == "poisson") dbglm_saturated_loglik(conn, base) else 0
  weighted_loglik = (saturated_loglik - fit[["deviance"]] / 2) * weight_scale
  reported_loglik = weighted_loglik
  if (family[["family"]] == "binomial" && isTRUE(snap_info[["fractional_w"]] > 0)) {
    working = if (length(fe)) fit[["working"]] else paste0("(",
      dbglm_working_sql(base, x, fit[["beta"]], family[["link"]]), "SELECT * FROM working)")
    reported_loglik = dbglm_binomial_loglik(conn, working, weight_scale)
  }
  result = list(coeftable = coeftable, vcov = V, collin.coef = full,
    collin.var = setdiff(all_names, public), coef_names = public, coef_names_all = all_names,
    family = family, fml = saved_formula, model_spec = spec, has_intercept = has_intercept,
    yvar = parsed[["yvar"]], xvars = public, offset = offset, weights = weights,
    nobs_orig = n, nobs_input = as.numeric(count_input), nobs_dropped = as.numeric(count_input) - n,
    n_params_model = counts[["model"]], n_params_vcov = counts[["covariance"]], n_slopes = sum(public != "(Intercept)"),
    df_residual = n - counts[["model"]], vcov_type = covariance[["vcov_type"]],
    cluster_var = covariance[["cluster_var"]], n_clusters = covariance_result[["clusters"]], ssc = ssc,
    loglik = reported_loglik, weighted_loglik = weighted_loglik, weight_scale = weight_scale,
    deviance = fit[["deviance"]] * weight_scale, converged = !isFALSE(fit[["converged"]]),
    finite_mle = !isFALSE(fit[["finite_mle"]]), separation = separation,
    iterations = fit[["iterations"]], score_max = fit[["score_max"]],
    deviance_history = fit[["deviance_history"]] * weight_scale, step_halvings = fit[["step_halvings"]],
    control = control, strategy = "irls", fe = if (length(fe)) display(parsed[["fe_labels"]]) else NULL,
    n_fe_levels = stats::setNames(counts[["sizes"]], display(parsed[["fe_labels"]])), fe_nested = counts[["nested"]],
    design_center = stats::setNames(unname(design_center[keep]), public),
    nobs_eligible = n_eligible, nobs_pruned = n_eligible - n,
    pruning_rounds = pruning[["rounds"]], exclusions = pruning[["counts"]], fixef_rm = fixef_rm,
    factor_sources = normalized[["factors"]], design_names = c(if (has_intercept) "(Intercept)", spec[["columns"]][["name"]])[keep],
    varying_slopes = normalized[["slopes"]], row_id = row_id, sample_table = sample_table, ap_cycles = fit[["ap_cycles"]], ap_residual = fit[["ap_residual"]],
    projection_iterations = fit[["ap_cycles"]], projection_residual = fit[["ap_residual"]],
    projection_method = if (!is.null(fit[["projection_method"]])) fit[["projection_method"]] else if (length(fe)) "group" else "none", call = call)
  if (retain == "prediction" && length(fe)) {
    result[["fitted_state"]] = dbglm_prediction_state(conn, fit[["working"]], base,
      fe, parsed[["fe_components"]], all_x, spec[["columns"]][["name"]],
      design_center, trend_aliases, counts[["sizes"]], create, drop, n,
      dictionary_components = prediction_keys[["columns"]])
  }
  if (!is.null(sample_table)) {
    keys = if (length(row_id)) row_id else "__dbglm_key"
    export_retained_sample(conn, sample_table, keys, from)
  }
  class(result) = c("dbglm", "dbreg")
  result
}
