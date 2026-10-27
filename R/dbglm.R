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

#' Fit a binomial, Poisson or Gaussian GLM on a database relation
#'
#' Fits binomial logit, probit, complementary-log-log and log-link Poisson models
#' using SQL-backed IRLS. Gaussian identity models use one certified WLS update.
#' Only bounded sufficient statistics are transferred to R.
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
#'   disabled. `data` or `table` may also be a `dbplyr` lazy table, whose
#'   connection is inferred. Caller connections retain their settings and remain open.
#' @param family A supported family name (`"gaussian"`, `"binomial"`, or
#'   `"poisson"`), a family function, or a family object. The default is
#'   `gaussian("identity")`. Supported alternatives are `binomial("cloglog")`,
#'   `binomial("logit")`, `binomial("probit")`, and `poisson("log")`. Poisson responses must be
#'   nonnegative integer counts; binomial responses must be zero or one.
#' @param weights Optional column name or one-sided formula such as `~w`,
#'   giving nonnegative observation weights, including
#'   fractional weights. Zero weights are omitted. Weights multiply each row's
#'   log likelihood and score; observation counts are counts of retained rows,
#'   not sums of weights. See the likelihood conventions below.
#' @param offset Optional numeric column name or one-sided formula such as
#'   `~o`, entering the predictor with coefficient one.
#' @param factor_refs Optional named categorical reference levels, as in [dbreg()].
#' @param vcov `"iid"`, `"hc1"`, or a one-sided formula naming one cluster column.
#' @param ssc Small-sample policy, as in [dbreg()]. Accepts `"full"`,
#'   `"nested"`, `"none"`, or a named correction-policy list. `NULL` preserves
#'   the GLM default of nonnested FE intercept counting. `K.adj=TRUE` by default;
#'   use `list(K.adj=FALSE)` for base `glm` covariance comparisons. Binomial and
#'   Poisson inference uses the normal distribution.
#' @param start Optional starting slope vector, including the intercept if present.
#'   Named vectors must identify the retained design coefficients.
#' @param control A list from [dbglm_control()].
#' @param strategy `"auto"` selects exact cell compression when its cardinality
#'   satisfies `compress_ratio` and `compress_nmax`. `"compress"` always uses
#'   cells; `"moments"` uses original-row IRLS sufficient statistics. `"ap"`
#'   requires fixed effects and selects SQL alternating projections.
#' @param compress_ratio,compress_nmax Compression thresholds, as in [dbreg()].
#'   Compression groups the expanded regressors, fixed effects, offsets and
#'   initial cluster. It preserves original observation counts, likelihoods,
#'   IID and sandwich covariance, and requested original-row retention.
#' @param verbose Print iteration diagnostics?
#' @param fixef_rm Recursive FE exclusion policy, using fixest conventions.
#'   If perfect-prediction groups remain under the selected policy, fitting errors.
#'   Gaussian `"perfect_fit"` removes singletons; `"infinite_coef"` removes none.
#' @param subset Optional logical column or restricted one-sided predicate
#'   formula, as in [dbreg()]. Original row identities precede this filter.
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
#'   input rows. `options(dbreg.prediction_max_entries = 2e6)` bounds matrix and
#'   transfer entries alongside the existing FE moment limits. One or two FEs
#'   without varying slopes can keep oversized dictionaries in SQL instead.
#'   That state needs its live owning session, cannot be restored with
#'   `saveRDS()`, and can be released with [close()]. Other oversized designs
#'   error instead of collecting the panel. Two-FE recovery uses the fitted
#'   projection iteration limit and a tolerance at least as strict as `1e-10`.
#'   `"all"` also saves original rows, fitted predictors, residuals and scores
#'   in an owned database relation, enabling default/sample prediction, [stats::case.names()],
#'   [residuals()] and post-fit covariance via [summary()]. It requires a finite
#'   fit, adds row-sized backend storage/SQL work and needs its live owning session.
#'   SQL inputs return keyed relations; R inputs preserve original positions.
#'   `"inference"` saves rows/scores without learned-effect dictionaries. It
#'   supports fitted-sample prediction and post-fit inference; learned-effect
#'   newdata prediction needs `"prediction"` or `"all"`.
#' @return A compact object of class `c("dbglm", "dbreg")`, containing slopes,
#'   covariance, log likelihood, sample counts and convergence diagnostics. Compact
#'   modes do not retain input rows. Database-backed retention keeps its owned
#'   dictionary and live connection until [close()] or garbage collection.
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
#'   Gaussian models estimate dispersion using weighted RSS divided by model
#'   residual df. Their IID covariance and t tests follow `fixest::feglm()`,
#'   including its additional IID SSC factor. Gaussian default likelihood uses
#'   R's Gaussian AIC convention (profile likelihood minus one); `type="weighted"`
#'   returns the profile likelihood. Likelihood df excludes the dispersion,
#'   matching fixest. These conventions differ from OLS [dbreg()].
#'   If all retained weights are below one, fitting rescales them by a power
#'   of two to avoid premature convergence and loss of precision for very small
#'   weights. The divisor is recorded in `weight_scale`. Likelihood, deviance
#'   and IID covariance are returned on the original weight scale; sandwich
#'   covariances are invariant to this rescaling. Score diagnostics use the
#'   internal scale. All families normalize uniformly large weights. Gaussian models
#'   combine dispersion and bread on the internal scale before restoring units.
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
                 family = stats::gaussian(), weights = NULL, offset = NULL,
                 factor_refs = NULL, vcov = "iid", ssc = NULL, start = NULL,
                 control = dbglm_control(), verbose = FALSE,
                 fixef_rm = c("perfect_fit", "singletons", "infinite_coef", "none"),
                 row_id = NULL, sample_table = NULL, separation = c("error", "prediction"),
                 retain = c("none", "prediction", "inference", "all"), subset = NULL,
                 strategy = c("auto", "compress", "moments", "ap"),
                 compress_ratio = NULL, compress_nmax = 1e6) {
  inputs = new.env(parent = emptyenv())
  env2env(environment(), inputs, names(formals(dbglm)))
  inputs[["resources"]] = fit_resources()
  on.exit(release_fit_resources(inputs[["resources"]]), add = TRUE)
  process_dbglm_inputs(inputs)
  execute_dbglm_fit(inputs)
}

process_dbglm_inputs = function(inputs) {
  env2env(inputs, environment(), c(names(formals(dbglm)), "resources"))
  weights = formula_column(weights, "weights")
  offset = formula_column(offset, "offset")
  family = dbglm_family(family)
  subset_variables(subset)
  if (is.character(fml) && length(fml) == 1 && !is.na(fml)) fml = stats::as.formula(fml, env = baseenv())
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  if (!is.list(control) || is.null(names(control)) || anyDuplicated(names(control))) {
    stop("`control` must be a named list from dbglm_control().", call. = FALSE)
  }
  control = do.call(dbglm_control, control)
  ssc = dbglm_ssc(ssc)
  covariance = parse_vcov_args(vcov)
  fixef_rm = match.arg(fixef_rm, c("perfect_fit", "singletons", "infinite_coef", "none"))
  separation = match.arg(separation, c("error", "prediction"))
  retain = match.arg(retain, c("none", "prediction", "inference", "all"))
  strategy = match.arg(strategy, c("auto", "compress", "moments", "ap"))
  if (!is.numeric(compress_nmax) || length(compress_nmax) != 1L ||
    is.na(compress_nmax) || compress_nmax < 1 || compress_nmax != floor(compress_nmax)) {
    stop("`compress_nmax` must be one positive whole number or Inf.", call. = FALSE)
  }
  if (retain_prediction(retain)) prediction_state_limit()
  if (family[["family"]] != "binomial" && separation == "prediction") {
    stop("Boundary prediction mode supports binomial families only.", call. = FALSE)
  }
  normalized = dbglm_formula(fml)
  # Replace only display names; fitted SQL/R design continues to use aliases.
  display = function(names) {
    names = as.character(names)
    i_columns = rep(FALSE, length(names))
    for (a in normalized[["i_aliases"]]) {
      selected = startsWith(names, paste0(a, "::"))
      i_columns = i_columns | selected
      names[selected] = paste0(normalized[["factor_display"]][[a]], substring(names[selected], nchar(a) + 1L))
    }
    for (a in setdiff(names(normalized[["factors"]]), normalized[["i_aliases"]])) {
      names[!i_columns] = gsub(a, normalized[["factor_display"]][[a]], names[!i_columns], fixed = TRUE)
    }
    names
  }
  parsed = parse_regression_formula(normalized[["fml"]], allow_empty = TRUE, allow_i = TRUE)
  unique_fe = !duplicated(parsed[["fe_labels"]])
  parsed[["fe"]] = parsed[["fe"]][unique_fe]
  parsed[["fe_labels"]] = parsed[["fe_labels"]][unique_fe]
  parsed[["fe_components"]] = parsed[["fe_components"]][unique_fe]
  fe = paste0("__dbglm_fe", seq_along(parsed[["fe"]]))
  if (!length(parsed[["fe"]])) fe = character()
  if (length(fe) && separation == "prediction") stop("Boundary prediction mode currently requires a no-FE model.", call. = FALSE)
  if (strategy == "ap" && !length(fe)) stop("The ap strategy requires fixed effects.", call. = FALSE)
  if (strategy == "ap") control[["projection_method"]] = "sql"
  compress_ratio = resolve_compress_ratio(compress_ratio, fe)
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

  retained_data = prepare_retained_data(data, table, retain)
  data = retained_data[["data"]]
  setup = setup_db_connection(conn, table, data, path,
    caller = "dbglm", resources = resources,
    duckdb_config = list(temp_directory = "")
  )
  conn = setup[["conn"]]
  data = setup[["data"]]
  table = setup[["table"]]
  retained_source = if (retain_sample(retain)) materialize_retained_source(
    conn, setup[["from_statement"]], resources,
    retained_data[["key"]], retained_data[["r_input"]]
  ) else NULL
  source_from = if (is.null(retained_source)) setup[["from_statement"]] else paste("FROM", DBI::dbQuoteIdentifier(conn, retained_source[["table"]]))
  resources[["retain_scores"]] = retain_sample(retain)
  backend_info = detect_backend(conn)
  if (backend_info[["name"]] != "duckdb" && !isTRUE(backend_info[["is_databricks"]])) {
    stop("dbglm() supports native DuckDB or Databricks connections.", call. = FALSE)
  }
  if (isTRUE(backend_info[["is_databricks"]])) validate_databricks_ap_backend(backend_info)
  if (!is.null(sample_table) && sample_table_exists(conn, sample_table)) stop("Sample export table already exists.", call. = FALSE)
  quote_id = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  tables = fit_table_manager(conn, resources, backend_info)
  create = tables[["create"]]
  drop = tables[["drop"]]
  variables = unique(c(parsed[["yvar"]], parsed[["xvars"]], unlist(parsed[["fe_components"]]), weights, offset, covariance[["cluster_var"]], unlist(normalized[["slopes"]])))
  real_vars = unique(c(setdiff(variables, names(normalized[["factors"]])), unlist(normalized[["factors"]]), row_id, subset_variables(subset)))
  if (any(startsWith(real_vars, "__dbglm_")) || any(startsWith(real_vars, "dbglm_factor_"))) stop("Source names use a reserved GLM prefix.", call. = FALSE)
  snapshot = create(paste("SELECT", paste(c(
    quote_id(real_vars),
    vapply(names(normalized[["factors"]]), function(a) paste(quote_id(normalized[["factors"]][[a]]), "AS", a), character(1)),
    if (!is.null(retained_source)) paste(quote_id(retained_source[["key"]]), "AS __dbglm_key") else if (isTRUE(backend_info[["is_databricks"]])) "monotonically_increasing_id() AS __dbglm_key" else "ROW_NUMBER() OVER () AS __dbglm_key"
  ), collapse = ", "), source_from), materialize = !is.null(sample_table) && is.null(row_id))
  validate_source_row_id(conn, paste("FROM", snapshot), row_id)
  schema = backend_get_query(conn, paste("SELECT * FROM", snapshot, "LIMIT 0"))
  for (name in c(parsed[["yvar"]], weights, offset, unlist(normalized[["slopes"]]))) {
    if (!is.numeric(schema[[name]]) && !(name == parsed[["yvar"]] && is.logical(schema[[name]]))) {
      stop("Outcome, weights and offset must be numeric (the outcome may be logical).", call. = FALSE)
    }
  }
  eligible = paste0(quote_id(variables), " IS NOT NULL")
  numeric_vars = variables[vapply(schema[variables], is.numeric, logical(1))]
  eligible = c(eligible, sql_is_finite(conn, quote_id(numeric_vars), backend_info))
  if (!length(numeric_vars)) eligible = paste0(quote_id(variables), " IS NOT NULL")
  if (!is.null(subset)) eligible = c(eligible, sql_subset_predicate(conn, subset, schema))
  if (!is.null(weights)) eligible = c(eligible, paste(quote_id(weights), "> 0"))
  elig_cond = paste(eligible, collapse = " AND ")
  outcome = quote_id(parsed[["yvar"]])
  invalid_y = if (family[["family"]] == "poisson") paste0("(", outcome, "<0 OR ", outcome, "!=FLOOR(", outcome, "))") else if (family[["family"]] == "gaussian") "FALSE" else paste(outcome, "NOT IN (0,1)")
  invalid_message = if (family[["family"]] == "poisson") "Poisson outcomes must be nonnegative integers." else "Binomial outcomes must be zero or one."
  w_bad_sql = if (!is.null(weights)) {
    w = quote_id(weights)
    paste0(
      "SUM(CASE WHEN ", w, " < 0 THEN 1 ELSE 0 END) AS bad_w, ",
      "SUM(CASE WHEN ", elig_cond, " AND ", w, " != FLOOR(", w, ") THEN 1 ELSE 0 END) AS fractional_w, "
    )
  } else "0 AS bad_w, 0 AS fractional_w, "
  snap_info = backend_get_query(conn, paste0(
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
  if (family[["family"]] == "gaussian") prune_mode = if (fixef_rm %in% c("perfect_fit", "singletons")) "singletons" else "none"
  if (length(fe) && prune_mode != "none") {
    groups = lapply(parsed[["fe_components"]], function(g) quote_id(unique(g)))
    if (!dbglm_pruning_needed(conn, retained, groups, outcome, prune_mode, family[["family"]])) prune_mode = "none"
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
    sample = list(
      n = n_eligible, ymin = snap_info[["ymin"]], ymax = snap_info[["ymax"]],
      weight_max = snap_info[["weight_max"]], bad = 0
    )
  } else {
    sample = backend_get_query(conn, paste(
      "SELECT COUNT(*) AS n, MIN(", quote_id(parsed[["yvar"]]),
      ") AS ymin, MAX(", quote_id(parsed[["yvar"]]), ") AS ymax, SUM(CASE WHEN",
      invalid_y, "THEN 1 ELSE 0 END) AS bad, MAX(",
      if (is.null(weights)) "1.0" else quote_id(weights), ") AS weight_max", from
    ))
  }
  n = as.numeric(sample[["n"]])
  if (n == 0) {
    stop("All observations were removed by fixed-effect singleton/separation pruning.", call. = FALSE)
  }
  # A binary power preserves representable input weights exactly, including
  # the half-integer ties used only for binomial likelihood reporting. Scaling
  # before the first WLS solve also protects SQL aggregate/FE certification.
  weight_max = as.numeric(sample[["weight_max"]])
  weight_scale = if (weight_max < 1 || weight_max > 2) 2^floor(log2(weight_max)) else 1
  if (sample[["bad"]] > 0) stop(invalid_message, call. = FALSE)
  if (has_intercept && (family[["family"]] == "poisson" && sample[["ymax"]] == 0 ||
    family[["family"]] == "binomial" && sample[["ymin"]] == sample[["ymax"]])) {
    stop("Constant outcome with an intercept has no finite ", family[["family"]], " solution.", call. = FALSE)
  }
  for (a in setdiff(names(normalized[["factors"]]), normalized[["i_aliases"]])) {
    source = normalized[["factors"]][[a]]
    ref = if (source %in% names(factor_refs)) factor_refs[[source]] else backend_get_query(
      conn,
      paste("SELECT", a, from, "GROUP BY", a, "ORDER BY", a, "LIMIT 1")
    )[[1]][1]
    factor_refs[[a]] = ref
  }
  if (family[["family"]] != "gaussian" && length(fe) && fixef_rm != "perfect_fit" && fixef_rm != "infinite_coef") {
    pure_condition = if (family[["family"]] == "poisson") paste0("MAX(", outcome, ")=0") else
      paste0("MIN(", outcome, ")=MAX(", outcome, ")")
    pure_queries = vapply(fe, function(f) {
      paste0("SELECT 1 AS p ", from, " GROUP BY ", f, " HAVING ", pure_condition)
    }, character(1))
    pure = backend_get_query(conn, paste0("SELECT COUNT(*) AS n FROM (", paste(pure_queries, collapse = " UNION ALL "), ") t"))[["n"]]
    if (pure) stop("GLM fixed-effect separation remains under fixef_rm='", fixef_rm, "'; no finite solution.", call. = FALSE)
  }
  spec = build_model_spec(parsed[["fml"]], conn, from, factor_refs = factor_refs, reserved_names = variables, absorbed = length(fe) > 0L)
  if (!has_intercept && length(setdiff(names(spec[["factor_levels"]]), names(spec[["i_specs"]])))) {
    stop("Categorical no-intercept designs are not implemented yet.", call. = FALSE)
  }
  has_intercept = has_intercept && !length(fe)
  public = c(if (has_intercept) "(Intercept)", spec[["columns"]][["name"]])
  if (anyDuplicated(display(public))) stop("Formula terms produce duplicate public coefficient names.", call. = FALSE)
  x = c(if (has_intercept) "__dbglm_intercept", spec[["columns"]][["alias"]])
  if (!length(x) && !length(fe)) stop("No GLM design columns were specified.", call. = FALSE)
  prediction_keys = if (retain_prediction(retain) && length(fe)) prediction_key_fields(conn, parsed[["fe_components"]]) else
    list(fields = character(), columns = list())
  trend_fields = vapply(trend_sources, function(v) paste0("CAST(", quote_id(v), " AS DOUBLE) AS ", trend_aliases[[v]]), character(1))
  fields = c(
    "__dbglm_key", fe, prediction_keys[["fields"]], trend_fields, sprintf("CAST((%s) AS DOUBLE) AS %s", quote_id(parsed[["yvar"]]), "__dbglm_y"),
    sprintf("CAST((%s) AS DOUBLE)/(%s) AS __dbglm_a", if (is.null(weights)) "1.0" else quote_id(weights), dbglm_number(weight_scale)),
    sprintf("CAST((%s) AS DOUBLE) AS __dbglm_o", if (is.null(offset)) "0.0" else quote_id(offset)),
    if (has_intercept) "CAST(1.0 AS DOUBLE) AS __dbglm_intercept",
    sprintf("CAST((%s) AS DOUBLE) AS %s", spec[["columns"]][["sql"]], spec[["columns"]][["alias"]]),
    if (!is.null(covariance[["cluster_var"]])) paste(quote_id(covariance[["cluster_var"]]), "AS __dbglm_cluster")
  )
  base = create(paste("SELECT", paste(fields, collapse = ", "), from))
  if (is.null(sample_table)) drop(retained)
  base_terms = c(
    if (family[["family"]] == "gaussian") "SUM(LN(__dbglm_a)) AS log_weights",
    if (length(x)) paste0("SUM(CASE WHEN ", paste(paste0("NOT ", sql_is_finite(conn, x, backend_info)), collapse = " OR "), " THEN 1 ELSE 0 END) AS nonfinite"),
    if (length(fe)) sprintf("COUNT(DISTINCT %s) AS fe_size_%d", fe, seq_along(fe)),
    if (length(fe) && length(x)) sprintf("AVG(%s) AS center_%d", x, seq_along(x)),
    dbglm_moment_terms(x, "__dbglm_a"),
    "SUM(__dbglm_a*__dbglm_y)/SUM(__dbglm_a) AS ybar"
  )
  row = backend_get_query(conn, paste("SELECT", paste(base_terms, collapse = ", "), "FROM", base))
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
    centered_sql = paste0(
      "SELECT ", sql_select_except(conn, x, backend_info), ",",
      paste(sprintf("(%s-(%s)) AS %s", x, dbglm_number(design_center), x), collapse = ","), " FROM ", base
    )
    base = if (isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates", FALSE))) create(centered_sql) else
      paste0("(", backend_sql(conn, centered_sql), ")")
    centered_row = backend_get_query(conn, paste("SELECT", paste(dbglm_moment_terms(x, "__dbglm_a"), collapse = ","), "FROM", base))
    for (name in names(centered_row)) row[[name]] = centered_row[[name]]
  }
  original_base = base
  cell_groups = unique(c(
    x, fe, unlist(attr(fe, "slopes")),
    unlist(prediction_keys[["columns"]]), "__dbglm_o",
    if (!is.null(covariance[["cluster_var"]])) "__dbglm_cluster"
  ))
  compression = dbglm_cell_source(
    conn, base, cell_groups, n,
    strategy,
    family[["family"]], as.numeric(row[["ybar"]]), create, drop, weight_scale,
    compress_ratio, compress_nmax
  )
  cell = compression[["enabled"]]
  base = compression[["table"]]
  if (length(fe)) {
    original = if (!isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates", FALSE))) paste0("(SELECT *, __dbglm_a AS w FROM ", base, ")") else
      create(paste("SELECT *, __dbglm_a AS w FROM", base), materialize = FALSE)
    raw_H = if (length(x)) dbglm_read_information(row, length(x)) else matrix(0, 0, 0)
    raw_scales = if (length(x)) sqrt(pmax(diag(raw_H), .Machine[["double.xmin"]])) else numeric()
    projection = solve_partial_wls(conn, original, x, fe, control, create, drop, response = FALSE, cache = resources, materialize = FALSE, scales = raw_scales, need_table = FALSE)
    if (length(x)) {
      H = if (!is.null(projection[["cross_reduction"]])) raw_H - projection[["cross_reduction"]] else
        dbglm_read_information(backend_get_query(conn, paste("SELECT", paste(dbglm_moment_terms(x, "w"), collapse = ","), "FROM", projection[["table"]])), length(x))
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
  env2env(environment(), inputs, c(
    "public", "all_names", "beta", "display", "all_x",
    "base", "cell", "compress_nmax", "compress_ratio", "compression",
    "conn", "control", "count_input", "covariance", "create",
    "data", "design_center", "drop", "family", "fe",
    "fixef_rm", "fml", "from", "has_intercept", "keep",
    "n", "n_eligible", "normalized", "offset", "original_base",
    "parsed", "path", "precomputed_fe_sizes", "prediction_keys", "pruning",
    "resources", "retain", "retained_source", "row", "row_id",
    "sample_table", "separation", "setup", "snap_info", "spec",
    "ssc", "strategy", "subset", "table", "trend_aliases",
    "verbose", "weight_scale", "weights", "x"
  ))
  invisible(NULL)
}

execute_dbglm_fit = function(inputs) {
  env2env(inputs, environment(), c(
    "public", "all_names", "beta", "display", "all_x",
    "base", "cell", "compress_nmax", "compress_ratio", "compression",
    "conn", "control", "count_input", "covariance", "create",
    "data", "design_center", "drop", "family", "fe",
    "fixef_rm", "fml", "from", "has_intercept", "keep",
    "n", "n_eligible", "normalized", "offset", "original_base",
    "parsed", "path", "precomputed_fe_sizes", "prediction_keys", "pruning",
    "resources", "retain", "retained_source", "row", "row_id",
    "sample_table", "separation", "setup", "snap_info", "spec",
    "ssc", "strategy", "subset", "table", "trend_aliases",
    "verbose", "weight_scale", "weights", "x"
  ))
  counts = dbglm_fe_counts(conn, base, fe, length(x), covariance[["vcov_type"]], ssc, parsed[["fe_labels"]], covariance[["cluster_var"]], sizes = precomputed_fe_sizes)
  if (length(fe)) {
    fit = dbglm_fe_fit(conn, base, x, fe, beta, family, control, create, drop, verbose, cache = resources, ybar = row[["ybar"]], vcov_type = covariance[["vcov_type"]], cell = cell)
    covariance_result = dbglm_fe_covariance(conn, x, fit, covariance[["vcov_type"]], ssc, n, counts, cell = cell)
  } else {
    fit = dbglm_irls(conn, base, x, beta, family, control, verbose, separation, cell = cell)
    covariance_result = if (isFALSE(fit[["finite_mle"]])) list(matrix = matrix(NA_real_, length(x), length(x)), clusters = NULL) else
      dbglm_covariance(conn, base, x, fit, family, covariance[["vcov_type"]], ssc, n, cell = cell)
  }
  public = display(public)
  all_names = display(all_names)
  V = covariance_result[["matrix"]]
  dispersion_scaled = if (family[["family"]] == "gaussian") fit[["deviance"]] / max(n - counts[["model"]], 1) else 1
  dispersion = if (family[["family"]] == "gaussian") dispersion_scaled * weight_scale else 1
  if (covariance[["vcov_type"]] == "iid") V = V * (if (family[["family"]] == "gaussian") dispersion_scaled else 1 / weight_scale)
  if (!isFALSE(fit[["finite_mle"]]) && (any(!is.finite(V)) || any(diag(V) < 0))) {
    stop("GLM covariance is invalid; check identification and weight scale.", call. = FALSE)
  }
  dimnames(V) = list(public, public)
  beta = stats::setNames(fit[["beta"]], public)
  test_df = if (is.numeric(ssc[["t.df"]])) ssc[["t.df"]] else if (!is.null(covariance_result[["clusters"]]) && ssc[["t.df"]] == "min")
    max(covariance_result[["clusters"]] - 1, 1) else max(n - counts[["covariance"]], 1)
  full = stats::setNames(rep(NA_real_, length(all_names)), all_names)
  full[public] = beta
  saved_formula = Formula::Formula(stats::as.formula(paste(deparse(Formula::Formula(fml), width.cutoff = 500L), collapse = ""), env = baseenv()))
  # A do.call() call can contain literal data frames, connections and formulas
  # whose environments capture rows. Record resolved compact metadata instead.
  # Build arguments as a list: removing an absent optional argument with
  # [[<- NULL is safe for lists, but not for language objects.
  call = as.list(call("dbglm", fml = saved_formula))
  call[["family"]] = if (family[["family"]] == "gaussian") quote(stats::gaussian("identity")) else if (family[["family"]] == "poisson") quote(stats::poisson("log")) else
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
  call[["strategy"]] = strategy
  call[["compress_ratio"]] = compress_ratio
  call[["compress_nmax"]] = compress_nmax
  if (!is.null(subset)) call[["subset"]] = if (inherits(subset, "formula"))
    stats::as.formula(paste(deparse(subset), collapse = ""), env = baseenv()) else subset
  call = as.call(call)
  saturated_loglik = if (family[["family"]] == "poisson") dbglm_saturated_loglik(conn, original_base) else 0
  weighted_loglik = (saturated_loglik - fit[["deviance"]] / 2) * weight_scale
  reported_loglik = weighted_loglik
  if (family[["family"]] == "gaussian") {
    log_weights = if (is.null(weights)) 0 else as.numeric(row[["log_weights"]]) + n * log(weight_scale)
    weighted_loglik = -n / 2 * (log(2 * pi) + log(fit[["deviance"]] / n) + log(weight_scale) + 1) + log_weights / 2
    reported_loglik = weighted_loglik - 1
  }
  if (family[["family"]] == "binomial" && isTRUE(snap_info[["fractional_w"]] > 0)) {
    working = if (length(fe)) fit[["working"]] else paste0(
      "(",
      dbglm_working_sql(base, x, fit[["beta"]], family[["link"]], cell = cell, conn = conn), "SELECT * FROM working)"
    )
    reported_loglik = dbglm_binomial_loglik(conn, working, weight_scale, cell = cell)
  }
  result = list(
    collin.coef = full,
    projection_diagnostics = resources[["projection_diagnostics"]],
    collin.var = setdiff(all_names, public), coef_names = public, coef_names_all = all_names,
    family = family, fml = saved_formula, model_spec = spec, has_intercept = has_intercept,
    yvar = parsed[["yvar"]], xvars = public, offset = offset, weights = weights,
    nobs_orig = n, nobs_input = as.numeric(count_input), nobs_dropped = as.numeric(count_input) - n,
    n_params_model = counts[["model"]], n_slopes = sum(public != "(Intercept)"),
    df_residual = n - counts[["model"]],
    loglik = reported_loglik, weighted_loglik = weighted_loglik, weight_scale = weight_scale,
    dispersion = dispersion, dispersion_scaled = dispersion_scaled,
    deviance = fit[["deviance"]] * weight_scale, converged = !isFALSE(fit[["converged"]]),
    finite_mle = !isFALSE(fit[["finite_mle"]]), separation = separation,
    iterations = fit[["iterations"]], score_max = fit[["score_max"]],
    deviance_history = fit[["deviance_history"]] * weight_scale, step_halvings = fit[["step_halvings"]],
    control = control, strategy = if (cell) "compress" else if (strategy == "ap") "ap" else "moments",
    n_cells = compression[["cells"]], compression_ratio = compression[["cells"]] / n, fe = if (length(fe)) display(parsed[["fe_labels"]]) else NULL,
    n_fe_levels = stats::setNames(counts[["sizes"]], display(parsed[["fe_labels"]])), fe_nested = counts[["nested"]],
    design_center = stats::setNames(unname(design_center[keep]), public),
    nobs_eligible = n_eligible, nobs_pruned = n_eligible - n,
    pruning_rounds = pruning[["rounds"]], exclusions = pruning[["counts"]], fixef_rm = fixef_rm,
    factor_sources = normalized[["factors"]], design_names = c(if (has_intercept) "(Intercept)", spec[["columns"]][["name"]])[keep],
    varying_slopes = normalized[["slopes"]], row_id = row_id, sample_table = sample_table, ap_cycles = fit[["ap_cycles"]], ap_residual = fit[["ap_residual"]],
    projection_iterations = fit[["ap_cycles"]], projection_residual = fit[["ap_residual"]],
    projection_method = if (!is.null(fit[["projection_method"]])) fit[["projection_method"]] else if (length(fe)) "group" else "none", call = call
  )
  class(result) = c("dbglm", "dbreg")
  unavailable = isFALSE(fit[["finite_mle"]])
  inference = new_inference_result(V, covariance[["vcov_type"]], public,
    if (covariance[["vcov_type"]] == "cluster") list(covariance[["cluster_var"]]) else list(),
    covariance_result[["clusters"]], ssc, counts[["covariance"]], test_df,
    if (family[["family"]] == "gaussian") "t" else "normal",
    status = if (unavailable) "unavailable" else "available",
    reason = if (unavailable) "nonfinite_mle" else NULL)
  result = install_inference_result(result, inference, beta)
  if (retain_prediction(retain) && length(fe)) {
    prediction_beta = numeric(length(all_x))
    prediction_beta[keep] = fit[["beta"]]
    result[["fitted_state"]] = dbglm_prediction_state(conn, fit[["working"]], original_base,
      fe, parsed[["fe_components"]], all_x, spec[["columns"]][["name"]],
      design_center, trend_aliases, counts[["sizes"]], create, drop, n,
      dictionary_components = prediction_keys[["columns"]], fitted_beta = prediction_beta,
      recovery_control = list(tol = control[["fe_tol"]], max_iter = control[["fe_max_iter"]])
    )
  }
  if (retain_sample(retain)) {
    if (cell) fit = dbglm_expand_cell_fit(
      conn, original_base, fit, x, fe,
      compression[["groups"]], family[["link"]], create
    )
    result = retain_glm_score_frame(
      result, fit, conn, original_base, x, fe, retained_source, resources,
      parsed[["fe_components"]], row_id, weights, weight_scale, counts, ssc
    )
  }
  if (!is.null(sample_table)) {
    keys = if (length(row_id)) row_id else "__dbglm_key"
    export_retained_sample(conn, sample_table, keys, from)
  }
  result = finalize_retained_resources(result, resources)
  result
}

# GLM-only formula normalization. Keep OLS parsing and contrasts unchanged.
dbglm_formula = function(fml) {
  fml = Formula::Formula(fml)
  if (length(fml)[1] != 1 || length(fml)[2] > 2) stop("Use one outcome and at most one fixed-effect block.", call. = FALSE)
  state = new.env(parent = emptyenv())
  state[["factors"]] = list()
  state[["factor_display"]] = character()
  state[["i_aliases"]] = character()
  state[["i_terms"]] = character()
  state[["slopes"]] = list()
  replace = function(e) {
    if (!is.call(e)) return(e)
    declaration = parse_sql_i_term(e, env = environment(fml))
    if (!is.null(declaration)) {
      index = match(declaration[["term"]], state[["i_terms"]])
      if (is.na(index)) {
        index = length(state[["i_terms"]]) + 1L
        state[["i_terms"]] = c(state[["i_terms"]], declaration[["term"]])
      }
      name = paste0("dbglm_factor_i", index, "__")
      state[["factors"]][[name]] = declaration[["first"]]
      state[["factor_display"]][name] = declaration[["first"]]
      state[["i_aliases"]] = unique(c(state[["i_aliases"]], name))
      result = declaration[["call"]]
      result[["factor_var"]] = as.name(name)
      return(result)
    }
    if ((identical(e[[1L]], as.name(":")) || identical(e[[1L]], as.name("*"))) && contains_i_term(e)) {
      stop("i() terms must be additive; interactions with i() are unsupported.", call. = FALSE)
    }
    if (identical(e[[1]], as.name("factor"))) {
      if (length(e) != 2 || !is.symbol(e[[2]])) stop("Use factor(column) without additional arguments.", call. = FALSE)
      source = as.character(e[[2]])
      name = paste0("dbglm_factor_", match(source, unique(c(unlist(state[["factors"]]), source))), "__")
      state[["factors"]][[name]] = source
      state[["factor_display"]][name] = paste0("factor(", source, ")")
      return(as.name(name))
    }
    as.call(lapply(as.list(e), replace))
  }
  rhs = stats::formula(fml, lhs = 0, rhs = 1)[[2]]
  rhs = replace(rhs)
  fe = if (length(fml)[2] > 1) stats::formula(fml, lhs = 0, rhs = 2)[[2]] else NULL
  if (contains_i_term(fe)) stop("i() is supported in additive regressor terms, not the fixed-effect block.", call. = FALSE)
  extract_slopes = function(e) {
    if (!is.call(e)) return(e)
    if (identical(e[[1]], as.name("["))) {
      args = as.list(e)[-1]
      if (length(args) < 2 || length(args) > 3 || !all(vapply(args, is.symbol, logical(1)))) stop("Use group[trend] or group[trend,trend2].", call. = FALSE)
      group = as.character(args[[1]])
      state[["slopes"]][[group]] = unique(c(state[["slopes"]][[group]], vapply(args[-1], as.character, character(1))))
      return(args[[1]])
    }
    as.call(lapply(as.list(e), extract_slopes))
  }
  if (!is.null(fe)) fe = replace(extract_slopes(fe))
  lhs = stats::formula(fml, lhs = 1, rhs = 0)[[2]]
  normalized = stats::as.formula(call("~", lhs, if (is.null(fe)) rhs else call("|", rhs, fe)), env = baseenv())
  list(
    fml = normalized, factors = state[["factors"]], slopes = state[["slopes"]],
    factor_display = state[["factor_display"]], i_aliases = state[["i_aliases"]]
  )
}

# A raw-key preflight avoids ranking and joining every FE dictionary just to
# discover that no group needs pruning. Only a proven empty removal set skips
# the existing recursive pruner; any removable group keeps that exact path.
dbglm_pruning_needed = function(conn, table, groups, y, mode, family) {
  if (!length(groups) || mode == "none") return(FALSE)
  if (length(unique(unlist(groups))) > 30L) return(TRUE)
  condition = switch(mode,
    singletons = "__n = 1",
    infinite_coef = if (family == "poisson") "__hi = 0" else "__lo = __hi",
    perfect_fit = if (family == "poisson") "(__n = 1 OR __hi = 0)" else "(__n = 1 OR __lo = __hi)"
  )
  sql = paste0(
    "SELECT COUNT(*) AS n FROM (SELECT COUNT(*) AS __n, MIN(", y,
    ") AS __lo, MAX(", y, ") AS __hi FROM ", table, " GROUP BY GROUPING SETS (",
    paste(vapply(groups, function(g) paste0("(", paste(g, collapse = ","), ")"), character(1)), collapse = ","),
    ")) __dbglm_prune WHERE ", condition
  )
  as.numeric(backend_get_query(conn, sql)[["n"]]) > 0
}

# Recursion operates on a snapshot and never rediscovers keys.
dbglm_prune = function(
  conn, table, fe, y, mode, create, drop, family = "binomial",
  flags = c("__dbglm_single", "__dbglm_infinite")
) {
  rounds = 0L
  counts = c(singleton = 0, infinite_coef = 0)
  if (!length(fe) || mode == "none") return(list(table = table, rounds = rounds, counts = counts))
  backend_info = detect_backend(conn)
  repeat {
    # Per-FE group summaries are tiny; joining them avoids N-row window sorts.
    singleton = infinite = joins = character()
    for (k in seq_along(fe)) {
      f = fe[k]
      s = paste0("__dbglm_p", k)
      joins = c(joins, paste0(
        " JOIN (SELECT ", f, " AS __dbglm_g, COUNT(*) AS n, MIN(", y, ") AS lo, MAX(", y,
        ") AS hi FROM ", table, " GROUP BY ", f, ") ", s, " ON __dbglm_t.", f, " = ", s, ".__dbglm_g"
      ))
      singleton = c(singleton, paste0(s, ".n = 1"))
      infinite = c(infinite, if (family == "poisson") paste0(s, ".hi = 0") else paste0(s, ".lo = ", s, ".hi"))
    }
    flagged = paste0(
      "(SELECT __dbglm_t.*, (", paste(singleton, collapse = " OR "), ") AS ", flags[[1L]], ", (",
      paste(infinite, collapse = " OR "), ") AS ", flags[[2L]], " FROM ", table, " __dbglm_t",
      paste(joins, collapse = ""), ")"
    )
    condition = switch(mode,
      perfect_fit = paste0("(", flags[[1L]], " OR ", flags[[2L]], ")"),
      singletons = flags[[1L]],
      infinite_coef = flags[[2L]]
    )
    info = backend_get_query(conn, paste0(
      "SELECT COUNT(*) FILTER (WHERE ", condition, ") AS removed, ",
      "COUNT(*) FILTER (WHERE ", condition, " AND ", flags[[1L]], ") AS single FROM ", flagged, " __dbglm_f"
    ))
    if (info[["removed"]] == 0) break
    counts = counts + c(singleton = as.numeric(info[["single"]]), infinite_coef = as.numeric(info[["removed"]] - info[["single"]]))
    # Each round refers to its input once per FE. Inlining successive rounds
    # on Databricks grows the plan exponentially, even on a tiny sample.
    # Materialize only after a removal; the common no-pruning case adds no write.
    next_table = create(paste0("SELECT ", sql_select_except(conn, flags, backend_info), " FROM ", flagged, " __dbglm_f WHERE NOT ", condition), materialize = TRUE)
    drop(table)
    table = next_table
    rounds = rounds + 1L
  }
  list(table = table, rounds = rounds, counts = counts)
}

# fixest's default parameter convention is intentionally not exact dummy rank.
dbglm_fe_counts = function(conn, table, fe, p, kind, ssc, fe_names = NULL, cluster_var = NULL, sizes = NULL) {
  if (!length(fe)) return(list(model = p, covariance = p, sizes = integer(), nested = logical()))
  if (is.null(sizes)) {
    size_row = backend_get_query(conn, paste(
      "SELECT",
      paste(sprintf("COUNT(DISTINCT %s) AS n%d", fe, seq_along(fe)), collapse = ", "),
      "FROM", table
    ))
    sizes = vapply(size_row, as.numeric, numeric(1))
  }
  slope_df = sum(sizes * lengths(attr(fe, "slopes")))
  model = p + slope_df + sum(sizes) - length(fe) + 1
  nested = rep(FALSE, length(fe))
  explicit = if (length(fe_names) && length(cluster_var)) vapply(strsplit(fe_names, "^", fixed = TRUE), function(parts) cluster_var %in% parts, logical(1)) else rep(FALSE, length(fe))
  if (kind == "cluster" && any(explicit)) nested = explicit
  if (kind == "cluster" && !any(explicit)) {
    nest_queries = vapply(seq_along(fe), function(i) {
      paste0(
        "SELECT ", i, " AS idx, MAX(n) AS max_n FROM (SELECT COUNT(DISTINCT __dbglm_cluster) AS n FROM ",
        table, " GROUP BY ", fe[i], ") t"
      )
    }, character(1))
    nest_df = backend_get_query(conn, paste(nest_queries, collapse = " UNION ALL "))
    nested[nest_df[["idx"]]] = (nest_df[["max_n"]] == 1)
  }
  K = model
  if (ssc[["K.fixef"]] == "none") K = p + slope_df
  if (kind == "cluster" && ssc[["K.fixef"]] == "nonnested" && any(nested)) {
    K = if (all(nested)) p + slope_df + 1 else model - sum(sizes[nested] - 1)
    K = max(K, p + 1)
  }
  list(model = model, covariance = K, sizes = sizes, nested = nested)
}

# Cell summaries use centered within-cell outcome moments. They preserve the
# original likelihood and individual-row sandwich meat at a common predictor.
dbglm_compress_base = function(conn, base, groups, family, center, create, drop,
                               weight_scale = 1) {
  groups = unique(groups)
  if (family == "binomial") return(dbglm_compress_binary(conn, base, groups, create, weight_scale))
  sum = function(expr) sql_stable_sum(conn, expr)
  centered = paste0("(__dbglm_y-(", dbglm_number(center), "))")
  mean = if (family == "gaussian") paste0(
    "(", dbglm_number(center), ")+",
    sum(paste0("__dbglm_a*", centered)), "/", sum("__dbglm_a")
  ) else
    paste0(sum("__dbglm_a*__dbglm_y"), "/", sum("__dbglm_a"))
  totals = paste0(
    "SELECT ", paste(c(
      groups, "MIN(__dbglm_key) AS __dbglm_key",
      paste(sum("__dbglm_a"), "AS __dbglm_a"),
      paste(mean, "AS __dbglm_y")
    ), collapse = ","),
    " FROM ", base, if (length(groups)) paste(" GROUP BY", paste(groups, collapse = ","))
  )
  means = create(totals)
  on.exit(drop(means), add = TRUE)
  delta = "(b.__dbglm_y-c.__dbglm_y)"
  within = switch(family,
    gaussian = paste0("b.__dbglm_a*", delta, "*", delta),
    poisson = paste0(
      "2.0*b.__dbglm_a*(CASE WHEN b.__dbglm_y=0 THEN c.__dbglm_y ELSE ",
      "b.__dbglm_y*LN(b.__dbglm_y/c.__dbglm_y)-b.__dbglm_y+c.__dbglm_y END)"
    )
  )
  fields = c(paste0("c.", groups), "c.__dbglm_key", "c.__dbglm_a", "c.__dbglm_y")
  aggregate = c(
    "COUNT(*) AS __dbglm_cell_n",
    paste(sum("b.__dbglm_a*b.__dbglm_a"), "AS __dbglm_a2"),
    paste(sum(paste0("b.__dbglm_a*b.__dbglm_a*", delta)), "AS __dbglm_a2delta"),
    paste(sum(paste0("b.__dbglm_a*b.__dbglm_a*", delta, "*", delta)), "AS __dbglm_a2delta2"),
    paste(sum(within), "AS __dbglm_dev_const")
  )
  join = if (length(groups)) paste0(
    " JOIN ", means, " c ON ",
    paste(paste0("b.", groups, "=c.", groups), collapse = " AND ")
  ) else paste0(" CROSS JOIN ", means, " c")
  create(paste0(
    "SELECT ", paste(c(fields, aggregate), collapse = ","), " FROM ", base, " b", join,
    " GROUP BY ", paste(fields, collapse = ",")
  ))
}

# A binary outcome needs only positive sums for its two outcome categories.
# Their squared-weight totals give exact HC1 moments without a second row pass.
dbglm_compress_binary = function(conn, base, groups, create, weight_scale) {
  sum = function(expr) sql_stable_sum(conn, expr)
  rounded = sql_round_even(paste0("__dbglm_a*(", dbglm_number(weight_scale), ")"))
  fields = c(
    groups, "MIN(__dbglm_key) AS __dbglm_key", "COUNT(*) AS __dbglm_cell_n",
    paste(sum("__dbglm_a"), "AS __dbglm_a"),
    paste0(sum("__dbglm_a*__dbglm_y"), "/", sum("__dbglm_a"), " AS __dbglm_y"),
    paste(sum("CASE WHEN __dbglm_y=1 THEN __dbglm_a*__dbglm_a ELSE 0.0 END"), "AS __dbglm_cell_s1"),
    paste(sum("CASE WHEN __dbglm_y=0 THEN __dbglm_a*__dbglm_a ELSE 0.0 END"), "AS __dbglm_cell_s0"),
    paste(sum(rounded), "AS __dbglm_round_a"),
    paste(sum(paste0("(", rounded, ")*__dbglm_y")), "AS __dbglm_round_ay")
  )
  grouped = paste0(
    "SELECT ", paste(fields, collapse = ","), " FROM ", base,
    if (length(groups)) paste(" GROUP BY", paste(groups, collapse = ","))
  )
  fields = c(
    groups, "__dbglm_key", "__dbglm_cell_n", "__dbglm_a", "__dbglm_y",
    "__dbglm_cell_s1+__dbglm_cell_s0 AS __dbglm_a2",
    "__dbglm_cell_s1*(1.0-__dbglm_y)-__dbglm_cell_s0*__dbglm_y AS __dbglm_a2delta",
    "__dbglm_cell_s1*POWER(1.0-__dbglm_y,2)+__dbglm_cell_s0*POWER(__dbglm_y,2) AS __dbglm_a2delta2",
    "0.0 AS __dbglm_dev_const", "__dbglm_round_a", "__dbglm_round_ay"
  )
  create(paste0("SELECT ", paste(fields, collapse = ","), " FROM (", grouped, ") binary_cells"))
}

sql_round_even = function(value) {
  paste0(
    "(CASE WHEN ", value, "-FLOOR(", value, ")<0.5 THEN FLOOR(", value, ") ",
    "WHEN ", value, "-FLOOR(", value, ")>0.5 THEN CEIL(", value, ") ",
    "WHEN MOD(FLOOR(", value, "),2.0)=0 THEN FLOOR(", value, ") ELSE CEIL(", value, ") END)"
  )
}

dbglm_cell_score_square = function(link) {
  variance = "GREATEST(0.0,__dbglm_a2delta2+2.0*(__dbglm_y-mu)*__dbglm_a2delta+POWER(__dbglm_y-mu,2)*__dbglm_a2)"
  if (link %in% c("identity", "log")) return(variance)
  paste0("(", variance, ")*POWER(d/(mu*(1.0-mu)),2)")
}

dbglm_cell_source = function(conn, base, groups, n, strategy, family, center,
                             create, drop, weight_scale, ratio, nmax) {
  if (!strategy %in% c("auto", "compress")) {
    return(list(table = base, enabled = FALSE, cells = n, groups = groups))
  }
  if (strategy == "compress") {
    table = dbglm_compress_base(conn, base, groups, family, center, create, drop, weight_scale)
    size = as.numeric(backend_get_query(conn, paste("SELECT COUNT(*) AS n FROM", table))[["n"]])
    return(list(table = table, enabled = TRUE, cells = size, groups = groups))
  }
  count = if (length(groups)) paste0(
    "SELECT COUNT(*) AS n FROM (SELECT DISTINCT ",
    paste(groups, collapse = ","), " FROM ", base, ") cell_keys"
  ) else "SELECT 1 AS n"
  size = as.numeric(backend_get_query(conn, count)[["n"]])
  if (strategy == "auto" && (size / n > ratio || size > nmax)) {
    return(list(table = base, enabled = FALSE, cells = n, groups = groups))
  }
  table = dbglm_compress_base(conn, base, groups, family, center, create, drop, weight_scale)
  list(table = table, enabled = TRUE, cells = size, groups = groups)
}

# Only requested original-row inference needs this relation. Join the compact
# final projector to original rows; never refit or repeat weighted projection.
dbglm_expand_cell_fit = function(conn, original, fit, x, fe, groups, link, create) {
  if (!length(fe)) return(fit)
  rows = fit[["raw_projected"]]
  if (is.null(rows)) {
    fields = c("b.*", if (length(x)) paste0("p.", x, " AS ", x, "_tilde"))
    rows = paste0(
      "(SELECT ", paste(fields, collapse = ","), " FROM ", fit[["working"]], " b",
      if (length(x)) paste(" JOIN", fit[["projected"]], "p USING (__dbglm_key)"), ")"
    )
  }
  fields = c(
    groups, "eta AS dbglm_cell_eta__", "mu AS dbglm_cell_mu__",
    "d AS dbglm_cell_derivative__", if (length(x)) paste0(x, "_tilde")
  )
  lookup = create(paste("SELECT", paste(fields, collapse = ","), "FROM", rows))
  multiplier = if (link %in% c("identity", "log")) "1.0" else
    "c.dbglm_cell_derivative__/(c.dbglm_cell_mu__*(1.0-c.dbglm_cell_mu__))"
  score = paste0("b.__dbglm_a*(b.__dbglm_y-c.dbglm_cell_mu__)*(", multiplier, ") AS score")
  fields = c(
    "b.*", "c.dbglm_cell_eta__ AS eta", "c.dbglm_cell_mu__ AS mu", score,
    if (length(x)) paste0("c.", x, "_tilde AS ", x, "_tilde")
  )
  fit[["raw_projected"]] = paste0(
    "(SELECT ", paste(fields, collapse = ","), " FROM ", original,
    " b JOIN ", lookup, " c ON ", paste(paste0("b.", groups, "=c.", groups), collapse = " AND "), ")"
  )
  fit
}
