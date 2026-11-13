#' Control database GLM convergence
#'
#' @param glm_tol Relative deviance tolerance.
#' @param score_tol Tolerance for the score scaled by expected information.
#' @param glm_max_iter Maximum number of scoring iterations.
#' @param step_max_halvings Maximum number of halvings during initial predictor
#'   recovery and during each likelihood iteration.
#' @param collin_tol Tolerance for identifying a scaled design basis.
#' @param projection_method Use bounded summaries with SQL fallback, or
#'   keep fixed-effect projection state in SQL from the start. Numerically weak
#'   multi-effect systems may transfer bounded exact working cells for QR recovery. Native Databricks
#'   single additive effects use a certified SQL window projection automatically.
#'   `"duckdb_bounded"` opts into materializing working rows and attempting
#'   bounded projection before weak-information recovery for DuckDB cloglog
#'   fits with multiple additive effects. Projection certificates and fallback
#'   remain active. Other fits handle this setting like `"auto"`.
#' @param irls Scoring rule for absorbed cloglog fits. Newton uses observed
#'   information for the step. Covariance uses expected information.
#' @param tolerate_max_iter Return the last accepted state with convergence
#'   set to false when the iteration limit or line search is exhausted.
#' @param compute_inference Compute covariance and retain score coordinates.
#'   Set to `FALSE` for intermediate fits that only need fitted predictors.
#' @param fe_tol Absolute and relative group-equation tolerance for projection.
#'   Defaults to \code{getOption("dbreg.ap_tol", 1e-10)}.
#' @param fe_max_iter Maximum total inner SQL CG iterations and bounded-moment
#'   refinement budget per projection. Residual refinements also have an eight-step cap.
#'   Defaults to \code{getOption("dbreg.cg_max_iter", 10000L)}.
#' @param tighten_fe_tol Tighten projection tolerance to \code{score_tol / 100}.
#'   Set to \code{FALSE} to preserve a relaxed projection tolerance for intermediate
#'   fits. Reporting convergence always requires the original strict score gates.
#' @return A validated list of controls.
#' @export
dbglm_control = function(glm_tol = 1e-8, score_tol = 1e-8, glm_max_iter = 100L,
                         step_max_halvings = 25L, collin_tol = 1e-10,
                         tolerate_max_iter = FALSE, irls = c("newton", "fisher"),
                         projection_method = c("auto", "sql", "duckdb_bounded"), compute_inference = TRUE,
                         fe_tol = getOption("dbreg.ap_tol", 1e-10),
                         fe_max_iter = getOption("dbreg.cg_max_iter", 10000L),
                         tighten_fe_tol = TRUE) {
  controls = list(glm_tol = glm_tol, score_tol = score_tol,
    glm_max_iter = glm_max_iter, step_max_halvings = step_max_halvings,
    collin_tol = collin_tol, fe_tol = fe_tol, fe_max_iter = fe_max_iter)
  for (name in names(controls)) {
    value = controls[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0) {
      stop("`", name, "` must be one finite positive number.", call. = FALSE)
    }
  }
  for (name in c("glm_max_iter", "step_max_halvings", "fe_max_iter")) {
    value = controls[[name]]
    if (value != floor(value) || value > .Machine[["integer.max"]]) {
      stop("`", name, "` must be a positive integer.", call. = FALSE)
    }
    controls[[name]] = as.integer(value)
  }
  if (!is.logical(tolerate_max_iter) || length(tolerate_max_iter) != 1L ||
      is.na(tolerate_max_iter)) {
    stop("`tolerate_max_iter` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(compute_inference) || length(compute_inference) != 1L || is.na(compute_inference)) {
    stop("`compute_inference` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(tighten_fe_tol) || length(tighten_fe_tol) != 1L || is.na(tighten_fe_tol)) {
    stop("`tighten_fe_tol` must be TRUE or FALSE.", call. = FALSE)
  }
  controls[["tighten_fe_tol"]] = tighten_fe_tol
  controls[["compute_inference"]] = compute_inference
  controls[["tolerate_max_iter"]] = tolerate_max_iter
  controls[["irls"]] = match.arg(irls)
  controls[["projection_method"]] = match.arg(projection_method)
  controls
}

#' Fit a generalized linear model on a database
#'
#' @param fml Formula with a single outcome column, structural regressors and
#'   optional additive fixed effects after a pipe. Join components with `^`
#'   for a native tuple key. Use `group[trend]` for a group-specific slope.
#'   Varying slopes use SQL projection and can include multiple numeric trends.
#' @param conn,table,data,path Input connection and source, as in [dbreg()].
#' @param family Gaussian identity, Poisson log, or binomial logit, probit or
#'   complementary log-log. Supply a built-in family, function, or family name.
#' @param weights,offset Optional names or one-sided formulas naming columns.
#' @param factor_refs Named reference levels for categorical regressors.
#' @param vcov `"iid"`, `"hc1"` or a one-sided clustering formula.
#' @param ssc Small-sample correction policy, as in [dbreg()].
#' @param start Finite starting values, one per identified coefficient. Numerically
#'   unusable starting predictors are damped toward the default predictor, up to
#'   \code{control$step_max_halvings} times. Known offsets remain unchanged.
#' @param etastart Optional numeric column name or one-sided formula giving
#'   initial linear predictors, including offsets. The initial predictor is
#'   projected onto the fitted design using prior weights and takes precedence
#'   over `start`. Missing or nonfinite starting values omit their rows.
#' @param control Convergence controls from [dbglm_control()].
#' @param verbose Print iteration diagnostics.
#' @param separation Use `"prediction"` for binomial boundary probabilities
#'   without inference when a finite estimate is unavailable.
#' @param subset Optional row filter in the same form as [dbreg()].
#' @param null_model Fit a structural constant or offset-only model on the
#'   same rows and store its likelihood and pseudo-R-squared diagnostics.
#' @param fixef_rm Recursive fixed-effect exclusion policy. Remove singleton
#'   groups, groups with infinite coefficients, both, or neither.
#' @param strategy Use bounded design-cell compression, original-row moments,
#'   or select between them based on estimated cell count.
#' @param compress_ratio Maximum fraction of original rows in automatic mode.
#' @param compress_nmax Maximum allowed number of compressed cells.
#' @param retain Choose \code{"prediction"} for learned fixed effects,
#'   \code{"inference"} for fitted rows and post-fit score inference, or
#'   \code{"all"} for both. Prediction alone keeps dictionaries and releases
#'   temporary rows before returning. Inference supports IID, HC1 and one-way
#'   clustered covariance. Use \code{close()} to release owned database tables.
#' @param row_id Optional vector of source columns that uniquely identify
#'   retained rows.
#' @return A fitted `dbglm` object.
#' @details The selected sample is snapshotted on DuckDB or native Databricks.
#'   The solver transfers bounded information, score and deviance aggregates.
#'   Weights are nonnegative and zero-weight rows are excluded. Binomial
#'   outcomes must be zero or one, and Poisson outcomes must be nonnegative
#'   (non-integer values such as trade flows are allowed, as in PPML).
#'
#'   IID covariance uses expected information. Gaussian dispersion is the
#'   weighted residual sum of squares divided by full-model residual degrees
#'   of freedom; Poisson and binomial fix dispersion at one. Independently,
#'   \code{ssc = list(K.adj = TRUE)} multiplies IID covariance by
#'   \code{(n - 1) / (n - K)}, including absorbed fits. Here \code{n} counts
#'   original fitted rows, and \code{K} follows the selected FE counting policy.
#'   Set \code{K.adj = FALSE} to disable this additional correction. The same
#'   policy applies to retained post-fit inference.
#'
#'   Initial recovery also applies to projected row-level warm starts. Its
#'   halvings are included in \code{step_halvings}. Recovery never relaxes the
#'   final score or finite-solution checks; exhausted recovery raises an error.
#'   Initial numerical projection failures may also trigger recovery. Driver,
#'   capacity and projection iteration-limit errors do not trigger initial recovery.
#'
#'   Convergence requires independent score certificates and finite identified
#'   information. A damped step may converge only when both the previous and
#'   current structural and fixed-effect scores pass their strict tolerance,
#'   with the same deviance and boundary checks as a full step. A small
#'   deviance change alone does not establish convergence. Absorbed fits also
#'   check the full Fisher predictor correction, including when no structural
#'   coefficients remain. Weak multi-effect information uses bounded working
#'   design cells and Householder QR recovery; cell summaries and residuals must
#'   fit the transfer and memory budgets. Unsupported recovery raises an error.
#'
#'   Fixed-effect groups with small information must also pass their
#'   Newton-scale group equation; changing a varying axis's units does not weaken
#'   this guard.
#' @param sample_table Optional new temporary table receiving fitted source keys.
#' @param cluster Optional cluster variable for clustered standard errors, as a
#'   character column name or one-sided formula.
#' @export
dbglm = function(fml, conn = NULL, table = NULL, data = NULL, path = NULL,
                 family = stats::gaussian(), weights = NULL, offset = NULL,
                 factor_refs = NULL, vcov = "iid", ssc = "full", start = NULL,
                 control = dbglm_control(), verbose = getOption("dbreg.verbose", FALSE),
                 separation = c("error", "prediction"), subset = NULL,
                 null_model = FALSE,
                 fixef_rm = c("perfect_fit", "singletons", "infinite_coef", "none"),
                 strategy = c("auto", "compress", "moments"),
                 compress_ratio = NULL, compress_nmax = 1e6,
                 retain = c("none", "prediction", "inference", "all"), row_id = NULL, etastart = NULL, sample_table = NULL,
                 cluster = NULL) {
  retain = match.arg(retain)
  source = normalize_data_source(table, data, path)
  table = source[["table"]]
  data = source[["data"]]
  path = source[["path"]]
  validate_sample_args(row_id, NULL, conn, table, FALSE, FALSE)
  validate_glm_sample_export(sample_table, conn, table, row_id)
  keep_rows = retain != "none"
  keep_source = retain %in% c("inference", "all")
  strategy = match.arg(strategy)
  if (!is.numeric(compress_nmax) || length(compress_nmax) != 1L ||
      is.na(compress_nmax) || compress_nmax < 1 || compress_nmax != floor(compress_nmax)) {
    stop("`compress_nmax` must be one positive whole number or Inf.", call. = FALSE)
  }
  fixef_rm = match.arg(fixef_rm)
  if (!is.logical(null_model) || length(null_model) != 1L || is.na(null_model)) {
    stop("`null_model` must be TRUE or FALSE.", call. = FALSE)
  }
  family = dbglm_family(family)
  separation = match.arg(separation)
  if (separation == "prediction" && family[["family"]] != "binomial") {
    stop("Boundary prediction requires a binomial family.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.list(control) || is.null(names(control)) || anyDuplicated(names(control)) ||
      any(!names(control) %in% names(formals(dbglm_control)))) {
    stop("`control` must be a named list of dbglm_control arguments.", call. = FALSE)
  }
  control = do.call(dbglm_control, control)
  weights = formula_column(weights, "weights")
  offset = formula_column(offset, "offset")
  etastart = formula_column(etastart, "etastart")
  for (value in list(weights, offset, etastart)) {
    if (!is.null(value) && (!is.character(value) || length(value) != 1L ||
        is.na(value) || !nzchar(value))) {
      stop("Weights, offset and etastart must each name one column.", call. = FALSE)
    }
  }
  if (is.null(vcov) || (is.character(vcov) && length(vcov) != 1L)) {
    stop("`vcov` must name one covariance type or cluster column.", call. = FALSE)
  }
  if (inherits(vcov, "formula")) formula_column(vcov, "vcov")
  covariance = parse_vcov_args(vcov, cluster)
  ssc = normalize_ssc(ssc)
  spec = dbglm_formula_spec(fml, factor_refs)
  fml = spec[["formula"]]
  if (length(spec[["groups"]]) && separation == "prediction") {
    stop("Boundary prediction requires a structural-only model.", call. = FALSE)
  }
  resources = fit_resources()
  on.exit(release_fit_resources(resources), add = TRUE)
  r_positions = keep_rows && is.null(table) && is.data.frame(data) && is.null(row_id)
  source_id = row_id
  if (r_positions) {
    requested = c(model_vars(fml), weights, offset, etastart, covariance[["cluster_var"]], subset_variables(subset))
    source_id = safe_identifiers("dbreg_input_position", taken = c(names(data), requested))
    data = as.data.frame(data)
    data[[source_id]] = seq_len(nrow(data))
  }
  setup = setup_db_connection(conn, table, data, path, caller = "dbglm",
    resources = resources)
  conn = setup[["conn"]]
  inherited = attr(conn, "dbreg_resource_lease", exact = TRUE)
  if (!is.null(inherited)) {
    retained_fit_connection(inherited)
    resources[["parent_session"]] = inherited[["session"]]
  }
  backend = detect_backend(conn)
  if (backend[["name"]] != "duckdb" && !isTRUE(backend[["is_databricks"]])) {
    stop("Database GLMs require DuckDB or native Databricks.", call. = FALSE)
  }
  validate_native_temp_query(conn, backend)
  if (!is.null(sample_table)) {
    if (isTRUE(resources[["own_conn"]]) || isTRUE(inherited[["session"]][["own_conn"]])) {
      stop("Sample exports require a caller-owned connection.", call. = FALSE)
    }
    if (sample_table_exists(conn, sample_table)) stop("Sample export table already exists.", call. = FALSE)
  }
  sample = dbglm_sample(conn, resources, setup[["from_statement"]], fml,
    family, weights, offset, etastart, covariance[["cluster_var"]], subset, spec[["factor_refs"]],
    spec[["intercept"]], spec[["groups"]], fixef_rm, source_id, keep_rows || !is.null(sample_table), keep_source)
  if (r_positions) {
    public = setdiff(names(sample[["cluster_fields"]]), source_id)
    sample[["cluster_fields"]] = sample[["cluster_fields"]][public]
    sample[["cluster_classes"]] = sample[["cluster_classes"]][public]
  }
  dbglm_run_sample(conn, resources, sample, spec, family, weights, offset, etastart,
    covariance, ssc, start, control, separation, strategy, compress_ratio, compress_nmax,
    retain, row_id, r_positions, null_model, verbose, sample_table)
}

dbglm_formula_spec = function(fml, factor_refs) {
  if (is.character(fml)) fml = stats::as.formula(fml, env = baseenv())
  fml = Formula::Formula(fml)
  shape = as.integer(length(fml))
  if (shape[1L] != 1L || !shape[2L] %in% 1:2) {
    stop("Expected one outcome and one structural formula part.", call. = FALSE)
  }
  structural = stats::formula(fml, lhs = 1, rhs = 1)
  parsed = parse_regression_formula(structural, allow_empty = TRUE)
  groups = if (shape[2L] == 2L) parse_fe_components(stats::formula(fml, lhs = 0, rhs = 2)[[2L]]) else list()
  intercept = parsed[["has_intercept"]]
  factor_refs = validate_factor_refs(factor_refs, parsed[["xvars"]])
  list(formula = fml, groups = groups, intercept = intercept,
    yvar = parsed[["yvar"]], factor_refs = factor_refs)
}

dbglm_run_sample = function(conn, resources, sample, spec, family, weights, offset, etastart,
                             covariance, ssc, start, control, separation, strategy, compress_ratio,
                             compress_nmax, retain, row_id, r_positions = FALSE,
                             null_model = FALSE, verbose = FALSE, sample_table = NULL) {
  keep_source = retain %in% c("inference", "all")
  keep_scores = keep_source && control[["compute_inference"]]
  sample[["original_base"]] = sample[["base"]]
  manager = fe_cg_manager(conn, resources)
  cluster = covariance[["vcov_type"]] == "cluster"
  cell = dbglm_cell_source(conn, sample[["base"]],
    c(sample[["columns"]][["alias"]], sample[["fe"]], unname(sample[["trend_columns"]]), "__dbglm_o",
      if (cluster) "__dbglm_cluster", if (length(etastart)) "__dbglm_etastart"),
    sample[["n"]], strategy, family[["family"]], sample[["moments"]][["ybar"]],
    manager[["create"]], manager[["drop"]], sample[["weight_scale"]],
    resolve_compress_ratio(compress_ratio, sample[["fe"]]), compress_nmax, sample[["sizes"]])
  sample[["base"]] = cell[["table"]]
  sample[["cell"]] = cell[["enabled"]]
  if (keep_scores) {
    sample[["score_keys"]] = if (cell[["enabled"]]) cell[["groups"]] else sample[["row_keys"]]
  }
  sample[["fit_rows"]] = cell[["cells"]]
  if (length(sample[["trend_columns"]])) {
    basis = dbglm_varying_basis(conn, sample, manager)
    sample[["base"]] = basis[["table"]]
    sample[["ops"]] = basis[["ops"]]
  }
  estimated = dbglm_estimate(conn, sample, family, control, start, covariance,
    ssc, separation, verbose, resources)
  result = estimated[["result"]]
  saved = stats::formula(spec[["formula"]])
  environment(saved) = baseenv()
  result[["fml"]] = Formula::Formula(saved)
  result[["offset"]] = offset
  result[["etastart"]] = etastart
  result[["weights"]] = weights
  result[["yvar"]] = spec[["yvar"]]
  result[["retain"]] = retain
  result[["row_id"]] = row_id
  result[["sample_table"]] = sample_table
  if (null_model) {
    result[["null_diagnostics"]] = dbglm_null_fit(conn, sample, result, spec[["intercept"]] || length(spec[["groups"]]) > 0L)
  }
  if (retain %in% c("prediction", "all") && length(sample[["fe"]])) {
    fit = estimated[["fit"]]
    if (sample[["y_center"]] != 0 && !is.null(fit[["working"]])) {
      fit[["working"]] = paste0("(SELECT ", sql_select_except(conn, "eta"), ", eta+(",
        dbglm_number(sample[["y_center"]]), ") AS eta FROM ", fit[["working"]], ") shifted_working")
    }
    state = learned_effect_attempt(dbglm_learned_state(conn, sample, result, fit),
      if (!length(sample[["trend_columns"]])) {
        beta = unname(result[["collin.coef"]])
        beta[is.na(beta)] = 0
        database_additive_effects(conn, fit[["working"]], sample[["key_source"]], sample[["fe"]],
          sample[["key_groups"]], sample[["public_groups"]],
          sample[["columns"]][["alias"]], beta, sample[["centers"]], sample[["fit_rows"]], resources)
      })
    result[["fitted_state"]] = state[["result"]]
    result[["prediction_diagnostics"]] = state[setdiff(names(state), "result")]
  }
  if (keep_source) {
    retained = dbglm_retain_rows(conn, resources, sample, cell[["groups"]],
      estimated[["fit"]], estimated[["x"]], family, row_id, keep_source, keep_scores,
      result[["fitted_state"]][["tables"]], transform = estimated[["transform"]])
    result[["retained_inference"]] = retained[["inference"]]
    retained[["rows"]][["r_input"]] = r_positions
    result[["retained_rows"]] = retained[["rows"]]
    result[["resource_state"]] = retained[["lease"]]
  } else if (retain == "prediction" && length(result[["fitted_state"]][["tables"]])) {
    result[["resource_state"]] = retain_fit_resources(resources,
      result[["fitted_state"]][["tables"]], session = resources[["parent_session"]])
  }
  if (!is.null(sample_table)) {
    tryCatch(export_glm_sample(conn, sample_table, sample, row_id), error = function(error) {
      release_fit_resources(result[["resource_state"]])
      stop(error)
    })
  }
  result
}

dbglm_sample = function(conn, resources, from, fml, family, weights, offset, etastart,
                         cluster, subset, factor_refs, intercept, groups, fixef_rm, row_id, keep_rows,
                         keep_source) {
  variables = unique(c(model_vars(fml), weights, offset, etastart, cluster))
  fields = c("__dbglm_y", "__dbglm_a", "__dbglm_o", "__dbglm_intercept",
    "__dbglm_cluster", "__dbglm_etastart")
  fe = if (length(groups)) paste0("__dbglm_fe", seq_along(groups)) else character()
  slopes = attr(groups, "slopes", exact = TRUE)
  trends = unique(unlist(slopes, use.names = FALSE))
  trend_columns = if (length(trends)) setNames(paste0("__dbglm_t", seq_along(trends)), trends) else character()
  needed = unique(c(variables, subset_variables(subset), row_id))
  original_columns = if (keep_source) names(DBI::dbGetQuery(conn, paste("SELECT *", from, "LIMIT 0")))
  src = safe_source(conn, from, c(needed, original_columns),
    reserved = c(fields, fe, unname(trend_columns), "__dbglm_row"))
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  safe = src[["to_safe"]](variables)
  snapshot = fe_projection_tables(conn, resources)
  source_columns = src[["to_safe"]](needed)
  keys = src[["to_safe"]](row_id)
  row_key = if (keep_rows && !length(keys)) {
    keys = "__dbglm_row"
    if (isTRUE(detect_backend(conn)[["is_databricks"]])) "monotonically_increasing_id() AS __dbglm_row"
    else "ROW_NUMBER() OVER () AS __dbglm_row"
  }
  snapshot_columns = src[["to_safe"]](unique(c(needed, original_columns)))
  source_fields = c(q(snapshot_columns), row_key)
  raw = snapshot(paste("SELECT", paste(source_fields, collapse = ","),
    src[["from_statement"]]), unique(c(snapshot_columns, keys)))
  source_columns = unique(c(source_columns, keys))
  from = paste("FROM", raw)
  if (length(row_id)) validate_source_row_id(conn, from, keys)
  schema = DBI::dbGetQuery(conn, paste("SELECT *", from, "LIMIT 0"))
  cluster_fields = setNames(src[["to_safe"]](original_columns), original_columns)
  cluster_classes = setNames(vapply(schema[cluster_fields], function(x) class(x)[1L], character(1)),
    original_columns)
  safe_fml = src[["to_safe"]](fml)
  y = all.vars(stats::formula(Formula::Formula(safe_fml), lhs = 1, rhs = 0))[[1L]]
  w = src[["to_safe"]](weights)
  o = src[["to_safe"]](offset)
  eta = src[["to_safe"]](etastart)
  cl = src[["to_safe"]](cluster)
  if (!(is.numeric(schema[[y]]) || is.logical(schema[[y]])) ||
      any(!vapply(schema[c(w, o, eta)], is.numeric, logical(1)))) {
    stop("Outcome, weights, offset and etastart must be numeric. Logical outcomes are allowed.",
      call. = FALSE)
  }
  if (length(trends) && !all(vapply(schema[src[["to_safe"]](trends)], is.numeric, logical(1)))) {
    stop("Varying-slope GLMs require numeric trend columns.", call. = FALSE)
  }
  validate_weights(w, NULL, conn, from)
  input_n = as.numeric(DBI::dbGetQuery(conn, paste("SELECT COUNT(*) AS n", from))[["n"]])
  from = apply_source_subset(conn, from, subset, src[["to_safe"]])
  numeric = safe[vapply(schema[safe], function(x) is.numeric(x) && !inherits(x, "integer64"), logical(1))]
  eligible = c(paste(q(safe), "IS NOT NULL"), sql_is_finite(conn, q(numeric)),
    if (length(w)) paste(q(w), "> 0"))
  from = filter_source_predicate(from, paste(eligible, collapse = " AND "))
  y_sql = sql_cast_double(q(y), conn)
  w_sql = if (length(w)) sql_cast_double(q(w), conn) else "1.0"
  invalid = switch(family[["family"]], gaussian = "FALSE",
    poisson = paste0(y_sql, " < 0"),
    binomial = paste(y_sql, "NOT IN (0, 1)"))
  sample_stats = function(from) DBI::dbGetQuery(conn, paste0(
    "SELECT COUNT(*) AS n, MAX(", w_sql, ") AS weight_max, ",
    "MIN(", y_sql, ") AS ymin, MAX(", y_sql, ") AS ymax, ",
    "SUM(CASE WHEN ", invalid, " THEN 1 ELSE 0 END) AS bad_y ", from))
  stats = sample_stats(from)
  n = as.numeric(stats[["n"]])
  eligible_n = n
  if (!n) stop("No observations in the GLM estimation sample.", call. = FALSE)
  if (stats[["bad_y"]] > 0) {
    stop(if (family[["family"]] == "poisson") "Poisson outcomes must be nonnegative."
      else "Binomial outcomes must be zero or one.", call. = FALSE)
  }
  key_source = NULL
  key_groups = list()
  pruning = list(rounds = 0L, counts = c(singleton = 0, infinite_coef = 0))
  if (length(fe)) {
    key_groups = lapply(groups, src[["to_safe"]])
    ranked = dbglm_ranked_sample(conn, resources, from, source_columns,
      key_groups, fe, y_sql, fixef_rm, family[["family"]])
    from = ranked[["from"]]
    key_source = sub("^FROM ", "", from)
    pruning = ranked[["pruning"]]
    stats = sample_stats(from)
    n = as.numeric(stats[["n"]])
    if (!n) stop("All observations were removed by fixed-effect singleton/separation pruning.", call. = FALSE)
  }
  if ((intercept || length(fe)) && (family[["family"]] == "binomial" && stats[["ymin"]] == stats[["ymax"]] ||
      family[["family"]] == "poisson" && stats[["ymax"]] == 0)) {
    stop("Constant outcome with an intercept has no finite GLM solution.", call. = FALSE)
  }
  weight_max = as.numeric(stats[["weight_max"]])
  weight_scale = if (weight_max < 1 || weight_max > 2) 2^floor(log2(weight_max)) else 1
  names(factor_refs) = src[["to_safe"]](names(factor_refs))
  design = sql_model_matrix(Formula::Formula(safe_fml), conn, paste0("(SELECT * ", from, ") glm_selected"),
    expand = "all", intercept = intercept || length(fe) > 0L, factor_refs = factor_refs)
  columns = design[["columns"]]
  aliases = if (nrow(columns)) paste0("__dbglm_x", seq_len(nrow(columns))) else character()
  columns[["alias"]] = aliases
  columns[["recipe"]] = lapply(columns[["recipe"]], restore_recipe_sources, map = src[["map"]])
  columns[["label"]] = vapply(columns[["recipe"]], linear_recipe_label, character(1))
  intercept = intercept && !length(fe)
  trend_fields = if (length(trends)) paste(sql_cast_double(q(src[["to_safe"]](trends)), conn),
    "AS", unname(trend_columns))
  key_aliases = if (length(keys)) paste0("__dbglm_key", seq_along(keys)) else character()
  select = c(fe, trend_fields, if (length(keys)) paste(q(keys), "AS", key_aliases), paste(y_sql, "AS __dbglm_y"),
    paste0("(", w_sql, ")/(", dbglm_number(weight_scale), ") AS __dbglm_a"),
    paste(if (length(o)) sql_cast_double(q(o), conn) else "0.0", "AS __dbglm_o"),
    if (length(cl)) paste(q(cl), "AS __dbglm_cluster"),
    if (length(eta)) paste(sql_cast_double(q(eta), conn), "AS __dbglm_etastart"),
    if (intercept) "CAST(1.0 AS DOUBLE) AS __dbglm_intercept",
    if (length(aliases)) paste(sql_cast_double(design[["select_exprs"]], conn), "AS", aliases))
  if (intercept) {
    constant = data.frame(id = "intercept", term = 0L, basis = 1L,
      alias = "__dbglm_intercept", label = "(Intercept)", role = "slope")
    constant[["recipe"]] = list(list(kind = "constant"))
    columns = rbind(constant, columns)
  }
  x = columns[["alias"]]
  if (!length(x) && !length(fe)) stop("No GLM design columns were specified.", call. = FALSE)
  limits = projection_limits()
  terms = length(x) * (as.double(length(x)) + 1) / 2
  if (terms > min(limits[["design_entries"]], limits[["transfer_cells"]] - length(x) - 8) ||
      64 * length(x)^2 > limits[["memory_bytes"]]) {
    stop("GLM information exceeds the design or memory budget.", call. = FALSE)
  }
  if (anyDuplicated(columns[["label"]])) {
    stop("Formula terms produce duplicate public coefficient names.", call. = FALSE)
  }
  base = snapshot(paste("SELECT", paste(select, collapse = ", "), from),
    c(x, fe, fields[1:3], if (length(cl)) "__dbglm_cluster"))
  is_gaussian = family[["family"]] == "gaussian"
  has_offset = length(o) > 0L
  partition_cands = list()
  if (!intercept && !length(fe) && nrow(columns)) {
    for (t in unique(columns[["term"]])) {
      cand = which(columns[["term"]] == t)
      kinds = vapply(columns[["recipe"]][cand], function(r) r[["kind"]], character(1))
      if (length(cand) >= 1L && all(kinds %in% c("indicator", "i", "constant"))) {
        partition_cands[[length(partition_cands) + 1L]] = cand
      }
    }
  }
  part_sums = unlist(lapply(seq_along(partition_cands), function(i) {
    expr = paste(x[partition_cands[[i]]], collapse = "+")
    c(sprintf("MIN(%s) AS __part_min_%d", expr, i),
      sprintf("MAX(%s) AS __part_max_%d", expr, i))
  }), use.names = FALSE)
  finite = paste(sql_is_finite(conn, x), collapse = " AND ")
  x_non_intercept = x[x != "__dbglm_intercept"]
  y_adj_expr = if (is_gaussian && has_offset) "(__dbglm_y-__dbglm_o)" else "__dbglm_y"
  sums = c(
    if (length(x)) paste0("SUM(CASE WHEN ", finite, " THEN 0 ELSE 1 END) AS bad") else "0 AS bad",
    "SUM(LN(__dbglm_a)) AS log_weights",
    if (is_gaussian && has_offset) "0.0 AS offset_mean" else "SUM(__dbglm_a*__dbglm_o)/SUM(__dbglm_a) AS offset_mean",
    paste0("SUM(__dbglm_a*", y_adj_expr, ")/SUM(__dbglm_a) AS ybar"),
    part_sums,
    if (length(x_non_intercept)) c(
      sprintf("MIN(%s) AS __min_%d", x_non_intercept, seq_along(x_non_intercept)),
      sprintf("MAX(%s) AS __max_%d", x_non_intercept, seq_along(x_non_intercept))
    ))
  moments = DBI::dbGetQuery(conn, paste("SELECT", paste(sums, collapse = ", "), "FROM", base))
  if (moments[["bad"]] > 0) stop("GLM design contains nonfinite transformed values.", call. = FALSE)
  partition_idx = integer()
  for (i in seq_along(partition_cands)) {
    if (isTRUE(moments[[paste0("__part_min_", i)]] == 1) &&
        isTRUE(moments[[paste0("__part_max_", i)]] == 1)) {
      partition_idx = partition_cands[[i]]
      break
    }
  }
  centers = numeric(length(x))
  scales = rep(1, length(x))
  step_scale = rep(1, length(x))
  y_center = 0
  if (is_gaussian && (intercept || length(fe) > 0L || length(partition_idx) > 0L)) {
    yc = as.numeric(moments[["ybar"]])
    if (is.finite(yc) && yc != 0) {
      y_center = yc
      moments[["ybar"]] = 0
    }
  }
  if (length(x_non_intercept)) {
    mins = vapply(seq_along(x_non_intercept), function(i) as.numeric(moments[[paste0("__min_", i)]]), numeric(1))
    maxs = vapply(seq_along(x_non_intercept), function(i) as.numeric(moments[[paste0("__max_", i)]]), numeric(1))
    idx = match(x_non_intercept, x)
    if (length(fe)) {
      centers[idx] = mins / 2 + maxs / 2
    } else {
      if (intercept) {
        centers[idx] = mins / 2 + maxs / 2
      } else if (length(partition_idx)) {
        center_pos = which(!(idx %in% partition_idx))
        if (length(center_pos)) {
          centers[idx[center_pos]] = mins[center_pos] / 2 + maxs[center_pos] / 2
        }
      }
      spans = (maxs - mins) / 2
      scales[idx] = ifelse(is.finite(spans) & spans > 0, spans, 1)
    }
    # Measure predictor changes in fixed design units, independently of
    # family weights or how close a supplied start is to a probability boundary.
    span = pmax(abs(mins - centers[idx]), abs(maxs - centers[idx])) / scales[idx]
    step_scale[idx] = ifelse(is.finite(span) & span > 0, span, 1)
  }
  transforms = character()
  drop_cols = character()
  if (is_gaussian && (has_offset || y_center != 0)) {
    if (has_offset) {
      transforms = c(transforms, "__dbglm_o AS __dbglm_saved_o", "CAST(0.0 AS DOUBLE) AS __dbglm_o")
      drop_cols = c(drop_cols, "__dbglm_o")
    }
    transforms = c(transforms, sprintf("(%s-(%s)) AS __dbglm_y", y_adj_expr, dbglm_number(y_center)))
    drop_cols = c(drop_cols, "__dbglm_y")
    if (length(eta)) {
      eta_adj_expr = if (has_offset) "(__dbglm_etastart-__dbglm_o)" else "__dbglm_etastart"
      transforms = c(transforms, sprintf("(%s-(%s)) AS __dbglm_etastart", eta_adj_expr, dbglm_number(y_center)))
      drop_cols = c(drop_cols, "__dbglm_etastart")
    }
  }
  changed_x = which(centers != 0 | scales != 1)
  if (length(changed_x)) {
    x_exprs = ifelse(scales[changed_x] != 1,
      sprintf("((%s-(%s))/(%s)) AS %s", x[changed_x], dbglm_number(centers[changed_x]),
        dbglm_number(scales[changed_x]), x[changed_x]),
      sprintf("(%s-(%s)) AS %s", x[changed_x], dbglm_number(centers[changed_x]), x[changed_x]))
    transforms = c(transforms, x_exprs)
    drop_cols = c(drop_cols, x[changed_x])
  }
  if (length(transforms)) {
    old_base = base
    base = snapshot(paste("SELECT", sql_select_except(conn, drop_cols), ",",
      paste(transforms, collapse = ","), "FROM", old_base), c(x, fe))
    fe_cg_manager(conn, resources)[["drop"]](old_base)
  }
  if (!length(fe) && length(x)) {
    h_row = DBI::dbGetQuery(conn, paste("SELECT", paste(dbglm_moment_terms(x, "__dbglm_a"), collapse = ", "), "FROM", base))
    for (nm in names(h_row)) moments[[nm]] = h_row[[nm]]
  }
  sizes = if (length(fe)) {
    row = DBI::dbGetQuery(conn, paste("SELECT",
      paste(sprintf("COUNT(DISTINCT %s) AS n%d", fe, seq_along(fe)), collapse = ","), "FROM", base))
    unname(vapply(row, as.numeric, numeric(1)))
  } else numeric()
  list(base = base, columns = columns, moments = moments, n = n, input_n = input_n,
    eligible_n = eligible_n, pruning = pruning, fe = fe, sizes = sizes,
    fe_labels = fe_component_labels(groups), weight_scale = weight_scale,
    slopes = slopes, trend_columns = trend_columns, row_keys = key_aliases,
    original_keys = keys, raw = raw, centers = centers, scales = scales,
    y_center = y_center, step_scale = step_scale, has_offset = is_gaussian && has_offset,
    offset_present = has_offset, partition_idx = partition_idx,
    key_source = key_source, key_groups = key_groups, public_groups = groups,
    cluster_fields = cluster_fields, cluster_classes = cluster_classes,
    default_cluster = cluster, warm = length(eta) > 0L)
}

dbglm_estimate = function(conn, sample, family, control, start, covariance,
                            ssc, separation, verbose, resources) {
  columns = sample[["columns"]]
  moments = sample[["moments"]]
  fe = sample[["fe"]]
  base = sample[["base"]]
  n = sample[["n"]]
  weight_scale = sample[["weight_scale"]]
  cell = sample[["cell"]]
  fit_rows = sample[["fit_rows"]]
  cluster = covariance[["vcov_type"]] == "cluster"
  alpha = dbglm_initial_intercept(conn, base, family, moments,
    offset_present = isTRUE(sample[["offset_present"]]))
  if (length(fe)) {
    sizes = sample[["sizes"]]
    ops = sample[["ops"]]
    all_x = columns[["alias"]]
    lookup_manager = fe_cg_manager(conn, resources)
    fit_absorbed = function(manager = NULL) {
      reused = NULL
      if (sample[["warm"]]) {
        warm = dbglm_warm_absorbed(conn, base, all_x, fe, sizes, fit_rows, control, manager, ops,
          lookup_manager = lookup_manager)
        effects = warm[["effects"]]
        keep = warm[["keep"]]
        beta = warm[["beta"]]
      } else {
        effects = if (is.null(manager)) {
          limits = projection_limits()
          width = length(all_x) + 1L
          if (sum(sizes) * width > limits[["design_entries"]] ||
              !fe_dense_feasible(sizes, width, limits)) {
            projection_abort("system_limit", "GLM fixed-effect information exceeds its projection budget.")
          }
          values = lapply(sizes, numeric)
          values[[1L]][] = alpha
          values
        } else {
          dbglm_fe_sql_start(base, fe, alpha, manager, ops)
        }
        # Carry the main fit's row columns so a full basis can reuse this
        # initial projection rather than solving it again.
        initial = dbglm_fe_irls(conn, base, all_x, fe, sizes, fit_rows,
          numeric(length(all_x)), effects, family, control, identify_only = TRUE,
          cluster = cluster && control[["compute_inference"]], cell = cell, manager = manager, ops = ops,
          score_keys = sample[["score_keys"]], initial_intercept = alpha, lookup_manager = lookup_manager)
        effects = initial[["effects"]]
        reused = initial[["initial"]]
        keep = initial[["restart"]]
        beta = dbglm_start(start, columns[["label"]][keep], alpha)
      }
      x = all_x[keep]
      labels = columns[["label"]][keep]
      fit = dbglm_fe_irls(conn, base, x, fe, sizes, fit_rows, beta, effects,
        family, control, verbose, cluster = cluster && control[["compute_inference"]], cell = cell, manager = manager, ops = ops,
        score_keys = sample[["score_keys"]], initial_intercept = alpha, lookup_manager = lookup_manager,
        step_scale = sample[["step_scale"]][keep], initial = reused)
      if (!is.null(fit[["restart"]])) {
        stop("GLM starting weights do not identify the selected structural basis.", call. = FALSE)
      }
      list(fit = fit, keep = keep, x = x, labels = labels)
    }
    # One plain fixed effect is an exact window projection; bounded moments
    # would inline every level into each query.
    window = length(fe) == 1L && is.null(ops)
    attempt = if (control[["projection_method"]] != "sql" && is.null(ops) && !window) {
      projection_attempt("bounded_moments", fit_absorbed())
    } else NULL
    result = attempt[["result"]]
    if (is.null(result)) {
      result = fit_absorbed(fe_cg_manager(conn, resources))
    }
    fit = result[["fit"]]
    fit[["projection_diagnostics"]] = attempt[setdiff(names(attempt), "result")]
    keep = result[["keep"]]
    x = result[["x"]]
    labels = result[["labels"]]
  } else {
    information = dbglm_read_information(moments, nrow(columns))
    # Centering uses the dummy partition as the constant coordinate. Keep that
    # partition before selecting other columns, so the retained raw design spans
    # the same space even when an earlier numeric term duplicates a dummy.
    order = unique(c(sample[["partition_idx"]], seq_len(nrow(columns))))
    keep = sort(order[dbglm_basis(information[order, order, drop = FALSE], control[["collin_tol"]])])
    x = columns[["alias"]][keep]
    labels = columns[["label"]][keep]
    # A complete dummy partition spans the same constant as an intercept.
    # Its fitted columns have been scaled, so seed that direction in the
    # solver's coordinates, including when recovering an unusable start.
    constant = numeric(length(keep))
    intercept_pos = match("(Intercept)", labels)
    if (!is.na(intercept_pos)) {
      constant[intercept_pos] = 1
    } else if (length(sample[["partition_idx"]]) && all(sample[["partition_idx"]] %in% keep)) {
      part_pos = which(keep %in% sample[["partition_idx"]])
      constant[part_pos] = sample[["scales"]][keep][part_pos]
    }
    if (sample[["warm"]]) {
      beta = dbglm_warm_structural(conn, base, x)
    } else {
      beta = dbglm_start(start, labels, alpha, constant)
      if (!is.null(start)) {
        scales_k = sample[["scales"]][keep]
        centers_k = sample[["centers"]][keep]
        intercept_pos = match("(Intercept)", labels)
        if (!is.na(intercept_pos)) {
          if (sample[["y_center"]] != 0) beta[intercept_pos] = beta[intercept_pos] - sample[["y_center"]]
          slope_pos = setdiff(seq_along(labels), intercept_pos)
          if (length(slope_pos)) beta[intercept_pos] = beta[intercept_pos] + sum(centers_k[slope_pos] * beta[slope_pos])
        } else if (length(sample[["partition_idx"]]) && all(sample[["partition_idx"]] %in% keep)) {
          part_pos = which(keep %in% sample[["partition_idx"]])
          if (sample[["y_center"]] != 0) beta[part_pos] = beta[part_pos] - sample[["y_center"]]
          slope_pos = setdiff(seq_along(labels), part_pos)
          if (length(slope_pos)) beta[part_pos] = beta[part_pos] + sum(centers_k[slope_pos] * beta[slope_pos])
        }
        beta = beta * scales_k
      }
    }
    fit = dbglm_irls(conn, base, x, beta, family, control, verbose, separation, cell = cell,
      anchor = dbglm_start(NULL, labels, alpha, constant), step_scale = sample[["step_scale"]][keep])
  }
  p = length(keep)
  requested = control[["compute_inference"]]
  counts = dbglm_parameter_counts(conn, sample, p, cluster && requested, ssc)
  df = n - counts[["model"]]
  gaussian = family[["family"]] == "gaussian"
  if (gaussian && df <= 0) stop("Insufficient residual degrees of freedom.", call. = FALSE)
  finite_mle = !isFALSE(fit[["finite_mle"]])
  kind = covariance[["vcov_type"]]
  available = finite_mle && requested
  cov = if (available) dbglm_covariance(conn, base, x, fit, family, kind, ssc, n, K = counts[["covariance"]], cell = cell)
    else list(matrix = matrix(NA_real_, p, p), clusters = NULL)
  dispersion_scaled = if (gaussian) fit[["deviance"]] / df else 1
  V = cov[["matrix"]]
  if (kind == "iid") V = V * (if (gaussian) dispersion_scaled else 1 / weight_scale)
  transform = NULL
  prediction_coordinates = NULL
  beta_out = fit[["beta"]]
  if (!length(fe) && p) {
    scales_k = sample[["scales"]][keep]
    centers_k = sample[["centers"]][keep]
    B = diag(1 / scales_k, nrow = p, ncol = p)
    intercept_pos = match("(Intercept)", labels)
    if (!is.na(intercept_pos)) {
      for (m in setdiff(seq_len(p), intercept_pos)) {
        B[intercept_pos, m] = -centers_k[m] / scales_k[m]
      }
    } else if (length(sample[["partition_idx"]]) && all(sample[["partition_idx"]] %in% keep)) {
      part_pos = which(keep %in% sample[["partition_idx"]])
      for (j in part_pos) {
        for (m in setdiff(seq_len(p), part_pos)) {
          B[j, m] = -centers_k[m] / scales_k[m]
        }
      }
    }
    transform = B
    prediction_coordinates = list(centers = centers_k, scales = scales_k,
      beta = fit[["beta"]], y_center = sample[["y_center"]], o_center = 0)
    beta_out = as.numeric(B %*% beta_out)
    if (sample[["y_center"]] != 0) {
      if (!is.na(intercept_pos)) {
        beta_out[intercept_pos] = beta_out[intercept_pos] + sample[["y_center"]]
      } else if (length(sample[["partition_idx"]])) {
        part_pos = which(keep %in% sample[["partition_idx"]])
        beta_out[part_pos] = beta_out[part_pos] + sample[["y_center"]]
      }
    }
    if (available) {
      V = B %*% V %*% t(B)
      V = (V + t(V)) / 2
    }
  }
  if (available && (any(!is.finite(V)) || any(diag(V) < 0))) {
    stop("GLM covariance is invalid; check identification and weight scale.", call. = FALSE)
  }
  dimnames(V) = list(labels, labels)
  clusters = cov[["clusters"]]
  test_df = if (is.numeric(ssc[["t.df"]])) ssc[["t.df"]] else
    if (length(clusters) && ssc[["t.df"]] == "min") clusters - 1 else n - counts[["covariance"]]
  attr(V, "type") = kind
  attr(V, "df.t") = test_df
  attr(V, "n_params") = attr(V, "df.K") = counts[["covariance"]]
  attr(V, "ssc") = ssc
  attr(V, "n_clusters") = clusters
  coeftable = inference_table(stats::setNames(beta_out, labels), V, test_df, normal = !gaussian)
  likelihood = dbglm_likelihood(conn, sample, x, fit, family)
  compact = compact_recipe_dictionaries(columns[["recipe"]])
  columns[["recipe"]] = compact[["recipes"]]
  full = stats::setNames(rep(NA_real_, nrow(columns)), columns[["label"]])
  full[keep] = beta_out
  inference = list(status = if (available) "available" else "unavailable",
    reason = if (!finite_mle) "nonfinite_mle" else if (!requested) "not_requested" else NULL,
    distribution = if (gaussian) "t" else "normal",
    conditional_on_iterate = isFALSE(fit[["converged"]]))
  result = structure(list(coeftable = coeftable, vcov = V, inference_result = inference, coef_names = labels,
    xvars = setdiff(labels, "(Intercept)"), collin.coef = full,
    collin.var = columns[["label"]][setdiff(seq_len(nrow(columns)), keep)], family = family,
    linear_design = list(version = 1L, columns = columns,
      fitted_ids = columns[["id"]][keep], dictionaries = compact[["dictionaries"]]),
    prediction_coordinates = prediction_coordinates,
    nobs_orig = n, nobs_input = sample[["input_n"]],
    nobs_dropped = sample[["input_n"]] - n, n_params_model = counts[["model"]], df_residual = df,
    fe = if (length(fe)) sample[["fe_labels"]],
    n_fe_levels = setNames(sample[["sizes"]], sample[["fe_labels"]]),
    nobs_eligible = sample[["eligible_n"]], nobs_pruned = sample[["eligible_n"]] - n,
    pruning_rounds = sample[["pruning"]][["rounds"]], exclusions = sample[["pruning"]][["counts"]],
    fe_score_max = fit[["fe_score_max"]], irls = fit[["irls"]],
    projection_method = fit[["projection_method"]],
    projection_tolerance = fit[["projection_tolerance"]],
    projection_max_iter = fit[["projection_max_iter"]], varying_slopes = sample[["slopes"]],
    projection_diagnostics = fit[["projection_diagnostics"]],
    dispersion = if (gaussian) dispersion_scaled * weight_scale else 1,
    loglik = likelihood[["reported"]], weighted_loglik = likelihood[["weighted"]],
    deviance = fit[["deviance"]] * weight_scale, finite_mle = finite_mle,
    converged = !isFALSE(fit[["converged"]]), iterations = fit[["iterations"]],
    score_max = fit[["score_max"]], stop_reason = if (!finite_mle) "boundary" else
      if (is.null(fit[["stop_reason"]])) "converged" else fit[["stop_reason"]],
    deviance_history = fit[["deviance_history"]] * weight_scale,
    step_halvings = fit[["step_halvings"]], strategy = if (cell) "compress" else "moments", n_cells = fit_rows,
    compression_ratio = fit_rows / n, control = control),
    class = c("dbglm", "dbreg"))
  if (!result[["converged"]] && finite_mle) {
    warning("dbglm did not converge after ", result[["iterations"]], " iterations (",
      result[["stop_reason"]], "); returning the last accepted state.", call. = FALSE)
  }
  list(result = result, fit = fit, x = x, transform = transform)
}

#' @rdname predict.dbreg
#' @param type Linear predictor (default) or response mean.
#' @export
predict.dbglm = function(object, newdata = NULL, type = c("link", "response"), ...,
                         conn = NULL, table = NULL, path = NULL, row_id = NULL, prediction_table = NULL,
                         sample = c("estimation", "original")) {
  check_unused_args(list(...))
  type = match.arg(type)
  sample = match.arg(sample)
  if (missing(newdata)) newdata = NULL
  sql = sql_prediction_requested(conn, table, path, row_id, prediction_table)
  if (sample == "original" && (!is.null(newdata) || sql)) {
    stop("Original-sample prediction uses the retained fitting source.", call. = FALSE)
  }
  if (sql) return(predict_keyed_sql(object, newdata, conn, table, path, row_id, prediction_table, type))
  if (sample == "original") require_inference_rows(object)
  if (is.null(newdata) && !is.null(object[["retained_rows"]])) {
    field = if (type == "link") "__dbglm_saved_eta" else "__dbglm_saved_mu"
    return(retained_row_values(object, field, if (type == "link") ".link" else ".fitted",
      na.rm = sample == "estimation"))
  }
  if (is.null(newdata)) stop("newdata is required for GLM predictions.", call. = FALSE)
  newdata = as.data.frame(newdata)
  if (length(object[["fe"]])) {
    if (is.null(object[["fitted_state"]])) {
      stop("Absorbed GLM prediction requires retained fixed effects.", call. = FALSE)
    }
    eta = predict_learned_effects(object, newdata)
  } else {
    matrix = linear_prediction_matrix(object, newdata)
    invalid = attr(matrix, "dbreg_invalid_rows")
    sc = object[["prediction_coordinates"]]
    if (!is.null(sc)) {
      if (ncol(matrix)) {
        for (j in seq_len(ncol(matrix))) {
          if (sc[["centers"]][j] != 0 || sc[["scales"]][j] != 1) {
            matrix[, j] = (matrix[, j] - sc[["centers"]][j]) / sc[["scales"]][j]
          }
        }
      }
      eta = as.numeric(matrix %*% sc[["beta"]]) + sc[["y_center"]] + linear_prediction_offset(object, newdata)
    } else {
      eta = as.numeric(matrix %*% coef(object)) + linear_prediction_offset(object, newdata)
    }
    eta[invalid | !is.finite(eta)] = NA_real_
  }
  if (type == "link") eta else if (identical(object[["family"]][["link"]], "log")) exp(eta) else
    object[["family"]][["linkinv"]](eta)
}

#' @rdname logLik.dbreg
#' @param type Use the fixest reporting convention or the unrounded weighted objective.
#' @export
logLik.dbglm = function(object, type = c("fixest", "weighted"), ...) {
  check_unused_args(list(...))
  type = match.arg(type)
  value = object[[if (type == "weighted") "weighted_loglik" else "loglik"]]
  extra_df = if (type == "weighted" && identical(object[["family"]][["family"]], "gaussian")) 1L else 0L
  structure(value, class = "logLik",
    df = if (object[["finite_mle"]]) object[["n_params_model"]] + extra_df else NA_integer_,
    nobs = object[["nobs_orig"]], finite_mle = object[["finite_mle"]])
}

dbglm_start = function(start, labels, alpha, constant = NULL) {
  beta = numeric(length(labels))
  if (is.null(start)) {
    if (is.null(constant)) constant = as.numeric(labels == "(Intercept)")
    positions = which(constant != 0)
    beta[positions] = alpha * constant[positions]
  } else {
    if (!is.numeric(start) || length(start) != length(beta) || any(!is.finite(start))) {
      stop("`start` must contain one finite value per retained coefficient.", call. = FALSE)
    }
    if (!is.null(names(start))) {
      if (anyDuplicated(names(start)) || !setequal(names(start), labels)) {
        stop("Named starts must match retained coefficients.", call. = FALSE)
      }
      start = start[labels]
    }
    beta[] = start
  }
  beta
}
