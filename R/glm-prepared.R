#' Prepare an immutable DuckDB GLM design
#'
#' A prepared context owns a snapshot of the original source and its encoded
#' design. Subsequent refits supply keyed weights and optional linear predictors.
#' Close the context when finished.
#' @param fml,family,offset,factor_refs,fixef_rm,control As in \code{\link{dbglm}}.
#' @param conn A live DuckDB connection supplied by the caller.
#' @param table Name of the source table.
#' @param row_id One unique, nonmissing source key column.
#' @param weights,etastart Distinct names of the mutable columns in a refit.
#' @param vcov The prepared covariance type or one-way cluster formula.
#' @return A prepared context for \code{\link{dbglm_refit}}.
#' @export
dbglm_prepare = function(fml, conn, table, row_id, family = stats::gaussian(),
                         weights = "w", etastart = "eta", vcov = "iid", offset = NULL,
                         factor_refs = NULL, fixef_rm = "perfect_fit", control = dbglm_control()) {
  if (!resource_connection_valid(conn) || detect_backend(conn)[["name"]] != "duckdb") {
    stop("Prepared GLM contexts require a live DuckDB connection.", call. = FALSE)
  }
  inherited = attr(conn, "dbreg_resource_lease", exact = TRUE)
  if (!is.null(inherited)) retained_fit_connection(inherited)
  validate_sample_args(row_id, NULL, conn, table, FALSE, FALSE)
  if (length(row_id) != 1L) stop("Prepared GLMs require one row_id column.", call. = FALSE)
  weights = formula_column(weights, "weights")
  etastart = formula_column(etastart, "etastart")
  offset = formula_column(offset, "offset")
  columns = list(weights, etastart, table, if (length(offset)) offset)
  if (any(!vapply(columns[!vapply(columns, is.null, logical(1))], function(x)
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x), logical(1))) ||
    !length(weights) || !length(etastart) || !length(table) ||
    identical(tolower(weights), tolower(etastart))) {
    stop("Provide a table and distinct weight and predictor column names.", call. = FALSE)
  }
  if (inherits(vcov, "formula")) formula_column(vcov, "vcov")
  covariance = parse_vcov_args(vcov)
  cluster = covariance[["cluster_var"]]
  family = dbglm_family(family)
  spec = dbglm_formula_spec(fml, factor_refs)
  fixef_rm = match.arg(fixef_rm, c("perfect_fit", "singletons", "infinite_coef", "none"))
  control = do.call(dbglm_control, control)
  immutable = unique(c(model_vars(spec[["formula"]]), row_id, offset, cluster))
  mutable = c(weights, etastart)
  if (any(tolower(mutable) %in% tolower(immutable))) {
    stop("Mutable columns cannot be immutable model variables.", call. = FALSE)
  }
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  resources = fit_resources()
  resources[["conn"]] = conn
  on.exit(release_fit_resources(resources))
  create = fe_projection_tables(conn, resources)
  schema = DBI::dbGetQuery(conn, paste("SELECT * FROM", q(table), "LIMIT 0"))
  kept = names(schema)[!tolower(names(schema)) %in% tolower(mutable)]
  if (!all(immutable %in% names(schema))) stop("Prepared source columns are missing.", call. = FALSE)
  source = create(paste("SELECT", paste(c(q(kept), paste("1.0::DOUBLE AS", q(weights)),
    paste("0.0::DOUBLE AS", q(etastart))), collapse = ","), "FROM", q(table)), c(kept, mutable))
  if (length(cluster)) {
    column = q(cluster)
    bad = paste0(column, " IS NULL", if (is.numeric(schema[[cluster]])) paste0(" OR NOT ", sql_is_finite(conn, column)))
    invalid = DBI::dbGetQuery(conn, paste("SELECT COUNT(*) AS n FROM", source, "WHERE", bad))[["n"]]
    if (invalid > 0) stop("Prepared cluster keys must be nonmissing and finite.", call. = FALSE)
  }
  sample = dbglm_sample(conn, resources, paste("FROM", source), spec[["formula"]], family,
    weights, offset, etastart, cluster, NULL, spec[["factor_refs"]], spec[["intercept"]],
    spec[["groups"]], fixef_rm, row_id, TRUE, TRUE)
  owned = resources[["tables"]][match(c(sample[["raw"]], sample[["base"]]), q(resources[["tables"]]))]
  formula = stats::formula(spec[["formula"]])
  environment(formula) = baseenv()
  spec[["formula"]] = Formula::Formula(formula)
  lease = retain_fit_resources(resources, owned, session = inherited[["session"]])
  structure(list(sample = sample, spec = spec, family = family, row_id = row_id,
    weights = weights, etastart = etastart, offset = offset, covariance = covariance,
    fixef_rm = fixef_rm, control = control, resource_state = lease), class = "dbglm_prepared")
}

#' Refit a prepared DuckDB GLM
#'
#' Positive finite weights and finite warm predictors reuse the encoded
#' absorbed design. Changes to eligibility and structural-only models use the
#' normal fitting path on the immutable original source.
#' @param object A context from \code{\link{dbglm_prepare}}.
#' @param state_table Name of a table with exactly one row per original source
#'   key, including observations omitted at fitting.
#' @param etastart Use the supplied warm predictor column.
#' @param vcov Optional covariance type or the prepared cluster formula.
#' @param control Controls as in \code{\link{dbglm}}.
#' @param compact Keep keys, weights and fitted diagnostics. Requires deferred
#'   inference and an unchanged prepared sample.
#' @return A `dbglm` fit with flags recording prepared reuse and compaction.
#' @export
dbglm_refit = function(object, state_table, etastart = TRUE, vcov = NULL,
                        control = object[["control"]], compact = FALSE) {
  if (!inherits(object, "dbglm_prepared")) stop("Expected a dbglm_prepared context.", call. = FALSE)
  parent = object[["resource_state"]]
  conn = retained_fit_connection(parent)
  if (!is.character(state_table) || length(state_table) != 1L || is.na(state_table) || !nzchar(state_table)) {
    stop("state_table must name one table.", call. = FALSE)
  }
  for (value in list(etastart, compact)) {
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop("etastart and compact must be TRUE or FALSE.", call. = FALSE)
    }
  }
  control = do.call(dbglm_control, control)
  if (compact && control[["compute_inference"]]) stop("Compact refits require compute_inference=FALSE.", call. = FALSE)
  original_covariance = object[["covariance"]]
  if (inherits(vcov, "formula")) formula_column(vcov, "vcov")
  covariance = if (is.null(vcov)) original_covariance else parse_vcov_args(vcov)
  cluster = covariance[["cluster_var"]]
  if (length(cluster) && !identical(cluster, original_covariance[["cluster_var"]])) {
    stop("Cluster must match the prepared context.", call. = FALSE)
  }
  resources = fit_resources()
  resources[["conn"]] = conn
  resources[["parent_session"]] = parent[["session"]]
  on.exit(release_fit_resources(resources))
  create = fe_projection_tables(conn, resources)
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  key = object[["row_id"]]
  weight = object[["weights"]]
  eta = object[["etastart"]]
  columns = c(key, weight, if (etastart) eta)
  schema = DBI::dbGetQuery(conn, paste("SELECT * FROM", q(state_table), "LIMIT 0"))
  if (!all(columns %in% names(schema))) stop("State table is missing required columns.", call. = FALSE)
  if (!all(vapply(schema[setdiff(columns, key)], is.numeric, logical(1)))) {
    stop("State weights and predictors must be numeric.", call. = FALSE)
  }
  state = create(paste("SELECT", paste(q(columns), collapse = ","), "FROM", q(state_table)), columns)
  sample = object[["sample"]]
  old = sample[["raw"]]
  raw_key = sample[["original_keys"]]
  raw_fields = sample[["cluster_fields"]]
  join = paste0("s.", q(key), "=b.", q(raw_key))
  w = sql_cast_double(paste0("s.", q(weight)), conn)
  start = if (etastart) sql_cast_double(paste0("s.", q(eta)), conn) else "0.0"
  validate_source_row_id(conn, paste("FROM", state), key)
  bad = DBI::dbGetQuery(conn, paste("SELECT COUNT(*) AS n, COUNT(", paste0("b.", q(raw_key)),
    ") AS matched, MIN(", w, ") AS min_weight,",
    "SUM(CASE WHEN", w, "IS NULL OR", w, "<=0 OR NOT", sql_is_finite(conn, w),
    "OR", start, "IS NULL OR NOT", sql_is_finite(conn, start), "THEN 1 ELSE 0 END) AS changed",
    "FROM", state, "s LEFT JOIN", old, "b ON", join))
  if (!isTRUE(as.numeric(bad[["n"]]) == sample[["input_n"]]) ||
      !isTRUE(as.numeric(bad[["matched"]]) == sample[["input_n"]])) {
    stop("State keys must match the prepared source exactly, once each.", call. = FALSE)
  }
  # Verify the reverse cardinality too. A coercion between different SQL key
  # types must never match two source keys to one state key.
  reverse = DBI::dbGetQuery(conn, paste("SELECT COUNT(*) AS n FROM", old, "b WHERE NOT EXISTS",
    "(SELECT 1 FROM", state, "s WHERE", join, ")"))[["n"]]
  if (reverse != 0) stop("State keys must match the prepared source exactly, once each.", call. = FALSE)
  if (isTRUE(bad[["min_weight"]] < 0)) stop("State weights must be nonnegative.", call. = FALSE)
  original_fields = setdiff(names(raw_fields), c(weight, eta))
  eligible = as.numeric(bad[["changed"]]) == 0
  reuse = eligible && length(sample[["fe"]]) > 0L
  selected = if (compact && reuse) key else original_fields
  target = if (reuse) unname(raw_fields) else names(raw_fields)
  names(target) = names(raw_fields)
  kept = c(selected, weight, if (!(compact && reuse)) eta)
  expressions = c(paste0("b.", q(raw_fields[selected])), w, if (!(compact && reuse)) start)
  raw = create(paste("SELECT", paste(paste(expressions, "AS", q(target[kept])), collapse = ","),
    "FROM", old, "b JOIN", state, "s ON", join), unname(target[kept]))
  spec = object[["spec"]]
  if (!reuse) {
    attr(conn, "dbreg_resource_lease") = parent
    result = dbglm(spec[["formula"]], conn = conn, table = raw,
      family = object[["family"]], weights = weight, offset = object[["offset"]],
      etastart = if (etastart) eta, factor_refs = spec[["factor_refs"]],
      vcov = if (length(cluster)) stats::as.formula(call("~", as.name(cluster)), env = baseenv()) else covariance[["vcov_type"]],
      fixef_rm = object[["fixef_rm"]], control = control, row_id = key, retain = "inference", strategy = "moments")
  } else {
    sample = dbglm_prepared_sample(conn, create, sample, state, raw, key, w, start, etastart, kept, cluster)
    result = dbglm_run_sample(conn, resources, sample, spec, object[["family"]],
      weight, object[["offset"]], if (etastart) eta, covariance, normalize_ssc("full"),
      NULL, control, "error", "moments", NULL, 1e6, "inference", key)
  }
  result[["prepared_reused"]] = reuse
  result[["prepared_compact"]] = compact && reuse
  result
}

dbglm_prepared_sample = function(conn, create, sample, state, raw, key, weight, eta, warm, fields, cluster) {
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  base = sample[["base"]]
  from = paste("FROM", base, "b JOIN", state, "s ON", paste0("b.", sample[["row_keys"]], "=s.", q(key)))
  maximum = as.numeric(DBI::dbGetQuery(conn, paste("SELECT MAX(", weight, ") AS maximum", from))[["maximum"]])
  scale = if (maximum < 1 || maximum > 2) 2^floor(log2(maximum)) else 1
  selected = DBI::dbGetQuery(conn, paste("SELECT * FROM", base, "LIMIT 0"))
  keep = setdiff(names(selected), c("__dbglm_a", "__dbglm_etastart"))
  eta_expr = eta
  if (warm) {
    if (isTRUE(sample[["has_offset"]])) eta_expr = paste0("((", eta_expr, ")-b.__dbglm_saved_o)")
    if (length(sample[["y_center"]]) && sample[["y_center"]] != 0) {
      eta_expr = paste0("((", eta_expr, ")-(", dbglm_number(sample[["y_center"]]), "))")
    }
  }
  expressions = c(paste0("b.", q(keep)),
    paste0(weight, "/(", dbglm_number(scale), ") AS __dbglm_a"),
    paste(eta_expr, "AS __dbglm_etastart"))
  sample[["base"]] = create(paste("SELECT", paste(expressions, collapse = ","), from),
    c(keep, "__dbglm_a", "__dbglm_etastart"))
  moments = DBI::dbGetQuery(conn, paste("SELECT SUM(LN(__dbglm_a)) AS log_weights,",
    "SUM(__dbglm_a*__dbglm_y)/SUM(__dbglm_a) AS ybar,",
    "SUM(__dbglm_a*__dbglm_o)/SUM(__dbglm_a) AS offset_mean FROM", sample[["base"]]))
  sample[["moments"]][names(moments)] = moments
  sample[["weight_scale"]] = scale
  sample[["warm"]] = warm
  sample[["raw"]] = raw
  sample[["cluster_fields"]] = sample[["cluster_fields"]][fields]
  sample[["cluster_classes"]] = sample[["cluster_classes"]][fields]
  sample[["default_cluster"]] = cluster
  sample
}

#' @rdname dbglm_prepare
#' @param con A prepared context.
#' @param ... Unused.
#' @export
close.dbglm_prepared = function(con, ...) close_retained_lease(con[["resource_state"]])
