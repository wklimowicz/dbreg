prepare_linear_retention = function(inputs, from, variables) {
  conn = inputs[["conn"]]
  backend = detect_backend(conn)
  if (backend[["name"]] != "duckdb" && !isTRUE(backend[["is_databricks"]])) {
    stop("Retained linear fits require DuckDB or native Databricks.", call. = FALSE)
  }
  validate_native_temp_query(conn, backend)
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  keys = inputs[["row_id"]]
  fields = unique(c(variables, keys))
  if (linear_keeps_scores(inputs)) {
    schema = DBI::dbGetQuery(conn, sql_limit(paste("SELECT *", from), 0, backend[["name"]]))
    fields = names(schema)
    public = setdiff(names(schema), inputs[["generated_row_id"]])
    inputs[["retained_source_classes"]] = vapply(schema[public], function(x) class(x)[1L], character(1))
  }
  expressions = q(fields)
  if (!length(keys)) {
    keys = safe_identifiers("dbreg_input_position", taken = fields)
    row_key = if (isTRUE(backend[["is_databricks"]])) "monotonically_increasing_id()" else
      "ROW_NUMBER() OVER ()"
    expressions = c(expressions, paste(row_key, "AS", q(keys)))
    fields = c(fields, keys)
  }
  create = fe_projection_tables(conn, inputs[["resources"]])
  original = create(paste("SELECT", paste(expressions, collapse = ","), from), fields)
  inputs[["row_id"]] = keys
  inputs[["retained_original"]] = original
  paste("FROM", original)
}

retain_linear_rows = function(inputs, frame, fit, n, fe_sizes = numeric()) {
  if (is.null(inputs[["retained_original"]])) return(invisible(NULL))
  conn = inputs[["conn"]]
  aliases = paste0("__dbreg_saved_key", seq_along(inputs[["row_id"]]))
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  safe_keys = inputs[["row_id_safe"]]
  residual = fit_residual_sql(fit, frame, inputs[["yvar"]], conn)
  score_keys = paste(q(safe_keys), "AS", aliases)
  score_names = linear_score_names(inputs, fit)
  score_terms = if (length(score_names)) {
    terms = vapply(fit[["x"]], function(name) {
      if (name == "(Intercept)") "1.0" else fit_column_sql(fit, frame, name, conn)
    }, character(1))
    paste(terms, "AS", score_names)
  }
  selected = paste(frame[["cte_sql"]], "SELECT",
    paste(c(score_keys, paste(residual, "AS saved_residual"), score_terms), collapse = ","),
    "FROM", frame[["cte_name"]])
  finalize_retained_linear_rows(inputs, selected, aliases, n, fit[["rss"]],
    fit, score_names, fe_sizes)
}

finalize_retained_linear_rows = function(inputs, selected, aliases, n, fitted_rss, fit, score_names, fe_sizes) {
  original = inputs[["retained_original"]]
  if (is.null(original)) return(invisible(NULL))
  conn = inputs[["conn"]]
  resources = inputs[["resources"]]
  keys = inputs[["row_id"]]
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  y = inputs[["offset_outcome"]]
  if (is.null(y)) y = restore_names(inputs[["yvar"]], inputs[["source_map"]])
  w = restore_names(inputs[["weights"]], inputs[["source_map"]])
  outcome = sql_cast_double(paste0("b.", q(y)), conn)
  weight = if (length(w)) sql_cast_double(paste0("b.", q(w)), conn) else "1.0"
  names = c(aliases, "__dbreg_saved_fit", "__dbreg_saved_residual",
    "__dbreg_saved_weight", "__dbreg_estimated", score_names)
  fields = c(paste0("b.", q(keys), " AS ", aliases),
    paste(outcome, "-a.saved_residual AS __dbreg_saved_fit"),
    "a.saved_residual AS __dbreg_saved_residual",
    paste(weight, "AS __dbreg_saved_weight"),
    paste0("a.", aliases[1L], " IS NOT NULL AS __dbreg_estimated"),
    if (length(score_names)) paste0("a.", score_names))
  joins = paste0("b.", q(keys), "=a.", aliases)
  sql = paste0("SELECT ", paste(fields, collapse = ","), " FROM ", original,
    " b LEFT JOIN (", selected, ") a ON ", paste(joins, collapse = " AND "))
  create = fe_projection_tables(conn, resources)
  table = create(sql, names)
  rss = sql_stable_sum(conn, "__dbreg_saved_weight*POWER(__dbreg_saved_residual,2)")
  checks = DBI::dbGetQuery(conn, paste("SELECT COUNT(*) AS n, COUNT(__dbreg_saved_residual) AS estimated,",
    rss, "AS rss,", sql_stable_sum(conn, "POWER(__dbreg_saved_residual,2)"), "AS rss_unweighted FROM", table))
  total = DBI::dbGetQuery(conn, paste("SELECT COUNT(*) AS n FROM", original))[["n"]]
  if (!isTRUE(as.numeric(checks[["n"]]) == as.numeric(total)) ||
      !isTRUE(as.numeric(checks[["estimated"]]) == n) || !is.finite(checks[["rss"]]) ||
      abs(checks[["rss"]] - fitted_rss) > 1e-7 * max(1, abs(fitted_rss))) {
    stop("Retained linear rows do not reproduce the fitted sample and RSS.", call. = FALSE)
  }
  names = resources[["tables"]]
  name = names[match(table, q(names))]
  r_input = is.data.frame(inputs[["data"]]) && !length(inputs[["public_row_id"]])
  row_id = if (length(inputs[["public_row_id"]])) inputs[["public_row_id"]] else "row_id"
  inference = NULL
  if (linear_keeps_scores(inputs)) {
    original_name = names[match(original, q(names))]
    name = c(name, original_name)
    inference = list(original = original, source_keys = keys,
      source_classes = inputs[["retained_source_classes"]],
      design = score_names, bread = fit[["bread"]],
      sandwich_bread = fit[["sandwich_bread"]], weight_scale = fit[["weight_scale"]],
      fe_sizes = fe_sizes, default_cluster = restore_names(inputs[["cluster_var"]], inputs[["source_map"]]))
  }
  inputs[["linear_retained_state"]] = list(name = name, inference = inference,
    rss_unweighted = as.numeric(checks[["rss_unweighted"]]),
    rows = list(table = table, keys = aliases, row_id = row_id, r_input = r_input, input_n = as.numeric(total)))
  invisible(NULL)
}

#' Fitted values and residuals from retained linear rows
#' @param object A database linear fit with retained rows.
#' @param ... Unused.
#' @return A numeric vector for R inputs or a keyed result for SQL inputs.
#' @export
fitted.dbreg = function(object, ...) {
  check_unused_args(list(...))
  require_inference_rows(object)
  retained_row_values(object, "__dbreg_saved_fit", ".fitted")
}

#' @rdname fitted.dbreg
#' @param type Raw response or weighted Pearson or deviance residuals.
#' @param na.rm Omit excluded input rows.
#' @export
residuals.dbreg = function(object, type = c("response", "deviance", "pearson"),
                            na.rm = TRUE, ...) {
  check_unused_args(list(...))
  require_inference_rows(object)
  type = match.arg(type)
  if (!is.logical(na.rm) || length(na.rm) != 1L || is.na(na.rm)) {
    stop("`na.rm` must be TRUE or FALSE.", call. = FALSE)
  }
  value = "__dbreg_saved_residual"
  if (type != "response") value = paste0(value, "*SQRT(__dbreg_saved_weight)")
  retained_row_values(object, value, ".residual", na.rm)
}

#' @rdname fitted.dbreg
#' @importFrom stats case.names
#' @export
case.names.dbreg = function(object, ...) {
  check_unused_args(list(...))
  require_inference_rows(object)
  retained_row_values(object, NULL, NULL)
}

retain_compressed_linear_rows = function(inputs, from, groups, ids, fitted,
                                          outcome, center, rss, n, keys, fit, design, fe_sizes) {
  if (is.null(inputs[["retained_original"]])) return(invisible(NULL))
  conn = inputs[["conn"]]
  limits = projection_limits()
  if (length(ids) > limits[["design_entries"]]) {
    projection_abort("system_limit", "Compressed fitted values exceed the retained lookup budget.")
  }
  order = order(ids)
  if (!identical(as.numeric(ids[order]), as.numeric(seq_along(ids)))) {
    stop("Invalid compressed fitted cell identities.", call. = FALSE)
  }
  budget = list2env(list(upload_bytes = 0), parent = emptyenv())
  score_names = linear_score_names(inputs, fit)
  if (length(ids) * (1 + length(score_names)) > limits[["design_entries"]]) {
    projection_abort("system_limit", "Compressed score coordinates exceed the retained lookup budget.")
  }
  values = matrix(fitted[order], ncol = 1L)
  if (length(score_names)) values = cbind(values, as.matrix(design[order, , drop = FALSE]))
  lookup = fe_coefficient_sql(conn, values, budget, limits)
  aliases = paste0("__dbreg_saved_key", seq_along(keys))
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  group_sql = q(groups)
  cells = compressed_cell_index_sql(from, group_sql)
  join = if (length(groups)) {
    paste0(" JOIN (", cells, ") c ON ",
      paste(paste0("r.", group_sql, "=c.", group_sql), collapse = " AND "))
  } else {
    paste0(" CROSS JOIN (", cells, ") c")
  }
  selected = paste0("SELECT ", paste(paste0("r.", q(keys), " AS ", aliases), collapse = ","),
    ",(", sql_cast_double(paste0("r.", q(outcome)), conn), "-(",
    dbglm_number(center), "))-f.c1 AS saved_residual",
    if (length(score_names)) paste0(",", paste(paste0("f.c", seq_along(score_names) + 1L,
      " AS ", score_names), collapse = ",")),
    " FROM (SELECT * ", from, ") r", join, " JOIN (", lookup, ") f ON c.dbreg_cell=f.id")
  finalize_retained_linear_rows(inputs, selected, aliases, n, rss, fit, score_names, fe_sizes)
}

linear_keeps_scores = function(inputs) {
  isTRUE(inputs[["retain"]] %in% c("inference", "all"))
}

linear_score_names = function(inputs, fit) {
  if (linear_keeps_scores(inputs) && length(fit[["x"]])) paste0("__dbreg_s", seq_along(fit[["x"]]))
  else character()
}
