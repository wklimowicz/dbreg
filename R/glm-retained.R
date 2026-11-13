dbglm_retain_rows = function(conn, resources, sample, groups, fit, x, family, row_id, keep_source, keep_scores, prediction_tables = NULL, transform = NULL) {
  working = fit[["working"]]
  if (is.null(working)) {
    sql = dbglm_working_sql(sample[["base"]], x, fit[["beta"]], family[["link"]],
      conn = conn, cell = sample[["cell"]])
    working = paste0("(", sql, "SELECT * FROM working) retained_working")
  }
  keys = sample[["row_keys"]]
  values = c("__dbglm_saved_y", "__dbglm_saved_mu", "__dbglm_saved_eta", "__dbglm_saved_d", "__dbglm_saved_w", "__dbglm_saved_res")
  projected = fit[["projected"]]
  design = if (keep_scores && length(x)) {
    if (is.null(projected)) x else paste0("v", seq_along(x) + 1L)
  }
  if (!is.null(projected) && length(design)) {
    ids = sample[["score_keys"]]
    on = paste0("w.", ids, "=p.", ids, collapse = " AND ")
    working = paste0("(SELECT w.*,", paste0("p.", design, collapse = ","),
      " FROM (SELECT * FROM ", working, ") w JOIN (SELECT * FROM ", projected, ") p ON ", on, ") retained_projected")
  }
  saved_design = if (length(design)) paste0("__dbglm_saved_x", seq_along(design))
  outputs = c(keys, values, saved_design, "__dbglm_estimated")
  if (sample[["cell"]]) {
    join = paste(paste0("b.", groups, "=a.", groups), collapse = " AND ")
    original_fields = if (!is.null(projected) && length(design)) setdiff(design, groups)
    source = paste0("(SELECT b.*,a.mu,a.eta,a.d",
      if (length(original_fields)) paste0(",a.", original_fields, collapse = ""), " FROM ", sample[["original_base"]],
      " b JOIN (SELECT ", paste(c(groups, "mu", "eta", "d", original_fields), collapse = ","),
      " FROM ", working, ") a ON ", join, ") retained_original")
  } else {
    source = working
  }
  create = fe_projection_tables(conn, resources)
  y_shift = paste0(
    if (length(sample[["y_center"]]) && sample[["y_center"]] != 0) paste0("+(", dbglm_number(sample[["y_center"]]), ")") else "",
    if (isTRUE(sample[["has_offset"]])) "+__dbglm_saved_o" else ""
  )
  fitted = paste("SELECT", paste(c(keys, paste0("__dbglm_y", y_shift, " AS __dbglm_saved_y"),
    paste0("mu", y_shift, " AS __dbglm_saved_mu"), paste0("eta", y_shift, " AS __dbglm_saved_eta"),
    "d AS __dbglm_saved_d",
    paste0("__dbglm_a*(", dbglm_number(sample[["weight_scale"]]), ") AS __dbglm_saved_w"),
    "(__dbglm_y-mu) AS __dbglm_saved_res",
    if (length(design)) paste(design, "AS", saved_design)),
    collapse = ","), "FROM", source)
  original = as.character(DBI::dbQuoteIdentifier(conn, sample[["original_keys"]]))
  match_keys = paste0("b.", original, "=a.", keys)
  select = c(paste("b.", original, " AS ", keys, sep = ""), paste0("a.", c(values, saved_design)),
    paste0("a.", keys[1L], " IS NOT NULL AS __dbglm_estimated"))
  quoted = create(paste0("SELECT ", paste(select, collapse = ","), " FROM ", sample[["raw"]],
    " b LEFT JOIN (", fitted, ") a ON ", paste(match_keys, collapse = " AND ")), outputs)
  counts = DBI::dbGetQuery(conn, paste("SELECT COUNT(*) AS n,",
    "SUM(CASE WHEN __dbglm_estimated THEN 1 ELSE 0 END) AS estimated FROM", quoted))
  if (!isTRUE(as.numeric(counts[["n"]]) == sample[["input_n"]]) ||
      !isTRUE(as.numeric(counts[["estimated"]]) == sample[["n"]])) {
    stop("Retained GLM rows do not match the fitted sample.", call. = FALSE)
  }
  names = resources[["tables"]]
  name = names[match(quoted, as.character(DBI::dbQuoteIdentifier(conn, names)))]
  inference = if (keep_source) list(design = saved_design, H = if (keep_scores) fit[["H"]],
    transform = transform,
    weight_scale = sample[["weight_scale"]], source = sample[["raw"]],
    source_keys = sample[["original_keys"]], cluster_fields = sample[["cluster_fields"]],
    cluster_classes = sample[["cluster_classes"]], default_cluster = sample[["default_cluster"]],
    families = sample[["public_groups"]], sizes = sample[["sizes"]],
    slope_df = sum(sample[["sizes"]] * lengths(sample[["slopes"]])))
  owned = if (keep_source) c(name, names[match(sample[["raw"]],
    as.character(DBI::dbQuoteIdentifier(conn, names)))]) else name
  lease = retain_fit_resources(resources, c(owned, prediction_tables), session = resources[["parent_session"]])
  list(rows = list(table = quoted, keys = keys, row_id = if (length(row_id)) row_id else "row_id",
    input_n = sample[["input_n"]]),
    lease = lease, inference = inference)
}

retained_row_values = function(object, field, label, na.rm = TRUE) {
  state = object[["retained_rows"]]
  lease = object[["resource_state"]]
  conn = retained_fit_connection(lease)
  keys = state[["keys"]]
  fields = c(keys, if (length(field)) paste(field, "AS dbreg_value"))
  estimated_col = if (inherits(object, "dbglm")) "__dbglm_estimated" else "__dbreg_estimated"
  sql = paste("SELECT", paste(fields, collapse = ","),
    "FROM", state[["table"]], if (na.rm) paste("WHERE", estimated_col), "ORDER BY", paste(keys, collapse = ","))
  if (isTRUE(state[["r_input"]])) {
    rows = DBI::dbGetQuery(conn, sql)
    return(if (length(field)) rows[["dbreg_value"]] else as.integer(rows[[keys]]))
  }
  materialize_glm_rows(conn, lease[["session"]], sql, keys,
    state[["row_id"]], label)
}

materialize_glm_rows = function(conn, session, sql, keys, row_id, label) {
  resources = fit_resources()
  resources[["conn"]] = conn
  on.exit(release_fit_resources(resources))
  columns = c(keys, if (length(label)) "dbreg_value")
  table = fe_projection_tables(conn, resources)(sql, columns)
  name = resources[["tables"]][[1L]]
  lease = retain_fit_resources(resources, name, session = session)
  structure(list(table = table, keys = keys, columns = columns,
    row_id = row_id, label = label, resource_state = lease), class = "dbreg_glm_rows")
}

#' @rdname fitted.dbreg
#' @export
fitted.dbglm = function(object, ...) {
  check_unused_args(list(...))
  require_inference_rows(object)
  retained_row_values(object, "__dbglm_saved_mu", ".fitted")
}

#' Inspect keyed fitted values from a database regression
#' @param x A retained result from `fitted()` or `predict()`.
#' @param ... Unused.
#' @return A data frame with source keys and fitted values from collection, the
#'   values in key order from `as.numeric()`, or invisibly the input from printing.
#' @export
as.data.frame.dbreg_glm_rows = function(x, ...) collect_glm_rows(x)

#' @rdname as.data.frame.dbreg_glm_rows
#' @export
as.double.dbreg_glm_rows = function(x, ...) {
  if (is.null(x[["label"]])) stop("These keyed rows carry no values; use as.data.frame().", call. = FALSE)
  data = collect_glm_rows(x)
  as.double(data[[ncol(data)]])
}

collect_glm_rows = function(x, n = NULL) {
  conn = retained_fit_connection(x[["resource_state"]])
  keys = x[["keys"]]
  columns = paste0("s.", x[["columns"]])
  backend = detect_backend(conn)
  if (isTRUE(backend[["is_databricks"]])) {
    fields = paste(sprintf("typeof(%s) AS %s", keys, keys), collapse = ",")
    types = DBI::dbGetQuery(conn, paste("SELECT", fields, "FROM", x[["table"]], "LIMIT 1"))
    exact = which(vapply(types, function(type) length(type) && tolower(type[1L]) == "bigint", logical(1)))
    columns[exact] = paste(sql_cast_text(conn, columns[exact]), "AS", keys[exact])
  }
  sql = paste("SELECT", paste(columns, collapse = ","), "FROM", x[["table"]],
    "s ORDER BY", paste(paste0("s.", keys), collapse = ","))
  if (!is.null(n)) sql = sql_limit(sql, n, backend[["name"]])
  data = read_keyed_query(conn, sql, keys, exact_keys = TRUE)
  names(data) = make.unique(c(x[["row_id"]], x[["label"]]))
  data
}

#' @rdname as.data.frame.dbreg_glm_rows
#' @param n Maximum number of rows to collect.
#' @importFrom utils head
#' @export
head.dbreg_glm_rows = function(x, n = 6L, ...) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 0 || n != floor(n)) {
    stop("`n` must be one nonnegative whole number.", call. = FALSE)
  }
  collect_glm_rows(x, n)
}

#' @rdname as.data.frame.dbreg_glm_rows
#' @param con A keyed result.
#' @export
close.dbreg_glm_rows = function(con, ...) close_retained_lease(con[["resource_state"]])

#' @rdname fitted.dbreg
#' @export
case.names.dbglm = function(object, ...) {
  check_unused_args(list(...))
  require_inference_rows(object)
  retained_row_values(object, NULL, NULL)
}

#' @rdname fitted.dbreg
#' @export
residuals.dbglm = function(object, type = c("response", "deviance", "pearson", "working"),
                            na.rm = TRUE, ...) {
  check_unused_args(list(...))
  type = match.arg(type)
  require_inference_rows(object)
  if (!is.logical(na.rm) || length(na.rm) != 1L || is.na(na.rm)) {
    stop("`na.rm` must be TRUE or FALSE.", call. = FALSE)
  }
  y = "__dbglm_saved_y"
  mu = "__dbglm_saved_mu"
  diff = "__dbglm_saved_res"
  family = object[["family"]][["family"]]
  variance = switch(family, gaussian = "1.0", poisson = mu,
    binomial = paste0(mu, "*(1.0-", mu, ")"))
  deviance = switch(family, gaussian = paste0("POWER(", diff, ",2)"),
    poisson = paste0("2.0*(", sql_poisson_bd0(y, mu), ")"),
    binomial = paste0("-2.0*(CASE WHEN ", y, "=1 THEN LN(", mu, ") ELSE LN(1.0-", mu, ") END)"))
  value = switch(type, response = diff,
    deviance = paste0("SIGN(", diff, ")*SQRT(GREATEST(0.0,__dbglm_saved_w*(", deviance, ")))"),
    pearson = paste0(diff, "*SQRT(__dbglm_saved_w)/SQRT(", variance, ")"),
    working = paste0(diff, "/__dbglm_saved_d"))
  retained_row_values(object, paste0("CASE WHEN __dbglm_estimated THEN ", value, " END"),
    ".residual", na.rm = na.rm)
}

# Stable Poisson deviance around y=mu. The odd-power series avoids subtracting
# nearly equal O(y) values when the desired deviance is O((y-mu)^2/y).
sql_poisson_bd0 = function(y, mu) {
  v = paste0("((", y, "-", mu, ")/(", y, "+", mu, "))")
  polynomial = "0.0"
  for (j in 12:1) polynomial = paste0("(", dbglm_number(1 / (2 * j + 1)), "+POWER(", v, ",2)*", polynomial, ")")
  series = paste0("(", y, "-", mu, ")*", v, "+2*", y, "*POWER(", v, ",3)*", polynomial)
  paste0(
    "CASE WHEN ", y, "=0 THEN ", mu, " WHEN ABS(", y, "-", mu, ")<0.1*(", y, "+", mu,
    ") THEN (", series, ") ELSE ", y, "*(LN(", y, ")-LN(", mu, "))-(", y, "-", mu, ") END"
  )
}

require_inference_rows = function(object) {
  if (!isTRUE(object[["retain"]] %in% c("inference", "all"))) {
    stop("Fit with retain = 'inference' or 'all' to access estimation rows.", call. = FALSE)
  }
}

#' @rdname as.data.frame.dbreg_glm_rows
#' @export
print.dbreg_glm_rows = function(x, ...) {
  cat("Database rows; keys: ", paste(x[["row_id"]], collapse = ", "),
    ". Use head() to preview.\n", sep = "")
  invisible(x)
}
