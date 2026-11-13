# Restricted predicate AST: never evaluate calls or resolve objects in the
# caller's environment. The database evaluates the compiled expression.
subset_literals = function(conn, values) {
  if (is.character(values)) return(ifelse(is.na(values), "NULL", category_string_literals(conn, values)))
  if (is.logical(values)) return(ifelse(is.na(values), "NULL", ifelse(values, "TRUE", "FALSE")))
  if (!is.numeric(values)) stop("Unsupported subset literal.", call. = FALSE)
  vapply(values, function(value) {
    if (is.nan(value)) return("CAST('NaN' AS DOUBLE)")
    if (is.na(value)) return("NULL")
    if (is.infinite(value)) return(paste0("CAST('", if (value < 0) "-" else "", "Infinity' AS DOUBLE)"))
    sql_cast_double(sprintf("%.17e", value), conn)
  }, character(1))
}

subset_variables = function(subset) {
  if (is.null(subset)) return(character())
  if (is.character(subset) && length(subset) == 1L && !is.na(subset) && nzchar(subset)) return(subset)
  if (!inherits(subset, "formula") || length(subset) != 2L) stop("`subset` must name a logical column or be a one-sided predicate formula.", call. = FALSE)
  all.vars(subset)
}

sql_subset_predicate = function(conn, subset, schema, to_safe = identity) {
  if (is.null(subset)) return(NULL)
  backend = detect_backend(conn)
  if (backend[["name"]] != "duckdb" && !isTRUE(backend[["is_databricks"]])) {
    stop("Subset predicates support DuckDB and native Databricks.", call. = FALSE)
  }
  subset_variables(subset)
  expr = if (is.character(subset)) as.name(subset) else subset[[2L]]
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  missing = function(value) paste0(
    "(", value[["sql"]], " IS NULL",
    if (value[["type"]] == "numeric") paste0(" OR isnan(CAST(", value[["sql"]], " AS DOUBLE))"), ")"
  )
  node = function(sql, type) list(sql = sql, type = type)
  literal = function(value) {
    if (length(value) != 1L || !(is.numeric(value) || is.character(value) || is.logical(value))) stop("Unsupported subset literal.", call. = FALSE)
    node(subset_literals(conn, value), if (is.character(value)) "character" else if (is.logical(value)) "logical" else "numeric")
  }
  constant = function(expr) {
    if (is.call(expr)) {
      op = as.character(expr[[1L]])
      args = as.list(expr)[-1L]
      if (op == "c") return(unlist(lapply(args, constant)))
      if (op == ":" && length(args) == 2L) {
        ends = c(constant(args[[1L]]), constant(args[[2L]]))
        if (is.numeric(ends) && length(ends) == 2L && all(is.finite(ends)) &&
            all(ends == floor(ends)) && abs(diff(ends)) < 1e4) return(ends[1L]:ends[2L])
      }
      if (op %in% c("+", "-") && length(args) == 1L) {
        value = constant(args[[1L]])
        if (length(value) == 1L && is.numeric(value)) return(if (op == "-") -value else value)
      }
      stop("Subset %in% requires literal values.", call. = FALSE)
    }
    if (is.symbol(expr)) stop("Subset %in% requires literal values.", call. = FALSE)
    expr
  }
  visit = function(expr) {
    if (is.symbol(expr)) {
      column = to_safe(as.character(expr))
      if (!column %in% names(schema)) stop("Unknown subset column: ", as.character(expr), call. = FALSE)
      value = schema[[column]]
      type = if (is.logical(value)) "logical" else if (is.numeric(value) || inherits(value, "integer64")) "numeric" else if (is.character(value) || is.factor(value)) "character" else "other"
      return(node(q(column), type))
    }
    if (!is.call(expr)) return(literal(expr))
    op = as.character(expr[[1L]])
    args = as.list(expr)[-1L]
    if (op == "(" && length(args) == 1L) return(visit(args[[1L]]))
    if (op == "is.na" && length(args) == 1L) return(node(missing(visit(args[[1L]])), "logical"))
    if (op == "!" && length(args) == 1L) {
      value = visit(args[[1L]])
      if (value[["type"]] != "logical") stop("Subset ! requires a logical operand.", call. = FALSE)
      return(node(paste0("(NOT ", value[["sql"]], ")"), "logical"))
    }
    if (op == "%in%" && length(args) == 2L) {
      value = visit(args[[1L]])
      constants = constant(args[[2L]])
      if (!length(constants)) return(node("FALSE", "logical"))
      present = constants[!is.na(constants)]
      type = if (is.character(constants)) "character" else if (is.logical(constants)) "logical" else "numeric"
      if (length(present) && value[["type"]] != type) stop("Subset %in% requires matching column and literal types.", call. = FALSE)
      membership = if (length(present)) paste0("COALESCE(", value[["sql"]], " IN (", paste(subset_literals(conn, present), collapse = ","), "),FALSE)") else "FALSE"
      has_nan = is.numeric(constants) && any(is.nan(constants))
      has_na = any(is.na(constants) & !(if (is.numeric(constants)) is.nan(constants) else FALSE))
      return(node(paste0(
        "(CASE WHEN ", value[["sql"]], " IS NULL THEN ", if (has_na) "TRUE" else "FALSE",
        if (value[["type"]] == "numeric") paste0(" WHEN isnan(CAST(", value[["sql"]], " AS DOUBLE)) THEN ", if (has_nan) "TRUE" else "FALSE"),
        " ELSE ", membership, " END)"
      ), "logical"))
    }
    if (op %in% c("+", "-") && length(args) == 1L) {
      value = visit(args[[1L]])
      if (value[["type"]] != "numeric") stop("Subset arithmetic requires numeric operands.", call. = FALSE)
      return(node(paste0("(", op, value[["sql"]], ")"), "numeric"))
    }
    if (length(args) == 2L && op %in% c("+", "-", "*", "/", "^", "==", "!=", ">", "<", ">=", "<=", "&", "|")) {
      left = visit(args[[1L]])
      right = visit(args[[2L]])
      if (op %in% c("&", "|")) {
        if (left[["type"]] != "logical" || right[["type"]] != "logical") stop("Subset boolean operators require logical operands.", call. = FALSE)
        return(node(paste0("(", left[["sql"]], if (op == "&") " AND " else " OR ", right[["sql"]], ")"), "logical"))
      }
      if (op %in% c("+", "-", "*", "/", "^")) {
        if (left[["type"]] != "numeric" || right[["type"]] != "numeric") stop("Subset arithmetic requires numeric operands.", call. = FALSE)
        return(node(if (op == "^") paste0("POWER(", sql_cast_double(left[["sql"]], conn), ",", sql_cast_double(right[["sql"]], conn), ")") else paste0("(", sql_cast_double(left[["sql"]], conn), op, sql_cast_double(right[["sql"]], conn), ")"), "numeric"))
      }
      if (left[["type"]] != right[["type"]] && !left[["sql"]] %in% c("NULL", "NA") && !right[["sql"]] %in% c("NULL", "NA")) stop("Subset comparisons require matching operand types.", call. = FALSE)
      comparison = paste0("(", left[["sql"]], if (op == "==") "=" else op, right[["sql"]], ")")
      return(node(paste0(
        "(CASE WHEN ", missing(left), " OR ", missing(right),
        " THEN CAST(NULL AS BOOLEAN) ELSE ", comparison, " END)"
      ), "logical"))
    }
    stop("Unsupported subset expression: ", paste(deparse(expr), collapse = ""), call. = FALSE)
  }
  result = visit(expr)
  if (result[["type"]] != "logical") stop("Subset formulas must return a logical predicate.", call. = FALSE)
  result[["sql"]]
}

filter_source_predicate = function(from_statement, predicate) {
  if (is.null(predicate)) return(from_statement)
  paste0("FROM (SELECT * ", from_statement, ") dbreg_subset_source WHERE (", predicate, ")")
}

apply_source_subset = function(conn, from_statement, subset, to_safe = identity) {
  if (is.null(subset)) return(from_statement)
  schema = DBI::dbGetQuery(conn, sql_limit(paste("SELECT *", from_statement), 0, detect_backend(conn)[["name"]]))
  filter_source_predicate(from_statement, sql_subset_predicate(conn, subset, schema, to_safe))
}

# Singleton removal precedes strategy selection and compression. Keep the
# native tuples and outcome/offset together in one owned temp table.
linear_prune_singletons = function(inputs, from, fe, y) {
  conn = inputs[["conn"]]
  backend = detect_backend(conn)
  if (backend[["name"]] != "duckdb" && !isTRUE(backend[["is_databricks"]])) {
    stop("OLS singleton removal requires DuckDB or native Databricks.", call. = FALSE)
  }
  fields = names(DBI::dbGetQuery(conn, paste("SELECT *", from, "LIMIT 0")))
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  cols_sql = paste(q(fields), collapse = ", ")
  flags = safe_identifiers(c("__dbreg_single"), taken = fields,
    fold_case = folds_identifier_case(conn))
  flag = flags[[1L]]
  make_temp = function(sql) {
    tbl = sprintf("dbreg_prune_%s_%06d", format(Sys.time(), "%H%M%S"), sample.int(999999L, 1L))
    q_tbl = as.character(DBI::dbQuoteIdentifier(conn, tbl))
    DBI::dbExecute(conn, paste("CREATE TEMP TABLE", q_tbl, "AS", sql))
    inputs[["temp_tables"]] = unique(c(inputs[["temp_tables"]], q_tbl))
    q_tbl
  }
  drop_temp = function(q_tbl) {
    drop_kw = if (isTRUE(backend[["is_databricks"]])) "DROP TEMP TABLE IF EXISTS" else "DROP TABLE IF EXISTS"
    try(DBI::dbExecute(conn, paste(drop_kw, q_tbl)), silent = TRUE)
    inputs[["temp_tables"]] = setdiff(inputs[["temp_tables"]], q_tbl)
    invisible(NULL)
  }
  table = make_temp(paste("SELECT *", from))
  rounds = 0L
  removed_total = 0
  if (length(fe)) {
    repeat {
      singleton = joins = character()
      for (k in seq_along(fe)) {
        f = fe[k]
        s = paste0("__dbreg_p", k)
        joins = c(joins, paste0(
          " JOIN (SELECT ", f, " AS __dbreg_g, COUNT(*) AS n FROM ", table,
          " GROUP BY ", f, ") ", s, " ON __dbreg_t.", f, " = ", s, ".__dbreg_g"
        ))
        singleton = c(singleton, paste0(s, ".n = 1"))
      }
      flagged = paste0(
        "(SELECT __dbreg_t.*, (", paste(singleton, collapse = " OR "), ") AS ", flag,
        " FROM ", table, " __dbreg_t", paste(joins, collapse = ""), ") __dbreg_f"
      )
      n_rem = as.numeric(DBI::dbGetQuery(conn, paste0(
        "SELECT COUNT(*) FILTER (WHERE ", flag, ") AS removed FROM ", flagged
      ))[["removed"]])
      if (n_rem == 0) break
      removed_total = removed_total + n_rem
      next_table = make_temp(paste0("SELECT ", cols_sql, " FROM ", flagged, " WHERE NOT ", flag))
      drop_temp(table)
      table = next_table
      rounds = rounds + 1L
    }
  }
  n = as.numeric(DBI::dbGetQuery(conn, paste("SELECT COUNT(*) AS n FROM", table))[["n"]])
  if (!n) stop("No observations remain after singleton removal.", call. = FALSE)
  list(table = table, rounds = rounds, counts = c(singleton = removed_total, infinite_coef = 0))
}
