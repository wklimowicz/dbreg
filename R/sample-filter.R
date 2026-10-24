# Restricted, typed sample predicates shared by the native and adapter APIs.
subset_variables = function(subset) {
  if (is.null(subset)) return(character())
  if (is.character(subset) && length(subset) == 1L && !is.na(subset) && nzchar(subset)) return(subset)
  if (!inherits(subset, 'formula') || length(subset) != 2L) stop('`subset` must name a logical column or be a one-sided predicate formula.', call. = FALSE)
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
  missing = function(value) paste0('(', value[['sql']], ' IS NULL',
    if (value[['type']] == 'numeric') paste0(' OR isnan(CAST(', value[['sql']], ' AS DOUBLE))'), ')')
  node = function(sql, type) list(sql = sql, type = type)
  literal = function(value) {
    if (length(value) != 1L || !(is.numeric(value) || is.character(value) || is.logical(value))) stop('Unsupported subset literal.', call. = FALSE)
    node(sql_value_literals(conn, value), if (is.character(value)) 'character' else if (is.logical(value)) 'logical' else 'numeric')
  }
  constant = function(expr) {
    if (is.call(expr)) {
      op = as.character(expr[[1L]])
      args = as.list(expr)[-1L]
      if (op == 'c') return(unlist(lapply(args, constant)))
      if (op %in% c('+','-') && length(args) == 1L) {
        value = constant(args[[1L]])
        if (length(value) == 1L && is.numeric(value)) return(if (op == '-') -value else value)
      }
      stop('Subset %in% requires literal values.', call. = FALSE)
    }
    if (is.symbol(expr)) stop('Subset %in% requires literal values.', call. = FALSE)
    expr
  }
  visit = function(expr) {
    if (is.symbol(expr)) {
      column = to_safe(as.character(expr))
      if (!column %in% names(schema)) stop('Unknown subset column: ', as.character(expr), call. = FALSE)
      value = schema[[column]]
      type = if (is.logical(value)) 'logical' else if (is.numeric(value) || inherits(value, 'integer64')) 'numeric' else if (is.character(value) || is.factor(value)) 'character' else 'other'
      return(node(q(column), type))
    }
    if (!is.call(expr)) return(literal(expr))
    op = as.character(expr[[1L]])
    args = as.list(expr)[-1L]
    if (op == '(' && length(args) == 1L) return(visit(args[[1L]]))
    if (op == 'is.na' && length(args) == 1L) return(node(missing(visit(args[[1L]])), 'logical'))
    if (op == '!' && length(args) == 1L) {
      value = visit(args[[1L]])
      if (value[['type']] != 'logical') stop('Subset ! requires a logical operand.', call. = FALSE)
      return(node(paste0('(NOT ', value[['sql']], ')'), 'logical'))
    }
    if (op == '%in%' && length(args) == 2L) {
      value = visit(args[[1L]])
      constants = constant(args[[2L]])
      if (!length(constants)) return(node('FALSE', 'logical'))
      present = constants[!is.na(constants)]
      type = if (is.character(constants)) 'character' else if (is.logical(constants)) 'logical' else 'numeric'
      if (length(present) && value[['type']] != type) stop('Subset %in% requires matching column and literal types.', call. = FALSE)
      membership = if (length(present)) paste0('COALESCE(', value[['sql']], ' IN (', paste(sql_value_literals(conn, present), collapse = ','), '),FALSE)') else 'FALSE'
      has_nan = is.numeric(constants) && any(is.nan(constants))
      has_na = any(is.na(constants) & !(if (is.numeric(constants)) is.nan(constants) else FALSE))
      return(node(paste0('(CASE WHEN ', value[['sql']], ' IS NULL THEN ', if (has_na) 'TRUE' else 'FALSE',
        if (value[['type']] == 'numeric') paste0(' WHEN isnan(CAST(', value[['sql']], ' AS DOUBLE)) THEN ', if (has_nan) 'TRUE' else 'FALSE'),
        ' ELSE ', membership, ' END)'), 'logical'))
    }
    if (op %in% c('+','-') && length(args) == 1L) {
      value = visit(args[[1L]])
      if (value[['type']] != 'numeric') stop('Subset arithmetic requires numeric operands.', call. = FALSE)
      return(node(paste0('(', op, value[['sql']], ')'), 'numeric'))
    }
    if (length(args) == 2L && op %in% c('+','-','*','/','^','==','!=','>','<','>=','<=','&','|')) {
      left = visit(args[[1L]]); right = visit(args[[2L]])
      if (op %in% c('&','|')) {
        if (left[['type']] != 'logical' || right[['type']] != 'logical') stop('Subset boolean operators require logical operands.', call. = FALSE)
        return(node(paste0('(',left[['sql']],if(op=='&') ' AND ' else ' OR ',right[['sql']],')'),'logical'))
      }
      if (op %in% c('+','-','*','/','^')) {
        if (left[['type']] != 'numeric' || right[['type']] != 'numeric') stop('Subset arithmetic requires numeric operands.', call. = FALSE)
        return(node(if(op=='^') paste0('POWER(',left[['sql']],',',right[['sql']],')') else paste0('(',left[['sql']],op,right[['sql']],')'),'numeric'))
      }
      if (left[['type']] != right[['type']] && !left[['sql']] %in% c('NULL', 'NA') && !right[['sql']] %in% c('NULL', 'NA')) stop('Subset comparisons require matching operand types.', call. = FALSE)
      comparison = paste0('(',left[['sql']],if(op=='==') '=' else op,right[['sql']],')')
      return(node(paste0('(CASE WHEN ', missing(left), ' OR ', missing(right),
        ' THEN CAST(NULL AS BOOLEAN) ELSE ', comparison, ' END)'), 'logical'))
    }
    stop('Unsupported subset expression: ', paste(deparse(expr), collapse=''), call. = FALSE)
  }
  result = visit(expr)
  if (result[['type']] != 'logical') stop('Subset formulas must return a logical predicate.', call. = FALSE)
  result[['sql']]
}

filter_source_predicate = function(from_statement, predicate) {
  if (is.null(predicate)) return(from_statement)
  paste0("FROM (SELECT * ", from_statement, ") dbreg_subset_source WHERE (", predicate, ")")
}

linear_prune_source = function(inputs, from_statement, fe, y) {
  conn = inputs[["conn"]]
  backend = detect_backend(conn)
  if (backend[["name"]] != "duckdb" && !isTRUE(backend[["is_databricks"]])) {
    stop("OLS singleton pruning supports DuckDB and native Databricks.", call. = FALSE)
  }
  resources = inputs[["sample_resources"]]
  resources[["conn"]] = conn
  schema = DBI::dbGetQuery(conn, sql_limit(paste("SELECT *", from_statement), 0, backend[["name"]]))
  flags = safe_identifiers(c("__dbglm_single", "__dbglm_infinite"), taken = names(schema), fold_case = folds_identifier_case(conn))
  manager = dbglm_table_manager(conn, resources, backend)
  source = paste0("(SELECT * ", from_statement, ")")
  pruning = dbglm_prune(conn, source, fe, y, "singletons", manager[["create"]], manager[["drop"]], flags = flags)
  if (pruning[["rounds"]] > 0 && as.numeric(dbglm_get_query(conn,
    paste("SELECT COUNT(*) AS n FROM", pruning[["table"]]))[["n"]]) == 0) {
    stop("No observations remain after singleton removal.", call. = FALSE)
  }
  pruning
}
