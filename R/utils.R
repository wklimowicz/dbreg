# Internal utility functions for dbreg
# These are general-purpose helpers used across multiple strategies

#' Copy named values between environments
#'
#' @param source Source environment to read from
#' @param target Target environment to write to
#' @param keys Character vector of names to copy. If NULL, copies all.
#' @keywords internal
env2env = function(source, target, keys = NULL) {
  if (is.null(keys)) {
    keys = ls(source, all.names = TRUE)
  }
  for (nm in keys) {
    assign(nm, source[[nm]], envir = target)
  }
}

#' Convert internal interaction notation to standard R notation
#' @keywords internal
standardize_coef_names = function(x) gsub("_x_", ":", x, fixed = TRUE)

#' Generate coefficient table from estimates and vcov matrix
#' @keywords internal
gen_coeftable = function(betahat, vcov_mat, df_residual) {
  coefs = as.numeric(betahat)
  names(coefs) = standardize_coef_names(rownames(betahat))
  ses = sqrt(Matrix::diag(vcov_mat))
  tstats = coefs / ses
  pvals = 2 * pt(-abs(tstats), df_residual)
  cbind(estimate = coefs, std.error = ses, statistic = tstats, p.value = pvals)
}

#' Detect and handle collinearity in design matrix
#' 
#' Uses QR decomposition with pivoting to identify rank deficiency.
#' Returns reduced XtX/Xty matrices and lists of kept/dropped variables.
#' 
#' @keywords internal
detect_collinearity = function(XtX, Xty, tol = 1e-10, verbose = FALSE) {
  p = ncol(XtX)
  var_names = colnames(XtX)
  
  qr_decomp = qr(XtX, tol = tol)
  rank = qr_decomp[["rank"]]
  
  if (rank < p) {
    keep_idx = qr_decomp[["pivot"]][seq_len(rank)]
    drop_idx = qr_decomp[["pivot"]][(rank + 1):p]
    drop_names = var_names[drop_idx]
    keep_names = var_names[keep_idx]
    
    if (verbose) {
      message(sprintf(
        "[dbreg] %d variable(s) removed due to collinearity: %s",
        length(drop_names),
        paste(drop_names, collapse = ", ")
      ))
    }
    
    list(
      XtX = XtX[keep_idx, keep_idx, drop = FALSE],
      Xty = Xty[keep_idx, , drop = FALSE],
      keep_names = keep_names,
      drop_names = drop_names,
      collinear = TRUE
    )
  } else {
    list(
      XtX = XtX,
      Xty = Xty,
      keep_names = var_names,
      drop_names = character(0),
      collinear = FALSE
    )
  }
}

#' Solve linear system using Cholesky with QR fallback
#' @keywords internal
solve_with_fallback = function(XtX, Xty) {
  Rch = tryCatch(chol(XtX), error = function(e) NULL)
  if (is.null(Rch)) {
    # Cholesky failed, use QR fallback
    qr_decomp = qr(XtX)
    betahat = qr.solve(qr_decomp, Xty)
    XtX_inv = qr.solve(qr_decomp, diag(ncol(XtX)))
  } else {
    # Cholesky succeeded
    betahat = backsolve(Rch, forwardsolve(Matrix::t(Rch), Xty))
    XtX_inv = chol2inv(Rch)
  }
  dimnames(XtX_inv) = dimnames(XtX)
  list(betahat = betahat, XtX_inv = XtX_inv)
}


#' Set up database connection and data source
#'
#' Handles the precedence logic for data sources (table > data > path),
#' connection inference from tbl_lazy objects, and DuckDB registration.
#'
#' @param conn Database connection, or NULL to create a DuckDB connection

#' @param table Character table name or tbl_lazy object
#' @param data R data.frame to register with DuckDB
#' @param path Character path to file(s) on disk
#' @param caller Character string identifying the calling function for error messages
#'   and temp table naming (e.g., "dbreg" or "dbbinsreg")
#'
#' @return A list containing:
#'   \item{conn}{The database connection (created or provided)}
#'   \item{table_name}{The table name to query from}
#'   \item{from_statement}{SQL FROM clause (for dbreg-style queries)}
#'   \item{own_conn}{Logical, TRUE if we created the connection}
#'   \item{registered_table}{Name of registered table if data.frame was provided, else NULL}
#'   \item{is_duckdb}{Logical, TRUE if the connection is DuckDB}
#'
#' @keywords internal
setup_db_connection = function(conn, table, data, path, caller = "dbreg") {
  own_conn = FALSE
  registered_table = NULL
  is_duckdb = FALSE
  table_name = NULL
  from_statement = NULL


  # If table is tbl_lazy and conn is NULL, try to infer conn from the table
  if (is.null(conn) && inherits(table, "tbl_lazy")) {
    inferred_con = tryCatch(dbplyr::remote_con(table), error = function(e) NULL)
    if (!is.null(inferred_con) && dbIsValid(inferred_con)) {
      conn = inferred_con
    } else {
      stop(
        "Could not extract a valid database connection from the provided tbl_lazy. ",
        "The connection may be closed or invalid. ",
        "Either provide `conn` explicitly or ensure the tbl_lazy has an active connection."
      )
    }
  }

  # Create default connection if still NULL (for data.frame or path inputs)
  if (is.null(conn)) {
    conn = dbConnect(duckdb(), shutdown = TRUE)
    own_conn = TRUE
    is_duckdb = TRUE
  } else {
    # Check if user-provided connection is DuckDB
    backend_info = detect_backend(conn)
    is_duckdb = (backend_info[["name"]] == "duckdb")
  }

  # Process data source with precedence: table > data > path
  if (!is.null(table)) {
    if (is.character(table)) {
      # Table name as string
      table_name = table
      from_statement = glue("FROM {table}")
    } else if (inherits(table, "tbl_lazy")) {
      # Lazy table: render SQL as subquery
      rendered_sql = tryCatch(dbplyr::sql_render(table), error = function(e) NULL)
      if (is.null(rendered_sql)) {
        stop("Failed to render SQL for provided tbl_lazy.")
      }
      from_statement = paste0("FROM (", rendered_sql, ") AS lazy_subquery")
      table_name = paste0("(", rendered_sql, ") AS lazy_subquery")
      # Verify connection is valid
      if (!dbIsValid(conn)) {
        stop(
          "Could not obtain a valid database connection. ",
          "Either provide `conn` explicitly or ensure the tbl_lazy has an active connection."
        )
      }
    } else {
      stop("`table` must be character or tbl_lazy object.")
    }
  } else if (!is.null(data)) {
    if (!inherits(data, "data.frame")) {
      stop("`data` must be data.frame.")
    }
    # Coerce to base data.frame (handles tibbles, data.tables, etc.)
    data = as.data.frame(data)
    # Register with DuckDB (requires DuckDB connection)
    if (!is_duckdb) {
      stop("In-memory data frames are only supported with DuckDB connections. ",
           "Use `table` or `path` for other backends.")
    }
    temp_name = sprintf("__%s_%s",
                        caller,
                        gsub("[^0-9]", "", format(Sys.time(), "%Y%m%d_%H%M%S_%OS3")))
    duckdb_register(conn, temp_name, data)
    registered_table = temp_name
    table_name = temp_name
    from_statement = paste("FROM", temp_name)
  } else if (!is.null(path)) {
    if (!is.character(path)) {
      stop("`path` must be character.")
    }
    # Handle path: check if it's already a read_* function call or raw path
    if (!(grepl("^read|^scan", path) && grepl("'", path))) {
      path = gsub('"', "'", path)
      from_statement = glue("FROM '{path}'")
      table_name = glue("'{path}'")
    } else {
      from_statement = glue("FROM {path}")
      table_name = path
    }
  } else {
    stop("Provide one of `table`, `data`, or `path`.")
  }

  list(
    conn = conn,
    table_name = table_name,
    from_statement = from_statement,
    own_conn = own_conn,
    registered_table = registered_table,
    is_duckdb = is_duckdb
  )
}


#' Parse regression formula
#'
#' Parses a formula using the Formula package and extracts components needed
#' for regression: outcome variable, RHS variables, term labels, and fixed effects.
#'
#' @param fml A formula object (will be converted to Formula)
#'
#' @return A list containing:
#'   \item{fml}{The Formula object}
#'   \item{yvar}{Character, the outcome variable name}
#'   \item{xvars}{Character vector of unique variable names on RHS part 1}
#'   \item{term_labels}{Character vector of term labels (preserves interactions like "a:b")}
#'   \item{has_interactions}{Logical, TRUE if any interaction terms present}
#'   \item{fe}{Character vector of fixed effect variable names, or NULL if none}
#'
#' @keywords internal
parse_regression_formula = function(fml) {
  fml = Formula(fml)

  validate_expression = function(expr, allowed_calls, component) {
    if (is.symbol(expr) || is.numeric(expr)) {
      return(invisible(NULL))
    }
    if (!is.call(expr)) {
      stop("Unsupported ", component, " expression: ", deparse(expr), call. = FALSE)
    }
    call_name = as.character(expr[[1]])
    if (!call_name %in% allowed_calls) {
      stop(
        "Unsupported ", component, " expression `", deparse(expr), "`. ",
        "Supported terms are column names and the operators ",
        paste(sprintf("`%s`", allowed_calls), collapse = ", "), ".",
        call. = FALSE
      )
    }
    for (arg in as.list(expr)[-1]) {
      validate_expression(arg, allowed_calls, component)
    }
    invisible(NULL)
  }

  # Extract outcome variable (LHS)
  yvar = all.vars(formula(fml, lhs = 1, rhs = 0))
  if (length(yvar) != 1) {
    stop("Exactly one outcome variable required.")
  }
  
  # Extract RHS part 1: variables and term structure
  rhs1 = formula(fml, lhs = 0, rhs = 1)
  validate_expression(rhs1[[2]], c("+", ":", "*", "("), "regressor")
  tt = terms(rhs1)
  term_labels = attr(tt, "term.labels")
  xvars = all.vars(rhs1)  # unique variable names
  has_interactions = any(grepl(":", term_labels))
  
  if (!length(xvars)) {
    stop("No regressors on RHS.")
  }
  
  # Extract fixed effects (RHS part 2, after |)
  fe = if (length(fml)[2] > 1) {
    rhs2 = formula(fml, lhs = 0, rhs = 2)
    validate_expression(rhs2[[2]], c("+", "("), "fixed-effect")
    all.vars(rhs2)
  } else {
    NULL
  }
  
  list(
    fml = fml,
    yvar = yvar,
    xvars = xvars,
    term_labels = term_labels,
    has_interactions = has_interactions,
    fe = fe
  )
}

#' Compute model-frame sample diagnostics
#' @keywords internal
compute_sample_diagnostics = function(conn, from_statement, yvar, xvars, fe,
                                       weights = NULL, cluster_var = NULL) {
  variables = c(yvar, xvars, fe, weights, cluster_var)
  roles = c(
    "outcome",
    rep("regressor", length(xvars)),
    rep("fixed_effect", length(fe)),
    rep("weight", length(weights)),
    rep("cluster", length(cluster_var))
  )
  aliases = paste0("n_missing_", seq_along(variables))
  terms = c(
    sql_count(conn, "n_input"),
    sprintf(
      "SUM(CASE WHEN %s IS NULL THEN 1 ELSE 0 END) AS %s",
      variables,
      aliases
    )
  )
  if (!is.null(weights)) {
    terms = c(
      terms,
      sprintf("SUM(CASE WHEN %s <= 0 THEN 1 ELSE 0 END) AS n_nonpositive_weight", weights)
    )
  }
  audit_sql = paste0(
    "SELECT\n  ", paste(terms, collapse = ",\n  "),
    "\nFROM (SELECT * ", from_statement, ") AS dbreg_input"
  )
  counts = dbGetQuery(conn, audit_sql)
  detail = data.frame(
    role = roles,
    variable = variables,
    n_missing = as.numeric(counts[1, aliases, drop = TRUE]),
    row.names = NULL
  )
  list(
    n_input = as.numeric(counts[["n_input"]]),
    n_nonpositive_weight = if (!is.null(weights)) {
      as.numeric(counts[["n_nonpositive_weight"]])
    } else {
      0
    },
    detail = detail,
    query_string = audit_sql
  )
}


#' Parse variance-covariance and clustering arguments
#'
#' Parses the vcov and cluster arguments to determine the variance-covariance
#' type and extract any clustering variable.
#'
#' @param vcov Character string ("iid", "hc1") or formula for clustering (e.g., ~firm)
#' @param cluster Optional separate cluster argument (formula or character)
#' @param valid_types Character vector of valid vcov type strings. Default is
#'   c("iid", "hc1"). Case-insensitive matching is used.
#'
#' @return A list containing:
#'   \item{vcov_type}{Character string: "iid", "hc1", or "cluster"}
#'   \item{cluster_var}{Character name of cluster variable, or NULL if not clustered}
#'
#' @keywords internal
parse_vcov_args = function(vcov, cluster = NULL, valid_types = c("iid", "hc1")) {
  vcov_type = NULL
  cluster_var = NULL

  validate_single_cluster = function(x) {
    if (length(x) != 1 || is.na(x) || !nzchar(x)) {
      stop("Only single-variable clustering is currently supported")
    }
    x
  }

  # Handle vcov: can be string or formula (for clustering)
  if (inherits(vcov, "formula")) {
    cluster_var = validate_single_cluster(all.vars(vcov))
    vcov_type = "cluster"
  } else if (is.character(vcov)) {
    vcov_type = tolower(vcov[1])
    vcov_type = match.arg(vcov_type, valid_types)
  } else if (!is.null(vcov)) {
    stop("vcov must be a character string or a formula for clustering")
  }
  
  # Handle separate cluster argument (overrides vcov if provided)
  if (!is.null(cluster)) {
    if (inherits(cluster, "formula")) {
      cluster_var = validate_single_cluster(all.vars(cluster))
    } else if (is.character(cluster)) {
      cluster_var = validate_single_cluster(cluster)
    } else {
      stop("cluster must be a formula (e.g., ~firm) or character string")
    }
    vcov_type = "cluster"
  }
  
  list(
    vcov_type = vcov_type,
    cluster_var = cluster_var
  )
}

# =============================================================================
# SQL dialect helpers
# =============================================================================

#' Check if backend supports COUNT_BIG
#' 
#' SQL Server and Azure SQL use COUNT_BIG for large tables to avoid overflow.
#' 
#' @param conn A DBI database connection object.
#' @return Logical: TRUE if backend supports COUNT_BIG, FALSE otherwise.
#' @keywords internal
backend_supports_count_big = function(conn) {
  info = try(dbGetInfo(conn), silent = TRUE)
  if (inherits(info, "try-error")) {
    return(FALSE)
  }
  dbms = tolower(paste(info[["dbms.name"]], collapse = " "))
  grepl("sql server|azure sql|microsoft sql server", dbms)
}

#' Interpret backend metadata
#'
#' Kept separate from `detect_backend()` so capability detection can be tested
#' without constructing connector-specific DBI connection objects.
#'
#' @param info A list returned by `DBI::dbGetInfo()`.
#' @param connection_string Optional connector connection-string metadata.
#' @param connection_class Character vector of connection classes.
#' @return List of backend name and capabilities.
#' @keywords internal
backend_info_from_metadata = function(
  info = list(),
  connection_string = NULL,
  connection_class = character()
) {
  flatten_metadata = function(x) {
    if (is.null(x)) {
      return(character())
    }
    if (is.list(x)) {
      return(unlist(x, recursive = TRUE, use.names = TRUE))
    }
    as.character(x)
  }

  info_flat = flatten_metadata(info)
  connection_flat = flatten_metadata(connection_string)
  metadata = tolower(paste(c(info_flat, connection_flat, connection_class), collapse = " "))
  dbms = tolower(paste(info[["dbms.name"]], collapse = " "))
  is_databricks = grepl("databricks", metadata)

  native_query = NA
  named_connection = connection_flat[!is.na(names(connection_flat))]
  if (length(named_connection)) {
    key_names = tolower(sub("^.*\\.", "", names(named_connection)))
    native_value = named_connection[key_names == "usenativequery"]
    if (length(native_value)) {
      native_value = tolower(as.character(native_value[[1]]))
      if (native_value %in% c("1", "true", "yes")) {
        native_query = TRUE
      } else if (native_value %in% c("0", "false", "no")) {
        native_query = FALSE
      }
    }
  }
  if (is.na(native_query) && length(connection_flat)) {
    connection_text = paste(connection_flat, collapse = ";")
    native_match = regmatches(
      connection_text,
      regexpr("(?i)useNativeQuery\\s*=\\s*([012]|true|false)", connection_text, perl = TRUE)
    )
    if (length(native_match) && nzchar(native_match)) {
      native_value = tolower(sub("^.*=\\s*", "", native_match))
      if (native_value %in% c("1", "true")) native_query = TRUE
      if (native_value %in% c("0", "false")) native_query = FALSE
    }
  }

  name = if ("duckdb_connection" %in% connection_class || grepl("duckdb", dbms)) {
    "duckdb"
  } else if (grepl("sql server|azure sql|microsoft sql server", dbms)) {
    "sqlserver"
  } else if (grepl("postgres", dbms)) {
    "postgres"
  } else if (grepl("sqlite", dbms)) {
    "sqlite"
  } else if (grepl("mysql|mariadb", dbms)) {
    "mysql"
  } else if (is_databricks || grepl("spark", dbms)) {
    "spark"
  } else if (grepl("athena|presto|trino", dbms)) {
    "athena"
  } else {
    "other"
  }

  list(
    name = name,
    supports_count_big = name == "sqlserver",
    is_databricks = is_databricks,
    supports_physical_temp_tables = is_databricks,
    native_query = native_query
  )
}

#' Detect SQL backend from connection
#' 
#' Returns backend name and capabilities for SQL dialect handling.
#' 
#' @param conn A DBI database connection object.
#' @return List with backend name and capabilities.
#' @keywords internal
detect_backend = function(conn) {
  info = try(dbGetInfo(conn), silent = TRUE)
  if (inherits(info, "try-error")) {
    info = list()
  }
  backend_info_from_metadata(
    info = info,
    connection_string = attr(conn, "connectionString"),
    connection_class = class(conn)
  )
}

#' Generate SQL COUNT expression with proper casting
#' 
#' Returns an aliased COUNT expression that handles:
#' - COUNT_BIG for SQL Server (avoids int overflow on large tables)
#' - CAST to BIGINT for other backends
#' - DISTINCT counting when needed
#' 
#' @param conn Database connection (used to detect backend)
#' @param alias Column alias for the count result
#' @param expr Expression to count (default "*")
#' @param distinct Logical; if TRUE, count distinct values
#' @return SQL expression string like "CAST(COUNT(*) AS BIGINT) AS n"
#' @keywords internal
sql_count = function(conn, alias, expr = "*", distinct = FALSE) {
  bd = detect_backend(conn)
  use_count_big = bd[["supports_count_big"]]
  if (distinct) {
    glue(
      "{if (use_count_big) paste0('COUNT_BIG(DISTINCT ', expr, ')') else paste0('CAST(COUNT(DISTINCT ', expr, ') AS BIGINT)')} AS {alias}"
    )
  } else {
    if (use_count_big) {
      glue("COUNT_BIG({expr}) AS {alias}")
    } else {
      glue("CAST(COUNT({expr}) AS BIGINT) AS {alias}")
    }
  }
}

#' Get simple COUNT expression for backend
#' 
#' Returns just the count expression without alias (for use in aggregations).
#' 
#' @param backend Backend name from detect_backend()
#' @return SQL expression: "COUNT_BIG(*)" for SQL Server, "COUNT(*)" otherwise
#' @keywords internal
sql_count_expr = function(backend) {
  if (backend == "sqlserver") "COUNT_BIG(*)" else "COUNT(*)"
}

#' Get NTILE window function expression
#' @param x_name Column name to partition by
#' @param n_bins Number of bins
#' @return SQL NTILE expression
#' @keywords internal
sql_ntile = function(x_name, n_bins) {
  glue("NTILE({n_bins}) OVER (ORDER BY {x_name})")
}

#' Get backend-specific RANDOM() expression
#' 
#' Different SQL backends use different syntax for random number generation.
#' 
#' @param backend Backend name from detect_backend()
#' @return SQL expression for random value (e.g., "RANDOM()", "NEWID()", "RAND()")
#' @keywords internal
sql_random = function(backend) {
  switch(backend,
    "duckdb" = "RANDOM()",
    "sqlserver" = "NEWID()",
    "postgres" = "RANDOM()",
    "sqlite" = "RANDOM()",
    "mysql" = "RAND()",
    "RANDOM()"  # default fallback
  )
}

#' Build efficient sample query (backend-specific)
#' 
#' Uses native SAMPLE/TABLESAMPLE when available, falls back to ORDER BY RANDOM().
#' 
#' @param select_clause The SELECT clause (e.g., "*" or "col1, col2")
#' @param table_name Table to sample from
#' @param where_clause WHERE clause without "WHERE" keyword (or NULL)
#' @param n Number of rows to sample
#' @param backend Backend name from detect_backend()
#' @return Complete SQL query for sampling
#' @keywords internal
sql_sample = function(select_clause, table_name, where_clause, n, backend) {
  where_sql = if (!is.null(where_clause)) paste("WHERE", where_clause) else ""
  
  if (backend == "duckdb") {
    # DuckDB: efficient reservoir sampling
    glue("SELECT {select_clause} FROM {table_name} {where_sql} USING SAMPLE {n} ROWS")
  } else if (backend == "spark") {
    # Spark SQL: TABLESAMPLE with ROWS
    glue("SELECT {select_clause} FROM {table_name} TABLESAMPLE ({n} ROWS) {where_sql}")
  } else if (backend %in% c("athena", "presto", "trino")) {
    # Athena/Presto/Trino: TABLESAMPLE BERNOULLI only supports percentages
    # Fall back to ORDER BY RANDOM() for exact row counts
    glue("SELECT {select_clause} FROM {table_name} {where_sql} ORDER BY RANDOM() LIMIT {n}")
  } else if (backend == "postgres") {
    # PostgreSQL: TABLESAMPLE BERNOULLI only supports percentages
    random_expr = sql_random(backend)
    glue("SELECT {select_clause} FROM {table_name} {where_sql} ORDER BY {random_expr} LIMIT {n}")
  } else if (backend == "sqlserver") {
    # SQL Server: TABLESAMPLE is block-level, not row-level; use TOP with NEWID()
    glue("SELECT TOP {n} {select_clause} FROM {table_name} {where_sql} ORDER BY NEWID()")
  } else {
    # SQLite, MySQL, others: ORDER BY RANDOM() LIMIT n
    random_expr = sql_random(backend)
    glue("SELECT {select_clause} FROM {table_name} {where_sql} ORDER BY {random_expr} LIMIT {n}")
  }
}

#' Apply row limit to a SQL query (backend-specific)
#' 
#' SQL Server uses TOP after SELECT, others use LIMIT at end.
#' 
#' @param query The SQL query string
#' @param n Number of rows to limit to
#' @param backend Backend name from detect_backend()
#' @return SQL query with appropriate LIMIT/TOP clause
#' @keywords internal
sql_limit = function(query, n, backend) {
  if (backend == "sqlserver") {
    # SQL Server uses TOP n after SELECT
    sub("^SELECT", paste("SELECT TOP", n), query, ignore.case = TRUE)
  } else {
    # Others use LIMIT at the end
    paste(query, "LIMIT", n)
  }
}

#' Generate a temp table name with optional SQL Server prefix
#' 
#' SQL Server uses # prefix for local temp tables.
#' 
#' @param base_name Base name for the temp table
#' @param backend Backend name from detect_backend()
#' @return Temp table name (with # prefix for SQL Server)
#' @keywords internal
temp_table_name = function(base_name, backend) {
  if (backend == "sqlserver") {
    paste0("#", base_name)
  } else {
    base_name
  }
}

#' Create a temp table from a SELECT query
#' 
#' Handles SQL dialect differences:
#' - SQL Server: SELECT ... INTO #table FROM (SELECT ...) AS __subq
#' - Others: CREATE TEMPORARY TABLE name AS SELECT ...
#' 
#' @param conn Database connection
#' @param table_name Name for the temp table
#' @param select_sql The SELECT statement (without CREATE TABLE)
#' @param backend Backend name from detect_backend()
#' @keywords internal
create_temp_table_as = function(conn, table_name, select_sql, backend) {
  if (backend == "sqlserver") {
    # SQL Server: SELECT ... INTO #table FROM ...
    sql = sub("^SELECT", paste0("SELECT * INTO ", table_name, " FROM (SELECT"), 
               select_sql, ignore.case = TRUE)
    sql = paste0(sql, ") AS __subq")
    dbExecute(conn, sql)
  } else {
    # Standard SQL: CREATE TEMPORARY TABLE name AS SELECT ...
    sql = glue("CREATE TEMPORARY TABLE {table_name} AS {select_sql}")
    dbExecute(conn, sql)
  }
}

#' Drop a table if it exists
#'
#' Handles SQL dialect differences for temporary-table cleanup.
#'
#' @param conn Database connection
#' @param table_name Name of the table to drop
#' @param backend Backend name from detect_backend()
#' @keywords internal
drop_table_if_exists = function(conn, table_name, backend) {
  if (backend == "sqlserver") {
    sql = glue("IF OBJECT_ID('tempdb..{table_name}') IS NOT NULL DROP TABLE {table_name}")
  } else {
    sql = glue("DROP TABLE IF EXISTS {table_name}")
  }
  tryCatch(dbExecute(conn, sql), error = function(e) NULL)
}

#' Validate a backend for Databricks AP
#' @keywords internal
validate_databricks_ap_backend = function(backend_info) {
  if (!isTRUE(backend_info[["is_databricks"]]) ||
      !isTRUE(backend_info[["supports_physical_temp_tables"]])) {
    stop(
      "[dbreg] Alternating-projection demeaning on Spark requires Databricks ",
      "session-scoped temporary tables (Databricks SQL or Runtime 18.1+).",
      call. = FALSE
    )
  }
  if (identical(backend_info[["native_query"]], FALSE)) {
    stop(
      "[dbreg] Databricks temporary tables require native query mode. ",
      "Reconnect with `useNativeQuery = TRUE`.",
      call. = FALSE
    )
  }
  invisible(backend_info)
}

#' Materialize an AP state table
#'
#' DuckDB support is intentional: it executes the production window AP path in
#' local tests while normal DuckDB regressions retain the generic AP algorithm.
#'
#' @keywords internal
create_ap_table_as = function(conn, table_name, select_sql, backend_info) {
  backend = backend_info[["name"]]
  if (isTRUE(backend_info[["is_databricks"]])) {
    sql = glue("CREATE OR REPLACE TEMP TABLE {table_name} AS {select_sql}")
  } else if (backend == "duckdb") {
    sql = glue("CREATE OR REPLACE TEMP TABLE {table_name} AS {select_sql}")
  } else {
    validate_databricks_ap_backend(backend_info)
  }
  tryCatch(
    dbExecute(conn, sql),
    error = function(e) {
      requirement = if (isTRUE(backend_info[["is_databricks"]])) {
        " Databricks SQL or Runtime 18.1+ is required."
      } else {
        ""
      }
      stop(
        "[dbreg] Failed to materialize a session-scoped AP table.",
        requirement,
        " Database error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
}

#' Drop a materialized AP state table
#' @keywords internal
drop_ap_table_if_exists = function(conn, table_name, backend_info) {
  backend = backend_info[["name"]]
  sql = if (isTRUE(backend_info[["is_databricks"]])) {
    glue("DROP TEMP TABLE IF EXISTS {table_name}")
  } else if (backend == "duckdb") {
    glue("DROP TABLE IF EXISTS {table_name}")
  } else {
    return(invisible(FALSE))
  }
  tryCatch(
    {
      dbExecute(conn, sql)
      invisible(TRUE)
    },
    error = function(e) invisible(FALSE)
  )
}

#' Generate unique pairs of variables (preserves original nested loop order)
#' @keywords internal
gen_xvar_pairs = function(xvars) {
  pairs = list()
  if (length(xvars) > 1) {
    for (i in seq_along(xvars)) {
      if (i == 1) {
        next
      }
      for (j in seq_len(i - 1)) {
        pairs = c(pairs, list(c(xvars[i], xvars[j])))
      }
    }
  }
  pairs
}
