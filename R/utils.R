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

#' Format variance-covariance type labels for printing
#' @keywords internal
format_vcov_label = function(vcov_type, n_clusters = NULL) {
  switch(vcov_type,
    "iid" = "IID",
    "hc1" = "Heteroskedasticity-robust",
    "cluster" = if (length(n_clusters) > 1L) paste0("Clustered (", paste(n_clusters, collapse = " and "), " clusters)") else if (!is.null(n_clusters)) sprintf("Clustered (%d clusters)", n_clusters) else "Clustered",
    vcov_type
  )
}

#' Generate coefficient table from estimates and vcov matrix
#' @keywords internal
gen_coeftable = function(betahat, vcov_mat, df_residual) {
  coefs = as.numeric(betahat)
  names(coefs) = standardize_coef_names(rownames(betahat))
  ses = sqrt(Matrix::diag(vcov_mat))
  tstats = coefs / ses
  df = attr(vcov_mat, "df.t")
  if (is.null(df)) df = df_residual
  pvals = 2 * pt(-abs(tstats), df)
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

  scales = sqrt(pmax(as.numeric(Matrix::diag(XtX)), .Machine[["double.xmin"]]))
  qr_decomp = qr(XtX / outer(scales, scales), tol = tol)
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
  scales = sqrt(pmax(as.numeric(Matrix::diag(XtX)), .Machine[["double.xmin"]]))
  XtX = XtX / outer(scales, scales)
  Xty = Xty / scales
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
  list(betahat = betahat / scales, XtX_inv = XtX_inv / outer(scales, scales))
}


formula_column = function(value, label) {
  if (!inherits(value, "formula")) return(value)
  if (length(value) != 2L || !is.symbol(value[[2L]])) {
    stop("`", label, "` formulas must name one column; calculate other expressions in the source.", call. = FALSE)
  }
  as.character(value[[2L]])
}

#' Set up database connection and data source
#'
#' Handles the precedence logic for data sources (table > data > path),
#' connection inference from tbl_lazy objects, and DuckDB registration.
#'
#' @param conn Database connection, or NULL to create a DuckDB connection

#' @param table Character table name or tbl_lazy object
#' @param data R data.frame to register with DuckDB, or tbl_lazy object
#' @param path Character path to file(s) on disk
#' @param resources Optional resource ledger to record connection ownership.
#' @param duckdb_config Configuration for a newly created DuckDB connection.
#' @param caller Character string identifying the calling function for error messages
#'   and temp table naming (e.g., "dbreg" or "dbbinsreg")
#'
#' @return A list containing:
#'   \item{conn}{The database connection (created or provided)}
#'   \item{data,table}{Normalized data sources, with lazy data moved to table}
#'   \item{table_name}{The table name to query from}
#'   \item{from_statement}{SQL FROM clause (for dbreg-style queries)}
#'   \item{own_conn}{Logical, TRUE if we created the connection}
#'   \item{registered_table}{Name of registered table if data.frame was provided, else NULL}
#'   \item{is_duckdb}{Logical, TRUE if the connection is DuckDB}
#'
#' @keywords internal
setup_db_connection = function(conn, table, data, path, caller = "dbreg",
                               resources = NULL, duckdb_config = list()) {
  if (inherits(data, "tbl_lazy") && is.null(table)) table = data
  if (!is.null(table)) data = NULL
  connection_lease = attr(conn, "dbreg_resource_lease", exact = TRUE)
  if (!is.null(connection_lease)) fitted_resource_connection(connection_lease)
  own_conn = FALSE
  registered_table = NULL
  is_duckdb = FALSE
  table_name = NULL
  from_statement = NULL
  setup_complete = FALSE
  on.exit(if (!setup_complete) release_db_setup(conn, own_conn, registered_table), add = TRUE)


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
    integer64_data = is.null(table) && is.data.frame(data) &&
      any(vapply(data, inherits, logical(1), "integer64"))
    conn = dbConnect(duckdb(
      config = duckdb_config,
      bigint = if (integer64_data) "integer64" else "numeric"
    ), shutdown = TRUE)
    own_conn = TRUE
    is_duckdb = TRUE
  } else {
    # Check if user-provided connection is DuckDB
    backend_info = detect_backend(conn)
    is_duckdb = (backend_info[["name"]] == "duckdb")
  }

  flush_retired_resources(conn)

  if (!is.null(resources)) {
    resources[["conn"]] = conn
    resources[["own_conn"]] = own_conn
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
      stop(
        "In-memory data frames are only supported with DuckDB connections. ",
        "Use `table` or `path` for other backends."
      )
    }
    temp_name = sprintf(
      "__%s_%s",
      caller,
      gsub("[^0-9]", "", format(Sys.time(), "%Y%m%d_%H%M%S_%OS3"))
    )
    if (!is.null(resources)) {
      temp_name = basename(tempfile(paste0("__", caller, "_")))
      resources[["registered_table"]] = temp_name
    }
    duckdb_register(conn, temp_name, data)
    registered_table = temp_name
    integer64_columns = names(data)[vapply(data, inherits, logical(1), "integer64")]
    if (length(integer64_columns)) {
      # Some DuckDB versions interpret integer64 storage as raw doubles when
      # the supplied connection uses bigint="numeric". Check SQL types before
      # fitting so keys cannot silently change or collapse above 2^53.
      schema = dbGetQuery(conn, paste("DESCRIBE", DBI::dbQuoteIdentifier(conn, temp_name)))
      types = schema[["column_type"]][match(integer64_columns, schema[["column_name"]])]
      if (anyNA(types) || any(types != "BIGINT")) {
        stop("DuckDB did not preserve integer64 input columns. Use a connection created with ",
          "duckdb::duckdb(bigint = 'integer64'), or convert those columns explicitly before fitting.",
          call. = FALSE
        )
      }
    }
    table_name = temp_name
    from_statement = paste("FROM", temp_name)
  } else if (!is.null(path)) {
    if (!is.character(path)) {
      stop("`path` must be character.")
    }
    # Handle path: check if it's already a read_* function call or raw path
    if (!(grepl("^read|^scan", path) && grepl("'", path))) {
      path = gsub("\"", "'", path)
      from_statement = glue("FROM '{path}'")
      table_name = glue("'{path}'")
    } else {
      from_statement = glue("FROM {path}")
      table_name = path
    }
  } else {
    stop("Provide one of `table`, `data`, or `path`.")
  }

  setup_complete = TRUE
  list(
    conn = conn,
    table_name = table_name,
    data = data,
    table = table,
    from_statement = from_statement,
    own_conn = own_conn,
    registered_table = registered_table,
    is_duckdb = is_duckdb
  )
}


#' Variables used by a formula
#'
#' Like [all.vars()], but only the variable arguments of `i()` terms count
#' (not, e.g., a `ref` value held in a variable).
#' @param x Formula or expression.
#' @param skip_factors Exclude variables used only inside explicit factor() calls.
#' @keywords internal
model_vars = function(x, skip_factors = FALSE) {
  if (inherits(x, "Formula")) x = formula(x)
  walk = function(e) {
    if (is.name(e)) return(as.character(e))
    if (!is.call(e)) return(character())
    if (skip_factors && identical(e[[1L]], as.name("factor"))) return(character())
    args = if (identical(e[[1]], as.name("i"))) {
      as.list(match.call(i_template, e))[c("factor_var", "var")]
    } else {
      as.list(e)[-1]
    }
    unlist(lapply(args, walk))
  }
  unique(walk(x))
}

#' Parse regression formula
#'
#' Parses a formula using the Formula package and extracts components needed
#' for regression: outcome variable, RHS variables, term labels, and fixed effects.
#'
#' @param fml A formula object (will be converted to Formula)
#' @param allow_empty Allow an intercept-only RHS (GLM).
#' @param allow_extended Allow restricted powers and composite fixed effects.
#' @param allow_i Allow restricted additive factor indicators/group slopes.
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
parse_regression_formula = function(fml, allow_empty = FALSE, allow_extended = FALSE, allow_i = FALSE) {
  fml = Formula(fml)
  # Extract outcome variable (LHS)
  lhs = formula(fml, lhs = 1, rhs = 0)
  yvar = all.vars(lhs)
  if (length(yvar) != 1) {
    stop("Exactly one outcome variable required.")
  }
  check_simple_terms(paste(deparse(lhs[[2]], backtick = TRUE), collapse = " "), "outcome", allow_interactions = FALSE)

  # Extract RHS part 1: variables and term structure
  rhs1 = formula(fml, lhs = 0, rhs = 1)
  tt = terms(rhs1)
  term_labels = attr(tt, "term.labels")
  has_intercept = attr(tt, "intercept") == 1L
  xvars = model_vars(rhs1) # unique variable names
  # i() terms need SQL expansion, like interactions
  has_i = vapply(term_labels, function(t) {
    !is.null(parse_sql_i_term(t))
  }, logical(1))
  has_interactions = any(grepl(":", term_labels)) || any(has_i)
  has_transformations = any(grepl("^I\\(", term_labels)) || length(explicit_factor_sources(term_labels)) > 0L

  if (isTRUE(allow_empty) || isTRUE(allow_extended)) {
    validate_regressor_expression = function(expr) {
      if (is.symbol(expr) || is.numeric(expr)) return(invisible(NULL))
      if (!is.call(expr)) {
        stop("Unsupported regressor expression: ", deparse(expr), call. = FALSE)
      }
      if (identical(expr[[1L]], as.name("factor"))) {
        explicit_factor_sources(paste(deparse(expr, width.cutoff = 500L), collapse = ""))
        return(invisible(NULL))
      }
      if (allow_i && !is.null(parse_sql_i_term(expr))) return(invisible(NULL))
      if (!allow_i && !is.null(parse_sql_i_term(paste(deparse(expr), collapse = "")))) return(invisible(NULL))
      call_name = as.character(expr[[1]])
      if (!allow_i && call_name %in% c(":", "*") && contains_i_term(expr)) {
        stop("i() terms can't be interacted with other variables yet.", call. = FALSE)
      }
      if (call_name == "-") {
        args = as.list(expr)[-1]
        last = args[[length(args)]]
        if (is.numeric(last) && length(last) == 1 && identical(as.numeric(last), 1)) {
          if (length(args) == 2) validate_regressor_expression(args[[1]])
          return(invisible(NULL))
        }
      }
      if (call_name == "I") {
        if (length(expr) != 2 || !is.call(expr[[2]]) ||
          !identical(as.character(expr[[2]][[1]]), "^") ||
          !is.symbol(expr[[2]][[2]]) || !is.numeric(expr[[2]][[3]]) ||
          !as.numeric(expr[[2]][[3]]) %in% c(2, 3)) {
          stop(
            "Unsupported regressor expression `", deparse(expr),
            "`. Only `I(column^2)` and `I(column^3)` are supported.",
            call. = FALSE
          )
        }
        return(invisible(NULL))
      }
      if (call_name == "^" || !call_name %in% c("+", ":", "*", "(")) {
        stop(
          "Unsupported regressor expression `", deparse(expr),
          "`. Supported terms are column names, `I(column^2)`, `I(column^3)`, ",
          "and the operators `+`, `:`, and `*`.",
          call. = FALSE
        )
      }
      for (arg in as.list(expr)[-1]) validate_regressor_expression(arg)
      invisible(NULL)
    }
    if (length(rhs1) > 1) validate_regressor_expression(rhs1[[2]])
  } else {
    check_simple_terms(term_labels, "regressors")
  }

  if (!length(xvars) && !isTRUE(allow_empty)) {
    stop("No regressors on RHS.")
  }

  # Extract fixed effects (RHS part 2, after |)
  fe_components = NULL
  fe_labels = NULL
  fe = NULL
  fe_slopes = list()
  if (length(fml)[2] > 1) {
    fe_rhs = formula(fml, lhs = 0, rhs = 2)
    if (isTRUE(allow_empty) || isTRUE(allow_extended)) {
      flatten_call = function(expr, operator) {
        if (is.call(expr) && identical(as.character(expr[[1]]), "(")) {
          return(flatten_call(expr[[2]], operator))
        }
        if (is.call(expr) && identical(as.character(expr[[1]]), operator)) {
          return(c(flatten_call(expr[[2]], operator), flatten_call(expr[[3]], operator)))
        }
        list(expr)
      }
      fe_terms = flatten_call(fe_rhs[[2]], "+")
      term_slopes = lapply(fe_terms, function(term) {
        if (!is.call(term) || !identical(term[[1L]], as.name("["))) return(character())
        args = as.list(term)[-1L]
        if (length(args) < 2L || !all(vapply(args, is.symbol, logical(1)))) {
          stop("Use group[trend] or group[trend1, trend2] with source columns.", call. = FALSE)
        }
        vapply(args[-1L], as.character, character(1))
      })
      fe_components = lapply(fe_terms, function(term) {
        if (is.call(term) && identical(term[[1L]], as.name("["))) term = term[[2L]]
        components = flatten_call(term, "^")
        if (!all(vapply(components, is.symbol, logical(1)))) {
          stop("Fixed-effect interactions must contain column names only.", call. = FALSE)
        }
        vapply(components, as.character, character(1))
      })
      fe_labels = vapply(fe_components, paste, collapse = "^", character(1))
      unique_labels = unique(fe_labels)
      fe_slopes = lapply(unique_labels, function(label) unique(unlist(term_slopes[fe_labels == label], use.names = FALSE)))
      fe_components = fe_components[match(unique_labels, fe_labels)]
      fe_labels = unique_labels
      fe = vapply(seq_along(fe_components), function(i) {
        if (length(fe_components[[i]]) == 1) fe_components[[i]] else paste0("dbreg_fe_", i, "__")
      }, character(1))
    } else {
      check_simple_terms(attr(terms(fe_rhs), "term.labels"), "fixed effects", allow_interactions = FALSE)
      fe = all.vars(fe_rhs)
      fe_labels = fe
      fe_components = as.list(fe)
    }
  }

  list(
    fml = fml,
    yvar = yvar,
    xvars = xvars,
    term_labels = term_labels,
    has_interactions = has_interactions,
    has_intercept = has_intercept,
    has_transformations = has_transformations,
    fe = fe,
    fe_labels = fe_labels,
    fe_components = fe_components,
    fe_slopes = fe_slopes
  )
}


#' Release a connection set up by setup_db_connection()
#'
#' Unregisters any registered in-memory table and disconnects connections
#' that we created.
#' @keywords internal
release_db_setup = function(conn, own_conn, registered_table) {
  # Nothing of ours to release (e.g. a user connection, or setup never got
  # that far). Never throw: this runs in on.exit() and must not mask errors.
  if (!isTRUE(own_conn) && is.null(registered_table)) {
    return(invisible(NULL))
  }
  if (is.null(conn) || !isTRUE(tryCatch(dbIsValid(conn), error = function(e) FALSE))) {
    return(invisible(NULL))
  }
  if (!is.null(registered_table)) {
    try(duckdb_unregister(conn, registered_table), silent = TRUE)
  }
  if (isTRUE(own_conn)) {
    dbDisconnect(conn, shutdown = TRUE)
  }
  invisible(NULL)
}

# R's matrix builders do not decode integer64 storage. Evaluate factor()
# first, preserving exact categorical labels, then convert numeric columns to
# the same double precision used by the SQL regression arithmetic.
model_matrix_data = function(formula, data) {
  if (!any(vapply(data, inherits, logical(1), "integer64"))) return(data)
  data = stats::model.frame(formula, data = data, na.action = stats::na.pass)
  for (name in names(data)) if (inherits(data[[name]], "integer64")) {
    data[[name]] = as.numeric(data[[name]])
  }
  data
}

# ---- SQL-safe identifiers ---------------------------------------------------
# Internally, all generated SQL uses identifiers that are valid unquoted on any
# backend. Model variables with other names (e.g. `Sepal.Length`, names with
# spaces, or reserved words) are renamed once at the source (see
# safe_source()) and mapped back to their display names in the final result
# (see restore_result_names()).

# Common reserved words across backends (lower case)
sql_reserved = c(
  "all", "and", "as", "asc", "between", "by", "case", "cast", "check",
  "column", "create", "cross", "default", "delete", "desc", "distinct",
  "drop", "else", "end", "except", "exists", "false", "for", "foreign", "from",
  "full", "group", "having", "in", "inner", "insert", "intersect", "into", "is",
  "join", "left", "like", "limit", "natural", "not", "null", "offset", "on",
  "or", "order", "outer", "primary", "references", "right", "select", "table",
  "then", "to", "true", "union", "unique", "update", "user", "using", "values",
  "when", "where", "window", "with"
)

#' Does the backend fold the case of identifiers?
#'
#' Athena/Trino stores and returns unquoted identifiers in lower case.
#' On these backends, safe identifiers must be
#' lower case (both to reference columns and to read results back by name).
#' @keywords internal
folds_identifier_case = function(conn) {
  identical(detect_backend(conn)[["name"]], "athena")
}

#' Is `x` usable as an unquoted SQL identifier?
#' @param fold_case Logical. Require lower case (see [folds_identifier_case()]).
#' @keywords internal
is_safe_identifier = function(x, fold_case = FALSE) {
  pattern = if (fold_case) "^[a-z_][a-z0-9_]*$" else "^[A-Za-z_][A-Za-z0-9_]*$"
  grepl(pattern, x) & !(tolower(x) %in% sql_reserved)
}

#' Make SQL-safe identifiers
#'
#' Leaves safe names unchanged; otherwise replaces invalid characters with `_`
#' (lower-casing on case-folding backends), prefixes names that would start
#' with a digit or are reserved, and de-duplicates against each other and
#' against `taken`. De-duplication ignores case, since most backends match
#' identifiers case-insensitively.
#' @inheritParams is_safe_identifier
#' @keywords internal
safe_identifiers = function(x, taken = character(), fold_case = FALSE) {
  keep = is_safe_identifier(x, fold_case) & !tolower(x) %in% tolower(taken) & !duplicated(tolower(x))
  out = gsub("[^A-Za-z0-9_]", "_", x)
  if (fold_case) out = tolower(out)
  bad_start = !grepl("^[A-Za-z_]", out) | tolower(out) %in% sql_reserved
  out[bad_start] = paste0("v_", out[bad_start])
  out[keep] = x[keep]
  gen = which(!keep)
  if (length(gen)) {
    # De-duplicate generated names (case-insensitively) against `taken`, the
    # kept names and each other in one pass; make.unique() leaves the earlier
    # (fixed) entries alone and suffixes later duplicates with _1, _2, ...
    fixed = tolower(c(taken, out[keep]))
    lc = tolower(out[gen])
    uniq = utils::tail(make.unique(c(fixed, lc), sep = "_"), length(gen))
    out[gen] = paste0(out[gen], substring(uniq, nchar(lc) + 1L))
  }
  out
}

# Rename variable symbols in an expression or formula (`map`: display -> safe).
# Function names (the head of a call) are left alone.
rename_expr_vars = function(e, map) {
  if (is.name(e)) {
    nm = as.character(e)
    return(if (nm %in% names(map)) as.name(map[[nm]]) else e)
  }
  if (is.call(e) && length(e) > 1L) {
    for (i in seq_along(e)[-1L]) {
      if (!is.null(e[[i]])) e[[i]] = rename_expr_vars(e[[i]], map)
    }
  }
  e
}

#' Give model variables SQL-safe names at the source
#'
#' Wraps `from_statement` in a subquery that selects `vars`, quoting (via
#' [glue::glue_sql()]) and aliasing any that are not safe identifiers.
#' If every variable is already safe, the source is returned unchanged.
#'
#' @return List with:
#'   - `from_statement`: the (possibly wrapped) FROM clause, and `table`, the
#'     same without the leading `FROM`;
#'   - `map`: display names indexed by their safe names (empty if nothing was
#'     renamed);
#'   - `to_safe()`, `to_display()`: map names (character vectors, `NULL`, or
#'     formulas for `to_safe()`) between the two. Identity if nothing was
#'     renamed.
#' @keywords internal
safe_source = function(conn, from_statement, vars) {
  vars = unique(vars[!is.na(vars) & nzchar(vars)])
  # Cell summaries use these aliases independently of source variables.
  # Reserve them before strategy selection, including automatic compression.
  reserved = c(
    "n", "sum_w", "sum_wy", "sum_wy_sq", "sum_w2", "sum_w2y",
    "sum_w2y_sq", "sum_w_gc", "sum_wy_gc", "mean_Y", "wts",
    "dbreg_cell_center__", "dbreg_cell_scale__"
  )
  safe = safe_identifiers(vars, taken = reserved, fold_case = folds_identifier_case(conn))
  renamed = safe != vars
  map = setNames(vars[renamed], safe[renamed])
  if (any(renamed)) {
    sel = vars
    # Quote the original names for this backend; keep the safe aliases unquoted
    # (case-folding backends need consistent aliases)
    sel[renamed] = glue_sql("{`vars[renamed]`} AS {DBI::SQL(safe[renamed])}", .con = conn)
    from_statement = paste0("FROM (SELECT ", paste(sel, collapse = ", "), " ", from_statement, ") AS dbreg_src")
  }
  lookup = function(x, from, to) {
    if (is.null(x) || !length(from)) return(x)
    i = match(x, from)
    ifelse(is.na(i), x, to[i])
  }
  list(
    from_statement = from_statement,
    table = sub("^FROM\\s+", "", from_statement),
    map = map,
    to_safe = function(x) {
      if (inherits(x, "formula")) {
        if (!length(map)) return(x)
        return(rename_expr_vars(formula(x), setNames(names(map), map)))
      }
      lookup(x, unname(map), names(map))
    },
    to_display = function(x) lookup(x, names(map), unname(map))
  )
}

#' Register display names for expanded factor columns
#'
#' Adds entries to `inputs$name_map` so that factor dummies named after their
#' (safe) variable and level map back to the display variable name plus the
#' original level text (e.g. `grptreated group`). Covers both R-side model
#' matrix names (`paste0(var, level)`) and the SQL-safe aliases generated by
#' [sql_model_matrix()].
#' @param factor_levels Named list of levels by (safe) variable name.
#' @param level_aliases Optional data frame of SQL aliases with their `var`
#'   and `level` (the `level_aliases` element of [sql_model_matrix()]).
#' @param i_columns Optional data frame of `i()` column aliases (the
#'   `i_columns` element of [sql_model_matrix()]). Their display names, e.g.
#'   `period::5:treat`, are also recorded in `inputs$i_columns` for
#'   `predict()`.
#' @param factor_classes Optional named fitted source classes for factors.
#' @param i_specs Typed fitted i() dictionaries shared by R and SQL prediction.
#' @keywords internal
register_factor_display = function(inputs, factor_levels = list(), level_aliases = NULL, col_names = NULL, aliases = NULL, factor_classes = NULL, i_columns = NULL, i_specs = list()) {
  saved = inputs[["factor_levels"]]
  if (is.null(saved)) saved = list()
  for (nm in names(factor_levels)) saved[[nm]] = factor_levels[[nm]]
  inputs[["factor_levels"]] = saved
  classes = inputs[["factor_classes"]]
  if (is.null(classes)) classes = character()
  classes[names(factor_classes)] = factor_classes
  inputs[["factor_classes"]] = classes
  map = inputs[["name_map"]]
  if (is.null(map)) map = setNames(character(), character())
  disp = function(v) vapply(v, function(term) {
    if (length(explicit_factor_sources(term))) restore_term_names(term, inputs[["source_map"]]) else
      coef_var_name(ifelse(term %in% names(map), map[term], term))
  }, character(1), USE.NAMES = FALSE)
  vars = rep(names(factor_levels), lengths(factor_levels))
  lvls = unlist(lapply(factor_levels, as.character), use.names = FALSE)
  factor_display = setNames(paste0(disp(vars), lvls), paste0(vars, lvls))
  col_names = vapply(col_names, function(nm) {
    parts = strsplit(nm, "_x_", fixed = TRUE)[[1]]
    parts = vapply(parts, function(part) {
      if (part %in% names(factor_display)) return(factor_display[[part]])
      if (part %in% names(map)) return(disp(part))
      if (is.null(parse_power_term(part))) return(part)
      paste(deparse(rename_expr_vars(str2lang(part), inputs[["source_map"]]), width.cutoff = 500L), collapse = "")
    }, character(1))
    paste(parts, collapse = "_x_")
  }, character(1))
  internal = c(paste0(vars, lvls), level_aliases[["alias"]], aliases)
  display = c(paste0(disp(vars), lvls), paste0(disp(level_aliases[["var"]]), level_aliases[["level"]]), col_names)
  new = setNames(display, internal)[internal != display]
  new = new[!duplicated(names(new)) & !names(new) %in% names(map)]
  if (NROW(i_columns)) {
    i_columns = i_columns[!duplicated(i_columns[["alias"]]), , drop = FALSE]
    raw = function(v) ifelse(is.na(v), NA_character_, ifelse(v %in% names(map), map[v], v))
    i_display = data.frame(
      name = paste0(raw(i_columns[["var"]]), "::", i_columns[["level"]], ifelse(is.na(i_columns[["by"]]), "", paste0(":", raw(i_columns[["by"]])))),
      var = raw(i_columns[["var"]]),
      level = i_columns[["level"]],
      by = raw(i_columns[["by"]])
    )
    if (anyDuplicated(i_display[["name"]])) stop("Design terms produce duplicate public coefficient names.", call. = FALSE)
    new = new[!names(new) %in% i_columns[["alias"]]]
    new = c(new, setNames(i_display[["name"]], i_columns[["alias"]]))
    i_display = rbind(inputs[["i_columns"]], i_display)
    inputs[["i_columns"]] = i_display[!duplicated(i_display[["name"]]), , drop = FALSE]
  }
  if (length(i_specs)) {
    public = function(value) restore_names(value, inputs[["source_map"]])
    specs = lapply(i_specs, function(spec) {
      spec[["first"]] = public(spec[["first"]])
      spec[["second"]] = public(spec[["second"]])
      spec[["term"]] = restore_term_names(spec[["term"]], inputs[["source_map"]])
      spec[["names"]] = if (!length(spec[["selected"]])) character() else paste0(
        spec[["first"]], "::", spec[["labels"]][spec[["selected"]]],
        if (!is.null(spec[["second"]])) paste0(":", spec[["second"]])
      )
      spec
    })
    specs = c(inputs[["i_specs"]], specs)
    inputs[["i_specs"]] = specs[!duplicated(vapply(specs, `[[`, character(1), "term"))]
  }
  inputs[["name_map"]] = c(map, new)
  invisible(inputs)
}

#' Map internal (safe) names back to display names
#'
#' Maps each name via `map` (display names indexed by safe names), or else
#' splits it on the interaction separators (`:` and `_x_`) and maps each token.
#' @keywords internal
restore_names = function(x, map) {
  if (!length(map) || is.null(x) || !length(x)) {
    return(x)
  }
  out = vapply(x, function(s) {
    if (is.na(s)) return(s)
    if (s %in% names(map)) return(map[[s]])
    m = gregexpr(":|_x_", s)
    tokens = regmatches(s, m, invert = TRUE)[[1]]
    seps = regmatches(s, m)[[1]]
    tokens = ifelse(tokens %in% names(map), map[tokens], tokens)
    paste0(c(rbind(tokens, c(seps, ""))), collapse = "")
  }, character(1), USE.NAMES = FALSE)
  if (!is.null(names(x))) names(out) = names(x)
  out
}

# Variable names as they appear inside coefficient names: like model.matrix(),
# backtick names that aren't syntactic in R (e.g. `my grp`)
coef_var_name = function(v) ifelse(make.names(v) == v, v, paste0("`", v, "`"))

#' Map a result's internal names back to display names
#'
#' Applies `inputs$name_map` (see [safe_source()] and
#' [register_factor_display()]) to coefficient, variable and data names, and
#' restores the user's formula.
#' @keywords internal
restore_result_names = function(result, inputs) {
  map = inputs[["name_map"]]
  if (!length(map)) {
    return(result)
  }
  # Coefficient names backtick non-syntactic variable names (as in lm())
  coef_map = map
  is_var = names(coef_map) %in% names(inputs[["source_map"]])
  coef_map[is_var] = coef_var_name(coef_map[is_var])
  if (is.data.frame(result)) {
    names(result) = restore_names(names(result), map)
    return(result)
  }
  if (!is.null(result[["coeftable"]])) {
    rownames(result[["coeftable"]]) = restore_names(rownames(result[["coeftable"]]), coef_map)
  }
  if (!is.null(result[["vcov"]])) {
    dimnames(result[["vcov"]]) = lapply(dimnames(result[["vcov"]]), restore_names, map = coef_map)
  }
  for (f in c("xvars", "collin.var", "coef_names", "coef_names_all", "endogenous", "instruments")) {
    if (!is.null(result[[f]])) result[[f]] = restore_names(result[[f]], coef_map)
  }
  for (f in c("yvar", "fe", "weights")) {
    if (!is.null(result[[f]])) result[[f]] = restore_names(result[[f]], map)
  }
  if (!is.null(result[["n_fe_levels"]]) && !is.null(names(result[["n_fe_levels"]]))) {
    names(result[["n_fe_levels"]]) = restore_names(names(result[["n_fe_levels"]]), map)
  }
  if (is.data.frame(result[["data"]])) {
    names(result[["data"]]) = restore_names(names(result[["data"]]), map)
  }
  for (d in c("first_stage_f", "first_stage_wald")) {
    if (!is.null(names(result[["diagnostics"]][[d]]))) {
      names(result[["diagnostics"]][[d]]) = restore_names(names(result[["diagnostics"]][[d]]), coef_map)
    }
  }
  # The user's own formula (IV results store the structural formula)
  if (!is.null(result[["fml"]])) {
    fml = if (inherits(result[["fml"]], "Formula")) formula(result[["fml"]]) else result[["fml"]]
    result[["fml"]] = Formula(rename_expr_vars(fml, map))
  }
  result
}

#' Reject formula terms that cannot be expressed in SQL (yet)
#'
#' Terms must be plain variables, interactions of them (`a:b`) or (where
#' interactions are allowed) `i()` terms. Anything else, e.g. `log(x)`, `I(x^2)` or `factor(x)`, would otherwise be silently
#' reduced to its underlying variable by `all.vars()`.
#'
#' @param term_labels Term labels, e.g. from `attr(terms(f), "term.labels")`.
#' @param what Label for the error message.
#' @param allow_interactions Logical. Allow `a:b` terms.
#' @keywords internal
check_simple_terms = function(term_labels, what = "regressors", allow_interactions = TRUE) {
  is_simple = function(term) {
    # fixest-style i() terms are expanded by sql_model_matrix()
    if (allow_interactions && !is.null(parse_sql_i_term(term))) return(TRUE)
    parts = if (allow_interactions) strsplit(term, ":", fixed = TRUE)[[1]] else term
    all(vapply(parts, function(x) {
      isTRUE(tryCatch(is.name(str2lang(x)), error = function(e) FALSE))
    }, logical(1)))
  }
  bad = term_labels[!vapply(term_labels, is_simple, logical(1))]
  if (allow_interactions && any(grepl("(^|:)i\\(", bad))) {
    stop(
      "i() terms can't be interacted with other variables yet: ",
      paste(bad[grepl("(^|:)i\\(", bad)], collapse = ", "),
      call. = FALSE
    )
  }
  if (length(bad)) {
    stop(
      "Transformed or function terms are not supported in the ", what, ": ",
      paste(bad, collapse = ", "), ". ",
      "Please create these variables in the data (or a database view) first.",
      call. = FALSE
    )
  }
  invisible(TRUE)
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

#' Backend metadata parser
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
  if ("BigQueryConnection" %in% connection_class) {
    return(list(
      name = "bigquery", supports_count_big = FALSE, is_databricks = FALSE,
      supports_physical_temp_tables = FALSE, supports_persistent_ap_tables = FALSE, native_query = NA
    ))
  }
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
  if (grepl("sqlite|postgres", dbms) ||
    any(grepl("^(SQLite|PostgreSQL|Pq)Connection$", connection_class, ignore.case = TRUE))) {
    stop("SQLite and PostgreSQL connections are not supported by dbreg. Use DuckDB for local analytics.", call. = FALSE)
  }
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
  } else if (grepl("mysql|mariadb", dbms)) {
    "mysql"
  } else if (is_databricks || grepl("spark", dbms) || grepl("spark", metadata)) {
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
    # OSS Spark uses managed Parquet tables to bound iterative AP plan depth.
    supports_persistent_ap_tables = name == "spark" && !is_databricks,
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
  if (inherits(conn, "duckdb_connection")) {
    return(list(
      name = "duckdb",
      supports_count_big = FALSE,
      is_databricks = FALSE,
      supports_physical_temp_tables = FALSE,
      supports_persistent_ap_tables = FALSE,
      native_query = NA
    ))
  }
  if (inherits(conn, "BigQueryConnection")) {
    return(list(
      name = "bigquery",
      supports_count_big = FALSE,
      is_databricks = FALSE,
      supports_physical_temp_tables = FALSE,
      supports_persistent_ap_tables = FALSE,
      native_query = NA
    ))
  }

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

#' Validate a backend for Spark AP
#' @keywords internal
validate_databricks_ap_backend = function(backend_info) {
  is_databricks = isTRUE(backend_info[["is_databricks"]]) &&
    isTRUE(backend_info[["supports_physical_temp_tables"]])
  is_oss_spark = identical(backend_info[["name"]], "spark") &&
    isTRUE(backend_info[["supports_persistent_ap_tables"]])
  if (!is_databricks && !is_oss_spark) {
    stop(
      "[dbreg] Alternating-projection demeaning on Spark requires either ",
      "Databricks session-scoped temporary tables or OSS Spark managed tables.",
      call. = FALSE
    )
  }
  if (is_databricks && identical(backend_info[["native_query"]], FALSE)) {
    stop(
      "[dbreg] Databricks temporary tables require native query mode. ",
      "Reconnect with `useNativeQuery = TRUE`.",
      call. = FALSE
    )
  }
  invisible(backend_info)
}

#' Cast to double precision
#'
#' `FLOAT` is single precision on DuckDB, MySQL, Spark and Athena (`REAL`),
#' which limits moment accuracy to ~1e-7. Use each backend's 8-byte type.
#' @keywords internal
sql_cast_double = function(sql, conn) {
  type = switch(detect_backend(conn)[["name"]],
    sqlserver = "FLOAT",
    bigquery = "FLOAT64",
    "DOUBLE"
  )
  sprintf("CAST(%s AS %s)", sql, type)
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
    "mysql" = "RAND()",
    "RANDOM()" # default fallback
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
  } else if (backend == "sqlserver") {
    # SQL Server: TABLESAMPLE is block-level, not row-level; use TOP with NEWID()
    glue("SELECT TOP {n} {select_clause} FROM {table_name} {where_sql} ORDER BY NEWID()")
  } else {
    # Other supported backends use their random-order fallback
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
      select_sql,
      ignore.case = TRUE
    )
    sql = paste0(sql, ") AS __subq")
    dbExecute(conn, sql)
  } else if (backend == "spark") {
    # Databricks has session temporary tables. OSS Spark uses owned managed
    # Parquet tables so repeated AP sweeps do not accumulate a lazy view plan.
    sql = if (isTRUE(detect_backend(conn)[["is_databricks"]]))
      glue("CREATE OR REPLACE TEMPORARY TABLE {table_name} AS {select_sql}") else
      glue("CREATE TABLE {table_name} USING PARQUET AS {select_sql}")
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
  } else if (backend == "spark") {
    sql = if (isTRUE(detect_backend(conn)[["is_databricks"]])) glue("DROP TEMP TABLE IF EXISTS {table_name}") else
      glue("DROP TABLE IF EXISTS {table_name}")
  } else {
    sql = glue("DROP TABLE IF EXISTS {table_name}")
  }
  tryCatch(dbExecute(conn, sql), error = function(e) NULL)
}

# Silently ignoring a requested offset or covariance specification changes the
# statistical model. Keep intentional legacy aliases explicit at each caller.
check_unused_args = function(dots, allowed = character()) {
  if (!length(dots)) return(invisible(NULL))
  labels = names(dots)
  if (is.null(labels)) labels = rep("", length(dots))
  labels[is.na(labels)] = ""
  unused = !nzchar(labels) | !labels %in% allowed
  if (any(unused)) {
    stop("Unused arguments: ", paste(ifelse(nzchar(labels[unused]),
      paste0("`", labels[unused], "`"), "<unnamed>"
    ), collapse = ", "), call. = FALSE)
  }
  invisible(NULL)
}

#' Drop a materialized AP state table
#' @keywords internal
drop_ap_table_if_exists = function(conn, table_name, backend_info) {
  backend = if (is.character(backend_info)) backend_info else backend_info[["name"]]
  if (is.list(backend_info) && isTRUE(backend_info[["is_databricks"]])) {
    sql = glue("DROP TEMP TABLE IF EXISTS {table_name}")
    tryCatch(
      {
        dbExecute(conn, sql)
        invisible(TRUE)
      },
      error = function(e) invisible(FALSE)
    )
  } else {
    drop_table_if_exists(conn, table_name, backend)
    invisible(TRUE)
  }
}
