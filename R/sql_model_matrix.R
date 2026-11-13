#' Construct SQL expressions for a design matrix
#'
#' Expands formula terms into SQL SELECT expressions, handling factor one-hot
#' encoding and interaction terms. Only the right-hand side of the formula is
#' processed; explicit `factor(column)`, `I(column^2)` and `I(column^3)` terms
#' and their interactions are supported. Other transformations are rejected.
#' The response is not expanded, but its source name is reserved
#' beside generated columns.
#'
#' Terms can also use `fixest`-style `i()` syntax: `i(f)` expands `f` like a
#' factor (numeric `f` included), and `i(f, x)` interacts each level of `f`
#' with a numeric or logical `x`. `ref` drops the given level(s) and `keep`
#' keeps only the given levels. A non-symbol second argument is a positional
#' reference selector. `ref = TRUE` drops the first level; `ref = NULL` retains
#' the full indicator basis before collinearity removal. Character selectors
#' support unambiguous case-insensitive partial matches, `@regex` and `@!regex`.
#' Numeric selectors match exact typed values. Columns are named as in `fixest`, e.g.
#' `period::5:treat`. One difference: in a model without an intercept, a plain
#' `i(f)` drops its first level if another factor comes first, whereas
#' `fixest` drops its last level (via collinearity) instead.
#'
#' @param formula A formula (or Formula) object. Only RHS terms are expanded;
#'   the LHS is reserved against generated target collisions.
#' @param conn Database connection
#' @param table Table name or FROM clause
#' @param expand Character: `"all"` expands factors and interactions, 
#'   `"interactions"` only expands interaction terms (factors in main effects 
#'   kept as-is for grouping)
#' @param sep Character separator for interaction term names. Default is
#'   `"_x_"`. Use `":"` for standard R naming convention.
#' @param fe_vars Character vector of fixed effect variable names. These are
#'   treated as "in the model" for determining whether to drop reference levels
#'   in interactions.
#' @param intercept Logical. Whether the model has an intercept. As in
#'   [stats::model.matrix()], without an intercept the first factor main effect
#'   keeps all of its levels (no reference level is dropped).
#' @param reserved_names Additional source-role names projected beside the design
#'   (for example, weights or cluster columns). Generated targets cannot use them.
#' @param factor_refs Optional named reference levels for ordinary categorical
#'   regressors. Numeric-coded columns named here expand categorically; i()
#'   references remain controlled by their own ref/keep selectors.
#' @section Typed categories and resource limits:
#' Integer keys are read as exact decimal text; DOUBLE, logical, character and
#' Date keys retain their native identities. Single-precision floating and
#' fractional DECIMAL categories are rejected. Unsupported native metadata
#' fails before expansion. Set `options(dbreg.i_max_columns = 1000L)` to bound
#' emitted i() columns and level discovery (at most the bound plus one level
#' is retained). The bound applies across IV blocks. Fitted category coding
#' and evaluated selectors are retained for prediction.
#'
#' @return List with:
#'   - `select_exprs`: character vector of SQL expressions
#'   - `col_names`: corresponding column names
#'   - `aliases`: SQL-safe versions of `col_names` (identical unless, e.g.,
#'     factor levels contain spaces or other special characters)
#'   - `level_aliases`: data frame mapping each factor dummy `alias` to its
#'     `var` and `level`
#'   - `factor_levels`: list of factor levels by variable (for reference)
#'   - `columns`: ordered producer metadata: unique IDs, SQL targets, display
#'     labels and construction recipes. SQL targets are not coefficient identities.
#'   - `i_columns`: data frame mapping each `i()` column `alias` to its factor
#'     `var`, `level` and interacted variable `by` (`NA` if none)
#' @importFrom stats terms
#' @export
#' @examples
#' library(DBI)
#' library(duckdb)
#' con = dbConnect(duckdb())
#' duckdb_register(con, "test", data.frame(x1 = 1:3, x2 = c("a", "b", "c")))
#' sql_model_matrix(~ x1 + x2, con, "test")
#' sql_model_matrix(~ x1:x2, con, "test")
#' sql_model_matrix(~ x1:x2, con, "test", sep = ":")
#' sql_model_matrix(~ i(x2, x1, ref = "a"), con, "test")
#' dbDisconnect(con)
sql_model_matrix = function(
  formula,
  conn,
  table,
  expand = c("all", "interactions"),
  sep = "_x_",
  fe_vars = character(),
  intercept = TRUE,
  reserved_names = character(),
  factor_refs = NULL
) {
  expand = match.arg(expand)
  formula_vars = model_vars(formula)
  
  # Extract terms from formula (handles Formula objects too)
  if (inherits(formula, "Formula")) {
    formula = formula(formula, lhs = 0, rhs = 1)
  }
  tt = terms(formula)
  term_labels = attr(tt, "term.labels")
  check_simple_terms(term_labels, "regressors")
  factor_refs = validate_factor_refs(factor_refs, unique(unlist(lapply(term_labels, function(term) model_vars(str2lang(term))))))
  
  if (!length(term_labels)) {
    return(list(
      select_exprs = character(), col_names = character(), aliases = character(),
      factor_levels = list(), level_aliases = data.frame(alias = character(), var = character(), level = character()),
      i_columns = empty_i_columns(), columns = empty_linear_columns()
    ))
  }
  env = environment(formula)
  i_specs = lapply(term_labels, parse_i_term, env = if (is.null(env)) baseenv() else env)

  # Get column types and factor levels from database (an i() factor is
  # expanded by level whatever its type)
  term_vars = Map(function(t, spec) if (is.null(spec)) model_vars(str2lang(t)) else c(spec[["f"]], spec[["var"]]), term_labels, i_specs)
  i_factors = unlist(lapply(i_specs, `[[`, "f"))
  nodes = linear_term_nodes(term_labels[vapply(i_specs, is.null, logical(1))])
  factor_nodes = Filter(function(node) node[["kind"]] == "factor", nodes)
  factor_vars = unique(vapply(factor_nodes, `[[`, character(1), "variable"))
  col_info = get_column_info(conn, table, unique(unlist(term_vars)), level_vars = i_factors, i_declarations = i_specs, intercept = intercept, factor_refs = factor_refs, factor_vars = factor_vars)

  source_info = col_info
  col_info[["term_nodes"]] = nodes
  col_info[["factor_sources"]] = list()
  for (term in names(nodes)) {
    node = nodes[[term]]
    source = node[["variable"]]
    if (node[["kind"]] == "power") next
    col_info[["types"]][[term]] = if (node[["kind"]] == "factor") "factor" else source_info[["types"]][[source]]
    if (!is.null(source_info[["levels"]][[source]])) col_info[["levels"]][[term]] = source_info[["levels"]][[source]]
    if (!is.null(source_info[["key_types"]][[source]])) col_info[["key_types"]][[term]] = source_info[["key_types"]][[source]]
    if (node[["kind"]] == "factor") col_info[["factor_sources"]][[term]] = source
    if (node[["kind"]] == "factor" || source %in% names(factor_refs)) {
      col_info[["forced"]] = unique(c(col_info[["forced"]], term))
    }
  }
  active_levels = intersect(names(col_info[["levels"]]), unique(c(names(nodes), i_factors)))
  col_info[["levels"]] = col_info[["levels"]][active_levels]
  col_info[["key_types"]] = col_info[["key_types"]][intersect(names(col_info[["key_types"]]), active_levels)]
  power_nodes = names(Filter(function(node) node[["kind"]] == "power", nodes))
  col_info[["power_aliases"]] = if (length(power_nodes)) stats::setNames(safe_identifiers(
    paste0("dbreg_power_", seq_along(power_nodes), "__"),
    taken = unique(c(formula_vars, reserved_names, fe_vars)), fold_case = folds_identifier_case(conn)), power_nodes) else character()

  # SQL-safe aliases for factor dummies, one per (variable, level), computed
  # once so that every term uses the same alias. They avoid the model's
  # variable names and each other, and follow the backend's case rules.
  fold_case = folds_identifier_case(conn)
  source_names = unique(c(names(col_info[["types"]]), formula_vars, fe_vars, reserved_names))
  taken = c(source_names, col_info[["power_aliases"]])
  level_aliases = list()
  for (v in names(col_info[["levels"]])) {
    labels = as.character(col_info[["levels"]][[v]])
    if (identical(col_info[["key_types"]][[v]], "numeric") && anyDuplicated(labels)) {
      labels = sprintf("%.17g", col_info[["levels"]][[v]])
    }
    lvls = c(labels, "_constant")
    a = safe_identifiers(paste0(v, lvls), taken = taken, fold_case = fold_case)
    taken = c(taken, a)
    level_aliases[[v]] = data.frame(alias = a, var = v, level = lvls)
  }
  level_aliases = do.call(rbind, c(list(data.frame(alias = character(), var = character(), level = character())), level_aliases))
  col_info[["level_aliases"]] = level_aliases

  result = list(select_exprs = character(), col_names = character(), aliases = character(), columns = empty_linear_columns())
  
  # Without an intercept, the first factor main effect (incl. a plain i(f))
  # keeps all its levels
  full_factor = if (!intercept) {
    is_factor_main = vapply(seq_along(term_labels), function(k) {
      t = term_labels[k]
      spec = i_specs[[k]]
      if (!is.null(spec)) return(is.null(spec[["var"]]) && is.null(spec[["ref"]]))
      length(split_interaction_term(t)) == 1L && identical(col_info[["types"]][[linear_term_key(t)]], "factor")
    }, logical(1))
    term_labels[is_factor_main][1]
  }
  i_columns = empty_i_columns()
  i_positions = integer()
  level_rows = integer()
  for (k in seq_along(term_labels)) {
    term = term_labels[k]
    if (!is.null(i_specs[[k]])) {
      expanded = expand_i_term(i_specs[[k]], col_info, k, keep_all_levels = identical(term, full_factor), conn = conn)
      i_positions = c(i_positions, length(result[["aliases"]])+seq_along(expanded[["alias"]]))
      i_columns = rbind(i_columns, expanded[["columns"]])
    } else {
      expanded = expand_term(
        term, col_info, expand, sep,
        all_terms = unique(c(term_labels, vapply(fe_vars, function(name)
          paste(deparse(as.name(name), backtick = TRUE), collapse = ""), character(1)))),
        fe_vars = fe_vars,
        keep_all_levels = identical(term, full_factor), conn = conn
      )
    }
    count = length(expanded[["sql"]])
    level_rows = c(level_rows, if (is.null(expanded[["level_row"]])) rep(NA_integer_, count) else expanded[["level_row"]])
    if (count) {
      columns = data.frame(id = paste0("t", k, "b", seq_len(count)), term = k,
        basis = seq_len(count), alias = expanded[["alias"]], label = expanded[["name"]], role = "slope")
      columns[["recipe"]] = expanded[["recipe"]]
      result[["columns"]] = rbind(result[["columns"]], columns)
    }
    result[["select_exprs"]] = c(result[["select_exprs"]], expanded[["sql"]])
    result[["col_names"]] = c(result[["col_names"]], expanded[["name"]])
    result[["aliases"]] = c(result[["aliases"]], expanded[["alias"]])
  }
  
  generated = vapply(result[["columns"]][["recipe"]], function(recipe)
    !recipe[["kind"]] %in% c("raw", "numeric"), logical(1))
  result[["aliases"]][!generated] = safe_identifiers(vapply(result[["columns"]][["recipe"]][!generated], `[[`, character(1), "variable"), fold_case = fold_case)
  result[["aliases"]] = allocate_sql_aliases(result[["aliases"]], generated, source_names, fold_case)
  result[["columns"]][["alias"]] = result[["aliases"]]
  result[["columns"]][["label"]] = vapply(result[["columns"]][["recipe"]], linear_recipe_label, character(1))
  i_columns[["alias"]] = result[["aliases"]][i_positions]
  present = which(!is.na(level_rows))
  level_aliases[["alias"]][level_rows[present]] = result[["aliases"]][present]
  result[["level_rows"]] = level_rows
  result[["factor_levels"]] = col_info[["levels"]]
  result[["level_aliases"]] = level_aliases
  result[["i_columns"]] = i_columns
  result
}

#' Get column types and factor levels from database
#' @keywords internal
get_column_info = function(conn, table, vars, level_vars = character(), i_declarations = list(), intercept = TRUE, factor_refs = NULL, factor_vars = character()) {
  # Build FROM clause
  from_sql = if (grepl("^FROM", table, ignore.case = TRUE)) table else paste("FROM", table)
  
  # Query to get types and sample values
  # Fetch zero rows to get the schema, then separate queries for factor levels
  # (sql_limit() handles dialects without LIMIT, e.g. SQL Server's TOP)
  backend = detect_backend(conn)[["name"]]
  schema_query = sql_limit(paste("SELECT", paste(sql_column_name(vars, conn), collapse = ","), from_sql), 0, backend)
  declarations = Filter(Negate(is.null), i_declarations)
  dictionary_vars = unique(vapply(declarations, `[[`, character(1), "f"))
  category_vars = unique(c(names(factor_refs), factor_vars))
  native_schema = NULL
  if ((length(declarations) || length(category_vars)) && backend %in% c("sqlserver", "mysql")) {
    native_schema = read_i_schema(conn, schema_query)
    schema = native_schema[["schema"]]
    if (backend == "sqlserver" && native_schema[["driver"]] == "odbc") {
      columns = native_schema[["columns"]]
      if (is.data.frame(columns) && all(c("name", "type") %in% names(columns))) {
        bits = columns[["name"]][!is.na(columns[["type"]]) & columns[["type"]] == "-7"]
        for (v in intersect(vars, bits)) schema[[v]] = logical()
      }
    }
  } else schema = dbGetQuery(conn, schema_query)
  
  types = list()
  levels = list()
  classes = list()
  key_types = list()
  
  for (v in vars) {
    col_class = class(schema[[v]])[1]
    if (v %in% level_vars && inherits(schema[[v]], "POSIXt")) {
      stop("Date-time variables cannot be i() factors: ", v, ". Convert to a date or character first.", call. = FALSE)
    }
    classes[[v]] = col_class
    # Treat character, factor, logical as categorical
    is_cat = col_class %in% c("character", "factor", "logical") || v %in% names(factor_refs)
    types[[v]] = if (is_cat) "factor" else "numeric"
    if ((is_cat || v %in% level_vars || v %in% factor_vars) && (!v %in% dictionary_vars || v %in% category_vars)) {
      # All distinct levels (a sample could miss rare levels)
      column = sql_column_name(v, conn)
      lvl_query = glue("SELECT DISTINCT {column} FROM (SELECT {column} {from_sql}) lvl WHERE {column} IS NOT NULL ORDER BY {column}")
      if (v %in% category_vars) {
        dictionary = read_category_dictionary(conn, from_sql, v, native_schema)
        keys = dictionary[["values"]]
        type = dictionary[["type"]]
        order = seq_along(keys)
        if (v %in% names(factor_refs)) {
          reference = factor_reference_index(factor_refs[[v]], keys, type)
          order = c(reference, setdiff(order, reference))
        }
        levels[[v]] = keys[order]
        key_types[[v]] = type
      } else {
        levels[[v]] = read_keyed_query(conn, lvl_query, v)[[1]]
        key_types[[v]] = if (identical(col_class, "logical")) "logical" else "character"
      }
    }
  }
  
  dictionaries = i_dictionaries(conn, from_sql, declarations, classes, intercept, native_schema)
  for (dictionary in dictionaries) {
    variable = dictionary[["first"]]
    if (is.null(levels[[variable]])) levels[[variable]] = dictionary[["values"]]
    if (is.null(key_types[[variable]])) key_types[[variable]] = dictionary[["type"]]
  }
  logical_sql = NULL
  logical_columns = names(classes)[vapply(classes, identical, logical(1), "logical")]
  if (backend == "sqlserver" && length(logical_columns)) {
    logical_sql = stats::setNames(paste0(as.character(DBI::dbQuoteIdentifier(conn, logical_columns)), " = 1"), logical_columns)
  }
  list(types = types, levels = levels, classes = classes, i_specs = dictionaries, logical_sql = logical_sql, key_types = key_types, forced = names(factor_refs))
}

#' Expand a single term into SQL expressions
#'
#' Expansions are lists of parallel character vectors: `sql` (expressions),
#' `name` (column names) and `alias` (SQL-safe column names).
#' @keywords internal
expand_term = function(term, col_info, expand, sep = "_x_", all_terms = character(), fe_vars = character(), keep_all_levels = FALSE, conn) {
  vars = vapply(split_interaction_term(term), linear_term_key, character(1))
  is_interaction = length(vars) > 1
  
  if (!is_interaction && expand == "interactions" && !vars[[1L]] %in% col_info[["forced"]] && !identical(col_info[["term_nodes"]][[vars[[1L]]]][["kind"]], "power")) {
    # Main effect, only expanding interactions: return as-is
    var = col_info[["term_nodes"]][[vars[[1L]]]][["variable"]]
    if (is.null(var)) var = vars[[1L]]
    return(list(sql = sql_column_name(var, conn), name = var, alias = safe_identifiers(var, fold_case = folds_identifier_case(conn)),
      recipe = list(list(kind = "raw", variable = var))))
  }
  
  canonical_terms = vapply(all_terms, function(value) paste(sort(vapply(split_interaction_term(value), linear_term_key, character(1))), collapse = ":"), character(1))
  # Expand each variable
  expansions = lapply(vars, function(v) {
    drop_ref = !keep_all_levels
    if (is_interaction && identical(col_info[["types"]][[v]], "factor")) {
      other_term = paste(sort(setdiff(vars, v)), collapse = ":")
      # A fixed effect absorbs its own intercepts, not slopes on its levels: drop
      # the reference only when the rest of the interaction enters on its own.
      drop_ref = other_term %in% canonical_terms || other_term %in% fe_vars
    }
    expand_variable(v, col_info, expand, is_interaction, drop_ref, conn = conn)
  })
  
  # Cross product for interactions
  Reduce(function(a, b) cross_product(a, b, sep), expansions, right = TRUE)
}

#' Expand a single variable into SQL expression(s)
#' @keywords internal
expand_variable = function(var, col_info, expand, in_interaction, drop_ref = TRUE, conn) {
  node = col_info[["term_nodes"]][[var]]
  if (!is.null(node) && node[["kind"]] == "power") {
    source = node[["variable"]]
    if (!col_info[["classes"]][[source]] %in% c("numeric", "integer", "integer64", "logical")) {
      stop("Polynomial sources must be numeric or logical: ", source, call. = FALSE)
    }
    value = sql_cast_double(sql_column_name(source, conn), conn)
    return(list(sql = paste0("POWER(", value, ",", node[["exponent"]], ")"), name = var,
      alias = unname(col_info[["power_aliases"]][var]),
      recipe = list(list(kind = "power", variable = source, exponent = node[["exponent"]]))))
  }
  source = col_info[["factor_sources"]][[var]]
  forced_factor = !is.null(source)
  if (!forced_factor) source = if (is.null(node)) var else node[["variable"]]
  is_factor = identical(col_info[["types"]][[var]], "factor")
  
  # Only expand factors if: expand="all", OR variable is part of interaction
  if (!(is_factor && (expand == "all" || in_interaction || var %in% col_info[["forced"]]))) {
    # Numeric products must be floating point before multiplication.
    sql = if (!is_factor) sql_cast_double(sql_column_name(source, conn), conn) else sql_column_name(source, conn)
    return(list(sql = sql, name = source, alias = safe_identifiers(source, fold_case = folds_identifier_case(conn)),
      recipe = list(list(kind = if (is_factor) "raw" else "numeric", variable = source))))
  }
  lvls = col_info[["levels"]][[var]]
  type = col_info[["key_types"]][[var]]
  la = col_info[["level_aliases"]]
  rows = which(la[["var"]] == var)
  if (length(lvls) < 2) {
    return(list(sql = "1", name = paste0(var, "_constant"), alias = la[["alias"]][utils::tail(rows, 1L)],
      recipe = list(list(kind = "constant", variable = source, levels = lvls, key_type = type, factor = forced_factor)), level_row = utils::tail(rows, 1L)))
  }
  all_levels = lvls
  selected = if (drop_ref) seq_along(lvls)[-1L] else seq_along(lvls)
  lvls = lvls[selected]
  labels = la[["level"]][rows[selected]]
  list(
    sql = paste0("CASE WHEN ", sql_column_name(source, conn), " = ", category_sql_literals(conn, lvls, type), " THEN 1.0 ELSE 0.0 END"),
    name = paste0(var, labels),
    alias = la[["alias"]][rows[selected]],
    level_row = rows[selected],
    recipe = lapply(seq_along(lvls), function(k) list(kind = "indicator", variable = source, value = lvls[k], levels = all_levels, key_type = type, label = labels[k], factor = forced_factor))
  )
}

#' Cartesian product of two term expansions
#'
#' Columns of `a` vary slowest (as in a nested loop over `a`, then `b`).
#' @keywords internal
cross_product = function(a, b, sep = "_x_") {
  i = rep(seq_along(a[["sql"]]), each = length(b[["sql"]]))
  j = rep(seq_along(b[["sql"]]), times = length(a[["sql"]]))
  list(
    sql = paste0("(", a[["sql"]][i], ") * (", b[["sql"]][j], ")"),
    name = paste0(a[["name"]][i], sep, b[["name"]][j]),
    alias = paste0(a[["alias"]][i], "_x_", b[["alias"]][j]),
    recipe = Map(function(left, right) list(kind = "product", left = left, right = right), a[["recipe"]][i], b[["recipe"]][j])
  )
}

empty_i_columns = function() {
  data.frame(alias = character(), var = character(), level = character(), by = character())
}

# Argument template for i() calls (as in fixest::i())
i_template = function(factor_var, var, ref, keep, bin, ref2, keep2, bin2) NULL

#' Expand an `i()` term into SQL expressions
#'
#' Without `var` or `ref`, `i(f)` drops the first level, like a factor main
#' effect (unless `keep_all_levels`). Columns get opaque aliases, since level
#' text is arbitrary.
#' @param spec Result of [parse_i_term()].
#' @param k Term index (to make aliases unique).
#' @keywords internal
expand_i_term = function(spec, col_info, k, keep_all_levels = FALSE, conn) {
  f = spec[["f"]]
  var = spec[["var"]]
  dictionaries = col_info[["i_specs"]]
  matches = vapply(dictionaries, function(item) identical(item[["term"]], spec[["term"]]), logical(1))
  if (!any(matches)) stop("Missing fitted i() dictionary.", call. = FALSE)
  dictionary = dictionaries[[which(matches)[[1L]]]]
  all_levels = dictionary[["values"]]
  selected = dictionary[["selected"]]
  # The first plain no-intercept factor gets a full basis. Other factors retain
  # the previous compiler's treatment-reference convention.
  if (is.null(var) && !isTRUE(spec[["explicit"]]) && !keep_all_levels &&
      identical(selected, seq_along(all_levels))) selected = setdiff(selected, 1L)
  lvls = all_levels[selected]
  lvl_chr = dictionary[["labels"]][selected]
  if (!length(lvls)) {
    return(list(sql = character(), name = character(), alias = character(), columns = empty_i_columns(), recipe = list()))
  }
  literal = category_sql_literals(conn, lvls, dictionary[["type"]])
  cond = paste0(sql_column_name(f, conn), " = ", literal)
  if (is.null(var)) {
    sql = paste0("CASE WHEN ", cond, " THEN 1.0 ELSE 0.0 END")
  } else if (identical(col_info[["classes"]][[var]], "logical")) {
    logical = col_info[["logical_sql"]][[var]]
    if (is.null(logical)) logical = sql_column_name(var, conn)
    sql = paste0("CASE WHEN ", cond, " AND (", logical, ") THEN 1.0 ELSE 0.0 END")
  } else if (col_info[["classes"]][[var]] %in% c("numeric", "integer", "integer64")) {
    sql = paste0("CASE WHEN ", cond, " THEN ", sql_cast_double(sql_column_name(var, conn), conn), " ELSE 0.0 END")
  } else {
    stop("The second argument of i() must be numeric or logical, but '", var, "' is not.", call. = FALSE)
  }
  alias = sprintf("dbreg_i%d_%d", k, seq_along(lvls))
  name = paste0(f, "::", lvl_chr, if (!is.null(var)) paste0(":", var))
  list(
    sql = sql,
    name = name,
    alias = alias,
    columns = data.frame(alias = alias, var = f, level = lvl_chr, by = if (is.null(var)) NA_character_ else var),
    recipe = lapply(seq_along(lvls), function(k) list(kind = "i", variable = f, value = lvls[k], levels = all_levels, key_type = dictionary[["type"]], label = lvl_chr[k], second = var))
  )
}
