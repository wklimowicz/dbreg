#' Construct SQL expressions for a design matrix
#'
#' Expands formula terms into SQL SELECT expressions, handling factor one-hot
#' encoding and interaction terms. Only the right-hand side of the formula is
#' processed; the response variable (LHS) is ignored.
#'
#' Terms can also use `fixest`-style `i()` syntax: `i(f)` expands `f` like a
#' factor (numeric `f` included), and `i(f, x)` interacts each level of `f`
#' with a numeric or logical `x`. `ref` drops the given level(s) and `keep`
#' keeps only the given levels. Columns are named as in `fixest`, e.g.
#' `period::5:treat`. One difference: in a model without an intercept, a plain
#' `i(f)` drops its first level if another factor comes first, whereas
#' `fixest` drops its last level (via collinearity) instead.
#'
#' @param formula A formula (or Formula) object. Only RHS terms are expanded;
#'   the LHS (response variable) is ignored.
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
#' @param factor_refs Optional named reference levels for categorical regressors.
#' @return List with:
#'   - `select_exprs`: character vector of SQL expressions
#'   - `col_names`: corresponding column names
#'   - `aliases`: SQL-safe versions of `col_names` (identical unless, e.g.,
#'     factor levels contain spaces or other special characters)
#'   - `level_aliases`: data frame mapping each factor dummy `alias` to its
#'     `var` and `level`
#'   - `factor_levels`: list of factor levels by variable (for reference)
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
  factor_refs = NULL
) {
  expand = match.arg(expand)
  
  # Extract terms from formula (handles Formula objects too)
  if (inherits(formula, "Formula")) {
    formula = formula(formula, lhs = 0, rhs = 1)
  }
  tt = terms(formula)
  term_labels = attr(tt, "term.labels")
  
  if (!length(term_labels)) {
    return(list(
      select_exprs = character(), col_names = character(), aliases = character(),
      factor_levels = list(), level_aliases = data.frame(alias = character(), var = character(), level = character()),
      i_columns = empty_i_columns(), coef_names = character(), terms = character(),
      source_classes = character(), i_specs = list(), ordinary_rhs = formula
    ))
  }
  env = environment(formula)
  i_specs = lapply(term_labels, parse_sql_i_term, env = if (is.null(env)) baseenv() else env)

  # Get column types and factor levels from database (an i() factor is
  # expanded by level whatever its type)
  term_vars = Map(function(t, spec) if (is.null(spec)) all.vars(str2lang(t)) else c(spec[["f"]], spec[["var"]]), term_labels, i_specs)
  i_factors = unlist(lapply(i_specs, `[[`, "f"))
  forced = explicit_factor_sources(term_labels)
  col_info = get_column_info(conn, table, unique(unlist(term_vars)), level_vars = i_factors,
    factor_refs = factor_refs, i_declarations = i_specs, intercept = intercept,
    factor_vars = unname(forced))
  col_info = explicit_factor_info(term_labels, col_info)
  # Forced categories do not change the same source's numeric main effect.
  for (source in setdiff(unname(forced), i_factors)) {
    if (!identical(col_info[["types"]][[source]], "factor")) {
      col_info[["levels"]][[source]] = NULL
      col_info[["level_sql"]][[source]] = NULL
      col_info[["factor_classes"]] = col_info[["factor_classes"]][setdiff(names(col_info[["factor_classes"]]), source)]
    }
  }

  result = sql_model_matrix_from_info(term_labels, col_info, expand, sep, fe_vars, intercept,
    fold_case = folds_identifier_case(conn), i_specs = i_specs)
  ordinary = term_labels[vapply(i_specs, is.null, logical(1))]
  result[["ordinary_rhs"]] = stats::reformulate(if (length(ordinary)) ordinary else "1",
    intercept = attr(tt, "intercept") == 1L, env = baseenv())
  result
}

# The fitted-design path supplies saved levels instead of discovering them from
# newdata. Both paths use exactly the same expansion and reference decisions.
sql_model_matrix_from_info = function(term_labels, col_info, expand = "all", sep = "_x_",
                                      fe_vars = character(), intercept = TRUE, fold_case = FALSE, i_specs = lapply(term_labels, parse_sql_i_term, env = baseenv())) {
  col_info = explicit_factor_info(term_labels, col_info)
  # SQL-safe aliases for factor dummies, one per (variable, level), computed
  # once so that every term uses the same alias. They avoid the model's
  # variable names and each other, and follow the backend's case rules.
  taken = names(col_info[["types"]])
  level_aliases = list()
  for (v in names(col_info[["levels"]])) {
    lvls = c(as.character(col_info[["levels"]][[v]]), "_constant")
    a = safe_identifiers(paste0(v, lvls), taken = taken, fold_case = fold_case)
    taken = c(taken, a)
    level_aliases[[v]] = data.frame(alias = a, var = v, level = lvls)
  }
  level_aliases = do.call(rbind, c(list(data.frame(alias = character(), var = character(), level = character())), level_aliases))
  col_info[["level_aliases"]] = level_aliases
  components = unique(unlist(lapply(term_labels, split_interaction_term)))
  powers = components[vapply(components, function(t) !is.null(parse_power_term(t)), logical(1))]
  col_info[["power_aliases"]] = if (length(powers)) setNames(safe_identifiers(
    paste0("dbreg_power_", seq_along(powers), "__"), taken = taken, fold_case = fold_case), powers) else character()

  fe_terms = vapply(fe_vars, function(name) paste(deparse(as.name(name), backtick = TRUE), collapse = ""), character(1))
  result = list(select_exprs = character(), col_names = character(), aliases = character())
  
  # Without an intercept, the first factor main effect (incl. a plain i(f))
  # keeps all its levels
  full_factor = if (!intercept) {
    is_factor_main = vapply(seq_along(term_labels), function(k) {
      t = term_labels[k]
      spec = i_specs[[k]]
      if (!is.null(spec)) return(is.null(spec[["var"]]) && is.null(spec[["ref"]]))
      !grepl(":", t, fixed = TRUE) && identical(col_info[["types"]][[t]], "factor")
    }, logical(1))
    term_labels[is_factor_main][1]
  }
  i_columns = empty_i_columns()
  for (k in seq_along(term_labels)) {
    term = term_labels[k]
    if (!is.null(i_specs[[k]])) {
      expanded = expand_i_term(i_specs[[k]], col_info, k, keep_all_levels = identical(term, full_factor))
      i_columns = rbind(i_columns, expanded[["columns"]])
    } else {
      expanded = expand_term(
        term, col_info, expand, sep,
        all_terms = unique(c(term_labels, fe_terms)),
        fe_vars = fe_vars,
        keep_all_levels = identical(term, full_factor)
      )
    }
    result[["select_exprs"]] = c(result[["select_exprs"]], expanded[["sql"]])
    result[["col_names"]] = c(result[["col_names"]], expanded[["name"]])
    result[["aliases"]] = c(result[["aliases"]], expanded[["alias"]])
    result[["terms"]] = c(result[["terms"]], rep(term, length(expanded[["sql"]])))
  }
  
  result[["factor_levels"]] = col_info[["levels"]]
  result[["factor_classes"]] = col_info[["factor_classes"]]
  result[["level_aliases"]] = level_aliases
  result[["i_columns"]] = i_columns
  result[["coef_names"]] = standardize_coef_names(result[["col_names"]])
  result[["source_classes"]] = col_info[["classes"]]
  result[["i_specs"]] = lapply(col_info[["i_specs"]], function(spec) {
    spec[["source_type"]] = col_info[["types"]][[spec[["first"]]]]
    spec
  })
  result
}

#' Get column types and factor levels from database
#' @keywords internal
get_column_info = function(conn, table, vars, level_vars = character(), factor_refs = NULL,
                           i_declarations = list(), intercept = TRUE, factor_vars = character()) {
  # Build FROM clause
  from_sql = if (grepl("^FROM", table, ignore.case = TRUE)) table else paste("FROM", table)
  
  # Query to get types and sample values
  # Fetch zero rows to get the schema, then separate queries for factor levels
  # (sql_limit() handles dialects without LIMIT, e.g. SQL Server's TOP)
  backend = detect_backend(conn)[["name"]]
  # SQLite has no static type for computed expressions; a zero-row CASE or
  # polynomial projection is reported as logical. Inspect one bounded row.
  schema_query = sql_limit(glue("SELECT {paste(vars, collapse = ', ')} {from_sql}"), if (backend == "sqlite") 1 else 0, backend)
  schema = dbGetQuery(conn, schema_query)
  
  types = list()
  levels = list()
  classes = list()
  
  for (v in vars) {
    col_class = class(schema[[v]])[1]
    if (backend == "sqlite" && col_class == "logical") {
      value = dbGetQuery(conn, sql_limit(glue("SELECT {v} FROM (SELECT {v} {from_sql}) dbreg_type_sample WHERE {v} IS NOT NULL"), 1, backend))[[1L]]
      if (length(value)) col_class = class(value)[[1L]]
    }
    classes[[v]] = col_class
    if (v %in% level_vars && inherits(schema[[v]], "POSIXt")) {
      stop("Date-time variables cannot be i() factors: ", v, ". Convert to a date or character first.", call. = FALSE)
    }
    # Treat character, factor, logical as categorical
    is_cat = col_class %in% c("character", "factor", "logical") || v %in% names(factor_refs)
    types[[v]] = if (is_cat) "factor" else "numeric"
    if ((is_cat || v %in% factor_vars) && !v %in% level_vars) {
      # All distinct levels (a sample could miss rare levels)
      lvl_query = glue("SELECT DISTINCT {v} FROM (SELECT {v} {from_sql}) lvl WHERE {v} IS NOT NULL ORDER BY {v}")
      levels[[v]] = dbGetQuery(conn, lvl_query)[[1]]
      if (v %in% names(factor_refs)) {
        ref = match(as.character(factor_refs[[v]]), as.character(levels[[v]]))
        if (is.na(ref)) stop("Reference level was not found in regressor `", v, "`.", call. = FALSE)
        levels[[v]] = levels[[v]][c(ref, setdiff(seq_along(levels[[v]]), ref))]
      }
    }
  }
  
  classes = vapply(schema, function(column) class(column)[[1L]], character(1))
  declarations = Filter(Negate(is.null), i_declarations)
  specifications = i_dictionaries(conn, from_sql, declarations, classes, intercept)
  for (spec in specifications) {
    first = spec[["first"]]
    levels[[first]] = spec[["values"]]
    classes[[first]] = switch(spec[["type"]], integer = "integer64", numeric = "numeric",
      logical = "logical", date = "Date", character = "character")
  }
  list(types = types, levels = levels, classes = classes, factor_classes = classes[names(levels)],
    numeric_sql = setNames(sql_cast_double(vars, conn), vars), i_specs = specifications,
    level_sql = sql_factor_literals(conn, levels, classes))
}

#' Expand a single term into SQL expressions
#'
#' Expansions are lists of parallel character vectors: `sql` (expressions),
#' `name` (column names) and `alias` (SQL-safe column names).
#' @keywords internal
expand_term = function(term, col_info, expand, sep = "_x_", all_terms = character(), fe_vars = character(), keep_all_levels = FALSE) {
  vars = split_interaction_term(term)
  is_interaction = length(vars) > 1
  
  if (!is_interaction && expand == "interactions" && is.null(parse_power_term(term)) &&
      is.null(col_info[["factor_sources"]][[term]])) {
    # Main effect, only expanding interactions: return as-is
    return(list(sql = term, name = term, alias = term))
  }
  
  # Formula term ordering is not a canonical interaction ordering. Compare
  # sets of variables so, e.g., z:history also matches history:z.
  canonical_terms = vapply(all_terms, function(term) {
    paste(sort(split_interaction_term(term)), collapse = ":")
  }, character(1))
  expansions = lapply(vars, function(v) {
    drop_ref = !keep_all_levels
    if (is_interaction && identical(col_info[["types"]][[v]], "factor")) {
      other_term = paste(sort(setdiff(vars, v)), collapse = ":")
      drop_ref = v %in% fe_vars || other_term %in% canonical_terms
    }
    expand_variable(v, col_info, expand, is_interaction, drop_ref)
  })
  
  # Cross product for interactions
  Reduce(function(a, b) cross_product(a, b, sep), expansions, right = TRUE)
}

#' Expand a single variable into SQL expression(s)
#' @keywords internal
expand_variable = function(var, col_info, expand, in_interaction, drop_ref = TRUE) {
  power = parse_power_term(var)
  if (!is.null(power)) {
    source = col_info[["numeric_sql"]][[power[["variable"]]]]
    if (is.null(source)) source = power[["variable"]]
    return(list(sql = sprintf("POWER(%s, %d)", source, power[["exponent"]]), name = var, alias = unname(col_info[["power_aliases"]][var])))
  }
  is_factor = identical(col_info[["types"]][[var]], "factor")
  
  # Only expand factors if: expand="all", OR variable is part of interaction
  forced = !is.null(col_info[["factor_sources"]][[var]])
  if (!(is_factor && (expand == "all" || in_interaction || forced))) {
    # Numeric or unexpanded factor: pass through
    sql = if (!is_factor) col_info[["numeric_sql"]][[var]] else NULL
    if (is.null(sql)) sql = var
    return(list(sql = sql, name = var, alias = var))
  }
  lvls = col_info[["levels"]][[var]]
  level_sql = col_info[["level_sql"]][[var]]
  if (is.null(level_sql)) level_sql = paste0("'", gsub("'", "''", lvls, fixed = TRUE), "'")
  la = col_info[["level_aliases"]]
  la = la[la[["var"]] == var, , drop = FALSE]
  if (length(lvls) < 2) {
    return(list(sql = "1", name = paste0(var, "_constant"), alias = la[["alias"]][la[["level"]] == "_constant"]))
  }
  # Drop first level (reference) only if appropriate
  if (drop_ref) {
    lvls = lvls[-1]
    level_sql = level_sql[-1]
  }
  list(
    sql = paste0("CASE WHEN ", if (forced) col_info[["factor_sources"]][[var]] else var,
      " = ", level_sql, " THEN 1.0 ELSE 0.0 END"),
    name = paste0(var, lvls),
    alias = la[["alias"]][match(lvls, la[["level"]])]
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
    alias = paste0(a[["alias"]][i], "_x_", b[["alias"]][j])
  )
}

empty_i_columns = function() {
  data.frame(alias = character(), var = character(), level = character(), by = character())
}

# Argument template for i() calls (as in fixest::i())
i_template = function(factor_var, var, ref, keep, bin, ref2, keep2, bin2) NULL

#' Parse a `fixest`-style `i()` term
#'
#' @param term A term label.
#' @param env Environment in which to evaluate `ref` and `keep` (e.g. the
#'   formula's), or `NULL` to only validate the term.
#' @return `NULL` if `term` is not an `i()` call, else a list with the factor
#'   `f`, interacted variable `var` and the `ref` and `keep` levels (each
#'   `NULL` if absent or not evaluated).
#' @keywords internal
parse_sql_i_term = function(term, env = NULL) {
  e = if (is.call(term)) term else tryCatch(str2lang(term), error = function(e) NULL)
  if (!is.call(e)) return(NULL)
  head = e[[1L]]
  qualified = is.call(head) && identical(head[[1L]], as.name("::")) &&
    as.character(head[[2L]]) %in% c("dbreg", "fixest") && identical(head[[3L]], as.name("i"))
  if (!identical(head, as.name("i")) && !qualified) return(NULL)
  e[[1L]] = as.name("i")
  label = paste(deparse(e, width.cutoff = 500L), collapse = "")
  e = tryCatch(match.call(i_template, e), error = function(err) {
    stop("Invalid i() term: ", label, " (", conditionMessage(err), ")", call. = FALSE)
  })
  unsupported = intersect(names(e), c("bin", "ref2", "keep2", "bin2"))
  if (length(unsupported)) {
    stop("i() arguments not supported yet: ", paste(unsupported, collapse = ", "), call. = FALSE)
  }
  # An expression in the second position denotes reference levels; a symbol
  # denotes a numeric/logical source column. Named ref/keep use the formula env.
  if (!is.null(e[["var"]]) && !is.symbol(e[["var"]])) {
    if ("ref" %in% names(e)) stop("Specify positional or named i() references, not both.", call. = FALSE)
    e[["ref"]] = e[["var"]]
    e[["var"]] = NULL
  }
  args = list(e[["factor_var"]], e[["var"]])
  if (any(vapply(args, function(x) is.symbol(x) && startsWith(as.character(x), "i."), logical(1)))) {
    stop("The i. prefix in i() terms is not supported yet: ", label, call. = FALSE)
  }
  if (!is.symbol(args[[1L]]) || !(is.null(args[[2L]]) || is.symbol(args[[2L]]))) {
    stop("i() terms must take plain variable names, e.g. i(period, treat, ref = 5): ", label, call. = FALSE)
  }
  value = function(arg) {
    expression = e[[arg]]
    if (is.null(expression) || is.null(env)) return(NULL)
    result = tryCatch(eval(expression, env), error = function(err) {
      stop("Could not evaluate `", arg, "` in ", label, ": ", conditionMessage(err), call. = FALSE)
    })
    if (!is.atomic(result) || anyNA(result) ||
        (is.numeric(result) && any(!is.finite(result)))) {
      stop("i() references/keep values must be finite, nonmissing vectors.", call. = FALSE)
    }
    result
  }
  ref = value("ref")
  keep = value("keep")
  if (!is.null(env)) {
    # Freeze evaluated selectors before compacting the formula environment.
    for (arg in intersect(c("ref", "keep"), names(e))) {
      e[arg] = list(if (arg == "ref") ref else keep)
    }
  }
  first = as.character(e[["factor_var"]])
  second = if (!is.null(e[["var"]])) as.character(e[["var"]])
  list(f = first, var = second, first = first, second = second,
    ref = ref, keep = keep, explicit = any(c("ref", "keep") %in% names(e)),
    call = e, term = paste(deparse(e, width.cutoff = 500L), collapse = ""))
}

# Selection is shared by native SQL expansion and retained typed dictionaries.
i_level_indices = function(value, keys, labels = as.character(keys)) {
  if (is.character(value)) {
    found = integer()
    for (item in value) {
      indices = if (startsWith(item, "@!")) {
        which(!grepl(substring(item, 3L), labels))
      } else if (startsWith(item, "@")) {
        grep(substring(item, 2L), labels)
      } else {
        exact = match(item, labels)
        if (is.na(exact)) pmatch(tolower(item), tolower(labels), duplicates.ok = TRUE) else exact
      }
      found = union(found, indices[!is.na(indices)])
    }
  } else {
    found = unique(match(value, keys))
    found = found[!is.na(found)]
  }
  if (!length(found)) stop("The i() reference/keep value(s) match no level.", call. = FALSE)
  found
}

#' Expand an `i()` term into SQL expressions
#'
#' Without `var` or `ref`, `i(f)` drops the first level, like a factor main
#' effect (unless `keep_all_levels`). Columns get opaque aliases, since level
#' text is arbitrary.
#' @param spec Result of [parse_sql_i_term()].
#' @param k Term index (to make aliases unique).
#' @keywords internal
expand_i_term = function(spec, col_info, k, keep_all_levels = FALSE) {
  f = spec[["f"]]
  var = spec[["var"]]
  lvls = col_info[["levels"]][[f]]
  if (inherits(lvls, "POSIXt")) {
    stop("Date-time variables can't be i() factors yet: ", f, ". Convert to a date or character first.", call. = FALSE)
  }
  dictionaries = col_info[["i_specs"]]
  matches = vapply(dictionaries, function(item) identical(item[["term"]], spec[["term"]]), logical(1))
  dictionary = if (any(matches)) dictionaries[[which(matches)[1L]]]
  lvl_chr = if (is.null(dictionary)) as.character(lvls) else dictionary[["labels"]]
  # As in fixest: ref values must match some level; with neither ref nor keep,
  # a plain i(f) drops its first level like a factor main effect (any other
  # collinearity is left to the fit)
  keep = rep(TRUE, length(lvls))
  if (!is.null(spec[["ref"]])) {
    reference = if (isTRUE(spec[["ref"]])) 1L else i_level_indices(spec[["ref"]], lvls, lvl_chr)
    keep[reference] = FALSE
  } else if (is.null(var) && !isTRUE(spec[["explicit"]]) && is.null(spec[["keep"]]) && !keep_all_levels) {
    keep[1] = FALSE
  }
  if (!is.null(spec[["keep"]])) {
    keep = keep & seq_along(lvls) %in% i_level_indices(spec[["keep"]], lvls, lvl_chr)
  }
  lvls = lvls[keep]
  lvl_chr = lvl_chr[keep]
  if (!length(lvls)) {
    return(list(sql = character(), name = character(), alias = character(), columns = empty_i_columns()))
  }
  literal = col_info[["level_sql"]][[f]][keep]
  if (is.null(literal)) literal = sql_level_literals(lvls)
  cond = paste0(f, " = ", literal)
  if (is.null(var)) {
    sql = paste0("CASE WHEN ", cond, " THEN 1.0 ELSE 0.0 END")
  } else if (identical(col_info[["classes"]][[var]], "logical")) {
    logical = col_info[["logical_sql"]][[var]]
    if (is.null(logical)) logical = var
    sql = paste0("CASE WHEN ", cond, " AND (", logical, ") THEN 1.0 ELSE 0.0 END")
  } else if (identical(col_info[["types"]][[var]], "numeric")) {
    sql = paste0("CASE WHEN ", cond, " THEN ", var, " ELSE 0 END")
  } else {
    stop("The second argument of i() must be numeric or logical, but '", var, "' is not.", call. = FALSE)
  }
  alias = sprintf("dbreg_i%d_%d", k, seq_along(lvls))
  name = paste0(f, "::", lvl_chr, if (!is.null(var)) paste0(":", var))
  list(
    sql = sql,
    name = name,
    alias = alias,
    columns = data.frame(alias = alias, var = f, level = lvl_chr, by = if (is.null(var)) NA_character_ else var)
  )
}

# Restricted i() declarations compile to ordinary SQL design columns. The same
# saved level dictionary drives R and SQL prediction; newdata discovers no levels.

contains_i_term = function(expr) {
  if (!is.call(expr)) return(FALSE)
  !is.null(parse_sql_i_term(expr)) || any(vapply(as.list(expr)[-1L],contains_i_term,logical(1)))
}

i_column_limit = function() {
  value = getOption("dbreg.i_max_columns",1000L)
  if (!is.numeric(value) || length(value)!=1L || !is.finite(value) || value<1 || value!=floor(value)) stop("dbreg.i_max_columns must be a positive integer.",call.=FALSE)
  value
}

i_selected_items = function(value, keys, labels, type) {
  if (is.logical(value) && type %in% c("integer", "numeric")) value = as.numeric(value)
  if (is.numeric(value) && type == "logical") value = as.logical(value)
  if (type == "integer" && !is.character(value)) value = as.character(value)
  i_level_indices(value, keys, labels)
}

make_i_spec = function(declaration,keys,type,second_class=NULL,intercept=TRUE) {
  labels = as.character(keys)
  if (type=="numeric" && anyDuplicated(labels)) labels = dbglm_number(keys)
  encoded = prediction_key_part(keys,type)
  if (anyNA(encoded) || anyDuplicated(encoded) || anyDuplicated(labels)) stop("Could not preserve i() levels exactly.",call.=FALSE)
  selected = seq_along(keys)
  if (!is.null(declaration[["second"]]) && !second_class %in% c("numeric","integer","integer64","logical")) stop("The second i() variable must be numeric or logical; categorical interactions use ordinary factor syntax.",call.=FALSE)
  ref = declaration[["ref"]]
  if (isTRUE(ref) || !declaration[["explicit"]] && intercept && is.null(declaration[["second"]])) {
    selected = setdiff(selected,1L)
  } else if (!is.null(ref)) selected = setdiff(selected,i_selected_items(ref,keys,labels,type))
  if (!is.null(declaration[["keep"]])) selected = intersect(selected,i_selected_items(declaration[["keep"]],keys,labels,type))
  names = if(length(selected)) paste0(declaration[["first"]],"::",labels[selected],
    if(!is.null(declaration[["second"]])) paste0(":",declaration[["second"]])) else character()
  list(first=declaration[["first"]],second=declaration[["second"]],second_class=second_class,
    values=keys,keys=encoded,type=type,labels=labels,selected=selected,names=names,term=declaration[["term"]])
}

i_match_values = function(spec,value,strict=TRUE) {
  valid = !is.na(value)
  if (is.numeric(value) && !inherits(value,"integer64")) valid = valid & is.finite(value)
  keys = rep(NA_character_,length(value))
  keys[valid] = prediction_key_part(value[valid],spec[["type"]])
  index = match(keys,spec[["keys"]])
  if (strict && any(valid & is.na(index))) stop("newdata contains unseen i() levels.",call.=FALSE)
  index
}

evaluate_i_spec = function(spec,first,second=NULL) {
  index = i_match_values(spec,first)
  columns = spec[["selected"]]
  if (!length(columns)) return(matrix(numeric(),length(first),0L))
  result = outer(index,columns,`==`)*1.0
  if (!is.null(spec[["second"]])) {
    if (is.factor(second)) second = as.character(second)
    if (spec[["second_class"]]=="logical" && !is.logical(second)) second = as.logical(second)
    result = result*as.numeric(second)
  }
  colnames(result) = spec[["names"]]
  result
}

# Only bounded, typed dictionaries cross the SQL/R boundary. Integer levels are
# transferred as decimal strings, never through the driver's double conversion.
i_dictionaries = function(conn,from_sql,declarations,classes,intercept) {
  result = list()
  if (!length(declarations)) return(result)
  limit = i_column_limit()
  source = sub("^FROM\\s+","",from_sql,ignore.case=TRUE)
  for (declaration in declarations) {
    first = declaration[["first"]]
    metadata = prediction_key_metadata(conn,source,first)
    type = metadata[["types"]][[1L]]
    if (type=="float") stop("i() floating-point categories require a DOUBLE source column.",call.=FALSE)
    value = if(type=="integer") paste0("CAST(",first," AS STRING)") else first
    # Group by the original typed value; textual ordering would put 10 before 2.
    keys = backend_get_query(conn,paste0("SELECT ",value," AS value ",from_sql," GROUP BY ",first,
      " ORDER BY ",first," LIMIT ",sprintf("%.0f",limit+2L)))[[1L]]
    if (length(keys)>limit+1L) stop("i() level dictionary exceeds dbreg.i_max_columns.",call.=FALSE)
    # Match R/fixest character collation on this bounded dictionary. Preserve
    # native enum order, numeric order and exact decimal integer text.
    if (type=="character" && !startsWith(tolower(metadata[["sql_types"]][[1L]]),"enum")) keys = sort(keys)
    if (type=="integer" && prediction_integer_sql_type(keys,conn)!="BIGINT") stop("i() integer levels must fit signed BIGINT.",call.=FALSE)
    second = declaration[["second"]]
    spec = make_i_spec(declaration,keys,type,if(!is.null(second)) classes[[second]],intercept || length(result)>0L)
    key = make.unique(c(names(result), first))[[length(result) + 1L]]
    result[[key]] = spec
    if (sum(vapply(result,function(x) length(x[["selected"]]),integer(1)))>limit) stop("i() columns exceed dbreg.i_max_columns.",call.=FALSE)
  }
  result
}

sql_level_literals = function(values) {
  if (is.numeric(values) && !inherits(values, "Date")) return(sprintf("%.17g", values))
  paste0("'", gsub("'", "''", as.character(values), fixed = TRUE), "'")
}

freeze_i_selectors = function(fml) {
  env = environment(fml)
  walk = function(expr) {
    if (!is.call(expr)) return(expr)
    declaration = parse_sql_i_term(expr, env = env)
    if (!is.null(declaration)) return(declaration[["call"]])
    as.call(lapply(as.list(expr), walk))
  }
  result = stats::as.formula(walk(as.call(fml)), env = baseenv())
  if (inherits(fml, "Formula")) Formula::Formula(result) else result
}

explicit_factor_sources = function(term_labels) {
  components = unique(unlist(lapply(term_labels, split_interaction_term)))
  result = character()
  for (term in components) {
    expr = tryCatch(str2lang(term), error = function(e) NULL)
    if (!is.call(expr) || !identical(expr[[1L]], as.name("factor"))) next
    if (length(expr) != 2L || !is.symbol(expr[[2L]])) {
      stop("Use factor(column) without additional arguments.", call. = FALSE)
    }
    result[[term]] = as.character(expr[[2L]])
  }
  result
}

explicit_factor_info = function(term_labels, info) {
  sources = explicit_factor_sources(term_labels)
  for (term in names(sources)) {
    source = sources[[term]]
    if (is.null(info[["levels"]][[term]])) info[["levels"]][[term]] = info[["levels"]][[source]]
    if (is.null(info[["level_sql"]][[term]])) info[["level_sql"]][[term]] = info[["level_sql"]][[source]]
    info[["types"]][[term]] = "factor"
    class = info[["factor_classes"]][term]
    if (!length(class) || is.na(class)) class = info[["classes"]][source]
    info[["factor_classes"]][term] = class
    if (is.null(info[["factor_sources"]][[term]])) info[["factor_sources"]][[term]] = source
  }
  info
}

# Names used as metadata keys can themselves be factor() expressions.
restore_term_names = function(names, map) {
  vapply(names, function(name) {
    expr = tryCatch(str2lang(name), error = function(e) NULL)
    if (!is.call(expr)) return(restore_names(name, map))
    paste(deparse(rename_expr_vars(expr, map), width.cutoff = 500L), collapse = "")
  }, character(1), USE.NAMES = FALSE)
}
