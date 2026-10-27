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
  result = sql_model_matrix_internal(formula, conn, table, expand, sep, fe_vars,
    intercept, factor_refs)
  attr(result, "sql_level_rows") = NULL
  result
}

# Internal callers additionally reserve source roles outside the RHS.
sql_model_matrix_internal = function(
  formula,
  conn,
  table,
  expand = c("all", "interactions"),
  sep = "_x_",
  fe_vars = character(),
  intercept = TRUE,
  factor_refs = NULL,
  reserved_names = character()
) {
  expand = match.arg(expand)
  reserved_names = unique(c(reserved_names, model_vars(formula), fe_vars))

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
      columns = empty_linear_columns(),
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
  col_info = get_column_info(conn, table, unique(unlist(term_vars)),
    level_vars = i_factors,
    factor_refs = factor_refs, i_declarations = i_specs, intercept = intercept,
    factor_vars = unname(forced)
  )
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
    fold_case = folds_identifier_case(conn), i_specs = i_specs, reserved_names = reserved_names
  )
  ordinary = term_labels[vapply(i_specs, is.null, logical(1))]
  result[["ordinary_rhs"]] = stats::reformulate(if (length(ordinary)) ordinary else "1",
    intercept = attr(tt, "intercept") == 1L, env = baseenv()
  )
  result
}

# The fitted-design path supplies saved levels instead of discovering them from
# newdata. Both paths use exactly the same expansion and reference decisions.
sql_model_matrix_from_info = function(term_labels, col_info, expand = "all", sep = "_x_",
                                      fe_vars = character(), intercept = TRUE, fold_case = FALSE, i_specs = lapply(term_labels, parse_sql_i_term, env = baseenv()), reserved_names = character()) {
  col_info = explicit_factor_info(term_labels, col_info)
  # SQL-safe aliases for factor dummies, one per (variable, level), computed
  # once so that every term uses the same alias. They avoid the model's
  # variable names and each other, and follow the backend's case rules.
  source_names = unique(c(names(col_info[["types"]]), fe_vars, reserved_names))
  taken = source_names
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
    paste0("dbreg_power_", seq_along(powers), "__"),
    taken = taken, fold_case = fold_case
  ), powers) else character()

  fe_terms = vapply(fe_vars, function(name) paste(deparse(as.name(name), backtick = TRUE), collapse = ""), character(1))
  result = list(select_exprs = character(), col_names = character(), aliases = character(),
    columns = empty_linear_columns())

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
  i_positions = integer()
  level_rows = integer()
  for (k in seq_along(term_labels)) {
    term = term_labels[k]
    if (!is.null(i_specs[[k]])) {
      expanded = expand_i_term(i_specs[[k]], col_info, k, keep_all_levels = identical(term, full_factor))
      i_positions = c(i_positions, length(result[["aliases"]]) + seq_along(expanded[["alias"]]))
      i_columns = rbind(i_columns, expanded[["columns"]])
    } else {
      expanded = expand_term(
        term, col_info, expand, sep,
        all_terms = unique(c(term_labels, fe_terms)),
        fe_vars = fe_vars,
        keep_all_levels = identical(term, full_factor)
      )
    }
    level_rows = c(level_rows, if (is.null(expanded[["level_row"]])) rep(NA_integer_, length(expanded[["alias"]])) else expanded[["level_row"]])
    result[["select_exprs"]] = c(result[["select_exprs"]], expanded[["sql"]])
    result[["col_names"]] = c(result[["col_names"]], expanded[["name"]])
    result[["aliases"]] = c(result[["aliases"]], expanded[["alias"]])
    result[["terms"]] = c(result[["terms"]], rep(term, length(expanded[["sql"]])))
    count = length(expanded[["sql"]])
    if (count) {
      columns = data.frame(id = paste0("t", k, "b", seq_len(count)), term = k,
        basis = seq_len(count), alias = expanded[["alias"]],
        label = expanded[["name"]], role = "slope", stringsAsFactors = FALSE)
      columns[["recipe"]] = expanded[["recipe"]]
      columns[["cell_recipe"]] = rep(list(NULL), count)
      result[["columns"]] = rbind(result[["columns"]], columns)
    }
  }

  # Allocate targets only after the whole expansion is known, before any SQL
  # projection. Source expressions and display labels remain untouched.
  generated = vapply(result[["columns"]][["recipe"]], function(recipe) {
    !recipe[["kind"]] %in% c("raw", "numeric")
  }, logical(1))
  result[["aliases"]] = allocate_sql_aliases(result[["aliases"]], generated,
    source_names, fold_case)
  result[["columns"]][["alias"]] = result[["aliases"]]
  i_columns[["alias"]] = result[["aliases"]][i_positions]
  present = which(!is.na(level_rows))
  level_aliases[["alias"]][level_rows[present]] = result[["aliases"]][present]

  result[["factor_levels"]] = col_info[["levels"]]
  result[["factor_classes"]] = col_info[["factor_classes"]]
  result[["level_aliases"]] = level_aliases
  result[["i_columns"]] = i_columns
  result[["coef_names"]] = standardize_coef_names(result[["col_names"]])
  result[["source_classes"]] = col_info[["classes"]]
  attr(result, "sql_level_rows") = level_rows
  result[["i_specs"]] = lapply(col_info[["i_specs"]], function(spec) {
    spec[["source_type"]] = col_info[["types"]][[spec[["first"]]]]
    spec
  })
  result
}

# One allocator for generated SQL targets. Raw passthroughs retain their names;
# generated targets reserve every source role and every other generated target.
allocate_sql_aliases = function(aliases, generated, reserved_names = character(), fold_case = FALSE) {
  reserved = unique(c(reserved_names, aliases[!generated]))
  positions = which(generated)
  conflicts = tolower(aliases[positions]) %in% tolower(reserved) |
    duplicated(tolower(aliases[positions]))
  changed = positions[conflicts]
  aliases[changed] = safe_identifiers(aliases[changed],
    taken = c(reserved, aliases[setdiff(seq_along(aliases), changed)]), fold_case = fold_case)
  aliases
}

# These names may be projected beside the design, even when absent from its RHS.
sql_design_reserved = function(inputs) {
  unique(c(model_vars(inputs[["fml"]]), inputs[["yvar"]], inputs[["xvars"]],
    inputs[["fe"]], unlist(inputs[["fe_slopes"]], use.names = FALSE),
    inputs[["weights"]], inputs[["offset"]], inputs[["cluster_var"]],
    inputs[["retained_key_safe"]], inputs[["row_id_safe"]], inputs[["row_id"]],
    unlist(inputs[["composite_fes"]], use.names = FALSE)))
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
  schema_query = sql_limit(glue("SELECT {paste(vars, collapse = ', ')} {from_sql}"), 0, backend)
  schema = dbGetQuery(conn, schema_query)

  types = list()
  levels = list()
  classes = list()

  for (v in vars) {
    col_class = class(schema[[v]])[1]
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
    classes[[first]] = switch(spec[["type"]],
      integer = "integer64",
      numeric = "numeric",
      logical = "logical",
      date = "Date",
      character = "character"
    )
  }
  list(
    types = types, levels = levels, classes = classes, factor_classes = classes[names(levels)],
    numeric_sql = setNames(sql_cast_double(vars, conn), vars), i_specs = specifications,
    level_sql = sql_factor_literals(conn, levels, classes)
  )
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
    return(list(sql = term, name = term, alias = term, recipe = list(list(kind = "raw", variable = term))))
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
    return(list(sql = sprintf("POWER(%s, %d)", source, power[["exponent"]]), name = var, alias = unname(col_info[["power_aliases"]][var]),
      recipe = list(list(kind = "power", variable = power[["variable"]], exponent = power[["exponent"]]))))
  }
  is_factor = identical(col_info[["types"]][[var]], "factor")

  # Only expand factors if: expand="all", OR variable is part of interaction
  forced = !is.null(col_info[["factor_sources"]][[var]])
  if (!(is_factor && (expand == "all" || in_interaction || forced))) {
    # Numeric or unexpanded factor: pass through
    sql = if (!is_factor) col_info[["numeric_sql"]][[var]] else NULL
    if (is.null(sql)) sql = var
    return(list(sql = sql, name = var, alias = var,
      recipe = list(list(kind = if (is_factor) "raw" else "numeric", variable = var))))
  }
  lvls = col_info[["levels"]][[var]]
  level_sql = col_info[["level_sql"]][[var]]
  if (is.null(level_sql)) level_sql = paste0("'", gsub("'", "''", lvls, fixed = TRUE), "'")
  la = col_info[["level_aliases"]]
  la = la[la[["var"]] == var, , drop = FALSE]
  if (length(lvls) < 2) {
    return(list(sql = "1", name = paste0(var, "_constant"), alias = la[["alias"]][la[["level"]] == "_constant"],
      recipe = list(list(kind = "constant")),
      level_row = which(col_info[["level_aliases"]][["var"]] == var & col_info[["level_aliases"]][["level"]] == "_constant")))
  }
  # Drop first level (reference) only if appropriate
  if (drop_ref) {
    lvls = lvls[-1]
    level_sql = level_sql[-1]
  }
  list(
    sql = paste0(
      "CASE WHEN ", if (forced) col_info[["factor_sources"]][[var]] else var,
      " = ", level_sql, " THEN 1.0 ELSE 0.0 END"
    ),
    name = paste0(var, lvls),
    alias = la[["alias"]][match(lvls, la[["level"]])],
    level_row = which(col_info[["level_aliases"]][["var"]] == var)[match(lvls, la[["level"]])],
    recipe = lapply(seq_along(lvls), function(i) list(kind = "indicator",
      variable = if (forced) unname(explicit_factor_sources(var)) else var,
      value = lvls[i]))
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
    recipe = Map(function(left, right) list(kind = "product", left = left, right = right),
      a[["recipe"]][i], b[["recipe"]][j])
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
  list(
    f = first, var = second, first = first, second = second,
    ref = ref, keep = keep, explicit = any(c("ref", "keep") %in% names(e)),
    call = e, term = paste(deparse(e, width.cutoff = 500L), collapse = "")
  )
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
    columns = data.frame(alias = alias, var = f, level = lvl_chr, by = if (is.null(var)) NA_character_ else var),
    recipe = lapply(seq_along(lvls), function(i) list(kind = "i", variable = f,
      second = var, term = spec[["term"]],
      level = if (is.null(dictionary)) match(lvl_chr[i], as.character(col_info[["levels"]][[f]])) else
        match(lvl_chr[i], dictionary[["labels"]])))
  )
}

# Restricted i() declarations compile to ordinary SQL design columns. The same
# saved level dictionary drives R and SQL prediction; newdata discovers no levels.

contains_i_term = function(expr) {
  if (!is.call(expr)) return(FALSE)
  !is.null(parse_sql_i_term(expr)) || any(vapply(as.list(expr)[-1L], contains_i_term, logical(1)))
}

i_column_limit = function() {
  value = getOption("dbreg.i_max_columns", 1000L)
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < 1 || value != floor(value)) stop("dbreg.i_max_columns must be a positive integer.", call. = FALSE)
  value
}

i_selected_items = function(value, keys, labels, type) {
  if (is.logical(value) && type %in% c("integer", "numeric")) value = as.numeric(value)
  if (is.numeric(value) && type == "logical") value = as.logical(value)
  if (type == "integer" && !is.character(value)) value = as.character(value)
  i_level_indices(value, keys, labels)
}

make_i_spec = function(declaration, keys, type, second_class = NULL, intercept = TRUE) {
  labels = as.character(keys)
  if (type == "numeric" && anyDuplicated(labels)) labels = dbglm_number(keys)
  encoded = prediction_key_part(keys, type)
  if (anyNA(encoded) || anyDuplicated(encoded) || anyDuplicated(labels)) stop("Could not preserve i() levels exactly.", call. = FALSE)
  selected = seq_along(keys)
  if (!is.null(declaration[["second"]]) && !second_class %in% c("numeric", "integer", "integer64", "logical")) stop("The second i() variable must be numeric or logical; categorical interactions use ordinary factor syntax.", call. = FALSE)
  ref = declaration[["ref"]]
  if (isTRUE(ref) || !declaration[["explicit"]] && intercept && is.null(declaration[["second"]])) {
    selected = setdiff(selected, 1L)
  } else if (!is.null(ref)) selected = setdiff(selected, i_selected_items(ref, keys, labels, type))
  if (!is.null(declaration[["keep"]])) selected = intersect(selected, i_selected_items(declaration[["keep"]], keys, labels, type))
  names = if (length(selected)) paste0(
    declaration[["first"]], "::", labels[selected],
    if (!is.null(declaration[["second"]])) paste0(":", declaration[["second"]])
  ) else character()
  list(
    first = declaration[["first"]], second = declaration[["second"]], second_class = second_class,
    values = keys, keys = encoded, type = type, labels = labels, selected = selected, names = names, term = declaration[["term"]]
  )
}

i_match_values = function(spec, value, strict = TRUE) {
  valid = !is.na(value)
  if (is.numeric(value) && !inherits(value, "integer64")) valid = valid & is.finite(value)
  keys = rep(NA_character_, length(value))
  keys[valid] = prediction_key_part(value[valid], spec[["type"]])
  index = match(keys, spec[["keys"]])
  if (strict && any(valid & is.na(index))) stop("newdata contains unseen i() levels.", call. = FALSE)
  index
}

evaluate_i_spec = function(spec, first, second = NULL) {
  index = i_match_values(spec, first)
  columns = spec[["selected"]]
  if (!length(columns)) return(matrix(numeric(), length(first), 0L))
  result = outer(index, columns, `==`) * 1.0
  if (!is.null(spec[["second"]])) {
    if (is.factor(second)) second = as.character(second)
    if (spec[["second_class"]] == "logical" && !is.logical(second)) second = as.logical(second)
    result = result * as.numeric(second)
  }
  colnames(result) = spec[["names"]]
  result
}

# Only bounded, typed dictionaries cross the SQL/R boundary. Integer levels are
# transferred as decimal strings, never through the driver's double conversion.
i_dictionaries = function(conn, from_sql, declarations, classes, intercept) {
  result = list()
  if (!length(declarations)) return(result)
  backend_info = detect_backend(conn)
  limit = i_column_limit()
  source = sub("^FROM\\s+", "", from_sql, ignore.case = TRUE)
  for (declaration in declarations) {
    first = declaration[["first"]]
    metadata = prediction_key_metadata(conn, source, first)
    type = metadata[["types"]][[1L]]
    if (type == "float") stop("i() floating-point categories require a DOUBLE source column.", call. = FALSE)
    value = if (type == "integer") sql_cast_text(conn, first, backend_info) else first
    # Group by the original typed value; textual ordering would put 10 before 2.
    keys = DBI::dbGetQuery(conn, paste0(
      "SELECT ", value, " AS value ", from_sql, " GROUP BY ", first,
      " ORDER BY ", first, " LIMIT ", sprintf("%.0f", limit + 2L)
    ))[[1L]]
    if (length(keys) > limit + 1L) stop("i() level dictionary exceeds dbreg.i_max_columns.", call. = FALSE)
    # Match R/fixest character collation on this bounded dictionary. Preserve
    # native enum order, numeric order and exact decimal integer text.
    if (type == "character" && !startsWith(tolower(metadata[["sql_types"]][[1L]]), "enum")) keys = sort(keys)
    if (type == "integer" && prediction_integer_sql_type(keys, conn) != "BIGINT") stop("i() integer levels must fit signed BIGINT.", call. = FALSE)
    second = declaration[["second"]]
    spec = make_i_spec(declaration, keys, type, if (!is.null(second)) classes[[second]], intercept || length(result) > 0L)
    key = make.unique(c(names(result), first))[[length(result) + 1L]]
    result[[key]] = spec
    if (sum(vapply(result, function(x) length(x[["selected"]]), integer(1))) > limit) stop("i() columns exceed dbreg.i_max_columns.", call. = FALSE)
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

# Formula and model-matrix helpers shared by OLS and GLM.

# Keep the fast matrix product for ordinary inputs. At extreme magnitudes,
# fused BLAS operations can hide Inf-Inf as Inf (and depend on batch shape).
# Re-evaluate potentially overflowing rows with the same staged arithmetic
# as the SQL predictor, using row vectors rather than another dense matrix.
prediction_matrix_product = function(x, beta) {
  value = as.vector(x %*% beta)
  if (!length(beta) || !nrow(x)) return(value)
  bound = .Machine[["double.xmax"]] / max(1, max(abs(beta))) / max(1, ncol(x)) / 2
  if (isTRUE(Matrix::norm(x, "M") <= bound)) return(value)
  risky = !is.finite(value)
  for (j in seq_along(beta)) {
    limit = .Machine[["double.xmax"]] / max(1, abs(beta[[j]])) / max(1, ncol(x)) / 2
    column = x[, j, drop = FALSE]
    if (isTRUE(Matrix::norm(column, "M") <= limit)) next
    column = as.numeric(column)
    risky = risky | !is.finite(column) | abs(column) > limit
  }
  if (any(risky)) {
    total = rep(0, sum(risky))
    for (j in seq_along(beta)) total = total + as.numeric(x[risky, j]) * beta[[j]]
    value[risky] = total
  }
  value
}

parse_power_term = function(term) {
  expr = tryCatch(parse(text = term, keep.source = FALSE)[[1]], error = function(e) NULL)
  if (is.null(expr) || !is.call(expr) || !identical(as.character(expr[[1]]), "I")) {
    return(NULL)
  }
  power = expr[[2]]
  if (!is.call(power) || !identical(as.character(power[[1]]), "^") ||
    !is.symbol(power[[2]]) || !is.numeric(power[[3]]) ||
    !as.numeric(power[[3]]) %in% c(2, 3)) {
    return(NULL)
  }
  variable = as.character(power[[2]])
  exponent = as.integer(power[[3]])
  list(
    variable = variable,
    exponent = exponent,
    sql = sprintf("POWER(%s, %d)", variable, exponent),
    name = sprintf("%s__pow%d", variable, exponent)
  )
}

split_interaction_term = function(term) {
  expr = parse(text = term, keep.source = FALSE)[[1]]
  flatten = function(x) {
    if (is.call(x) && identical(as.character(x[[1]]), ":")) {
      return(c(flatten(x[[2]]), flatten(x[[3]])))
    }
    paste(deparse(x, width.cutoff = 500L), collapse = "")
  }
  flatten(expr)
}

# SQL factor expansions use treatment coding independently of R's global options.
fitted_model_matrix = function(formula, data) {
  data = model_matrix_data(formula, data)
  categorical = names(data)[vapply(data, function(x) is.factor(x) || is.character(x) || is.logical(x), logical(1))]
  categorical = intersect(categorical, all.vars(formula))
  contrasts = lapply(data[categorical], function(x) {
    lev = if (is.factor(x)) levels(x) else levels(factor(x))
    stats::contr.treatment(lev)
  })
  if (!length(contrasts)) contrasts = NULL
  # Tiny compressed cell designs are faster as ordinary matrices. Bound an
  # upper estimate including factor interactions before allowing dense storage.
  labels = attr(stats::terms(formula, data = data), "term.labels")
  width = (length(labels) + 1L) *
    prod(vapply(contrasts, nrow, integer(1)))
  if (!any(grepl("(", labels, fixed = TRUE)) && nrow(data) * width <= 1e6) {
    return(stats::model.matrix(formula, data = data, contrasts.arg = contrasts))
  }
  tryCatch(Matrix::sparse.model.matrix(formula, data = data, contrasts.arg = contrasts),
    error = function(e) stats::model.matrix(formula, data = data, contrasts.arg = contrasts)
  )
}

# A single GLM design specification shared by SQL fitting and R prediction.
build_model_spec = function(formula, conn, table, fe_vars = character(),
                            factor_refs = NULL, reserved_names = character(), absorbed = FALSE) {
  design = sql_model_matrix_internal(
    formula = formula,
    conn = conn,
    table = table,
    expand = "all",
    sep = "_x_",
    fe_vars = fe_vars,
    factor_refs = factor_refs,
    reserved_names = reserved_names,
    intercept = attr(stats::terms(stats::formula(formula, lhs = 0, rhs = 1)), "intercept") == 1L || absorbed
  )
  # Prediction labels are distinct from typed SQL comparison values. Adjacent
  # doubles can share as.character() output, so retain the dictionary labels.
  for (spec in design[["i_specs"]]) {
    design[["factor_levels"]][[spec[["first"]]]] = spec[["labels"]]
  }
  public_names = design[["coef_names"]]
  if (anyDuplicated(public_names)) {
    duplicates = unique(public_names[duplicated(public_names)])
    stop(
      "Formula terms produce duplicate coefficient names: ",
      paste(duplicates, collapse = ", "),
      call. = FALSE
    )
  }

  aliases = character(length(public_names))
  reserved = unique(reserved_names)
  for (i in seq_along(public_names)) {
    candidate = sprintf("dbreg_internal_x_%04d__", i)
    while (candidate %in% reserved) candidate = paste0(candidate, "_")
    aliases[[i]] = candidate
    reserved = c(reserved, candidate)
  }
  list(
    rhs = if (inherits(formula, "Formula")) {
      formula(formula, lhs = 0, rhs = 1)
    } else {
      formula
    },
    source_vars = model_vars(if (inherits(formula, "Formula")) {
      formula(formula, lhs = 0, rhs = 1)
    } else {
      formula
    }),
    columns = data.frame(
      alias = aliases,
      name = public_names,
      sql = design[["select_exprs"]],
      term = design[["terms"]],
      stringsAsFactors = FALSE
    ),
    factor_levels = design[["factor_levels"]],
    factor_classes = design[["factor_classes"]],
    source_classes = design[["source_classes"]],
    i_specs = design[["i_specs"]],
    ordinary_rhs = design[["ordinary_rhs"]],
    factor_refs = factor_refs
  )
}

evaluate_model_spec = function(model_spec, newdata) {
  # Public newdata always contains source variables. Extra columns must never
  # bypass the fitted formula merely by matching generated SQL aliases.
  missing_vars = setdiff(model_spec[["source_vars"]], names(newdata))
  if (length(missing_vars)) {
    stop(
      "newdata is missing required regressor columns: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }

  for (variable in names(model_spec[["factor_levels"]])) {
    if (variable %in% names(model_spec[["i_specs"]])) next
    training_levels = model_spec[["factor_levels"]][[variable]]
    values = as.character(newdata[[variable]])
    unseen = setdiff(unique(values[!is.na(values)]), as.character(training_levels))
    if (length(unseen)) {
      stop(
        "newdata contains unseen levels for `", variable, "`: ",
        paste(unseen, collapse = ", "),
        call. = FALSE
      )
    }
    newdata[[variable]] = factor(values, levels = as.character(training_levels))
  }

  rhs = if (length(model_spec[["i_specs"]])) model_spec[["ordinary_rhs"]] else model_spec[["rhs"]]
  design_data = model_matrix_data(rhs, newdata)
  matrix = fitted_model_matrix(rhs, design_data)
  if ("(Intercept)" %in% colnames(matrix)) {
    matrix = matrix[, colnames(matrix) != "(Intercept)", drop = FALSE]
  }
  for (spec in model_spec[["i_specs"]]) {
    matrix = cbind(matrix, evaluate_i_spec(
      spec, newdata[[spec[["first"]]]],
      if (!is.null(spec[["second"]])) newdata[[spec[["second"]]]]
    ))
  }
  expected = model_spec[["columns"]][["name"]]
  missing_columns = setdiff(expected, colnames(matrix))
  if (length(missing_columns)) {
    stop(
      "Could not reconstruct fitted design columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  matrix = matrix[, expected, drop = FALSE]
  colnames(matrix) = model_spec[["columns"]][["alias"]]
  matrix
}

# Composite FEs are tuples, never delimiter-concatenated SQL strings.
linear_composite_frame = function(conn, from_statement, parsed) {
  components = parsed[["fe_components"]]
  composite = which(lengths(components) > 1L)
  if (!length(composite)) return(list(from_statement = from_statement, parsed = parsed, composites = list()))
  columns = names(DBI::dbGetQuery(conn, sql_limit(paste("SELECT *", from_statement), 0, detect_backend(conn)[["name"]])))
  aliases = safe_identifiers(paste0("dbreg_joint_fe_", composite, "__"), taken = columns, fold_case = folds_identifier_case(conn))
  parsed[["fe"]][composite] = aliases
  expressions = vapply(seq_along(composite), function(j) {
    paste0("DENSE_RANK() OVER (ORDER BY ", paste(components[[composite[[j]]]], collapse = ", "), ") AS ", aliases[[j]])
  }, character(1))
  from_statement = paste0(
    "FROM (SELECT dbreg_source.*, ", paste(expressions, collapse = ", "),
    " FROM (SELECT * ", from_statement, ") dbreg_source) dbreg_composite"
  )
  rhs = stats::formula(parsed[["fml"]], lhs = 0, rhs = 1)[[2]]
  lhs = stats::formula(parsed[["fml"]], lhs = 1, rhs = 0)[[2]]
  fe_rhs = Reduce(function(a, b) call("+", a, b), lapply(parsed[["fe"]], as.name))
  parsed[["fml"]] = Formula::Formula(stats::as.formula(call("~", lhs, call("|", rhs, fe_rhs))))
  specs = setNames(lapply(composite, function(i) components[[i]]), aliases)
  list(from_statement = from_statement, parsed = parsed, composites = specs)
}

validate_factor_refs = function(refs, variables) {
  if (is.null(refs)) return(NULL)
  refs = as.list(refs)
  if (is.null(names(refs)) || anyNA(names(refs)) || any(!nzchar(names(refs))) ||
    anyDuplicated(names(refs)) || any(!names(refs) %in% variables) || any(lengths(refs) != 1L) ||
    any(vapply(refs, anyNA, logical(1)))) {
    stop("`factor_refs` must name regressor columns with one nonmissing reference each.", call. = FALSE)
  }
  refs
}

# Base factor() stringifies doubles before matching and can merge adjacent
# representable keys. Match numeric values first, keeping ordinary labels when
# they are unambiguous and using round-trip labels only when necessary.
exact_numeric_factor = function(x) {
  if (inherits(x, "integer64")) {
    values = sort(unique(x))
    result = factor(as.character(x), levels = as.character(values))
    attr(result, "dbreg_key_values") = values
    attr(result, "dbreg_key_type") = "integer"
    return(result)
  }
  if (!is.numeric(x)) return(factor(x))
  values = sort(unique(x))
  labels = as.character(values)
  if (anyDuplicated(labels)) labels = sprintf("%.17g", values)
  codes = match(x, values)
  names(codes) = names(x)
  structure(codes,
    levels = labels, class = "factor", dbreg_key_values = values,
    dbreg_key_type = if (is.integer(x)) "integer" else "numeric"
  )
}

# Exact decimal integer labels can be sorted without a new runtime dependency.
exact_integer_factor = function(x) {
  values = unique(as.character(x[!is.na(x)]))
  negative = startsWith(values, "-")
  lo = values[negative]
  hi = values[!negative]
  levels = c(lo[order(nchar(lo), lo, decreasing = TRUE)], hi[order(nchar(hi), hi)])
  result = factor(as.character(x), levels = levels)
  attr(result, "dbreg_key_values") = levels
  attr(result, "dbreg_key_type") = "integer"
  result
}

# Re-read only an exceptional lossy BIGINT result. Ordinary compression keeps
# its single query. Cast after SQL grouping, preserving both exact group keys
# and numeric level ordering, without changing the caller's DBI configuration.
read_compressed_frame = function(inputs, sql, from_statement = inputs[["from_statement"]], integer_keys = inputs[["fe"]]) {
  conn = inputs[["conn"]]
  data = DBI::dbGetQuery(conn, sql)
  risky = integer_keys[vapply(integer_keys, function(f) {
    value = data[[f]]
    is.numeric(value) && !inherits(value, "integer64") && any(abs(value) >= 2^53, na.rm = TRUE)
  }, logical(1))]
  integers = character()
  if (length(risky)) {
    backend = detect_backend(conn)
    if (backend[["name"]] == "duckdb" || isTRUE(backend[["is_databricks"]])) {
      q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
      types = DBI::dbGetQuery(conn, paste0(
        "SELECT ", paste0("typeof(", q(risky), ") AS t", seq_along(risky), collapse = ","),
        " ", from_statement, " LIMIT 1"
      ))
      integers = risky[vapply(types, function(type) grepl("^(u?(bigint|hugeint)|decimal\\([0-9]+,0\\))$", tolower(type[[1]])), logical(1))]
      if (length(integers)) {
        cast = sql_text_type(conn, backend)
        fields = c(q(setdiff(names(data), integers)), paste0("CAST(", q(integers), " AS ", cast, ") AS ", q(integers)))
        sql = paste0("SELECT ", paste(fields, collapse = ","), " FROM (", sql, ") dbreg_exact_cells")
        data = DBI::dbGetQuery(conn, sql)
      }
    }
  }
  # Keep typed identifiers intact, but cross the same numeric boundary used by
  # cross_moments() for sufficient statistics. integer64 arithmetic otherwise
  # truncates fractional residuals, even when cell means/coefficients are doubles.
  moments = intersect(c(
    "sum_w", "sum_wy", "sum_wy_sq", "sum_w2", "sum_w2y", "sum_w2y_sq",
    "sum_w_gc", "sum_wy_gc", "mean_Y", "wts"
  ), names(data))
  for (name in moments) data[[name]] = as.numeric(data[[name]])
  list(data = data, sql = sql, integers = integers)
}

# Only compressed fits retain explicit FE coefficients. Their tuple-to-rank
# mapping comes from cells already collected for estimation, with no extra scan.
finalize_linear_design = function(result, inputs) {
  if (!is.null(result[["linear_design"]])) {
    result[["linear_design"]] = restore_native_linear_layout(result[["linear_design"]],
      inputs[["source_map"]], inputs[["name_map"]], rownames(result[["coeftable"]]))
  }
  result[["i_specs"]] = inputs[["i_specs"]]
  result[["offset"]] = inputs[["offset"]]
  result[["fixef_rm"]] = inputs[["fixef_rm"]]
  result[["pruning_rounds"]] = if (is.null(inputs[["pruning"]])) 0L else inputs[["pruning"]][["rounds"]]
  result[["exclusions"]] = if (is.null(inputs[["pruning"]])) c(singleton = 0, infinite_coef = 0) else inputs[["pruning"]][["counts"]]
  levels = inputs[["factor_levels"]]
  if (length(levels)) names(levels) = restore_term_names(names(levels), inputs[["source_map"]])
  result[["factor_levels"]] = levels
  classes = inputs[["factor_classes"]]
  if (length(classes)) names(classes) = restore_term_names(names(classes), inputs[["source_map"]])
  result[["factor_classes"]] = classes
  composites = inputs[["composite_fes"]]
  result[["composite_fes"]] = lapply(composites, function(vars) restore_names(vars, inputs[["source_map"]]))
  result[["composite_levels"]] = inputs[["composite_levels"]]
  result[["fe_labels"]] = inputs[["fe_labels"]]
  result[["varying_slopes"]] = lapply(inputs[["fe_slopes"]], restore_names, map = inputs[["source_map"]])
  fml = if (!is.null(inputs[["original_fml"]])) inputs[["original_fml"]] else result[["fml"]]
  # The supported design uses columns and built-in transformations only.
  # Retaining the caller/strategy frame can capture the entire original panel
  # (and a live connection) through an otherwise compact formula object.
  result[["fml"]] = if (!is.null(inputs[["original_fml"]])) fml else freeze_i_selectors(fml)
  design = result[["compressed_design"]]
  if (!is.null(design)) {
    environment(design[["formula"]]) = baseenv()
    design[["source_names"]] = setNames(design[["source_names"]], restore_names(design[["source_names"]], inputs[["name_map"]]))
    design[["columns"]] = setNames(rownames(result[["coeftable"]]), design[["columns"]])
    if (length(design[["trend_sources"]])) {
      design[["prediction_trends"]] = stats::setNames(
        restore_names(unname(design[["trend_sources"]]), inputs[["source_map"]]), names(design[["trend_sources"]])
      )
    }
    result[["compressed_design"]] = design
  }
  result
}

weighted_group_mean = function(x, group, weights = NULL) {
  if (is.null(weights)) return(ave(x, group, FUN = mean))
  weights = as.numeric(weights)
  ave(seq_along(x), group, FUN = function(i) sum(x[i] * weights[i]) / sum(weights[i]))
}

demean_prediction_data = function(matrix, groups, weights = NULL, tol = 1e-10, max_iter = 1000L, slopes = NULL) {
  out = as.matrix(matrix)
  if (!length(slopes)) slopes = rep(list(NULL), length(groups))
  for (iteration in seq_len(max_iter)) {
    change = 0
    for (k in seq_along(groups)) {
      group = groups[[k]]
      if (length(slopes[[k]])) {
        basis = cbind(1, as.matrix(slopes[[k]]))
        means = matrix(0, nrow(out), ncol(out))
        for (rows in split(seq_len(nrow(out)), group)) {
          w = if (is.null(weights)) rep(1, length(rows)) else weights[rows]
          local = basis[rows, , drop = FALSE]
          local[, -1L] = sweep(
            local[, -1L, drop = FALSE], 2L,
            colSums(w * local[, -1L, drop = FALSE]) / sum(w)
          )
          scale = sqrt(colSums(w * local^2))
          scale[scale == 0] = 1
          local = sweep(local, 2L, scale, `/`)
          qr = qr(sqrt(w) * local, tol = 1e-11)
          coefficients = qr.coef(qr, sqrt(w) * out[rows, , drop = FALSE])
          coefficients[is.na(coefficients)] = 0
          means[rows, ] = local %*% coefficients
        }
      } else {
        means = vapply(seq_len(ncol(out)), function(j) weighted_group_mean(out[, j], group, weights), numeric(nrow(out)))
        means = matrix(means, nrow = nrow(out), ncol = ncol(out))
      }
      out = out - means
      change = max(change, abs(means))
    }
    if (is.finite(change) && change < tol) return(out)
  }
  stop("Prediction demeaning did not converge.", call. = FALSE)
}

# Length prefixes distinguish tuples containing delimiters, empty strings and
# numeric/integer64 keys without enumerating the Cartesian product of levels.
tuple_key = function(data, columns) {
  values = lapply(columns, function(column) {
    value = as.character(data[[column]])
    paste0(nchar(value, type = "bytes"), ":", value)
  })
  do.call(paste, c(values, sep = "|"))
}

prediction_factors = function(object, data) {
  for (fe in names(object[["composite_fes"]])) {
    if (fe %in% names(data)) next
    columns = object[["composite_fes"]][[fe]]
    if (!all(columns %in% names(data))) stop("Missing composite fixed-effect columns: ", paste(setdiff(columns, names(data)), collapse = ", "), call. = FALSE)
    key = tuple_key(data, columns)
    lookup = object[["composite_levels"]][[fe]]
    if (!is.null(lookup)) {
      data[[fe]] = unname(lookup[key])
      if (anyNA(data[[fe]])) stop("newdata contains unseen composite fixed-effect levels.", call. = FALSE)
    } else data[[fe]] = key
  }
  typed = unique(unlist(lapply(object[["i_specs"]], function(spec) c(spec[["first"]], spec[["second"]])), use.names = FALSE))
  for (variable in setdiff(names(object[["factor_levels"]]), typed)) {
    if (!variable %in% names(data)) next
    levels = as.character(object[["factor_levels"]][[variable]])
    values = as.character(data[[variable]])
    unseen = setdiff(unique(values[!is.na(values)]), levels)
    if (length(unseen)) stop("newdata contains unseen levels for `", variable, "`: ", paste(unseen, collapse = ", "), call. = FALSE)
    data[[variable]] = factor(values, levels = levels)
  }
  data
}

# Saved i() columns define the fitted coding without evaluating formula calls.
expand_newdata_i = function(fml, data, i_columns, i_specs = list()) {
  labels = attr(stats::terms(fml), "term.labels")
  if (NROW(i_columns) || length(i_specs)) {
    labels = labels[vapply(labels, function(term) is.null(parse_sql_i_term(term)), logical(1))]
    typed = list()
    for (spec in i_specs) {
      variable = spec[["first"]]
      if (!variable %in% names(data) || variable %in% names(typed)) next
      typed[[variable]] = spec[["labels"]][i_match_values(spec, data[[variable]])]
    }
    for (k in seq_len(NROW(i_columns))) {
      variable = i_columns[["var"]][k]
      name = i_columns[["name"]][k]
      if (!variable %in% names(data) && name %in% names(data)) next
      values = if (variable %in% names(typed)) typed[[variable]] else as.character(data[[variable]])
      value = as.numeric(values == i_columns[["level"]][k])
      by = i_columns[["by"]][k]
      if (!is.na(by)) value = value * data[[by]]
      data[[i_columns[["name"]][k]]] = value
    }
    if (NROW(i_columns)) labels = c(labels, vapply(i_columns[["name"]], function(name) paste(deparse(as.name(name), backtick = TRUE), collapse = ""), character(1)))
  }
  list(data = data, terms = labels)
}

build_newdata_design = function(fml, data, factor_levels = list(), i_columns = NULL, i_specs = list()) {
  expanded = expand_newdata_i(fml, data, i_columns, i_specs)
  firsts = vapply(i_specs, `[[`, character(1), "first")
  ordinary_levels = factor_levels[setdiff(names(factor_levels), firsts)]
  data = prediction_factors(list(factor_levels = ordinary_levels), expanded[["data"]])
  for (spec in i_specs[!duplicated(firsts)]) {
    first = spec[["first"]]
    if (!first %in% names(data) || !identical(spec[["source_type"]], "factor")) next
    data[[first]] = factor(spec[["labels"]][i_match_values(spec, expanded[["data"]][[first]])], levels = spec[["labels"]])
  }
  expanded[["data"]] = data
  labels = expanded[["terms"]]
  if (!length(labels)) labels = "1"
  fml = stats::reformulate(labels,
    intercept = attr(stats::terms(fml), "intercept") == 1L, env = baseenv()
  )
  prepared = explicit_factor_newdata(fml, expanded[["data"]], factor_levels)
  design = fitted_model_matrix(prepared[["formula"]], prepared[["data"]])
  for (term in names(prepared[["names"]])) {
    colnames(design) = gsub(term, prepared[["names"]][[term]], colnames(design), fixed = TRUE)
  }
  if (NROW(i_columns)) {
    quoted = vapply(i_columns[["name"]], function(name) paste(deparse(as.name(name), backtick = TRUE), collapse = ""), character(1))
    matched = match(colnames(design), quoted)
    present = !is.na(matched)
    colnames(design)[present] = i_columns[["name"]][matched[present]]
  }
  design
}

explicit_factor_newdata = function(fml, data, factor_levels) {
  sources = explicit_factor_sources(attr(stats::terms(fml), "term.labels"))
  names = character()
  for (term in names(sources)) {
    value = data[[sources[[term]]]]
    levels = explicit_prediction_levels(term, value, factor_levels)
    data[[term]] = if (length(levels) == 1L) rep(1, length(value)) else factor(value, levels = levels)
    key = paste(deparse(as.name(term), width.cutoff = 500L, backtick = TRUE), collapse = "")
    names[[key]] = if (length(levels) == 1L) paste0(term, "_constant") else term
  }
  replace = function(expr) {
    if (!is.call(expr)) return(expr)
    if (identical(expr[[1L]], as.name("factor"))) {
      return(as.name(paste(deparse(expr, width.cutoff = 500L), collapse = "")))
    }
    as.call(lapply(as.list(expr), replace))
  }
  list(formula = stats::as.formula(replace(as.call(fml)), env = baseenv()), data = data, names = names)
}

linear_fe_df = function(inputs, sizes) {
  if (!length(sizes)) return(0)
  slopes = inputs[["fe_slopes"]]
  if (!length(slopes)) slopes = rep(list(NULL), length(sizes))
  sum(sizes * (1L + lengths(slopes))) - length(sizes) + 1L
}

# Check rank after removing nuisance effects on the compressed cells, preserving
# slope semantics when an ordinary regressor is also in the absorbed FE span.
compressed_slope_rank = function(X, dat, inputs, columns, trend_sources = NULL) {
  if (!length(columns) || !length(unlist(inputs[["fe_slopes"]]))) return(character())
  values = as.matrix(X[, columns, drop = FALSE])
  groups = dat[inputs[["fe"]]]
  slopes = lapply(inputs[["fe_slopes"]], function(vars) {
    if (length(trend_sources)) vars = names(trend_sources)[match(vars, trend_sources)]
    dat[vars]
  })
  w = dat[["sum_w"]]
  residual = demean_prediction_data(values, groups, w,
    tol = inputs[["fe_tol"]], max_iter = max(1000L, inputs[["fe_max_iter"]]), slopes = slopes
  )
  S = crossprod(residual, w * residual)
  floor = sum(w) * (128 * .Machine[["double.eps"]] * apply(abs(values), 2L, max))^2
  absorbed = which(diag(S) <= floor)
  S[absorbed, ] = 0
  S[, absorbed] = 0
  detect_collinearity(S, matrix(0, length(columns), 1L), verbose = inputs[["verbose"]])[["drop_names"]]
}

compressed_trend_data = function(data, sources, centers = NULL) {
  for (alias in names(sources)) {
    value = data[[sources[[alias]]]]
    if (is.factor(value)) {
      keys = attr(value, "dbreg_key_values")
      value = if (is.null(keys)) as.character(value) else keys[as.integer(value)]
    }
    value = as.numeric(value)
    data[[alias]] = if (is.null(centers)) value else value - centers[[alias]]
  }
  data
}

# A column identity is compiler metadata, never a public coefficient label.
empty_linear_columns = function() {
  columns = data.frame(id = character(), term = integer(), basis = integer(),
    alias = character(), label = character(), role = character(), stringsAsFactors = FALSE)
  columns[["recipe"]] = list()
  columns[["cell_recipe"]] = list()
  columns
}

linear_intercept_column = function() {
  column = data.frame(id = "intercept", term = 0L, basis = 1L,
    alias = "(Intercept)", label = "(Intercept)", role = "intercept", stringsAsFactors = FALSE)
  column[["recipe"]] = list(list(kind = "constant"))
  column[["cell_recipe"]] = list(NULL)
  column
}

native_linear_layout = function(design, fitted, intercept) {
  columns = design[["columns"]]
  if (is.null(columns)) columns = empty_linear_columns()
  if (intercept) columns = rbind(linear_intercept_column(), columns)
  positions = match(fitted, columns[["alias"]])
  if (anyNA(positions) || anyDuplicated(columns[["id"]])) {
    stop("Fitted native coordinates do not match their compiled columns.", call. = FALSE)
  }
  ids = columns[["id"]][positions]
  list(version = 1L, columns = columns, fitted_ids = ids,
    slope_ids = ids[columns[["role"]][positions] == "slope"])
}

# Capture the expansion actually emitted by the compressed matrix producer.
# Its term assignment disambiguates a factor dummy from a same-labelled source.
compressed_linear_layout = function(X, formula, data, design, original_formula, inputs) {
  assignment = attr(X, "assign")
  terms = attr(stats::terms(formula), "term.labels")
  original = attr(stats::terms(stats::formula(original_formula, lhs = 0, rhs = 1)), "term.labels")
  compiled = design[["columns"]]
  columns = empty_linear_columns()
  for (j in seq_len(ncol(X))) {
    term_id = assignment[[j]]
    if (!term_id) {
      column = linear_intercept_column()
    } else {
      variable = terms[[term_id]]
      emitted = if (!is.null(compiled)) match(variable, compiled[["alias"]]) else NA_integer_
      if (!is.na(emitted) && !identical(compiled[["recipe"]][[emitted]][["kind"]], "raw")) {
        column = compiled[emitted, , drop = FALSE]
      } else {
        basis = match(j, which(assignment == term_id))
        value = data[[variable]]
        categorical = is.factor(value) || is.character(value) || is.logical(value)
        if (categorical) {
          levels = if (is.factor(value)) levels(value) else levels(factor(value))
          width = sum(assignment == term_id)
          selected = if (width == length(levels)) seq_along(levels) else seq.int(2L, length(levels))
          if (length(selected) != width) stop("Compressed factor coding does not match its emitted basis.", call. = FALSE)
          recipe = list(kind = "indicator", variable = variable, value = levels[selected[[basis]]])
        } else recipe = list(kind = "numeric", variable = variable)
        source_term = match(variable, original)
        if (is.na(source_term)) stop("Compressed term has no fitted source identity.", call. = FALSE)
        column = data.frame(id = paste0("t", source_term, "b", basis), term = source_term,
          basis = basis, alias = colnames(X)[[j]], label = colnames(X)[[j]], role = "slope", stringsAsFactors = FALSE)
        column[["recipe"]] = list(recipe)
      }
    }
    if (!term_id) {
      cell_recipe = list(kind = "constant")
    } else {
      position = match(variable, names(data))
      if (is.na(position)) stop("Compressed cell term has no prepared column.", call. = FALSE)
      cell_recipe = if (!is.na(emitted) && !identical(compiled[["recipe"]][[emitted]][["kind"]], "raw")) {
        list(kind = "numeric", variable = position)
      } else {
        recipe = column[["recipe"]][[1L]]
        recipe[["variable"]] = position
        recipe
      }
    }
    column[["cell_recipe"]] = list(cell_recipe)
    columns = rbind(columns, column)
  }
  aliases = columns[["alias"]]
  positions = which(columns[["role"]] != "intercept")
  aliases[positions] = safe_identifiers(paste0("dbreg_linear_column", seq_along(positions), "__"),
    taken = c(names(data), names(inputs[["name_map"]])))
  display = restore_names(standardize_coef_names(colnames(X)),
    coefficient_name_map(inputs[["name_map"]], inputs[["source_map"]]))
  additions = stats::setNames(display[positions], aliases[positions])
  inputs[["name_map"]] = c(inputs[["name_map"]], additions)
  columns[["alias"]] = aliases
  list(columns = columns, aliases = aliases)
}

restore_linear_recipe = function(recipe, map) {
  if (!is.null(recipe[["variable"]])) recipe[["variable"]] = restore_names(recipe[["variable"]], map)
  if (!is.null(recipe[["second"]])) recipe[["second"]] = restore_names(recipe[["second"]], map)
  if (!is.null(recipe[["term"]])) recipe[["term"]] = restore_term_names(recipe[["term"]], map)
  if (!is.null(recipe[["left"]])) recipe[["left"]] = restore_linear_recipe(recipe[["left"]], map)
  if (!is.null(recipe[["right"]])) recipe[["right"]] = restore_linear_recipe(recipe[["right"]], map)
  recipe
}

linear_recipe_value = function(recipe, data, i_specs) {
  kind = recipe[["kind"]]
  if (kind == "constant") return(rep(1, nrow(data)))
  if (kind == "product") return(linear_recipe_value(recipe[["left"]], data, i_specs) *
    linear_recipe_value(recipe[["right"]], data, i_specs))
  value = data[[recipe[["variable"]]]]
  if (kind == "numeric") return(as.numeric(value))
  if (kind == "power") return(as.numeric(value^recipe[["exponent"]]))
  if (kind == "indicator") {
    expected = recipe[["value"]]
    if (inherits(expected, "integer64") || !is.numeric(expected)) {
      return(as.numeric(as.character(value) == as.character(expected)))
    }
    if (is.factor(value) || is.character(value)) return(as.numeric(as.character(value) == as.character(expected)))
    return(as.numeric(value == expected))
  }
  if (kind == "i") {
    matching = which(vapply(i_specs, function(spec) identical(spec[["term"]], recipe[["term"]]), logical(1)))
    if (length(matching) != 1L) stop("Missing fitted indicator dictionary.", call. = FALSE)
    out = as.numeric(i_match_values(i_specs[[matching]], value, strict = TRUE) == recipe[["level"]])
    if (!is.null(recipe[["second"]])) out = out * as.numeric(data[[recipe[["second"]]]])
    return(out)
  }
  stop("Unsupported internal native column recipe.", call. = FALSE)
}

native_linear_matrix = function(object, data, stored = FALSE) {
  layout = object[["linear_design"]]
  columns = layout[["columns"]]
  positions = match(layout[["fitted_ids"]], columns[["id"]])
  if (!stored) {
    rhs = stats::formula(object[["fml"]], lhs = 0, rhs = 1)
    forced = explicit_factor_sources(attr(stats::terms(rhs), "term.labels"))
    for (term in names(forced)) explicit_prediction_levels(term, data[[forced[[term]]]], object[["factor_levels"]])
    validate_native_numeric_types(object, data, columns[["recipe"]][positions])
  }
  recipes = if (stored) columns[["cell_recipe"]] else columns[["recipe"]]
  values = lapply(recipes[positions], linear_recipe_value, data = data,
    i_specs = object[["i_specs"]])
  out = if (length(values)) do.call(cbind, values) else matrix(0, nrow(data), 0L)
  colnames(out) = columns[["label"]][positions]
  out
}

restore_native_linear_layout = function(layout, source_map, name_map, fitted_labels) {
  layout[["columns"]][["recipe"]] = lapply(layout[["columns"]][["recipe"]],
    restore_linear_recipe, map = source_map)
  layout[["columns"]][["label"]] = standardize_coef_names(
    restore_names(layout[["columns"]][["label"]], name_map))
  positions = match(layout[["fitted_ids"]], layout[["columns"]][["id"]])
  layout[["columns"]][["label"]][positions] = fitted_labels
  layout
}

# The raw source determines forced-factor validation even when a column with
# the generated factor expression's name is present in arbitrary newdata.
explicit_prediction_levels = function(term, value, factor_levels) {
  fitted = factor_levels[[term]]
  if (is.null(fitted)) fitted = levels(factor(value))
  unseen = setdiff(unique(as.character(value[!is.na(value)])), as.character(fitted))
  if (length(unseen)) stop("newdata contains unseen levels for `", term, "`.", call. = FALSE)
  fitted
}

# Preserve the existing R expression/type contract for exceptional inputs.
# Numeric vectors take the recipe fast path; categorical reinterpretations do
# not become factor codes merely because an output column has the same label.
validate_native_numeric_types = function(object, data, recipes) {
  inspect = function(recipe) {
    kind = recipe[["kind"]]
    if (kind == "product") return(c(inspect(recipe[["left"]]), inspect(recipe[["right"]])))
    if (kind %in% c("numeric", "power")) {
      value = data[[recipe[["variable"]]]]
      regular = is.numeric(value) && !is.factor(value) && !is.complex(value) && is.null(dim(value))
      return(list(list(exceptional = !regular,
        invalid = is.factor(value) || is.character(value) || is.complex(value) ||
          kind == "numeric" && is.logical(value))))
    }
    if (kind == "i" && !is.null(recipe[["second"]])) {
      value = data[[recipe[["second"]]]]
      regular = (is.numeric(value) || is.logical(value)) && !is.factor(value) && !is.complex(value)
      return(list(list(exceptional = !regular, invalid = !regular)))
    }
    list()
  }
  inputs = unlist(lapply(recipes, inspect), recursive = FALSE)
  if (!length(inputs) || !any(vapply(inputs, `[[`, logical(1), "exceptional"))) return(invisible(NULL))
  # This adapter validates/evaluates the existing formula semantics only. Its
  # columns are not used to align coefficients or replace fitted coordinates.
  rhs = stats::formula(object[["fml"]], lhs = 0, rhs = 1)
  legacy = build_newdata_design(rhs, data, object[["factor_levels"]],
    object[["i_columns"]], object[["i_specs"]])
  missing = setdiff(names(coef.dbreg(object, fe = TRUE)), colnames(legacy))
  if (length(missing)) stop("predict() could not rebuild the design for coefficient(s): ",
    paste(missing, collapse = ", "), call. = FALSE)
  if (any(vapply(inputs, `[[`, logical(1), "invalid"))) {
    stop("newdata changes a fitted numeric term to an incompatible input type.", call. = FALSE)
  }
  invisible(NULL)
}
