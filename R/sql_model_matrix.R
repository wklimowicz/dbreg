#' Construct SQL expressions for a design matrix
#'
#' Expands formula terms into SQL SELECT expressions, handling factor one-hot
#' encoding and interaction terms. Only the right-hand side of the formula is
#' processed; the response variable (LHS) is ignored.
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
#' @param factor_refs Optional named list of categorical reference levels.
#' @return List with:
#'   - `select_exprs`: character vector of SQL expressions
#'   - `col_names`: corresponding column names
#'   - `factor_levels`: list of factor levels by variable (for reference)
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
#' dbDisconnect(con)
sql_model_matrix = function(
  formula,
  conn,
  table,
  expand = c("all", "interactions"),
  sep = "_x_",
  fe_vars = character(),
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
    return(list(select_exprs = character(), col_names = character(), factor_levels = list()))
  }
  
  # Get column types and factor levels from database
  col_info = get_column_info(conn, table, term_labels, factor_refs)
  
  result = list(select_exprs = character(), col_names = character())
  
  for (term in term_labels) {
    expanded = expand_term(term, col_info, expand, sep, all_terms = unique(c(term_labels, fe_vars)), fe_vars = fe_vars)
    result[["select_exprs"]] = c(result[["select_exprs"]], vapply(expanded, `[[`, character(1), "sql"))
    result[["col_names"]] = c(result[["col_names"]], vapply(expanded, `[[`, character(1), "name"))
  }
  
  result[["factor_levels"]] = col_info[["levels"]]
  result
}

#' Get column types and factor levels from database
#' @keywords internal
get_column_info = function(conn, table, term_labels, factor_refs = NULL) {
  # Extract unique variable names from all terms
  vars = unique(unlist(lapply(term_labels, function(term) {
    all.vars(parse(text = term, keep.source = FALSE)[[1]])
  })))
  
  # Build FROM clause
  from_sql = if (grepl("^FROM", table, ignore.case = TRUE)) table else paste("FROM", table)
  
  # Query to get types and sample values
  # Use LIMIT 0 to get schema, then separate queries for factor levels
  schema_query = glue("SELECT {paste(vars, collapse = ', ')} {from_sql} LIMIT 0")
  schema = dbGetQuery(conn, schema_query)
  
  types = list()
  levels = list()
  
  for (v in vars) {
    col_class = class(schema[[v]])[1]
    # Treat character, factor, logical as categorical
    if (col_class %in% c("character", "factor", "logical") || v %in% names(factor_refs)) {
      types[[v]] = "factor"
      # Get distinct levels
      lvl_query = glue(
        "SELECT DISTINCT {v} FROM (SELECT * {from_sql}) AS dbreg_levels ",
        "WHERE {v} IS NOT NULL ORDER BY {v}"
      )
      levels[[v]] = dbGetQuery(conn, lvl_query)[[1]]
      if (v %in% names(factor_refs)) {
        ref_idx = match(as.character(factor_refs[[v]]), as.character(levels[[v]]))
        if (is.na(ref_idx)) {
          stop(
            "Reference level `", factor_refs[[v]], "` was not found in regressor `", v, "`.",
            call. = FALSE
          )
        }
        levels[[v]] = levels[[v]][c(ref_idx, setdiff(seq_along(levels[[v]]), ref_idx))]
      }
    } else {
      types[[v]] = "numeric"
    }
  }
  
  level_sql = lapply(levels, function(x) {
    vapply(x, function(value) as.character(DBI::dbQuoteLiteral(conn, value)), character(1))
  })
  list(types = types, levels = levels, level_sql = level_sql)
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

#' Expand a single term into SQL expressions
#' @keywords internal
expand_term = function(term, col_info, expand, sep = "_x_", all_terms = character(), fe_vars = character()) {
  vars = split_interaction_term(term)
  is_interaction = length(vars) > 1
  
  if (!is_interaction && expand == "interactions" && is.null(parse_power_term(term))) {
    # Main effect, only expanding interactions: return as-is
    return(list(list(sql = term, name = term)))
  }
  
  # Expand each variable
  expansions = lapply(vars, function(v) {
    drop_ref = TRUE
    if (is_interaction && identical(col_info[["types"]][[v]], "factor")) {
      other_term = paste(sort(setdiff(vars, v)), collapse = ":")
      drop_ref = v %in% fe_vars || other_term %in% all_terms
    }
    expand_variable(v, col_info, expand, is_interaction, drop_ref)
  })
  
  # Cross product for interactions
  cross_product(expansions, sep)
}

#' Expand a single variable into SQL expression(s)
#' @keywords internal
expand_variable = function(var, col_info, expand, in_interaction, drop_ref = TRUE) {
  power = parse_power_term(var)
  source_var = if (is.null(power)) var else power[["variable"]]
  is_factor = is.null(power) && identical(col_info[["types"]][[source_var]], "factor")
  
  # Only expand factors if: expand="all", OR variable is part of interaction
  if (is_factor && (expand == "all" || in_interaction)) {
    lvls = col_info[["levels"]][[source_var]]
    lvl_sql = col_info[["level_sql"]][[source_var]]
    if (length(lvls) < 2) {
      return(list(list(sql = "1", name = paste0(var, "_constant"))))
    }
    # Drop first level (reference) only if appropriate
    if (drop_ref) {
      lvls = lvls[-1]
      lvl_sql = lvl_sql[-1]
    }
    lapply(seq_along(lvls), function(i) {
      lvl = lvls[[i]]
      list(
        sql = glue("CASE WHEN {source_var} = {lvl_sql[[i]]} THEN 1.0 ELSE 0.0 END"),
        name = paste0(source_var, lvl)
      )
    })
  } else if (!is.null(power)) {
    list(list(sql = power[["sql"]], name = power[["name"]]))
  } else {
    # Numeric or unexpanded factor: pass through
    list(list(sql = var, name = var))
  }
}

#' Cartesian product of term expansions
#' @keywords internal
cross_product = function(expansions, sep = "_x_") {
  if (length(expansions) == 1) return(expansions[[1]])
  
  # Recursive: cross first with rest
  rest = cross_product(expansions[-1], sep)
  
  result = list()
  for (e1 in expansions[[1]]) {
    for (e2 in rest) {
      e1_sql = e1[["sql"]]
      e2_sql = e2[["sql"]]
      result = c(result, list(list(
        sql = glue("({e1_sql}) * ({e2_sql})"),
        name = paste0(e1[["name"]], sep, e2[["name"]])
      )))
    }
  }
  result
}
