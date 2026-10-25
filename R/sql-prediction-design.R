# Recompile from fitted coding, never category discovery in prediction data.
# The same term expanders serve estimation and this SQL prediction path.
sql_fitted_prediction_design = function(object, conn, factor_labels = FALSE) {
  glm = inherits(object, "dbglm")
  state = object[["fitted_state"]]
  if (glm) {
    specification = object[["model_spec"]]
    rhs = specification[["rhs"]]
    variables = specification[["source_vars"]]
    if (any(!is_safe_identifier(variables))) {
      stop("SQL prediction requires SQL-safe fitted GLM regressor names.", call. = FALSE)
    }
    aliases = variables
    sources = stats::setNames(variables, aliases)
    for (name in intersect(names(object[["factor_sources"]]), aliases)) sources[name] = object[["factor_sources"]][[name]]
    levels = specification[["factor_levels"]]
    classes = specification[["factor_classes"]]
    if (factor_labels) classes = stats::setNames(rep("character", length(levels)), names(levels))
    info = list(i_specs = specification[["i_specs"]],types = as.list(stats::setNames(ifelse(aliases %in% names(levels), "factor", "numeric"), aliases)),
      levels = levels, factor_classes = classes, level_sql = if (factor_labels)
        lapply(levels, function(value) sql_value_literals(conn, as.character(value))) else sql_factor_literals(conn, levels, classes))
    if (!factor_labels) for (name in names(specification[["i_specs"]])) {
      values = stats::setNames(list(specification[["i_specs"]][[name]][["values"]]),name)
      info[["level_sql"]][[name]] = sql_factor_literals(conn,values,classes)[[name]]
    }
    info[["classes"]] = specification[["source_classes"]]
    info[["numeric_sql"]] = stats::setNames(sql_cast_double(aliases, conn), aliases)
    design = sql_model_matrix_from_info(attr(stats::terms(rhs), "term.labels"), info,
      intercept = attr(stats::terms(rhs), "intercept") == 1L || length(object[["fe"]]) > 0L)
    names = design[["coef_names"]]
    expressions = design[["select_exprs"]]
    if (!identical(names, specification[["columns"]][["name"]])) {
      stop("Could not reconstruct the fitted GLM SQL design.", call. = FALSE)
    }
    intercept = isTRUE(object[["has_intercept"]])
  } else {
    rhs = stats::formula(object[["fml"]], lhs = 0, rhs = 1)
    variables = all.vars(rhs)
    aliases = safe_identifiers(variables, fold_case = folds_identifier_case(conn))
    sources = stats::setNames(variables, aliases)
    renamed = aliases != variables
    source_map = stats::setNames(variables[renamed], aliases[renamed])
    rhs = rename_expr_vars(rhs, stats::setNames(aliases, variables))
    source_aliases = stats::setNames(aliases, variables)
    raw_rhs = stats::formula(object[["fml"]], lhs = 0, rhs = 1)
    forced = explicit_factor_sources(attr(stats::terms(raw_rhs), "term.labels"))
    relevant = c(variables, names(forced))
    levels = object[["factor_levels"]][intersect(names(object[["factor_levels"]]), relevant)]
    if (length(levels)) names(levels) = restore_term_names(names(levels), source_aliases)
    classes = object[["factor_classes"]][intersect(names(object[["factor_classes"]]), relevant)]
    if (length(classes)) names(classes) = restore_term_names(names(classes), source_aliases)
    forced_safe = explicit_factor_sources(attr(stats::terms(rhs), "term.labels"))
    for (term in names(forced_safe)) sources[[term]] = variables[match(forced_safe[[term]], aliases)]
    if (length(forced)) {
      unused = setdiff(unname(forced), model_vars(raw_rhs, skip_factors = TRUE))
      unused = aliases[match(unused, variables)]
      sources = sources[!names(sources) %in% unused]
      levels = levels[!names(levels) %in% unused]
      classes = classes[!names(classes) %in% unused]
    }
    if (factor_labels) classes = stats::setNames(rep("character", length(levels)), names(levels))
    # Scalar FE terms affect the native interaction reference rule; composite
    # FE ranks are distinct from their component variables in that rule.
    scalar_fe = unlist(lapply(state[["specification"]][["dictionaries"]], function(dictionary) {
      if (length(dictionary[["columns"]]) == 1L) dictionary[["columns"]] else character()
    }), use.names = FALSE)
    fe = aliases[match(intersect(scalar_fe, variables), variables)]
    info = list(types = as.list(stats::setNames(ifelse(aliases %in% names(levels), "factor", "numeric"), aliases)),
      levels = levels, factor_classes = classes, level_sql = if (factor_labels)
        lapply(levels, function(value) sql_value_literals(conn, as.character(value))) else sql_factor_literals(conn, levels, classes))
    info[["classes"]] = classes
    info[["factor_sources"]] = as.list(stats::setNames(as.character(DBI::dbQuoteIdentifier(conn, as.character(names(forced_safe)))), names(forced_safe)))
    info[["numeric_sql"]] = stats::setNames(sql_cast_double(aliases, conn), aliases)
    intercept = attr(stats::terms(rhs), "intercept") == 1L
    design = sql_model_matrix_from_info(attr(stats::terms(rhs), "term.labels"), info,
      fe_vars = fe, intercept = intercept, fold_case = folds_identifier_case(conn))
    inputs = new.env(parent = emptyenv())
    inputs[["source_map"]] = source_map
    inputs[["name_map"]] = source_map
    register_factor_display(inputs, design[["factor_levels"]], design[["level_aliases"]],
      design[["col_names"]], design[["aliases"]], i_columns = design[["i_columns"]])
    map = inputs[["name_map"]]
    variable = names(map) %in% names(source_map)
    map[variable] = coef_var_name(map[variable])
    names = standardize_coef_names(restore_names(design[["aliases"]], map))
    expressions = design[["select_exprs"]]
    intercept = intercept && !length(object[["fe"]])
  }
  if (intercept) {
    expressions = c("1.0", expressions)
    names = c("(Intercept)", names)
  }
  list(expressions = expressions, names = names, sources = sources,
    factor_levels = levels, factor_classes = classes,
    level_sql = info[["level_sql"]],
    factor_labels = factor_labels,
    numeric_sources = setdiff(aliases, names(levels)))
}
