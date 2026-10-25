# A single GLM design specification shared by SQL fitting and R prediction.
build_model_spec = function(formula, conn, table, fe_vars = character(),
                             factor_refs = NULL, reserved_names = character(), absorbed = FALSE) {
  design = sql_model_matrix(
    formula = formula,
    conn = conn,
    table = table,
    expand = "all",
    sep = "_x_",
    fe_vars = fe_vars,
    factor_refs = factor_refs,
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
  aliases = model_spec[["columns"]][["alias"]]
  if (all(aliases %in% names(newdata))) {
    alias_data = model_matrix_data(~ ., newdata[, aliases, drop = FALSE])
    matrix = Matrix::Matrix(as.matrix(alias_data), sparse = TRUE)
    colnames(matrix) = aliases
    return(matrix)
  }
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

  rhs = if(length(model_spec[["i_specs"]])) model_spec[["ordinary_rhs"]] else model_spec[["rhs"]]
  design_data = model_matrix_data(rhs, newdata)
  matrix = fitted_model_matrix(rhs, design_data)
  if ("(Intercept)" %in% colnames(matrix)) {
    matrix = matrix[, colnames(matrix) != "(Intercept)", drop = FALSE]
  }
  for (spec in model_spec[["i_specs"]]) {
    matrix = cbind(matrix,evaluate_i_spec(spec,newdata[[spec[["first"]]]],
      if(!is.null(spec[["second"]])) newdata[[spec[["second"]]]]))
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
