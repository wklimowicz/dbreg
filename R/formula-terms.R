# Formula and model-matrix helpers shared by OLS and GLM.
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
  tryCatch(Matrix::sparse.model.matrix(formula, data = data, contrasts.arg = contrasts),
    error = function(e) stats::model.matrix(formula, data = data, contrasts.arg = contrasts))
}

