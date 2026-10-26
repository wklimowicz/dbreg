# Formula and model-matrix helpers shared by OLS and GLM.

# Keep the fast matrix product for ordinary inputs. At extreme magnitudes,
# fused BLAS operations can hide Inf-Inf as Inf (and depend on batch shape).
# Re-evaluate potentially overflowing rows with the same staged arithmetic
# as the SQL predictor, using row vectors rather than another dense matrix.
prediction_matrix_product = function(x, beta) {
  value = as.vector(x %*% beta)
  if (!length(beta) || !nrow(x)) return(value)
  bound = .Machine[["double.xmax"]]/max(1,max(abs(beta)))/max(1,ncol(x))/2
  if (isTRUE(Matrix::norm(x,"M")<=bound)) return(value)
  risky = !is.finite(value)
  for (j in seq_along(beta)) {
    limit = .Machine[["double.xmax"]]/max(1,abs(beta[[j]]))/max(1,ncol(x))/2
    column = x[,j,drop=FALSE]
    if (isTRUE(Matrix::norm(column,"M")<=limit)) next
    column = as.numeric(column)
    risky = risky | !is.finite(column) | abs(column)>limit
  }
  if (any(risky)) {
    total = rep(0,sum(risky))
    for (j in seq_along(beta)) total = total+as.numeric(x[risky,j])*beta[[j]]
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
    error = function(e) stats::model.matrix(formula, data = data, contrasts.arg = contrasts))
}

