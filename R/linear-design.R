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
  from_statement = paste0("FROM (SELECT dbreg_source.*, ", paste(expressions, collapse = ", "),
    " FROM (SELECT * ", from_statement, ") dbreg_source) dbreg_composite")
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

# Only compressed fits retain explicit FE coefficients. Their tuple-to-rank
# mapping comes from cells already collected for estimation, with no extra scan.
finalize_linear_design = function(result, inputs) {
  levels = inputs[["factor_levels"]]
  if (length(levels)) names(levels) = restore_names(names(levels), inputs[["source_map"]])
  result[["factor_levels"]] = levels
  composites = inputs[["composite_fes"]]
  result[["composite_fes"]] = lapply(composites, function(vars) restore_names(vars, inputs[["source_map"]]))
  result[["composite_levels"]] = inputs[["composite_levels"]]
  result[["fe_labels"]] = inputs[["fe_labels"]]
  fml = if (!is.null(inputs[["original_fml"]])) inputs[["original_fml"]] else result[["fml"]]
  # The supported design uses columns and built-in transformations only.
  # Retaining the caller/strategy frame can capture the entire original panel
  # (and a live connection) through an otherwise compact formula object.
  environment(fml) = baseenv()
  result[["fml"]] = fml
  design = result[["compressed_design"]]
  if (!is.null(design)) {
    environment(design[["formula"]]) = baseenv()
    design[["source_names"]] = setNames(design[["source_names"]], restore_names(design[["source_names"]], inputs[["name_map"]]))
    design[["columns"]] = setNames(rownames(result[["coeftable"]]), design[["columns"]])
    result[["compressed_design"]] = design
  }
  result
}

weighted_group_mean = function(x, group, weights = NULL) {
  if (is.null(weights)) return(ave(x, group, FUN = mean))
  weights = as.numeric(weights)
  ave(seq_along(x), group, FUN = function(i) sum(x[i] * weights[i]) / sum(weights[i]))
}

demean_prediction_data = function(matrix, groups, weights = NULL, tol = 1e-10, max_iter = 1000L) {
  out = as.matrix(matrix)
  for (iteration in seq_len(max_iter)) {
    change = 0
    for (group in groups) {
      means = vapply(seq_len(ncol(out)), function(j) weighted_group_mean(out[, j], group, weights), numeric(nrow(out)))
      means = matrix(means, nrow = nrow(out), ncol = ncol(out))
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
  for (variable in names(object[["factor_levels"]])) {
    if (!variable %in% names(data)) next
    levels = as.character(object[["factor_levels"]][[variable]])
    values = as.character(data[[variable]])
    unseen = setdiff(unique(values[!is.na(values)]), levels)
    if (length(unseen)) stop("newdata contains unseen levels for `", variable, "`: ", paste(unseen, collapse = ", "), call. = FALSE)
    data[[variable]] = factor(values, levels = levels)
  }
  data
}
