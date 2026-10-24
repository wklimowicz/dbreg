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
  structure(codes, levels = labels, class = "factor", dbreg_key_values = values,
    dbreg_key_type = if (is.integer(x)) "integer" else "numeric")
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
      types = DBI::dbGetQuery(conn, paste0("SELECT ", paste0("typeof(", q(risky), ") AS t", seq_along(risky), collapse = ","),
        " ", from_statement, " LIMIT 1"))
      integers = risky[vapply(types, function(type) grepl("^(u?(bigint|hugeint)|decimal\\([0-9]+,0\\))$", tolower(type[[1]])), logical(1))]
      if (length(integers)) {
        cast = if (backend[["name"]] == "duckdb") "VARCHAR" else "STRING"
        fields = c(q(setdiff(names(data), integers)), paste0("CAST(", q(integers), " AS ", cast, ") AS ", q(integers)))
        sql = paste0("SELECT ", paste(fields, collapse = ","), " FROM (", sql, ") dbreg_exact_cells")
        data = DBI::dbGetQuery(conn, sql)
      }
    }
  }
  # Keep typed identifiers intact, but cross the same numeric boundary used by
  # cross_moments() for sufficient statistics. integer64 arithmetic otherwise
  # truncates fractional residuals, even when cell means/coefficients are doubles.
  moments = intersect(c("sum_w", "sum_wy", "sum_wy_sq", "sum_w2", "sum_w2y", "sum_w2y_sq",
    "sum_w_gc", "sum_wy_gc", "mean_Y", "wts"), names(data))
  for (name in moments) data[[name]] = as.numeric(data[[name]])
  list(data = data, sql = sql, integers = integers)
}

# Only compressed fits retain explicit FE coefficients. Their tuple-to-rank
# mapping comes from cells already collected for estimation, with no extra scan.
finalize_linear_design = function(result, inputs) {
  result[["offset"]] = inputs[["offset"]]
  result[["fixef_rm"]] = inputs[["fixef_rm"]]
  result[["pruning_rounds"]] = if (is.null(inputs[["pruning"]])) 0L else inputs[["pruning"]][["rounds"]]
  result[["exclusions"]] = if (is.null(inputs[["pruning"]])) c(singleton = 0, infinite_coef = 0) else inputs[["pruning"]][["counts"]]
  levels = inputs[["factor_levels"]]
  if (length(levels)) names(levels) = restore_names(names(levels), inputs[["source_map"]])
  result[["factor_levels"]] = levels
  classes = inputs[["factor_classes"]]
  if (length(classes)) names(classes) = restore_names(names(classes), inputs[["source_map"]])
  result[["factor_classes"]] = classes
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
