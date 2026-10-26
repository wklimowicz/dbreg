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
  result[["fml"]] = freeze_i_selectors(fml)
  design = result[["compressed_design"]]
  if (!is.null(design)) {
    environment(design[["formula"]]) = baseenv()
    design[["source_names"]] = setNames(design[["source_names"]], restore_names(design[["source_names"]], inputs[["name_map"]]))
    design[["columns"]] = setNames(rownames(result[["coeftable"]]), design[["columns"]])
    if (length(design[["trend_sources"]])) {
      design[["prediction_trends"]] = stats::setNames(
        restore_names(unname(design[["trend_sources"]]), inputs[["source_map"]]), names(design[["trend_sources"]]))
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
          local[, -1L] = sweep(local[, -1L, drop = FALSE], 2L,
            colSums(w * local[, -1L, drop = FALSE]) / sum(w))
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
    intercept = attr(stats::terms(fml), "intercept") == 1L, env = baseenv())
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
    levels = factor_levels[[term]]
    if (is.null(levels)) levels = levels(factor(value))
    unseen = setdiff(unique(as.character(value[!is.na(value)])), as.character(levels))
    if (length(unseen)) stop("newdata contains unseen levels for `", term, "`.", call. = FALSE)
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
    tol = inputs[["fe_tol"]], max_iter = max(1000L, inputs[["fe_max_iter"]]), slopes = slopes)
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
