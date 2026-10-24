# Portable prediction state contains values only. Observation-sized relations
# remain fit-owned and are released by the usual fitting resource manager.
prediction_state_limit = function() {
  limit = getOption("dbreg.prediction_max_entries", 2e6)
  dbglm_validate_moment_limits(limit, limit)
  limit
}

prediction_state_overflow = function() {
  stop(structure(list(message = paste0("Portable prediction state exceeds its transfer or matrix limit; ",
    "database-backed recovery is unavailable for this design or budget. Use retain='none' for a compact fit."),
    call = NULL), class = c("dbreg_prediction_overflow", "error", "condition")))
}

prediction_key_part = function(value, type) {
  if (type == "integer") {
    if (is.numeric(value) && !inherits(value, "integer64")) {
      bad = !is.na(value) & (!is.finite(value) | value != floor(value) | abs(value) > 2^53)
      if (any(bad)) stop("Integer FE keys must be exact; use integer64 or decimal strings for large keys.", call. = FALSE)
      value = format(value, scientific = FALSE, trim = TRUE, digits = 22)
    } else value = as.character(value)
    valid = !is.na(value) & grepl("^[+-]?[0-9]+$", value)
    value[!valid] = NA_character_
    value = sub("^\\+", "", value)
    value = sub("^(-?)0+([0-9])", "\\1\\2", value)
    value[value %in% c("-0", "+0")] = "0"
  } else if (type %in% c("numeric", "float")) {
    if (!is.numeric(value) || inherits(value, "integer64")) {
      stop("Numeric FE keys require numeric newdata columns.", call. = FALSE)
    }
    if (type == "float") value = readBin(writeBin(as.double(value), raw(), size = 4L), "double", n = length(value), size = 4L)
    value[value == 0 & !is.na(value)] = 0
    value = ifelse(is.finite(value), sprintf("%.17g", value), NA_character_)
  } else if (type == "logical") {
    if (!is.logical(value)) stop("Logical FE keys require logical newdata columns.", call. = FALSE)
    value = ifelse(is.na(value), NA_character_, ifelse(value, "TRUE", "FALSE"))
  } else value = as.character(value)
  ifelse(is.na(value), NA_character_, paste0(nchar(enc2utf8(value), type = "bytes"), ":", enc2utf8(value)))
}

prediction_tuple_key = function(data, types) {
  parts = Map(prediction_key_part, data, types)
  key = do.call(paste0, unname(parts))
  missing = Reduce(`|`, lapply(parts, is.na))
  key[missing] = NA_character_
  key
}

prediction_key_fields = function(conn, components) {
  components = lapply(components, unique)
  columns = lapply(seq_along(components), function(i) paste0("__dbglm_key_fe", i, "_", seq_along(components[[i]])))
  fields = unlist(lapply(seq_along(components), function(i) {
    paste(as.character(DBI::dbQuoteIdentifier(conn, components[[i]])), "AS", columns[[i]])
  }), use.names = FALSE)
  list(fields = fields, columns = columns)
}

prediction_key_metadata = function(conn, retained, components) {
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  sql_types = dbglm_get_query(conn, paste0("SELECT ", paste0("typeof(", q(components), ") AS t", seq_along(components), collapse = ","),
    " FROM ", retained, " LIMIT 1"))
  types = vapply(sql_types, function(type) {
    type = tolower(type[[1]])
    if (grepl("^(tinyint|smallint|int|integer|bigint|hugeint|u(tinyint|smallint|integer|bigint|hugeint)|decimal\\([0-9]+,0\\))$", type)) "integer" else
    if (type %in% c("double", "double precision")) "numeric" else
    if (type %in% c("float", "real")) "float" else
    if (type %in% c("boolean", "bool")) "logical" else
    if (grepl("^(varchar|char|string|enum)", type)) "character" else
    if (type == "date") "date" else stop("Portable FE prediction does not support key SQL type '", type, "'.", call. = FALSE)
  }, character(1))
  list(types = unname(types), sql_types = vapply(sql_types, function(type) type[[1]], character(1)))
}

prediction_fe_dictionary = function(conn, retained, fe, components, limit, public_components = components) {
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  types = prediction_key_metadata(conn, retained, components)[["types"]]
  fields = vapply(seq_along(components), function(i) {
    value = q(components[i])
    if (types[i] == "integer") value = paste0("CAST(", value, " AS STRING)")
    paste(value, "AS", paste0("k", i))
  }, character(1))
  max_rows = floor(limit / (length(components) + 1L))
  if (max_rows < 1L) prediction_state_overflow()
  dictionary = dbglm_get_query(conn, paste0("SELECT DISTINCT ", fe, " AS id,", paste(fields, collapse = ","),
    " FROM ", retained, " ORDER BY id LIMIT ", sprintf("%.0f", max_rows + 1)))
  if (nrow(dictionary) > max_rows) prediction_state_overflow()
  keys = prediction_tuple_key(dictionary[paste0("k", seq_along(components))], types)
  if (anyNA(keys) || anyDuplicated(keys) || anyDuplicated(dictionary[["id"]])) {
    stop("Could not preserve the fitted FE key dictionary exactly.", call. = FALSE)
  }
  list(columns = public_components, types = unname(types), key = keys, id = dictionary[["id"]])
}

dbglm_prediction_state = function(...) build_prediction_state(...)

build_prediction_state = function(conn, working, retained, fe, components, x, design_names,
                                  center, trend_aliases, sizes, create, drop, n,
                                  dictionary_components = components, beta = NULL,
                                  expected_rss = NULL, constant = 0, fitted_beta = beta) {
  tryCatch(build_portable_prediction_state(conn, working, retained, fe, components, x, design_names,
    center, trend_aliases, sizes, create, drop, n, dictionary_components, beta, expected_rss, constant),
    dbreg_prediction_overflow = function(error) {
      if (length(fe) != 1L || length(unlist(attr(fe, "slopes"))) || is.null(fitted_beta)) stop(error)
      build_database_prediction_state(conn, working, retained, fe, components, x, design_names,
        center, create, n, dictionary_components, fitted_beta, expected_rss, constant)
    })
}

build_portable_prediction_state = function(conn, working, retained, fe, components, x, design_names,
                                  center, trend_aliases, sizes, create, drop, n,
                                  dictionary_components = components, beta = NULL,
                                  expected_rss = NULL, constant = 0) {
  limit = prediction_state_limit()
  max_rows = getOption("dbreg.duckdb_ap_max_edges", 1e6)
  dbglm_validate_moment_limits(max_rows, limit)
  used = 0
  bounded_query = function(sql, width) {
    remaining = min(max_rows - used, floor(limit / width))
    if (remaining < 1) prediction_state_overflow()
    result = dbglm_get_query(conn, paste0(sql, " LIMIT ", sprintf("%.0f", remaining + 1)))
    if (nrow(result) > remaining) prediction_state_overflow()
    used = used + nrow(result)
    result
  }
  ops = lapply(seq_along(fe), function(i) list(f = fe[i], v = "1.0", group = i, variable = NULL))
  for (i in seq_along(fe)) for (variable in attr(fe, "slopes")[[i]]) {
    ops[[length(ops) + 1L]] = list(f = fe[i], v = variable, group = i,
      variable = names(trend_aliases)[match(variable, trend_aliases)])
  }
  source = paste0("(", dbglm_sql(conn, paste0("SELECT * EXCLUDE(w), __dbglm_a AS w,",
    "(eta-__dbglm_o) AS __dbglm_prediction_target FROM ", working)), ")")
  columns = c(x, "__dbglm_prediction_target")
  op_sizes = sizes[vapply(ops, `[[`, integer(1), "group")]
  total = sum(op_sizes) + length(x)
  rest = total - max(op_sizes)
  if (!dbglm_fe_dense_feasible(c(op_sizes, rep(1, length(x)))) || total * max(1, rest, length(columns)) > limit) {
    prediction_state_overflow()
  }
  if (length(fe) >= 2L && length(ops) == length(fe)) {
    sql = dbglm_fe_moments_sql(conn, source, fe, dbglm_fe_value_terms(conn, columns),
      limit = min(max_rows, floor(limit / (length(columns) + 4L))) + 1)
    res = dbglm_read_fe_moments(conn, sql)
    if (nrow(res) > min(max_rows, floor(limit / (length(columns) + 4L)))) prediction_state_overflow()
    used = nrow(res)
    moments = dbglm_unpack_fe_moments(res, length(fe), length(columns), 0L)
    if (is.null(moments)) prediction_state_overflow()
    info = moments[["info"]]
    sys = moments[["sys"]]
    B = moments[["B"]]
  } else {
    info = lapply(ops, function(op) bounded_query(paste0("SELECT ", op[["f"]], " AS id,",
      "FSUM(w*(", op[["v"]], ")*(", op[["v"]], ")) AS h,",
      paste(sprintf("FSUM(w*(%s)*(%s)) AS b%d", op[["v"]], columns, seq_along(columns)), collapse = ","),
      " FROM ", source, " GROUP BY ", op[["f"]], " ORDER BY ", op[["f"]]), length(columns) + 2L))
    pairs = list()
    for (i in seq_along(ops)) if (i > 1L) for (j in seq_len(i - 1L)) {
      oi = ops[[i]]
      oj = ops[[j]]
      edge = bounded_query(paste0("SELECT ", oi[["f"]], " AS i,", oj[["f"]], " AS j,",
        "FSUM(w*(", oi[["v"]], ")*(", oj[["v"]], ")) AS h FROM ", source,
        " GROUP BY ", paste(unique(c(oi[["f"]], oj[["f"]])), collapse = ",")), 3L)
      pairs[[length(pairs) + 1L]] = list(i = i, j = j, h = as.numeric(edge[["h"]]),
        ir = match(edge[["i"]], info[[i]][["id"]]), jr = match(edge[["j"]], info[[j]][["id"]]))
    }
    sys = dbglm_fe_system(info, pairs)
    if (is.null(sys)) prediction_state_overflow()
    B = do.call(rbind, lapply(info, function(group) as.matrix(group[paste0("b", seq_along(columns))])))
  }
  global = dbglm_get_query(conn, paste("SELECT", paste(dbglm_moment_terms(columns, "w"), collapse = ","), "FROM", source))
  cross = dbglm_read_information(global, length(columns))
  p = length(x)
  fe_n = sys[["n"]]
  fe_system = sys
  if (p) {
    old_rest = sys[["rest"]]
    bx = B[, seq_len(p), drop = FALSE]
    sys[["Cm"]] = cbind(sys[["Cm"]], bx[sys[["big"]], , drop = FALSE])
    sys[["E"]] = rbind(cbind(sys[["E"]], bx[old_rest, , drop = FALSE]),
      cbind(t(bx[old_rest, , drop = FALSE]), cross[seq_len(p), seq_len(p), drop = FALSE]))
    sys[["rest"]] = c(old_rest, fe_n + seq_len(p))
    sys[["n"]] = fe_n + p
    sys[["diag"]] = c(sys[["diag"]], diag(cross)[seq_len(p)])
  }
  rhs = matrix(c(B[, p + 1L], if (p) cross[seq_len(p), p + 1L]), ncol = 1L)
  solver = dbglm_fe_solver(sys, new.env(parent = emptyenv()), retain_null = TRUE)
  if (is.null(beta)) {
    rhs = solver[["project"]](rhs)
    coefficients = solver[["solve"]](rhs)
    # Recover the fitted GLM predictor, independently of its fitting projector.
    for (i in seq_len(3L)) coefficients = coefficients + solver[["solve"]](rhs - dbglm_fe_multiply(sys, coefficients))
  } else {
    # OLS already supplies its fitted slopes. Solve only the nuisance equations
    # for y-Xb, preserving those slopes and their original estimation strategy.
    if (length(beta) != p || any(!is.finite(beta)) || is.null(expected_rss)) {
      stop("OLS prediction recovery requires fitted slopes and residual sum of squares.", call. = FALSE)
    }
    fe_solver = dbglm_fe_solver(fe_system, new.env(parent = emptyenv()))
    residual_rhs = matrix(B[, p + 1L] - B[, seq_len(p), drop = FALSE] %*% beta, ncol = 1L)
    residual_rhs = fe_solver[["project"]](residual_rhs)
    gamma = fe_solver[["solve"]](residual_rhs)
    for (i in seq_len(3L)) gamma = gamma + fe_solver[["solve"]](residual_rhs - dbglm_fe_multiply(fe_system, gamma))
    coefficients = matrix(c(gamma, beta), ncol = 1L)
  }
  if (any(!is.finite(coefficients))) stop("Prediction coefficient recovery is nonfinite.", call. = FALSE)
  frame = dbglm_coefficient_frame(info, coefficients[seq_len(fe_n), , drop = FALSE])
  lookup = dbglm_coefficient_lookup(conn, frame, create)
  ctes = dbglm_coefficient_ctes(lookup[["source"]], length(ops), "pred")
  joins = vapply(seq_along(ops), function(i) paste0(" LEFT JOIN pred_g", i, " ON s.", ops[[i]][["f"]], "=pred_g", i, ".__glm_group"), character(1))
  terms = c(if (p) sprintf("(%s)*s.%s", dbglm_number(coefficients[fe_n + seq_len(p)]), x),
    vapply(seq_along(ops), function(i) sprintf("(%s)*pred_g%d.c1", if (ops[[i]][["v"]] == "1.0") "1.0" else paste0("s.", ops[[i]][["v"]]), i), character(1)))
  residual = paste0("((", paste(terms, collapse = "+"), ")-s.__dbglm_prediction_target)")
  diagnostics = if (is.null(expected_rss)) paste0("MAX(ABS(", residual,
    ")) AS error,MAX(ABS(s.__dbglm_prediction_target)) AS scale") else
    paste0("FSUM(s.w*POWER(", residual, ",2)) AS rss")
  error = dbglm_get_query(conn, paste0("WITH ", paste(ctes, collapse = ","), " SELECT ", diagnostics, " FROM ",
    source, " s", paste(joins, collapse = "")))
  for (table in lookup[["tables"]]) drop(table)
  if (is.null(expected_rss)) {
    recovery_error = error[["error"]]
    if (!is.finite(recovery_error) || recovery_error > 1e-7 * max(1, error[["scale"]])) {
      stop("Recovered prediction state does not reproduce the fitted linear predictor.", call. = FALSE)
    }
  } else {
    recovery_error = abs(error[["rss"]] - expected_rss)
    if (!is.finite(recovery_error) || recovery_error > 1e-6 * max(1, abs(expected_rss), cross[p + 1L, p + 1L])) {
      stop("Recovered OLS effects do not reproduce the fitted residual sum of squares (recovered=",
        signif(error[["rss"]], 9), ", fitted=", signif(expected_rss, 9), ").", call. = FALSE)
    }
  }
  dictionaries = lapply(seq_along(fe), function(i) {
    prediction_fe_dictionary(conn, retained, fe[i], unique(dictionary_components[[i]]),
      floor(limit / length(fe)), public_components = unique(components[[i]]))
  })
  list(version = 1L, specification = list(design_names = design_names, design_center = unname(center),
      dictionaries = dictionaries, operations = lapply(ops, function(op) op[c("group", "variable")])),
    coefficients = list(joint = as.numeric(coefficients), starts = sys[["starts"]],
      ids = lapply(info, `[[`, "id"), fe_n = fe_n),
    sample = list(nobs = n), inference = NULL, resources = NULL,
    prediction = list(null_basis = solver[["null_basis"]], tolerance = 1e-7,
      recovery_error = recovery_error, recovery_kind = if (is.null(expected_rss)) "linear_predictor" else "rss", constant = constant))
}

prediction_match_groups = function(state, data) {
  if (!identical(state[["version"]], 1L)) stop("Unsupported fitted-state version.", call. = FALSE)
  lapply(state[["specification"]][["dictionaries"]], function(dictionary) {
    match(prediction_tuple_key(data[dictionary[["columns"]]], dictionary[["types"]]), dictionary[["key"]])
  })
}

predict_from_fitted_state = function(state, design, data, groups = NULL, lookup = NULL) {
  if (!identical(state[["version"]], 1L)) stop("Unsupported fitted-state version.", call. = FALSE)
  if (identical(state[["prediction"]][["storage"]], "database")) {
    return(predict_from_database_state(state, design, lookup))
  }
  specification = state[["specification"]]
  coef = state[["coefficients"]]
  x = sweep(design[, specification[["design_names"]], drop = FALSE], 2L, specification[["design_center"]], "-")
  if (is.null(groups)) groups = prediction_match_groups(state, data)
  ops = specification[["operations"]]
  indices = values = vector("list", length(ops))
  for (i in seq_along(ops)) {
    group = ops[[i]][["group"]]
    dictionary = specification[["dictionaries"]][[group]]
    indices[[i]] = coef[["starts"]][i] + match(dictionary[["id"]][groups[[group]]], coef[["ids"]][[i]])
    values[[i]] = if (is.null(ops[[i]][["variable"]])) rep(1, nrow(data)) else as.numeric(data[[ops[[i]][["variable"]]]])
  }
  prediction = as.vector(x %*% coef[["joint"]][coef[["fe_n"]] + seq_len(ncol(x))])
  if (!is.null(state[["prediction"]][["constant"]])) prediction = prediction + state[["prediction"]][["constant"]]
  norm = rowSums(x * x)
  for (i in seq_along(ops)) {
    prediction = prediction + values[[i]] * coef[["joint"]][indices[[i]]]
    norm = norm + values[[i]]^2
  }
  null = state[["prediction"]][["null_basis"]]
  if (!is.null(null)) {
    chunk = max(1, min(4096, floor(prediction_state_limit() / ncol(null))))
    for (start in seq.int(1L, nrow(data), by = chunk)) {
      rows = seq.int(start, min(nrow(data), start + chunk - 1L))
      score = x[rows, , drop = FALSE] %*% null[coef[["fe_n"]] + seq_len(ncol(x)), , drop = FALSE]
      for (i in seq_along(ops)) score = score + values[[i]][rows] * null[indices[[i]][rows], , drop = FALSE]
      non_estimable = !is.finite(rowSums(score)) | sqrt(rowSums(score^2)) > state[["prediction"]][["tolerance"]] * sqrt(norm[rows])
      prediction[rows[non_estimable]] = NA_real_
    }
  }
  prediction[!is.finite(prediction)] = NA_real_
  prediction
}
