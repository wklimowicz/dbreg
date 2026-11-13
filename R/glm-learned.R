# Preserve training-key identities before the owned fitting tables disappear.
learned_key_dictionary = function(conn, table, fe, columns, public, n, byte_budget) {
  quoted = as.character(DBI::dbQuoteIdentifier(conn, columns))
  types = category_key_metadata(conn, table, columns)[["types"]]
  values = quoted
  integer = types == "integer"
  values[integer] = sql_cast_text(conn, quoted[integer])
  aliases = paste0("k", seq_along(columns))
  keys = paste(values, "AS", aliases)
  text = sql_cast_text(conn, quoted)
  if (detect_backend(conn)[["name"]] == "duckdb") text = paste0("encode(", text, ")")
  sizes = paste0("OCTET_LENGTH(", text, ")")
  largest = if (length(sizes) == 1L) sizes else paste0("GREATEST(", paste(sizes, collapse = ","), ")")
  bounds = DBI::dbGetQuery(conn, paste0("SELECT MAX(", largest, ") AS largest, MAX(",
    paste(sizes, collapse = "+"), ") AS combined FROM ", table))
  if (any(!is.finite(unlist(bounds))) || bounds[["largest"]] >= 65535 ||
      n * bounds[["combined"]] > byte_budget) {
    projection_abort("transfer_limit", "A fitted key exceeds the transfer byte limit.")
  }
  limit = format(n + 1, scientific = FALSE, trim = TRUE)
  rows = DBI::dbGetQuery(conn, paste("SELECT DISTINCT", paste(c(paste(fe, "AS id"), keys),
    collapse = ","), "FROM", table, "ORDER BY id LIMIT", limit))
  if (nrow(rows) != n || !identical(as.numeric(rows[["id"]]), as.numeric(seq_len(n)))) {
    stop("The fitted key dictionary does not match its coefficient ranks.", call. = FALSE)
  }
  domain = composite_key(rows, aliases, types)
  if (anyNA(domain) || anyDuplicated(domain)) stop("Invalid fitted key domain.", call. = FALSE)
  list(columns = public, types = types, keys = domain)
}

check_learned_budget = function(sizes, slopes, p, groups, limits) {
  axis_sizes = c(sizes, if (length(slopes)) rep(sizes, lengths(slopes)))
  total = sum(axis_sizes) + p
  entries = getOption("dbreg.prediction_max_entries", 2e6)
  if (!is.numeric(entries) || length(entries) != 1L || !is.finite(entries) || entries < 1) {
    stop("Invalid dbreg.prediction_max_entries.", call. = FALSE)
  }
  # The option budgets dictionary entries; the dense recovery system is bounded
  # by the projection limits, as at fit time.
  if (total > entries ||
      sum(sizes * (1 + lengths(groups))) > limits[["transfer_cells"]] ||
      !fe_dense_feasible(axis_sizes, max(1, p), limits)) {
    projection_abort("system_limit", "The portable prediction state exceeds its resource budget.")
  }
}

dbglm_learned_state = function(conn, sample, object, fit) {
  fe = sample[["fe"]]
  sizes = sample[["sizes"]]
  x = sample[["columns"]][["alias"]]
  p = length(x)
  limits = projection_limits()
  check_learned_budget(sizes, sample[["slopes"]], p, sample[["key_groups"]], limits)
  dictionaries = lapply(seq_along(fe), function(i) {
    learned_key_dictionary(conn, sample[["key_source"]], fe[i],
      sample[["key_groups"]][[i]], sample[["public_groups"]][[i]], sizes[i],
      limits[["memory_bytes"]] / (8 * length(fe)))
  })
  budget = list2env(list(rows = 0, cells = 0, summary_bytes = 0, upload_bytes = 0), parent = emptyenv())
  if (length(sample[["trend_columns"]])) {
    beta = unname(object[["collin.coef"]])
    beta[is.na(beta)] = 0
    state = recover_learned_effects(conn, fit[["working"]], fe, sample[["slopes"]],
      sample[["trend_columns"]], sizes, x, beta, sample[["centers"]], sample[["fit_rows"]],
      budget, limits)
    state[["dictionaries"]] = dictionaries
    return(state)
  }
  source = paste0("(SELECT ", paste(c(fe, x, "__dbglm_a AS w", "1.0 AS check_axis"),
    collapse = ","), " FROM ", sample[["original_base"]], ") prediction_source")
  columns = c(x, "check_axis")
  summaries = query_fe_summaries(conn, source, fe, columns, sizes, budget, limits)
  sys = fe_system(summaries[["info"]], summaries[["pairs"]], limits)
  if (is.null(sys)) projection_abort("system_limit", "The prediction fixed-effect system exceeds its resource budget.")
  starts = sys[["starts"]]
  null = learned_design_null(conn, sys, summaries[["info"]], source, x, budget, limits)
  effects = fit[["effects"]]
  if (is.character(effects)) {
    rows = DBI::dbGetQuery(conn, paste("SELECT c1 FROM", effects,
      "ORDER BY __cg_op,__cg_group LIMIT", format(sum(sizes) + 1, scientific = FALSE)))
    effects = list(as.numeric(rows[["c1"]]))
  }
  alpha = unlist(effects, use.names = FALSE)
  if (length(alpha) != sum(sizes) || any(!is.finite(alpha))) {
    stop("Invalid fitted GLM fixed-effect coefficients.", call. = FALSE)
  }
  beta = unname(object[["collin.coef"]])
  beta[is.na(beta)] = 0
  list(dictionaries = dictionaries, alpha = alpha, beta = beta, starts = starts,
    center = sample[["centers"]], null = null)
}

learned_design_null = function(conn, sys, info, source, x, budget, limits) {
  if (length(x)) {
    X = fe_information_query(conn, source, x, budget, limits)[["cross"]]
    B = do.call(rbind, lapply(info, `[[`, "b"))[, seq_along(x), drop = FALSE]
    sys[["Cm"]] = cbind(sys[["Cm"]], B[sys[["big"]], , drop = FALSE])
    BX = B[sys[["rest"]], , drop = FALSE]
    sys[["E"]] = rbind(cbind(sys[["E"]], BX), cbind(t(BX), X))
    sys[["rest"]] = c(sys[["rest"]], sys[["n"]] + seq_along(x))
    sys[["n"]] = sys[["n"]] + length(x)
    sys[["diag"]] = c(sys[["diag"]], diag(X))
  }
  solver = fe_solver(sys, list2env(list(projection_limits = limits), parent = emptyenv()),
    retain_null = TRUE)
  solver[["null_basis"]]
}

learned_effect_operations = function(state) {
  ops = state[["ops"]]
  if (!is.null(ops)) return(ops)
  dictionaries = state[["dictionaries"]]
  lapply(seq_along(dictionaries), function(i) {
    list(group = i, index = state[["starts"]][i] + seq_along(dictionaries[[i]][["keys"]]))
  })
}

predict_learned_effects = function(object, data) {
  state = object[["fitted_state"]]
  design = linear_prediction_matrix(object, data, all_columns = TRUE, ids = state[["design_ids"]])
  invalid = attr(design, "dbreg_invalid_rows")
  offset = linear_prediction_offset(object, data)
  n = nrow(data)
  if (!n) return(numeric())
  groups = lapply(seq_along(state[["dictionaries"]]), function(i) {
    dict = state[["dictionaries"]][[i]]
    match(composite_key(data, dict[["columns"]], dict[["types"]]), dict[["keys"]])
  })
  eta = offset + if (is.null(state[["constant"]])) 0 else state[["constant"]]
  ops = learned_effect_operations(state)
  values = indices = list()
  for (i in seq_along(ops)) {
    op = ops[[i]]
    g = groups[[op[["group"]]]]
    if (is.null(op[["trend"]])) {
      value = rep(1, n)
    } else {
      trend = data[[op[["trend"]]]]
      if (is.null(trend) || !is.numeric(trend)) {
        stop("Prediction requires numeric trend column: ", op[["trend"]], call. = FALSE)
      }
      value = (trend - op[["center"]][g]) / op[["scale"]][g]
    }
    index = op[["index"]][g]
    invalid = invalid | is.na(g) | !is.finite(value) | (is.na(index) & abs(value) > 1e-10)
    values[[i]] = value
    indices[[i]] = index
    eta = eta + ifelse(is.na(index), 0, state[["alpha"]][index] * value)
  }
  null = state[["null"]]
  width = ncol(design)
  block = max(1L, min(4096L, floor(2e6 / max(1L, width, ncol(null)))))
  for (start in seq.int(1L, n, by = block)) {
    rows = seq.int(start, min(n, start + block - 1L))
    centered = sweep(as.matrix(design[rows, , drop = FALSE]), 2L, state[["center"]], "-")
    if (width) eta[rows] = eta[rows] + as.numeric(centered %*% state[["beta"]])
    if (!is.null(null)) {
      q = null[length(state[["alpha"]]) + seq_len(width), , drop = FALSE]
      score = centered %*% q
      norm = rowSums(centered^2)
      for (i in seq_along(indices)) {
        idx = indices[[i]][rows]
        value = values[[i]][rows]
        known = !is.na(idx) & is.finite(value)
        score[known, ] = score[known, , drop = FALSE] +
          null[idx[known], , drop = FALSE] * value[known]
        norm[known] = norm[known] + value[known]^2
      }
      invalid[rows] = invalid[rows] | !is.finite(rowSums(score)) |
        rowSums(score^2) > 1e-14 * norm
    }
  }
  eta[invalid | !is.finite(eta)] = NA_real_
  eta
}
