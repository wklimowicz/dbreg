# A single FE block has a diagonal nuisance system. Keep its dictionary and
# within-group regressor means in SQL; only a p-by-p estimability system moves
# into R. Multiple FE blocks/varying slopes still use the portable solver.
build_database_prediction_state = function(conn, working, retained, fe, components, x,
  design_names, center, create, n, dictionary_components, beta, expected_rss, constant) {
  p = length(x)
  if ((p + 1)^2 > prediction_state_limit()) prediction_state_overflow()
  if (length(beta) != p || any(!is.finite(beta))) stop("Database prediction requires finite fitted slopes.", call. = FALSE)
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  metadata = prediction_key_metadata(conn, retained, dictionary_components[[1L]])
  mx = paste0("m", seq_len(p))
  slope = if (p) paste(sprintf("(%s)*s.%s", dbglm_number(beta), x), collapse = "+") else "0.0"
  target = "(s.eta-s.__dbglm_o)"
  fields = c(paste0("s.", fe, " AS id"),
    paste0("FSUM(s.__dbglm_a*(", target, "-(", slope, ")))/FSUM(s.__dbglm_a) AS gamma"),
    sprintf("FSUM(s.__dbglm_a*s.%s)/FSUM(s.__dbglm_a) AS %s", x, mx))
  aggregate = paste("SELECT", paste(fields, collapse = ","), "FROM", working, "s GROUP BY", paste0("s.", fe))
  keys = paste0("k", seq_along(dictionary_components[[1L]]))
  dictionary_keys = paste("SELECT DISTINCT", paste(c(paste(fe, "AS id"),
    paste(q(dictionary_components[[1L]]), "AS", keys)), collapse = ","), "FROM", retained)
  dictionary = create(paste0("SELECT a.*,", paste0("d.", keys, collapse = ","),
    " FROM (", aggregate, ") a JOIN (", dictionary_keys, ") d ON a.id=d.id"), materialize = TRUE)
  join = paste("FROM", working, "s JOIN", dictionary, "g ON", paste0("s.", fe, "=g.id"))
  residual = paste0("((", slope, ")+g.gamma-", target, ")")
  diagnostics = if (is.null(expected_rss)) paste0("MAX(ABS(", residual,
    ")) AS error,MAX(ABS(", target, ")) AS scale") else
    paste0("FSUM(s.__dbglm_a*POWER(", residual, ",2)) AS rss,FSUM(s.__dbglm_a*POWER(", target, ",2)) AS scale")
  diagnostic = dbglm_get_query(conn, paste("SELECT", diagnostics, join))
  recovery_error = if (is.null(expected_rss)) as.numeric(diagnostic[["error"]]) else abs(as.numeric(diagnostic[["rss"]]) - expected_rss)
  tolerance = if (is.null(expected_rss)) 1e-7 else 1e-6
  scale = max(1, as.numeric(diagnostic[["scale"]]), if (!is.null(expected_rss)) abs(expected_rss))
  if (!is.finite(recovery_error) || recovery_error > tolerance * scale) {
    stop("Database prediction recovery does not reproduce the fitted predictor or RSS.", call. = FALSE)
  }
  null = NULL
  if (p) {
    within = paste0("(", paste("SELECT s.__dbglm_a AS w,", paste(sprintf("(s.%s-g.%s) AS %s", x, mx, x), collapse = ","), join), ")")
    H = dbglm_read_information(dbglm_get_query(conn, paste("SELECT", paste(dbglm_moment_terms(x, "w"), collapse = ","), "FROM", within)), p)
    raw = dbglm_read_information(dbglm_get_query(conn, paste("SELECT", paste(dbglm_moment_terms(x, "__dbglm_a"), collapse = ","), "FROM", working)), p)
    absorbed = diag(H) <= 1e-12 * pmax(diag(raw), .Machine[["double.xmin"]])
    H[absorbed, ] = 0
    H[, absorbed] = 0
    system = list(n = p + 1L, big = 1L, rest = 1L + seq_len(p), Dg = 1,
      Cm = matrix(0, 1L, p), E = H, diag = c(1, diag(raw)))
    solver = dbglm_fe_solver(system, new.env(parent = emptyenv()), retain_null = TRUE)
    if (!is.null(solver[["null_basis"]])) null = solver[["null_basis"]][-1L, , drop = FALSE]
  }
  list(version = 1L, specification = list(design_names = design_names, design_center = unname(center),
      dictionaries = list(list(columns = unique(components[[1L]]), types = metadata[["types"]],
        sql_types = unname(metadata[["sql_types"]]), table = dictionary, keys = keys, means = mx))),
    coefficients = list(slopes = unname(beta)), sample = list(nobs = n), inference = NULL,
    resources = list(tables = gsub('["`]', "", dictionary), storage = "database", version = 1L),
    prediction = list(storage = "database", null_basis = null, tolerance = 1e-7,
      recovery_error = recovery_error, recovery_kind = if (is.null(expected_rss)) "linear_predictor" else "rss", constant = constant))
}

finalize_retained_resources = function(result, resources) {
  retained = result[["fitted_state"]][["resources"]]
  if (identical(retained[["storage"]], "database")) {
    result[["resource_state"]] = retain_fit_resources(resources, retained[["tables"]])
  }
  result
}

# Match only the requested rows, in bounded batches. The SQL dictionary keeps
# original key types; incoming integers use exact decimal text, never doubles.
prediction_database_lookup = function(state, lease, data) {
  conn = fitted_resource_connection(lease)
  dictionary = state[["specification"]][["dictionaries"]][[1L]]
  keys = dictionary[["keys"]]
  values = Map(function(name, type) {
    encoded = prediction_key_part(data[[name]], type)
    sub("^[0-9]+:", "", encoded)
  }, dictionary[["columns"]], dictionary[["types"]])
  fields = c("gamma", dictionary[["means"]])
  result = as.data.frame(stats::setNames(rep(list(rep(NA_real_, nrow(data))), length(fields)), fields))
  if (!nrow(data)) return(result)
  batch = max(1L, min(4096L, floor(prediction_state_limit() / (length(fields) + length(keys) + 1L))))
  casts = sprintf("TRY_CAST(v.%s AS %s)", keys, dictionary[["sql_types"]])
  join = paste(sprintf("g.%s=%s", keys, casts), collapse = " AND ")
  format = if (isTRUE(detect_backend(conn)[["is_databricks"]])) "json" else "values"
  for (first in seq.int(1L, nrow(data), by = batch)) {
    rows = seq.int(first, min(nrow(data), first + batch - 1L))
    requested = data.frame(i = as.numeric(rows))
    requested[keys] = lapply(values, function(value) value[rows])
    values_sql = dbglm_values_sql(conn, requested, max_bytes = 2 * 1024^2, format = format)
    if (is.null(values_sql)) {
      stop("Prediction key batch exceeds the SQL transfer limit.", call. = FALSE)
    }
    source = paste0("(", values_sql, ") v")
    sql = paste0("SELECT v.i,", paste0("g.", fields, collapse = ","), " FROM ", source,
      " LEFT JOIN ", dictionary[["table"]], " g ON ", join, " ORDER BY v.i")
    found = dbglm_get_query(conn, sql)
    if (nrow(found) != length(rows) || !identical(as.numeric(found[["i"]]), as.numeric(rows))) {
      stop("Prediction dictionary no longer has unique fitted keys.", call. = FALSE)
    }
    result[rows, ] = lapply(found[fields], as.numeric)
  }
  result
}

predict_from_database_state = function(state, design, lookup) {
  if (is.null(lookup)) stop("Database prediction requires its fitted key lookup.", call. = FALSE)
  specification = state[["specification"]]
  x = sweep(design[, specification[["design_names"]], drop = FALSE], 2L,
    specification[["design_center"]], "-")
  prediction = as.vector(x %*% state[["coefficients"]][["slopes"]]) + lookup[["gamma"]] + state[["prediction"]][["constant"]]
  null = state[["prediction"]][["null_basis"]]
  if (!is.null(null)) {
    means = as.matrix(lookup[specification[["dictionaries"]][[1L]][["means"]]])
    score = (x - means) %*% null
    norm = sqrt(1 + rowSums(x * x) + rowSums(means * means))
    prediction[sqrt(rowSums(score * score)) > state[["prediction"]][["tolerance"]] * norm] = NA_real_
  }
  prediction[!is.finite(prediction)] = NA_real_
  prediction
}
