build_database_two_fe_state = function(conn, working, retained, fe, components, x,
  design_names, center, create, drop, n, dictionary_components, beta, expected_rss,
  constant, recovery_control) {
  p = length(x)
  if ((p + 1)^2 > prediction_state_limit()) prediction_state_overflow()
  if (length(beta) != p || any(!is.finite(beta))) stop("Database prediction requires finite fitted slopes.", call. = FALSE)
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  slope = if (p) paste(sprintf("(%s)*s.%s", dbglm_number(beta), x), collapse = "+") else "0.0"
  target = "__dbglm_prediction_target"
  fields = c(paste0("s.", fe), "s.__dbglm_a AS w", paste0("(s.eta-s.__dbglm_o-(", slope, ")) AS ", target), if (p) paste0("s.", x))
  source = create(paste("SELECT", paste(fields, collapse = ","), "FROM", working, "s"))
  projection = prediction_two_fe_projection(conn, source, fe, c(target, x), create, drop,
    tol = min(recovery_control[["tol"]], 1e-10), max_iter = recovery_control[["max_iter"]])
  scale = projection[["scales"]]
  within = projection[["residual"]]
  columns = projection[["columns"]]
  diagnostics = if (is.null(expected_rss)) paste0("MAX(ABS(", columns[[1L]], "))*(", dbglm_number(scale[[1L]]), ") AS error") else
    paste0("FSUM(w*POWER(", columns[[1L]], "*(", dbglm_number(scale[[1L]]), "),2)) AS rss")
  diagnostic = backend_get_query(conn, paste("SELECT", diagnostics, "FROM", within))
  recovery_error = if (is.null(expected_rss)) as.numeric(diagnostic[["error"]]) else abs(as.numeric(diagnostic[["rss"]]) - expected_rss)
  if (is.null(expected_rss)) {
    bound = 1e-7 * scale[[1L]]
  } else {
    norm = backend_get_query(conn, paste("SELECT FSUM(w*POWER(", target, ",2)) AS norm FROM", source))[["norm"]]
    bound = 1e-6 * max(1, abs(expected_rss), as.numeric(norm))
  }
  if (!is.finite(recovery_error) || recovery_error > bound) {
    stop("Database two-FE recovery does not reproduce the fitted predictor or RSS.", call. = FALSE)
  }
  null = NULL
  if (p) {
    H = dbglm_read_information(backend_get_query(conn, paste("SELECT", paste(dbglm_moment_terms(columns[-1L], "w"), collapse = ","), "FROM", within)), p)
    normalized = paste0("(SELECT w,", paste(sprintf("(%s)/(%s) AS %s", x, dbglm_number(scale[-1L]), columns[-1L]), collapse = ","), " FROM ", source, ")")
    raw = dbglm_read_information(backend_get_query(conn, paste("SELECT", paste(dbglm_moment_terms(columns[-1L], "w"), collapse = ","), "FROM", normalized)), p)
    null = prediction_within_null(H, raw)
  }
  labels = prediction_two_fe_components(conn, source, fe, create, drop,
    max_iter = recovery_control[["max_iter"]])
  component_tables = list(labels[["left"]], labels[["right"]])
  mx = if (p) paste0("m", seq_len(p)) else character()
  dictionaries = lapply(seq_along(fe), function(i) {
    metadata = prediction_key_metadata(conn, retained, dictionary_components[[i]])
    keys = paste0("k", seq_along(dictionary_components[[i]]))
    dictionary_keys = paste("SELECT DISTINCT", paste(c(paste(fe[[i]], "AS id"),
      paste(q(dictionary_components[[i]]), "AS", keys)), collapse = ","), "FROM", retained)
    fields = c("a.id", paste0("a.", columns[[1L]], "*(", dbglm_number(scale[[1L]]), ") AS gamma"),
      sprintf("a.%s*(%s) AS %s", columns[-1L], dbglm_number(scale[-1L]), mx),
      "l.component", paste0("d.", keys))
    dictionary = create(paste("SELECT", paste(fields, collapse = ","), "FROM", projection[["dictionaries"]][[i]],
      "a JOIN", component_tables[[i]], "l ON a.id=l.id JOIN (", dictionary_keys, ") d ON a.id=d.id"))
    list(columns = unique(components[[i]]), types = metadata[["types"]], sql_types = unname(metadata[["sql_types"]]),
      table = dictionary, keys = keys, means = mx, component = "component")
  })
  list(version = 1L, specification = list(design_names = design_names, design_center = unname(center), dictionaries = dictionaries),
    coefficients = list(slopes = unname(beta)), sample = list(nobs = n), inference = NULL,
    resources = list(tables = gsub('["`]', "", vapply(dictionaries, `[[`, character(1), "table")), storage = "database", version = 1L),
    prediction = list(storage = "database", null_basis = null, tolerance = 1e-7, design_scale = scale[-1L],
      recovery_error = recovery_error, recovery_kind = if (is.null(expected_rss)) "linear_predictor" else "rss", constant = constant,
      iterations = projection[["iterations"]], projection_error = projection[["error"]], component_iterations = labels[["iterations"]]))
}
