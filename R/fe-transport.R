# Bounded FE moment and coefficient transport.
# SQL/list/JSON schemas live here; system assembly and solving stay in fe-solver.R.

# One statement returns the global row (__set = 0), per-FE marginals
# (__set = i) and pairwise FE information (__set = K + pair index). Pairs need
# only the weight sum, so each is a narrow GROUP BY over one shared source;
# GROUPING SETS would evaluate every aggregate for every set.
fe_moments_sql = function(conn, source, fe, value_terms, global_terms = character(), limit = NULL, compact = FALSE, text_keys = FALSE) {
  if (identical(detect_backend(conn)[["name"]], "spark")) {
    sql = dbglm_databricks_moments_sql(source, fe, value_terms, global_terms, limit, text_keys = text_keys, conn = conn)
    if (!is.null(sql)) return(if (compact) pack_fe_moments_sql(conn, sql, length(fe), value_terms, global_terms, globals_packed = TRUE) else sql)
  }
  K = length(fe)
  pairs = utils::combn(K, 2L, simplify = FALSE)
  alias = function(terms) sub(".* AS ([A-Za-z_0-9]+)$", "\\1", terms)
  null = function(terms) if (length(terms)) paste0("CAST(NULL AS DOUBLE) AS ", alias(terms)) else character()
  id_type = if (text_keys) sql_text_type(conn) else "BIGINT"
  id = function(f) if (is.null(f)) paste0("CAST(NULL AS ", id_type, ")") else paste0("CAST(", f, " AS ", id_type, ")")
  branch = function(set, i, j, columns, group) paste0(
    "SELECT ", set, " AS __set,", id(i), " AS __i,", id(j), " AS __j,",
    paste(columns, collapse = ","), " FROM __dbglm_src", if (length(group)) paste(" GROUP BY", paste(group, collapse = ","))
  )
  branches = c(
    if (length(global_terms)) branch(0L, NULL, NULL, c(global_terms, null(value_terms)), character()),
    vapply(seq_len(K), function(i) branch(i, fe[i], NULL, c(null(global_terms), value_terms), fe[i]), character(1)),
    vapply(seq_along(pairs), function(k) {
      p = pairs[[k]]
      branch(K + k, fe[p[1]], fe[p[2]], c(null(global_terms), value_terms[1L], null(value_terms[-1L])), fe[p])
    }, character(1))
  )
  materialized = if (isTRUE(detect_backend(conn)[["is_databricks"]])) "" else "MATERIALIZED "
  sql = paste0("WITH __dbglm_src AS ", materialized, "(SELECT * FROM ", source, ") ", paste(branches, collapse = " UNION ALL "))
  if (!is.null(limit)) sql = paste0("SELECT * FROM (", sql, ") __dbglm_m ORDER BY __set LIMIT ", sprintf("%.0f", limit))
  if (compact) pack_fe_moments_sql(conn, sql, K, value_terms, global_terms) else sql
}

# Pair rows carry only indices and weight. Only marginal rows carry RHS
# vectors; only the global row carries global moments. Packing before DBI
# transfer avoids recreating a wide mostly-NULL rectangle in client memory.
pack_fe_moments_sql = function(conn, sql, K, value_terms, global_terms = character(), globals_packed = FALSE) {
  alias = function(terms) sub(".* AS ([A-Za-z_0-9]+)$", "\\1", terms)
  values = alias(value_terms[-1L])
  native = identical(detect_backend(conn)[["name"]], "spark")
  payload = if (length(values)) paste0(
    "CASE WHEN __set BETWEEN 1 AND ", K,
    " THEN ", if (native) "to_json(array(" else "list_value(", paste(values, collapse = ","),
    if (native) ")) END" else ") END"
  ) else if (native) sql_cast_text(conn, "NULL") else "CAST(NULL AS DOUBLE[])"
  fields = c("__set", "__i", "__j", "__h", paste(payload, "AS __values"))
  if (length(global_terms)) {
    globals = alias(global_terms)
    encoded = if (globals_packed) "__globals" else paste0(
      "CASE WHEN __set=0 THEN list_value(",
      paste(globals, collapse = ","), ") END"
    )
    fields = c(fields, paste(encoded, "AS __globals"))
  }
  result = paste0("SELECT ", paste(fields, collapse = ","), " FROM (", sql, ") compact_fe_moments")
  if (length(global_terms)) attr(result, "fe_global_names") = alias(global_terms)
  result
}

read_fe_moments = function(conn, sql, max_rows = Inf) {
  result = DBI::dbGetQuery(conn, sql)
  if (nrow(result) > max_rows) return(NULL)
  if ("__values" %in% names(result)) {
    payload = result[["__values"]]
    rows = if (is.list(payload)) which(lengths(payload) > 0L) else which(!is.na(payload))
    values = lapply(payload[rows], function(x) if (is.list(payload)) as.numeric(x) else as.numeric(read_numeric_json(x)))
    if (length(unique(lengths(values))) > 1L) stop("Inconsistent fixed-effect marginal widths.", call. = FALSE)
    result[["__values"]] = NULL
    attr(result, "fe_values") = list(rows = rows, values = if (length(rows)) do.call(rbind, values) else matrix(numeric(), 0L, 0L))
  }
  if (!"__globals" %in% names(result)) return(result)
  global = which(as.integer(result[["__set"]]) == 0L)
  payload = result[["__globals"]]
  if (length(global) != 1L || (is.list(payload) && !length(payload[[global]])) ||
    (!is.list(payload) && is.na(payload[global]))) {
    stop("GLM moment query did not return its global summary.", call. = FALSE)
  }
  values = if (is.list(payload)) stats::setNames(as.list(as.numeric(payload[[global]])), attr(sql, "fe_global_names")) else
    read_numeric_json(payload[global])
  # Do not recreate mostly-NULL global columns in R: their number grows as
  # p^2, and padding them across every FE-pair row can exhaust client memory.
  globals = as.data.frame(lapply(values, as.numeric), check.names = FALSE)
  result[["__globals"]] = NULL
  attr(result, "dbglm_globals") = globals
  result
}

decode_fe_moments = function(res, K, n_cols, n_extra) {
  sets = as.integer(res[["__set"]])
  b_cols = if (n_cols) paste0("__b", seq_len(n_cols)) else character()
  ex_cols = if (n_extra) paste0("__ex", seq_len(n_extra)) else character()
  compact = attr(res, "fe_values", exact = TRUE)
  info = vector("list", K)
  for (i in seq_len(K)) {
    rows = which(sets == i)
    if (!length(rows)) return(NULL)
    rows = rows[order(res[["__i"]][rows])]
    h = as.numeric(res[["__h"]][rows])
    if (is.null(compact)) {
      b = as.matrix(res[rows, b_cols, drop = FALSE])
      ex = as.matrix(res[rows, ex_cols, drop = FALSE])
    } else {
      positions = match(rows, compact[["rows"]])
      if ((n_cols + n_extra) && (anyNA(positions) || ncol(compact[["values"]]) != n_cols + n_extra)) return(NULL)
      b = if (n_cols) compact[["values"]][positions, seq_len(n_cols), drop = FALSE] else matrix(numeric(), length(rows), 0L)
      ex = if (n_extra) compact[["values"]][positions, n_cols + seq_len(n_extra), drop = FALSE] else matrix(numeric(), length(rows), 0L)
    }
    if (anyNA(res[["__i"]][rows]) || anyDuplicated(res[["__i"]][rows]) || any(!is.finite(h)) || any(h <= 0) || any(!is.finite(b)) || any(!is.finite(ex))) return(NULL)
    info[[i]] = list(id = res[["__i"]][rows], h = h, b = unname(b), ex = unname(ex))
  }
  combos = utils::combn(K, 2L, simplify = FALSE)
  pairs = vector("list", length(combos))
  for (k in seq_along(combos)) {
    p = combos[[k]]
    rows = which(sets == K + k)
    if (!length(rows)) return(NULL)
    h = as.numeric(res[["__h"]][rows])
    ir = match(res[["__j"]][rows], info[[p[2]]][["id"]])
    jr = match(res[["__i"]][rows], info[[p[1]]][["id"]])
    if (any(!is.finite(h)) || anyNA(ir) || anyNA(jr)) return(NULL)
    pairs[[k]] = list(i = p[2], j = p[1], ir = ir, jr = jr, h = h)
  }
  list(info = info, pairs = pairs, rows = sum(sets != 0L), n_extra = n_extra)
}

fe_value_terms = function(conn, columns, extra_sums = character(), backend_info = detect_backend(conn)) {
  sum_sql = function(expr) sql_sum(conn, expr, backend_info)
  c(
    paste(sum_sql("w"), "AS __h"),
    if (length(columns)) paste(sum_sql(sprintf("w*(%s)", columns)), "AS", paste0("__b", seq_along(columns))),
    if (length(extra_sums)) sprintf("%s AS __ex%d", vapply(extra_sums, function(e) sql_stable_sum(conn, e, backend_info), character(1)), seq_along(extra_sums))
  )
}

# Account across the separate marginal/pair statements needed by varying slopes.
# Reserve an overflow probe before each DBI fetch and propagate typed capacity
# rejection through the same solver policy as packed intercept-only summaries.
projection_fetch = function(conn, sql, budget, width, limits) {
  rows = min(
    limits[["summary_rows"]] - budget[["rows"]],
    floor((limits[["transfer_cells"]] - budget[["cells"]]) / width) - 1L
  )
  if (rows < 1L) projection_abort("transfer_limit", "Fixed-effect summaries exceed the transfer budget.")
  result = DBI::dbGetQuery(conn, paste(sql, "LIMIT", sprintf("%.0f", rows + 1L)))
  if (nrow(result) > rows) projection_abort("transfer_limit", "Fixed-effect summaries exceed the transfer budget.")
  budget[["rows"]] = budget[["rows"]] + nrow(result)
  budget[["cells"]] = budget[["cells"]] + width * (nrow(result) + 1L)
  result
}

# Databricks inlines CTEs, so a UNION branch per FE pair repeatedly plans the
# predictor. Keep three aggregation branches: globals, FE marginals and pairs.
# Pairs need only their weight sum, not every diagnostic/score aggregate.
dbglm_databricks_moments_sql = function(source, fe, value_terms, global_terms = character(), limit = NULL, text_keys = FALSE, conn) {
  K = length(fe)
  # GROUPING_ID uses a signed BIGINT. Keep masks exactly representable in R.
  if (K < 2L || K > 30L) return(NULL)
  aliases = function(terms) sub(".* AS ([A-Za-z_0-9]+)$", "\\1", terms)
  id_type = if (text_keys) sql_text_type(conn) else "BIGINT"
  value_names = aliases(value_terms)
  global_names = aliases(global_terms)
  nulls = function(names) if (length(names)) paste0("CAST(NULL AS DOUBLE) AS ", names) else character()
  grouped = function(groups, sets, terms) {
    masks = vapply(groups, function(g) sum(2^(K - setdiff(seq_len(K), g))), numeric(1))
    key = function(j) {
      cases = vapply(seq_along(groups), function(i) {
        g = groups[[i]]
        if (length(g) < j) return("")
        sprintf("WHEN %.0f THEN CAST(%s AS %s)", masks[[i]], fe[g[[j]]], id_type)
      }, character(1))
      if (!any(nzchar(cases))) return(paste0("CAST(NULL AS ", id_type, ")"))
      paste("CASE __gid", paste(cases, collapse = " "), "END")
    }
    fields = c(
      paste("CASE __gid", paste(sprintf("WHEN %.0f THEN %d", masks, sets), collapse = " "), "END AS __set"),
      paste(key(1L), "AS __i"), paste(key(2L), "AS __j"),
      nulls(global_names), aliases(terms), nulls(setdiff(value_names, aliases(terms)))
    )
    # Keep the UNION schema in the same order even when pairs omit values.
    paste0(
      "SELECT ", paste(fields, collapse = ","), " FROM (SELECT GROUPING_ID(", paste(fe, collapse = ","),
      ") AS __gid,", paste(c(fe, terms), collapse = ","), " FROM ", source, " GROUP BY GROUPING SETS (",
      paste(vapply(groups, function(g) paste0("(", paste(fe[g], collapse = ","), ")"), character(1)), collapse = ","), ")) grouped"
    )
  }
  pairs = utils::combn(K, 2L, simplify = FALSE)
  branches = c(
    if (length(global_terms)) paste0(
      "SELECT 0 AS __set,CAST(NULL AS ", id_type, ") AS __i,CAST(NULL AS ", id_type, ") AS __j,",
      paste(c(global_terms, nulls(value_names)), collapse = ","), " FROM ", source
    ),
    grouped(as.list(seq_len(K)), seq_len(K), value_terms),
    grouped(pairs, K + seq_along(pairs), value_terms[1L])
  )
  sql = paste(branches, collapse = " UNION ALL ")
  # Keep stable set ordering for bounded preflight queries. The immutable FE
  # support lets IRLS reuse this bound without sorting every evaluation.
  if (!is.null(limit)) sql = paste0("SELECT * FROM (", sql, ") limited ORDER BY __set LIMIT ", sprintf("%.0f", limit))
  if (!length(global_terms)) return(sql)
  # Only the global row has these scalars. Avoid sending NULL copies on every
  # pair row; the reader keeps globals separate from the numeric FE moments.
  fields = paste(sprintf("'%s',%s", global_names, global_names), collapse = ",")
  paste0(
    "SELECT __set,__i,__j,", paste(value_names, collapse = ","),
    ",CASE WHEN __set=0 THEN to_json(named_struct(", fields, ")) END AS __globals FROM (", sql, ") packed"
  )
}

# Small, dense FE rank dictionaries can use constant array indexing instead
# of coefficient joins. Large/sparse dictionaries keep the bounded join path.
# These are FE-level values only; observation data never enter SQL literals.
dbglm_databricks_array_terms = function(
  ops, info, coefficients,
  max_entries = getOption("dbreg.glm_databricks_array_max_entries", 65536L), max_bytes = 512 * 1024, conn
) {
  backend_info = detect_backend(conn)
  if (!is.numeric(max_entries) || length(max_entries) != 1L || !is.finite(max_entries) ||
    max_entries < 0 || max_entries != floor(max_entries)) {
    stop("GLM coefficient array limit must be a nonnegative whole number.", call. = FALSE)
  }
  if (max_entries == 0 || !length(ops)) return(NULL)
  ids = lapply(info, function(inf) as.numeric(inf[["id"]]))
  sizes = lengths(ids)
  if (length(ids) != length(ops) || any(!sizes) || length(coefficients) != sum(sizes)) {
    stop("Invalid FE coefficient dictionary.", call. = FALSE)
  }
  if (any(!is.finite(coefficients))) stop("Nonfinite projection coefficients.", call. = FALSE)
  if (any(vapply(ids, function(id) any(!is.finite(id) | id < 1 | id != floor(id)) || anyDuplicated(id) > 0L, logical(1)))) return(NULL)
  maxima = vapply(ids, max, numeric(1))
  if (any(maxima > 4 * sizes) || sum(maxima) > max_entries || any(maxima > .Machine[["integer.max"]])) return(NULL)
  starts = cumsum(c(0L, sizes))
  terms = character(length(ops))
  bytes = 0
  for (i in seq_along(ops)) {
    values = numeric(maxima[[i]])
    values[ids[[i]]] = coefficients[starts[[i]] + seq_len(sizes[[i]])]
    payload = paste0("[", paste(dbglm_number(values), collapse = ","), "]")
    encoded = raw_to_base64(charToRaw(enc2utf8(payload)))
    effect = paste0(
      "COALESCE(element_at(from_json(", sql_cast_text(conn, paste0("unbase64('", encoded, "')"), backend_info),
      ",'ARRAY<DOUBLE>'),CAST(s.", ops[[i]][["f"]], " AS INT)),0.0)"
    )
    v = ops[[i]][["v"]]
    terms[[i]] = if (v == "1.0") effect else paste0("s.", v, "*(", effect, ")")
    bytes = bytes + nchar(terms[[i]], type = "bytes")
    if (bytes > max_bytes) return(NULL)
  }
  terms
}

# One coefficient-frame/lookup contract serves prediction and projection.
# Decode integer64 keys before unlist() can strip their class.
dbglm_coefficient_frame = function(info, coefficients) {
  coefficients = as.matrix(coefficients)
  frame = data.frame(op = rep(seq_along(info), vapply(info, function(group) length(group[["id"]]), integer(1))))
  names(frame) = "__glm_op"
  ids = lapply(info, function(group) group[["id"]])
  frame[["__glm_group"]] = if (any(vapply(ids, is.character, logical(1)))) unlist(lapply(ids, as.character), use.names = FALSE) else
    unlist(lapply(ids, as.numeric), use.names = FALSE)
  frame[paste0("c", seq_len(ncol(coefficients)))] = as.data.frame(coefficients)
  frame
}

# `tables` lists physical objects owned by the fit; inline SQL owns nothing.
dbglm_coefficient_lookup = function(conn, frame, create, inline = FALSE) {
  if (isTRUE(detect_backend(conn)[["is_databricks"]])) {
    sql = sql_values_frame(conn, frame)
    if (is.null(sql)) projection_abort("upload_limit", "FE coefficient lookup exceeds the SQL transfer limit.")
    if (inline) return(list(source = paste0("(", sql, ")"), tables = character()))
    lookup = create(sql, materialize = TRUE)
  } else {
    registered = basename(tempfile("dbglm_coef_"))
    tryCatch(
      {
        duckdb::duckdb_register(conn, registered, frame)
        lookup = create(paste("SELECT * FROM", DBI::dbQuoteIdentifier(conn, registered)), materialize = TRUE)
      },
      finally = try(duckdb::duckdb_unregister(conn, registered), silent = TRUE)
    )
  }
  list(source = lookup, tables = lookup)
}

# Define a coefficient lookup once, then filter it per FE operation. Repeating
# an inline JSON lookup in every operation can exceed the warehouse query-text
# limit even when the lookup itself is within its transfer bound.
dbglm_coefficient_ctes = function(source, n_ops, prefix) {
  inline = startsWith(source, "(")
  lookup = if (inline) paste0(prefix, "_lookup") else source
  c(
    if (inline) paste0(lookup, " AS (SELECT * FROM ", source, ")"),
    vapply(seq_len(n_ops), function(i) paste0(
      prefix, "_g", i,
      " AS (SELECT * FROM ", lookup, " WHERE __glm_op=", i, ")"
    ), character(1))
  )
}
