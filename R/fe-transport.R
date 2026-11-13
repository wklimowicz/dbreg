# One closed protocol: marginal (id,h,RHS) and pair (id,id,h) summaries.
# Typed source keys remain in SQL; only bounded dense ordinal keys cross DBI.
fe_summary_sql = function(conn, source, fe, columns, max_rows) {
  K = length(fe)
  pairs = if (K > 1L) utils::combn(K, 2L, simplify = FALSE) else list()
  backend = detect_backend(conn)
  sum_sql = function(x) sql_stable_sum(conn, x, backend_info = backend)
  values = paste0("b", seq_along(columns))
  terms = c(paste(sum_sql("w"), "AS h"), paste(sum_sql(paste0("w*", columns)), "AS", values))
  if (isTRUE(detect_backend(conn)[["is_databricks"]])) {
    grouped = function(groups, sets, rhs) {
      masks = vapply(groups, function(g) sum(2^(K - setdiff(seq_len(K), g))), numeric(1))
      key = function(j) {
        cases = vapply(seq_along(groups), function(i) {
          g = groups[[i]]
          if (length(g) < j) return("")
          sprintf("WHEN %.0f THEN CAST(%s AS BIGINT)", masks[i], fe[g[j]])
        }, character(1))
        if (!any(nzchar(cases))) return("CAST(NULL AS BIGINT)")
        paste("CASE gid", paste(cases, collapse = " "), "END")
      }
      selected = c(paste("CASE gid", paste(sprintf("WHEN %.0f THEN %d", masks, sets), collapse = " "), "END AS kind"),
        paste(key(1L), "AS i"), paste(key(2L), "AS j"), "h",
        if (rhs) values else paste("CAST(NULL AS DOUBLE) AS", values))
      paste0("SELECT ", paste(selected, collapse = ","), " FROM (SELECT GROUPING_ID(", paste(fe, collapse = ","), ") AS gid,",
        paste(c(fe, if (rhs) terms else terms[1L]), collapse = ","), " FROM fe_source GROUP BY GROUPING SETS (",
        paste(vapply(groups, function(g) paste0("(", paste(fe[g], collapse = ","), ")"), character(1)), collapse = ","), ")) grouped")
    }
    branches = c(grouped(as.list(seq_len(K)), seq_len(K), TRUE), if (length(pairs)) grouped(pairs, K + seq_along(pairs), FALSE))
    body = paste0("WITH fe_source AS (SELECT ", paste(c(fe, columns, "w"), collapse = ","), " FROM ", source, ") ", paste(branches, collapse = " UNION ALL "))
    limited = paste0("SELECT * FROM (", body, ") bounded LIMIT ", sprintf("%.0f", max_rows + 1))
    return(fe_packet_sql(conn, limited, K, length(columns)))
  }
  branch = function(kind, i, j = NULL) {
    groups = c(i, j)
    paste0("SELECT ", kind, " AS kind,CAST(", i, " AS DOUBLE) AS i,",
      if (is.null(j)) "CAST(NULL AS DOUBLE)" else paste0("CAST(", j, " AS DOUBLE)"), " AS j,", sum_sql("w"), " AS h,",
      if (is.null(j)) paste0("list_value(", paste(sum_sql(paste0("w*", columns)), collapse = ","), ")") else "CAST(NULL AS DOUBLE[])",
      " AS payload FROM fe_source GROUP BY ", paste(groups, collapse = ","))
  }
  branches = c(vapply(seq_len(K), function(i) branch(i, fe[i]), character(1)),
    vapply(seq_along(pairs), function(k) branch(K + k, fe[pairs[[k]][1L]], fe[pairs[[k]][2L]]), character(1)))
  paste0("WITH fe_source AS MATERIALIZED (SELECT ", paste(c(fe, columns, "w"), collapse = ","), " FROM ", source, ") ",
    "SELECT * FROM (", paste(branches, collapse = " UNION ALL "), ") bounded LIMIT ", sprintf("%.0f", max_rows + 1))
}

fe_packet_sql = function(conn, sql, K, width, packet_bytes = getOption("dbreg.databricks_packet_bytes", 20000L)) {
  text = function(x) paste0("COALESCE(CAST(", x, " AS STRING),'NaN')")
  marginal = c("CAST(i AS STRING)", text("h"), text(paste0("b", seq_len(width))))
  pair = c("CAST(i AS STRING)", "CAST(j AS STRING)", text("h"))
  entries = paste0("SELECT kind,i,j,CASE WHEN kind<=", K, " THEN concat_ws(',',", paste(marginal, collapse = ","),
    ") ELSE concat_ws(',',", paste(pair, collapse = ","), ") END AS entry FROM (", sql, ") summaries")
  fe_packetize_sql(entries, packet_bytes)
}

# Packetization is separate from numeric field formatting so byte-boundary
# checks can exercise wide valid vectors without a huge SQL projection schema.
fe_packetize_sql = function(entries, packet_bytes = getOption("dbreg.databricks_packet_bytes", 20000L)) {
  if (!is.numeric(packet_bytes) || length(packet_bytes) != 1L || !is.finite(packet_bytes) || packet_bytes != floor(packet_bytes) || packet_bytes < 1 || packet_bytes > 32000) {
    stop("dbreg.databricks_packet_bytes must be a whole number from 1 to 32000.", call. = FALSE)
  }
  # Entries larger than the target get their own packet. Smaller entries use
  # running byte buckets whose payload is <2*target, hence below ODBC's 65,535
  # byte ceiling. Even a single oversized entry is flagged, never truncated.
  costed = paste0("SELECT *,OCTET_LENGTH(entry) AS bytes FROM (", entries, ") entries")
  buckets = paste0("SELECT *,FLOOR(SUM(CASE WHEN bytes<=", packet_bytes, " THEN bytes+1 ELSE 0 END) OVER (PARTITION BY kind ORDER BY i,j ",
    "ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW)/", packet_bytes, ") AS bucket,",
    "CASE WHEN bytes>", packet_bytes, " THEN i ELSE 0 END AS single_i,CASE WHEN bytes>", packet_bytes, " THEN COALESCE(j,0) ELSE 0 END AS single_j FROM (", costed, ") costed")
  paste0("SELECT kind,COUNT(*) AS n,SUM(bytes+1)-1 AS bytes,MAX(CASE WHEN bytes>=65535 THEN 1 ELSE 0 END) AS overflow,",
    "array_join(collect_list(CASE WHEN bytes<65535 THEN entry END),';') AS payload FROM (", buckets,
    ") packets GROUP BY kind,bucket,single_i,single_j")
}

# Decode each bounded packet into preallocated numeric arrays. No giant split
# of every payload, wide NULL-padded rectangle or source-key string dictionary.
fe_decode_summaries = function(result, sizes, width, budget, limits, packets = FALSE) {
  required = if (packets) c("kind", "n", "bytes", "overflow", "payload") else c("kind", "i", "j", "h", "payload")
  if (!all(required %in% names(result))) stop("Incomplete fixed-effect summary schema.", call. = FALSE)
  K = length(sizes)
  combinations = if (K > 1L) utils::combn(K, 2L, simplify = FALSE) else list()
  count = if (packets) as.numeric(result[["n"]]) else rep(1, nrow(result))
  sets = as.numeric(result[["kind"]])
  if (length(sets) != length(count) || any(!is.finite(count) | count < 1 | count != floor(count)) ||
      any(!is.finite(sets) | sets < 1 | sets > K + length(combinations) | sets != floor(sets))) stop("Invalid fixed-effect summary metadata.", call. = FALSE)
  cells = sum(count[sets<=K]) * (2 + width) + sum(count[sets>K]) * 3 + max(3, 2 + width)
  if (sum(count) + budget[["rows"]] > limits[["summary_rows"]] || cells + budget[["cells"]] > limits[["transfer_cells"]]) {
    projection_abort("transfer_limit", "Fixed-effect summaries exceed the transfer budget.", list(rows = sum(count), cells = cells, limits = limits))
  }
  if (packets) {
    flags = as.numeric(result[["overflow"]])
    if (any(!is.finite(flags) | !flags %in% c(0,1))) stop("Invalid fixed-effect packet overflow metadata.", call. = FALSE)
    if (any(flags != 0)) projection_abort("transfer_limit", "A fixed-effect summary exceeds the ODBC cell byte limit.")
  }
  sizes_by_set = vapply(seq_len(K + length(combinations)), function(set) sum(count[sets==set]), numeric(1))
  if (!identical(as.numeric(sizes_by_set[seq_len(K)]), as.numeric(sizes)) || any(sizes_by_set <= 0)) stop("Fixed-effect summary support changed.", call. = FALSE)
  state = list2env(list(
    info = lapply(sizes, function(size) list(id = seq_len(size), h = numeric(size), b = matrix(0, size, width))),
    seen = lapply(sizes, function(size) logical(size)),
    pairs = lapply(seq_along(combinations), function(k) {
      pair = combinations[[k]]; n = sizes_by_set[K + k]
      list(i = pair[2L], j = pair[1L], ir = integer(n), jr = integer(n), h = numeric(n))
    }),
    positions = integer(length(combinations))
  ), parent = emptyenv())
  number = function(x) {
    valid = grepl("^([+-]?([0-9]+(\\.[0-9]*)?|\\.[0-9]+)([eE][+-]?[0-9]+)?|NaN|Infinity|-Infinity)$", x)
    if (!all(valid)) stop("Invalid numeric fixed-effect packet.", call. = FALSE)
    suppressWarnings(as.numeric(x))
  }
  absorb = function(set, values) {
    keys = values[, if (set<=K) 1L else 1:2, drop = FALSE]
    bounds = if (set<=K) sizes[set] else sizes[combinations[[set-K]]]
    if (any(!is.finite(keys) | keys < 1 | keys != floor(keys)) || any(sweep(keys, 2L, bounds, ">"))) stop("Invalid fixed-effect ordinal keys.", call. = FALSE)
    if (any(!is.finite(values)) || any(values[, if (set<=K) 2L else 3L] <= 0)) projection_abort("singular", "Nonfinite or nonpositive fixed-effect moments.")
    if (set<=K) {
      key = as.integer(keys[,1L])
      if (anyDuplicated(key) || any(state[["seen"]][[set]][key])) stop("Duplicate fixed-effect marginal keys.", call. = FALSE)
      state[["seen"]][[set]][key] = TRUE
      state[["info"]][[set]][["h"]][key] = values[,2L]
      state[["info"]][[set]][["b"]][key,] = values[,-(1:2),drop=FALSE]
    } else {
      k = set-K
      where = state[["positions"]][k] + seq_len(nrow(values))
      state[["pairs"]][[k]][["ir"]][where] = as.integer(keys[,2L])
      state[["pairs"]][[k]][["jr"]][where] = as.integer(keys[,1L])
      state[["pairs"]][[k]][["h"]][where] = values[,3L]
      state[["positions"]][k] = state[["positions"]][k] + nrow(values)
    }
  }
  if (packets) {
    for (row in seq_len(nrow(result))) {
      set = sets[row]
      payload = result[["payload"]][row]
      bytes = as.numeric(result[["bytes"]][row])
      if (is.na(payload) || !is.finite(bytes) || bytes < 1 || bytes >= 65535 || nchar(payload, type = "bytes") != bytes) stop("Fixed-effect packet byte count is incomplete.", call. = FALSE)
      fields = strsplit(payload, "[;,]")[[1L]]
      fields_per_row = if (set<=K) 2 + width else 3L
      if (length(fields) != count[row] * fields_per_row) stop("Fixed-effect packet field count is incomplete.", call. = FALSE)
      absorb(set, matrix(number(fields), ncol = fields_per_row, byrow = TRUE))
    }
  } else {
    # One block per summary set; each row is one marginal or pair cell.
    for (set in unique(sets)) {
      rows = which(sets == set)
      values = if (set<=K) {
        payload = lapply(result[["payload"]][rows], as.numeric)
        if (any(lengths(payload) != width)) stop("Inconsistent fixed-effect marginal width.", call. = FALSE)
        cbind(as.numeric(result[["i"]][rows]), as.numeric(result[["h"]][rows]),
          matrix(unlist(payload, use.names = FALSE), ncol = width, byrow = TRUE))
      } else cbind(as.numeric(result[["i"]][rows]), as.numeric(result[["j"]][rows]), as.numeric(result[["h"]][rows]))
      absorb(set, values)
    }
  }
  info = state[["info"]]
  pairs = state[["pairs"]]
  for (pair in pairs) {
    ordering = order(pair[["ir"]], pair[["jr"]])
    if (any(diff(pair[["ir"]][ordering]) == 0 & diff(pair[["jr"]][ordering]) == 0)) stop("Duplicate fixed-effect pair keys.", call. = FALSE)
  }
  budget[["rows"]] = budget[["rows"]] + sum(count)
  budget[["cells"]] = budget[["cells"]] + cells
  list(info = info, pairs = pairs, wire_rows = nrow(result), wire_bytes = if (packets) sum(as.numeric(result[["bytes"]])) else NULL)
}

query_fe_summaries = function(conn, source, fe, columns, sizes, budget, limits) {
  marginal = sum(sizes) * (2 + length(columns))
  room = limits[["transfer_cells"]] - budget[["cells"]] - marginal - max(3, 2 + length(columns))
  max_rows = min(limits[["summary_rows"]] - budget[["rows"]], sum(sizes) + floor(room/3), .Machine[["integer.max"]] - 1L)
  if (room < 0 || max_rows < sum(sizes)) projection_abort("transfer_limit", "Fixed-effect summaries exceed the transfer budget.")
  sql = fe_summary_sql(conn, source, fe, columns, max_rows)
  result = DBI::dbGetQuery(conn, sql)
  counts = if (isTRUE(detect_backend(conn)[["is_databricks"]])) as.numeric(result[["n"]]) else rep(1, nrow(result))
  if (any(!is.finite(counts) | counts < 1) || sum(counts) > max_rows) projection_abort("transfer_limit", "Fixed-effect summary row cap exceeded.")
  decoded = fe_decode_summaries(result, sizes, length(columns), budget, limits, packets = isTRUE(detect_backend(conn)[["is_databricks"]]))
  decoded[["statements"]] = 1L
  decoded[["transport"]] = if (isTRUE(detect_backend(conn)[["is_databricks"]])) "packets" else "lists"
  decoded
}
