# Databricks inlines CTEs, so a UNION branch per FE pair repeatedly plans the
# predictor. Keep three aggregation branches: globals, FE marginals and pairs.
# Pairs need only their weight sum, not every diagnostic/score aggregate.
dbglm_databricks_moments_sql = function(source, fe, value_terms, global_terms = character(), limit = NULL) {
  K = length(fe)
  # GROUPING_ID uses a signed BIGINT. Keep masks exactly representable in R.
  if (K < 2L || K > 30L) return(NULL)
  aliases = function(terms) sub(".* AS ([A-Za-z_0-9]+)$", "\\1", terms)
  value_names = aliases(value_terms)
  global_names = aliases(global_terms)
  nulls = function(names) if (length(names)) paste0("CAST(NULL AS DOUBLE) AS ", names) else character()
  grouped = function(groups, sets, terms) {
    masks = vapply(groups, function(g) sum(2^(K - setdiff(seq_len(K), g))), numeric(1))
    key = function(j) {
      cases = vapply(seq_along(groups), function(i) {
        g = groups[[i]]
        if (length(g) < j) return("")
        sprintf("WHEN %.0f THEN CAST(%s AS BIGINT)", masks[[i]], fe[g[[j]]])
      }, character(1))
      if (!any(nzchar(cases))) return("CAST(NULL AS BIGINT)")
      paste("CASE __gid", paste(cases, collapse = " "), "END")
    }
    fields = c(paste("CASE __gid", paste(sprintf("WHEN %.0f THEN %d", masks, sets), collapse = " "), "END AS __set"),
      paste(key(1L), "AS __i"), paste(key(2L), "AS __j"),
      nulls(global_names), aliases(terms), nulls(setdiff(value_names, aliases(terms))))
    # Keep the UNION schema in the same order even when pairs omit values.
    paste0("SELECT ", paste(fields, collapse = ","), " FROM (SELECT GROUPING_ID(", paste(fe, collapse = ","),
      ") AS __gid,", paste(c(fe, terms), collapse = ","), " FROM ", source, " GROUP BY GROUPING SETS (",
      paste(vapply(groups, function(g) paste0("(", paste(fe[g], collapse = ","), ")"), character(1)), collapse = ","), ")) grouped")
  }
  pairs = utils::combn(K, 2L, simplify = FALSE)
  branches = c(
    if (length(global_terms)) paste0("SELECT 0 AS __set,CAST(NULL AS BIGINT) AS __i,CAST(NULL AS BIGINT) AS __j,",
      paste(c(global_terms, nulls(value_names)), collapse = ","), " FROM ", source),
    grouped(as.list(seq_len(K)), seq_len(K), value_terms),
    grouped(pairs, K + seq_along(pairs), value_terms[1L]))
  sql = paste(branches, collapse = " UNION ALL ")
  # Keep stable set ordering for bounded preflight queries. The immutable FE
  # support lets IRLS reuse this bound without sorting every evaluation.
  if (!is.null(limit)) sql = paste0("SELECT * FROM (", sql, ") limited ORDER BY __set LIMIT ", sprintf("%.0f", limit))
  if (!length(global_terms)) return(sql)
  # Only the global row has these scalars. Avoid sending NULL copies on every
  # pair row; the reader keeps globals separate from the numeric FE moments.
  fields = paste(sprintf("'%s',%s", global_names, global_names), collapse = ",")
  paste0("SELECT __set,__i,__j,", paste(value_names, collapse = ","),
    ",CASE WHEN __set=0 THEN to_json(named_struct(", fields, ")) END AS __globals FROM (", sql, ") packed")
}

# Small, dense FE rank dictionaries can use constant array indexing instead
# of coefficient joins. Large/sparse dictionaries keep the bounded join path.
# These are FE-level values only; observation data never enter SQL literals.
dbglm_databricks_array_terms = function(ops, info, coefficients,
  max_entries = getOption("dbreg.glm_databricks_array_max_entries", 65536L), max_bytes = 512 * 1024) {
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
    payload = jsonlite::toJSON(values, digits = 17, auto_unbox = FALSE)
    encoded = gsub("[\r\n]", "", jsonlite::base64_enc(charToRaw(enc2utf8(payload))))
    effect = paste0("COALESCE(element_at(from_json(CAST(unbase64('", encoded,
      "') AS STRING),'ARRAY<DOUBLE>'),CAST(s.", ops[[i]][["f"]], " AS INT)),0.0)")
    v = ops[[i]][["v"]]
    terms[[i]] = if (v == "1.0") effect else paste0("s.", v, "*(", effect, ")")
    bytes = bytes + nchar(terms[[i]], type = "bytes")
    if (bytes > max_bytes) return(NULL)
  }
  terms
}
