inference_ssc = function(ssc) {
  numeric_df = !is.null(ssc) && is.numeric(ssc[["t.df"]])
  if (numeric_df) {
    df = ssc[["t.df"]]
    if (length(df) != 1L || !is.finite(df) || df < 1) stop("Numeric t.df must be at least one.", call. = FALSE)
    ssc[["t.df"]] = "min"
  }
  value = dbglm_ssc(ssc)
  if (numeric_df) value[["t.df"]] = df
  value
}

inference_cluster_terms = function(value) {
  if (is.character(value)) {
    if (!length(value) || anyNA(value) || any(!nzchar(value))) stop("Cluster columns must be nonempty names.", call. = FALSE)
    return(lapply(value, function(name) name))
  }
  if (!inherits(value, "formula")) stop("Clusters require column names or a formula.", call. = FALSE)
  if (length(value) == 3L && (!is.symbol(value[[2L]]) || !as.character(value[[2L]]) %in% c("cluster", "twoway"))) {
    stop("Unsupported covariance formula type.", call. = FALSE)
  }
  node = value[[length(value)]]
  split = function(node) {
    if (is.call(node) && identical(node[[1L]], as.name("+"))) return(c(split(node[[2L]]), split(node[[3L]])))
    list(node)
  }
  components = function(node) {
    if (is.symbol(node)) return(as.character(node))
    if (is.call(node) && identical(node[[1L]], as.name("^"))) return(unique(c(components(node[[2L]]), components(node[[3L]]))))
    stop("Cluster terms must be columns or typed column combinations with ^.", call. = FALSE)
  }
  lapply(split(node), components)
}

inference_vcov_request = function(object, vcov, cluster = NULL) {
  if (!is.null(cluster) && !is.null(vcov)) stop("Specify either vcov or cluster.", call. = FALSE)
  if (!is.null(cluster)) return(list(kind = "cluster", clusters = inference_cluster_terms(cluster)))
  if (inherits(vcov, "formula")) return(list(kind = "cluster", clusters = inference_cluster_terms(vcov)))
  if (!is.character(vcov) || length(vcov) != 1L || is.na(vcov)) stop("Unsupported covariance request.", call. = FALSE)
  kind = tolower(vcov)
  if (kind %in% c("iid", "standard")) return(list(kind = "iid", clusters = list()))
  if (kind %in% c("hc1", "hetero", "white")) return(list(kind = "hc1", clusters = list()))
  if (!kind %in% c("cluster", "twoway")) stop("Unsupported covariance type: ", vcov, ".", call. = FALSE)
  groups = object[["fitted_state"]][["inference"]][["fe_components"]]
  count = if (kind == "cluster") 1L else 2L
  if (length(groups) < count) stop("This covariance type requires explicit cluster columns.", call. = FALSE)
  list(kind = "cluster", clusters = groups[seq_len(count)])
}

inference_fe_nesting = function(conn, state, clusters) {
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  nested = rep(FALSE, length(state[["fe_components"]]))
  if (!length(nested)) return(nested)
  labels = vapply(state[["fe_components"]], paste, collapse = "^", character(1))
  cluster_labels = vapply(clusters, paste, collapse = "^", character(1))
  nested = vapply(seq_along(nested), function(i) labels[[i]] %in% cluster_labels ||
    any(state[["fe_components"]][[i]] %in% cluster_labels), logical(1))
  # Match fixest's K.exact=FALSE shortcut when every cluster is itself an FE.
  # Otherwise check the remaining functional dependencies on the fitted rows.
  if (all(nested) || all(cluster_labels %in% labels)) return(nested)
  queries = unlist(lapply(which(!nested), function(i) {
    vapply(clusters, function(cluster) paste0("SELECT ", i, " AS fe,MAX(n) AS maximum FROM (SELECT ",
      sql_distinct_tuple_count(conn, cluster), " AS n FROM ", q(state[["table"]]), " WHERE ", q(state[["estimated"]]),
      " GROUP BY ", paste(q(state[["fe_components"]][[i]]), collapse = ","), ") g"), character(1))
  }), use.names = FALSE)
  rows = dbglm_get_query(conn, paste(queries, collapse = " UNION ALL "))
  nested[as.integer(rows[["fe"]][as.numeric(rows[["maximum"]]) == 1])] = TRUE
  nested
}

inference_parameter_count = function(conn, state, request, ssc) {
  K = state[["K_model"]]
  sizes = state[["fe_sizes"]]
  fe_df = if (length(sizes)) sum(sizes) - length(sizes) + 1 else 0
  if (ssc[["K.fixef"]] == "none") return(K - fe_df)
  if (request[["kind"]] == "cluster" && ssc[["K.fixef"]] == "nonnested" && length(sizes)) {
    nested = inference_fe_nesting(conn, state, request[["clusters"]])
    K = K - sum(sizes[nested] - 1)
    K = max(K, length(state[["coefficient_names"]]) + 1L)
  }
  K
}

inference_repair_vcov = function(V, fix = TRUE) {
  if (!length(V)) return(V)
  if (any(!is.finite(V))) stop("Post-fit covariance is nonfinite.", call. = FALSE)
  e = eigen(V, symmetric = TRUE)
  if (all(e[["values"]] > 0)) return(V)
  repaired = tcrossprod(e[["vectors"]] %*% diag(pmax(e[["values"]], 1e-16), nrow(V)), e[["vectors"]])
  dimnames(repaired) = dimnames(V)
  if (max(abs(V - repaired)) > 1e-8) warning(if (fix)
    "The VCOV matrix is not positive definite and was 'fixed'." else "The VCOV matrix is not positive definite.", call. = FALSE)
  if (fix) repaired else V
}

retained_covariance = function(object, vcov = "iid", cluster = NULL, ssc = NULL, vcov_fix = TRUE) {
  state = object[["fitted_state"]][["inference"]]
  if (is.null(state)) stop("Post-fit covariance requires fitting with retain='all'.", call. = FALSE)
  if (!identical(state[["version"]], 1L)) stop("Unsupported retained-inference version.", call. = FALSE)
  if (!is.logical(vcov_fix) || length(vcov_fix) != 1L || is.na(vcov_fix)) stop("vcov_fix must be TRUE or FALSE.", call. = FALSE)
  conn = fitted_resource_connection(object[["resource_state"]])
  request = inference_vcov_request(object, if (is.null(cluster)) vcov else NULL, cluster)
  clusters = request[["clusters"]]
  if (length(clusters) > 2L || (inherits(object, "dbglm") && length(clusters) > 1L)) {
    stop("Supported clustering is one-/two-way OLS or one-way GLM.", call. = FALSE)
  }
  if (is.null(ssc)) ssc = state[["default_ssc"]]
  ssc = inference_ssc(ssc)
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  from = paste("FROM", q(state[["table"]]), "WHERE", q(state[["estimated"]]))
  if (length(clusters)) {
    columns = unique(unlist(clusters))
    classes = state[["column_classes"]]
    if (any(!columns %in% names(classes))) stop("Cluster columns are absent from the retained source.", call. = FALSE)
    missing = paste0(q(columns), " IS NULL")
    numeric = columns[classes[columns] %in% c("numeric", "integer")]
    if (length(numeric)) missing = c(missing, paste0("isnan(CAST(", q(numeric), " AS DOUBLE))"))
    bad = dbglm_get_query(conn, paste("SELECT COUNT(*) AS n", from, "AND (", paste(missing, collapse = " OR "), ")"))
    if (as.numeric(bad[["n"]]) > 0) stop("Post-fit clusters have missing values on the fitted sample; the sample cannot change.", call. = FALSE)
  }
  n = state[["nobs"]]
  K = inference_parameter_count(conn, state, request, ssc)
  if (ssc[["K.adj"]] && n <= K) stop("Insufficient residual degrees of freedom.", call. = FALSE)
  bread = state[["bread"]]
  G = numeric()
  if (request[["kind"]] == "iid") {
    if (state[["kind"]] == "ols") V = bread*state[["rss"]]/(if (ssc[["K.adj"]]) n-K else n-1) else if (state[["kind"]] == "gaussian")
      V = bread*state[["dispersion_scaled"]]*(if (ssc[["K.adj"]]) (n-1)/(n-K) else 1) else
      V = bread/state[["weight_scale"]]*(if (ssc[["K.adj"]]) (n-1)/(n-K) else 1)
  } else if (request[["kind"]] == "hc1") {
    V = score_moments(conn, from, q(state[["scores"]]), projection = bread)[["matrix"]]
    if (ssc[["K.adj"]]) V = V*n/(n-K)
  } else {
    components = c(clusters, if (length(clusters) == 2L) list(unique(unlist(clusters))))
    moments = lapply(components, function(columns) score_moments(conn, from, q(state[["scores"]]), q(columns), projection = bread))
    G = vapply(moments, `[[`, numeric(1), "groups")
    if (any(G < 2)) stop("Clustered covariance requires at least two clusters.", call. = FALSE)
    meat = bread*0
    for (i in seq_along(moments)) {
      adjustment = if (!ssc[["G.adj"]]) 1 else if (ssc[["G.df"]] == "min") min(G[seq_along(clusters)])/(min(G[seq_along(clusters)])-1) else G[[i]]/(G[[i]]-1)
      meat = meat + (if (i > length(clusters)) -1 else 1)*adjustment*moments[[i]][["matrix"]]
    }
    V = meat
    if (ssc[["K.adj"]]) V = V*(n-1)/(n-K)
    G = G[seq_along(clusters)]
  }
  V = inference_repair_vcov((V+t(V))/2, vcov_fix)
  dimnames(V) = rep(list(state[["coefficient_names"]]), 2L)
  attr(V, "type") = request[["kind"]]
  attr(V, "n_params") = attr(V, "df.K") = K
  attr(V, "ssc") = ssc
  attr(V, "rss") = state[["rss"]]
  if (!is.null(attr(object[["vcov"]], "tss"))) attr(V, "tss") = attr(object[["vcov"]], "tss")
  attr(V, "df.t") = if (is.numeric(ssc[["t.df"]])) ssc[["t.df"]] else if (length(G) && ssc[["t.df"]] == "min") max(min(G)-1, 1) else max(n-K, 1)
  if (length(G)) attr(V, "n_clusters") = G
  V
}
