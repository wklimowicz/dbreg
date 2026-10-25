# Fisher information, score meat and fixest-compatible GLM corrections.
dbglm_ssc = function(ssc) {
  defaults = list(K.adj = TRUE, K.fixef = "nonnested", K.exact = FALSE,
    G.adj = TRUE, G.df = "min", t.df = "min")
  if (!is.null(ssc)) {
    if (!is.list(ssc) || is.null(names(ssc)) || anyDuplicated(names(ssc)) ||
      any(!names(ssc) %in% names(defaults))) stop("Invalid `ssc` list.", call. = FALSE)
    defaults[names(ssc)] = ssc
  }
  for (name in c("K.adj", "G.adj", "K.exact")) {
    if (!is.logical(defaults[[name]]) || length(defaults[[name]]) != 1 || is.na(defaults[[name]])) {
      stop("`ssc` adjustment flags must be TRUE or FALSE.", call. = FALSE)
    }
  }
  for (name in c("K.fixef", "G.df", "t.df")) {
    if (!is.character(defaults[[name]]) || length(defaults[[name]]) != 1 || is.na(defaults[[name]])) {
      stop("`ssc` policies must be single character strings.", call. = FALSE)
    }
  }
  if (defaults[["K.exact"]] || !defaults[["K.fixef"]] %in% c("full", "none", "nonnested") ||
    !defaults[["G.df"]] %in% c("min", "conventional") || !defaults[["t.df"]] %in% c("min", "conventional")) {
    stop("Unsupported `ssc` settings.", call. = FALSE)
  }
  defaults
}

dbglm_covariance = function(conn, table, x, fit, family, kind, ssc, n) {
  inverse = dbglm_solve(fit[["H"]])[["inverse"]]
  p = length(x)
  if (ssc[["K.adj"]] && n <= p) stop("Insufficient residual degrees of freedom.", call. = FALSE)
  if (kind == "iid") {
    V = inverse
    if (ssc[["K.adj"]]) V = V * (n - 1) / (n - p)
    return(list(matrix = V, clusters = NULL))
  }
  cte = dbglm_working_sql(table, x, fit[["beta"]], family[["link"]])
  if (kind == "hc1") {
    row = backend_get_query(conn, paste0(cte, "SELECT ",
      paste(dbglm_moment_terms(x, "score*score"), collapse = ", "), " FROM working"))
    G = NULL
  } else {
    scores = sprintf("SUM((%s)*score) AS s_%d", x, seq_along(x))
    terms = dbglm_moment_terms(paste0("s_", seq_along(x)), "1.0")
    row = backend_get_query(conn, paste0(cte,
      ", cluster_scores AS (SELECT __dbglm_cluster, ", paste(scores, collapse = ", "),
      " FROM working GROUP BY __dbglm_cluster) SELECT COUNT(*) AS n_clusters, ",
      paste(terms, collapse = ", "), " FROM cluster_scores"))
    G = as.numeric(row[["n_clusters"]])
    if (G < 2) stop("Clustered covariance requires at least two clusters.", call. = FALSE)
  }
  meat = dbglm_read_information(row, p)
  V = inverse %*% meat %*% inverse
  if (ssc[["K.adj"]]) V = V * (if (kind == "hc1") n else n - 1) / (n - p)
  if (kind == "cluster" && ssc[["G.adj"]]) V = V * G / (G - 1)
  list(matrix = V, clusters = G)
}

dbglm_fe_covariance = function(conn, x, fit, kind, ssc, n, counts) {
  p = length(x)
  K = if (kind == "cluster") counts[["covariance"]] else counts[["covariance"]]
  if (ssc[["K.adj"]] && n <= K) stop("Insufficient residual degrees of freedom.", call. = FALSE)
  inverse = if (p) dbglm_solve(fit[["H"]])[["inverse"]] else matrix(numeric(), 0, 0)
  G = NULL
  if (kind == "cluster" && !p) {
    G = as.numeric(backend_get_query(conn, paste("SELECT COUNT(DISTINCT __dbglm_cluster) AS n FROM", fit[["working"]]))[["n"]])
    if (G < 2) stop("Clustered covariance requires at least two clusters.", call. = FALSE)
  }
  if (!p) return(list(matrix = inverse, clusters = G, K = K))
  if (kind == "iid") V = inverse else {
    if (!is.null(fit[["raw_projected"]])) {
      from = paste0(" FROM ", fit[["raw_projected"]], " b")
      x_cols = paste0("b.", x, "_tilde")
    } else {
      from = paste0(" FROM ", fit[["working"]], " b JOIN ", fit[["projected"]], " p USING (__dbglm_key)")
      x_cols = paste0("p.", x)
    }
    if (kind == "hc1") {
      sql = paste0("SELECT ", paste(dbglm_moment_terms(x_cols, "b.score*b.score"), collapse = ","), from)
      row = backend_get_query(conn, sql)
    } else {
      scores = sprintf("SUM(%s*b.score) AS s%d", x_cols, seq_along(x))
      sql = paste0("SELECT COUNT(*) AS __n_clusters,", paste(dbglm_moment_terms(paste0("s", seq_along(x)), "1.0"), collapse = ","), " FROM (SELECT ",
        paste(scores, collapse = ","), from, " GROUP BY __dbglm_cluster) t")
      row = backend_get_query(conn, sql)
      G = as.numeric(row[["__n_clusters"]])
      if (G < 2) stop("Clustered covariance requires at least two clusters.", call. = FALSE)
    }
    meat = dbglm_read_information(row, p)
    V = inverse %*% meat %*% inverse
  }
  if (ssc[["K.adj"]]) V = V * (if (kind == "hc1") n else n - 1) / (n - K)
  if (kind == "cluster" && ssc[["G.adj"]]) V = V * G / (G - 1)
  list(matrix = V, clusters = G, K = K)
}
