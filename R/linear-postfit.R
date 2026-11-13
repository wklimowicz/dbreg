linear_cluster_terms = function(value) {
  if (is.character(value) && length(value) && !anyNA(value) && all(nzchar(value))) {
    return(as.list(value))
  }
  if (!inherits(value, "formula")) stop("Clusters require column names or a formula.", call. = FALSE)
  if (length(value) == 3L &&
      (!is.symbol(value[[2L]]) || !as.character(value[[2L]]) %in% c("cluster", "twoway"))) {
    stop("Unsupported covariance formula type.", call. = FALSE)
  }
  split = function(node) {
    if (is.call(node) && identical(node[[1L]], as.name("+"))) return(c(split(node[[2L]]), split(node[[3L]])))
    list(node)
  }
  components = function(node) {
    if (is.symbol(node)) return(as.character(node))
    if (is.call(node) && identical(node[[1L]], as.name("^"))) {
      return(unique(c(components(node[[2L]]), components(node[[3L]]))))
    }
    stop("Cluster terms must be columns or native tuples joined with ^.", call. = FALSE)
  }
  lapply(split(value[[length(value)]]), components)
}

prepare_retained_covariance = function(object, vcov, cluster, ssc, source_table, source_keys,
                                       source_classes, default_cluster, families, sizes,
                                       estimated_col, allow_twoway = TRUE, native_map = NULL) {
  conn = retained_fit_connection(object[["resource_state"]])
  if (!is.null(vcov) && !is.null(cluster)) stop("Specify either vcov or cluster.", call. = FALSE)
  V0 = object[["vcov"]]
  groups = if (!is.null(cluster)) linear_cluster_terms(cluster) else if (inherits(vcov, "formula")) {
    if (!allow_twoway && length(vcov) == 3L && identical(vcov[[2L]], as.name("twoway"))) {
      stop("Post-fit GLM covariance accepts one cluster term.", call. = FALSE)
    }
    linear_cluster_terms(vcov)
  } else list()
  allowed = if (allow_twoway) c("iid", "hc1", "cluster", "twoway") else c("iid", "hc1", "cluster")
  kind = if (length(groups)) "cluster" else if (is.null(vcov)) attr(V0, "type") else match.arg(vcov, allowed)
  if (is.null(vcov) && kind == "cluster" && !length(groups) && length(default_cluster)) {
    groups = if (is.list(default_cluster)) default_cluster else as.list(default_cluster)
  }
  if (allow_twoway && kind %in% c("cluster", "twoway") && !length(groups)) {
    needed = if (kind == "cluster") 1L else 2L
    if (length(families) < needed) stop("This covariance type requires cluster columns.", call. = FALSE)
    groups = families[seq_len(needed)]
    kind = "cluster"
  } else if (!allow_twoway && kind == "cluster" && !length(groups) && length(families)) {
    groups = families[1L]
  }
  max_groups = if (allow_twoway) 2L else 1L
  if (allow_twoway && (length(groups) > 2L || anyDuplicated(groups))) {
    stop("Post-fit clustering accepts one or two distinct cluster terms.", call. = FALSE)
  }
  if (!allow_twoway && kind == "cluster" && length(groups) != 1L) {
    stop("Post-fit GLM covariance requires one cluster term.", call. = FALSE)
  }
  if (is.null(ssc)) ssc = attr(V0, "ssc")
  ssc = normalize_ssc(ssc)
  rows = object[["retained_rows"]]
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  cluster_columns = unique(unlist(groups, use.names = FALSE))
  if (length(cluster_columns) && !all(cluster_columns %in% names(source_classes))) {
    stop("Cluster columns are absent from the retained source.", call. = FALSE)
  }
  from = paste("FROM", rows[["table"]], "s")
  if (length(groups)) {
    from = paste(from, "JOIN", source_table, "b ON", paste(
      paste0("s.", rows[["keys"]], "=b.", q(source_keys)), collapse = " AND "))
  }
  from = paste(from, "WHERE", paste0("s.", estimated_col))
  native = function(x) if (is.null(native_map)) x else unname(native_map[x])
  if (length(cluster_columns)) {
    cols = paste0("b.", q(native(cluster_columns)))
    numeric = source_classes[cluster_columns] %in% c("integer", "numeric")
    bad = c(paste(cols, "IS NULL"), if (any(numeric)) paste0("NOT ", sql_is_finite(conn, cols[numeric])))
    bad_sql = paste("SELECT COUNT(*) AS n", from, "AND (", paste(bad, collapse = " OR "), ")")
    if (DBI::dbGetQuery(conn, bad_sql)[["n"]] > 0) {
      stop("Post-fit clusters have missing or nonfinite values on the fitted sample.", call. = FALSE)
    }
  }
  nested_df = 0
  if (length(sizes) && length(groups) && ssc[["K.fixef"]] == "nonnested") {
    native_families = if (is.null(native_map)) families else lapply(families, native)
    native_groups = if (is.null(native_map)) groups else lapply(groups, native)
    nested_df = retained_nested_fe_df(conn, from, native_families, sizes, native_groups)
  }
  list(conn = conn, V0 = V0, groups = groups, kind = kind, ssc = ssc, q = q,
    cluster_columns = cluster_columns, from = from, nested_df = nested_df, native = native)
}

retained_linear_covariance = function(object, vcov, cluster, ssc, vcov_fix) {
  state = object[["retained_inference"]]
  if (is.null(state)) stop("Post-fit covariance requires retain = 'inference' or 'all'.", call. = FALSE)
  if (!is.logical(vcov_fix) || length(vcov_fix) != 1L || is.na(vcov_fix)) {
    stop("`vcov_fix` must be TRUE or FALSE.", call. = FALSE)
  }
  prep = prepare_retained_covariance(object, vcov, cluster, ssc,
    state[["original"]], state[["source_keys"]], state[["source_classes"]],
    state[["default_cluster"]], object[["fe_components"]], state[["fe_sizes"]],
    "__dbreg_estimated", allow_twoway = TRUE)
  conn = prep[["conn"]]
  V0 = prep[["V0"]]
  groups = prep[["groups"]]
  kind = prep[["kind"]]
  ssc = prep[["ssc"]]
  q = prep[["q"]]
  cluster_columns = prep[["cluster_columns"]]
  from = prep[["from"]]
  n = nobs(object)
  df = object[["df_residual"]]
  sizes = state[["fe_sizes"]]
  K = n - df
  if (ssc[["K.fixef"]] == "none") K = K - fe_intercept_df(sizes)
  if (length(sizes) && length(groups) && ssc[["K.fixef"]] == "nonnested") {
    K = K - prep[["nested_df"]]
  }
  rss = attr(V0, "rss")
  public_bread = state[["bread"]]
  coordinates = object[["prediction_coordinates"]]
  centered = !is.null(coordinates[["bread"]]) && !is.null(coordinates[["transform"]])
  if (centered) {
    # Retained scores already use fitting coordinates. Form their covariance
    # directly; undoing a huge public-origin transform loses interval precision.
    state[["bread"]] = coordinates[["bread"]]
    state[["sandwich_bread"]] = coordinates[["bread"]]
  }
  bread = state[["bread"]]
  if (kind == "iid") {
    V = compute_vcov("iid", "moments", bread, rss, df, n, n_params = K, ssc = ssc)
  } else {
    scores = if (length(state[["design"]])) paste0("s.", state[["design"]], "*s.__dbreg_saved_residual*",
      "(s.__dbreg_saved_weight/(", dbglm_number(state[["weight_scale"]]), "))") else character()
    components = c(groups, if (length(groups) == 2L) list(cluster_columns))
    meats = if (!length(components)) list(retained_score_cross(conn, from, scores)) else
      lapply(components, function(keys) retained_score_cross(conn, from, scores, paste0("b.", q(keys))))
    if (length(groups) == 2L) {
      G = vapply(meats, `[[`, numeric(1), "groups")
      if (any(G < 2)) stop("Clustered covariance requires at least two clusters.", call. = FALSE)
      adjustment = if (!ssc[["G.adj"]]) rep(1, 3L) else if (ssc[["G.df"]] == "min") {
        rep(min(G[1:2]) / (min(G[1:2]) - 1), 3L)
      } else G / (G - 1)
      meat = meats[[1L]][["matrix"]] * adjustment[1L] +
        meats[[2L]][["matrix"]] * adjustment[2L] - meats[[3L]][["matrix"]] * adjustment[3L]
      attr(meat, "n_clusters") = min(G[1:2])
      covariance_ssc = ssc
      covariance_ssc[["G.adj"]] = FALSE
      V = retained_score_vcov(kind, meat, state, rss, df, n, K, covariance_ssc)
      attr(V, "ssc") = ssc
      attr(V, "n_clusters") = G[1:2]
      attr(V, "type") = "twoway"
    } else {
      meat = meats[[1L]][["matrix"]]
      if (length(groups)) attr(meat, "n_clusters") = meats[[1L]][["groups"]]
      V = retained_score_vcov(kind, meat, state, rss, df, n, K, ssc)
    }
  }
  attr(V, "rss") = rss
  attr(V, "tss") = attr(V0, "tss")
  if (centered) {
    prediction = V
    internal_names = dimnames(V)
    transform = coordinates[["transform"]]
    fields = c("type", "n_params", "df.K", "ssc", "df.t", "n_clusters", "repair", "rss", "tss")
    metadata = lapply(fields, function(name) attr(V, name))
    V = transform %*% V %*% t(transform)
    V = (V + t(V)) / 2
    dimnames(V) = dimnames(public_bread)
    if (length(groups) == 2L) {
      # Repair the reported coefficients, as fixest does, and give intervals
      # the same covariance back in the fitting coordinates.
      visible = linear_reported_positions(object, nrow(V))
      block = retained_vcov_psd(V[visible, visible, drop = FALSE], vcov_fix)
      metadata[match("repair", fields)] = list(attr(block, "repair"))
      # Intervals combine slopes and fixed effects, so they use the repair of
      # the full public covariance (fixest with FE dummies) in fitting coordinates.
      full = if (vcov_fix) suppressWarnings(retained_vcov_psd(V, TRUE)) else V
      inverse = if (nrow(transform) == ncol(transform)) tryCatch(solve(transform), error = function(e) NULL)
      if (!identical(full, V)) {
        prediction = if (!is.null(inverse)) {
          value = inverse %*% full %*% t(inverse)
          value = (value + t(value)) / 2
          # A far-from-origin transform can lose the repair; repair in place then.
          values = eigen(value, symmetric = TRUE, only.values = TRUE)[["values"]]
          if (!all(is.finite(values)) || min(values) < -1e-8 * max(abs(values))) {
            suppressWarnings(retained_vcov_psd(prediction, TRUE))
          } else value
        } else suppressWarnings(retained_vcov_psd(prediction, TRUE))
        dimnames(prediction) = internal_names
      }
      if (isTRUE(attr(block, "repair")[["repaired"]])) V[visible, visible] = block
    }
    for (i in seq_along(fields)) if (!is.null(metadata[[i]])) attr(V, fields[[i]]) = metadata[[i]]
    attr(V, "prediction_covariance") = prediction
  }
  V
}

retained_score_cross = function(conn, from, scores, clusters = character()) {
  names = if (length(scores)) paste0("q", seq_along(scores)) else character()
  values = if (length(clusters)) sql_stable_sum(conn, scores) else scores
  terms = c("COUNT(*) AS n_groups", dbglm_moment_terms(names, "1.0"))
  projection = if (length(names)) paste(paste(values, "AS", names), collapse = ",") else "1 AS dummy"
  sql = paste0("SELECT ", paste(terms, collapse = ","),
    " FROM (SELECT ", projection, " ", from,
    if (length(clusters)) paste(" GROUP BY", paste(clusters, collapse = ",")),
    ") retained_score_groups")
  row = DBI::dbGetQuery(conn, sql)
  list(matrix = dbglm_read_information(row, length(scores)), groups = as.numeric(row[["n_groups"]]))
}

retained_score_vcov = function(kind, meat, state, rss, df, n, K, ssc) {
  attr(meat, "weight_scale") = state[["weight_scale"]]
  compute_vcov(kind, "moments", state[["bread"]], rss, df, n, n_params = K, meat = meat,
    sandwich_bread = state[["sandwich_bread"]], ssc = ssc)
}

retained_nested_fe_df = function(conn, from, families, sizes, groups) {
  q = function(x) paste0("b.", as.character(DBI::dbQuoteIdentifier(conn, x)))
  nested = vapply(families, function(fe) any(vapply(groups,
    function(cl) all(cl %in% fe), logical(1))), logical(1))
  if (!all(nested) && !all(vapply(groups, function(cl) any(vapply(families,
      function(fe) identical(cl, fe), logical(1))), logical(1)))) {
    for (i in which(!nested)) for (cl in groups) {
      fields = paste(q(unique(c(families[[i]], cl))), collapse = ",")
      family = paste(as.character(DBI::dbQuoteIdentifier(conn, families[[i]])), collapse = ",")
      sql = paste0("SELECT MAX(n) AS n FROM (SELECT COUNT(*) AS n FROM (SELECT DISTINCT ",
        fields, " ", from, ") tuples GROUP BY ", family, ") nesting")
      nested[i] = nested[i] || DBI::dbGetQuery(conn, sql)[["n"]] == 1
    }
  }
  sum(sizes[nested] - 1)
}

retained_vcov_psd = function(V, fix) {
  if (!length(V)) return(V)
  e = eigen(as.matrix(V), symmetric = TRUE)
  negative = min(e[["values"]])
  if (negative >= 0) return(V)
  corrected = tcrossprod(sweep(e[["vectors"]], 2L, sqrt(pmax(e[["values"]], 0)), "*"))
  difference = max(abs(V - corrected))
  if (difference > 1e-8) warning(if (fix) "The covariance matrix was repaired to be positive semidefinite."
    else "The covariance matrix is not positive semidefinite.", call. = FALSE)
  if (fix) {
    old = attributes(V)
    for (name in setdiff(names(old), c("class", "Dim", "Dimnames", "dim"))) {
      attr(corrected, name) = old[[name]]
    }
    V = corrected
  }
  attr(V, "repair") = list(min_eigenvalue = negative, repaired = fix, max_adjustment = difference)
  V
}
