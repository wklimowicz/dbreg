retained_glm_covariance = function(object, vcov, cluster, ssc) {
  if (identical(object[["inference_result"]][["reason"]], "not_requested")) {
    stop("The model used compute_inference=FALSE. Refit with inference enabled to compute covariance.", call. = FALSE)
  }
  state = object[["retained_inference"]]
  if (is.null(state)) stop("Post-fit covariance requires retain = 'inference' or 'all'.", call. = FALSE)
  if (!isTRUE(object[["finite_mle"]])) {
    stop("Post-fit covariance requires a finite GLM estimate.", call. = FALSE)
  }
  prep = prepare_retained_covariance(object, vcov, cluster, ssc,
    state[["source"]], state[["source_keys"]], state[["cluster_classes"]],
    state[["default_cluster"]], state[["families"]], state[["sizes"]],
    "__dbglm_estimated", allow_twoway = FALSE, native_map = state[["cluster_fields"]])
  conn = prep[["conn"]]
  V0 = prep[["V0"]]
  groups = prep[["groups"]]
  kind = prep[["kind"]]
  ssc = prep[["ssc"]]
  q = prep[["q"]]
  from = prep[["from"]]
  columns = if (length(groups)) groups[[1L]]
  cols = if (length(columns)) paste0("b.", q(prep[["native"]](columns))) else character()
  n = nobs(object)
  p = length(state[["design"]])
  K = object[["n_params_model"]]
  if (ssc[["K.fixef"]] == "none") K = p + state[["slope_df"]]
  if (length(state[["sizes"]]) && length(cols) && ssc[["K.fixef"]] == "nonnested") {
    K = max(p + state[["slope_df"]] + 1, K - prep[["nested_df"]])
  }
  if (ssc[["K.adj"]] && n <= K) stop("Insufficient residual degrees of freedom.", call. = FALSE)
  H = state[["H"]]
  V = if (p) dbglm_solve(H)[["inverse"]] else H
  weight_scale = state[["weight_scale"]]
  family = object[["family"]][["family"]]
  G = NULL
  if (kind == "iid") {
    V = V * (if (family == "gaussian") object[["dispersion"]] else 1) / weight_scale
    if (ssc[["K.adj"]]) V = V * (n - 1) / (n - K)
  } else {
    mu = "s.__dbglm_saved_mu"
    variance = switch(family, gaussian = "1.0", poisson = mu,
      binomial = paste0(mu, "*(1.0-", mu, ")"))
    scores = if (p) paste0("s.", state[["design"]], "*(s.__dbglm_saved_w/(",
      dbglm_number(weight_scale), "))*s.__dbglm_saved_d*s.__dbglm_saved_res/(",
      variance, ")") else character()
    meat = retained_score_cross(conn, from, scores, cols)
    if (p) V = V %*% meat[["matrix"]] %*% V
    if (length(cols)) {
      G = meat[["groups"]]
      if (G < 2) stop("Clustered covariance requires at least two clusters.", call. = FALSE)
      if (ssc[["G.adj"]]) V = V * G / (G - 1)
    }
    if (ssc[["K.adj"]]) V = V * (if (kind == "hc1") n else n - 1) / (n - K)
  }
  if (!is.null(state[["transform"]]) && p) {
    B = state[["transform"]]
    V = B %*% V %*% t(B)
    V = (V + t(V)) / 2
    dimnames(V) = dimnames(V0)
  }
  attr(V, "type") = kind
  attr(V, "n_params") = attr(V, "df.K") = K
  attr(V, "df.t") = if (is.numeric(ssc[["t.df"]])) ssc[["t.df"]] else
    if (length(G) && ssc[["t.df"]] == "min") G - 1 else n - K
  attr(V, "ssc") = ssc
  attr(V, "n_clusters") = G
  V
}
