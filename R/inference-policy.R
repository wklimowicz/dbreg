# Correction policy is shared by fitting and post-fit inference consumers.
normalize_ssc = function(ssc) {
  if (is.character(ssc)) {
    if (length(ssc) != 1L || is.na(ssc) || !ssc %in% c("full", "nested", "none")) {
      stop("`ssc` must be full, nested, none, or a correction-policy list.", call. = FALSE)
    }
    ssc = list(K.fixef = switch(ssc,
      full = "full",
      nested = "nonnested",
      none = "none"
    ))
  }
  defaults = list(
    K.adj = TRUE, K.fixef = "nonnested", K.exact = FALSE,
    G.adj = TRUE, G.df = "min", t.df = "min"
  )
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
  numeric_df = is.numeric(defaults[["t.df"]])
  if (numeric_df && (length(defaults[["t.df"]]) != 1L || !is.finite(defaults[["t.df"]]) || defaults[["t.df"]] < 1)) {
    stop("Numeric t.df must be at least one.", call. = FALSE)
  }
  for (name in c("K.fixef", "G.df", if (!numeric_df) "t.df")) {
    if (!is.character(defaults[[name]]) || length(defaults[[name]]) != 1 || is.na(defaults[[name]])) {
      stop("`ssc` policies must be single character strings.", call. = FALSE)
    }
  }
  if (defaults[["K.exact"]] || !defaults[["K.fixef"]] %in% c("full", "none", "nonnested") ||
    !defaults[["G.df"]] %in% c("min", "conventional") || !numeric_df && !defaults[["t.df"]] %in% c("min", "conventional")) {
    stop("Unsupported `ssc` settings.", call. = FALSE)
  }
  defaults
}

fe_intercept_df = function(sizes) {
  if (length(sizes)) sum(sizes) - length(sizes) + 1 else 0
}

linear_parameter_count = function(inputs, model, fe_df = 0, nested = NULL) {
  policy = inputs[["ssc"]][["K.fixef"]]
  if (policy == "none") return(model - fe_df)
  if (inputs[["vcov_type_req"]] == "cluster" && policy == "nonnested" && fe_df > 0) {
    if (is.null(nested)) nested = nested_fe_df(inputs)
    return(model - nested)
  }
  model
}

linear_test_df = function(object) {
  df = attr(object[["vcov"]], "df.t")
  if (is.null(df)) df = object[["df_residual"]]
  df
}
