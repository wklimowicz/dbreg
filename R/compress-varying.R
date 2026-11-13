# Varying compressed fits bind independent numeric, category and payload roles.
# SQL keeps observations; only bounded cells enter the existing moment estimator.
varying_compress_frame = function(inputs, resources) {
  conn = inputs[["conn"]]
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  from = inputs[["from_statement"]]
  fe = inputs[["fe"]]
  trends = unique(unlist(inputs[["fe_slopes"]], use.names = FALSE))
  raw = unique(c(
    inputs[["xvars"]],
    unlist(inputs[["fe_components"]], use.names = FALSE),
    trends
  ))
  design = sql_model_matrix(
    inputs[["fml"]],
    conn,
    sub("^FROM\\s+", "", from),
    expand = "all",
    fe_vars = fe,
    intercept = TRUE,
    factor_refs = inputs[["factor_refs"]],
    reserved_names = sql_design_reserved(inputs)
  )
  register_sql_design(inputs, design)
  x = if (length(design[["aliases"]])) {
    paste0("cx", seq_along(design[["aliases"]]))
  } else {
    character()
  }
  f = paste0("cf", seq_along(fe))
  t = paste0("ct", seq_along(trends))
  payload = paste0("cr", seq_along(raw))
  reserved = unique(c(
    raw,
    fe,
    trends,
    sql_design_reserved(inputs),
    names(inputs[["source_map"]]),
    names(inputs[["name_map"]]),
    x,
    payload,
    "cy",
    "cw",
    "cc"
  ))
  f = safe_identifiers(
    f,
    taken = c(reserved, t),
    fold_case = folds_identifier_case(conn)
  )
  t = safe_identifiers(
    t,
    taken = c(reserved, f),
    fold_case = folds_identifier_case(conn)
  )
  keys = if (!is.null(inputs[["retained_original"]])) inputs[["row_id_safe"]]
  key_aliases = if (length(keys)) paste0("ck", seq_along(keys)) else character()
  fields = c(
    if (length(keys)) paste(q(keys), "AS", key_aliases),
    if (length(x)) {
      paste(sql_cast_double(design[["select_exprs"]], conn), "AS", x)
    },
    paste(q(fe), "AS", f),
    paste(sql_cast_double(q(trends), conn), "AS", t),
    paste(q(raw), "AS", payload),
    paste(sql_cast_double(q(inputs[["yvar"]]), conn), "AS cy"),
    if (is.null(inputs[["weights"]])) {
      "1.0 AS cw"
    } else {
      paste(sql_cast_double(q(inputs[["weights"]]), conn), "AS cw")
    },
    if (length(inputs[["cluster_var"]])) {
      paste(q(inputs[["cluster_var"]]), "AS cc")
    }
  )
  groups = c(x, f, t)
  select = paste("SELECT", paste(fields, collapse = ","), from)
  context = list2env(as.list(inputs), parent = emptyenv())
  context[["retained_compression_keys"]] = key_aliases
  context[["yvar"]] = "cy"
  context[["weights"]] = if (is.null(inputs[["weights"]])) NULL else "cw"
  context[["cluster_var"]] = if (length(inputs[["cluster_var"]])) {
    "cc"
  } else {
    NULL
  }
  if (isTRUE(inputs[["sql_only"]])) {
    context[["from_statement"]] = paste0(
      "FROM (",
      select,
      ") dbreg_coordinate_source"
    )
    return(list(context = context, groups = groups, raw = payload))
  }
  backend = detect_backend(conn)
  if (!(backend[["name"]] == "duckdb" || isTRUE(backend[["is_databricks"]]))) {
    stop(
      "Compressed varying slopes require DuckDB or native Databricks.",
      call. = FALSE
    )
  }
  validate_native_temp_query(conn, backend)
  snapshot = fe_projection_tables(conn, resources)
  table = snapshot(
    select,
    c(groups, payload, "cy", "cw", context[["cluster_var"]], key_aliases)
  )
  key_fields = c(groups, payload)
  # Payloads also carry numeric regressors (e.g. DECIMAL money columns), which
  # are not category keys; only fall back per column when one is rejected.
  types = tryCatch(category_key_metadata(conn, table, key_fields)[["types"]], error = function(e) {
    vapply(key_fields, function(field) tryCatch(category_key_metadata(conn, table, field)[["types"]][[1L]],
      error = function(inner) if (field %in% payload) "float" else stop(inner)), character(1))
  })
  types = setNames(types, key_fields)
  limits = projection_limits()
  stats = DBI::dbGetQuery(
    conn,
    paste0(
      "SELECT COUNT(*) AS n FROM (SELECT DISTINCT ",
      paste(groups, collapse = ","),
      " FROM ",
      table,
      ") cells"
    )
  )
  n = as.numeric(stats[["n"]])
  sizes = as.numeric(DBI::dbGetQuery(
    conn,
    paste(
      "SELECT",
      paste(
        paste0("COUNT(DISTINCT ", f, ") AS n", seq_along(f)),
        collapse = ","
      ),
      "FROM",
      table
    )
  )[1, ])
  dimension = length(x) + sum(sizes * (1 + lengths(inputs[["fe_slopes"]]))) + 1
  cells = n * (length(key_fields) + 14)
  bytes = 64 * cells + 48 * dimension^2
  if (
    !is.finite(n) ||
      n < 1 ||
      n > limits[["summary_rows"]] ||
      cells > limits[["transfer_cells"]] ||
      dimension > limits[["design_entries"]] ||
      !is.finite(bytes) ||
      bytes > limits[["memory_bytes"]]
  ) {
    projection_abort(
      "system_limit",
      "Compressed varying cells/design exceed their projection row, cell, entry or memory budget."
    )
  }
  # Native string keys cannot be allowed to hit ODBC's 65,535-byte ceiling.
  row_bytes = 0
  text = key_fields[types %in% c("character", "date", "integer")]
  if (length(text)) {
    values = sql_cast_text(conn, text)
    if (detect_backend(conn)[["name"]] == "duckdb") {
      values = paste0("encode(", values, ")")
    }
    widths = paste0("OCTET_LENGTH(", values, ")")
    largest_expr = if (length(widths) == 1L) widths[[1L]] else paste0("GREATEST(", paste(widths, collapse = ","), ")")
    cost = DBI::dbGetQuery(
      conn,
      paste0(
        "SELECT MAX(",
        largest_expr,
        ") AS largest,MAX(",
        paste(widths, collapse = "+"),
        ") AS bytes FROM ",
        table
      )
    )
    if (
      !is.finite(cost[["largest"]]) ||
        cost[["largest"]] >= 65535 ||
        !is.finite(cost[["bytes"]]) ||
        bytes + 8 * n * cost[["bytes"]] > limits[["memory_bytes"]]
    ) {
      projection_abort(
        "transfer_limit",
        "Compressed varying keys exceed their byte budget."
      )
    }
    row_bytes = as.numeric(cost[["bytes"]])
    bytes = bytes + 8 * n * row_bytes
  }
  context[["from_statement"]] = paste("FROM", table)
  context[["score_upload_bytes"]] = min(
    limits[["upload_bytes"]],
    limits[["memory_bytes"]] / 8
  )
  context[["composite_integer_keys"]] = c(
    key_fields[types == "integer"],
    if (length(inputs[["cluster_var"]])) {
      cluster_type = category_key_metadata(conn, table, "cc")[["types"]][[1L]]
      if (cluster_type == "integer") "cc"
    }
  )
  list(
    context = context,
    groups = groups,
    types = types,
    sizes = sizes,
    x = x,
    fe = f,
    trends = t,
    raw = payload,
    fe_sources = fe,
    trend_sources = trends,
    raw_sources = raw,
    design = design,
    n = n,
    bytes = bytes,
    row_bytes = row_bytes,
    transfer_cells = cells,
    score_width = length(x) +
      length(f) +
      sum(lengths(inputs[["fe_slopes"]])) +
      1L,
    limits = limits
  )
}

varying_key_factor = function(value, type) {
  keys = category_key_part(value, type)
  if (anyNA(keys)) {
    stop("Compressed FE key receipt is incomplete.", call. = FALSE)
  }
  levels = sort(unique(keys))
  original = value[match(levels, keys)]
  labels = if (type %in% c("numeric", "float")) {
    sprintf("%.17g", original)
  } else {
    as.character(original)
  }
  if (anyDuplicated(labels)) {
    stop(
      "Compressed FE labels do not identify their native keys.",
      call. = FALSE
    )
  }
  factor(labels[match(keys, levels)], levels = labels)
}

# Apply typed FE domains and stable, group-specific trend origins to the matrix
# producer. Recipes save every affine transform for independent raw prediction.
varying_compressed_layout = function(X, formula, data, design, frame) {
  producer_data = data
  for (name in frame[["fe"]]) {
    if (nlevels(data[[name]]) == 1L) producer_data[[name]] = rep(1, nrow(data))
  }
  columns = matrix_column_layout(
    X,
    formula,
    producer_data,
    design[["columns"]],
    slope_terms = seq_along(frame[["x"]])
  )
  internal = columns
  key_receipts = lapply(frame[["fe"]], function(name) {
    values = attr(data[[name]], "dbreg_keys")
    list(
      values = values,
      type = frame[["types"]][[name]],
      labels = levels(data[[name]])
    )
  })
  names(key_receipts) = frame[["fe"]]
  decorate = function(recipe) {
    if (
      recipe[["kind"]] == "numeric" &&
        recipe[["variable"]] %in% names(key_receipts)
    ) {
      receipt = key_receipts[[recipe[["variable"]]]]
      return(list(
        kind = "constant",
        variable = recipe[["variable"]],
        levels = receipt[["values"]],
        key_type = receipt[["type"]]
      ))
    }
    if (
      recipe[["kind"]] == "indicator" &&
        recipe[["variable"]] %in% names(key_receipts)
    ) {
      receipt = key_receipts[[recipe[["variable"]]]]
      recipe[["value"]] = receipt[["values"]][match(
        recipe[["value"]],
        receipt[["labels"]]
      )]
      recipe[["levels"]] = receipt[["values"]]
      recipe[["key_type"]] = receipt[["type"]]
    }
    for (child in c("left", "right")) {
      if (!is.null(recipe[[child]])) {
        recipe[[child]] = decorate(recipe[[child]])
      }
    }
    recipe
  }
  for (j in seq_len(nrow(columns))) {
    recipe = columns[["recipe"]][[j]]
    if (recipe[["kind"]] == "product") {
      parts = list(recipe[["left"]], recipe[["right"]])
      key = which(vapply(
        parts,
        function(p) {
          p[["kind"]] %in%
            c("indicator", "numeric") &&
            p[["variable"]] %in% frame[["fe"]]
        },
        logical(1)
      ))
      trend = which(vapply(
        parts,
        function(p) {
          identical(p[["kind"]], "numeric") &&
            p[["variable"]] %in% frame[["trends"]]
        },
        logical(1)
      ))
      if (length(key) == 1L && length(trend) == 1L) {
        rows = if (parts[[key]][["kind"]] == "numeric") {
          rep(TRUE, nrow(data))
        } else {
          as.character(data[[parts[[key]][["variable"]]]]) ==
            parts[[key]][["value"]]
        }
        value = data[[parts[[trend]][["variable"]]]]
        center = min(value[rows]) / 2 + max(value[rows]) / 2
        X[, j] = as.numeric(rows) * (value - center)
        parts[[trend]][["center"]] = center
        recipe[["left"]] = parts[[1L]]
        recipe[["right"]] = parts[[2L]]
      }
    }
    columns[["recipe"]][[j]] = decorate(recipe)
    if (columns[["role"]][j] == "slope" && columns[["term"]][j] > 0L) {
      internal[["recipe"]][[j]] = list(
        kind = "numeric",
        variable = frame[["x"]][columns[["term"]][j]]
      )
    } else {
      internal[["recipe"]][[j]] = columns[["recipe"]][[j]]
    }
  }
  # Opaque matrix targets never depend on public coefficient labels.
  columns[["alias"]] = internal[["alias"]] = paste0(
    "dbreg_c_basis",
    seq_len(nrow(columns))
  )
  columns[["alias"]][columns[["term"]] == 0L] = internal[["alias"]][
    columns[["term"]] == 0L
  ] = "(Intercept)"
  colnames(X) = columns[["alias"]]
  list(
    X = X,
    columns = columns,
    internal = internal,
    key_receipts = key_receipts
  )
}

execute_varying_compress = function(inputs) {
  conn = inputs[["conn"]]
  resources = fit_resources()
  resources[["conn"]] = conn
  on.exit(release_fit_resources(resources), add = TRUE)
  frame = varying_compress_frame(inputs, resources)
  ctx = frame[["context"]]
  packet = compressed_moment_query(
    conn,
    ctx[["from_statement"]],
    frame[["groups"]],
    "cy",
    ctx[["weights"]],
    TRUE,
    payload = frame[["raw"]]
  )
  if (isTRUE(inputs[["sql_only"]])) {
    return(packet[["sql"]])
  }
  packet[["sql"]] = paste0(
    "SELECT * FROM (",
    sql_limit(
      packet[["sql"]],
      frame[["n"]] + 1,
      detect_backend(conn)[["name"]]
    ),
    ") bounded_cells"
  )
  dat = read_compressed_moments(
    conn,
    packet,
    integer_text = ctx[["composite_integer_keys"]]
  )
  if (nrow(dat) != frame[["n"]]) {
    stop("Compressed coordinate cell support changed.", call. = FALSE)
  }
  token = "raw-v1"
  attr(dat, "dbreg_varying_payload") = list(
    columns = setNames(frame[["raw"]], frame[["raw_sources"]]),
    token = token
  )
  if (isTRUE(inputs[["data_only"]])) {
    output = varying_cell_payload(dat, frame)
    output[["dbreg_cell"]] = NULL
    attr(output, "dbreg_varying_payload") = attr(dat, "dbreg_varying_payload")
    return(output)
  }
  # Cell payloads preserve native values; categorical copies alone become factors.
  matrix_data = dat
  for (name in frame[["fe"]]) {
    value = dat[[name]]
    type = frame[["types"]][[name]]
    matrix_data[[name]] = varying_key_factor(value, type)
    # Recover values by the exact display labels used for this factor.
    labels = if (type %in% c("numeric", "float")) {
      sprintf("%.17g", value)
    } else {
      as.character(value)
    }
    attr(matrix_data[[name]], "dbreg_keys") = value[match(
      levels(matrix_data[[name]]),
      labels
    )]
  }
  nuisance = unlist(
    lapply(seq_along(frame[["fe"]]), function(i) {
      if (length(inputs[["fe_slopes"]][[i]])) {
        paste(
          frame[["fe"]][i],
          frame[["trends"]][match(
            inputs[["fe_slopes"]][[i]],
            frame[["trend_sources"]]
          )],
          sep = ":"
        )
      }
    }),
    use.names = FALSE
  )
  formula = reformulate(
    c(frame[["x"]], frame[["fe"]], nuisance),
    intercept = TRUE
  )
  design = frame[["design"]]
  design[["columns"]][["alias"]] = frame[["x"]]
  producer_data = matrix_data
  one_level = frame[["fe"]][
    vapply(matrix_data[frame[["fe"]]], nlevels, integer(1)) == 1L
  ]
  for (name in one_level) {
    producer_data[[name]] = rep(1, nrow(dat))
  }
  contrasts = lapply(
    matrix_data[setdiff(frame[["fe"]], one_level)],
    function(x) stats::contr.treatment(levels(x))
  )
  if (!length(contrasts)) {
    contrasts = NULL
  }
  X = Matrix::sparse.model.matrix(
    formula,
    producer_data,
    contrasts.arg = contrasts
  )
  layout = varying_compressed_layout(X, formula, matrix_data, design, frame)
  X = layout[["X"]]
  columns = layout[["columns"]]
  S = compressed_cross_moments(
    X,
    dat,
    "cy",
    center_columns = which(
      columns[["role"]] == "slope" & columns[["term"]] > 0L
    )
  )
  slopes = columns[["alias"]][
    columns[["role"]] == "slope" & columns[["term"]] > 0L
  ]
  # Rank priority belongs to the absorbed span. Ordinary regressors already in
  # that span must not survive because of roundoff in a large normal equation.
  residual = compressed_coordinate_projection(
    as.matrix(X[, slopes, drop = FALSE]),
    matrix_data[frame[["fe"]]],
    lapply(inputs[["fe_slopes"]], function(v) {
      matrix_data[frame[["trends"]][match(v, frame[["trend_sources"]])]]
    }),
    dat[["sum_w"]],
    tol = inputs[["fe_tol"]],
    max_iter = max(1000L, inputs[["fe_max_iter"]])
  )
  original_scale = apply(as.matrix(X[, slopes, drop = FALSE]), 2L,
    function(value) max(value)/2 - min(value)/2)
  floor = sum(dat[["sum_w"]]) *
    (128 *
      .Machine[["double.eps"]] *
      original_scale)^2
  absorbed = slopes[colSums(residual^2 * dat[["sum_w"]]) <= floor]
  eligible = setdiff(colnames(X), absorbed)
  order = c(setdiff(eligible, slopes), intersect(eligible, slopes))
  chosen = detect_collinearity(
    S[order, order, drop = FALSE],
    S[order, "cy", drop = FALSE],
    verbose = inputs[["verbose"]]
  )[["keep_names"]]
  kept = eligible[eligible %in% chosen]
  fit = fit_ols(S, "cy", kept)
  collin = setdiff(colnames(X), kept)
  score_X = center_compressed_design(
    X[, kept, drop = FALSE],
    c("(Intercept)" = 0, fit[["centers"]])[kept]
  )
  yhat = as.numeric(score_X %*% fit[["moment_coef"]])
  raw_residual = matrix(dat[["mean_Y"]] - yhat, ncol = 1L)
  within_residual = compressed_coordinate_projection(
    raw_residual,
    matrix_data[frame[["fe"]]],
    lapply(inputs[["fe_slopes"]], function(v) {
      matrix_data[frame[["trends"]][match(v, frame[["trend_sources"]])]]
    }),
    dat[["sum_w"]],
    tol = inputs[["fe_tol"]],
    max_iter = max(1000L, inputs[["fe_max_iter"]])
  )
  drift = max(abs(raw_residual - within_residual))
  rms = sqrt(sum(dat[["sum_w"]] * raw_residual^2) / sum(dat[["sum_w"]]))
  # Drift is judged in units of max(1, residual RMS, residual half-range), the
  # scale the coordinate projection converges in (absolute below unit scale).
  scale = max(1, rms, diff(range(raw_residual)) / 2)
  if (!is.finite(drift) || !is.finite(rms) || drift > inputs[["fe_tol"]] * scale) {
    stop(
      "Compressed coordinate fit failed its independent row certificate.",
      call. = FALSE
    )
  }
  rss = fit[["rss"]]
  df_intercepts = fe_intercept_df(frame[["sizes"]])
  K = length(intersect(kept, slopes)) +
    df_intercepts +
    sum(frame[["sizes"]] * lengths(inputs[["fe_slopes"]]))
  n = sum(dat[["n"]])
  df = n - K
  meat = NULL
  if (inputs[["vcov_type_req"]] == "hc1") {
    outcome = dat[["sum_w2y_sq"]]
    cross = 2 * yhat * dat[["sum_w2y"]]
    fitted = dat[["sum_w2"]] * yhat^2
    rss_w2 = outcome - cross + fitted
    rounding = 64 *
      .Machine[["double.eps"]] *
      pmax(abs(outcome), abs(cross), abs(fitted))
    if (any(!is.finite(rss_w2)) || any(rss_w2 < -rounding)) {
      stop("Invalid compressed HC1 residual moments.", call. = FALSE)
    }
    rss_w2 = pmax(rss_w2, 0)
    meat = crossprod(
      score_X,
      Matrix::Diagonal(x = as.numeric(rss_w2)) %*% score_X
    )
  }
  if (inputs[["vcov_type_req"]] == "cluster") {
    meat = varying_compressed_cluster(
      frame,
      dat,
      score_X,
      yhat,
      fit[["weight_scale"]]
    )
  }
  if (!is.null(meat)) {
    attr(meat, "weight_scale") = fit[["weight_scale"]]
  }
  V = compute_vcov(
    inputs[["vcov_type_req"]],
    "compress",
    fit[["bread"]],
    rss,
    df,
    n,
    n_params = linear_parameter_count(inputs, K, df_intercepts),
    meat = meat,
    sandwich_bread = fit[["sandwich_bread"]],
    ssc = inputs[["ssc"]]
  )
  retain_compressed_linear_rows(inputs, ctx[["from_statement"]], frame[["groups"]],
    dat[["dbreg_cell"]], yhat, "cy", attr(dat, "moment_center"), rss, n,
    ctx[["retained_compression_keys"]], fit, score_X, frame[["sizes"]])
  attr(V, "rss") = rss
  attr(V, "tss") = moment_tss(S, "cy")
  internal = list(
    version = 1L,
    columns = layout[["internal"]],
    fitted_ids = columns[["id"]][match(kept, columns[["alias"]])]
  )
  compact = compact_recipe_dictionaries(internal[["columns"]][["recipe"]])
  internal[["columns"]][["recipe"]] = compact[["recipes"]]
  internal[["dictionaries"]] = compact[["dictionaries"]]
  payload = varying_cell_payload(dat, frame)
  attr(payload, "dbreg_varying_payload") = attr(dat, "dbreg_varying_payload")
  source_map = inputs[["source_map"]]
  bindings = list(
    x = frame[["x"]],
    fe = setNames(
      unname(restore_names(frame[["fe_sources"]], source_map)),
      frame[["fe"]]
    ),
    trends = setNames(
      unname(restore_names(frame[["trend_sources"]], source_map)),
      frame[["trends"]]
    ),
    structural = lapply(
      design[["columns"]][["recipe"]],
      restore_recipe_sources,
      map = source_map
    ),
    internal = internal,
    token = token,
    span = compressed_coordinate_span(S, X, kept, dat[["sum_w"]])
  )
  inputs[["source_map"]] = c(
    source_map,
    setNames(
      unname(restore_names(frame[["fe_sources"]], source_map)),
      frame[["fe"]]
    ),
    setNames(
      unname(restore_names(frame[["trend_sources"]], source_map)),
      frame[["trends"]]
    )
  )
  inputs[["composite_dictionaries"]] = composite_prediction_dictionaries(
    inputs
  )
  list(
    coeftable = gen_coeftable(fit[["coef"]], V, df),
    vcov = V,
    linear_design = list(
      version = 1L,
      columns = columns,
      fitted_ids = internal[["fitted_ids"]]
    ),
    compressed_design = list(bindings = bindings),
    data = payload,
    fml = inputs[["fml"]],
    yvar = inputs[["yvar"]],
    xvars = inputs[["xvars"]],
    fe = inputs[["fe"]],
    weights = inputs[["weights"]],
    collin.var = collin,
    coef_names = intersect(kept, slopes),
    query_string = packet[["sql"]],
    nobs = nrow(dat),
    nobs_orig = n,
    strategy = "compress",
    prediction_coordinates = fit_prediction_state(fit, "cy", V, rss, df, n, meat),
    compression_ratio = nrow(dat) / n,
    compression_ratio_est = inputs[["compression_ratio_est"]],
    compression_diagnostics = list(
      cells = nrow(dat),
      cell_transfer_cells = frame[["transfer_cells"]],
      cell_plan_bytes = frame[["bytes"]],
      cluster_score = attr(meat, "compressed_score_method", exact = TRUE),
      max_abs = drift,
      max_relative = if (
        rms > 128 * .Machine[["double.eps"]] * max(abs(dat[["mean_Y"]]))
      ) {
        drift / rms
      } else {
        0
      }
    ),
    n_fe_levels = setNames(frame[["sizes"]], inputs[["fe"]]),
    df_residual = df
  )
}

# Public cell tables carry an explicit raw-payload receipt. Counts and source
# columns with the same spelling never compete for the same producer.
restore_varying_cells = function(data, inputs) {
  spec = attr(data, "dbreg_varying_payload", exact = TRUE)
  if (is.null(spec)) {
    return(data)
  }
  mapping = character()
  for (source in names(spec[["columns"]])) {
    alias = spec[["columns"]][[source]]
    target = unname(restore_names(source, inputs[["source_map"]]))
    taken = names(data)[names(data) != alias]
    if (target %in% taken) {
      target = utils::tail(
        make.unique(c(taken, paste0(target, "_source"))),
        1L
      )
    }
    names(data)[match(alias, names(data))] = target
    mapping[unname(restore_names(source, inputs[["source_map"]]))] = target
  }
  attr(data, "dbreg_varying_payload") = NULL
  attr(data, "dbreg_cell_sources") = mapping
  attr(data, "dbreg_cell_token") = spec[["token"]]
  data
}

varying_compressed_prediction_matrix = function(object, data) {
  binding = object[["compressed_design"]][["bindings"]]
  prepared_cells = !is.null(attr(data, "dbreg_cell_token", exact = TRUE))
  if (prepared_cells) {
    if (
      !identical(
        attr(data, "dbreg_cell_token", exact = TRUE),
        binding[["token"]]
      )
    ) {
      stop("Unsupported prepared raw-cell format.", call. = FALSE)
    }
    source = attr(data, "dbreg_cell_sources", exact = TRUE)
    if (any(!source %in% names(data))) {
      stop(
        "Prepared cells are missing their raw source columns.",
        call. = FALSE
      )
    }
    data = data[, unname(source), drop = FALSE]
    names(data) = names(source)
  }
  data = prepare_composite_prediction(object, data)
  # Numeric source roles may be exact BIGINT payload strings. Convert only
  # their numeric copy; typed category producers retain their original keys.
  numeric_data = data
  numeric_sources = function(recipe) {
    own = if (recipe[["kind"]] %in% c("numeric", "raw", "power")) {
      recipe[["variable"]]
    } else {
      recipe[["second"]]
    }
    unique(c(
      own,
      unlist(
        lapply(c("left", "right", "source"), function(child) {
          if (is.list(recipe[[child]])) numeric_sources(recipe[[child]])
        }),
        use.names = FALSE
      )
    ))
  }
  sources = unique(unlist(
    lapply(binding[["structural"]], numeric_sources),
    use.names = FALSE
  ))
  for (name in sources) {
    if (
      prepared_cells &&
        is.character(numeric_data[[name]]) &&
        all(
          is.na(numeric_data[[name]]) |
            grepl("^[+-]?[0-9]+$", numeric_data[[name]])
        )
    ) {
      numeric_data[[name]] = as.numeric(numeric_data[[name]])
    }
  }
  invalid = linear_category_invalid(object[["linear_design"]], data)
  prepared = data.frame(row.names = seq_len(nrow(data)))
  for (i in seq_along(binding[["x"]])) {
    prepared[[binding[["x"]][i]]] = linear_recipe_value(
      binding[["structural"]][[i]],
      numeric_data
    )
  }
  for (name in names(binding[["fe"]])) {
    source = binding[["fe"]][[name]]
    if (!source %in% names(data)) {
      stop("Prediction data is missing FE source: ", source, call. = FALSE)
    }
    prepared[[name]] = data[[source]]
  }
  for (name in names(binding[["trends"]])) {
    value = data[[binding[["trends"]][[name]]]]
    if (
      prepared_cells &&
        is.character(value) &&
        all(is.na(value) | grepl("^[+-]?[0-9]+$", value))
    ) {
      value = as.numeric(value)
    }
    if (!is.numeric(value) && !inherits(value, "integer64")) {
      stop("Prediction trends must be numeric.", call. = FALSE)
    }
    prepared[[name]] = as.numeric(value)
  }
  adapted = object
  adapted[["linear_design"]] = binding[["internal"]]
  matrix = linear_prediction_matrix(adapted, prepared)
  span = binding[["span"]]
  current = center_compressed_design(matrix, span[["centers"]][span[["kept"]]])
  norm = 1 +
    sqrt(as.numeric(Matrix::rowSums(
      (current %*% Matrix::Diagonal(x = 1 / span[["scales"]][span[["kept"]]]))^2
    )))
  for (i in seq_along(span[["dropped"]])) {
    index = span[["dropped"]][i]
    value = linear_recipe_value(
      binding[["internal"]][["columns"]][["recipe"]][[index]],
      prepared,
      binding[["internal"]][["dictionaries"]]
    )
    error = (value -
      span[["centers"]][index] -
      as.numeric(current %*% span[["coefficients"]][, i])) /
      span[["scales"]][index]
    invalid = invalid |
      !is.finite(error) |
      abs(error) > span[["tolerance"]] * norm
  }
  for (name in c(binding[["x"]], names(binding[["trends"]]))) {
    invalid = invalid | !is.finite(prepared[[name]])
  }
  attr(matrix, "dbreg_invalid_rows") = attr(
    matrix,
    "dbreg_invalid_rows",
    exact = TRUE
  ) |
    invalid
  matrix
}

# Save dependencies in centered score coordinates. This makes prediction
# estimability independent of a large raw origin or an arbitrary FE reference.
compressed_coordinate_span = function(S, X, kept, weights) {
  dropped = setdiff(colnames(X), kept)
  centers = attr(S, "centers")[colnames(X)]
  centers[is.na(centers)] = 0
  scales = sqrt(pmax(diag(S)[seq_len(ncol(X))] / sum(weights), 0))
  scales[scales == 0] = 1
  names(centers) = names(scales) = colnames(X)
  coefficients = if (length(dropped)) {
    solve_with_fallback(
      S[kept, kept, drop = FALSE],
      S[kept, dropped, drop = FALSE]
    )[["betahat"]]
  } else {
    matrix(0, length(kept), 0)
  }
  current = center_compressed_design(X[, kept, drop = FALSE], centers[kept])
  norm = 1 +
    sqrt(as.numeric(Matrix::rowSums(
      (current %*% Matrix::Diagonal(x = 1 / scales[kept]))^2
    )))
  for (i in seq_along(dropped)) {
    error = (as.numeric(X[, dropped[i]]) -
      centers[[dropped[i]]] -
      as.numeric(current %*% coefficients[, i])) /
      scales[[dropped[i]]]
    if (any(!is.finite(error)) || any(abs(error) > 1e-8 * norm)) {
      stop(
        "Compressed rank dependencies failed their row certificate.",
        call. = FALSE
      )
    }
  }
  list(
    kept = match(kept, colnames(X)),
    dropped = match(dropped, colnames(X)),
    coefficients = coefficients,
    centers = unname(centers),
    scales = unname(scales),
    tolerance = 1e-8
  )
}


varying_cell_payload = function(data, frame) {
  payload = data[c(
    frame[["raw"]],
    setdiff(names(data), c(frame[["groups"]], frame[["raw"]]))
  )]
  attr(payload, "moment_center") = attr(data, "moment_center", exact = TRUE)
  attr(payload, "weight_scale") = attr(data, "weight_scale", exact = TRUE)
  restore_compressed_data(payload)
}


varying_compressed_cluster = function(frame, data, X, yhat, weight_scale) {
  inputs = frame[["context"]]
  conn = inputs[["conn"]]
  limits = frame[["limits"]]
  groups = frame[["groups"]]
  from = inputs[["from_statement"]]
  center = attr(data, "moment_center", exact = TRUE)
  rank = ncol(X)
  cells_left = limits[["transfer_cells"]] - frame[["transfer_cells"]]
  if (
    use_cluster_sql(nrow(data), rank, inputs) &&
      1 + rank * (rank + 1) / 2 <= cells_left &&
      nrow(data) < limits[["summary_rows"]]
  ) {
    result = tryCatch(
      compute_meat_cluster_compress_sql(
        inputs,
        from,
        groups,
        data[["dbreg_cell"]],
        X,
        list(yhat),
        y_center = center,
        weight_scale = weight_scale
      )[[1L]],
      dbreg_projection_capacity = function(e) {
        if (!identical(e[["reason"]], "upload_limit")) {
          stop(e)
        }
        NULL
      }
    )
    if (!is.null(result)) {
      attr(result, "compressed_score_method") = "sql"
      return(result)
    }
  }
  keys = c("cc", groups)
  support = DBI::dbGetQuery(
    conn,
    paste0(
      "SELECT COUNT(*) AS n, COUNT(DISTINCT cc) AS g FROM (SELECT DISTINCT ",
      paste(keys, collapse = ","),
      " ",
      from,
      ") cluster_cells"
    )
  )
  n = as.numeric(support[["n"]])
  clusters = as.numeric(support[["g"]])
  width = length(keys) + 2L
  bytes = n *
    (64 * (width + frame[["score_width"]]) + 8 * frame[["row_bytes"]]) +
    16 * clusters * rank
  cc = sql_cast_text(conn, "cc")
  if (detect_backend(conn)[["name"]] == "duckdb") {
    cc = paste0("encode(", cc, ")")
  }
  longest = as.numeric(DBI::dbGetQuery(
    conn,
    paste0("SELECT MAX(OCTET_LENGTH(", cc, ")) AS n ", from)
  )[["n"]])
  if (
    !is.finite(n) ||
      n < 1 ||
      n + nrow(data) > limits[["summary_rows"]] ||
      n * width > cells_left ||
      !is.finite(longest) ||
      longest >= 65535 ||
      !is.finite(bytes) ||
      bytes + 8 * n * longest + frame[["bytes"]] > limits[["memory_bytes"]]
  ) {
    projection_abort(
      "transfer_limit",
      "Compressed varying cluster cells exceed their transfer or byte budget."
    )
  }
  inputs[["cluster_cell_count"]] = n
  result = compute_meat_cluster_compress(
    inputs,
    from,
    groups,
    data,
    X,
    yhat,
    y_center = center,
    weight_scale = weight_scale
  )
  attr(result, "compressed_score_method") = "cluster_cells"
  result
}

# The compressed varying-slope route projects cells by plain alternating
# projections in R; weakly connected FE graphs can exceed its sweep cap.
compressed_coordinate_projection = function(...) {
  tryCatch(demean_prediction_coordinates(...), error = function(e) {
    if (!grepl("did not converge", conditionMessage(e))) stop(e)
    stop("Compressed varying-slope projection did not converge (weakly connected fixed effects). ",
      "Consider strategy = 'demean', which strategy = 'auto' uses for varying slopes.", call. = FALSE)
  })
}
