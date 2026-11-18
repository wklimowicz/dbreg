# Build scaled within-group coordinates without mixing key/RHS/trend roles.
# Local reorthogonalisation makes redundant trends exact zero axes; they stay
# in the specification and count convention while solvers use the actual span.
fe_varying_basis = function(
  conn,
  table,
  fe,
  slopes,
  trends,
  scales,
  fields,
  create,
  drop
) {
  basis = character()
  ops = lapply(fe, function(key) list(f = key, v = "1.0"))
  for (i in seq_along(fe)) {
    previous = character()
    partition = paste0(" OVER (PARTITION BY ", fe[i], ")")
    for (trend in slopes[[i]]) {
      column = paste0("axis", length(basis) + 1L)
      raw = trends[[trend]]
      # Subtract a per-group midpoint BEFORE scaling. Dividing a trend near
      # 1e12 by its global maximum first destroys representable within-group
      # differences, even when a certificate on that altered basis looks good.
      anchor = paste0("MIN(", raw, ")", partition, "/2+MAX(", raw, ")", partition, "/2")
      delta = paste0("(", raw, "-(", anchor, "))/(", sprintf("%.17e", scales[[trend]]), ")")
      next_table = create(paste("SELECT *,", delta, "AS", column, "FROM", table))
      if (length(basis)) drop(table)
      table = next_table
      fields = unique(c(fields, column))
      centered = paste0(column, "-SUM(w*", column, ")", partition, "/SUM(w)", partition)
      next_table = create(paste("SELECT", paste(c(setdiff(fields, column), paste(centered, "AS", column)), collapse = ","), "FROM", table))
      drop(table); table = next_table
      reference = paste0(column, "_scale")
      next_table = create(paste(
        "SELECT *,",
        paste0("MAX(ABS(", column, "))", partition),
        "AS",
        reference,
        "FROM",
        table
      ))
      drop(table)
      table = next_table
      fields = c(fields, reference)
      for (pass in 1:2) {
        for (v in previous) {
          value = paste0(
            column,
            "-",
            v,
            "*COALESCE(SUM(w*",
            v,
            "*",
            column,
            ")",
            partition,
            "/NULLIF(SUM(w*",
            v,
            "*",
            v,
            ")",
            partition,
            ",0),0)"
          )
          next_table = create(paste(
            "SELECT",
            paste(
              c(setdiff(fields, column), paste(value, "AS", column)),
              collapse = ","
            ),
            "FROM",
            table
          ))
          drop(table)
          table = next_table
        }
      }
      # Only axes at local rounding precision relative to their raw trend are
      # null. Avoid a fixed absolute cutoff that discards small real trends.
      local_scale = paste0("MAX(ABS(", column, "))", partition)
      raw_scale = reference
      norm = paste0(
        "SQRT(SUM(w*",
        column,
        "*",
        column,
        ")",
        partition,
        "/SUM(w)",
        partition,
        ")"
      )
      value = paste0(
        "CASE WHEN ",
        local_scale,
        "<=32*2.220446049250313e-16*",
        raw_scale,
        " THEN 0.0 ELSE ",
        column,
        "/NULLIF(",
        norm,
        ",0) END AS ",
        column
      )
      next_table = create(paste(
        "SELECT",
        paste(c(setdiff(fields, column), value), collapse = ","),
        "FROM",
        table
      ))
      drop(table)
      table = next_table
      previous = c(previous, column)
      basis = c(basis, column)
      ops[[length(ops) + 1L]] = list(f = fe[i], v = column)
    }
  }
  list(table = table, basis = basis, ops = ops)
}

# Bounded scalar rows for arbitrary group-coordinate products. SQL LIMIT
# probes precede transfer; zero diagonal axes never get promoted by scaling.
fe_coordinate_fetch = function(conn, sql, width, budget, limits) {
  room = min(
    limits[["summary_rows"]] - budget[["rows"]],
    floor((limits[["transfer_cells"]] - budget[["cells"]]) / width) - 1,
    floor(
      (limits[["memory_bytes"]] - budget[["summary_bytes"]]) / (64 * width)
    ) -
      1
  )
  if (room < 1) {
    projection_abort(
      "transfer_limit",
      "Coordinate summaries exceed the transfer budget."
    )
  }
  value = DBI::dbGetQuery(
    conn,
    sql_limit(sql, room + 1, detect_backend(conn)[["name"]])
  )
  if (nrow(value) > room || ncol(value) != width) {
    projection_abort(
      "transfer_limit",
      "Coordinate summary overflow probe exceeds its budget."
    )
  }
  budget[["summary_bytes"]] = budget[["summary_bytes"]] +
    64 * (nrow(value) + 1) * width
  budget[["rows"]] = budget[["rows"]] + nrow(value)
  budget[["cells"]] = budget[["cells"]] + (nrow(value) + 1) * width
  if (budget[["cells"]] > limits[["transfer_cells"]]) {
    projection_abort(
      "transfer_limit",
      "Coordinate summaries exceed the transfer budget."
    )
  }
  value
}

fe_coordinate_system = function(
  conn,
  table,
  ops,
  columns,
  budget,
  limits,
  max_iter = 8L
) {
  sum_sql = function(x) sql_sum(conn, x)
  info = lapply(ops, function(op) {
    value = fe_coordinate_fetch(
      conn,
      paste0(
        "SELECT ",
        op[["f"]],
        " AS id,",
        sum_sql(paste0("w*", op[["v"]], "*", op[["v"]])),
        " AS h,",
        paste(
          paste(
            sum_sql(paste0("w*", op[["v"]], "*", columns)),
            "AS",
            paste0("b", seq_along(columns))
          ),
          collapse = ","
        ),
        " FROM ",
        table,
        " GROUP BY ",
        op[["f"]],
        " ORDER BY ",
        op[["f"]]
      ),
      2 + length(columns),
      budget,
      limits
    )
    if (any(!is.finite(as.matrix(value))) || any(value[["h"]] < 0)) {
      projection_abort("singular", "Invalid coordinate information.")
    }
    value = value[value[["h"]] > 0, , drop = FALSE]
    list(
      id = as.numeric(value[["id"]]),
      h = as.numeric(value[["h"]]),
      b = as.matrix(value[paste0("b", seq_along(columns))])
    )
  })
  sizes = lengths(lapply(info, `[[`, "id"))
  active = which(sizes > 0)
  info = info[active]
  ops = ops[active]
  sizes = sizes[active]
  # Signed general coordinate blocks cannot use the additive component PCG.
  if (
    !length(sizes) ||
      sum(sizes) - max(sizes) > limits[["dense_levels"]] ||
      sum(sizes) * length(columns) > limits[["design_entries"]] ||
      !fe_dense_feasible(sizes, length(columns), limits)
  ) {
    projection_abort(
      "system_limit",
      "General coordinate system exceeds its dense Schur budget."
    )
  }
  pairs = list()
  if (length(ops) > 1L) {
    for (i in seq_along(ops)) {
      if (i > 1L) {
        for (j in seq_len(i - 1L)) {
          oi = ops[[i]]
          oj = ops[[j]]
          value = fe_coordinate_fetch(
            conn,
            paste0(
              "SELECT ",
              oi[["f"]],
              " AS i,",
              oj[["f"]],
              " AS j,",
              sum_sql(paste0("w*", oi[["v"]], "*", oj[["v"]])),
              " AS h FROM ",
              table,
              " GROUP BY ",
              paste(unique(c(oi[["f"]], oj[["f"]])), collapse = ",")
            ),
            3L,
            budget,
            limits
          )
          if (any(!is.finite(as.matrix(value)))) {
            projection_abort("singular", "Nonfinite coordinate cross products.")
          }
          ir = match(value[["i"]], info[[i]][["id"]])
          jr = match(value[["j"]], info[[j]][["id"]])
          invalid = is.na(ir) | is.na(jr)
          if (any(invalid & value[["h"]] != 0)) {
            projection_abort(
              "singular",
              "Zero coordinate has nonzero cross information."
            )
          }
          keep = !invalid & value[["h"]] != 0
          pairs[[length(pairs) + 1L]] = list(
            i = i,
            j = j,
            ir = ir[keep],
            jr = jr[keep],
            h = value[["h"]][keep]
          )
        }
      }
    }
  }
  if (
    projection_memory(sizes, length(columns), limits = limits) +
      budget[["summary_bytes"]] >
      limits[["memory_bytes"]]
  ) {
    projection_abort(
      "system_limit",
      "Coordinate summaries and factor buffers exceed the memory budget."
    )
  }
  sys = fe_system(info, pairs, limits)
  if (is.null(sys)) {
    projection_abort(
      "system_limit",
      "General coordinate cross blocks exceed their budget."
    )
  }
  cache = new.env(parent = emptyenv())
  cache[["projection_limits"]] = limits
  solver = fe_solver(sys, cache)
  if (is.null(solver)) {
    projection_abort("singular", "General coordinate solver is unavailable.")
  }
  B = solver[["project"]](do.call(rbind, lapply(info, `[[`, "b")))
  coefficients = solver[["solve"]](B)
  # Refine to 1e-12 relative to the right-hand side (absolute below unit scale).
  target = 1e-12 * max(1, abs(B / sqrt(sys[["diag"]])))
  for (cycle in seq_len(min(max_iter, 8L))) {
    residual = B - fe_multiply(sys, coefficients)
    err = max(abs(residual / sqrt(sys[["diag"]])))
    if (is.finite(err) && err <= target) {
      break
    }
    if (cycle == min(max_iter, 8L)) {
      projection_abort(
        "nonconverged",
        "Coordinate refinement did not converge."
      )
    }
    coefficients = coefficients + solver[["solve"]](residual)
  }
  list(
    info = info,
    ops = ops,
    sys = sys,
    coefficients = coefficients,
    iterations = cycle,
    solver = solver[["method"]]
  )
}

# Certify every original operator, including zero/redundant axes. The norm
# reduces to the usual group-mean certificate for intercept operators.
fe_coordinate_certificate = function(
  conn,
  table,
  ops,
  columns,
  nobs,
  original_scale
) {
  sum_sql = function(x) sql_sum(conn, x)
  branches = vapply(
    ops,
    function(op) {
      v = op[["v"]]
      denom = paste0(
        "SQRT(",
        sum_sql("w"),
        "*",
        sum_sql(paste0("w*", v, "*", v)),
        ")"
      )
      means = paste0(
        "COALESCE(ABS(",
        sum_sql(paste0("w*", v, "*", columns)),
        ")/NULLIF(",
        denom,
        ",0),0) AS r",
        seq_along(columns)
      )
      checked = unique(c("w", v, columns))
      invalid = paste(
        c(
          "w<=0",
          paste0(checked, " IS NULL"),
          paste0("NOT ", sql_is_finite(conn, checked))
        ),
        collapse = " OR "
      )
      paste0(
        "SELECT ",
        paste(
          c(
            means,
            "COUNT(*) AS n",
            paste0("SUM(CASE WHEN ", invalid, " THEN 1 ELSE 0 END) AS bad")
          ),
          collapse = ","
        ),
        " FROM ",
        table,
        " GROUP BY ",
        op[["f"]]
      )
    },
    character(1)
  )
  row = DBI::dbGetQuery(
    conn,
    paste0(
      "SELECT ",
      paste(
        c(
          paste0("MAX(r", seq_along(columns), ") AS r", seq_along(columns)),
          "MAX(bad) AS bad",
          "SUM(n) AS n"
        ),
        collapse = ","
      ),
      " FROM (",
      paste(branches, collapse = " UNION ALL "),
      ") certified"
    )
  )
  if (
    !isTRUE(as.numeric(row[["bad"]]) == 0) ||
      !isTRUE(as.numeric(row[["n"]]) == nobs * length(ops))
  ) {
    return(list(max_abs = Inf, max_level = Inf, max_relative = Inf))
  }
  residual = vapply(
    row[paste0("r", seq_along(columns))],
    as.numeric,
    numeric(1)
  )
  units = sprintf("%.17e", pmax(original_scale, .Machine[["double.xmin"]]))
  terms = paste0(
    "(",
    units,
    ")*SQRT(",
    sum_sql(paste0("w*POWER(", columns, "/(", units, "),2)")),
    "/",
    sum_sql("w"),
    ") AS r",
    seq_along(columns)
  )
  rms = vapply(
    DBI::dbGetQuery(
      conn,
      paste("SELECT", paste(terms, collapse = ","), "FROM", table)
    ),
    as.numeric,
    numeric(1)
  )
  if (any(!is.finite(c(residual, rms)))) {
    return(list(max_abs = Inf, max_level = Inf, max_relative = Inf))
  }
  identified = rms > 128 * .Machine[["double.eps"]] * original_scale
  list(
    max_abs = max(residual),
    max_level = max(residual / pmax(1, original_scale)),
    max_relative = if (any(identified)) {
      max(residual[identified] / rms[identified])
    } else {
      0
    }
  )
}

# OLS varying slopes share the native snapshot and coefficient kernel. Signed
# bounded operator products are preferred; a typed decline tries SQL CG.
dbreg_varying_projections = function(inputs, x_sql, x_names, max_iter, tol) {
  policy = getOption("dbreg.ap_solver", "auto")
  if (
    !is.character(policy) ||
      length(policy) != 1L ||
      is.na(policy) ||
      !policy %in% c("auto", "moments", "iterative")
  ) {
    stop("dbreg.ap_solver must be auto, moments or iterative.", call. = FALSE)
  }
  prepared = prepare_fe_projection(
    inputs,
    x_sql,
    x_names,
    materialize = TRUE,
    coefficient_cg = TRUE
  )
  if (is.null(prepared)) {
    stop(
      "Varying-slope projection requires DuckDB or native Databricks.",
      call. = FALSE
    )
  }
  env2env(prepared, environment(), names(prepared))
  success = FALSE
  on.exit(if (!success) release_fit_resources(resources))
  manager = fe_cg_manager(conn, resources)
  create = manager[["create"]]
  drop = manager[["drop"]]
  operator_sizes = c(sizes, rep(sizes, lengths(inputs[["fe_slopes"]])))
  if (
    sum(operator_sizes) * length(columns) > limits[["design_entries"]] ||
      128 * sum(operator_sizes) * length(columns) > limits[["memory_bytes"]]
  ) {
    projection_abort(
      "system_limit",
      "Varying coordinates exceed their state budget."
    )
  }
  identifiers = c(id_columns, "raw_w", fe_columns, unname(trend_columns))
  fields = c(identifiers, "w", columns)
  basis = fe_varying_basis(
    conn,
    centered,
    fe_columns,
    inputs[["fe_slopes"]],
    trend_columns,
    trend_scale,
    fields,
    create,
    drop
  )
  table = basis[["table"]]
  ops = basis[["ops"]]
  identifiers = c(identifiers, basis[["basis"]])
  budget = list2env(
    list(rows = 0, cells = 0, upload_bytes = 0, summary_bytes = 0),
    parent = emptyenv()
  )
  projected = NULL
  iterations = 0L
  restarts = 0L
  retries = 0L
  solver = NULL
  if (policy != "iterative") {
    owned = resources[["tables"]]
    attempt = projection_attempt("coordinate_moments", {
      system = fe_coordinate_system(
        conn,
        table,
        ops,
        columns,
        budget,
        limits,
        max_iter
      )
      lookup = vapply(
        seq_along(system[["ops"]]),
        function(i) {
          n = max(system[["info"]][[i]][["id"]])
          coeff = matrix(0, n, length(columns))
          positions = seq.int(
            system[["sys"]][["starts"]][i] + 1L,
            system[["sys"]][["starts"]][i + 1L]
          )
          coeff[system[["info"]][[i]][["id"]], ] = system[["coefficients"]][
            positions,
            ,
            drop = FALSE
          ]
          fe_lookup_source(conn, coeff, budget, limits, create)
        },
        character(1)
      )
      active_ops = system[["ops"]]
      joins = vapply(
        seq_along(lookup),
        function(i) {
          paste0(
            " LEFT JOIN ",
            lookup[i],
            " p",
            i,
            " ON s.",
            active_ops[[i]][["f"]],
            "=p",
            i,
            ".id"
          )
        },
        character(1)
      )
      subtract = vapply(
        seq_along(columns),
        function(j) {
          terms = vapply(
            seq_along(lookup),
            function(i) {
              paste0(
                "(",
                if (active_ops[[i]][["v"]] == "1.0") {
                  "1.0"
                } else {
                  paste0("s.", active_ops[[i]][["v"]])
                },
                ")*COALESCE(p",
                i,
                ".c",
                j,
                ",0)"
              )
            },
            character(1)
          )
          paste0(
            "(s.",
            columns[j],
            "-(",
            paste(terms, collapse = "+"),
            "))*(",
            sprintf("%.17e", scale[j]),
            ") AS ",
            columns[j]
          )
        },
        character(1)
      )
      projected = create(paste0(
        "SELECT ",
        paste(c(paste0("s.", c(identifiers, "w")), subtract), collapse = ","),
        " FROM ",
        table,
        " s",
        paste(joins, collapse = "")
      ))
      certificate = fe_coordinate_certificate(
        conn,
        projected,
        ops,
        columns,
        nobs,
        original_scale
      )
      if (
        fe_certificate_shortfall(certificate, tol) > 1
      ) {
        projection_abort(
          "certificate_failed",
          "Coordinate moment projection failed its independent row certificate.",
          certificate
        )
      }
      list(
        table = projected,
        iterations = system[["iterations"]],
        solver = system[["solver"]],
        certificate = certificate
      )
    })
    accepted = record_projection_attempt(inputs, attempt)
    if (!is.null(accepted)) {
      projected = accepted[["table"]]
      iterations = accepted[["iterations"]]
      solver = accepted[["solver"]]
      certificate = accepted[["certificate"]]
    } else {
      for (name in rev(setdiff(resources[["tables"]], owned))) {
        drop(as.character(DBI::dbQuoteIdentifier(conn, name)))
      }
      projected = NULL
    }
  }
  if (is.null(projected)) {
    budget = list2env(
      list(rows = 0, cells = 0, upload_bytes = 0, summary_bytes = 0),
      parent = emptyenv()
    )
    control = fe_cg_control()
    if (isTRUE(inputs[["fe_max_iter_explicit"]])) control[["fe_max_iter"]] = max_iter
    weight_groups = vapply(
      fe_columns,
      function(key) {
        paste0("SELECT SUM(w) AS sw FROM ", table, " GROUP BY ", key)
      },
      character(1)
    )
    totals = fe_cg_fetch(
      conn,
      paste0(
        "SELECT MIN(sw) AS smallest,SUM(sw)/",
        length(fe_columns),
        " AS total FROM (",
        paste(weight_groups, collapse = " UNION ALL "),
        ") groups"
      ),
      2L,
      budget,
      limits
    )
    share = as.numeric(totals[["smallest"]]) / as.numeric(totals[["total"]])
    if (!is.finite(share) || share <= 0) {
      projection_abort("singular", "Invalid varying-coordinate group weights.")
    }
    target = max(4e-15, min(tol, 1e-12) * min(1, sqrt(share)))
    for (retry in 1:2) {
      cg = fe_sql_cg(
        conn,
        table,
        columns,
        identifiers,
        ops,
        rep(1, length(columns)),
        target,
        control,
        create,
        drop,
        budget,
        limits,
        allow_zero = TRUE,
        replace = manager[["replace"]],
        output_scales = scale
      )
      iterations = iterations + cg[["cycles"]]
      restarts = restarts + cg[["restarts"]]
      retries = retry - 1L
      projected = cg[["table"]]
      certificate = fe_coordinate_certificate(
        conn,
        projected,
        ops,
        columns,
        nobs,
        original_scale
      )
      shortfall = fe_certificate_shortfall(certificate, tol)
      if (is.finite(shortfall) && shortfall <= 1) {
        break
      }
      drop(projected)
      target = target / (4 * shortfall)
      if (retry == 2L || !is.finite(target) || target < 5e-16) {
        projection_abort(
          "certificate_failed",
          "Coordinate CG failed its independent row certificate.",
          certificate
        )
      }
    }
    solver = "coefficient_cg"
    record_projection_attempt(
      inputs,
      list(
        method = "coordinate_cg",
        status = "completed",
        reason = NULL,
        details = list(iterations = iterations),
        result = list(table = projected)
      )
    )
  }
  result = create(paste(
    "SELECT",
    paste(
      c(
        paste(id_columns, "AS", q(ids)),
        "raw_w AS __w",
        paste(columns, "AS", q(paste0(vars, "_tilde")))
      ),
      collapse = ","
    ),
    "FROM",
    projected
  ))
  # Only the final physical row table survives into the estimation frame.
  final_name = gsub('["`]', '', result)
  for (name in rev(setdiff(resources[["tables"]], final_name))) {
    drop(as.character(DBI::dbQuoteIdentifier(conn, name)))
  }
  success = TRUE
  list(
    table = result,
    resources = resources,
    original_scale = setNames(original_scale, vars),
    iterations = iterations,
    accelerated_steps = 0L,
    max_abs = certificate[["max_abs"]],
    max_relative = certificate[["max_relative"]],
    method = if (solver == "coefficient_cg") {
      "coordinate_cg"
    } else {
      "coordinate_moments"
    },
    solver = solver,
    solver_iterations = iterations,
    transfer_cells = budget[["cells"]],
    summary_rows = budget[["rows"]],
    upload_bytes = budget[["upload_bytes"]],
    estimated_bytes = if (solver == "coefficient_cg") 128 * sum(operator_sizes) * length(columns) else
      projection_memory(system[["sys"]][["sizes"]], length(columns), limits = limits) + budget[["summary_bytes"]],
    moment_method = "rows",
    observation_tables = 1L,
    transformed_observation_tables = 1L,
    peak_tables = resources[["peak_tables"]],
    cg_restarts = restarts, cg_retries = retries
  )
}
