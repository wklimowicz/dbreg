# Bounded direct nuisance solves for moderate FE dimensions. The database sends
# only FE cross-products and FE-by-design moments; no observation vectors enter R.
# The largest FE block is diagonal and is absorbed analytically, so only the
# remaining levels form a dense system (Schur complement).
dbglm_fe_dense_feasible = function(op_sizes) {
  op_sizes = as.numeric(op_sizes)
  if (!length(op_sizes) || any(!is.finite(op_sizes))) return(FALSE)
  rest = sum(op_sizes) - max(op_sizes)
  rest <= getOption("dbreg.glm_fe_dense_max", 4096) &&
    max(op_sizes) * rest <= getOption("dbreg.glm_fe_cross_max", 5e7)
}

# info: per-op list(id, h); pairs: list of list(i, j, ir, jr, h) with local
# level indices into info[[i]]$id and info[[j]]$id.
dbglm_fe_system = function(info, pairs) {
  sizes = vapply(info, function(x) length(x$id), integer(1))
  if (!dbglm_fe_dense_feasible(sizes)) return(NULL)
  starts = cumsum(c(0L, sizes))
  n = starts[length(info) + 1L]
  a = which.max(sizes)
  big = seq.int(starts[a] + 1L, starts[a + 1L])
  rest = setdiff(seq_len(n), big)
  pos = integer(n)
  pos[rest] = seq_along(rest)
  diagonal = as.numeric(unlist(lapply(info, `[[`, "h"), use.names = FALSE))
  Cm = matrix(0, length(big), length(rest))
  E = matrix(0, length(rest), length(rest))
  E[cbind(seq_along(rest), seq_along(rest))] = diagonal[rest]
  for (p in pairs) {
    gi = starts[p$i] + p$ir
    gj = starts[p$j] + p$jr
    if (p$i == a) {
      Cm[cbind(gi - starts[a], pos[gj])] = p$h
    } else if (p$j == a) {
      Cm[cbind(gj - starts[a], pos[gi])] = p$h
    } else {
      E[cbind(pos[gi], pos[gj])] = p$h
      E[cbind(pos[gj], pos[gi])] = p$h
    }
  }
  list(n = n, starts = starts, big = big, rest = rest, Dg = diagonal[big], Cm = Cm, E = E, diag = diagonal)
}

dbglm_fe_multiply = function(sys, X) {
  x1 = X[sys$big, , drop = FALSE]
  x2 = X[sys$rest, , drop = FALSE]
  out = matrix(0, sys$n, ncol(X))
  out[sys$big, ] = sys$Dg * x1 + sys$Cm %*% x2
  out[sys$rest, ] = crossprod(sys$Cm, x1) + sys$E %*% x2
  out
}

# Factor the Schur complement on the non-absorbed levels; the identified basis
# is chosen once per fit so FE normalisations stay fixed across iterations.
dbglm_fe_solver = function(sys, cache) {
  if (any(!is.finite(sys$Dg)) || any(sys$Dg <= 0)) stop("GLM fixed-effect information is singular.", call. = FALSE)
  scaled = sys$Cm / sys$Dg
  S = sys$E - crossprod(sys$Cm, scaled)
  sdiag = sqrt(pmax(diag(S), 0))
  if (is.null(cache$moment_basis)) {
    positive = which(sdiag > 0 & is.finite(sdiag))
    keep = integer()
    if (length(positive)) {
      C = S[positive, positive, drop = FALSE] / outer(sdiag[positive], sdiag[positive])
      pivoted = suppressWarnings(chol(C, pivot = TRUE, tol = 1e-11))
      keep = positive[attr(pivoted, "pivot")[seq_len(attr(pivoted, "rank"))]]
    }
    if (!length(keep) && !length(sys$big)) stop("No identified fixed-effect span remains.", call. = FALSE)
    cache$moment_basis = keep
  }
  keep = cache$moment_basis
  scale = sdiag[keep]
  factor = if (length(keep)) tryCatch(chol(S[keep, keep, drop = FALSE] / outer(scale, scale)), error = function(e) NULL) else matrix(0, 0, 0)
  if (is.null(factor)) stop("GLM fixed-effect information is singular.", call. = FALSE)
  # B lies in range(H) exactly, but plain SQL sums (e.g. Spark) leave roundoff
  # in the redundant FE directions that no solve can remove. Null vectors of
  # H follow from the dropped Schur pivots: v2[keep] = -S_kk^{-1} S_kd.
  dropped = setdiff(seq_along(sys$rest), keep)
  null_q = NULL
  if (length(dropped)) {
    v2 = matrix(0, length(sys$rest), length(dropped))
    v2[cbind(dropped, seq_along(dropped))] = 1
    if (length(keep)) {
      v2[keep, ] = -backsolve(factor, backsolve(factor, S[keep, dropped, drop = FALSE] / scale, transpose = TRUE)) / scale
    }
    V = matrix(0, sys$n, length(dropped))
    V[sys$rest, ] = v2
    V[sys$big, ] = -scaled %*% v2
    null_q = qr.Q(qr(V))
  }
  solve = function(rhs) {
    r1 = rhs[sys$big, , drop = FALSE]
    r2 = rhs[sys$rest, , drop = FALSE] - crossprod(scaled, r1)
    x2 = matrix(0, length(sys$rest), ncol(rhs))
    if (length(keep)) {
      x2[keep, ] = backsolve(factor, backsolve(factor, r2[keep, , drop = FALSE] / scale, transpose = TRUE)) / scale
    }
    out = matrix(0, sys$n, ncol(rhs))
    out[sys$rest, ] = x2
    out[sys$big, ] = r1 / sys$Dg - scaled %*% x2
    out
  }
  project = function(rhs) if (is.null(null_q)) rhs else rhs - null_q %*% crossprod(null_q, rhs)
  list(solve = solve, project = project)
}

# One statement returns the global row (__set = 0), per-FE marginals
# (__set = i) and pairwise FE information (__set = K + pair index). Pairs need
# only the weight sum, so each is a narrow GROUP BY over one shared source;
# GROUPING SETS would evaluate every aggregate for every set.
dbglm_fe_moments_sql = function(conn, source, fe, value_terms, global_terms = character(), limit = NULL) {
  K = length(fe)
  pairs = utils::combn(K, 2L, simplify = FALSE)
  alias = function(terms) sub(".* AS ([A-Za-z_0-9]+)$", "\\1", terms)
  null = function(terms) if (length(terms)) paste0("CAST(NULL AS DOUBLE) AS ", alias(terms)) else character()
  id = function(f) if (is.null(f)) "CAST(NULL AS BIGINT)" else paste0("CAST(", f, " AS BIGINT)")
  branch = function(set, i, j, columns, group) paste0("SELECT ", set, " AS __set,", id(i), " AS __i,", id(j), " AS __j,",
    paste(columns, collapse = ","), " FROM __dbglm_src", if (length(group)) paste(" GROUP BY", paste(group, collapse = ",")))
  branches = c(
    if (length(global_terms)) branch(0L, NULL, NULL, c(global_terms, null(value_terms)), character()),
    vapply(seq_len(K), function(i) branch(i, fe[i], NULL, c(null(global_terms), value_terms), fe[i]), character(1)),
    vapply(seq_along(pairs), function(k) {
      p = pairs[[k]]
      branch(K + k, fe[p[1]], fe[p[2]], c(null(global_terms), value_terms[1L], null(value_terms[-1L])), fe[p])
    }, character(1)))
  materialized = if (isTRUE(detect_backend(conn)$is_databricks)) "" else "MATERIALIZED "
  sql = paste0("WITH __dbglm_src AS ", materialized, "(SELECT * FROM ", source, ") ", paste(branches, collapse = " UNION ALL "))
  if (!is.null(limit)) sql = paste0("SELECT * FROM (", sql, ") __dbglm_m ORDER BY __set LIMIT ", sprintf("%.0f", limit))
  sql
}

# Value columns are __h, __b1.., __ex1..; see dbglm_fe_moments_sql().
dbglm_unpack_fe_moments = function(res, K, n_cols, n_extra) {
  sets = as.integer(res[["__set"]])
  b_cols = if (n_cols) paste0("__b", seq_len(n_cols)) else character()
  ex_cols = if (n_extra) paste0("__ex", seq_len(n_extra)) else character()
  info = vector("list", K)
  for (i in seq_len(K)) {
    rows = which(sets == i)
    if (!length(rows)) return(NULL)
    rows = rows[order(res[["__i"]][rows])]
    h = as.numeric(res[["__h"]][rows])
    b = as.matrix(res[rows, b_cols, drop = FALSE])
    ex = as.matrix(res[rows, ex_cols, drop = FALSE])
    if (anyNA(res[["__i"]][rows]) || any(!is.finite(h)) || any(h <= 0) || any(!is.finite(b)) || any(!is.finite(ex))) return(NULL)
    info[[i]] = list(id = res[["__i"]][rows], h = h, b = unname(b), ex = unname(ex))
  }
  combos = utils::combn(K, 2L, simplify = FALSE)
  pairs = vector("list", length(combos))
  for (k in seq_along(combos)) {
    p = combos[[k]]
    rows = which(sets == K + k)
    if (!length(rows)) return(NULL)
    h = as.numeric(res[["__h"]][rows])
    ir = match(res[["__j"]][rows], info[[p[2]]]$id)
    jr = match(res[["__i"]][rows], info[[p[1]]]$id)
    if (any(!is.finite(h)) || anyNA(ir) || anyNA(jr)) return(NULL)
    pairs[[k]] = list(i = p[2], j = p[1], ir = ir, jr = jr, h = h)
  }
  sys = dbglm_fe_system(info, pairs)
  if (is.null(sys)) return(NULL)
  B = do.call(rbind, lapply(info, `[[`, "b"))
  B_extra = do.call(rbind, lapply(info, `[[`, "ex"))
  fe_score_max = if (n_extra >= 1L) max(abs(B_extra[, 1L]) / sqrt(pmax(sys$diag, .Machine$double.xmin))) else 0
  list(info = info, sys = sys, B = B, B_extra = B_extra, starts = sys$starts, n = sys$n, fe_score_max = fe_score_max)
}

dbglm_fe_value_terms = function(conn, columns, extra_sums = character()) {
  c("FSUM(w) AS __h",
    if (length(columns)) sprintf("FSUM(w*(%s)) AS __b%d", columns, seq_along(columns)),
    if (length(extra_sums)) sprintf("%s AS __ex%d", vapply(extra_sums, function(e) dbglm_exact_sum(conn, e), character(1)), seq_along(extra_sums)))
}

dbglm_query_fe_moments = function(conn, table, fe, columns, extra_sums = character(),
                                  max_edges = getOption("dbreg.duckdb_ap_max_edges", 1e6),
                                  max_entries = getOption("dbreg.duckdb_ap_max_entries", 2e6)) {
  valid_cap = function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 1 && x == floor(x)
  }
  if (!valid_cap(max_edges) || !valid_cap(max_entries)) {
    stop("DuckDB AP transfer limits must be positive whole numbers.", call. = FALSE)
  }
  sql = dbglm_fe_moments_sql(conn, table, fe, dbglm_fe_value_terms(conn, columns, extra_sums), limit = max_edges + 1)
  res = dbglm_get_query(conn, sql)
  if (!nrow(res) || nrow(res) > max_edges) return(NULL)
  unpacked = dbglm_unpack_fe_moments(res, length(fe), length(columns), length(extra_sums))
  if (is.null(unpacked)) return(NULL)
  if (unpacked$n * max(1L, length(columns)) > max_entries) return(NULL)
  unpacked
}

dbglm_project_moments = function(conn, table, columns, identifiers, ops, scales, target, control, create, drop, cache,
                                 materialize = TRUE, extra_sums = character(), need_table = TRUE,
                                 precomputed = NULL) {
  all_intercepts = all(vapply(ops, function(op) identical(op$v, "1.0"), logical(1)))
  n_cols = length(columns)
  n_extra = length(extra_sums)
  if (!is.null(precomputed)) {
    info = lapply(precomputed$info, function(x) list(id = x$id, h = x$h))
    sys = precomputed$sys
    B = if (n_cols) precomputed$B[, seq_len(n_cols), drop = FALSE] else matrix(0, precomputed$n, 0L)
    B_extra = if (n_extra) precomputed$B_extra[, seq_len(n_extra), drop = FALSE] else matrix(0, precomputed$n, 0L)
    n = precomputed$n
  } else if (all_intercepts && length(ops) >= 2L) {
    fe = vapply(ops, `[[`, character(1), "f")
    unpacked = dbglm_query_fe_moments(conn, table, fe, columns, extra_sums)
    if (is.null(unpacked)) return(NULL)
    info = lapply(unpacked$info, function(x) list(id = x$id, h = x$h))
    sys = unpacked$sys
    B = unpacked$B
    B_extra = unpacked$B_extra
    n = unpacked$n
  } else {
    info = lapply(ops, function(op) dbglm_get_query(conn, paste0("SELECT ", op$f, " AS id,FSUM(w*(", op$v, ")*(", op$v, ")) AS h",
      if (n_cols) paste0(",", paste(sprintf("FSUM(w*(%s)*%s) AS b%d", op$v, columns, seq_len(n_cols)), collapse = ",")),
      if (n_extra) paste0(",", paste(sprintf("FSUM((%s)*%s) AS ex%d", op$v, extra_sums, seq_len(n_extra)), collapse = ",")),
      " FROM ", table, " GROUP BY ", op$f, " ORDER BY ", op$f)))
    pairs = list()
    for (i in seq_along(ops)) if (i > 1) for (j in seq_len(i - 1L)) {
      oi = ops[[i]]; oj = ops[[j]]
      edge = dbglm_get_query(conn, paste0("SELECT ", oi$f, " AS i,", oj$f, " AS j,FSUM(w*(", oi$v, ")*(", oj$v, ")) AS h FROM ",
        table, " GROUP BY ", oi$f, ",", oj$f))
      pairs[[length(pairs) + 1L]] = list(i = i, j = j, h = as.numeric(edge$h),
        ir = match(edge$i, info[[i]]$id), jr = match(edge$j, info[[j]]$id))
    }
    sys = dbglm_fe_system(lapply(info, function(x) list(id = x$id, h = as.numeric(x$h))), pairs)
    if (is.null(sys)) return(NULL)
    n = sys$n
    B = if (n_cols) do.call(rbind, lapply(info, function(x) as.matrix(x[, paste0("b", seq_len(n_cols)), drop = FALSE]))) else matrix(0, n, 0L)
    B_extra = if (n_extra) do.call(rbind, lapply(info, function(x) as.matrix(x[, paste0("ex", seq_len(n_extra)), drop = FALSE]))) else matrix(0, n, 0L)
  }
  ids = unlist(lapply(seq_along(ops), function(i) paste(i, info[[i]]$id, sep = "/")), use.names = FALSE)
  if (!is.null(cache$moment_keys) && !identical(ids, cache$moment_keys)) stop("The fixed-effect sample changed within fitting.", call. = FALSE)
  cache$moment_keys = ids
  diagonal = sqrt(pmax(sys$diag, 0))
  solver = dbglm_fe_solver(sys, cache)
  solve_rhs = solver$solve
  if (n_cols) B = solver$project(B)
  coefficients = matrix(0, n, n_cols)
  if (n_cols) coefficients = solve_rhs(B)
  residual = 0
  residual_rhs = matrix(0, n, n_cols)
  if (n_cols) {
    max_cycles = min(control$fe_max_iter, 8L)
    for (cycle in seq_len(max_cycles)) {
      residual_rhs = B - dbglm_fe_multiply(sys, coefficients)
      normalized = residual_rhs / pmax(diagonal, .Machine$double.xmin)
      residual = max(abs(sweep(normalized, 2, scales, "/")))
      if (is.finite(residual) && residual <= target) break
      if (cycle == max_cycles) {
        stop("GLM fixed-effect projection did not converge to its group-equation tolerance.", call. = FALSE)
      }
      coefficients = coefficients + solve_rhs(residual_rhs)
    }
  } else {
    cycle = 1L
  }
  cross_reduction = if (n_cols) {
    cr = crossprod(B, coefficients) + crossprod(residual_rhs, coefficients)
    0.5 * (cr + t(cr))
  } else matrix(0, 0L, 0L)
  extra_reduction = if (n_extra && n_cols) crossprod(coefficients, B_extra) else NULL
  if (!need_table) {
    return(list(table = NULL, raw_projected = NULL, cycles = cycle, residual = residual, method = "moments",
      cross_reduction = cross_reduction, extra_reduction = extra_reduction,
      coefficients = coefficients, info = info, ops = ops, sub_tables = character()))
  }
  if (!n_cols) {
    projected = create(paste("SELECT __dbglm_key,", paste(c(identifiers, "w"), collapse = ","), "FROM", table), materialize = materialize)
    return(list(table = projected, raw_projected = NULL, cycles = cycle, residual = 0, method = "moments",
      cross_reduction = cross_reduction, extra_reduction = extra_reduction,
      coefficients = coefficients, info = info, ops = ops, sub_tables = character()))
  }
  frame = data.frame(op = rep(seq_along(ops), vapply(info, function(x) length(x$id), integer(1))), check.names = FALSE)
  names(frame) = "__glm_op"
  frame[["__glm_group"]] = as.numeric(unlist(lapply(info, `[[`, "id"), use.names = FALSE))
  frame[paste0("c", seq_len(n_cols))] = as.data.frame(coefficients)
  is_dbx = isTRUE(detect_backend(conn)$is_databricks)
  inline_intermediate = !materialize && !isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates", FALSE))
  lookup_sub = character()
  if (is_dbx) {
    lookup_sql = dbglm_values_sql(conn, frame)
    if (is.null(lookup_sql)) stop("FE coefficient lookup exceeds the SQL transfer limit.", call. = FALSE)
    if (inline_intermediate) {
      lookup = paste0("(", lookup_sql, ")")
    } else {
      lookup = create(lookup_sql, materialize = TRUE)
      lookup_sub = lookup
    }
  } else {
    registered = basename(tempfile("dbglm_coef_"))
    tryCatch({
      duckdb::duckdb_register(conn, registered, frame)
      lookup = create(paste("SELECT * FROM", DBI::dbQuoteIdentifier(conn, registered)), materialize = TRUE)
    }, finally = try(duckdb::duckdb_unregister(conn, registered), silent = TRUE))
    lookup_sub = lookup
  }
  ctes = vapply(seq_along(ops), function(i) paste0("proj_g", i, " AS (SELECT * FROM ", lookup, " WHERE __glm_op=", i, ")"), character(1))
  joins = vapply(seq_along(ops), function(i) paste0(" LEFT JOIN proj_g", i, " ON s.", ops[[i]]$f, "=proj_g", i, ".__glm_group"), character(1))
  broadcast = if (is_dbx) paste0("/*+ BROADCAST(", paste0("proj_g", seq_along(ops), collapse = ","), ") */ ") else ""
  if (inline_intermediate) {
    tilde_exprs = vapply(seq_len(n_cols), function(j) {
      terms = vapply(seq_along(ops), function(i) {
        v = if (ops[[i]]$v == "1.0") "1.0" else paste0("s.", ops[[i]]$v)
        sprintf("(%s)*COALESCE(proj_g%d.c%d,0)", v, i, j)
      }, character(1))
      paste0("CAST((s.", columns[j], "-(", paste(terms, collapse = " + "), ")) AS DOUBLE) AS ", columns[j], "_tilde")
    }, character(1))
    raw_sql = paste0("WITH ", paste(ctes, collapse = ","), " SELECT ", broadcast, "s.*, s.w AS __w, ",
      paste(tilde_exprs, collapse = ","), " FROM ", table, " s", paste(joins, collapse = ""))
    raw_projected = paste0("(", dbglm_sql(conn, raw_sql), ")")
    return(list(table = NULL, raw_projected = raw_projected, cycles = cycle, residual = residual, method = "moments",
      cross_reduction = cross_reduction, extra_reduction = extra_reduction,
      coefficients = coefficients, info = info, ops = ops, sub_tables = lookup_sub))
  }
  expressions = vapply(seq_len(n_cols), function(j) {
    terms = vapply(seq_along(ops), function(i) {
      v = if (ops[[i]]$v == "1.0") "1.0" else paste0("s.", ops[[i]]$v)
      sprintf("(%s)*COALESCE(proj_g%d.c%d,0)", v, i, j)
    }, character(1))
    paste0("CAST((s.", columns[j], "-(", paste(terms, collapse = " + "), ")) AS DOUBLE) AS ", columns[j])
  }, character(1))
  projected = create(paste0("WITH ", paste(ctes, collapse = ","), " SELECT ", broadcast,
    paste(paste0("s.", c("__dbglm_key", identifiers, "w")), collapse = ","), ",",
    paste(expressions, collapse = ","), " FROM ", table, " s", paste(joins, collapse = "")), materialize = TRUE)
  for (st in lookup_sub) drop(st)
  list(table = projected, raw_projected = NULL, cycles = cycle, residual = residual, method = "moments",
    cross_reduction = cross_reduction, extra_reduction = extra_reduction,
    coefficients = coefficients, info = info, ops = ops, sub_tables = character())
}
