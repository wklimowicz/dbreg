# SQL owns row state; bounded FE moments or FE-level CG determine projections.
solve_partial_wls = function(conn, table, x, fe, control, create, drop, response = TRUE, cache = NULL, materialize = TRUE, scales = NULL, extra_sums = character(), need_table = TRUE, precomputed = NULL) {
  columns = c(x, if (response) "z")
  backend_info = detect_backend(conn)
  slopes = attr(fe,"slopes")
  if(is.null(cache)) cache = new.env(parent=emptyenv())
  if(length(fe)>1) {
    dbglm_moments_enabled(conn, table, fe, control, cache)
    if (!cache$moment_enabled && length(fe) == 2L && !length(unlist(slopes)) &&
        identical(control$projection_method, "auto") && length(columns)) {
      identifiers = unique(c(fe, unlist(slopes)))
      if (is.null(scales)) {
        norm = as.numeric(dbglm_get_query(conn,paste("SELECT",paste(sprintf(
          "FSUM(w*POWER(%s,2)) AS n%d",columns,seq_along(columns)),collapse=","),"FROM",table))[1,])
        scales = sqrt(pmax(norm,.Machine$double.xmin))
      }
      target = control$fe_tol
      if(response) target = min(target,control$score_tol/(100*max(1,utils::tail(scales,1))))
      first = utils::tail(columns,1)
      rest = utils::head(columns,-1)
      order = match(c(first,rest),columns)
      projected = dbreg_duckdb_cg_projections(conn,paste("FROM",table),fe,first,rest,rest,"w",
        cluster_var="__dbglm_key",max_iter=control$fe_max_iter,tol=target,scales=scales[order],
        materialize_result=FALSE, inline_result=!materialize, extra_sums=extra_sums)
      if(!is.null(projected)) {
        on.exit(for(t in rev(projected$all_tables)) {
          if (!materialize && isTRUE(backend_info$is_databricks) && !isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates",FALSE))) next
          if(t %in% projected$view_tables) try(DBI::dbExecute(conn,paste("DROP VIEW IF EXISTS",t)),silent=TRUE) else
            drop_ap_table_if_exists(conn,t,backend_info)
        },add=TRUE)
        if (!materialize && isTRUE(backend_info$is_databricks) && !isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates",FALSE))) {
          cache$views = unique(c(cache$views, projected$view_tables))
          cache$tables = unique(c(cache$tables, setdiff(projected$all_tables, projected$view_tables)))
        }
        result = create(paste("SELECT __dbglm_key,",paste(c(identifiers,"__w AS w",
          sprintf("%s_tilde AS %s",columns,columns)),collapse=","),"FROM",projected$table), materialize=materialize)
        inv_order = match(columns, c(first, rest))
        reduction = if (!is.null(projected$cross_reduction)) projected$cross_reduction[inv_order, inv_order, drop = FALSE] else NULL
        extra_red = if (!is.null(projected$extra_reduction)) projected$extra_reduction[inv_order, , drop = FALSE] else NULL
        raw_proj = if (!materialize && isTRUE(backend_info$is_databricks) &&
          !isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates",FALSE))) projected$table else NULL
        return(list(table=result,raw_projected=raw_proj,cycles=projected$iterations,residual=projected$max_abs,method="bounded_pcg",
          cross_reduction=reduction,extra_reduction=extra_red,sub_tables=if (!materialize) projected$all_tables else character(),
          s1=projected$s1[, inv_order, drop=FALSE], s2=projected$s2[, inv_order, drop=FALSE],
          a=projected$a[, inv_order, drop=FALSE], b=projected$b[, inv_order, drop=FALSE]))
      }
    }
    if(cache$moment_enabled) {
      ops = lapply(fe,function(f) list(f=f,v="1.0"))
      for(i in seq_along(fe)) for(v in slopes[[i]]) ops[[length(ops)+1]] = list(f=fe[i],v=v)
      identifiers = unique(c(fe,unlist(slopes)))
      if (is.null(scales)) {
        norm = if(length(columns)) as.numeric(dbglm_get_query(conn,paste("SELECT",paste(sprintf("FSUM(w*POWER(%s,2)) AS n%d",columns,seq_along(columns)),collapse=","),"FROM",table))[1,]) else numeric()
        scales = sqrt(pmax(norm,.Machine$double.xmin))
      }
      target = control$fe_tol
      if(response) target = min(target,control$score_tol/(100*max(1,utils::tail(scales,1))))
      ans = dbglm_project_moments(conn,table,columns,identifiers,ops,scales,target,control,create,drop,cache,
        materialize=materialize,extra_sums=extra_sums,need_table=need_table,precomputed=precomputed)
      if(!is.null(ans)) return(ans)
      cache$moment_enabled = FALSE
    }
    if(isTRUE(backend_info$is_databricks) && !length(unlist(slopes))) {
      return(dbglm_project_window(conn,table,columns,fe,control,create,drop,response))
    }
  }
  projection = create(paste("SELECT __dbglm_key,", paste(unique(c(fe, unlist(slopes), "w", columns)), collapse = ","), "FROM", table))
  ops = lapply(fe,function(f) list(f=f,v="1.0"))
  basis = character()
  for(i in seq_along(fe)) {
    local_basis = character()
    for(j in seq_along(slopes[[i]])) {
      v = slopes[[i]][j]
      name = paste0("__dbglm_v",i,"_",j)
      partition = paste0(" OVER (PARTITION BY ",fe[i],")")
      centered = paste0(v," - SUM(w*",v,")",partition,"/SUM(w)",partition)
      next_projection = create(paste("SELECT *,",centered,"AS",name,"FROM",projection))
      drop(projection); projection = next_projection
      for(prev in local_basis) {
        adjustment = paste0(name," - ",prev,"*COALESCE(SUM(w*",prev,"*",name,")",partition,
          "/NULLIF(SUM(w*",prev,"*",prev,")",partition,",0),0)")
        next_projection = create(paste("SELECT * EXCLUDE(",name,"),",adjustment,"AS",name,"FROM",projection))
        drop(projection); projection = next_projection
      }
      local_basis = c(local_basis,name)
      basis = c(basis,name)
      ops[[length(ops)+1]] = list(f=fe[i],v=name)
    }
  }
  identifiers = unique(c(fe,unlist(slopes),basis))
  if (!length(columns)) return(list(table = projection, cycles = 0L, residual = 0, method = "group"))
  norm = dbglm_get_query(conn, paste("SELECT", paste(sprintf("SUM(w*POWER(%s,2)) AS n%d", columns, seq_along(columns)), collapse=","), "FROM", projection))
  scales = sqrt(pmax(as.numeric(norm[1, ]), .Machine$double.xmin))
  target = control$fe_tol
  if(response) target = min(target,control$score_tol/(100*max(1,utils::tail(scales,1))))
  if(length(fe)>1) return(dbglm_project_cg(conn,projection,columns,identifiers,ops,scales,target,control,create,drop))
  residual = 0
  for (cycle in seq_len(control$fe_max_iter)) {
    for (op in ops) {
      f = op$f; v = op$v
      values = sprintf("(%1$s - (%3$s)*COALESCE(SUM(w*(%3$s)*%1$s) OVER (PARTITION BY %2$s)/NULLIF(SUM(w*(%3$s)*(%3$s)) OVER (PARTITION BY %2$s),0),0)) AS %1$s", columns, f, v)
      next_projection = create(paste("SELECT __dbglm_key,", paste(c(identifiers, "w", values), collapse=","), "FROM", projection))
      drop(projection)
      projection = next_projection
    }
    residual = 0
    for (op in ops) {
      f = op$f; v = op$v
      terms = sprintf("COALESCE(ABS(SUM(w*(%3$s)*%1$s))/NULLIF(SQRT(SUM(w*(%3$s)*(%3$s))),0),0)/%2$s AS r%4$d", columns, dbglm_number(scales), v, seq_along(columns))
      row = dbglm_get_query(conn, paste0("SELECT ", paste(sprintf("MAX(r%d) AS r%d", seq_along(columns), seq_along(columns)), collapse=","),
        " FROM (SELECT ", paste(terms, collapse=","), " FROM ", projection, " GROUP BY ", f, ") t"))
      residual = max(residual, as.numeric(row[1, ]))
    }
    if (residual <= target) return(list(table = projection, cycles = cycle, residual = residual, method = "group"))
  }
  stop("GLM fixed-effect projection did not converge.", call. = FALSE)
}

# Direct FE moments need the non-absorbed levels to fit a dense solve in R;
# see dbglm_fe_dense_feasible(). The decision is cached for the whole fit.
dbglm_moments_enabled = function(conn, table, fe, control, cache) {
  if (!is.null(cache$moment_enabled)) return(cache$moment_enabled)
  if (identical(control$projection_method, "sql")) {
    cache$moment_enabled = FALSE
    return(FALSE)
  }
  sizes = cache$fe_sizes
  if (is.null(sizes)) {
    count = dbglm_get_query(conn, paste("SELECT", paste(sprintf("COUNT(DISTINCT %s) AS n%d", fe, seq_along(fe)), collapse = ","), "FROM", table))
    sizes = vapply(count, as.numeric, numeric(1))
  }
  cache$moment_enabled = dbglm_fe_dense_feasible(rep(sizes, 1 + lengths(attr(fe, "slopes"))))
  cache$moment_enabled
}

dbglm_fe_working = function(conn, table, link, create, target = NULL) {
  create(paste0(dbglm_working_sql(table, character(), numeric(), link, eta = "__dbglm_eta"), "SELECT * FROM working"), target = target)
}

dbglm_fe_evaluate = function(conn, working, x, fe, collect_fe_moments = FALSE) {
  cols = c(x, "z")
  terms = c(dbglm_deviance_term(conn), "MIN(mu) AS mu_min", "MAX(mu) AS mu_max",
    "SUM(CASE WHEN NOT isfinite(w) OR w<=0 OR NOT isfinite(eta) THEN 1 ELSE 0 END) AS bad",
    dbglm_moment_terms(cols, "w"))
  if (length(x)) {
    x_score = vapply(paste0(x, "*score"), function(e) dbglm_exact_sum(conn, e), character(1))
    terms = c(terms, sprintf("%s AS raw_s%d", x_score, seq_along(x)),
      sprintf("ABS(%s)/SQRT(SUM(w*%s*%s)) AS s%d", x_score, x, x, seq_along(x)))
  }
  if (collect_fe_moments && length(fe) >= 2L && !length(unlist(attr(fe, "slopes")))) {
    max_edges = getOption("dbreg.duckdb_ap_max_edges", 1e6)
    max_entries = getOption("dbreg.duckdb_ap_max_entries", 2e6)
    K = length(fe)
    res = dbglm_get_query(conn, dbglm_fe_moments_sql(conn, working, fe, dbglm_fe_value_terms(conn, cols, "score"), terms))
    sets = as.integer(res[["__set"]])
    row = res[sets == 0L, , drop = FALSE]
    score = if (length(x)) max(as.numeric(row[1, paste0("s", seq_along(x)), drop = FALSE])) else 0
    raw_score = if (length(x)) as.numeric(row[1, paste0("raw_s", seq_along(x)), drop = FALSE]) else numeric()
    raw_moments = dbglm_read_information(row, length(cols))
    fe_moments = NULL
    if (!is.na(row$bad) && row$bad == 0 && sum(sets != 0L) <= max_edges) {
      unpacked = dbglm_unpack_fe_moments(res, K, length(cols), 1L)
      if (!is.null(unpacked) && unpacked$n * max(1L, length(cols)) <= max_entries) fe_moments = unpacked
    }
    singles = sets >= 1L & sets <= K
    h = as.numeric(res[["__h"]][singles])
    s_val = ifelse(is.finite(h) & h > 0 & is.finite(res[["__ex1"]][singles]), abs(res[["__ex1"]][singles]) / sqrt(h), 0)
    if (length(s_val)) score = max(score, s_val)
    return(list(deviance = row$deviance, mu_min = row$mu_min, mu_max = row$mu_max, bad = row$bad,
                score_max = score, raw_score = raw_score, raw_moments = raw_moments,
                scales = sqrt(pmax(diag(raw_moments), .Machine$double.xmin)),
                fe_moments = fe_moments))
  }
  fe_queries = character()
  for (i in seq_along(fe)) for(v in c("1.0",attr(fe,"slopes")[[i]])) {
    fe_queries = c(fe_queries, paste0("SELECT COALESCE(ABS(",dbglm_exact_sum(conn,paste0("score*(",v,")")),")/NULLIF(SQRT(SUM(w*(",v,")*(",v,"))),0),0) AS s FROM ", working, " GROUP BY ", fe[i]))
  }
  if (length(fe_queries)) {
    terms = c(terms, paste0("(SELECT MAX(s) FROM (", paste(fe_queries, collapse = " UNION ALL "), ") t) AS fe_score_max"))
  }
  row = dbglm_get_query(conn, paste("SELECT", paste(terms, collapse=","), "FROM", working))
  score = if (length(x)) max(as.numeric(row[1, paste0("s", seq_along(x)), drop=FALSE])) else 0
  raw_score = if (length(x)) as.numeric(row[1, paste0("raw_s", seq_along(x)), drop=FALSE]) else numeric()
  if (length(fe_queries) && !is.na(row$fe_score_max)) {
    score = max(score, as.numeric(row$fe_score_max))
  }
  raw_moments = dbglm_read_information(row, length(cols))
  list(deviance = row$deviance, mu_min = row$mu_min, mu_max = row$mu_max, bad = row$bad,
       score_max = score, raw_score = raw_score, raw_moments = raw_moments,
       scales = sqrt(pmax(diag(raw_moments), .Machine$double.xmin)))
}

dbglm_fe_fit = function(conn, base, x, fe, beta, family, control, create, drop, verbose, cache = NULL, ybar = NULL, vcov_type = "iid") {
  if (is.null(ybar)) ybar = dbglm_get_query(conn, paste("SELECT SUM(__dbglm_a*__dbglm_y)/SUM(__dbglm_a) AS y FROM", base))$y
  intercept = family$linkfun(ybar)
  backend_info = detect_backend(conn)
  slopes = attr(fe, "slopes")
  if (is.null(cache)) cache = new.env(parent = emptyenv())
  if (length(fe) > 1L) dbglm_moments_enabled(conn, base, fe, control, cache)
  use_moments = length(fe) > 1L && isTRUE(cache$moment_enabled)
  collect_fe_moments = use_moments && !length(unlist(slopes))
  inline_sql = function(sql) {
    if (isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates", FALSE))) {
      create(sql, materialize = TRUE)
    } else {
      paste0("(", dbglm_sql(conn, sql), ")")
    }
  }
  make_working = function(beta_val, c_eta_val, gamma_val = NULL, ops = NULL, info = NULL) {
    lookup_sub = NULL
    if (is.null(gamma_val) || all(gamma_val == 0)) {
      eta = paste(c("__dbglm_o", if (c_eta_val != 0) dbglm_number(c_eta_val), sprintf("(%s)*%s", dbglm_number(beta_val), x)), collapse = " + ")
      st = if (use_moments) inline_sql(paste("SELECT *,", eta, "AS __dbglm_eta FROM", base)) else
        create(paste("SELECT *,", eta, "AS __dbglm_eta FROM", base), materialize = FALSE)
    } else {
      frame = data.frame(op = rep(seq_along(ops), vapply(info, function(inf) length(inf$id), integer(1))), check.names = FALSE)
      names(frame) = "__glm_op"
      frame[["__glm_group"]] = as.numeric(unlist(lapply(info, `[[`, "id"), use.names = FALSE))
      frame[["c1"]] = as.numeric(gamma_val)
      if (isTRUE(backend_info$is_databricks)) {
        lsql = dbglm_values_sql(conn, frame)
        if (is.null(lsql)) stop("FE coefficient lookup exceeds the SQL transfer limit.", call. = FALSE)
        if (!isTRUE(getOption("dbreg.glm_databricks_materialize_intermediates", FALSE))) {
          lookup = paste0("(", lsql, ")")
        } else {
          lookup = create(lsql, materialize = TRUE)
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
      ctes = vapply(seq_along(ops), function(i) paste0("eta_g", i, " AS (SELECT __glm_group, c1 FROM ", lookup, " WHERE __glm_op=", i, ")"), character(1))
      fe_terms = vapply(seq_along(ops), function(i) {
        v = if (ops[[i]]$v == "1.0") "1.0" else paste0("s.", ops[[i]]$v)
        if (v == "1.0") sprintf("COALESCE(eta_g%d.c1,0.0)", i) else sprintf("(%s)*COALESCE(eta_g%d.c1,0.0)", v, i)
      }, character(1))
      base_terms = c("s.__dbglm_o", if (c_eta_val != 0) dbglm_number(c_eta_val), sprintf("(%s)*s.%s", dbglm_number(beta_val), x), fe_terms)
      joins = vapply(seq_along(ops), function(i) paste0(" LEFT JOIN eta_g", i, " ON s.", ops[[i]]$f, "=eta_g", i, ".__glm_group"), character(1))
      broadcast = if (isTRUE(backend_info$is_databricks)) paste0("/*+ BROADCAST(", paste0("eta_g", seq_along(ops), collapse = ","), ") */ ") else ""
      st = inline_sql(paste0("WITH ", paste(ctes, collapse = ","), " SELECT ", broadcast, "s.*, CAST((", paste(base_terms, collapse = " + "), ") AS DOUBLE) AS __dbglm_eta FROM ", base, " s", paste(joins, collapse = "")))
    }
    wk = if (use_moments) inline_sql(paste0(dbglm_working_sql(st, character(), numeric(), family$link, eta = "__dbglm_eta"), "SELECT * FROM working")) else
      dbglm_fe_working(conn, st, family$link, create)
    list(state = st, working = wk, lookup = lookup_sub)
  }
  c_eta = intercept
  gamma = NULL
  wk_obj = make_working(beta, c_eta)
  state = wk_obj$state; working = wk_obj$working; cur_lookup = wk_obj$lookup
  current = dbglm_fe_evaluate(conn, working, x, fe, collect_fe_moments = collect_fe_moments)
  initial_halvings = 0L
  while(current$bad>0 || dbglm_boundary(current,family$link)) {
    if(initial_halvings>=control$step_max_halvings) stop("Initial GLM predictor reaches a numerical boundary; choose smaller starts or offsets.",call.=FALSE)
    beta = beta/2
    drop(working);drop(state);drop(cur_lookup)
    wk_obj = make_working(beta, c_eta)
    state = wk_obj$state; working = wk_obj$working; cur_lookup = wk_obj$lookup
    current = dbglm_fe_evaluate(conn,working,x,fe,collect_fe_moments=collect_fe_moments)
    initial_halvings = initial_halvings+1L
  }
  if (collect_fe_moments && is.null(current$fe_moments)) {
    cache$moment_enabled = FALSE
    use_moments = FALSE
    collect_fe_moments = FALSE
    drop(working);drop(state);drop(cur_lookup)
    wk_obj = make_working(beta, c_eta)
    state = wk_obj$state; working = wk_obj$working; cur_lookup = wk_obj$lookup
  }
  history = current$deviance
  total_halvings = initial_halvings
  total_cycles = 0L
  alt_working = NULL
  for (iter in seq_len(control$glm_max_iter)) {
    projected = solve_partial_wls(conn, working, x, fe, control, create, drop, cache=cache,
      materialize=FALSE, scales=current$scales, need_table=!use_moments, precomputed=current$fe_moments)
    total_cycles = total_cycles + projected$cycles
    p = length(x)
    if (p) {
      if (!is.null(projected$cross_reduction)) {
        M_tilde = current$raw_moments - projected$cross_reduction
        H_step = M_tilde[seq_len(p), seq_len(p), drop=FALSE]
        b_step = as.numeric(M_tilde[seq_len(p), p + 1L])
        candidate = dbglm_solve(H_step, b_step)$beta
      } else {
        row = dbglm_get_query(conn, paste("SELECT", paste(c(dbglm_moment_terms(x,"w"), sprintf("SUM(w*%s*z) AS b_%d",x,seq_along(x))),collapse=","), "FROM", projected$table))
        candidate = dbglm_solve(dbglm_read_information(row,p), vapply(seq_len(p),function(i) as.numeric(row[[paste0("b_",i)]]),numeric(1)))$beta
      }
    } else candidate = numeric()
    if (!is.null(projected$coefficients)) {
      alpha_z = projected$coefficients[, p + 1L]
      gamma_cand = if (p) as.numeric(alpha_z - projected$coefficients[, seq_len(p), drop = FALSE] %*% candidate) else as.numeric(alpha_z)
      if (is.null(gamma)) gamma = numeric(length(gamma_cand))
      step = 1
      accepted = FALSE
      for (halvings in 0:control$step_max_halvings) {
        next_beta = beta + step * (candidate - beta)
        next_gamma = gamma + step * (gamma_cand - gamma)
        next_c_eta = (1 - step) * c_eta
        next_wk = make_working(next_beta, next_c_eta, next_gamma, projected$ops, projected$info)
        next_state = next_wk$state; next_working = next_wk$working; next_lookup = next_wk$lookup
        next_eval = dbglm_fe_evaluate(conn, next_working, x, fe, collect_fe_moments = (halvings == 0L && collect_fe_moments))
        if (next_eval$score_max<=max(2*current$score_max,control$score_tol) && next_eval$bad==0 && !dbglm_boundary(next_eval,family$link) && is.finite(next_eval$deviance) && next_eval$deviance <= current$deviance + 64*.Machine$double.eps*max(1,abs(current$deviance))) {
          if (halvings > 0L && collect_fe_moments && is.null(next_eval$fe_moments)) {
            next_eval = dbglm_fe_evaluate(conn, next_working, x, fe, collect_fe_moments = TRUE)
          }
          accepted = TRUE
          break
        }
        drop(next_working);drop(next_state);drop(next_lookup)
        step = step/2
      }
      drop(projected$table);for(st in projected$sub_tables) drop(st)
      if (!accepted) {
        if(isTRUE(verbose)) message("[dbglm] line search exhausted: deviance change=",dbglm_number(next_eval$deviance-current$deviance),
          "; allowed=",dbglm_number(64*.Machine$double.eps*max(1,abs(current$deviance))),
          "; current score=",dbglm_number(current$score_max),"; trial score=",dbglm_number(next_eval$score_max),
          "; mu range=",dbglm_number(next_eval$mu_min),",",dbglm_number(next_eval$mu_max))
        stop("GLM step halving failed to find a valid improving step.",call.=FALSE)
      }
      relative = abs(current$deviance-next_eval$deviance)/(0.1+abs(next_eval$deviance))
      beta = next_beta; gamma = next_gamma; c_eta = next_c_eta
      drop(working);drop(state);drop(cur_lookup)
      state = next_state; working = next_working; cur_lookup = next_lookup; current = next_eval
    } else {
      if (!is.null(projected$raw_projected)) {
        residual = paste(c("z_tilde", sprintf("-(%s)*%s_tilde",dbglm_number(candidate),x)),collapse=" + ")
        tilde_cols = paste(c("__w", paste0(c(x, "z"), "_tilde")), collapse=",")
        candidate_state = create(paste0("SELECT * EXCLUDE(",tilde_cols,"), (z - (",residual,") + __dbglm_o) AS __dbglm_candidate FROM ",projected$raw_projected),materialize=FALSE)
      } else {
        residual = paste(c("p.z", sprintf("-(%s)*p.%s",dbglm_number(candidate),x)),collapse=" + ")
        candidate_state = create(paste0("SELECT b.*, (b.z - (",residual,") + b.__dbglm_o) AS __dbglm_candidate FROM ",working," b JOIN ",projected$table," p USING (__dbglm_key)"),materialize=FALSE)
      }
      step = 1
      accepted = FALSE
      for (halvings in 0:control$step_max_halvings) {
        next_state = create(paste0("SELECT * EXCLUDE(eta,mu,d,w,z,score,dev,__dbglm_eta,__dbglm_candidate), ",
          "__dbglm_eta + ",dbglm_number(step),"*(__dbglm_candidate-__dbglm_eta) AS __dbglm_eta FROM ",candidate_state),materialize=FALSE)
        next_working = dbglm_fe_working(conn,next_state,family$link,create,target=alt_working)
        next_eval = dbglm_fe_evaluate(conn,next_working,x,fe)
        # Reject Fisher steps that collapse information before the next WLS.
        if (next_eval$score_max<=max(2*current$score_max,control$score_tol) && next_eval$bad==0 && !dbglm_boundary(next_eval,family$link) && is.finite(next_eval$deviance) && next_eval$deviance <= current$deviance + 64*.Machine$double.eps*max(1,abs(current$deviance))) {
          accepted = TRUE
          break
        }
        if (is.null(alt_working)) drop(next_working)
        drop(next_state)
        step = step/2
      }
      # working is materialized before releasing its view dependencies.
      if(accepted) {drop(next_state);next_state=NULL}
      drop(candidate_state);drop(projected$table);for(st in projected$sub_tables) drop(st)
      if (!accepted) {
        if(isTRUE(verbose)) message("[dbglm] line search exhausted: deviance change=",dbglm_number(next_eval$deviance-current$deviance),
          "; allowed=",dbglm_number(64*.Machine$double.eps*max(1,abs(current$deviance))),
          "; current score=",dbglm_number(current$score_max),"; trial score=",dbglm_number(next_eval$score_max),
          "; mu range=",dbglm_number(next_eval$mu_min),",",dbglm_number(next_eval$mu_max))
        stop("GLM step halving failed to find a valid improving step.",call.=FALSE)
      }
      relative = abs(current$deviance-next_eval$deviance)/(0.1+abs(next_eval$deviance))
      beta = beta + step*(candidate-beta)
      if (isTRUE(detect_backend(conn)$is_databricks)) {
        alt_working = working
      } else {
        drop(working)
      }
      drop(state)
      state = next_state; working = next_working; current = next_eval
    }
    history = c(history,current$deviance)
    total_halvings = total_halvings+halvings
    if (isTRUE(verbose)) message("[dbglm] iteration ",iter,": deviance=",signif(current$deviance,9),", score=",signif(current$score_max,3))
    if (step==1 && relative<=control$glm_tol && current$score_max<=control$score_tol) {
      if (dbglm_boundary(current,family$link)) stop("GLM reached the numerical mean boundary.",call.=FALSE)
      final = solve_partial_wls(conn,working,x,fe,control,create,drop,response=FALSE,cache=cache,
        materialize=FALSE,scales=if(p) current$scales[seq_len(p)] else numeric(),
        extra_sums=if(p && !is.null(current$raw_score)) "score" else character(),
        need_table=!use_moments || (p>0L && vcov_type!="iid"),precomputed=current$fe_moments)
      H = if (p) {
        if (!is.null(final$cross_reduction)) current$raw_moments[seq_len(p), seq_len(p), drop=FALSE] - final$cross_reduction else
          dbglm_read_information(dbglm_get_query(conn,paste("SELECT",paste(dbglm_moment_terms(x,"w"),collapse=","),"FROM",final$table)),p)
      } else matrix(numeric(),0,0)
      total_cycles = total_cycles+final$cycles
      profiled_score = 0
      if(p) {
        if(any(!is.finite(diag(H))) || any(diag(H)<=0)) stop("GLM slope information lost identification.",call.=FALSE)
        if (!is.null(final$extra_reduction)) {
          s_tilde = current$raw_score - final$extra_reduction[, 1L]
        } else {
          score_row = dbglm_get_query(conn,paste0("SELECT ",paste(sprintf("%s AS s%d",vapply(sprintf("p.%s*b.score",x),function(e) dbglm_exact_sum(conn,e),character(1)),seq_along(x)),collapse=","),
            " FROM ",working," b JOIN ",final$table," p USING (__dbglm_key)"))
          s_tilde = as.numeric(score_row[1,])
        }
        profiled_score = max(abs(s_tilde)/sqrt(diag(H)))
      }
      if(profiled_score>control$score_tol) {drop(final$table);for(st in final$sub_tables) drop(st);next}
      current$score_max = max(current$score_max,profiled_score)
      current$score_profile = profiled_score
      return(c(current,list(beta=beta,H=H,working=working,projected=final$table,raw_projected=final$raw_projected,iterations=iter,deviance_history=history,
        step_halvings=total_halvings,ap_cycles=total_cycles,ap_residual=final$residual,projection_method=if(is.null(final$method)) "cg" else final$method)))
    }
  }
  stop("GLM did not converge within ",control$glm_max_iter," iterations.",call.=FALSE)
}

dbglm_fe_covariance = function(conn,x,fit,kind,ssc,n,counts) {
  p = length(x)
  K = if (kind=="cluster") counts$covariance else counts$covariance
  if (ssc$K.adj && n<=K) stop("Insufficient residual degrees of freedom.",call.=FALSE)
  inverse = if (p) dbglm_solve(fit$H)$inverse else matrix(numeric(),0,0)
  G = NULL
  if (kind=="cluster" && !p) {
    G = as.numeric(dbglm_get_query(conn,paste("SELECT COUNT(DISTINCT __dbglm_cluster) AS n FROM",fit$working))$n)
    if (G<2) stop("Clustered covariance requires at least two clusters.",call.=FALSE)
  }
  if (!p) return(list(matrix=inverse,clusters=G,K=K))
  if (kind=="iid") V = inverse else {
    if (!is.null(fit$raw_projected)) {
      from = paste0(" FROM ", fit$raw_projected, " b")
      x_cols = paste0("b.", x, "_tilde")
    } else {
      from = paste0(" FROM ",fit$working," b JOIN ",fit$projected," p USING (__dbglm_key)")
      x_cols = paste0("p.", x)
    }
    if (kind=="hc1") {
      sql = paste0("SELECT ",paste(dbglm_moment_terms(x_cols,"b.score*b.score"),collapse=","),from)
      row = dbglm_get_query(conn,sql)
    } else {
      scores = sprintf("SUM(%s*b.score) AS s%d",x_cols,seq_along(x))
      sql = paste0("SELECT COUNT(*) AS __n_clusters,",paste(dbglm_moment_terms(paste0("s",seq_along(x)),"1.0"),collapse=",")," FROM (SELECT ",
        paste(scores,collapse=","),from," GROUP BY __dbglm_cluster) t")
      row = dbglm_get_query(conn,sql)
      G = as.numeric(row[["__n_clusters"]])
      if (G<2) stop("Clustered covariance requires at least two clusters.",call.=FALSE)
    }
    meat = dbglm_read_information(row,p)
    V = inverse %*% meat %*% inverse
  }
  if (ssc$K.adj) V = V * (if(kind=="hc1") n else n-1)/(n-K)
  if (kind=="cluster" && ssc$G.adj) V = V*G/(G-1)
  list(matrix=V,clusters=G,K=K)
}

# Preconditioned CG on FE coefficients. Every row correction is reconstructed
# from group lookups, so roundoff cannot introduce arbitrary row-space effects.
# State size is O(number of FE coefficients * number of projected columns).
dbglm_project_cg = function(conn, table, columns, identifiers, ops, scales, target, control, create, drop) {
  index = seq_along(columns)
  r = paste0("__cg_r",index); d = paste0("__cg_d",index)
  g = paste0("__cg_g",index); a = paste0("__cg_a",index)
  direction_names = paste0("__cg_value",index)
  key = "__cg_op,__cg_group"
  moments = function(source, inputs, outputs, info = NULL) {
    pieces = vapply(seq_along(ops),function(i) {
      op = ops[[i]]
      terms = c(if(is.null(info)) sprintf("SQRT(FSUM(w*(%1$s)*(%1$s))) AS __cg_scale",op$v),
        sprintf("FSUM(w*(%s)*(%s)) AS u%d",op$v,inputs,index))
      paste0("SELECT ",i," AS __cg_op,",op$f," AS __cg_group,",paste(terms,collapse=","),
        " FROM ",source," GROUP BY ",op$f)
    },character(1))
    denominator = if(is.null(info)) "u.__cg_scale" else "i.__cg_scale"
    fields = sprintf("COALESCE(u.u%d/NULLIF(%s,0),0) AS %s",index,denominator,outputs)
    create(paste0("SELECT u.__cg_op,u.__cg_group,",denominator," AS __cg_scale,",paste(fields,collapse=","),
      " FROM (",paste(pieces,collapse=" UNION ALL "),") u",
      if(!is.null(info)) paste0(" JOIN ",info," i USING (",key,")")))
  }
  rows = function(state, coefficients, outputs, subtract = FALSE) {
    ctes = vapply(seq_along(ops),function(i) paste0("cg",i," AS (SELECT __cg_group,__cg_scale,",
      paste(coefficients,collapse=",")," FROM ",state," WHERE __cg_op=",i,")"),character(1))
    fields = vapply(index,function(j) {
      terms = vapply(seq_along(ops),function(i) {
        v = if(ops[[i]]$v=="1.0") "1.0" else paste0("s.",ops[[i]]$v)
        sprintf("(%s)*COALESCE(cg%d.%s/NULLIF(cg%d.__cg_scale,0),0)",v,i,coefficients[j],i)
      },character(1))
      value = paste(terms,collapse=" + ")
      if(subtract) value = paste0("s.",columns[j],"-(",value,")")
      paste0("(",value,") AS ",outputs[j])
    },character(1))
    joins = vapply(seq_along(ops),function(i) paste0(" LEFT JOIN cg",i," ON s.",ops[[i]]$f,"=cg",i,".__cg_group"),character(1))
    create(paste0("WITH ",paste(ctes,collapse=",")," SELECT ",
      paste(paste0("s.",c("__dbglm_key",identifiers,"w")),collapse=","),",",paste(fields,collapse=","),
      " FROM ",table," s",paste(joins,collapse="")))
  }
  info = moments(table,columns,r)
  state = create(paste("SELECT *,",paste(c(paste(r,"AS",d),paste("0.0 AS",g)),collapse=","),"FROM",info))
  active = rep(TRUE,length(columns))
  floor_rr = (32*.Machine$double.eps*scales)^2
  for(cycle in seq_len(control$fe_max_iter)) {
    direction = rows(state,d,direction_names)
    applied = moments(direction,direction_names,a,info)
    drop(direction)
    joined = create(paste0("SELECT s.*,",paste(paste0("a.",a),collapse=",")," FROM ",state," s JOIN ",applied," a USING (",key,")"))
    moment = dbglm_get_query(conn,paste("SELECT",paste(c(sprintf("FSUM(%1$s*%1$s) AS rr%2$d",r,index),
      sprintf("FSUM(%s*%s) AS da%d",d,a,index)),collapse=","),"FROM",joined))
    rr = vapply(index,function(j) as.numeric(moment[[paste0("rr",j)]]),numeric(1))
    da = vapply(index,function(j) as.numeric(moment[[paste0("da",j)]]),numeric(1))
    if(any(!is.finite(rr)) || any(!is.finite(da)) || any(active & da<=0 & rr>pmax((target*scales)^2,floor_rr))) {
      stop("GLM fixed-effect projection lost numerical information.",call.=FALSE)
    }
    alpha = ifelse(active & da>0 & rr>floor_rr,rr/da,0)
    updated = create(paste("SELECT",paste(c(key,"__cg_scale",
      sprintf("(%s+(%s)*%s) AS %s",g,dbglm_number(alpha),d,g),
      sprintf("(%s-(%s)*%s) AS %s",r,dbglm_number(alpha),a,r),d),collapse=","),"FROM",joined))
    drop(joined);drop(applied);drop(state)
    state = updated
    restart = cycle %% 20L==0L || all(alpha==0)
    if(cycle %% 4L==0L || cycle==1L || cycle==control$fe_max_iter || restart) {
      projected = rows(state,g,columns,subtract=TRUE)
      certificate = moments(projected,columns,r,info)
      row = dbglm_get_query(conn,paste("SELECT",paste(sprintf("MAX(ABS(%s))/%s AS r%d",r,dbglm_number(scales),index),collapse=","),"FROM",certificate))
      residual_columns = as.numeric(row[1,])
      residual = max(residual_columns)
      active = residual_columns>target
      if(isTRUE(getOption("dbreg.glm_trace")) && cycle %% 100L==0L) message("[projection] cycle ",cycle," residual ",residual," target ",target)
      if(residual<=target) {
        drop(certificate);drop(state);drop(info);drop(table)
        return(list(table=projected,cycles=cycle,residual=residual,method="cg"))
      }
      drop(projected)
      if(restart) {
        restarted = create(paste0("SELECT s.__cg_op,s.__cg_group,s.__cg_scale,",paste(paste0("s.",g),collapse=","),",",
          paste(c(paste0("c.",r," AS ",r),paste0("c.",r," AS ",d)),collapse=","),
          " FROM ",state," s JOIN ",certificate," c USING (",key,")"))
        drop(state);state = restarted
      }
      drop(certificate)
    }
    if(!restart) {
      new_rr = as.numeric(dbglm_get_query(conn,paste("SELECT",paste(sprintf("FSUM(%1$s*%1$s) AS r%2$d",r,index),collapse=","),"FROM",state))[1,])
      beta = ifelse(active & rr>floor_rr,new_rr/rr,0)
      next_state = create(paste("SELECT * EXCLUDE(",paste(d,collapse=","),"),",
        paste(sprintf("((%s)*%s+(%s)*%s) AS %s",dbglm_number(as.numeric(active)),r,dbglm_number(beta),d,d),collapse=","),"FROM",state))
      drop(state);state = next_state
    }
  }
  stop("GLM fixed-effect projection did not converge.",call.=FALSE)
}
