# Interpret FE syntax directly: R's formula power operator has other semantics.
parse_fe_components = function(expr) {
  flatten = function(node, op) {
    if (is.call(node) && identical(node[[1L]], as.name("("))) return(flatten(node[[2L]], op))
    if (is.call(node) && length(node) == 3L && identical(node[[1L]], as.name(op))) {
      return(c(flatten(node[[2L]], op), flatten(node[[3L]], op)))
    }
    list(node)
  }
  components = lapply(flatten(expr, "+"), function(term) {
    columns = flatten(term, "^")
    if (!length(columns) || !all(vapply(columns, is.symbol, logical(1)))) {
      stop("Unsupported fixed effects: use simple additive variable names or column-name tuples joined by ^.", call. = FALSE)
    }
    unique(vapply(columns, as.character, character(1)))
  })
  components[!duplicated(vapply(components, function(x) paste(nchar(x), x, collapse = "|"), character(1)))]
}

fe_component_labels = function(components) {
  vapply(components, function(columns) paste(vapply(columns, function(x) paste(deparse(as.name(x), backtick = TRUE), collapse = ""), character(1)), collapse = "^"), character(1))
}

replace_fe_formula = function(fml, fe) {
  lhs = stats::formula(fml, lhs = 1, rhs = 0)[[2L]]
  rhs = stats::formula(fml, lhs = 0, rhs = 1)[[2L]]
  fixed = Reduce(function(a, b) call("+", a, b), lapply(fe, as.name))
  Formula::Formula(stats::as.formula(call("~", lhs, call("|", rhs, fixed)), env = environment(fml)))
}

# A native tuple rank has no delimiter, string conversion or floating key cast.
# Freeze it once so later score, dictionary and sample-export reads agree.
prepare_composite_fe = function(inputs, from, components) {
  conn = inputs[["conn"]]
  joint = which(lengths(components) > 1L)
  fe = vapply(components, `[[`, character(1), 1L)
  if (!length(joint)) return(list(from = from, fe = fe, specs = list(), types = list()))
  schema = DBI::dbGetQuery(conn, sql_limit(paste("SELECT *", from), 0, detect_backend(conn)[["name"]]))
  aliases = safe_identifiers(paste0("dbreg_joint_fe_", joint, "__"), taken = names(schema), fold_case = folds_identifier_case(conn))
  fe[joint] = aliases
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  # Deduplicate tuple keys before ranking so Spark/DuckDB never sorts all
  # observation rows on a single partition.
  joins = vapply(seq_along(joint), function(i) {
    cols = q(components[[joint[i]]])
    col_list = paste(cols, collapse = ",")
    on_clause = paste(paste0("s.", cols, "=d", i, ".", cols), collapse = " AND ")
    paste0(" JOIN (SELECT ", col_list, ", ROW_NUMBER() OVER (ORDER BY ", col_list, ") AS ", q(aliases[i]),
      " FROM (SELECT DISTINCT ", col_list, " FROM __dbreg_comp_src) k", i, ") d", i, " ON ", on_clause)
  }, character(1))
  sql = paste0("WITH __dbreg_comp_src AS (SELECT * ", from, ") SELECT s.*,",
    paste(paste0("d", seq_along(joint), ".", q(aliases)), collapse = ","),
    " FROM __dbreg_comp_src s", paste(joins, collapse = ""))
  if (isTRUE(inputs[["sql_only"]])) {
    from = paste0("FROM (", sql, ") dbreg_composite_source")
  } else {
    backend = detect_backend(conn)
    if (!(backend[["name"]] == "duckdb" || isTRUE(backend[["is_databricks"]]))) {
      stop("Composite FE snapshots currently support DuckDB and native Databricks.", call. = FALSE)
    }
    validate_native_temp_query(conn, backend)
    create = fe_projection_tables(conn, inputs[["resources"]])
    from = paste("FROM", create(sql, c(names(schema), aliases)))
  }
  specs = setNames(components[joint], aliases)
  types = if (isTRUE(inputs[["sql_only"]])) list() else lapply(specs, function(columns) category_key_metadata(conn, sub("^FROM\\s+", "", from), columns)[["types"]])
  list(from = from, fe = fe, specs = specs, types = types)
}

composite_key = function(data, columns, types = NULL) {
  if (any(!columns %in% names(data))) stop("Prediction data is missing composite FE components.", call. = FALSE)
  if (is.null(types)) types = vapply(data[columns], function(x) {
    if (inherits(x, "integer64") || is.integer(x)) "integer" else if (is.numeric(x) && !inherits(x, "Date")) "numeric" else if (is.logical(x)) "logical" else if (inherits(x, "Date")) "date" else "character"
  }, character(1))
  if (length(types) != length(columns) || anyNA(types) || any(!types %in% c("integer", "numeric", "float", "logical", "date", "character"))) {
    stop("Composite key types must be complete and supported.", call. = FALSE)
  }
  pieces = Map(category_key_part, data[columns], types)
  out = do.call(paste0, unname(pieces))
  out[Reduce(`|`, lapply(pieces, is.na))] = NA_character_
  out
}

# Raw components are also present in compressed cells/data_only. Validate
# their wire size before either those cells or a dictionary can cross DBI.
validate_composite_key_transport = function(inputs) {
  columns = unique(unlist(inputs[["composite_fes"]], use.names = FALSE))
  if (!length(columns)) return(invisible(NULL))
  conn = inputs[["conn"]]
  quoted = as.character(DBI::dbQuoteIdentifier(conn, columns))
  text = sql_cast_text(conn, quoted)
  if (detect_backend(conn)[["name"]] == "duckdb") text = paste0("encode(", text, ")")
  costs = paste0("OCTET_LENGTH(", text, ")")
  row = DBI::dbGetQuery(conn, paste0("SELECT MAX(GREATEST(", paste(costs, collapse = ","), ")) AS largest ", inputs[["from_statement"]]))
  largest = as.numeric(row[["largest"]])
  if (length(largest) != 1L || !is.finite(largest) || largest >= 65535) {
    stop("Composite FE key strings exceed the ODBC cell byte limit.", call. = FALSE)
  }
  invisible(NULL)
}

# Only compressed FE prediction needs the training tuple -> fitted-rank map.
# Collect bounded distinct levels, never observations, before releasing SQL.
composite_prediction_dictionaries = function(inputs) {
  specs = inputs[["composite_fes"]]
  if (!length(specs)) return(list())
  conn = inputs[["conn"]]
  from = inputs[["from_statement"]]
  limits = projection_limits()
  q = function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
  budget = list2env(list(cells = 0, bytes = 0, rows = 0), parent = emptyenv())
  lapply(setNames(names(specs), names(specs)), function(alias) {
    columns = specs[[alias]]
    types = inputs[["composite_types"]][[alias]]
    if (is.null(types)) types = category_key_metadata(conn, sub("^FROM\\s+", "", from), columns)[["types"]]
    # Text transports exact native integers; numeric components stay doubles.
    values = vapply(seq_along(columns), function(i) if (types[i] %in% c("integer", "character", "date")) sql_cast_text(conn, q(columns[i])) else q(columns[i]), character(1))
    cost = vapply(values, function(value) if (detect_backend(conn)[["name"]] == "duckdb") paste0("OCTET_LENGTH(encode(CAST(", value, " AS VARCHAR)))") else paste0("OCTET_LENGTH(CAST(", value, " AS STRING))"), character(1))
    source = paste0("(SELECT DISTINCT ", paste(c(q(alias), q(columns)), collapse = ","), " ", from, ") levels")
    stats = DBI::dbGetQuery(conn, paste0("SELECT COUNT(*) AS n,SUM(", paste(cost, collapse = "+"), ") AS bytes,MAX(GREATEST(", paste(cost, collapse = ","), ")) AS largest FROM ", source))
    n = as.numeric(stats[["n"]]); size = as.numeric(stats[["bytes"]]); largest = as.numeric(stats[["largest"]])
    budget[["rows"]] = budget[["rows"]] + n
    budget[["cells"]] = budget[["cells"]] + n * (length(columns) + 1)
    budget[["bytes"]] = budget[["bytes"]] + size
    if (!is.finite(n) || n < 1 || budget[["rows"]] > limits[["summary_rows"]] || budget[["cells"]] > limits[["transfer_cells"]] || !is.finite(size) || 8 * budget[["bytes"]] + 64 * budget[["cells"]] > limits[["memory_bytes"]] || !is.finite(largest) || largest >= 65535) {
      stop("Composite FE prediction dictionaries exceed their row/cell/byte budget.", call. = FALSE)
    }
    sql = paste0("SELECT ", paste(c(paste(q(alias), "AS id"), paste(values, "AS", paste0("v", seq_along(columns)))), collapse = ","), " FROM ", source)
    data = DBI::dbGetQuery(conn, sql_limit(sql, n + 1, detect_backend(conn)[["name"]]))
    if (nrow(data) != n) stop("Composite FE dictionary support changed.", call. = FALSE)
    keys = composite_key(data, paste0("v", seq_along(columns)), types)
    if (anyNA(keys) || anyDuplicated(keys) || anyDuplicated(as.numeric(data[["id"]]))) stop("Invalid composite FE dictionary.", call. = FALSE)
    list(columns = unname(restore_names(columns, inputs[["source_map"]])), types = types,
      keys = keys, ids = sprintf("%.0f", as.numeric(data[["id"]])))
  })
}

restore_composite_formula = function(fml, specs) {
  walk = function(node) {
    if (is.symbol(node) && as.character(node) %in% names(specs)) return(Reduce(function(a, b) call("^", a, b), lapply(specs[[as.character(node)]], as.name)))
    if (is.call(node)) for (i in seq_along(node)[-1L]) node[[i]] = walk(node[[i]])
    node
  }
  value = stats::formula(fml)
  value[[2L]] = walk(value[[2L]])
  value[[3L]] = walk(value[[3L]])
  Formula::Formula(value)
}

prediction_fe_groups = function(object, data) {
  components = object[["fe_components"]]
  if (is.null(components)) return(lapply(object[["fe"]], function(name) data[[name]]))
  Map(function(columns, types) if (length(columns) == 1L) data[[columns]] else composite_key(data, columns, types), components, if (is.null(object[["fe_component_types"]])) rep(list(NULL), length(components)) else object[["fe_component_types"]])
}


prepare_composite_prediction = function(object, data) {
  dictionaries = object[["composite_dictionaries"]]
  if (!length(dictionaries)) return(data)
  source_columns = attr(data, "dbreg_source_columns", exact = TRUE)
  data = as.data.frame(data)
  for (source in names(source_columns)) data[[source]] = data[[source_columns[[source]]]]
  for (alias in names(dictionaries)) {
    dictionary = dictionaries[[alias]]
    keys = composite_key(data, dictionary[["columns"]], dictionary[["types"]])
    data[[alias]] = dictionary[["ids"]][match(keys, dictionary[["keys"]])]
  }
  data
}

# Count/moment fields keep their public contract. Conflicting raw component
# columns get distinct output names, with an explicit source-to-column receipt.
restore_composite_cell_names = function(data, inputs) {
  sources = inputs[["composite_cell_sources"]]
  present = intersect(names(sources), names(data))
  if (!length(present)) return(data)
  mapping = character()
  for (alias in present) {
    source = unname(restore_names(sources[[alias]], inputs[["source_map"]]))
    target = source
    taken = names(data)[names(data) != alias]
    if (target %in% taken) target = utils::tail(make.unique(c(taken, paste0(source, "_source"))), 1L)
    names(data)[match(alias, names(data))] = target
    mapping[source] = target
  }
  attr(data, "dbreg_source_columns") = mapping
  data
}

composite_recipe_metadata = function(recipe, specs, map, types = list()) {
  variable = recipe[["variable"]]
  if (length(variable) && variable %in% names(specs)) {
    components = unname(restore_names(specs[[variable]], map))
    recipe[["display_variable"]] = fe_component_labels(list(components))[[1L]]
    if (identical(recipe[["kind"]], "mean")) {
      recipe[["components"]] = components
      recipe[["key_types"]] = types[[variable]]
    }
  }
  for (child in c("left", "right", "source")) {
    if (!is.null(recipe[[child]])) recipe[[child]] = composite_recipe_metadata(recipe[[child]], specs, map, types)
  }
  recipe
}


# Any term built entirely from one tuple's components is exactly absorbed.
# Encode that structural fact before floating normal equations can invent rank.
composite_absorbed_columns = function(columns, specs) {
  variables = function(recipe) {
    unique(c(recipe[["variable"]], recipe[["second"]], unlist(lapply(c("left", "right", "source"), function(child) if (!is.null(recipe[[child]])) variables(recipe[[child]])), use.names = FALSE)))
  }
  which(vapply(seq_len(nrow(columns)), function(i) {
    sources = variables(columns[["recipe"]][[i]])
    columns[["role"]][i] == "slope" && length(sources) > 0 && any(vapply(specs, function(components) all(sources %in% components), logical(1)))
  }, logical(1)))
}
