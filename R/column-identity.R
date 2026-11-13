# SQL targets, coefficient labels and source variables are separate identities.
# Recipes carry their construction across SQL fitting and R prediction.
empty_linear_columns = function() {
  columns = data.frame(
    id = character(),
    term = integer(),
    basis = integer(),
    alias = character(),
    label = character(),
    role = character()
  )
  columns[["recipe"]] = list()
  columns
}

allocate_sql_aliases = function(
  aliases,
  generated,
  reserved_names = character(),
  fold_case = FALSE
) {
  reserved = unique(c(reserved_names, aliases[!generated]))
  positions = which(generated)
  conflicts = tolower(aliases[positions]) %in%
    tolower(reserved) |
    duplicated(tolower(aliases[positions]))
  changed = positions[conflicts]
  aliases[changed] = safe_identifiers(
    aliases[changed],
    taken = c(reserved, aliases[setdiff(seq_along(aliases), changed)]),
    fold_case = fold_case
  )
  aliases
}

sql_design_reserved = function(inputs) {
  unique(c(
    model_vars(inputs[["fml"]]),
    inputs[["yvar"]],
    inputs[["xvars"]],
    inputs[["fe"]],
    inputs[["weights"]],
    inputs[["offset_safe"]],
    inputs[["cluster_var"]]
  ))
}

linear_recipe_label = function(recipe) {
  variable = recipe[["variable"]]
  kind = recipe[["kind"]]
  source = if (is.null(variable)) NULL else coef_var_name(variable)
  if (isTRUE(recipe[["factor"]])) source = paste0("factor(", source, ")")
  if (kind == "power") return(paste0("I(", coef_var_name(variable), "^", recipe[["exponent"]], ")"))
  if (kind %in% c("raw", "numeric")) {
    return(coef_var_name(variable))
  }
  if (kind == "constant") {
    return(
      if (is.null(variable)) {
        "(Intercept)"
      } else {
        paste0(source, "_constant")
      }
    )
  }
  if (kind == "indicator") {
    return(paste0(source, if (is.null(recipe[["label"]])) as.character(recipe[["value"]]) else recipe[["label"]]))
  }
  if (kind == "i") {
    return(paste0(
      variable,
      "::",
      if (is.null(recipe[["label"]])) as.character(recipe[["value"]]) else recipe[["label"]],
      if (!is.null(recipe[["second"]])) paste0(":", recipe[["second"]])
    ))
  }
  if (kind == "product") {
    return(paste0(
      linear_recipe_label(recipe[["left"]]),
      ":",
      linear_recipe_label(recipe[["right"]])
    ))
  }
  if (kind == "mean") {
    return(paste0(linear_recipe_label(recipe[["source"]]), "_bar_", variable))
  }
  stop("Unknown fitted column recipe.", call. = FALSE)
}

restore_recipe_sources = function(recipe, map) {
  for (name in c("variable", "second")) {
    value = recipe[[name]]
    if (length(value) && value %in% names(map)) {
      recipe[[name]] = unname(map[[value]])
    }
  }
  for (name in c("left", "right", "source")) {
    if (!is.null(recipe[[name]])) {
      recipe[[name]] = restore_recipe_sources(recipe[[name]], map)
    }
  }
  recipe
}

register_sql_design = function(inputs, design) {
  register_factor_display(
    inputs,
    design[["factor_levels"]],
    design[["level_aliases"]],
    design[["i_columns"]]
  )
  inputs[["design_columns"]] = design[["columns"]]
  invisible(inputs)
}

linear_fit_layout = function(columns, fitted) {
  if ("(Intercept)" %in% fitted) {
    intercept = data.frame(
      id = "intercept",
      term = 0L,
      basis = 1L,
      alias = "(Intercept)",
      label = "(Intercept)",
      role = "slope"
    )
    intercept[["recipe"]] = list(list(kind = "constant"))
    columns = rbind(intercept, columns)
  }
  positions = match(fitted, columns[["alias"]])
  if (anyNA(positions) || anyDuplicated(columns[["alias"]])) {
    stop("Fitted columns do not match their producer.", call. = FALSE)
  }
  list(version = 1L, columns = columns, fitted_ids = columns[["id"]][positions])
}

# Follow the R matrix producer's term/basis order, including its no-intercept
# factor promotion. Generated SQL columns retain their original recipes.
matrix_column_layout = function(
  X,
  formula,
  data,
  design = NULL,
  slope_terms = integer()
) {
  tt = stats::terms(formula, data = data)
  factors = attr(tt, "factors")
  variable_name = function(variable) as.character(str2lang(variable))
  if (identical(attr(tt, "intercept"), 0L)) {
    promoted = FALSE
    for (term in seq_len(ncol(factors))) {
      for (variable in rownames(factors)) {
        value = data[[variable_name(variable)]]
        categorical = is.factor(value) ||
          is.character(value) ||
          is.logical(value)
        if (
          factors[variable, term] != 0L &&
            categorical &&
            (is.logical(value) ||
              length(if (is.factor(value)) levels(value) else unique(value)) >
                1L)
        ) {
          factors[variable, term] = 2L
          promoted = TRUE
          break
        }
      }
      if (promoted) break
    }
  }
  assignment = attr(X, "assign")
  recipes = rep(list(list(kind = "constant")), ncol(X))
  for (term in seq_along(attr(tt, "term.labels"))) {
    positions = which(assignment == term)
    variables = rownames(factors)[factors[, term] != 0L]
    parts = lapply(variables, function(variable) {
      actual = variable_name(variable)
      emitted = match(actual, design[["alias"]])
      if (!is.na(emitted)) {
        recipe = design[["recipe"]][[emitted]]
        if (!identical(recipe[["kind"]], "raw")) return(list(recipe))
      }
      value = data[[actual]]
      if (is.factor(value) || is.character(value) || is.logical(value)) {
        levels = if (is.logical(value)) {
          c("FALSE", "TRUE")
        } else if (is.factor(value)) {
          levels(value)
        } else {
          levels(factor(value))
        }
        chosen = if (factors[variable, term] == 1L) levels[-1L] else levels
        return(lapply(chosen, function(level) {
          list(
            kind = "indicator",
            variable = actual,
            value = level,
            levels = levels
          )
        }))
      }
      list(list(kind = "numeric", variable = actual))
    })
    combinations = expand.grid(
      lapply(parts, seq_along),
      KEEP.OUT.ATTRS = FALSE
    )
    if (nrow(combinations) != length(positions)) {
      stop("Matrix basis does not match its term assignment.", call. = FALSE)
    }
    recipes[positions] = lapply(seq_len(nrow(combinations)), function(j) {
      operands = lapply(seq_along(parts), function(k) {
        parts[[k]][[combinations[j, k]]]
      })
      Reduce(
        function(left, right) {
          list(kind = "product", left = left, right = right)
        },
        operands
      )
    })
  }
  columns = data.frame(
    id = paste0("m", seq_len(ncol(X))),
    term = assignment,
    basis = seq_len(ncol(X)),
    alias = colnames(X),
    label = colnames(X),
    role = ifelse(
      assignment == 0L | assignment %in% slope_terms,
      "slope",
      "fixed_effect"
    )
  )
  columns[["recipe"]] = recipes
  columns
}

finalize_linear_layout = function(result, inputs) {
  layout = result[["linear_design"]]
  if (is.null(layout)) {
    return(result)
  }
  columns = layout[["columns"]]
  columns[["recipe"]] = lapply(
    columns[["recipe"]],
    restore_recipe_sources,
    map = inputs[["source_map"]]
  )
  columns[["label"]] = vapply(
    columns[["recipe"]],
    linear_recipe_label,
    character(1)
  )
  positions = match(layout[["fitted_ids"]], columns[["id"]])
  if (anyNA(positions) || length(positions) != nrow(result[["coeftable"]])) {
    stop("Invalid fitted column identity.", call. = FALSE)
  }
  labels = columns[["label"]][positions]
  rownames(result[["coeftable"]]) = labels
  dimnames(result[["vcov"]]) = rep(list(labels), 2L)
  visible = columns[["role"]][positions] == "slope"
  if (length(result[["fe"]])) {
    visible = visible & columns[["alias"]][positions] != "(Intercept)"
  }
  result[["coef_names"]] = labels[visible]
  result[["xvars"]] = labels[
    visible & columns[["alias"]][positions] != "(Intercept)"
  ]
  if (length(result[["collin.var"]])) {
    hit = match(result[["collin.var"]], columns[["alias"]])
    result[["collin.var"]][!is.na(hit)] = columns[["label"]][hit[!is.na(hit)]]
  }
  compact = compact_recipe_dictionaries(columns[["recipe"]])
  columns[["recipe"]] = compact[["recipes"]]
  layout[["dictionaries"]] = compact[["dictionaries"]]
  layout[["columns"]] = columns
  result[["linear_design"]] = layout
  means = which(vapply(
    columns[["recipe"]],
    function(recipe) identical(recipe[["kind"]], "mean"),
    logical(1)
  ))
  if (length(means)) {
    result[["mundlak_means"]] = data.frame(
      term = vapply(
        columns[["recipe"]][means],
        function(recipe) linear_recipe_label(recipe[["source"]]),
        character(1)
      ),
      fe = vapply(columns[["recipe"]][means], `[[`, character(1), "variable"),
      name = columns[["label"]][means]
    )
  }
  iv = inputs[["iv_label_columns"]]
  for (name in names(iv)) {
    recipes = lapply(
      iv[[name]][["recipe"]],
      restore_recipe_sources,
      map = inputs[["source_map"]]
    )
    result[[name]] = vapply(recipes, linear_recipe_label, character(1))
    if (name == "endogenous") {
      for (diagnostic in c("first_stage_f", "first_stage_wald")) {
        if (length(result[["diagnostics"]][[diagnostic]])) {
          names(result[["diagnostics"]][[diagnostic]]) = result[[name]]
        }
      }
    }
  }
  validate_linear_layout(result)
  result
}

visible_linear_positions = function(object, fe = FALSE) {
  layout = object[["linear_design"]]
  positions = seq_along(layout[["fitted_ids"]])
  if (!isTRUE(fe) && length(object[["fe"]])) {
    roles = layout[["columns"]][["role"]][match(
      layout[["fitted_ids"]],
      layout[["columns"]][["id"]]
    )]
    aliases = layout[["columns"]][["alias"]][match(
      layout[["fitted_ids"]],
      layout[["columns"]][["id"]]
    )]
    positions = positions[roles == "slope" & aliases != "(Intercept)"]
  }
  positions
}

linear_recipe_value = function(recipe, data, dictionaries = list(), weights = NULL) {
  levels = recipe[["levels"]]
  if (!is.null(recipe[["dictionary"]])) {
    levels = dictionaries[[recipe[["dictionary"]]]]
  }
  kind = recipe[["kind"]]
  variable = recipe[["variable"]]
  column = function(name) {
    if (!name %in% names(data)) {
      stop("Prediction data is missing variable: ", name, call. = FALSE)
    }
    data[[name]]
  }
  numeric_column = function(name) {
    value = column(name)
    if (
      !(is.numeric(value) || is.logical(value)) ||
        !is.null(dim(value)) ||
        inherits(value, c("Date", "POSIXt"))
    ) {
      stop("Prediction variable must be numeric: ", name, call. = FALSE)
    }
    if (
      inherits(value, "integer64") && !requireNamespace("bit64", quietly = TRUE)
    ) {
      stop("integer64 prediction requires bit64.", call. = FALSE)
    }
    as.numeric(value)
  }
  if (kind == "constant") {
    value = rep(1, nrow(data))
    if (!is.null(variable)) {
      original = column(variable)
      unknown = if (is.null(levels)) {
        rep(FALSE, length(original))
      } else {
        if (!is.null(recipe[["key_type"]])) !category_key_part(original, recipe[["key_type"]]) %in% category_key_part(levels, recipe[["key_type"]]) else !as.character(original) %in% as.character(levels)
      }
      value[is.na(original) | unknown] = NA_real_
    }
    return(value)
  }
  if (kind %in% c("numeric", "raw")) {
    return(numeric_column(variable))
  }
  if (kind == "power") return(numeric_column(variable)^recipe[["exponent"]])
  if (kind %in% c("i", "indicator") && !is.null(recipe[["key_type"]])) {
    value = category_key_part(column(variable), recipe[["key_type"]])
    target = category_key_part(recipe[["value"]], recipe[["key_type"]])
    domain = category_key_part(levels, recipe[["key_type"]])
    out = as.numeric(value == target)
    out[is.na(value) | !value %in% domain] = NA_real_
    if (!is.null(recipe[["second"]])) out = out * numeric_column(recipe[["second"]])
    return(out)
  }
  if (kind %in% c("indicator", "i")) {
    value = column(variable)
    target = recipe[["value"]]
    numeric_keys = is.numeric(levels) &&
      !inherits(levels, c("integer64", "Date", "POSIXt"))
    value = if (numeric_keys) {
      if (is.factor(value) || is.character(value)) {
        as.numeric(as.character(value))
      } else {
        as.numeric(value)
      }
    } else {
      as.character(value)
    }
    levels = if (numeric_keys) as.numeric(levels) else as.character(levels)
    target = if (numeric_keys) as.numeric(target) else as.character(target)
    out = as.numeric(value == target)
    out[is.na(value) | !value %in% levels] = NA_real_
    if (kind == "i" && !is.null(recipe[["second"]])) {
      out = out * numeric_column(recipe[["second"]])
    }
    return(out)
  }
  if (kind == "product") {
    return(
      linear_recipe_value(recipe[["left"]], data, dictionaries, weights) *
        linear_recipe_value(recipe[["right"]], data, dictionaries, weights)
    )
  }
  if (kind == "mean") {
    return(weighted_group_mean(
      linear_recipe_value(recipe[["source"]], data, dictionaries, weights),
      column(variable), weights
    ))
  }
  stop("Unknown fitted column recipe.", call. = FALSE)
}

# Category validity belongs to the fitted design, including columns dropped for
# collinearity. Check each saved domain once without constructing indicator cells.
linear_category_invalid = function(layout, data) {
  domains = function(recipe) {
    result = list()
    if (recipe[["kind"]] %in% c("constant", "indicator", "i") && !is.null(recipe[["variable"]])) {
      levels = recipe[["levels"]]
      if (!is.null(recipe[["dictionary"]])) levels = layout[["dictionaries"]][[recipe[["dictionary"]]]]
      if (!is.null(levels)) {
        type = recipe[["key_type"]]
        if (is.null(type)) type = if (inherits(levels, "integer64")) "integer" else
          if (is.numeric(levels) && !inherits(levels, "Date")) "numeric" else "character"
        result = list(list(variable = recipe[["variable"]], type = type, levels = levels))
      }
    }
    for (child in c("left", "right", "source")) {
      if (!is.null(recipe[[child]])) result = c(result, domains(recipe[[child]]))
    }
    result
  }
  records = unique(unlist(lapply(layout[["columns"]][["recipe"]], domains), recursive = FALSE))
  invalid = rep(FALSE, nrow(data))
  for (record in records) {
    variable = record[["variable"]]
    if (!variable %in% names(data)) stop("Prediction data is missing variable: ", variable, call. = FALSE)
    value = data[[variable]]
    # SQLite stores logical and date columns as numbers; such levels cannot
    # match R logical or Date values, which would otherwise miss every level.
    levels = as.character(record[["levels"]])
    if ((is.logical(value) && record[["type"]] != "logical" && !any(levels %in% c("TRUE", "FALSE"))) ||
        (inherits(value, "Date") && record[["type"]] != "date" &&
          !any(grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", levels)))) {
      stop("Prediction column `", variable, "` is ", class(value)[1L],
        ", but the fitted levels were stored as ", record[["type"]],
        " values; convert it to the stored values.", call. = FALSE)
    }
    value = category_key_part(value, record[["type"]])
    keys = category_key_part(record[["levels"]], record[["type"]])
    invalid = invalid | is.na(value) | !value %in% keys
  }
  invalid
}

linear_prediction_matrix = function(object, data) {
  validate_linear_layout(object)
  layout = object[["linear_design"]]
  columns = layout[["columns"]]
  recipes = columns[["recipe"]][match(layout[["fitted_ids"]], columns[["id"]])]
  weights = if (identical(object[["strategy"]], "mundlak")) prediction_weights(data, object[["weights"]]) else NULL
  rows = positions = values = vector("list", length(recipes))
  invalid = linear_category_invalid(layout, data)
  for (k in seq_along(recipes)) {
    value = linear_recipe_value(recipes[[k]], data, layout[["dictionaries"]], weights)
    if (length(value) != nrow(data)) {
      stop("Prediction column does not match its source rows.", call. = FALSE)
    }
    invalid = invalid | is.na(value)
    present = which(!is.na(value) & value != 0)
    rows[[k]] = present
    positions[[k]] = rep.int(k, length(present))
    values[[k]] = value[present]
  }
  matrix = Matrix::sparseMatrix(
    i = as.integer(unlist(rows, use.names = FALSE)),
    j = as.integer(unlist(positions, use.names = FALSE)),
    x = as.numeric(unlist(values, use.names = FALSE)),
    dims = c(nrow(data), length(recipes)),
    dimnames = list(NULL, layout[["fitted_ids"]])
  )
  attr(matrix, "dbreg_invalid_rows") = invalid
  matrix
}

# Validate fitted axes before consuming serialized producer metadata. Public
# labels may repeat; producer IDs and SQL targets may not.
validate_linear_layout = function(object) {
  layout = object[["linear_design"]]
  columns = layout[["columns"]]
  invalid = !identical(layout[["version"]], 1L) ||
    !is.data.frame(columns) ||
    !all(c("id", "alias", "label", "role", "recipe") %in% names(columns)) ||
    !is.character(columns[["id"]]) ||
    anyNA(columns[["id"]]) ||
    any(!nzchar(columns[["id"]])) ||
    anyDuplicated(columns[["id"]]) ||
    !is.character(columns[["alias"]]) ||
    anyNA(columns[["alias"]]) ||
    any(!nzchar(columns[["alias"]])) ||
    anyDuplicated(columns[["alias"]]) ||
    !is.character(columns[["label"]]) ||
    anyNA(columns[["label"]]) ||
    !is.list(layout[["dictionaries"]]) ||
    any(!columns[["role"]] %in% c("slope", "fixed_effect", "nuisance")) ||
    !is.character(layout[["fitted_ids"]]) ||
    anyNA(layout[["fitted_ids"]]) ||
    anyDuplicated(layout[["fitted_ids"]]) ||
    any(!layout[["fitted_ids"]] %in% columns[["id"]]) ||
    !is.list(columns[["recipe"]]) ||
    length(layout[["fitted_ids"]]) != nrow(object[["coeftable"]]) ||
    !identical(
      as.integer(dim(object[["vcov"]])),
      rep.int(as.integer(length(layout[["fitted_ids"]])), 2L)
    )
  if (invalid) {
    stop("Invalid fitted column metadata.", call. = FALSE)
  }
  invisible(TRUE)
}

# Split formula interactions as syntax, preserving quoted source names that
# themselves contain a colon or other operator characters.
split_interaction_term = function(term) {
  split = function(expression) {
    if (is.call(expression) && identical(expression[[1L]], as.name(":"))) {
      return(c(split(expression[[2L]]), split(expression[[3L]])))
    }
    paste(deparse(expression, backtick = TRUE), collapse = "")
  }
  split(str2lang(term))
}

sql_column_name = function(name, conn) {
  safe = is_safe_identifier(name, fold_case = folds_identifier_case(conn))
  out = name
  out[!safe] = as.character(DBI::dbQuoteIdentifier(conn, name[!safe]))
  out
}

# Share immutable category domains once in a serialized fitted model, rather
# than copying a G-level dictionary into each of its G indicator recipes.
compact_recipe_dictionaries = function(recipes) {
  state = new.env(parent = emptyenv())
  state[["dictionaries"]] = list()
  compact = function(recipe) {
    levels = recipe[["levels"]]
    if (!is.null(levels)) {
      matches = vapply(state[["dictionaries"]], identical, logical(1), levels)
      if (any(matches)) {
        index = which(matches)[[1L]]
      } else {
        index = length(state[["dictionaries"]]) + 1L
        state[["dictionaries"]][[index]] = levels
      }
      recipe[["levels"]] = NULL
      recipe[["dictionary"]] = index
    }
    for (name in c("left", "right", "source")) {
      if (!is.null(recipe[[name]])) recipe[[name]] = compact(recipe[[name]])
    }
    recipe
  }
  list(
    recipes = lapply(recipes, compact),
    dictionaries = state[["dictionaries"]]
  )
}

# R matrix constructors otherwise interpret integer64's underlying double bits.
# Category columns must be factored first so their exact decimal labels survive.
linear_matrix_data = function(data) {
  for (name in names(data)) if (inherits(data[[name]], "integer64")) {
    data[[name]] = as.numeric(data[[name]])
  }
  data
}
