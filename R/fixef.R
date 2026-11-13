# Keys are concatenated UTF-8 components prefixed with their byte lengths.
# Decoding preserves large integer and adjacent floating-point identities.
decode_learned_key = function(key, width) {
  bytes = charToRaw(key)
  colon = as.raw(58)
  offset = 1L
  values = character(width)
  for (i in seq_len(width)) {
    delimiter = match(colon, bytes[seq.int(offset, length(bytes))]) + offset - 1L
    count = as.integer(rawToChar(bytes[seq.int(offset, delimiter - 1L)]))
    offset = delimiter + 1L
    if (count) {
      values[i] = rawToChar(bytes[seq.int(offset, length.out = count)])
      Encoding(values[i]) = "UTF-8"
    }
    offset = offset + count
  }
  values
}

learned_effect_labels = function(dictionary) {
  width = length(dictionary[["columns"]])
  vapply(dictionary[["keys"]], function(key) {
    values = decode_learned_key(key, width)
    if (width == 1L) values else paste0("(", paste(encodeString(values, quote = '"'), collapse = ", "), ")")
  }, character(1), USE.NAMES = FALSE)
}

#' Extract learned fixed effects
#'
#' Returns the saved fixed-effect intercepts and varying slopes on the
#' original trend scale. The identifying normalization assigns the model
#' constant to the first fixed-effect family.
#' @param object A database OLS or GLM fit with learned effects.
#' @param ... Unused.
#' @return A list of named coefficient vectors. Composite groups use quoted
#'   tuples of the native key components. Varying slopes use names of the form
#'   `group[[trend]]`.
#' @exportS3Method fixest::fixef
fixef.dbreg = function(object, ...) {
  if (length(list(...))) stop("Unused arguments to fixef().", call. = FALSE)
  if (!length(object[["fe"]])) return(list())
  state = object[["fitted_state"]]
  if (is.null(state)) {
    stop("Fixed-effect extraction requires a portable learned-effect state from retain = 'prediction' or 'all'.",
      call. = FALSE)
  }
  if (is_database_effect_state(state)) return(database_fixed_effects(object, state))
  dictionaries = state[["dictionaries"]]
  group_names = vapply(dictionaries, function(dict) paste(dict[["columns"]], collapse = "^"), character(1))
  labels = lapply(dictionaries, learned_effect_labels)
  ops = learned_effect_operations(state)
  intercepts = seq_along(dictionaries)
  output = vector("list", length(ops))
  names(output) = vapply(ops, function(op) {
    group = group_names[[op[["group"]]]]
    if (is.null(op[["trend"]])) group else paste0(group, "[[", op[["trend"]], "]]")
  }, character(1))
  for (i in seq_along(ops)) {
    op = ops[[i]]
    group = op[["group"]]
    indices = op[["index"]]
    values = state[["alpha"]][indices]
    values[is.na(indices)] = 0
    group_labels = labels[[group]]
    if (is.null(op[["trend"]])) {
      output[[i]] = stats::setNames(values, group_labels)
    } else {
      slope = values / op[["scale"]]
      output[[i]] = stats::setNames(slope, group_labels)
      intercept = intercepts[group]
      output[[intercept]] = output[[intercept]] - op[["center"]] * slope
    }
  }
  constant = state[["constant"]]
  if (is.null(constant)) constant = 0
  output[[1L]] = output[[1L]] + constant - sum(state[["beta"]] * state[["center"]])
  output
}
