# Restricted i() declarations compile to ordinary SQL design columns. The same
# saved level dictionary drives R and SQL prediction; newdata discovers no levels.

#' Factor indicators and numeric group slopes
#'
#' Use additive `i(group)` or `i(group, numeric_column)` terms in [feols()],
#' [feglm()] and [dbglm()] formulas. Indicators omit the first fitted level
#' when the model has an intercept or absorbed effects. Group slopes keep all
#' levels unless `ref` or `keep` is supplied. Formula controls must be literals;
#' categorical second variables, binning and interactions around `i()` are
#' unsupported. Fitted dictionaries are bounded by
#' `options(dbreg.i_max_columns = 1000)` and reused for prediction.
#'
#' @param factor_var Categorical vector, or a source column in a formula.
#' @param var Optional numeric/logical vector or column for group slopes.
#' @param ref Levels to omit; `TRUE` omits the first level. `NULL` keeps all.
#' @param keep Levels to retain. Character selectors use case-insensitive
#'   partial matching; a leading `@` selects levels by regular expression.
#' @param ... Unsupported controls produce an error.
#' @return Outside a regression formula, a numeric indicator/group-slope matrix
#'   with all levels by default. Formula expansion happens in SQL.
#' @md
#' @export
i = function(factor_var,var,ref=NULL,keep=NULL,...) {
  check_unused_args(list(...))
  first = factor_var
  if (!is.atomic(first) || is.matrix(first) || !length(first)) stop("factor_var must be a nonempty vector.",call.=FALSE)
  second = if(!missing(var)) var else NULL
  if (!is.null(second) && length(second)!=length(first)) {
    if (missing(ref) && length(second) %in% 1:2) { ref = second; second = NULL } else stop("var must have the same length as factor_var.",call.=FALSE)
  }
  type = if(inherits(first,"integer64") || is.integer(first)) "integer" else if(inherits(first,"Date")) "date" else
    if(is.numeric(first)) "numeric" else if(is.logical(first)) "logical" else "character"
  values = if(is.factor(first)) levels(first)[levels(first) %in% as.character(first)] else sort(unique(first[!is.na(first)]))
  declaration = list(first="group",second=if(!is.null(second)) "var",ref=ref,keep=keep,explicit=TRUE,term="i()")
  spec = make_i_spec(declaration,values,type,if(!is.null(second)) class(second)[[1L]],intercept=FALSE)
  result = evaluate_i_spec(spec,first,second)
  colnames(result) = spec[["labels"]][spec[["selected"]]]
  result
}

i_literal = function(expr) {
  if (is.null(expr)) return(NULL)
  if (is.call(expr)) {
    op = as.character(expr[[1L]])
    args = as.list(expr)[-1L]
    if (op == "c") return(unlist(lapply(args,i_literal),use.names=FALSE))
    if (op %in% c("+","-") && length(args)==1L) {
      value = i_literal(args[[1L]])
      if (is.numeric(value)) return(if(op=="-") -value else value)
    }
  } else if ((is.numeric(expr) || is.character(expr) || is.logical(expr)) && !anyNA(expr) &&
    (!is.numeric(expr) || all(is.finite(expr)))) return(expr)
  stop("i() references/keep values must be finite, nonmissing literals or literal c() vectors.",call.=FALSE)
}

parse_i_term = function(term) {
  expr = if(is.character(term)) tryCatch(str2lang(term),error=function(e) NULL) else term
  if (!is.call(expr)) return(NULL)
  head = expr[[1L]]
  recognized = identical(head,as.name("i")) || is.call(head) && identical(head[[1L]],as.name("::")) &&
    as.character(head[[2L]]) %in% c("dbreg","fixest") && identical(head[[3L]],as.name("i"))
  if (!recognized) return(NULL)
  expr[[1L]] = as.name("i")
  definition = function(factor_var,var,ref=NULL,keep=NULL,...) NULL
  matched = match.call(definition,expr,expand.dots=FALSE)
  args = as.list(matched)[-1L]
  if (length(args[["..."]])) stop("Supported i() controls are ref and keep; categorical second variables and binning are unsupported.",call.=FALSE)
  if (is.null(args[["factor_var"]]) || !is.symbol(args[["factor_var"]])) stop("i() requires a categorical source column.",call.=FALSE)
  var = args[["var"]]
  if (!is.null(var) && !is.symbol(var)) {
    shorthand = i_literal(var)
    if (!"ref" %in% names(args) && length(shorthand) %in% 1:2) {
      matched[["ref"]] = var
      matched[["var"]] = NULL
      args = as.list(matched)[-1L]
      var = NULL
    } else stop("The second i() variable must be a numeric/logical column.",call.=FALSE)
  }
  if (!is.null(var) && startsWith(as.character(var),"i.")) stop("Categorical second-variable i. forcing is unsupported; use ordinary factor interactions.",call.=FALSE)
  list(first=as.character(args[["factor_var"]]),second=if(!is.null(var)) as.character(var),
    ref=if("ref" %in% names(args)) i_literal(args[["ref"]]),
    keep=if("keep" %in% names(args)) i_literal(args[["keep"]]),
    explicit=any(c("ref","keep") %in% names(args)),call=matched,
    term=paste(deparse(matched,width.cutoff=500L),collapse=""))
}

contains_i_term = function(expr) {
  if (!is.call(expr)) return(FALSE)
  !is.null(parse_i_term(expr)) || any(vapply(as.list(expr)[-1L],contains_i_term,logical(1)))
}

i_column_limit = function() {
  value = getOption("dbreg.i_max_columns",1000L)
  if (!is.numeric(value) || length(value)!=1L || !is.finite(value) || value<1 || value!=floor(value)) stop("dbreg.i_max_columns must be a positive integer.",call.=FALSE)
  value
}

i_selected_items = function(value,keys,labels,type) {
  if (is.character(value)) {
    found = integer()
    for (item in value) {
      indices = if(startsWith(item,"@")) grep(substring(item,2L),labels) else pmatch(tolower(item),tolower(labels),duplicates.ok=TRUE)
      if (!length(indices) || anyNA(indices)) stop("An i() reference/keep value does not uniquely match fitted levels: ",item,call.=FALSE)
      found = union(found,indices)
    }
    return(found)
  }
  if (is.logical(value) && type %in% c("integer","numeric")) value = as.numeric(value)
  if (is.numeric(value) && type=="logical") {
    if (any(!value %in% c(0,1))) stop("i() reference/keep values do not match logical levels.",call.=FALSE)
    value = as.logical(value)
  }
  encoded = prediction_key_part(value,type)
  found = match(encoded,prediction_key_part(keys,type))
  if (anyNA(found)) stop("An i() reference/keep value does not match fitted levels.",call.=FALSE)
  unique(found)
}

make_i_spec = function(declaration,keys,type,second_class=NULL,intercept=TRUE) {
  labels = as.character(keys)
  if (type=="numeric" && anyDuplicated(labels)) labels = dbglm_number(keys)
  encoded = prediction_key_part(keys,type)
  if (anyNA(encoded) || anyDuplicated(encoded) || anyDuplicated(labels)) stop("Could not preserve i() levels exactly.",call.=FALSE)
  selected = seq_along(keys)
  if (!is.null(declaration[["second"]]) && !second_class %in% c("numeric","integer","integer64","logical")) stop("The second i() variable must be numeric/logical; categorical interactions use ordinary factor syntax.",call.=FALSE)
  ref = declaration[["ref"]]
  if (isTRUE(ref) || !declaration[["explicit"]] && intercept && is.null(declaration[["second"]])) {
    selected = setdiff(selected,1L)
  } else if (!is.null(ref)) selected = setdiff(selected,i_selected_items(ref,keys,labels,type))
  if (!is.null(declaration[["keep"]])) selected = intersect(selected,i_selected_items(declaration[["keep"]],keys,labels,type))
  names = if(length(selected)) paste0(declaration[["first"]],"::",labels[selected],
    if(!is.null(declaration[["second"]])) paste0(":",declaration[["second"]])) else character()
  list(first=declaration[["first"]],second=declaration[["second"]],second_class=second_class,
    values=keys,keys=encoded,type=type,labels=labels,selected=selected,names=names,term=declaration[["term"]])
}

i_match_values = function(spec,value,strict=TRUE) {
  valid = !is.na(value)
  if (is.numeric(value) && !inherits(value,"integer64")) valid = valid & is.finite(value)
  keys = rep(NA_character_,length(value))
  keys[valid] = prediction_key_part(value[valid],spec[["type"]])
  index = match(keys,spec[["keys"]])
  if (strict && any(valid & is.na(index))) stop("newdata contains unseen i() levels.",call.=FALSE)
  index
}

evaluate_i_spec = function(spec,first,second=NULL) {
  index = i_match_values(spec,first)
  columns = spec[["selected"]]
  if (!length(columns)) return(matrix(numeric(),length(first),0L))
  result = outer(index,columns,`==`)*1.0
  if (!is.null(spec[["second"]])) {
    if (is.factor(second)) second = as.character(second)
    if (spec[["second_class"]]=="logical" && !is.logical(second)) second = as.logical(second)
    result = result*as.numeric(second)
  }
  colnames(result) = spec[["names"]]
  result
}

expand_i_spec = function(spec,info) {
  selected = spec[["selected"]]
  if (!length(selected)) return(list())
  literals = info[["level_sql"]][[spec[["first"]]]][selected]
  expression = paste0("CASE WHEN ",spec[["first"]],"=",literals," THEN 1.0 ELSE 0.0 END")
  if (!is.null(spec[["second"]])) {
    value = if(spec[["second_class"]]=="logical") paste0("CAST(CAST(",spec[["second"]]," AS BOOLEAN) AS DOUBLE)") else paste0("CAST(",spec[["second"]]," AS DOUBLE)")
    expression = paste0("(",expression,")*(",value,")")
  }
  lapply(seq_along(expression),function(j) list(sql=expression[[j]],name=spec[["names"]][[j]],display=spec[["names"]][[j]]))
}

# Only bounded, typed dictionaries cross the SQL/R boundary. Integer levels are
# transferred as decimal strings, never through the driver's double conversion.
glm_i_dictionaries = function(conn,from_sql,declarations,classes,intercept) {
  result = list()
  if (!length(declarations)) return(result)
  limit = i_column_limit()
  source = sub("^FROM\\s+","",from_sql,ignore.case=TRUE)
  for (declaration in declarations) {
    first = declaration[["first"]]
    metadata = prediction_key_metadata(conn,source,first)
    type = metadata[["types"]][[1L]]
    if (type=="float") stop("i() floating-point categories require a DOUBLE source column.",call.=FALSE)
    value = if(type=="integer") paste0("CAST(",first," AS STRING)") else first
    # Group by the original typed value; textual ordering would put 10 before 2.
    keys = dbglm_get_query(conn,paste0("SELECT ",value," AS value ",from_sql," GROUP BY ",first,
      " ORDER BY ",first," LIMIT ",sprintf("%.0f",limit+2L)))[[1L]]
    if (length(keys)>limit+1L) stop("i() level dictionary exceeds dbreg.i_max_columns.",call.=FALSE)
    # Match R/fixest character collation on this bounded dictionary. Preserve
    # native enum order, numeric order and exact decimal integer text.
    if (type=="character" && !startsWith(tolower(metadata[["sql_types"]][[1L]]),"enum")) keys = sort(keys)
    if (type=="integer" && prediction_integer_sql_type(keys,conn)!="BIGINT") stop("i() integer levels must fit signed BIGINT.",call.=FALSE)
    second = declaration[["second"]]
    spec = make_i_spec(declaration,keys,type,if(!is.null(second)) classes[[second]],intercept || length(result)>0L)
    result[[first]] = spec
    if (sum(vapply(result,function(x) length(x[["selected"]]),integer(1)))>limit) stop("i() columns exceed dbreg.i_max_columns.",call.=FALSE)
  }
  result
}
