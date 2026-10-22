# Independent oracle: refine only from fixest's own estimates. Its deviance
# stopping rule and previous-WLS covariance can stop too early for strict gates.
dbglm_reference_feglm <- function(fml,data,...,max_refine=60L) {
 args <- c(list(fml=fml,data=data),list(...))
 defaults <- list(glm.tol=3e-13,glm.iter=100,fixef.tol=1e-10,fixef.iter=10000,notes=FALSE)
 for(n in names(defaults)) if(is.null(args[[n]])) args[[n]] <- defaults[[n]]
 fit <- do.call(fixest::feglm,args)
 if(!length(stats::coef(fit))) {
  if(!isTRUE(fit$convStatus)) stop('The FE-only fixest reference did not converge.')
  return(fit)
 }
 for(i in seq_len(max_refine)) {
  previous_coef <- stats::coef(fit);previous_vcov <- stats::vcov(fit)
  # A slope-only start makes fixest reinitialize the FE component at loose
  # tolerances. Warm the entire predictor from fixest's OWN previous fit.
  args$start <- NULL
  args$mustart <- NULL
  args$etastart <- fit$linear.predictors
  fit <- do.call(fixest::feglm,args)
  current_coef <- stats::coef(fit);current_vcov <- stats::vcov(fit)
  if(!identical(names(current_coef),names(previous_coef))) stop('The fixest reference basis changed.')
  coef_stable <- all(abs(current_coef-previous_coef)<=1e-10+1e-10*abs(current_coef))
  vcov_stable <- all(abs(current_vcov-previous_vcov)<=1e-10+1e-10*abs(current_vcov))
  if(i>=2L && coef_stable && vcov_stable && isTRUE(fit$convStatus)) {
   attr(fit,'reference_refinements') <- i
   return(fit)
  }
 }
 stop('The independent fixest reference did not stabilize within the refinement cap.')
}
