#' @export
get_refmodel_poc.brmsfit <- function(fit, newdata = NULL, resp = NULL, 
                                     folds = NULL, ...) {
  resp <- validate_resp(resp, fit, multiple = FALSE)
  formula <- formula(fit)
  if (!is.null(resp)) {
    formula <- formula$forms[[resp]]
  }
  # TODO: check if 'formula' is valid
  bterms <- parse_bf(formula)
  # only use the raw formula for selection of terms
  formula <- formula$formula

  family <- family(fit, resp = resp)
  if (family$family == "bernoulli") {
    family$family <- "binomial"
  }
  family <- get(family$family, mode = "function")(link = family$link)
  family <- projpred::extend_family(family)
  
  if (is.null(newdata)) {
    newdata <- model.frame(fit)
    y <- get_y(fit)
  } else {
    newdata <- validate_newdata(newdata, fit, resp = resp, check_response = TRUE)
    y <- get_y(fit, newdata = newdata)
  }
  
  weights <- NULL
  if (is.formula(bterms$adforms$trials)) {
    if (is.formula(bterms$adforms$weights)) {
      stop2("Cannot handle 'trials' and 'weights' at the same time in projpred.")
    }
    weights <- get_ad_values(bterms, "trials", "trials", newdata)
  } else if (is.formula(bterms$adforms$weights)) {
    weights <- get_ad_values(bterms, "weights", "weights", newdata)
  }
  # TODO: remove as soon as handled in projpred
  if (!is.null(weights)) {
    y <- y / weights
  }
  
  # using default prediction functions from projpred is fine
  args <- nlist(
    fit, data = newdata, y, formula, family, folds, weights,
    predfun = NULL, proj_predfun = NULL, mle = NULL, ...
  )
  do_call(projpred::init_refmodel_poc, args)
}
