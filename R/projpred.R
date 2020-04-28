#' @export
get_refmodel.brmsfit <- function(object, newdata = NULL, resp = NULL, 
                                 folds = NULL, ...) {
  resp <- validate_resp(resp, object, multiple = FALSE)
  formula <- formula(object)
  if (!is.null(resp)) {
    formula <- formula$forms[[resp]]
  }
  # TODO: check if 'formula' is valid
  bterms <- parse_bf(formula)
  # only use the raw formula for selection of terms
  formula <- formula$formula

  family <- family(object, resp = resp)
  if (family$family == "bernoulli") {
    family$family <- "binomial"
  }
  family <- get(family$family, mode = "function")(link = family$link)
  family <- projpred::extend_family(family)
  
  if (is.null(newdata)) {
    newdata <- model.frame(object)
    y <- get_y(object)
  } else {
    newdata <- validate_newdata(
      newdata, object, resp = resp, check_response = TRUE
    )
    y <- get_y(object, newdata = newdata)
  }
  
  # allows to handle additional arguments implicitely
  extract_aux_data <- function(object, newdata = NULL, ...) {
    .extract_aux_data(object, newdata = newdata, resp = resp, ...)
  }
  
  # using default prediction functions from projpred is fine
  args <- nlist(
    object, data = newdata, y, formula, family, folds, weights,
    ref_predfun = NULL, proj_predfun = NULL, div_minimizer = NULL, 
    extract_aux_data = extract_aux_data, ...
  )
  do_call(projpred::init_refmodel, args)
}

# auxiliary data required in predictions via projpred
# @return a named list with slots 'weights' and 'offset'
.extract_aux_data <- function(object, newdata = NULL, resp = NULL, ...) {
  stopifnot(is.brmsfit(object))
  resp <- validate_resp(resp, object, multiple = FALSE)
  
  # call standata to ensure the correct format of the data
  args <- nlist(
    object, newdata, resp, re_formula = NA,
    check_response = TRUE, internal = TRUE
  )
  sdata <- do_call(standata, args)
  
  # extract relevant auxiliary data
  usc_resp <- usc(resp)
  weights <- sdata[[paste0("weights", usc_resp)]]
  trials <- sdata[[paste0("trials", usc_resp)]]
  offset <- sdata[[paste0("offset", usc_resp)]]
  if (!is.null(trials)) {
    if (!is.null(weights)) {
      stop2("Cannot handle 'trials' and 'weights' at the same time in projpred.") 
    }
    weights <- trials
  }
  nlist(weights, offset)
}
