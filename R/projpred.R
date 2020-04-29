#' @importFrom projpred get_refmodel
#' @export get_refmodel
#' @export
get_refmodel.brmsfit <- function(object, newdata = NULL, resp = NULL, 
                                 folds = NULL, ...) {
  resp <- validate_resp(resp, object, multiple = FALSE)
  formula <- formula(object)
  if (!is.null(resp)) {
    formula <- formula$forms[[resp]]
  }
  
  # check if the model is supported by projpred
  bterms <- parse_bf(formula)
  if (length(bterms$dpars) > 1L) {
    stop2("Projpred does not support distributional models.")
  }
  if (length(bterms$nlpars) > 0L) {
    stop2("Projpred does not support non-linear models.")
  }
  not_ok_term_types <- setdiff(all_term_types(), c("fe", "re", "offset"))
  if (any(not_ok_term_types %in% names(bterms$dpars$mu))) {
    stop2("Projpred only supports standard multilevel terms and offsets.")
  }
  
  # only use the raw formula for selection of terms
  formula <- formula$formula
  # LHS should only contain the response variable
  formula[[2]] <- bterms$respform[[2]]

  # prepare the family object for use in projpred
  family <- family(object, resp = resp)
  if (family$family == "bernoulli") {
    family$family <- "binomial"
  }
  family <- get(family$family, mode = "function")(link = family$link)
  family <- projpred::extend_family(family)
  
  # projpred requires the dispersion parameter if present
  dis <- NULL
  if (family$family == "gaussian") {
    dis <- paste0("sigma", usc(resp))
    dis <- as.data.frame(object, pars = dis, fixed = TRUE)[[dis]]
  }
  
  # allows to handle additional arguments implicitely
  extract_model_data <- function(object, newdata = NULL, ...) {
    .extract_model_data(object, newdata = newdata, resp = resp, ...)
  }
  
  # using default prediction functions from projpred is fine
  args <- nlist(
    object, data = newdata, formula, family, folds, dis,
    ref_predfun = NULL, proj_predfun = NULL, div_minimizer = NULL, 
    extract_model_data = extract_model_data, ...
  )
  do_call(projpred::init_refmodel, args)
}

# auxiliary data required in predictions via projpred
# @return a named list with slots 'weights' and 'offset'
.extract_model_data <- function(object, newdata = NULL, resp = NULL, ...) {
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
  y <- sdata[[paste0("Y", usc_resp)]]
  weights <- sdata[[paste0("weights", usc_resp)]]
  trials <- sdata[[paste0("trials", usc_resp)]]
  offset <- sdata[[paste0("offset", usc_resp)]]
  if (!is.null(trials)) {
    if (!is.null(weights)) {
      stop2("Projpred cannot handle 'trials' and 'weights' at the same time.") 
    }
    weights <- trials
  }
  nlist(y, weights, offset)
}
