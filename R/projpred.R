#' Projection Predictive Variable Selection
#' 
#' Perform projection predictive variable selection with the \pkg{projpred}
#' package. See \code{\link[projpred:varsel]{varsel}} and
#' \code{\link[projpred:cv_varsel]{cv_varsel}} for more details.
#' 
#' @aliases varsel cv_varsel
#' 
#' @param object A \code{brmsfit} object.
#' @param ... Further arguments passed to \code{\link{get_refmodel.brmsfit}}
#' as well as \code{\link[projpred:varsel]{varsel.refmodel}} or 
#' \code{\link[projpred:cv_varsel]{cv_varsel.refmodel}}.
#' 
#' @return A \code{vsel} object for which several methods are available
#' in the \pkg{projpred} package.
#' 
#' @examples 
#' \dontrun{
#' # fit a simple model
#' fit <- brm(count ~ zAge + zBase * Trt,
#'            data = epilepsy, family = poisson())
#' summary(fit)
#' 
#' # perform variable selection without cross-validation
#' vs <- varsel(fit)
#' summary(vs)
#' plot(vs)
#' 
#' # perform variable selection with cross-validation
#' cv_vs <- cv_varsel(fit)
#' summary(cv_vs)
#' plot(cv_vs)
#' }
#' 
#' @importFrom projpred varsel
#' @export varsel
#' @export
varsel.brmsfit <- function(object, ...) {
  refmodel <- get_refmodel(object, ...)
  varsel(refmodel, ...)
}

#' @rdname varsel.brmsfit
#' @importFrom projpred cv_varsel
#' @export cv_varsel
#' @export
cv_varsel.brmsfit <- function(object, ...) {
  refmodel <- get_refmodel(object, ...)
  cv_varsel(refmodel, ...)
}

#' Get Reference Models
#' 
#' Get reference model structure from \code{brmsfit} objects for use in
#' \code{\link[projpred:varsel]{varsel}} and related variable selection methods.
#' This method is called automatically when performing variable selection via
#' \code{\link{varsel.brmsfit}} and so you will rarely need to call it manually
#' yourself.
#' 
#' @aliases get_refmodel
#' 
#' @inheritParams posterior_predict.brmsfit
#' @param folds Only used for k-fold variable selection. A vector of fold
#' indices for each data point in data.
#' @param cvfun Optional cross-validation function
#' (see \code{\link[projpred:get-refmodel]{get_refmodel}} for details).
#' If \code{NULL} (the default), \code{cvfun} is defined internally
#' based on \code{\link{kfold.brmsfit}}.
#' @param ... Further arguments passed to 
#' \code{\link[projpred:get-refmodel]{init_refmodel}}.
#' 
#' @return A \code{refmodel} object to be used in
#'   \code{\link[projpred:varsel]{varsel}} and related variable selection
#'   methods.
#' 
#' @importFrom projpred get_refmodel
#' @export get_refmodel
#' @export
get_refmodel.brmsfit <- function(object, newdata = NULL, resp = NULL, 
                                 folds = NULL, cvfun = NULL, ...) {
  resp <- validate_resp(resp, object, multiple = FALSE)
  formula <- formula(object)
  if (!is.null(resp)) {
    formula <- formula$forms[[resp]]
  }
  
  # prepare the family object for use in projpred
  family <- family(object, resp = resp)
  if (family$family == "bernoulli") {
    family$family <- "binomial"
  }
  dot_args <- list(...)
  # For the augmented-data approach, do not re-define ordinal or categorical
  # families to preserve their family-specific extra arguments ("extra" meaning
  # "additionally to `link`") like `refcat` and `thresholds` (see ?brmsfamily):
  if (!(isTRUE(dot_args$aug_data) &&
        (is_ordinal(family$family) || is_categorical(family$family)))) {
    family <- get(family$family, mode = "function")(link = family$link)
  }
  family <- projpred::extend_family(family)
  
  # check if the model is supported by projpred
  bterms <- brmsterms(formula)
  if (length(bterms$dpars) > 1L && !is_categorical(family$family)) {
    stop2("Projpred does not support distributional models.")
  }
  if (length(bterms$nlpars) > 0L) {
    stop2("Projpred does not support non-linear models.")
  }
  not_ok_term_types <- setdiff(all_term_types(), c("fe", "re", "offset", "sm"))
  if (any(not_ok_term_types %in% names(bterms$dpars$mu))) {
    stop2("Projpred only supports standard multilevel terms and offsets.")
  }
  
  # only use the raw formula for selection of terms
  formula <- formula$formula
  # LHS should only contain the response variable
  formula[[2]] <- bterms$respform[[2]]
  
  # projpred requires the dispersion parameter if present
  dis <- NULL
  if (family$family == "gaussian") {
    dis <- paste0("sigma", usc(resp))
    dis <- as.data.frame(object, pars = dis, fixed = TRUE)[[dis]]
  }
  
  ref_predfun <- if (isTRUE(dot_args$aug_data) && is_ordinal(family$family)) {
    # Use argument `incl_thres` of posterior_linpred() (and convert the
    # 3-dimensional array to an "augmented-rows" matrix):
    function(fit, newdata = NULL) {
      # Note: `transform = FALSE` is not needed, but included here for
      # consistency with projpred's default ref_predfun():
      linpred_out <- posterior_linpred(
        fit, transform = FALSE, newdata = newdata, incl_thres = TRUE
      )
      stopifnot(length(dim(linpred_out)) == 3)
      linpred_out <- projpred:::arr2augmat(linpred_out, margin_draws = 1)
      return(linpred_out)
    }
  } else {
    # Using the default prediction function from projpred is fine:
    NULL
  }
  
  # prepare data passed to projpred
  data <- current_data(object, newdata, resp = resp, check_response = TRUE)
  attr(data, "terms") <- NULL
  
  # allows to handle additional arguments implicitly
  extract_model_data <- function(object, newdata = NULL, ...) {
    .extract_model_data(object, newdata = newdata, resp = resp, ...)
  }
  
  # extract a list of K-fold sub-models
  if (is.null(cvfun)) {
    cvfun <- function(folds, ...) {
      cvres <- kfold(
        object, K = max(folds),
        save_fits = TRUE, folds = folds,
        ...
      )
      fits <- cvres$fits[, "fit"]
      return(fits)
    }
  } else {
    if (!is.function(cvfun)) {
      stop2("'cvfun' should be a function.")
    }
  }
  
  args <- nlist(
    object, data, formula, family, folds, dis,
    ref_predfun = ref_predfun, proj_predfun = NULL, div_minimizer = NULL, 
    cvfun = cvfun, extract_model_data = extract_model_data, ...
  )
  do_call(projpred::init_refmodel, args)
}

# auxiliary data required in predictions via projpred
# @return a named list with slots 'weights' and 'offset'
.extract_model_data <- function(object, newdata = NULL, resp = NULL, ...) {
  stopifnot(is.brmsfit(object))
  resp <- validate_resp(resp, object, multiple = FALSE)
  family <- family(object, resp = resp)
  
  # call standata to ensure the correct format of the data
  args <- nlist(
    object, newdata, resp,
    check_response = TRUE, 
    internal = TRUE
  )
  sdata <- do_call(standata, args)
  
  # extract relevant auxiliary data
  usc_resp <- usc(resp)
  y <- as.vector(sdata[[paste0("Y", usc_resp)]])
  offset <- as.vector(sdata[[paste0("offset", usc_resp)]])
  weights <- as.vector(sdata[[paste0("weights", usc_resp)]])
  trials <- as.vector(sdata[[paste0("trials", usc_resp)]])
  if (is_binary(family)) {
    trials <- rep(1, length(y))
  }
  if (!is.null(trials)) {
    if (!is.null(weights)) {
      stop2("Projpred cannot handle 'trials' and 'weights' at the same time.") 
    }
    weights <- trials
  }
  nlist(y, weights, offset)
}
