#' Projection Predictive Variable Selection: Get Reference Model
#' 
#' The \code{get_refmodel.brmsfit} method can be used to create the reference
#' model structure which is needed by the \pkg{projpred} package for performing
#' a projection predictive variable selection. This method is called
#' automatically when performing variable selection via
#' \code{\link[projpred:varsel]{varsel}} or
#' \code{\link[projpred:cv_varsel]{cv_varsel}}, so you will rarely need to call
#' it manually yourself.
#' 
#' @inheritParams posterior_predict.brmsfit
#' @param cvfun Optional cross-validation function
#' (see \code{\link[projpred:get-refmodel]{get_refmodel}} for details).
#' If \code{NULL} (the default), \code{cvfun} is defined internally
#' based on \code{\link{kfold.brmsfit}}.
#' @param kfold_seed A seed passed to \code{\link{kfold.brmsfit}}.
#' @param ... Further arguments passed to 
#' \code{\link[projpred:get-refmodel]{init_refmodel}}.
#' 
#' @details Note that the \code{extract_model_data} function used internally by
#'   \code{get_refmodel.brmsfit} ignores arguments \code{wrhs}, \code{orhs}, and
#'   \code{extract_y}. This is relevant for
#'   \code{\link[projpred:predict.refmodel]{predict.refmodel}}, for example.
#' 
#' @return A \code{refmodel} object to be used in conjunction with the
#'   \pkg{projpred} package.
#' 
#' @examples 
#' \dontrun{
#' # fit a simple model
#' fit <- brm(count ~ zAge + zBase * Trt,
#'            data = epilepsy, family = poisson())
#' summary(fit)
#' 
#' # The following code requires the 'projpred' package to be installed:
#' library(projpred)
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
get_refmodel.brmsfit <- function(object, newdata = NULL, resp = NULL, 
                                 cvfun = NULL, kfold_seed = NULL, ...) {
  require_package("projpred")
  dots <- list(...)
  resp <- validate_resp(resp, object, multiple = FALSE)
  formula <- formula(object)
  if (!is.null(resp)) {
    formula <- formula$forms[[resp]]
  }
  
  # prepare the family object for use in projpred
  family <- family(object, resp = resp)
  if (family$family == "bernoulli") {
    family$family <- "binomial"
  } else if (family$family == "gamma") {
    family$family <- "Gamma"
  }
  # For the augmented-data approach, do not re-define ordinal or categorical
  # families to preserve their family-specific extra arguments ("extra" meaning
  # "additionally to `link`") like `refcat` and `thresholds` (see ?brmsfamily):
  if (!(isTRUE(dots$aug_data) && is_polytomous(family))) {
    family <- get(family$family, mode = "function")(link = family$link)
  } else {
    # TODO: uncomment the lines below as soon as the
    # `extend_family_<family_name>` exist (in brms):
    # family <- get(paste0("extend_family_", family$family, mode = "function"))(
    #   family
    # )
  }
  
  # check if the model is supported by projpred
  bterms <- brmsterms(formula)
  if (length(bterms$dpars) > 1L && !conv_cats_dpars(family$family)) {
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
    dis <- as.data.frame(object, variable = dis)[[dis]]
  } else if (family$family == "Gamma") {
    dis <- paste0("shape", usc(resp))
    dis <- as.data.frame(object, variable = dis)[[dis]]
  }
  
  # prepare data passed to projpred
  data <- current_data(object, newdata, resp = resp, check_response = TRUE)
  attr(data, "terms") <- NULL
  
  # allows to handle additional arguments implicitly
  extract_model_data <- function(object, newdata = NULL, ...) {
    .extract_model_data(object, newdata = newdata, resp = resp, ...)
  }
  
  # Using the default prediction function from projpred is usually fine
  ref_predfun <- NULL
  if (isTRUE(dots$aug_data) && is_ordinal(family$family)) {
    stop2("This case is not yet supported.")
    # Use argument `incl_thres` of posterior_linpred() (and convert the
    # 3-dimensional array to an "augmented-rows" matrix)
    # TODO: uncomment the lines below as soon as arr2augmat() is exported
    # ref_predfun <- function(fit, newdata = NULL) {
    #   # Note: `transform = FALSE` is not needed, but included here for
    #   # consistency with projpred's default ref_predfun():
    #   linpred_out <- posterior_linpred(
    #     fit, transform = FALSE, newdata = newdata, incl_thres = TRUE
    #   )
    #   stopifnot(length(dim(linpred_out)) == 3L)
    #   # Since posterior_linpred() is supposed to include the offsets in its
    #   # result, subtract them here:
    #   # Observation weights are not needed here, so use `wrhs = NULL` to avoid
    #   # potential conflicts for a non-`NULL` default `wrhs`:
    #   offs <- extract_model_data(fit, newdata = newdata, wrhs = NULL)$offset
    #   if (length(offs)) {
    #     stopifnot(length(offs) %in% c(1L, dim(linpred_out)[2]))
    #     linpred_out <- sweep(linpred_out, 2, offs)
    #   }
    #   linpred_out <- projpred:::arr2augmat(linpred_out, margin_draws = 1)
    #   return(linpred_out)
    # }
  }
  
  # extract a list of K-fold sub-models
  if (is.null(cvfun)) {
    cvfun <- function(folds, ...) {
      if (is.null(kfold_seed)) {
        # Since kfold() doesn't seem to accept `seed = NULL`, set a random seed:
        kfold_seed <- sample.int(.Machine$integer.max, 1)
      }
      kfold(
        object, K = max(folds), save_fits = TRUE, folds = folds,
        seed = kfold_seed, ...
      )$fits[, "fit"]
    }
  } else {
    if (!is.function(cvfun)) {
      stop2("'cvfun' should be a function.")
    }
  }
  
  args <- nlist(
    object, data, formula, family, dis, ref_predfun = ref_predfun,
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
  offset <- as.vector(sdata[[paste0("offsets", usc_resp)]])
  weights <- as.vector(sdata[[paste0("weights", usc_resp)]])
  trials <- as.vector(sdata[[paste0("trials", usc_resp)]])
  stopifnot(!is.null(y))
  if (is_binary(family)) {
    trials <- rep(1, length(y))
  }
  if (!is.null(trials)) {
    if (!is.null(weights)) {
      stop2("Projpred cannot handle 'trials' and 'weights' at the same time.") 
    }
    weights <- trials
  }
  if (is.null(weights)) {
    weights <- rep(1, length(y))
  }
  if (is.null(offset)) {
    offset <- rep(0, length(y))
  }
  nlist(y, weights, offset)
}
