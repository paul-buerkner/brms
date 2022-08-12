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
#' (see \code{\link[projpred:get_refmodel]{get_refmodel}} for details).
#' If \code{NULL} (the default), \code{cvfun} is defined internally
#' based on \code{\link{kfold.brmsfit}}.
#' @param brms_seed A seed used to infer seeds for \code{\link{kfold.brmsfit}}
#'   and for sampling group-level effects for new levels (in multilevel models).
#' @param ... Further arguments passed to
#' \code{\link[projpred:init_refmodel]{init_refmodel}}.
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
                                 cvfun = NULL, brms_seed = NULL, ...) {
  require_package("projpred")
  object <- restructure(object)
  resp <- validate_resp(resp, object, multiple = FALSE)
  formula <- formula(object)
  if (!is.null(resp)) {
    formula <- formula$forms[[resp]]
  }

  # Infer "sub-seeds":
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    rng_state_old <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  }
  if (!is.null(brms_seed)) {
    set.seed(brms_seed)
  }
  kfold_seed <- sample.int(.Machine$integer.max, 1)
  refprd_seed <- sample.int(.Machine$integer.max, 1)

  # prepare the family object for use in projpred
  family <- family(object, resp = resp)
  if (family$family == "bernoulli") {
    family$family <- "binomial"
  } else if (family$family == "gamma") {
    family$family <- "Gamma"
  } else if (family$family == "beta") {
    family$family <- "Beta"
  }
  # For the augmented-data approach, do not re-define ordinal or categorical
  # families to preserve their family-specific extra arguments ("extra" meaning
  # "additionally to `link`") like `refcat` and `thresholds` (see ?brmsfamily):
  aug_data <- is_categorical(family) || is_ordinal(family)
  if (!aug_data) {
    family <- get(family$family, mode = "function")(link = family$link)
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

  # allows to handle additional arguments implicitly
  extract_model_data <- function(object, newdata = NULL, ...) {
    .extract_model_data(object, newdata = newdata, resp = resp, ...)
  }

  # The default `ref_predfun` from projpred does not set `allow_new_levels`, so
  # use a customized `ref_predfun` which also handles some preparations for the
  # augmented-data projection:
  ref_predfun <- function(fit, newdata = NULL) {
    # Setting a seed is necessary for reproducible sampling of group-level
    # effects for new levels:
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      rng_state_old <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
    }
    set.seed(refprd_seed)
    lprd_args <- nlist(
      object = fit, newdata, resp, allow_new_levels = TRUE,
      sample_new_levels = "gaussian"
    )
    if (is_ordinal(family)) {
      c(lprd_args) <- list(incl_thres = TRUE)
    }
    out <- do_call(posterior_linpred, lprd_args)
    if (length(dim(out)) == 2) {
      out <- t(out)
    }
    out
  }

  if (utils::packageVersion("projpred") <= "2.0.2" && NROW(object$ranef)) {
    warning2("In projpred versions <= 2.0.2, projpred's K-fold CV results may ",
             "not be reproducible for multilevel brms reference models.")
  }

  # extract a list of K-fold sub-models
  if (is.null(cvfun)) {
    cvfun <- function(folds, ...) {
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

  cvrefbuilder <- function(cvfit) {
    # For `brms_seed` in fold `cvfit$projpred_k` (= k) of K, choose a new seed
    # which is based on the original `brms_seed`:
    if (is.null(brms_seed)) {
      brms_seed_k <- NULL
    } else {
      brms_seed_k <- brms_seed + cvfit$projpred_k
    }
    projpred::get_refmodel(cvfit, resp = resp, brms_seed = brms_seed_k, ...)
  }

  # prepare data passed to projpred
  data <- current_data(
    object, newdata, resp = resp, check_response = TRUE,
    allow_new_levels = TRUE
  )
  attr(data, "terms") <- NULL
  args <- nlist(
    object, data, formula, family, dis, ref_predfun,
    cvfun, extract_model_data, cvrefbuilder, ...
  )
  if (aug_data) {
    c(args) <- list(
      augdat_link = get(paste0("link_", family$family), mode = "function"),
      augdat_ilink = get(paste0("inv_link_", family$family), mode = "function")
    )
    if (is_ordinal(family)) {
      c(args) <- list(
        augdat_args_link = list(link = family$link),
        augdat_args_ilink = list(link = family$link)
      )
    }
  }
  do_call(projpred::init_refmodel, args)
}

# auxiliary data required in predictions via projpred
# @return a named list with slots 'y', 'weights', and 'offset'
.extract_model_data <- function(object, newdata = NULL, resp = NULL, ...) {
  stopifnot(is.brmsfit(object))
  resp <- validate_resp(resp, object, multiple = FALSE)

  # extract the response variable manually instead of from make_standata
  # so that it passes input checks of validate_newdata later on (#1314)
  formula <- formula(object)
  if (!is.null(resp)) {
    formula <- formula$forms[[resp]]
  }
  respform <- brmsterms(formula)$respform
  data <- current_data(
    object, newdata, resp = resp, check_response = TRUE,
    allow_new_levels = TRUE
  )
  y <- unname(model.response(model.frame(respform, data, na.action = na.pass)))

  # extract relevant auxiliary data
  # call standata to ensure the correct format of the data
  args <- nlist(
    object, newdata, resp,
    allow_new_levels = TRUE,
    check_response = TRUE,
    internal = TRUE
  )
  sdata <- do_call(standata, args)

  usc_resp <- usc(resp)
  weights <- as.vector(sdata[[paste0("weights", usc_resp)]])
  trials <- as.vector(sdata[[paste0("trials", usc_resp)]])
  if (is_binary(formula)) {
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
  offset <- as.vector(sdata[[paste0("offsets", usc_resp)]])
  if (is.null(offset)) {
    offset <- rep(0, length(y))
  }
  nlist(y, weights, offset)
}
