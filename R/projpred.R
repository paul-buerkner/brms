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
#' @param dis Passed to argument \code{dis} of
#'   \code{\link[projpred:init_refmodel]{init_refmodel}}, but leave this at
#'   \code{NULL} unless \pkg{projpred} complains about it.
#' @param latent See argument \code{latent} of
#'   \code{\link[projpred:extend_family]{extend_family}}. Setting this to
#'   \code{TRUE} requires a \pkg{projpred} version >= 2.4.0.
#' @param brms_seed A seed used to infer seeds for \code{\link{kfold.brmsfit}}
#'   and for sampling group-level effects for new levels (in multilevel models).
#'   If \code{NULL}, then \code{\link{set.seed}} is not called at all. If not
#'   \code{NULL}, then the pseudorandom number generator (PRNG) state is reset
#'   (to the state before calling this function) upon exiting this function.
#' @param ... Further arguments passed to
#' \code{\link[projpred:init_refmodel]{init_refmodel}}.
#'
#' @details The \code{extract_model_data} function used internally by
#'   \code{get_refmodel.brmsfit} ignores arguments \code{wrhs} and \code{orhs}
#'   (a warning is thrown if these are non-\code{NULL}). For example, arguments
#'   \code{weightsnew} and \code{offsetnew} of
#'   \code{\link[projpred:proj_linpred]{proj_linpred}},
#'   \code{\link[projpred:proj_predict]{proj_predict}}, and
#'   \code{\link[projpred:predict.refmodel]{predict.refmodel}} are passed to
#'   \code{wrhs} and \code{orhs}, respectively.
#'
#' @return A \code{refmodel} object to be used in conjunction with the
#'   \pkg{projpred} package.
#'
#' @examples
#' \dontrun{
#' # fit a simple model
#' fit <- brm(count ~ zAge + zBase * Trt,
#'   data = epilepsy, family = poisson()
#' )
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
#' @exportS3Method projpred::get_refmodel brmsfit
get_refmodel.brmsfit <- function(object, newdata = NULL, resp = NULL,
                                 cvfun = NULL, dis = NULL, latent = FALSE,
                                 brms_seed = NULL, ...) {
  require_package("projpred")
  object <- restructure(object)
  stopifnot_resp(object, resp)
  resp <- validate_resp(resp, object, multiple = FALSE)
  formula <- formula(object)
  if (!is.null(resp)) {
    formula <- formula$forms[[resp]]
  }

  # Infer "sub-seeds":
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    rng_state_old <- get(".Random.seed", envir = .GlobalEnv)
  }
  if (!is.null(brms_seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
    }
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
  aug_data <- (is_categorical(family) || is_ordinal(family)) && !latent
  # For the augmented-data and the latent approach, do not re-define the family
  # to preserve family-specific extra arguments ("extra" meaning "additionally
  # to `link`") like `refcat` and `thresholds` (see ?brmsfamily):
  if (!aug_data && !latent) {
    family <- get(family$family, mode = "function")(link = family$link)
  }

  # check if the model is supported by projpred
  bterms <- brmsterms(formula)
  if (length(bterms$dpars) > 1L && !conv_cats_dpars(family)) {
    stop2("Projpred does not support distributional models.")
  }
  if (conv_cats_dpars(family) && length(formula$pforms)) {
    stop2("Projpred does not support category-specific formulas.")
  }
  if (length(bterms$nlpars) > 0L) {
    stop2("Projpred does not support non-linear models.")
  }
  not_ok_term_types <- setdiff(all_term_types(), c("fe", "re", "offset", "sm"))
  if (any(not_ok_term_types %in% names(bterms$dpars$mu))) {
    stop2(
      "Projpred only supports standard multilevel and smoothing terms as ",
      "well as offsets."
    )
  }

  # only use the raw formula for selection of terms
  formula <- formula$formula
  # LHS should only contain the response variable
  formula[[2]] <- bterms$respform[[2]]

  # projpred requires the dispersion parameter if present
  if (is.null(dis) && !latent) {
    if (family$family == "gaussian") {
      dis <- paste0("sigma", usc(resp))
      dis <- as.data.frame(object, variable = dis)[[dis]]
    } else if (family$family == "Gamma") {
      dis <- paste0("shape", usc(resp))
      dis <- as.data.frame(object, variable = dis)[[dis]]
    }
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
    if (is_ordinal(family) && !latent) {
      c(lprd_args) <- list(incl_thres = TRUE)
    }
    out <- do_call(posterior_linpred, lprd_args)
    if (length(dim(out)) == 2) {
      out <- t(out)
    }
    out
  }

  if (utils::packageVersion("projpred") <= "2.0.2" && NROW(object$ranef)) {
    warning2(
      "In projpred versions <= 2.0.2, projpred's K-fold CV results may ",
      "not be reproducible for multilevel brms reference models."
    )
  }

  # extract a list of K-fold sub-models
  if (is.null(cvfun)) {
    cvfun <- function(folds, ...) {
      kfold(
        object,
        K = max(folds), save_fits = TRUE, folds = folds,
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
    projpred::get_refmodel(cvfit,
      resp = resp, dis = dis, latent = latent,
      brms_seed = brms_seed_k,
      called_from_cvrefbuilder = TRUE, ...
    )
  }

  # prepare data passed to projpred
  if (!is.null(newdata)) {
    warning2(
      "Argument 'newdata' of get_refmodel.brmsfit() is deprecated and ",
      "will be removed in the future."
    )
  }
  data <- current_data(
    object, newdata,
    resp = resp, check_response = TRUE,
    allow_new_levels = TRUE
  )
  attr(data, "terms") <- NULL
  args <- nlist(
    object, data, formula, family, dis, ref_predfun,
    cvfun, extract_model_data, cvrefbuilder, latent, ...
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
  } else if (latent) {
    require_package("projpred", "2.4.0")
    if (family$family == "cumulative") {
      args$latent_ilink <- latent_ilink_cumulative(
        object = object, family = family, bterms = bterms, resp = resp
      )
    }
    # TODO: If requested by users, add response-scale support for more families:
    # For response-scale support, they all need a specific `latent_ilink`
    # function; some families (those for which the response can be numeric) also
    # require specific `latent_ll_oscale` and `latent_ppd_oscale` functions. The
    # binomial family (and thereby also the brms::bernoulli() family) has
    # response-scale support implemented natively in projpred.
  }
  do_call(projpred::init_refmodel, args)
}

# auxiliary data required in predictions via projpred
# @return a named list with slots 'y', 'weights', and 'offset'
.extract_model_data <- function(object, newdata = NULL, resp = NULL,
                                extract_y = TRUE, wrhs = NULL, orhs = NULL,
                                ...) {
  stopifnot(is.brmsfit(object))
  resp <- validate_resp(resp, object, multiple = FALSE)
  if (utils::packageVersion("projpred") >= "2.8.0") {
    if (!is.null(wrhs)) warn_wrhs_orhs("wrhs")
    if (!is.null(orhs)) warn_wrhs_orhs("orhs")
  }

  # extract the response variable manually instead of from standata
  # so that it passes input checks of validate_newdata later on (#1314)
  formula <- formula(object)
  if (!is.null(resp)) {
    formula <- formula$forms[[resp]]
  }
  bterms <- brmsterms(formula)
  y <- NULL
  if (extract_y) {
    data <- current_data(
      object, newdata,
      resp = resp, check_response = TRUE,
      allow_new_levels = TRUE, req_vars = all.vars(bterms$respform)
    )
    y <- model.response(model.frame(bterms$respform, data, na.action = na.pass))
    y <- unname(y)
  }

  # extract relevant auxiliary data (offsets and weights (or numbers of trials))
  # call standata to ensure the correct format of the data
  # For this, we use `check_response = FALSE` and only include offsets and
  # weights (or numbers of trials) in `req_vars` because of issue #1457 (note
  # that all.vars(NULL) gives character(0), as desired).
  req_vars <- unlist(lapply(bterms$dpars, function(x) all.vars(x[["offset"]])))
  req_vars <- unique(req_vars)
  c(req_vars) <- all.vars(bterms$adforms$weights)
  c(req_vars) <- all.vars(bterms$adforms$trials)
  args <- nlist(
    object, newdata, resp,
    allow_new_levels = TRUE,
    check_response = FALSE,
    internal = TRUE,
    req_vars = req_vars
  )
  # NOTE: Missing weights don't cause an error here (see #1459)
  sdata <- do_call(standata, args)

  usc_resp <- usc(resp)
  N <- sdata[[paste0("N", usc_resp)]]
  weights <- as.vector(sdata[[paste0("weights", usc_resp)]])
  trials <- as.vector(sdata[[paste0("trials", usc_resp)]])
  if (is_binary(formula)) {
    trials <- rep(1, N)
  }
  if (!is.null(trials)) {
    if (!is.null(weights)) {
      stop2("Projpred cannot handle 'trials' and 'weights' at the same time.")
    }
    weights <- trials
  }
  if (is.null(weights)) {
    weights <- rep(1, N)
  }
  offset <- as.vector(sdata[[paste0("offsets", usc_resp)]])
  if (is.null(offset)) {
    offset <- rep(0, N)
  }
  nlist(y, weights, offset)
}

# Helper function for throwing a warning if argument `wrhs` or `orhs` is
# non-`NULL`.
warn_wrhs_orhs <- function(arg_nm) {
  warning2(
    "Argument `", arg_nm, "` is currently ignored. See section ",
    "'Details' of `?brms:::get_refmodel.brmsfit` for details."
  )
}

# Construct the inverse-link function required for the latent projection in case
# of the cumulative family.
#
# @param object See argument `object` of get_refmodel.brmsfit(), but here, the
#   `object` as modified inside of get_refmodel.brmsfit() is required.
# @param family The `family` object corresponding to `object` (taking `resp`
#   into account). Could be re-inferred from `object` and `resp`, but for
#   computational efficiency, this is avoided.
# @param bterms The `brmsterms` object corresponding to `object` (or rather
#   `object`'s `formula`, taking `resp` into account). Could be re-inferred from
#   `object` and `resp`, but for computational efficiency, this is avoided.
# @param resp See argument `resp` of get_refmodel.brmsfit(), but here, the
#   `resp` as modified inside of get_refmodel.brmsfit() is required.
#
# @return A function to be supplied to projpred::extend_family()'s argument
#   `latent_ilink`.
latent_ilink_cumulative <- function(object, family, bterms, resp) {
  stopifnot(!is.null(family$cats))
  draws_mat <- as_draws_matrix(object)
  thres_regex <- paste0("^b", usc(combine_prefix(bterms)), "_Intercept\\[")
  thres_draws <- prepare_draws(draws_mat, variable = thres_regex, regex = TRUE)
  if (ncol(thres_draws) > length(family$cats) - 1L) {
    stop2(
      "Currently, projpred does not support group-specific thresholds ",
      "(argument `gr` of resp_thres())."
    )
  }
  # Note: Currently, `disc` should always be constantly 1 because
  # distributional models are not allowed here.
  disc_regex <- paste0("^", "disc", resp, "$")
  disc_draws <- prepare_draws(draws_mat, variable = disc_regex, regex = TRUE)

  out <- function(lpreds, cl_ref, wdraws_ref = rep(1, length(cl_ref))) {
    thres_agg <- projpred::cl_agg(thres_draws, cl = cl_ref, wdraws = wdraws_ref)
    disc_agg <- projpred::cl_agg(disc_draws, cl = cl_ref, wdraws = wdraws_ref)
    disc_agg <- as.vector(disc_agg)
    lpreds_thres <- apply(thres_agg, 2, function(thres_agg_c) {
      # Notes on dimensionalities (with S_agg = `nrow(lpreds)`):
      # * `disc_agg` is a vector of length S_agg (because `disc` is not
      #   predicted here),
      # * `thres_agg` is S_agg x C_lat (with C_lat = `ncats - 1L` =
      #   `nthres`) and thus `thres_agg_c` is a vector of length S_agg,
      # * `lpreds` is S_agg x N (with N denoting the number of (possibly
      #   new) observations (not necessarily the original number of
      #   observations)).
      disc_agg * (thres_agg_c - lpreds)
    }, simplify = FALSE)
    # Coerce to an S_agg x N x C_lat array:
    lpreds_thres <- do.call(abind, c(lpreds_thres, rev.along = 0))
    # Transform to response space, yielding an S_agg x N x C_cat array:
    return(inv_link_cumulative(lpreds_thres, link = family$link))
  }
  # Free up some memory (keeping `draws_mat` would lead to unnecessary memory
  # usage because `draws_mat` would continue to live in the environment of the
  # returned function):
  rm(draws_mat)
  out
}
