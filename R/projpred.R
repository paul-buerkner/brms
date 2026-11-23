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
#'   \code{TRUE} requires a \pkg{projpred} version >= 2.10.0. For the families
#'   \code{bernoulli}, \code{poisson}, \code{negbinomial}, \code{weibull}, and
#'   \code{lognormal}, \pkg{brms} defines the functions for response-scale
#'   analyses (see the latent-projection vignette of \pkg{projpred}).
#' @param brms_seed A seed used to infer seeds for \code{\link{kfold.brmsfit}}
#'   and for sampling group-level effects for new levels (in multilevel models).
#'   If \code{NULL}, then \code{\link{set.seed}} is not called at all. If not
#'   \code{NULL}, then the pseudorandom number generator (PRNG) state is reset
#'   (to the state before calling this function) upon exiting this function.
#' @param cvrefbuilder Optional function that creates the reference model for
#'   one fold of \eqn{K}-fold cross-validation
#'   (see \code{\link[projpred:init_refmodel]{init_refmodel}} for details).
#'   If \code{NULL} (the default), \code{cvrefbuilder} is defined internally
#'   and calls \code{get_refmodel} with the arguments of the original call.
#'   A custom function is needed, for example, if the reference model uses
#'   functions that depend on the posterior draws (see
#'   \code{\link[projpred:extend_family]{extend_family}}). It gets only the
#'   fold-wise fit, so set the other arguments of \code{get_refmodel} inside
#'   it.
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
#' @exportS3Method projpred::get_refmodel brmsfit
get_refmodel.brmsfit <- function(object, newdata = NULL, resp = NULL,
                                 cvfun = NULL, dis = NULL, latent = FALSE,
                                 brms_seed = NULL, cvrefbuilder = NULL, ...) {
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
  family_name <- family$family
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
  # get_effect() checks the terms of all dpars: if conv_cats_dpars(family),
  # the dpars are 'mu<category>' (no 'mu'); other dpars stopped above
  not_ok_term_types <- setdiff(
    all_term_types(), c("fe", "re", "offset", "sm", "sp")
  )
  has_not_ok <- ulapply(not_ok_term_types, function(t) {
    length(get_effect(bterms, target = t)) > 0L
  })
  if (any(has_not_ok)) {
    stop2("Projpred only supports standard multilevel, smoothing, and ",
          "monotonic terms as well as offsets.")
  }
  sp_terms <- get_effect(bterms, target = "sp")
  if (length(sp_terms)) {
    # 'sp' terms are mo(), me(), and mi() terms (see terms_sp())
    uni_me <- ulapply(sp_terms, attr, "uni_me")
    uni_mi <- ulapply(sp_terms, attr, "uni_mi")
    if (length(uni_me) || length(uni_mi)) {
      stop2("Projpred does not support me() or mi() terms.")
    }
    # check for the feature instead of a version number, so that development
    # versions of projpred with mo() support also work
    if (!"mo_effects" %in% getNamespaceExports("projpred")) {
      stop2("Monotonic terms need a version of projpred that supports them. ",
            "Please update package 'projpred'.")
    }
  }

  if (latent) {
    # the latent projection of projpred does not support observation weights
    mdata <- .extract_model_data(object, resp = resp, extract_y = FALSE)
    if (!isTRUE(all(mdata$weights == 1))) {
      stop2("Projpred does not support 'weights' or 'trials' other than 1 ",
            "in the latent projection.")
    }
  }

  # only use the raw formula for selection of terms
  formula <- formula$formula
  # LHS should only contain the response variable
  formula[[2]] <- bterms$respform[[2]]

  # keep the user's `dis` (possibly `NULL`) for the K-fold refits, so that
  # each refit computes its default `dis` from its own draws
  dis_user <- dis
  # projpred requires the dispersion parameter if present
  if (is.null(dis) && !latent) {
    if (family$family == "gaussian") {
      dis <- paste0("sigma", usc(resp))
      dis <- as.data.frame(object, variable = dis)[[dis]]
    } else if (family$family == "Gamma") {
      dis <- paste0("shape", usc(resp))
      dis <- as.data.frame(object, variable = dis)[[dis]]
    }
  } else if (is.null(dis) && latent) {
    # the SD of the latent Gaussian distribution, if the family has a dpar
    # for it (see latent_oscale_families())
    dis_dpar <- latent_dis_dpar(family_name, family$link)
    if (!is.null(dis_dpar)) {
      dis <- paste0(dis_dpar, usc(resp))
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

  if (is.null(cvrefbuilder)) {
    cvrefbuilder <- function(cvfit) {
      # For `brms_seed` in fold `cvfit$projpred_k` (= k) of K, choose a new
      # seed which is based on the original `brms_seed`:
      if (is.null(brms_seed)) {
        brms_seed_k <- NULL
      } else {
        brms_seed_k <- brms_seed + cvfit$projpred_k
      }
      projpred::get_refmodel(cvfit, resp = resp, dis = dis_user,
                             latent = latent, brms_seed = brms_seed_k,
                             called_from_cvrefbuilder = TRUE, ...)
    }
  } else {
    if (!is.function(cvrefbuilder)) {
      stop2("'cvrefbuilder' should be a function.")
    }
  }

  # prepare data passed to projpred
  if (!is.null(newdata)) {
    warning2("Argument 'newdata' of get_refmodel.brmsfit() is deprecated and ",
             "will be removed in the future.")
  }
  data <- current_data(
    object, newdata, resp = resp, check_response = TRUE,
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
    # the response-scale functions need argument 'cens' (projpred 2.10.0)
    require_package("projpred", "2.10.0")
    # brms adds its functions only if the user passes none of these arguments
    oscale_args <- c("latent_ilink", "latent_ll_oscale", "latent_ppd_oscale",
                     "latent_epred", "latent_y_unqs")
    if (family$family == "cumulative") {
      args$latent_ilink <- latent_ilink_cumulative(
        object = object, family = family, bterms = bterms, resp = resp
      )
    } else if (!any(oscale_args %in% names(args))) {
      c(args) <- latent_oscale_funs(
        object = object, family = family, family_name = family_name,
        bterms = bterms, resp = resp,
        warn = !isTRUE(args[["called_from_cvrefbuilder"]])
      )
    }
    # TODO: If requested by users, add response-scale support for more families.
    # Ordinal families other than cumulative() need a `latent_ilink` function
    # that gives the category probabilities (projpred has the other functions
    # for categorical responses). Families with a numeric response can be added
    # in latent_oscale_families().
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
      object, newdata, resp = resp, check_response = TRUE,
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
  warning2("Argument `", arg_nm, "` is currently ignored. See section ",
           "'Details' of `?brms:::get_refmodel.brmsfit` for details.")
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
    stop2("Currently, projpred does not support group-specific thresholds ",
          "(argument `gr` of resp_thres()).")
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

# Families for which brms defines the functions required for response-scale
# analyses with the latent projection. The functions call log_lik_<family>(),
# posterior_predict_<family>(), and posterior_epred_<family>() with `mu`
# replaced by the (projected) predictions. A family can be added here if these
# functions use only `mu`, dpars that are not predicted, and the data elements
# set in latent_oscale_funs() (`Y`, `cens`, `lb`, and `ub`), for example not
# `trials`.
# Element `dis` names a dpar that is the SD of the latent Gaussian distribution
# if the link is the identity (see latent_dis_dpar()).
latent_oscale_families <- function() {
  list(
    bernoulli = list(),
    poisson = list(),
    negbinomial = list(),
    weibull = list(),
    lognormal = list(dis = "sigma")
  )
}

# The dpar that is the SD of the latent Gaussian distribution. Its draws are the
# default of argument `dis` of get_refmodel.brmsfit(), and the functions of
# latent_oscale_funs() take it from their argument `dis`, so that submodels use
# their projected SD. With other links than the identity, the dpar is not the
# SD of the latent Gaussian distribution, and its draws are aggregated as those
# of other dpars.
# @return The name of the dpar or NULL.
latent_dis_dpar <- function(family_name, link) {
  if (link != "identity") {
    return(NULL)
  }
  latent_oscale_families()[[family_name]]$dis
}

# Construct the functions required for response-scale analyses with the latent
# projection for the families in latent_oscale_families().
#
# @param object, family, bterms, resp See latent_ilink_cumulative().
# @param family_name The brms name of the family.
# @param warn Warn if the model is not supported?
#
# @return NULL if projpred defines the functions for the family or if the
#   model is not supported. Otherwise, a named list with elements
#   `latent_ilink`, `latent_ll_oscale`, `latent_ppd_oscale`, and (if projpred
#   supports it) `latent_epred` to be passed to projpred::extend_family().
latent_oscale_funs <- function(object, family, family_name, bterms, resp,
                               warn = TRUE) {
  if (family_name == "binomial" || is_categorical(family) ||
      is_ordinal(family)) {
    # projpred defines the functions for these families
    return(NULL)
  }
  unsupported <- NULL
  if (!family_name %in% names(latent_oscale_families())) {
    unsupported <- paste0("family '", family_name, "'")
  } else if (has_interval_cens(bterms)) {
    unsupported <- "interval censoring"
  } else if (is.formula(bterms$adforms$rate)) {
    unsupported <- "addition argument 'rate'"
  } else {
    # use the model data instead of prepare_predictions(), so that this check
    # also works for fits without draws ('empty = TRUE')
    bounds <- trunc_bounds(bterms, data = object$data)
    if (length(unique(bounds$lb)) > 1L || length(unique(bounds$ub)) > 1L) {
      unsupported <- "truncation bounds that vary between observations"
    }
  }
  if (!is.null(unsupported)) {
    if (warn) {
      warning2("The response-scale functions of the latent projection are ",
               "not yet implemented in brms for ", unsupported, ", so the ",
               "response-scale results are NA or not valid. Pass these ",
               "functions via '...' (see ?projpred::extend_family).")
    }
    return(NULL)
  }

  prep <- prepare_predictions(object, resp = resp)
  # Free up some memory (keeping `object` would lead to unnecessary memory
  # usage because `object` would continue to live in the environment of the
  # returned functions, so a saved reference model or one sent to parallel
  # workers would contain the fit twice):
  rm(object)
  dis_dpar <- latent_dis_dpar(family_name, family$link)
  # other dpars are not predicted (no distributional models) and their draws
  # are aggregated in the same way as the projected draws; rep_len() is
  # required for fixed dpars of older models, which are single values
  prep$dpars$mu <- NULL
  aux_draws <- lapply(prep$dpars[setdiff(names(prep$dpars), dis_dpar)],
                      function(x) as.matrix(rep_len(x, prep$ndraws)))
  # truncation bounds are the same for all observations (see above)
  lb <- prep$data[["lb"]][1]
  ub <- prep$data[["ub"]][1]
  prep$data <- list()
  prep$old_order <- NULL
  log_lik_fun <- paste0("log_lik_", prep$family$fun)
  log_lik_fun <- get(log_lik_fun, asNamespace("brms"))
  pp_fun <- paste0("posterior_predict_", prep$family$fun)
  pp_fun <- get(pp_fun, asNamespace("brms"))

  # prepare a brmsprep object for (projected) predictions 'ilpreds'
  latent_prep <- function(ilpreds, dis, cl_ref, wdraws_ref, idxs = NULL) {
    out <- prep
    out$ndraws <- nrow(ilpreds)
    out$nobs <- ncol(ilpreds)
    out$dpars$mu <- ilpreds
    # cl_agg() is slow and not needed if each draw is its own cluster
    is_cl <- !identical(as.integer(cl_ref), seq_len(prep$ndraws))
    for (dp in names(aux_draws)) {
      draws <- aux_draws[[dp]]
      if (is_cl) {
        draws <- projpred::cl_agg(draws, cl = cl_ref, wdraws = wdraws_ref)
      }
      if (!is.null(idxs)) {
        draws <- draws[idxs, , drop = FALSE]
      }
      out$dpars[[dp]] <- as.vector(draws)
    }
    if (!is.null(dis_dpar)) {
      out$dpars[[dis_dpar]] <- dis
    }
    if (!is.null(lb)) {
      out$data$lb <- rep(lb, out$nobs)
    }
    if (!is.null(ub)) {
      out$data$ub <- rep(ub, out$nobs)
    }
    out
  }

  latent_ilink <- function(lpreds, cl_ref,
                           wdraws_ref = rep(1, length(cl_ref))) {
    inv_link(lpreds, family$link)
  }
  latent_ll_oscale <- function(ilpreds, dis, y_oscale,
                               wobs = rep(1, ncol(ilpreds)), cens, cl_ref,
                               wdraws_ref = rep(1, length(cl_ref))) {
    prep_ll <- latent_prep(
      ilpreds, dis, cl_ref = cl_ref, wdraws_ref = wdraws_ref
    )
    prep_ll$data$Y <- y_oscale
    if (!is.null(cens)) {
      prep_ll$data$cens <- prepare_cens(cens)
    }
    out <- lapply(seq_cols(ilpreds), log_lik_fun, prep = prep_ll)
    out <- do_call(cbind, out)
    out * rep(wobs, each = nrow(out))
  }
  cens_expr <- get_ad_expr(bterms, "cens", "cens")
  if (!is.null(cens_expr)) {
    # projpred evaluates this formula in the data and passes the result to
    # argument 'cens' of latent_ll_oscale()
    cens_var <- str2formula(
      cens_expr, env = environment(bterms$adforms$cens)
    )
    attr(latent_ll_oscale, "cens_var") <- cens_var
  }
  # predictions do not take censoring into account
  latent_ppd_oscale <- function(ilpreds_resamp, dis_resamp,
                                wobs = rep(1, ncol(ilpreds_resamp)), cl_ref,
                                wdraws_ref = rep(1, length(cl_ref)),
                                idxs_prjdraws) {
    prep_pp <- latent_prep(
      ilpreds_resamp, dis_resamp, cl_ref = cl_ref, wdraws_ref = wdraws_ref,
      idxs = idxs_prjdraws
    )
    out <- lapply(seq_cols(ilpreds_resamp), pp_fun, prep = prep_pp)
    do_call(cbind, out)
  }
  out <- nlist(latent_ilink, latent_ll_oscale, latent_ppd_oscale)
  # check for the feature instead of a version number, so that development
  # versions of projpred with 'latent_epred' also work
  if ("latent_epred" %in% names(formals(projpred::extend_family))) {
    posterior_epred_fun <- paste0("posterior_epred_", prep$family$family)
    posterior_epred_fun <- get(posterior_epred_fun, asNamespace("brms"))
    # the same as posterior_epred.brmsprep() for these families
    out$latent_epred <- function(ilpreds, dis, cl_ref,
                                 wdraws_ref = rep(1, length(cl_ref))) {
      prep_ep <- latent_prep(
        ilpreds, dis, cl_ref = cl_ref, wdraws_ref = wdraws_ref
      )
      if (is_trunc(prep_ep)) {
        return(posterior_epred_trunc(prep_ep))
      }
      posterior_epred_fun(prep_ep)
    }
  }
  out
}
