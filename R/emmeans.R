#' Support Functions for \pkg{emmeans}
#' 
#' Functions required for compatibility of \pkg{brms} with \pkg{emmeans}.
#' Users are not required to call these functions themselves. Instead,
#' they will be called automatically by the \code{emmeans} function
#' of the \pkg{emmeans} package.
#' 
#' @name emmeans-brms-helpers
#' 
#' @inheritParams posterior_epred.brmsfit
#' @param re_formula Optional formula containing group-level effects to be
#'   considered in the prediction. If \code{NULL}, include all group-level
#'   effects; if \code{NA} (default), include no group-level effects.
#' @param epred Logical. If \code{TRUE} compute predictions of
#'   the posterior predictive distribution's mean
#'   (see \code{\link{posterior_epred.brmsfit}}) while ignoring
#'   arguments \code{dpar} and \code{nlpar}. Defaults to \code{FALSE}.
#' @param data,trms,xlev,grid,vcov. Arguments required by \pkg{emmeans}.
#' @param ... Additional arguments passed to \pkg{emmeans}.
#' 
#' @details 
#' In order to ensure compatibility of most \pkg{brms} models with
#' \pkg{emmeans}, predictions are not generated 'manually' via a design matrix
#' and coefficient vector, but rather via \code{\link{posterior_linpred.brmsfit}}. 
#' This appears to generally work well, but note that it produces an `.@linfct`
#' slot that contains the computed predictions as columns instead of the
#' coefficients.
#' 
#' @examples 
#' \dontrun{
#' fit <- brm(time | cens(censored) ~ age * sex + disease + (1|patient),
#'             data = kidney, family = lognormal())
#' summary(fit)           
#'
#' # summarize via 'emmeans'
#' library(emmeans)
#' rg <- ref_grid(fit)
#' em <- emmeans(rg, "disease")
#' summary(em, point.est = mean)
#' 
#' # obtain estimates for the posterior predictive distribution's mean
#' epred <- emmeans(fit, "disease", epred = TRUE)
#' summary(epred, point.est = mean)
#' }
NULL

# recover the variables used in the model predictions
# @param data only added to prevent it from being passed further via ...
#' @rdname emmeans-brms-helpers
recover_data.brmsfit <- function (object, data, resp = NULL, dpar = NULL, 
                                  nlpar = NULL, re_formula = NA,
                                  epred = FALSE, ...) {
  bterms <- .extract_par_terms(object, resp, dpar, nlpar, re_formula, epred)
  trms <- attr(model.frame(bterms$allvars, data = object$data), "terms")
  # brms has no call component so the call is just a dummy
  emmeans::recover_data(call("brms"), trms, "na.omit", data = object$data, ...)
}

# Calculate the basis for making predictions. In some sense, this is
# similar to the fitted() function with new data on the link scale. 
# Transforming to response scale, if desired, is handled by emmeans.
#' @rdname emmeans-brms-helpers
emm_basis.brmsfit <- function (object, trms, xlev, grid, vcov., resp = NULL, 
                               dpar = NULL, nlpar = NULL, re_formula = NA,
                               epred = FALSE, ...) {
  if (is_equal(dpar, "mean")) {
    # deprecated as of version 2.15.9
    warning2("dpar = 'mean' is deprecated. Please use epred = TRUE instead.")
    epred <- TRUE
    dpar <- NULL
  }
  epred <- as_one_logical(epred)
  if (epred) {
    post.beta <- posterior_epred(
      object, newdata = grid, re_formula = re_formula,
      incl_autocor = FALSE, ...
    )
    misc <- list()
  } else {
    post.beta <- posterior_linpred(
      object, newdata = grid, re_formula = re_formula, 
      resp = resp, dpar = dpar, nlpar = nlpar, 
      incl_autocor = FALSE, ...
    )
    bterms <- .extract_par_terms(object, resp, dpar, nlpar)
    misc <- emmeans::.std.link.labels(bterms$family, list())
  }
  attr(post.beta, "n.chains") <- object$fit@sim$chains
  bhat <- apply(post.beta, 2, mean)
  V <- cov(post.beta)
  X <- diag(length(bhat))
  nbasis <- matrix(NA)
  dfargs <- list()
  dffun <- function(k, dfargs) Inf
  environment(dffun) <- baseenv()
  nlist(X, bhat, nbasis, V, dffun, dfargs, misc, post.beta)
}

# extract terms of a specific predicted parameter in the model
.extract_par_terms <- function(object, resp = NULL, dpar = NULL, nlpar = NULL,
                               re_formula = NA, epred = FALSE) {
  if (is_equal(dpar, "mean")) {
    # deprecation warning already provided in emm_basis.brmsfit
    epred <- TRUE
    dpar <- NULL
  }
  epred <- as_one_logical(epred)
  resp <- validate_resp(resp, object, multiple = FALSE)
  stopifnot_resp(object, resp)
  new_formula <- update_re_terms(formula(object), re_formula)
  bterms <- brmsterms(new_formula, resp_rhs_all = FALSE)
  if (is.mvbrmsterms(bterms)) {
    bterms <- bterms$terms[[resp]]
  }
  if (epred) {
    return(bterms)
  }
  all_dpars <- names(bterms$dpars)
  all_nlpars <- names(bterms$nlpars)
  if (!is.null(nlpar)) {
    if (!is.null(dpar)) {
      stop2("'dpar' and 'nlpar' cannot be specified at the same time.")
    }
    nlpar <- as_one_character(nlpar)
    if (!nlpar %in% all_nlpars) {
      stop2(
        "Non-linear parameter '", nlpar, "' is not part of the model.",
        "\nSupported parameters are: ", collapse_comma(all_nlpars)
      )
    }
    out <- bterms$nlpars[[nlpar]]
  } else if (!is.null(dpar)) {
    dpar <- as_one_character(dpar)
    if (!dpar %in% all_dpars) {
      stop2(
        "Distributional parameter '", dpar, "' is not part of the model.",
        "\nSupported parameters are: ", collapse_comma(all_dpars)
      )
    }
    out <- bterms$dpars[[dpar]]
  } else {
    # neither dpar nor nlpar specified
    out <- bterms$dpars[["mu"]]
  }
  if (!is.btl(out)) {
    # TODO: can non-linear formulas be supported by emmeans as well?
    btl_dpars <- all_dpars[ulapply(bterms$dpars, is.btl)]
    btl_nlpars <- all_nlpars[ulapply(bterms$nlpars, is.btl)]
    stop2(
      "The select parameter is not predicted by a linear formula. ",
      "Use the 'dpar' and 'nlpar' arguments to select the ",
      "parameter for which marginal means should be computed.",
      "\nPredicted distributional parameters are: ", collapse_comma(btl_dpars),
      "\nPredicted non-linear parameters are: ", collapse_comma(btl_nlpars)
    )
  }
  out
}
