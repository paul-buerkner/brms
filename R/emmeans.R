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
#'   If you have specified a response transformation within the formula,
#'   you need to set \code{epred} to \code{TRUE} for \pkg{emmeans} to
#'   detect this transformation.
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
#' fit1 <- brm(time | cens(censored) ~ age * sex + disease + (1|patient),
#'             data = kidney, family = lognormal())
#' summary(fit1)
#'
#' # summarize via 'emmeans'
#' library(emmeans)
#' rg <- ref_grid(fit1)
#' em <- emmeans(rg, "disease")
#' summary(em, point.est = mean)
#'
#' # obtain estimates for the posterior predictive distribution's mean
#' epred <- emmeans(fit1, "disease", epred = TRUE)
#' summary(epred, point.est = mean)
#'
#'
#' # model with transformed response variable
#' fit2 <- brm(log(mpg) ~ factor(cyl), data = mtcars)
#' summary(fit2)
#'
#' # results will be on the log scale by default
#' emmeans(fit2, ~ cyl)
#' # log transform is detected and can be adjusted automatically
#' emmeans(fit2, ~ cyl, epred = TRUE, type = "response")
#' }
NULL

# recover the variables used in the model predictions
# @param data only added to prevent it from being passed further via ...
#' @rdname emmeans-brms-helpers
recover_data.brmsfit <- function(object, data, resp = NULL, dpar = NULL,
                                 nlpar = NULL, re_formula = NA,
                                 epred = FALSE, ...) {
  bterms <- .extract_par_terms(
    object, resp = resp, dpar = dpar, nlpar = nlpar,
    re_formula = re_formula, epred = epred
  )
  trms <- terms(bterms$allvars, data = object$data)
  # brms has no call component so the call is just a dummy for the most part
  cl <- call("brms")
  if (epred) {
    # fixes issue #1360 for in-formula response transformations
    cl$formula <- bterms$respform
  }
  emmeans::recover_data(cl, trms, "na.omit", data = object$data, ...)
}

# Calculate the basis for making predictions. In some sense, this is
# similar to the fitted() function with new data on the link scale.
# Transforming to response scale, if desired, is handled by emmeans.
#' @rdname emmeans-brms-helpers
emm_basis.brmsfit <- function(object, trms, xlev, grid, vcov., resp = NULL,
                              dpar = NULL, nlpar = NULL, re_formula = NA,
                              epred = FALSE, ...) {
  if (is_equal(dpar, "mean")) {
    # deprecated as of version 2.15.9
    warning2("dpar = 'mean' is deprecated. Please use epred = TRUE instead.")
    epred <- TRUE
    dpar <- NULL
  }
  epred <- as_one_logical(epred)
  bterms <- .extract_par_terms(
    object, resp = resp, dpar = dpar, nlpar = nlpar,
    re_formula = re_formula, epred = epred
  )
  if (epred) {
    post.beta <- posterior_epred(
      object, newdata = grid, re_formula = re_formula,
      resp = resp, incl_autocor = FALSE, ...
    )
  } else {
    req_vars <- all_vars(bterms$allvars)
    post.beta <- posterior_linpred(
      object, newdata = grid, re_formula = re_formula,
      resp = resp, dpar = dpar, nlpar = nlpar,
      incl_autocor = FALSE, req_vars = req_vars,
      # offsets are handled by emmeans (#1096)
      transform = FALSE, offset = FALSE, ...
    )
  }
  if (anyNA(post.beta)) {
    stop2("emm_basis.brmsfit created NAs. Please check your reference grid.")
  }
  misc <- bterms$.misc
  if (length(dim(post.beta)) == 3L) {
    # reshape to a 2D matrix, for example, in multivariate models
    ynames <- dimnames(post.beta)[[3]]
    if (is.null(ynames)) {
      ynames <- as.character(seq_len(dim(post.beta)[3]))
    }
    dims <- dim(post.beta)
    post.beta <- matrix(post.beta, ncol = prod(dims[2:3]))
    misc$ylevs = list(rep.meas = ynames)
  }
  attr(post.beta, "n.chains") <- object$fit@sim$chains
  X <- diag(ncol(post.beta))
  bhat <- apply(post.beta, 2, mean)
  V <- cov(post.beta)
  nbasis <- matrix(NA)
  dfargs <- list()
  dffun <- function(k, dfargs) Inf
  environment(dffun) <- baseenv()
  nlist(X, bhat, nbasis, V, dffun, dfargs, misc, post.beta)
}

# extract terms of specific predicted parameter(s) in the model
# currently, the only slots that matter in the returned object are
# allvars: formula with all required variables on the right-hand side
# .misc: a named list with additional info to be interpreted by emmeans
.extract_par_terms <- function(x, ...) {
  UseMethod(".extract_par_terms")
}

#' @export
.extract_par_terms.brmsfit <- function(x, resp = NULL, re_formula = NA,
                                       dpar = NULL, epred = FALSE, ...) {
  if (is_equal(dpar, "mean")) {
    # deprecation warning already provided in emm_basis.brmsfit
    epred <- TRUE
    dpar <- NULL
  }
  resp <- validate_resp(resp, x)
  new_formula <- update_re_terms(formula(x), re_formula)
  # autocorrelation terms are always excluded for emmeans predictions (#1424)
  new_formula <- exclude_terms(new_formula, incl_autocor = FALSE)
  bterms <- brmsterms(new_formula, resp_rhs_all = FALSE)
  if (is_ordinal(bterms)) {
    warning2("brms' emmeans support for ordinal models is experimental ",
             "and currently ignores the threshold parameters.")
  }
  .extract_par_terms(bterms, resp = resp, dpar = dpar, epred = epred, ...)
}

#' @export
.extract_par_terms.mvbrmsterms <- function(x, resp, epred, ...) {
  stopifnot(is.character(resp))
  epred <- as_one_logical(epred)
  out <- x
  # only use selected univariate models
  out$terms <- out$terms[resp]
  if (epred) {
    out$allvars <- allvars_formula(lapply(out$terms, get_allvars))
    out$.misc <- list()
    return(out)
  }
  for (i in seq_along(out$terms)) {
    out$terms[[i]] <- .extract_par_terms(out$terms[[i]], epred = epred, ...)
  }
  out$allvars <- allvars_formula(lapply(out$terms, get_allvars))
  misc_list <- unique(from_list(out$terms, ".misc"))
  if (length(misc_list) > 1L){
    stop2("brms' emmeans support for multivariate models is limited ",
          "to cases where all univariate models have the same family.")
  }
  out$.misc <- misc_list[[1]]
  out
}

#' @export
.extract_par_terms.brmsterms <- function(x, dpar, nlpar, epred, ...) {
  epred <- as_one_logical(epred)
  all_dpars <- names(x$dpars)
  all_nlpars <- names(x$nlpars)
  out <- x
  if (epred) {
    out$.misc <- list()
    return(out)
  }
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
    out <- x$nlpars[[nlpar]]
  } else if (!is.null(dpar)) {
    dpar <- as_one_character(dpar)
    if (!dpar %in% all_dpars) {
      stop2(
        "Distributional parameter '", dpar, "' is not part of the model.",
        "\nSupported parameters are: ", collapse_comma(all_dpars)
      )
    }
    out <- x$dpars[[dpar]]
  } else {
    # extract 'mu' parameter by default
    if (!"mu" %in% names(x$dpars)) {
      # concerns categorical-like and mixture models
      stop2("emmeans is not yet supported for this brms model.")
    }
    out <- x$dpars[["mu"]]
  }
  if (!is.null(out$offset)) {
    # ensure that offsets are detected by emmeans (#1096)
    out$allvars <- allvars_formula(out$allvars, out$offset)
  }
  out$.misc <- emmeans::.std.link.labels(out$family, list())
  out
}
