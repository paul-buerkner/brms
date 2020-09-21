#' Support Functions for \pkg{emmeans}
#' 
#' Functions required for compatibility of \pkg{brms} with \pkg{emmeans}.
#' Users are not required to call these functions themselves. Instead,
#' they will be called automatically by the \code{emmeans} function
#' of the \pkg{emmeans} package.
#' 
#' In addition to the usual choices for \code{dpar}, the special value
#' \code{dpar = "mean"} requests that we use the expected values of the posterior 
#' predictive distribution, obtained via \code{\link{posterior_epred.brmsfit}}.
#' 
#' @name emmeans-brms-helpers
#' 
#' @inheritParams posterior_epred.brmsfit
#' @param data,trms,xlev,grid,vcov. Arguments required by \pkg{emmeans}.
#' @param ... Additional arguments passed to \pkg{emmeans}.
#' 
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
#' epred <- emmeans(fit, "disease", dpar = "mean")
#' summary(epred, point.est = mean)
#' }
NULL

# recover the predictors used in the population-level part of the model
#' @rdname emmeans-brms-helpers
recover_data.brmsfit <- function (object, data, resp = NULL, dpar = NULL, 
                                  nlpar = NULL, ...) {
  bterms <- .extract_par_terms(object, resp, dpar, nlpar)
  trms <- attr(model.frame(bterms$fe, data = object$data), "terms")
  # brms has no call component so the call is just a dummy
  emmeans::recover_data(call("brms"), trms, "na.omit", object$data, ...)
}

# Calculate the basis for making predictions. This is essentially the
# inside of the predict() function with new data on the link scale. 
# Transforming to response scale, if desired, is handled by emmeans.
#' @rdname emmeans-brms-helpers
emm_basis.brmsfit <- function (object, trms, xlev, grid, vcov., resp = NULL, 
                               dpar = NULL, nlpar = NULL, ...) {
  if (is_equal(dpar, "mean")) {
    post.beta <- posterior_epred(object, newdata = grid, re_formula = NA, ...)
    bhat <- apply(post.beta, 2, mean)
    V <- cov(post.beta)
    X <- diag(length(bhat))
    misc <- list()
  } else {
    bterms <- .extract_par_terms(object, resp, dpar, nlpar)
    m <- model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    contr <- lapply(object$data, function(x) attr(x, "contrasts"))
    contr <- contr[!ulapply(contr, is.null)]
    p <- combine_prefix(bterms)
    cols2remove <- str_if(is_ordinal(bterms), "(Intercept)")
    X <- get_model_matrix(trms, m, cols2remove, contrasts.arg = contr)
    nm <- paste0(usc(p, "suffix"), colnames(X))
    V <- vcov(object)[nm, nm, drop = FALSE]
    misc <- emmeans::.std.link.labels(bterms$family, list())
    post.beta <- as.matrix(object, pars = paste0("b_", nm), fixed = TRUE)
    bhat <- apply(post.beta, 2, mean)
  }
  attr(post.beta, "n.chains") <- object$fit@sim$chains
  nbasis <- matrix(NA)
  dfargs <- list()
  dffun <- function(k, dfargs) Inf
  nlist(X, bhat, nbasis, V, dffun, dfargs, misc, post.beta)
}

# extract terms of a specific predicted parameter in the model
.extract_par_terms <- function(object, resp = NULL, dpar = NULL, nlpar = NULL) {
  resp <- validate_resp(resp, object, multiple = FALSE)
  stopifnot_resp(object, resp)
  bterms <- brmsterms(formula(object))
  if (is.mvbrmsterms(bterms)) {
    bterms <- bterms$terms[[resp]]
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
    if (dpar == "mean") {
      # prepare posterior mean predictions as a special case
      # TODO: decide which variables to actually include in 'all_vars'
      all_vars <- NULL
      for (dp in all_dpars) {
        if (is.btl(bterms$dpars[[dp]])) {
          vars <- bterms$dpars[[dp]]$fe
        } else {
          vars <- bterms$dpars[[dp]]$covars
        }
        all_vars <- union(all_vars, vars)
      }
      for (nlp in all_nlpars) {
        if (is.btl(bterms$nlpars[[nlp]])) {
          vars <- bterms$nlpars[[nlp]]$fe
        } else {
          vars <- bterms$nlpars[[nlp]]$covars
        }
        all_vars <- union(all_vars, vars)
      } 
      out <- list(fe = terms_fe(str2formula(all_vars)))
      class(out) <- "btl"
    } else {
      if (!dpar %in% all_dpars) {
        stop2(
          "Distributional parameter '", dpar, "' is not part of the model.",
          "\nSupported parameters are: ", collapse_comma(all_dpars)
        )
      }
      out <- bterms$dpars[[dpar]]
    }
  } else {
    out <- bterms$dpars[["mu"]]
  }
  if (!is.btl(out)) {
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
