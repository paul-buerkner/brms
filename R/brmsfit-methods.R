# This file contains several extractor methods for brmsfit objects.
# A lot of other brmsfit methods have their own dedicated files.

#' Extract Population-Level Estimates
#'
#' Extract the population-level ('fixed') effects
#' from a \code{brmsfit} object.
#'
#' @aliases fixef
#'
#' @inheritParams predict.brmsfit
#' @param pars Optional names of coefficients to extract.
#'   By default, all coefficients are extracted.
#' @param ... Currently ignored.
#'
#' @return If \code{summary} is \code{TRUE}, a matrix returned
#'   by \code{\link{posterior_summary}} for the population-level effects.
#'   If \code{summary} is \code{FALSE}, a matrix with one row per
#'   posterior draw and one column per population-level effect.
#'
#' @examples
#' \dontrun{
#' fit <- brm(time | cens(censored) ~ age + sex + disease,
#'            data = kidney, family = "exponential")
#' fixef(fit)
#' # extract only some coefficients
#' fixef(fit, pars = c("age", "sex"))
#' }
#'
#' @method fixef brmsfit
#' @export
#' @export fixef
#' @importFrom nlme fixef
fixef.brmsfit <-  function(object, summary = TRUE, robust = FALSE,
                           probs = c(0.025, 0.975), pars = NULL, ...) {
  contains_draws(object)
  all_pars <- variables(object)
  fpars <- all_pars[grepl(fixef_pars(), all_pars)]
  if (!is.null(pars)) {
    pars <- as.character(pars)
    fpars <- fpars[sub("^[^_]+_", "", fpars) %in% pars]
  }
  if (!length(fpars)) {
    return(NULL)
  }
  out <- as.matrix(object, variable = fpars)
  colnames(out) <- gsub(fixef_pars(), "", fpars)
  if (summary) {
    out <- posterior_summary(out, probs, robust)
  }
  out
}

#' Covariance and Correlation Matrix of Population-Level Effects
#'
#' Get a point estimate of the covariance or
#' correlation matrix of population-level parameters
#'
#' @inheritParams fixef.brmsfit
#' @param correlation Logical; if \code{FALSE} (the default), compute
#'   the covariance matrix, if \code{TRUE}, compute the correlation matrix.
#'
#' @return covariance or correlation matrix of population-level parameters
#'
#' @details Estimates are obtained by calculating the maximum likelihood
#'   covariances (correlations) of the posterior draws.
#'
#' @examples
#' \dontrun{
#' fit <- brm(count ~ zAge + zBase * Trt + (1+Trt|visit),
#'            data = epilepsy, family = gaussian(), chains = 2)
#' vcov(fit)
#' }
#'
#' @export
vcov.brmsfit <- function(object, correlation = FALSE, pars = NULL, ...) {
  contains_draws(object)
  all_pars <- variables(object)
  fpars <- all_pars[grepl(fixef_pars(), all_pars)]
  if (!is.null(pars)) {
    pars <- as.character(pars)
    fpars <- intersect(fpars, paste0("b_", pars))
  }
  if (!length(fpars)) {
    return(NULL)
  }
  draws <- as.data.frame(object, variable = fpars)
  names(draws) <- sub(fixef_pars(), "", names(draws))
  if (correlation) {
    out <- cor(draws)
  } else {
    out <- cov(draws)
  }
  out
}

#' Extract Group-Level Estimates
#'
#' Extract the group-level ('random') effects of each level
#' from a \code{brmsfit} object.
#'
#' @aliases ranef
#'
#' @inheritParams fixef.brmsfit
#' @param groups Optional names of grouping variables
#'   for which to extract effects.
#' @param ... Currently ignored.
#'
#' @return A list of 3D arrays (one per grouping factor).
#'  If \code{summary} is \code{TRUE},
#'  the 1st dimension contains the factor levels,
#'  the 2nd dimension contains the summary statistics
#'  (see \code{\link{posterior_summary}}), and
#'  the 3rd dimension contains the group-level effects.
#'  If \code{summary} is \code{FALSE}, the 1st dimension contains
#'  the posterior draws, the 2nd dimension contains the factor levels,
#'  and the 3rd dimension contains the group-level effects.
#'
#' @examples
#' \dontrun{
#' fit <- brm(count ~ zAge + zBase * Trt + (1+Trt|visit),
#'            data = epilepsy, family = gaussian(), chains = 2)
#' ranef(fit)
#' }
#'
#' @method ranef brmsfit
#' @export
#' @export ranef
#' @importFrom nlme ranef
ranef.brmsfit <- function(object, summary = TRUE, robust = FALSE,
                          probs = c(0.025, 0.975), pars = NULL,
                          groups = NULL, ...) {
  contains_draws(object)
  object <- restructure(object)
  if (!nrow(object$ranef)) {
    stop2("The model does not contain group-level effects.")
  }
  all_pars <- variables(object)
  if (!is.null(pars)) {
    pars <- as.character(pars)
  }
  ranef <- object$ranef
  all_groups <- unique(ranef$group)
  if (!is.null(groups)) {
    groups <- as.character(groups)
    all_groups <- intersect(all_groups, groups)
  }
  out <- named_list(all_groups)
  for (g in all_groups) {
    r <- subset2(ranef, group = g)
    coefs <- paste0(usc(combine_prefix(r), "suffix"), r$coef)
    rpars <- all_pars[grepl(paste0("^r_", g, "(__.+\\[|\\[)"), all_pars)]
    if (!is.null(pars)) {
      coefs <- coefs[r$coef %in% pars]
      if (!length(coefs)) {
        next
      }
      regex <- paste0("(", escape_all(coefs), ")", collapse = "|")
      regex <- paste0(",", regex, "\\]$")
      rpars <- rpars[grepl(regex, rpars)]
    }
    levels <- attr(ranef, "levels")[[g]]
    if (length(rpars)) {
      # draws of varying coefficients were saved
      out[[g]] <- as.matrix(object, variable = rpars)
      dim(out[[g]]) <- c(nrow(out[[g]]), length(levels), length(coefs))
    } else {
      # draws of varying coefficients were not saved
      out[[g]] <- array(dim = c(ndraws(object), length(levels), length(coefs)))
    }
    dimnames(out[[g]])[2:3] <- list(levels, coefs)
    if (summary) {
      out[[g]] <- posterior_summary(out[[g]], probs, robust)
    }
  }
  rmNULL(out, recursive = FALSE)
}

#' Extract Model Coefficients
#'
#' Extract model coefficients, which are the sum of population-level
#' effects and corresponding group-level effects
#'
#' @inheritParams ranef.brmsfit
#' @param ... Further arguments passed to \code{\link{fixef.brmsfit}}
#'   and \code{\link{ranef.brmsfit}}.
#'
#' @return A list of 3D arrays (one per grouping factor).
#'  If \code{summary} is \code{TRUE},
#'  the 1st dimension contains the factor levels,
#'  the 2nd dimension contains the summary statistics
#'  (see \code{\link{posterior_summary}}), and
#'  the 3rd dimension contains the group-level effects.
#'  If \code{summary} is \code{FALSE}, the 1st dimension contains
#'  the posterior draws, the 2nd dimension contains the factor levels,
#'  and the 3rd dimension contains the group-level effects.
#'
#' @examples
#' \dontrun{
#' fit <- brm(count ~ zAge + zBase * Trt + (1+Trt|visit),
#'            data = epilepsy, family = gaussian(), chains = 2)
#' ## extract population and group-level coefficients separately
#' fixef(fit)
#' ranef(fit)
#' ## extract combined coefficients
#' coef(fit)
#' }
#'
#' @export
coef.brmsfit <- function(object, summary = TRUE, robust = FALSE,
                         probs = c(0.025, 0.975), ...) {
  contains_draws(object)
  object <- restructure(object)
  if (!nrow(object$ranef)) {
    stop2("No group-level effects detected. Call method ",
          "'fixef' to access population-level effects.")
  }
  fixef <- fixef(object, summary = FALSE, ...)
  coef <- ranef(object, summary = FALSE, ...)
  # add missing coefficients to fixef
  all_ranef_names <- unique(ulapply(coef, function(x) dimnames(x)[[3]]))
  fixef_names <- colnames(fixef)
  fixef_no_digits <- get_matches("^[^\\[]+", fixef_names)
  miss_fixef <- setdiff(all_ranef_names, fixef_names)
  miss_fixef_no_digits <- get_matches("^[^\\[]+", miss_fixef)
  new_fixef <- named_list(miss_fixef)
  for (k in seq_along(miss_fixef)) {
    # digits occur in ordinal models with category specific effects
    match_fixef <- match(miss_fixef_no_digits[k], fixef_names)
    if (!is.na(match_fixef)) {
      new_fixef[[k]] <- fixef[, match_fixef]
    } else if (!miss_fixef[k] %in% fixef_no_digits) {
      new_fixef[[k]] <- 0
    }
  }
  rm_fixef <- fixef_names %in% miss_fixef_no_digits
  fixef <- fixef[, !rm_fixef, drop = FALSE]
  fixef <- do_call(cbind, c(list(fixef), rmNULL(new_fixef)))

  for (g in names(coef)) {
    # add missing coefficients to ranef
    ranef_names <- dimnames(coef[[g]])[[3]]
    ranef_no_digits <- get_matches("^[^\\[]+", ranef_names)
    miss_ranef <- setdiff(fixef_names, ranef_names)
    miss_ranef_no_digits <- get_matches("^[^\\[]+", miss_ranef)
    new_ranef <- named_list(miss_ranef)
    for (k in seq_along(miss_ranef)) {
      # digits occur in ordinal models with category specific effects
      match_ranef <- match(miss_ranef_no_digits[k], ranef_names)
      if (!is.na(match_ranef)) {
        new_ranef[[k]] <- coef[[g]][, , match_ranef]
      } else if (!miss_ranef[k] %in% ranef_no_digits) {
        new_ranef[[k]] <- array(0, dim = dim(coef[[g]])[1:2])
      }
    }
    rm_ranef <- ranef_names %in% miss_ranef_no_digits
    coef[[g]] <- coef[[g]][, , !rm_ranef, drop = FALSE]
    coef[[g]] <- abind(c(list(coef[[g]]), rmNULL(new_ranef)))
    for (nm in dimnames(coef[[g]])[[3]]) {
      is_ord_intercept <- grepl("(^|_)Intercept\\[[[:digit:]]+\\]$", nm)
      if (is_ord_intercept) {
        # correct the sign of thresholds in ordinal models
        resp <- if (is_mv(object)) get_matches("^[^_]+", nm)
        family <- family(object, resp = resp)$family
        if (has_thres_minus_eta(family)) {
          coef[[g]][, , nm] <- fixef[, nm] - coef[[g]][, , nm]
        } else if (has_eta_minus_thres(family)) {
          coef[[g]][, , nm] <- coef[[g]][, , nm] - fixef[, nm]
        } else {
          coef[[g]][, , nm] <- fixef[, nm] + coef[[g]][, , nm]
        }
      } else {
        coef[[g]][, , nm] <- fixef[, nm] + coef[[g]][, , nm]
      }
    }
    if (summary) {
      coef[[g]] <- posterior_summary(coef[[g]], probs, robust)
    }
  }
  coef
}

#' Extract Variance and Correlation Components
#'
#' This function calculates the estimated standard deviations,
#' correlations and covariances of the group-level terms
#' in a multilevel model of class \code{brmsfit}.
#' For linear models, the residual standard deviations,
#' correlations and covariances are also returned.
#'
#' @aliases VarCorr
#'
#' @param x An object of class \code{brmsfit}.
#' @inheritParams fixef.brmsfit
#' @param sigma Ignored (included for compatibility with
#'  \code{\link[nlme:VarCorr]{VarCorr}}).
#' @param ... Currently ignored.
#'
#' @return A list of lists (one per grouping factor), each with
#' three elements: a matrix containing the standard deviations,
#' an array containing the correlation matrix, and an array
#' containing the covariance matrix with variances on the diagonal.
#'
#' @examples
#' \dontrun{
#' fit <- brm(count ~ zAge + zBase * Trt + (1+Trt|visit),
#'            data = epilepsy, family = gaussian(), chains = 2)
#' VarCorr(fit)
#' }
#'
#' @method VarCorr brmsfit
#' @import abind abind
#' @importFrom nlme VarCorr
#' @export VarCorr
#' @export
VarCorr.brmsfit <- function(x, sigma = 1, summary = TRUE, robust = FALSE,
                            probs = c(0.025, 0.975), ...) {
  contains_draws(x)
  x <- restructure(x)
  if (!(nrow(x$ranef) || any(grepl("^sigma($|_)", variables(x))))) {
    stop2("The model does not contain covariance matrices.")
  }
  .VarCorr <- function(y) {
    # extract draws for sd, cor and cov
    out <- list(sd = as.matrix(x, variable = y$sd_pars))
    colnames(out$sd) <- y$rnames
    # compute correlation and covariance matrices
    found_cor_pars <- intersect(y$cor_pars, variables(x))
    if (length(found_cor_pars)) {
      cor <- as.matrix(x, variable = found_cor_pars)
      if (length(found_cor_pars) < length(y$cor_pars)) {
        # some correlations are missing and will be replaced by 0
        cor_all <- matrix(0, nrow = nrow(cor), ncol = length(y$cor_pars))
        names(cor_all) <- y$cor_pars
        for (i in seq_len(ncol(cor_all))) {
          found <- match(names(cor_all)[i], colnames(cor))
          if (!is.na(found)) {
            cor_all[, i] <- cor[, found]
          }
        }
        cor <- cor_all
      }
      out$cor <- get_cor_matrix(cor = cor)
      out$cov <- get_cov_matrix(sd = out$sd, cor = cor)
      dimnames(out$cor)[2:3] <- list(y$rnames, y$rnames)
      dimnames(out$cov)[2:3] <- list(y$rnames, y$rnames)
      if (summary) {
        out$cor <- posterior_summary(out$cor, probs, robust)
        out$cov <- posterior_summary(out$cov, probs, robust)
      }
    }
    if (summary) {
      out$sd <- posterior_summary(out$sd, probs, robust)
    }
    return(out)
  }

  if (nrow(x$ranef)) {
    get_names <- function(group) {
      # get names of group-level parameters
      r <- subset2(x$ranef, group = group)
      rnames <- as.vector(get_rnames(r))
      cor_type <- paste0("cor_", group)
      sd_pars <- paste0("sd_", group, "__", rnames)
      cor_pars <- get_cornames(rnames, cor_type, brackets = FALSE)
      nlist(rnames, sd_pars, cor_pars)
    }
    group <- unique(x$ranef$group)
    tmp <- lapply(group, get_names)
    names(tmp) <- group
  } else {
    tmp <- list()
  }
  # include residual variances in the output as well
  bterms <- brmsterms(x$formula)
  if (is.brmsterms(bterms)) {
    if (simple_sigma(bterms) && !is.mixfamily(x$family)) {
      tmp_resid <- list(rnames = bterms$resp, sd_pars = "sigma")
      tmp <- c(tmp, residual__ = list(tmp_resid))
    }
  } else if (is.mvbrmsterms(bterms)) {
    simple_sigma <- ulapply(bterms$terms, simple_sigma)
    pred_sigma <- ulapply(bterms$terms, pred_sigma)
    is_mix <- ulapply(x$family, is.mixfamily)
    if (any(simple_sigma) && !any(pred_sigma) && !any(is_mix)) {
      resps <- bterms$responses[simple_sigma]
      sd_pars <- paste0("sigma_", resps)
      if (bterms$rescor) {
        cor_pars <- get_cornames(resps, type = "rescor", brackets = FALSE)
      } else {
        cor_pars <- character(0)
      }
      tmp_resid <- nlist(rnames = resps, sd_pars, cor_pars)
      tmp <- c(tmp, residual__ = list(tmp_resid))
    }
  }
  lapply(tmp, .VarCorr)
}

#' @export
model.frame.brmsfit <- function(formula, ...) {
  formula$data
}

#' (Deprecated) Number of Posterior Samples
#'
#' Extract the number of posterior samples (draws) stored in a fitted Bayesian
#' model. Method \code{nsamples} is deprecated. Please use \code{ndraws}
#' instead.
#'
#' @aliases nsamples
#'
#' @param object An object of class \code{brmsfit}.
#' @param subset An optional integer vector defining a subset of samples
#'   to be considered.
#' @param incl_warmup A flag indicating whether to also count warmup / burn-in
#'   samples.
#' @param ... Currently ignored.
#'
#' @method nsamples brmsfit
#' @export
#' @export nsamples
#' @importFrom rstantools nsamples
nsamples.brmsfit <- function(object, subset = NULL,
                             incl_warmup = FALSE, ...) {
  warning2("'nsamples.brmsfit' is deprecated. Please use 'ndraws' instead.")
  if (!is(object$fit, "stanfit") || !length(object$fit@sim)) {
    out <- 0
  } else {
    ntsamples <- object$fit@sim$n_save[1]
    if (!incl_warmup) {
      ntsamples <- ntsamples - object$fit@sim$warmup2[1]
    }
    ntsamples <- ntsamples * object$fit@sim$chains
    if (length(subset)) {
      out <- length(subset)
      if (out > ntsamples || max(subset) > ntsamples) {
        stop2("Argument 'subset' is invalid.")
      }
    } else {
      out <- ntsamples
    }
  }
  out
}

#' @export
nobs.brmsfit <- function(object, resp = NULL, ...) {
  if (is_mv(object) && length(resp)) {
    resp <- validate_resp(resp, object, multiple = FALSE)
    bterms <- brmsterms(object$formula$forms[[resp]])
    out <- nrow(subset_data(model.frame(object), bterms))
  } else {
    out <- nrow(model.frame(object))
  }
  out
}

#' Number of Grouping Factor Levels
#'
#' Extract the number of levels of one or more grouping factors.
#'
#' @aliases ngrps.brmsfit
#'
#' @param object An \R object.
#' @param ... Currently ignored.
#'
#' @return A named list containing the number of levels per
#'   grouping factor.
#'
#' @export
ngrps.brmsfit <- function(object, ...) {
  object <- restructure(object)
  if (nrow(object$ranef)) {
    out <- lapply(attr(object$ranef, "levels"), length)
  } else {
    out <- NULL
  }
  out
}

#' @rdname ngrps.brmsfit
#' @export
ngrps <- function(object, ...) {
  UseMethod("ngrps")
}

#' @export
formula.brmsfit <- function(x, ...) {
  x$formula
}

#' @export
getCall.brmsfit <- function(x, ...) {
  x$formula
}

#' Extract Model Family Objects
#'
#' @inheritParams posterior_predict.brmsfit
#' @param ... Currently unused.
#'
#' @return A \code{brmsfamily} object
#' or a list of such objects for multivariate models.
#'
#' @export
family.brmsfit <- function(object, resp = NULL, ...) {
  resp <- validate_resp(resp, object)
  if (!is.null(resp)) {
    # multivariate model
    family <- from_list(object$formula$forms[resp], "family")
    if (length(resp) == 1L) {
      family <- family[[1]]
    }
  } else {
    # univariate model
    family <- object$formula$family
    if (is.null(family)) {
      family <- object$family
    }
  }
  family
}

#' Expose user-defined \pkg{Stan} functions
#'
#' Export user-defined \pkg{Stan} function and
#' optionally vectorize them. For more details see
#' \code{\link[rstan:expose_stan_functions]{expose_stan_functions}}.
#'
#' @param x An object of class \code{brmsfit}.
#' @param vectorize Logical; Indicates if the exposed functions
#'   should be vectorized via \code{\link{Vectorize}}.
#'   Defaults to \code{FALSE}.
#' @param env Environment where the functions should be made
#'   available. Defaults to the global environment.
#' @param ... Further arguments passed to
#'   \code{\link[rstan:expose_stan_functions]{expose_stan_functions}}.
#'
#' @export
expose_functions.brmsfit <- function(x, vectorize = FALSE,
                                     env = globalenv(), ...) {
  vectorize <- as_one_logical(vectorize)
  stanmodel <- compiled_model(x)
  if (x$backend == "cmdstanr") {
    if ("expose_functions" %in% names(stanmodel) & !vectorize) {
      funs <- .expose_functions_cmdstanr(
        stanmodel, vectorize = vectorize, env = env, ...
      )
    } else {
      # older versions of cmdstanr cannot export stan functions (#1176)
      scode  <- strsplit(stancode(x), "\n")[[1]]
      data_line <- grep("^data[ ]+\\{$", scode)
      scode <- paste0(c(scode[seq_len(data_line - 1)], "\n"), collapse = "\n")
      stanmodel <- tempfile(fileext = ".stan")
      cat(scode, file = stanmodel)
      funs <- .expose_functions_rstan(
        stanmodel, vectorize = vectorize, env = env, ...
      )
    }
  } else {
    funs <- .expose_functions_rstan(
      stanmodel, vectorize = vectorize, env = env, ...
    )
  }
  invisible(funs)
}

# expose stan functions via rstan
.expose_functions_rstan <- function(stanmodel, vectorize = FALSE,
                                    env = globalenv(), ...) {
  if (vectorize) {
    funs <- rstan::expose_stan_functions(stanmodel, env = environment(), ...)
    for (i in seq_along(funs)) {
      FUN <- Vectorize(get(funs[i], mode = "function"))
      assign(funs[i], FUN, pos = env)
    }
  } else {
    funs <- rstan::expose_stan_functions(stanmodel, env = env, ...)
  }
  funs
}

# expose stan functions via cmdstanr
.expose_functions_cmdstanr <- function(stanmodel, vectorize = FALSE,
                                       env = globalenv(), ...) {
  # TODO: support argument 'vectorize'
  # TODO: support argument 'env'
  stopifnot(!vectorize)
  funs <- stanmodel$expose_functions(global = TRUE)
  funs
}

#' @rdname expose_functions.brmsfit
#' @export
expose_functions <- function(x, ...) {
  UseMethod("expose_functions")
}
