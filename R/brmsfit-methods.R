#' @export
parnames.brmsfit <- function(x, ...) {
  out <- dimnames(x$fit)
  if (is.list(out)) {
    out <- out$parameters
  }
  out
}

#' Extract Population-Level Estimates
#' 
#' Extract the population-level ('fixed') effects 
#' from a \code{brmsfit} object. 
#' 
#' @aliases fixef
#' 
#' @param object An object of class \code{brmsfit}.
#' @param ... Currently ignored.
#' @inheritParams predict.brmsfit
#' 
#' @return If \code{summary} is \code{TRUE}, a matrix with one row per 
#'   population-level effect and one column per calculated estimate. 
#'   If \code{summary} is \code{FALSE}, a matrix with one row per 
#'   posterior sample and one column per population-level effect.
#'   
#' @author Paul-Christian Buerkner \email{paul.buerkner@gmail.com}
#' 
#' @examples
#' \dontrun{
#' fit <- brm(time | cens(censored) ~ age + sex + disease, 
#'            data = kidney, family = "exponential")
#' fixef(fit)
#' }
#' 
#' @method fixef brmsfit
#' @export
#' @export fixef
#' @importFrom nlme fixef
fixef.brmsfit <-  function(object, summary = TRUE, robust = FALSE,
                           probs = c(0.025, 0.975), ...) {
  contains_samples(object)
  pars <- parnames(object)
  fpars <- pars[grepl(fixef_pars(), pars)]
  if (!length(fpars)) {
    stop2("The model does not contain population-level effects.")
  }
  out <- as.matrix(object, pars = fpars, exact_match = TRUE)
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
#' @param object An object of class \code{brmsfit}
#' @param correlation logical; if \code{FALSE} (the default), 
#'   compute the covariance matrix,
#'   if \code{TRUE}, compute the correlation matrix
#' @param ... Currently ignored
#' 
#' @return covariance or correlation matrix of population-level parameters
#' 
#' @details Estimates are obtained by calculating the maximum likelihood 
#'   covariances (correlations) of the posterior samples. 
#'   
#' @examples
#' \dontrun{
#' fit <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1+Trt_c|visit), 
#'            data = epilepsy, family = gaussian(), chains = 2)
#' vcov(fit)
#' }
#'
#' @export
vcov.brmsfit <- function(object, correlation = FALSE, ...) {
  contains_samples(object)
  pars <- parnames(object)
  fpars <- pars[grepl(fixef_pars(), pars)]
  if (!length(fpars)) {
    stop2("The model does not contain population-level effects.")
  }
  samples <- posterior_samples(object, pars = fpars, exact_match = TRUE)
  names(samples) <- sub(fixef_pars(), "", names(samples))
  if (correlation) {
    out <- cor(samples) 
  } else {
    out <- cov(samples)
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
#' @param object An object of class \code{brmsfit}.
#' @inheritParams fixef.brmsfit
#' @param ... Currently ignored.
#'
#' @return If \code{old} is \code{FALSE}: A list of arrays 
#'  (one per grouping factor). If \code{summary} is \code{TRUE},
#'  names of the first dimension are the factor levels and names
#'  of the third dimension are the group-level effects. 
#'  If \code{summary} is \code{FALSE}, names of the second dimension
#'  are the factor levels and names of the third dimension are the 
#'  group-level effects.
#'  
#' @author Paul-Christian Buerkner \email{paul.buerkner@gmail.com}   
#'   
#' @examples
#' \dontrun{
#' fit <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1+Trt_c|visit), 
#'            data = epilepsy, family = gaussian(), chains = 2)
#' ranef(fit)
#' }
#' 
#' @method ranef brmsfit
#' @export
#' @export ranef
#' @importFrom nlme ranef
ranef.brmsfit <- function(object, summary = TRUE, robust = FALSE,
                          probs = c(0.025, 0.975), ...) {
  contains_samples(object)
  object <- restructure(object)
  if (!nrow(object$ranef)) {
    stop2("The model does not contain group-level effects.")
  }
  pars <- parnames(object)
  ranef <- object$ranef
  groups <- unique(ranef$group)
  out <- named_list(groups)
  for (g in groups) {
    r <- subset2(ranef, group = g)
    coefs <- paste0(usc(combine_prefix(r), "suffix"), r$coef)
    levels <- attr(ranef, "levels")[[g]]
    rpars <- pars[grepl(paste0("^r_", g, "(__.+\\[|\\[)"), pars)]
    out[[g]] <- as.matrix(object, rpars, exact_match = TRUE)
    dim(out[[g]]) <- c(nrow(out[[g]]), length(levels), length(coefs))
    dimnames(out[[g]])[2:3] <- list(levels, coefs)
    if (summary) {
      out[[g]] <- posterior_summary(out[[g]], probs, robust)
    }
  }
  out
} 

#' Extract Model Coefficients
#'
#' Extract model coefficients, which are the sum of population-level 
#' effects and corresponding group-level effects
#' 
#' @param object An object of class \code{brmsfit}
#' @inheritParams ranef.brmsfit
#'
#' @return If \code{old} is \code{FALSE}: A list of arrays 
#'  (one per grouping factor). If \code{summary} is \code{TRUE},
#'  names of the first dimension are the factor levels and names
#'  of the third dimension are the group-level effects. 
#'  If \code{summary} is \code{FALSE}, names of the second dimension
#'  are the factor levels and names of the third dimension are the 
#'  group-level effects.
#'  
#' @author Paul-Christian Buerkner \email{paul.buerkner@gmail.com}
#'  
#' @examples
#' \dontrun{
#' fit <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1+Trt_c|visit), 
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
  contains_samples(object)
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
  fixef <- do.call(cbind, c(list(fixef), rmNULL(new_fixef)))
  
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
    coef[[g]] <- do.call(abind, c(list(coef[[g]]), rmNULL(new_ranef)))
    for (nm in dimnames(coef[[g]])[[3]]) {
      is_ord_intercept <- grepl("(^|_)Intercept\\[[[:digit:]]+\\]$", nm)
      if (is_ord_intercept) {
        # correct the sign of thresholds in ordinal models
        resp <- if (is_mv(object)) get_matches("^[^_]+", nm)
        family <- family(object, resp = resp)$family
        if (family %in% c("cumulative", "sratio")) {
          # threshold - mu
          coef[[g]][, , nm] <- fixef[, nm] - coef[[g]][, , nm] 
        } else {
          # mu - threshold
          coef[[g]][, , nm] <- coef[[g]][, , nm] - fixef[, nm]
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
#' @author Paul-Christian Buerkner \email{paul.buerkner@gmail.com}
#' 
#' @examples
#' \dontrun{
#' fit <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1+Trt_c|visit), 
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
  contains_samples(x)
  x <- restructure(x)
  if (!(nrow(x$ranef) || any(grepl("^sigma($|_)", parnames(x))))) {
    stop2("The model does not contain covariance matrices.")
  }
  .VarCorr <- function(y) {
    # extract samples for sd, cor and cov
    out <- list(sd = as.matrix(x, pars = y$sd_pars, exact_match = TRUE))
    colnames(out$sd) <- y$rnames
    # compute correlation and covariance matrices
    found_cor_pars <- intersect(y$cor_pars, parnames(x))
    if (length(found_cor_pars)) {
      cor <- as.matrix(x, pars = found_cor_pars, exact_match = TRUE)
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
  bterms <- parse_bf(x$formula)
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

#' @rdname posterior_samples
#' @export
posterior_samples.brmsfit <- function(x, pars = NA, exact_match = FALSE, 
                                      add_chain = FALSE, subset = NULL, 
                                      as.matrix = FALSE, as.array = FALSE,
                                      ...) {
  if (as.matrix && as.array) {
    stop2("Cannot use 'as.matrix' and 'as.array' at the same time.")
  }
  if (add_chain && as.array) {
    stop2("Cannot use 'add_chain' and 'as.array' at the same time.")
  }
  contains_samples(x)
  pars <- extract_pars(pars, parnames(x), exact_match = exact_match, ...)
  
  # get basic information on the samples 
  iter <- x$fit@sim$iter
  warmup <- x$fit@sim$warmup
  thin <- x$fit@sim$thin
  chains <- x$fit@sim$chains
  final_iter <- ceiling((iter - warmup) / thin)
  samples_taken <- seq(warmup + 1, iter, thin)
  
  if (length(pars)) {
    if (as.matrix) {
      samples <- as.matrix(x$fit, pars = pars)
    } else if (as.array) {
      samples <- as.array(x$fit, pars = pars)
    } else {
      samples <- as.data.frame(x$fit, pars = pars) 
    }
    if (add_chain) {
      # name the column 'chain' not 'chains' (#32)
      samples <- cbind(samples,
        chain = factor(rep(1:chains, each = final_iter)),
        iter = rep(samples_taken, chains)           
      )
    }
    if (!is.null(subset)) {
      if (as.array) {
        samples <- samples[subset, , , drop = FALSE] 
      } else {
        samples <- samples[subset, , drop = FALSE] 
      }
    }
  } else {
    samples <- NULL 
  }
  samples
}

#' @rdname posterior_samples
#' @export
as.data.frame.brmsfit <- function(x, row.names = NULL, optional = FALSE, ...) {
  out <- posterior_samples(x, ..., as.matrix = FALSE)
  data.frame(out, row.names = row.names, check.names = !optional)
}

#' @rdname posterior_samples
#' @export
as.matrix.brmsfit <- function(x, ...) {
  posterior_samples(x, ..., as.matrix = TRUE)
}

#' @rdname posterior_samples
#' @export
as.array.brmsfit <- function(x, ...) {
  posterior_samples(x, ..., as.array = TRUE)
}

#' Compute posterior uncertainty intervals 
#' 
#' Compute posterior uncertainty intervals for \code{brmsfit} objects.
#' 
#' @inheritParams summary.brmsfit
#' @param pars Names of parameters for which posterior samples should be 
#'   returned, as given by a character vector or regular expressions. 
#'   By default, all posterior samples of all parameters are extracted.
#' @param ... More arguments passed to 
#'   \code{\link[brms:as.matrix.brmsfit]{as.matrix.brmsfit}}.
#' 
#' @return A \code{matrix} with lower and upper interval bounds
#'   as columns and as many rows as selected parameters.
#'   
#' @examples 
#' \dontrun{
#' fit <- brm(count ~ log_Age_c + log_Base4_c * Trt_c,
#'            data = epilepsy, family = negbinomial())
#' posterior_interval(fit)
#' }
#' 
#' @aliases posterior_interval
#' @method posterior_interval brmsfit
#' @export
#' @export posterior_interval
#' @importFrom rstantools posterior_interval
posterior_interval.brmsfit <- function(
  object, pars = NA, prob = 0.95, ...
) {
  ps <- as.matrix(object, pars = pars, ...)
  rstantools::posterior_interval(ps, prob = prob)
}

#' @rdname posterior_summary
#' @export
posterior_summary.brmsfit <- function(x, pars = NA, 
                                      probs = c(0.025, 0.975), 
                                      robust = FALSE, ...) {
  out <- as.matrix(x, pars = pars, ...)
  posterior_summary(out, probs = probs, robust = robust, ...)
}

#' Extract posterior samples for use with the \pkg{coda} package
#' 
#' @aliases as.mcmc
#' 
#' @inheritParams posterior_samples
#' @param ... currently unused
#' @param combine_chains Indicates whether chains should be combined.
#' @param inc_warmup Indicates if the warmup samples should be included.
#'   Default is \code{FALSE}. Warmup samples are used to tune the 
#'   parameters of the sampling algorithm and should not be analyzed.
#'   
#' @return If \code{combine_chains = TRUE} an \code{mcmc} object is returned.
#'   If \code{combine_chains = FALSE} an \code{mcmc.list} object is returned.
#' 
#' @method as.mcmc brmsfit
#' @export
#' @export as.mcmc
#' @importFrom coda as.mcmc
as.mcmc.brmsfit <- function(x, pars = NA, exact_match = FALSE,
                            combine_chains = FALSE, inc_warmup = FALSE,
                            ...) {
  contains_samples(x)
  pars <- extract_pars(
    pars, all_pars = parnames(x),
    exact_match = exact_match, ...
  )
  if (combine_chains) {
    if (inc_warmup) {
      stop2("Cannot include warmup samples when 'combine_chains' is TRUE.")
    }
    out <- as.matrix(x$fit, pars)
    mcpar <- c(
      x$fit@sim$warmup * x$fit@sim$chain + 1, 
      x$fit@sim$iter * x$fit@sim$chain, x$fit@sim$thin
    )
    attr(out, "mcpar") <- mcpar
    class(out) <- "mcmc"
  } else {
    ps <- extract(x$fit, pars, permuted = FALSE, inc_warmup = inc_warmup)
    mcpar <- c(
      if (inc_warmup) 1 else x$fit@sim$warmup + 1, 
      x$fit@sim$iter, x$fit@sim$thin
    )
    out <- vector("list", length = dim(ps)[2])
    for (i in seq_along(out)) {
      out[[i]] <- ps[, i, ]
      attr(out[[i]], "mcpar") <- mcpar
      class(out[[i]]) <- "mcmc" 
    }
    class(out) <- "mcmc.list"
  }
  out
}

#' Extract Priors of a Bayesian Model Fitted with \pkg{brms}
#' 
#' @aliases prior_summary
#' 
#' @param object A \code{brmsfit} object
#' @param all Logical; Show all parameters in the model which may have 
#'   priors (\code{TRUE}) or only those with proper priors (\code{FALSE})?
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return For \code{brmsfit} objects, an object of class \code{brmsprior}.
#' 
#' @examples 
#' \dontrun{
#' fit <- brm(count ~ log_Age_c + log_Base4_c * Trt_c  
#'              + (1|patient) + (1|obs), 
#'            data = epilepsy, family = poisson(), 
#'            prior = c(prior(student_t(5,0,10), class = b),
#'                      prior(cauchy(0,2), class = sd)))
#'                    
#' prior_summary(fit)
#' prior_summary(fit, all = FALSE)
#' print(prior_summary(fit, all = FALSE), show_df = FALSE)
#' }
#' 
#' @method prior_summary brmsfit
#' @export
#' @export prior_summary
#' @importFrom rstantools prior_summary
prior_summary.brmsfit <- function(object, all = TRUE, ...) {
  object <- restructure(object)
  prior <- object$prior
  if (!all) {
    prior <- prior[nzchar(prior$prior), ]
  }
  prior
}

#' @rdname prior_samples
#' @export
prior_samples.brmsfit <- function(x, pars = NA, ...) {
  if (!anyNA(pars) && !is.character(pars)) {
    stop2("Argument 'pars' must be a character vector.")
  }
  par_names <- parnames(x)
  prior_names <- unique(par_names[grepl("^prior_", par_names)])
  if (length(prior_names)) {
    samples <- posterior_samples(x, prior_names, exact_match = TRUE)
    names(samples) <- sub("^prior_", "", prior_names)
    if (!anyNA(pars)) {
      .prior_samples <- function(par) {
        # get prior samples for parameter par 
        matches <- lapply(paste0("^", names(samples)), regexpr, text = par)
        matches <- ulapply(matches, attr, which = "match.length")
        if (max(matches) == -1) {
          out <- NULL
        } else {
          take <- match(max(matches), matches)
          # order samples randomly to avoid artifical dependencies
          # between parameters using the same prior samples
          samples <- list(samples[sample(nsamples(x)), take])
          out <- structure(samples, names = par)
        }
        return(out)
      }
      samples <- data.frame(
        rmNULL(lapply(pars, .prior_samples)), check.names = FALSE
      )
    }
  } else {
    samples <- NULL
  }
  samples
}

#' Print a summary for a fitted model represented by a \code{brmsfit} object
#' 
#' @aliases print.brmssummary
#' 
#' @param x An object of class \code{brmsfit}
#' @param digits The number of significant digits for printing out the summary; 
#'  defaults to 2. The effective sample size is always rounded to integers.
#' @param ... Additional arguments that would be passed 
#'  to method \code{summary} of \code{brmsfit}.
#'
#' @author Paul-Christian Buerkner \email{paul.buerkner@gmail.com}
#' 
#' @export
print.brmsfit <- function(x, digits = 2, ...) {
  print(summary(x, ...), digits = digits, ...)
}  

#' Create a summary of a fitted model represented by a \code{brmsfit} object
#'
#' @param object An object of class \code{brmsfit}
#' @param waic,loo Logical; Indicating if the LOO or WAIC information
#'   criteria should be computed and shown in the summary. 
#'   Defaults to \code{FALSE}.
#' @param R2 Logical; Indicating if the Bayesian R-squared
#'   should be computed and shown in the summary. 
#'   Defaults to \code{FALSE}.
#' @param priors Logical; Indicating if priors should be included 
#'   in the summary. Default is \code{FALSE}.
#' @param prob A value between 0 and 1 indicating the desired probability 
#'   to be covered by the uncertainty intervals. The default is 0.95.
#' @param mc_se Logical; Indicating if the uncertainty caused by the 
#'   MCMC sampling should be shown in the summary. Defaults to \code{FALSE}.
#' @param use_cache Logical; Indicating if summary results should
#'   be cached for future use by \pkg{rstan}. Defaults to \code{TRUE}.
#'   For models fitted with earlier versions of \pkg{brms},
#'   it may be necessary to set \code{use_cache} to
#'   \code{FALSE} in order to get the \code{summary} 
#'   method working correctly.
#' @param ... Other potential arguments
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@gmail.com}
#' 
#' @method summary brmsfit
#' @export
summary.brmsfit <- function(object, waic = FALSE, loo = FALSE, 
                            R2 = FALSE, priors = FALSE, prob = 0.95,
                            mc_se = FALSE, use_cache = TRUE, ...) {
  object <- restructure(object, rstr_summary = use_cache)
  bterms <- parse_bf(object$formula)
  out <- list(
    formula = object$formula,
    data.name = object$data.name, 
    group = unique(object$ranef$group), 
    nobs = nobs(object), 
    ngrps = ngrps(object), 
    autocor = object$autocor,
    prior = empty_brmsprior(),
    algorithm = algorithm(object),
    waic = NA, loo = NA, R2 = NA
  )
  class(out) <- "brmssummary"
  if (!length(object$fit@sim)) {
    # the model does not contain posterior samples
    return(out)
  }
  out$chains <- object$fit@sim$chains
  out$iter <- object$fit@sim$iter
  out$warmup <- object$fit@sim$warmup
  out$thin <- object$fit@sim$thin
  stan_args <- object$fit@stan_args[[1]]
  out$sampler <- paste0(stan_args$method, "(", stan_args$algorithm, ")")
  if (length(prob) != 1L || prob < 0 || prob > 1) {
    stop2("'prob' must be a single numeric value in [0, 1].")
  }
  if (priors) {
    out$prior <- prior_summary(object, all = FALSE)
  }
  if (loo || is.ic(object[["loo"]])) {
    out$loo <- SW(LOO(object)$looic)
  }
  if (waic || is.ic(object[["waic"]])) {
    out$waic <- SW(WAIC(object)$waic)
  }
  if (R2 || is.matrix(object[["R2"]])) {
    out$R2 <- mean(bayes_R2(object, summary = FALSE))
  }
  
  pars <- parnames(object)
  meta_pars <- object$fit@sim$pars_oi
  meta_pars <- meta_pars[!grepl("^(r|s|zgp|Xme|prior|lp)_", meta_pars)]
  probs <- c((1 - prob) / 2, 1 - (1 - prob) / 2)
  fit_summary <- summary(
    object$fit, pars = meta_pars, 
    probs = probs, use_cache = use_cache
  )
  fit_summary <- fit_summary$summary
  if (!mc_se) {
    fit_summary <- fit_summary[, -2, drop = FALSE] 
  }
  CIs <- paste0(c("l-", "u-"), prob * 100, "% CI")
  colnames(fit_summary) <- c(
    "Estimate", if (mc_se) "MC.Error", 
    "Est.Error", CIs, "Eff.Sample", "Rhat"
  )
  if (algorithm(object) == "sampling") {
    Rhats <- fit_summary[, "Rhat"]
    if (any(Rhats > 1.1, na.rm = TRUE) || anyNA(Rhats)) {
      warning2(
        "The model has not converged (some Rhats are > 1.1). ",
        "Do not analyse the results! \nWe recommend running ", 
        "more iterations and/or setting stronger priors."
      )
    }
    # nuts_params may not work for some models fitted with brms < 1.0.0
    div_trans <- try(
      sum(nuts_params(object, pars = "divergent__")$Value), silent = TRUE
    )
    if (is(div_trans, "try-error")) {
      warning2("Could not extract information about divergent transitions.")
    } else {
      adapt_delta <- control_params(object)$adapt_delta
      if (div_trans > 0) {
        warning2(
          "There were ", div_trans, " divergent transitions after warmup. ", 
          "Increasing adapt_delta above ", adapt_delta, " may help.\nSee ",
          "http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup"
        )
      }
    }
  }
  
  # summary of population-level effects
  fe_pars <- pars[grepl(fixef_pars(), pars)]
  out$fixed <- fit_summary[fe_pars, , drop = FALSE]
  rownames(out$fixed) <- gsub(fixef_pars(), "", fe_pars)
  
  # summary of family specific parameters
  spec_pars <- c(dpars(), "delta", "theta", "rescor")
  spec_pars <- paste0(spec_pars, collapse = "|")
  spec_pars <- paste0("^(", spec_pars, ")($|_|[[:digit:]])")
  spec_pars <- pars[grepl(spec_pars, pars)]
  out$spec_pars <- fit_summary[spec_pars, , drop = FALSE]
  is_rescor <- grepl("^rescor_", spec_pars)
  if (any(is_rescor)) {
    rescor_pars <- spec_pars[is_rescor]
    rescor_names <- sub("__", ",", sub("__", "(", rescor_pars))
    spec_pars[is_rescor] <- paste0(rescor_names, ")")
  }    
  rownames(out$spec_pars) <- spec_pars
  
  # summary of autocorrelation effects
  cor_pars <- pars[grepl(regex_cor_pars(), pars)]
  out$cor_pars <- fit_summary[cor_pars, , drop = FALSE]
  rownames(out$cor_pars) <- cor_pars
  
  # summary of group-level effects
  for (g in out$group) {
    gregex <- escape_dot(g)
    sd_prefix <- paste0("^sd_", gregex, "__")
    sd_pars <- pars[grepl(sd_prefix, pars)]
    cor_prefix <- paste0("^cor_", gregex, "__")
    cor_pars <- pars[grepl(cor_prefix, pars)]
    out$random[[g]] <- fit_summary[c(sd_pars, cor_pars), , drop = FALSE]
    if (nrow(out$random[[g]])) {
      sd_names <- sub(sd_prefix, "sd(", sd_pars)
      cor_names <- sub(cor_prefix, "cor(", cor_pars)
      cor_names <- sub("__", ",", cor_names)
      rownames(out$random[[g]]) <- paste0(c(sd_names, cor_names), ")")
    }
  }
  # summary of smooths
  sm_pars <- pars[grepl("^sds_", pars)]
  if (length(sm_pars)) {
    out$splines <- fit_summary[sm_pars, , drop = FALSE]
    rownames(out$splines) <- paste0(gsub("^sds_", "sds(", sm_pars), ")")
  }
  # summary of monotonic parameters
  mo_pars <- pars[grepl("^simo_", pars)]
  if (length(mo_pars)) {
    out$mo <- fit_summary[mo_pars, , drop = FALSE]
    rownames(out$mo) <- gsub("^simo_", "", mo_pars)
  }
  # summary of gaussian processes
  gp_pars <- pars[grepl("^(sdgp|lscale)_", pars)]
  if (length(gp_pars)) {
    out$gp <- fit_summary[gp_pars, , drop = FALSE]
    rownames(out$gp) <- gsub("^sdgp_", "sdgp(", rownames(out$gp))
    rownames(out$gp) <- gsub("^lscale_", "lscale(", rownames(out$gp))
    rownames(out$gp) <- paste0(rownames(out$gp), ")")
  }
  out
}

#' @rdname nsamples
#' @export
nsamples.brmsfit <- function(x, subset = NULL, 
                             incl_warmup = FALSE, ...) {
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) {
    out <- 0
  } else {
    ntsamples <- x$fit@sim$n_save[1]
    if (!incl_warmup) {
      ntsamples <- ntsamples - x$fit@sim$warmup2[1]
    }
    ntsamples <- ntsamples * x$fit@sim$chains
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
nobs.brmsfit <- function(object, ...) {
  nrow(model.frame(object))
}

#' @rdname ngrps
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

#' @export
formula.brmsfit <- function(x, ...) {
  x$formula
}

#' @export
getCall.brmsfit <- function(x, ...) {
  x$formula
}

#' @export
family.brmsfit <- function(object, resp = NULL, ...) {
  if (!is.null(resp)) {
    stopifnot(is_mv(object))
    resp <- as_one_character(resp)
    resp <- validate_resp(resp, object$formula$responses)
    family <- object$formula$forms[[resp]]$family
  } else {
    family <- object$family
  }
  family
}

#' @rdname stancode
#' @export
stancode.brmsfit <- function(object, version = TRUE, ...) {
  out <- object$model
  if (!version) {
    out <- sub("^[^\n]+[[:digit:]]\\.[^\n]+\n", "", out) 
  }
  out
}

#' @rdname standata
#' @export
standata.brmsfit <- function(object, newdata = NULL, re_formula = NULL, 
                             incl_autocor = TRUE, new_objects = list(),
                             internal = FALSE, control = list(), ...) {
  dots <- list(...)
  object <- restructure(object)
  if (!incl_autocor) {
    object <- remove_autocor(object)
  }
  if (internal) {
    control[c("not4stan", "save_order")] <- TRUE
  }
  new_formula <- update_re_terms(object$formula, re_formula)
  is_original_data <- isTRUE(attr(newdata, "original"))
  if (is.null(newdata)) {
    newdata <- object$data
  } else if (!is_original_data) {
    if (!isTRUE(attr(newdata, "valid"))) {
      newdata <- validate_newdata(
        newdata, object, re_formula = re_formula, ...
      )
    }
    object <- add_new_objects(object, newdata, new_objects)
    control$new <- TRUE
    # ensure correct handling of functions like poly or scale
    bterms <- parse_bf(new_formula)
    old_terms <- attr(object$data, "terms")
    terms_attr <- c("variables", "predvars")
    control$terms_attr <- attributes(old_terms)[terms_attr]
    control$old_sdata <- extract_old_standata(bterms, object$data)
    control$old_levels <- get_levels(
      tidy_ranef(bterms, object$data),
      tidy_meef(bterms, object$data)
    )
  }
  sample_prior <- attr(object$prior, "sample_prior")
  sample_prior <- ifelse(is.null(sample_prior), "no", sample_prior)
  args <- list(
    formula = new_formula, data = newdata, 
    prior = object$prior, cov_ranef = object$cov_ranef, 
    sample_prior = sample_prior, stan_vars = object$stan_vars, 
    knots = attr(object$data, "knots"), control = control
  )
  do.call(make_standata, c(args, dots))
}

#' Interface to \pkg{shinystan}
#' 
#' Provide an interface to \pkg{shinystan} for models fitted with \pkg{brms}
#' 
#' @aliases launch_shinystan
#' 
#' @param object A fitted model object typically of class \code{brmsfit}. 
#' @param rstudio Only relevant for RStudio users. 
#' The default (\code{rstudio=FALSE}) is to launch the app 
#' in the default web browser rather than RStudio's pop-up Viewer. 
#' Users can change the default to \code{TRUE} 
#' by setting the global option \cr \code{options(shinystan.rstudio = TRUE)}.
#' @param ... Optional arguments to pass to \code{\link[shiny:runApp]{runApp}}
#' 
#' @return An S4 shinystan object
#' 
#' @examples
#' \dontrun{
#' fit <- brm(rating ~ treat + period + carry + (1|subject),
#'            data = inhaler, family = "gaussian")
#' launch_shinystan(fit)                         
#' }
#' 
#' @seealso \code{\link[shinystan:launch_shinystan]{launch_shinystan}}
#' 
#' @method launch_shinystan brmsfit
#' @importFrom shinystan launch_shinystan
#' @export launch_shinystan
#' @export
launch_shinystan.brmsfit <- function(
  object, rstudio = getOption("shinystan.rstudio"), ...
) {
  contains_samples(object)
  launch_shinystan(object$fit, rstudio = rstudio, ...)
}

#' Trace and Density Plots for MCMC Samples
#' 
#' @param x An object of class \code{brmsfit}.
#' @param pars Names of the parameters to plot, as given by a character vector 
#'   or a regular expression. By default, all parameters except 
#'   for group-level and smooth effects are plotted. 
#' @param combo A character vector with at least two elements. 
#'   Each element of \code{combo} corresponds to a column in the resulting 
#'   graphic and should be the name of one of the available 
#'   \code{\link[bayesplot:MCMC-overview]{MCMC}} functions 
#'   (omitting the \code{mcmc_} prefix).
#' @param N The number of parameters plotted per page.
#' @param theme A \code{\link[ggplot2:theme]{theme}} object 
#'   modifying the appearance of the plots. 
#'   For some basic themes see \code{\link[ggplot2:ggtheme]{ggtheme}}
#'   and \code{\link[bayesplot:theme_default]{theme_default}}.
#' @param exact_match Indicates whether parameter names 
#'   should be matched exactly or treated as regular expression. 
#'   Default is \code{FALSE}.
#' @param plot logical; indicates if plots should be
#'   plotted directly in the active graphic device.
#'   Defaults to \code{TRUE}.
#' @param ask logical; indicates if the user is prompted 
#'   before a new page is plotted. 
#'   Only used if \code{plot} is \code{TRUE}.
#' @param newpage logical; indicates if the first set of plots
#'   should be plotted to a new page. 
#'   Only used if \code{plot} is \code{TRUE}.
#' @param ... Further arguments passed to 
#'   \code{\link[bayesplot:mcmc_combo]{mcmc_combo}}.
#' 
#' @return An invisible list of 
#'   \code{\link[gtable:gtable]{gtable}} objects.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@gmail.com}
#' 
#' @examples
#' \dontrun{ 
#' fit <- brm(count ~ log_Age_c + log_Base4_c * Trt_c 
#'            + (1|patient) + (1|visit), 
#'            data = epilepsy, family = "poisson")
#' plot(fit)
#' ## plot population-level effects only
#' plot(fit, pars = "^b_") 
#' }
#' 
#' @method plot brmsfit
#' @import ggplot2
#' @importFrom grDevices devAskNewPage
#' @export
plot.brmsfit <- function(x, pars = NA, combo = c("dens", "trace"), 
                         N = 5, exact_match = FALSE, theme = NULL, 
                         plot = TRUE, ask = TRUE, 
                         newpage = TRUE, ...) {
  contains_samples(x)
  if (!is_wholenumber(N) || N < 1) {
    stop2("Argument 'N' must be a positive integer.")
  }
  if (!is.character(pars)) {
    pars <- default_plot_pars()
    exact_match <- FALSE
  }
  samples <- posterior_samples(
    x, pars = pars, add_chain = TRUE, exact_match = exact_match
  )
  pars <- names(samples)[!names(samples) %in% c("chain", "iter")] 
  if (!length(pars)) {
    stop2("No valid parameters selected.")
  }
  
  if (plot) {
    default_ask <- devAskNewPage()
    on.exit(devAskNewPage(default_ask))
    devAskNewPage(ask = FALSE)
  }
  n_plots <- ceiling(length(pars) / N)
  plots <- vector(mode = "list", length = n_plots)
  for (i in seq_len(n_plots)) {
    sub_pars <- pars[((i - 1) * N + 1):min(i * N, length(pars))]
    sub_samples <- samples[, c(sub_pars, "chain"), drop = FALSE]
    plots[[i]] <- bayesplot::mcmc_combo(
      sub_samples, combo = combo, gg_theme = theme, ...
    )
    if (plot) {
      plot(plots[[i]], newpage = newpage || i > 1)
      if (i == 1) {
        devAskNewPage(ask = ask)
      }
    }
  }
  invisible(plots) 
}

#' @rdname stanplot
#' @export
stanplot.brmsfit <- function(object, pars = NA, type = "intervals", 
                             exact_match = FALSE, ...) {
  contains_samples(object)
  object <- restructure(object)
  if (length(type) != 1L) {
    stop2("Argument 'type' must be of length 1.")
  }
  if (!is.character(pars)) {
    pars <- default_plot_pars()
    exact_match <- FALSE
  }
  nuts_types <- c(
    "acceptance", "divergence", "stepsize", "treedepth", "energy"
  )
  valid_types <- c(
    "hist", "dens", "hist_by_chain", "dens_overlay", 
    "violin", "intervals", "areas", "trace", "trace_highlight", 
    "scatter", "hex", "pairs", "rhat", "rhat_hist", "neff", 
    "neff_hist", "acf", "acf_bar", "recover_intervals",
    paste0("nuts_", nuts_types)
  )
  if (!type %in% valid_types) {
    stop2("Invalid plot type. Valid plot types are: \n",
          collapse_comma(valid_types))
  }
  mcmc_fun <- get(paste0("mcmc_", type), pos = asNamespace("bayesplot"))
  mcmc_arg_names <- names(formals(mcmc_fun))
  mcmc_args <- list(...)
  if ("x" %in% mcmc_arg_names) {
    if (grepl("^nuts_", type)) {
      # x refers to a molten data.frame of NUTS parameters
      mcmc_args[["x"]] <- nuts_params(object)
    } else {
      # x refers to a data.frame of samples
      samples <- posterior_samples(
        object, pars, add_chain = TRUE, exact_match = exact_match
      )
      samples$iter <- NULL
      sel_pars <- names(samples)[!names(samples) %in% "chain"]
      if (type %in% c("scatter", "hex") && length(sel_pars) != 2L) {
        stop2("Exactly 2 parameters must be selected for this type.",
              "\nParameters selected: ", collapse_comma(sel_pars))
      }
      mcmc_args[["x"]] <- samples
    }
  }
  if ("lp" %in% mcmc_arg_names) {
    mcmc_args[["lp"]] <- log_posterior(object)
  }
  use_nuts <- isTRUE(object$algorithm == "sampling")
  if ("np" %in% mcmc_arg_names && use_nuts) {
    mcmc_args[["np"]] <- nuts_params(object)
  }
  interval_type <- type %in% c("intervals", "areas")
  if ("rhat" %in% mcmc_arg_names && !interval_type) {
    mcmc_args[["rhat"]] <- rhat(object)
  }
  if ("ratio" %in% mcmc_arg_names) {
    mcmc_args[["ratio"]] <- neff_ratio(object)
  }
  do.call(mcmc_fun, mcmc_args)
}

#' Posterior Predictive Checks for \code{brmsfit} Objects
#' 
#' Perform posterior predictive checks with the help
#' of the \pkg{bayesplot} package.
#' 
#' @aliases pp_check
#' 
#' @param object An object of class \code{brmsfit}.
#' @param type Type of the ppc plot as given by a character string.
#'   See \code{\link[bayesplot:PPC-overview]{PPC}} for an overview
#'   of currently supported types. You may also use an invalid
#'   type (e.g. \code{type = "xyz"}) to get a list of supported 
#'   types in the resulting error message.
#' @param nsamples Positive integer indicating how many
#'  posterior samples should be used.
#'  If \code{NULL} all samples are used. If not specified,
#'  the number of posterior samples is chosen automatically.
#'  Ignored if \code{subset} is not \code{NULL}.
#' @param ntrys Parameter used in rejection sampling
#'  for truncated discrete models only
#'  (defaults to \code{5}). For more details see
#'  \code{\link[brms:predict.brmsfit]{predict.brmsfit}}.
#' @param group Optional name of a factor variable in the model
#'  by which to stratify the ppc plot. This argument is required for
#'  ppc \code{*_grouped} types and ignored otherwise.
#' @param x Optional name of a variable in the model. 
#'  Only used for ppc types having an \code{x} argument 
#'  and ignored otherwise.
#' @param loo_args An optional list of additional arguments 
#'  passed to \code{\link[loo:psislw]{psislw}}. 
#'  Ignored for non \code{loo_*} ppc types.
#' @param ... Further arguments passed to the ppc functions
#'   of \pkg{\link[bayesplot:bayesplot]{bayesplot}}.
#' @inheritParams predict.brmsfit
#' 
#' @return A ggplot object that can be further
#'  customized using the \pkg{ggplot2} package.
#' 
#' @details For a detailed explanation of each of the ppc functions, 
#' see the \code{\link[bayesplot:PPC-overview]{PPC}} 
#' documentation of the \pkg{\link[bayesplot:bayesplot]{bayesplot}} 
#' package.
#' 
#' @examples
#' \dontrun{
#' fit <-  brm(count ~ log_Age_c + log_Base4_c * Trt_c
#'             + (1|patient) + (1|obs),
#'             data = epilepsy, family = poisson())
#' 
#' pp_check(fit)  # shows dens_overlay plot by default
#' pp_check(fit, type = "error_hist", nsamples = 11)
#' pp_check(fit, type = "scatter_avg", nsamples = 100)
#' pp_check(fit, type = "stat_2d")
#' pp_check(fit, type = "rootogram")
#' pp_check(fit, type = "loo_pit")
#' 
#' ## get an overview of all valid types
#' pp_check(fit, type = "xyz")
#' }
#' 
#' @importFrom bayesplot pp_check
#' @export pp_check
#' @export
pp_check.brmsfit <- function(object, type, nsamples, group = NULL,
                             x = NULL, resp = NULL, newdata = NULL,
                             re_formula = NULL, allow_new_levels = FALSE,
                             sample_new_levels = "uncertainty", 
                             new_objects = list(), incl_autocor = TRUE, 
                             subset = NULL, nug = NULL, ntrys = 5, 
                             loo_args = list(), ...) {
  if (missing(type)) {
    type <- "dens_overlay"
  }
  type <- as_one_character(type)
  if (!is.null(group)) {
    group <- as_one_character(group)
  }
  if (!is.null(x)) {
    x <- as_one_character(x)
  }
  if (!is.null(resp)) {
    resp <- as_one_character(resp)
  }
  ppc_funs <- as.character(bayesplot::available_ppc(""))
  valid_ppc_types <- sub("^ppc_", "", ppc_funs)
  if (!type %in% valid_ppc_types) {
    stop2("Type '", type, "' is not a valid ppc type. Valid types are: \n", 
          collapse_comma(valid_ppc_types))
  }
  ppc_fun <- get(paste0("ppc_", type), pos = asNamespace("bayesplot"))
  # validate arguments 'resp', 'group', and 'x'
  object <- restructure(object)
  stopifnot_resp(object, resp)
  valid_groups <- get_cat_vars(object)
  if (!is.null(group) && !group %in% valid_groups) {
    stop2("Group '", group, "' is not a valid grouping factor. ",
          "Valid groups are: \n", collapse_comma(valid_groups))
  }
  is_group_type <- "group" %in% names(formals(ppc_fun))
  if (is.null(group) && is_group_type) {
    stop2("Argument 'group' is required for ppc type '", type, "'.")
  }
  bterms <- parse_bf(object$formula)
  if ("x" %in% names(formals(ppc_fun)) && !is.null(x)) {
    ae_coll <- ulapply(get_all_effects(bterms), paste, collapse = ":")
    if (!x %in% ae_coll) {
      stop2("Variable '", x, "' is not a valid variable for this model.",
            "\nValid variables are: ", collapse_comma(ae_coll))
    }
  }
  if (type == "error_binned") {
    if (is_ordinal(object$family)) {
      stop2("Type '", type, "' is not available for ordinal models.")
    }
    method <- "fitted"
  } else {
    method <- "predict"
  }
  if (missing(nsamples)) {
    aps_types <- c(
      "error_scatter_avg", "error_scatter_avg_vs_x",
      "intervals", "intervals_grouped", "loo_pit", 
      "loo_intervals", "loo_ribbon", "ribbon", 
      "ribbon_grouped", "rootogram", "scatter_avg", 
      "scatter_avg_grouped", "stat", "stat_2d", 
      "stat_freqpoly_grouped", "stat_grouped", 
      "violin_grouped"
    )
    if (!is.null(subset)) {
      nsamples <- NULL
    } else if (type %in% aps_types) {
      nsamples <- NULL
      message("Using all posterior samples for ppc type '", 
              type, "' by default.")
    } else {
      nsamples <- 10
      message("Using 10 posterior samples for ppc type '",
              type, "' by default.")
    }
  }
  newd_args <- nlist(
    object, newdata, re_formula = NA, incl_autocor, 
    new_objects, check_response = TRUE
  )
  newdata <- do.call(validate_newdata, newd_args)
  newd_args$newdata <- newdata
  newd_args[c("internal", "only_response")] <- TRUE
  sdata <- do.call(standata, newd_args)
  if (any(grepl("^cens_", names(sdata)))) {
    warning2("'pp_check' may not be meaningful for censored models.")
  }
  y <- as.vector(sdata[[paste0("Y", usc(resp))]])
  pred_args <- nlist(
    object, newdata, re_formula, allow_new_levels, resp,
    sample_new_levels, new_objects, incl_autocor, nsamples, 
    subset, nug, ntrys, sort = FALSE, summary = FALSE
  )
  yrep <- do.call(method, pred_args)
  if (has_trials(family(object, resp = resp))) {
    # use success proportions following Gelman and Hill (2006)
    y <- y / sdata$trials
    yrep <- yrep / as_draws_matrix(sdata$trials, dim = dim(yrep))
  }
  ppc_args <- list(y, yrep, ...)
  if ("lw" %in% names(formals(ppc_fun)) && !"lw" %in% names(ppc_args)) {
    # required for 'loo' types only
    pred_args$loo_args <- loo_args
    ppc_args$lw <- do.call(loo_weights, c(pred_args, log = TRUE))
  }
  # allow using arguments 'group' and 'x' for new data
  mf <- update_data(newdata, bterms, na.action = na.pass)
  if (!is.null(group)) {
    ppc_args$group <- mf[[group]]
  }
  if (!is.null(x)) {
    ppc_args$x <- mf[[x]]
    if (!is_like_factor(ppc_args$x)) {
      # fixes issue #229
      ppc_args$x <- as.numeric(ppc_args$x)
    }
  }
  do.call(ppc_fun, ppc_args)
}

#' Create a matrix of output plots from a \code{brmsfit} object
#'
#' A \code{\link[graphics:pairs]{pairs}} 
#' method that is customized for MCMC output.
#' 
#' @param x An object of class \code{brmsfit}
#' @inheritParams posterior_samples
#' @param ... Further arguments to be passed to 
#'   \code{\link[bayesplot:mcmc_pairs]{mcmc_pairs}}.
#'  
#' @details For a detailed description see  
#'   \code{\link[bayesplot:mcmc_pairs]{mcmc_pairs}}.
#'  
#' @examples 
#' \dontrun{
#' fit <- brm(count ~ log_Age_c + log_Base4_c * Trt_c 
#'            + (1|patient) + (1|visit), 
#'            data = epilepsy, family = "poisson")  
#' pairs(fit, pars = parnames(fit)[1:3], exact_match = TRUE)
#' pairs(fit, pars = "^sd_")
#' }
#'
#' @export
pairs.brmsfit <- function(x, pars = NA, exact_match = FALSE, ...) {
  if (!is.character(pars)) {
    pars <- default_plot_pars()
    exact_match <- FALSE
  }
  samples <- posterior_samples(
    x, pars, add_chain = TRUE, exact_match = exact_match
  )
  samples$iter <- NULL
  bayesplot::mcmc_pairs(samples, ...)
}

#' @rdname marginal_effects
#' @export
marginal_effects.brmsfit <- function(x, effects = NULL, conditions = NULL, 
                                     int_conditions = NULL, re_formula = NA, 
                                     robust = TRUE, probs = c(0.025, 0.975),
                                     method = c("fitted", "predict"), 
                                     spaghetti = FALSE, surface = FALSE,
                                     ordinal = FALSE, transform = NULL, 
                                     resolution = 100, select_points = 0, 
                                     too_far = 0, ...) {
  method <- match.arg(method)
  spaghetti <- as_one_logical(spaghetti)
  surface <- as_one_logical(surface)
  ordinal <- as_one_logical(ordinal)
  contains_samples(x)
  x <- restructure(x)
  new_formula <- update_re_terms(x$formula, re_formula = re_formula)
  bterms <- parse_bf(new_formula)
  if (length(probs) != 2L) {
    stop2("Arguments 'probs' must be of length 2.")
  }
  if (!is.null(transform) && method != "predict") {
    stop2("'transform' is only allowed when 'method' is set to 'predict'.")
  }
  if (any(is_ordinal(family_names(x))) && !ordinal) {
    warning2(
      "Predictions are treated as continuous variables in ",
      "'marginal_effects' by default, which is likely invalid ", 
      "for ordinal families. Consider setting 'ordinal' to TRUE."
    )
  }
  rsv_vars <- rsv_vars(bterms)
  use_def_effects <- is.null(effects)
  if (use_def_effects) {
    effects <- get_all_effects(bterms, rsv_vars = rsv_vars)
  } else {
    # allow to define interactions in any order
    effects <- strsplit(as.character(effects), split = ":")
    if (any(unique(unlist(effects)) %in% rsv_vars)) {
      stop2("Variables ", collapse_comma(rsv_vars),
            " should not be used as effects for this model")
    }
    if (any(lengths(effects) > 2L)) {
      stop2("To display interactions of order higher than 2 ",
            "please use the 'conditions' argument.")
    }
    all_effects <- get_all_effects(
      bterms, rsv_vars = rsv_vars, comb_all = TRUE
    )
    ae_coll <- all_effects[lengths(all_effects) == 1L]
    ae_coll <- ulapply(ae_coll, paste, collapse = ":")
    matches <- match(lapply(all_effects, sort), lapply(effects, sort), 0L)
    if (sum(matches) > 0 && sum(matches > 0) < length(effects)) {
      invalid <- effects[setdiff(seq_along(effects), sort(matches))]  
      invalid <- ulapply(invalid, paste, collapse = ":")
      warning2(
        "Some specified effects are invalid for this model: ",
        collapse_comma(invalid), "\nValid effects are ", 
        "(combinations of): ", collapse_comma(ae_coll)
      )
    }
    effects <- unique(effects[sort(matches)])
    if (!length(effects)) {
      stop2(
        "All specified effects are invalid for this model.\n", 
        "Valid effects are (combinations of): ", 
        collapse_comma(ae_coll)
      )
    }
  }
  if (ordinal) {
    int_effs <- lengths(effects) == 2L
    if (any(int_effs)) {
      effects <- effects[!int_effs]
      warning2("Interactions cannot be plotted if 'ordinal' is TRUE.")
    }
  }
  if (!length(effects)) {
    stop2("No valid effects detected.")
  }
  mf <- model.frame(x)
  conditions <- prepare_conditions(
    x, conditions = conditions, effects = effects, 
    re_formula = re_formula, rsv_vars = rsv_vars
  )
  int_vars <- get_int_vars(bterms)
  int_conditions <- lapply(int_conditions, 
    function(x) if (is.numeric(x)) sort(x, TRUE) else x
  )
  out <- list()
  for (i in seq_along(effects)) {
    eff <- effects[[i]]
    marg_args <- nlist(
      data = mf[, eff, drop = FALSE], conditions, int_conditions,
      int_vars, surface, resolution, reorder = use_def_effects
    )
    marg_data <- do.call(prepare_marg_data, marg_args)
    if (surface && length(eff) == 2L && too_far > 0) {
      # exclude prediction grid points too far from data
      ex_too_far <- mgcv::exclude.too.far(
        g1 = marg_data[[eff[1]]], 
        g2 = marg_data[[eff[2]]], 
        d1 = mf[, eff[1]],
        d2 = mf[, eff[2]],
        dist = too_far)
      marg_data <- marg_data[!ex_too_far, ]  
    }
    me_args <- nlist(
      x = bterms, fit = x, marg_data, method, surface, 
      spaghetti, ordinal, re_formula, transform, conditions,
      int_conditions, select_points, probs, robust, ...
    )
    out <- c(out, do.call(marginal_effects_internal, me_args))
  }
  structure(out, class = "brmsMarginalEffects")
}

#' @rdname marginal_smooths
#' @export
marginal_smooths.brmsfit <- function(x, smooths = NULL,
                                     int_conditions = NULL,
                                     probs = c(0.025, 0.975),
                                     spaghetti = FALSE,
                                     resolution = 100, too_far = 0,
                                     subset = NULL, nsamples = NULL,
                                     ...) {
  spaghetti <- as_one_logical(spaghetti)
  contains_samples(x)
  x <- restructure(x)
  x <- remove_autocor(x)
  smooths <- rm_wsp(as.character(smooths))
  bterms <- parse_bf(x$formula)
  conditions <- prepare_conditions(x)
  subset <- subset_samples(x, subset, nsamples)
  # call as.matrix only once to save time and memory
  samples <- as.matrix(x, subset = subset)
  ms_args <- nlist(
    x = bterms, fit = x, samples, smooths, conditions, 
    int_conditions, too_far, resolution, probs, spaghetti
  )
  out <- do.call(marginal_smooths_internal, ms_args)
  if (!length(out)) {
    stop2("No valid smooth terms found in the model.")
  }
  too_many_covars <- any(ulapply(out, attr, "too_many_covars"))
  if (too_many_covars) {
    warning2("Smooth terms with more than two covariates ",
             "are not yet supported by 'marginal_smooths'.")
  }
  structure(out, class = "brmsMarginalEffects", smooths_only = TRUE)
}

#' Model Predictions of \code{brmsfit} Objects
#' 
#' Predict responses based on the fitted model.
#' Can be performed for the data used to fit the model 
#' (posterior predictive checks) or for new data.
#' By definition, these predictions have higher variance than 
#' predictions of the fitted values (i.e., the 'regression line')
#' performed by the \code{\link[brms:fitted.brmsfit]{fitted}}
#' method. This is because the measurement error is incorporated.
#' The estimated means of both methods should, however, be very similar.
#' 
#' @param object An object of class \code{brmsfit}
#' @param newdata An optional data.frame for which to evaluate predictions.
#'   If \code{NULL} (default), the original data of the model is used.
#' @param re_formula formula containing group-level effects 
#'   to be considered in the prediction. 
#'   If \code{NULL} (default), include all group-level effects; 
#'   if \code{NA}, include no group-level effects.
#' @param re.form Alias of \code{re_formula}.
#' @param transform A function or a character string naming 
#'   a function to be applied on the predicted responses
#'   before summary statistics are computed.
#' @param allow_new_levels A flag indicating if new
#'   levels of group-level effects are allowed 
#'   (defaults to \code{FALSE}). 
#'   Only relevant if \code{newdata} is provided.
#' @param sample_new_levels Indicates how to sample new levels 
#'   for grouping factors specified in \code{re_formula}.
#'   This argument is only relevant if \code{newdata} is provided and 
#'   \code{allow_new_levels} is set to \code{TRUE}.
#'   If \code{"uncertainty"} (default), include group-level uncertainty
#'   in the predictions based on the variation of the existing levels. 
#'   If \code{"gaussian"}, sample new levels from the (multivariate) 
#'   normal distribution implied by the group-level standard deviations 
#'   and correlations. This options may be useful for conducting 
#'   Bayesian power analysis. 
#'   If \code{"old_levels"}, directly sample new levels from the
#'   existing levels. 
#' @param new_objects A named \code{list} of objects containing 
#'   new data, which cannot be passed via argument \code{newdata}.
#'   Currently, only required for objects passed to 
#'   \code{\link[brms:cor_sar]{cor_sar}} and 
#'   \code{\link[brms:cor_fixed]{cor_fixed}}. 
#' @param incl_autocor A flag indicating if ARMA autocorrelation
#'  parameters should be included in the predictions. Defaults to 
#'  \code{TRUE}. Setting it to \code{FALSE} will not affect other 
#'  correlation structures such as \code{\link[brms:cor_bsts]{cor_bsts}},
#'  or \code{\link[brms:cor_fixed]{cor_fixed}}.
#' @param negative_rt Only relevant for Wiener diffusion models. 
#'   A flag indicating whether response times of responses
#'   on the lower boundary should be returned as negative values.
#'   This allows to distinguish responses on the upper and
#'   lower boundary. Defaults to \code{FALSE}.
#' @param resp Optional names of response variables.
#'  If specified, fitted values of these response variables are returned.
#' @param summary Should summary statistics 
#'   (i.e. means, sds, and 95\% intervals) be returned
#'  instead of the raw values? Default is \code{TRUE}.
#' @param robust If \code{FALSE} (the default) the mean is used as 
#'  the measure of central tendency and the standard deviation as 
#'  the measure of variability. If \code{TRUE}, the median and the 
#'  median absolute deviation (MAD) are applied instead.
#'  Only used if \code{summary} is \code{TRUE}.
#' @param probs  The percentiles to be computed by the \code{quantile} 
#'  function. Only used if \code{summary} is \code{TRUE}. 
#' @param subset A numeric vector specifying
#'  the posterior samples to be used. 
#'  If \code{NULL} (the default), all samples are used.
#' @param nsamples Positive integer indicating how many 
#'  posterior samples should be used. 
#'  If \code{NULL} (the default) all samples are used.
#'  Ignored if \code{subset} is not \code{NULL}.
#' @param sort Logical. Only relevant for time series models. 
#'  Indicating whether to return predicted values in the original 
#'  order (\code{FALSE}; default) or in the order of the 
#'  time series (\code{TRUE}). 
#' @param ntrys Parameter used in rejection sampling 
#'   for truncated discrete models only 
#'   (defaults to \code{5}). See Details for more information.
#' @param nug Small positive number for Gaussian process terms only. 
#'   For numerical reasons, the covariance matrix of a Gaussian 
#'   process might not be positive definite. Adding a very small 
#'   number to the matrix's diagonal often solves this problem. 
#'   If \code{NULL} (the default), \code{nug} is chosen internally.
#' @param ... Currently ignored.
#' 
#' @return Predicted values of the response variable. 
#'   If \code{summary = TRUE} the output depends on the family:
#'   For categorical and ordinal families, it is a N x C matrix, 
#'   where N is the number of observations and
#'   C is the number of categories. 
#'   For all other families, it is a N x E matrix where E is equal 
#'   to \code{length(probs) + 2}.
#'   If \code{summary = FALSE}, the output is as a S x N matrix, 
#'   where S is the number of samples.
#'   In multivariate models, the output is an array of 3 dimensions, 
#'   where the third dimension indicates the response variables.
#' 
#' @details \code{NA} values within factors in \code{newdata}, 
#'   are interpreted as if all dummy variables of this factor are 
#'   zero. This allows, for instance, to make predictions of the grand mean 
#'   when using sum coding.  
#'   
#'   Method \code{posterior_predict.brmsfit} is an alias of 
#'   \code{predict.brmsfit} with \code{summary = FALSE}. 
#' 
#'   For truncated discrete models only:
#'   In the absence of any general algorithm to sample 
#'   from truncated discrete distributions,
#'   rejection sampling is applied in this special case. 
#'   This means that values are sampled until 
#'   a value lies within the defined truncation boundaries. 
#'   In practice, this procedure may be rather slow (especially in \R). 
#'   Thus, we try to do approximate rejection sampling 
#'   by sampling each value \code{ntrys} times and then select a valid value. 
#'   If all values are invalid, the closest boundary is used, instead. 
#'   If there are more than a few of these pathological cases, 
#'   a warning will occur suggesting to increase argument \code{ntrys}.
#' 
#' @examples 
#' \dontrun{
#' ## fit a model
#' fit <- brm(time | cens(censored) ~ age + sex + (1+age||patient), 
#'            data = kidney, family = "exponential", inits = "0")
#' 
#' ## predicted responses
#' pp <- predict(fit)
#' head(pp)
#' 
#' ## predicted responses excluding the group-level effect of age
#' pp2 <- predict(fit, re_formula = ~ (1|patient))
#' head(pp2)
#' 
#' ## predicted responses of patient 1 for new data
#' newdata <- data.frame(sex = factor(c("male", "female")),
#'                       age = c(20, 50),
#'                       patient = c(1, 1))
#' predict(fit, newdata = newdata)
#' }
#' 
#' @export 
predict.brmsfit <- function(object, newdata = NULL, re_formula = NULL,
                            transform = NULL, allow_new_levels = FALSE,
                            sample_new_levels = "uncertainty", 
                            new_objects = list(), incl_autocor = TRUE, 
                            resp = NULL, negative_rt = FALSE, subset = NULL, 
                            nsamples = NULL, sort = FALSE, nug = NULL, 
                            ntrys = 5, summary = TRUE, robust = FALSE,
                            probs = c(0.025, 0.975), ...) {
  contains_samples(object)
  object <- restructure(object)
  draws_args <- nlist(
    x = object, newdata, re_formula, incl_autocor, allow_new_levels,
    sample_new_levels, new_objects, resp, subset, nsamples, nug,
    check_response = FALSE
  )
  draws <- do.call(extract_draws, draws_args)
  predict_args <- nlist(
    draws, transform, ntrys, negative_rt, 
    summary, robust, probs, sort
  )
  do.call(predict_internal, predict_args)
}

#' @rdname predict.brmsfit
#' @aliases posterior_predict
#' @method posterior_predict brmsfit
#' @export
#' @export posterior_predict
#' @importFrom rstantools posterior_predict
posterior_predict.brmsfit <- function(object, newdata = NULL, re_formula = NULL,
                                      re.form = NULL, transform = NULL, 
                                      allow_new_levels = FALSE,
                                      sample_new_levels = "uncertainty", 
                                      new_objects = list(), incl_autocor = TRUE, 
                                      resp = NULL, negative_rt = FALSE, subset = NULL, 
                                      nsamples = NULL, sort = FALSE, nug = NULL, 
                                      ntrys = 5, robust = FALSE,
                                      probs = c(0.025, 0.975), ...) {
  cl <- match.call()
  cl[[1]] <- quote(predict)
  cl$summary <- FALSE
  if ("re.form" %in% names(cl)) {
    cl$re_formula <- cl$re.form
    cl$re.form <- NULL
  }
  eval(cl, parent.frame())
}

#' Extract Model Fitted Values of \code{brmsfit} Objects
#' 
#' Predict fitted values (i.e., the 'regression line') of a fitted model.
#' Can be performed for the data used to fit the model 
#' (posterior predictive checks) or for new data.
#' By definition, these predictions have smaller variance
#' than the response predictions performed by
#' the \code{\link[brms:predict.brmsfit]{predict}} method. 
#' This is because the measurement error is not incorporated.
#' The estimated means of both methods should, however, be very similar.
#' 
#' @inheritParams predict.brmsfit
#' @param scale Either \code{"response"} or \code{"linear"}. 
#'  If \code{"response"} results are returned on the scale 
#'  of the response variable. If \code{"linear"} 
#'  fitted values are returned on the scale of the linear predictor.
#' @param transform Logical; alias of \code{scale}.
#'  If \code{TRUE}, \code{scale} is set to \code{"response"}.
#'  If \code{FALSE}, \code{scale} is set to \code{"linear"}. 
#'  Only implemented for compatibility with the 
#'  \code{\link[rstantools:posterior_linpred]{posterior_linpred}}
#'  generic. 
#' @param dpar Optional name of a predicted distributional parameter.
#'  If specified, fitted values of this parameters are returned.
#'
#' @return Fitted values extracted from \code{object}. 
#'  The output depends on the family:
#'  If \code{summary = TRUE} it is a N x E x C array 
#'  for categorical and ordinal models and a N x E matrix else.
#'  If \code{summary = FALSE} it is a S x N x C array 
#'  for categorical and ordinal models and a S x N matrix else.
#'  N is the number of observations, S is the number of samples, 
#'  C is the number of categories, and E is equal to \code{length(probs) + 2}.
#'  In multivariate models, the output is an array of 3 dimensions, 
#'  where the third dimension indicates the response variables.
#'   
#' @details \code{NA} values within factors in \code{newdata}, 
#'   are interpreted as if all dummy variables of this factor are 
#'   zero. This allows, for instance, to make predictions of the grand mean 
#'   when using sum coding.
#'   
#'   Method \code{posterior_linpred.brmsfit} is an alias of 
#'   \code{fitted.brmsfit} with \code{scale = "linear"} and
#'   \code{summary = FALSE}. 
#'
#' @examples 
#' \dontrun{
#' ## fit a model
#' fit <- brm(rating ~ treat + period + carry + (1|subject), 
#'            data = inhaler)
#' 
#' ## extract fitted values
#' fitted_values <- fitted(fit)
#' head(fitted_values)
#' 
#' ## plot fitted means against actual response
#' dat <- as.data.frame(cbind(Y = standata(fit)$Y, fitted_values))
#' ggplot(dat) + geom_point(aes(x = Estimate, y = Y))
#' }
#' 
#' @export 
fitted.brmsfit <- function(object, newdata = NULL, re_formula = NULL,
                           scale = c("response", "linear"),
                           allow_new_levels = FALSE, 
                           sample_new_levels = "uncertainty", 
                           new_objects = list(), incl_autocor = TRUE, 
                           dpar = NULL, resp = NULL, subset = NULL, 
                           nsamples = NULL, sort = FALSE, nug = NULL,
                           summary = TRUE, robust = FALSE, 
                           probs = c(0.025, 0.975), ...) {
  scale <- match.arg(scale)
  contains_samples(object)
  object <- restructure(object)
  draws_args <- nlist(
    x = object, newdata, re_formula, incl_autocor, allow_new_levels, 
    sample_new_levels, new_objects, resp, subset, nsamples, nug,
    check_response = FALSE
  )
  draws <- do.call(extract_draws, draws_args)
  fitted_args <- nlist(draws, scale, dpar, summary, robust, probs, sort)
  do.call(fitted_internal, fitted_args)
}

#' @rdname fitted.brmsfit
#' @aliases posterior_linpred
#' @method posterior_linpred brmsfit
#' @export
#' @export posterior_linpred
#' @importFrom rstantools posterior_linpred
posterior_linpred.brmsfit <- function(
  object, transform = FALSE, newdata = NULL, re_formula = NULL,
  re.form = NULL, allow_new_levels = FALSE, sample_new_levels = "uncertainty",
  new_objects = list(), incl_autocor = TRUE, dpar = NULL, resp = NULL,
  subset = NULL, nsamples = NULL, sort = FALSE, nug = NULL, 
  robust = FALSE, probs = c(0.025, 0.975), ...
) {
  cl <- match.call()
  cl[[1]] <- quote(fitted)
  cl$summary <- FALSE
  cl$scale <- if (transform) "response" else "linear"
  if ("re.form" %in% names(cl)) {
    cl$re_formula <- cl$re.form
    cl$re.form <- NULL
  }
  eval(cl, parent.frame())
}

#' Extract Model Residuals from brmsfit Objects
#' 
#' @inheritParams predict.brmsfit
#' @param type The type of the residuals, 
#'  either \code{"ordinary"} or \code{"pearson"}. 
#'  More information is provided under 'Details'. 
#' @param method Indicates the method to compute
#'  model implied values. Either \code{"fitted"}
#'  (predicted values of the regression curve) or
#'  \code{"predict"} (predicted response values). 
#'  Using \code{"predict"} is recommended
#'  but \code{"fitted"} is the current default for 
#'  reasons of backwards compatibility.
#' 
#' @details Residuals of type \code{ordinary} 
#'  are of the form \eqn{R = Y - Yp}, where \eqn{Y} is the observed 
#'  and \eqn{Yp} is the predicted response.
#'  Residuals of type \code{pearson} are 
#'  of the form \eqn{R = (Y - Yp) / SD(Y)},
#'  where \eqn{SD(Y)} is an estimation of the standard deviation 
#'  of \eqn{Y}. 
#'   
#'  Currently, \code{residuals.brmsfit} does not support 
#'  \code{categorical} or ordinal models. 
#'  
#'  Method \code{predictive_error.brmsfit} is an alias of 
#'  \code{residuals.brmsfit} with \code{method = "predict"} and
#'  \code{summary = FALSE}. 
#' 
#' @return Model residuals. If \code{summary = TRUE} 
#'  this is a N x C matrix and if \code{summary = FALSE} 
#'  a S x N matrix, where S is the number of samples, 
#'  N is the number of observations, and C is equal to 
#'  \code{length(probs) + 2}.  
#' 
#' @examples 
#' \dontrun{
#' ## fit a model
#' fit <- brm(rating ~ treat + period + carry + (1|subject), 
#'            data = inhaler, cores = 2)
#' 
#' ## extract residuals 
#' res <- residuals(fit, summary = TRUE)
#' head(res)
#' }
#' 
#' @export
residuals.brmsfit <- function(object, newdata = NULL, re_formula = NULL, 
                              type = c("ordinary", "pearson"),
                              method = c("fitted", "predict"),
                              allow_new_levels = FALSE, 
                              sample_new_levels = "uncertainty",
                              new_objects = list(), incl_autocor = TRUE, 
                              resp = NULL, subset = NULL, nsamples = NULL, 
                              sort = FALSE, nug = NULL, summary = TRUE, 
                              robust = FALSE, probs = c(0.025, 0.975), ...) {
  type <- match.arg(type)
  method <- match.arg(method)
  contains_samples(object)
  object <- restructure(object)
  family_names <- family_names(object)
  if (is_ordinal(family_names) || is_categorical(family_names)) {
    stop2("Residuals are not defined for ordinal or categorical models.")
  }
  newd_args <- nlist(
    object, newdata, re_formula, allow_new_levels,
    new_objects, check_response = TRUE, internal = TRUE
  )
  sdata <- do.call(standata, newd_args)
  if (any(grepl("^cens_", names(sdata)))) {
    warning2("'residuals' may not be meaningful for censored models.")
  }
  if (is.null(subset) && !is.null(nsamples)) {
    subset <- sample(nsamples(object), nsamples)
  }
  pred_args <- nlist(
    object, newdata, re_formula, allow_new_levels,
    sample_new_levels, new_objects, incl_autocor, 
    resp, subset, nug, summary = FALSE, sort = TRUE
  )
  mu <- do.call(method, pred_args)
  if (length(dim(mu)) == 3L) {
    # multivariate model
    Y <- sdata[paste0("Y_", dimnames(mu)[[3]])]
    names(Y) <- dimnames(mu)[[3]]
    Y <- lapply(Y, as_draws_matrix, dim = dim(mu)[1:2])
    Y <- do.call(abind, c(Y, along = 3))
  } else {
    Y <- sdata[[paste0("Y", usc(resp))]]
    Y <- as_draws_matrix(Y, dim = dim(mu))
  }
  res <- Y - mu
  remove(Y, mu)
  if (type == "pearson") {
    # get predicted standard deviation for each observation
    pred_args$summary <- TRUE
    pred <- do.call(predict, pred_args)
    if (length(dim(pred)) == 3L) {
      sd_pred <- array2list(pred[, 2, ])
      sd_pred <- lapply(sd_pred, as_draws_matrix, dim = dim(res)[1:2])
      sd_pred <- do.call(abind, c(sd_pred, along = 3))
    } else {
      sd_pred <- as_draws_matrix(pred[, 2], dim = dim(res))
    }
    res <- res / sd_pred
  }
  res <- reorder_obs(res, attr(sdata, "old_order"), sort = sort)
  if (summary) {
    res <- posterior_summary(res, probs = probs, robust = robust)
  }
  res
}

#' @rdname residuals.brmsfit
#' @aliases predictive_error
#' @method predictive_error brmsfit
#' @export
#' @export predictive_error
#' @importFrom rstantools predictive_error
predictive_error.brmsfit <- function(object, newdata = NULL, re_formula = NULL,
                                     re.form = NULL, allow_new_levels = FALSE, 
                                     sample_new_levels = "uncertainty",
                                     new_objects = list(), incl_autocor = TRUE, 
                                     resp = NULL, subset = NULL, nsamples = NULL, 
                                     sort = FALSE, nug = NULL, robust = FALSE, 
                                     probs = c(0.025, 0.975), ...) {
  cl <- match.call()
  cl[[1]] <- quote(residuals)
  cl$method <- "predict"
  cl$type <- "ordinary"
  cl$summary <- FALSE
  if ("re.form" %in% names(cl)) {
    cl$re_formula <- cl$re.form
    cl$re.form <- NULL
  }
  eval(cl, parent.frame())
}

#' @rdname pp_average
#' @export
pp_average.brmsfit <- function(
  x, ..., weights = "loo",
  method = c("predict", "fitted", "residuals"),
  newdata = NULL, re_formula = NULL,
  allow_new_levels = FALSE,
  sample_new_levels = "uncertainty", 
  new_objects = list(), incl_autocor = TRUE, 
  resp = NULL, nsamples = NULL, sort = FALSE, nug = NULL, 
  summary = TRUE, robust = FALSE, probs = c(0.025, 0.975),
  more_args = NULL, control = NULL
) {
  models <- list(x, ...)
  if (!match_response(models)) {
    stop2("Can only average models predicting the same response.")
  }
  method <- match.arg(method)
  if (is.null(nsamples)) {
    nsamples <- nsamples(x)
  }
  model_names <- ulapply(substitute(list(...))[-1], deparse_combine)
  model_names <- c(deparse_combine(substitute(x)), model_names)
  if (!is.numeric(weights)) {
    fun <- tolower(weights)
    fun <- match.arg(fun, c("loo", "waic", "kfold", "bridge"))
    args <- c(models, control)
    if (fun %in% c("loo", "waic", "kfold")) {
      args$compare <- FALSE
      ename <- switch(fun, loo = "looic", waic = "waic", kfold = "kfoldic")
      ics <- ulapply(SW(do.call(fun, args)), "[[", ename)
      ic_diffs <- ics - min(ics)
      weights <- exp(- ic_diffs / 2)
    } else if (weights %in% "bridge") {
      weights <- do.call("post_prob", args)
    }
  } else {
    if (length(weights) != length(models)) {
      stop2("If numeric, 'weights' must have the same length ",
            "as the number of models.")
    }
    if (any(weights < 0)) {
      stop2("If numeric, 'weights' must be positive.")
    }
  }
  weights <- weights / sum(weights)
  nsamples <- round(weights * nsamples)
  names(weights) <- names(nsamples) <- model_names
  method_args <- nlist(
    newdata, re_formula, allow_new_levels, sample_new_levels,
    new_objects, incl_autocor, resp, sort, nug, summary = FALSE
  )
  method_args <- c(method_args, more_args)
  out <- named_list(model_names)
  for (i in seq_along(out)) {
    if (nsamples[i] > 0) {
      method_args$object <- models[[i]]
      method_args$nsamples <- nsamples[i]
      out[[i]] <- do.call(method, method_args)
    }
  }
  out <- do.call(rbind, out)
  if (summary) {
    out <- posterior_summary(out, probs = probs, robust = robust) 
  }
  attr(out, "weights") <- weights
  attr(out, "nsamples") <- nsamples
  out
}

#' Compute a Bayesian version of R-squared for regression models
#' 
#' @aliases bayes_R2
#' 
#' @inheritParams predict.brmsfit
#' 
#' @return If \code{summary = TRUE} a 1 x C matrix is returned
#'  (\code{C = length(probs) + 2}) containing summary statistics
#'  of Bayesian R-squared values.
#'  If \code{summary = FALSE} the posterior samples of the R-squared values
#'  are returned in a S x 1 matrix (S is the number of samples).
#'  
#' @details For an introduction to the approach, see
#'   \url{https://github.com/jgabry/bayes_R2/blob/master/bayes_R2.pdf}.
#'  
#' @examples 
#' \dontrun{
#' fit <- brm(mpg ~ wt + cyl, data = mtcars)
#' summary(fit)
#' bayes_R2(fit)
#' 
#' # compute R2 with new data
#' nd <- data.frame(mpg = c(10, 20, 30), wt = c(4, 3, 2), cyl = c(8, 6, 4))
#' bayes_R2(fit, newdata = nd)
#' }
#' 
#' @method bayes_R2 brmsfit
#' @importFrom rstantools bayes_R2
#' @export bayes_R2
#' @export
bayes_R2.brmsfit <- function(object, newdata = NULL, re_formula = NULL, 
                             allow_new_levels = FALSE, 
                             sample_new_levels = "uncertainty",
                             new_objects = list(), incl_autocor = TRUE, 
                             subset = NULL, nsamples = NULL, resp = NULL,
                             nug = NULL, summary = TRUE, robust = FALSE, 
                             probs = c(0.025, 0.975), ...) {
  # do it like residuals.brmsfit
  contains_samples(object)
  object <- restructure(object)
  family_names <- family_names(object)
  if (is_ordinal(family_names) || is_categorical(family_names)) {
    stop2("Residuals are not defined for ordinal or categorical models.")
  }
  use_stored_ic <- !length(
    intersect(names(match.call()), args_not_for_reloo())
  )
  if (use_stored_ic && is.matrix(object[["R2"]])) {
    R2 <- object[["R2"]]
  } else {
    newd_args <- nlist(
      object, newdata, re_formula, allow_new_levels,
      new_objects, check_response = TRUE, internal = TRUE
    )
    sdata <- do.call(standata, newd_args)
    if (any(grepl("^cens_", names(sdata)))) {
      warning2("'bayes_R2' may not be meaningful for censored models.")
    }
    pred_args <- nlist(
      object, newdata, re_formula, allow_new_levels,
      sample_new_levels, new_objects, incl_autocor, resp,
      subset, nsamples, nug, summary = FALSE, sort = TRUE
    )
    ypred <- do.call(fitted, pred_args)
    # see https://github.com/jgabry/bayes_R2/blob/master/bayes_R2.pdf
    .bayes_R2 <- function(y, ypred) {
      y <- as.numeric(y)
      e <- - 1 * sweep(ypred, 2, y)
      var_ypred <- matrixStats::rowVars(ypred)
      var_e <- matrixStats::rowVars(e)
      return(as.matrix(var_ypred / (var_ypred + var_e)))
    }
    if (is.matrix(ypred)) {
      if (!length(resp)) resp <- ""
      resp <- usc(resp)
      y <- sdata[[paste0("Y", resp)]]
      R2 <- .bayes_R2(y, ypred)
    } else {
      resp <- usc(dimnames(ypred)[[3]])
      R2 <- named_list(resp)
      for (i in seq_along(R2)) {
        y <- as.numeric(sdata[[paste0("Y", resp[i])]])
        R2[[i]] <- .bayes_R2(y, ypred[, , i])
      }
      R2 <- do.call(cbind, R2)
    }
    colnames(R2) <- paste0("R2", resp)
  }
  if (summary) {
    R2 <- posterior_summary(R2, probs = probs, robust = robust)
  }
  R2
}

#' Update \pkg{brms} models
#' 
#' This method allows to update an existing \code{brmsfit} object
#' 
#' @param object An object of class \code{brmsfit}.
#' @param formula. Changes to the formula; for details see 
#'   \code{\link{update.formula}} and \code{\link{brmsformula}}.
#' @param newdata Optional \code{data.frame} 
#'   to update the model with new data.
#' @param recompile Logical, indicating whether the Stan model should 
#'   be recompiled. If \code{NULL} (the default), \code{update} tries
#'   to figure out internally, if recompilation is necessary. 
#'   Setting it to \code{FALSE} will cause all Stan code changing 
#'   arguments to be ignored. 
#' @param ... Other arguments passed to \code{\link{brm}}.
#'  
#' @details Sometimes, when updating the model formula, 
#'  it may happen that \R complains about a mismatch
#'  between \code{model frame} and \code{formula}.
#'  This error can be avoided by supplying your original data
#'  again via argument \code{newdata}.
#'  
#' @examples 
#' \dontrun{
#' fit1 <- brm(time | cens(censored) ~ age * sex + disease + (1|patient), 
#'             data = kidney, family = gaussian("log"))
#' summary(fit1)
#' 
#' ## remove effects of 'disease'
#' fit2 <- update(fit1, formula. = ~ . - disease)
#' summary(fit2)
#' 
#' ## remove the group specific term of 'patient' and
#' ## change the data (just take a subset in this example)
#' fit3 <- update(fit1, formula. = ~ . - (1|patient), 
#'                newdata = kidney[1:38, ])
#' summary(fit3)
#' 
#' ## use another family and add population-level priors
#' fit4 <- update(fit1, family = weibull(), inits = "0",
#'                prior = set_prior("normal(0,5)"))
#' summary(fit4)
#' }
#'
#' @export
update.brmsfit <- function(object, formula., newdata = NULL, 
                           recompile = NULL, ...) {
  dots <- list(...)
  testmode <- isTRUE(dots[["testmode"]])
  dots$testmode <- NULL
  if ("data" %in% names(dots)) {
    # otherwise the data name cannot be found by substitute 
    stop2("Please use argument 'newdata' to update the data.")
  }
  object <- restructure(object)
  if (isTRUE(object$version$brms < "2.0.0")) {
    warning2("Updating models fitted with older versions of brms may fail.")
  }
  if (missing(formula.)) {
    dots$formula <- object$formula
    if (!is.null(dots[["family"]])) {
      dots$formula <- dots$formula + check_family(dots$family)
    } 
    if (!is.null(dots[["autocor"]])) {
      dots$formula <- dots$formula + check_autocor(dots$autocor)
    }
  } else {
    # TODO: restructure updating of the model formula
    if (is.mvbrmsformula(formula.) || is.mvbrmsformula(object$formula)) {
      stop2("Updating formulas of multivariate models is not yet possible.")
    }
    if (is.brmsformula(formula.)) {
      nl <- get_nl(formula.)
    } else {
      formula. <- as.formula(formula.) 
      nl <- get_nl(formula(object))
    }
    family <- get_arg("family", formula., dots, object)
    autocor <- get_arg("autocor", formula., dots, object)
    dots$formula <- bf(formula., family = family, autocor = autocor, nl = nl)
    if (is_nonlinear(object)) {
      if (length(setdiff(all.vars(dots$formula$formula), ".")) == 0L) {
        dots$formula <- update(object$formula, dots$formula, mode = "keep")
      } else {
        dots$formula <- update(object$formula, dots$formula, mode = "replace")
        message("Argument 'formula.' will completely replace the ", 
                "original formula in non-linear models.")
      }
    } else {
      mvars <- all.vars(dots$formula$formula)
      mvars <- setdiff(mvars, c(names(object$data), "."))
      if (length(mvars) && is.null(newdata)) {
        stop2("New variables found: ", collapse_comma(mvars),
              "\nPlease supply your data again via argument 'newdata'.")
      }
      dots$formula <- update(formula(object), dots$formula)
    }
  }
  
  arg_names <- c("prior", "cov_ranef", "stan_vars", "stan_funs")
  old_args <- setdiff(arg_names, names(dots))
  dots[old_args] <- object[old_args]
  if (!is.null(newdata)) {
    dots$data <- newdata
  } else {
    dots$data <- rm_attr(object$data, c("terms", "brmsframe"))
  }
  if (!is.null(dots$prior)) {
    if (!is.brmsprior(dots$prior)) { 
      stop2("Invalid 'prior' argument.")
    }
    dots$prior <- rbind(dots$prior, object$prior)
    dupl_priors <- duplicated(dots$prior[, rcols_prior()])
    dots$prior <- dots$prior[!dupl_priors, ]
  }
  if (is.null(dots$sample_prior)) {
    dots$sample_prior <- attr(object$prior, "sample_prior")
    if (is.null(dots$sample_prior)) {
      has_prior_pars <- any(grepl("^prior_", parnames(object)))
      dots$sample_prior <- if (has_prior_pars) "yes" else "no"
    }
  }
  if (is.null(dots$save_ranef)) {
    dots$save_ranef <- isTRUE(attr(object$exclude, "save_ranef"))
  }
  if (is.null(dots$save_mevars)) {
    dots$save_mevars <- isTRUE(attr(object$exclude, "save_mevars"))
  }
  if (is.null(dots$save_all_pars)) {
    dots$save_all_pars <- isTRUE(attr(object$exclude, "save_all_pars"))
  }
  if (is.null(dots$sparse)) {
    dots$sparse <- grepl("sparse matrix", stancode(object))
  }
  dots$iter <- first_not_null(dots$iter, object$fit@sim$iter)
  # brm computes warmup automatically based on iter 
  dots$chains <- first_not_null(dots$chains, object$fit@sim$chains)
  dots$thin <- first_not_null(dots$thin, object$fit@sim$thin)
  control <- attr(object$fit@sim$samples[[1]], "args")$control
  control <- control[setdiff(names(control), names(dots$control))]
  dots$control[names(control)] <- control
  
  if (is.null(recompile)) {
    # only recompile if new and old stan code do not match
    new_stancode <- suppressMessages(do.call(make_stancode, dots))
    # stan code may differ just because of the version number (#288)
    new_stancode <- sub("^[^\n]+\n", "", new_stancode)
    old_stancode <- stancode(object, version = FALSE)
    recompile <- !is_equal(new_stancode, old_stancode)
    if (recompile) {
      message("The desired updates require recompiling the model") 
    }
  }
  recompile <- as_one_logical(recompile)
  if (recompile) {
    # recompliation is necessary
    dots$fit <- NA
    if (!is.null(newdata)) {
      dots$data.name <- Reduce(paste, deparse(substitute(newdata)))
      dots$data.name <- substr(dots$data.name, 1, 50)
    } else {
      dots$data.name <- object$data.name
    }
    if (!testmode) {
      object <- do.call(brm, dots)
    }
  } else {
    # refit the model without compiling it again
    if (!is.null(dots$formula)) {
      object$formula <- dots$formula
      dots$formula <- NULL
    }
    bterms <- parse_bf(object$formula)
    object$data <- update_data(dots$data, bterms = bterms)
    object$family <- get_element(object$formula, "family")
    object$autocor <- get_element(object$formula, "autocor")
    object$ranef <- tidy_ranef(bterms, data = object$data)
    object$stan_vars <- validate_stanvars(dots$stan_vars)
    if (!is.null(newdata)) {
      object$data.name <- Reduce(paste, deparse(substitute(newdata)))
    }
    if (!is.null(dots$sample_prior)) {
      dots$sample_prior <- check_sample_prior(dots$sample_prior)
      attr(object$prior, "sample_prior") <- dots$sample_prior
    }
    object$exclude <- exclude_pars(
      bterms, data = object$data, ranef = object$ranef, 
      save_ranef = dots$save_ranef, save_mevars = dots$save_mevars,
      save_all_pars = dots$save_all_pars
    )
    if (!is.null(dots$algorithm)) {
      aopts <- c("sampling", "meanfield", "fullrank")
      algorithm <- match.arg(dots$algorithm, aopts)
      dots$algorithm <- object$algorithm <- algorithm
    } else if (!is.null(object$algorithm)) {
      dots$algorithm <- object$algorithm
    }
    if (!testmode) {
      object <- do.call(brm, c(list(fit = object), dots))
    }
  }
  object
}

#' @export
#' @describeIn WAIC \code{WAIC} method for \code{brmsfit} objects
WAIC.brmsfit <- function(x, ..., compare = TRUE, newdata = NULL, 
                         re_formula = NULL, allow_new_levels = FALSE,
                         sample_new_levels = "uncertainty", resp = NULL,
                         new_objects = list(), subset = NULL,
                         nsamples = NULL, pointwise = NULL, nug = NULL) {
  models <- list(x, ...)
  model_names <- ulapply(substitute(list(...))[-1], deparse_combine)
  model_names <- c(deparse_combine(substitute(x)), model_names)
  if (is.null(subset) && !is.null(nsamples)) {
    subset <- sample(nsamples(x), nsamples)
  }
  pointwise <- set_pointwise(x, pointwise, newdata, subset)
  use_stored_ic <- !any(names(match.call()) %in% args_not_for_reloo())
  args <- nlist(
    models, model_names, ic = "waic", use_stored_ic,
    newdata, re_formula, resp, subset, nug, allow_new_levels, 
    sample_new_levels, new_objects, pointwise, compare
  )
  do.call(compute_ics, args)
}

#' @importFrom loo waic
#' @export waic
#' @export
waic.brmsfit <- function(x, ..., compare = TRUE, newdata = NULL,
                         re_formula = NULL, allow_new_levels = FALSE,
                         sample_new_levels = "uncertainty",
                         resp = NULL, new_objects = list(), subset = NULL, 
                         nsamples = NULL, pointwise = NULL, nug = NULL) {
  cl <- match.call()
  cl[[1]] <- quote(WAIC)
  eval(cl, parent.frame())
}

#' @export
#' @describeIn LOO \code{LOO} method for \code{brmsfit} objects
LOO.brmsfit <- function(x, ..., compare = TRUE, reloo = FALSE, 
                        newdata = NULL, re_formula = NULL, 
                        allow_new_levels = FALSE, 
                        sample_new_levels = "uncertainty", 
                        resp = NULL, new_objects = list(), 
                        subset = NULL, nsamples = NULL, pointwise = NULL,
                        nug = NULL, k_threshold = 0.7, 
                        update_args = list(), cores = 1, 
                        wcp = 0.2, wtrunc = 3/4) {
  models <- list(x, ...)
  model_names <- ulapply(substitute(list(...))[-1], deparse_combine)
  model_names <- c(deparse_combine(substitute(x)), model_names)
  if (is.null(subset) && !is.null(nsamples)) {
    subset <- sample(nsamples(x), nsamples)
  }
  pointwise <- set_pointwise(x, pointwise, newdata, subset)
  loo_args <- nlist(wcp, wtrunc, cores)
  not_for_reloo <- intersect(names(match.call()), args_not_for_reloo())
  if (reloo && length(not_for_reloo)) {
    not_for_reloo <- collapse_comma(not_for_reloo)
    stop2("Cannot use 'reloo' with arguments ", not_for_reloo, ".")
  }
  use_stored_ic <- !length(not_for_reloo)
  args <- nlist(
    models, model_names, ic = "loo", use_stored_ic, loo_args,
    newdata, re_formula, resp, subset, nug, allow_new_levels, 
    sample_new_levels, new_objects, pointwise, compare,
    reloo, k_threshold, update_args
  )
  do.call(compute_ics, args)
}

#' @importFrom loo loo
#' @export loo
#' @export
loo.brmsfit <-  function(x, ..., compare = TRUE, reloo = FALSE, 
                         newdata = NULL, re_formula = NULL, 
                         allow_new_levels = FALSE, 
                         sample_new_levels = "uncertainty", 
                         new_objects = list(), resp = NULL, 
                         subset = NULL, nsamples = NULL, pointwise = NULL, 
                         nug = NULL, k_threshold = 0.7, 
                         update_args = list(), cores = 1, 
                         wcp = 0.2, wtrunc = 3/4) {
  cl <- match.call()
  cl[[1]] <- quote(LOO)
  eval(cl, parent.frame())
}

#' @export
#' @describeIn kfold \code{kfold} method for \code{brmsfit} objects
kfold.brmsfit <- function(x, ..., compare = TRUE,
                          K = 10, Ksub = NULL, exact_loo = FALSE, 
                          group = NULL, newdata = NULL, resp = NULL,
                          save_fits = FALSE, update_args = list()) {
  models <- list(x, ...)
  model_names <- ulapply(substitute(list(...))[-1], deparse_combine)
  model_names <- c(deparse_combine(substitute(x)), model_names)
  use_stored_ic <- ulapply(models, 
    function(x) is.brmsfit(x) && is_equal(x$kfold$K, K)
  )
  args <- nlist(
    models, model_names, ic = "kfold", K, save_fits, 
    use_stored_ic, compare, update_args, newdata, 
    resp, Ksub, exact_loo, group
  )
  do.call(compute_ics, args)
}

#' Compute Weighted Expectations Using LOO
#' 
#' These functions are wrappers around the \code{\link[loo]{E_loo}} function 
#' of the \pkg{loo} package.
#'
#' @aliases loo_predict loo_linpred loo_predictive_interval
#' 
#' @param object An object of class \code{brmsfit}.
#' @param type The statistic to be computed on the results. 
#'   Can by either \code{"mean"} (default), \code{"var"}, or
#'   \code{"quantile"}.
#' @param probs A vector of quantiles to compute. 
#'   Only used if \code{type = quantile}.
#' @param scale Passed to \code{\link[brms:fitted.brmsfit]{fitted}}.
#' @param prob For \code{loo_predictive_interval}, a scalar in \eqn{(0,1)}
#'   indicating the desired probability mass to include in the intervals. The
#'   default is \code{prob = 0.9} (\eqn{90}\% intervals).
#' @param lw An optional matrix of (smoothed) log-weights. If \code{lw} is 
#'   missing then \code{\link[loo]{psislw}} is executed internally, which may be
#'   time consuming for models fit to very large datasets. 
#'   If \code{lw} is specified, arguments passed via \code{...} may be ignored.
#' @param ... Optional arguments passed to the underlying methods that is 
#'   \code{\link[brms:log_lik.brmsfit]{log_lik}}, as well as
#'   \code{\link[brms:predict.brmsfit]{predict}} or
#'   \code{\link[brms:fitted.brmsfit]{fitted}}. 
#' @inheritParams LOO.brmsfit
#'   
#' @return \code{loo_predict} and \code{loo_linpred} return a vector with one 
#'   element per observation. The only exception is if \code{type = "quantile"} 
#'   and \code{length(probs) >= 2}, in which case a separate vector for each 
#'   element of \code{probs} is computed and they are returned in a matrix with 
#'   \code{length(probs)} rows and one column per observation.
#'   
#'   \code{loo_predictive_interval} returns a matrix with one row per 
#'   observation and two columns. 
#'   \code{loo_predictive_interval(..., prob = p)} is equivalent to 
#'   \code{loo_predict(..., type = "quantile", probs = c(a, 1-a))} with 
#'   \code{a = (1 - p)/2}, except it transposes the result and adds informative 
#'   column names.
#'   
#' @examples
#' \dontrun{
#' ## data from help("lm")
#' ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
#' trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
#' d <- data.frame(
#'   weight = c(ctl, trt), 
#'   group = gl(2, 10, 20, labels = c("Ctl", "Trt"))
#' ) 
#' fit <- brm(weight ~ group, data = d)
#' loo_predictive_interval(fit, prob = 0.8)
#' 
#' ## optionally log-weights can be pre-computed and reused
#' psis <- loo::psislw(-log_lik(fit), cores = 2)
#' loo_predictive_interval(fit, prob = 0.8, lw = psis$lw_smooth)
#' loo_predict(fit, type = "var", lw = psis$lw_smooth)
#' }
#' 
#' @method loo_predict brmsfit
#' @importFrom rstantools loo_predict
#' @export loo_predict
#' @export 
loo_predict.brmsfit <- function(object, type = c("mean", "var", "quantile"), 
                                resp = NULL, probs = 0.5, lw = NULL, cores = 1, 
                                wcp = 0.2, wtrunc = 3/4, ...) {
  type <- match.arg(type)
  stopifnot_resp(object, resp)
  loo_args <- nlist(cores, wcp, wtrunc)
  lw <- loo_weights(
    object, lw = lw, log = TRUE, 
    loo_args = loo_args, resp = resp, ...
  )
  preds <- predict(object, summary = FALSE, resp = resp, ...)
  loo::E_loo(x = preds, lw = lw, type = type, probs = probs)
}

#' @rdname loo_predict.brmsfit
#' @method loo_linpred brmsfit
#' @importFrom rstantools loo_linpred
#' @export loo_linpred
#' @export 
loo_linpred.brmsfit <- function(object, type = c("mean", "var", "quantile"), 
                                resp = NULL, probs = 0.5, scale = "linear", lw = NULL, 
                                cores = 1, wcp = 0.2, wtrunc = 3/4, ...) {
  type <- match.arg(type)
  stopifnot_resp(object, resp)
  if (is_ordinal(object$family) || is_categorical(object$family)) {
    stop2("Method 'loo_linpred' is not yet working ", 
          "for categorical or ordinal models")
  }
  loo_args <- nlist(cores, wcp, wtrunc)
  lw <- loo_weights(
    object, lw = lw, log = TRUE, 
    loo_args = loo_args, resp = resp, ...
  )
  preds <- fitted(object, scale = scale, summary = FALSE, resp = resp, ...)
  loo::E_loo(x = preds, lw = lw, type = type, probs = probs)
}

#' @rdname loo_predict.brmsfit
#' @method loo_predictive_interval brmsfit
#' @importFrom rstantools loo_predictive_interval
#' @export loo_predictive_interval
#' @export
loo_predictive_interval.brmsfit <- function(object, prob = 0.9,
                                            lw = NULL, ...) {
  if (length(prob) != 1L) {
    stop2("Argument 'prob' should be of length 1.")
  }
  alpha <- (1 - prob) / 2
  probs <- c(alpha, 1 - alpha)
  labs <- paste0(100 * probs, "%")
  intervals <- loo_predict(
    object, type = "quantile", probs = probs, lw = lw, ...
  )
  rownames(intervals) <- labs
  t(intervals)
}
 
#' Compute the Pointwise Log-Likelihood
#' 
#' @aliases log_lik logLik.brmsfit
#' 
#' @param object A fitted model object of class \code{brmsfit}. 
#' @inheritParams predict.brmsfit
#' @param combine Only relevant in multivariate models.
#'   Indicates if the log-likelihoods of the submodels should
#'   be combined per observation (i.e. added together; the default) 
#'   or if the log-likelihoods should be returned separately.
#' @param pointwise A flag indicating whether to compute the full
#'   log-likelihood matrix at once (the default), or just return
#'   the likelihood function along with all data and samples
#'   required to compute the log-likelihood separately for each
#'   observation. The latter option is rarely useful when
#'   calling \code{log_lik} directly, but rather when computing
#'   \code{\link[brms:WAIC]{WAIC}} or \code{\link[brms:LOO]{LOO}}.
#' @param ... Currently ignored
#' 
#' @return Usually, an S x N matrix containing the pointwise log-likelihood
#'  samples, where S is the number of samples and N is the number 
#'  of observations in the data. For multivariate models and if 
#'  \code{combine} is \code{FALSE}, an S x N x R array is returned, 
#'  where R is the number of response variables.
#'  If \code{pointwise = TRUE}, the output is a function
#'  with a \code{draws} attribute containing all relevant
#'  data and posterior samples.
#' 
#' @aliases log_lik
#' @method log_lik brmsfit
#' @export
#' @export log_lik
#' @importFrom rstantools log_lik
log_lik.brmsfit <- function(object, newdata = NULL, re_formula = NULL,
                            allow_new_levels = FALSE, 
                            sample_new_levels = "uncertainty", 
                            new_objects = list(), incl_autocor = TRUE, 
                            resp = NULL, subset = NULL, nsamples = NULL, 
                            pointwise = FALSE, nug = NULL, combine = TRUE, 
                            ...) {
  contains_samples(object)
  object <- restructure(object)
  draws_args <- nlist(
    x = object, newdata, re_formula, allow_new_levels, incl_autocor,
    sample_new_levels, new_objects, resp, subset, nsamples, nug, 
    check_response = TRUE
  )
  draws <- do.call(extract_draws, draws_args)
  if (pointwise) {
    stopifnot(combine)
    log_lik <- log_lik_pointwise
    attr(log_lik, "args") <- list(
      draws = draws, N = choose_N(draws), 
      S = draws$nsamples, data = data.frame()
    )
  } else {
    log_lik <- log_lik_internal(draws, combine = combine)
    if (anyNA(log_lik)) {
      warning2(
        "NAs were found in the log-likelihood. Possibly this is because ",
        "some of your predictors contain NAs. If you use 'mi' terms, try ", 
        "setting 'resp' to those response variables without missing values."
      )
    }
  }
  log_lik
}

#' @export
logLik.brmsfit <- function(object, newdata = NULL, re_formula = NULL,
                           allow_new_levels = FALSE,
                           sample_new_levels = "uncertainty",
                           new_objects = list(), incl_autocor = TRUE, 
                           resp = NULL, subset = NULL, nsamples = NULL, 
                           pointwise = FALSE, nug = NULL, ...) {
  cl <- match.call()
  cl[[1]] <- quote(log_lik)
  eval(cl, parent.frame())
}

#' @rdname pp_mixture
#' @export
pp_mixture.brmsfit <- function(x, newdata = NULL, re_formula = NULL,
                               allow_new_levels = FALSE, 
                               sample_new_levels = "uncertainty", 
                               new_objects = list(), incl_autocor = TRUE, 
                               resp = NULL, subset = NULL, nsamples = NULL, 
                               nug = NULL, summary = TRUE, 
                               robust = FALSE, probs = c(0.025, 0.975), 
                               log = FALSE, ...) {
  stopifnot_resp(x, resp)
  contains_samples(x)
  x <- restructure(x)
  if (is_mv(x)) {
    resp <- validate_resp(resp, x$formula$responses, multiple = FALSE)
    family <- x$family[[resp]]
  } else {
    family <- x$family
  }
  if (!is.mixfamily(family)) {
    stop2("Method 'pp_mixture' can only be applied on mixture models.")
  }
  draws_args <- nlist(
    x, newdata, re_formula, allow_new_levels, incl_autocor, 
    sample_new_levels, new_objects, subset, nsamples, nug,
    resp, check_response = TRUE
  )
  draws <- do.call(extract_draws, draws_args)
  stopifnot(is.brmsdraws(draws))
  draws$pp_mixture <- TRUE
  for (dp in names(draws$dpars)) {
    draws$dpars[[dp]] <- get_dpar(draws, dpar = dp)
  }
  N <- choose_N(draws)
  log_lik <- lapply(seq_len(N), log_lik_mixture, draws = draws)
  log_lik <- do.call(abind, c(log_lik, along = 3))
  log_lik <- aperm(log_lik, c(1, 3, 2))
  old_order <- draws$old_order
  sort <- isTRUE(ncol(log_lik) != length(old_order))
  log_lik <- reorder_obs(log_lik, old_order, sort = sort)
  if (!log) {
    log_lik <- exp(log_lik)
  }
  if (summary) {
    log_lik <- posterior_summary(log_lik, probs = probs, robust = robust)
    dimnames(log_lik) <- list(
      seq_len(nrow(log_lik)), colnames(log_lik),
      paste0("P(K = ", seq_len(dim(log_lik)[3]), " | Y)")
    )
  }
  log_lik
}

#' @rdname hypothesis
#' @export
hypothesis.brmsfit <- function(x, hypothesis, class = "b", group = "",
                               scope = c("standard", "ranef", "coef"),
                               alpha = 0.05, seed = NULL, ...) {
  # use a seed as prior_samples.brmsfit randomly permutes samples
  if (!is.null(seed)) {
    set.seed(seed) 
  }
  contains_samples(x)
  x <- restructure(x)
  group <- as_one_character(group)
  scope <- match.arg(scope)
  if (scope == "standard") {
    if (!length(class)) {
      class <- "" 
    }
    class <- as_one_character(class)
    if (class %in% c("sd", "cor") && nzchar(group)) {
      class <- paste0(class, "_", group, "__")
    } else if (nzchar(class)) {
      class <- paste0(class, "_")
    }
    out <- hypothesis_internal(
      x, hypothesis, class = class, alpha = alpha, ...
    )
  } else {
    co <- do.call(scope, list(x, summary = FALSE))
    if (!group %in% names(co)) {
      stop2("'group' should be one of ", collapse_comma(names(co)))
    }
    out <- hypothesis_coef(co[[group]], hypothesis, alpha = alpha, ...)
  }
  out
}

#' @rdname expose_functions
#' @export
expose_functions.brmsfit <- function(x, vectorize = FALSE, 
                                     env = globalenv(), ...) {
  vectorize <- as_one_logical(vectorize)
  if (vectorize) {
    funs <- expose_stan_functions(x$fit, env = environment(), ...)
    for (i in seq_along(funs)) {
      FUN <- Vectorize(get(funs[i], mode = "function"))
      assign(funs[i], FUN, pos = env) 
    }
  } else {
    funs <- expose_stan_functions(x$fit, env = env, ...)
  }
  invisible(funs)
}

#' @rdname diagnostic-quantities
#' @importFrom bayesplot log_posterior
#' @export log_posterior
#' @export
log_posterior.brmsfit <- function(object, ...) {
  contains_samples(object)
  bayesplot::log_posterior(object$fit, ...)
}

#' @rdname diagnostic-quantities
#' @importFrom bayesplot nuts_params
#' @export nuts_params
#' @export
nuts_params.brmsfit <- function(object, pars = NULL, ...) {
  contains_samples(object)
  bayesplot::nuts_params(object$fit, pars = pars, ...)
}

#' @rdname diagnostic-quantities
#' @importFrom bayesplot rhat
#' @export rhat
#' @export
rhat.brmsfit <- function(object, pars = NULL, ...) {
  contains_samples(object)
  bayesplot::rhat(object$fit, pars = pars, ...)
}

#' @rdname diagnostic-quantities
#' @importFrom bayesplot neff_ratio
#' @export neff_ratio
#' @export
neff_ratio.brmsfit <- function(object, pars = NULL, ...) {
  contains_samples(object)
  bayesplot::neff_ratio(object$fit, pars = pars, ...)
}

#' @rdname control_params
#' @export
control_params.brmsfit <- function(x, pars = NULL, ...) {
  contains_samples(x)
  out <- attr(x$fit@sim$samples[[1]], "args")$control
  if (!is.null(pars)) {
    out <- out[pars]
  }
  out
}

#' Log Marginal Likelihood via Bridge Sampling
#' 
#' Computes log marginal likelihood via bridge sampling,
#' which can be used in the computation of bayes factors
#' and posterior model probabilities.
#' The \code{brmsfit} method is just a thin wrapper around
#' the corresponding method for \code{stanfit} objects.
#' 
#' @aliases bridge_sampler
#' 
#' @param samples A \code{brmsfit} object.
#' @param ... Additional arguments passed to 
#'   \code{\link[bridgesampling:bridge_sampler]{bridge_sampler.stanfit}}.
#' 
#' @details Computing the marginal likelihood requires samples 
#'   of all variables defined in Stan's \code{parameters} block
#'   to be saved. Otherwise \code{bridge_sampler} cannot be computed.
#'   Thus, please set \code{save_all_pars = TRUE} in the call to \code{brm},
#'   if you are planning to apply \code{bridge_sampler} to your models.
#' 
#'   More details are provided under
#'   \code{\link[bridgesampling:bridge_sampler]{bridge_sampler}}.
#'   
#' @seealso \code{
#'   \link[brms:bayes_factor]{bayes_factor},
#'   \link[brms:post_prob]{post_prob}
#' }
#' 
#' @examples 
#' \dontrun{
#' # model with the treatment effect
#' fit1 <- brm(
#'   count ~ log_Age_c + log_Base4_c + Trt_c,
#'   data = epilepsy, family = negbinomial(), 
#'   prior = prior(normal(0, 1), class = b),
#'   save_all_pars = TRUE
#' )
#' summary(fit1)
#' bridge_sampler(fit1)
#' 
#' # model without the treatment effect
#' fit2 <- brm(
#'   count ~ log_Age_c + log_Base4_c,
#'   data = epilepsy, family = negbinomial(), 
#'   prior = prior(normal(0, 1), class = b),
#'   save_all_pars = TRUE
#' )
#' summary(fit2)
#' bridge_sampler(fit2)
#' }
#' 
#' @method bridge_sampler brmsfit
#' @importFrom bridgesampling bridge_sampler 
#' @export bridge_sampler 
#' @export
bridge_sampler.brmsfit <- function(samples, ...) {
  if (inherits(samples[["bridge"]], "bridge")) {
    if (!is.na(samples[["bridge"]]$logml)) {
      return(samples[["bridge"]]) 
    }
  }
  samples <- restructure(samples)
  if (samples$version$brms <= "1.8.0") {
    stop2(
      "Models fitted with brms 1.8.0 or lower are not ",
      "usable in method 'bridge_sampler'."
    )
  }
  if (is.cor_car(samples$autocor)) {
    warning2(
      "Some constants were dropped from the log-posterior ",
      "of CAR models so that the output of 'bridge_sampler' ",
      "may not be valid."
    )
  }
  sample_prior <- attr(samples$prior, "sample_prior")
  if (isTRUE(sample_prior %in% c("yes", "only"))) {
    stop2(
      "Models including prior samples are not usable ",
      "in method 'bridge_sampler'."
    )
  }
  # otherwise bridge_sampler might not work in a new R session
  stanfit_tmp <- suppressMessages(brm(fit = samples, chains = 0))$fit
  out <- try(
    bridge_sampler(samples$fit, stanfit_model = stanfit_tmp, ...),
    silent = TRUE
  )
  if (is(out, "try-error")) {
    stop2(
      "Bridgesampling failed. Did you set 'save_all_pars' ",
      "to TRUE when fitting your model?"
    )
  }
  out
}

#' Bayes Factors from Marginal Likelihoods
#' 
#' Compute Bayes factors from marginal likelihoods.
#' 
#' @aliases bayes_factor
#' 
#' @param x1 A \code{brmsfit} object
#' @param x2 Another \code{brmsfit} object based on the same responses.
#' @param log Report Bayes factors on the log-scale?
#' @param ... Additional arguments passed to 
#'   \code{\link[brms:bridge_sampler]{bridge_sampler}}.
#' 
#' @details Computing the marginal likelihood requires samples 
#'   of all variables defined in Stan's \code{parameters} block
#'   to be saved. Otherwise \code{bayes_factor} cannot be computed.
#'   Thus, please set \code{save_all_pars = TRUE} in the call to \code{brm},
#'   if you are planning to apply \code{bayes_factor} to your models.
#' 
#'   More details are provided under 
#'   \code{\link[bridgesampling:bayes_factor]{bayes_factor}}.
#'  
#' @seealso \code{
#'   \link[brms:bridge_sampler]{bridge_sampler},
#'   \link[brms:post_prob]{post_prob}
#' }
#' 
#' @examples 
#' \dontrun{
#' # model with the treatment effect
#' fit1 <- brm(
#'   count ~ log_Age_c + log_Base4_c + Trt_c,
#'   data = epilepsy, family = negbinomial(), 
#'   prior = prior(normal(0, 1), class = b),
#'   save_all_pars = TRUE
#' )
#' summary(fit1)
#' 
#' # model without the treatment effect
#' fit2 <- brm(
#'   count ~ log_Age_c + log_Base4_c,
#'   data = epilepsy, family = negbinomial(), 
#'   prior = prior(normal(0, 1), class = b),
#'   save_all_pars = TRUE
#' )
#' summary(fit2)
#' 
#' # compute the bayes factor
#' bayes_factor(fit1, fit2)
#' }
#' 
#' @method bayes_factor brmsfit
#' @importFrom bridgesampling bayes_factor
#' @export bayes_factor
#' @export
bayes_factor.brmsfit <- function(x1, x2, log = FALSE, ...) {
  match_response(list(x1, x2))
  bridge1 <- bridge_sampler(x1, ...)
  bridge2 <- bridge_sampler(x2, ...)
  bridgesampling::bf(bridge1, bridge2, log = log)
}

#' Posterior Model Probabilities from Marginal Likelihoods
#' 
#' Compute posterior model probabilities from marginal likelihoods.
#' The \code{brmsfit} method is just a thin wrapper around
#' the corresponding method for \code{bridge} objects.
#' 
#' @aliases post_prob
#' 
#' @param x A \code{brmsfit} object.
#' @param ... More \code{brmsfit} objects.
#' @param prior_prob Numeric vector with prior model probabilities. 
#'   If omitted, a uniform prior is used (i.e., all models are equally 
#'   likely a priori). The default \code{NULL} corresponds to equal 
#'   prior model weights.
#' @param model_names If \code{NULL} (the default) will use model names 
#'   derived from deparsing the call. Otherwise will use the passed 
#'   values as model names.
#' @param bs_args A list of additional arguments passed to 
#'   \code{\link[brms:bridge_sampler]{bridge_sampler}}.
#'   
#' @details Computing the marginal likelihood requires samples 
#'   of all variables defined in Stan's \code{parameters} block
#'   to be saved. Otherwise \code{post_prob} cannot be computed.
#'   Thus, please set \code{save_all_pars = TRUE} in the call to \code{brm},
#'   if you are planning to apply \code{post_prob} to your models.
#' 
#'   More details are provided under 
#'   \code{\link[bridgesampling:post_prob]{post_prob}}. 
#'   
#' @seealso \code{
#'   \link[brms:bridge_sampler]{bridge_sampler},
#'   \link[brms:bayes_factor]{bayes_factor}
#' }
#' 
#' @examples 
#' \dontrun{
#' # model with the treatment effect
#' fit1 <- brm(
#'   count ~ log_Age_c + log_Base4_c + Trt_c,
#'   data = epilepsy, family = negbinomial(), 
#'   prior = prior(normal(0, 1), class = b),
#'   save_all_pars = TRUE
#' )
#' summary(fit1)
#' 
#' # model without the treatent effect
#' fit2 <- brm(
#'   count ~ log_Age_c + log_Base4_c,
#'   data = epilepsy, family = negbinomial(), 
#'   prior = prior(normal(0, 1), class = b),
#'   save_all_pars = TRUE
#' )
#' summary(fit2)
#' 
#' # compute the posterior model probabilities
#' post_prob(fit1, fit2)
#' 
#' # specify prior model probabilities
#' post_prob(fit1, fit2, prior_prob = c(0.8, 0.2))
#' }
#' 
#' @method post_prob brmsfit
#' @importFrom bridgesampling post_prob
#' @export post_prob 
#' @export
post_prob.brmsfit <- function(x, ..., prior_prob = NULL, 
                              model_names = NULL,
                              bs_args = list()) {
  models <- list(x, ...)
  if (is.null(model_names)) {
    model_names <- ulapply(substitute(list(...))[-1], deparse_combine)
    model_names <- c(deparse_combine(substitute(x)), model_names)
  } else if (length(model_names) != length(models)) {
    stop2("Number of model names is not equal to the number of models.") 
  }
  for (i in seq_along(models)) {
    if (!is.brmsfit(models[[i]])) {
      stop2("Object '", model_names[i], "' is not of class 'brmsfit'.")
    }
  }
  match_response(models)
  bs <- vector("list", length(models))
  for (i in seq_along(models)) {
    bs[[i]] <- do.call(bridge_sampler, c(list(models[[i]]), bs_args))
  }
  do.call(post_prob, c(bs, nlist(prior_prob, model_names)))
}
