#' @export
parnames.brmsfit <- function(x, ...) {
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  dimnames(x$fit)$parameters
}

#' Extract Population-Level Estimates
#' 
#' Extract the population-level ('fixed') effects 
#' from a \code{brmsfit} object. 
#' 
#' @aliases fixef
#' 
#' @param object An object of class \code{brmsfit}
#' @param estimate A character vector specifying which coefficients 
#'  (e.g., "mean", "median", "sd", or "quantile") 
#'  should be calculated for the population-level effects.
#' @param ... Further arguments to be passed to the functions 
#'  specified in \code{estimate}
#' 
#' @return A matrix with one row per population-level effect 
#'   and one column per calculated estimate.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' fit <- brm(time | cens(censored) ~ age + sex + disease, 
#'            data = kidney, family = "exponential")
#' fixef(fit, estimate = c("mean", "sd"))
#' }
#' 
#' @export
#' @export fixef
#' @importFrom lme4 fixef
fixef.brmsfit <-  function(object, estimate = "mean", ...) {
  if (!is(object$fit, "stanfit") || !length(object$fit@sim)) 
    stop("The model does not contain posterior samples")
  pars <- parnames(object)
  fpars <- pars[grepl("^b_", pars)]
  if (!length(fpars)) 
    stop("The model does not contain fixed effects", call. = FALSE) 
  out <- posterior_samples(object, pars = fpars, exact_match = TRUE)
  out <- do.call(cbind, lapply(estimate, get_estimate, samples = out, ...))
  rownames(out) <- gsub("^b_", "", fpars)
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
#' @export
vcov.brmsfit <- function(object, correlation = FALSE, ...) {
  if (!is(object$fit, "stanfit") || !length(object$fit@sim)) 
    stop("The model does not contain posterior samples")
  pars <- parnames(object)
  fpars <- pars[grepl("^b_", pars)]
  if (!length(fpars)) 
    stop("The model does not contain fixed effects", call. = FALSE) 
  samples <- posterior_samples(object, pars = fpars, exact_match = TRUE)
  names(samples) <- sub("^b_", "", names(samples))
  if (correlation) {
    cor(samples) 
  } else {
    cov(samples)
  }
}

#' Extract Group-Level Estimates
#' 
#' Extract the group-level ('random') effects of each level 
#' from a \code{brmsfit} object. 
#' 
#' @aliases ranef
#' 
#' @param object An object of class \code{brmsfit}.
#' @param estimate The point estimate to be calculated 
#'  for the random effects, either "mean" or "median".
#' @param var logical; indicating if the covariance matrix 
#'  for each random effects should be computed.
#' @param ... Further arguments to be passed to the function 
#'  specified in \code{estimate}
#'
#' @return A list of matrices (one per grouping factor), 
#'  with factor levels as row names and 
#'  group-level effects as column names.
#'     
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}   
#'   
#' @examples
#' \dontrun{
#' fit <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1+Trt_c|visit), 
#'              data = epilepsy, family = "poisson", chains = 1)
#' ## group-level means with corresponding covariances
#' rf <- ranef(fit, var = TRUE)
#' attr(rf, "var")
#' ## group-level medians
#' ranef(fit, estimate = "median")                                                        
#' }
#' 
#' @export
#' @export ranef
#' @importFrom lme4 ranef
ranef.brmsfit <- function(object, estimate = "mean", var = FALSE, ...) {
  if (!is(object$fit, "stanfit") || !length(object$fit@sim)) 
    stop("The model does not contain posterior samples")
  if (!estimate %in% c("mean","median"))
    stop("Argument estimate must be either 'mean' or 'median'", call. = FALSE)
  if (!length(object$ranef))
    stop("The model does not contain random effects", call. = FALSE)
  group <- names(object$ranef)
  pars <- parnames(object)
  
  get_ranef <- function(i) {
    # get random effects of a grouping factor
    #
    # Args:
    #   i: index of a grouping factor
    g <- group[i]
    rnames <- object$ranef[[i]]
    nlpar <- get_nlpar(object$ranef[[i]], suffix = "_")
    rpars <- pars[grepl(paste0("^r_", nlpar, group[i],"\\["), pars)]
    if (!length(rpars))
      stop(paste0("The model does not contain random effects for group '",g,"'\n",
                  "You should use argument ranef = TRUE in function brm."),
           call. = FALSE)
    rdims <- object$fit@sim$dims_oi[[paste0("r_", nlpar, group[i])]]
    levels <- attr(object$ranef[[i]], "levels")
    if (is.null(levels)) {
      # avoid error in dimnames if levels are NULL 
      # for backwards compatibility with brms < 0.5.0 
      levels <- 1:rdims[1]
    }
    rs <- posterior_samples(object, pars = rpars, exact_match = TRUE)
    ncol <- ifelse(is.na(rdims[2]), 1, rdims[2])
    rs_array <- array(dim = c(rdims[1], ncol, nrow(rs)))
    k <- 0
    for (j in 1:ncol) {
      for (l in 1:rdims[1]) {
        k <- k + 1
        rs_array[l, j, ] <- rs[, k]
      }
    }
    out <- get_estimate(estimate, samples = rs_array, margin = 1:2, ...)
    colnames(out) <- rnames
    if(var) {
      Var <- array(dim = c(rep(ncol, 2), rdims[1]), 
                   dimnames = list(rnames, rnames, 1:rdims[1]))
      for (j in 1:rdims[1]) {
        if (is.na(rdims[2])) Var[, , j] <- var(rs_array[j, 1, ]) 
        else Var[, , j] <- cov(t(rs_array[j, , ])) 
      }
      dimnames(Var)[[3]] <- levels
      attr(out, "var") <- Var
    }
    rownames(out) <- levels
    if (nchar(nlpar)) 
      attr(out, "nlpar") <- get_nlpar(object$ranef[[i]])
    return(out)
  }
  
  ranef <- lapply(seq_along(group), get_ranef)
  names(ranef) <- group
  ranef 
} 

#' Extract model coefficients
#'
#' Extract model coefficients, which are the sum of population-level 
#' effects and corresponding group-level effects
#' 
#' @param object An object of class \code{brmsfit}
#' @inheritParams ranef.brmsfit
#'
#' @return A list of matrices (one per grouping factor), 
#'  with factor levels as row names and 
#'  coefficients as column names 
#'  
#' @examples
#' \dontrun{
#' fit <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1+Trt_c|visit), 
#'            data = epilepsy, family = "poisson", chains = 1)
#' ## extract fixed and random effects coefficients seperately
#' fixef(fit)
#' ranef(fit)
#' ## extract combined coefficients     
#' coef(fit)
#' }
#' 
#' @export
coef.brmsfit <- function(object, estimate = "mean", ...) {
  if (!estimate %in% c("mean","median"))
    stop("Argument estimate must be either 'mean' or 'median'", call. = FALSE)
  fixef <- fixef(object, estimate = estimate, ...)
  if (!length(object$ranef)) {
    return(fixef)  # no random effects present
  }
  coef <- ranef(object, estimate = estimate, ...)
  ranef_names <- unique(ulapply(coef, colnames))
  missing_fixef <- setdiff(ranef_names, rownames(fixef))
  if (length(missing_fixef)) {
    zero_mat <- matrix(0, nrow = length(missing_fixef))
    rownames(zero_mat) <- missing_fixef
    fixef <- rbind(fixef, zero_mat)
  }
  for (i in seq_along(coef)) {
    missing_ranef <-  setdiff(rownames(fixef), colnames(coef[[i]]))
    if (length(missing_ranef)) {
      zero_mat <- matrix(0, nrow = nrow(coef[[i]]), 
                         ncol = length(missing_ranef))
      colnames(zero_mat) <- missing_ranef
      coef[[i]] <- cbind(coef[[i]], zero_mat)
    }
    for (nm in colnames(coef[[i]])) {
      coef[[i]][, nm] <- coef[[i]][, nm] + fixef[nm, 1]
    }
  }
  coef
}

#' Extract variance and correlation components
#' 
#' This function calculates the estimated standard deviations, 
#' correlations and covariances of the random-effects terms 
#' in a mixed-effects model of class \code{brmsfit}. 
#' For linear models, the residual standard deviations, 
#' correlations and covariances are also returned. 
#' 
#' @aliases VarCorr
#' 
#' @param x An object of class \code{brmsift}. 
#' @param estimate A character vector specifying which summary
#'  statistics (e.g., "mean", "median", "sd", or "quantile")
#'  should be calculated for the extracted parameters.
#' @param as.list logical; Indicates if covariance 
#'  and correlation matrices should be returned as 
#'  lists of matrices (the default), or as 3-dimensional arrays.
#'  We recommend not to set \code{as.list} to \code{FALSE}.
#' @param sigma Ignored (included for compatibility with 
#'  \code{\link[lme4:VarCorr]{VarCorr}}).
#' @param ... Further arguments to be passed to the functions 
#'  specified in \code{estimate}
#' 
#' @return An object of class \code{brmsVarCorr}, 
#' which is a list of lists (one per grouping factor), 
#' each containing 3 elements: a matrix containing the standard deviations, 
#' a list of correlation matrices, and a list of covariance matrices. 
#' Can be coerced to a \code{data.frame} by using the \code{as.data.frame} method.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' fit <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1+Trt_c|visit), 
#'              data = epilepsy, family = "poisson", chains = 1)
#' ## return the means of random effects covariances
#' (vc <- VarCorr(fit))
#' as.data.frame(vc)
#' 
#' ## return 2.5% and 97.5% quantiles of random effects covariances
#' VarCorr(fit, estimate = "quantile", probs = c(0.025, 0.975))
#' }
#' 
#' @import abind abind
#' @importFrom lme4 VarCorr
#' @export VarCorr
#' @export
VarCorr.brmsfit <- function(x, sigma = 1, estimate = "mean", 
                            as.list = TRUE, ...) {
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  if (!(length(x$ranef) || any(grepl("^sigma_", parnames(x)))))
    stop("The model does not contain covariance matrices", call. = FALSE)
  if (!is_equal(sigma, 1))
    warning("argument 'sigma' is unused")
  
  # extracts samples for sd, cor and cov
  extract <- function(p) {
    nr <- length(p$sd_pars)
    sd <- posterior_samples(x, pars = p$sd_pars, exact_match = TRUE)
    nsamples <- nrow(sd)
    out <- list(sd = do.call(cbind, lapply(estimate, get_estimate, 
                                           samples = sd, ...)))
    rownames(out$sd) <- p$rnames 
    
    # calculate correlation and covariance matrices
    found_cor_pars <- intersect(p$cor_pars, parnames(x))
    if (length(found_cor_pars)) {
      cor <- posterior_samples(x, pars = paste0("^",found_cor_pars,"$"))
      if (length(found_cor_pars) < length(p$cor_pars)) { 
        # some correlations are missing and will be replaced by 0
        cor_all <- as.data.frame(matrix(0, nrow = nrow(cor), 
                                        ncol = length(p$cor_pars)))
        names(cor_all) <- p$cor_pars
        for (i in 1:ncol(cor_all)) {
          found <- match(names(cor_all)[i], names(cor))
          if (!is.na(found))  # correlation was estimated
            cor_all[, i] <- cor[, found]
        }
        cor <- cor_all
      }
    } else cor <- NULL
    
    # get_cov_matrix and array2list can be found in brsmfit-helpers.R
    matrices <- get_cov_matrix(sd = sd, cor = cor) 
    out$cor <- abind(lapply(estimate, get_estimate, samples = matrices$cor, 
                            margin = c(2,3), to.array = TRUE, ...))
    out$cov <- abind(lapply(estimate, get_estimate, samples = matrices$cov, 
                            margin = c(2,3), to.array = TRUE, ...)) 
    if (length(p$rnames) > 1) {
      dimnames(out$cor) <- list(p$rnames, p$rnames, dimnames(out$cor)[[3]])
      dimnames(out$cov) <- dimnames(out$cor)   
    }
    if (as.list) {
      out$cor <- lapply(array2list(out$cor), function(x)
        if (is.null(dim(x))) 
          structure(matrix(x), dimnames = list(p$rnames, p$rnames)) 
        else x)
      out$cov <- lapply(array2list(out$cov), function(x)
        if (is.null(dim(x))) 
          structure(matrix(x), dimnames = list(p$rnames, p$rnames)) 
        else x)
    }
    out
  }
  
  family <- family(x)
  ee <- extract_effects(x$formula, family = family, nonlinear = x$nonlinear)
  if (length(x$ranef)) {
    gather_names <- function(i) {
      # gather names of random effects parameters
      cor_type <- paste0("cor_", group[i])
      sd_pars <- paste0("sd_", group[i], "_", x$ranef[[i]])
      cor_pars <- get_cornames(x$ranef[[i]], type = cor_type, brackets = FALSE)
      nlist(rnames = x$ranef[[i]], type = cor_type, sd_pars, cor_pars)
    }
    group <- paste0(ulapply(x$ranef, get_nlpar, suffix = "_"), names(x$ranef))
    p <- lapply(seq_along(group), gather_names)
  } else {
    p <- group <- NULL
  } 
  # special treatment of residuals variances in linear models
  if (has_sigma(family, se = ee$se, autocor = x$autocor)) {
    cor_pars <- get_cornames(ee$response, type = "rescor", 
                             brackets = FALSE)
    p <- lc(p, list(rnames = ee$response, 
                    sd_pars = paste0("sigma_", ee$response),
                    cor_pars = cor_pars))
    group <- c(group, "RESIDUAL")
  } 
  VarCorr <- lapply(p, extract)
  names(VarCorr) <- group
  if (as.list) class(VarCorr) <- "brmsVarCorr"
  VarCorr
}

#' @export
model.frame.brmsfit <- function(formula, ...) {
  if (!is.data.frame(formula$data)) {
    stop("Cannot extract model.frame for models fitted with brms <= 0.5.0.",
         call. = FALSE)
  }
  formula$data 
}

#' @rdname posterior_samples
#' @export
posterior_samples.brmsfit <- function(x, pars = NA, parameters = NA,  
                                      exact_match = FALSE, 
                                      add_chain = FALSE,
                                      add_chains = FALSE, 
                                      subset = NULL, as.matrix = FALSE, 
                                      ...) {
  if (is.na(pars[1])) 
    pars <- parameters  
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  pars <- extract_pars(pars, all_pars = parnames(x), 
                       exact_match = exact_match, ...)
  
  # get basic information on the samples 
  iter <- x$fit@sim$iter
  warmup <- x$fit@sim$warmup
  thin <- x$fit@sim$thin
  chains <- x$fit@sim$chains
  final_iter <- (iter - warmup) / thin
  samples_taken <- seq(warmup + 1, iter, thin)
  
  if (length(pars)) {
    samples <- as.data.frame(x$fit, pars = pars)
    if (add_chain) {
      # name the column 'chain' not 'chains' (#32)
      samples$chain <- factor(rep(1:chains, each = final_iter))
      samples$iter <- rep(samples_taken, chains)
    }
    if (add_chains) {
      warning(paste("Argument 'add_chains' is deprecated.",
                    "Please use argument 'add_chain' instead."))
      if (!add_chain) {
        samples$chains <- factor(rep(1:chains, each = final_iter))
        samples$iter <- rep(samples_taken, chains)
      }
    }
    if (!is.null(subset)) {
      samples <- samples[subset, , drop = FALSE]
    }
    if (as.matrix) {
      samples <- as.matrix(samples)
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

#' Extract posterior samples for use with the \pkg{coda} package
#' 
#' @aliases as.mcmc
#' 
#' @inheritParams posterior_samples
#' @param ... currently unused
#' @param inc_warmup Indicates if the warmup samples should be included.
#'   Default is \code{FALSE}. Warmup samples are used to tune the 
#'   parameters of the sampling algorithm and should not be analyzed.
#'   
#' @return A \code{list} of \code{mcmc} objects (not an \code{mcmc} object itself).
#' 
#' @method as.mcmc brmsfit
#' @export
#' @export as.mcmc
#' @importFrom coda as.mcmc
as.mcmc.brmsfit <- function(x, pars = NA, exact_match = FALSE,
                                 inc_warmup = FALSE, ...) {
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  pars <- extract_pars(pars, all_pars = parnames(x),
                       exact_match = exact_match, ...)
  ps <- extract(x$fit, pars = pars, permuted = FALSE, 
                inc_warmup = inc_warmup)
  first <- if (inc_warmup) 1 else x$fit@sim$warmup + 1
  out <- vector("list", length = dim(ps)[2])
  for (i in seq_along(out)) {
    out[[i]] <- ps[, i, ]
    attr(out[[i]], "mcpar") <- c(first, x$fit@sim$iter, x$fit@sim$thin)
    class(out[[i]]) <- "mcmc" 
  }
  class(out) <- "mcmc.list"
  out
}

#' @rdname prior_samples
#' @export
prior_samples.brmsfit <- function(x, pars = NA, parameters = NA, ...) {
  if (is.na(pars[1])) 
    pars <- parameters 
  if (!anyNA(pars) && !is.character(pars)) 
    stop("pars must be a character vector", call. = FALSE)
  par_names <- parnames(x)
  prior_names <- par_names[grepl("^prior_", par_names)]
  if (length(prior_names)) {
    samples <- posterior_samples(x, pars = prior_names, 
                                 exact_match = TRUE)
    names(samples) <- sub("^prior_", "", prior_names)
    if (!anyNA(pars)) {
      get_samples <- function(par) {
        # get prior samples for parameter par 
        is_cse <- grepl("^b_", par) && grepl("\\[[[:digit:]]+\\]", par)
        # ensures correct parameter to prior mapping for cse effects
        par_internal <- ifelse(is_cse, sub("^b_", "bp_", par), par)
        matches <- lapply(paste0("^",sub("^prior_", "", prior_names)), 
                          regexpr, text = par_internal)
        matches <- ulapply(matches, attr, which = "match.length")
        if (max(matches) == -1) {
          return(NULL)
        } else {
          take <- match(max(matches), matches)
          # order samples randomly to avoid artifical dependencies
          # between parameters using the same prior samples
          samples <- list(samples[sample(Nsamples(x)), take])
          return(structure(samples, names = par))
        }
      }
      samples <- data.frame(rmNULL(lapply(pars, get_samples)), 
                            check.names = FALSE)
    }
  } else {
    samples <- NULL
  }
  samples
}

#' Print a summary for a fitted model represented by a \code{brmsfit} object
#'
#' Print basic information regarding the fitted model and a summary 
#' for the fixed and random effects
#' estimated by the samples included in a \code{brmsfit} object.
#' 
#' @aliases print.brmssummary
#' 
#' @param x An object of class \code{brmsfit}
#' @param digits The number of significant digits for printing out the summary; 
#'  defaults to 2. The effective sample size is always rounded to integers.
#' @param ... Additional arguments that would be passed 
#'  to method \code{summary} of \code{brmsfit}.
#'
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @export
print.brmsfit <- function(x, digits = 2, ...) {
  print(summary(x), digits = digits, ...)
}  

#' Create a summary of a fitted model represented by a \code{brmsfit} object
#' 
#' Summarize estimated fixed and random effects as well as other useful
#' results included in a \code{brmsfit} object.
#'
#' @param object An object of class \code{brmsfit}
#' @param waic logical; indicating if the WAIC should be computed
#'   (this will take some time for larger models). 
#'   Default is \code{FALSE}.
#' @param ... Other potential arguments
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @method summary brmsfit
#' @export
summary.brmsfit <- function(object, waic = FALSE, ...) {
  family <- family(object)
  ee <- extract_effects(object$formula, family = family,
                        nonlinear = object$nonlinear)
  formula <- SW(update_formula(object$formula, partial = object$partial))
  out <- brmssummary(formula = formula,
                     family = family, 
                     data.name = object$data.name, 
                     group = names(object$ranef), 
                     nobs = nobs(object), 
                     ngrps = ngrps(object), 
                     autocor = object$autocor,
                     nonlinear = object$nonlinear,
                     algorithm = algorithm(object))
  
  if (length(object$fit@sim)) {
    out$chains <- object$fit@sim$chains
    out$iter <- object$fit@sim$iter
    out$warmup <- object$fit@sim$warmup
    out$thin <- object$fit@sim$thin
    stan_args <- object$fit@stan_args[[1]]
    out$sampler <- paste0(stan_args$method, "(", stan_args$algorithm, ")")
    allow_waic <- !length(object$ranef) || any(grepl("^r_", parnames(object)))
    if (waic && allow_waic) out$WAIC <- WAIC(object)$waic
    
    pars <- parnames(object)
    meta_pars <- object$fit@sim$pars_oi
    meta_pars <- meta_pars[!apply(sapply(paste0("^", c("r_", "prior_")), 
                                  grepl, x = meta_pars, ...), 1, any)]
    fit_summary <- rstan::summary(object$fit, pars = meta_pars,
                                  probs = c(0.025, 0.975))$summary
    algorithm <- algorithm(object)
    if (algorithm == "sampling") {
      fit_summary <- fit_summary[, -2]
      colnames(fit_summary) <- c("Estimate", "Est.Error", "l-95% CI", 
                                 "u-95% CI", "Eff.Sample", "Rhat")
      Rhats <- fit_summary[, "Rhat"]
      if (any(Rhats > 1.1, na.rm = TRUE) || anyNA(Rhats)) {
        msg <- paste("The model has not converged (some Rhats are > 1.1).",
                     "Do not analyse the results! \nWe recommend running", 
                     "more iterations and/or setting stronger priors.")
        warning(msg, call. = FALSE)
      }
    } else {
      colnames(fit_summary) <- c("Estimate", "Est.Error", "l-95% CI", 
                                 "u-95% CI")
    }
    
    # fixed effects summary
    fix_pars <- pars[grepl("^b_", pars)]
    out$fixed <- fit_summary[fix_pars, , drop = FALSE]
    rownames(out$fixed) <- gsub("^b_", "", fix_pars)
    
    # summary of family specific parameters
    spec_pars <- pars[pars %in% c("nu", "shape", "delta", "phi") | 
      apply(sapply(c("^sigma_", "^rescor_"), grepl, x = pars), 1, any)]
    out$spec_pars <- fit_summary[spec_pars, , drop = FALSE]
    if (is.linear(family)) {
      sigma_names <- paste0("sigma(",ee$response,")")
      rescor_names <- get_cornames(ee$response, type = "rescor")   
      spec_pars[grepl("^sigma_", spec_pars)] <- sigma_names
      spec_pars[grepl("^rescor_", spec_pars)] <- rescor_names 
    }    
    rownames(out$spec_pars) <- spec_pars
    
    # summary of autocorrelation effects
    cor_pars <- pars[grepl("^ar|^ma|^sigmaLL$", pars)]
    out$cor_pars <- fit_summary[cor_pars, , drop = FALSE]
    rownames(out$cor_pars) <- cor_pars
    
    # summary of random effects
    for (i in seq_along(out$group)) {
      nlp <- get_nlpar(object$ranef[[i]])
      nlp_ <- ifelse(nchar(nlp), paste0(nlp, "_"), nlp)
      rnames <- object$ranef[[i]]
      sd_pars <- paste0("sd_", nlp_, out$group[i], "_", rnames)
      sd_names <- paste0("sd", "(", nlp_, rnames,")")
      # construct correlation names
      full_type <- paste0("cor_", nlp_, out$group[i])
      all_cor_pars <- get_cornames(rnames, brackets = FALSE, type = full_type)
      take <- all_cor_pars %in% parnames(object)
      cor_pars <- all_cor_pars[take]
      cor_names <- get_cornames(paste0(nlp_, rnames))[take]
      # extract sd and cor parameters from the summary
      new_random <- fit_summary[c(sd_pars, cor_pars), , drop = FALSE]
      rownames(new_random) <- c(sd_names, cor_names)
      out$random[[out$group[i]]] <- 
        rbind(out$random[[out$group[i]]], new_random)
    }
    
    # summary of splines
    spline_pars <- pars[grepl("^sds_", pars)]
    if (length(spline_pars)) {
      out$splines <- fit_summary[spline_pars, , drop = FALSE]
      rownames(out$splines) <- paste0(gsub("^sds_", "sds(", spline_pars), ")")
    }
  }  
  out
}

#' @export
nobs.brmsfit <- function(object, ...) {
  length(standata(object)$Y)
}

#' Number of levels
#' 
#' Number of levels of one or more grouping factors
#' 
#' @aliases ngrps
#' 
#' @param object An object of class \code{brmsfit}.
#' @param ... Currently ignored.
#' 
#' @return A named list containing the number of levels per
#'   grouping factor
#' 
#' @export
#' @export ngrps
#' @importFrom lme4 ngrps
ngrps.brmsfit <- function(object, ...) {
  if (length(object$ranef)) {
    out <- named_list(names(object$ranef))
    for (i in seq_along(out)) {
      out[[i]] <- length(attr(object$ranef[[i]], "levels"))
    } 
    out <- out[!duplicated(names(out))]
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
family.brmsfit <- function(object, ...) {
  if (is(object$family, "family")) {
    # brms > 0.6.0
    family <- object$family
  } else {
    family <- family(object$family, link = object$link) 
  }
  family
}

#' @export
stancode.brmsfit <- function(object, ...) {
  object$model
}

#' @export
standata.brmsfit <- function(object, ...) {
  dots <- list(...)
  if (is.data.frame(object$data)) {
    # brms > 0.5.0 stores the original model.frame
    new_formula <- update_re_terms(object$formula, dots$re_formula)
    new_formula <- SW(update_formula(new_formula, partial = object$partial))
    dots$control$old_cat <- is.old_categorical(object)
    prior_only <- attr(object$prior, "prior_only")
    sample_prior <- ifelse(isTRUE(prior_only), "only", FALSE)
    args <- list(formula = new_formula, data = model.frame(object), 
                 family = object$family, prior = object$prior, 
                 nonlinear = object$nonlinear, autocor = object$autocor, 
                 cov_ranef = object$cov_ranef, sample_prior = sample_prior,
                 knots = attr(model.frame(object), "knots"))
    standata <- do.call(make_standata, c(args, dots))
  } else {
    # brms <= 0.5.0 only stores the data passed to Stan 
    standata <- object$data
    # for a short period in 0.4.1.9000, "lev" was used instead of "J"
    names(standata) <- sub("^lev_", "J_", names(standata))
  }
  standata
}
  
#' @export
launch_shiny.brmsfit <- function(x, rstudio = getOption("shinystan.rstudio"), 
                                 ...) {
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  shinystan::launch_shinystan(x$fit, rstudio = rstudio, ...)
}

#' Trace and Density Plots for MCMC Samples
#' 
#' @param x An object of class \code{brmsfit}.
#' @param pars Names of the parameters to plot, as given by a character vector 
#'   or a regular expression. By default, all parameters except for random effects 
#'   are plotted. 
#' @param parameters A deprecated alias of \code{pars}   
#' @param N The number of parameters plotted per page.
#' @param theme A \code{\link[ggplot2:theme]{theme}} object 
#'   modifying the appearance of the plots. 
#'   For some basic themes see \code{\link[ggplot2:ggtheme]{ggtheme}}. 
#'   Can be defined globally for the current session, via
#'   \code{\link[ggplot2:theme_update]{theme_set}}.
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
#'   \code{\link[gridExtra:arrangeGrob]{arrangeGrob}}.
#' 
#' @return An invisible list of 
#'   \code{\link[gtable:gtable]{gtable}} objects.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' \dontrun{ 
#' fit <- brm(count ~ log_Age_c + log_Base4_c * Trt_c 
#'            + (1|patient) + (1|visit), 
#'            data = epilepsy, family = "poisson")
#' ## plot fixed effects as well as standard devations of the random effects
#' plot(fit)
#' ## plot fixed effects only
#' plot(fit, pars = "^b_") 
#' }
#' 
#' @method plot brmsfit
#' @import ggplot2
#' @importFrom gridExtra arrangeGrob
#' @importFrom grDevices devAskNewPage
#' @importFrom grid grid.draw grid.newpage
#' @importFrom bayesplot mcmc_trace mcmc_dens_overlay
#' @export
plot.brmsfit <- function(x, pars = NA, parameters = NA, N = 5, 
                         theme = bayesplot::theme_ppc(), ask = TRUE, 
                         plot = TRUE, newpage = TRUE, ...) {
  dots <- list(...)
  plot <- use_alias(plot, dots$do_plot)
  dots$do_plot <- NULL
  if (is.na(pars[1])) 
    pars <- parameters 
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  if (!is.wholenumber(N) || N < 1)
    stop("N must be a positive integer", call. = FALSE)
  if (!is.character(pars)) {
    pars <- c("^b_", "^bm_", "^sd_", "^cor_", "^sigma", "^rescor", 
              "^nu$", "^shape$", "^delta$", "^phi$", "^ar", "^ma", 
              "^arr", "^simplex_", "^sds_")
  }
  samples <- posterior_samples(x, pars = pars, add_chain = TRUE)
  pars <- names(samples)[!names(samples) %in% c("chain", "iter")] 
  if (length(pars) == 0) {
    stop("No valid parameters selected", call. = FALSE)
  }
  
  if (plot) {
    default_ask <- devAskNewPage()
    on.exit(devAskNewPage(default_ask))
    devAskNewPage(ask = FALSE)
  }
  n_plots <- ceiling(length(pars) / N)
  plots <- vector(mode = "list", length = n_plots)
  for (i in 1:n_plots) {
    sub_pars <- pars[((i - 1) * N + 1):min(i * N, length(pars))]
    sub_samples <- samples[, c(sub_pars, "chain"), drop = FALSE]
    bp_args <- list(sub_samples, facet_args = list(ncol = 1), plot = FALSE)
    trace <- do.call(mcmc_trace, bp_args) + xlab("") + ylab("") + 
      theme + theme(legend.position = "none")
    dens <- do.call(mcmc_dens_overlay, bp_args) + xlab("") + ylab("") + 
      theme + theme(legend.position = "right")
    ge_args <- list(trace, dens, nrow = 1, ncol = 2, widths = c(1, 1.2))
    plots[[i]] <- do.call(arrangeGrob, c(ge_args, dots))
    if (plot) {
      if (newpage || i > 1) grid.newpage()
      grid.draw(plots[[i]])
      if (i == 1) devAskNewPage(ask = ask)
    }
  }
  invisible(plots) 
}

#' @rdname stanplot
#' @export
stanplot.brmsfit <- function(object, pars = NA, type = "plot", 
                             exact_match = FALSE, quiet = FALSE, ...) {
  
  # check validity of type first
  basic_types <- c("plot", "trace", "scat", "hist", "dens", "ac")
  diag_types <- c("diag", "par", "rhat", "ess", "mcse")
  if (!type %in% c(basic_types, diag_types)) {
    stop(paste("Invalid plot type. Valid plot types are: \n",
               paste(c(basic_types, diag_types), collapse = ", ")),
         call. = FALSE)
  }
  dots <- list(...)
  args <- c(object = object$fit, dots)
  plot_fun <- get(paste0("stan_", type), mode = "function")
  
  # ensure that only desired parameters are plotted
  if (!is(object$fit, "stanfit") || !length(object$fit@sim)) 
    stop("The model does not contain posterior samples")
  pars <- extract_pars(pars, all_pars = parnames(object),
                       exact_match = exact_match, 
                       na_value = NA)
  
  if (type %in% basic_types) {
    if (!anyNA(pars)) {
      args <- c(args, list(pars = pars))
    }
  } else if (type %in% diag_types) {
    if (type %in% c("rhat", "ess", "msce") && !anyNA(pars)) {
      args <- c(args, list(pars = pars))
    } else if (type == "par") {
      if (length(pars) > 1L) {
        warning(paste("stan_par expects a single parameter name",
                      "so that only the first one will be used"))
      }
      args <- c(args, list(par = pars[1]))
    } 
  }
  # make the plot
  if (quiet) {
    suppressMessages(do.call(plot_fun, args))
  } else {
    do.call(plot_fun, args)
  }
}

#' Posterior Predictive Checks for \code{brmsfit} Objects
#' 
#' Perform posterior predictive checks with the help
#' of the \pkg{bayesplot} package.
#' 
#' @param object An object of class \code{brmsfit}.
#' @param type Type of the ppc plot as given by a character string. 
#'   Currently, the following plots (as names) are implemented:
#'   \code{dens} \code{dens_overlay}, \code{hist}, \code{resid}, 
#'   \code{resid_binned}, \code{scatter}, \code{scatter_avg}, 
#'   \code{scatter_avg_grouped}, \code{stat}, \code{stat_2d},
#'   \code{stat_grouped}, \code{ts}, \code{ts_grouped}, and
#'   \code{violin_grouped}.
#' @param nsamples Positive integer indicating how many 
#'  posterior samples should be used. 
#'  If \code{NULL} all samples are used. If not specified, 
#'  the number of posterior samples is chosen automatically.
#'  Ignored if \code{subset} is not \code{NULL}.
#' @param ntrys Parameter used in rejection sampling 
#'  for truncated discrete models only 
#'  (defaults to \code{5}). For more details see
#'  \code{\link[brms:predict.brmsfit]{predict.brmsfit}}.
#' @param group Optional name of a grouping factor in the model
#'  by which to stratify the ppc plot. This argument is required for 
#'  ppc \code{*_grouped} types and ignored otherwise.
#' @param time Optional name of a time variable in the model by
#'  which to order time-series plots. Only used for
#'  ppc \code{ts*} types and ignored otherwise.
#' @param ... Further arguments passed to the ppc functions
#'   of the \pkg{bayesplot} package.
#' @inheritParams predict.brmsfit
#' 
#' @return A ggplot object that can be further 
#'  customized using the \pkg{ggplot2} package.
#'  
#' @details For a detailed explanation of each of the
#' ppc functions, see \code{\link[bayesplot:PPC-overview]{PPC-overview}}.
#' 
#' @examples
#' \dontrun{
#' fit <-  brm(count ~ log_Age_c + log_Base4_c * Trt_c 
#'             + (1|patient) + (1|visit),
#'             data = epilepsy, family = poisson())
#' 
#' pp_check(fit) # shows dens_overlay plot by default             
#' pp_check(fit, type = "resid", nsamples = 12)
#' pp_check(fit, type = "scatter_average", nsamples = 100)  
#' pp_check(fit, type = "stat_2d")
#' } 
#' 
#' @export
pp_check.brmsfit <- function(object, type, nsamples, group = NULL,
                             time = NULL, re_formula = NULL, 
                             subset = NULL, ntrys = 5, ...) {
  if (!requireNamespace("bayesplot", quietly = TRUE)) {
    # remove check as soon as bayesplot is on CRAN
    stop(paste0("please install the bayesplot package via\n",
                "devtools::install_github('jgabry/bayesplot')"),
         call. = FALSE)
  }
  if (missing(type)) {
    type <- "dens_overlay"
  }
  if (length(type) != 1L) {
    stop("argument 'type' must be of length 1", call. = FALSE)
  }
  if (!is.null(group) && length(group) != 1L) {
    stop("argument 'group' must be of length 1", call. = FALSE)
    
  }
  if (!is.null(time) && length(time) != 1L) {
    stop("argument 'time' must be of length 1", call. = FALSE)
    
  }
  ppc_funs <- lsp("bayesplot", what = "exports", pattern = "^ppc_")
  valid_ppc_types <- sub("^ppc_", "", ppc_funs)
  if (!type %in% valid_ppc_types) {
    stop(paste0("Type '", type, "' is not a valid ppc type. ",
                "Valid types are: \n", 
                paste(valid_ppc_types, collapse = ", ")),
         call. = FALSE)
  }
  ppc_fun <- get(paste0("ppc_", type), pos = asNamespace("bayesplot"))
  # validate argument "group"
  valid_groups <- names(object$ranef)
  time_group <- extract_time(object$autocor$formula)$group
  if (!is.null(time_group) && nchar(time_group)) {
    valid_groups <- unique(c(valid_groups, time_group)) 
  }
  if (!is.null(group) && !group %in% valid_groups) {
    stop(paste0("Group '", group, "' is not a valid grouping factor. ",
                "Valid groups are: \n", paste(valid_groups, collapse = ", ")),
         call. = FALSE)
  }
  is_group_type <- "group" %in% names(formals(ppc_fun))
  if (is.null(group) && is_group_type) {
    stop(paste0("Argument 'group' is required for ppc type '", type, "'."), 
         call. = FALSE)
  }
  if (!is.null(group) && !is_group_type) {
    warning(paste0("Argument 'group' is ignored for ppc type '", type, "'."), 
            call. = FALSE)
  }
  is_ts_type <- "time" %in% names(formals(ppc_fun))
  if (!is.null(time) && !is_ts_type) {
    warning(paste0("Argument 'time' is ignored for ppc type '", type, "'."), 
            call. = FALSE)
  }
  if (names(formals(ppc_fun))[2] == "Ey") {
    if (is.ordinal(object$family)) {
      stop(paste0("ppc type '", type, "' is not available", 
                  "for ordinal models"), call. = FALSE)
    }
    method <- "fitted"
  } else {
    method <- "predict"
  }
  if (missing(nsamples)) {
    aps_types <- c("scatter_avg", "scatter_avg_grouped", "stat", "stat_2d", 
                   "stat_grouped", "ts", "ts_grouped", "violin_grouped")
    if (!is.null(subset)) {
      nsamples <- NULL
    } else if (type %in% aps_types) {
      nsamples <- NULL
      message(paste0("Using all posterior samples for ppc type '", 
                     type, "'"))
    } else {
      nsamples <- 10
      message(paste0("Using 10 posterior samples for ppc type '", 
                     type, "'"))
    }
  }
  args <- nlist(object, nsamples, subset, re_formula, 
                ntrys, sort = TRUE, summary = FALSE)
  yrep <- as.matrix(do.call(method, args))
  standata <- standata(object, control = list(save_order = TRUE))
  y <- as.vector(standata$Y)
  if (family(object)$family %in% "binomial") {
    # use success proportions following Gelman and Hill (2006)
    y <- y / standata$trials
    yrep <- yrep / matrix(standata$trials, nrow = nrow(yrep),
                          ncol = ncol(yrep), byrow = TRUE)
  }
  ppc_args <- list(y, yrep, ...)
  if (!is.null(group)) {
    group <- model.frame(object)[[group]]
    if (!is.null(attr(standata, "old_order"))) {
      group <- group[order(attr(standata, "old_order"))]
    }
    ppc_args$group <- group
  }
  if (!is.null(time)) {
    time <- model.frame(object)[[time]]
    if (!is.null(attr(standata, "old_order"))) {
      time <- time[order(attr(standata, "old_order"))]
    }
    ppc_args$time <- as.numeric(time)
  }
  rstan::quietgg(do.call(ppc_fun, ppc_args))
}

#' Create a matrix of output plots from a \code{brmsfit} object
#'
#' A \code{\link[graphics:pairs]{pairs}} 
#' method that is customized for MCMC output
#' 
#' @param x An object of class \code{stanfit}
#' @param pars Names of the parameters to plot, as given by 
#'  a character vector or a regular expression. 
#'  By default, all parameters are plotted. 
#' @param exact_match Indicates whether parameter names 
#'   should be matched exactly or treated as regular expression. 
#'   Default is \code{FALSE}.
#' @param ... Further arguments to be passed to 
#'  \code{\link[rstan:pairs.stanfit]{pairs.stanfit}}.
#'  
#' @details For a detailed description see 
#'  \code{\link[rstan:pairs.stanfit]{pairs.stanfit}}.
#'  
#' @examples 
#' \dontrun{
#' fit <- brm(count ~ log_Age_c + log_Base4_c * Trt_c 
#'            + (1|patient) + (1|visit), 
#'            data = epilepsy, family = "poisson")  
#' pairs(fit, pars = parnames(fit)[1:3], exact_match = TRUE)
#' pairs(fit, pars = "^sd")
#' }
#'
#' @export
pairs.brmsfit <- function(x, pars = NA, exact_match = FALSE, ...) {
  pars <- extract_pars(pars, all_pars = parnames(x),
                       exact_match = exact_match, na_value = NULL)
  graphics::pairs(x$fit, pars = pars, ...)
}

#' @rdname marginal_effects
#' @export
marginal_effects.brmsfit <- function(x, effects = NULL, conditions = NULL, 
                                     re_formula = NA, robust = FALSE, 
                                     probs = c(0.025, 0.975),
                                     method = c("fitted", "predict"), ...) {
  method <- match.arg(method)
  dots <- list(...)
  conditions <- use_alias(conditions, dots$data)
  x$formula <- SW(update_formula(x$formula, partial = x$partial))
  new_formula <- update_re_terms(x$formula, re_formula = re_formula)
  new_nonlinear <- lapply(x$nonlinear, update_re_terms, re_formula = re_formula)
  ee <- extract_effects(new_formula, family = x$family, 
                        nonlinear = new_nonlinear)
  if (is.linear(x$family) && length(ee$response) > 1L) {
    stop("Marginal plots are not yet implemented for multivariate models.",
         call. = FALSE)
  } else if (is.categorical(x$family)) {
    stop("Marginal plots are not yet implemented for categorical models.",
         call. = FALSE)
  } else if (is.ordinal(x$family)) {
    warning(paste0("Predictions are treated as continuous variables ", 
                   "in marginal plots, \nwhich is likely an invalid ", 
                   "assumption for family ", x$family$family, "."),
            call. = FALSE)
  }
  rsv_vars <- rsv_vars(x$family, nresp = length(ee$response),
                       rsv_intercept = attr(ee$fixed, "rsv_intercept"))
  if (length(ee$nonlinear)) {
    # allow covariates as well as fixed effects of non-linear parameters
    covars <- setdiff(all.vars(rhs(ee$fixed)), names(ee$nonlinear))
    nlpar_effects <- unlist(lapply(ee$nonlinear, function(nl)
      get_var_combs(nl$fixed, nl$mono, nl$gam)), recursive = FALSE)
    all_effects <- unique(c(list(covars), nlpar_effects))
  } else {
    all_effects <- get_var_combs(ee$fixed, ee$cse, ee$mono, ee$gam)
  }
  all_effects <- rmNULL(lapply(all_effects, setdiff, y = rsv_vars))
  ae_collapsed <- ulapply(all_effects, function(e) paste(e, collapse = ":"))
  if (is.null(effects)) {
    effects <- all_effects[ulapply(all_effects, length) < 3]
  } else {
    # allow to define interactions in any order
    effects <- strsplit(effects, split = ":")
    if (any(unique(unlist(effects)) %in% rsv_vars)) {
      stop(paste("Variables", paste0(rsv_vars, collapse = ", "),
                 "should not be used as effects for this model"),
           call. = FALSE)
    }
    matches <- match(lapply(all_effects, sort), lapply(effects, sort), 0L)
    if (sum(matches) > 0 && sum(matches > 0) < length(effects)) {
      invalid <- effects[setdiff(1:length(effects), sort(matches))]  
      invalid <- ulapply(invalid, function(e) paste(e, collapse = ":"))
      warning(paste0("Some specified effects are invalid for this model: ",
                     paste(invalid, collapse = ", "), "\nValid effects are: ", 
                     paste(ae_collapsed, collapse = ", ")),
              call. = FALSE)
    }
    effects <- unique(effects[sort(matches)])
  }
  if (!length(unlist(effects))) {
    stop(paste0("All specified effects are invalid for this model.\n", 
                "Valid effects are: ", paste(ae_collapsed, collapse = ", ")), 
         call. = FALSE)
  }
  if (any(ulapply(effects, length) > 2L)) {
    stop("Interactions of order higher than 2 are currently not supported.",
         call. = FALSE)
  }
  if (length(probs) != 2L) {
    stop("Arguments 'probs' must be of length 2.", call. = FALSE)
  }
  
  # prepare marginal conditions
  mf <- model.frame(x)
  mono_vars <- unique(ulapply(get_effect(ee, "mono"), all.vars))
  if (is.null(conditions)) {
    if (length(rmNULL(ee[c("trials", "cat")]))) {
      stop("Please specify argument 'conditions' manually for this model.", 
           call. = FALSE)
    }
    # list all required variables
    req_vars <- c(lapply(get_effect(ee), rhs), get_random(ee)$form, 
                  get_effect(ee, "mono"), get_effect(ee, "gam"), ee$cse,
                  ee$se, ee$disp, extract_time(x$autocor$formula)$all)
    req_vars <- unique(ulapply(req_vars, all.vars))
    req_vars <- setdiff(req_vars, c(rsv_vars, names(ee$nonlinear)))
    conditions <- as.data.frame(as.list(rep(NA, length(req_vars))))
    names(conditions) <- req_vars
    for (v in req_vars) {
      if (is.numeric(mf[[v]])) {
        if (v %in% mono_vars) {
          conditions[[v]] <- round(median(mf[[v]]))
        } else {
          conditions[[v]] <- mean(mf[[v]])
        }
      } else {
        # use reference category
        lev <- attr(as.factor(mf[[v]]), "levels")
        conditions[[v]] <- factor(lev[1], levels = lev, 
                                  ordered = is.ordered(mf[[v]]))
      }
    }
  } else if (is.data.frame(conditions)) {
    if (!nrow(conditions)) {
      stop("'conditions' must have a least one row", call. = FALSE)
    }
    if (any(duplicated(rownames(conditions)))) {
      stop("Row names of 'conditions' should be unique.", call. = FALSE)
    }
    conditions <- unique(conditions)
    eff_vars <- lapply(effects, function(e) all.vars(parse(text = e)))
    uni_eff_vars <- unique(unlist(eff_vars))
    is_everywhere <- ulapply(uni_eff_vars, function(uv)
      all(ulapply(eff_vars, function(vs) uv %in% vs)))
    # variables that are present in every effect term
    # do not need to be defined in conditions
    missing_vars <- setdiff(uni_eff_vars[is_everywhere], names(conditions))
    for (v in missing_vars) {
      conditions[, v] <- mf[[v]][1]
    }
  } else {
    stop("conditions must be a data.frame or NULL")
  }
  conditions <- amend_newdata(conditions, fit = x, re_formula = re_formula,
                              allow_new_levels = TRUE, return_standata = FALSE)

  results <- list()
  for (i in seq_along(effects)) {
    marg_data <- mf[, effects[[i]], drop = FALSE]
    pred_types <- ifelse(ulapply(marg_data, is.numeric), "numeric", "factor")
    # numeric effects should come first
    new_order <- order(pred_types, decreasing = TRUE)
    effects[[i]] <- effects[[i]][new_order]
    pred_types <- pred_types[new_order]
    is_mono <- effects[[i]] %in% mono_vars
    if (pred_types[1] == "numeric") {
      min1 <- min(marg_data[, effects[[i]][1]])
      max1 <- max(marg_data[, effects[[i]][1]])
      if (is_mono[1]) {
        values <- seq(min1, max1, by = 1)
      } else {
        values <- seq(min1, max1, length.out = 100)
      }
    }
    if (length(effects[[i]]) == 2L) {
      if (pred_types[1] == "numeric") {
        values <- setNames(list(values, NA), effects[[i]])
        if (pred_types[2] == "numeric") {
          if (is_mono[2]) {
            median2 <- median(marg_data[, effects[[i]][2]])
            mad2 <- mad(marg_data[, effects[[i]][2]])
            values[[2]] <- round((-1:1) * mad2 + median2)
          } else {
            mean2 <- mean(marg_data[, effects[[i]][2]])
            sd2 <- sd(marg_data[, effects[[i]][2]])
            values[[2]] <- (-1:1) * sd2 + mean2
          }
        } else {
          values[[2]] <- unique(marg_data[, effects[[i]][2]])
        }
        marg_data <- do.call(expand.grid, values)
      }
    } else if (pred_types == "numeric") {
      # just a single numeric predictor
      marg_data <- structure(data.frame(values), names = effects[[i]])
    }
    # no need to have the same value combination more than once
    marg_data <- unique(marg_data)
    marg_data <- replicate(nrow(conditions), simplify = FALSE,
     expr = marg_data[do.call(order, as.list(marg_data)), , drop = FALSE])
    marg_vars <- setdiff(names(conditions), effects[[i]])
    for (j in 1:nrow(conditions)) {
      marg_data[[j]][, marg_vars] <- conditions[j, marg_vars]
      marg_data[[j]][["MargCond"]] <- rownames(conditions)[j]
    }
    marg_data <- do.call(rbind, marg_data)
    marg_data$MargCond <- factor(marg_data$MargCond, rownames(conditions))
    args <- list(x, newdata = marg_data, re_formula = re_formula,
                 allow_new_levels = TRUE, incl_autocor = FALSE,
                 probs = probs, robust = robust)
    if (is.ordinal(x$family) || is.categorical(x$family)) {
      args$summary <- FALSE 
      marg_res <- do.call(method, args)
      if (method == "fitted") {
        for (k in 1:dim(marg_res)[3]) {
          marg_res[, , k] <- marg_res[, , k] * k
        }
        marg_res <- do.call(cbind, lapply(1:dim(marg_res)[2], 
                            function(s) rowSums(marg_res[, s, ])))
      } 
      marg_res <- get_summary(marg_res, probs = probs, robust = robust)
    } else {
      marg_res <- do.call(method, args)
    }
    colnames(marg_res)[3:4] <- c("lowerCI", "upperCI")
     
    if (length(effects[[i]]) == 2L && all(pred_types == "numeric")) {
      # can only be converted to factor after having called method
      if (is_mono[2]) {
        labels <- c("Median - MAD", "Median", "Median + MAD")
      } else {
        labels <- c("Mean - SD", "Mean", "Mean + SD") 
      }
      marg_data[[effects[[i]][2]]] <- 
        factor(marg_data[[effects[[i]][2]]], labels = labels)
    }
    marg_res = cbind(marg_data, marg_res)
    attr(marg_res, "response") <- as.character(x$formula[2])
    attr(marg_res, "effects") <- effects[[i]]
    point_args <- nlist(mf, effects = effects[[i]], conditions,
                        groups = get_random(ee)$group, family = x$family)
    # see brmsfit-helpers.R
    attr(marg_res, "points") <- do.call(make_point_frame, point_args)
    results[[paste0(effects[[i]], collapse = ":")]] <- marg_res
  }
  class(results) <- "brmsMarginalEffects"
  results
}

#' Model Predictions of \code{brmsfit} Objects
#' 
#' Predict responses based on the fitted model.
#' Can be performed for the data used to fit the model 
#' (posterior predictive checks) or for new data.
#' By definition, these predictions have higher variance than 
#' predictions of the fitted values (i.e. the 'regression line')
#' performed by the \code{\link[brms:fitted.brmsfit]{fitted}}
#' method. This is because the measurement error is incorporated.
#' The estimated means of both methods should, however, be very similar.
#' 
#' @param object An object of class \code{brmsfit}
#' @param newdata An optional data.frame for which to evaluate predictions.
#'   If \code{NULL} (default), the orginal data of the model is used.
#' @param re_formula formula containing random effects 
#'   to be considered in the prediction. 
#'   If \code{NULL} (default), include all random effects; 
#'   if \code{NA}, include no random effects.
#' @param transform A function or a character string naming 
#'   a function to be applied on the predicted responses
#'   before summary statistics are computed.
#' @param allow_new_levels A flag indicating if new
#'   levels of random effects are allowed (defaults to \code{FALSE}). 
#'   Only relevant if \code{newdata} is provided.
#' @param incl_autocor A flag indicating if autocorrelation
#'  parameters should be included in the predictions. 
#'  Defaults to \code{TRUE}.
#' @param summary Should summary statistics 
#'   (i.e. means, sds, and 95\% intervals) be returned
#'  instead of the raw values? Default is \code{TRUE}.
#' @param robust If \code{FALSE} (the default) the mean is used as 
#'  the measure of central tendency and the standard deviation as 
#'  the measure of variability. If \code{TRUE}, the median and the 
#'  median absolute deivation (MAD) are applied instead.
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
#' @param ... Currently ignored
#' 
#' @return Predicted values of the response variable. 
#'   If \code{summary = TRUE} the output depends on the family:
#'   For catagorical and ordinal families, it is a N x C matrix, 
#'   where N is the number of observations and
#'   C is the number of categories. 
#'   For all other families, it is a N x E matrix where E is equal 
#'   to \code{length(probs) + 2}.
#'   If \code{summary = FALSE}, the output is as a S x N matrix, 
#'   where S is the number of samples.
#' 
#' @details \code{NA} values within factors in \code{newdata}, 
#'   are interpreted as if all dummy variables of this factor are 
#'   zero. This allows, for instance, to make predictions of the grand mean 
#'   when using sum coding.  
#' 
#'   For truncated discrete models only:
#'   In the absence of any general algorithm to sample 
#'   from truncated discrete distributions,
#'   rejection sampling is applied in this special case. 
#'   This means that values are sampled until 
#'   a value lies within the defined truncation boundaries. 
#'   In practice, this procedure may be rather slow (especially in R). 
#'   Thus, we try to do approximate rejection sampling 
#'   by sampling each value \code{ntrys} times and then select a valid value. 
#'   If all values are invalid, the closest boundary is used, instead. 
#'   If there are more than a few of these pathological cases, 
#'   a warning will occure suggesting to increase argument \code{ntrys}.
#'   
#'   For models fitted with \pkg{brms} <= 0.5.0 only: 
#'   Be careful when using \code{newdata} with factors 
#'   in fixed or random effects. The predicted results are only valid 
#'   if all factor levels present in the initial 
#'   data are also defined and ordered correctly 
#'   for the factors in \code{newdata}.
#'   Grouping factors may contain fewer levels than in the 
#'   inital data without causing problems.
#'   When using higher versions of \pkg{brms}, 
#'   all factors are automatically checked 
#'   for correctness and amended if necessary.
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
#' ## predicted responses excluding the random effect of age
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
#' @importFrom statmod rinvgauss pinvgauss qinvgauss
#' @export 
predict.brmsfit <- function(object, newdata = NULL, re_formula = NULL,
                            transform = NULL, allow_new_levels = FALSE,
                            incl_autocor = TRUE, subset = NULL, 
                            nsamples = NULL, sort = FALSE,
                            ntrys = 5, summary = TRUE, robust = FALSE,
                            probs = c(0.025, 0.975), ...) {
  if (!is(object$fit, "stanfit") || !length(object$fit@sim)) 
    stop("The model does not contain posterior samples")
  draws_args <- nlist(x = object, newdata, re_formula, incl_autocor, 
                      allow_new_levels, subset, nsamples)
  draws <- do.call(extract_draws, draws_args)
  if (length(object$nonlinear)) {
    draws$eta <- nonlinear_predictor(draws)
  } else {
    draws$eta <- linear_predictor(draws)
  }
  # see predict.R
  predict_fun <- get(paste0("predict_", draws$f$family), mode = "function")
  N <- if (!is.null(draws$data$N_trait)) draws$data$N_trait
       else if (!is.null(draws$data$N_tg)) draws$data$N_tg
       else if (is(object$autocor, "cor_fixed")) 1
       else draws$data$N
  out <- do.call(cbind, lapply(1:N, predict_fun, draws = draws, ntrys = ntrys))
  # percentage of invalid samples for truncated discrete models
  # should always be zero for all other models
  pct_invalid <- get_pct_invalid(out, data = draws$data)  # see predict.R
  if (pct_invalid >= 0.01) {
    warning(paste0(round(pct_invalid * 100), "% of all predicted values ", 
                   "were invalid. Increasing argument ntrys may help."))
  }
  # reorder predicted responses in case of multivariate models
  # as they are sorted after units first not after traits
  if (grepl("^multi_", draws$f$family)) {
    reorder <- with(draws$data, ulapply(1:K_trait, seq, to = N, by = K_trait))
    out <- out[, reorder, drop = FALSE]
    colnames(out) <- 1:ncol(out) 
  }
  # reorder predicted responses to be in the initial user defined order
  # currently only relevant for autocorrelation models 
  old_order <- attr(draws$data, "old_order")
  if (!is.null(old_order) && !sort) {
    out <- out[, old_order, drop = FALSE]  
    colnames(out) <- NULL
  }
  # transform predicted response samples before summarizing them 
  is_catordinal <- is.ordinal(object$family) || is.categorical(object$family)
  if (!is.null(transform) && !is_catordinal) {
    out <- do.call(transform, list(out))
  }
  if (summary) {
    if (is_catordinal) {
      # compute frequencies of categories 
      out <- get_table(out, levels = 1:max(draws$data$max_obs)) 
    } else {
      out <- get_summary(out, probs = probs, robust = robust)
    }
    rownames(out) <- 1:nrow(out)
  }
  out
}

#' Extract Model Fitted Values of \code{brmsfit} Objects
#' 
#' Predict fitted values (i.e. the 'regression line') of a fitted model.
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
#'
#' @return Fitted values extracted from \code{object}. 
#'  The output depends on the family:
#'  If \code{summary = TRUE} it is a N x E x C array 
#'  for categorical and ordinal models and a N x E matrix else.
#'  If \code{summary = FALSE} it is a S x N x C array 
#'  for categorical and ordinal models and a S x N matrix else.
#'  N is the number of observations, S is the number of samples, 
#'  C is the number of categories, and E is equal to \code{length(probs) + 2}.
#'   
#' @details \code{NA} values within factors in \code{newdata}, 
#'   are interpreted as if all dummy variables of this factor are 
#'   zero. This allows, for instance, to make predictions of the grand mean 
#'   when using sum coding.  
#' 
#'   For models fitted with \pkg{brms} <= 0.5.0 only: 
#'   Be careful when using \code{newdata} with factors 
#'   in fixed or random effects. The predicted results are only valid 
#'   if all factor levels present in the initial 
#'   data are also defined and ordered correctly 
#'   for the factors in \code{newdata}.
#'   Grouping factors may contain fewer levels than in the 
#'   inital data without causing problems.
#'   When using higher versions of \pkg{brms}, 
#'   all factors are automatically checked 
#'   for correctness and amended if necessary.
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
                           allow_new_levels = FALSE, incl_autocor = TRUE,
                           subset = NULL, nsamples = NULL, sort = FALSE,
                           summary = TRUE, robust = FALSE,
                           probs = c(0.025, 0.975), ...) {
  scale <- match.arg(scale)
  if (!is(object$fit, "stanfit") || !length(object$fit@sim)) 
    stop("The model does not contain posterior samples")
  draws_args <- nlist(x = object, newdata, re_formula, incl_autocor, 
                      allow_new_levels, subset, nsamples)
  draws <- do.call(extract_draws, draws_args)
  # get mu and scale it appropriately
  if (length(object$nonlinear)) {
    mu <- nonlinear_predictor(draws)
  } else {
    mu <- linear_predictor(draws)
  }
  if (scale == "response") {
    mu <- fitted_response(draws = draws, mu = mu)  # see fitted.R
  }
  # reorder fitted values to be in the initial user defined order
  # currently only relevant for autocorrelation models 
  old_order <- attr(draws$data, "old_order")
  if (!is.null(old_order) && !sort) {
    mu <- mu[, old_order, drop = FALSE]  
    colnames(mu) <- NULL
  }
  if (summary) {
    mu <- get_summary(mu, probs = probs, robust = robust)
    rownames(mu) <- 1:nrow(mu)
  }
  mu
}

#' Extract Model Residuals from brmsfit Objects
#' 
#' @inheritParams predict.brmsfit
#' @param type The type of the residuals, 
#'  either \code{"ordinary"} or \code{"pearson"}. 
#'  More information is provided under 'Details'. 
#' 
#' @details Residuals of type \code{ordinary} 
#'  are of the form \eqn{R = Y - Yp}, where \eqn{Y} is the observed 
#'  and \eqn{Yp} is the predicted response.
#'  Residuals of type \code{pearson} are 
#'  of the form \eqn{R = (Y - Yp) / SD(Y)},
#'  where \eqn{SD(Y)} is an estimation of the standard deviation 
#'  of \eqn{Y}. \cr
#'   
#'  Currently, \code{residuals.brmsfit} does not support 
#'  \code{categorical} or ordinal models. 
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
#'            data = inhaler, cluster = 2)
#' 
#' ## extract residuals 
#' res <- residuals(fit, summary = TRUE)
#' head(res)
#' }
#' 
#' @export
residuals.brmsfit <- function(object, newdata = NULL, re_formula = NULL, 
                              type = c("ordinary", "pearson"), 
                              allow_new_levels = FALSE, 
                              incl_autocor = TRUE, subset = NULL, 
                              nsamples = NULL, sort = FALSE,
                              summary = TRUE, robust = FALSE, 
                              probs = c(0.025, 0.975), ...) {
  type <- match.arg(type)
  family <- family(object)
  if (!is(object$fit, "stanfit") || !length(object$fit@sim)) 
    stop("The model does not contain posterior samples")
  if (is.ordinal(family) || is.categorical(family))
    stop(paste("residuals not yet implemented for family", family$family),
         call. = FALSE)
  
  standata <- amend_newdata(newdata, fit = object, re_formula = re_formula,
                            allow_new_levels = allow_new_levels, 
                            check_response = TRUE)
  if (is.null(subset) && !is.null(nsamples)) {
    subset <- sample(Nsamples(object), nsamples)
  }
  pred_args <- nlist(object, newdata, re_formula, allow_new_levels,
                     incl_autocor, subset, sort, summary = FALSE)
  mu <- do.call(fitted, pred_args)
  Y <- matrix(rep(as.numeric(standata$Y), nrow(mu)), 
              nrow = nrow(mu), byrow = TRUE)
  res <- Y - mu
  colnames(res) <- NULL
  if (type == "pearson") {
    # get predicted standard deviation for each observation
    pred_args$summary <- TRUE
    sd <- do.call(predict, pred_args)[, 2]
    sd <- matrix(rep(sd, nrow(mu)), nrow = nrow(mu), byrow = TRUE)
    res <- res / sd
  }
  if (summary) {
    res <- get_summary(res, probs = probs, robust = robust)
  }
  res
}

#' Update \pkg{brms} models
#' 
#' This method allows to update an existing \code{brmsfit} object
#' 
#' @param object object of class \code{brmsfit}
#' @param formula. changes to the formula; for details see 
#'   \code{\link[stats:update.formula]{update.formula}}
#' @param newdata optional \code{data.frame} 
#'  to update the model with new data
#' @param ... other arguments passed to 
#'  \code{\link[brms:brm]{brm}}
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
#' ## use another family and add fixed effects priors
#' fit4 <- update(fit1, family = weibull(), inits = "0",
#'                prior = set_prior("normal(0,5)"))
#' summary(fit4)
#' }
#'
#' @export
update.brmsfit <- function(object, formula., newdata = NULL, ...) {
  dots <- list(...)
  if ("data" %in% names(dots)) {
    # otherwise the data name cannot be found by substitute 
    stop("Please use argument 'newdata' to update the data", 
         call. = FALSE)
  }
  recompile <- FALSE
  if (missing(formula.)) {
    dots$formula <- object$formula
  } else {
    dots$formula <- as.formula(formula.)
    if (length(object$nonlinear)) {
      warning(paste("Argument 'formula.' will completely replace the", 
                    "original formula in non-linear models.", call. = FALSE))
      recompile <- TRUE
    } else {
      dots$formula <- update.formula(object$formula, dots$formula)
      ee_old <- extract_effects(object$formula, family = object$family)
      family <- get_arg("family", dots, object)
      ee_new <- extract_effects(dots$formula, family = family)
      # no need to recompile the model when changing fixed effects only
      recompile <- !(is_equal(names(ee_old), names(ee_new)) && 
        is_equal(ee_old$random, ee_new$random) &&
        is_equal(length(ee_old$response), length(ee_new$response)))
    }
    if (recompile) {
      message("The desired formula changes require recompling the model")
    }
  }
  # update gaussian("log") to lognormal() family
  resp <- extract_effects(object$formula, family = object$family)$response
  if (is.old_lognormal(object$family, nresp = length(resp),
                       version = object$version)) {
    object$family <- lognormal()
  }
  
  dots$iter <- first_not_null(dots$iter, object$fit@sim$iter)
  # brm computes warmup automatically based on iter 
  dots$chains <- first_not_null(dots$chains, object$fit@sim$chains)
  dots$thin <- first_not_null(dots$thin, object$fit@sim$thin)
  object$prior <- update_prior_frame(object$prior, ranef = object$ranef)
  rc_args <- c("family", "prior", "autocor", "nonlinear", "partial", 
               "threshold", "cov_ranef", "sparse", "sample_prior")
  new_args <- intersect(rc_args, names(dots))
  recompile <- recompile || length(new_args)
  if (recompile) {
    if (length(new_args)) {
      message(paste("Changing argument(s)", 
                    paste0("'", new_args, "'", collapse = ", "),
                    "requires recompiling the model"))
    }
    old_args <- setdiff(rc_args, new_args)
    dots[old_args] <- object[old_args]
    if (!is.null(newdata)) {
      dots$data <- newdata
      dots$data.name <- Reduce(paste, deparse(substitute(newdata)))
      dots$data.name <- substr(dots$data.name, 1, 50)
    } else  {
      dots$data <- object$data
      dots$data.name <- object$data.name
    }
    if (is.null(dots$threshold)) {
      # for backwards compatibility with brms <= 0.8.0
      if (grepl("(k - 1.0) * delta", object$model, fixed = TRUE)) {
        dots$threshold <- "equidistant"
      } else dots$threshold <- "flexible"
    }
    if ("prior" %in% new_args) {
      if (is(dots$prior, "brmsprior")) { 
        dots$prior <- c(dots$prior)
      } else if (!is(dots$prior, "prior_frame")) {
        stop("invalid prior argument")
      }
      dots$prior <- rbind(dots$prior, object$prior)
      dots$prior <- dots$prior[!duplicated(dots$prior[, 2:5]), ]
    }
    pnames <- parnames(object)
    if (is.null(dots$sample_prior)) {
      dots$sample_prior <- any(grepl("^prior_", pnames))
    }
    if (is.null(dots$ranef)) {
      dots$ranef <- any(grepl("^r_", pnames)) || !length(object$ranef)
    }
    if (is.null(dots$sparse)) {
      dots$sparse <- grepl("sparse matrix", stancode(object))
    }
    if (!isTRUE(dots$testmode)) {
      object <- do.call(brm, dots)
    }
  } else {
    # refit the model without compiling it again
    if (!is.null(dots$formula)) {
      object$formula <- dots$formula
      dots$formula <- NULL
    }
    ee <- extract_effects(object$formula, family = object$family, 
                          nonlinear = object$nonlinear)
    if (!is.null(newdata)) {
      object$data <- update_data(newdata, family = object$family, effects = ee)
      object$data.name <- Reduce(paste, deparse(substitute(newdata)))
      object$ranef <- gather_ranef(ee, data = object$data, 
                                   forked = is.forked(object$family))
      dots$is_newdata <- TRUE
    }
    if (!is.null(dots$ranef)) {
      object$exclude <- exclude_pars(ee, ranef = object$ranef, 
                                     save_ranef = dots$ranef)
    }
    if (!is.null(dots$algorithm)) {
      object$algorithm <- match.arg(dots$algorithm, 
                                    c("sampling", "meanfield", "fullrank"))
    }
    if (!isTRUE(dots$testmode)) {
      object <- do.call(brm, c(list(fit = object), dots))
    }
  }
  object
}

#' @export
#' @describeIn WAIC method for class \code{brmsfit}
WAIC.brmsfit <- function(x, ..., compare = TRUE, newdata = NULL, 
                         re_formula = NULL, allow_new_levels = FALSE, 
                         subset = NULL, nsamples = NULL, pointwise = NULL) {
  models <- list(x, ...)
  names <- c(deparse(substitute(x)), sapply(substitute(list(...))[-1], 
                                            deparse))
  if (is.null(subset) && !is.null(nsamples)) {
    subset <- sample(Nsamples(x), nsamples)
  }
  if (is.null(pointwise)) {
    pointwise <- set_pointwise(x, subset = subset, newdata = newdata)
  }
  ll_args = nlist(newdata, re_formula, allow_new_levels, subset, pointwise)
  if (length(models) > 1) {
    args <- nlist(X = models, FUN = compute_ic, ic = "waic", ll_args)
    out <- setNames(do.call(lapply, args), names)
    class(out) <- c("iclist", "list")
    if (compare && match_response(models)) {
      comp <- compare_ic(out, ic = "waic")
      attr(out, "compare") <- comp$ic_diffs
      attr(out, "weights") <- comp$weights
    }
  } else { 
    out <- do.call(compute_ic, nlist(x, ic = "waic", ll_args))
  }
  out
}

#' @export
#' @describeIn LOO method for class \code{brmsfit}
LOO.brmsfit <- function(x, ..., compare = TRUE, newdata = NULL, 
                        re_formula = NULL, allow_new_levels = FALSE, 
                        subset = NULL, nsamples = NULL, pointwise = NULL,
                        cores = 1, wcp = 0.2, wtrunc = 3/4) {
  models <- list(x, ...)
  names <- c(deparse(substitute(x)), sapply(substitute(list(...))[-1], 
                                            deparse))
  if (is.null(subset) && !is.null(nsamples)) {
    subset <- sample(Nsamples(x), nsamples)
  }
  if (is.null(pointwise)) {
    pointwise <- set_pointwise(x, subset = subset, newdata = newdata)
  }
  ll_args = nlist(newdata, re_formula, allow_new_levels, subset, pointwise)
  if (length(models) > 1) {
    args <- nlist(X = models, FUN = compute_ic, ic = "loo", 
                  ll_args, wcp, wtrunc, cores)
    out <- setNames(do.call(lapply, args), names)
    class(out) <- c("iclist", "list")
    if (compare && match_response(models)) {
      comp <- compare_ic(out, ic = "loo")
      attr(out, "compare") <- comp$ic_diffs
      attr(out, "weights") <- comp$weights
    }
  } else {
    out <- do.call(compute_ic, nlist(x, ic = "loo", ll_args, 
                                     wcp, wtrunc, cores))
  }
  out
}

#' Compute the pointwise log-likelihood
#' 
#' @param object A fitted model object of class \code{brmsfit}. 
#' @inheritParams predict.brmsfit
#' @param pointwise A flag indicating whether to compute the full
#'   log-likelihood matrix at once (the default), or just return
#'   the likelihood function along with all data and samples
#'   required to compute the log-likelihood separately for each
#'   observation. The latter option is rarely useful when
#'   calling \code{logLik} directly, but rather when computing
#'   \code{\link[brms:WAIC]{WAIC}} or \code{\link[brms:LOO]{LOO}}.
#' @param ... Currently ignored
#' 
#' @return Usually, an S x N matrix containing 
#'  the pointwise log-likelihood samples, 
#'  where S is the number of samples and N is the number 
#'  of observations in the data. 
#'  If \code{pointwise = TRUE}, the output is a function
#'  with a \code{draws} attribute containing all relevant
#'  data and posterior samples.
#' 
#' @importFrom statmod dinvgauss
#' @export
logLik.brmsfit <- function(object, newdata = NULL, re_formula = NULL,
                           allow_new_levels = FALSE, subset = NULL,
                           nsamples = NULL, pointwise = FALSE, ...) {
  if (!is(object$fit, "stanfit") || !length(object$fit@sim)) 
    stop("The model does not contain posterior samples")
  draws <- extract_draws(x = object, newdata = newdata, 
                         re_formula = re_formula, subset = subset,
                         allow_new_levels = allow_new_levels,
                         nsamples = nsamples, check_response = TRUE)
  N <- if (!is.null(draws$data$N_tg)) draws$data$N_tg
       else if (is(object$autocor, "cor_fixed")) 1
       else nrow(as.matrix(draws$data$Y))
  loglik_fun <- get(paste0("loglik_", draws$f$family), mode = "function")
  if (pointwise) {
    loglik <- structure(loglik_fun, draws = draws, N = N)
  } else {
    if (length(object$nonlinear)) {
      draws$eta <- nonlinear_predictor(draws)
    } else {
      draws$eta <- linear_predictor(draws)
    }
    loglik <- do.call(cbind, lapply(1:N, loglik_fun, draws = draws))
    # reorder loglik values to be in the initial user defined order
    # currently only relevant for autocorrelation models
    # that are not using covariance formulation
    old_order <- attr(draws$data, "old_order")
    if (!is.null(old_order) && !isTRUE(object$autocor$cov)) {
      loglik <- loglik[, old_order[1:N]]  
    }
    colnames(loglik) <- NULL
  }
  loglik
}

#' @rdname hypothesis
#' @export
hypothesis.brmsfit <- function(x, hypothesis, class = "b", group = "",
                               alpha = 0.05, ...) {
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  if (!is.character(hypothesis)) 
    stop("Argument hypothesis must be a character vector", call. = FALSE)
  if (alpha < 0 || alpha > 1)
    stop("Argument alpha must be in [0,1]", call. = FALSE)
  
  # process class and group arguments
  if (is.null(class)) class <- ""
  valid_classes <- c("", "b", "r", "sd", "cor", "ar", "ma", "arr", 
                     "sigma", "rescor", "nu", "shape", "delta")
  if (!class %in% valid_classes)
    stop(paste(class, "is not a valid paramter class"), call. = FALSE)
  if (class %in% c("b", "r", "sd", "cor", "sigma", "rescor")) {
    if (class %in% c("sd", "cor") && nchar(group)) {
      class <- paste0(class, "_", group, "_")
    } else {
      class <- paste0(class, "_")
    }
  } else {
    # no class required
    class <- ""  
  }

  hyp_fun <- function(h) {
    # internal function to evaluate hypotheses
    # 
    # Args:
    #   h: A string containing a hypothesis
    h <- rename(h, c("[ \t\r\n]", ":"), c("", "__"), fixed = FALSE)
    sign <- unlist(regmatches(h, gregexpr("=|<|>", h)))
    lr <- unlist(regmatches(h, gregexpr("[^=<>]+", h)))
    if (length(sign) != 1 || length(lr) != 2)
      stop("Every hypothesis must be of the form 'left (= OR < OR >) right'",
           call. = FALSE)
    h <- paste0(lr[1], ifelse(lr[2] != "0", paste0("-(",lr[2],")"), ""))
    varsH <- unique(find_names(h))
    parsH <- paste0(class, varsH)
    if (!all(parsH %in% pars)) 
      stop(paste("The following parameters cannot be found in the model:", 
                 paste0(gsub("__", ":", parsH[which(!parsH %in% pars)]), 
                        collapse = ", ")), call. = FALSE)
    
    # prepare for renaming of parameters so that h can be evaluated
    parsH <- rename(parsH, "__", ":")
    h_renamed <- rename(h, c("[", "]", ","), c(".", ".", ".."))
    symbols <- c(paste0("^",class), ":", "\\[", "\\]", ",")
    subs <- c("", "__", ".", ".", "..")
    
    # get posterior samples
    samples <- posterior_samples(x, pars = parsH, exact_match = TRUE)
    names(samples) <- rename(names(samples), symbols = symbols, 
                             subs = subs, fixed = FALSE)
    samples <- matrix(with(samples, eval(parse(text = h_renamed))), ncol = 1)
    
    # get prior samples
    prior_samples <- prior_samples(x, pars = parsH, fixed = TRUE)
    if (!is.null(prior_samples) && ncol(prior_samples) == length(varsH)) {
      names(prior_samples) <- rename(names(prior_samples), symbols = symbols, 
                                     subs = subs, fixed = FALSE)
      prior_samples <- matrix(with(prior_samples, eval(parse(text = h_renamed))), 
                              ncol = 1)
    } else prior_samples <- NULL
    
    # evaluate hypothesis
    wsign <- ifelse(sign == "=", "equal", ifelse(sign == "<", "less", "greater"))
    probs <- switch(wsign, equal = c(alpha / 2, 1 - alpha / 2), 
                    less = c(0, 1 - alpha), greater = c(alpha, 1))
    sm <- lapply(c("mean", "sd", "quantile", "evidence_ratio"), 
                 get_estimate, samples = samples, probs = probs, 
                 wsign = wsign, prior_samples = prior_samples, ...)
    sm <- as.data.frame(matrix(unlist(sm), nrow = 1))
    if (sign == "<") {
      sm[1, 3] <- -Inf
    } else if (sign == ">") {
      sm[1, 4] <- Inf
    }
    sm <- cbind(sm, ifelse(!(sm[1, 3] <= 0 && 0 <= sm[1, 4]), '*', ''))
    h <- paste0(lr[1], ifelse(lr[2] != "0", paste0("-(", lr[2], ")"), ""))
    rownames(sm) <- paste(rename(h, "__", ":"), sign, "0")
    cl <- (1 - alpha) * 100
    colnames(sm) <- c("Estimate", "Est.Error", paste0("l-",cl,"% CI"), 
                      paste0("u-",cl,"% CI"), "Evid.Ratio", "")
    if (!is.null(prior_samples)) {
      samples <- c(samples, prior_samples)
    } else {
      samples <- c(samples, rep(NA, nrow(samples)))
    }
    nlist(summary = sm, samples = samples)
  }

  pars <- rename(parnames(x)[grepl("^", class, parnames(x))],
                 symbols = ":", subs = "__")
  hlist <- lapply(hypothesis, hyp_fun)
  hs <- do.call(rbind, lapply(hlist, function(h) h$summary))
  samples <- do.call(cbind, lapply(hlist, function(h) h$samples))
  samples <- as.data.frame(samples) 
  names(samples) <- paste0("H", seq_along(hlist))
  samples$Type <- factor(rep(c("posterior", "prior"), each = Nsamples(x)))
  class <- substr(class, 1, nchar(class) - 1)
  out <- nlist(hypothesis = hs, samples, class, alpha)
  class(out) <- "brmshypothesis"
  out
}

#' @rdname expose_functions
#' @export
expose_functions.brmsfit <- function(x, ...) {
  expose_stan_functions(x$fit)
}
