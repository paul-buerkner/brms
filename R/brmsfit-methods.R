#' @export
parnames.brmsfit <- function(x, ...) {
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  dimnames(x$fit)$parameters
}

#' Extract Fixed Effects Estimates
#' 
#' Extract the fixed effects from a \code{brmsfit} object. 
#' 
#' @aliases fixef
#' 
#' @param object An object of class \code{brmsfit}
#' @param estimate A character vector specifying which coefficients 
#'  (e.g., "mean", "median", "sd", or "quantile") 
#'  should be calculated for the fixed effects.
#' @param ... Further arguments to be passed to the functions 
#'  specified in \code{estimate}
#' 
#' @return A matrix with one row per fixed effect 
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

#' Covariance and Correlation Matrix of Fixed Effects
#' 
#' Get a point estimate of the covariance or 
#' correlation matrix of fixed effects parameters
#' 
#' @param object An object of class \code{brmsfit}
#' @param correlation logical; if \code{FALSE} (the default), 
#'   compute the covariance matrix,
#'   if \code{TRUE}, compute the correlation matrix
#' @param ... Currently ignored
#' 
#' @return covariance or correlation matrix of fixed effects parameters
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

#' Extract Random Effects Estimates
#' 
#' Extract the random effects of each level from \code{brmsfit} object. 
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
#'  random effects as column names 
#'     
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}   
#'   
#' @examples
#' \dontrun{
#' fit <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1+Trt_c|visit), 
#'              data = epilepsy, family = "poisson", chains = 1)
#' ## random effects means with corresponding covariances
#' rf <- ranef(fit, var = TRUE)
#' attr(rf, "var")
#' ## random effects medians
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
#' Extract model coefficients, which are the sum of fixed effects
#' and corresponding random effects
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
#' @param sigma,rdig Ignored (included for compatibility with 
#'  \code{\link[nlme:VarCorr]{VarCorr}}).
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
#' @importFrom nlme VarCorr
#' @export VarCorr
#' @export
VarCorr.brmsfit <- function(x, sigma = 1, rdig = 3, estimate = "mean", 
                            as.list = TRUE, ...) {
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  if (!(length(x$ranef) || any(grepl("^sigma_", parnames(x)))))
    stop("The model does not contain covariance matrices", call. = FALSE)
  if (!isTRUE(all.equal(sigma, 1)))
    warning("argument 'sigma' is unused")
  if (!isTRUE(all.equal(rdig, 3)))
    warning("argument 'rdig' is unused")
  
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
  samples_taken <- seq((warmup + 1), iter, thin)
  
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
        is_partial <- grepl("^b_", par) && grepl("\\[[[:digit:]]+\\]", par)
        # ensures correct parameter to prior mapping for partial effects
        par_internal <- ifelse(is_partial, sub("^b_", "bp_", par), par)
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
  formula <- update_formula(object$formula, partial = object$partial)
  out <- brmssummary(formula = formula,
                     family = family, 
                     data.name = object$data.name, 
                     group = names(object$ranef), 
                     nobs = nobs(object), 
                     ngrps = ngrps(object), 
                     autocor = object$autocor,
                     algorithm = algorithm(object))
  
  if (length(object$fit@sim)) {
    out$chains <- object$fit@sim$chains
    out$iter <- object$fit@sim$iter
    out$warmup <- object$fit@sim$warmup
    out$thin <- object$fit@sim$thin
    stan_args <- object$fit@stan_args[[1]]
    out$sampler <- paste0(stan_args$method, "(", stan_args$algorithm, ")")
    if (length(object$ranef) && !any(grepl("^r_", parnames(object)))
        || length(ee$response) > 1 && is.linear(family)) {
      # if brm(..., ranef = FALSE) or model is multivariate
      waic <- FALSE
    }
    if (waic) out$WAIC <- WAIC(object)$waic
    
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
    } else {
      colnames(fit_summary) <- c("Estimate", "Est.Error", "l-95% CI", 
                                 "u-95% CI")
    }
    
    # fixed effects summary
    fix_pars <- pars[grepl("^b_", pars)]
    out$fixed <- fit_summary[fix_pars, , drop = FALSE]
    rownames(out$fixed) <- gsub("^b_", "", fix_pars)
    
    # summary of family specific parameters
    spec_pars <- pars[pars %in% c("nu","shape","delta", "phi") | 
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
    cor_pars <- pars[grepl("^ar|^ma", pars)]
    out$cor_pars <- fit_summary[cor_pars, , drop = FALSE]
    rownames(out$cor_pars) <- cor_pars
    
    for (i in seq_along(out$group)) {
      nlp <- get_nlpar(object$ranef[[i]])
      nlp_ <- ifelse(nchar(nlp), paste0(nlp, "_"), nlp)
      rnames <- object$ranef[[i]]
      sd_pars <- paste0("sd_", nlp_, out$group[i], "_", rnames)
      sd_type <- ifelse(nchar(nlp), paste0("sd_", nlp), "sd")
      sd_names <- paste0(sd_type, "(", rnames,")")
      all_cor_pars <- get_cornames(rnames, brackets = FALSE,
                                   type = paste0("cor_", nlp_, out$group[i]))
      cor_pars <- intersect(all_cor_pars, parnames(object))
      cor_type <- ifelse(nchar(nlp), paste0("cor_", nlp), "cor") 
      cor_names <- get_cornames(rnames, type = cor_type,
                                subset = cor_pars, 
                                subtype = out$group[i])
      out$random[[out$group[i]]] <- 
        fit_summary[c(sd_pars, cor_pars), , drop = FALSE]
      rownames(out$random[[out$group[i]]]) <- c(sd_names, cor_names)
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
  standata <- standata(object)
  ee <- extract_effects(object$formula, family = object$family,
                        nonlinear = object$nonlinear)
  rand <- get_random(ee)
  group <- rand$group
  if (length(group)) {
    nlp <- if (length(ee$nonlinear)) paste0(rownames(rand), "_") 
           else rep("", nrow(rand))
    .fun <- function(i) {
      n <- standata[[paste0("N_", nlp[i], i)]]
      if (is.null(n)) {
        n <- standata[[paste0("N_", group[[i]])]]
      }
      return(n)
    }
    out <- setNames(lapply(seq_along(group), .fun), group)
    out <- out[!duplicated(group)]
  } else out <- NULL
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
stancode.brmsfit <- function(object, ...)
  object$model

#' @export
standata.brmsfit <- function(object, ...) {
  dots <- list(...)
  if (is.data.frame(object$data)) {
    # brms > 0.5.0 stores the original model.frame
    new_formula <- update_re_terms(object$formula, dots$re_formula)
    standata <- make_standata(new_formula,  
                              data = object$data, 
                              family = object$family, 
                              autocor = object$autocor, 
                              nonlinear = object$nonlinear,
                              cov_ranef = object$cov_ranef, 
                              partial = object$partial, ...)
  } else {
    # brms <= 0.5.0 only stores the data passed to Stan 
    standata <- object$data
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
#' @param theme The ggplot theme to use. For details see
#'  \code{\link[ggplot2:ggtheme]{ggtheme}}.
#' @param do_plot logical; indicates if plots should be
#'   plotted directly in the active graphic device.
#'   Defaults to \code{TRUE}.
#' @param ask logical; indicates if the user is prompted 
#'   before a new page is plotted. 
#'   Only used if \code{do_plot} is \code{TRUE}.
#' @param newpage logical; indicates if the first set of plots
#'   should be plotted to a new page. 
#'   Only used if \code{do_plot} is \code{TRUE}.
#' @param ... Further arguments passed to 
#'   \code{\link[gridExtra:arrangeGrob]{arrangeGrob}}.
#' 
#' @return A (possibly invisible) list of 
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
#' @export
plot.brmsfit <- function(x, pars = NA, parameters = NA, N = 5, 
                         theme = "classic", ask = TRUE, 
                         do_plot = TRUE, newpage = TRUE, ...) {
  dots <- list(...)
  if (is.na(pars[1])) 
    pars <- parameters 
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  if (!is.wholenumber(N) || N < 1) 
    stop("N must be a positive integer", call. = FALSE)
  if (!is.character(pars)) {
    pars <- c("^b_", "^bm_", "^sd_", "^cor_", "^sigma", "^rescor", 
              "^nu$", "^shape$", "^delta$", "^phi$", "^ar", "^ma", "^arr")
  }
  samples <- posterior_samples(x, pars = pars, add_chain = TRUE)
  pars <- names(samples)[!names(samples) %in% c("chain", "iter")] 
  if (length(pars) == 0) {
    stop("No valid parameters selected", call. = FALSE)
  }
  
  if (do_plot) {
    default_ask <- devAskNewPage()
    on.exit(devAskNewPage(default_ask))
    devAskNewPage(ask = FALSE)
  }
  n_plots <- ceiling(length(pars) / N)
  plots <- vector(mode = "list", length = n_plots)
  for (i in 1:n_plots) {
    temp_plot <- lapply(pars[((i - 1) * N + 1):min(i * N, length(pars))], 
                        trace_density_plot, x = samples, theme = theme)
    plots[[i]] <- arrangeGrob(grobs = unlist(temp_plot, recursive = FALSE), 
                              nrow = length(temp_plot), ncol = 2, ...)
    if (do_plot) {
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
      if (length(pars) > 1) {
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
marginal_effects.brmsfit <- function(x, effects = NULL, data = NULL, 
                                     re_formula = NA, probs = c(0.025, 0.975),
                                     method = c("fitted", "predict"), ...) {
  method <- match.arg(method)
  ee <- extract_effects(x$formula, family = x$family)
  if (is.linear(x$family) && length(ee$response) > 1) {
    stop("Marginal plots are not yet implemented for multivariate models.",
         call. = FALSE)
  }
  rsv_vars <- rsv_vars(x$family, nresp = length(ee$response))
  all_effects <- strsplit(attr(terms(ee$fixed), "term.labels"), split = ":")
  all_effects <- rmNULL(lapply(all_effects, setdiff, y = rsv_vars))
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
    matches <- match(lapply(all_effects, sort), 
                     lapply(effects, sort), 0L)
    effects <- unique(effects[sort(matches)])
  }
  if (!length(effects)) {
    stop("No valid effects specified.", call. = FALSE)
  }
  if (any(ulapply(effects, length) > 2)) {
    stop("Interactions of order higher than 2 are currently not supported.",
         call. = FALSE)
  }
  if (length(probs) != 2L) {
    stop("Arguments 'probs' must be of length 2.", call. = FALSE)
  }
  if (is.ordinal(x$family) || is.categorical(x$family)) {
    warning(paste0("Predictions are treated as continuous variables ", 
                   "in marginal plots, \nwhich is likely an invalid ", 
                   "assumption for family ", x$family$family, "."),
            call. = FALSE)
  }
  
  # prepare marginal data
  if (is.null(data)) {
    if (!isTRUE(all.equal(x$autocor, cor_arma())) || 
        length(rmNULL(ee[c("se", "trials", "cat")]))) {
      stop("Please specify argument 'data' manually for this model.", 
           call. = FALSE)
    }
    vars <- unique(ulapply(c(ee$fixed[[3]], ee$random$form), all.vars))
    vars <- setdiff(vars, rsv_vars)
    data <- as.data.frame(as.list(rep(NA, length(vars))))
    names(data) <- vars
    for (v in vars) {
      if (is.numeric(x$data[[v]])) {
        data[[v]] <- mean(x$data[[v]])
      } else {
        # use reference category
        data[[v]] <- attr(as.factor(x$data[[v]]), "levels")[1]
      }
    }
  } else if (is.data.frame(data)) {
    used_effects <- unique(unlist(effects))
    is_everywhere <- ulapply(used_effects, function(up)
      all(ulapply(effects, function(pred) up %in% pred)))
    non_marg_effects <- used_effects[is_everywhere]
    # effects that are present in every effect term
    # do not need to be defined in data
    missing_effects <- setdiff(non_marg_effects, names(data)) 
    data[, missing_effects] <- x$data[1, missing_effects] 
  } else {
    stop("data must be a data.frame or NULL")
  }
  data <- amend_newdata(data, fit = x, re_formula = re_formula,
                        allow_new_levels = TRUE, return_standata = FALSE)

  results <- list()
  for (i in seq_along(effects)) {
    marg_data <- x$data[, effects[[i]], drop = FALSE]
    pred_types <- ifelse(ulapply(marg_data, is.numeric), "numeric", "factor")
    if (length(effects[[i]]) == 2L) {
      # numeric effects should come first
      new_order <- order(pred_types, decreasing = TRUE)
      effects[[i]] <- effects[[i]][new_order]
      pred_types <- pred_types[new_order]
      if (pred_types[1] == "numeric") {
        values <- setNames(vector("list", length = 2), effects[[i]])
        values[[1]] <- unique(marg_data[, effects[[i]][1]])
        if (pred_types[2] == "numeric") {
          mean2 <- mean(marg_data[, effects[[i]][2]])
          sd2 <- sd(marg_data[, effects[[i]][2]])
          values[[2]] <- (-1:1) * sd2 + mean2
        } else {
          values[[2]] <- unique(marg_data[, effects[[i]][2]])
        }
        marg_data <- do.call(expand.grid, values)
      }
    } else {
      marg_data <- unique(marg_data)
    }
    marg_data <- replicate(nrow(data), simplify = FALSE,
     expr = marg_data[do.call(order, as.list(marg_data)), , drop = FALSE])
    marg_vars <- setdiff(names(data), effects[[i]])
    for (j in 1:nrow(data)) {
      marg_data[[j]][, marg_vars] <- data[j, marg_vars]
      marg_data[[j]][["MargRow"]] <- j
    }
    marg_data <- do.call(rbind, marg_data)
    args <- list(x, newdata = marg_data, re_formula = re_formula,
                 allow_new_levels = TRUE, probs = probs)
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
      marg_res <- get_summary(marg_res, probs = probs)
    } else {
      marg_res <- do.call(method, args)
    }
    colnames(marg_res)[3:4] <- c("lowerCI", "upperCI")
     
    if (length(effects[[i]]) == 2L && all(pred_types == "numeric")) {
      # can only be converted to factor after having called method
      labels <- c("Mean - SD", "Mean", "Mean + SD")
      marg_data[[effects[[i]][2]]] <- 
        factor(marg_data[[effects[[i]][2]]], labels = labels)
    }
    marg_res = cbind(marg_data, marg_res)
    attr(marg_res, "response") <- as.character(x$formula[2])
    attr(marg_res, "effects") <- effects[[i]]
    attr(marg_res, "rug") <- marg_data[, effects[[i]], drop = FALSE]
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
#' @param summary Should summary statistics 
#'   (i.e. means, sds, and 95\% intervals) be returned
#'  instead of the raw values? Default is \code{TRUE}.
#' @param probs The percentiles to be computed 
#'  by the \code{quantile} function. 
#'  Only used if \code{summary} is \code{TRUE}.
#' @param subset A numeric vector specifying
#'  the posterior samples to be used. 
#'  If \code{NULL} (the default), all samples are used.
#' @param nsamples Positive integer indicating how many 
#'  posterior samples should be used. 
#'  If \code{NULL} (the default) all samples are used.
#'  Ignored if \code{subset} is not \code{NULL}.
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
#' @details For truncated discrete models only:
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
                            subset = NULL, nsamples = NULL, ntrys = 5, 
                            summary = TRUE, probs = c(0.025, 0.975), ...) {
  if (!is(object$fit, "stanfit") || !length(object$fit@sim)) 
    stop("The model does not contain posterior samples")
  family <- family(object)
  ee <- extract_effects(object$formula, family = family)
  standata <- amend_newdata(newdata, fit = object, re_formula = re_formula,
                            allow_new_levels = allow_new_levels)

  # compute all necessary samples
  if (is.null(subset) && !is.null(nsamples)) {
    subset <- sample(Nsamples(object), nsamples)
  }
  nresp <- length(ee$response)
  samples <- list(eta = linear_predictor(object, newdata = standata, 
                                         re_formula = re_formula,
                                         subset = subset))
  args <- list(x = object, as.matrix = TRUE, subset = subset) 
  if (has_sigma(family, se = ee$se, autocor = object$autocor))
    samples$sigma <- do.call(posterior_samples, c(args, pars = "^sigma_"))
  if (family$family == "student") 
    samples$nu <- do.call(posterior_samples, c(args, pars = "^nu$"))
  if (family$family %in% c("beta", "zero_inflated_beta"))
    samples$phi <- do.call(posterior_samples, c(args, pars = "^phi$"))
  if (has_shape(family)) 
    samples$shape <- do.call(posterior_samples, c(args, pars = "^shape$"))
  if (is.linear(family) && nresp > 1) {
    samples$rescor <- do.call(posterior_samples, c(args, pars = "^rescor_"))
    samples$Sigma <- get_cov_matrix(sd = samples$sigma, cor = samples$rescor)$cov
    message(paste("Computing predicted values of a multivariate model. \n", 
                  "This may take a while."))
  }
  
  # call predict functions
  autocor <- object$autocor
  if (is.lognormal(family, nresp = nresp)) {
    family$family <- "lognormal"
    family$link <- "identity"
  } else if (is.linear(family) && nresp > 1) {
    family$family <- paste0("multi_", family$family)
  } else if (use_cov(autocor) && (get_ar(autocor) || get_ma(autocor))) {
    # special family for ARMA models using residual covariance matrices
    family$family <- paste0(family$family, "_cov")
    samples$ar <- do.call(posterior_samples, c(args, pars = "^ar\\["))
    samples$ma <- do.call(posterior_samples, c(args, pars = "^ma\\["))
  } 
  
  is_catordinal <- is.ordinal(family) || is.categorical(family)
  # see predict.R
  predict_fun <- get(paste0("predict_", family$family), mode = "function")
  call_predict_fun <- function(n) {
    do.call(predict_fun, list(n = n, data = standata, samples = samples, 
                              link = family$link, ntrys = ntrys))
  }
  N <- if (!is.null(standata$N_trait)) standata$N_trait
       else if (!is.null(standata$N_tg)) standata$N_tg
       else standata$N
  out <- do.call(cbind, lapply(1:N, call_predict_fun))
  rm(samples)
  
  # percentage of invalid samples for truncated discrete models
  # should always be zero for all other models
  pct_invalid <- get_pct_invalid(out, data = standata)  # see predict.R
  if (pct_invalid >= 0.01) {
    warning(paste0(round(pct_invalid * 100), "% of all predicted values ", 
                   "were invalid. Increasing argument ntrys may help."))
  }
  
  # reorder predicted responses in case of multivariate models
  # as they are sorted after units first not after traits
  if (grepl("^multi_", family$family)) {
    reorder <- with(standata, ulapply(1:K_trait, seq, to = N, by = K_trait))
    # observations in columns
    out <- out[, reorder, drop = FALSE]  
    colnames(out) <- 1:ncol(out) 
  }
  # reorder predicted responses to be in the initial user defined order
  # currently only relevant for autocorrelation models 
  old_order <- attr(standata, "old_order")
  if (!isTRUE(all.equal(old_order, 1:ncol(out)))) {
    out <- out[, old_order, drop = FALSE]  
    colnames(out) <- 1:ncol(out) 
  }
  # transform predicted response samples before summarizing them 
  if (!is.null(transform) && !is_catordinal) {
    out <- do.call(transform, list(out))
  }
  if (summary && !is_catordinal) {
    out <- get_summary(out, probs = probs)
  } else if (summary && is_catordinal) { 
    # compute frequencies of categories for categorical and ordinal models
    out <- get_table(out, levels = 1:max(standata$max_obs)) 
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
#' @details For models fitted with \pkg{brms} <= 0.5.0 only: 
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
                           allow_new_levels = FALSE,
                           subset = NULL, nsamples = NULL, 
                           summary = TRUE, probs = c(0.025, 0.975), ...) {
  scale <- match.arg(scale)
  if (!is(object$fit, "stanfit") || !length(object$fit@sim)) 
    stop("The model does not contain posterior samples")
  family <- family(object)
  ee <- extract_effects(object$formula, family = family)
  standata <- amend_newdata(newdata, fit = object, re_formula = re_formula,
                            allow_new_levels = allow_new_levels)
  
  if (is.null(subset) && !is.null(nsamples)) {
    subset <- sample(Nsamples(object), nsamples)
  }
  # get mu and scale it appropriately
  mu <- linear_predictor(object, newdata = standata, subset = subset, 
                         re_formula = re_formula)
  if (scale == "response") {
    # see fitted.R
    mu <- fitted_response(object, eta = mu, data = standata)
  }
  # reorder fitted values to be in the initial user defined order
  # currently only relevant for autocorrelation models 
  old_order <- attr(standata, "old_order")
  if (!isTRUE(all.equal(old_order, 1:ncol(mu)))) {
    mu <- mu[, old_order, drop = FALSE]  
    colnames(mu) <- 1:ncol(mu) 
  }
  if (summary) {
    mu <- get_summary(mu, probs = probs)
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
                              subset = NULL, nsamples = NULL,
                              summary = TRUE, probs = c(0.025, 0.975), ...) {
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
  mu <- fitted(object, newdata = standata, re_formula = re_formula, 
               allow_new_levels = allow_new_levels, 
               summary = FALSE, subset = subset)
  Y <- matrix(rep(as.numeric(standata$Y), nrow(mu)), 
              nrow = nrow(mu), byrow = TRUE)
  res <- Y - mu
  colnames(res) <- NULL
  if (type == "pearson") {
    # get predicted standard deviation for each observation
    sd <- predict(object, newdata = standata, re_formula = re_formula, 
                  allow_new_levels = allow_new_levels, 
                  summary = TRUE, subset = subset)[, 2]
    sd <- matrix(rep(sd, nrow(mu)), nrow = nrow(mu), byrow = TRUE)
    res <- res / sd
  }
  if (summary) {
    res <- get_summary(res, probs = probs)
  }
  res
}

#' Update \pkg{brms} models
#' 
#' This method allows to update an existing \code{brmsfit} object 
#' with new data as well as changed configuration of the chains.
#' 
#' @param object object of class \code{brmsfit}
#' @param newdata optional \code{data.frame} 
#'  to update the model with new data
#' @param ... other arguments passed to 
#'  \code{\link[brms:brm]{brm}} such as
#'  \code{iter} or \code{chains}.
#'
#' @export
update.brmsfit <- function(object, newdata = NULL, ...) {
  dots <- list(...)
  invalid_args <- c("formula", "family", "prior", "autocor", 
                    "partial", "threshold", "cov_ranef", 
                    "sample_prior", "save_model")
  z <- which(names(dots) %in% invalid_args)
  if (length(z)) {
    stop(paste("Argument(s)", paste(names(dots)[z], collapse = ", "),
               "cannot be updated"), call. = FALSE)
  }
  if ("data" %in% names(dots)) {
    stop("Please use argument 'newdata' to update your data", call. = FALSE)
  }
  # update arguments if required
  ee <- extract_effects(object$formula)
  if (!is.null(newdata)) {
    object$data <- amend_newdata(newdata, fit = object, 
                                 return_standata = FALSE)
    object$data.name <- Reduce(paste, deparse(substitute(newdata)))
    object$ranef <- gather_ranef(ee, data = object$data, 
                                 is_forked = is.forked(object$family))
    dots$is_newdata <- TRUE
  }
  if (!is.null(dots$ranef)) {
    object$exclude <- exclude_pars(object$formula, ranef = dots$ranef)
  }
  if (!isFALSE(dots$refit)) {
    # allows test 'update' without having to fit a Stan model
    dots$refit <- NULL
    object <- do.call(brm, c(list(fit = object), dots))
  }
  object
}

#' @export
WAIC.brmsfit <- function(x, ..., compare = TRUE) {
  models <- list(x, ...)
  names <- c(deparse(substitute(x)), sapply(substitute(list(...))[-1], 
                                            deparse))
  if (length(models) > 1) {
    out <- setNames(lapply(models, compute_ic, ic = "waic"), names)
    class(out) <- c("iclist", "list")
    if (compare && match_response(models)) {
      comp <- compare_ic(out, ic = "waic")
      attr(out, "compare") <- comp$ic_diffs
      attr(out, "weights") <- comp$weights
    }
  } else { 
    out <- compute_ic(x, ic = "waic")
  }
  out
}

#' @export
#' @describeIn LOO method for class \code{brmsfit}
LOO.brmsfit <- function(x, ..., compare = TRUE,
                        cores = getOption("loo.cores", parallel::detectCores()),
                        wcp = 0.2, wtrunc = 3/4) {
  models <- list(x, ...)
  names <- c(deparse(substitute(x)), sapply(substitute(list(...))[-1], 
                                            deparse))
  if (length(models) > 1) {
    out <- setNames(lapply(models, compute_ic, ic = "loo", wcp = wcp, 
                           wtrunc = wtrunc, cores = cores), names)
    class(out) <- c("iclist", "list")
    if (compare && match_response(models)) {
      comp <- compare_ic(out, ic = "loo")
      attr(out, "compare") <- comp$ic_diffs
      attr(out, "weights") <- comp$weights
    }
  } else {
    out <- compute_ic(x, ic = "loo", wcp = wcp, wtrunc = wtrunc, 
                      cores = cores)
  }
  out
}

#' Compute the pointwise log-likelihood
#' 
#' @param object A fitted model object of class \code{brmsfit}. 
#' @param ... Currently ignored
#' 
#' @return Usually, an S x N matrix containing 
#'  the pointwise log-likelihood samples, 
#'  where S is the number of samples and N is the number 
#'  of observations in the data. 
#' 
#' @importFrom statmod dinvgauss
#' @export
logLik.brmsfit <- function(object, ...) {
  if (!is(object$fit, "stanfit") || !length(object$fit@sim)) 
    stop("The model does not contain posterior samples")
  family <- family(object)
  ee <- extract_effects(object$formula, family = family)
  nresp <- length(ee$response)
  data <- standata(object, control = list(save_order = TRUE))
  N <- ifelse(is.null(data$N_tg), nrow(as.matrix(data$Y)), data$N_tg)
  
  # extract relevant samples
  samples <- list(eta = linear_predictor(object))
  if (has_sigma(family, se = ee$se, autocor = object$autocor))
    samples$sigma <- as.matrix(posterior_samples(object, pars = "^sigma_"))
  if (family$family == "student") 
    samples$nu <- as.matrix(posterior_samples(object, pars = "^nu$"))
  if (family$family %in% c("beta", "zero_inflated_beta"))
    samples$phi <- as.matrix(posterior_samples(object, pars = "^phi$"))
  if (has_shape(family)) 
    samples$shape <- as.matrix(posterior_samples(object, pars = "^shape$"))
  if (is.linear(family) && nresp > 1) {
    samples$rescor <- as.matrix(posterior_samples(object, pars = "^rescor_"))
    samples$Sigma <- get_cov_matrix(sd = samples$sigma, cor = samples$rescor)$cov
    message(paste("Computing pointwise log-likelihood of a", 
                  "multivariate model.\n This may take a while."))
  }
  
  # prepare for calling family specific loglik functions
  autocor <- object$autocor
  if (is.lognormal(family, nresp = nresp)) {
    family$family <- "lognormal"
    family$link <- "identity"
  } else if (is.linear(family) && nresp > 1) {
    family$family <- paste0("multi_", family$family)
  } else if (use_cov(autocor) && (get_ar(autocor) || get_ma(autocor))) {
    # special family for ARMA models using residual covariance matrices
    family$family <- paste0(family$family, "_cov")
    samples$ar <- posterior_samples(object, pars = "^ar\\[", as.matrix = TRUE)
    samples$ma <- posterior_samples(object, pars = "^ma\\[", as.matrix = TRUE)
  } 

  loglik_fun <- get(paste0("loglik_", family$family), mode = "function")
  call_loglik_fun <- function(n) {
    do.call(loglik_fun, list(n = n, data = data, samples = samples, 
                             link = family$link)) 
  }
  loglik <- do.call(cbind, lapply(1:N, call_loglik_fun))
  # reorder loglik values to be in the initial user defined order
  # currently only relevant for autocorrelation models
  # that are not using covariance formulation
  old_order <- attr(data, "old_order")
  if (!isTRUE(all.equal(old_order[1:N], 1:N)) && !isTRUE(autocor$cov)) {
    loglik <- loglik[, old_order[1:N]]  
  }
  colnames(loglik) <- 1:ncol(loglik)
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
    h <- rename(h, symbols = c(" ", ":"), subs = c("", "__"))
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