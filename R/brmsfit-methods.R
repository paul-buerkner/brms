#' @export
parnames.brmsfit <- function(x, ...) {
  contains_samples(x)
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
#'  (e.g., \code{"mean"}, \code{"median"}, \code{"sd"}, or \code{"quantile"})
#'  should be calculated for the population-level effects.
#' @param ... Further arguments to be passed to the functions 
#'  specified in \code{estimate}
#' 
#' @return A matrix with one row per population-level effect 
#'   and one column per calculated estimate.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@gmail.com}
#' 
#' @examples
#' \dontrun{
#' fit <- brm(time | cens(censored) ~ age + sex + disease, 
#'            data = kidney, family = "exponential")
#' fixef(fit, estimate = c("mean", "sd"))
#' }
#' 
#' @method fixef brmsfit
#' @export
#' @export fixef
#' @importFrom nlme fixef
fixef.brmsfit <-  function(object, estimate = "mean", ...) {
  contains_samples(object)
  pars <- parnames(object)
  fpars <- pars[grepl(fixef_pars(), pars)]
  if (!length(fpars)) {
    stop2("The model does not contain population-level effects.")
  }
  out <- posterior_samples(object, pars = fpars, exact_match = TRUE)
  out <- do.call(cbind, lapply(estimate, get_estimate, samples = out, ...))
  rownames(out) <- gsub(fixef_pars(), "", fpars)
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
  contains_samples(object)
  pars <- parnames(object)
  fpars <- pars[grepl(fixef_pars(), pars)]
  if (!length(fpars)) {
    stop2("The model does not contain population-level effects.")
  }
  samples <- posterior_samples(object, pars = fpars, exact_match = TRUE)
  names(samples) <- sub(fixef_pars(), "", names(samples))
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
#'  for the group-level effects, either \code{"mean"} or \code{"median"}.
#' @param var logical; indicating if the covariance matrix 
#'  for each group-level effects should be computed.
#' @param ... Further arguments to be passed to the function 
#'  specified in \code{estimate}.
#'
#' @return A list of matrices (one per grouping factor), 
#'  with factor levels as row names and 
#'  group-level effects as column names.
#'     
#' @author Paul-Christian Buerkner \email{paul.buerkner@gmail.com}   
#'   
#' @examples
#' \dontrun{
#' fit <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1+Trt_c|visit), 
#'            data = epilepsy, family = "poisson", chains = 1)
#' ## group-level means with corresponding covariances
#' rf <- ranef(fit, var = TRUE)
#' attr(rf, "var")
#' ## group-level medians
#' ranef(fit, estimate = "median")                                                        
#' }
#' 
#' @method ranef brmsfit
#' @export
#' @export ranef
#' @importFrom nlme ranef
ranef.brmsfit <- function(object, estimate = c("mean", "median"), 
                          var = FALSE, ...) {
  contains_samples(object)
  object <- restructure(object)
  estimate <- match.arg(estimate)
  if (!nrow(object$ranef)) {
    stop2("The model does not contain group-level effects.")
  }
  pars <- parnames(object)
  
  .ranef <- function(group, nlpar = "") {
    # get group-level effects of a grouping factor
    # Args:
    #   group: name of a grouping factor
    #   nlpar: name of a non-linear parameter
    rnames <- object$ranef[object$ranef$group == group & 
                           object$ranef$nlpar == nlpar, "coef"]
    usc_nlpar <- usc(usc(nlpar))
    rpars <- pars[grepl(paste0("^r_", group, usc_nlpar, "\\["), pars)]
    if (!length(rpars)) {
      return(NULL)
    }
    rdims <- object$fit@sim$dims_oi[[paste0("r_", group, usc_nlpar)]]
    if (is.na(rdims[2])) rdims[2] <- 1
    levels <- attr(object$ranef, "levels")[[group]]
    if (is.null(levels)) {
      # avoid error in dimnames if levels are NULL 
      # for backwards compatibility with brms < 0.5.0 
      levels <- seq_len(rdims[1])
    }
    rs <- posterior_samples(object, pars = rpars, exact_match = TRUE)
    rs_array <- array(dim = c(rdims[1], rdims[2], nrow(rs)))
    k <- 0
    for (j in seq_len(rdims[2])) {
      for (l in seq_len(rdims[1])) {
        k <- k + 1
        rs_array[l, j, ] <- rs[, k]
      }
    }
    out <- get_estimate(estimate, samples = rs_array, margin = 1:2, ...)
    colnames(out) <- rnames
    if (var) {
      Var <- array(dim = c(rep(rdims[2], 2), rdims[1]), 
                   dimnames = list(rnames, rnames, seq_len(rdims[1])))
      if (rdims[2] == 1L) {
        for (j in seq_len(rdims[1])) {
          Var[, , j] <- var(rs_array[j, 1, ]) 
        }
      } else {
        for (j in seq_len(rdims[1])) {
          Var[, , j] <- cov(t(rs_array[j, , ]))
        }
      }
      dimnames(Var)[[3]] <- levels
      attr(out, "var") <- Var
    }
    rownames(out) <- levels
    if (nchar(nlpar)) {
      attr(out, "nlpar") <- nlpar
    }
    return(out)
  }
  
  group_nlpar <- unique(object$ranef[, c("group", "nlpar")])
  ranef <- named_list(group_nlpar$group)
  for (i in seq_along(ranef)) {
    ranef[[i]] <- .ranef(group = group_nlpar$group[i], 
                         nlpar = group_nlpar$nlpar[i])
  }
  rmNULL(ranef) 
} 

#' Extract model coefficients
#'
#' Extract model coefficients, which are the sum of population-level 
#' effects and corresponding group-level effects
#' 
#' @param object An object of class \code{brmsfit}
#' @inheritParams ranef.brmsfit
#'
#' @return A list of matrices (one per grouping factor and 
#'  non-linear parameter) with factor levels as row names and 
#'  coefficients as column names.
#'  
#' @author Paul-Christian Buerkner \email{paul.buerkner@gmail.com}
#'  
#' @examples
#' \dontrun{
#' fit <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1+Trt_c|visit), 
#'            data = epilepsy, family = "poisson", chains = 1)
#' ## extract population and group-level coefficients separately
#' fixef(fit)
#' ranef(fit)
#' ## extract combined coefficients 
#' coef(fit)
#' }
#' 
#' @export
coef.brmsfit <- function(object, estimate = c("mean", "median"), ...) {
  contains_samples(object)
  object <- restructure(object)
  estimate <- match.arg(estimate)
  if (!nrow(object$ranef)) {
    stop2("No group-level effects detected. Call method ", 
          "'fixef' to access population-level effects.")
  }
  
  .coef <- function(ranef, fixef) {
    # helper function to combine group and population-level effects
    all_ranef_names <- unique(ulapply(ranef, colnames))
    fixef_no_digits <- get_matches("^[^\\[]+", rownames(fixef))
    miss_fixef <- setdiff(all_ranef_names, rownames(fixef))
    miss_fixef_no_digits <- get_matches("^[^\\[]+", miss_fixef)
    new_fixef <- named_list(miss_fixef)
    for (k in seq_along(miss_fixef)) {
      # digits occur in ordinal models with category specific effects
      match_fixef <- match(miss_fixef_no_digits[k], rownames(fixef))
      if (!is.na(match_fixef)) {
        new_fixef[[k]] <- fixef[match_fixef, 1]
      } else if (!miss_fixef[k] %in% fixef_no_digits) {
        new_fixef[[k]] <- 0
      }
    }
    rm_fixef <- rownames(fixef) %in% miss_fixef_no_digits
    fixef <- fixef[!rm_fixef, , drop = FALSE]
    fixef <- do.call(rbind, c(list(fixef), rmNULL(new_fixef)))
    coef <- ranef
    for (i in seq_along(coef)) {
      ranef_names <- colnames(coef[[i]])
      ranef_no_digits <- get_matches("^[^\\[]+", ranef_names)
      miss_ranef <- setdiff(rownames(fixef), ranef_names)
      miss_ranef_no_digits <- get_matches("^[^\\[]+", miss_ranef)
      new_ranef <- named_list(miss_ranef)
      for (k in seq_along(miss_ranef)) {
        # digits occur in ordinal models with category specific effects
        match_ranef <- match(miss_ranef_no_digits[k], ranef_names)
        if (!is.na(match_ranef)) {
          new_ranef[[k]] <- coef[[i]][, match_ranef]
        } else if (!miss_ranef[k] %in% ranef_no_digits) {
          new_ranef[[k]] <- 0
        }
      }
      rm_ranef <- ranef_names %in% miss_ranef_no_digits
      coef[[i]] <- coef[[i]][, !rm_ranef, drop = FALSE]
      coef[[i]] <- do.call(cbind, c(list(coef[[i]]), rmNULL(new_ranef)))
      for (nm in colnames(coef[[i]])) {
        # correct the sign of thresholds in ordinal models
        sign <- ifelse(grepl("^Intercept\\[[[:digit:]]+\\]$", nm), -1, 1)
        coef[[i]][, nm] <- coef[[i]][, nm] + sign * fixef[nm, 1]
      }
    }
    return(coef)
  }
  
  fixef <- fixef(object, estimate = estimate, ...)
  ranef <- ranef(object, estimate = estimate, ...)
  nlpars <- ulapply(ranef, attr, "nlpar")
  if (length(nlpars)) {
    # do not combine effects of different nlpars
    unique_nlpars <- unique(nlpars)
    coef <- named_list(unique_nlpars)
    for (p in unique_nlpars) {
      ranef_temp <- ranef[nlpars %in% p]
      rx <- paste0("^", p, "_")
      take_rows <- grepl(rx, rownames(fixef))
      fixef_temp <- fixef[take_rows, , drop = FALSE]
      rownames(fixef_temp) <- sub(rx, "", rownames(fixef_temp))
      coef[[p]] <- .coef(ranef_temp, fixef_temp)
      for (i in seq_along(coef[[p]])) {
        attr(coef[[p]][[i]], "nlpar") <- p
      }
    }
    coef <- unlist(unname(coef), recursive = FALSE)
  } else {
    coef <- .coef(ranef, fixef)
  }
  coef
}

#' Extract variance and correlation components
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
#' @param estimate A character vector specifying which coefficients 
#'  (e.g. \code{"mean"}, \code{"median"}, \code{"sd"}, or \code{"quantile"})
#'  should be calculated for the extracted parameters.
#' @param as.list logical; Indicates if covariance 
#'  and correlation matrices should be returned as 
#'  lists of matrices (\code{TRUE}; the default)
#'  or as 3-dimensional arrays (\code{FALSE}).
#'  This option is kept for backwards compatibility and 
#'  we recommend not to change it.
#' @param sigma Ignored (included for compatibility with 
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
#' @author Paul-Christian Buerkner \email{paul.buerkner@gmail.com}
#' 
#' @examples
#' \dontrun{
#' fit <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1+Trt_c|visit), 
#'              data = epilepsy, family = "poisson", chains = 1)
#' ## return the means of group-level covariances
#' (vc <- VarCorr(fit))
#' as.data.frame(vc)
#' 
#' ## return 2.5% and 97.5% quantiles of group-level covariances
#' VarCorr(fit, estimate = "quantile", probs = c(0.025, 0.975))
#' }
#' 
#' @method VarCorr brmsfit
#' @import abind abind
#' @importFrom nlme VarCorr
#' @export VarCorr
#' @export
VarCorr.brmsfit <- function(x, sigma = 1, estimate = "mean", 
                            as.list = TRUE, ...) {
  contains_samples(x)
  x <- restructure(x)
  if (!(nrow(x$ranef) || any(grepl("^sigma($|_)", parnames(x))))) {
    stop2("The model does not contain covariance matrices.")
  }
  x <- restructure(x)
  
  # extracts samples for sd, cor and cov
  extract <- function(p) {
    nr <- length(p$sd_pars)
    sd <- posterior_samples(x, pars = p$sd_pars, exact_match = TRUE)
    nsamples <- nrow(sd)
    out <- list(
      sd = do.call(cbind, lapply(estimate, get_estimate, samples = sd, ...))
    )
    rownames(out$sd) <- p$rnames 
    # calculate correlation and covariance matrices
    found_cor_pars <- intersect(p$cor_pars, parnames(x))
    if (length(found_cor_pars)) {
      cor <- posterior_samples(x, pars = paste0("^", found_cor_pars, "$"))
      if (length(found_cor_pars) < length(p$cor_pars)) { 
        # some correlations are missing and will be replaced by 0
        cor_all <- as.data.frame(
          matrix(0, nrow = nrow(cor), ncol = length(p$cor_pars))
        )
        names(cor_all) <- p$cor_pars
        for (i in seq_len(ncol(cor_all))) {
          found <- match(names(cor_all)[i], names(cor))
          if (!is.na(found)) {
            # correlation was estimated
            cor_all[, i] <- cor[, found]
          }
        }
        cor <- cor_all
      }
    } else {
      cor <- NULL
    } 
    # get_cov_matrix and array2list can be found in brmsfit-helpers.R
    matrices <- get_cov_matrix(sd = sd, cor = cor) 
    out$cor <- abind(lapply(
      estimate, get_estimate, samples = matrices$cor, 
      margin = c(2, 3), to.array = TRUE, ...
    ))
    out$cov <- abind(lapply(
      estimate, get_estimate, samples = matrices$cov, 
      margin = c(2,3), to.array = TRUE, ...
    )) 
    if (length(p$rnames) > 1) {
      dimnames(out$cor) <- list(p$rnames, p$rnames, dimnames(out$cor)[[3]])
      dimnames(out$cov) <- dimnames(out$cor)   
    }
    if (as.list) {
      out$cor <- lapply(array2list(out$cor), function(x)
        if (is.null(dim(x))) 
          structure(matrix(x), dimnames = list(p$rnames, p$rnames)) 
        else x
      )
      out$cov <- lapply(array2list(out$cov), function(x)
        if (is.null(dim(x))) 
          structure(matrix(x), dimnames = list(p$rnames, p$rnames)) 
        else x
      )
    }
    out
  }
  
  family <- family(x)
  bterms <- parse_bf(x$formula, family = family)
  if (nrow(x$ranef)) {
    get_names <- function(group) {
      # get names of group-level parameters
      r <- x$ranef[x$ranef$group == group, ]
      rnames <- paste0(usc(r$nlpar, "suffix"), r$coef)
      cor_type <- paste0("cor_", group)
      sd_pars <- paste0("sd_", group, "__", rnames)
      cor_pars <- get_cornames(rnames, type = cor_type, brackets = FALSE)
      nlist(rnames, type = cor_type, sd_pars, cor_pars)
    }
    group <- unique(x$ranef$group)
    p <- lapply(group, get_names)
  } else {
    p <- group <- NULL
  } 
  # special treatment of residuals variances in linear models
  has_sigma <- has_sigma(family, bterms, incmv = TRUE)
  if (has_sigma && !"sigma" %in% names(bterms$auxpars)) {
    cor_pars <- get_cornames(
      bterms$response, type = "rescor", brackets = FALSE
    )
    p <- lc(p, 
      list(rnames = bterms$response, cor_pars = cor_pars,
           sd_pars = c("sigma", paste0("sigma_", bterms$response)))
    )
    group <- c(group, "RESIDUAL")
  } 
  VarCorr <- lapply(p, extract)
  names(VarCorr) <- group
  if (as.list) {
    class(VarCorr) <- "brmsVarCorr" 
  }
  VarCorr
}

#' @export
model.frame.brmsfit <- function(formula, ...) {
  if (!is.data.frame(formula$data)) {
    stop2("Cannot extract the model.frame for this model.")
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
  pars <- use_alias(pars, parameters, default = NA)
  add_chain <- use_alias(add_chain, add_chains, default = FALSE)
  contains_samples(x)
  pars <- extract_pars(pars, all_pars = parnames(x), 
                       exact_match = exact_match, ...)
  
  # get basic information on the samples 
  iter <- x$fit@sim$iter
  warmup <- x$fit@sim$warmup
  thin <- x$fit@sim$thin
  chains <- x$fit@sim$chains
  final_iter <- ceiling((iter - warmup) / thin)
  samples_taken <- seq(warmup + 1, iter, thin)
  
  if (length(pars)) {
    samples <- as.data.frame(x$fit, pars = pars)
    if (add_chain) {
      # name the column 'chain' not 'chains' (#32)
      samples$chain <- factor(rep(1:chains, each = final_iter))
      samples$iter <- rep(samples_taken, chains)
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
  pars <- extract_pars(pars, all_pars = parnames(x),
                       exact_match = exact_match, ...)
  if (combine_chains) {
    if (inc_warmup) {
      stop2("Cannot include warmup samples when 'combine_chains' is TRUE.")
    }
    out <- as.matrix(x$fit, pars)
    mcpar <- c(x$fit@sim$warmup * x$fit@sim$chain + 1, 
               x$fit@sim$iter * x$fit@sim$chain, x$fit@sim$thin)
    attr(out, "mcpar") <- mcpar
    class(out) <- "mcmc"
  } else {
    ps <- extract(x$fit, pars, permuted = FALSE, 
                  inc_warmup = inc_warmup)
    mcpar <- c(if (inc_warmup) 1 else x$fit@sim$warmup + 1, 
               x$fit@sim$iter, x$fit@sim$thin)
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
prior_samples.brmsfit <- function(x, pars = NA, parameters = NA, ...) {
  pars <- use_alias(pars, parameters, default = NA)
  if (!anyNA(pars) && !is.character(pars)) {
    stop2("Argument 'pars' must be a character vector.")
  }
  par_names <- parnames(x)
  prior_names <- unique(par_names[grepl("^prior_", par_names)])
  if (length(prior_names)) {
    samples <- posterior_samples(x, pars = prior_names, exact_match = TRUE)
    names(samples) <- sub("^prior_", "", prior_names)
    if (!anyNA(pars)) {
      .prior_samples <- function(par) {
        # get prior samples for parameter par 
        if (grepl("^b_Intercept", par) && !par %in% names(samples)) {
          # population-level intercepts do not inherit priors
          out  <- NULL
        } else {
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
        }
        return(out)
      }
      samples <- data.frame(rmNULL(lapply(pars, .prior_samples)), 
                            check.names = FALSE)
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
#' @param waic Logical; Indicating if the WAIC should be computed
#'   (this will take some time for larger models). 
#'   Default is \code{FALSE}.
#' @param priors Logical; Indicating if priors should be included 
#'   in the summary. Default is \code{FALSE}.
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
summary.brmsfit <- function(object, waic = FALSE, priors = FALSE,
                            use_cache = TRUE, ...) {
  object <- restructure(object, rstr_summary = use_cache)
  bterms <- parse_bf(formula(object), family = family(object))
  out <- brmssummary(formula = formula(object), 
                     family = family(object), 
                     data.name = object$data.name, 
                     group = unique(object$ranef$group), 
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
    allow_waic <- !nrow(object$ranef) || any(grepl("^r_", parnames(object)))
    if (waic && allow_waic) {
      out$WAIC <- SW(WAIC(object)$waic)
    }
    if (priors) {
      out$prior <- prior_summary(object, all = FALSE)
    }
    
    pars <- parnames(object)
    meta_pars <- object$fit@sim$pars_oi
    meta_pars <- meta_pars[!grepl("^(r|s|Xme|prior)_", meta_pars)]
    fit_summary <- summary(object$fit, pars = meta_pars,
                           probs = c(0.025, 0.975),
                           use_cache = use_cache)$summary
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
        warning2(msg)
      }
    } else {
      colnames(fit_summary) <- c("Estimate", "Est.Error", 
                                 "l-95% CI", "u-95% CI")
    }
    
    # fixed effects summary
    fe_pars <- pars[grepl(fixef_pars(), pars)]
    out$fixed <- fit_summary[fe_pars, , drop = FALSE]
    rownames(out$fixed) <- gsub(fixef_pars(), "", fe_pars)
    
    # summary of family specific parameters
    is_mv_par <- apply(sapply(c("^sigma_", "^rescor_"), grepl, pars), 1, any)
    spec_pars <- pars[pars %in% c(auxpars(), "delta") | is_mv_par]
    out$spec_pars <- fit_summary[spec_pars, , drop = FALSE]
    if (is_linear(family(object)) && length(bterms$response) > 1L) {
      sigma_names <- paste0("sigma(", bterms$response, ")")
      rescor_names <- get_cornames(bterms$response, type = "rescor")   
      spec_pars[grepl("^sigma_", spec_pars)] <- sigma_names
      spec_pars[grepl("^rescor_", spec_pars)] <- rescor_names 
    }    
    rownames(out$spec_pars) <- spec_pars
    
    # summary of autocorrelation effects
    cor_pars <- pars[grepl("^ar|^ma|^sigmaLL$", pars)]
    out$cor_pars <- fit_summary[cor_pars, , drop = FALSE]
    rownames(out$cor_pars) <- cor_pars
    
    # summary of group-level effects
    for (g in out$group) {
      r <- object$ranef[object$ranef$group == g, ]
      nlpar_usc <- usc(r$nlpar, "suffix")
      rnames <- paste0(nlpar_usc, r$coef)
      sd_pars <- paste0("sd_", g, "__", rnames)
      sd_names <- paste0("sd", "(", rnames ,")")
      # construct correlation names
      type <- paste0("cor_", g)
      all_cor_pars <- get_cornames(rnames, brackets = FALSE, type = type)
      take <- all_cor_pars %in% parnames(object)
      cor_pars <- all_cor_pars[take]
      cor_names <- get_cornames(rnames)[take]
      # extract sd and cor parameters from the summary
      out$random[[g]] <- fit_summary[c(sd_pars, cor_pars), , drop = FALSE]
      rownames(out$random[[g]]) <- c(sd_names, cor_names)
    }
    
    # summary of smooths
    spline_pars <- pars[grepl("^sds_", pars)]
    if (length(spline_pars)) {
      out$splines <- fit_summary[spline_pars, , drop = FALSE]
      rownames(out$splines) <- paste0(gsub("^sds_", "sds(", spline_pars), ")")
    }
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
  bf(x$formula, ...)
}

#' @export
family.brmsfit <- function(object, ...) {
  if (is.character(object$family)) {
    # brms <= 0.6.0
    family <- brmsfamily(object$family, link = object$link) 
  } else {
    family <- object$family
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
    object <- restructure(object)
    new_formula <- update_re_terms(object$formula, dots$re_formula)
    dots$control$old_cat <- is_old_categorical(object)
    prior_only <- attr(object$prior, "prior_only")
    sample_prior <- ifelse(isTRUE(prior_only), "only", FALSE)
    args <- list(formula = new_formula, data = model.frame(object), 
                 family = object$family, prior = object$prior, 
                 autocor = object$autocor, cov_ranef = object$cov_ranef, 
                 knots = attr(model.frame(object), "knots"),
                 sample_prior = sample_prior)
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
  contains_samples(x)
  shinystan::launch_shinystan(x$fit, rstudio = rstudio, ...)
}

#' Trace and Density Plots for MCMC Samples
#' 
#' @param x An object of class \code{brmsfit}.
#' @param pars Names of the parameters to plot, as given by a character vector 
#'   or a regular expression. By default, all parameters except 
#'   for group-level and smooth effects are plotted. 
#' @param parameters A deprecated alias of \code{pars}
#' @param combo A character vector with at least two elements. 
#'   Each element of \code{combo} corresponds to a column in the resulting 
#'   graphic and should be the name of one of the available 
#'   \code{link[bayesplot:MCMC-overview]{MCMC}} functions 
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
plot.brmsfit <- function(x, pars = NA, parameters = NA, 
                         combo = c("dens", "trace"), N = 5, 
                         exact_match = FALSE, theme = NULL, 
                         plot = TRUE, ask = TRUE, 
                         newpage = TRUE, ...) {
  dots <- list(...)
  pars <- use_alias(pars, parameters, default = NA)
  plot <- use_alias(plot, dots$do_plot)
  contains_samples(x)
  if (!is_wholenumber(N) || N < 1) {
    stop2("Argument 'N' must be a positive integer.")
  }
  if (!is.character(pars)) {
    pars <- default_plot_pars()
    exact_match <- FALSE
  }
  samples <- posterior_samples(x, pars = pars, add_chain = TRUE,
                               exact_match = exact_match)
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
    plots[[i]] <- bayesplot::mcmc_combo(sub_samples, combo = combo, 
                                        gg_theme = theme, ...)
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
      samples <- posterior_samples(object, pars, add_chain = TRUE,
                                   exact_match = exact_match)
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
#'   Currently, the following plots (as names) are implemented
#'   in \pkg{bayesplot} (v1.1.0):
#'   \code{boxplot}, \code{dens} \code{dens_overlay} (default), 
#'   \code{ecdf_overlay}, \code{error_binned}, 
#'   \code{error_hist}, \code{error_scatter}, 
#'   \code{error_scatter_avg}, \code{error_scatter_avg_vs_x},
#'   \code{freqpoly}, \code{freqpoly_grouped}, \code{hist},
#'   \code{intervals}, \code{intervals_grouped}, 
#'   \code{ribbon}, \code{ribbon_grouped}, \code{scatter},
#'   \code{scatter_avg}, \code{scatter_avg_grouped}, 
#'   \code{stat}, \code{stat_2d}, \code{stat_freqpoly_grouped},
#'   \code{stat_grouped}, and \code{violin_grouped}.
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
#' @param x Optional name of a variable in the model. 
#'  Only used for ppc types having an \code{x} argument 
#'  and ignored otherwise.
#' @param ... Further arguments passed to the ppc functions
#'   of the \pkg{\link[bayesplot:bayesplot]{bayesplot}} package.
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
#'             + (1|patient) + (1|visit),
#'             data = epilepsy, family = poisson())
#' 
#' pp_check(fit)  # shows dens_overlay plot by default
#' pp_check(fit, type = "error_hist", nsamples = 11)
#' pp_check(fit, type = "scatter_avg", nsamples = 100)
#' pp_check(fit, type = "stat_2d")
#' }
#' 
#' @importFrom bayesplot pp_check
#' @export pp_check
#' @export
pp_check.brmsfit <- function(object, type, nsamples, group = NULL,
                             x = NULL, newdata = NULL, 
                             re_formula = NULL, allow_new_levels = FALSE,
                             sample_new_levels = "uncertainty", 
                             incl_autocor = TRUE, subset = NULL, 
                             ntrys = 5, ...) {
  if (missing(type)) {
    type <- "dens_overlay"
  }
  if (length(type) != 1L) {
    stop2("Argument 'type' must be of length 1.")
  }
  if (!is.null(group) && length(group) != 1L) {
    stop2("Argument 'group' must be of length 1.")
  }
  if (!is.null(x) && length(x) != 1L) {
    stop2("Argument 'x' must be of length 1.")
  }
  ppc_funs <- as.character(bayesplot::available_ppc(""))
  valid_ppc_types <- sub("^ppc_", "", ppc_funs)
  if (!type %in% valid_ppc_types) {
    stop2("Type '", type, "' is not a valid ppc type. Valid types are: \n", 
          collapse_comma(valid_ppc_types))
  }
  ppc_fun <- get(paste0("ppc_", type), pos = asNamespace("bayesplot"))
  # validate argument 'group'
  object <- restructure(object)
  not_num <- !sapply(model.frame(object), is.numeric)
  valid_groups <- c(
    names(model.frame(object))[not_num],
    parse_time(object$autocor$formula)$group,
    object$ranef$group
  )
  valid_groups <- unique(valid_groups[nzchar(valid_groups)])
  if (!is.null(group) && !group %in% valid_groups) {
    stop2("Group '", group, "' is not a valid grouping factor. ",
          "Valid groups are: \n", collapse_comma(valid_groups))
  }
  is_group_type <- "group" %in% names(formals(ppc_fun))
  if (is.null(group) && is_group_type) {
    stop2("Argument 'group' is required for ppc type '", type, "'.")
  }
  is_vs_x_type <- "x" %in% names(formals(ppc_fun))
  if (is_vs_x_type) {
    if (is.null(x)) {
      stop2("Argument 'x' is required for ppc type '", type, "'.")
    }
    bterms <- parse_bf(formula(object), family = family(object))
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
      "intervals", "intervals_grouped", "ribbon", 
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
    newdata, fit = object, re_formula, 
    allow_new_levels, incl_autocor,
    check_response = TRUE
  )
  standata <- do.call(amend_newdata, newd_args)
  y <- as.vector(standata$Y)
  if (!is.null(standata$cens)) {
    warning2("Posterior predictive checks may not be ", 
             "meaningful for censored models.")
  }
  pred_args <- nlist(
    object, newdata, re_formula, allow_new_levels, 
    sample_new_levels, incl_autocor, nsamples, subset, 
    ntrys, sort = TRUE, summary = FALSE
  )
  yrep <- as.matrix(do.call(method, pred_args))
  if (family(object)$family %in% "binomial") {
    # use success proportions following Gelman and Hill (2006)
    y <- y / standata$trials
    yrep <- yrep / as_draws_matrix(standata$trials, dim = dim(yrep))
  }
  ppc_args <- list(y, yrep, ...)
  old_order <- attr(standata, "old_order")
  if (!is.null(group)) {
    group_var <- model.frame(object)[[group]]
    if (!is.null(old_order)) {
      group_var <- group_var[order(old_order)]
    }
    ppc_args$group <- group_var
  }
  if (!is.null(x)) {
    x_var <- model.frame(object)[[x]]
    if (!is.null(old_order)) {
      x_var <- x_var[order(old_order)]
    }
    ppc_args$x <- x_var
  }
  do.call(ppc_fun, ppc_args)
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
                                     re_formula = NA, robust = TRUE, 
                                     probs = c(0.025, 0.975),
                                     method = c("fitted", "predict"), 
                                     surface = FALSE, resolution = 100,
                                     too_far = 0, ...) {
  method <- match.arg(method)
  dots <- list(...)
  conditions <- use_alias(conditions, dots[["data"]])
  surface <- use_alias(surface, dots[["contour"]])
  dots[["data"]] <- dots[["contour"]] <- NULL
  contains_samples(x)
  x <- restructure(x)
  new_formula <- update_re_terms(x$formula, re_formula = re_formula)
  bterms <- parse_bf(new_formula, family = x$family)
  if (is_linear(x$family) && length(bterms$response) > 1L) {
    stop2("Marginal plots are not yet implemented for multivariate models.")
  } else if (is_categorical(x$family)) {
    stop2("Marginal plots are not yet implemented for categorical models.")
  } else if (is_ordinal(x$family)) {
    warning2("Predictions are treated as continuous variables ", 
             "in marginal plots, \nwhich is likely an invalid ", 
             "assumption for family ", x$family$family, ".")
  }
  rsv_vars <- rsv_vars(
    x$family, nresp = length(bterms$response),
    rsv_intercept = attr(bterms$fe, "rsv_intercept"),
    old_mv = attr(bterms$formula, "old_mv")
  )
  if (is.null(effects)) {
    effects <- get_all_effects(bterms, rsv_vars = rsv_vars)
    if (!length(effects)) {
      stop2("No valid effects detected.")
    }
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
      warning2("Some specified effects are invalid for this model: ",
               collapse_comma(invalid), "\nValid effects are ", 
               "(combinations of): ", collapse_comma(ae_coll))
    }
    effects <- unique(effects[sort(matches)])
    if (!length(effects)) {
      stop2("All specified effects are invalid for this model.\n", 
            "Valid effects are (combinations of): ", 
            collapse_comma(ae_coll))
    }
  }
  if (length(probs) != 2L) {
    stop2("Arguments 'probs' must be of length 2.")
  }
  
  conditions <- prepare_conditions(
    x, conditions = conditions, effects = effects, 
    re_formula = re_formula, rsv_vars = rsv_vars
  )
  int_effects <- c(get_effect(bterms, "mo"), 
                   rmNULL(bterms[c("trials", "cat")]))
  int_vars <- unique(ulapply(int_effects, all.vars))
  mf <- model.frame(x)
  results <- list()
  for (i in seq_along(effects)) {
    marg_data <- mf[, effects[[i]], drop = FALSE]
    marg_args <- nlist(data = marg_data, conditions, 
                       int_vars, surface, resolution)
    marg_data <- do.call(prepare_marg_data, marg_args)
    if (surface && length(effects[[i]]) == 2L && too_far > 0) {
      # exclude prediction grid points too far from data
      ex_too_far <- mgcv::exclude.too.far(
        g1 = marg_data[[effects[[i]][1]]], 
        g2 = marg_data[[effects[[i]][2]]], 
        d1 = x$data[, effects[[i]][1]],
        d2 = x$data[, effects[[i]][2]],
        dist = too_far)
      marg_data <- marg_data[!ex_too_far, ]  
    }
    # make sure numeric variables come first
    effects[[i]] <- attr(marg_data, "effects")
    args <- list(
      x, newdata = marg_data, re_formula = re_formula,
      allow_new_levels = TRUE, incl_autocor = FALSE, 
      probs = probs, robust = robust
    )
    args <- c(args, dots)
    if (is_ordinal(x$family) || is_categorical(x$family)) {
      args$summary <- FALSE 
      marg_res <- do.call(method, args)
      if (method == "fitted") {
        for (k in seq_len(dim(marg_res)[3])) {
          marg_res[, , k] <- marg_res[, , k] * k
        }
        marg_res <- do.call(cbind, lapply(seq_len(dim(marg_res)[2]), 
                            function(s) rowSums(marg_res[, s, ])))
      } 
      marg_res <- get_summary(marg_res, probs = probs, robust = robust)
    } else {
      marg_res <- do.call(method, args)
    }
    colnames(marg_res) <- c("estimate__", "se__", "lower__", "upper__")
    
    types <- attr(marg_data, "types")
    both_numeric <- length(types) == 2L && all(types == "numeric")
    if (both_numeric && !surface) {
      # can only be converted to factor after having called method
      if (isTRUE(attr(marg_data, "mono")[2])) {
        labels <- c("Median - MAD", "Median", "Median + MAD")
      } else {
        labels <- c("Mean - SD", "Mean", "Mean + SD") 
      }
      marg_data[[effects[[i]][2]]] <- 
        factor(marg_data[[effects[[i]][2]]], labels = labels)
    }
    marg_res = cbind(marg_data, marg_res)
    attr(marg_res, "response") <- as.character(x$formula$formula[2])
    attr(marg_res, "effects") <- effects[[i]]
    attr(marg_res, "surface") <- both_numeric && surface
    point_args <- nlist(mf, effects = effects[[i]], conditions,
                        groups = get_re(bterms)$group, family = x$family)
    attr(marg_res, "points") <- do.call(make_point_frame, point_args)
    results[[paste0(effects[[i]], collapse = ":")]] <- marg_res
  }
  class(results) <- "brmsMarginalEffects"
  results
}

#' @rdname marginal_smooths
#' @export
marginal_smooths.brmsfit <- function(x, smooths = NULL,
                                     probs = c(0.025, 0.975),
                                     resolution = 100, too_far = 0,
                                     subset = NULL, nsamples = NULL,
                                     ...) {
  contains_samples(x)
  x <- restructure(x)
  mf <- model.frame(x)
  conditions <- prepare_conditions(x)
  smooths <- rename(as.character(smooths), " ", "")
  bterms <- parse_bf(formula(x), family = family(x))
  lee <- list()
  if (length(bterms$response) > 1L) {
    for (r in bterms$response) {
      lee <- c(lee, setNames(bterms$auxpars["mu"], r))
    }
    bterms$auxpars[["mu"]] <- NULL
  }
  for (ap in names(bterms$auxpars)) {
    bt <- bterms$auxpars[ap]
    if (is.btnl(bt[[1]])) {
      lee <- c(lee, bt[[1]]$nlpars)
    } else {
      lee <- c(lee, bt)
    }
  }
  subset <- subset_samples(x, subset, nsamples)
  nsamples <- nsamples(x, subset = subset)
  args <- nlist(
    fit = x, allow_new_levels = TRUE,
    subset, nsamples, incl_autocor = FALSE, 
    smooths_only = TRUE 
  )
  too_many_covars <- FALSE
  results <- list()
  for (k in seq_along(lee)) {
    # loop over elements that may contain smooth terms
    sm_labels <- get_sm_labels(lee[[k]])
    sm_labels_by <- get_sm_labels(lee[[k]], data = mf)
    covars <- get_sm_labels(lee[[k]], covars = TRUE, combine = FALSE)
    for (i in seq_along(sm_labels)) {
      # loop over smooth terms and compute their predictions
      covars_no_by_factor <- covars[[i]]
      byvars <- attr(covars, "byvars")[[i]]
      byfactors <- !sapply(mf[, byvars, drop = FALSE], is.numeric)
      byfactors <- byvars[byfactors]
      covars_no_byfactor <- setdiff(covars[[i]], byfactors)
      ncovars <- length(covars_no_byfactor)
      if (ncovars > 2L) {
        too_many_covars <- TRUE
      }
      include_smooth <- !length(smooths) || sm_labels[[i]] %in% smooths
      if (include_smooth && ncovars <= 2L) {
        values <- named_list(covars[[i]])
        for (cv in names(values)) {
          if (is.numeric(mf[[cv]])) {
            values[[cv]] <- seq(min(mf[[cv]]), max(mf[[cv]]), 
                                length.out = resolution)
          } else {
            values[[cv]] <- levels(factor(mf[[cv]]))
          }
        }
        newdata <- expand.grid(values)
        if (ncovars == 2L && too_far > 0) {
          # exclude prediction grid points too far from data
          ex_too_far <- mgcv::exclude.too.far(
            g1 = newdata[[covars_no_byfactor[1]]], 
            g2 = newdata[[covars_no_byfactor[2]]], 
            d1 = x$data[, covars_no_byfactor[1]],
            d2 = x$data[, covars_no_byfactor[2]],
            dist = too_far
          )
          newdata <- newdata[!ex_too_far, ]  
        }
        other_vars <- setdiff(names(conditions), covars[[i]])
        newdata[, other_vars] <- conditions[1, other_vars]
        # prepare draws for linear_predictor
        more_args <- nlist(x = lee[[k]], newdata, nlpar = names(lee)[k])
        draws <- do.call(extract_draws, c(args, more_args))
        J <- which(attr(sm_labels_by, "termnum") == i)
        scs <- unlist(attr(draws$data[["X"]], "smooth_cols")[J])
        draws$data[["X"]] <- draws$data[["X"]][, scs, drop = FALSE]
        draws[["b"]] <- draws[["b"]][, scs, drop = FALSE]
        draws[["Zs"]] <- draws[["Zs"]][J] 
        draws[["s"]] <- draws[["s"]][J]
        eta <- get_eta(draws = draws, i = NULL)
        eta <- get_summary(eta, robust = TRUE, probs = probs)
        colnames(eta) <- c("estimate__", "se__", "lower__", "upper__")
        res <- cbind(newdata[, covars[[i]], drop = FALSE], eta)
        if (length(byfactors)) {
          res$cond__ <- Reduce(paste_colon, res[, byfactors, drop = FALSE]) 
        }
        response <- sm_labels[[i]]
        if (isTRUE(nzchar(names(lee)[k]))) {
          response <- paste0(names(lee)[k], ": ", response)
        }
        attr(res, "response") <- response
        attr(res, "effects") <- covars_no_byfactor
        attr(res, "surface") <- ncovars == 2L
        attr(res, "points") <- mf[, covars[[i]], drop = FALSE]
        results[[response]] <- res
      }
    }
  }
  if (!length(results)) {
    stop2("No valid smooth terms found in the model.")
  }
  if (too_many_covars) {
    warning2("Smooth terms with more than two covariates ",  
             "are not yet supported by 'marginal_smooths'.")
  }
  attr(results, "smooths_only") <- TRUE
  class(results) <- "brmsMarginalEffects"
  results
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
#'   If \code{NULL} (default), the orginal data of the model is used.
#' @param re_formula formula containing group-level effects 
#'   to be considered in the prediction. 
#'   If \code{NULL} (default), include all group-level effects; 
#'   if \code{NA}, include no group-level effects.
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
#'   a warning will occure suggesting to increase argument \code{ntrys}.
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
                            incl_autocor = TRUE, subset = NULL, 
                            nsamples = NULL, sort = FALSE,
                            ntrys = 5, summary = TRUE, robust = FALSE,
                            probs = c(0.025, 0.975), ...) {
  contains_samples(object)
  object <- restructure(object)
  draws_args <- nlist(
    x = object, newdata, re_formula, incl_autocor, 
    allow_new_levels, sample_new_levels, subset, nsamples
  )
  draws <- do.call(extract_draws, draws_args)
  if (is.list(draws$mu[["mv"]])) {
    draws$mu <- get_eta(draws$mu)
  }
  for (ap in intersect(auxpars(), names(draws))) {
    if (is.list(draws[[ap]])) {
      draws[[ap]] <- get_auxpar(draws[[ap]])
    }
  }
  # see predict.R
  predict_fun <- paste0("predict_", draws$f$family)
  predict_fun <- get(predict_fun, asNamespace("brms"))
  N <- if (!is.null(draws$data$N_trait)) draws$data$N_trait
       else if (!is.null(draws$data$N_tg)) draws$data$N_tg
       else if (is.cor_fixed(draws$autocor)) 1
       else draws$data$N
  out <- do.call(cbind, 
    lapply(seq_len(N), predict_fun, draws = draws, ntrys = ntrys)
  )
  # percentage of invalid samples for truncated discrete models
  # should always be zero for all other models; see predict.R
  pct_invalid <- get_pct_invalid(out, lb = draws$data$lb, ub = draws$data$ub) 
  if (pct_invalid >= 0.01) {
    warning2(round(pct_invalid * 100), "% of all predicted values ", 
             "were invalid. Increasing argument 'ntrys' may help.")
  }

  # reorder predicted responses in case of multivariate models
  # as they are sorted after units first not after traits
  if (grepl("_mv$", draws$f$family)) {
    nresp <- draws$data$nresp
    reorder <- ulapply(seq_len(nresp), seq, to = N * nresp, by = nresp)
    out <- out[, reorder, drop = FALSE]
  }
  old_order <- attr(draws$data, "old_order")
  out <- reorder_obs(out, old_order, sort = sort)
  colnames(out) <- NULL
  # transform predicted response samples before summarizing them 
  is_catordinal <- is_ordinal(draws$f) || is_categorical(draws$f)
  if (!is.null(transform) && !is_catordinal) {
    out <- do.call(transform, list(out))
  }
  if (summary) {
    if (is_catordinal) {
      # compute frequencies of categories 
      out <- get_table(out, levels = seq_len(max(draws$data$ncat)))
    } else {
      out <- get_summary(out, probs = probs, robust = robust)
    }
    rownames(out) <- seq_len(nrow(out))
  }
  out
}

#' @rdname predict.brmsfit
#' @aliases posterior_predict
#' @method posterior_predict brmsfit
#' @export
#' @export posterior_predict
#' @importFrom rstantools posterior_predict
posterior_predict.brmsfit <- function(object, newdata = NULL, re_formula = NULL,
                                      transform = NULL, allow_new_levels = FALSE,
                                      sample_new_levels = "uncertainty", 
                                      incl_autocor = TRUE, subset = NULL, 
                                      nsamples = NULL, sort = FALSE,
                                      ntrys = 5, robust = FALSE,
                                      probs = c(0.025, 0.975), ...) {
  cl <- match.call()
  cl[[1]] <- quote(predict)
  cl[["summary"]] <- quote(FALSE)
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
#' @param auxpar Optional name of a predicted auxiliary parameter.
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
#'   
#' @details \code{NA} values within factors in \code{newdata}, 
#'   are interpreted as if all dummy variables of this factor are 
#'   zero. This allows, for instance, to make predictions of the grand mean 
#'   when using sum coding.
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
                           incl_autocor = TRUE, auxpar = NULL, subset = NULL, 
                           nsamples = NULL, sort = FALSE, summary = TRUE, 
                           robust = FALSE, probs = c(0.025, 0.975), ...) {
  scale <- match.arg(scale)
  contains_samples(object)
  object <- restructure(object)
  draws_args <- nlist(
    x = object, newdata, re_formula, incl_autocor, 
    allow_new_levels, sample_new_levels, subset, nsamples
  )
  draws <- do.call(extract_draws, draws_args)
  auxpars <- intersect(auxpars(), names(draws))
  if (!length(auxpar)) {
    if (is.list(draws$mu[["mv"]])) {
      draws$mu <- get_eta(draws$mu)
    }
    for (ap in auxpars) {
      if (is.list(draws[[ap]])) {
        draws[[ap]] <- get_auxpar(draws[[ap]])
      }
    }
    if (grepl("_mv$", draws$f$family) && length(dim(draws$mu)) == 3L) {
      # collapse over responses in linear MV models
      dim(draws$mu) <- c(dim(draws$mu)[1], prod(dim(draws$mu)[2:3]))
    }
    if (scale == "response") {
      # original families are required for fitted helper functions 
      draws$f <- family(object)
      fitted_fun <- paste0("fitted_", draws$f$family)
      fitted_fun <- get(fitted_fun, asNamespace("brms"))
      draws$mu <- fitted_fun(draws)
    }
  } else {
    auxpars <- setdiff(auxpars, "mu")
    if (length(auxpar) != 1L || !auxpar %in% auxpars) {
      stop2("Invalid argument 'auxpar'. Valid auxiliary ",
            "parameters are: ", collapse_comma(auxpars))
    }
    if (!is.list(draws[[auxpar]])) {
      stop2("Auxiliary parameter '", auxpar, "' was not predicted.")
    }
    if (scale == "linear") {
      draws[[auxpar]]$f$link <- "identity"
    }
    draws$mu <- get_auxpar(draws[[auxpar]])
  }
  old_order <- attr(draws$data, "old_order")
  draws$mu <- reorder_obs(draws$mu, old_order, sort = sort)
  colnames(draws$mu) <- NULL
  if (summary) {
    draws$mu <- get_summary(draws$mu, probs = probs, robust = robust)
    rownames(draws$mu) <- seq_len(nrow(draws$mu))
  }
  draws$mu
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
                              incl_autocor = TRUE, subset = NULL, 
                              nsamples = NULL, sort = FALSE,
                              summary = TRUE, robust = FALSE, 
                              probs = c(0.025, 0.975), ...) {
  type <- match.arg(type)
  method <- match.arg(method)
  contains_samples(object)
  object <- restructure(object)
  family <- family(object)
  if (is_ordinal(family) || is_categorical(family)) {
    stop2("Residuals not implemented for family '", family$family, "'.")
  }
  
  standata <- amend_newdata(
    newdata, fit = object, re_formula = re_formula,
    allow_new_levels = allow_new_levels, check_response = TRUE
  )
  if (!is.null(standata$cens)) {
    warning2("Residuals may not be meaningful for censored models.")
  }
  if (is.null(subset) && !is.null(nsamples)) {
    subset <- sample(nsamples(object), nsamples)
  }
  pred_args <- nlist(
    object, newdata, re_formula, allow_new_levels,
    sample_new_levels, incl_autocor, subset, sort, 
    summary = FALSE
  )
  mu <- do.call(method, pred_args)
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

#' @rdname residuals.brmsfit
#' @aliases predictive_error
#' @method predictive_error brmsfit
#' @export
#' @export predictive_error
#' @importFrom rstantools predictive_error
predictive_error.brmsfit <- function(object, newdata = NULL, re_formula = NULL, 
                                     type = c("ordinary", "pearson"),
                                     allow_new_levels = FALSE, 
                                     sample_new_levels = "uncertainty",
                                     incl_autocor = TRUE, subset = NULL, 
                                     nsamples = NULL, sort = FALSE,
                                     robust = FALSE, probs = c(0.025, 0.975),
                                     ...) {
  cl <- match.call()
  cl[[1]] <- quote(residuals)
  cl[c("method", "summary")] <- list(quote("predict"), quote(FALSE))
  eval(cl, parent.frame())
}

#' Update \pkg{brms} models
#' 
#' This method allows to update an existing \code{brmsfit} object
#' 
#' @param object An object of class \code{brmsfit}.
#' @param formula. Changes to the formula; for details see 
#'   \code{\link[stats:update.formula]{update.formula}} and
#'   \code{\link[brms:brmsformula]{brmsformula}}.
#' @param newdata Optional \code{data.frame} 
#'   to update the model with new data.
#' @param recompile Logical, indicating whether the Stan model should 
#'  be recompiled. If \code{FALSE} (the default), the model is only 
#'  recompiled when necessary.
#' @param ... Other arguments passed to \code{\link[brms:brm]{brm}}.
#'  
#' @details Sometimes, when updating the model formula, 
#'  it may happen that \R complains about a mismatch
#'  between \code{model frame} and \code{formula}.
#'  This error can be avoided by supplying your orginal data
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
                           recompile = FALSE, ...) {
  dots <- list(...)
  if ("data" %in% names(dots)) {
    # otherwise the data name cannot be found by substitute 
    stop2("Please use argument 'newdata' to update the data.")
  }
  object <- restructure(object)
  if (isTRUE(object$version$brms < utils::packageVersion("brms"))) {
    warning2("Updating models fitted with older versions of brms may fail.")
  }
  if (missing(formula.)) {
    dots$formula <- object$formula
  } else {
    family <- get_arg("family", dots, formula., object)
    nl <- get_arg("nl", formula., formula(object))
    dots$formula <- bf(formula., family = family, nl = nl)
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
  
  arg_names <- c("prior", "autocor", "nonlinear", "threshold", 
                 "cov_ranef", "sparse", "sample_prior")
  new_args <- intersect(arg_names, names(dots))
  old_args <- setdiff(arg_names, new_args)
  dots[old_args] <- object[old_args]
  if (!is.null(newdata)) {
    dots$data <- newdata
  } else  {
    dots$data <- rm_attr(object$data, c("terms", "brmsframe"))
  }
  if (is.null(dots$threshold)) {
    # for backwards compatibility with brms <= 0.8.0
    if (grepl("(k - 1.0) * delta", object$model, fixed = TRUE)) {
      dots$threshold <- "equidistant"
    } else {
      dots$threshold <- "flexible"
    }
  }
  if ("prior" %in% new_args) {
    if (!is.brmsprior(dots$prior)) { 
      stop2("Invalid 'prior' argument.")
    }
    dots$prior <- rbind(dots$prior, object$prior)
    dots$prior <- dots$prior[!duplicated(dots$prior[, 2:5]), ]
  }
  pnames <- parnames(object)
  if (is.null(dots$sample_prior)) {
    dots$sample_prior <- any(grepl("^prior_", pnames))
  }
  if (is.null(dots$save_ranef)) {
    dots$save_ranef <- any(grepl("^r_", pnames)) || !nrow(object$ranef)
  }
  if (is.null(dots$save_mevars)) {
    dots$save_mevars <- any(grepl("^Xme_", pnames))
  }
  if (is.null(dots$sparse)) {
    dots$sparse <- grepl("sparse matrix", stancode(object))
  }
  dots$iter <- first_not_null(dots$iter, object$fit@sim$iter)
  # brm computes warmup automatically based on iter 
  dots$chains <- first_not_null(dots$chains, object$fit@sim$chains)
  dots$thin <- first_not_null(dots$thin, object$fit@sim$thin)
  
  new_stancode <- suppressMessages(
    do.call(make_stancode, dots[!names(dots) %in% "testmode"])
  )
  # only recompile if new and old stan code do not match
  if (recompile || !identical(new_stancode, stancode(object))) {
    # recompliation is necessary
    message("The desired updates require recompling the model")
    dots$fit <- NA
    if (!is.null(newdata)) {
      dots$data.name <- Reduce(paste, deparse(substitute(newdata)))
      dots$data.name <- substr(dots$data.name, 1, 50)
    } else {
      dots$data.name <- object$data.name
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
    bterms <- parse_bf(object$formula, family = object$family)
    object$data <- update_data(
      dots$data, bterms = bterms, family = object$family
    )
    if (!is.null(newdata)) {
      object$data.name <- Reduce(paste, deparse(substitute(newdata)))
      object$ranef <- tidy_ranef(bterms, data = object$data)
      dots$is_newdata <- TRUE
    }
    if (!is.null(dots$sample_prior)) {
      prior_only <- identical(dots$sample_prior, "only")
      attr(object$prior, "prior_only") <- prior_only
    }
    if (!is.null(dots$save_ranef) || !is.null(dots$save_mevars)) {
      if (is.null(dots$save_ranef)) {
        dots$save_ranef <- any(grepl("^r_", pnames)) || !nrow(object$ranef)
      }
      if (is.null(dots$save_mevars)) {
        dots$save_mevars <- any(grepl("^Xme_", pnames))
      }
      object$exclude <- exclude_pars(
        bterms, data = object$data, ranef = object$ranef, 
        save_ranef = dots$save_ranef, save_mevars = dots$save_mevars
      )
    }
    if (!is.null(dots$algorithm)) {
      aopts <- c("sampling", "meanfield", "fullrank")
      algorithm <- match.arg(dots$algorithm, aopts)
      dots$algorithm <- object$algorithm <- algorithm
    } else if (!is.null(object$algorithm)) {
      dots$algorithm <- object$algorithm
    }
    if (!isTRUE(dots$testmode)) {
      object <- do.call(brm, c(list(fit = object), dots))
    }
  }
  object
}

#' @export
#' @describeIn WAIC \code{WAIC} method for \code{brmsfit} objects
WAIC.brmsfit <- function(x, ..., compare = TRUE, newdata = NULL, 
                         re_formula = NULL, allow_new_levels = FALSE,
                         sample_new_levels = "uncertainty", subset = NULL,
                         nsamples = NULL, pointwise = NULL) {
  models <- list(x, ...)
  mnames <- deparse(substitute(x))
  mnames <- c(mnames, sapply(substitute(list(...))[-1], deparse))
  if (is.null(subset) && !is.null(nsamples)) {
    subset <- sample(nsamples(x), nsamples)
  }
  if (is.null(pointwise)) {
    pointwise <- set_pointwise(x, subset = subset, newdata = newdata)
  }
  ll_args = nlist(
    newdata, re_formula, allow_new_levels, 
    sample_new_levels, subset, pointwise
  )
  args <- nlist(ic = "waic", ll_args)
  if (length(models) > 1L) {
    out <- named_list(mnames)
    for (i in seq_along(models)) {
      args[["x"]] <- models[[i]]
      out[[i]] <- do.call(compute_ic, args)
      out[[i]]$model_name <- mnames[i]
    }
    if (compare) {
      match_response(models)
      out <- compare_ic(x = out)
    }
    class(out) <- c("iclist", "list")
  } else {
    out <- do.call(compute_ic, c(nlist(x), args))
    out$model_name <- mnames
  }
  out
}

#' @importFrom loo waic
#' @export waic
#' @export
waic.brmsfit <- function(x, ..., compare = TRUE, newdata = NULL,
                         re_formula = NULL, allow_new_levels = FALSE,
                         sample_new_levels = "uncertainty", subset = NULL, 
                         nsamples = NULL, pointwise = NULL) {
  cl <- match.call()
  cl[[1]] <- quote(WAIC)
  eval(cl, parent.frame())
}

#' @export
#' @describeIn LOO \code{LOO} method for \code{brmsfit} objects
LOO.brmsfit <- function(x, ..., compare = TRUE, newdata = NULL, 
                        re_formula = NULL, allow_new_levels = FALSE, 
                        sample_new_levels = "uncertainty", subset = NULL, 
                        nsamples = NULL, pointwise = NULL,
                        cores = 1, wcp = 0.2, wtrunc = 3/4) {
  models <- list(x, ...)
  mnames <- deparse(substitute(x))
  mnames <- c(mnames, sapply(substitute(list(...))[-1], deparse))
  if (is.null(subset) && !is.null(nsamples)) {
    subset <- sample(nsamples(x), nsamples)
  }
  if (is.null(pointwise)) {
    pointwise <- set_pointwise(x, subset = subset, newdata = newdata)
  }
  ll_args = nlist(
    newdata, re_formula, allow_new_levels, 
    sample_new_levels, subset, pointwise
  )
  args <- nlist(ic = "loo", ll_args, wcp, wtrunc, cores)
  if (length(models) > 1L) {
    out <- named_list(mnames)
    for (i in seq_along(models)) {
      args[["x"]] <- models[[i]]
      out[[i]] <- do.call(compute_ic, args)
      out[[i]]$model_name <- mnames[i]
    }
    if (compare) {
      match_response(models)
      out <- compare_ic(x = out)
    }
    class(out) <- c("iclist", "list")
  } else {
    out <- do.call(compute_ic, c(nlist(x), args))
    out$model_name <- mnames
  }
  out
}

#' @importFrom loo loo
#' @export loo
#' @export
loo.brmsfit <- function(x, ..., compare = TRUE, newdata = NULL,
                        re_formula = NULL, allow_new_levels = FALSE,
                        sample_new_levels = "uncertainty", subset = NULL,
                        nsamples = NULL, pointwise = NULL,
                        cores = 1, wcp = 0.2, wtrunc = 3/4) {
  cl <- match.call()
  cl[[1]] <- quote(LOO)
  eval(cl, parent.frame())
}

#' Compute the Pointwise Log-Likelihood
#' 
#' @aliases log_lik logLik.brmsfit
#' 
#' @param object A fitted model object of class \code{brmsfit}. 
#' @inheritParams predict.brmsfit
#' @param pointwise A flag indicating whether to compute the full
#'   log-likelihood matrix at once (the default), or just return
#'   the likelihood function along with all data and samples
#'   required to compute the log-likelihood separately for each
#'   observation. The latter option is rarely useful when
#'   calling \code{log_lik} directly, but rather when computing
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
#' @aliases log_lik
#' @method log_lik brmsfit
#' @export
#' @export log_lik
#' @importFrom rstantools log_lik
log_lik.brmsfit <- function(object, newdata = NULL, re_formula = NULL,
                            allow_new_levels = FALSE, 
                            sample_new_levels = "uncertainty", subset = NULL,
                            nsamples = NULL, pointwise = FALSE, ...) {
  contains_samples(object)
  object <- restructure(object)
  draws_args <- nlist(
    x = object, newdata, re_formula, allow_new_levels,
    sample_new_levels, subset, nsamples, check_response = TRUE
  )
  draws <- do.call(extract_draws, draws_args)
  
  N <- if (!is.null(draws$data$N_tg)) draws$data$N_tg
       else if (is.cor_fixed(object$autocor)) 1
       else nrow(as.matrix(draws$data$Y))
  loglik_fun <- paste0("loglik_", draws$f$family)
  loglik_fun <- get(loglik_fun, asNamespace("brms"))
  if (pointwise) {
    loglik <- structure(loglik_fun, draws = draws, N = N)
  } else {
    if (is.list(draws$mu[["mv"]])) {
      draws$mu <- get_eta(draws$mu)
    }
    for (ap in intersect(auxpars(), names(draws))) {
      if (is.list(draws[[ap]])) {
        draws[[ap]] <- get_auxpar(draws[[ap]])
      }
    }
    loglik <- do.call(cbind, lapply(seq_len(N), loglik_fun, draws = draws))
    old_order <- attr(draws$data, "old_order")
    # do not loglik reorder for ARMA covariance models
    sort <- use_cov(object$autocor)
    loglik <- reorder_obs(loglik, old_order, sort = sort)
    colnames(loglik) <- NULL
  }
  loglik
}

#' @export
logLik.brmsfit <- function(object, newdata = NULL, re_formula = NULL,
                           allow_new_levels = FALSE,
                           sample_new_levels = "uncertainty", subset = NULL,
                           nsamples = NULL, pointwise = FALSE, ...) {
  cl <- match.call()
  cl[[1]] <- quote(log_lik)
  eval(cl, parent.frame())
}

#' @rdname hypothesis
#' @export
hypothesis.brmsfit <- function(x, hypothesis, class = "b", group = "",
                               alpha = 0.05, seed = 1234, ...) {
  set.seed(1234)
  if (!is.character(hypothesis)) {
    stop2("Argument 'hypothesis' must be a character vector.")
  }
  if (length(alpha) != 1L || alpha < 0 || alpha > 1) {
    stop2("Argument 'alpha' must be a single value in [0,1].")
  }
  contains_samples(x)
  x <- restructure(x)
  
  if (!length(class)) {
    class <- "" 
  }
  if (length(class) != 1L || length(group) != 1L) {
    stop2("Arguments 'class' and 'group' must be of length one.")
  }
  valid_classes <- c(
    "", "b", "bcs", "bmo", "bme", "bm", "sd", "cor", 
    "r", "sds", "s", "simplex", "sigma", "rescor"
  )
  if (!class %in% valid_classes) {
    stop2("'", class, "' is not a valid parameter class.")
  }
  if (class %in% c("sd", "cor", "r") && nzchar(group)) {
    class <- paste0(class, "_", group, "__")
  } else if (nzchar(class)) {
    class <- paste0(class, "_")
  }

  hyp_fun <- function(h) {
    # internal function to evaluate hypotheses
    # Args:
    #   h: A string containing a hypothesis
    h <- rename(h, c("[ \t\r\n]", ":"), c("", "___"), fixed = FALSE)
    sign <- unlist(regmatches(h, gregexpr("=|<|>", h)))
    lr <- unlist(regmatches(h, gregexpr("[^=<>]+", h)))
    if (length(sign) != 1 || length(lr) != 2) {
      stop2("Every hypothesis must be of the form 'left (= OR < OR >) right'.")
    }
    h <- paste0("(", lr[1], ")")
    h <- paste0(h, ifelse(lr[2] != "0", paste0("-(", lr[2], ")"), ""))
    varsH <- unique(find_names(h))
    parsH <- paste0(class, varsH)
    missing_pars <- setdiff(parsH, pars)
    if (length(missing_pars)) {
      stop2("The following parameters cannot be found in the model: \n", 
            collapse_comma(gsub("___", ":", missing_pars)))
    }
    # prepare for renaming of parameters so that h can be evaluated
    parsH <- rename(parsH, "___", ":")
    h_renamed <- rename(h, c("[", "]", ","), c(".", ".", ".."))
    symbols <- c(paste0("^", class), ":", "\\[", "\\]", ",")
    subs <- c("", "___", ".", ".", "..")
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
    } else {
      prior_samples <- NULL
    }
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
    rownames(sm) <- paste(rename(h, "___", ":"), sign, "0")
    cl <- (1 - alpha) * 100
    colnames(sm) <- c("Estimate", "Est.Error", paste0("l-", cl, "% CI"),
                      paste0("u-", cl, "% CI"), "Evid.Ratio", "Star")
    if (!is.null(prior_samples)) {
      samples <- c(samples, prior_samples)
    } else {
      samples <- c(samples, rep(NA, nrow(samples)))
    }
    nlist(summary = sm, samples = samples)
  }

  pars <- rename(parnames(x)[grepl("^", class, parnames(x))],
                 symbols = ":", subs = "___")
  hlist <- lapply(hypothesis, hyp_fun)
  hs <- do.call(rbind, lapply(hlist, function(h) h$summary))
  samples <- do.call(cbind, lapply(hlist, function(h) h$samples))
  samples <- as.data.frame(samples) 
  names(samples) <- paste0("H", seq_along(hlist))
  samples$Type <- factor(rep(c("posterior", "prior"), each = nsamples(x)))
  class <- sub("_+$", "", class)
  out <- nlist(hypothesis = hs, samples, class, alpha)
  class(out) <- "brmshypothesis"
  out
}

#' @rdname expose_functions
#' @export
expose_functions.brmsfit <- function(x, ...) {
  expose_stan_functions(x$fit)
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
