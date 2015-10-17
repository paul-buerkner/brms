#' @export
parnames.brmsfit <- function(x, ...) {
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  dimnames(x$fit)$parameters
}

#' @export
fixef.brmsfit <-  function(x, estimate = "mean", ...) {
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  pars <- parnames(x)
  fpars <- pars[grepl("^b_", pars)]
  if (!length(fpars)) 
    stop(paste("The model does not contain fixed effects")) 
  out <- posterior_samples(x, pars = fpars, exact_match = TRUE)
  out <- do.call(cbind, lapply(estimate, get_estimate, samples = out, ...))
  rownames(out) <- gsub("^b_", "", fpars)
  out
}

#' Covariance and Correlation Matrix of Fixed Effects
#' 
#' Get a point estimate of the covariance or correlation matrix of fixed effects parameters
#' 
#' @param object An object of class \code{brmsfit}
#' @param correlation logical; if \code{FALSE} (the default), compute the covariance matrix,
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
    stop(paste("The model does not contain fixed effects")) 
  samples <- posterior_samples(object, pars = fpars, exact_match = TRUE)
  names(samples) <- sub("^b_", "", names(samples))
  if (correlation) {
    cor(samples) 
  } else {
    cov(samples)
  }
}

#' @export
ranef.brmsfit <- function(x, estimate = "mean", var = FALSE, ...) {
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  if (!estimate %in% c("mean","median"))
    stop("Argument estimate must be either 'mean' or 'median'")
  if (!length(x$ranef))
    stop("The model does not contain random effects")
  group <- names(x$ranef)
  pars <- parnames(x)
  
  get_ranef <- function(i) {
    # get random effects of a grouping factor
    #
    # Args:
    #   i: index of a grouping factor
    g <- group[i]
    rnames <- x$ranef[[i]]
    rpars <- pars[grepl(paste0("^r_",group[i],"\\["), pars)]
    if (!length(rpars))
      stop(paste0("The model does not contain random effects for group '",g,"'\n",
                  "You should use argument ranef = TRUE in function brm."))
    rdims <- x$fit@sim$dims_oi[[paste0("r_",group[i])]]
    levels <- attr(x$ranef[[i]], "levels")
    if (is.null(levels)) {
      # avoid error in dimnames if levels are NULL 
      # for backwards compatibility with brms < 0.5.0 
      levels <- 1:rdims[1]
    }
    rs <- posterior_samples(x, pars = rpars, exact_match = TRUE)
    ncol <- ifelse(is.na(rdims[2]), 1, rdims[2])
    rs_array <- array(dim = c(rdims[1], ncol, nrow(rs)))
    k <- 0
    for (j in 1:ncol) {
      for (i in 1:rdims[1]) {
        k <- k + 1
        rs_array[i, j, ] <- rs[, k]
      }
    }
    out <- get_estimate(estimate, samples = rs_array, margin = 1:2, ...)
    colnames(out) <- rnames
    if(var) {
      Var <- array(dim = c(rep(ncol, 2), rdims[1]), 
                   dimnames = list(rnames, rnames, 1:rdims[1]))
      for (i in 1:rdims[1]) {
        if (is.na(rdims[2])) Var[, , i] <- var(rs_array[i, 1, ]) 
        else Var[, , i] <- cov(t(rs_array[i, , ])) 
      }
      dimnames(Var)[[3]] <- levels
      attr(out, "var") <- Var
    }
    rownames(out) <- levels
    out
  }
  
  ranef <- lapply(1:length(group), get_ranef)
  names(ranef) <- group
  ranef 
} 

#' @import abind
#' @export
VarCorr.brmsfit <- function(x, estimate = "mean", as.list = TRUE, ...) {
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  if (!(length(x$ranef) || (indicate_linear(x$family) && 
                            any(grepl("^sigma_", parnames(x))))))
    stop("The model does not contain covariance matrices")

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
  
  ee <- extract_effects(x$formula, family = x$family)
  if (length(x$ranef)) {
    gather_names <- function(i) {
      # gather names of random effects parameters
      cor_type <- paste0("cor_",group[i])
      list(rnames = x$ranef[[i]],
           type = paste0("cor_",group[i]),
           sd_pars = paste0("sd_",group[i],"_",x$ranef[[i]]),
           cor_pars = get_cornames(x$ranef[[i]], type = cor_type, 
                                   brackets = FALSE))
    }
    group <- names(x$ranef)
    p <- lapply(1:length(group), gather_names)
  } else {
    p <- group <- NULL
  } 
  
  # special treatment of residuals variances in linear models
  if (indicate_linear(x$family) && !is.formula(ee$se)) {
    cor_pars <- get_cornames(ee$response, type = "rescor", brackets = FALSE)
    p <- lc(p, list(rnames = ee$response, sd_pars = paste0("sigma_",ee$response),
                    cor_pars = cor_pars))
    group <- c(group, "RESIDUAL")
  } 
  VarCorr <- lapply(p, extract)
  names(VarCorr) <- group
  class(VarCorr) <- "VarCorr_brmsfit"
  VarCorr
}

#' @export
posterior_samples.brmsfit <- function(x, pars = NA, parameters = NA,  
                                      exact_match = FALSE, add_chains = FALSE, 
                                      as.matrix = FALSE, ...) {
  if (is.na(pars[1])) 
    pars <- parameters  
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  all_pars <- parnames(x)
  if (!(anyNA(pars) || is.character(pars))) 
    stop("Argument pars must be NA or a character vector")
  if (!anyNA(pars)) {
    if (exact_match) {
      pars <- all_pars[all_pars %in% pars]
    } else {
      pars <- all_pars[apply(sapply(pars, grepl, x = all_pars, ...), 1, any)]
    }
  }
  else {
    pars <- all_pars
  }
  
  # get basic information on the samples 
  iter <- attr(x$fit@sim$samples[[1]],"args")$iter
  warmup <- attr(x$fit@sim$samples[[1]],"args")$warmup
  thin <- attr(x$fit@sim$samples[[1]],"args")$thin
  chains <- length(x$fit@sim$samples) 
  final_iter <- (iter - warmup) / thin
  samples_taken <- seq((warmup + 1), iter, thin)
  
  if (length(pars)) {
    samples <- as.data.frame(x$fit, pars = pars)
    if (add_chains) {
      samples$chains <- factor(rep(1:chains, each = final_iter))
      samples$iter <- rep(samples_taken, chains)
    }
    if (as.matrix) {
      samples <- as.matrix(samples)
    }
  } else {
    samples <- NULL 
  }
  samples
}

#' @export
prior_samples.brmsfit <- function(x, pars = NA, parameters = NA, ...) {
  if (is.na(pars[1])) 
    pars <- parameters 
  if (!anyNA(pars) && !is.character(pars)) 
    stop("pars must be a character vector")
  par_names <- parnames(x)
  prior_names <- par_names[grepl("^prior_", par_names)]
  if (length(prior_names)) {
    samples <- posterior_samples(x, pars = prior_names, exact_match = TRUE)
    names(samples) <- sub("^prior_", "", prior_names)
    if (!anyNA(pars)) {
      get_samples <- function(par) {
        # get prior samples for parameter par 
        matches <- lapply(paste0("^",sub("^prior_", "", prior_names)), 
                          regexpr, text = par)
        matches <- ulapply(matches, attr, which = "match.length")
        if (max(matches) == -1) {
          return(NULL)
        } else {
          take <- match(max(matches), matches)
          return(structure(list(samples[, take]), names = par))
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
#'   defaults to 2. The effective sample size is always rounded to integers.
#' @param ... Additional arguments that would be passed to method \code{summary} of \code{brmsfit}.
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
#'   (this will take some time for larger models)
#' @param ... Other potential arguments
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @method summary brmsfit
#' @export
summary.brmsfit <- function(object, waic = TRUE, ...) {
  ee <- extract_effects(object$formula, family = object$family)
  formula <- update_formula(object$formula, partial = object$partial)
  out <- brmssummary(formula = formula,
                     family = object$family, 
                     link = object$link, 
                     data.name = object$data.name, 
                     group = names(object$ranef), 
                     nobs = nobs(object), 
                     ngrps = brms::ngrps(object), 
                     autocor = object$autocor)
  if (length(object$fit@sim)) {
    out$n.chains <- length(object$fit@sim$samples)
    out$n.iter <- attr(object$fit@sim$samples[[1]],"args")$iter
    out$n.warmup <- attr(object$fit@sim$samples[[1]],"args")$warmup
    out$n.thin <- attr(object$fit@sim$samples[[1]],"args")$thin
    out$sampler <- attr(object$fit@sim$samples[[1]],"args")$sampler_t
    if (waic) out$WAIC <- WAIC(object)$waic
    
    pars <- parnames(object)
    meta_pars <- object$fit@sim$pars_oi
    meta_pars <- meta_pars[!apply(sapply(paste0("^",c("r_","prior_")), 
                                  grepl, x = meta_pars, ...), 1, any)]
    fit_summary <- rstan::summary(object$fit, pars = meta_pars,
                                  probs = c(0.025, 0.975))
    col_names <- c("Estimate", "Est.Error", "l-95% CI", 
                   "u-95% CI", "Eff.Sample", "Rhat")
    
    fix_pars <- pars[grepl("^b_", pars)]
    out$fixed <- matrix(fit_summary$summary[fix_pars,-c(2)], ncol = 6)
    colnames(out$fixed) <- col_names
    rownames(out$fixed) <- gsub("^b_", "", fix_pars)
    
    spec_pars <- pars[pars %in% c("nu","shape","delta") | 
      apply(sapply(c("^sigma_", "^rescor_"), grepl, x = pars), 1, any)]
    out$spec_pars <- matrix(fit_summary$summary[spec_pars,-c(2)], ncol = 6)
    if (indicate_linear(object$family)) {
      spec_pars[grepl("^sigma_", spec_pars)] <- paste0("sigma(",ee$response,")")
      spec_pars[grepl("^rescor_", spec_pars)] <- get_cornames(ee$response, type = "rescor")   
    }    
    colnames(out$spec_pars) <- col_names
    rownames(out$spec_pars) <- spec_pars
    
    cor_pars <- pars[grepl("^ar|^ma", pars)]
    out$cor_pars <- matrix(fit_summary$summary[cor_pars,-c(2)], ncol = 6)
    colnames(out$cor_pars) <- col_names
    rownames(out$cor_pars) <- cor_pars
    
    if (length(out$group)) {
      for (i in 1:length(out$group)) {
        rnames <- object$ranef[[i]]
        sd_pars <- paste0("sd_", out$group[i],"_",rnames)
        cor_pars <- intersect(get_cornames(rnames, type = paste0("cor_",out$group[i]),
                                           brackets = FALSE), parnames(object))
        sd_names <- paste0("sd(",rnames,")")
        cor_names <- get_cornames(rnames, subset = cor_pars, subtype = out$group[i])
        out$random[[out$group[i]]] <- 
          matrix(fit_summary$summary[c(sd_pars, cor_pars), -c(2)], ncol = 6)
        colnames(out$random[[out$group[i]]]) <- col_names
        rownames(out$random[[out$group[i]]]) <- c(sd_names, cor_names)
      }
    }
  }  
  out
}

#' @export
nobs.brmsfit <- function(object, ...) {
  length(standata(object)$Y)
}
  
#' @export
ngrps.brmsfit <- function(object, ...) {
  standata <- standata(object)
  group <- extract_effects(object$formula, family = object$family)$group
  if (length(group)) {
    out <- setNames(lapply(1:length(group), function(i) 
      standata[[paste0("N_",i)]]), group)
    out <- out[!duplicated(group)]
  } else out <- NULL
  out
}

#' @export
formula.brmsfit <- function(x, ...) 
  x$formula

#' @export
family.brmsfit <- function(object, ...) 
  family(object$family, link = object$link)

#' @export
stancode.brmsfit <- function(object, ...)
  object$model

#' @export
standata.brmsfit <- function(object, ...) {
  if (is.data.frame(object$data)) {
    # brms > 0.5.0 stores the original model.frame 
    standata <- brmdata(object$formula, data = object$data, 
                        family = object$family, autocor = object$autocor, 
                        cov.ranef = object$cov.ranef, partial = object$partial,
                        ...)
  } else {
    # brms <= 0.5.0 only stores the data passed to Stan 
    standata <- object$data
  }
  standata
}
  
#' @export
launch_shiny.brmsfit <- function(x, rstudio = getOption("shinystan.rstudio"), ...) {
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
#' @param ask logical; Indicates if the user is prompted before a new page is plotted.   
#' @param ... Further arguments passed to \code{\link[gridExtra:arrangeGrob]{arrangeGrob}}.
#' 
#' @return NULL
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' \dontrun{ 
#' fit_e <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1|patient) + (1|visit), 
#'              data = epilepsy, family = "poisson")
#' ## plot fixed effects as well as standard devations of the random effects
#' plot(fit_e)
#' ## plot fixed effects only and combine the chains into one posterior
#' plot(fit_e, pars = "^b_", combine = TRUE) 
#' }
#' 
#' @method plot brmsfit
#' @import ggplot2
#' @importFrom gridExtra arrangeGrob
#' @importFrom grDevices devAskNewPage
#' @export
plot.brmsfit <- function(x, pars = NA, parameters = NA, N = 5, ask = TRUE, ...) {
  if (is.na(pars[1])) 
    pars <- parameters 
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  if (!is.wholenumber(N) || N < 1) 
    stop("N must be a positive integer")
  if (!is.character(pars)) 
    pars <- c("^b_", "^sd_", "^cor_", "^sigma", "^rescor", "^nu$", 
              "^shape$", "^delta$", "^ar", "^ma")
  samples <- posterior_samples(x, pars = pars, add_chains = TRUE)
  pars <- names(samples)[which(!names(samples) %in% c("chains", "iter"))] 
  
  default_ask <- devAskNewPage()
  on.exit(devAskNewPage(default_ask))
  devAskNewPage(ask = FALSE)
  for (i in 1:ceiling(length(pars) / N)) {
    plots <- lapply(pars[((i - 1) * N + 1):min(i * N, length(pars))], 
                    td_plot, x = samples)
    plot(arrangeGrob(grobs = unlist(plots, recursive = FALSE), 
                     nrow = length(plots), ncol = 2, ...))
    if (i == 1) devAskNewPage(ask = ask)
  }
}

#' @describeIn stanplot A fitted model object returned by function
#'                      \code{\link[brms:brm]{brm}} 
#' @export
stanplot.brmsfit <- function(object, pars = NA, type = "plot", 
                             exact_match = FALSE, quiet = FALSE, ...) {
  
  # check validity of type first
  basic_types <- c("plot", "trace", "scat", "hist", "dens", "ac")
  diag_types <- c("diag", "par", "rhat", "ess", "mcse")
  if (!type %in% c(basic_types, diag_types))
    stop(paste("Invalid plot type. Valid plot types are: \n",
               paste(c(basic_types, diag_types), collapse = ", ")))
  plot_fun <- get(paste0("stan_", type), mode = "function")
  
  # ensure that only desired parameters are plotted
  if (!is(object$fit, "stanfit") || !length(object$fit@sim)) 
    stop("The model does not contain posterior samples")
  all_pars <- parnames(object)
  if (!(anyNA(pars) || is.character(pars))) 
    stop("Argument pars must be NA or a character vector")
  if (!anyNA(pars)) {
    if (exact_match) {
      pars <- all_pars[all_pars %in% pars]
    } else {
      pars <- all_pars[apply(sapply(pars, grepl, x = all_pars), 1, any)]
    }
  } 
  
  dots <- list(...)
  args <- c(object = object$fit, dots)
  basic_types <- c("plot", "trace", "scat", "hist", "dens", "ac")
  diag_types <- c("diag", "par", "rhat", "ess", "mcse")
  if (type %in% basic_types) {
    if (!anyNA(pars)) {
      args <- c(args, list(pars = pars))
    }
  } else if (type %in% diag_types) {
    if (type %in% c("rhat", "ess", "msce") && !anyNA(pars)) {
      args <- c(args, list(pars = pars))
    } else if (type == "par") {
      args <- c(args, list(pars = pars))
    } 
  }
  # make the plot
  if (quiet) {
    quietgg(do.call(plot_fun, args))
  } else {
    do.call(plot_fun, args)
  }
}

#' Model Predictions of \code{brmsfit} Objects
#' 
#' Make predictions based on the fitted model parameters. 
#' Can be performed for the data used to fit the model (posterior predictive checks) or for new data.
#' 
#' @param object An object of class \code{brmsfit}
#' @param newdata An optional data.frame containing new data to make predictions for.
#'   If \code{NULL} (default), the data used to fit the model is applied.
#' @param re_formula formula containing random effects to be considered in the prediction. 
#'   If \code{NULL} (default), include all random effects; if \code{NA}, include no random effects.
#'   Other options will be implemented in the future.
#' @param transform A function or a character string naming a function to be applied on the predicted responses
#'   before summary statistics are computed.
#' @param allow_new_levels Currenly, \code{FALSE} (no new levels allowed) is the only option. 
#'   This will change in future versions of the package.
#' @param summary logical. Should summary statistics (i.e. means, sds, and 95\% intervals) be returned
#'  instead of the raw values. Default is \code{TRUE}
#' @param probs The percentiles to be computed by the \code{quantile} function. 
#'  Only used if \code{summary} is \code{TRUE}.
#' @param ntrys Parameter used in rejection sampling for truncated discrete models only 
#'   (defaults to \code{5}). See Details for more information.
#' @param ... Currently ignored
#' 
#' @return Predicted values of the response variable. If \code{summary = TRUE} the output depends on the family:
#'   For catagorical and ordinal families, it is a N x C matrix where N is the number of observations and
#'   C is the number of categories. For all other families, it is a N x E matrix where E is equal 
#'   to \code{length(probs) + 2}.
#'   If \code{summary = FALSE}, the output is as a S x N matrix, where S is the number of samples.
#' 
#' @details For truncated discrete models only:
#'   In the absence of any general algorithm to sample from truncated discrete distributions,
#'   rejection sampling is applied in this special case. This means that values are sampled until 
#'   a value lies within the defined truncation boundaries. In practice, this procedure may be 
#'   rather slow (especially in R). Thus, we try to do approximate rejection sampling 
#'   by sampling each value \code{ntrys} times and then select a valid value. 
#'   If all values are invalid, the closest boundary is used, instead. 
#'   If there are more than a few of these pathological cases, 
#'   a warning will occure suggesting to increase argument \code{ntrys}.
#'   
#'   For \pkg{brms} <= 0.5.0 only: 
#'   Be careful when using \code{newdata} with factors in fixed or random effects. 
#'   The predicted results are only valid if all factor levels present in the initial 
#'   data are also defined and ordered correctly for the factors in \code{newdata}.
#'   Grouping factors may contain fewer levels than in the inital data without causing problems.
#'   When using higher versions of \pkg{brms}, all factors are automatically checked 
#'   for correctness and amended if necessary.
#' 
#' @examples 
#' \dontrun{
#' ## fit a model
#' fit <- brm(time | cens(censored) ~ age + sex, data = kidney,
#'            family = "exponential", silent = TRUE)
#' 
#' ## posterior predictive checks
#' pp <- predict(fit)
#' head(pp)
#' 
#' ## predict response for new data
#' newdata <- data.frame(sex = factor(c("male", "female")),
#'                       age = c(20,50))
#' predict(fit, newdata = newdata)
#' }
#' 
#' @importFrom statmod rinvgauss pinvgauss qinvgauss
#' @export 
predict.brmsfit <- function(object, newdata = NULL, re_formula = NULL,
                            transform = NULL, allow_new_levels = FALSE,
                            summary = TRUE, probs = c(0.025, 0.975), 
                            ntrys = 5, ...) {
  if (!is(object$fit, "stanfit") || !length(object$fit@sim)) 
    stop("The model does not contain posterior samples")
  ee <- extract_effects(object$formula, family = object$family)
  # use newdata if defined
  if (is.null(newdata)) {
    data <- standata(object, keep_intercept = TRUE)
  } else {
    data <- amend_newdata(newdata, fit = object, re_formula = re_formula,
                          allow_new_levels = allow_new_levels)
  }
  
  # compute all necessary samples
  nresp <- length(ee$response)
  samples <- list(eta = linear_predictor(object, newdata = data, 
                                         re_formula = re_formula))
  if (indicate_linear(object$family) && !is.formula(ee$se))
    samples$sigma <- as.matrix(posterior_samples(object, pars = "^sigma_"))
  if (object$family == "student") 
    samples$nu <- as.matrix(posterior_samples(object, pars = "^nu$"))
  if (indicate_shape(object$family)) 
    samples$shape <- as.matrix(posterior_samples(object, pars = "^shape$"))
  if (object$family == "gaussian" && nresp > 1) {
    samples$rescor <- as.matrix(posterior_samples(object, pars = "^rescor_"))
    samples$Sigma <- get_cov_matrix(sd = samples$sigma, cor = samples$rescor)$cov
    message(paste("Computing posterior predictive samples of multinormal distribution. \n", 
                  "This may take a while."))
  }
  
  # call predict functions
  family <- object$family
  if (object$link == "log" && family == "gaussian" && nresp == 1) {
    family <- "lognormal"
  } else if (family == "gaussian" && nresp > 1) {
    family <- "multinormal"
  }
  is_catordinal <- indicate_ordinal(family) || family == "categorical"
  # see predict.R
  predict_fun <- get(paste0("predict_", family), mode = "function")
  call_predict_fun <- function(n) {
    do.call(predict_fun, list(n = n, data = data, samples = samples, 
                              link = object$link, ntrys = ntrys))
  }
  N <- ifelse(is.null(data$N_trait), data$N, data$N_trait)
  out <- do.call(cbind, lapply(1:N, call_predict_fun))
  
  # percentage of invalid samples for truncated discrete models
  # should always be zero for all other models
  pct_invalid <- get_pct_invalid(out, data = data)  # see predict.R
  if (pct_invalid >= 0.01) {
    warning(paste0(round(pct_invalid * 100), "% of all predicted values ", 
                   "were invalid. Increasing argument ntrys may help."))
  }
  
  if (!is.null(transform) && !is_catordinal) {
    out <- do.call(transform, list(out))
  }
  if (summary && !is_catordinal) {
    out <- get_summary(out, probs = probs)
  } else if (summary && is_catordinal) { 
    # compute frequencies of categories for categorical and ordinal models
    out <- get_table(out, levels = 1:max(data$max_obs)) 
  }
  
  # sort predicted responses in case of multinormal models
  if (family == "multinormal") {
    to_order <- ulapply(1:data$K_trait, seq, to = data$N, by = data$K_trait)
    if (summary) {
      # observations in rows
      out <- out[to_order, ]
      rownames(out) <- 1:nrow(out)
    } else {
      # observations in columns
      out <- out[, to_order]  
      colnames(out) <- 1:ncol(out) 
    }
  }
  out
}

#' Extract Model Fitted Values of \code{brmsfit} Objects
#' 
#' @inheritParams predict.brmsfit
#' @param scale Either \code{"response"} or \code{"linear"}. If \code{scale = "response"} 
#'   results are returned on the scale of the response variable. If \code{scale = "linear"} 
#'   fitted values are returned on the scale of the linear predictor.
#'
#' @return Fitted values extracted from \code{object}. The output depends on the family:
#'   If \code{summary = TRUE} it is a N x E x C array for categorical and ordinal models 
#'   and a N x E matrix else.
#'   If \code{summary = FALSE} it is a S x N x C array for categorical and ordinal models 
#'   and a S x N matrix else.
#'   N is the number of observations, S is the number of samples, C is the number of categories,
#'   and E is equal to \code{length(probs) + 2}.
#'   
#' @details For \pkg{brms} <= 0.5.0 only: 
#'  Be careful when using \code{newdata} with factors in fixed or random effects. 
#'  The predicted results are only valid if all factor levels present in the initial 
#'  data are also defined and ordered correctly for the factors in \code{newdata}.
#'  Grouping factors may contain fewer levels than in the inital data without causing problems.
#'  When using higher versions of \pkg{brms}, all factors are automatically checked 
#'  for correctness and amended if necessary.
#'
#' @examples 
#' \dontrun{
#' ## fit a model
#' fit <- brm(rating ~ treat + period + carry + (1|subject), data = inhaler,
#'            n.cluster = 2)
#' 
#' ## extract fitted values
#' fitted_values <- fitted(fit)
#' head(fitted_values)
#' 
#' ## plot fitted means against actual response
#' dat <- as.data.frame(cbind(Y = fit$data$Y, fitted_values))
#' ggplot(dat) + geom_point(aes(x = Estimate, y = Y))
#' }
#' 
#' @export 
fitted.brmsfit <- function(object, newdata = NULL, re_formula = NULL,
                           scale = c("response", "linear"),
                           allow_new_levels = FALSE,
                           summary = TRUE, probs = c(0.025, 0.975), ...) {
  scale <- match.arg(scale)
  if (!is(object$fit, "stanfit") || !length(object$fit@sim)) 
    stop("The model does not contain posterior samples")
  ee <- extract_effects(object$formula, family = object$family)
  # use newdata if defined
  if (is.null(newdata)) {
    data <- standata(object, keep_intercept = TRUE)
  } else {
    data <- amend_newdata(newdata, fit = object, re_formula = re_formula,
                          allow_new_levels = allow_new_levels)
  }
  # get mu and scale it appropriately
  mu <- linear_predictor(object, newdata = data, re_formula = re_formula)
  if (scale == "response") {
    # see fitted.R
    mu <- fitted_response(object, eta = mu, data = data)
  }
  if (summary) {
    mu <- get_summary(mu, probs = probs)
  }
  mu
}

#' Extract Model Residuals from brmsfit Objects
#' 
#' @inheritParams predict.brmsfit
#' @param type The type of the residuals, either \code{"ordinary"} or \code{"pearson"}. 
#'   More information is provided under 'Details'. 
#' 
#' @details Residuals of type \code{ordinary} are of the form \eqn{R = Y - Yp},
#'   where \eqn{Y} is the observed and \eqn{Yp} is the predicted response.
#'   Residuals of type \code{pearson} are of the form \eqn{R = (Y - Yp) / Var(Y)},
#'   where \eqn{Var(Y)} is an estimation of the variance of \eqn{Y}. \cr
#'   
#'   Currently, \code{residuals.brmsfit} does not support \code{categorical} or ordinal models. 
#' 
#' @return Model residuals. If \code{summary = TRUE} this is a N x C matrix 
#'   and if \code{summary = FALSE} a S x N matrix, where S is the number of samples, 
#'   N is the number of observations, and C is equal to \code{length(probs) + 2}.  
#' 
#' @examples 
#' \dontrun{
#' ## fit a model
#' fit <- brm(rating ~ treat + period + carry + (1|subject), data = inhaler,
#'            n.cluster = 2)
#' 
#' ## extract residuals 
#' res <- residuals(fit, summary = TRUE)
#' head(res)
#' }
#' 
#' @export
residuals.brmsfit <- function(object, re_formula = NULL, type = c("ordinary", "pearson"), 
                              summary = TRUE, probs = c(0.025, 0.975), ...) {
  type <- match.arg(type)
  if (!is(object$fit, "stanfit") || !length(object$fit@sim)) 
    stop("The model does not contain posterior samples")
  if (indicate_ordinal(object$family) || object$family == "categorical")
    stop(paste("residuals not yet implemented for family", object$family))
  
  standata <- standata(object)
  mu <- fitted(object, re_formula = re_formula, summary = FALSE)
  Y <- matrix(rep(as.numeric(standata$Y), nrow(mu)), 
              nrow = nrow(mu), byrow = TRUE)
  res <- Y - mu
  colnames(res) <- NULL
  if (type == "pearson") {
    # get predicted standard deviation for each observation
    sd <- matrix(rep(predict(object, re_formula = re_formula, 
                             summary = TRUE)[, 2], nrow(mu)), 
                 nrow = nrow(mu), byrow = TRUE)
    res <- res / sd
  }
  # for compatibility with the macf function (see correlations.R)
  # so that the colnames of the output correspond to the levels of the autocor grouping factor
  if (has_arma(object$autocor)) {
    tgroup <- extract_time(object$autocor$formula)$group
    if (nchar(tgroup)) colnames(res) <- standata[[tgroup]]
  }
  if (summary) {
    res <- get_summary(res, probs = probs)
  }
  res
}

#' Update \pkg{brms} models
#' 
#' This method allows to update an existing \code{brmsfit} object with new data
#' as well as changed configuration of the chains.
#' 
#' @param object object of class \code{brmsfit}
#' @param newdata optional \code{data.frame} to update the model with new data
#' @param ... other arguments passed to \code{\link[brms:brm]{brm}} such as
#'    \code{n.iter}, \code{n.warmup}, \code{n.chains}, \code{n.thin}, 
#'    \code{n.cluster}, and \code{inits}.
#'
#' @export
update.brmsfit <- function(object, newdata = NULL, ...) {
  dots <- list(...)
  z <- which(names(dots) %in% c("formula", "family", "prior", "autocor", 
                                "partial", "threshold", "cov.ranef", 
                                "sample.prior"))
  if (length(z)) {
    stop(paste("Argument(s)", paste(names(dots)[z], collapse = ", "),
               "cannot be updated"))
  }
  if (!is.null(newdata)) {
    object$data <- newdata
  }
  do.call(brm, c(list(fit = object), dots))
}

#' @export
WAIC.brmsfit <- function(x, ..., compare = TRUE) {
  models <- list(x, ...)
  names <- c(deparse(substitute(x)), sapply(substitute(list(...))[-1], deparse))
  if (length(models) > 1) {
    out <- setNames(lapply(models, calculate_ic, ic = "waic"), names)
    class(out) <- c("iclist", "list")
    if (compare) 
      attr(out, "compare") <- compare_ic(out, ic = "waic")
  } else { 
    out <- calculate_ic(x, ic = "waic")
  }
  out
}

#' @export
#' @describeIn LOO method for class \code{brmsfit}
LOO.brmsfit <- function(x, ..., compare = TRUE,
                        cores = getOption("loo.cores", parallel::detectCores()),
                        wcp = 0.2, wtrunc = 3/4) {
  models <- list(x, ...)
  names <- c(deparse(substitute(x)), sapply(substitute(list(...))[-1], deparse))
  if (length(models) > 1) {
    out <- setNames(lapply(models, calculate_ic, ic = "loo", wcp = wcp, 
                           wtrunc = wtrunc, cores = cores), names)
    class(out) <- c("iclist", "list")
    if (compare) 
      attr(out, "compare") <- compare_ic(out, ic = "loo")
  } else {
    out <- calculate_ic(x, ic = "loo", wcp = wcp, wtrunc = wtrunc, cores = cores)
  }
  out
}

#' Compute the pointwise log-likelihood
#' 
#' @param object A fitted model object of class \code{brmsfit}. 
#' @param ... Currently ignored
#' 
#' @return Usually, an S x N matrix containing the pointwise log-likelihood samples, 
#'         where S is the number of samples and N is the number of observations in the data. 
#' 
#' @importFrom statmod dinvgauss
#' @export
logLik.brmsfit <- function(object, ...) {
  if (!is(object$fit, "stanfit") || !length(object$fit@sim)) 
    stop("The model does not contain posterior samples")
  ee <- extract_effects(object$formula, family = object$family)
  nresp <- length(ee$response)
  
  # extract relevant samples
  samples <- list(eta = linear_predictor(object))
  if (indicate_linear(object$family) && !is.formula(ee$se))
    samples$sigma <- as.matrix(posterior_samples(object, pars = "^sigma_"))
  if (object$family == "student") 
    samples$nu <- as.matrix(posterior_samples(object, pars = "^nu$"))
  if (indicate_shape(object$family)) 
    samples$shape <- as.matrix(posterior_samples(object, pars = "^shape$"))
  if (object$family == "gaussian" && nresp > 1) {
    samples$rescor <- as.matrix(posterior_samples(object, pars = "^rescor_"))
    samples$Sigma <- get_cov_matrix(sd = samples$sigma, cor = samples$rescor)$cov
    message(paste("Computing pointwise log-likelihood of multinormal distribution. \n",
                  "This may take a while."))
  }
  
  # prepare for calling family specific loglik functions
  family <- object$family
  if (object$link == "log" && family == "gaussian" && nresp == 1) {
    family <- "lognormal"
  } else if (family == "gaussian" && nresp > 1) {
    family <- "multinormal"
  }
  standata <- standata(object)
  loglik_fun <- get(paste0("loglik_", family), mode = "function")
  call_loglik_fun <- function(n) {
    do.call(loglik_fun, list(n = n, data = standata, samples = samples, 
                             link = object$link)) 
  }
  N <- nrow(as.matrix(standata$Y))
  loglik <- do.call(cbind, lapply(1:N, call_loglik_fun))
  colnames(loglik) <- 1:ncol(loglik)
  loglik
}

#' @export
hypothesis.brmsfit <- function(x, hypothesis, class = "b", group = "",
                               alpha = 0.05, ...) {
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  if (!is.character(hypothesis)) 
    stop("Argument hypothesis must be a character vector")
  if (alpha < 0 || alpha > 1)
    stop("Argument alpha must be in [0,1]")
  
  # process class and group arguments
  if (is.null(class)) class <- ""
  valid_classes <- c("", "b", "sd", "cor", "ar", "ma", "sigma", 
                     "rescor", "nu", "shape", "delta")
  if (!class %in% valid_classes)
    stop(paste(class, "is not a valid paramter class"))
  if (class %in% c("b", "sd", "cor", "sigma", "rescor")) {
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
      stop("Every hypothesis must be of the form 'left (= OR < OR >) right'")
    h <- paste0(lr[1], ifelse(lr[2] != "0", paste0("-(",lr[2],")"), ""))
    varsH <- unique(find_names(h))
    parsH <- paste0(class, varsH)
    if (!all(parsH %in% pars)) 
      stop(paste("The following parameters cannot be found in the model:", 
                 paste0(gsub("__", ":", parsH[which(!parsH %in% pars)]), 
                        collapse = ", ")))
    
    # prepare for renaming of parameters so that h can be evaluated
    parsH <- rename(parsH, "__", ":")
    h_renamed <- rename(h, c("[", "]"), c(".", "."))
    symbols <- c(paste0("^",class), ":", "\\[", "\\]")
    subs <- c("", "__", ".", ".")
    
    # get posterior samples
    samples <- posterior_samples(x, pars = parsH, exact_match = TRUE)
    names(samples) <- rename(names(samples), symbols = symbols, 
                             subs = subs, fixed = FALSE)
    samples <- matrix(with(samples, eval(parse(text = h_renamed))), ncol=1)
    
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
    out <- lapply(c("mean", "sd", "quantile", "evidence_ratio"), 
                  get_estimate, samples = samples, probs = probs, 
                  wsign = wsign, prior_samples = prior_samples, ...)
    out <- as.data.frame(matrix(unlist(out), nrow = 1))
    if (sign == "<") {
      out[1, 3] <- -Inf
    } else if (sign == ">") {
      out[1, 4] <- Inf
    }
    out <- cbind(out, ifelse(!(out[1, 3] <= 0 && 0 <= out[1, 4]), '*', ''))
    rownames(out) <- paste(rename(h, "__", ":"), sign, "0")
    cl <- (1 - alpha) * 100
    colnames(out) <- c("Estimate", "Est.Error", paste0("l-",cl,"% CI"), 
                       paste0("u-",cl,"% CI"), "Evid.Ratio", "")
    out
  }
  
  pars <- rename(parnames(x)[grepl("^", class, parnames(x))],
                 symbols = ":", subs = "__")
  out <- do.call(rbind, lapply(hypothesis, hyp_fun))
  out <- list(hypothesis = out, class = substr(class, 1, nchar(class) - 1), 
              alpha = alpha)
  class(out) <- "brmshypothesis"
  out
}