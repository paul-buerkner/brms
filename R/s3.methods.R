#' @export
fixef.brmsfit <-  function(x, estimate = "mean", ...) {
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  pars <- parnames(x)
  fpars <- pars[grepl("^b_", pars)]
  if (!length(fpars)) stop(paste("No fixed effect present in argument x")) 
  out <- posterior_samples(x, parameters = paste0("^",fpars,"$"))
  out <- do.call(cbind, lapply(estimate, get_estimate, samples = out, ...))
  rownames(out) <- gsub("^b_", "", fpars)
  out
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
  
  ranef <- lapply(group, function(g) {
    rpars <- pars[grepl(paste0("^r_",g,"\\["), pars)]
    rnames <- x$ranef[[match(g, names(x$ranef))]]
    if (!length(rpars))
      stop(paste0("The model does not contain random effects for group '",g,"'\n",
                  "You should use argument ranef = TRUE in function brm."))
    rdims <- x$fit@par_dims[[paste0("r_",gsub(":", "__", g))]]
    rs <- posterior_samples(x, parameters = rpars, fixed = TRUE)
    ncol <- ifelse(is.na(rdims[2]), 1, rdims[2])
    rs_array <- array(dim = c(rdims[1], ncol, nrow(rs)))
    k <- 0
    for (j in 1:ncol) {
      for (i in 1:rdims[1]) {
        k <- k + 1
        rs_array[i,j,] <- rs[,k]
      }
    }
    out <- get_estimate(estimate, samples = rs_array, margin = 1:2, ...)
    colnames(out) <- rnames
    if(var) {
      Var <- array(dim = c(rep(ncol, 2), rdims[1]), 
                   dimnames = list(rnames, rnames, 1:rdims[1]))
      for (i in 1:rdims[1]) {
        if (is.na(rdims[2])) Var[,,i] <- var(rs_array[i,1,]) 
        else Var[,,i] <- cov(t(rs_array[i,,])) 
      }
      attr(out, "var") <- Var
    }
    out
  })
  names(ranef) <- group
  ranef 
} 

#' @export
VarCorr.brmsfit <- function(x, estimate = "mean", as.list = TRUE, ...) {
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  if (!(length(x$ranef) || x$family %in% c("gaussian", "student", "cauchy")))
    stop("The model does not contain covariance matrices")

  # extracts samples for sd, cor and cov
  extract <- function(p) {
    nr <- length(p$sd_pars)
    sd <- posterior_samples(x, parameters = paste0("^",p$sd_pars,"$"))
    nsamples <- nrow(sd)
    out <- list(sd = do.call(cbind, lapply(estimate, get_estimate, samples = sd, ...)))
    rownames(out$sd) <- p$rnames 
    
    # calculate correlation and covariance matrices
    if (length(p$cor_pars))
      cor <- posterior_samples(x, parameters = paste0("^",p$cor_pars,"$"))
    else cor <- NULL
    matrices <- cov_matrix(sd = sd, cor = cor)
    out$cor <- list2array(lapply(estimate, get_estimate, samples = matrices$cor, 
                            margin = c(2,3), to.array = TRUE, ...))
    out$cov <- list2array(lapply(estimate, get_estimate, samples = matrices$cov, 
                            margin = c(2,3), to.array = TRUE, ...))
    if (length(p$rnames) > 1)
      dimnames(out$cov) <- dimnames(out$cor) <- list(p$rnames, p$rnames, dimnames(out$cor)[[3]])
    else dimnames(out$cov) <- dimnames(out$cor) <- list(dimnames(out$cor)[[1]])
    if (as.list) {
      out$cor <- lapply(array2list(out$cor), function(x)
        if (is.null(dim(x))) structure(matrix(x), dimnames = list(p$rnames, p$rnames)) else x)
      out$cov <- lapply(array2list(out$cov), function(x)
        if (is.null(dim(x))) structure(matrix(x), dimnames = list(p$rnames, p$rnames)) else x)
    }
    out
  }
  
  ee <- extract_effects(x$formula, family = x$family)
  if (length(x$ranef)) {
    group <- names(x$ranef)
    p <- lapply(1:length(group), function(i)
    list(rnames = x$ranef[[i]],
         sd_pars = paste0("sd_",group[i],"_",x$ranef[[i]]),
         cor_pars = get_cornames(x$ranef[[i]], type = paste0("cor_",group[i]),
                                  eval = length(x$ranef[[i]]) > 1 && ee$cor[[i]], 
                                  brackets = FALSE)))
  } else p <- group <- NULL
  if (x$family %in% c("gaussian", "student", "cauchy") && !is.formula(ee$se)) {
    p[[length(p)+1]] <- list(rnames = ee$response, 
                             sd_pars = paste0("sigma_",ee$response),
                             cor_pars = get_cornames(ee$response, type = "rescor", 
                                                      eval = length(ee$response) > 1, 
                                                      brackets = FALSE))
    group <- c(group, "RESIDUAL")
  } 
  VarCorr <- lapply(p, extract)
  names(VarCorr) <- group
  VarCorr
}

#' @export
posterior_samples.brmsfit <- function(x, parameters = NA, add.chains = FALSE, ...) {
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  pars <- parnames(x)
  if (!(anyNA(parameters) || is.character(parameters))) 
    stop("Argument parameters must be NA or a character vector")
  if (!anyNA(parameters)) pars <- pars[apply(sapply(parameters, grepl, x = pars, ...), 1, any)]
  
  iter <- attr(x$fit@sim$samples[[1]],"args")$iter
  warmup <- attr(x$fit@sim$samples[[1]],"args")$warmup
  thin <- attr(x$fit@sim$samples[[1]],"args")$thin
  chains <- length(x$fit@sim$samples) 
  
  if (length(pars)) {
    samples <- data.frame(sapply(1:length(pars), function(i)
      unlist(lapply(1:chains, function(j) 
        x$fit@sim$samples[[j]][[pars[i]]][(warmup/thin+1):(iter/thin)]))))
    names(samples) <- pars
    if (add.chains) {
      samples$chains <- factor(do.call(c, lapply(1:chains, rep, times = nrow(samples)/chains)))
      samples$iter <- rep((warmup+1):(nrow(samples)/chains+warmup), chains)
    }
  } else samples <- NULL  
  samples
}

#' @export
prior_samples.brmsfit <- function(x, parameters = NA, ...) {
  if (!anyNA(parameters) && !is.character(parameters)) 
    stop("pars must be a character vector")
  par_names <- parnames(x)
  prior_names <- par_names[grepl("^prior_", par_names)]
  if (length(prior_names)) {
    samples <- posterior_samples(x, parameters = prior_names, fixed = TRUE)
    names(samples) <- sub("^prior_", "", prior_names)
    if (!anyNA(parameters)) {
      samples <- data.frame(rmNULL(lapply(parameters, function(par) {
        matches <- lapply(paste0("^",sub("^prior_", "", prior_names)), regexpr, 
                          text = par)
        matches <- unlist(lapply(matches, attr, which = "match.length"))
        if (max(matches) == -1) NULL
        else structure(list(samples[,match(max(matches), matches)]), names = par)
      })))
    }
  }
  else samples <- NULL
  samples
}

#' Create a summary of a fitted model represented by a \code{brmsfit} object
#' 
#' Summarize estimated fixed and random effects as well as other useful
#' results included in a \code{brmsfit} object.
#'
#' @param object An object of class \code{brmsfit}
#' @param ... Other potential arguments
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @method summary brmsfit
#' @export
summary.brmsfit <- function(object, ...) {
  ee <- extract_effects(object$formula, family = object$family)
  out <- brmssummary(formula = update_formula(object$formula, partial = object$partial),
             family = object$family, link = object$link, data.name = object$data.name, 
             group = names(object$ranef), nobs = nobs(object), ngrps = brms::ngrps(object), 
             autocor = object$autocor)
  if (length(object$fit@sim)) {
    out$n.chains <- length(object$fit@sim$samples)
    out$n.iter = attr(object$fit@sim$samples[[1]],"args")$iter
    out$n.warmup = attr(object$fit@sim$samples[[1]],"args")$warmup
    out$n.thin = attr(object$fit@sim$samples[[1]],"args")$thin
    out$sampler = attr(object$fit@sim$samples[[1]],"args")$sampler_t
    if (length(ee$response) == 1) out$WAIC <- WAIC(object)$waic
    
    pars <- parnames(object)
    meta_pars <- object$fit@sim$pars_oi
    meta_pars <- meta_pars[!apply(sapply(paste0("^",c("r_","prior_")), 
                                  grepl, x = meta_pars, ...), 1, any)]
    fit_summary <- rstan::summary(object$fit, pars = meta_pars, probs = c(0.025, 0.975))
    col_names <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Eff.Sample", "Rhat")
    
    fix_pars <- pars[grepl("^b_", pars)]
    out$fixed <- matrix(fit_summary$summary[fix_pars,-c(2)], ncol = 6)
    colnames(out$fixed) <- col_names
    rownames(out$fixed) <- gsub("^b_", "", fix_pars)
    
    spec_pars <- pars[pars %in% c("nu","shape","delta") | 
      apply(sapply(c("^sigma_", "^rescor_"), grepl, x = pars), 1, any)]
    out$spec_pars <- matrix(fit_summary$summary[spec_pars,-c(2)], ncol = 6)
    if (object$family %in% c("gaussian", "student", "cauchy")) {
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
        has_cor <- length(rnames) > 1 && ee$cor[[i]]
        sd_pars <- paste0("sd_", out$group[i],"_",rnames)
        cor_pars <- get_cornames(rnames, type = paste0("cor_",out$group[i]),
                                  eval = has_cor, brackets = FALSE)
        sd_names <- paste0("sd(",rnames,")")
        cor_names <- get_cornames(rnames, eval = has_cor)
        out$random[[out$group[i]]] <- matrix(fit_summary$summary[c(sd_pars, cor_pars),-c(2)], ncol = 6)
        colnames(out$random[[out$group[i]]]) <- col_names
        rownames(out$random[[out$group[i]]]) <- c(sd_names, cor_names)
      }
    }
  }  
  out
}

#' @export
nobs.brmsfit <- function(object, ...) length(object$data$Y)

#' @export
ngrps.brmsfit <- function(object, ...) {
  group <- names(object$ranef)
  setNames(lapply(group, function(g) object$data[[paste0("N_",gsub(":","__",g))]]), group)
}

#' @export
formula.brmsfit <- function(x, ...) x$formula

#' Extract Model Fitted Values of \code{brmsfit} Objects
#' 
#' @inheritParams residuals.brmsfit
#' 
#' @details Currently, the method does not support \code{categorical} or ordinal models. 
#'
#' @return Fitted values extracted from \code{object}.
#'
#' @examples 
#' \dontrun{
#' ## fit a model
#' fit <- brm(rating ~ treat + period + carry + (1|subject), data = inhaler,
#'            n.cluster = 2)
#' 
#' ## extract fitted values
#' fitted_values <- fitted(fit, summary = TRUE)
#' head(fitted_values)
#' }
#' 
#' @export 
fitted.brmsfit <- function(object, summary = FALSE, ...) {
  if (!is(object$fit, "stanfit") || !length(object$fit@sim)) 
    stop("The model does not contain posterior samples")
  if (object$family %in% c("categorical", "cumulative", "sratio", "cratio", "acat"))
    stop(paste("fitted values not yet implemented for family", object$family))
  mu <- ilink(linear_predictor(object), object$link)
  if (summary) {
    mu <- do.call(cbind, lapply(c("mean", "sd", "quantile"), get_estimate, 
                                 samples = mu, probs = c(0.025, 0.975)))
    colnames(mu) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI")
  }
  mu
}

#' Extract Model Residuals from brmsfit Objects
#' 
#' @param object An object of class \code{brmsfit}
#' @param type The type of the residuals, either \code{ordinary} or \code{pearson}. Currently, the latter 
#'   is only supported for families \code{binomial} and \code{bernoulli}
#' @param summary logical. Should summary statistics (i.e. means, sds, and 95\% intervals) be returned
#'  instead of the raw values. Default is \code{FALSE}
#' @param ... Currently ignored
#' 
#' @details Currently, the method does not support \code{categorical} or ordinal models. 
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
residuals.brmsfit <- function(object, type = c("ordinary", "pearson"), summary = FALSE, ...) {
  type <- match.arg(type)
  if (!is(object$fit, "stanfit") || !length(object$fit@sim)) 
    stop("The model does not contain posterior samples")
  if (object$family %in% c("categorical", "cumulative", "sratio", "cratio", "acat"))
    stop(paste("residuals not yet implemented for family", object$family))
  if (type == "pearson" && !family %in% c("binomial", "bernoulli"))
    stop("Pearson residuals are currently only available for families binomial and bernoulli")
  
  mu <- fitted(object)
  if (object$family %in% c("binomial", "bernoulli")) {
    if (object$family == "binomial")
      max_obs <- matrix(rep(object$data$max_obs, nrow(mu)), nrow = nrow(mu), byrow = TRUE)
    else max_obs <- 1
    mu <- mu * max_obs
  }
  Y <- matrix(rep(as.numeric(object$data$Y), nrow(mu)), nrow = nrow(mu), byrow = TRUE)
  res <- Y - mu
  colnames(res) <- NULL
  
  if (type == "pearson")
    res <- res / sqrt(mu * (max_obs - mu) / max_obs)  
  
  # for compatibility with the macf function
  if (is(object$autocor, "cor_arma") && sum(object$autocor$p, object$autocor$q) > 0) {
    tgroup <- extract_time(object$autocor$formula)$group
    if (nchar(tgroup)) colnames(res) <- object$data[[tgroup]]
  }
  
  if (summary) {
    res <- do.call(cbind, lapply(c("mean", "sd", "quantile"), get_estimate, 
                                 samples = res, probs = c(0.025, 0.975)))
    colnames(res) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI")
  }
  res
}

#' @export
parnames.brmsfit <- function(x, ...) {
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  dimnames(x$fit)$parameters
}

#' @export
hypothesis.brmsfit <- function(x, hypothesis, class = "b", alpha = 0.05, ...) {
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  if (!is.character(hypothesis)) 
    stop("Argument hypothesis must be a character vector")
  if (alpha < 0 || alpha > 1)
    stop("Argument alpha must be in [0,1]")
  if (!is.null(class) && nchar(class)) class <- paste0(class,"_")
  else class <- ""
  pars <- rename(parnames(x)[grepl("^",class, parnames(x))],
                 symbols = ":", subs = "__")
  
  out <- do.call(rbind, lapply(hypothesis, function(h) {
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
                 paste0(gsub("__", ":", parsH[which(!parsH %in% pars)]), collapse = ", ")))
    
    #get posterior samples
    samples <- posterior_samples(x, parameters = rename(parsH, "__", ":"), fixed = TRUE)
    names(samples) <- rename(names(samples), symbols = c(paste0("^",class), ":"), 
                             subs = c("", "__"), fixed = FALSE)
    samples <- matrix(with(samples, eval(parse(text = rename(h, c("[", "]"), c("OB", "CB"))))), ncol=1)
    
    #get prior samples
    prior_samples <- prior_samples(x, parameters = rename(parsH, "__", ":"), fixed = TRUE)
    if (!is.null(prior_samples) && ncol(prior_samples) == length(varsH)) {
      names(prior_samples) <- rename(names(prior_samples), symbols = c(paste0("^",class), ":"), 
                                     subs = c("", "__"), fixed = FALSE)
      prior_samples <- matrix(with(prior_samples, eval(parse(text = rename(h, c("[", "]"), c("OB", "CB"))))), ncol=1)
    } else prior_samples <- NULL

    #evaluate hypothesis
    wsign <- ifelse(sign == "=", "equal", ifelse(sign == "<", "less", "greater"))
    probs <- switch(wsign, equal = c(alpha/2, 1-alpha/2), less = c(0, 1-alpha), greater = c(alpha, 1))
    out <- as.data.frame(matrix(unlist(lapply(c("mean", "sd", "quantile", "eratio"), get_estimate, 
      samples = samples, probs = probs, wsign = wsign, prior_samples = prior_samples, ...)), nrow = 1))
    if (sign == "<") out[1,3] <- -Inf
    else if (sign == ">") out[1,4] <- Inf
    out <- cbind(out, ifelse(!(out[1,3] <= 0 && 0 <= out[1,4]), '*', ''))
    rownames(out) <- paste(rename(h, "__", ":"), sign, "0")
    cl <- (1-alpha)*100
    colnames(out) <- c("Estimate", "Est.Error", paste0("l-",cl,"% CI"), paste0("u-",cl,"% CI"), "Evid.Ratio", "")
    out
  }))
  out <- list(hypothesis = out, class = substr(class, 1, nchar(class)-1), alpha = alpha)
  class(out) <- "brmshypothesis"
  out
}

#' @export
launch_shiny.brmsfit <- function(x, rstudio = getOption("shinystan.rstudio"), ...) {
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  shinystan::launch_shinystan(x$fit, rstudio = rstudio, ...)
}