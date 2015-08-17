#' @export
fixef.brmsfit <-  function(x, estimate = "mean", ...) {
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  pars <- par.names(x)
  f.pars <- pars[grepl("^b_", pars)]
  if (!length(f.pars)) stop(paste("No fixed effect present in argument x")) 
  out <- posterior.samples(x, parameters = paste0("^",f.pars,"$"))
  out <- do.call(cbind, lapply(estimate, get.estimate, samples = out, ...))
  rownames(out) <- gsub("^b_", "", f.pars)
  out
}

#' @export
ranef.brmsfit <- function(x, estimate = "mean", var = FALSE, ...) {
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  if (!estimate %in% c("mean","median"))
    stop("Argument estimate must be either 'mean' or 'median'")
  pars <- par.names(x)
  group <- names(x$ranef)
  
  ranef <- lapply(group, function(g) {
    r.pars <- pars[grepl(paste0("^r_",g,"\\["), pars)]
    r.names <- gsub(paste0("^sd_",g,"_"), "", pars[grepl(paste0("^sd_",g,"_"), pars)])
    if (!length(r.pars))
      stop(paste0("The model does not contain random effects for group '",g,"'\n",
                  "You should use argument ranef = TRUE in function brm."))
    r_dims <- x$fit@par_dims[[paste0("r_",gsub(":", "__", g))]]
    rs <- posterior.samples(x, parameters = r.pars, fixed = TRUE)
    n.samples <- nrow(rs)
    n.col <- ifelse(is.na(r_dims[2]), 1, r_dims[2])
    rs.array <- array(dim = c(r_dims[1], n.col, n.samples))
    k <- 0
    for (j in 1:n.col) {
      for (i in 1:r_dims[1]) {
        k <- k + 1
        rs.array[i,j,] <- rs[,k]
      }
    }
    out <- get.estimate(estimate, samples = rs.array, margin = 1:2, ...)
    colnames(out) <- r.names
    if(var) {
      Var <- array(dim = c(rep(n.col, 2), r_dims[1]), 
                   dimnames = list(r.names, r.names, 1:r_dims[1]))
      for (i in 1:r_dims[1]) {
        if (is.na(r_dims[2])) Var[,,i] <- var(rs.array[i,1,]) 
        else Var[,,i] <- cov(t(rs.array[i,,])) 
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
  pars <- par.names(x)
  group <- names(x$ranef)
  
  # extracts samples for sd, cor and cov
  extract <- function(pattern) {
    if (length(pattern) != 2) stop("pattern must be of length 2")
    sd.pars <- pars[grepl(pattern[1], pars)]
    cor.pars <- pars[grepl(pattern[2], pars)]
    r.names <- gsub(pattern[1], "", sd.pars)
    nr <- length(sd.pars)
    sds <- posterior.samples(x, parameters = paste0("^",sd.pars,"$"))
    n.samples <- nrow(sds)
    out <- list(sd = do.call(cbind, lapply(estimate, get.estimate, samples = sds, ...)))
    rownames(out$sd) <- r.names 
    if (length(cor.pars)) {
      cors <- posterior.samples(x, parameters = paste0("^",cor.pars,"$"))
      out$cor <- array(diag(1,nr), dim = c(nr, nr, n.samples))
      out$cov <- out$cor
      k <- 0 
      for (i in 1:nr) {
        for (j in 1:i) {
          if (i == j) out$cov[i,j,] <- sds[,i]^2
          else {
            k = k + 1
            out$cor[i,j,] <- cors[,k]
            out$cor[j,i,] <- out$cor[i,j,]
            out$cov[i,j,] <- out$cor[i,j,] * sds[,i] * sds[,j]
            out$cov[j,i,] <- out$cov[i,j,]
          }}}
      out$cor <- abind(lapply(estimate, get.estimate, samples = out$cor, 
                              margin=  c(1,2), to.array=TRUE, ...))
      out$cov <- abind(lapply(estimate, get.estimate, samples = out$cov, 
                              margin = c(1,2), to.array=TRUE, ...))
      dimnames(out$cor) <- list(r.names, r.names, dimnames(out$cor)[[3]])
      dimnames(out$cov) <- dimnames(out$cor)
      if (as.list) {
        out$cor <- array2list(out$cor)
        out$cov <- array2list(out$cov)
      }
    }
    out
  }
  
  pattern <- lapply(group, function(g) c(paste0("^sd_",g,"_"), paste0("^cor_",g,"_")))
  if (x$family %in% c("gaussian", "student", "cauchy")) {
    pattern <- c(pattern, list(c("^sigma_", "^rescor_")))
    group <- c(group, "RESIDUAL")
  } 
  VarCorr <- lapply(pattern, extract)
  names(VarCorr) <- group
  VarCorr
}

#' @export
posterior.samples.brmsfit <- function(x, parameters = NA, add.chains = FALSE, ...) {
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  pars <- par.names(x)
  if (!(anyNA(parameters) || is.character(parameters))) 
    stop("Argument parameters must be NA or a character vector")
  if (!anyNA(parameters)) pars <- pars[apply(sapply(parameters, grepl, x = pars, ...), 1, any)]
  
  iter <- attr(x$fit@sim$samples[[1]],"args")$iter
  warmup <- attr(x$fit@sim$samples[[1]],"args")$warmup
  thin <- attr(x$fit@sim$samples[[1]],"args")$thin
  chains <- length(x$fit@sim$samples) 
  
  samples <- data.frame(sapply(1:length(pars), function(i)
    unlist(lapply(1:chains, function(j) 
      x$fit@sim$samples[[j]][[pars[i]]][(warmup/thin+1):(iter/thin)]))))
  names(samples) <- pars
  if (add.chains) {
    samples$chains <- factor(do.call(c, lapply(1:chains, rep, times = nrow(samples)/chains)))
    samples$iter <- rep((warmup+1):(nrow(samples)/chains+warmup), chains)
  }  
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
  ee <- extract.effects(object$formula, add.ignore = TRUE)
  out <- brmssummary(formula = brm.update.formula(object$formula, partial = object$partial),
             family = object$family, link = object$link, data.name = object$data.name, 
             group = names(object$ranef), nobs = nobs(object), ngrps = brms::ngrps(object), 
             autocor = object$autocor)
  if (length(object$fit@sim)) {
    out$n.chains <- length(object$fit@sim$samples)
    out$n.iter = attr(object$fit@sim$samples[[1]],"args")$iter
    out$n.warmup = attr(object$fit@sim$samples[[1]],"args")$warmup
    out$n.thin = attr(object$fit@sim$samples[[1]],"args")$thin
    out$sampler = attr(object$fit@sim$samples[[1]],"args")$sampler_t
    if ("log_llh" %in% object$fit@model_pars) out$WAIC <- WAIC(object)$waic
    
    pars <- par.names(object)
    meta_pars <- object$fit@sim$pars_oi
    meta_pars <- meta_pars[!apply(sapply(paste0("^",c("r_","log_llh","prior_","Y_pred")), 
                                  grepl, x = meta_pars, ...), 1, any)]
    fit.summary <- rstan::summary(object$fit, pars = meta_pars, probs = c(0.025, 0.975))
    col.names <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Eff.Sample", "Rhat")
    
    fix.pars <- pars[grepl("^b_", pars)]
    out$fixed <- matrix(fit.summary$summary[fix.pars,-c(2)], ncol = 6)
    colnames(out$fixed) <- col.names
    rownames(out$fixed) <- gsub("^b_", "", fix.pars)
    
    spec.pars <- pars[pars %in% c("nu","shape","delta") | 
      apply(sapply(c("^sigma_", "^rescor_"), grepl, x = pars), 1, any)]
    out$spec.pars <- matrix(fit.summary$summary[spec.pars,-c(2)], ncol = 6)
    if (object$family %in% c("gaussian", "student", "cauchy")) {
      spec.pars[grepl("^sigma_", spec.pars)] <- paste0("sigma(",ee$response,")")
      spec.pars[grepl("^rescor_", spec.pars)] <- get.cor.names(ee$response, type = "rescor")   
    }    
    colnames(out$spec.pars) <- col.names
    rownames(out$spec.pars) <- spec.pars
    
    cor.pars <- pars[grepl("^ar|^ma", pars)]
    out$cor.pars <- matrix(fit.summary$summary[cor.pars,-c(2)], ncol = 6)
    colnames(out$cor.pars) <- col.names
    rownames(out$cor.pars) <- cor.pars
    
    if (length(out$group)) {
      for (i in 1:length(out$group)) {
        sd.pars <- pars[grepl(paste0("^sd_", out$group[i],"_"), pars)]
        cor.pars <- pars[grepl(paste0("^cor_", out$group[i],"_"), pars)]
        r.names <- gsub(paste0("^sd_",out$group[i],"_"), "", sd.pars)
        sd.names <- paste0("sd(",r.names,")")
        cor.names <- get.cor.names(r.names, eval = length(cor.pars))
        out$random[[out$group[i]]] <- matrix(fit.summary$summary[c(sd.pars, cor.pars),-c(2)], ncol = 6)
        colnames(out$random[[out$group[i]]]) <- col.names
        rownames(out$random[[out$group[i]]]) <- c(sd.names,cor.names)
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

#' @export 
predict.brmsfit <- function(object, ...) {
  if (!"Y_pred" %in% object$fit@model_pars) 
    stop(paste0("The model does not contain predicted values. \n",
                "You should use argument predict = TRUE in function brm."))
  if (!is(object$fit, "stanfit") || !length(object$fit@sim)) 
    stop("The model does not contain posterior samples")
  else {
    ee <- extract.effects(object$formula, add.ignore = TRUE)
    pars <- par.names(object)
    fit.summary <- rstan::summary(object$fit, probs = c(0.025, 0.975))
    pred.pars <- pars[grepl("^Y_pred\\[", pars)]
    out <- fit.summary$summary[pred.pars,-c(2,6,7)]
    colnames(out) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI")
  } 
  out
}

#' @export
WAIC.brmsfit <- function(x, ..., compare = TRUE) {
  models <- list(x, ...)
  names <- c(deparse(substitute(x)), sapply(substitute(list(...))[-1], deparse))
  if (length(models) > 1) {
    out <- setNames(lapply(models, calculate_ic, ic = "waic"), names)
    class(out) <- c("iclist", "list")
    if (compare) attr(out, "compare") <- compare_ic(out, ic = "waic")
  }  
  else out <- calculate_ic(x, ic = "waic")
  out
}

#' @export
LOO.brmsfit <- function(x, ..., compare = TRUE) {
  models <- list(x, ...)
  names <- c(deparse(substitute(x)), sapply(substitute(list(...))[-1], deparse))
  if (length(models) > 1) {
    out <- setNames(lapply(models, calculate_ic, ic = "loo"), names)
    class(out) <- c("iclist", "list")
    if (compare) attr(out, "compare") <- compare_ic(out, ic = "loo")
  }  
  else out <- calculate_ic(x, ic = "loo")
  out
}

#' @export
par.names.brmsfit <- function(x, ...) {
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
  pars <- rename(par.names(x)[grepl("^",class, par.names(x))],
                 symbols = ":", subs = "__")
  
  out <- do.call(rbind, lapply(hypothesis, function(h) {
    h <- rename(h, symbols = c(" ", ":"), subs = c("", "__"))
    sign <- unlist(regmatches(h, gregexpr("=|<|>", h)))
    lr <- unlist(regmatches(h, gregexpr("[^=<>]+", h)))
    if (length(sign) != 1 || length(lr) != 2)
      stop("Every hypothesis must be of the form 'left (= OR < OR >) right'")
    h <- paste0(lr[1], ifelse(lr[2] != "0", paste0("-(",lr[2],")"), ""))
    varsH <- find.names(h)
    parsH <- paste0(class, varsH)
    if (!all(parsH %in% pars)) 
      stop(paste("The following parameters cannot be found in the model:", 
                 paste0(gsub("__", ":", parsH[which(!parsH %in% pars)]), collapse = ", ")))
    samples <- structure(posterior.samples(x, parameters = rename(parsH, "__", ":"), fixed = TRUE),
                         names = rename(varsH, c("[", "]"), c("OB", "CB")))
    samples <- matrix(with(samples, eval(parse(text = rename(h, c("[", "]"), c("OB", "CB"))))), ncol=1)
    prior_samples <- get_prior_samples(x, pars = rename(parsH, "__", ":"))
    if (!is.null(prior_samples)) {
      names(prior_samples) <- rename(varsH, c("[", "]"), c("OB", "CB"))
      prior_samples <- matrix(with(prior_samples, eval(parse(text = rename(h, c("[", "]"), c("OB", "CB"))))), ncol=1)
    } 

    #evaluate hypothesis
    wsign <- ifelse(sign == "=", "equal", ifelse(sign == "<", "less", "greater"))
    probs <- switch(wsign, equal = c(alpha/2, 1-alpha/2), less = c(0, 1-alpha), greater = c(alpha, 1))
    out <- as.data.frame(matrix(unlist(lapply(c("mean","sd","quantile", "eratio"), get.estimate, 
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