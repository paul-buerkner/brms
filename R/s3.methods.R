#' @export
fixef.brmsfit <-  function(x, estimate = "mean", ...) {
  if (!is(x$fit, "stanfit")) 
    stop("Argument x does not contain posterior samples")
  pars <- names(x$fit@sim$samples[[1]])
  iter <- attr(x$fit@sim$samples[[1]],"args")$iter
  warmup <- attr(x$fit@sim$samples[[1]],"args")$warmup
  thin <- attr(x$fit@sim$samples[[1]],"args")$thin
  chains <- length(x$fit@sim$samples) 
  
  f.pars <- pars[grepl("^b_", pars)]
  f.names <- gsub("^b_", "", f.pars)
  f.names <- gsub("__", ":", f.names)
  nf <- length(f.pars)
  if (!nf)
    stop(paste("No fixed effect present in argument x")) 
  out <- t(sapply(1:nf, function(i)
    unlist(lapply(1:chains, function(j) 
      x$fit@sim$samples[[j]][[f.pars[i]]][(warmup/thin+1):(iter/thin)]))))
  out <- do.call(cbind, lapply(estimate, get.estimate, samples = out, ...))
  rownames(out) <- f.names
  out
}

#' @export
ranef.brmsfit <- function(x, estimate = "mean", var = FALSE, center.zero = TRUE, ...) {
  if (!is(x$fit, "stanfit")) 
    stop("Argument x does not contain posterior samples")
  if (!estimate %in% c("mean","median"))
    stop("Argument estimate must be either 'mean' or 'median'")
  pars <- names(x$fit@sim$samples[[1]])
  iter <- attr(x$fit@sim$samples[[1]],"args")$iter
  warmup <- attr(x$fit@sim$samples[[1]],"args")$warmup
  thin <- attr(x$fit@sim$samples[[1]],"args")$thin
  chains <- length(x$fit@sim$samples)
  n.samples <- (iter-warmup)/thin*chains
  group <- unlist(extract.effects(x$formula, add.ignore = TRUE)$group)
  
  ranef <- lapply(group, function(g) {
    r.pars <- pars[grepl(paste0("^r_",g,"\\["), pars)]
    r.names <- gsub(paste0("^sd_",g,"_"), "", pars[grepl(paste0("^sd_",g,"_"), pars)])
    r.names <- gsub("__",":",r.names) 
    nr <- length(r.pars)
    if (!nr)
      stop(paste0("The model does not contain random effects for group '",g,"'\n",
                  "You should use argument ranef = TRUE in function brm."))
    r_dims <- x$fit@par_dims[[paste0("r_",g)]]
    rs <- t(sapply(1:nr, function(i)
      unlist(lapply(1:chains, function(j) 
        x$fit@sim$samples[[j]][[r.pars[i]]][(warmup/thin+1):(iter/thin)]))))
    n.col <- ifelse(is.na(r_dims[2]), 1, r_dims[2])
    rs.array <- array(dim = c(r_dims[1], n.col, n.samples))
    k <- 0
    for (j in 1:n.col) {
      for (i in 1:r_dims[1]) {
        k <- k + 1
        rs.array[i,j,] <- rs[k,]
      }}
    if (center.zero) {
      center <- t(sapply(1:dim(rs.array)[2], function(i)
        unlist(lapply(1:n.samples, function(k) mean(rs.array[,i,k])))))
      for (j in 1:n.col) 
        rs.array[,j,] <- rs.array[,j,] - 
        matrix(center[j,], nrow = r_dims[1], ncol = n.samples, byrow = TRUE)
    }
    out <- get.estimate(estimate, samples = rs.array, margin = 1:2, ...)
    colnames(out) <- r.names
    if(var) {
      Var <- array(dim = c(rep(n.col, 2), r_dims[1]), 
                   dimnames = list(r.names, r.names, 1:r_dims[1]))
      for (i in 1:r_dims[1])
        if (is.na(r_dims[2])) Var[,,i] <- var(rs.array[i,]) 
      else Var[,,i] <- cov(t(rs.array[i,,])) 
      attr(out, "var") <- Var
    }
    out
  })
  names(ranef) <- group
  ranef 
} 

#' @export
VarCorr.brmsfit <- function(x, estimate = "mean", as.list = TRUE, ...) {
  if (!is(x$fit, "stanfit")) 
    stop("Argument x does not contain posterior samples")
  pars <- names(x$fit@sim$samples[[1]])
  iter <- attr(x$fit@sim$samples[[1]],"args")$iter
  warmup <- attr(x$fit@sim$samples[[1]],"args")$warmup
  thin <- attr(x$fit@sim$samples[[1]],"args")$thin
  chains <- length(x$fit@sim$samples) 
  group <- unlist(extract.effects(x$formula, add.ignore = TRUE)$group)
  
  VarCorr <- lapply(group, function(g) {
    sd.pars <- pars[grepl(paste0("^sd_",g,"_"), pars)]
    cor.pars <- pars[grepl(paste0("^cor_",g,"_"), pars)]
    r.names <- gsub(paste0("^sd_",g,"_"), "", sd.pars)
    r.names <- gsub("__",":",r.names)
    nr <- length(sd.pars)
    if (!nr)
      stop(paste("Grouping variable",g,"is not present in argument x"))
    out <- list() 
    sds <- t(sapply(1:nr, function(i)
      unlist(lapply(1:chains, function(j) 
        x$fit@sim$samples[[j]][[sd.pars[i]]][(warmup/thin+1):(iter/thin)]))))
    out$sd <- do.call(cbind, lapply(estimate, get.estimate, samples = sds, ...))
    rownames(out$sd) <- r.names 
    if (length(cor.pars)) {
      cors <- t(sapply(1:length(cor.pars), function(i)
        unlist(lapply(1:chains, function(j) 
          x$fit@sim$samples[[j]][[cor.pars[i]]][(warmup/thin+1):(iter/thin)])))) 
      out$cor <- array(diag(1,nr), dim = c(nr, nr, (iter-warmup)/thin*chains))
      out$cov <- out$cor
      k <- 0 
      for (i in 1:nr) {
        for (j in 1:i) {
          if (i == j) out$cov[i,j,] <- sds[i,]^2
          else {
            k = k + 1
            out$cor[i,j,] <- cors[k,]
            out$cor[j,i,] <- out$cor[i,j,]
            out$cov[i,j,] <- out$cor[i,j,] * sds[i,] * sds[j,]
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
  })
  names(VarCorr) <- group
  VarCorr
}

#' Create a summary of a fitted model represented by a \code{brmsfit} object
#' 
#' Summarize estimated fixed and random effects as well as other useful
#' results included in a \code{brmsfit} object.
#'
#' @param object An object of class \code{brmsfit}
#' @param ... Other potential arguments
#' 
#' @author Paul-Christian Buerkner <\email{paul.buerkner@@gmail.com}>
#' 
#' @export
summary.brmsfit <- function(object, ...) {
  if (!is(object$fit, "stanfit")) 
    out <- brmssummary(formula = object$formula, family = object$family, link = object$link,
             data.name = object$data.name, 
             group = unlist(extract.effects(object$formula, add.ignore = TRUE)$group),
             nobs = nobs(object), ngrps = brms::ngrps(object), autocor = object$autocor)
  else {
    out <- brmssummary(formula = object$formula, family = object$family, link = object$link,
             data.name = object$data.name, 
             group = unlist(extract.effects(object$formula, add.ignore = TRUE)$group),
             nobs = nobs(object), ngrps <- ngrps(object), autocor = object$autocor,
             n.chain = length(object$fit@sim$samples),
             n.iter = attr(object$fit@sim$samples[[1]],"args")$iter,
             n.warmup = attr(object$fit@sim$samples[[1]],"args")$warmup,
             n.thin = attr(object$fit@sim$samples[[1]],"args")$thin,
             sampler = attr(object$fit@sim$samples[[1]],"args")$sampler_t) 
    pars <- names(object$fit@sim$samples[[1]])
    fit.summary <- rstan::summary(object$fit, probs = c(0.025, 0.975))
    col.names <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Eff.Sample", "Rhat")
    
    fix.pars <- pars[grepl("^b_", pars)]
    out$fixed <- matrix(fit.summary$summary[fix.pars,-c(2)], ncol = 6)
    colnames(out$fixed) <- col.names
    rownames(out$fixed) <- gsub("__",":",gsub("^b_","",fix.pars))
    
    spec.pars <- pars[pars %in% c("sigma","nu","shape","delta")]
    out$spec.pars <- matrix(fit.summary$summary[spec.pars,-c(2)], ncol = 6)
    colnames(out$spec.pars) <- col.names
    rownames(out$spec.pars) <- spec.pars
    
    cor.pars <- pars[grepl("^ar|^ma", pars)]
    out$cor.pars <- matrix(fit.summary$summary[cor.pars,-c(2)], ncol = 6)
    colnames(out$cor.pars) <- col.names
    rownames(out$cor.pars) <- cor.pars
    
    if (length(out$group)) {
      for (i in 1:length(out$group)) {
        sd.pars <- pars[grepl(paste0("^sd_", out$group[i]), pars)]
        cor.pars <- pars[grepl(paste0("^cor_", out$group[i]), pars)]
        r.names <- gsub(paste0("^sd_",out$group[i],"_"), "", sd.pars)
        sd.names <- paste0("sd(",r.names,")")
        cor.names <- get.cor.names(r.names)
        out$random[[out$group[i]]] <- matrix(fit.summary$summary[c(sd.pars, cor.pars),-c(2)], ncol = 6)
        colnames(out$random[[out$group[i]]]) <- col.names
        rownames(out$random[[out$group[i]]]) <- gsub("__",":",c(sd.names,cor.names)) 
      }
    }
  }  
  out
}

#' @export
print.brmssummary <- function(x, digits = 2, ...) {
  cat(paste0(" Family: ", x$family, " (", x$link, ") \n"))
  cat(paste("Formula:", gsub(" {1,}", " ", Reduce(paste, deparse(x$formula))), "\n"))
  cat(paste0("   Data: ", x$data.name, " (Number of observations: ",x$nobs,") \n"))
  if (x$sampler == "") {
    cat(paste("\nArgument x does not contain posterior samples. Most likely, this is \n",
        "because the package rstan was not installed when the model was fitted. \n",
        "Please see https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started \n",
        "for instructions on how to install rstan."))
  }
  else {
    cat(paste0("Samples: ", x$n.chain, " chains, each with n.iter = ", x$n.iter, 
               "; n.warmup = ", x$n.warmup, "; n.thin = ", x$n.thin, "; \n",
      "         total post-warmup samples = ", (x$n.iter-x$n.warmup)/x$n.thin*x$n.chain, "\n \n"))  
    
    if (length(x$group)) {
      cat("Random Effects: \n")
      for (i in 1:length(x$group)) {
        cat(paste0("~",x$group[i], " (Number of levels: ",x$ngrps[[x$group[i]]],") \n"))
        x$random[[x$group[i]]][,"Eff.Sample"] <- round(x$random[[x$group[i]]][,"Eff.Sample"], digits = 0)
        print(round(x$random[[x$group[i]]], digits = digits))
        cat("\n")
      }
    }
    
    if (nrow(x$cor.pars)) {
      cat("Correlation Structure: "); print(x$autocor); cat("\n")
      x$cor.pars[,"Eff.Sample"] <- round(x$cor.pars[,"Eff.Sample"], digits = 0)
      print(round(x$cor.pars, digits = digits))
      cat("\n")
    }
    
    cat("Fixed Effects: \n")
    x$fixed[,"Eff.Sample"] <- round(x$fixed[,"Eff.Sample"], digits = 0)
    print(round(x$fixed, digits = digits)) 
    cat("\n")
    
    if (nrow(x$spec.pars)) {
      cat("Family Specific Parameters: \n")
      x$spec.pars[,"Eff.Sample"] <- round(x$spec.pars[,"Eff.Sample"], digits = 0)
      print(round(x$spec.pars, digits = digits))
      cat("\n")
    }
    
    cat(paste0("Samples were drawn using ",x$sampler,". For each parameter, Eff.Sample is a \n",
        "crude measure of effective sample size, and Rhat is the potential scale \n",
        "reduction factor on split chains (at convergence, Rhat = 1)."))
  }
}  

#' Print a summary for a fitted model represented by a \code{brmsfit} object
#'
#' Print basic information regarding the fitted model and a summary for the fixed and random effects
#' estimated by the samples included in a \code{brmsfit} object.
#' 
#' @aliases print.brmssummary
#' 
#' @param x An object of class \code{brmsfit}
#' @param digits The number of significant digits for printing out the summary; defaults to 2. 
#'   The effective sample size is always rounded to integers.
#' @param ... Additional arguments that would be passed to method \code{summary} of \code{brmsfit}.
#'
#' @author Paul-Christian Buerkner <\email{paul.buerkner@@gmail.com}>
#' 
#' @export
print.brmsfit <- function(x, digits = 2, ...) {
  print(summary(x), digits = digits, ...)
}  

#' @export
nobs.brmsfit <- function(object, ...) length(object$data$Y)

#' @export
ngrps.brmsfit <- function(object, ...) {
  group <- unlist(extract.effects(object$formula, add.ignore = TRUE)$group)
  setNames(lapply(group, function(g) object$data[[paste0("N_",g)]]), group)
}

#' @export
formula.brmsfit <- function(x, ...) x$formula

#' @export 
predict.brmsfit <- function(object, ...) {
  if (!"Y_pred" %in% object$fit@model_pars) 
    stop(paste0("The model does not contain predicted values. \n",
         "You should use argument predict = TRUE in function brm."))
  else {
    pars <- names(object$fit@sim$samples[[1]])
    fit.summary <- rstan::summary(object$fit, probs = c(0.025, 0.975))
    pred.pars <- pars[grepl("^Y_pred\\[", pars)]
    out <- matrix(fit.summary$summary[pred.pars,-c(2,6,7)], ncol = 4)
    rownames(out) <- NULL
    colnames(out) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI")
  } 
  out
}

#' @export
print.brmsmodel <- function(x, ...) {
  cat(x)
}

#' Trace and density plots for MCMC samples
#' 
#' Trace and density plots for MCMC samples using the \code{ggmcmc} package
#' 
#' @param x An object of class \code{brmsfit}.
#' @param parameters Name of the parameters to plot, as given by a character vector or a regular expression.
#'   By default, all parameters except for random effects and posterior predictives are plotted. 
#' @param combine logical; Indicates if the samples of all chains should be combined into one posterior distribution. 
#' @param N The number of parameters plotted per page.
#' @param ask logical; Indicates if the user is prompted before a new page is plotted.   
#' @param ... Currently ignored.
#' 
#' @return NULL
#' 
#' @author Paul-Christian Buerkner <\email{paul.buerkner@@gmail.com}>
#' 
#' @examples
#' \dontrun{ 
#' fit_e <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1|patient) + (1|visit), 
#'              data = epilepsy, family = "poisson")
#' ## plot fixed effects as well as standard devations and correlations (if present) of random effects
#' plot(fit_e)
#' ## plot fixed effects only and combine the chains into one posterior
#' plot(fit_e, parameters = "^b_", combine = TRUE) 
#' }
#' 
#' @import ggplot2
#' @export
plot.brmsfit <- function(x, parameters = NA, combine = FALSE, N = 5, ask = TRUE, ...) {
  if (!is(x$fit, "stanfit")) 
    stop("Argument x does not contain posterior samples")
  if (is.na(parameters)) 
    parameters <- c("^b_", "^sd_", "^cor_", "^sigma$", "^nu$", 
                    "^shape$", "^delta$", "^ar", "^ma")
  
  pars <- sort(names(x$fit@sim$samples[[1]]))
  pars <- gsub("__", ":", pars[apply(sapply(parameters, grepl, x = pars), 1, any)])
  pfit <- ggmcmc::ggs(x$fit)
  pfit$Parameter <- gsub("__", ":", pfit$Parameter)
  att <- attributes(pfit)
  rel.att <- c("class", "nChains", "nIterations", "nBurnin", "nThin", "description")
  pfit <- pfit[which(pfit$Parameter %in% pars),]
  
  default.ask <- devAskNewPage()
  devAskNewPage(ask = FALSE)
  for (i in 1:ceiling(length(pars)/N)) {
    pfit.sub1 <- pfit[which(pfit$Parameter %in% pars[((i-1)*N+1):min(i*N,length(pars))]),]
    for (j in 1:length(rel.att)) 
      attr(pfit.sub1, rel.att[j]) <- att[[rel.att[j]]]
    pfit.sub2 <- pfit.sub1
    if (combine) pfit.sub2$Chain <- 1
    gridExtra::grid.arrange(ggmcmc::ggs_traceplot(pfit.sub1) + 
        ggplot2::theme(legend.position = "none"), 
        ggmcmc::ggs_density(pfit.sub2), ncol = 2, nrow = 1)
    if (i == 1) devAskNewPage(ask = ask)
  }
  devAskNewPage(default.ask)
}