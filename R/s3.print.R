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
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @export
print.brmsfit <- function(x, digits = 2, ...) {
  print(summary(x), digits = digits, ...)
}  

#' @export
print.brmssummary <- function(x, digits = 2, ...) {
  cat(paste0(" Family: ", x$family, " (", x$link, ") \n"))
  cat(paste("Formula:", gsub(" {1,}", " ", Reduce(paste, deparse(x$formula))), "\n"))
  cat(paste0("   Data: ", x$data.name, " (Number of observations: ",x$nobs,") \n"))
  if (x$sampler == "") {
    cat(paste("\nThe model does not contain posterior samples."))
  }
  else {
    cat(paste0("Samples: ", x$n.chains, " chains, each with n.iter = ", x$n.iter, 
               "; n.warmup = ", x$n.warmup, "; n.thin = ", x$n.thin, "; \n",
               "         total post-warmup samples = ", (x$n.iter-x$n.warmup)/x$n.thin*x$n.chains, "\n"))
    cat(paste0("   WAIC: ", ifelse(is.numeric(x$WAIC), round(x$WAIC, digits = digits), x$WAIC), "\n \n"))
    
    if (length(x$group)) {
      cat("Random Effects: \n")
      for (i in 1:length(x$group)) {
        cat(paste0("~",x$group[i], " (Number of levels: ",x$ngrps[[x$group[i]]],") \n"))
        x$random[[x$group[i]]][,"Eff.Sample"] <- round(x$random[[x$group[i]]][,"Eff.Sample"], digits = 0)
        print(round(x$random[[x$group[i]]], digits = digits))
        cat("\n")
      }
    }
    
    if (nrow(x$cor_pars)) {
      cat("Correlation Structure: "); print(x$autocor); cat("\n")
      x$cor_pars[,"Eff.Sample"] <- round(x$cor_pars[,"Eff.Sample"], digits = 0)
      print(round(x$cor_pars, digits = digits))
      cat("\n")
    }
    
    cat("Fixed Effects: \n")
    x$fixed[,"Eff.Sample"] <- round(x$fixed[,"Eff.Sample"], digits = 0)
    print(round(x$fixed, digits = digits)) 
    cat("\n")
    
    if (nrow(x$spec_pars)) {
      cat("Family Specific Parameters: \n")
      x$spec_pars[,"Eff.Sample"] <- round(x$spec_pars[,"Eff.Sample"], digits = 0)
      print(round(x$spec_pars, digits = digits))
      cat("\n")
    }
    
    cat(paste0("Samples were drawn using ",x$sampler,". For each parameter, Eff.Sample is a \n",
               "crude measure of effective sample size, and Rhat is the potential scale \n",
               "reduction factor on split chains (at convergence, Rhat = 1)."))
  }
}  

#' @export
print.brmshypothesis <- function(x, digits = 2, ...) {
  cat(paste0("Hypothesis Tests for class ", x$class, ":\n"))
  x$hypothesis[,1:5] <- round(x$hypothesis[,1:5], digits = digits)
  print(x$hypothesis, quote = FALSE)
  cat(paste0("---\n'*': The expected value under the hypothesis lies outside the ",(1-x$alpha)*100,"% CI."))
}

#' @export
print.brmsmodel <- function(x, ...) cat(x)

#' @export
print.ic <- function(x, digits = 2, ...) {
  ic <- names(x)[3]
  mat <- matrix(c(x[[ic]], x[[paste0("se_",ic)]]), ncol = 2, 
                dimnames = list("", c(toupper(ic), "SE")))
  print(round(mat, digits = digits))
}

#' @export
print.iclist <- function(x, digits = 2, ...) {
  ic <- names(x[[1]])[3]
  mat <- matrix(0, nrow = length(x), ncol = 2, 
                dimnames = list(names(x), c(toupper(ic), "SE")))
  for (i in 1:length(x))
    mat[i, ] <- c(x[[i]][[ic]], x[[i]][[paste0("se_",ic)]])
  if (is.matrix(attr(x, "compare")))
    mat <- rbind(mat, attr(x, "compare"))
  print(round(mat, digits = digits))
}