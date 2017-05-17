#' @export
print.brmssummary <- function(x, digits = 2, ...) {
  if (!is.null(x[["WAIC"]])) {
    # deprecated as of 1.7.0
    x[["waic"]] <- x[["WAIC"]]
  }
  cat(" Family: ")
  if (is.family(x$family)) {
    cat(summary(x$family), "\n")
  } else {
    cat(paste0(x$family, " (", x$link, ") \n"))  
  }
  cat("Formula: ")
  print(x$formula, wsp = 9)
  cat(paste0("   Data: ", x$data.name, 
             " (Number of observations: ", x$nobs, ") \n"))
  if (!nzchar(x$sampler)) {
    cat("\nThe model does not contain posterior samples.\n")
  } else {
    if (!is.null(x$n.iter)) {
      # deprecated names are used
      args <- c("iter", "warmup", "thin", "chains")
      x[args] <- x[paste0("n.", args)]
    }
    final_samples <- ceiling((x$iter - x$warmup) / x$thin * x$chains)
    if (is.numeric(x$loo)) {
      x$loo <- round(x$loo, digits = digits)
    }
    if (is.numeric(x$waic)) {
      x$waic <- round(x$waic, digits = digits)
    }
    cat(paste0("Samples: ", x$chains, " chains, each with iter = ", x$iter, 
               "; warmup = ", x$warmup, "; thin = ", x$thin, "; \n",
               "         total post-warmup samples = ", final_samples, "\n"))
    cat(paste0("    ICs: LOO = ", x$loo, "; WAIC = ", x$waic, "\n \n"))
    
    if (nrow(x$prior)) {
      cat("Priors: \n")
      print(x$prior, show_df = FALSE)
      cat("\n")
    }
    
    if (length(x$splines)) {
      cat("Smooth Terms: \n")
      if (x$algorithm == "sampling") {
        x$splines[, "Eff.Sample"] <- 
          round(x$splines[, "Eff.Sample"], digits = 0)
      }
      print(round(x$splines, digits = digits)) 
      cat("\n")
    }
    
    if (length(x$gp)) {
      cat("Gaussian Process Terms: \n")
      if (x$algorithm == "sampling") {
        x$gp[, "Eff.Sample"] <- 
          round(x$gp[, "Eff.Sample"], digits = 0)
      }
      print(round(x$gp, digits = digits)) 
      cat("\n")
    }
    
    if (length(x$random)) {
      cat("Group-Level Effects: \n")
      for (i in seq_along(x$random)) {
        g <- names(x$random)[i]
        cat(paste0("~", g, " (Number of levels: ", x$ngrps[[g]], ") \n"))
        if (x$algorithm == "sampling") {
          x$random[[g]][, "Eff.Sample"] <- 
            round(x$random[[g]][, "Eff.Sample"], digits = 0)
        }
        print(round(x$random[[g]], digits = digits))
        cat("\n")
      }
    }
    
    if (nrow(x$cor_pars)) {
      cat("Correlation Structure: ")
      print(x$autocor)
      cat("\n")
      if (x$algorithm == "sampling") {
        x$cor_pars[, "Eff.Sample"] <- 
          round(x$cor_pars[, "Eff.Sample"], digits = 0)
      }
      print(round(x$cor_pars, digits = digits))
      cat("\n")
    }
    
    if (nrow(x$fixed)) {
      cat("Population-Level Effects: \n")
      if (x$algorithm == "sampling") {
        x$fixed[, "Eff.Sample"] <- 
          round(x$fixed[, "Eff.Sample"], digits = 0)
      }
      print(round(x$fixed, digits = digits)) 
      cat("\n")
    }
    
    if (length(x$mo)) {
      cat("Simplex Parameters: \n")
      if (x$algorithm == "sampling") {
        x$mo[, "Eff.Sample"] <- 
          round(x$mo[, "Eff.Sample"], digits = 0)
      }
      print(round(x$mo, digits = digits)) 
      cat("\n")
    }
    
    if (nrow(x$spec_pars)) {
      cat("Family Specific Parameters: \n")
      if (x$algorithm == "sampling") {
        x$spec_pars[, "Eff.Sample"] <- 
          round(x$spec_pars[, "Eff.Sample"], digits = 0)
      }
      print(round(x$spec_pars, digits = digits))
      cat("\n")
    }
    
    cat(paste0("Samples were drawn using ", x$sampler, ". "))
    if (x$algorithm == "sampling") {
      cat(paste0("For each parameter, Eff.Sample \n",
          "is a crude measure of effective sample size, ", 
          "and Rhat is the potential \n",
          "scale reduction factor on split chains ",
          "(at convergence, Rhat = 1)."))
    }
    cat("\n")
  }
  invisible(x)
}

#' @export
print.brmsmodel <- function(x, ...) {
  cat(x)
  invisible(x) 
}

#' @export
print.ic <- function(x, digits = 2, ...) {
  # print the output of LOO(x) and WAIC(x)
  ic <- names(x)[3]
  mat <- matrix(c(x[[ic]], x[[paste0("se_",ic)]]), ncol = 2, 
                dimnames = list("", c(toupper(ic), "SE")))
  print(round(mat, digits = digits))
  invisible(x)
}

#' @export
print.iclist <- function(x, digits = 2, ...) {
  # print the output of LOO and WAIC with multiple models
  m <- x
  m$ic_diffs__ <- NULL
  if (length(m)) {
    ic <- names(m[[1]])[3]
    mat <- matrix(0, nrow = length(m), ncol = 2)
    dimnames(mat) <- list(names(m), c(toupper(ic), "SE"))
    for (i in seq_along(m)) { 
      mat[i, ] <- c(m[[i]][[ic]], m[[i]][[paste0("se_", ic)]])
    }
  } else {
    mat <- NULL
  }
  ic_diffs <- x$ic_diffs__
  if (is.matrix(attr(x, "compare"))) {
    # deprecated as of brms 1.4.0
    ic_diffs <- attr(x, "compare")
  }
  if (is.matrix(ic_diffs)) {
    # models were compared using the compare_ic function
    mat <- rbind(mat, ic_diffs)
  }
  print(round(mat, digits = digits), na.print = "")
  invisible(x)
}

#' @rdname brmshypothesis
#' @export
print.brmshypothesis <- function(x, digits = 2, chars = 20, ...) {
  # make sure rownames are not too long
  rnames <- limit_chars(rownames(x$hypothesis), chars = chars)
  rownames(x$hypothesis) <- make.unique(rnames, sep = " #")
  cat(paste0("Hypothesis Tests for class ", x$class, ":\n"))
  x$hypothesis[, 1:5] <- round(x$hypothesis[, 1:5], digits = digits)
  print(x$hypothesis, quote = FALSE)
  cat(paste0("---\n'*': The expected value under the hypothesis ", 
             "lies outside the ", (1 - x$alpha) * 100, "% CI."))
  invisible(x)
}

#' @export
print.brmsMarginalEffects <- function(x, ...) {
  plot(x, ...)
}

#' @export
parnames.default <- function(x, ...) {
  names(x)
}

#' @export
posterior_samples.default <- function(x, pars = NA, exact_match = FALSE, ...) {
  x <- as.data.frame(x)
  if (!anyNA(pars)) {
    pars <- extract_pars(
      pars, all_pars = names(x), exact_match = exact_match, ...
    )
    x <- x[, pars, drop = FALSE]
  }
  if (!ncol(x)) {
    x <- NULL
  }
  x
}

#' @export
prior_samples.default <- function(x, pars = NA, exact_match = FALSE, ...) {
  if (anyNA(pars)) {
    pars <- "^prior_"
    exact_match <- FALSE
  } else {
    if (exact_match) {
      pars <- paste0("prior_", pars) 
    } else {
      hat <- substr(pars, 1, 1) == "^"
      pars <- ifelse(hat, substr(pars, 2, nchar(pars)), pars)
      pars <- paste0("^prior_", pars)  
    }
  }
  posterior_samples(x, pars = pars, exact_match = exact_match, ...)
}

#' @rdname hypothesis
#' @export
hypothesis.default <- function(x, hypothesis, alpha = 0.05, ...) {
  hypothesis_internal(x, hypothesis, alpha = alpha, ...)
}
