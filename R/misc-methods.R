#' @export
print.brmssummary <- function(x, digits = 2, ...) {
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
    waic <- ifelse(is.numeric(x$WAIC), round(x$WAIC, digits = digits), x$WAIC)
    cat(paste0("Samples: ", x$chains, " chains, each with iter = ", x$iter, 
               "; warmup = ", x$warmup, "; thin = ", x$thin, "; \n",
               "         total post-warmup samples = ", final_samples, "\n"))
    cat(paste0("   WAIC: ", waic, "\n \n"))
    
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

#' @rdname VarCorr.brmsfit
#' @export
as.data.frame.brmsVarCorr <- function(x, ...) {
  estimates <- colnames(x[[1]]$sd)
  groups <- names(x)
  n_groups <- length(groups)
  names_coef <- lapply(x, function(y) rownames(y$sd))
  groups_col <- ulapply(1:n_groups, function(i) 
    c(groups[i], rep("", length(names_coef[[i]]) - 1)))
  max_cor <- max(ulapply(names_coef, length)) - 1
  # basic data.frame to be used in fill_base_frame
  base_frame <- as.data.frame(matrix(NA, nrow = length(groups_col),
                                     ncol = 4 + 2 * max_cor))
  names(base_frame) <- c("Group", "Name", "Std.Dev", rep("Cor", max_cor),
                         rep("Cov", max_cor + 1))
  base_frame[, 1:2] <- cbind(groups_col, unlist(names_coef))
  
  fill_base_frame <- function(estimate) {
    # fills the base_frame with SD and COR estimates
    # Args:
    #   estimate: The estimate being applied on the SD and COR parameters
    out <- base_frame
    pos <- 1
    for (i in 1:n_groups) {
      len <- length(names_coef[[i]])
      rows <- pos:(pos + len - 1)
      out[rows, "Std.Dev"] <- x[[i]]$sd[, estimate]
      if (len > 1) {
        # covariances and correlations present; add correlations
        cor_pos <- 4:(2 + len)
        cormat <- x[[i]]$cor[[estimate]][2:len, 1:(len-1), drop = FALSE]
        lt <- lower.tri(cormat, diag = TRUE)
        out[rows[2:length(rows)], cor_pos][lt] <- cormat[lt]
      }
      # add covariances
      cov_pos <- (4 + max_cor):(3 + max_cor + len)
      covmat <- x[[i]]$cov[[estimate]]
      lt <- lower.tri(covmat, diag = TRUE)
      out[rows, cov_pos][lt] <- covmat[lt]
      pos <- pos + len
    }
    out
  }
  
  out <- do.call(rbind, lapply(estimates, fill_base_frame))
  estimates_col <- ulapply(estimates, function(e)
    c(e, rep("", length(groups_col) - 1)))
  cbind(Estimate = estimates_col, out)
}

#' @export
print.brmsVarCorr <- function(x, digits = 2, ...) {
  dat <- as.data.frame(x)
  dat[, 4:ncol(dat)] <- round(as.matrix(dat[, 4:ncol(dat)]), digits = digits)
  dat[is.na(dat)] <- ""
  print(dat, row.names = FALSE, ...)
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

#' @rdname hypothesis
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
  