#' @export
print.brmssummary <- function(x, digits = 2, ...) {
  cat(" Family: ")
  if (is(x$family, "family")) {
    type <- ifelse(is.null(x$family$type), "", paste(",", x$family$type))
    cat(paste0(x$family$family, " (", x$family$link, type, ") \n"))
  } else {
    cat(paste0(x$family, " (", x$link, ") \n"))  
  }
  cat(paste("Formula:", 
            gsub(" {1,}", " ", Reduce(paste, deparse(x$formula))), "\n"))
  cat(paste0("   Data: ", x$data.name, 
             " (Number of observations: ",x$nobs,") \n"))
  if (x$sampler == "") {
    cat(paste("\nThe model does not contain posterior samples."))
  }
  else {
    final_samples <- (x$n.iter - x$n.warmup) / x$n.thin * x$n.chains
    waic <- ifelse(is.numeric(x$WAIC), round(x$WAIC, digits = digits), x$WAIC)
    cat(paste0("Samples: ", x$n.chains, " chains, each with n.iter = ", x$n.iter, 
               "; n.warmup = ", x$n.warmup, "; n.thin = ", x$n.thin, "; \n",
               "         total post-warmup samples = ", final_samples, "\n"))
    cat(paste0("   WAIC: ", waic, "\n \n"))
    
    if (length(x$group)) {
      cat("Random Effects: \n")
      for (i in 1:length(x$group)) {
        g <- x$group[i]
        cat(paste0("~",g," (Number of levels: ",x$ngrps[[g]],") \n"))
        x$random[[g]][, "Eff.Sample"] <- 
          round(x$random[[g]][, "Eff.Sample"], digits = 0)
        print(round(x$random[[g]], digits = digits))
        cat("\n")
      }
    }
    
    if (nrow(x$cor_pars)) {
      cat("Correlation Structure: ")
      print(x$autocor)
      cat("\n")
      x$cor_pars[, "Eff.Sample"] <- round(x$cor_pars[, "Eff.Sample"], 
                                          digits = 0)
      print(round(x$cor_pars, digits = digits))
      cat("\n")
    }
    
    if (nrow(x$fixed)) {
      cat("Fixed Effects: \n")
      x$fixed[, "Eff.Sample"] <- round(x$fixed[, "Eff.Sample"], 
                                       digits = 0)
      print(round(x$fixed, digits = digits)) 
      cat("\n")
    }
    
    if (nrow(x$spec_pars)) {
      cat("Family Specific Parameters: \n")
      x$spec_pars[, "Eff.Sample"] <- 
        round(x$spec_pars[, "Eff.Sample"], digits = 0)
      print(round(x$spec_pars, digits = digits))
      cat("\n")
    }
    
    cat(paste0("Samples were drawn using ",x$sampler,". ", 
               "For each parameter, Eff.Sample is a \n",
               "crude measure of effective sample size, ", 
               "and Rhat is the potential scale \n",
               "reduction factor on split chains ",
               "(at convergence, Rhat = 1)."))
  }
  invisible(x)
}

#' @rdname VarCorr
#' @export
as.data.frame.VarCorr_brmsfit <- function(x, ...) {
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
        # covariances and correlations present
        # add correlations
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
print.VarCorr_brmsfit <- function(x, digits = 2, ...) {
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
  # print the output of LOO(x1, x2, ...) and WAIC(x1, x2, ...)
  ic <- names(x[[1]])[3]
  mat <- matrix(0, nrow = length(x), ncol = 2, 
                dimnames = list(names(x), c(toupper(ic), "SE")))
  for (i in 1:length(x)) { 
    mat[i, ] <- c(x[[i]][[ic]], x[[i]][[paste0("se_",ic)]])
  }
  if (is.matrix(attr(x, "compare"))) {
    # models were compared using the compare_ic function
    mat <- rbind(mat, attr(x, "compare"))
    weights <- c(attr(x, "weights"), rep(NA, nrow(attr(x, "compare")))) 
    mat <- cbind(mat, Weights = weights)
  }
  print(round(mat, digits = digits), na.print = "")
  invisible(x)
}

#' @export
print.brmshypothesis <- function(x, digits = 2, ...) {
  cat(paste0("Hypothesis Tests for class ", x$class, ":\n"))
  x$hypothesis[, 1:5] <- round(x$hypothesis[, 1:5], digits = digits)
  print(x$hypothesis, quote = FALSE)
  cat(paste0("---\n'*': The expected value under the hypothesis lies outside the ",
             (1 - x$alpha) * 100, "% CI."))
  invisible(x)
}

#' @rdname hypothesis
#' @method plot brmshypothesis
#' @export
plot.brmshypothesis <- function(x, N = 5, ignore_prior = FALSE, 
                                theme = "classic", ask = TRUE, 
                                do_plot = TRUE, newpage = TRUE, ...) {
  if (!is.data.frame(x$samples)) {
    stop("No posterior samples found")
  }
  .plot_fun <- function(i) {
    # create the ggplot object for each hypothesis
    # Args: i: index variable
    ggplot(x$samples, aes_string(x = hypnames[i])) + 
      geom_density(aes_string(fill = "Type"), alpha = 0.5, na.rm = TRUE) + 
      scale_fill_manual(values = c("red", "blue")) + 
      xlab("") + ylab("") + ggtitle(hyps[i]) + 
      do.call(paste0("theme_", theme), args = list())
  }
  if (ignore_prior) {
    x$samples <- subset(x$samples, x$samples$Type == "posterior")
  }
  if (do_plot) {
    default_ask <- devAskNewPage()
    on.exit(devAskNewPage(default_ask))
    devAskNewPage(ask = FALSE)
  }
  hyps <- rownames(x$hypothesis)
  hypnames <- names(x$samples)[seq_along(hyps)]
  n_plots <- ceiling(length(hyps) / N)
  plots <- vector(mode = "list", length = n_plots)
  for (i in 1:n_plots) {
    I <- ((i - 1) * N + 1):min(i * N, length(hyps))
    temp_plot <- lapply(I, .plot_fun)
    plots[[i]] <- arrangeGrob(grobs = temp_plot, ncol = 1, 
                              nrow = length(temp_plot), ...)
    if (do_plot) {
      if (newpage || i > 1) grid.newpage()
      grid.draw(plots[[i]])
      if (i == 1) devAskNewPage(ask = ask)
    }
  }
  if (do_plot) {
    invisible(plots) 
  } else {
    plots
  }
}