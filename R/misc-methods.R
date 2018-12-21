#' @export
print.brmssummary <- function(x, digits = 2, ...) {
  cat(" Family: ")
  cat(summarise_families(x$formula), "\n")
  cat("  Links: ")
  cat(summarise_links(x$formula, wsp = 9), "\n")
  cat("Formula: ")
  print(x$formula, wsp = 9)
  cat(paste0(
    "   Data: ", x$data.name, 
    " (Number of observations: ", x$nobs, ") \n"
  ))
  if (!isTRUE(nzchar(x$sampler))) {
    cat("\nThe model does not contain posterior samples.\n")
  } else {
    final_samples <- ceiling((x$iter - x$warmup) / x$thin * x$chains)
    cat(paste0(
      "Samples: ", x$chains, " chains, each with iter = ", x$iter, 
      "; warmup = ", x$warmup, "; thin = ", x$thin, ";\n",
      "         total post-warmup samples = ", final_samples, "\n\n"
    ))
    if (nrow(x$prior)) {
      cat("Priors: \n")
      print(x$prior, show_df = FALSE)
      cat("\n")
    }
    if (length(x$splines)) {
      cat("Smooth Terms: \n")
      print_format(x$splines, digits)
      cat("\n")
    }
    if (length(x$gp)) {
      cat("Gaussian Process Terms: \n")
      print_format(x$gp, digits)
      cat("\n")
    }
    if (nrow(x$cor_pars)) {
      cat("Correlation Structures:\n")
      # TODO: better printing for correlation structures?
      print_format(x$cor_pars, digits)
      cat("\n")
    }
    if (length(x$random)) {
      cat("Group-Level Effects: \n")
      for (i in seq_along(x$random)) {
        g <- names(x$random)[i]
        cat(paste0("~", g, " (Number of levels: ", x$ngrps[[g]], ") \n"))
        print_format(x$random[[g]], digits)
        cat("\n")
      }
    }
    if (nrow(x$fixed)) {
      cat("Population-Level Effects: \n")
      print_format(x$fixed, digits)
      cat("\n")
    }
    if (length(x$mo)) {
      cat("Simplex Parameters: \n")
      print_format(x$mo, digits)
      cat("\n")
    }
    if (nrow(x$spec_pars)) {
      cat("Family Specific Parameters: \n")
      print_format(x$spec_pars, digits)
      cat("\n")
    }
    if (length(x$rescor_pars)) {
      cat("Residual Correlations: \n")
      print_format(x$rescor, digits)
      cat("\n")
    }
    cat(paste0("Samples were drawn using ", x$sampler, ". "))
    if (x$algorithm == "sampling") {
      cat(paste0(
        "For each parameter, Eff.Sample \n",
        "is a crude measure of effective sample size, ", 
        "and Rhat is the potential \n",
        "scale reduction factor on split chains ",
        "(at convergence, Rhat = 1)."
      ))
    }
    cat("\n")
  }
  invisible(x)
}

print_format <- function(x, digits = 2, no_digits = "Eff.Sample") {
  # helper function to print summary matrices
  # in nice format, also showing -0.00 (#263)
  # Args:
  #   x: object to be printed; coerced to matrix
  #   digits: number of digits to show
  #   no_digits: names of columns for which no digits should be shown
  x <- as.matrix(x)
  digits <- as.numeric(digits)
  if (length(digits) != 1L) {
    stop2("'digits' should be a single numeric value.")
  }
  out <- x
  fmt <- paste0("%.", digits, "f")
  for (i in seq_cols(x)) {
    if (isTRUE(colnames(x)[i] %in% no_digits)) {
      out[, i] <- sprintf("%.0f", x[, i])
    } else {
      out[, i] <- sprintf(fmt, x[, i])
    }
  }
  print(out, quote = FALSE, right = TRUE)
  invisible(x)
}

#' @export
print.brmsmodel <- function(x, ...) {
  cat(x)
  invisible(x) 
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

#' @rdname posterior_summary
#' @export
posterior_summary.default <- function(x, probs = c(0.025, 0.975), 
                                      robust = FALSE, ...) {
  if (robust) {
    coefs <- c("median", "mad", "quantile")
  } else {
    coefs <- c("mean", "sd", "quantile")
  }
  .posterior_summary <- function(x) {
    run(cbind, lapply(
      coefs, get_estimate, samples = x, 
      probs = probs, na.rm = TRUE
    ))
  }
  if (length(dim(x)) <= 2L) {
    # data.frames cause trouble in as.array
    x <- as.matrix(x)
  } else {
    x <- as.array(x) 
  }
  if (length(dim(x)) == 2L) {
    out <- .posterior_summary(x)
    rownames(out) <- colnames(x)
  } else if (length(dim(x)) == 3L) {
    out <- lapply(array2list(x), .posterior_summary)
    out <- abind(out, along = 3)
    dnx <- dimnames(x)
    dimnames(out) <- list(dnx[[2]], dimnames(out)[[2]], dnx[[3]])
  } else {
    stop("'x' must be of dimension 2 or 3.")
  }
  colnames(out) <- c("Estimate", "Est.Error", paste0("Q", probs * 100))
  out  
}
