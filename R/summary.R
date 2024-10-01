#' Create a summary of a fitted model represented by a \code{brmsfit} object
#'
#' @param object An object of class \code{brmsfit}.
#' @param priors Logical; Indicating if priors should be included
#'   in the summary. Default is \code{FALSE}.
#' @param prob A value between 0 and 1 indicating the desired probability
#'   to be covered by the uncertainty intervals. The default is 0.95.
#' @param mc_se Logical; Indicating if the uncertainty in \code{Estimate}
#'   caused by the MCMC sampling should be shown in the summary. Defaults to
#'   \code{FALSE}.
#' @param ... Other potential arguments
#' @inheritParams posterior_summary
#'
#' @details The convergence diagnostics \code{Rhat}, \code{Bulk_ESS}, and
#' \code{Tail_ESS} are described in detail in Vehtari et al. (2020).
#'
#' @references
#' Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
#' Paul-Christian Bürkner (2020). Rank-normalization, folding, and
#' localization: An improved R-hat for assessing convergence of
#' MCMC. *Bayesian Analysis*. 1–28. dpi:10.1214/20-BA1221
#'
#' @method summary brmsfit
#' @importFrom posterior subset_draws summarize_draws
#' @export
summary.brmsfit <- function(object, priors = FALSE, prob = 0.95,
                            robust = FALSE, mc_se = FALSE, ...) {
  priors <- as_one_logical(priors)
  probs <- validate_ci_bounds(prob)
  robust <- as_one_logical(robust)
  mc_se <- as_one_logical(mc_se)
  object <- restructure(object)
  bterms <- brmsterms(object$formula)
  out <- list(
    formula = object$formula,
    data_name = get_data_name(object$data),
    group = unique(object$ranef$group),
    nobs = nobs(object),
    ngrps = ngrps(object),
    autocor = object$autocor,
    prior = empty_prior(),
    algorithm = algorithm(object)
  )
  class(out) <- "brmssummary"

  # check if the model contains any posterior draws
  model_is_empty <- !length(object$fit@sim) ||
    isTRUE(object$fit@sim$iter <= object$fit@sim$warmup)
  if (model_is_empty) {
    return(out)
  }

  stan_args <- object$fit@stan_args[[1]]
  out$sampler <- paste0(stan_args$method, "(", stan_args$algorithm, ")")
  if (priors) {
    out$prior <- prior_summary(object, all = FALSE)
  }

  variables <- variables(object)
  incl_classes <- c(
    "b", "bs", "bcs", "bsp", "bmo", "bme", "bmi", "bm",
    valid_dpars(object), "delta", "lncor", "rescor", "ar", "ma", "sderr",
    "cosy", "cortime", "lagsar", "errorsar", "car", "sdcar", "rhocar",
    "sd", "cor", "df", "sds", "sdgp", "lscale", "simo"
  )
  incl_regex <- paste0("^", regex_or(incl_classes), "(_|$|\\[)")
  variables <- variables[grepl(incl_regex, variables)]
  draws <- as_draws_array(object, variable = variables)

  out$total_ndraws <- ndraws(draws)
  out$chains <- nchains(object)
  if (length(object$fit@sim$iter)) {
    # MCMC algorithms
    out$iter <- object$fit@sim$iter
    out$warmup <- object$fit@sim$warmup
  } else {
    # non-MCMC algorithms
    out$iter <- out$total_ndraws
    out$warmup <- 0
  }
  out$thin <- nthin(object)

  # compute a summary for given set of parameters
  # TODO: align names with summary outputs of other methods and packages
  .summary <- function(draws, variables, probs, robust) {
    # quantiles with appropriate names to retain backwards compatibility
    .quantile <- function(x, ...) {
      qs <- posterior::quantile2(x, probs = probs, ...)
      prob <- probs[2] - probs[1]
      names(qs) <- paste0(c("l-", "u-"), prob * 100, "% CI")
      return(qs)
    }
    draws <- subset_draws(draws, variable = variables)
    measures <- list()
    if (robust) {
      measures$Estimate <- median
      if (mc_se) {
        measures$MCSE <- posterior::mcse_median
      }
      measures$Est.Error <- mad
    } else {
      measures$Estimate <- mean
      if (mc_se) {
        measures$MCSE <- posterior::mcse_mean
      }
      measures$Est.Error <- sd
    }
    c(measures) <- list(quantiles = .quantile)
    if (!is.brmsfit_multiple(object)) {
      # TODO: add a viable post-processing solution for brm_multiple models
      c(measures) <- list(
        Rhat = posterior::rhat,
        Bulk_ESS = posterior::ess_bulk,
        Tail_ESS = posterior::ess_tail
      )
    }
    out <- do.call(summarize_draws, c(list(draws), measures))
    out <- as.data.frame(out)
    rownames(out) <- out$variable
    out$variable <- NULL
    return(out)
  }

  full_summary <- .summary(draws, variables, probs, robust)
  out$has_rhat <- "Rhat" %in% colnames(full_summary)
  if (algorithm(object) == "sampling") {
    if (out$has_rhat) {
      Rhats <- full_summary[, "Rhat"]
      if (any(Rhats > 1.05, na.rm = TRUE)) {
        warning2(
          "Parts of the model have not converged (some Rhats are > 1.05). ",
          "Be careful when analysing the results! We recommend running ",
          "more iterations and/or setting stronger priors."
        )
      }
    }
    div_trans <- sum(nuts_params(object, pars = "divergent__")$Value)
    adapt_delta <- control_params(object)$adapt_delta
    if (div_trans > 0) {
      warning2(
        "There were ", div_trans, " divergent transitions after warmup. ",
        "Increasing adapt_delta above ", adapt_delta, " may help. See ",
        "http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup"
      )
    }
  }

  # summary of population-level effects
  fe_pars <- variables[grepl(fixef_pars(), variables)]
  out$fixed <- full_summary[fe_pars, , drop = FALSE]
  rownames(out$fixed) <- gsub(fixef_pars(), "", fe_pars)

  # summary of family specific parameters
  spec_pars <- c(valid_dpars(object), "delta")
  spec_pars <- paste0(spec_pars, collapse = "|")
  spec_pars <- paste0("^(", spec_pars, ")($|_)")
  spec_pars <- variables[grepl(spec_pars, variables)]
  out$spec_pars <- full_summary[spec_pars, , drop = FALSE]
  # correlation parameters require renaming to look good in the summary
  lncor_pars <- variables[grepl("^lncor_", variables)]
  if (length(lncor_pars)) {
    lncor_summary <- full_summary[lncor_pars, , drop = FALSE]
    lncor_pars <- sub("__", ",", sub("__", "(", lncor_pars))
    rownames(lncor_summary) <- paste0(lncor_pars, ")")
    out$spec_pars <- rbind(out$spec_pars, lncor_summary)
  }

  # summary of residual correlations
  rescor_pars <- variables[grepl("^rescor_", variables)]
  if (length(rescor_pars)) {
    out$rescor_pars <- full_summary[rescor_pars, , drop = FALSE]
    rescor_pars <- sub("__", ",", sub("__", "(", rescor_pars))
    rownames(out$rescor_pars) <- paste0(rescor_pars, ")")
  }

  # summary of autocorrelation effects
  cor_pars <- variables[grepl(regex_autocor_pars(), variables)]
  out$cor_pars <- full_summary[cor_pars, , drop = FALSE]
  rownames(out$cor_pars) <- cor_pars
  cortime_pars <- variables[grepl("^cortime_", variables)]
  if (length(cortime_pars)) {
    tmp <- full_summary[cortime_pars, , drop = FALSE]
    cortime_pars <- sub("__", ",", sub("__", "(", cortime_pars))
    rownames(tmp) <- paste0(cortime_pars, ")")
    out$cor_pars <- rbind(out$cor_pars, tmp)
  }

  # summary of group-level effects
  for (g in out$group) {
    gregex <- escape_dot(g)
    sd_prefix <- paste0("^sd_", gregex, "__")
    sd_pars <- variables[grepl(sd_prefix, variables)]
    cor_prefix <- paste0("^cor_", gregex, "__")
    cor_pars <- variables[grepl(cor_prefix, variables)]
    df_prefix <- paste0("^df_", gregex, "$")
    df_pars <- variables[grepl(df_prefix, variables)]
    gpars <- c(df_pars, sd_pars, cor_pars)
    out$random[[g]] <- full_summary[gpars, , drop = FALSE]
    if (has_rows(out$random[[g]])) {
      sd_names <- sub(sd_prefix, "sd(", sd_pars)
      cor_names <- sub(cor_prefix, "cor(", cor_pars)
      cor_names <- sub("__", ",", cor_names)
      df_names <- sub(df_prefix, "df", df_pars)
      gnames <- c(df_names, paste0(c(sd_names, cor_names), ")"))
      rownames(out$random[[g]]) <- gnames
    }
  }
  # summary of smooths
  sm_pars <- variables[grepl("^sds_", variables)]
  if (length(sm_pars)) {
    out$splines <- full_summary[sm_pars, , drop = FALSE]
    rownames(out$splines) <- paste0(gsub("^sds_", "sds(", sm_pars), ")")
  }
  # summary of monotonic parameters
  mo_pars <- variables[grepl("^simo_", variables)]
  if (length(mo_pars)) {
    out$mo <- full_summary[mo_pars, , drop = FALSE]
    rownames(out$mo) <- gsub("^simo_", "", mo_pars)
  }
  # summary of gaussian processes
  gp_pars <- variables[grepl("^(sdgp|lscale)_", variables)]
  if (length(gp_pars)) {
    out$gp <- full_summary[gp_pars, , drop = FALSE]
    rownames(out$gp) <- gsub("^sdgp_", "sdgp(", rownames(out$gp))
    rownames(out$gp) <- gsub("^lscale_", "lscale(", rownames(out$gp))
    rownames(out$gp) <- paste0(rownames(out$gp), ")")
  }
  out
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
#' @seealso \code{\link{summary.brmsfit}}
#'
#' @export
print.brmsfit <- function(x, digits = 2, ...) {
  print(summary(x, ...), digits = digits, ...)
}

#' @export
print.brmssummary <- function(x, digits = 2, ...) {
  cat(" Family: ")
  cat(summarise_families(x$formula), "\n")
  cat("  Links: ")
  cat(summarise_links(x$formula, wsp = 9), "\n")
  cat("Formula: ")
  print(x$formula, wsp = 9)
  cat(paste0(
    "   Data: ", x$data_name,
    " (Number of observations: ", x$nobs, ") \n"
  ))
  if (!isTRUE(nzchar(x$sampler))) {
    cat("\nThe model does not contain posterior draws.\n")
    return(invisible(x))
  }
  # TODO: make this option a user-facing argument?
  short <- as_one_logical(getOption("brms.short_summary", FALSE))
  if (!short) {
    cat(paste0(
      "  Draws: ", x$chains, " chains, each with iter = ", x$iter,
      "; warmup = ", x$warmup, "; thin = ", x$thin, ";\n",
      "         total post-warmup draws = ", x$total_ndraws, "\n"
    ))
  }
  cat("\n")
  # TODO: change order of the displayed summaries?
  if (nrow(x$prior)) {
    cat("Priors:\n")
    print(x$prior, show_df = FALSE)
    cat("\n")
  }
  if (length(x$splines)) {
    cat("Smoothing Spline Hyperparameters:\n")
    print_format(x$splines, digits)
    cat("\n")
  }
  if (length(x$gp)) {
    cat("Gaussian Process Hyperparameters:\n")
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
    cat("Multilevel Hyperparameters:\n")
    for (i in seq_along(x$random)) {
      g <- names(x$random)[i]
      cat(paste0("~", g, " (Number of levels: ", x$ngrps[[g]], ") \n"))
      print_format(x$random[[g]], digits)
      cat("\n")
    }
  }
  if (nrow(x$fixed)) {
    cat("Regression Coefficients:\n")
    print_format(x$fixed, digits)
    cat("\n")
  }
  if (length(x$mo)) {
    cat("Monotonic Simplex Parameters:\n")
    print_format(x$mo, digits)
    cat("\n")
  }
  if (nrow(x$spec_pars)) {
    cat("Further Distributional Parameters:\n")
    print_format(x$spec_pars, digits)
    cat("\n")
  }
  if (length(x$rescor_pars)) {
    cat("Residual Correlations: \n")
    print_format(x$rescor, digits)
    cat("\n")
  }
  if (!short) {
    cat(paste0("Draws were sampled using ", x$sampler, ". "))
    if (x$algorithm == "sampling") {
      if (isTRUE(x$has_rhat)) {
        cat(paste0(
          "For each parameter, Bulk_ESS\n",
          "and Tail_ESS are effective sample size measures, ",
          "and Rhat is the potential\n",
          "scale reduction factor on split chains ",
          "(at convergence, Rhat = 1)."
        ))
      } else {
        cat(paste0(
          "Standard Rhat and ESS estimates\n",
          "should not be trusted for brm_multiple models and are hence not displayed.\n",
          "Please see ?brm_multiple for how to assess convergence of such models."
        ))
      }
    }
    cat("\n")
  }
  invisible(x)
}

# helper function to print summary matrices in nice format
# also displays -0.00 as a result of round negative values to zero (#263)
# @param x object to be printed; coerced to matrix
# @param digits number of digits to show
# @param no_digits names of columns for which no digits should be shown
print_format <- function(x, digits = 2, no_digits = c("Bulk_ESS", "Tail_ESS")) {
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

# regex to extract population-level coefficients
fixef_pars <- function() {
  types <- c("", "s", "cs", "sp", "mo", "me", "mi", "m")
  types <- paste0("(", types, ")", collapse = "|")
  paste0("^b(", types, ")_")
}

# algorithm used in the model fitting
algorithm <- function(x) {
  stopifnot(is.brmsfit(x))
  if (is.null(x$algorithm)) "sampling"
  else x$algorithm
}

#' Summarize Posterior draws
#'
#' Summarizes posterior draws based on point estimates (mean or median),
#' estimation errors (SD or MAD) and quantiles. This function mainly exists to
#' retain backwards compatibility. It will eventually be replaced by functions
#' of the \pkg{posterior} package (see examples below).
#'
#' @param x An \R object.
#' @inheritParams as.matrix.brmsfit
#' @param probs The percentiles to be computed by the
#'   \code{\link[stats:quantile]{quantile}} function.
#' @param robust If \code{FALSE} (the default) the mean is used as
#'  the measure of central tendency and the standard deviation as
#'  the measure of variability. If \code{TRUE}, the median and the
#'  median absolute deviation (MAD) are applied instead.
#' @param ... More arguments passed to or from other methods.
#'
#' @return A matrix where rows indicate variables
#' and columns indicate the summary estimates.
#'
#' @seealso \code{\link[posterior:summarize_draws]{summarize_draws}}
#'
#' @examples
#' \dontrun{
#' fit <- brm(time ~ age * sex, data = kidney)
#' posterior_summary(fit)
#'
#' # recommended workflow using posterior
#' library(posterior)
#' draws <- as_draws_array(fit)
#' summarise_draws(draws, default_summary_measures())
#' }
#'
#' @export
posterior_summary <- function(x, ...) {
  UseMethod("posterior_summary")
}

#' @rdname posterior_summary
#' @export
posterior_summary.default <- function(x, probs = c(0.025, 0.975),
                                      robust = FALSE, ...) {
  # TODO: replace with summary functions from posterior
  # TODO: find a way to represent 3D summaries as well
  if (!length(x)) {
    stop2("No posterior draws supplied.")
  }
  if (robust) {
    coefs <- c("median", "mad", "quantile")
  } else {
    coefs <- c("mean", "sd", "quantile")
  }
  .posterior_summary <- function(x) {
    do_call(cbind, lapply(
      coefs, get_estimate, draws = x,
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
  # TODO: align names with summary outputs of other methods and packages
  colnames(out) <- c("Estimate", "Est.Error", paste0("Q", probs * 100))
  out
}

#' @rdname posterior_summary
#' @export
posterior_summary.brmsfit <- function(x, pars = NA, variable = NULL,
                                      probs = c(0.025, 0.975),
                                      robust = FALSE, ...) {
  out <- as.matrix(x, pars = pars, variable = variable, ...)
  posterior_summary(out, probs = probs, robust = robust, ...)
}

# calculate estimates over posterior draws
# @param coef coefficient to be applied on the draws (e.g., "mean")
# @param draws the draws over which to apply coef
# @param margin see 'apply'
# @param ... additional arguments passed to get(coef)
# @return typically a matrix with colnames(draws) as colnames
get_estimate <- function(coef, draws, margin = 2, ...) {
  # TODO: replace with summary functions from posterior
  dots <- list(...)
  args <- list(X = draws, MARGIN = margin, FUN = coef)
  fun_args <- names(formals(coef))
  if (!"..." %in% fun_args) {
    dots <- dots[names(dots) %in% fun_args]
  }
  x <- do_call(apply, c(args, dots))
  if (is.null(dim(x))) {
    x <- matrix(x, dimnames = list(NULL, coef))
  } else if (coef == "quantile") {
    x <- aperm(x, length(dim(x)):1)
  }
  x
}

# validate bounds of credible intervals
# @return a numeric vector of length 2
validate_ci_bounds <- function(prob, probs = NULL) {
  if (!is.null(probs)) {
    # deprecated as of version 2.13.7
    warning2("Argument 'probs' is deprecated. Please use 'prob' instead.")
    if (length(probs) != 2L) {
      stop2("Arguments 'probs' must be of length 2.")
    }
    probs <- as.numeric(probs)
  } else {
    prob <- as_one_numeric(prob)
    if (prob < 0 || prob > 1) {
      stop2("'prob' must be a single numeric value in [0, 1].")
    }
    probs <- c((1 - prob) / 2, 1 - (1 - prob) / 2)
  }
  probs
}

#' Table Creation for Posterior Draws
#'
#' Create a table for unique values of posterior draws.
#' This is usually only useful when summarizing predictions
#' of ordinal models.
#'
#' @param x A matrix of posterior draws where rows
#'   indicate draws and columns indicate parameters.
#' @param levels Optional values of possible posterior values.
#'   Defaults to all unique values in \code{x}.
#'
#' @return A matrix where rows indicate parameters
#'  and columns indicate the unique values of
#'  posterior draws.
#'
#' @examples
#' \dontrun{
#' fit <- brm(rating ~ period + carry + treat,
#'            data = inhaler, family = cumulative())
#' pr <- predict(fit, summary = FALSE)
#' posterior_table(pr)
#' }
#'
#' @export
posterior_table <- function(x, levels = NULL) {
  x <- as.matrix(x)
  if (anyNA(x)) {
    warning2("NAs will be ignored in 'posterior_table'.")
  }
  if (is.null(levels)) {
    levels <- sort(unique(as.vector(x)))
  }
  xlevels <- attr(x, "levels")
  if (length(xlevels) != length(levels)) {
    xlevels <- levels
  }
  out <- lapply(seq_len(ncol(x)),
    function(n) table(factor(x[, n], levels = levels))
  )
  out <- do_call(rbind, out)
  # compute relative frequencies
  out <- out / rowSums(out)
  rownames(out) <- colnames(x)
  colnames(out) <- paste0("P(Y = ", xlevels, ")")
  out
}

#' Compute posterior uncertainty intervals
#'
#' Compute posterior uncertainty intervals for \code{brmsfit} objects.
#'
#' @param object An object of class \code{brmsfit}.
#' @param prob A value between 0 and 1 indicating the desired probability
#'   to be covered by the uncertainty intervals. The default is 0.95.
#' @inheritParams as.matrix.brmsfit
#' @param ... More arguments passed to \code{\link{as.matrix.brmsfit}}.
#'
#' @return A \code{matrix} with lower and upper interval bounds
#'   as columns and as many rows as selected variables.
#'
#' @examples
#' \dontrun{
#' fit <- brm(count ~ zAge + zBase * Trt,
#'            data = epilepsy, family = negbinomial())
#' posterior_interval(fit)
#' }
#'
#' @aliases posterior_interval
#' @method posterior_interval brmsfit
#' @export
#' @export posterior_interval
#' @importFrom rstantools posterior_interval
posterior_interval.brmsfit <- function(
  object, pars = NA, variable = NULL, prob = 0.95, ...
) {
  ps <- as.matrix(object, pars = pars, variable = variable, ...)
  rstantools::posterior_interval(ps, prob = prob)
}
