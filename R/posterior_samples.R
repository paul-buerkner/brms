#' (Deprecated) Extract Posterior Samples
#'
#' Extract posterior samples of specified parameters. The
#' \code{posterior_samples} method is deprecated. We recommend using the more
#' modern and consistent \code{\link[brms:draws-brms]{as_draws_*}} extractor
#' functions of the \pkg{posterior} package instead.
#'
#' @param x An \code{R} object typically of class \code{brmsfit}
#' @param pars Names of parameters for which posterior samples
#'   should be returned, as given by a character vector or regular expressions.
#'   By default, all posterior samples of all parameters are extracted.
#' @param fixed Indicates whether parameter names
#'   should be matched exactly (\code{TRUE}) or treated as
#'   regular expressions (\code{FALSE}). Default is \code{FALSE}.
#' @param add_chain A flag indicating if the returned \code{data.frame}
#'   should contain two additional columns. The \code{chain} column
#'   indicates the chain in which each sample was generated, the \code{iter}
#'   column indicates the iteration number within each chain.
#' @param subset A numeric vector indicating the rows
#'   (i.e., posterior samples) to be returned.
#'   If \code{NULL} (the default), all  posterior samples are returned.
#' @param as.matrix Should the output be a \code{matrix}
#'   instead of a \code{data.frame}? Defaults to \code{FALSE}.
#' @param as.array Should the output be an \code{array}
#'   instead of a \code{data.frame}? Defaults to \code{FALSE}.
#' @param ... Arguments passed to individual methods (if applicable).
#'
#' @return A data.frame (matrix or array) containing the posterior samples.
#'
#' @seealso \code{\link[brms:draws-brms]{as_draws}},
#'    \code{\link[brms:as.data.frame.brmsfit]{as.data.frame}}
#'
#' @examples
#' \dontrun{
#' fit <- brm(rating ~ treat + period + carry + (1|subject),
#'            data = inhaler, family = "cumulative")
#'
#' # extract posterior samples of population-level effects
#' samples1 <- posterior_samples(fit, pars = "^b")
#' head(samples1)
#'
#' # extract posterior samples of group-level standard deviations
#' samples2 <- posterior_samples(fit, pars = "^sd_")
#' head(samples2)
#' }
#'
#' @export
posterior_samples.brmsfit <- function(x, pars = NA, fixed = FALSE,
                                      add_chain = FALSE, subset = NULL,
                                      as.matrix = FALSE, as.array = FALSE,
                                      ...) {
  if (as.matrix && as.array) {
    stop2("Cannot use 'as.matrix' and 'as.array' at the same time.")
  }
  if (add_chain && as.array) {
    stop2("Cannot use 'add_chain' and 'as.array' at the same time.")
  }
  contains_draws(x)
  pars <- extract_pars(pars, variables(x), fixed = fixed, ...)

  # get basic information on the samples
  iter <- x$fit@sim$iter
  warmup <- x$fit@sim$warmup
  thin <- x$fit@sim$thin
  chains <- x$fit@sim$chains
  final_iter <- ceiling((iter - warmup) / thin)
  samples_taken <- seq(warmup + 1, iter, thin)

  samples <- NULL
  if (length(pars)) {
    if (as.matrix) {
      samples <- as.matrix(x$fit, pars = pars)
    } else if (as.array) {
      samples <- as.array(x$fit, pars = pars)
    } else {
      samples <- as.data.frame(x$fit, pars = pars)
    }
    if (add_chain) {
      # name the column 'chain' not 'chains' (#32)
      samples <- cbind(samples,
        chain = factor(rep(1:chains, each = final_iter)),
        iter = rep(samples_taken, chains)
      )
    }
    if (!is.null(subset)) {
      if (as.array) {
        samples <- samples[subset, , , drop = FALSE]
      } else {
        samples <- samples[subset, , drop = FALSE]
      }
    }
  }
  samples
}

#' @rdname posterior_samples.brmsfit
#' @export
posterior_samples <- function(x, pars = NA, ...) {
  warning2("Method 'posterior_samples' is deprecated. ",
           "Please see ?as_draws for recommended alternatives.")
  UseMethod("posterior_samples")
}

#' @export
posterior_samples.default <- function(x, pars = NA, fixed = FALSE, ...) {
  x <- as.data.frame(x)
  if (!anyNA(pars)) {
    pars <- extract_pars(pars, all_pars = names(x), fixed = fixed, ...)
    x <- x[, pars, drop = FALSE]
  }
  if (!ncol(x)) {
    x <- NULL
  }
  x
}

#' Extract Parameter Names
#'
#' Extract all parameter names of a given model.
#'
#' @aliases parnames.brmsfit
#'
#' @param x An \R object
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A character vector containing the parameter names of the model.
#'
#' @export
parnames <- function(x, ...) {
  warning2("'parnames' is deprecated. Please use 'variables' instead.")
  UseMethod("parnames")
}

#' @export
parnames.default <- function(x, ...) {
  names(x)
}

#' @export
parnames.brmsfit <- function(x, ...) {
  out <- dimnames(x$fit)
  if (is.list(out)) {
    out <- out$parameters
  }
  out
}

# extract all valid parameter names that match pars
# @param pars A character vector or regular expression
# @param all_pars all parameter names of the fitted model
# @param fixed should parameter names be matched exactly?
# @param exact_match deprecated alias of fixed
# @param na_value: what should be returned if pars is NA?
# @param ... Further arguments to be passed to grepl
# @return A character vector of parameter names
extract_pars <- function(pars, all_pars, fixed = FALSE,
                         exact_match = FALSE,
                         na_value = all_pars, ...) {
  if (!(anyNA(pars) || is.character(pars))) {
    stop2("Argument 'pars' must be NA or a character vector.")
  }
  fixed <- check_deprecated_fixed(fixed, exact_match)
  if (!anyNA(pars)) {
    fixed <- as_one_logical(fixed)
    if (fixed) {
      out <- intersect(pars, all_pars)
    } else {
      out <- vector("list", length(pars))
      for (i in seq_along(pars)) {
        out[[i]] <- all_pars[grepl(pars[i], all_pars, ...)]
      }
      out <- unique(unlist(out))
    }
  } else {
    out <- na_value
  }
  out
}

# check deprecated alias of argument 'fixed'
check_deprecated_fixed <- function(fixed, exact_match) {
  if (!isFALSE(exact_match)) {
    # deprecated as of brms 2.10.6; remove in brms 3.0
    warning2("Argument 'exact_match' is deprecated. ",
             "Please use 'fixed' instead.")
    fixed <- exact_match
  }
  fixed
}

#' Extract posterior samples for use with the \pkg{coda} package
#'
#' @aliases as.mcmc
#'
#' @inheritParams posterior_samples.brmsfit
#' @param ... currently unused
#' @param combine_chains Indicates whether chains should be combined.
#' @param inc_warmup Indicates if the warmup samples should be included.
#'   Default is \code{FALSE}. Warmup samples are used to tune the
#'   parameters of the sampling algorithm and should not be analyzed.
#'
#' @return If \code{combine_chains = TRUE} an \code{mcmc} object is returned.
#'   If \code{combine_chains = FALSE} an \code{mcmc.list} object is returned.
#'
#' @method as.mcmc brmsfit
#' @export
#' @export as.mcmc
#' @importFrom coda as.mcmc
as.mcmc.brmsfit <- function(x, pars = NA, fixed = FALSE,
                            combine_chains = FALSE, inc_warmup = FALSE,
                            ...) {
  warning2("as.mcmc.brmsfit is deprecated and will eventually be removed.")
  contains_draws(x)
  pars <- extract_pars(pars, all_pars = variables(x), fixed = fixed, ...)
  combine_chains <- as_one_logical(combine_chains)
  inc_warmup <- as_one_logical(inc_warmup)
  if (combine_chains) {
    if (inc_warmup) {
      stop2("Cannot include warmup samples when 'combine_chains' is TRUE.")
    }
    out <- as.matrix(x$fit, pars)
    ndraws <- nrow(out)
    end <- x$fit@sim$iter * x$fit@sim$chains
    thin <- x$fit@sim$thin
    start <- end - (ndraws - 1) * thin
    mcpar <- c(start, end, thin)
    attr(out, "mcpar") <- mcpar
    class(out) <- "mcmc"
  } else {
    thin <- x$fit@sim$thin
    if (inc_warmup && thin >= 2) {
      stop2("Cannot include warmup samples when 'thin' >= 2.")
    }
    ps <- rstan::extract(x$fit, pars, permuted = FALSE, inc_warmup = inc_warmup)
    ndraws <- dim(ps)[1]
    end <- x$fit@sim$iter
    start <- end - (ndraws - 1) * thin
    mcpar <- c(start, end, thin)
    out <- vector("list", length = dim(ps)[2])
    for (i in seq_along(out)) {
      out[[i]] <- ps[, i, ]
      attr(out[[i]], "mcpar") <- mcpar
      class(out[[i]]) <- "mcmc"
    }
    class(out) <- "mcmc.list"
  }
  out
}
