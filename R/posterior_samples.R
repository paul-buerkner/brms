#' Extract Posterior Samples
#' 
#' Extract posterior samples of specified parameters.
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
#' @param row.names,optional See \code{\link[base:as.data.frame]{as.data.frame}}.
#' @param ... For \code{as.data.frame}, \code{as.matrix}, and \code{as.array}:
#'   Further arguments to be passed to \code{posterior_samples}.
#'   
#' @details Currently there are methods for \code{brmsfit} objects.
#'   \code{as.data.frame.brmsfit}, \code{as.matrix.brmsfit}, and
#'   \code{as.array.brmsfit} are basically aliases of 
#'   \code{posterior_samples.brmsfit} and differ from
#'   each other only in type of the returned object.
#'   
#' @return A data.frame (matrix or array) containing the posterior samples, 
#'   with one column per parameter. In case an array is returned,
#'   it contains one additional dimension for the chains.
#' 
#' @examples
#' \dontrun{
#' fit <- brm(rating ~ treat + period + carry + (1|subject), 
#'            data = inhaler, family = "cumulative")
#' 
#' # extract posterior samples of population-level effects 
#' samples1 <- posterior_samples(fit, "^b")
#' head(samples1)
#' 
#' # extract posterior samples of group-level standard deviations
#' samples2 <- posterior_samples(fit, "^sd_")
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
  contains_samples(x)
  pars <- extract_pars(pars, parnames(x), fixed = fixed, ...)
  
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

#' @rdname posterior_samples.brmsfit
#' @export
as.data.frame.brmsfit <- function(x, row.names = NULL, optional = TRUE, ...) {
  out <- posterior_samples(x, ..., as.matrix = FALSE)
  data.frame(out, row.names = row.names, check.names = !optional)
}

#' @rdname posterior_samples.brmsfit
#' @export
as.matrix.brmsfit <- function(x, ...) {
  posterior_samples(x, ..., as.matrix = TRUE)
}

#' @rdname posterior_samples.brmsfit
#' @export
as.array.brmsfit <- function(x, ...) {
  posterior_samples(x, ..., as.array = TRUE)
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
# @param fixed should parnames be matched exactly?
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

#' Extract prior samples
#' 
#' Extract prior samples of specified parameters 
#' 
#' @aliases prior_samples.brmsfit
#' 
#' @param x An \code{R} object typically of class \code{brmsfit}
#' @param pars Names of parameters for which prior samples should be returned, 
#'   as given by a character vector or regular expressions.
#'   By default, all prior samples are extracted
#' @param ... Currently ignored
#'   
#' @details To make use of this function, the model must contain samples of
#'   prior distributions. This can be ensured by setting \code{sample_prior =
#'   TRUE} in function \code{brm}. Priors of certain parameters cannot be saved
#'   for technical reasons. For instance, this is the case for the
#'   population-level intercept, which is only computed after fitting the model
#'   by default. If you want to treat the intercept as part of all the other
#'   regression coefficients, so that sampling from its prior becomes possible,
#'   use \code{... ~ 0 + Intercept + ...} in the formulas.
#'   
#' @return A data frame containing the prior samples.
#' 
#' @examples
#' \dontrun{
#' fit <- brm(rating ~ treat + period + carry + (1|subject), 
#'            data = inhaler, family = "cumulative", 
#'            prior = set_prior("normal(0,2)", class = "b"), 
#'            sample_prior = TRUE)
#' 
#' # extract all prior samples
#' samples1 <- prior_samples(fit)
#' head(samples1)
#' 
#' # extract prior samples for the population-level effects of 'treat'
#' samples2 <- prior_samples(fit, "b_treat")
#' head(samples2)
#' }
#' 
#' @export
prior_samples.brmsfit <- function(x, pars = NA, ...) {
  if (!anyNA(pars) && !is.character(pars)) {
    stop2("Argument 'pars' must be a character vector.")
  }
  par_names <- parnames(x)
  prior_names <- unique(par_names[grepl("^prior_", par_names)])
  if (length(prior_names)) {
    samples <- posterior_samples(x, prior_names, fixed = TRUE)
    names(samples) <- sub("^prior_", "", prior_names)
    if (!anyNA(pars)) {
      .prior_samples <- function(par) {
        # get prior samples for parameter par 
        matches <- paste0("^", escape_all(names(samples)))
        matches <- lapply(matches, regexpr, text = par)
        matches <- ulapply(matches, attr, which = "match.length")
        if (max(matches) == -1 || ignore_prior(x, par)) {
          out <- NULL
        } else {
          take <- match(max(matches), matches)
          # order samples randomly to avoid artifical dependencies
          # between parameters using the same prior samples
          samples <- list(samples[sample(nsamples(x)), take])
          out <- structure(samples, names = par)
        }
        return(out)
      }
      samples <- data.frame(
        rmNULL(lapply(pars, .prior_samples)), check.names = FALSE
      )
    }
  } else {
    samples <- NULL
  }
  samples
}

#' @rdname prior_samples.brmsfit
#' @export 
prior_samples <- function(x, pars = NA, ...) {
  UseMethod("prior_samples")
}

#' @export
prior_samples.default <- function(x, pars = NA, fixed = FALSE, ...) {
  if (anyNA(pars)) {
    pars <- "^prior_"
    fixed <- FALSE
  } else {
    if (fixed) {
      pars <- paste0("prior_", pars) 
    } else {
      hat <- substr(pars, 1, 1) == "^"
      pars <- ifelse(hat, substr(pars, 2, nchar(pars)), pars)
      pars <- paste0("^prior_", pars)  
    }
  }
  posterior_samples(x, pars = pars, fixed = fixed, ...)
}

# ignore priors of certain parameters from whom we cannot obtain prior samples
# currently applies only to overall intercepts of centered design matrices
# fixes issue #696
# @param x a brmsfit object
# @param par name of the parameter
# @return TRUE (if the prior should be ignored) or FALSE
ignore_prior <- function(x, par) {
  stopifnot(is.brmsfit(x))
  par <- as_one_character(par)
  out <- FALSE
  if (grepl("^b_.*Intercept($|\\[)", par)) {
    # cannot sample from intercepts if 'center' was TRUE
    intercept_priors <- subset2(x$prior, class = "Intercept")
    if (NROW(intercept_priors)) {
      # prefixes of the model intercepts
      p_intercepts <- usc(combine_prefix(intercept_priors))
      # prefix of the parameter under question
      p_par <- sub("^b", "", par)
      p_par <- sub("_Intercept($|\\[)", "", p_par)
      out <- p_par %in% p_intercepts
      if (out) {
        warning2(
          "Sampling from the prior of an overall intercept is not ", 
          "possible by default. See the documentation of the ", 
          "'sample_prior' argument in help('brm')."
        )
      }
    }
  }
  out
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
  contains_samples(x)
  pars <- extract_pars(pars, all_pars = parnames(x), fixed = fixed, ...)
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
