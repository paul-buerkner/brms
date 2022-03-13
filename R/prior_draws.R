#' Extract Prior Draws
#'
#' Extract prior draws of specified parameters
#'
#' @aliases prior_draws.brmsfit prior_samples
#'
#' @param x An \code{R} object typically of class \code{brmsfit}.
#' @inheritParams as.data.frame.brmsfit
#' @param ... Arguments passed to individual methods (if applicable).
#'
#' @details To make use of this function, the model must contain draws of
#'   prior distributions. This can be ensured by setting \code{sample_prior =
#'   TRUE} in function \code{brm}. Priors of certain parameters cannot be saved
#'   for technical reasons. For instance, this is the case for the
#'   population-level intercept, which is only computed after fitting the model
#'   by default. If you want to treat the intercept as part of all the other
#'   regression coefficients, so that sampling from its prior becomes possible,
#'   use \code{... ~ 0 + Intercept + ...} in the formulas.
#'
#' @return A \code{data.frame} containing the prior draws.
#'
#' @examples
#' \dontrun{
#' fit <- brm(rating ~ treat + period + carry + (1|subject),
#'            data = inhaler, family = "cumulative",
#'            prior = set_prior("normal(0,2)", class = "b"),
#'            sample_prior = TRUE)
#'
#' # extract all prior draws
#' draws1 <- prior_draws(fit)
#' head(draws1)
#'
#' # extract prior draws for the coefficient of 'treat'
#' draws2 <- prior_draws(fit, "b_treat")
#' head(draws2)
#' }
#'
#' @export
prior_draws.brmsfit <- function(x, variable = NULL, pars = NULL, ...) {
  variable <- use_alias(variable, pars)
  if (!is.null(variable)) {
    variable <- as.character(variable)
  }
  all_names <- variables(x)
  prior_names <- unique(all_names[grepl("^prior_", all_names)])
  if (!length(prior_names)) {
    return(data.frame(NULL))
  }
  draws <- as.data.frame(x, variable = prior_names)
  names(draws) <- sub("^prior_", "", prior_names)
  if (is.null(variable)) {
    return(draws)
  }
  # get prior draws for a single variable
  .prior_draws <- function(variable) {
    matches <- paste0("^", escape_all(names(draws)))
    matches <- lapply(matches, regexpr, text = variable)
    matches <- ulapply(matches, attr, which = "match.length")
    if (max(matches) == -1 || ignore_prior(x, variable)) {
      out <- NULL
    } else {
      take <- match(max(matches), matches)
      # order draws randomly to avoid artificial dependencies
      # between parameters using the same prior draws
      draws <- list(draws[sample(ndraws(x)), take])
      out <- structure(draws, names = variable)
    }
    return(out)
  }
  draws <- rmNULL(lapply(variable, .prior_draws))
  draws <- data.frame(draws, check.names = FALSE)
  draws
}

#' @rdname prior_draws.brmsfit
#' @export
prior_draws <- function(x, ...) {
  UseMethod("prior_draws")
}

#' @export
prior_draws.default <- function(x, variable = NULL, pars = NULL,
                                regex = FALSE, fixed = FALSE, ...) {
  call <- match.call()
  if ("pars" %in% names(call)) {
    variable <- use_alias(variable, pars)
    regex <- !as_one_logical(fixed)
  }
  if (is.null(variable)) {
    variable <- "^prior_"
    regex <- TRUE
  } else {
    variable <- as.character(variable)
    regex <- as_one_logical(regex)
    if (regex) {
      hat <- substr(variable, 1, 1) == "^"
      variable <- ifelse(hat, substr(variable, 2, nchar(variable)), variable)
      variable <- paste0("^prior_", variable)
    } else {
      variable <- paste0("prior_", variable)
    }
  }
  x <- as_draws_df(as.data.frame(x))
  if (!regex) {
    # missing variables will leads to an error in posterior
    variable <- intersect(variable, variables(x))
    if (!length(variable)) {
      return(data.frame(NULL))
    }
  }
  x <- subset_draws(x, variable = variable, regex = regex, ...)
  unclass_draws(x)
}

#' @rdname prior_draws.brmsfit
#' @export
prior_samples <- function(x, ...) {
  warning2("'prior_samples' is deprecated. Please use 'prior_draws' instead.")
  UseMethod("prior_draws")
}

# ignore priors of certain parameters from whom we cannot obtain prior draws
# currently applies only to overall intercepts of centered design matrices
# fixes issue #696
# @param x a brmsfit object
# @param variable name of a single variable
# @return TRUE (if the prior should be ignored) or FALSE
ignore_prior <- function(x, variable) {
  stopifnot(is.brmsfit(x))
  variable <- as_one_character(variable)
  out <- FALSE
  if (grepl("^b_.*Intercept($|\\[)", variable)) {
    # cannot sample from intercepts if 'center' was TRUE
    intercept_priors <- subset2(x$prior, class = "Intercept")
    if (NROW(intercept_priors)) {
      # prefixes of the model intercepts
      p_intercepts <- usc(combine_prefix(intercept_priors))
      # prefix of the parameter under question
      p_par <- sub("^b", "", variable)
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
