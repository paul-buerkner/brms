#' Index \code{brmsfit} objects
#'
#' @aliases variables nvariables niterations nchains ndraws
#'
#' Index variables, iterations, chains, and draws.
#'
#' @param x A \code{brmsfit} object or another \R object for which
#' the methods are defined.
#' @param ... Arguments passed to individual methods (if applicable).
#'
#' @name draws-index-brms
NULL

#' @rdname draws-index-brms
#' @importFrom posterior variables
#' @method variables brmsfit
#' @export
#' @export variables
variables.brmsfit <- function(x, ...) {
  # TODO: simplify once rstan and cmdstanr support these methods
  out <- dimnames(x$fit)
  if (is.list(out)) {
    out <- out$parameters
  }
  out
}

#' @method variables data.frame
variables.data.frame <- function(x, ...) {
  names(x)
}

#' @rdname draws-index-brms
#' @importFrom posterior nvariables
#' @method nvariables brmsfit
#' @export
#' @export nvariables
nvariables.brmsfit <- function(x, ...) {
  length(variables(x, ...))
}

#' @rdname draws-index-brms
#' @importFrom posterior niterations
#' @method niterations brmsfit
#' @export
#' @export niterations
niterations.brmsfit <- function(x) {
  if (!is.stanfit(x$fit)) return(0)
  niterations <- x$fit@sim$n_save[1] %||% 0
  niterations - nwarmup(x)
}

#' @rdname draws-index-brms
#' @importFrom posterior nchains
#' @method nchains brmsfit
#' @export
#' @export nchains
nchains.brmsfit <- function(x) {
  if (!is.stanfit(x$fit)) return(0)
  x$fit@sim$chains %||% 0
}

#' @rdname draws-index-brms
#' @importFrom posterior ndraws
#' @method ndraws brmsfit
#' @export
#' @export ndraws
ndraws.brmsfit <- function(x) {
  niterations(x) * nchains(x)
}

nwarmup <- function(x) {
  if (!is.stanfit(x$fit)) return(0)
  x$fit@sim$warmup2[1] %||% 0
}

nthin <- function(x) {
  if (!is.stanfit(x$fit)) return(1)
  x$fit@sim$thin %||% 1
}

#' Transform \code{brmsfit} to \code{draws} objects
#'
#' Transform a \code{brmsfit} object to a format supported by the
#' \pkg{posterior} package.
#'
#' @aliases as_draws as_draws_matrix as_draws_array as_draws_df
#' @aliases as_draws_rvars as_draws_list
#'
#' @param x A \code{brmsfit} object or another \R object for which
#' the methods are defined.
#' @param variable A character vector providing the variables to extract.
#'   By default, all variables are extracted.
#' @param regex Logical; Should variable should be treated as a (vector of)
#'   regular expressions? Any variable in \code{x} matching at least one of the
#'   regular expressions will be selected. Defaults to \code{FALSE}.
#' @param inc_warmup Should warmup draws be included? Defaults to \code{FALSE}.
#' @param ... Arguments passed to individual methods (if applicable).
#'
#' @details To subset iterations, chains, or draws, use the
#'   \code{\link[posterior:subset_draws]{subset_draws}} method after
#'   transforming the \code{brmsfit} to a \code{draws} object.
#'
#' @seealso \code{\link[posterior:draws]{draws}}
#'   \code{\link[posterior:subset_draws]{subset_draws}}
#'
#' @examples
#' \dontrun{
#' fit <- brm(count ~ zAge + zBase * Trt + (1|patient),
#'            data = epilepsy, family = poisson())
#'
#' # extract posterior draws in an array format
#' (draws_fit <- as_draws_array(fit))
#' posterior::summarize_draws(draws_fit)
#'
#' # extract only certain variables
#' as_draws_array(fit, variable = "r_patient")
#' as_draws_array(fit, variable = "^b_", regex = TRUE)
#'
#' # extract posterior draws in a random variables format
#' as_draws_rvars(fit)
#' }
#'
#' @name draws-brms
NULL

#' @rdname draws-brms
#' @importFrom posterior as_draws
#' @method as_draws brmsfit
#' @export
#' @export as_draws
as_draws.brmsfit <- function(x, variable = NULL, regex = FALSE,
                             inc_warmup = FALSE, ...) {
  # draws_list is the fastest format to convert to at the moment
  as_draws_list(
    x, variable = variable, regex = regex,
    inc_warmup = inc_warmup, ...
  )
}

#' @rdname draws-brms
#' @importFrom posterior as_draws_matrix
#' @method as_draws_matrix brmsfit
#' @export
#' @export as_draws_matrix
as_draws_matrix.brmsfit <- function(x, variable = NULL, regex = FALSE,
                                    inc_warmup = FALSE, ...) {
  as_draws_matrix(as_draws_list(
    x, variable = variable, regex = regex,
    inc_warmup = inc_warmup, ...
  ))
}

#' @rdname draws-brms
#' @importFrom posterior as_draws_array
#' @method as_draws_array brmsfit
#' @export
#' @export as_draws_array
as_draws_array.brmsfit <- function(x, variable = NULL, regex = FALSE,
                                   inc_warmup = FALSE, ...) {
  as_draws_array(as_draws_list(
    x, variable = variable, regex = regex,
    inc_warmup = inc_warmup, ...
  ))
}

#' @rdname draws-brms
#' @importFrom posterior as_draws_df
#' @method as_draws_df brmsfit
#' @export
#' @export as_draws_df
as_draws_df.brmsfit <- function(x, variable = NULL, regex = FALSE,
                                inc_warmup = FALSE, ...) {
  as_draws_df(as_draws_list(
    x, variable = variable, regex = regex,
    inc_warmup = inc_warmup, ...
  ))
}

#' @rdname draws-brms
#' @importFrom posterior as_draws_list
#' @method as_draws_list brmsfit
#' @export
#' @export as_draws_list
as_draws_list.brmsfit <- function(x, variable = NULL, regex = FALSE,
                                  inc_warmup = FALSE, ...) {
  .as_draws_list(
    x$fit, variable = variable, regex = regex,
    inc_warmup = inc_warmup, ...
  )
}

#' @rdname draws-brms
#' @importFrom posterior as_draws_rvars
#' @method as_draws_rvars brmsfit
#' @export
#' @export as_draws_rvars
as_draws_rvars.brmsfit <- function(x, variable = NULL, regex = FALSE,
                                   inc_warmup = FALSE, ...) {
  as_draws_rvars(as_draws_list(
    x, variable = variable, regex = regex,
    inc_warmup = inc_warmup, ...
  ))
}

# in stanfit objects draws are stored in a draws_list-like format
# so converting from there will be most efficient
# may be removed once rstan supports posterior natively
.as_draws_list <- function(x, variable = NULL, regex = FALSE,
                           inc_warmup = FALSE, ...) {
  stopifnot(is.stanfit(x))
  inc_warmup <- as_one_logical(inc_warmup)
  if (!length(x@sim$samples)) {
    stop2("The model does not contain posterior draws.")
  }
  out <- as_draws_list(x@sim$samples)
  # first subset variables then remove warmup as removing warmup
  # will take a lot of time when extracting many variables
  out <- subset_draws(out, variable = variable, regex = regex)
  if (!inc_warmup) {
    nwarmup <- x@sim$warmup2[1] %||% 0
    warmup_ids <- seq_len(nwarmup)
    iteration_ids <- posterior::iteration_ids(out)
    if (length(warmup_ids)) {
      iteration_ids <- iteration_ids[-warmup_ids]
    }
    out <- subset_draws(out, iteration = iteration_ids)
  }
  out
}

#' Extract Posterior Draws
#'
#' Extract posterior draws in conventional formats
#' as data.frames, matrices, or arrays.
#'
#' @inheritParams as_draws.brmsfit
#' @param pars Deprecated alias of \code{variable}. For reasons of backwards
#'   compatibility, \code{pars} is interpreted as a vector of regular
#'   expressions by default unless \code{fixed = TRUE} is specified.
#' @param draw The draw indices to be select. Subsetting draw indices will lead
#'   to an automatic merging of chains.
#' @param subset Deprecated alias of \code{draw}.
#' @param row.names,optional Unused and only added for consistency with
#'   the \code{\link[base:as.data.frame]{as.data.frame}} generic.
#' @param ... Further arguments to be passed to the corresponding
#'   \code{\link[brms:draws-brms]{as_draws_*}} methods as well as to
#'   \code{\link[posterior:subset_draws]{subset_draws}}.
#'
#' @return A data.frame, matrix, or array containing the posterior draws.
#'
#' @seealso \code{\link[brms:draws-brms]{as_draws}},
#'   \code{\link[posterior:subset_draws]{subset_draws}}
#'
#' @export
as.data.frame.brmsfit <- function(x, row.names = NULL, optional = TRUE,
                                  pars = NA, variable = NULL, draw = NULL,
                                  subset = NULL, ...) {
  variable <- use_variable_alias(variable, x, pars = pars, ...)
  draw <- use_alias(draw, subset)
  out <- as_draws_df(x, variable = variable, ...)
  out <- suppressMessages(subset_draws(out, draw = draw, ...))
  unclass_draws(out)
}

#' @rdname as.data.frame.brmsfit
#' @export
as.matrix.brmsfit <- function(x, pars = NA, variable = NULL,
                              draw = NULL, subset = NULL, ...) {
  variable <- use_variable_alias(variable, x, pars = pars, ...)
  draw <- use_alias(draw, subset)
  out <- as_draws_matrix(x, variable = variable, ...)
  out <- suppressMessages(subset_draws(out, draw = draw, ...))
  unclass_draws(out)
}

#' @rdname as.data.frame.brmsfit
#' @export
as.array.brmsfit <- function(x, pars = NA, variable = NULL,
                             draw = NULL, subset = NULL, ...) {
  variable <- use_variable_alias(variable, x, pars = pars, ...)
  draw <- use_alias(draw, subset)
  out <- as_draws_array(x, variable = variable, ...)
  out <- suppressMessages(subset_draws(out, draw = draw, ...))
  unclass_draws(out)
}

# use the deprecated 'pars' alias to 'variable'
use_variable_alias <- function(variable, object, pars = NA, ...) {
  if (!anyNA(pars)) {
    warning2("Argument 'pars' is deprecated. Please use 'variable' instead.")
    variable <- extract_pars(pars, variables(object), ...)
  }
  variable
}

# remove the posterior draws format classes from objects
unclass_draws <- function(x, ...) {
  UseMethod("unclass_draws")
}

#' @export
unclass_draws.default <- function(x, ...) {
  unclass(x)
}

#' @export
unclass_draws.draws_df <- function(x, ...) {
  x <- as.data.frame(x)
  x$.chain <- x$.iteration <- x$.draw <- NULL
  x
}
