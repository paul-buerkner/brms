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
