#' Index \code{brmsfit} objects
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
nvariables.brmsfit <- function(x) {
  length(variables(x))
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
  niteration(x) * nchains(x)
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
#' @param x A \code{brmsfit} object or another \R object for which
#' the methods are defined.
#' @param ... Arguments passed to individual methods (if applicable).
#' 
#' @seealso \code{\link[posterior:draws]{draws}}
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
as_draws.brmsfit <- function(x, ...) {
  as_draws(x$fit, ...) 
}

#' @rdname draws-brms
#' @importFrom posterior as_draws_matrix
#' @method as_draws_matrix brmsfit
#' @export
#' @export as_draws_matrix
as_draws_matrix.brmsfit <- function(x, ...) {
  as_draws_matrix(x$fit, ...) 
}

#' @rdname draws-brms
#' @importFrom posterior as_draws_array
#' @method as_draws_array brmsfit
#' @export
#' @export as_draws_array
as_draws_array.brmsfit <- function(x, ...) {
  as_draws_array(x$fit, ...) 
}

#' @rdname draws-brms
#' @importFrom posterior as_draws_df
#' @method as_draws_df brmsfit
#' @export
#' @export as_draws_df
as_draws_df.brmsfit <- function(x, ...) {
  as_draws_df(x$fit, ...) 
}

#' @rdname draws-brms
#' @importFrom posterior as_draws_list
#' @method as_draws_list brmsfit
#' @export
#' @export as_draws_list
as_draws_list.brmsfit <- function(x, ...) {
  as_draws_list(x$fit, ...) 
}

#' @rdname draws-brms
#' @importFrom posterior as_draws_rvars
#' @method as_draws_rvars brmsfit
#' @export
#' @export as_draws_rvars
as_draws_rvars.brmsfit <- function(x, ...) {
  as_draws_rvars(x$fit, ...) 
}
