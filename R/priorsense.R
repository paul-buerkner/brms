#' Prior sensitivity: Create priorsense data
#'
#' The \code{create_priorsense_data.brmsfit} method can be used to
#' create the data structure needed by the \pkg{priorsense} package
#' for performing power-scaling sensitivity analysis. This method is
#' called automatically when performing powerscaling via
#' \code{\link[priorsense:powerscale]{powerscale}} or other related
#' functions, so you will rarely need to call it manually yourself.
#'
#' @param x A \code{brmsfit} object.
#' @param ... Additional arguments passed to \code{\link{log_lik}},
#'   for example \code{newdata}.
#'
#' @return A \code{priorsense_data} object to be used in conjunction
#' with the \pkg{priorsense} package.
#'
#' @examples
#' \dontrun{
#' # fit a model with non-uniform priors
#' fit <- brm(rating ~ treat + period + carry,
#'            data = inhaler, family = sratio(),
#'            prior = set_prior("normal(0, 0.5)"))
#' summary(fit)
#'
#' # The following code requires the 'priorsense' package to be installed:
#' library(priorsense)
#'
#' # perform power-scaling of the prior
#' powerscale(fit, alpha = 1.5, component = "prior")
#'
#' # perform power-scaling sensitivity checks
#' powerscale_sensitivity(fit)
#'
#' # create power-scaling sensitivity plots (for one variable)
#' powerscale_plot_dens(fit, variable = "b_treat")
#' }
#'
#' @exportS3Method priorsense::create_priorsense_data brmsfit
create_priorsense_data.brmsfit <- function(x, ...) {
  priorsense::create_priorsense_data(
    x = get_draws_ps(x),
    fit = x,
    log_prior = log_prior_draws.brmsfit(x, ...),
    log_lik = log_lik_draws.brmsfit(x, ...),
    log_prior_fn = log_prior_draws.brmsfit,
    log_lik_fn = log_lik_draws.brmsfit,
    log_ratio_fn = powerscale_log_ratio,
    ...
  )
}

#' @exportS3Method priorsense::log_lik_draws
log_lik_draws.brmsfit <- function(x, ...) {

  log_lik <- log_lik(x, ...)

  nchains <- nchains(x)
  niters <- niterations(x)

  # check if log-lik was subset, if so, merge the chains
  if (nrow(log_lik) < niters) {
    niters <- nrow(log_lik)
    nchains <- 1
  }

  nobs <- ncol(log_lik)

  dim(log_lik) <- c(niters, nchains, nobs)
  log_lik <- as_draws_array(log_lik)
  nvars <- nvariables(log_lik)
  posterior::variables(log_lik) <- paste0("log_lik[", seq_len(nvars), "]")

  log_lik
}


#' @exportS3Method priorsense::log_prior_draws
log_prior_draws.brmsfit <- function(x, log_prior_name = "lprior", ...) {
  stopifnot(length(log_prior_name) == 1)
  if (!log_prior_name %in% variables(x)) {
    warning2("Variable '", log_prior_name, "' was not found. ",
             "Perhaps you used normalize = FALSE?")
  }
  posterior::subset_draws(
    as_draws_array(x),
    variable = paste0("^", log_prior_name),
    regex = TRUE
  )
}

get_draws_ps <- function(x, variable = NULL, regex = FALSE,
                         log_prior_name = "lprior") {
  excluded_variables <- c(log_prior_name, "lp__")
  draws <- as_draws_df(x, regex = regex)
  if (is.null(variable)) {
    # remove unnecessary variables
    variable <- variables(x)
    variable <- variable[!(variable %in% excluded_variables)]
    draws <- posterior::subset_draws(draws, variable = variable)
  }
  draws
}

powerscale_log_ratio <- function(draws, fit, alpha, component_fn) {
  component_draws <- component_fn(fit)
  component_draws <- rowsums_draws(component_draws)
  component_draws * (alpha - 1)
}

rowsums_draws <- function(x) {
  posterior::draws_array(
    sum = rowSums(
      as_draws_array(x),
      dims = 2
    ),
    .nchains = nchains(x)
  )
}
