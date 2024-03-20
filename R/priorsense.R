##' @rdname create-powerscaling-data
##' @importFrom priorsense create_powerscaling_data
##' @export create_powerscaling_data
##' @export
create_powerscaling_data.brmsfit <- function(x, ...) {

  priorsense::create_powerscaling_data.default(
    x = x,
    log_prior = log_prior_brmsfit,
    log_lik = joint_log_lik_brmsfit,
    get_draws = get_draws_brmsfit,
    unconstrain_pars = unconstrain_pars,
    constrain_pars = rstan::constrain_pars,
    log_prob_upars = log_prob_upars,
    log_ratio_upars = log_ratio_upars,
    ...
  )
}

##' @rdname powerscale-overview
##' @importFrom priorsense powerscale
##' @export powerscale
##' @export
powerscale.brmsfit <- function(x,
                               component,
                               alpha,
                               ...
                               ) {
  psd <- create_powerscaling_data.brmsfit(x, ...)

  priorsense::powerscale.powerscaling_data(
    psd,
    component = component,
    alpha = alpha,
    ...
  )

}

##' @rdname powerscale-overview
##' @importFrom priorsense powerscale_sequence
##' @export powerscale_sequence
##' @export
powerscale_sequence.brmsfit <- function(x,
                                        ...
                                        ) {

  psd <- create_powerscaling_data.brmsfit(x, ...)

  priorsense::powerscale_sequence.powerscaling_data(psd, ...)

}

##' @rdname powerscale-sensitivity
##' @importFrom priorsense powerscale_sensitivity
##' @exrpot powerscale_sensitivity
##' @export
powerscale_sensitivity.brmsfit <- function(x,
                                           ...
                                           ) {

  psd <- create_powerscaling_data.brmsfit(x, ...)

  priorsense::powerscale_sensitivity.powerscaling_data(
    psd,
    ...
  )
}


##' @rdname joint_log_lik
##' @export
joint_log_lik_brmsfit <- function(x, ...) {

  log_lik <- rowSums(log_lik(x, ...))
  chains <- x$fit@sim$chains

  log_lik <- posterior::draws_array(
    log_lik = log_lik,
    .nchains = chains
  )

  return(log_lik)
}


##' @rdname log_prior
##' @export
log_prior_brmsfit <- function(x, ...) {

  log_prior <- posterior::subset_draws(posterior::as_draws_array(x), variable = "lprior")

  return(log_prior)
}

get_draws_brmsfit <- function(x, variable = NULL, regex = FALSE, ...) {

  excluded_variables <- c("lprior", "lp__")
  draws <- posterior::as_draws_df(x, variable = variable, regex = regex)

  if (is.null(variable)) {
    # remove unnecessary variables
    variable <- posterior::variables(x)
    variable <- variable[!(variable %in% excluded_variables)]

    draws <- posterior::subset_draws(draws, variable = variable)
  }

  return(draws)
}

moment_match.brmsfit <- function(x, psis, ...) {
  # ensure compatibility with objects not created in the current R session
  x$fit@.MISC <- suppressMessages(brm(fit = x, chains = 0))$fit@.MISC
  mm <- try(moment_match.default(
    x,
    psis = psis, post_draws = as.matrix,
    unconstrain_pars = unconstrain_pars.brmsfit,
    log_prob_upars = log_prob_upars.brmsfit,
    log_ratio_upars = log_ratio_upars.brmsfit,
    nchains = posterior::nchains(x),
    ...
  ))
  if (methods::is(out, "try-error")) {
    stop(
      "'moment_match' failed. Did you set 'save_all_pars' ",
      "to TRUE when fitting your brms model?"
    )
  }
  return(mm)
}

unconstrain_pars.brmsfit <- function(x, pars, ...) {
  unconstrain_pars.stanfit(x$fit, pars = pars, ...)
}

log_prob_upars.brmsfit <- function(x, upars, ...) {
  log_prob_upars.stanfit(x$fit, upars = upars, ...)
}

update_pars.brmsfit <- function(x, upars, ...) {
  x$fit <- update_pars(x$fit, upars = upars, save_old_pars = FALSE, ...)
  brms::rename_pars(x)
}

log_ratio_upars.brmsfit <- function(x, upars, component_fn, samples = NULL,
                                    subset = NULL, ...) {
  # do not pass subset or nsamples further to avoid subsetting twice
  x <- update_pars(x, upars = upars, ...)
  component_draws <- component_fn(x)
  scaled_log_ratio(component_draws, ...)
}

