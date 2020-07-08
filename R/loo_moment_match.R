#' Moment matching for efficient approximate leave-one-out cross-validation
#' 
#' Moment matching for efficient approximate leave-one-out cross-validation 
#' (LOO-CV). See \code{\link[loo:loo_moment_match]{loo_moment_match}} 
#' for more details.
#' 
#' @aliases loo_moment_match
#' 
#' @param x An object of class \code{brmsfit}.
#' @param loo An object of class \code{loo} originally created from \code{x}.
#' @param ... Further arguments passed to the underlying methods.
#' @return An object of the class \code{loo}.
#' 
#' @details The moment matching algorithm requires samples 
#'   of all variables defined in Stan's \code{parameters} block
#'   to be saved. Otherwise \code{loo_moment_match} cannot be computed.
#'   Thus, please set \code{save_all_pars = TRUE} in the call to \code{brm},
#'   if you are planning to apply \code{loo_moment_match} to your models.
#' 
#' @examples 
#' \dontrun{
#' fit1 <- brm(count ~ zAge + zBase * Trt + (1|patient),
#'             data = epilepsy, family = poisson(),
#'             save_all_pars = TRUE)
#' # throws warning about some pareto k estimates being too high
#' (loo1 <- loo(fit1))
#' (mmloo1 <- loo_moment_match(fit1, loo = loo1, k_thres = 0.7))
#' }
#' 
#' @importFrom loo loo_moment_match
#' @export loo_moment_match
#' @export
loo_moment_match.brmsfit <- function(x, loo, ...) {
  # TODO: support more arguments such as 'newdata' or 'resp'
  # ensure compatibility with objects not created in the current R session
  x$fit@.MISC <- suppressMessages(brm(fit = x, chains = 0))$fit@.MISC
  out <- try(loo::loo_moment_match.default(
    x, loo = loo, 
    post_draws = as.matrix, 
    log_lik_i = .log_lik_i, 
    unconstrain_pars = .unconstrain_pars,
    log_prob_upars = .log_prob_upars,
    log_lik_i_upars = .log_lik_i_upars,
    ...
  ))
  if (is(out, "try-error")) {
    stop2(
      "'loo_moment_match' failed. Did you set 'save_all_pars' ",
      "to TRUE when fitting your model?"
    )
  }
  out
}

# compute a vector of log-likelihood values for the ith observation
.log_lik_i <- function(x, i, ...) {
  as.vector(log_lik(x, newdata = x$data[i, , drop = FALSE], ...))
}

# transform parameters to the unconstraint space
.unconstrain_pars <- function(x, pars, ...) {
  unconstrain_pars_stanfit(x$fit, pars = pars, ...)
}

# compute log_prob for each posterior draws on the unconstrained space
.log_prob_upars <- function(x, upars, ...) {
  log_prob_upars_stanfit(x$fit, upars = upars, ...)
}

# transform parameters to the constraint space
.update_pars <- function(x, upars, ...) {
  # list with one element per posterior draw
  pars <- apply(upars, 1, rstan::constrain_pars, object = x$fit)
  # transform samples
  nsamples <- length(pars)
  pars <- unlist(pars)
  npars <- length(pars) / nsamples
  dim(pars) <- c(npars, nsamples)
  # add dummy 'lp__' samples
  pars <- rbind(pars, rep(0, nsamples))
  # bring samples into the right structure
  new_samples <- named_list(x$fit@sim$fnames_oi_old, list(numeric(nsamples)))
  stopifnot(length(new_samples) == nrow(pars))
  for (i in seq_len(npars)) {
    new_samples[[i]] <- pars[i, ]
  }
  # create new sim object to overwrite x$fit@sim
  x$fit@sim <- list(
    samples = list(new_samples),
    iter = nsamples,
    thin = 1,
    warmup = 0,
    chains = 1,
    n_save = nsamples,
    warmup2 = 0,
    permutation = list(seq_len(nsamples)),
    pars_oi = x$fit@sim$pars_oi_old,
    dims_oi = x$fit@sim$dims_oi_old,
    fnames_oi = x$fit@sim$fnames_oi_old,
    n_flatnames = length(x$fit@sim$fnames_oi_old)
  ) 
  x$fit@stan_args <- list(
    list(chain_id = 1, iter = nsamples, thin = 1, warmup = 0)
  )
  rename_pars(x)
}

# compute log_lik values based on the unconstrained parameters
.log_lik_i_upars <- function(x, upars, i, samples = NULL, 
                             subset = NULL, ...) {
  # do not pass subset or nsamples further to avoid subsetting twice
  x <- .update_pars(x, upars = upars, ...)
  .log_lik_i(x, i = i, ...)
}

# -------- will be imported from rstan at some point -------
# transform parameters to the unconstraint space
unconstrain_pars_stanfit <- function(x, pars, ...) {
  skeleton <- .create_skeleton(x@sim$pars_oi, x@par_dims[x@sim$pars_oi])
  upars <- apply(pars, 1, FUN = function(theta) {
    rstan::unconstrain_pars(x, pars = .rstan_relist(theta, skeleton))
  })
  # for one parameter models
  if (is.null(dim(upars))) {
    dim(upars) <- c(1, length(upars))
  }
  t(upars)
}

# compute log_prob for each posterior draws on the unconstrained space
log_prob_upars_stanfit <- function(x, upars, ...) {
  apply(upars, 1, rstan::log_prob, object = x,
        adjust_transform = TRUE, gradient = FALSE)
}

# create a named list of draws for use with rstan methods
.rstan_relist <- function (x, skeleton) {
  out <- utils::relist(x, skeleton)
  for (i in seq_along(skeleton)) {
    dim(out[[i]]) <- dim(skeleton[[i]])
  }
  out
}

# rstan helper function to get dims of parameters right
.create_skeleton <- function (pars, dims) {
  out <- lapply(seq_along(pars), function(i) {
    len_dims <- length(dims[[i]])
    if (len_dims < 1) return(0)
    return(array(0, dim = dims[[i]]))
  })
  names(out) <- pars
  out
}
