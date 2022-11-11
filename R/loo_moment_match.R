#' Moment matching for efficient approximate leave-one-out cross-validation
#'
#' Moment matching for efficient approximate leave-one-out cross-validation
#' (LOO-CV). See \code{\link[loo:loo_moment_match]{loo_moment_match}}
#' for more details.
#'
#' @aliases loo_moment_match
#'
#' @inheritParams predict.brmsfit
#' @param x An object of class \code{brmsfit}.
#' @param loo An object of class \code{loo} originally created from \code{x}.
#' @param k_threshold The threshold at which Pareto \eqn{k}
#'   estimates are treated as problematic. Defaults to \code{0.7}.
#'   See \code{\link[loo:pareto-k-diagnostic]{pareto_k_ids}}
#'   for more details.
#' @param check Logical; If \code{TRUE} (the default), some checks
#'   check are performed if the \code{loo} object was generated
#'   from the \code{brmsfit} object passed to argument \code{fit}.
#' @param recompile Logical, indicating whether the Stan model should be
#'   recompiled. This may be necessary if you are running moment matching on
#'   another machine than the one used to fit the model. No recompilation
#'   is done by default.
#' @param ... Further arguments passed to the underlying methods.
#'   Additional arguments initially passed to \code{\link{loo}},
#'   for example, \code{newdata} or \code{resp} need to be passed
#'   again to \code{loo_moment_match} in order for the latter
#'   to work correctly.
#' @return An updated object of class \code{loo}.
#'
#' @details The moment matching algorithm requires draws of all variables
#'   defined in Stan's \code{parameters} block to be saved. Otherwise
#'   \code{loo_moment_match} cannot be computed. Thus, please set
#'   \code{save_pars = save_pars(all = TRUE)} in the call to \code{\link{brm}},
#'   if you are planning to apply \code{loo_moment_match} to your models.
#'
#' @references
#'   Paananen, T., Piironen, J., Buerkner, P.-C., Vehtari, A. (2021).
#'   Implicitly Adaptive Importance Sampling. Statistics and Computing.
#'
#' @examples
#' \dontrun{
#' fit1 <- brm(count ~ zAge + zBase * Trt + (1|patient),
#'             data = epilepsy, family = poisson(),
#'             save_pars = save_pars(all = TRUE))
#'
#' # throws warning about some pareto k estimates being too high
#' (loo1 <- loo(fit1))
#' (mmloo1 <- loo_moment_match(fit1, loo = loo1))
#' }
#'
#' @importFrom loo loo_moment_match
#' @export loo_moment_match
#' @export
loo_moment_match.brmsfit <- function(x, loo, k_threshold = 0.7, newdata = NULL,
                                     resp = NULL, check = TRUE,
                                     recompile = FALSE, ...) {
  stopifnot(is.loo(loo), is.brmsfit(x))
  if (is.null(newdata)) {
    newdata <- model.frame(x)
  } else {
    newdata <- as.data.frame(newdata)
  }
  check <- as_one_logical(check)
  if (check) {
    yhash_loo <- attr(loo, "yhash")
    yhash_fit <- hash_response(x, newdata = newdata)
    if (!is_equal(yhash_loo, yhash_fit)) {
      stop2(
        "Response values used in 'loo' and 'x' do not match. ",
        "If this is a false positive, please set 'check' to FALSE."
      )
    }
  }
  # otherwise loo_moment_match may fail in a new R session or on another machine
  x <- update_misc_env(x, recompile = recompile)
  out <- try(loo::loo_moment_match.default(
    x, loo = loo,
    post_draws = as.matrix,
    log_lik_i = .log_lik_i,
    unconstrain_pars = .unconstrain_pars,
    log_prob_upars = .log_prob_upars,
    log_lik_i_upars = .log_lik_i_upars,
    k_threshold = k_threshold,
    newdata = newdata,
    resp = resp, ...
  ))
  if (is(out, "try-error")) {
    stop2(
      "Moment matching failed. Perhaps you did not set ",
      "'save_pars = save_pars(all = TRUE)' when fitting your model? ",
      "If you are running moment matching on another machine than the one ",
      "used to fit the model, you may need to set recompile = TRUE."
    )
  }
  out
}

# compute a vector of log-likelihood values for the ith observation
.log_lik_i <- function(x, i, newdata, ...) {
  as.vector(log_lik(x, newdata = newdata[i, , drop = FALSE], ...))
}

# transform parameters to the unconstrained space
.unconstrain_pars <- function(x, pars, ...) {
  unconstrain_pars_stanfit(x$fit, pars = pars, ...)
}

# compute log_prob for each posterior draws on the unconstrained space
.log_prob_upars <- function(x, upars, ...) {
  x <- update_misc_env(x, only_windows = TRUE)
  log_prob_upars_stanfit(x$fit, upars = upars, ...)
}

# transform parameters to the constraint space
.update_pars <- function(x, upars, ...) {
  # list with one element per posterior draw
  pars <- apply(upars, 1, .constrain_pars, x = x)
  # select required parameters only
  pars <- lapply(pars, "[", x$fit@sim$pars_oi_old)
  # transform draws
  ndraws <- length(pars)
  pars <- unlist(pars)
  npars <- length(pars) / ndraws
  dim(pars) <- c(npars, ndraws)
  # add dummy 'lp__' draws
  pars <- rbind(pars, rep(0, ndraws))
  # bring draws into the right structure
  new_draws <- named_list(x$fit@sim$fnames_oi_old, list(numeric(ndraws)))
  if (length(new_draws) != nrow(pars)) {
    stop2("Updating parameters in `loo_moment_match.brmsfit' failed. ",
          "Please report a bug at https://github.com/paul-buerkner/brms.")
  }
  for (i in seq_len(npars)) {
    new_draws[[i]] <- pars[i, ]
  }
  # create new sim object to overwrite x$fit@sim
  x$fit@sim <- list(
    samples = list(new_draws),
    iter = ndraws,
    thin = 1,
    warmup = 0,
    chains = 1,
    n_save = ndraws,
    warmup2 = 0,
    permutation = list(seq_len(ndraws)),
    pars_oi = x$fit@sim$pars_oi_old,
    dims_oi = x$fit@sim$dims_oi_old,
    fnames_oi = x$fit@sim$fnames_oi_old,
    n_flatnames = length(x$fit@sim$fnames_oi_old)
  )
  x$fit@stan_args <- list(
    list(chain_id = 1, iter = ndraws, thin = 1, warmup = 0)
  )
  rename_pars(x)
}

# wrapper around rstan::constrain_pars
# ensures that the right posterior draws are excluded
.constrain_pars <- function(upars, x) {
  out <- rstan::constrain_pars(upars, object = x$fit)
  out[x$exclude] <- NULL
  out
}

# compute log_lik values based on the unconstrained parameters
.log_lik_i_upars <- function(x, upars, i, ndraws = NULL,
                             draw_ids = NULL, ...) {
  # do not pass draw_ids or ndraws further to avoid subsetting twice
  x <- update_misc_env(x, only_windows = TRUE)
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
