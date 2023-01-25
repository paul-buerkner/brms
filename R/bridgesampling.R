#' Log Marginal Likelihood via Bridge Sampling
#'
#' Computes log marginal likelihood via bridge sampling,
#' which can be used in the computation of bayes factors
#' and posterior model probabilities.
#' The \code{brmsfit} method is just a thin wrapper around
#' the corresponding method for \code{stanfit} objects.
#'
#' @aliases bridge_sampler
#'
#' @param samples A \code{brmsfit} object.
#' @param recompile Logical, indicating whether the Stan model should be
#'   recompiled. This may be necessary if you are running bridge sampling on
#'   another machine than the one used to fit the model. No recompilation
#'   is done by default.
#' @param ... Additional arguments passed to
#'   \code{\link[bridgesampling:bridge_sampler]{bridge_sampler.stanfit}}.
#'
#' @details Computing the marginal likelihood requires samples of all variables
#'   defined in Stan's \code{parameters} block to be saved. Otherwise
#'   \code{bridge_sampler} cannot be computed. Thus, please set \code{save_pars
#'   = save_pars(all = TRUE)} in the call to \code{brm}, if you are planning to
#'   apply \code{bridge_sampler} to your models.
#'
#'   The computation of marginal likelihoods based on bridge sampling requires
#'   a lot more posterior draws than usual. A good conservative
#'   rule of thump is perhaps 10-fold more draws (read: the default of 4000
#'   draws may not be enough in many cases). If not enough posterior
#'   draws are provided, the bridge sampling algorithm tends to be
#'   unstable leading to considerably different results each time it is run.
#'   We thus recommend running \code{bridge_sampler}
#'   multiple times to check the stability of the results.
#'
#'   More details are provided under
#'   \code{\link[bridgesampling:bridge_sampler]{bridgesampling::bridge_sampler}}.
#'
#' @seealso \code{
#'   \link[brms:bayes_factor.brmsfit]{bayes_factor},
#'   \link[brms:post_prob.brmsfit]{post_prob}
#' }
#'
#' @examples
#' \dontrun{
#' # model with the treatment effect
#' fit1 <- brm(
#'   count ~ zAge + zBase + Trt,
#'   data = epilepsy, family = negbinomial(),
#'   prior = prior(normal(0, 1), class = b),
#'   save_pars = save_pars(all = TRUE)
#' )
#' summary(fit1)
#' bridge_sampler(fit1)
#'
#' # model without the treatment effect
#' fit2 <- brm(
#'   count ~ zAge + zBase,
#'   data = epilepsy, family = negbinomial(),
#'   prior = prior(normal(0, 1), class = b),
#'   save_pars = save_pars(all = TRUE)
#' )
#' summary(fit2)
#' bridge_sampler(fit2)
#' }
#'
#' @method bridge_sampler brmsfit
#' @importFrom bridgesampling bridge_sampler
#' @export bridge_sampler
#' @export
bridge_sampler.brmsfit <- function(samples, recompile = FALSE, ...) {
  out <- get_criterion(samples, "marglik")
  if (inherits(out, "bridge") && !is.na(out$logml)) {
    # return precomputed criterion
    return(out)
  }
  samples <- restructure(samples)
  if (samples$version$brms <= "1.8.0") {
    stop2(
      "Models fitted with brms 1.8.0 or lower are not ",
      "usable in method 'bridge_sampler'."
    )
  }
  if (!is_normalized(samples$model)) {
    stop2(
      "The Stan model has to be normalized to be ",
      "usable in method 'bridge_sampler'."
    )
  }
  # otherwise bridge_sampler may fail in a new R session or on another machine
  samples <- update_misc_env(samples, recompile = recompile)
  out <- try(bridge_sampler(samples$fit, ...))
  if (is(out, "try-error")) {
    stop2(
      "Bridgesampling failed. Perhaps you did not set ",
      "'save_pars = save_pars(all = TRUE)' when fitting your model? ",
      "If you are running bridge sampling on another machine than the one ",
      "used to fit the model, you may need to set recompile = TRUE."
    )
  }
  out
}

#' Bayes Factors from Marginal Likelihoods
#'
#' Compute Bayes factors from marginal likelihoods.
#'
#' @aliases bayes_factor
#'
#' @param x1 A \code{brmsfit} object
#' @param x2 Another \code{brmsfit} object based on the same responses.
#' @param log Report Bayes factors on the log-scale?
#' @param ... Additional arguments passed to
#'   \code{\link[brms:bridge_sampler.brmsfit]{bridge_sampler}}.
#'
#' @details Computing the marginal likelihood requires samples
#'   of all variables defined in Stan's \code{parameters} block
#'   to be saved. Otherwise \code{bayes_factor} cannot be computed.
#'   Thus, please set \code{save_all_pars = TRUE} in the call to \code{brm},
#'   if you are planning to apply \code{bayes_factor} to your models.
#'
#'   The computation of Bayes factors based on bridge sampling requires
#'   a lot more posterior samples than usual. A good conservative
#'   rule of thumb is perhaps 10-fold more samples (read: the default of 4000
#'   samples may not be enough in many cases). If not enough posterior
#'   samples are provided, the bridge sampling algorithm tends to be unstable,
#'   leading to considerably different results each time it is run.
#'   We thus recommend running \code{bayes_factor}
#'   multiple times to check the stability of the results.
#'
#'   More details are provided under
#'   \code{\link[bridgesampling:bf]{bridgesampling::bayes_factor}}.
#'
#' @seealso \code{
#'   \link[brms:bridge_sampler.brmsfit]{bridge_sampler},
#'   \link[brms:post_prob.brmsfit]{post_prob}
#' }
#'
#' @examples
#' \dontrun{
#' # model with the treatment effect
#' fit1 <- brm(
#'   count ~ zAge + zBase + Trt,
#'   data = epilepsy, family = negbinomial(),
#'   prior = prior(normal(0, 1), class = b),
#'   save_all_pars = TRUE
#' )
#' summary(fit1)
#'
#' # model without the treatment effect
#' fit2 <- brm(
#'   count ~ zAge + zBase,
#'   data = epilepsy, family = negbinomial(),
#'   prior = prior(normal(0, 1), class = b),
#'   save_all_pars = TRUE
#' )
#' summary(fit2)
#'
#' # compute the bayes factor
#' bayes_factor(fit1, fit2)
#' }
#'
#' @method bayes_factor brmsfit
#' @importFrom bridgesampling bayes_factor
#' @export bayes_factor
#' @export
bayes_factor.brmsfit <- function(x1, x2, log = FALSE, ...) {
  model_name_1 <- deparse0(substitute(x1))
  model_name_2 <- deparse0(substitute(x2))
  match_response(list(x1, x2))
  bridge1 <- bridge_sampler(x1, ...)
  bridge2 <- bridge_sampler(x2, ...)
  out <- bayes_factor(bridge1, bridge2, log = log)
  attr(out, "model_names") <- c(model_name_1, model_name_2)
  out
}

#' Posterior Model Probabilities from Marginal Likelihoods
#'
#' Compute posterior model probabilities from marginal likelihoods.
#' The \code{brmsfit} method is just a thin wrapper around
#' the corresponding method for \code{bridge} objects.
#'
#' @aliases post_prob
#'
#' @inheritParams loo.brmsfit
#' @param prior_prob Numeric vector with prior model probabilities.
#'   If omitted, a uniform prior is used (i.e., all models are equally
#'   likely a priori). The default \code{NULL} corresponds to equal
#'   prior model weights.
#'
#' @details Computing the marginal likelihood requires samples
#'   of all variables defined in Stan's \code{parameters} block
#'   to be saved. Otherwise \code{post_prob} cannot be computed.
#'   Thus, please set \code{save_all_pars = TRUE} in the call to \code{brm},
#'   if you are planning to apply \code{post_prob} to your models.
#'
#'   The computation of model probabilities based on bridge sampling requires
#'   a lot more posterior samples than usual. A good conservative
#'   rule of thump is perhaps 10-fold more samples (read: the default of 4000
#'   samples may not be enough in many cases). If not enough posterior
#'   samples are provided, the bridge sampling algorithm tends to be
#'   unstable leading to considerably different results each time it is run.
#'   We thus recommend running \code{post_prob}
#'   multiple times to check the stability of the results.
#'
#'   More details are provided under
#'   \code{\link[bridgesampling:post_prob]{bridgesampling::post_prob}}.
#'
#' @seealso \code{
#'   \link[brms:bridge_sampler.brmsfit]{bridge_sampler},
#'   \link[brms:bayes_factor.brmsfit]{bayes_factor}
#' }
#'
#' @examples
#' \dontrun{
#' # model with the treatment effect
#' fit1 <- brm(
#'   count ~ zAge + zBase + Trt,
#'   data = epilepsy, family = negbinomial(),
#'   prior = prior(normal(0, 1), class = b),
#'   save_all_pars = TRUE
#' )
#' summary(fit1)
#'
#' # model without the treatent effect
#' fit2 <- brm(
#'   count ~ zAge + zBase,
#'   data = epilepsy, family = negbinomial(),
#'   prior = prior(normal(0, 1), class = b),
#'   save_all_pars = TRUE
#' )
#' summary(fit2)
#'
#' # compute the posterior model probabilities
#' post_prob(fit1, fit2)
#'
#' # specify prior model probabilities
#' post_prob(fit1, fit2, prior_prob = c(0.8, 0.2))
#' }
#'
#' @method post_prob brmsfit
#' @importFrom bridgesampling post_prob
#' @export post_prob
#' @export
post_prob.brmsfit <- function(x, ..., prior_prob = NULL, model_names = NULL) {
  args <- split_dots(x, ..., model_names = model_names)
  models <- args$models
  args$models <- NULL
  bs <- vector("list", length(models))
  for (i in seq_along(models)) {
    bs[[i]] <- do_call(bridge_sampler, c(list(models[[i]]), args))
  }
  model_names <- names(models)
  do_call(post_prob, c(bs, nlist(prior_prob, model_names)))
}
