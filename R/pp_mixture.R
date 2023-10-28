#' Posterior Probabilities of Mixture Component Memberships
#'
#' Compute the posterior probabilities of mixture component
#' memberships for each observation including uncertainty
#' estimates.
#'
#' @inheritParams predict.brmsfit
#' @param x An \R object usually of class \code{brmsfit}.
#' @param log Logical; Indicates whether to return
#'   probabilities on the log-scale.
#' @param z_mix Logical. If TRUE, returns the mixture component
#'   memberships for the hurdle/zero-inflation state in hurdle
#'   and zero-inflated models. For an explicit mixture model 
#'   including a zero-inflated model in the mixture, probabilities
#'   are returned for each mixture component, treating the 
#'   zero-inflation state and the non-zero-inflation state as
#'   separate mixture components.
#'
#' @return
#' If \code{summary = TRUE}, an N x E x K array,
#' where N is the number of observations, K is the number
#' of mixture components, and E is equal to \code{length(probs) + 2}.
#' If \code{summary = FALSE}, an S x N x K array, where
#' S is the number of posterior draws.
#'
#' @details
#' The returned probabilities can be written as
#' \eqn{P(Kn = k | Yn)}, that is the posterior probability
#' that observation n originates from component k.
#' They are computed using Bayes' Theorem
#' \deqn{P(Kn = k | Yn) = P(Yn | Kn = k) P(Kn = k) / P(Yn),}
#' where \eqn{P(Yn | Kn = k)} is the (posterior) likelihood
#' of observation n for component k, \eqn{P(Kn = k)} is
#' the (posterior) mixing probability of component k
#' (i.e. parameter \code{theta<k>}), and
#' \deqn{P(Yn) = \sum (k=1,...,K) P(Yn | Kn = k) P(Kn = k)}
#' is a normalizing constant.
#'
#' @examples
#' \dontrun{
#' ## simulate some data
#' set.seed(1234)
#' dat <- data.frame(
#'   y = c(rnorm(100), rnorm(50, 2)),
#'   x = rnorm(150)
#' )
#' ## fit a simple normal mixture model
#' mix <- mixture(gaussian, nmix = 2)
#' prior <- c(
#'   prior(normal(0, 5), Intercept, nlpar = mu1),
#'   prior(normal(0, 5), Intercept, nlpar = mu2),
#'   prior(dirichlet(2, 2), theta)
#' )
#' fit1 <- brm(bf(y ~ x), dat, family = mix,
#'             prior = prior, chains = 2, init = 0)
#' summary(fit1)
#'
#' ## compute the membership probabilities
#' ppm <- pp_mixture(fit1)
#' str(ppm)
#'
#' ## extract point estimates for each observation
#' head(ppm[, 1, ])
#'
#' ## classify every observation according to
#' ## the most likely component
#' apply(ppm[, 1, ], 1, which.max)
#' }
#'
#' @export
pp_mixture.brmsfit <- function(x, newdata = NULL, re_formula = NULL,
                               resp = NULL, ndraws = NULL, draw_ids = NULL,
                               log = FALSE, summary = TRUE, robust = FALSE,
                               probs = c(0.025, 0.975), z_mix = FALSE, 
                               ...) {
  log <- as_one_logical(log)
  z_mix <- as_one_logical(z_mix)
  contains_draws(x)
  x <- restructure(x)
  stopifnot_resp(x, resp)
  if (is_mv(x)) {
    resp <- validate_resp(resp, x$formula$responses, multiple = FALSE)
    family <- x$family[[resp]]
  } else {
    family <- x$family
  }
  
  if (!is.mixfamily(family)) {
    all_families <- family$family
  } else {
    all_families <- lapply(fit$family$mix, get_family)
  }
  any_zihu <- any(unlist(lapply(all_families, is.zihufamily)))
  
  if (z_mix & !any_zihu) {
    stop2("Method 'pp_mixture(..., z_mix = TRUE)' can only be applied to zero-inflated and hurdle models.")
  }
  if (!is.mixfamily(family) & !z_mix) {
    stop2("Method 'pp_mixture' can only be applied to mixture models.")
  }
  
  prep <- prepare_predictions(
    x, newdata = newdata, re_formula = re_formula, resp = resp,
    draw_ids = draw_ids, ndraws = ndraws, check_response = TRUE, ...
  )
  stopifnot(is.brmsprep(prep))
  prep$pp_mixture <- TRUE
  for (dp in names(prep$dpars)) {
    prep$dpars[[dp]] <- get_dpar(prep, dpar = dp)
  }
  N <- choose_N(prep)
  if (is.mixfamily(family) & !z_mix) {
    out <- lapply(seq_len(N), log_lik_mixture, prep = prep)
  } else if (!is.mixfamily(family) & z_mix) {
    out <- lapply(seq_len(N), log_lik_mixture, prep = prep, mix = FALSE, zihu = TRUE)
  } else {
    out <- lapply(seq_len(N), log_lik_mixture, prep = prep, zihu = TRUE)
  }
  
  out <- abind(out, along = 3)
  out <- aperm(out, c(1, 3, 2))
  old_order <- prep$old_order
  sort <- isTRUE(ncol(out) != length(old_order))
  out <- reorder_obs(out, old_order, sort = sort)
  if (!log) {
    out <- exp(out)
  }
  if (summary) {
    out <- posterior_summary(out, probs = probs, robust = robust)
    dimnames(out) <- list(
      seq_len(nrow(out)), colnames(out),
      paste0("P(K = ", seq_len(dim(out)[3]), " | Y)")
    )
  }
  out
}

#' @rdname pp_mixture.brmsfit
#' @export
pp_mixture <- function(x, ...) {
  UseMethod("pp_mixture")
}

# 
# latent_zi_state_prob <- function (brmsfit, ...) {
#   fam <- brmsfit$family$family
#   if(!is.zifamily(fam)) {
#     stop2(paste0(
#       "latent zi probs available only for discrete zero-inflated families (",
#       paste(zi_families(), collapse = ", "),
#       "). Supplied brmsfit has family: ",
#       fam
#     ))
#   }
#   lik <- exp(log_lik(brmsfit, ...))
#   zi_lik_given_obs_zero <- posterior_linpred(
#     brmsfit, 
#     dpar = "zi",
#     transform = TRUE,
#     ...
#   )
#   resp_name <- brmsfit$formula$resp
#   resp_data <- brmsfit$data[[resp_name]]
#   resp_is_zero <- as.integer(resp_data == 0)
#   zi_lik <- sweep(zi_lik_given_obs_zero, MARGIN=2, resp_is_zero, `*`)
#   zi_lik / lik
# }
# 
# 
# latent_hu_state_prob <- function (brmsfit) {
#   fam <- brmsfit$family$family
#   if(!is.hufamily(fam)) {
#     stop2(paste0(
#       "latent hu probs available only for hurdle and zero-inflated continuous families (",
#       paste(hu_families(), collapse = ", "),
#       "). Supplied brmsfit has family: ",
#       fam
#     ))
#   }
#   
#   resp_name <- brmsfit$formula$resp 
#   
#   if (is.null(newdata)) {
#     the_data <- brmsfit$data
#   } else {
#     the_data <- newdata
#   }
#   
#   if (fam == "zero_one_inflated_beta") {
#     out <- data.frame(
#       ni = as.integer(!(the_data[[resp]] %in% c(0, 1))),
#       zi = as.integer(the_data[[resp]] == 0),
#       oi = as.integer(the_data[[resp]] == 1)
#   } else {
#     out <- data.frame(
#       ni = as.integer(the_data != 0),
#       zi = as.integer(the_data == 0)
#     )
#   }
#   out
# }
# 

#' List of all hurdle families, where "hurdle" means anything where observed 
#' zeros are guaranteed to be due to the zero-inflation component.
hu_families <- function() {
  c("hurdle_poisson", "hurdle_negbinomial", "hurdle_gamma", 
    "hurdle_lognormal", "hurdle_cumulative", "zero_inflated_beta",
    "zero_one_inflated_beta")
}

#' List of all true zero-inflated families, where observed zeros may arise
#' either from the zero-inflation component or from the other component.
zi_families <- function() {
  c("zero_inflated_poisson", "zero_inflated_negbinomial",
    "zero_inflated_binomial", "zero_inflated_beta_binomial")
}

is.zihufamily <- function(x) {
  x %in% c(hu_families(), zi_families())
}

is.zifamily <- function(x) {
  x %in% zi_families()
}

is.hufamily <- function(x) {
  x %in% hu_families()
}

get_family <- function(x) {
  x$family
}

