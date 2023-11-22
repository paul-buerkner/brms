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
                               probs = c(0.025, 0.975), ...) {
  log <- as_one_logical(log)
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
  out <- lapply(seq_len(N), log_lik_mixture, prep = prep)
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
