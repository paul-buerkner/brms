#' Posterior Predictions of Smooth Terms
#'
#' Compute posterior predictions of smooth \code{s} and \code{t2} terms of
#' models fitted with \pkg{brms}.
#'
#' @inheritParams posterior_epred.brmsfit
#' @param smooth Name of a single smooth term for which predictions should
#'   be computed.
#' @param newdata An optional \code{data.frame} for which to evaluate
#'   predictions. If \code{NULL} (default), the original data of the model is
#'   used. Only those variables appearing in the chosen \code{smooth} term are
#'   required.
#' @param ... Currently ignored.
#'
#' @return An S x N matrix, where S is the number of
#'   posterior draws and N is the number of observations.
#'
#' @examples
#' \dontrun{
#' set.seed(0)
#' dat <- mgcv::gamSim(1, n = 200, scale = 2)
#' fit <- brm(y ~ s(x0) + s(x1) + s(x2) + s(x3), data = dat)
#' summary(fit)
#'
#' newdata <- data.frame(x2 = seq(0, 1, 10))
#' str(posterior_smooths(fit, smooth = "s(x2)", newdata = newdata))
#' }
#'
#' @export
posterior_smooths.brmsfit <- function(object, smooth, newdata = NULL,
                                      resp = NULL, dpar = NULL, nlpar = NULL,
                                      ndraws = NULL, draw_ids = NULL, ...) {
  resp <- validate_resp(resp, object, multiple = FALSE)
  bterms <- brmsterms(exclude_terms(object$formula, smooths_only = TRUE))
  if (!is.null(resp)) {
    stopifnot(is.mvbrmsterms(bterms))
    bterms <- bterms$terms[[resp]]
  }
  if (!is.null(nlpar)) {
    if (length(dpar)) {
      stop2("Cannot use 'dpar' and 'nlpar' at the same time.")
    }
    nlpar <- as_one_character(nlpar)
    nlpars <- names(bterms$nlpars)
    if (!nlpar %in% nlpars) {
      stop2("Invalid argument 'nlpar'. Valid non-linear ",
            "parameters are: ", collapse_comma(nlpars))
    }
    bterms <- bterms$nlpars[[nlpar]]
  } else {
    dpar <- dpar %||% "mu"
    dpar <- as_one_character(dpar)
    dpars <- names(bterms$dpars)
    if (!dpar %in% dpars) {
      stop2("Invalid argument 'dpar'. Valid distributional ",
            "parameters are: ", collapse_comma(dpars))
    }
    bterms <- bterms$dpars[[dpar]]
  }
  posterior_smooths(
    bterms, fit = object, smooth = smooth, newdata = newdata,
    ndraws = ndraws, draw_ids = draw_ids, ...
  )
}

#' @export
posterior_smooths.btl <- function(object, fit, smooth, newdata = NULL,
                                  ndraws = NULL, draw_ids = NULL,
                                  nsamples = NULL, subset = NULL, ...) {
  smooth <- rm_wsp(as_one_character(smooth))
  ndraws <- use_alias(ndraws, nsamples)
  draw_ids <- use_alias(draw_ids, subset)
  object$frame$sm <- tidy_smef(object, fit$data)
  class(object) <- c("bfrl", class(object))
  smef <- object$frame$sm
  smef$term <- rm_wsp(smef$term)
  smterms <- unique(smef$term)
  if (!smooth %in% smterms) {
    stop2("Term '", smooth, "' cannot be found. Available ",
          "smooth terms are: ", collapse_comma(smterms))
  }
  # find relevant variables
  sub_smef <- subset2(smef, term = smooth)
  covars <- all_vars(sub_smef$covars[[1]])
  byvars <- all_vars(sub_smef$byvars[[1]])
  req_vars <- c(covars, byvars)
  # prepare predictions for splines
  sdata <- standata(
    fit, newdata, re_formula = NA, internal = TRUE,
    check_response = FALSE, req_vars = req_vars
  )
  draw_ids <- validate_draw_ids(fit, draw_ids, ndraws)
  draws <- as_draws_matrix(fit)
  draws <- suppressMessages(subset_draws(draws, draw = draw_ids))
  prep_args <- nlist(x = object, draws, sdata, data = fit$data)
  prep <- do_call(prepare_predictions, prep_args)
  # select subset of smooth parameters and design matrices
  i <- which(smterms %in% smooth)[1]
  J <- which(smef$termnum == i)
  scs <- unlist(attr(prep$sm$fe$Xs, "smcols")[J])
  prep$sm$fe$Xs <- prep$sm$fe$Xs[, scs, drop = FALSE]
  prep$sm$fe$bs <- prep$sm$fe$bs[, scs, drop = FALSE]
  prep$sm$re <- prep$sm$re[J]
  prep$family <- brmsfamily("gaussian")
  predictor(prep, i = NULL)
}

#' @export
posterior_smooths.btnl <- function(object, ...) {
  stop2("Non-linear formulas do not contain smooth terms.")
}

#' @rdname posterior_smooths.brmsfit
#' @export
posterior_smooths <- function(object, ...) {
  UseMethod("posterior_smooths")
}
