#' Draws from the Expected Value of the Posterior Predictive Distribution
#'
#' Compute posterior draws of the expected value of the posterior predictive
#' distribution. Can be performed for the data used to fit the model (posterior
#' predictive checks) or for new data. By definition, these predictions have
#' smaller variance than the posterior predictions performed by the
#' \code{\link{posterior_predict.brmsfit}} method. This is because only the
#' uncertainty in the expected value of the posterior predictive distribution is
#' incorporated in the draws computed by \code{posterior_epred} while the
#' residual error is ignored there. However, the estimated means of both methods
#' averaged across draws should be very similar.
#'
#' @aliases pp_expect
#'
#' @inheritParams posterior_predict.brmsfit
#' @param dpar Optional name of a predicted distributional parameter.
#'  If specified, expected predictions of this parameters are returned.
#' @param nlpar Optional name of a predicted non-linear parameter.
#'  If specified, expected predictions of this parameters are returned.
#'
#' @return An \code{array} of draws. For
#'   categorical and ordinal models, the output is an S x N x C array.
#'   Otherwise, the output is an S x N matrix, where S is the number of
#'   posterior draws, N is the number of observations, and C is the number of
#'   categories. In multivariate models, an additional dimension is added to the
#'   output which indexes along the different response variables.
#'
#' @template details-newdata-na
#' @template details-allow_new_levels
#'
#' @examples
#' \dontrun{
#' ## fit a model
#' fit <- brm(rating ~ treat + period + carry + (1|subject),
#'            data = inhaler)
#'
#' ## compute expected predictions
#' ppe <- posterior_epred(fit)
#' str(ppe)
#' }
#'
#' @aliases posterior_epred
#' @method posterior_epred brmsfit
#' @importFrom rstantools posterior_epred
#' @export posterior_epred
#' @export
posterior_epred.brmsfit <- function(object, newdata = NULL, re_formula = NULL,
                                    re.form = NULL, resp = NULL, dpar = NULL,
                                    nlpar = NULL, ndraws = NULL, draw_ids = NULL,
                                    sort = FALSE, ...) {
  cl <- match.call()
  if ("re.form" %in% names(cl) && !missing(re.form)) {
    re_formula <- re.form
  }
  contains_draws(object)
  object <- restructure(object)
  prep <- prepare_predictions(
    object, newdata = newdata, re_formula = re_formula, resp = resp,
    ndraws = ndraws, draw_ids = draw_ids, check_response = FALSE, ...
  )
  posterior_epred(
    prep,  dpar = dpar, nlpar = nlpar, sort = sort,
    scale = "response", summary = FALSE
  )
}

#' @export
posterior_epred.mvbrmsprep <- function(object, ...) {
  out <- lapply(object$resps, posterior_epred, ...)
  along <- ifelse(length(out) > 1L, 3, 2)
  do_call(abind, c(out, along = along))
}

#' @export
posterior_epred.brmsprep <- function(object, dpar, nlpar, sort,
                                     scale = "response", incl_thres = NULL,
                                     summary = FALSE, robust = FALSE,
                                     probs = c(0.025, 0.975), ...) {
  summary <- as_one_logical(summary)
  dpars <- names(object$dpars)
  nlpars <- names(object$nlpars)
  if (length(dpar)) {
    # predict a distributional parameter
    dpar <- as_one_character(dpar)
    if (!dpar %in% dpars) {
      stop2("Invalid argument 'dpar'. Valid distributional ",
            "parameters are: ", collapse_comma(dpars))
    }
    if (length(nlpar)) {
      stop2("Cannot use 'dpar' and 'nlpar' at the same time.")
    }
    predicted <- is.bprepl(object$dpars[[dpar]]) ||
      is.bprepnl(object$dpars[[dpar]])
    if (predicted) {
      # parameter varies across observations
      if (scale == "linear") {
        object$dpars[[dpar]]$family$link <- "identity"
      }
      if (is_ordinal(object$family)) {
        object$dpars[[dpar]]$cs <- NULL
        object$family <- object$dpars[[dpar]]$family <-
          .dpar_family(link = object$dpars[[dpar]]$family$link)
      }
      if (dpar_class(dpar) == "theta" && scale == "response") {
        ap_id <- as.numeric(dpar_id(dpar))
        out <- get_theta(object)[, , ap_id, drop = FALSE]
        dim(out) <- dim(out)[c(1, 2)]
      } else {
        out <- get_dpar(object, dpar = dpar, inv_link = TRUE)
      }
    } else {
      # parameter is constant across observations
      out <- object$dpars[[dpar]]
      out <- matrix(out, nrow = object$ndraws, ncol = object$nobs)
    }
  } else if (length(nlpar)) {
    # predict a non-linear parameter
    nlpar <- as_one_character(nlpar)
    if (!nlpar %in% nlpars) {
      stop2("Invalid argument 'nlpar'. Valid non-linear ",
            "parameters are: ", collapse_comma(nlpars))
    }
    out <- get_nlpar(object, nlpar = nlpar)
  } else {
    # no dpar or nlpar specified
    incl_thres <- as_one_logical(incl_thres %||% FALSE)
    incl_thres <- incl_thres && is_ordinal(object$family) && scale == "linear"
    if (incl_thres) {
      # extract linear predictor array with thresholds etc. included
      if (is.mixfamily(object$family)) {
        stop2("'incl_thres' is not supported for mixture models.")
      }
      object$family$link <- "identity"
    }
    if (scale == "response" || incl_thres) {
      # predict the mean of the response distribution
      for (nlp in nlpars) {
        object$nlpars[[nlp]] <- get_nlpar(object, nlpar = nlp)
      }
      for (dp in dpars) {
        object$dpars[[dp]] <- get_dpar(object, dpar = dp)
      }
      if (is_trunc(object)) {
        out <- posterior_epred_trunc(object)
      } else {
        posterior_epred_fun <- paste0("posterior_epred_", object$family$family)
        posterior_epred_fun <- get(posterior_epred_fun, asNamespace("brms"))
        out <- posterior_epred_fun(object)
      }
    } else {
      # return results on the linear scale
      # extract all 'mu' parameters
      if (conv_cats_dpars(object$family)) {
        out <- dpars[grepl("^mu", dpars)]
      } else {
        out <- dpars[dpar_class(dpars) %in% "mu"]
      }
      if (length(out) == 1L) {
        out <- get_dpar(object, dpar = out, inv_link = FALSE)
      } else {
        # multiple mu parameters in categorical or mixture models
        out <- lapply(out, get_dpar, prep = object, inv_link = FALSE)
        out <- abind::abind(out, along = 3)
      }
    }
  }
  if (is.null(dim(out))) {
    out <- as.matrix(out)
  }
  colnames(out) <- NULL
  out <- reorder_obs(out, object$old_order, sort = sort)
  if (scale == "response" && is_polytomous(object$family) &&
      length(dim(out)) == 3L && dim(out)[3] == length(object$cats)) {
    # for ordinal models with varying thresholds, dim[3] may not match cats
    dimnames(out)[[3]] <- object$cats
  }
  if (summary) {
    # only for compatibility with the 'fitted' method
    out <- posterior_summary(out, probs = probs, robust = robust)
    if (is_polytomous(object$family) && length(dim(out)) == 3L) {
      if (scale == "linear") {
        dimnames(out)[[3]] <- paste0("eta", seq_dim(out, 3))
      } else {
        dimnames(out)[[3]] <- paste0("P(Y = ", dimnames(out)[[3]], ")")
      }
    }
  }
  out
}

#' Expected Values of the Posterior Predictive Distribution
#'
#' This method is an alias of \code{\link{posterior_epred.brmsfit}}
#' with additional arguments for obtaining summaries of the computed draws.
#'
#' @inheritParams posterior_epred.brmsfit
#' @param object An object of class \code{brmsfit}.
#' @param scale Either \code{"response"} or \code{"linear"}.
#'  If \code{"response"}, results are returned on the scale
#'  of the response variable. If \code{"linear"},
#'  results are returned on the scale of the linear predictor term,
#'  that is without applying the inverse link function or
#'  other transformations.
#' @param summary Should summary statistics be returned
#'  instead of the raw values? Default is \code{TRUE}..
#' @param robust If \code{FALSE} (the default) the mean is used as
#'  the measure of central tendency and the standard deviation as
#'  the measure of variability. If \code{TRUE}, the median and the
#'  median absolute deviation (MAD) are applied instead.
#'  Only used if \code{summary} is \code{TRUE}.
#' @param probs The percentiles to be computed by the \code{quantile}
#'  function. Only used if \code{summary} is \code{TRUE}.
#'
#' @return An \code{array} of predicted \emph{mean} response values.
#'   If \code{summary = FALSE} the output resembles those of
#'   \code{\link{posterior_epred.brmsfit}}.
#'
#'   If \code{summary = TRUE} the output depends on the family: For categorical
#'   and ordinal families, the output is an N x E x C array, where N is the
#'   number of observations, E is the number of summary statistics, and C is the
#'   number of categories. For all other families, the output is an N x E
#'   matrix. The number of summary statistics E is equal to \code{2 +
#'   length(probs)}: The \code{Estimate} column contains point estimates (either
#'   mean or median depending on argument \code{robust}), while the
#'   \code{Est.Error} column contains uncertainty estimates (either standard
#'   deviation or median absolute deviation depending on argument
#'   \code{robust}). The remaining columns starting with \code{Q} contain
#'   quantile estimates as specified via argument \code{probs}.
#'
#'   In multivariate models, an additional dimension is added to the output
#'   which indexes along the different response variables.
#'
#' @seealso \code{\link{posterior_epred.brmsfit}}
#'
#' @examples
#' \dontrun{
#' ## fit a model
#' fit <- brm(rating ~ treat + period + carry + (1|subject),
#'            data = inhaler)
#'
#' ## compute expected predictions
#' fitted_values <- fitted(fit)
#' head(fitted_values)
#'
#' ## plot expected predictions against actual response
#' dat <- as.data.frame(cbind(Y = standata(fit)$Y, fitted_values))
#' ggplot(dat) + geom_point(aes(x = Estimate, y = Y))
#' }
#'
#' @export
fitted.brmsfit <- function(object, newdata = NULL, re_formula = NULL,
                           scale = c("response", "linear"),
                           resp = NULL, dpar = NULL, nlpar = NULL,
                           ndraws = NULL, draw_ids = NULL, sort = FALSE,
                           summary = TRUE, robust = FALSE,
                           probs = c(0.025, 0.975), ...) {
  scale <- match.arg(scale)
  summary <- as_one_logical(summary)
  contains_draws(object)
  object <- restructure(object)
  prep <- prepare_predictions(
    object, newdata = newdata, re_formula = re_formula, resp = resp,
    ndraws = ndraws, draw_ids = draw_ids, check_response = FALSE, ...
  )
  posterior_epred(
    prep, dpar = dpar, nlpar = nlpar, sort = sort, scale = scale,
    summary = summary, robust = robust, probs = probs
  )
}

#' Posterior Draws of the Linear Predictor
#'
#' Compute posterior draws of the linear predictor, that is draws before
#' applying any link functions or other transformations. Can be performed for
#' the data used to fit the model (posterior predictive checks) or for new data.
#'
#' @inheritParams posterior_epred.brmsfit
#' @param object An object of class \code{brmsfit}.
#' @param transform Logical; if \code{FALSE}
#'  (the default), draws of the linear predictor are returned.
#'  If \code{TRUE}, draws of the transformed linear predictor,
#'  that is, after applying the inverse link function are returned.
#' @param dpar Name of a predicted distributional parameter
#'  for which draws are to be returned. By default, draws
#'  of the main distributional parameter(s) \code{"mu"} are returned.
#' @param incl_thres Logical; only relevant for ordinal models when
#'   \code{transform} is \code{FALSE}, and ignored otherwise. Shall the
#'   thresholds and category-specific effects be included in the linear
#'   predictor? For backwards compatibility, the default is to not include them.
#'
#' @seealso \code{\link{posterior_epred.brmsfit}}
#'
#' @examples
#' \dontrun{
#' ## fit a model
#' fit <- brm(rating ~ treat + period + carry + (1|subject),
#'            data = inhaler)
#'
#' ## extract linear predictor values
#' pl <- posterior_linpred(fit)
#' str(pl)
#' }
#'
#' @aliases posterior_linpred
#' @method posterior_linpred brmsfit
#' @importFrom rstantools posterior_linpred
#' @export
#' @export posterior_linpred
posterior_linpred.brmsfit <- function(
  object, transform = FALSE, newdata = NULL, re_formula = NULL,
  re.form = NULL, resp = NULL, dpar = NULL, nlpar = NULL,
  incl_thres = NULL, ndraws = NULL, draw_ids = NULL, sort = FALSE, ...
) {
  cl <- match.call()
  if ("re.form" %in% names(cl) && !missing(re.form)) {
    re_formula <- re.form
  }
  scale <- "linear"
  transform <- as_one_logical(transform)
  if (transform) {
    scale <- "response"
    # if transform, return inv-link draws of only a single
    # distributional or non-linear parameter for consistency
    # of brms and rstanarm
    if (is.null(dpar) && is.null(nlpar)) {
      dpar <- "mu"
    }
  }
  contains_draws(object)
  object <- restructure(object)
  prep <- prepare_predictions(
    object, newdata = newdata, re_formula = re_formula, resp = resp,
    ndraws = ndraws, draw_ids = draw_ids, check_response = FALSE, ...
  )
  posterior_epred(
    prep, dpar = dpar, nlpar = nlpar, sort = sort,
    scale = scale, incl_thres = incl_thres, summary = FALSE
  )
}

# ------------------- family specific posterior_epred methods ---------------------
# All posterior_epred_<family> functions have the same arguments structure
# @param prep A named list returned by prepare_predictions containing
#   all required data and draws
# @return transformed linear predictor representing the mean
#   of the posterior predictive distribution
posterior_epred_gaussian <- function(prep) {
  if (!is.null(prep$ac$lagsar)) {
    prep$dpars$mu <- posterior_epred_lagsar(prep)
  }
  prep$dpars$mu
}

posterior_epred_student <- function(prep) {
  if (!is.null(prep$ac$lagsar)) {
    prep$dpars$mu <- posterior_epred_lagsar(prep)
  }
  prep$dpars$mu
}

posterior_epred_skew_normal <- function(prep) {
  prep$dpars$mu
}

posterior_epred_lognormal <- function(prep) {
  with(prep$dpars, exp(mu + sigma^2 / 2))
}

posterior_epred_shifted_lognormal <- function(prep) {
  with(prep$dpars, exp(mu + sigma^2 / 2) + ndt)
}

posterior_epred_binomial <- function(prep) {
  trials <- data2draws(prep$data$trials, dim_mu(prep))
  prep$dpars$mu * trials
}

posterior_epred_beta_binomial <- function(prep) {
  # beta part included in mu
  trials <- data2draws(prep$data$trials, dim_mu(prep))
  prep$dpars$mu * trials
}

posterior_epred_bernoulli <- function(prep) {
  prep$dpars$mu
}

posterior_epred_poisson <- function(prep) {
  multiply_dpar_rate_denom(prep$dpars$mu, prep)
}

posterior_epred_negbinomial <- function(prep) {
  multiply_dpar_rate_denom(prep$dpars$mu, prep)
}

posterior_epred_negbinomial2 <- function(prep) {
  multiply_dpar_rate_denom(prep$dpars$mu, prep)
}

posterior_epred_geometric <- function(prep) {
  multiply_dpar_rate_denom(prep$dpars$mu, prep)
}

posterior_epred_discrete_weibull <- function(prep) {
  mean_discrete_weibull(prep$dpars$mu, prep$dpars$shape)
}

posterior_epred_com_poisson <- function(prep) {
  mean_com_poisson(prep$dpars$mu, prep$dpars$shape)
}

posterior_epred_exponential <- function(prep) {
  prep$dpars$mu
}

posterior_epred_gamma <- function(prep) {
  prep$dpars$mu
}

posterior_epred_weibull <- function(prep) {
  prep$dpars$mu
}

posterior_epred_frechet <- function(prep) {
  prep$dpars$mu
}

posterior_epred_gen_extreme_value <- function(prep) {
  with(prep$dpars, mu + sigma * (gamma(1 - xi) - 1) / xi)
}

posterior_epred_loglogistic <- function(prep) {
  if (prep$dpar$shape <= 1) {
    stop2("The mean is undefined for shape parameter <= 1.")
  }
  with(prep$dpars, scale * (pi / shape) / sin(pi / shape))
}

posterior_epred_inverse.gaussian <- function(prep) {
  prep$dpars$mu
}

posterior_epred_exgaussian <- function(prep) {
  prep$dpars$mu
}

posterior_epred_wiener <- function(prep) {
  # obtained from https://doi.org/10.1016/j.jmp.2009.01.006
  # mu is the drift rate
  with(prep$dpars,
   ndt - bias / mu + bs / mu *
     (exp(-2 * mu * bias) - 1) / (exp(-2 * mu * bs) - 1)
  )
}

posterior_epred_beta <- function(prep) {
  prep$dpars$mu
}

posterior_epred_von_mises <- function(prep) {
  prep$dpars$mu
}

posterior_epred_asym_laplace <- function(prep) {
  with(prep$dpars,
    mu + sigma * (1 - 2 * quantile) / (quantile * (1 - quantile))
  )
}

posterior_epred_zero_inflated_asym_laplace <- function(prep) {
  posterior_epred_asym_laplace(prep) * (1 - prep$dpars$zi)
}

posterior_epred_cox <- function(prep) {
  stop2("Cannot compute expected values of the posterior predictive ",
        "distribution for family 'cox'.")
}

posterior_epred_hurdle_poisson <- function(prep) {
  with(prep$dpars, mu / (1 - exp(-mu)) * (1 - hu))
}

posterior_epred_hurdle_negbinomial <- function(prep) {
  with(prep$dpars, mu / (1 - (shape / (mu + shape))^shape) * (1 - hu))
}

posterior_epred_hurdle_gamma <- function(prep) {
  with(prep$dpars, mu * (1 - hu))
}

posterior_epred_hurdle_lognormal <- function(prep) {
  with(prep$dpars, exp(mu + sigma^2 / 2) * (1 - hu))
}

posterior_epred_mixcure_lognormal <- function(prep) {
  stop2("Cannot compute expected values of the posterior predictive ",
        "distribution for family 'mixcure_lognormal' because some values ",
        "of the response are predicted to be infinite.")
}

posterior_epred_mixcure_weibull <- function(prep) {
  stop2("Cannot compute expected values of the posterior predictive ",
        "distribution for family 'mixcure_weibull' because some values ",
        "of the response are predicted to be infinite.")
}

posterior_epred_mixcure_loglogistic <- function(prep) {
  stop2("Cannot compute expected values of the posterior predictive ",
        "distribution for family 'mixcure_loglogistic' because some values ",
        "of the response are predicted to be infinite.")
}

posterior_epred_hurdle_cumulative <- function(prep) {
  adjust <- ifelse(prep$family$link == "identity", 0, 1)
  ncat_max <- max(prep$data$nthres) + adjust
  nact_min <- min(prep$data$nthres) + adjust
  init_mat <- matrix(
    ifelse(prep$family$link == "identity", NA, 0),
    nrow = prep$ndraws, ncol = ncat_max - nact_min
  )
  args <- list(link = prep$family$link)
  out <- vector("list", prep$nobs)

  for (i in seq_along(out)) {
    args_i <- args
    args_i$eta <- slice_col(get_dpar(prep, "mu", i))
    args_i$disc <- slice_col(get_dpar(prep, "disc", i))
    args_i$thres <- subset_thres(prep, i)
    ncat_i <- NCOL(args_i$thres) + adjust
    args_i$x <- seq_len(ncat_i)
    out[[i]] <- do_call(dcumulative, args_i)
    if (ncat_i < ncat_max) {
      sel <- seq_len(ncat_max - ncat_i)
      out[[i]] <- cbind(out[[i]], init_mat[, sel])
    }
    hu <- get_dpar(prep, "hu", i)
    out[[i]] <- cbind(hu, out[[i]] * (1 - hu))
  }
  out <- abind(out, along = 3)
  out <- aperm(out, perm = c(1, 3, 2))
  dimnames(out)[[3]] <- c(paste0(0), seq_len(ncat_max))
  out
}

posterior_epred_zero_inflated_poisson <- function(prep) {
  with(prep$dpars, mu * (1 - zi))
}

posterior_epred_zero_inflated_negbinomial <- function(prep) {
  with(prep$dpars, mu * (1 - zi))
}

posterior_epred_zero_inflated_binomial <- function(prep) {
  trials <- data2draws(prep$data$trials, dim_mu(prep))
  prep$dpars$mu * trials * (1 - prep$dpars$zi)
}

posterior_epred_zero_inflated_beta_binomial <- function(prep) {
  # same as zero_inflated_binom as beta part is included in mu
  trials <- data2draws(prep$data$trials, dim_mu(prep))
  prep$dpars$mu * trials * (1 - prep$dpars$zi)
}

posterior_epred_zero_inflated_beta <- function(prep) {
  with(prep$dpars, mu * (1 - zi))
}

posterior_epred_zero_one_inflated_beta <- function(prep) {
  with(prep$dpars, zoi * coi + mu * (1 - zoi))
}

posterior_epred_categorical <- function(prep) {
  get_probs <- function(i) {
    eta <- insert_refcat(slice_col(eta, i), refcat = prep$refcat)
    dcategorical(cats, eta = eta)
  }
  eta <- get_Mu(prep)
  cats <- seq_len(prep$data$ncat)
  out <- abind(lapply(seq_cols(eta), get_probs), along = 3)
  out <- aperm(out, perm = c(1, 3, 2))
  dimnames(out)[[3]] <- prep$cats
  out
}

posterior_epred_multinomial <- function(prep) {
  get_counts <- function(i) {
    eta <- insert_refcat(slice_col(eta, i), refcat = prep$refcat)
    dcategorical(cats, eta = eta) * trials[i]
  }
  eta <- get_Mu(prep)
  cats <- seq_len(prep$data$ncat)
  trials <- prep$data$trials
  out <- abind(lapply(seq_cols(eta), get_counts), along = 3)
  out <- aperm(out, perm = c(1, 3, 2))
  dimnames(out)[[3]] <- prep$cats
  out
}

posterior_epred_dirichlet <- function(prep) {
  get_probs <- function(i) {
    eta <- insert_refcat(slice_col(eta, i), refcat = prep$refcat)
    dcategorical(cats, eta = eta)
  }
  eta <- get_Mu(prep)
  cats <- seq_len(prep$data$ncat)
  out <- abind(lapply(seq_cols(eta), get_probs), along = 3)
  out <- aperm(out, perm = c(1, 3, 2))
  dimnames(out)[[3]] <- prep$cats
  out
}

posterior_epred_dirichlet2 <- function(prep) {
  mu <- get_Mu(prep)
  sums_mu <- apply(mu, 1:2, sum)
  cats <- seq_len(prep$data$ncat)
  for (i in cats) {
    mu[, , i] <- mu[, , i] / sums_mu
  }
  dimnames(mu)[[3]] <- prep$cats
  mu
}

posterior_epred_logistic_normal <- function(prep) {
  stop2("Cannot compute expected values of the posterior predictive ",
        "distribution for family 'logistic_normal'.")
}

posterior_epred_cumulative <- function(prep) {
  posterior_epred_ordinal(prep)
}

posterior_epred_sratio <- function(prep) {
  posterior_epred_ordinal(prep)
}

posterior_epred_cratio <- function(prep) {
  posterior_epred_ordinal(prep)
}

posterior_epred_acat <- function(prep) {
  posterior_epred_ordinal(prep)
}

posterior_epred_custom <- function(prep) {
  custom_family_method(prep$family, "posterior_epred")(prep)
}

posterior_epred_mixture <- function(prep) {
  families <- family_names(prep$family)
  prep$dpars$theta <- get_theta(prep)
  out <- 0
  for (j in seq_along(families)) {
    posterior_epred_fun <- paste0("posterior_epred_", families[j])
    posterior_epred_fun <- get(posterior_epred_fun, asNamespace("brms"))
    tmp_prep <- pseudo_prep_for_mixture(prep, j)
    if (length(dim(prep$dpars$theta)) == 3L) {
      theta <- prep$dpars$theta[, , j]
    } else {
      theta <- prep$dpars$theta[, j]
    }
    out <- out + theta * posterior_epred_fun(tmp_prep)
  }
  out
}

# ------ posterior_epred helper functions ------
# compute 'posterior_epred' for ordinal models
posterior_epred_ordinal <- function(prep) {
  dens <- get(paste0("d", prep$family$family), mode = "function")
  # the linear scale has one column less than the response scale
  adjust <- ifelse(prep$family$link == "identity", 0, 1)
  ncat_max <- max(prep$data$nthres) + adjust
  nact_min <- min(prep$data$nthres) + adjust
  init_mat <- matrix(ifelse(prep$family$link == "identity", NA, 0),
                     nrow = prep$ndraws,
                     ncol = ncat_max - nact_min)
  args <- list(link = prep$family$link)
  out <- vector("list", prep$nobs)
  for (i in seq_along(out)) {
    args_i <- args
    args_i$eta <- slice_col(prep$dpars$mu, i)
    args_i$disc <- slice_col(prep$dpars$disc, i)
    args_i$thres <- subset_thres(prep, i)
    ncat_i <- NCOL(args_i$thres) + adjust
    args_i$x <- seq_len(ncat_i)
    out[[i]] <- do_call(dens, args_i)
    if (ncat_i < ncat_max) {
      sel <- seq_len(ncat_max - ncat_i)
      out[[i]] <- cbind(out[[i]], init_mat[, sel])
    }
  }
  out <- abind(out, along = 3)
  out <- aperm(out, perm = c(1, 3, 2))
  dimnames(out)[[3]] <- seq_len(ncat_max)
  out
}

# compute 'posterior_epred' for lagsar models
posterior_epred_lagsar <- function(prep) {
  stopifnot(!is.null(prep$ac$lagsar))
  I <- diag(prep$nobs)
  .posterior_epred <- function(s) {
    IB <- I - with(prep$ac, lagsar[s, ] * Msar)
    as.numeric(solve(IB, prep$dpars$mu[s, ]))
  }
  out <- rblapply(seq_len(prep$ndraws), .posterior_epred)
  rownames(out) <- NULL
  out
}

# expand data to dimension appropriate for
# vectorized multiplication with posterior draws
data2draws <- function(x, dim) {
  stopifnot(length(dim) == 2L, length(x) %in% c(1, dim[2]))
  matrix(x, nrow = dim[1], ncol = dim[2], byrow = TRUE)
}

# expected dimension of the main parameter 'mu'
dim_mu <- function(prep) {
  c(prep$ndraws, prep$nobs)
}

# is the model truncated?
is_trunc <- function(prep) {
  stopifnot(is.brmsprep(prep))
  any(prep$data[["lb"]] > -Inf) || any(prep$data[["ub"]] < Inf)
}

# prepares data required for truncation and calles the
# family specific truncation function for posterior_epred values
posterior_epred_trunc <- function(prep) {
  stopifnot(is_trunc(prep))
  lb <- data2draws(prep$data[["lb"]], dim_mu(prep))
  ub <- data2draws(prep$data[["ub"]], dim_mu(prep))
  posterior_epred_trunc_fun <- paste0("posterior_epred_trunc_", prep$family$family)
  posterior_epred_trunc_fun <- try(
    get(posterior_epred_trunc_fun, asNamespace("brms")),
    silent = TRUE
  )
  if (is(posterior_epred_trunc_fun, "try-error")) {
    stop2("posterior_epred values on the respone scale not yet implemented ",
          "for truncated '", prep$family$family, "' models.")
  }
  trunc_args <- nlist(prep, lb, ub)
  do_call(posterior_epred_trunc_fun, trunc_args)
}

# ----- family specific truncation functions -----
# @param prep output of 'prepare_predictions'
# @param lb lower truncation bound
# @param ub upper truncation bound
# @return draws of the truncated mean parameter
posterior_epred_trunc_gaussian <- function(prep, lb, ub) {
  zlb <- (lb - prep$dpars$mu) / prep$dpars$sigma
  zub <- (ub - prep$dpars$mu) / prep$dpars$sigma
  # truncated mean of standard normal; see Wikipedia
  trunc_zmean <- (dnorm(zlb) - dnorm(zub)) / (pnorm(zub) - pnorm(zlb))
  prep$dpars$mu + trunc_zmean * prep$dpars$sigma
}

posterior_epred_trunc_student <- function(prep, lb, ub) {
  zlb <- with(prep$dpars, (lb - mu) / sigma)
  zub <- with(prep$dpars, (ub - mu) / sigma)
  nu <- prep$dpars$nu
  # see Kim 2008: Moments of truncated Student-t distribution
  G1 <- gamma((nu - 1) / 2) * nu^(nu / 2) /
    (2 * (pt(zub, df = nu) - pt(zlb, df = nu))
     * gamma(nu / 2) * gamma(0.5))
  A <- (nu + zlb^2) ^ (-(nu - 1) / 2)
  B <- (nu + zub^2) ^ (-(nu - 1) / 2)
  trunc_zmean <- G1 * (A - B)
  prep$dpars$mu + trunc_zmean * prep$dpars$sigma
}

posterior_epred_trunc_lognormal <- function(prep, lb, ub) {
  lb <- ifelse(lb < 0, 0, lb)
  m1 <- with(prep$dpars,
    exp(mu + sigma^2 / 2) *
      (pnorm((log(ub) - mu) / sigma - sigma) -
       pnorm((log(lb) - mu) / sigma - sigma))
  )
  with(prep$dpars,
    m1 / (plnorm(ub, meanlog = mu, sdlog = sigma) -
          plnorm(lb, meanlog = mu, sdlog = sigma))
  )
}

posterior_epred_trunc_gamma <- function(prep, lb, ub) {
  lb <- ifelse(lb < 0, 0, lb)
  prep$dpars$scale <- prep$dpars$mu / prep$dpars$shape
  # see Jawitz 2004: Moments of truncated continuous univariate distributions
  m1 <- with(prep$dpars,
    scale / gamma(shape) *
      (incgamma(1 + shape, ub / scale) -
       incgamma(1 + shape, lb / scale))
  )
  with(prep$dpars,
    m1 / (pgamma(ub, shape, scale = scale) -
          pgamma(lb, shape, scale = scale))
  )
}

posterior_epred_trunc_exponential <- function(prep, lb, ub) {
  lb <- ifelse(lb < 0, 0, lb)
  inv_mu <- 1 / prep$dpars$mu
  # see Jawitz 2004: Moments of truncated continuous univariate distributions
  m1 <- with(prep$dpars, mu * (incgamma(2, ub / mu) - incgamma(2, lb / mu)))
  m1 / (pexp(ub, rate = inv_mu) - pexp(lb, rate = inv_mu))
}

posterior_epred_trunc_weibull <- function(prep, lb, ub) {
  lb <- ifelse(lb < 0, 0, lb)
  prep$dpars$a <- 1 + 1 / prep$dpars$shape
  prep$dpars$scale <- with(prep$dpars, mu / gamma(a))
  # see Jawitz 2004: Moments of truncated continuous univariate distributions
  m1 <- with(prep$dpars,
    scale * (incgamma(a, (ub / scale)^shape) -
             incgamma(a, (lb / scale)^shape))
  )
  with(prep$dpars,
    m1 / (pweibull(ub, shape, scale = scale) -
          pweibull(lb, shape, scale = scale))
  )
}

posterior_epred_trunc_binomial <- function(prep, lb, ub) {
  lb <- ifelse(lb < -1, -1, lb)
  max_value <- max(prep$data$trials)
  ub <- ifelse(ub > max_value, max_value, ub)
  trials <- prep$data$trials
  if (length(trials) > 1) {
    trials <- data2draws(trials, dim_mu(prep))
  }
  args <- list(size = trials, prob = prep$dpars$mu)
  posterior_epred_trunc_discrete(dist = "binom", args = args, lb = lb, ub = ub)
}

posterior_epred_trunc_poisson <- function(prep, lb, ub) {
  lb <- ifelse(lb < -1, -1, lb)
  mu <- multiply_dpar_rate_denom(prep$dpars$mu, prep)
  max_value <- 3 * max(mu)
  ub <- ifelse(ub > max_value, max_value, ub)
  args <- list(lambda = mu)
  posterior_epred_trunc_discrete(dist = "pois", args = args, lb = lb, ub = ub)
}

posterior_epred_trunc_negbinomial <- function(prep, lb, ub) {
  lb <- ifelse(lb < -1, -1, lb)
  mu <- multiply_dpar_rate_denom(prep$dpars$mu, prep)
  max_value <- 3 * max(mu)
  ub <- ifelse(ub > max_value, max_value, ub)
  shape <- multiply_dpar_rate_denom(prep$dpars$shape, prep)
  args <- list(mu = mu, size = shape)
  posterior_epred_trunc_discrete(dist = "nbinom", args = args, lb = lb, ub = ub)
}

posterior_epred_trunc_negbinomial2 <- function(prep, lb, ub) {
  lb <- ifelse(lb < -1, -1, lb)
  mu <- multiply_dpar_rate_denom(prep$dpars$mu, prep)
  max_value <- 3 * max(mu)
  ub <- ifelse(ub > max_value, max_value, ub)
  shape <- multiply_dpar_rate_denom(1 / prep$dpars$sigma, prep)
  args <- list(mu = mu, size = shape)
  posterior_epred_trunc_discrete(dist = "nbinom", args = args, lb = lb, ub = ub)
}

posterior_epred_trunc_geometric <- function(prep, lb, ub) {
  lb <- ifelse(lb < -1, -1, lb)
  mu <- multiply_dpar_rate_denom(prep$dpars$mu, prep)
  max_value <- 3 * max(mu)
  ub <- ifelse(ub > max_value, max_value, ub)
  shape <- multiply_dpar_rate_denom(1, prep)
  args <- list(mu = mu, size = shape)
  posterior_epred_trunc_discrete(dist = "nbinom", args = args, lb = lb, ub = ub)
}

# posterior_epred values for truncated discrete distributions
posterior_epred_trunc_discrete <- function(dist, args, lb, ub) {
  stopifnot(is.matrix(lb), is.matrix(ub))
  message(
    "Computing posterior_epred values for truncated ",
    "discrete models may take a while."
  )
  pdf <- get(paste0("d", dist), mode = "function")
  cdf <- get(paste0("p", dist), mode = "function")
  mean_kernel <- function(x, args) {
    # just x * density(x)
    x * do_call(pdf, c(x, args))
  }
  if (any(is.infinite(c(lb, ub)))) {
    stop("lb and ub must be finite")
  }
  # simplify lb and ub back to vector format
  vec_lb <- lb[1, ]
  vec_ub <- ub[1, ]
  min_lb <- min(vec_lb)
  # array of dimension S x N x length((lb+1):ub)
  mk <- lapply((min_lb + 1):max(vec_ub), mean_kernel, args = args)
  mk <- do_call(abind, c(mk, along = 3))
  m1 <- vector("list", ncol(mk))
  for (n in seq_along(m1)) {
    # summarize only over non-truncated values for this observation
    J <- (vec_lb[n] - min_lb + 1):(vec_ub[n] - min_lb)
    m1[[n]] <- rowSums(mk[, n, ][, J, drop = FALSE])
  }
  rm(mk)
  m1 <- do.call(cbind, m1)
  m1 / (do.call(cdf, c(list(ub), args)) - do.call(cdf, c(list(lb), args)))
}

#' @export
pp_expect <- function(object, ...) {
  warning2("Method 'pp_expect' is deprecated. ",
           "Please use 'posterior_epred' instead.")
  UseMethod("posterior_epred")
}
