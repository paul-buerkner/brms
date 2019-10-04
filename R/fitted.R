fitted_internal <- function(draws, ...) {
  UseMethod("fitted_internal")
}

#' @export
fitted_internal.mvbrmsdraws <- function(draws, ...) {
  out <- lapply(draws$resps, fitted_internal, ...)
  along <- ifelse(length(out) > 1L, 3, 2)
  do_call(abind, c(out, along = along))
}

#' @export
fitted_internal.brmsdraws <- function(
  draws, scale = "response", dpar = NULL, nlpar = NULL,
  summary = TRUE, sort = FALSE, robust = FALSE, 
  probs = c(0.025, 0.975), ...
) {
  dpars <- names(draws$dpars)
  nlpars <- names(draws$nlpars)
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
    predicted <- is.bdrawsl(draws$dpars[[dpar]]) ||
      is.bdrawsnl(draws$dpars[[dpar]])
    if (predicted) {
      # parameter varies across observations
      if (scale == "linear") {
        draws$dpars[[dpar]]$family$link <- "identity"
      }
      if (is_ordinal(draws$family)) {
        draws$dpars[[dpar]]$cs <- NULL
        draws$family <- draws$dpars[[dpar]]$family <- 
          .dpar_family(link = draws$dpars[[dpar]]$family$link)
      }
      if (dpar_class(dpar) == "theta" && scale == "response") {
        ap_id <- as.numeric(dpar_id(dpar))
        out <- get_theta(draws)[, , ap_id, drop = FALSE]
        dim(out) <- dim(out)[c(1, 2)]
      } else {
        out <- get_dpar(draws, dpar = dpar, ilink = TRUE)
      }
    } else {
      # parameter is constant across observations
      out <- draws$dpars[[dpar]]
      out <- matrix(out, nrow = draws$nsamples, ncol = draws$nobs)
    }
  } else if (length(nlpar)) {
    # predict a non-linear parameter
    nlpar <- as_one_character(nlpar)
    if (!nlpar %in% nlpars) {
      stop2("Invalid argument 'nlpar'. Valid non-linear ",
            "parameters are: ", collapse_comma(nlpars))
    }
    out <- get_nlpar(draws, nlpar = nlpar)
  } else {
    # predict the mean of the response distribution
    if (scale == "response") {
      for (nlp in nlpars) {
        draws$nlpars[[nlp]] <- get_nlpar(draws, nlpar = nlp)
      }
      for (dp in dpars) {
        draws$dpars[[dp]] <- get_dpar(draws, dpar = dp)
      }
      if (is_trunc(draws)) {
        out <- fitted_trunc(draws)
      } else {
        fitted_fun <- paste0("fitted_", draws$family$family)
        fitted_fun <- get(fitted_fun, asNamespace("brms"))
        out <- fitted_fun(draws)
      }
    } else {
      if (conv_cats_dpars(draws$family)) {
        mus <- dpars[grepl("^mu", dpars)] 
      } else {
        mus <- dpars[dpar_class(dpars) %in% "mu"]
      }
      if (length(mus) == 1L) {
        out <- get_dpar(draws, dpar = mus, ilink = FALSE)
      } else {
        # multiple mu parameters in categorical or mixture models
        out <- lapply(mus, get_dpar, draws = draws, ilink = FALSE)
        out <- abind::abind(out, along = 3)
      }
    }
  }
  if (is.null(dim(out))) {
    out <- as.matrix(out)
  }
  colnames(out) <- NULL
  out <- reorder_obs(out, draws$old_order, sort = sort)
  if (summary) {
    out <- posterior_summary(out, probs = probs, robust = robust)
    if (has_cat(draws$family) && length(dim(out)) == 3L) {
      if (scale == "linear") {
        dimnames(out)[[3]] <- paste0("eta", seq_dim(out, 3))
      } else {
        dimnames(out)[[3]] <- paste0("P(Y = ", dimnames(out)[[3]], ")")
      }
    }
  }
  out
}

# All fitted_<family> functions have the same arguments structure
# @param draws A named list returned by extract_draws containing 
#   all required data and samples
# @return transformed linear predictor representing the mean
#   of the response distribution
fitted_gaussian <- function(draws) {
  if (!is.null(draws$ac$lagsar)) {
    draws$dpars$mu <- fitted_lagsar(draws)
  }
  draws$dpars$mu
}

fitted_student <- function(draws) {
  if (!is.null(draws$ac$lagsar)) {
    draws$dpars$mu <- fitted_lagsar(draws)
  }
  draws$dpars$mu
}

fitted_skew_normal <- function(draws) {
  draws$dpars$mu
}

fitted_lognormal <- function(draws) {
  with(draws$dpars, exp(mu + sigma^2 / 2))
}

fitted_shifted_lognormal <- function(draws) {
  with(draws$dpars, exp(mu + sigma^2 / 2) + ndt)
}

fitted_binomial <- function(draws) {
  trials <- as_draws_matrix(draws$data$trials, dim_mu(draws))
  draws$dpars$mu * trials 
}

fitted_bernoulli <- function(draws) {
  draws$dpars$mu
}

fitted_poisson <- function(draws) {
  draws$dpars$mu
}

fitted_negbinomial <- function(draws) {
  draws$dpars$mu
}

fitted_geometric <- function(draws) {
  draws$dpars$mu
}

fitted_discrete_weibull <- function(draws) {
  mean_discrete_weibull(draws$dpars$mu, draws$dpars$shape)
}

fitted_com_poisson <- function(draws) {
  mean_com_poisson(draws$dpars$mu, draws$dpars$shape)
}

fitted_exponential <- function(draws) {
  draws$dpars$mu
}

fitted_gamma <- function(draws) {
  draws$dpars$mu
}

fitted_weibull <- function(draws) {
  draws$dpars$mu
}

fitted_frechet <- function(draws) {
  draws$dpars$mu
}

fitted_gen_extreme_value <- function(draws) {
  with(draws$dpars, mu + sigma * (gamma(1 - xi) - 1) / xi)
}

fitted_inverse.gaussian <- function(draws) {
  draws$dpars$mu
}

fitted_exgaussian <- function(draws) {
  draws$dpars$mu
}

fitted_wiener <- function(draws) {
  # mu is the drift rate
  with(draws$dpars,
   ndt - bias / mu + bs / mu * 
     (exp(-2 * mu * bias) - 1) / (exp(-2 * mu * bs) - 1)
  )
}

fitted_beta <- function(draws) {
  draws$dpars$mu
}

fitted_von_mises <- function(draws) {
  draws$dpars$mu
}

fitted_asym_laplace <- function(draws) {
  with(draws$dpars, 
    mu + sigma * (1 - 2 * quantile) / (quantile * (1 - quantile))
  )
}

fitted_zero_inflated_asym_laplace <- function(draws) {
  fitted_asym_laplace(draws) * (1 - draws$dpars$zi)
}

fitted_cox <- function(draws) {
  stop2("Cannot compute expected values of the posterior predictive ",
        "distribution for family 'cox'.")
}

fitted_hurdle_poisson <- function(draws) {
  with(draws$dpars, mu / (1 - exp(-mu)) * (1 - hu))
}

fitted_hurdle_negbinomial <- function(draws) {
  with(draws$dpars, mu / (1 - (shape / (mu + shape))^shape) * (1 - hu))
}

fitted_hurdle_gamma <- function(draws) {
  with(draws$dpars, mu * (1 - hu))
}

fitted_hurdle_lognormal <- function(draws) {
  with(draws$dpars, exp(mu + sigma^2 / 2) * (1 - hu))
}

fitted_zero_inflated_poisson <- function(draws) {
  with(draws$dpars, mu * (1 - zi))
}

fitted_zero_inflated_negbinomial <- function(draws) {
  with(draws$dpars, mu * (1 - zi))  
}

fitted_zero_inflated_binomial <- function(draws) {
  trials <- as_draws_matrix(draws$data$trials, dim_mu(draws))
  draws$dpars$mu * trials * (1 - draws$dpars$zi)
}

fitted_zero_inflated_beta <- function(draws) {
  with(draws$dpars, mu * (1 - zi)) 
}

fitted_zero_one_inflated_beta <- function(draws) {
  with(draws$dpars, zoi * coi + mu * (1 - zoi))
}

fitted_categorical <- function(draws) {
  get_probs <- function(i) {
    eta <- insert_refcat(extract_col(eta, i), family = draws$family)
    dcategorical(cats, eta = eta)
  }
  eta <- abind(draws$dpars, along = 3)
  cats <- seq_len(draws$data$ncat)
  out <- abind(lapply(seq_cols(eta), get_probs), along = 3)
  out <- aperm(out, perm = c(1, 3, 2))
  dimnames(out)[[3]] <- draws$data$cats
  out
}

fitted_multinomial <- function(draws) {
  get_counts <- function(i) {
    eta <- insert_refcat(extract_col(eta, i), family = draws$family)
    dcategorical(cats, eta = eta) * trials[i]
  }
  eta <- abind(draws$dpars, along = 3)
  cats <- seq_len(draws$data$ncat)
  trials <- draws$data$trials
  out <- abind(lapply(seq_cols(eta), get_counts), along = 3)
  out <- aperm(out, perm = c(1, 3, 2))
  dimnames(out)[[3]] <- draws$data$cats
  out
}

fitted_dirichlet <- function(draws) {
  get_probs <- function(i) {
    eta <- insert_refcat(extract_col(eta, i), family = draws$family)
    dcategorical(cats, eta = eta)
  }
  eta <- draws$dpars[grepl("^mu", names(draws$dpars))]
  eta <- abind(eta, along = 3)
  cats <- seq_len(draws$data$ncat)
  out <- abind(lapply(seq_cols(eta), get_probs), along = 3)
  out <- aperm(out, perm = c(1, 3, 2))
  dimnames(out)[[3]] <- draws$data$cats
  out
}

fitted_cumulative <- function(draws) {
  fitted_ordinal(draws)
}

fitted_sratio <- function(draws) {
  fitted_ordinal(draws)
}

fitted_cratio <- function(draws) {
  fitted_ordinal(draws)
}

fitted_acat <- function(draws) {
  fitted_ordinal(draws)
}

fitted_custom <- function(draws) {
  fitted_fun <- draws$family$fitted
  if (!is.function(fitted_fun)) {
    fitted_fun <- paste0("fitted_", draws$family$name)
    fitted_fun <- get(fitted_fun, draws$family$env)
  }
  fitted_fun(draws)
}

fitted_mixture <- function(draws) {
  families <- family_names(draws$family)
  draws$dpars$theta <- get_theta(draws)
  out <- 0
  for (j in seq_along(families)) {
    fitted_fun <- paste0("fitted_", families[j])
    fitted_fun <- get(fitted_fun, asNamespace("brms"))
    tmp_draws <- pseudo_draws_for_mixture(draws, j)
    if (length(dim(draws$dpars$theta)) == 3L) {
      theta <- draws$dpars$theta[, , j]
    } else {
      theta <- draws$dpars$theta[, j]
    }
    out <- out + theta * fitted_fun(tmp_draws)
  }
  out
}

# ------ fitted helper functions ------
# compute 'fitted' for ordinal models
fitted_ordinal <- function(draws) {
  dens <- get(paste0("d", draws$family$family), mode = "function")
  args <- list(
    seq_len(draws$data$ncat), 
    thres = draws$thres,
    link = draws$family$link
  )
  out <- vector("list", draws$nobs)
  for (i in seq_along(out)) {
    args_i <- args
    args_i$eta <- extract_col(draws$dpars$mu, i)
    args_i$disc <- extract_col(draws$dpars$disc, i)
    out[[i]] <- do_call(dens, args_i)
  }
  out <- abind(out, along = 3)
  out <- aperm(out, perm = c(1, 3, 2))
  dimnames(out)[[3]] <- draws$data$cats
  out
}

# compute 'fitted' for lagsar models
fitted_lagsar <- function(draws) {
  stopifnot(!is.null(draws$ac$lagsar))
  .fitted <- function(s) {
    W_new <- with(draws, diag(nobs) - ac$lagsar[s, ] * ac$W)
    as.numeric(solve(W_new) %*% draws$dpars$mu[s, ])
  }
  do_call(rbind, lapply(1:draws$nsamples, .fitted))
}

# expand data to dimension appropriate for
# vectorized multiplication with posterior samples
as_draws_matrix <- function(x, dim) {
  stopifnot(length(dim) == 2L, length(x) %in% c(1, dim[2]))
  matrix(x, nrow = dim[1], ncol = dim[2], byrow = TRUE)
}

# expected dimension of the main parameter 'mu'
dim_mu <- function(draws) {
  c(draws$nsamples, draws$nobs)
}

# is the model truncated?
is_trunc <- function(draws) {
  stopifnot(is.brmsdraws(draws))
  any(draws$data[["lb"]] > -Inf) || any(draws$data[["ub"]] < Inf)
}

# prepares data required for truncation and calles the 
# family specific truncation function for fitted values
fitted_trunc <- function(draws) {
  stopifnot(is_trunc(draws))
  lb <- as_draws_matrix(draws$data[["lb"]], dim_mu(draws))
  ub <- as_draws_matrix(draws$data[["ub"]], dim_mu(draws))
  fitted_trunc_fun <- paste0("fitted_trunc_", draws$family$family)
  fitted_trunc_fun <- try(
    get(fitted_trunc_fun, asNamespace("brms")), 
    silent = TRUE
  )
  if (is(fitted_trunc_fun, "try-error")) {
    stop2("Fitted values on the respone scale not yet implemented ",
          "for truncated '", draws$family$family, "' models.")
  }
  trunc_args <- nlist(draws, lb, ub)
  do_call(fitted_trunc_fun, trunc_args)
}

# ----- family specific truncation functions -----
# @param draws output of 'extract_draws'
# @param lb lower truncation bound
# @param ub upper truncation bound
# @return samples of the truncated mean parameter
fitted_trunc_gaussian <- function(draws, lb, ub) {
  zlb <- (lb - draws$dpars$mu) / draws$dpars$sigma
  zub <- (ub - draws$dpars$mu) / draws$dpars$sigma
  # truncated mean of standard normal; see Wikipedia
  trunc_zmean <- (dnorm(zlb) - dnorm(zub)) / (pnorm(zub) - pnorm(zlb))  
  draws$dpars$mu + trunc_zmean * draws$dpars$sigma  
}

fitted_trunc_student <- function(draws, lb, ub) {
  zlb <- with(draws$dpars, (lb - mu) / sigma)
  zub <- with(draws$dpars, (ub - mu) / sigma)
  nu <- draws$dpars$nu
  # see Kim 2008: Moments of truncated Student-t distribution
  G1 <- gamma((nu - 1) / 2) * nu^(nu / 2) / 
    (2 * (pt(zub, df = nu) - pt(zlb, df = nu))
     * gamma(nu / 2) * gamma(0.5))
  A <- (nu + zlb^2) ^ (-(nu - 1) / 2)
  B <- (nu + zub^2) ^ (-(nu - 1) / 2)
  trunc_zmean <- G1 * (A - B)
  draws$dpars$mu + trunc_zmean * draws$dpars$sigma 
}

fitted_trunc_lognormal <- function(draws, lb, ub) {
  lb <- ifelse(lb < 0, 0, lb)
  m1 <- with(draws$dpars, 
    exp(mu + sigma^2 / 2) * 
      (pnorm((log(ub) - mu) / sigma - sigma) - 
       pnorm((log(lb) - mu) / sigma - sigma))
  )
  with(draws$dpars, 
    m1 / (plnorm(ub, meanlog = mu, sdlog = sigma) - 
          plnorm(lb, meanlog = mu, sdlog = sigma))
  )
}

fitted_trunc_gamma <- function(draws, lb, ub) {
  lb <- ifelse(lb < 0, 0, lb)
  draws$dpars$scale <- draws$dpars$mu / draws$dpars$shape
  # see Jawitz 2004: Moments of truncated continuous univariate distributions
  m1 <- with(draws$dpars, 
    scale / gamma(shape) * 
      (incgamma(1 + shape, ub / scale) - 
       incgamma(1 + shape, lb / scale))
  )
  with(draws$dpars, 
    m1 / (pgamma(ub, shape, scale = scale) - 
          pgamma(lb, shape, scale = scale))
  )
}

fitted_trunc_exponential <- function(draws, lb, ub) {
  lb <- ifelse(lb < 0, 0, lb)
  inv_mu <- 1 / draws$dpars$mu
  # see Jawitz 2004: Moments of truncated continuous univariate distributions
  m1 <- with(draws$dpars, mu * (incgamma(2, ub / mu) - incgamma(2, lb / mu)))
  m1 / (pexp(ub, rate = inv_mu) - pexp(lb, rate = inv_mu))
}

fitted_trunc_weibull <- function(draws, lb, ub) {
  lb <- ifelse(lb < 0, 0, lb)
  draws$dpars$a <- 1 + 1 / draws$dpars$shape
  draws$dpars$scale <- with(draws$dpars, mu / gamma(a))
  # see Jawitz 2004: Moments of truncated continuous univariate distributions
  m1 <- with(draws$dpars,
    scale * (incgamma(a, (ub / scale)^shape) - 
             incgamma(a, (lb / scale)^shape))
  )
  with(draws$dpars,
    m1 / (pweibull(ub, shape, scale = scale) - 
          pweibull(lb, shape, scale = scale))
  )
}

fitted_trunc_binomial <- function(draws, lb, ub) {
  lb <- ifelse(lb < -1, -1, lb)
  max_value <- max(draws$data$trials)
  ub <- ifelse(ub > max_value, max_value, ub)
  trials <- draws$data$trials
  if (length(trials) > 1) {
    trials <- as_draws_matrix(trials, dim_mu(draws))
  }
  args <- list(size = trials, prob = draws$dpars$mu)
  fitted_trunc_discrete(dist = "binom", args = args, lb = lb, ub = ub)
}

fitted_trunc_poisson <- function(draws, lb, ub) {
  lb <- ifelse(lb < -1, -1, lb)
  max_value <- 3 * max(draws$dpars$mu)
  ub <- ifelse(ub > max_value, max_value, ub)
  args <- list(lambda = draws$dpars$mu)
  fitted_trunc_discrete(dist = "pois", args = args, lb = lb, ub = ub)
}

fitted_trunc_negbinomial <- function(draws, lb, ub) {
  lb <- ifelse(lb < -1, -1, lb)
  max_value <- 3 * max(draws$dpars$mu)
  ub <- ifelse(ub > max_value, max_value, ub)
  args <- list(mu = draws$dpars$mu, size = draws$dpars$shape)
  fitted_trunc_discrete(dist = "nbinom", args = args, lb = lb, ub = ub)
}

fitted_trunc_geometric <- function(draws, lb, ub) {
  lb <- ifelse(lb < -1, -1, lb)
  max_value <- 3 * max(draws$dpars$mu)
  ub <- ifelse(ub > max_value, max_value, ub)
  args <- list(mu = draws$dpars$mu, size = 1)
  fitted_trunc_discrete(dist = "nbinom", args = args, lb = lb, ub = ub)
}

# fitted values for truncated discrete distributions
fitted_trunc_discrete <- function(dist, args, lb, ub) {
  stopifnot(is.matrix(lb), is.matrix(ub))
  message(
    "Computing fitted values for truncated ", 
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
  m1 <- do_call(cbind, m1)
  m1 / (do_call(cdf, c(list(ub), args)) - do_call(cdf, c(list(lb), args)))
}
