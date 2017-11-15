stan_llh <- function(family, ...) {
  # Stan code for the model likelihood 
  UseMethod("stan_llh")
}

#' @export
stan_llh.default <- function(family, bterms, data, mix = "", 
                             ptheta = FALSE, ...) {
  # Likelihood in Stan language
  # Args:
  #   family: the model family
  #   bterms: object of class brmsterms
  #   data: data passed by the user
  #   autocor: object of classe cor_brms
  #   mix: optional mixture component ID
  #   ptheta: are mixing proportions predicted?
  stopifnot(is.family(family), is.brmsterms(bterms), length(mix) == 1L)
  link <- family$link
  family <- family$family
  autocor <- bterms$autocor
  mix <- as.character(mix)
  is_mix <- nzchar(mix)
  is_categorical <- is_categorical(family)
  is_ordinal <- is_ordinal(family)
  is_hurdle <- is_hurdle(family)
  is_zero_inflated <- is_zero_inflated(family)
  is_mv <- grepl("_mv$", family)
  has_sigma <- has_sigma(family, bterms)
  has_se <- is.formula(bterms$adforms$se)
  has_weights <- is.formula(bterms$adforms$weights)
  has_cens <- has_cens(bterms$adforms$cens, data = data)
  has_cs <- has_cs(bterms)
  bounds <- get_bounds(bterms$adforms$trunc, data = data)
  has_trunc <- any(bounds$lb > -Inf) || any(bounds$ub < Inf)
  llh_adj <- stan_llh_adj(bterms$adforms)
  
  dpars <- names(bterms$dpars)
  reqn <- llh_adj || !stan_is_vectorized(family) || 
    is_mix || is_mv && isTRUE(bterms$sigma_pred) ||
    has_sigma && has_se && !use_cov(autocor) ||
    any(c("phi", "kappa") %in% dpars)
  n <- ifelse(reqn, "[n]", "")
  # prepare distributional parameters
  resp <- usc(combine_prefix(bterms))
  p <- named_list(dpars())
  p$sigma <- stan_llh_sigma(family, bterms, resp = resp, mix = mix)
  for (ap in setdiff(dpars(), "sigma")) {
    is_ap_pred <- ap %in% dpar_class(dpars)
    p[[ap]] <- paste0(ap, mix, resp, if (reqn && is_ap_pred) "[n]")
  }
  if (family == "skew_normal") {
    # required because of CP parameterization of mu and sigma
    nomega <- if (reqn && any(c("sigma", "alpha") %in% dpars)) "[n]"
    p$omega <- paste0("omega", mix, resp, nomega)
  }
  ord_args <- sargs(
    p$mu, if (has_cs) paste0("mucs[n]", resp), 
    paste0("temp_Intercept", resp), p$disc
  )
  ord_family <- paste0(family, "_", link, if (has_cs) "_cs")
  usc_logit <- stan_llh_dpar_usc_logit(c("zi", "hu"), bterms)
  trials <- ifelse(llh_adj || is_zero_inflated, "trials[n]", "trials")
  nSigma <- ifelse(is_mv && isTRUE(bterms$sigma_pred), "[n]", "")
  
  if (use_cov(autocor)) {
    # ARMA effects in covariance matrix formulation
    if (llh_adj) {
      stop2("Invalid addition arguments for this model.")
    }
    str_add(family) <- "_cov"
  } else if (is.cor_sar(autocor)) {
    if (has_se || llh_adj) {
      stop2("Invalid addition arguments for this model.")
    }
    if (identical(autocor$type, "lag")) {
      str_add(family) <- "_lagsar"
    } else if (identical(autocor$type, "error")) {
      str_add(family) <- "_errorsar"
    }
  } else if (is.cor_fixed(autocor)) {
    if (has_se || llh_adj) {
      stop2("Invalid addition arguments for this model.")
    }
    str_add(family) <- "_fixed"
  }
  
  simplify <- stan_has_built_in_fun(nlist(family, link)) &&
    !has_trunc && !has_cens && !"disc" %in% dpars
  if (simplify) { 
    llh_pre <- switch(family,
      poisson = c(
        "poisson_log", 
        p$mu
      ), 
      negbinomial = c(
        "neg_binomial_2_log", 
        sargs(p$mu, p$shape)
      ),
      geometric = c(
        "neg_binomial_2_log",
        sargs(p$mu, "1")
      ),
      cumulative = c(
        "ordered_logistic", 
        sargs(p$mu, "temp_Intercept")
      ),
      categorical = c(
        "categorical_logit", 
        paste0("append_row(", sargs("zero", "Mu[n]"), ")")
      ),
      binomial = c(
        "binomial_logit", 
        sargs(trials, p$mu)
      ), 
      bernoulli = c(
        "bernoulli_logit", 
        p$mu
      ),
      hurdle_poisson = c(
        paste0("hurdle_poisson_log", usc_logit), 
        sargs(p$mu, p$hu)
      ),
      hurdle_negbinomial = c(
        paste0("hurdle_neg_binomial_log", usc_logit), 
        sargs(p$mu, p$shape, p$hu)
      ),
      zero_inflated_poisson = c(
        paste0("zero_inflated_poisson_log", usc_logit), 
        sargs(p$mu, p$zi)
      ),
      zero_inflated_negbinomial = c(
        paste0("zero_inflated_neg_binomial_log", usc_logit),
        sargs(p$mu, p$shape, p$zi)
      ),
      zero_inflated_binomial = c(
        paste0("zero_inflated_binomial_blogit", usc_logit), 
        sargs(trials, p$mu, p$zi)
      )
    )
  } else {
    llh_pre <- switch(family,
      gaussian = c(
        "normal", 
        sargs(p$mu, p$sigma)
      ),
      gaussian_cov = c(
        "normal_cov", 
        sargs(p$mu, "se2, N_tg, begin_tg, end_tg", 
              "nobs_tg, res_cov_matrix")
      ),
      gaussian_mv = c(
        "multi_normal_cholesky", 
        sargs(paste0("Mu", n), paste0("LSigma", nSigma))
      ),
      gaussian_fixed = c(
        "multi_normal_cholesky", 
        sargs(p$mu, "LV")
      ),
      gaussian_lagsar = c(
        "normal_lagsar",
        sargs(p$mu, p$sigma, "lagsar", "W")
      ),
      gaussian_errorsar = c(
        "normal_errorsar",
        sargs(p$mu, p$sigma, "errorsar", "W")
      ),
      student = c(
        "student_t", 
        sargs(p$nu, p$mu, p$sigma)
      ),
      student_cov = c(
        "student_t_cov", 
        sargs(p$nu, p$mu, "se2, N_tg, begin_tg", 
              "end_tg, nobs_tg, res_cov_matrix")
      ),
      student_mv = c(
        "multi_student_t", 
        sargs(p$nu, paste0("Mu", n), paste0("LSigma", nSigma))
      ),
      student_fixed = c(
        "multi_student_t", 
        sargs(p$nu, p$mu, "V")
      ),
      student_lagsar = c(
        "student_t_lagsar",
        sargs(p$nu, p$mu, p$sigma, "lagsar", "W")
      ),
      student_errorsar = c(
        "student_t_errorsar",
        sargs(p$nu, p$mu, p$sigma, "errorsar", "W")
      ),
      asym_laplace = c(
        "asym_laplace",
        sargs(p$mu, p$sigma, p$quantile)
      ),
      lognormal = c(
        "lognormal", 
        sargs(p$mu, p$sigma)
      ),
      skew_normal = c(
        "skew_normal",
        sargs(p$mu, p$omega, p$alpha)
      ),
      poisson = c(
        "poisson", 
        p$mu
      ),
      negbinomial = c(
        "neg_binomial_2", 
        sargs(p$mu, p$shape)
      ),
      geometric = c(
        "neg_binomial_2", 
        sargs(p$mu, "1")
      ),
      binomial = c(
        "binomial", 
        sargs(trials, p$mu)
      ),
      bernoulli = c(
        "bernoulli", 
        p$mu
      ),
      gamma = c(
        "gamma", 
        sargs(p$shape, p$mu)
      ), 
      exponential = c(
        "exponential", 
        p$mu
      ),
      weibull = c(
        "weibull", 
        sargs(p$shape, p$mu)
      ), 
      frechet = c(
        "frechet", 
        sargs(p$nu, p$mu)
      ),
      gen_extreme_value = c(
        "gen_extreme_value",
        sargs(p$mu, p$sigma, p$xi)
      ),
      exgaussian = c(
        "exp_mod_normal", 
        sargs(p$mu, p$sigma, paste0("inv(", p$beta, ")"))
      ),
      inverse.gaussian = c(
        paste0("inv_gaussian", if (!reqn) "_vector"),
        sargs(p$mu, p$shape, 
          paste0(if (!reqn) "sum_", "log_Y", resp, n),
          paste0("sqrt_Y", resp, n)
        )
      ),
      wiener = c(
        "wiener_diffusion", 
        sargs("dec[n]", p$bs, p$ndt, p$bias, p$mu)
      ),
      beta = c(
        "beta", 
        sargs(paste0(p$mu, " * ", p$phi), 
              paste0("(1 - ", p$mu, ") * ", p$phi))
      ),
      von_mises = c(
        paste0("von_mises_", ifelse(reqn, "real", "vector")), 
        sargs(p$mu, p$kappa)
      ),
      cumulative = c(
        ord_family,
        ord_args
      ),
      sratio = c(
        ord_family,
        ord_args
      ),
      cratio = c(
        ord_family,
        ord_args
      ),
      acat = c(
        ord_family,
        ord_args
      ),
      hurdle_poisson = c(
        paste0("hurdle_poisson", usc_logit), 
        sargs(p$mu, p$hu)
      ),
      hurdle_negbinomial = c(
        paste0("hurdle_neg_binomial", usc_logit), 
        sargs(p$mu, p$shape, p$hu)
      ),
      hurdle_gamma = c(
        paste0("hurdle_gamma", usc_logit), 
        sargs(p$shape, p$mu, p$hu)
      ),
      hurdle_lognormal = c(
        paste0("hurdle_lognormal", usc_logit), 
        sargs(p$mu, p$sigma, p$hu)
      ),
      zero_inflated_poisson = c(
        paste0("zero_inflated_poisson", usc_logit), 
        sargs(p$mu, p$zi)
      ),
      zero_inflated_negbinomial = c(
        paste0("zero_inflated_neg_binomial", usc_logit),
        sargs(p$mu, p$shape, p$zi)
      ),
      zero_inflated_binomial = c(
        paste0("zero_inflated_binomial", usc_logit), 
        sargs(trials, p$mu, p$zi)
      ),
      zero_inflated_beta = c(
        paste0("zero_inflated_beta", usc_logit), 
        sargs(p$mu, p$phi, p$zi)
      ),
      zero_one_inflated_beta = c(
        "zero_one_inflated_beta", 
        sargs(p$mu, p$phi, p$zoi, p$coi)
      )
    )
  }
  
  # write likelihood code
  interval <- isTRUE(attr(has_cens, "interval"))
  type <- match(TRUE, c(is_mix, has_cens, has_weights))
  type <- c("mix", "cens", "weights")[type]
  type <- ifelse(is.na(type), "general", type)
  args <- nlist(
    llh_pre, family, reqn, resp, bounds, mix, 
    ptheta, interval, weights = has_weights
  )
  llh <- switch(type, 
    mix = do.call(stan_llh_mix, args),
    cens = do.call(stan_llh_cens, args),
    weights = do.call(stan_llh_weights, args),
    general = do.call(stan_llh_general, args)
  ) 
  if (reqn && !is_mix) {
    # loop over likelihood if it cannot be vectorized
    llh <- paste0("  for (n in 1:N) { \n    ", llh, "    } \n")
  }
  llh
}

#' @export
stan_llh.mixfamily <- function(family, bterms, ...) {
  ap_ids <- dpar_id(names(bterms$dpars))
  fap_ids <- dpar_id(names(bterms$fdpars))
  ptheta <- any(dpar_class(names(bterms$dpars)) %in% "theta")
  llh <- rep(NA, length(family$mix))
  for (i in seq_along(family$mix)) {
    sbterms <- bterms
    sbterms$dpars <- sbterms$dpars[ap_ids == i]
    sbterms$fdpars <- sbterms$fdpars[fap_ids == i]
    llh[i] <- stan_llh(
      family$mix[[i]], sbterms, mix = i, ptheta = ptheta, ...
    )
  }
  resp <- usc(combine_prefix(bterms))
  has_weights <- is.formula(bterms$adforms$weights)  
  weights <- if (has_weights) paste0("weights", resp, "[n] * ")
  paste0(
    "  for (n in 1:N) { \n",
    "      real ps[", length(llh), "]; \n",
    collapse("    ", llh),
    "    ", tp(), weights, "log_sum_exp(ps); \n",
    "    } \n"
  )
}

#' @export
stan_llh.brmsterms <- function(family, ...) {
  paste0("  ", stan_llh(family$family, bterms = family, ...))
}

#' @export
stan_llh.mvbrmsterms <- function(family, ...) {
  if (family$rescor) {
    out <- stan_llh(as.brmsterms(family), ...)
  } else {
    out <- collapse(ulapply(family$terms, stan_llh, ...)) 
  }
  out
}

stan_llh_general <- function(llh_pre, family, reqn, resp = "", 
                             bounds = NULL, ...) {
  # default likelihood in Stan language
  # Args:
  #   reqn: does Y require the index 'n'?
  #   bounds: a list containing elements lb and ub
  stopifnot(length(llh_pre) == 2L)
  lpdf <- ifelse(use_int(family), "_lpmf", "_lpdf")
  tr <- stan_llh_trunc(llh_pre, bounds = bounds)
  paste0(
    tp(), llh_pre[1], lpdf, "(Y", resp, ifelse(reqn, "[n]", ""), 
    " | ", llh_pre[2], ")", tr, "; \n"
  )
}

stan_llh_cens <- function(llh_pre, family, interval, resp = "",
                          weights = FALSE, bounds = NULL, ...) {
  # censored likelihood in Stan language
  # Args: 
  #   interval: are there interval censored responses present?
  #   weights: is the model additionally weighted?
  stopifnot(length(llh_pre) == 2L)
  s <- collapse(rep(" ", 6))
  lpdf <- ifelse(use_int(family), "_lpmf", "_lpdf")
  w <- if (weights) paste0("weights", resp, "[n] * ")
  tr <- stan_llh_trunc(llh_pre, bounds = bounds)
  tp <- tp()
  if (interval) {
    int_cens <- paste0(
      s, "} else if (cens[n] == 2) { \n",
      s, tp, w, "log_diff_exp(", 
      llh_pre[1], "_lcdf(rcens[n] | ", llh_pre[2], "), ",
      llh_pre[1], "_lcdf(Y[n] | ", llh_pre[2], "))", tr, "; \n"
    )
  } else {
    int_cens <- ""
  }
  paste0(
    "  // special treatment of censored data \n",
    s, "if (cens[n] == 0) {\n", 
    s, tp, w, llh_pre[1], lpdf, "(Y", resp, "[n] | ", llh_pre[2], ")", tr, ";\n",
    s, "} else if (cens[n] == 1) {\n",         
    s, tp, w, llh_pre[1], "_lccdf(Y", resp, "[n] | ", llh_pre[2], ")", tr, ";\n",
    s, "} else if (cens[n] == -1) {\n",
    s, tp, w, llh_pre[1], "_lcdf(Y", resp, "[n] | ", llh_pre[2], ")", tr, ";\n",
    int_cens, s, "} \n"
  )
}

stan_llh_weights <- function(llh_pre, family, resp = "", 
                             bounds = NULL, ...) {
  # weighted likelihood in Stan language
  stopifnot(length(llh_pre) == 2L)
  tr <- stan_llh_trunc(llh_pre, bounds = bounds)
  lpdf <- ifelse(use_int(family), "lpmf", "lpdf")
  paste0(
    tp(), "weights", resp, "[n] * ", llh_pre[1], "_", lpdf, 
    "(Y", resp, "[n] | ", llh_pre[2],")", tr, "; \n"
  )
}

stan_llh_mix <- function(llh_pre, family, mix, ptheta, 
                         resp = "", bounds = NULL, ...) {
  # likelihood of a single mixture component
  stopifnot(length(llh_pre) == 2L)
  theta <- ifelse(ptheta,
    paste0("theta", mix, resp, "[n]"), 
    paste0("log(theta", mix, resp, ")")
  )
  tr <- stan_llh_trunc(llh_pre, bounds = bounds)
  lpdf <- ifelse(use_int(family), "lpmf", "lpdf")
  paste0(
    "  ps[", mix, "] = ", theta, " + ",
    llh_pre[1], "_", lpdf, "(Y", resp, "[n] | ", llh_pre[2], ")", tr, "; \n"
  )
}

stan_llh_trunc <- function(llh_pre, bounds, short = FALSE) {
  # truncated part of the likelihood
  # Args:
  #   short: use the T[, ] syntax?
  if (short) {
    if (any(bounds$lb > -Inf) || any(bounds$ub < Inf)) {
      # truncation using T[, ] syntax
      lb <- ifelse(any(bounds$lb > -Inf), "lb[n]", "")
      ub <- ifelse(any(bounds$ub < Inf), "ub[n]", "")
      tr <- paste0(" T[", lb, ", ", ub, "]")
    } else {
      tr <- ""
    }
  } else {
    # truncation making use of _lcdf functions
    ms <- paste0(" - \n", collapse(rep(" ", 8)))
    if (any(bounds$lb > -Inf) && !any(bounds$ub < Inf)) {
      tr <- paste0(ms, llh_pre[1], "_lccdf(lb[n] | ", llh_pre[2], ")")
    } else if (!any(bounds$lb > -Inf) && any(bounds$ub < Inf)) {
      tr <- paste0(ms, llh_pre[1], "_lcdf(ub[n] | ", llh_pre[2], ")")
    } else if (any(bounds$lb > -Inf) && any(bounds$ub < Inf)) {
      trr <- paste0(llh_pre[1], "_lcdf(ub[n] | ", llh_pre[2], ")")
      trl <- paste0(llh_pre[1], "_lcdf(lb[n] | ", llh_pre[2], ")")
      tr <- paste0(ms, "log_diff_exp(", trr, ", ", trl, ")")
    } else {
      tr <- ""
    }
  }
  tr
}

stan_llh_sigma <- function(family, bterms, resp = "", mix = "") {
  # prepare the code for 'sigma' in the likelihood statement
  has_sigma <- has_sigma(family, bterms)
  has_se <- is.formula(bterms$adforms$se)
  llh_adj <- stan_llh_adj(bterms$adforms)
  nsigma <- llh_adj || has_se || nzchar(mix) ||
    family %in% c("exgaussian", "gen_extreme_value", "asym_laplace")
  nsigma <- nsigma && paste0("sigma", mix) %in% names(bterms$dpars)
  nsigma <- if (nsigma) "[n]"
  nse <- if (llh_adj) "[n]"
  if (has_sigma) {
    if (has_se) {
      out <- paste0("sqrt(sigma", mix, resp, nsigma, "^2 + se2", resp, "[n])")
    } else {
      out <- paste0("sigma", mix, resp, nsigma)
    }
  } else {
    if (has_se) {
      out <- paste0("se", resp, nse) 
    } else {
      out <- ""
    }
  }
  out
}

stan_llh_dpar_usc_logit <- function(dpars, bterms) {
  # prepare _logit suffix for distributional parameters
  # currently only used in zero-inflated and hurdle models
  stopifnot(is.brmsterms(bterms))
  use_logit <- any(ulapply(dpars, function(dp) 
    isTRUE(bterms$dpars[[dp]]$family$link == "logit")
  ))
  ifelse(use_logit, "_logit", "")
}

stan_llh_adj <- function(adforms, adds = c("weights", "cens", "trunc")) {
  # checks if certain 'adds' are present so that the LL has to be adjusted
  # Args:
  #   adforms: named list of formulas
  #   adds: vector of addition argument names
  stopifnot(all(adds %in% c("weights", "cens", "trunc")))
  any(ulapply(adforms[adds], is.formula))
}

sargs <- function(...) {
  # prepare arguments of Stan likelihood statements
  paste0(c(...), collapse = ", ")
}

tp <- function(wsp = 2) {
  wsp <- collapse(rep(" ", wsp))
  paste0(wsp, "target += ")
}
