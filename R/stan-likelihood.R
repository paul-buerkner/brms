stan_llh <- function(family, bterms, data, autocor) {
  # Likelihood in Stan language
  # Args:
  #   family: the model family
  #   bterms: object of class brmsterms
  #   data: data passed by the user
  #   autocor: object of classe cor_brms
  stopifnot(is.family(family))
  stopifnot(is.brmsterms(bterms))
  link <- family$link
  family <- family$family
  is_categorical <- is_categorical(family)
  is_ordinal <- is_ordinal(family)
  is_hurdle <- is_hurdle(family)
  is_zero_inflated <- is_zero_inflated(family)
  is_forked <- is_forked(family)
  is_mv <- is_linear(family) && length(bterms$response) > 1L
  
  has_sigma <- has_sigma(family, bterms)
  has_se <- is.formula(bterms$adforms$se)
  has_weights <- is.formula(bterms$adforms$weights)
  has_cens <- has_cens(bterms$adforms$cens, data = data)
  has_disp <- is.formula(bterms$adforms$disp)
  has_cs <- has_cs(bterms)
  bounds <- get_bounds(bterms$adforms$trunc, data = data)
  has_trunc <- any(bounds$lb > -Inf) || any(bounds$ub < Inf)
  llh_adj <- stan_llh_adj(bterms$adforms)

  if (is_mv) {
    # prepare for use of a multivariate likelihood
    family <- paste0(family, "_mv")
  } else if (use_cov(autocor)) {
    # ARMA effects in covariance matrix formulation
    if (llh_adj) {
      stop2("Invalid addition arguments for this model.")
    }
    family <- paste0(family, "_cov")
  } else if (is.cor_fixed(autocor)) {
    if (has_se || llh_adj) {
      stop2("Invalid addition arguments for this model.")
    }
    family <- paste0(family, "_fixed")
  }
  
  auxpars <- names(bterms$auxpars)
  reqn <- llh_adj || is_categorical || is_ordinal || 
          is_hurdle || is_zero_inflated || 
          is_wiener(family) || is_exgaussian(family) || 
          is_asym_laplace(family) || is_gev(family) ||
          has_sigma && has_se && !use_cov(autocor) ||
          any(c("phi", "kappa") %in% auxpars)
  n <- ifelse(reqn, "[n]", "")
  # prepare auxiliary parameters
  p <- named_list(auxpars())
  p$sigma <- stan_llh_sigma(family, bterms)
  p$shape <- stan_llh_shape(family, bterms)
  for (ap in setdiff(auxpars(), c("sigma", "shape"))) {
    p[[ap]] <- paste0(ap, if (reqn && ap %in% auxpars) "[n]")
  }
  .logit <- ifelse(any(c("zi", "hu") %in% auxpars), "_logit", "")
  trials <- ifelse(llh_adj || is_zero_inflated, "trials[n]", "trials")

  simplify <- stan_has_built_in_fun(nlist(family, link)) &&
              !has_trunc && !has_cens && !"disc" %in% auxpars
  eta <- paste0(ifelse(is_mv, "Eta", "eta"), n)
  ord_args <- sargs("eta[n]", if (has_cs) "etacs[n]", 
                    "temp_Intercept", p$disc)
  
  if (simplify) { 
    llh_pre <- switch(family,
      poisson = c(
        "poisson_log", 
        eta
      ), 
      negbinomial = c(
        "neg_binomial_2_log", 
        sargs(eta, p$shape)
      ),
      geometric = c(
        "neg_binomial_2_log",
        sargs(eta, "1")
      ),
      cumulative = c(
        "ordered_logistic", 
        sargs("eta[n]", "temp_Intercept")
      ),
      categorical = c(
        "categorical_logit", 
        "append_row(zero, Eta[n])"
      ),
      binomial = c(
        "binomial_logit", 
        sargs(trials, eta)
      ), 
      bernoulli = c(
        "bernoulli_logit", 
        eta)
      )
  } else {
    llh_pre <- switch(family,
      gaussian = c(
        "normal", 
        sargs(eta, p$sigma)
      ),
      gaussian_cov = c(
        "normal_cov", 
        sargs(eta, "se2, N_tg, begin_tg, end_tg", 
              "nobs_tg, res_cov_matrix")
      ),
      gaussian_mv = c(
        "multi_normal_cholesky", 
        sargs(eta, "LSigma")
      ),
      gaussian_fixed = c(
        "multi_normal_cholesky", 
        sargs(eta, "LV")
      ),
      student = c(
        "student_t", 
        sargs(p$nu, eta, p$sigma)
      ),
      student_cov = c(
        "student_t_cov", 
        sargs(p$nu, eta, "se2, N_tg, begin_tg", 
              "end_tg, nobs_tg, res_cov_matrix")
      ),
      student_mv = c(
        "multi_student_t", 
        sargs(p$nu, eta, "Sigma")
      ),
      student_fixed = c(
        "multi_student_t", 
        sargs(p$nu, eta, "V")
      ),
      asym_laplace = c(
        "asym_laplace",
        sargs(eta, p$sigma, p$quantile)
      ),
      lognormal = c(
        "lognormal", 
        sargs(eta, p$sigma)
      ),
      poisson = c(
        "poisson", 
        eta
      ),
      negbinomial = c(
        "neg_binomial_2", 
        sargs(eta, p$shape)
      ),
      geometric = c(
        "neg_binomial_2", 
        sargs(eta, "1")
      ),
      binomial = c(
        "binomial", 
        sargs(trials, eta)
      ),
      bernoulli = c(
        "bernoulli", 
        eta
      ),
      gamma = c(
        "gamma", 
        sargs(p$shape, eta)
      ), 
      exponential = c(
        "exponential", 
        eta
      ),
      weibull = c(
        "weibull", 
        sargs(p$shape, eta)
      ), 
      frechet = c(
        "frechet", 
        sargs(p$nu, eta)
      ),
      gen_extreme_value = c(
        "gen_extreme_value",
        sargs(eta, p$sigma, p$xi)
      ),
      exgaussian = c(
        "exgaussian", 
        sargs(eta, p$sigma, p$beta)
      ),
      inverse.gaussian = c(
        paste0("inv_gaussian", if (!reqn) "_vector"),
        sargs(eta, "shape", paste0(if (!reqn) "sum_", "log_Y", n),
              paste0("sqrt_Y", n))
      ),
      wiener = c(
        "wiener_diffusion", 
        sargs("dec[n]", p$bs, p$ndt, p$bias, eta)
      ),
      beta = c(
        "beta", 
        sargs(paste0(eta, " * ", p$phi), 
              paste0("(1 - ", eta, ") * ", p$phi))
      ),
      von_mises = c(
        paste0("von_mises_", ifelse(reqn, "real", "vector")), 
        sargs(eta, p$kappa)
      ),
      cumulative = c(
        "cumulative", 
        ord_args
      ),
      sratio = c(
        "sratio",
        ord_args
      ),
      cratio = c(
        "cratio", 
        ord_args
      ),
      acat = c(
        "acat", 
        ord_args
      ),
      hurdle_poisson = c(
        paste0("hurdle_poisson", .logit), 
        sargs(eta, p$hu)
      ),
      hurdle_negbinomial = c(
        paste0("hurdle_neg_binomial", .logit), 
        sargs(eta, p$hu, p$shape)
      ),
      hurdle_gamma = c(
        paste0("hurdle_gamma", .logit), 
        sargs(p$shape, eta, p$hu)
      ),
      hurdle_lognormal = c(
        paste0("hurdle_lognormal", .logit), 
        sargs(eta, p$hu, p$sigma)
      ),
      zero_inflated_poisson = c(
        paste0("zero_inflated_poisson", .logit), 
        sargs(eta, p$zi)
      ),
      zero_inflated_negbinomial = c(
        paste0("zero_inflated_neg_binomial", .logit),
        sargs(eta, p$zi, p$shape)
      ),
      zero_inflated_binomial = c(
        paste0("zero_inflated_binomial", .logit), 
        sargs(trials, eta, p$zi)
      ),
      zero_inflated_beta = c(
        paste0("zero_inflated_beta", .logit), 
        sargs(eta, p$zi, p$phi)
      )
    )
  }
  
  # write likelihood code
  interval <- isTRUE(attr(has_cens, "interval"))
  type <- match(TRUE, c(has_cens, has_weights))
  type <- c("cens", "weights")[type]
  type <- ifelse(is.na(type), "general", type)
  llh <- switch(type, 
    cens = stan_llh_cens(llh_pre, family, interval, has_weights, bounds),
    weights = stan_llh_weights(llh_pre, family, bounds),
    general = stan_llh_general(llh_pre, reqn, bounds)
  ) 
  if (reqn) {
    # loop over likelihood if it cannot be vectorized
    llh <- paste0("  for (n in 1:N) { \n    ", llh, "    } \n")
  }
  llh
}

stan_llh_general <- function(llh_pre, reqn, bounds = NULL) {
  # default likelihood in Stan language
  # Args:
  #   reqn: does Y require the index 'n'?
  #   bounds: a list containing elements lb and ub
  stopifnot(length(llh_pre) == 2L)
  tr <- stan_llh_trunc(llh_pre, bounds = bounds)
  paste0(
    "  Y", ifelse(reqn, "[n]", ""), " ~ ", llh_pre[1], 
    "(", llh_pre[2], ")", tr, "; \n"
  )
}

stan_llh_cens <- function(llh_pre, family, interval, 
                          weights = FALSE, bounds = NULL) {
  # censored likelihood in Stan language
  # Args: 
  #   interval: are there interval censored responses present?
  #   weights: is the model additionally weighted?
  stopifnot(length(llh_pre) == 2L)
  s <- collapse(rep(" ", 6))
  tp <- "  target += "
  lpdf <- ifelse(use_int(family), "lpmf", "lpdf")
  w <- ifelse(weights, "weights[n] * ", "")
  tr <- stan_llh_trunc(llh_pre, bounds = bounds, general = FALSE)
  if (interval) {
    int_cens <- paste0(
      s, "} else if (cens[n] == 2) { \n",
      s, tp, w, "log_diff_exp(", 
      llh_pre[1], "_lcdf(rcens[n] | ", llh_pre[2], "), \n",
      collapse(rep(" ", 31)),
      llh_pre[1], "_lcdf(Y[n] | ", llh_pre[2], "))", tr, "; \n"
    )
  } else {
    int_cens <- ""
  }
  paste0(
    "  // special treatment of censored data \n",
    s, "if (cens[n] == 0) {\n", 
    s, tp, w, llh_pre[1], "_", lpdf, "(Y[n] | ", llh_pre[2], ")", tr, ";\n",
    s, "} else if (cens[n] == 1) {\n",         
    s, tp, w, llh_pre[1], "_lccdf(Y[n] | ", llh_pre[2], ")", tr, ";\n",
    s, "} else if (cens[n] == -1) {\n",
    s, tp, w, llh_pre[1], "_lcdf(Y[n] | ", llh_pre[2], ")", tr, ";\n",
    int_cens, s, "} \n"
  )
}

stan_llh_weights <- function(llh_pre, family, bounds = NULL) {
  # weighted likelihood in Stan language
  stopifnot(length(llh_pre) == 2L)
  tr <- stan_llh_trunc(llh_pre, bounds = bounds, general = FALSE)
  lpdf <- ifelse(use_int(family), "lpmf", "lpdf")
  paste0(
    "  lp_pre[n] = ", llh_pre[1], "_", lpdf, 
    "(Y[n] | ", llh_pre[2],")", tr, "; \n"
  )
}

stan_llh_trunc <- function(llh_pre, bounds, general = TRUE) {
  # truncated part of the likelihood
  # Args:
  #   general: use the T[, ] syntax?
  if (general) {
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
    ms <- paste0(" - \n", collapse(rep(" ", 18)))
    if (any(bounds$lb > -Inf) && !any(bounds$ub < Inf)) {
      tr <- paste0(ms, llh_pre[1], "_lccdf(lb[n] | ", llh_pre[2], ")")
    } else if (!any(bounds$lb > -Inf) && any(bounds$ub < Inf)) {
      tr <- paste0(ms, llh_pre[1], "_lcdf(ub[n] | ", llh_pre[2], ")")
    } else if (any(bounds$lb > -Inf) && any(bounds$ub < Inf)) {
      trr <- paste0(llh_pre[1], "_lcdf(ub[n] | ", llh_pre[2], ")")
      trl <- paste0(llh_pre[1], "_lcdf(lb[n] | ", llh_pre[2], ")")
      tr <- paste0(
        ms, "log_diff_exp(", trr, ", \n",
        collapse(rep(" ", 31)), trl, ")"
      )
    } else {
      tr <- ""
    }
  }
  tr
}

stan_llh_sigma <- function(family, bterms) {
  # prepare the code for 'sigma' in the likelihood statement
  has_sigma <- has_sigma(family, bterms)
  has_se <- is.formula(bterms$adforms$se)
  has_disp <- is.formula(bterms$adforms$disp)
  llh_adj <- stan_llh_adj(bterms$adforms)
  auxpars <- names(bterms$auxpars)
  nsigma <- llh_adj || has_se || is_exgaussian(family) || is_gev(family)
  nsigma <- nsigma && (has_disp || "sigma" %in% auxpars)
  nsigma <- if (nsigma) "[n]"
  nse <- if (llh_adj) "[n]"
  if (has_sigma) {
    if (has_se) {
      out <- paste0("sqrt(sigma", nsigma, "^2 + se2[n])")
    } else {
      out <- paste0(if (has_disp) "disp_", "sigma", nsigma) 
    }
  } else {
    if (has_se) {
      out <- paste0("se", nse) 
    } else {
      out <- NULL
    }
  }
  out
}

stan_llh_shape <- function(family, bterms) {
  # prepare the code for 'shape' in the likelihood statement
  has_disp <- is.formula(bterms$adforms$disp)
  llh_adj <- stan_llh_adj(bterms$adforms)
  auxpars <- names(bterms$auxpars)
  nshape <- (llh_adj || is_forked(family)) &&
            (has_disp || "shape" %in% auxpars)
  nshape <- if (nshape) "[n]"
  paste0(if (has_disp) "disp_", "shape", nshape)
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
