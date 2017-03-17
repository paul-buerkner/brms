#' @export
stan_llh.default <- function(family, bterms, data, autocor, mix = "", ...) {
  # Likelihood in Stan language
  # Args:
  #   family: the model family
  #   bterms: object of class brmsterms
  #   data: data passed by the user
  #   autocor: object of classe cor_brms
  stopifnot(is.family(family))
  stopifnot(is.brmsterms(bterms))
  stopifnot(length(mix) == 1L)
  link <- family$link
  family <- family$family
  mix <- as.character(mix)
  is_mix <- nzchar(mix)
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
          is_hurdle || is_zero_inflated || is_mix ||
          is_wiener(family) || is_exgaussian(family) || 
          is_asym_laplace(family) || is_gev(family) ||
          has_sigma && has_se && !use_cov(autocor) ||
          any(c("phi", "kappa") %in% auxpars)
  n <- ifelse(reqn, "[n]", "")
  # prepare auxiliary parameters
  p <- named_list(auxpars())
  p$mu <- paste0(ifelse(is_mv || is_categorical, "Mu", "mu"), mix, n)
  p$sigma <- stan_llh_sigma(family, bterms, mix)
  p$shape <- stan_llh_shape(family, bterms, mix)
  for (ap in setdiff(auxpars(), c("mu", "sigma", "shape"))) {
    p[[ap]] <- paste0(ap, mix, if (reqn && ap %in% auxpars) "[n]")
  }
  ord_args <- sargs(p$mu, if (has_cs) "mucs[n]", "temp_Intercept", p$disc)
  usc_logit <- ifelse(any(c("zi", "hu") %in% auxpars), "_logit", "")
  trials <- ifelse(llh_adj || is_zero_inflated, "trials[n]", "trials")
  simplify <- stan_has_built_in_fun(nlist(family, link)) &&
    !has_trunc && !has_cens && !"disc" %in% auxpars
  
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
        paste0("append_row(", sargs("zero", p$mu), ")")
      ),
      binomial = c(
        "binomial_logit", 
        sargs(trials, p$mu)
      ), 
      bernoulli = c(
        "bernoulli_logit", 
        p$mu)
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
        sargs(p$mu, "LSigma")
      ),
      gaussian_fixed = c(
        "multi_normal_cholesky", 
        sargs(p$mu, "LV")
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
        sargs(p$nu, p$mu, "Sigma")
      ),
      student_fixed = c(
        "multi_student_t", 
        sargs(p$nu, p$mu, "V")
      ),
      asym_laplace = c(
        "asym_laplace",
        sargs(p$mu, p$sigma, p$quantile)
      ),
      lognormal = c(
        "lognormal", 
        sargs(p$mu, p$sigma)
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
        "exgaussian", 
        sargs(p$mu, p$sigma, p$beta)
      ),
      inverse.gaussian = c(
        paste0("inv_gaussian", if (!reqn) "_vector"),
        sargs(p$mu, "shape", paste0(if (!reqn) "sum_", "log_Y", n),
              paste0("sqrt_Y", n))
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
        paste0("hurdle_poisson", usc_logit), 
        sargs(p$mu, p$hu)
      ),
      hurdle_negbinomial = c(
        paste0("hurdle_neg_binomial", usc_logit), 
        sargs(p$mu, p$hu, p$shape)
      ),
      hurdle_gamma = c(
        paste0("hurdle_gamma", usc_logit), 
        sargs(p$shape, p$mu, p$hu)
      ),
      hurdle_lognormal = c(
        paste0("hurdle_lognormal", usc_logit), 
        sargs(p$mu, p$hu, p$sigma)
      ),
      zero_inflated_poisson = c(
        paste0("zero_inflated_poisson", usc_logit), 
        sargs(p$mu, p$zi)
      ),
      zero_inflated_negbinomial = c(
        paste0("zero_inflated_neg_binomial", usc_logit),
        sargs(p$mu, p$zi, p$shape)
      ),
      zero_inflated_binomial = c(
        paste0("zero_inflated_binomial", usc_logit), 
        sargs(trials, p$mu, p$zi)
      ),
      zero_inflated_beta = c(
        paste0("zero_inflated_beta", usc_logit), 
        sargs(p$mu, p$zi, p$phi)
      )
    )
  }
  
  # write likelihood code
  interval <- isTRUE(attr(has_cens, "interval"))
  type <- match(TRUE, c(is_mix, has_cens, has_weights))
  type <- c("mix", "cens", "weights")[type]
  type <- ifelse(is.na(type), "general", type)
  llh <- switch(type, 
    mix = stan_llh_mix(llh_pre, family, mix, bounds),
    cens = stan_llh_cens(llh_pre, family, interval, has_weights, bounds),
    weights = stan_llh_weights(llh_pre, family, bounds),
    general = stan_llh_general(llh_pre, reqn, bounds)
  ) 
  if (reqn && !is_mix) {
    # loop over likelihood if it cannot be vectorized
    llh <- paste0("  for (n in 1:N) { \n    ", llh, "    } \n")
  }
  llh
}

#' @export
stan_llh.mixfamily <- function(family, bterms, ...) {
  ap_ids <- auxpar_id(names(bterms$auxpars))
  fap_ids <- auxpar_id(names(bterms$fauxpars))
  llh <- rep(NA, length(family$mix))
  for (i in seq_along(family$mix)) {
    sbterms <- bterms
    sbterms$auxpars <- sbterms$auxpars[ap_ids == i]
    sbterms$fauxpars <- sbterms$fauxpars[fap_ids == i]
    llh[i] <- stan_llh(family$mix[[i]], sbterms, mix = i, ...)
  }
  has_weights <- is.formula(bterms$adforms$weights)  
  lhs <- ifelse(has_weights, "lp_pre[n] = ", "target += ")
  paste0(
    "  for (n in 1:N) { \n",
    "      real ps[", length(llh), "]; \n",
    collapse("    ", llh),
    "      ", lhs, "log_sum_exp(ps); \n",
    "    } \n"
  )
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

stan_llh_mix <- function(llh_pre, family, mix, bounds = NULL) {
  stopifnot(length(llh_pre) == 2L)
  tr <- stan_llh_trunc(llh_pre, bounds = bounds, general = FALSE)
  lpdf <- ifelse(use_int(family), "lpmf", "lpdf")
  paste0(
    "  ps[", mix, "] = log(theta[", mix, "]) + ",
    llh_pre[1], "_", lpdf, "(Y[n] | ", llh_pre[2],")", tr, "; \n"
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

stan_llh_sigma <- function(family, bterms, mix = "") {
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
      out <- paste0("sqrt(sigma", mix, nsigma, "^2 + se2[n])")
    } else {
      out <- paste0(if (has_disp) "disp_", "sigma", mix, nsigma)
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

stan_llh_shape <- function(family, bterms, mix = "") {
  # prepare the code for 'shape' in the likelihood statement
  has_disp <- is.formula(bterms$adforms$disp)
  llh_adj <- stan_llh_adj(bterms$adforms)
  auxpars <- names(bterms$auxpars)
  nshape <- (llh_adj || is_forked(family)) &&
            (has_disp || "shape" %in% auxpars)
  nshape <- if (nshape) "[n]"
  paste0(if (has_disp) "disp_", "shape", mix, nshape)
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
