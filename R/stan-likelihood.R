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
  stopifnot(is.family(family))
  stopifnot(is.brmsterms(bterms))
  stopifnot(length(mix) == 1L)
  bterms$family <- family
  resp <- usc(combine_prefix(bterms))
  # prepare family part of the likelihood
  llh_args <- nlist(bterms, resp, mix)
  llh_fun <- paste0("stan_llh_", prepare_family(bterms)$fun)
  llh_pre <- do.call(llh_fun, llh_args)
  # incorporate other parts into the likelihood
  args <- nlist(llh_pre, bterms, data, resp, mix, ptheta)
  if (nzchar(mix)) {
    out <- do.call(stan_llh_mix, args)
  } else if (is.formula(bterms$adforms$cens)) {
    out <- do.call(stan_llh_cens, args)
  } else if (is.formula(bterms$adforms$weights)) {
    out <- do.call(stan_llh_weights, args)
  } else {
    out <- do.call(stan_llh_general, args)
  }
  if (grepl("\\[n\\]", out) && !nzchar(mix)) {
    # loop over likelihood if it cannot be vectorized
    out <- paste0("  for (n in 1:N) { \n    ", out, "    } \n")
  }
  out
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

stan_llh_general <- function(llh_pre, bterms, data, resp = "", ...) {
  # default likelihood in Stan language
  stopifnot(length(llh_pre) == 2L)
  n <- if (grepl("\\[n\\]", llh_pre[2])) "[n]"
  lpdf <- ifelse(use_int(bterms$family), "_lpmf", "_lpdf")
  tr <- stan_llh_trunc(llh_pre, bterms, data)
  paste0(
    tp(), llh_pre[1], lpdf, "(Y", resp, n, 
    " | ", llh_pre[2], ")", tr, "; \n"
  )
}

stan_llh_cens <- function(llh_pre, bterms, data, resp = "", ...) {
  # censored likelihood in Stan language
  stopifnot(length(llh_pre) == 2L)
  s <- collapse(rep(" ", 6))
  cens <- has_cens(bterms$adforms$cens, data = data)
  interval <- isTRUE(attr(cens, "interval"))
  lpdf <- ifelse(use_int(bterms$family), "_lpmf", "_lpdf")
  has_weights <- is.formula(bterms$adforms$weights)
  w <- if (has_weights) paste0("weights", resp, "[n] * ")
  tr <- stan_llh_trunc(llh_pre, bterms, data)
  tp <- tp()
  if (interval) {
    int_cens <- paste0(
      s, "} else if (cens", resp, "[n] == 2) { \n",
      s, tp, w, "log_diff_exp(", 
      llh_pre[1], "_lcdf(rcens", resp, "[n] | ", llh_pre[2], "), ",
      llh_pre[1], "_lcdf(Y", resp, "[n] | ", llh_pre[2], "))", tr, "; \n"
    )
  } else {
    int_cens <- ""
  }
  paste0(
    "  // special treatment of censored data \n",
    s, "if (cens", resp, "[n] == 0) {\n", 
    s, tp, w, llh_pre[1], lpdf, "(Y", resp, "[n] | ", llh_pre[2], ")", tr, ";\n",
    s, "} else if (cens", resp, "[n] == 1) {\n",         
    s, tp, w, llh_pre[1], "_lccdf(Y", resp, "[n] | ", llh_pre[2], ")", tr, ";\n",
    s, "} else if (cens", resp, "[n] == -1) {\n",
    s, tp, w, llh_pre[1], "_lcdf(Y", resp, "[n] | ", llh_pre[2], ")", tr, ";\n",
    int_cens, s, "} \n"
  )
}

stan_llh_weights <- function(llh_pre, bterms, data, resp = "", ...) {
  # weighted likelihood in Stan language
  stopifnot(length(llh_pre) == 2L)
  tr <- stan_llh_trunc(llh_pre, bterms, data)
  lpdf <- ifelse(use_int(bterms$family), "lpmf", "lpdf")
  paste0(
    tp(), "weights", resp, "[n] * ", llh_pre[1], "_", lpdf, 
    "(Y", resp, "[n] | ", llh_pre[2],")", tr, "; \n"
  )
}

stan_llh_mix <- function(llh_pre, bterms, data, mix, 
                         ptheta, resp = "", ...) {
  # likelihood of a single mixture component
  stopifnot(length(llh_pre) == 2L)
  theta <- ifelse(ptheta,
    paste0("theta", mix, resp, "[n]"), 
    paste0("log(theta", mix, resp, ")")
  )
  tr <- stan_llh_trunc(llh_pre, bterms, data)
  lpdf <- ifelse(use_int(bterms$family), "lpmf", "lpdf")
  paste0(
    "  ps[", mix, "] = ", theta, " + ",
    llh_pre[1], "_", lpdf, "(Y", resp, "[n] | ", llh_pre[2], ")", tr, "; \n"
  )
}

stan_llh_trunc <- function(llh_pre, bterms, data, short = FALSE) {
  # truncated part of the likelihood
  # Args:
  #   short: use the T[, ] syntax?
  bounds <- get_bounds(bterms$adforms$trunc, data = data)
  if (!any(bounds$lb > -Inf | bounds$ub < Inf)) {
    return("")
  }
  if (short) {
    # truncation using T[, ] syntax
    lb <- ifelse(any(bounds$lb > -Inf), "lb[n]", "")
    ub <- ifelse(any(bounds$ub < Inf), "ub[n]", "")
    out <- paste0(" T[", lb, ", ", ub, "]")
  } else {
    # truncation making use of _lcdf functions
    ms <- paste0(" - \n", collapse(rep(" ", 8)))
    if (any(bounds$lb > -Inf) && !any(bounds$ub < Inf)) {
      out <- paste0(ms, llh_pre[1], "_lccdf(lb[n] | ", llh_pre[2], ")")
    } else if (!any(bounds$lb > -Inf) && any(bounds$ub < Inf)) {
      out <- paste0(ms, llh_pre[1], "_lcdf(ub[n] | ", llh_pre[2], ")")
    } else if (any(bounds$lb > -Inf) && any(bounds$ub < Inf)) {
      trr <- paste0(llh_pre[1], "_lcdf(ub[n] | ", llh_pre[2], ")")
      trl <- paste0(llh_pre[1], "_lcdf(lb[n] | ", llh_pre[2], ")")
      out <- paste0(ms, "log_diff_exp(", trr, ", ", trl, ")")
    }
  }
  out
}

stan_llh_dpars <- function(bterms, reqn, resp = "", mix = "", dpars = NULL) {
  # prepare names of distributional parameters
  # Args:
  #   reqn: will the likelihood be wrapped in a loop over n?
  #   dpars: optinal names of distributional parameters to be prepared
  if (is.null(dpars)) {
    dpars <- valid_dpars(bterms) 
  }
  pred_dpars <- names(bterms$dpars)
  is_pred <- dpars %in% c("mu", dpar_class(pred_dpars))
  out <- paste0(dpars, mix, resp, ifelse(reqn & is_pred, "[n]", ""))
  named_list(dpars, out)
}

stan_llh_simple_lpdf <- function(lpdf, link, bterms) {
  # adjust lpdf name if a more efficient version is available
  # for a specific link. For instance poisson_log
  cens_or_trunc <- stan_llh_adj(bterms, c("cens", "trunc"))
  if (bterms$family$link == link && !cens_or_trunc) {
    lpdf <- paste0(lpdf, paste0("_", link))
  }
  lpdf
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

stan_llh_add_se <- function(sigma, bterms, reqn, resp = "") {
  # prepare the code for 'sigma' in the likelihood statement
  if (is.formula(bterms$adforms$se)) {
    nse <- if (reqn) "[n]"
    if (no_sigma(bterms)) {
      sigma <- paste0("se", resp, nse) 
    } else {
      sigma <- paste0("sqrt(square(", sigma, ") + se2", resp, nse, ")")
    }
  }
  sigma
}

stan_llh_adj <- function(x, adds = c("weights", "cens", "trunc")) {
  # checks if certain 'adds' are present so that the LL has to be adjusted
  # Args:
  #   x: named list of formulas or brmsterms object
  #   adds: vector of addition argument names
  stopifnot(all(adds %in% c("weights", "cens", "trunc")))
  if (is.brmsterms(x)) x <- x$adforms
  any(ulapply(x[adds], is.formula))
}

# one function per family
stan_llh_gaussian <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  p$sigma <- stan_llh_add_se(p$sigma, bterms, reqn, resp)
  c("normal", sargs(p$mu, p$sigma))
}

stan_llh_gaussian_mv <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix) || bterms$sigma_pred
  p <- list()
  p$Mu <- paste0("Mu", if (reqn) "[n]")
  p$LSigma <- paste0("LSigma", if (bterms$sigma_pred) "[n]")
  c("multi_normal_cholesky", sargs(p$Mu, p$LSigma))
}

stan_llh_gaussian_cov <- function(bterms, resp = "", mix = "") {
  if (stan_llh_adj(bterms)) {
    stop2("Invalid addition arguments for this model.")
  }
  p <- stan_llh_dpars(bterms, FALSE, resp, mix)
  l <- c("se2", "N_tg", "begin_tg", "end_tg", "nobs_tg", "res_cov_matrix")
  p[l] <- as.list(paste0(l, resp))
  c("normal_cov", sargs(
    p$mu, p$se2, p$N_tg, p$begin_tg,
    p$end_tg, p$nobs_tg, p$res_cov_matrix
  ))
}

stan_llh_gaussian_fixed <- function(bterms, resp = "", mix = "") {
  has_se <- is.formula(bterms$adforms$se)
  if (stan_llh_adj(bterms) || has_se) {
    stop2("Invalid addition arguments for this model.")
  }
  p <- stan_llh_dpars(bterms, FALSE, resp, mix)
  p$LV <- paste0("LV", resp)
  c("multi_normal_cholesky", sargs(p$mu, p$LV))
}

stan_llh_gaussian_lagsar <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, FALSE, resp, mix)
  p$sigma <- stan_llh_add_se(p$sigma, bterms, FALSE, resp)
  l <- c("lagsar", "W")
  p[l] <- as.list(paste0(l, resp))
  c("normal_lagsar", sargs(p$mu, p$sigma, p$lagsar, p$W))
}

stan_llh_gaussian_errorsar <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, FALSE, resp, mix)
  p$sigma <- stan_llh_add_se(p$sigma, bterms, FALSE, resp)
  l <- c("errorsar", "W")
  p[l] <- as.list(paste0(l, resp))
  c("normal_errorsar", sargs(p$mu, p$sigma, p$errorsar, p$W))
}

stan_llh_student <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  p$sigma <- stan_llh_add_se(p$sigma, bterms, reqn, resp)
  c("student_t", sargs(p$nu, p$mu, p$sigma))
}

stan_llh_student_mv <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix) || bterms$sigma_pred
  p <- stan_llh_dpars(bterms, reqn, resp, mix, dpars = "nu")
  p$Mu <- paste0("Mu", if (reqn) "[n]")
  p$Sigma <- paste0("Sigma", if (bterms$sigma_pred) "[n]")
  c("multi_student_t", sargs(p$nu, p$Mu, p$Sigma))
}

stan_llh_student_cov <- function(bterms, resp = "", mix = "") {
  if (stan_llh_adj(bterms)) {
    stop2("Invalid addition arguments for this model.")
  }
  p <- stan_llh_dpars(bterms, FALSE, resp, mix)
  l <- c("se2", "N_tg", "begin_tg", "end_tg", "nobs_tg", "res_cov_matrix")
  p[l] <- as.list(paste0(l, resp))
  c("student_t_cov", sargs(
    p$nu, p$mu, p$se2, p$N_tg, p$begin_tg,
    p$end_tg, p$nobs_tg, p$res_cov_matrix
  ))
}

stan_llh_student_fixed <- function(bterms, resp = "", mix = "") {
  has_se <- is.formula(bterms$adforms$se)
  if (stan_llh_adj(bterms) || has_se) {
    stop2("Invalid addition arguments for this model.")
  }
  p <- stan_llh_dpars(bterms, FALSE, resp, mix)
  p$V <- paste0("V", resp)
  c("multi_student_t", sargs(p$nu, p$mu, p$V))
}

stan_llh_student_lagsar <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, FALSE, resp, mix)
  p$sigma <- stan_llh_add_se(p$sigma, bterms, FALSE, resp)
  l <- c("lagsar", "W")
  p[l] <- as.list(paste0(l, resp))
  c("student_t_lagsar", sargs(p$nu, p$mu, p$sigma, p$lagsar, p$W))
}

stan_llh_student_errorsar <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, FALSE, resp, mix)
  p$sigma <- stan_llh_add_se(p$sigma, bterms, FALSE, resp)
  l <- c("errorsar", "W")
  p[l] <- as.list(paste0(l, resp))
  c("student_t_errorsar", sargs(p$nu, p$mu, p$sigma, p$errorsar, p$W))
}

stan_llh_lognormal <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  c("lognormal", sargs(p$mu, p$sigma))
}

stan_llh_asym_laplace <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  c("asym_laplace", sargs(p$mu, p$sigma, p$quantile))
}

stan_llh_skew_normal <- function(bterms, resp = "", mix = "", ...) {
  reqn <- stan_llh_adj(bterms) || nzchar(mix)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  p$sigma <- stan_llh_add_se(p$sigma, bterms, reqn, resp)
  # required because of CP parameterization of mu and sigma
  nomega <- any(grepl("\\[n\\]", c(p$sigma, p$alpha)))
  nomega <- if (reqn && nomega) "[n]"
  p$omega <- paste0("omega", mix, resp, nomega)
  c("skew_normal", sargs(p$mu, p$omega, p$alpha))
}

stan_llh_poisson <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  lpdf <- stan_llh_simple_lpdf("poisson", "log", bterms)
  c(lpdf, p$mu)
}

stan_llh_negbinomial <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  lpdf <- stan_llh_simple_lpdf("neg_binomial_2", "log", bterms)
  c(lpdf, sargs(p$mu, p$shape))
}

stan_llh_geometric <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  lpdf <- stan_llh_simple_lpdf("neg_binomial_2", "log", bterms)
  c(lpdf, sargs(p$mu, "1"))
}

stan_llh_binomial <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  p$trials <- paste0("trials", resp, if (reqn) "[n]")
  lpdf <- stan_llh_simple_lpdf("binomial", "logit", bterms)
  c(lpdf, sargs(p$trials, p$mu))
}

stan_llh_bernoulli <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  lpdf <- stan_llh_simple_lpdf("bernoulli", "logit", bterms)
  c(lpdf, p$mu)
}

stan_llh_gamma <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  c("gamma", sargs(p$shape, p$mu))
}

stan_llh_exponential <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  c("exponential", p$mu)
}

stan_llh_weibull <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  c("weibull", sargs(p$shape, p$mu))
}

stan_llh_frechet <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  c("frechet", sargs(p$nu, p$mu))
}

stan_llh_gen_extreme_value <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  c("gen_extreme_value", sargs(p$mu, p$sigma, p$xi))
}

stan_llh_exgaussian <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  c("exp_mod_normal", sargs(p$mu, p$sigma, paste0("inv(", p$beta, ")")))
}

stan_llh_inverse.gaussian <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  lpdf <- paste0("inv_gaussian", if (!reqn) "_vector")
  n <- if (reqn) "[n]"
  p$log_Y <- paste0("log_Y", resp, n)
  p$log_Y <- ifelse(reqn, p$log_Y, paste0("sum_", p$log_Y))
  p$sqrt_Y <- paste0("sqrt_Y", resp, n)
  c(lpdf, sargs(p$mu, p$shape, p$log_Y, p$sqrt_Y))
}

stan_llh_wiener <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  p$dec <- paste0("dec", resp, "[n]")
  c("wiener_diffusion", sargs(p$dec, p$bs, p$ndt, p$bias, p$mu))
}

stan_llh_beta <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  c("beta", sargs(
    paste0(p$mu, " * ", p$phi), 
    paste0("(1 - ", p$mu, ") * ", p$phi)
  ))
}

stan_llh_von_mises <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix) ||
    "kappa" %in% names(bterms$dpars)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  lpdf <- paste0("von_mises_", ifelse(reqn, "real", "vector"))
  c(lpdf, sargs(p$mu, p$kappa))
}

stan_llh_cumulative <- function(bterms, resp = "", mix = "") {
  simplify <- bterms$family$link == "logit" && 
    !"disc" %in% names(bterms$dpars) && !has_cs(bterms)
  if (simplify) {
    p <- stan_llh_dpars(bterms, TRUE, resp, mix)
    p$ord_intercept <- paste0("temp", resp, "_Intercept")
    out <- c("ordered_logistic", sargs(p$mu, p$ord_intercept))
  } else {
    out <- stan_llh_ordinal(bterms, resp, mix)
  }
  out
}

stan_llh_sratio <- function(bterms, resp = "", mix = "") {
  stan_llh_ordinal(bterms, resp, mix)
}

stan_llh_cratio <- function(bterms, resp = "", mix = "") {
  stan_llh_ordinal(bterms, resp, mix)
}

stan_llh_acat <- function(bterms, resp = "", mix = "") {
  stan_llh_ordinal(bterms, resp, mix)
}

stan_llh_categorical <- function(bterms, resp = "", mix = "") {
  stopifnot(bterms$family$link == "logit")
  p <- stan_llh_dpars(bterms, TRUE, resp, mix, dpars = "mu")
  c("categorical_logit", p$mu)
}

stan_llh_ordinal <- function(bterms, resp = "", mix = "") {
  # helper function for ordinal families
  has_cs <- has_cs(bterms)
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  p$ord_intercept <- paste0("temp", resp, "_Intercept")
  p$mucs <- if (has_cs) paste0("mucs", resp, "[n]")
  lpdf <- bterms$family$family
  lpdf <- paste0(lpdf, "_", bterms$family$link, if (has_cs) "_cs")
  c(lpdf, sargs(p$mu, p$mucs, p$ord_intercept, p$disc))
}

stan_llh_hurdle_poisson <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  usc_log <- if (bterms$family$link == "log") "_log"
  usc_logit <- stan_llh_dpar_usc_logit("hu", bterms)
  lpdf <- paste0("hurdle_poisson", usc_log, usc_logit)
  c(lpdf, sargs(p$mu, p$hu))
}

stan_llh_hurdle_negbinomial <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  usc_log <- if (bterms$family$link == "log") "_log"
  usc_logit <- stan_llh_dpar_usc_logit("hu", bterms)
  lpdf <- paste0("hurdle_neg_binomial", usc_log, usc_logit)
  c(lpdf, sargs(p$mu, p$shape, p$hu))
}

stan_llh_hurdle_gamma <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  usc_logit <- stan_llh_dpar_usc_logit("hu", bterms)
  lpdf <- paste0("hurdle_gamma", usc_logit)
  c(lpdf, sargs(p$shape, p$mu, p$hu))
}

stan_llh_hurdle_lognormal <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  usc_logit <- stan_llh_dpar_usc_logit("hu", bterms)
  lpdf <- paste0("hurdle_lognormal", usc_logit)
  c(lpdf, sargs(p$mu, p$sigma, p$hu))
}

stan_llh_zero_inflated_poisson <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  usc_log <- if (bterms$family$link == "log") "_log"
  usc_logit <- stan_llh_dpar_usc_logit("zi", bterms)
  lpdf <- paste0("zero_inflated_poisson", usc_log, usc_logit)
  c(lpdf, sargs(p$mu, p$zi))
}

stan_llh_zero_inflated_negbinomial <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  usc_log <- if (bterms$family$link == "log") "_log"
  usc_logit <- stan_llh_dpar_usc_logit("zi", bterms)
  lpdf <- paste0("zero_inflated_neg_binomial", usc_log, usc_logit)
  c(lpdf, sargs(p$mu, p$shape, p$zi))
}

stan_llh_zero_inflated_binomial <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  p$trials <- "trials[n]"
  usc_blogit <- if (bterms$family$link == "logit") "_blogit"
  usc_logit <- stan_llh_dpar_usc_logit("zi", bterms)
  lpdf <- paste0("zero_inflated_binomial", usc_blogit, usc_logit)
  c(lpdf, sargs(p$trials, p$mu, p$zi))
}

stan_llh_zero_inflated_beta <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  c("zero_inflated_beta", sargs(p$mu, p$phi, p$zi))
}

stan_llh_zero_one_inflated_beta <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  c("zero_one_inflated_beta", sargs(p$mu, p$phi, p$zoi, p$coi))
}

sargs <- function(...) {
  # prepare arguments of Stan likelihood statements
  paste0(c(...), collapse = ", ")
}

tp <- function(wsp = 2) {
  wsp <- collapse(rep(" ", wsp))
  paste0(wsp, "target += ")
}
