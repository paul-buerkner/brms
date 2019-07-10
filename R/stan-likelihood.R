# unless otherwise specified, functions return a single character 
# string defining the likelihood of the model in Stan language

stan_llh <- function(family, ...) {
  UseMethod("stan_llh")
}

# Stan code for the model likelihood 
# @param family the model family
# @param bterms object of class brmsterms
# @param data data passed by the user
# @param mix optional mixture component ID
# @param ptheta are mixing proportions predicted?
#' @export
stan_llh.default <- function(family, bterms, data, mix = "", 
                             ptheta = FALSE, ...) {
  stopifnot(is.family(family))
  stopifnot(is.brmsterms(bterms))
  stopifnot(length(mix) == 1L)
  bterms$family <- family
  resp <- usc(combine_prefix(bterms))
  # prepare family part of the likelihood
  llh_args <- nlist(bterms, resp, mix)
  llh_fun <- paste0("stan_llh_", prepare_family(bterms)$fun)
  llh <- do_call(llh_fun, llh_args)
  # incorporate other parts into the likelihood
  args <- nlist(llh, bterms, data, resp, mix, ptheta)
  if (nzchar(mix)) {
    out <- do_call(stan_llh_mix, args)
  } else if (is.formula(bterms$adforms$cens)) {
    out <- do_call(stan_llh_cens, args)
  } else if (is.formula(bterms$adforms$weights)) {
    out <- do_call(stan_llh_weights, args)
  } else {
    out <- do_call(stan_llh_general, args)
  }
  if (grepl("\\[n\\]", out) && !nzchar(mix)) {
    # loop over likelihood if it cannot be vectorized
    out <- paste0("  for (n in 1:N", resp, ") {\n    ", out, "    }\n")
  }
  out
}

#' @export
stan_llh.mixfamily <- function(family, bterms, ...) {
  dp_ids <- dpar_id(names(bterms$dpars))
  fdp_ids <- dpar_id(names(bterms$fdpars))
  resp <- usc(bterms$resp)
  ptheta <- any(dpar_class(names(bterms$dpars)) %in% "theta")
  llh <- rep(NA, length(family$mix))
  for (i in seq_along(family$mix)) {
    sbterms <- bterms
    sbterms$dpars <- sbterms$dpars[dp_ids == i]
    sbterms$fdpars <- sbterms$fdpars[fdp_ids == i]
    llh[i] <- stan_llh(
      family$mix[[i]], sbterms, mix = i, ptheta = ptheta, ...
    )
  }
  resp <- usc(combine_prefix(bterms))
  has_weights <- is.formula(bterms$adforms$weights)  
  weights <- str_if(has_weights, glue("weights{resp}[n] * "))
  out <- glue(
    "  // likelihood of the mixture model\n",
    "  for (n in 1:N{resp}) {{\n",
    "      real ps[{length(llh)}];\n"
  )
  str_add(out) <- collapse("    ", llh)
  str_add(out) <- glue(
    "    {tp()}{weights}log_sum_exp(ps);\n",
    "  }}\n"
  )
  out
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

# default likelihood in Stan language
stan_llh_general <- function(llh, bterms, data, resp = "", ...) {
  stopifnot(is.sdist(llh))
  n <- str_if(grepl("\\[n\\]", llh$args), "[n]")
  lpdf <- stan_llh_lpdf_name(bterms)
  Y <- stan_llh_Y_name(bterms)
  tr <- stan_llh_trunc(llh, bterms, data, resp = resp)
  glue("{tp()}{llh$dist}_{lpdf}({Y}{resp}{n}{llh$shift} | {llh$args}){tr};\n")
}

# censored likelihood in Stan language
stan_llh_cens <- function(llh, bterms, data, resp = "", ...) {
  stopifnot(is.sdist(llh))
  s <- wsp(nsp = 6)
  cens <- has_cens(bterms, data = data)
  lpdf <- stan_llh_lpdf_name(bterms)
  has_weights <- is.formula(bterms$adforms$weights)
  Y <- stan_llh_Y_name(bterms)
  w <- str_if(has_weights, glue("weights{resp}[n] * "))
  tr <- stan_llh_trunc(llh, bterms, data, resp = resp)
  tp <- tp()
  out <- glue(
    "  // special treatment of censored data\n",
    s, "if (cens{resp}[n] == 0) {{\n", 
    s, "{tp}{w}{llh$dist}_{lpdf}({Y}{resp}[n]{llh$shift} | {llh$args}){tr};\n",
    s, "}} else if (cens{resp}[n] == 1) {{\n",         
    s, "{tp}{w}{llh$dist}_lccdf({Y}{resp}[n]{llh$shift} | {llh$args}){tr};\n",
    s, "}} else if (cens{resp}[n] == -1) {{\n",
    s, "{tp}{w}{llh$dist}_lcdf({Y}{resp}[n]{llh$shift} | {llh$args}){tr};\n"
  )
  if (isTRUE(attr(cens, "interval"))) {
    str_add(out) <- glue(
      s, "}} else if (cens{resp}[n] == 2) {{\n",
      s, "{tp}{w}log_diff_exp(\n", 
      s, "    {llh$dist}_lcdf(rcens{resp}[n]{llh$shift} | {llh$args}),\n", 
      s, "    {llh$dist}_lcdf({Y}{resp}[n]{llh$shift} | {llh$args})\n", 
      s, "  ){tr};\n"
    )
  }
  str_add(out) <- glue(s, "}}\n")
  out
}

# weighted likelihood in Stan language
stan_llh_weights <- function(llh, bterms, data, resp = "", ...) {
  stopifnot(is.sdist(llh))
  tr <- stan_llh_trunc(llh, bterms, data, resp = resp)
  lpdf <- stan_llh_lpdf_name(bterms)
  Y <- stan_llh_Y_name(bterms)
  glue(
    "{tp()}weights{resp}[n] * {llh$dist}_{lpdf}", 
    "({Y}{resp}[n]{llh$shift} | {llh$args}){tr};\n"
  )
}

# likelihood of a single mixture component
stan_llh_mix <- function(llh, bterms, data, mix, ptheta, resp = "", ...) {
  stopifnot(is.sdist(llh))
  theta <- str_if(ptheta,
    glue("theta{mix}{resp}[n]"), 
    glue("log(theta{mix}{resp})")
  )
  tr <- stan_llh_trunc(llh, bterms, data, resp = resp)
  lpdf <- stan_llh_lpdf_name(bterms)
  Y <- stan_llh_Y_name(bterms)
  if (is.formula(bterms$adforms$cens)) {
    # mostly copied over from stan_llh_cens
    cens <- has_cens(bterms, data = data)
    s <- wsp(nsp = 6)
    out <- glue(
      "  // special treatment of censored data\n",
      s, "if (cens{resp}[n] == 0) {{\n", 
      s, "  ps[{mix}] = {theta} + ", 
      "{llh$dist}_{lpdf}({Y}{resp}[n]{llh$shift} | {llh$args}){tr};\n",
      s, "}} else if (cens{resp}[n] == 1) {{\n",         
      s, "  ps[{mix}] = {theta} + ",
      "{llh$dist}_lccdf({Y}{resp}[n]{llh$shift} | {llh$args}){tr};\n",
      s, "}} else if (cens{resp}[n] == -1) {{\n",
      s, "  ps[{mix}] = {theta} + ",
      "{llh$dist}_lcdf({Y}{resp}[n]{llh$shift} | {llh$args}){tr};\n"
    )
    if (isTRUE(attr(cens, "interval"))) {
      str_add(out) <- glue(
        s, "}} else if (cens{resp}[n] == 2) {{\n",
        s, "  ps[{mix}] = {theta} + log_diff_exp(\n", 
        s, "    {llh$dist}_lcdf(rcens{resp}[n]{llh$shift} | {llh$args}),\n", 
        s, "    {llh$dist}_lcdf({Y}{resp}[n]{llh$shift} | {llh$args})\n", 
        s, "  ){tr};\n"
      )
    }
    str_add(out) <- glue(s, "}}\n")
  } else {
    out <- glue(
      "  ps[{mix}] = {theta} + ", 
      "{llh$dist}_{lpdf}({Y}{resp}[n]{llh$shift} | {llh$args}){tr};\n"
    ) 
  }
  out
}

# truncated part of the likelihood
# @param short use the T[, ] syntax?
stan_llh_trunc <- function(llh, bterms, data, resp = "", short = FALSE) {
  stopifnot(is.sdist(llh))
  bounds <- trunc_bounds(bterms, data = data)
  if (!any(bounds$lb > -Inf | bounds$ub < Inf)) {
    return("")
  }
  m1 <- str_if(use_int(bterms), " - 1")
  lb <- str_if(any(bounds$lb > -Inf), glue("lb{resp}[n]{m1}"))
  ub <- str_if(any(bounds$ub < Inf), glue("ub{resp}[n]"))
  if (short) {
    # truncation using T[, ] syntax
    out <- glue(" T[{lb}, {ub}]")
  } else {
    # truncation making use of _lcdf functions
    ms <- paste0(" -\n", wsp(nsp = 8))
    if (any(bounds$lb > -Inf) && !any(bounds$ub < Inf)) {
      out <- glue("{ms}{llh$dist}_lccdf({lb} | {llh$args})")
    } else if (!any(bounds$lb > -Inf) && any(bounds$ub < Inf)) {
      out <- glue("{ms}{llh$dist}_lcdf({ub} | {llh$args})")
    } else if (any(bounds$lb > -Inf) && any(bounds$ub < Inf)) {
      trr <- glue("{llh$dist}_lcdf({ub} | {llh$args})")
      trl <- glue("{llh$dist}_lcdf({lb} | {llh$args})")
      out <- glue("{ms}log_diff_exp({trr}, {trl})")
    }
  }
  out
}

stan_llh_lpdf_name <- function(bterms) {
  ifelse(use_int(bterms$family), "lpmf", "lpdf")
}

stan_llh_Y_name <- function(bterms) {
  ifelse(is.formula(bterms$adforms$mi), "Yl", "Y")
}

# prepare names of distributional parameters
# @param reqn will the likelihood be wrapped in a loop over n?
# @param dpars optional names of distributional parameters to be prepared
#   if not specified will prepare all distributional parameters
stan_llh_dpars <- function(bterms, reqn, resp = "", mix = "", dpars = NULL) {
  if (is.null(dpars)) {
    dpars <- paste0(valid_dpars(bterms), mix)
  }
  is_pred <- dpars %in% c("mu", names(bterms$dpars))
  out <- paste0(dpars, resp, ifelse(reqn & is_pred, "[n]", ""))
  named_list(dpars, out)
}

# adjust lpdf name if a more efficient version is available
# for a specific link. For instance 'poisson_log'
stan_llh_simple_lpdf <- function(lpdf, link, bterms, sep = "_") {
  stopifnot(is.brmsterms(bterms))
  cens_or_trunc <- stan_llh_adj(bterms, c("cens", "trunc"))
  if (bterms$family$link == link && !cens_or_trunc) {
    lpdf <- paste0(lpdf, sep, link)
  }
  lpdf
}

# prepare _logit suffix for distributional parameters
# used in zero-inflated and hurdle models
stan_llh_dpar_usc_logit <- function(dpar, bterms) {
  stopifnot(dpar %in% c("zi", "hu"))
  stopifnot(is.brmsterms(bterms))
  cens_or_trunc <- stan_llh_adj(bterms, c("cens", "trunc"))
  usc_logit <- isTRUE(bterms$dpars[[dpar]]$family$link == "logit")
  str_if(usc_logit && !cens_or_trunc, "_logit")
}

# prepare the code for 'sigma' in the likelihood statement
stan_llh_add_se <- function(sigma, bterms, reqn, resp = "") {
  if (is.formula(bterms$adforms$se)) {
    nse <- str_if(reqn, "[n]")
    if (no_sigma(bterms)) {
      sigma <- glue("se{resp}{nse}") 
    } else {
      sigma <- glue("sqrt(square({sigma}) + se2{resp}{nse})")
    }
  }
  sigma
}

# check if the log-liklihood needs to be adjused
# @param x named list of formulas or brmsterms object
# @param adds vector of addition argument names
# @return a single logical value
stan_llh_adj <- function(x, adds = c("weights", "cens", "trunc")) {
  stopifnot(all(adds %in% c("weights", "cens", "trunc")))
  if (is.brmsterms(x)) x <- x$adforms
  any(ulapply(x[adds], is.formula))
}

# one function per family
stan_llh_gaussian <- function(bterms, resp = "", mix = "") {
  if (use_glm_primitive(bterms)) {
    p <- args_glm_primitive(bterms$dpars$mu, resp = resp)
    p$sigma <- paste0("sigma", resp)
    out <- sdist("normal_id_glm", p$x, p$alpha, p$beta, p$sigma)
  } else {
    reqn <- stan_llh_adj(bterms) || nzchar(mix)
    p <- stan_llh_dpars(bterms, reqn, resp, mix)
    p$sigma <- stan_llh_add_se(p$sigma, bterms, reqn, resp)
    out <- sdist("normal", p$mu, p$sigma) 
  }
  out
}

stan_llh_gaussian_mv <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix) || bterms$sigma_pred
  p <- list(Mu = paste0("Mu", if (reqn) "[n]"))
  p$LSigma <- paste0("LSigma", if (bterms$sigma_pred) "[n]")
  sdist("multi_normal_cholesky", p$Mu, p$LSigma)
}

stan_llh_gaussian_cov <- function(bterms, resp = "", mix = "") {
  if (stan_llh_adj(bterms)) {
    stop2("Invalid addition arguments for this model.")
  }
  p <- stan_llh_dpars(bterms, FALSE, resp, mix)
  v <- c("chol_cor", "se2", "nobs_tg", "begin_tg", "end_tg")
  p[v] <- as.list(paste0(v, resp))
  sfx <- str_if("sigma" %in% names(bterms$dpars), "het", "hom")
  sdist(glue("normal_cov_{sfx}"), 
    p$mu, p$sigma, p$chol_cor, p$se2,
    p$nobs_tg, p$begin_tg, p$end_tg
  )
}

stan_llh_gaussian_fixed <- function(bterms, resp = "", mix = "") {
  has_se <- is.formula(bterms$adforms$se)
  if (stan_llh_adj(bterms) || has_se) {
    stop2("Invalid addition arguments for this model.")
  }
  p <- stan_llh_dpars(bterms, FALSE, resp, mix)
  p$LV <- paste0("LV", resp)
  sdist("multi_normal_cholesky", p$mu, p$LV)
}

stan_llh_gaussian_lagsar <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, FALSE, resp, mix)
  p$sigma <- stan_llh_add_se(p$sigma, bterms, FALSE, resp)
  v <- c("lagsar", "W", "eigenW")
  p[v] <- as.list(paste0(v, resp))
  sdist("normal_lagsar", p$mu, p$sigma, p$lagsar, p$W, p$eigenW)
}

stan_llh_gaussian_errorsar <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, FALSE, resp, mix)
  p$sigma <- stan_llh_add_se(p$sigma, bterms, FALSE, resp)
  v <- c("errorsar", "W", "eigenW")
  p[v] <- as.list(paste0(v, resp))
  sdist("normal_errorsar", p$mu, p$sigma, p$errorsar, p$W, p$eigenW)
}

stan_llh_student <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  p$sigma <- stan_llh_add_se(p$sigma, bterms, reqn, resp)
  sdist("student_t", p$nu, p$mu, p$sigma)
}

stan_llh_student_mv <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix) || bterms$sigma_pred
  p <- stan_llh_dpars(bterms, reqn, resp, mix, dpars = "nu")
  p$Mu <- paste0("Mu", if (reqn) "[n]")
  p$Sigma <- paste0("Sigma", if (bterms$sigma_pred) "[n]")
  sdist("multi_student_t", p$nu, p$Mu, p$Sigma)
}

stan_llh_student_cov <- function(bterms, resp = "", mix = "") {
  if (stan_llh_adj(bterms)) {
    stop2("Invalid addition arguments for this model.")
  }
  p <- stan_llh_dpars(bterms, FALSE, resp, mix)
  v <- c("chol_cor", "se2", "nobs_tg", "begin_tg", "end_tg")
  p[v] <- as.list(paste0(v, resp))
  sfx <- str_if("sigma" %in% names(bterms$dpars), "het", "hom")
  sdist(glue("student_t_cov_{sfx}"), 
    p$nu, p$mu, p$sigma, p$chol_cor, p$se2,
    p$nobs_tg, p$begin_tg, p$end_tg
  )
}

stan_llh_student_fixed <- function(bterms, resp = "", mix = "") {
  has_se <- is.formula(bterms$adforms$se)
  if (stan_llh_adj(bterms) || has_se) {
    stop2("Invalid addition arguments for this model.")
  }
  p <- stan_llh_dpars(bterms, FALSE, resp, mix)
  p$V <- paste0("V", resp)
  sdist("multi_student_t", p$nu, p$mu, p$V)
}

stan_llh_student_lagsar <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, FALSE, resp, mix)
  p$sigma <- stan_llh_add_se(p$sigma, bterms, FALSE, resp)
  v <- c("lagsar", "W", "eigenW")
  p[v] <- as.list(paste0(v, resp))
  sdist("student_t_lagsar", p$nu, p$mu, p$sigma, p$lagsar, p$W, p$eigenW)
}

stan_llh_student_errorsar <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, FALSE, resp, mix)
  p$sigma <- stan_llh_add_se(p$sigma, bterms, FALSE, resp)
  v <- c("errorsar", "W", "eigenW")
  p[v] <- as.list(paste0(v, resp))
  sdist("student_t_errorsar", p$nu, p$mu, p$sigma, p$errorsar, p$W, p$eigenW)
}

stan_llh_lognormal <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  sdist("lognormal", p$mu, p$sigma)
}

stan_llh_shifted_lognormal <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  sdist("lognormal", p$mu, p$sigma, shift = paste0(" - ", p$ndt))
}

stan_llh_asym_laplace <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  sdist("asym_laplace", p$mu, p$sigma, p$quantile)
}

stan_llh_skew_normal <- function(bterms, resp = "", mix = "", ...) {
  reqn <- stan_llh_adj(bterms) || nzchar(mix)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  p$sigma <- stan_llh_add_se(p$sigma, bterms, reqn, resp)
  # required because of CP parameterization of mu and sigma
  nomega <- any(grepl("\\[n\\]", c(p$sigma, p$alpha)))
  nomega <- str_if(reqn && nomega, "[n]")
  p$omega <- paste0("omega", mix, resp, nomega)
  sdist("skew_normal", p$mu, p$omega, p$alpha)
}

stan_llh_poisson <- function(bterms, resp = "", mix = "") {
  if (use_glm_primitive(bterms)) {
    p <- args_glm_primitive(bterms$dpars$mu, resp = resp)
    out <- sdist("poisson_log_glm", p$x, p$alpha, p$beta)
  } else {
    reqn <- stan_llh_adj(bterms) || nzchar(mix)
    p <- stan_llh_dpars(bterms, reqn, resp, mix)
    lpdf <- stan_llh_simple_lpdf("poisson", "log", bterms)
    out <- sdist(lpdf, p$mu)
  }
  out
}

stan_llh_negbinomial <- function(bterms, resp = "", mix = "") {
  if (use_glm_primitive(bterms)) {
    p <- args_glm_primitive(bterms$dpars$mu, resp = resp)
    p$shape <- paste0("shape", resp)
    out <- sdist("neg_binomial_2_log_glm", p$x, p$alpha, p$beta, p$shape)
  } else {
    reqn <- stan_llh_adj(bterms) || nzchar(mix)
    p <- stan_llh_dpars(bterms, reqn, resp, mix)
    lpdf <- stan_llh_simple_lpdf("neg_binomial_2", "log", bterms)
    out <- sdist(lpdf, p$mu, p$shape)
  }
  out
}

stan_llh_geometric <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  lpdf <- stan_llh_simple_lpdf("neg_binomial_2", "log", bterms)
  sdist(lpdf, p$mu, "1")
}

stan_llh_binomial <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  p$trials <- paste0("trials", resp, if (reqn) "[n]")
  lpdf <- stan_llh_simple_lpdf("binomial", "logit", bterms)
  sdist(lpdf, p$trials, p$mu)
}

stan_llh_bernoulli <- function(bterms, resp = "", mix = "") {
  if (use_glm_primitive(bterms)) {
    p <- args_glm_primitive(bterms$dpars$mu, resp = resp)
    out <- sdist("bernoulli_logit_glm", p$x, p$alpha, p$beta)
  } else {
    reqn <- stan_llh_adj(bterms) || nzchar(mix)
    p <- stan_llh_dpars(bterms, reqn, resp, mix)
    lpdf <- stan_llh_simple_lpdf("bernoulli", "logit", bterms)
    out <- sdist(lpdf, p$mu)
  }
  out
}

stan_llh_discrete_weibull <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  sdist("discrete_weibull", p$mu, p$shape)
}

stan_llh_com_poisson <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  lpdf <- stan_llh_simple_lpdf("com_poisson", "log", bterms)
  sdist(lpdf, p$mu, p$shape)
}

stan_llh_gamma <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  sdist("gamma", p$shape, p$mu)
}

stan_llh_exponential <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  sdist("exponential", p$mu)
}

stan_llh_weibull <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  sdist("weibull", p$shape, p$mu)
}

stan_llh_frechet <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  sdist("frechet", p$nu, p$mu)
}

stan_llh_gen_extreme_value <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  sdist("gen_extreme_value", p$mu, p$sigma, p$xi)
}

stan_llh_exgaussian <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  sdist(
    "exp_mod_normal", paste0(p$mu, " - ", p$beta), 
    p$sigma, paste0("inv(", p$beta, ")")
  )
}

stan_llh_inverse.gaussian <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  lpdf <- paste0("inv_gaussian", if (!reqn) "_vector")
  n <- str_if(reqn, "[n]")
  sdist(lpdf, p$mu, p$shape)
}

stan_llh_wiener <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  p$dec <- paste0("dec", resp, "[n]")
  sdist("wiener_diffusion", p$dec, p$bs, p$ndt, p$bias, p$mu)
}

stan_llh_beta <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix) ||
    paste0("phi", mix) %in% names(bterms$dpars)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  sdist("beta",
    paste0(p$mu, " * ", p$phi), 
    paste0("(1 - ", p$mu, ") * ", p$phi)
  )
}

stan_llh_von_mises <- function(bterms, resp = "", mix = "") {
  reqn <- stan_llh_adj(bterms) || nzchar(mix) ||
    "kappa" %in% names(bterms$dpars)
  p <- stan_llh_dpars(bterms, reqn, resp, mix)
  lpdf <- paste0("von_mises_", str_if(reqn, "real", "vector"))
  sdist(lpdf, p$mu, p$kappa)
}

stan_llh_cox <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  p$bhaz <- paste0("bhaz", resp, "[n]")
  p$cbhaz <- paste0("cbhaz", resp, "[n]")
  lpdf <- "cox"
  if (bterms$family$link == "log") {
    str_add(lpdf) <- "_log"
  }
  sdist(lpdf, p$mu, p$bhaz, p$cbhaz)
}

stan_llh_cumulative <- function(bterms, resp = "", mix = "") {
  simplify <- bterms$family$link == "logit" && 
    !"disc" %in% names(bterms$dpars) && !has_cs(bterms)
  if (simplify) {
    prefix <- paste0(resp, if (nzchar(mix)) paste0("_mu", mix))
    p <- stan_llh_dpars(bterms, TRUE, resp, mix)
    p$ord_intercept <- paste0("temp", prefix, "_Intercept")
    out <- sdist("ordered_logistic", p$mu, p$ord_intercept)
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
  stopifnot(!isTRUE(nzchar(mix)))  # mixture models are not allowed
  p <- stan_llh_dpars(bterms, TRUE, resp, mix, dpars = "mu")
  sdist("categorical_logit", p$mu)
}

stan_llh_multinomial <- function(bterms, resp = "", mix = "") {
  stopifnot(bterms$family$link == "logit")
  stopifnot(!isTRUE(nzchar(mix)))  # mixture models are not allowed
  p <- stan_llh_dpars(bterms, TRUE, resp, mix, dpars = "mu")
  sdist("multinomial_logit", p$mu)
}

stan_llh_dirichlet <- function(bterms, resp = "", mix = "") {
  stopifnot(bterms$family$link == "logit")
  stopifnot(!isTRUE(nzchar(mix)))  # mixture models are not allowed
  mu <- stan_llh_dpars(bterms, TRUE, resp, mix, dpars = "mu")$mu
  reqn <- glue("phi{mix}") %in% names(bterms$dpars)
  phi <- stan_llh_dpars(bterms, reqn, resp, mix, dpars = "phi")$phi
  sdist("dirichlet_logit", mu, phi)
}

stan_llh_ordinal <- function(bterms, resp = "", mix = "") {
  prefix <- paste0(str_if(nzchar(mix), paste0("_mu", mix)), resp)
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  p$thres <- paste0("temp", prefix, "_Intercept")
  if (has_cs(bterms)) {
    str_add(p$thres) <- paste0(" - mucs", prefix, "[n]'")
  }
  lpdf <- paste0(bterms$family$family, "_", bterms$family$link)
  sdist(lpdf, p$mu, p$thres, p$disc)
}

stan_llh_hurdle_poisson <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  lpdf <- stan_llh_simple_lpdf("hurdle_poisson", "log", bterms)
  lpdf <- paste0(lpdf, stan_llh_dpar_usc_logit("hu", bterms))
  sdist(lpdf, p$mu, p$hu)
}

stan_llh_hurdle_negbinomial <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  lpdf <- stan_llh_simple_lpdf("hurdle_neg_binomial", "log", bterms)
  lpdf <- paste0(lpdf, stan_llh_dpar_usc_logit("hu", bterms))
  sdist(lpdf, p$mu, p$shape, p$hu)
}

stan_llh_hurdle_gamma <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  usc_logit <- stan_llh_dpar_usc_logit("hu", bterms)
  lpdf <- paste0("hurdle_gamma", usc_logit)
  sdist(lpdf, p$shape, p$mu, p$hu)
}

stan_llh_hurdle_lognormal <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  usc_logit <- stan_llh_dpar_usc_logit("hu", bterms)
  lpdf <- paste0("hurdle_lognormal", usc_logit)
  sdist(lpdf, p$mu, p$sigma, p$hu)
}

stan_llh_zero_inflated_poisson <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  lpdf <- stan_llh_simple_lpdf("zero_inflated_poisson", "log", bterms)
  lpdf <- paste0(lpdf, stan_llh_dpar_usc_logit("zi", bterms))
  sdist(lpdf, p$mu, p$zi)
}

stan_llh_zero_inflated_negbinomial <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  lpdf <- stan_llh_simple_lpdf("zero_inflated_neg_binomial", "log", bterms)
  lpdf <- paste0(lpdf, stan_llh_dpar_usc_logit("zi", bterms))
  sdist(lpdf, p$mu, p$shape, p$zi)
}

stan_llh_zero_inflated_binomial <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  p$trials <- "trials[n]"
  lpdf <- "zero_inflated_binomial"
  lpdf <- stan_llh_simple_lpdf(lpdf, "logit", bterms, sep = "_b")
  lpdf <- paste0(lpdf, stan_llh_dpar_usc_logit("zi", bterms))
  sdist(lpdf, p$trials, p$mu, p$zi)
}

stan_llh_zero_inflated_beta <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  usc_logit <- stan_llh_dpar_usc_logit("zi", bterms)
  lpdf <- paste0("zero_inflated_beta", usc_logit)
  sdist(lpdf, p$mu, p$phi, p$zi)
}

stan_llh_zero_one_inflated_beta <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  sdist("zero_one_inflated_beta", p$mu, p$phi, p$zoi, p$coi)
}

stan_llh_custom <- function(bterms, resp = "", mix = "") {
  p <- stan_llh_dpars(bterms, TRUE, resp, mix)
  family <- bterms$family
  dpars <- paste0(family$dpars, mix)
  if (is_ordinal(family)) {
    prefix <- paste0(resp, if (nzchar(mix)) paste0("_mu", mix))
    p$ord_intercept <- paste0("temp", prefix, "_Intercept")
  }
  sdist(family$name, p[dpars], p$ord_intercept, family$vars)
}

# use Stan GLM primitive functions?
# @param bterms a brmsterms object
# @return TRUE or FALSE
use_glm_primitive <- function(bterms) {
  stopifnot(is.brmsterms(bterms))
  # the model can only have a single predicted parameter
  # and no additional residual or autocorrelation structure
  mu <- bterms$dpars[["mu"]]
  if (!is.btl(mu) || length(bterms$dpars) > 1L ||
      isTRUE(bterms$rescor) || length(bterms$adforms) ||
      !is.cor_empty(bterms$autocor)) {
    return(FALSE)
  }
  # supported families and link functions
  glm_links <- list(
    gaussian = "identity", bernoulli = "logit",
    poisson = "log", negbinomial = "log"
  )
  if (!isTRUE(glm_links[[mu$family$family]] == mu$family$link)) {
    return(FALSE)
  }
  # can only use GLM primitives if solely 'fixed effects' are present
  special_term_names <- c("sp", "cs", "sm", "gp", "offset")
  length(all_terms(mu$fe)) && !is_sparse(mu$fe) &&
    !NROW(mu$re) && !any(lengths(mu[special_term_names]))
}

# standard arguments for primitive Stan GLM functions
# @param bterms a btl object
# @param resp optional name of the response variable
# @return a named list of Stan code snippets
args_glm_primitive <- function(bterms, resp = "") {
  stopifnot(is.btl(bterms))
  decomp <- get_decomp(bterms$fe)
  center_X <- stan_center_X(bterms)
  sfx_X <- sfx_b <- ""
  if (decomp == "QR") {
    sfx_X <- sfx_b <- "Q"
  } else if (center_X) {
    sfx_X <- "c"
  }
  if (center_X) {
    intercept <- paste0("temp", resp, "_Intercept")
  } else {
    intercept <- "0"
  }
  list(
    x = paste0("X", sfx_X, resp),
    alpha = intercept,
    beta = paste0("b", sfx_b, resp)
  )
}

# prepare distribution and arguments for use in Stan
sdist <- function(dist, ..., shift = "") {
  args <- sargs(...)
  structure(nlist(dist, args, shift), class = "sdist")
}

# prepare arguments for Stan likelihood statements
sargs <- function(...) {
  dots <- as.character(c(...))
  dots <- dots[nzchar(dots)]
  paste0(dots, collapse = ", ")
}

is.sdist <- function(x) {
  inherits(x, "sdist")
}

tp <- function(wsp = 2) {
  wsp <- collapse(rep(" ", wsp))
  paste0(wsp, "target += ")
}
