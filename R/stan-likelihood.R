# unless otherwise specified, functions return a single character 
# string defining the likelihood of the model in Stan language

# TODO: adjust white spaces
stan_log_lik <- function(x, ...) {
  UseMethod("stan_log_lik")
}

# Stan code for the model likelihood 
# @param bterms object of class brmsterms
# @param data data passed by the user
# @param mix optional mixture component ID
# @param ptheta are mixing proportions predicted?
#' @export
stan_log_lik.family <- function(x, bterms, data, threads, 
                                mix = "", ptheta = FALSE, ...) {
  stopifnot(is.brmsterms(bterms))
  stopifnot(length(mix) == 1L)
  bterms$family <- x
  resp <- usc(combine_prefix(bterms))
  # prepare family part of the likelihood
  log_lik_args <- nlist(bterms, resp, mix, threads)
  log_lik_fun <- paste0("stan_log_lik_", prepare_family(bterms)$fun)
  ll <- do_call(log_lik_fun, log_lik_args)
  # incorporate other parts into the likelihood
  args <- nlist(ll, bterms, data, resp, threads, mix, ptheta)
  if (nzchar(mix)) {
    out <- do_call(stan_log_lik_mix, args)
  } else if (is.formula(bterms$adforms$cens)) {
    out <- do_call(stan_log_lik_cens, args)
  } else if (is.formula(bterms$adforms$weights)) {
    out <- do_call(stan_log_lik_weights, args)
  } else {
    out <- do_call(stan_log_lik_general, args)
  }
  if (grepl("\\[n\\]", out) && !nzchar(mix)) {
    # loop over likelihood if it cannot be vectorized
    out <- paste0(
      "  for (n in 1:N", resp, ") {\n", 
      "  ", stan_nn_def(threads), 
      "  ", out, 
      "  }\n"
    )
  }
  out
}

#' @export
stan_log_lik.mixfamily <- function(x, bterms, threads, ...) {
  dp_ids <- dpar_id(names(bterms$dpars))
  fdp_ids <- dpar_id(names(bterms$fdpars))
  resp <- usc(bterms$resp)
  ptheta <- any(dpar_class(names(bterms$dpars)) %in% "theta")
  ll <- rep(NA, length(x$mix))
  for (i in seq_along(x$mix)) {
    sbterms <- bterms
    sbterms$dpars <- sbterms$dpars[dp_ids == i]
    sbterms$fdpars <- sbterms$fdpars[fdp_ids == i]
    ll[i] <- stan_log_lik(
      x$mix[[i]], sbterms, mix = i, ptheta = ptheta, 
      threads = threads, ...
    )
  }
  resp <- usc(combine_prefix(bterms))
  n <- stan_nn(threads)
  has_weights <- is.formula(bterms$adforms$weights)  
  weights <- str_if(has_weights, glue("weights{resp}{n} * "))
  out <- glue(
    "  // likelihood of the mixture model\n",
    "  for (n in 1:N{resp}) {{\n",
    stan_nn_def(threads),
    "    real ps[{length(ll)}];\n"
  )
  str_add(out) <- collapse("    ", ll)
  str_add(out) <- glue(
    "  {tp()}{weights}log_sum_exp(ps);\n",
    "  }}\n"
  )
  out
}

#' @export
stan_log_lik.brmsterms <- function(x, ...) {
  stan_log_lik(x$family, bterms = x, ...)
}

#' @export
stan_log_lik.mvbrmsterms <- function(x, ...) {
  if (x$rescor) {
    out <- stan_log_lik(as.brmsterms(x), ...)
  } else {
    out <- ulapply(x$terms, stan_log_lik, ...) 
  }
  out
}

# default likelihood in Stan language
stan_log_lik_general <- function(ll, bterms, data, threads, resp = "", ...) {
  stopifnot(is.sdist(ll))
  require_n <- grepl("\\[n\\]", ll$args)
  n <- str_if(require_n, stan_nn(threads), stan_slice(threads))
  lpdf <- stan_log_lik_lpdf_name(bterms)
  Y <- stan_log_lik_Y_name(bterms)
  tr <- stan_log_lik_trunc(ll, bterms, data, resp = resp, threads = threads)
  glue("{tp()}{ll$dist}_{lpdf}({Y}{resp}{n}{ll$shift} | {ll$args}){tr};\n")
}

# censored likelihood in Stan language
stan_log_lik_cens <- function(ll, bterms, data, threads, resp = "", ...) {
  stopifnot(is.sdist(ll))
  s <- wsp(nsp = 4)
  cens <- eval_rhs(bterms$adforms$cens)
  lpdf <- stan_log_lik_lpdf_name(bterms)
  has_weights <- is.formula(bterms$adforms$weights)
  Y <- stan_log_lik_Y_name(bterms)
  n <- stan_nn(threads)
  w <- str_if(has_weights, glue("weights{resp}{n} * "))
  tr <- stan_log_lik_trunc(ll, bterms, data, resp = resp, threads = threads)
  tp <- tp()
  out <- glue(
    "// special treatment of censored data\n",
    s, "if (cens{resp}{n} == 0) {{\n", 
    s, "{tp}{w}{ll$dist}_{lpdf}({Y}{resp}{n}{ll$shift} | {ll$args}){tr};\n",
    s, "}} else if (cens{resp}{n} == 1) {{\n",         
    s, "{tp}{w}{ll$dist}_lccdf({Y}{resp}{n}{ll$shift} | {ll$args}){tr};\n",
    s, "}} else if (cens{resp}{n} == -1) {{\n",
    s, "{tp}{w}{ll$dist}_lcdf({Y}{resp}{n}{ll$shift} | {ll$args}){tr};\n"
  )
  if (cens$vars$y2 != "NA") {
    # interval censoring is required
    str_add(out) <- glue(
      s, "}} else if (cens{resp}{n} == 2) {{\n",
      s, "{tp}{w}log_diff_exp(\n", 
      s, "    {ll$dist}_lcdf(rcens{resp}{n}{ll$shift} | {ll$args}),\n", 
      s, "    {ll$dist}_lcdf({Y}{resp}{n}{ll$shift} | {ll$args})\n", 
      s, "  ){tr};\n"
    )
  }
  str_add(out) <- glue(s, "}}\n")
  out
}

# weighted likelihood in Stan language
stan_log_lik_weights <- function(ll, bterms, data, threads, resp = "", ...) {
  stopifnot(is.sdist(ll))
  tr <- stan_log_lik_trunc(ll, bterms, data, resp = resp, threads = threads)
  lpdf <- stan_log_lik_lpdf_name(bterms)
  Y <- stan_log_lik_Y_name(bterms)
  n <- stan_nn(threads)
  glue(
    "{tp()}weights{resp}{n} * ({ll$dist}_{lpdf}", 
    "({Y}{resp}{n}{ll$shift} | {ll$args}){tr});\n"
  )
}

# likelihood of a single mixture component
stan_log_lik_mix <- function(ll, bterms, data, mix, ptheta, threads, 
                             resp = "", ...) {
  stopifnot(is.sdist(ll))
  theta <- str_if(ptheta,
    glue("theta{mix}{resp}[n]"), 
    glue("log(theta{mix}{resp})")
  )
  tr <- stan_log_lik_trunc(ll, bterms, data, resp = resp, threads = threads)
  lpdf <- stan_log_lik_lpdf_name(bterms)
  Y <- stan_log_lik_Y_name(bterms)
  n <- stan_nn(threads)
  if (is.formula(bterms$adforms$cens)) {
    # mostly copied over from stan_log_lik_cens
    cens <- eval_rhs(bterms$adforms$cens)
    s <- wsp(nsp = 4)
    out <- glue(
      "// special treatment of censored data\n",
      s, "if (cens{resp}{n} == 0) {{\n", 
      s, "  ps[{mix}] = {theta} + ", 
      "{ll$dist}_{lpdf}({Y}{resp}{n}{ll$shift} | {ll$args}){tr};\n",
      s, "}} else if (cens{resp}{n} == 1) {{\n",         
      s, "  ps[{mix}] = {theta} + ",
      "{ll$dist}_lccdf({Y}{resp}{n}{ll$shift} | {ll$args}){tr};\n",
      s, "}} else if (cens{resp}{n} == -1) {{\n",
      s, "  ps[{mix}] = {theta} + ",
      "{ll$dist}_lcdf({Y}{resp}{n}{ll$shift} | {ll$args}){tr};\n"
    )
    if (cens$vars$y2 != "NA") {
      # interval censoring is required
      str_add(out) <- glue(
        s, "}} else if (cens{resp}{n} == 2) {{\n",
        s, "  ps[{mix}] = {theta} + log_diff_exp(\n", 
        s, "    {ll$dist}_lcdf(rcens{resp}{n}{ll$shift} | {ll$args}),\n", 
        s, "    {ll$dist}_lcdf({Y}{resp}{n}{ll$shift} | {ll$args})\n", 
        s, "  ){tr};\n"
      )
    }
    str_add(out) <- glue(s, "}}\n")
  } else {
    out <- glue(
      "  ps[{mix}] = {theta} + ", 
      "{ll$dist}_{lpdf}({Y}{resp}{n}{ll$shift} | {ll$args}){tr};\n"
    ) 
  }
  out
}

# truncated part of the likelihood
# @param short use the T[, ] syntax?
stan_log_lik_trunc <- function(ll, bterms, data, threads,resp = "", 
                               short = FALSE) {
  stopifnot(is.sdist(ll))
  bounds <- trunc_bounds(bterms, data = data)
  if (!any(bounds$lb > -Inf | bounds$ub < Inf)) {
    return("")
  }
  n <- stan_nn(threads)
  m1 <- str_if(use_int(bterms), " - 1")
  lb <- str_if(any(bounds$lb > -Inf), glue("lb{resp}{n}{m1}"))
  ub <- str_if(any(bounds$ub < Inf), glue("ub{resp}{n}"))
  if (short) {
    # truncation using T[, ] syntax
    out <- glue(" T[{lb}, {ub}]")
  } else {
    # truncation making use of _lcdf functions
    ms <- paste0(" -\n", wsp(nsp = 6))
    if (any(bounds$lb > -Inf) && !any(bounds$ub < Inf)) {
      out <- glue("{ms}{ll$dist}_lccdf({lb} | {ll$args})")
    } else if (!any(bounds$lb > -Inf) && any(bounds$ub < Inf)) {
      out <- glue("{ms}{ll$dist}_lcdf({ub} | {ll$args})")
    } else if (any(bounds$lb > -Inf) && any(bounds$ub < Inf)) {
      trr <- glue("{ll$dist}_lcdf({ub} | {ll$args})")
      trl <- glue("{ll$dist}_lcdf({lb} | {ll$args})")
      out <- glue("{ms}log_diff_exp({trr}, {trl})")
    }
  }
  out
}

stan_log_lik_lpdf_name <- function(bterms) {
  ifelse(use_int(bterms$family), "lpmf", "lpdf")
}

stan_log_lik_Y_name <- function(bterms) {
  ifelse(is.formula(bterms$adforms$mi), "Yl", "Y")
}

# prepare names of distributional parameters
# @param reqn will the likelihood be wrapped in a loop over n?
# @param dpars optional names of distributional parameters to be prepared
#   if not specified will prepare all distributional parameters
stan_log_lik_dpars <- function(bterms, reqn, resp = "", mix = "", dpars = NULL) {
  if (is.null(dpars)) {
    dpars <- paste0(valid_dpars(bterms), mix)
  }
  is_pred <- dpars %in% c("mu", names(bterms$dpars))
  out <- paste0(dpars, resp, ifelse(reqn & is_pred, "[n]", ""))
  named_list(dpars, out)
}

# adjust lpdf name if a more efficient version is available
# for a specific link. For instance 'poisson_log'
stan_log_lik_simple_lpdf <- function(lpdf, link, bterms, sep = "_") {
  stopifnot(is.brmsterms(bterms))
  cens_or_trunc <- stan_log_lik_adj(bterms, c("cens", "trunc"))
  if (bterms$family$link == link && !cens_or_trunc) {
    lpdf <- paste0(lpdf, sep, link)
  }
  lpdf
}

# prepare _logit suffix for distributional parameters
# used in zero-inflated and hurdle models
stan_log_lik_dpar_usc_logit <- function(dpar, bterms) {
  stopifnot(dpar %in% c("zi", "hu"))
  stopifnot(is.brmsterms(bterms))
  cens_or_trunc <- stan_log_lik_adj(bterms, c("cens", "trunc"))
  usc_logit <- isTRUE(bterms$dpars[[dpar]]$family$link == "logit")
  str_if(usc_logit && !cens_or_trunc, "_logit")
}

# add 'se' to 'sigma' within the Stan likelihood
stan_log_lik_add_se <- function(sigma, bterms, reqn, resp = "") {
  if (!is.formula(bterms$adforms$se)) {
    return(sigma) 
  }
  nse <- str_if(reqn, "[n]")
  if (no_sigma(bterms)) {
    sigma <- glue("se{resp}{nse}") 
  } else {
    sigma <- glue("sqrt(square({sigma}) + se2{resp}{nse})")
  }
  sigma
}

# multiply 'dpar' by the 'rate' denominator within the Stan likelihood
# @param log add the rate denominator on the log scale if sensible?
stan_log_lik_multiply_rate_denom <- function(dpar, bterms, reqn, resp = "", 
                                         log = FALSE) {
  if (!is.formula(bterms$adforms$rate)) {
    return(dpar)
  }
  ndenom <- str_if(reqn, "[n]")
  denom <- glue("denom{resp}{ndenom}")
  cens_or_trunc <- stan_log_lik_adj(bterms, c("cens", "trunc"))
  if (log && bterms$family$link == "log" && !cens_or_trunc) {
    denom <- glue("log_{denom}")
    operator <- "+"
  } else {
    is_pred <- dpar %in% c("mu", names(bterms$dpars))
    operator <- str_if(reqn || !is_pred, "*", ".*")
  }
  glue("{dpar} {operator} {denom}")
}

# check if the log-liklihood needs to be adjused
# @param x named list of formulas or brmsterms object
# @param adds vector of addition argument names
# @return a single logical value
stan_log_lik_adj <- function(x, adds = c("weights", "cens", "trunc")) {
  stopifnot(all(adds %in% c("weights", "cens", "trunc")))
  if (is.brmsterms(x)) x <- x$adforms
  any(ulapply(x[adds], is.formula))
}

# one function per family
stan_log_lik_gaussian <- function(bterms, resp = "", mix = "", threads = 1,
                                  ...) {
  if (use_glm_primitive(bterms)) {
    p <- args_glm_primitive(bterms$dpars$mu, resp = resp, threads = threads)
    p$sigma <- paste0("sigma", resp)
    out <- sdist("normal_id_glm", p$x, p$alpha, p$beta, p$sigma)
  } else {
    reqn <- stan_log_lik_adj(bterms) || nzchar(mix)
    p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
    p$sigma <- stan_log_lik_add_se(p$sigma, bterms, reqn, resp)
    out <- sdist("normal", p$mu, p$sigma) 
  }
  out
}

stan_log_lik_gaussian_mv <- function(bterms, resp = "", mix = "", ...) {
  reqn <- stan_log_lik_adj(bterms) || nzchar(mix) || bterms$sigma_pred
  p <- list(Mu = paste0("Mu", if (reqn) "[n]"))
  p$LSigma <- paste0("LSigma", if (bterms$sigma_pred) "[n]")
  sdist("multi_normal_cholesky", p$Mu, p$LSigma)
}

stan_log_lik_gaussian_time <- function(bterms, resp = "", mix = "", ...) {
  if (stan_log_lik_adj(bterms)) {
    stop2("Invalid addition arguments for this model.")
  }
  p <- stan_log_lik_dpars(bterms, FALSE, resp, mix)
  v <- c("chol_cor", "se2", "nobs_tg", "begin_tg", "end_tg")
  p[v] <- as.list(paste0(v, resp))
  sfx <- str_if("sigma" %in% names(bterms$dpars), "het", "hom")
  sdist(glue("normal_time_{sfx}"), 
    p$mu, p$sigma, p$chol_cor, p$se2,
    p$nobs_tg, p$begin_tg, p$end_tg
  )
}

stan_log_lik_gaussian_fcor <- function(bterms, resp = "", mix = "", ...) {
  has_se <- is.formula(bterms$adforms$se)
  if (stan_log_lik_adj(bterms) || has_se) {
    stop2("Invalid addition arguments for this model.")
  }
  p <- stan_log_lik_dpars(bterms, FALSE, resp, mix)
  p$Lfcor <- paste0("Lfcor", resp)
  sfx <- str_if("sigma" %in% names(bterms$dpars), "het", "hom")
  sdist(glue("normal_fcor_{sfx}"), p$mu, p$sigma, p$Lfcor)
}

stan_log_lik_gaussian_lagsar <- function(bterms, resp = "", mix = "", ...) {
  p <- stan_log_lik_dpars(bterms, FALSE, resp, mix)
  p$sigma <- stan_log_lik_add_se(p$sigma, bterms, FALSE, resp)
  v <- c("lagsar", "Msar", "eigenMsar")
  p[v] <- as.list(paste0(v, resp))
  sdist("normal_lagsar", p$mu, p$sigma, p$lagsar, p$Msar, p$eigenMsar)
}

stan_log_lik_gaussian_errorsar <- function(bterms, resp = "", mix = "", ...) {
  p <- stan_log_lik_dpars(bterms, FALSE, resp, mix)
  p$sigma <- stan_log_lik_add_se(p$sigma, bterms, FALSE, resp)
  v <- c("errorsar", "Msar", "eigenMsar")
  p[v] <- as.list(paste0(v, resp))
  sdist("normal_errorsar", p$mu, p$sigma, p$errorsar, p$Msar, p$eigenMsar)
}

stan_log_lik_student <- function(bterms, resp = "", mix = "", ...) {
  reqn <- stan_log_lik_adj(bterms) || nzchar(mix)
  p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
  p$sigma <- stan_log_lik_add_se(p$sigma, bterms, reqn, resp)
  sdist("student_t", p$nu, p$mu, p$sigma)
}

stan_log_lik_student_mv <- function(bterms, resp = "", mix = "", ...) {
  reqn <- stan_log_lik_adj(bterms) || nzchar(mix) || bterms$sigma_pred
  p <- stan_log_lik_dpars(bterms, reqn, resp, mix, dpars = "nu")
  p$Mu <- paste0("Mu", if (reqn) "[n]")
  p$Sigma <- paste0("Sigma", if (bterms$sigma_pred) "[n]")
  sdist("multi_student_t", p$nu, p$Mu, p$Sigma)
}

stan_log_lik_student_time <- function(bterms, resp = "", mix = "", ...) {
  if (stan_log_lik_adj(bterms)) {
    stop2("Invalid addition arguments for this model.")
  }
  p <- stan_log_lik_dpars(bterms, FALSE, resp, mix)
  v <- c("chol_cor", "se2", "nobs_tg", "begin_tg", "end_tg")
  p[v] <- as.list(paste0(v, resp))
  sfx <- str_if("sigma" %in% names(bterms$dpars), "het", "hom")
  sdist(glue("student_t_time_{sfx}"), 
    p$nu, p$mu, p$sigma, p$chol_cor, p$se2,
    p$nobs_tg, p$begin_tg, p$end_tg
  )
}

stan_log_lik_student_fcor <- function(bterms, resp = "", mix = "", ...) {
  has_se <- is.formula(bterms$adforms$se)
  if (stan_log_lik_adj(bterms) || has_se) {
    stop2("Invalid addition arguments for this model.")
  }
  p <- stan_log_lik_dpars(bterms, FALSE, resp, mix)
  p$Lfcor <- paste0("Lfcor", resp)
  sfx <- str_if("sigma" %in% names(bterms$dpars), "het", "hom")
  sdist(glue("student_t_fcor_{sfx}"), p$nu, p$mu, p$sigma, p$Lfcor)
}

stan_log_lik_student_lagsar <- function(bterms, resp = "", mix = "", ...) {
  p <- stan_log_lik_dpars(bterms, FALSE, resp, mix)
  p$sigma <- stan_log_lik_add_se(p$sigma, bterms, FALSE, resp)
  v <- c("lagsar", "Msar", "eigenMsar")
  p[v] <- as.list(paste0(v, resp))
  sdist("student_t_lagsar", p$nu, p$mu, p$sigma, 
        p$lagsar, p$Msar, p$eigenMsar)
}

stan_log_lik_student_errorsar <- function(bterms, resp = "", mix = "", ...) {
  p <- stan_log_lik_dpars(bterms, FALSE, resp, mix)
  p$sigma <- stan_log_lik_add_se(p$sigma, bterms, FALSE, resp)
  v <- c("errorsar", "Msar", "eigenMsar")
  p[v] <- as.list(paste0(v, resp))
  sdist("student_t_errorsar", p$nu, p$mu, p$sigma, 
        p$errorsar, p$Msar, p$eigenMsar)
}

stan_log_lik_lognormal <- function(bterms, resp = "", mix = "", ...) {
  reqn <- stan_log_lik_adj(bterms) || nzchar(mix)
  p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
  sdist("lognormal", p$mu, p$sigma)
}

stan_log_lik_shifted_lognormal <- function(bterms, resp = "", mix = "", ...) {
  reqn <- stan_log_lik_adj(bterms) || nzchar(mix)
  p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
  sdist("lognormal", p$mu, p$sigma, shift = paste0(" - ", p$ndt))
}

stan_log_lik_asym_laplace <- function(bterms, resp = "", mix = "", ...) {
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix)
  sdist("asym_laplace", p$mu, p$sigma, p$quantile)
}

stan_log_lik_skew_normal <- function(bterms, resp = "", mix = "", ...) {
  reqn <- stan_log_lik_adj(bterms) || nzchar(mix)
  p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
  p$sigma <- stan_log_lik_add_se(p$sigma, bterms, reqn, resp)
  # required because of CP parameterization of mu and sigma
  nomega <- any(grepl("\\[n\\]", c(p$sigma, p$alpha)))
  nomega <- str_if(reqn && nomega, "[n]")
  p$omega <- paste0("omega", mix, resp, nomega)
  sdist("skew_normal", p$mu, p$omega, p$alpha)
}

stan_log_lik_poisson <- function(bterms, resp = "", mix = "", threads = 1,
                                 ...) {
  if (use_glm_primitive(bterms)) {
    p <- args_glm_primitive(bterms$dpars$mu, resp = resp, threads = threads)
    out <- sdist("poisson_log_glm", p$x, p$alpha, p$beta)
  } else {
    reqn <- stan_log_lik_adj(bterms) || nzchar(mix)
    p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
    lpdf <- stan_log_lik_simple_lpdf("poisson", "log", bterms)
    p$mu <- stan_log_lik_multiply_rate_denom(p$mu, bterms, reqn, resp, log = TRUE)
    out <- sdist(lpdf, p$mu)
  }
  out
}

stan_log_lik_negbinomial <- function(bterms, resp = "", mix = "", threads = 1,
                                     ...) {
  if (use_glm_primitive(bterms)) {
    p <- args_glm_primitive(bterms$dpars$mu, resp = resp, threads = threads)
    p$shape <- paste0("shape", resp)
    out <- sdist("neg_binomial_2_log_glm", p$x, p$alpha, p$beta, p$shape)
  } else {
    reqn <- stan_log_lik_adj(bterms) || nzchar(mix)
    p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
    p$mu <- stan_log_lik_multiply_rate_denom(p$mu, bterms, reqn, resp, log = TRUE)
    p$shape <- stan_log_lik_multiply_rate_denom(p$shape, bterms, reqn, resp)
    lpdf <- stan_log_lik_simple_lpdf("neg_binomial_2", "log", bterms)
    out <- sdist(lpdf, p$mu, p$shape)
  }
  out
}

stan_log_lik_geometric <- function(bterms, resp = "", mix = "", threads = 1, 
                                   ...) {
  if (use_glm_primitive(bterms)) {
    p <- args_glm_primitive(bterms$dpars$mu, resp = resp, threads = threads)
    p$shape <- "1"
    out <- sdist("neg_binomial_2_log_glm", p$x, p$alpha, p$beta, p$shape)
  } else {
    reqn <- stan_log_lik_adj(bterms) || nzchar(mix)
    p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
    p$shape <- "1"
    p$mu <- stan_log_lik_multiply_rate_denom(p$mu, bterms, reqn, resp, log = TRUE)
    p$shape <- stan_log_lik_multiply_rate_denom(p$shape, bterms, reqn, resp)
    lpdf <- stan_log_lik_simple_lpdf("neg_binomial_2", "log", bterms)
    out <- sdist(lpdf, p$mu, p$shape)
  }
}

stan_log_lik_binomial <- function(bterms, resp = "", mix = "", threads = 1, 
                                  ...) {
  reqn <- stan_log_lik_adj(bterms) || nzchar(mix)
  p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
  slice <- str_if(reqn, stan_nn(threads), stan_slice(threads))
  p$trials <- paste0("trials", resp, slice)
  lpdf <- stan_log_lik_simple_lpdf("binomial", "logit", bterms)
  sdist(lpdf, p$trials, p$mu)
}

stan_log_lik_bernoulli <- function(bterms, resp = "", mix = "", threads = 1, 
                                   ...) {
  if (use_glm_primitive(bterms)) {
    p <- args_glm_primitive(bterms$dpars$mu, resp = resp, threads = threads)
    out <- sdist("bernoulli_logit_glm", p$x, p$alpha, p$beta)
  } else {
    reqn <- stan_log_lik_adj(bterms) || nzchar(mix)
    p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
    lpdf <- stan_log_lik_simple_lpdf("bernoulli", "logit", bterms)
    out <- sdist(lpdf, p$mu)
  }
  out
}

stan_log_lik_discrete_weibull <- function(bterms, resp = "", mix = "", ...) {
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix)
  sdist("discrete_weibull", p$mu, p$shape)
}

stan_log_lik_com_poisson <- function(bterms, resp = "", mix = "", ...) {
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix)
  lpdf <- stan_log_lik_simple_lpdf("com_poisson", "log", bterms)
  sdist(lpdf, p$mu, p$shape)
}

stan_log_lik_gamma <- function(bterms, resp = "", mix = "", ...) {
  reqn <- stan_log_lik_adj(bterms) || nzchar(mix)
  p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
  sdist("gamma", p$shape, p$mu)
}

stan_log_lik_exponential <- function(bterms, resp = "", mix = "", ...) {
  reqn <- stan_log_lik_adj(bterms) || nzchar(mix)
  p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
  sdist("exponential", p$mu)
}

stan_log_lik_weibull <- function(bterms, resp = "", mix = "", ...) {
  reqn <- stan_log_lik_adj(bterms) || nzchar(mix)
  p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
  sdist("weibull", p$shape, p$mu)
}

stan_log_lik_frechet <- function(bterms, resp = "", mix = "", ...) {
  reqn <- stan_log_lik_adj(bterms) || nzchar(mix)
  p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
  sdist("frechet", p$nu, p$mu)
}

stan_log_lik_gen_extreme_value <- function(bterms, resp = "", mix = "", ...) {
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix)
  sdist("gen_extreme_value", p$mu, p$sigma, p$xi)
}

stan_log_lik_exgaussian <- function(bterms, resp = "", mix = "", ...) {
  reqn <- stan_log_lik_adj(bterms) || nzchar(mix)
  p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
  sdist(
    "exp_mod_normal", paste0(p$mu, " - ", p$beta), 
    p$sigma, paste0("inv(", p$beta, ")")
  )
}

stan_log_lik_inverse.gaussian <- function(bterms, resp = "", mix = "", ...) {
  reqn <- stan_log_lik_adj(bterms) || nzchar(mix) ||
    glue("shape{mix}") %in% names(bterms$dpars)
  p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
  lpdf <- paste0("inv_gaussian", if (!reqn) "vector")
  n <- str_if(reqn, "[n]")
  sdist(lpdf, p$mu, p$shape)
}

stan_log_lik_wiener <- function(bterms, resp = "", mix = "", ...) {
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix)
  p$dec <- paste0("dec", resp, "[n]")
  sdist("wiener_diffusion", p$dec, p$bs, p$ndt, p$bias, p$mu)
}

stan_log_lik_beta <- function(bterms, resp = "", mix = "", ...) {
  reqn <- stan_log_lik_adj(bterms) || nzchar(mix) ||
    paste0("phi", mix) %in% names(bterms$dpars)
  p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
  sdist("beta",
    paste0(p$mu, " * ", p$phi), 
    paste0("(1 - ", p$mu, ") * ", p$phi)
  )
}

stan_log_lik_von_mises <- function(bterms, resp = "", mix = "", ...) {
  reqn <- stan_log_lik_adj(bterms) || nzchar(mix) ||
    "kappa" %in% names(bterms$dpars)
  p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
  lpdf <- paste0("von_mises_", str_if(reqn, "real", "vector"))
  sdist(lpdf, p$mu, p$kappa)
}

stan_log_lik_cox <- function(bterms, resp = "", mix = "", threads = 1,
                             ...) {
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix)
  n <- stan_nn(threads)
  p$bhaz <- paste0("bhaz", resp, n)
  p$cbhaz <- paste0("cbhaz", resp, n)
  lpdf <- "cox"
  if (bterms$family$link == "log") {
    str_add(lpdf) <- "_log"
  }
  sdist(lpdf, p$mu, p$bhaz, p$cbhaz)
}

stan_log_lik_cumulative <- function(bterms, resp = "", mix = "", ...) {
  stan_log_lik_ordinal(bterms, resp, mix)
}

stan_log_lik_sratio <- function(bterms, resp = "", mix = "", ...) {
  stan_log_lik_ordinal(bterms, resp, mix)
}

stan_log_lik_cratio <- function(bterms, resp = "", mix = "", ...) {
  stan_log_lik_ordinal(bterms, resp, mix)
}

stan_log_lik_acat <- function(bterms, resp = "", mix = "", ...) {
  stan_log_lik_ordinal(bterms, resp, mix)
}

stan_log_lik_categorical <- function(bterms, resp = "", mix = "", ...) {
  stopifnot(bterms$family$link == "logit")
  stopifnot(!isTRUE(nzchar(mix)))  # mixture models are not allowed
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix, dpars = "mu")
  sdist("categorical_logit", p$mu)
}

stan_log_lik_multinomial <- function(bterms, resp = "", mix = "", ...) {
  stopifnot(bterms$family$link == "logit")
  stopifnot(!isTRUE(nzchar(mix)))  # mixture models are not allowed
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix, dpars = "mu")
  sdist("multinomial_logit", p$mu)
}

stan_log_lik_dirichlet <- function(bterms, resp = "", mix = "", ...) {
  stopifnot(bterms$family$link == "logit")
  stopifnot(!isTRUE(nzchar(mix)))  # mixture models are not allowed
  mu <- stan_log_lik_dpars(bterms, TRUE, resp, mix, dpars = "mu")$mu
  reqn <- glue("phi{mix}") %in% names(bterms$dpars)
  phi <- stan_log_lik_dpars(bterms, reqn, resp, mix, dpars = "phi")$phi
  sdist("dirichlet_logit", mu, phi)
}

stan_log_lik_ordinal <- function(bterms, resp = "", mix = "", ...) {
  prefix <- paste0(str_if(nzchar(mix), paste0("_mu", mix)), resp)
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix)
  if (use_ordered_logistic(bterms)) {
    lpdf <- "ordered_logistic"
    p[grepl("^disc", names(p))] <- NULL
  } else {
    lpdf <- paste0(bterms$family$family, "_", bterms$family$link) 
  }
  if (has_thres_groups(bterms)) {
    str_add(lpdf) <- "_merged"
    p$Jthres <- paste0("Jthres", resp, "[n]")
    p$thres <- "merged_Intercept"
  } else {
    p$thres <- "Intercept"
  }
  str_add(p$thres) <- prefix
  if (has_sum_to_zero_thres(bterms)) {
    str_add(p$thres) <- "_stz"
  }
  if (has_cs(bterms)) {
    if (has_thres_groups(bterms)) {
      stop2("Cannot use category specific effects ", 
            "in models with multiple thresholds.")
    }
    str_add(p$thres) <- paste0(" - mucs", prefix, "[n]'")
  }
  sdist(lpdf, p$mu, p$disc, p$thres, p$Jthres)
}

stan_log_lik_hurdle_poisson <- function(bterms, resp = "", mix = "", ...) {
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix)
  lpdf <- stan_log_lik_simple_lpdf("hurdle_poisson", "log", bterms)
  lpdf <- paste0(lpdf, stan_log_lik_dpar_usc_logit("hu", bterms))
  sdist(lpdf, p$mu, p$hu)
}

stan_log_lik_hurdle_negbinomial <- function(bterms, resp = "", mix = "", ...) {
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix)
  lpdf <- stan_log_lik_simple_lpdf("hurdle_neg_binomial", "log", bterms)
  lpdf <- paste0(lpdf, stan_log_lik_dpar_usc_logit("hu", bterms))
  sdist(lpdf, p$mu, p$shape, p$hu)
}

stan_log_lik_hurdle_gamma <- function(bterms, resp = "", mix = "", ...) {
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix)
  usc_logit <- stan_log_lik_dpar_usc_logit("hu", bterms)
  lpdf <- paste0("hurdle_gamma", usc_logit)
  sdist(lpdf, p$shape, p$mu, p$hu)
}

stan_log_lik_hurdle_lognormal <- function(bterms, resp = "", mix = "", ...) {
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix)
  usc_logit <- stan_log_lik_dpar_usc_logit("hu", bterms)
  lpdf <- paste0("hurdle_lognormal", usc_logit)
  sdist(lpdf, p$mu, p$sigma, p$hu)
}

stan_log_lik_zero_inflated_poisson <- function(bterms, resp = "", mix = "", 
                                               ...) {
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix)
  lpdf <- stan_log_lik_simple_lpdf("zero_inflated_poisson", "log", bterms)
  lpdf <- paste0(lpdf, stan_log_lik_dpar_usc_logit("zi", bterms))
  sdist(lpdf, p$mu, p$zi)
}

stan_log_lik_zero_inflated_negbinomial <- function(bterms, resp = "", mix = "",
                                                   ...) {
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix)
  lpdf <- stan_log_lik_simple_lpdf("zero_inflated_neg_binomial", "log", bterms)
  lpdf <- paste0(lpdf, stan_log_lik_dpar_usc_logit("zi", bterms))
  sdist(lpdf, p$mu, p$shape, p$zi)
}

stan_log_lik_zero_inflated_binomial <- function(bterms, resp = "", mix = "",
                                                threads = 1, ...) {
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix)
  n <- stan_nn(threads)
  p$trials <- paste0("trials", resp, n)
  lpdf <- "zero_inflated_binomial"
  lpdf <- stan_log_lik_simple_lpdf(lpdf, "logit", bterms, sep = "_b")
  lpdf <- paste0(lpdf, stan_log_lik_dpar_usc_logit("zi", bterms))
  sdist(lpdf, p$trials, p$mu, p$zi)
}

stan_log_lik_zero_inflated_beta <- function(bterms, resp = "", mix = "", ...) {
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix)
  usc_logit <- stan_log_lik_dpar_usc_logit("zi", bterms)
  lpdf <- paste0("zero_inflated_beta", usc_logit)
  sdist(lpdf, p$mu, p$phi, p$zi)
}

stan_log_lik_zero_one_inflated_beta <- function(bterms, resp = "", mix = "", 
                                                ...) {
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix)
  sdist("zero_one_inflated_beta", p$mu, p$phi, p$zoi, p$coi)
}

stan_log_lik_zero_inflated_asym_laplace <- function(bterms, resp = "", mix = "",
                                                    ...) {
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix)
  usc_logit <- stan_log_lik_dpar_usc_logit("zi", bterms)
  lpdf <- paste0("zero_inflated_asym_laplace", usc_logit)
  sdist(lpdf, p$mu, p$sigma, p$quantile, p$zi)
}

stan_log_lik_custom <- function(bterms, resp = "", mix = "", ...) {
  # TODO: support reduce_sum
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix)
  family <- bterms$family
  dpars <- paste0(family$dpars, mix)
  if (is_ordinal(family)) {
    prefix <- paste0(resp, if (nzchar(mix)) paste0("_mu", mix))
    p$thres <- paste0("Intercept", prefix)
  }
  # insert the response name into the 'vars' strings
  # addition terms contain the response in their variable name
  var_names <- sub("\\[.+$", "", family$vars)
  var_indices <- sub("^.+(?=\\[)", "", family$vars, perl = TRUE)
  is_var_adterms <- var_names %in% c("se", "trials", "dec") |
    grepl("^((vint)|(vreal))[[:digit:]]+$", var_names)
  var_resps <- ifelse(is_var_adterms, resp, "")
  vars <- paste0(var_names, var_resps, var_indices)
  sdist(family$name, p[dpars], p$thres, vars)
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
      is.formula(mu$ac)) {
    return(FALSE)
  }
  # supported families and link functions
  # TODO: support ordered_logistic and categorical_logit primitives
  glm_links <- list(
    gaussian = "identity", bernoulli = "logit",
    poisson = "log", negbinomial = "log"
  )
  if (!isTRUE(glm_links[[mu$family$family]] == mu$family$link)) {
    return(FALSE)
  }
  length(all_terms(mu$fe)) > 0 && !is_sparse(mu$fe)
}

# standard arguments for primitive Stan GLM functions
# @param bterms a btl object
# @param resp optional name of the response variable
# @return a named list of Stan code snippets
args_glm_primitive <- function(bterms, resp = "", threads = 1) {
  stopifnot(is.btl(bterms))
  decomp <- get_decomp(bterms$fe)
  center_X <- stan_center_X(bterms)
  slice <- stan_slice(threads)
  sfx_X <- sfx_b <- ""
  if (decomp == "QR") {
    sfx_X <- sfx_b <- "Q"
  } else if (center_X) {
    sfx_X <- "c"
  }
  x <- glue("X{sfx_X}{resp}{slice}")
  beta <- glue("b{sfx_b}{resp}")
  if (has_special_terms(bterms)) {
    # the intercept vector will contain all the remaining terms
    alpha <- glue("mu{resp}")
  } else {
    if (center_X) {
      alpha <- glue("Intercept{resp}")
    } else {
      alpha <- "0"
    } 
  }
  nlist(x, alpha, beta)
}

# use the ordered_logistic built-in functions
use_ordered_logistic <- function(bterms) {
  stopifnot(is.brmsterms(bterms))
  isTRUE(bterms$family$family == "cumulative") &&
    isTRUE(bterms$family$link == "logit") && 
    isTRUE(bterms$fdpars$disc$value == 1) &&
    !has_cs(bterms)
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
