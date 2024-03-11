# unless otherwise specified, functions return a single character
# string defining the likelihood of the model in Stan language

stan_log_lik <- function(x, ...) {
  UseMethod("stan_log_lik")
}

# Stan code for the model likelihood
# @param bterms object of class brmsterms
# @param data data passed by the user
# @param mix optional mixture component ID
# @param ptheta are mixing proportions predicted?
#' @export
stan_log_lik.family <- function(x, bterms, data, threads, normalize,
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
  args <- nlist(ll, bterms, data, resp, threads, normalize, mix, ptheta)
  if (nzchar(mix)) {
    out <- do_call(stan_log_lik_mix, args)
  } else if (is.formula(bterms$adforms$cens)) {
    out <- do_call(stan_log_lik_cens, args)
  } else if (is.formula(bterms$adforms$weights)) {
    out <- do_call(stan_log_lik_weights, args)
  } else {
    out <- do_call(stan_log_lik_general, args)
  }
  if (grepl(stan_nn_regex(), out) && !nzchar(mix)) {
    # loop over likelihood if it cannot be vectorized
    out <- paste0(
      "  for (n in 1:N", resp, ") {\n",
      stan_nn_def(threads),
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
    "    array[{length(ll)}] real ps;\n"
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
stan_log_lik_general <- function(ll, bterms, data, threads, normalize, resp = "", ...) {
  stopifnot(is.sdist(ll))
  require_n <- grepl(stan_nn_regex(), ll$args)
  n <- str_if(require_n, stan_nn(threads), stan_slice(threads))
  lpdf <- stan_log_lik_lpdf_name(bterms, normalize, dist = ll$dist)
  Y <- stan_log_lik_Y_name(bterms)
  tr <- stan_log_lik_trunc(ll, bterms, data, resp = resp, threads = threads)
  glue("{tp()}{ll$dist}_{lpdf}({Y}{resp}{n}{ll$shift} | {ll$args}){tr};\n")
}

# censored likelihood in Stan language
stan_log_lik_cens <- function(ll, bterms, data, threads, normalize, resp = "", ...) {
  stopifnot(is.sdist(ll))
  s <- wsp(nsp = 4)
  cens <- eval_rhs(bterms$adforms$cens)
  lpdf <- stan_log_lik_lpdf_name(bterms, normalize, dist = ll$dist)
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
stan_log_lik_weights <- function(ll, bterms, data, threads, normalize, resp = "", ...) {
  stopifnot(is.sdist(ll))
  tr <- stan_log_lik_trunc(ll, bterms, data, resp = resp, threads = threads)
  lpdf <- stan_log_lik_lpdf_name(bterms, normalize, dist = ll$dist)
  Y <- stan_log_lik_Y_name(bterms)
  n <- stan_nn(threads)
  glue(
    "{tp()}weights{resp}{n} * ({ll$dist}_{lpdf}",
    "({Y}{resp}{n}{ll$shift} | {ll$args}){tr});\n"
  )
}

# likelihood of a single mixture component
stan_log_lik_mix <- function(ll, bterms, data, mix, ptheta, threads,
                             normalize, resp = "", ...) {
  stopifnot(is.sdist(ll))
  theta <- str_if(ptheta,
    glue("theta{mix}{resp}[n]"),
    glue("log(theta{mix}{resp})")
  )
  tr <- stan_log_lik_trunc(ll, bterms, data, resp = resp, threads = threads)
  lpdf <- stan_log_lik_lpdf_name(bterms, normalize, dist = ll$dist)
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
      "ps[{mix}] = {theta} + ",
      "{ll$dist}_{lpdf}({Y}{resp}{n}{ll$shift} | {ll$args}){tr};\n"
    )
  }
  out
}

# truncated part of the likelihood
# @param short use the T[, ] syntax?
stan_log_lik_trunc <- function(ll, bterms, data, threads, resp = "",
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
      out <- glue("{ms}{ll$dist}_lccdf({lb}{ll$shift} | {ll$args})")
    } else if (!any(bounds$lb > -Inf) && any(bounds$ub < Inf)) {
      out <- glue("{ms}{ll$dist}_lcdf({ub}{ll$shift} | {ll$args})")
    } else if (any(bounds$lb > -Inf) && any(bounds$ub < Inf)) {
      trr <- glue("{ll$dist}_lcdf({ub}{ll$shift} | {ll$args})")
      trl <- glue("{ll$dist}_lcdf({lb}{ll$shift} | {ll$args})")
      out <- glue("{ms}log_diff_exp({trr}, {trl})")
    }
  }
  out
}

stan_log_lik_lpdf_name <- function(bterms, normalize, dist = NULL) {
  if (!is.null(dist) && !normalize) {
    # some Stan lpdfs or lpmfs only exist as normalized versions
    always_normalized <- always_normalized(bterms)
    if (length(always_normalized)) {
      always_normalized <- paste0(escape_all(always_normalized), "$")
      normalize <- any(ulapply(always_normalized, grepl, x = dist))
    }
  }
  if (normalize) {
    out <- ifelse(use_int(bterms$family), "lpmf", "lpdf")
  } else {
    out <- ifelse(use_int(bterms$family), "lupmf", "lupdf")
  }
  out
}

stan_log_lik_Y_name <- function(bterms) {
  ifelse(is.formula(bterms$adforms$mi), "Yl", "Y")
}

# prepare names of distributional parameters
# @param reqn will the likelihood be wrapped in a loop over n?
# @param dpars optional names of distributional parameters to be prepared
#   if not specified will prepare all distributional parameters
stan_log_lik_dpars <- function(bterms, reqn, resp = "", mix = "", dpars = NULL,
                               type = NULL) {
  if (is.null(dpars)) {
    dpars <- paste0(valid_dpars(bterms, type = type), mix)
  }
  pred_dpars <- names(bterms$dpars)
  if (is_equal(type, "multi")) {
    pred_dpars <- unique(dpar_class(pred_dpars, bterms))
  }
  is_pred <- dpars %in% pred_dpars
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
stan_log_lik_add_se <- function(sigma, bterms, reqn, resp = "",
                                threads = NULL) {
  if (!is.formula(bterms$adforms$se)) {
    return(sigma)
  }
  nse <- str_if(reqn, stan_nn(threads), stan_slice(threads))
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
                                             log = FALSE, transform = NULL,
                                             threads = NULL) {
  dpar_transform <- dpar
  if (!is.null(transform)) {
    dpar_transform <- glue("{transform}({dpar})")
  }
  if (!is.formula(bterms$adforms$rate)) {
    return(dpar_transform)
  }
  ndenom <- str_if(reqn, stan_nn(threads), stan_slice(threads))
  denom <- glue("denom{resp}{ndenom}")
  cens_or_trunc <- stan_log_lik_adj(bterms, c("cens", "trunc"))
  if (log && bterms$family$link == "log" && !cens_or_trunc) {
    denom <- glue("log_{denom}")
    operator <- "+"
  } else {
    # dpar without resp name or index
    dpar_clean <- sub("(_|\\[).*", "", dpar)
    is_pred <- dpar_clean %in% c("mu", names(bterms$dpars))
    operator <- str_if(reqn || !is_pred, "*", ".*")
  }
  glue("{dpar_transform} {operator} {denom}")
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
stan_log_lik_gaussian <- function(bterms, resp = "", mix = "", threads = NULL,
                                  ...) {
  if (use_glm_primitive(bterms)) {
    p <- args_glm_primitive(bterms$dpars$mu, resp = resp, threads = threads)
    p$sigma <- paste0("sigma", resp)
    out <- sdist("normal_id_glm", p$x, p$alpha, p$beta, p$sigma)
  } else {
    reqn <- stan_log_lik_adj(bterms) || nzchar(mix)
    p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
    p$sigma <- stan_log_lik_add_se(p$sigma, bterms, reqn, resp, threads)
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
  has_se <- is.formula(bterms$adforms$se)
  flex <- has_ac_class(tidy_acef(bterms), "unstr")
  p <- stan_log_lik_dpars(bterms, FALSE, resp, mix)
  v <- c("Lcortime", "nobs_tg", "begin_tg", "end_tg")
  if (has_se) {
    c(v) <- "se2"
  }
  if (flex) {
    c(v) <- "Jtime_tg"
  }
  p[v] <- as.list(paste0(v, resp))
  sfx <- str_if("sigma" %in% names(bterms$dpars), "het", "hom")
  sfx <- str_if(has_se, paste0(sfx, "_se"), sfx)
  sfx <- str_if(flex, paste0(sfx, "_flex"), sfx)
  sdist(glue("normal_time_{sfx}"),
    p$mu, p$sigma, p$se2, p$Lcortime,
    p$nobs_tg, p$begin_tg, p$end_tg, p$Jtime_tg
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

stan_log_lik_gaussian_lagsar <- function(bterms, resp = "", mix = "",
                                         threads = NULL, ...) {
  p <- stan_log_lik_dpars(bterms, FALSE, resp, mix)
  p$sigma <- stan_log_lik_add_se(p$sigma, bterms, FALSE, resp, threads)
  v <- c("lagsar", "Msar", "eigenMsar")
  p[v] <- as.list(paste0(v, resp))
  sdist("normal_lagsar", p$mu, p$sigma, p$lagsar, p$Msar, p$eigenMsar)
}

stan_log_lik_gaussian_errorsar <- function(bterms, resp = "", mix = "",
                                           threads = NULL, ...) {
  p <- stan_log_lik_dpars(bterms, FALSE, resp, mix)
  p$sigma <- stan_log_lik_add_se(p$sigma, bterms, FALSE, resp, threads)
  v <- c("errorsar", "Msar", "eigenMsar")
  p[v] <- as.list(paste0(v, resp))
  sdist("normal_errorsar", p$mu, p$sigma, p$errorsar, p$Msar, p$eigenMsar)
}

stan_log_lik_student <- function(bterms, resp = "", mix = "",
                                 threads = NULL, ...) {
  reqn <- stan_log_lik_adj(bterms) || nzchar(mix)
  p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
  p$sigma <- stan_log_lik_add_se(p$sigma, bterms, reqn, resp, threads)
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
  has_se <- is.formula(bterms$adforms$se)
  flex <- has_ac_class(tidy_acef(bterms), "unstr")
  p <- stan_log_lik_dpars(bterms, FALSE, resp, mix)
  v <- c("Lcortime", "nobs_tg", "begin_tg", "end_tg")
  if (has_se) {
    c(v) <- "se2"
  }
  if (flex) {
    c(v) <- "Jtime_tg"
  }
  p[v] <- as.list(paste0(v, resp))
  sfx <- str_if("sigma" %in% names(bterms$dpars), "het", "hom")
  sfx <- str_if(has_se, paste0(sfx, "_se"), sfx)
  sfx <- str_if(flex, paste0(sfx, "_flex"), sfx)
  sdist(glue("student_t_time_{sfx}"),
    p$nu, p$mu, p$sigma, p$se2, p$Lcortime,
    p$nobs_tg, p$begin_tg, p$end_tg, p$Jtime_tg
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

stan_log_lik_student_lagsar <- function(bterms, resp = "", mix = "",
                                        threads = NULL, ...) {
  p <- stan_log_lik_dpars(bterms, FALSE, resp, mix)
  p$sigma <- stan_log_lik_add_se(p$sigma, bterms, FALSE, resp, threads)
  v <- c("lagsar", "Msar", "eigenMsar")
  p[v] <- as.list(paste0(v, resp))
  sdist("student_t_lagsar", p$nu, p$mu, p$sigma,
        p$lagsar, p$Msar, p$eigenMsar)
}

stan_log_lik_student_errorsar <- function(bterms, resp = "", mix = "",
                                          threads = NULL, ...) {
  p <- stan_log_lik_dpars(bterms, FALSE, resp, mix)
  p$sigma <- stan_log_lik_add_se(p$sigma, bterms, FALSE, resp, threads)
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

stan_log_lik_skew_normal <- function(bterms, resp = "", mix = "",
                                     threads = NULL, ...) {
  reqn <- stan_log_lik_adj(bterms) || nzchar(mix)
  p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
  p$sigma <- stan_log_lik_add_se(p$sigma, bterms, reqn, resp, threads)
  # required because of CP parameterization of mu and sigma
  nomega <- any(grepl(stan_nn_regex(), c(p$sigma, p$alpha)))
  nomega <- str_if(reqn && nomega, "[n]")
  p$omega <- paste0("omega", mix, resp, nomega)
  sdist("skew_normal", p$mu, p$omega, p$alpha)
}

stan_log_lik_poisson <- function(bterms, resp = "", mix = "", threads = NULL,
                                 ...) {
  if (use_glm_primitive(bterms)) {
    p <- args_glm_primitive(bterms$dpars$mu, resp = resp, threads = threads)
    out <- sdist("poisson_log_glm", p$x, p$alpha, p$beta)
  } else {
    reqn <- stan_log_lik_adj(bterms) || nzchar(mix)
    p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
    lpdf <- stan_log_lik_simple_lpdf("poisson", "log", bterms)
    p$mu <- stan_log_lik_multiply_rate_denom(
      p$mu, bterms, reqn, resp, log = TRUE, threads = threads
    )
    out <- sdist(lpdf, p$mu)
  }
  out
}

stan_log_lik_negbinomial <- function(bterms, resp = "", mix = "", threads = NULL,
                                     ...) {
  if (use_glm_primitive(bterms)) {
    p <- args_glm_primitive(bterms$dpars$mu, resp = resp, threads = threads)
    p$shape <- paste0("shape", resp)
    out <- sdist("neg_binomial_2_log_glm", p$x, p$alpha, p$beta, p$shape)
  } else {
    reqn <- stan_log_lik_adj(bterms) || nzchar(mix)
    p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
    p$mu <- stan_log_lik_multiply_rate_denom(
      p$mu, bterms, reqn, resp, log = TRUE, threads = threads
    )
    p$shape <- stan_log_lik_multiply_rate_denom(
      p$shape, bterms, reqn, resp, threads = threads
    )
    lpdf <- stan_log_lik_simple_lpdf("neg_binomial_2", "log", bterms)
    out <- sdist(lpdf, p$mu, p$shape)
  }
  out
}

stan_log_lik_negbinomial2 <- function(bterms, resp = "", mix = "", threads = NULL,
                                      ...) {
  if (use_glm_primitive(bterms)) {
    p <- args_glm_primitive(bterms$dpars$mu, resp = resp, threads = threads)
    p$sigma <- paste0("sigma", resp)
    p$shape <- paste0("inv(", p$sigma, ")")
    out <- sdist("neg_binomial_2_log_glm", p$x, p$alpha, p$beta, p$shape)
  } else {
    reqn <- stan_log_lik_adj(bterms) || nzchar(mix)
    p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
    p$mu <- stan_log_lik_multiply_rate_denom(
      p$mu, bterms, reqn, resp, log = TRUE, threads = threads
    )
    p$shape <- stan_log_lik_multiply_rate_denom(
      p$sigma, bterms, reqn, resp, transform = "inv", threads = threads
    )
    lpdf <- stan_log_lik_simple_lpdf("neg_binomial_2", "log", bterms)
    out <- sdist(lpdf, p$mu, p$shape)
  }
  out
}

stan_log_lik_geometric <- function(bterms, resp = "", mix = "", threads = NULL,
                                   ...) {
  if (use_glm_primitive(bterms)) {
    p <- args_glm_primitive(bterms$dpars$mu, resp = resp, threads = threads)
    p$shape <- "1"
    out <- sdist("neg_binomial_2_log_glm", p$x, p$alpha, p$beta, p$shape)
  } else {
    reqn <- stan_log_lik_adj(bterms) || nzchar(mix)
    p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
    p$shape <- "1"
    p$mu <- stan_log_lik_multiply_rate_denom(
      p$mu, bterms, reqn, resp, log = TRUE, threads = threads
    )
    p$shape <- stan_log_lik_multiply_rate_denom(
      p$shape, bterms, reqn, resp, threads = threads
    )
    lpdf <- stan_log_lik_simple_lpdf("neg_binomial_2", "log", bterms)
    out <- sdist(lpdf, p$mu, p$shape)
  }
}

stan_log_lik_binomial <- function(bterms, resp = "", mix = "", threads = NULL,
                                  ...) {
  reqn <- stan_log_lik_adj(bterms) || nzchar(mix)
  p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
  slice <- str_if(reqn, stan_nn(threads), stan_slice(threads))
  p$trials <- paste0("trials", resp, slice)
  lpdf <- stan_log_lik_simple_lpdf("binomial", "logit", bterms)
  sdist(lpdf, p$trials, p$mu)
}

stan_log_lik_beta_binomial <- function(bterms, resp = "", mix = "",
                                       threads = NULL, ...) {
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix)
  n <- stan_nn(threads)
  sdist(
    "beta_binomial",
    paste0("trials", resp, n),
    paste0(p$mu, " * ", p$phi),
    paste0("(1 - ", p$mu, ") * ", p$phi)
  )
}

stan_log_lik_bernoulli <- function(bterms, resp = "", mix = "", threads = NULL,
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
  reqn <- stan_log_lik_adj(bterms) || nzchar(mix) ||
    paste0("shape", mix) %in% names(bterms$dpars)
  p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
  # Stan uses shape-rate parameterization with rate = shape / mean
  div_op <- str_if(reqn, " / ", " ./ ")
  sdist("gamma", p$shape, paste0(p$shape, div_op, p$mu))
}

stan_log_lik_exponential <- function(bterms, resp = "", mix = "", ...) {
  reqn <- stan_log_lik_adj(bterms) || nzchar(mix)
  p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
  # Stan uses rate parameterization with rate = 1 / mean
  sdist("exponential", paste0("inv(", p$mu, ")"))
}

stan_log_lik_weibull <- function(bterms, resp = "", mix = "", ...) {
  reqn <- stan_log_lik_adj(bterms) || nzchar(mix)
  p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
  # Stan uses shape-scale parameterization for weibull
  need_dot_div <- !reqn && paste0("shape", mix) %in% names(bterms$dpars)
  div_op <- str_if(need_dot_div, " ./ ", " / ")
  p$scale <- paste0(p$mu, div_op, "tgamma(1 + 1", div_op, p$shape, ")")
  sdist("weibull", p$shape, p$scale)
}

stan_log_lik_frechet <- function(bterms, resp = "", mix = "", ...) {
  reqn <- stan_log_lik_adj(bterms) || nzchar(mix)
  p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
  # Stan uses shape-scale parameterization for frechet
  need_dot_div <- !reqn && paste0("nu", mix) %in% names(bterms$dpars)
  div_op <- str_if(need_dot_div, " ./ ", " / ")
  p$scale <- paste0(p$mu, div_op, "tgamma(1 - 1", div_op, p$nu, ")")
  sdist("frechet", p$nu, p$scale)
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
  n <- str_if(reqn, "[n]")
  sdist("inv_gaussian", p$mu, p$shape)
}

stan_log_lik_wiener <- function(bterms, resp = "", mix = "", threads = NULL,
                                ...) {
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix)
  n <- stan_nn(threads)
  p$dec <- paste0("dec", resp, n)
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
  sdist("von_mises2", p$mu, p$kappa)
}

stan_log_lik_cox <- function(bterms, resp = "", mix = "", threads = NULL,
                             ...) {
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix)
  p$bhaz <- paste0("bhaz", resp, "[n]")
  p$cbhaz <- paste0("cbhaz", resp, "[n]")
  lpdf <- "cox"
  if (bterms$family$link == "log") {
    str_add(lpdf) <- "_log"
  }
  sdist(lpdf, p$mu, p$bhaz, p$cbhaz)
}

stan_log_lik_cumulative <- function(bterms, resp = "", mix = "",
                                    threads = NULL, ...) {
  if (use_glm_primitive(bterms, allow_special_terms = FALSE)) {
    p <- args_glm_primitive(bterms$dpars$mu, resp = resp, threads = threads)
    out <- sdist("ordered_logistic_glm", p$x, p$beta, p$alpha)
  } else {
    out <- stan_log_lik_ordinal(bterms, resp, mix, threads, ...)
  }
  out
}

stan_log_lik_sratio <- function(bterms, resp = "", mix = "",
                                threads = NULL, ...) {
  stan_log_lik_ordinal(bterms, resp, mix, threads, ...)
}

stan_log_lik_cratio <- function(bterms, resp = "", mix = "",
                                threads = NULL, ...) {
  stan_log_lik_ordinal(bterms, resp, mix, threads, ...)
}

stan_log_lik_acat <- function(bterms, resp = "", mix = "",
                              threads = NULL, ...) {
  stan_log_lik_ordinal(bterms, resp, mix, threads, ...)
}

stan_log_lik_categorical <- function(bterms, resp = "", mix = "",
                                     threads = NULL, ...) {
  stopifnot(bterms$family$link == "logit")
  stopifnot(!isTRUE(nzchar(mix)))  # mixture models are not allowed
  # if (use_glm_primitive_categorical(bterms)) {
  #   # TODO: support categorical_logit_glm
  # }
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix, dpars = "mu", type = "multi")
  sdist("categorical_logit", p$mu)
}

stan_log_lik_multinomial <- function(bterms, resp = "", mix = "", ...) {
  stopifnot(bterms$family$link == "logit")
  stopifnot(!isTRUE(nzchar(mix)))  # mixture models are not allowed
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix, dpars = "mu", type = "multi")
  sdist("multinomial_logit2", p$mu)
}

stan_log_lik_dirichlet <- function(bterms, resp = "", mix = "", ...) {
  stopifnot(bterms$family$link == "logit")
  stopifnot(!isTRUE(nzchar(mix)))  # mixture models are not allowed
  mu <- stan_log_lik_dpars(bterms, TRUE, resp, mix, dpars = "mu", type = "multi")$mu
  reqn <- glue("phi{mix}") %in% names(bterms$dpars)
  phi <- stan_log_lik_dpars(bterms, reqn, resp, mix, dpars = "phi")$phi
  sdist("dirichlet_logit", mu, phi)
}

stan_log_lik_dirichlet2 <- function(bterms, resp = "", mix = "", ...) {
  stopifnot(!isTRUE(nzchar(mix)))  # mixture models are not allowed
  mu <- stan_log_lik_dpars(bterms, TRUE, resp, mix, dpars = "mu", type = "multi")$mu
  sdist("dirichlet", mu)
}

stan_log_lik_logistic_normal <- function(bterms, resp = "", mix = "", ...) {
  stopifnot(bterms$family$link == "identity")
  stopifnot(!isTRUE(nzchar(mix)))  # mixture models are not allowed
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix, type = "multi")
  p$Llncor <- glue("Llncor{mix}{resp}")
  p$refcat <- get_refcat(bterms$family, int = TRUE)
  sdist("logistic_normal_cholesky_cor", p$mu, p$sigma, p$Llncor, p$refcat)
}

stan_log_lik_ordinal <- function(bterms, resp = "", mix = "",
                                 threads = NULL, ...) {
  prefix <- paste0(str_if(nzchar(mix), paste0("_mu", mix)), resp)
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix)
  if (use_ordered_builtin(bterms, "logit")) {
    lpdf <- "ordered_logistic"
    p[grepl("^disc", names(p))] <- NULL
  } else if (use_ordered_builtin(bterms, "probit")) {
    lpdf <- "ordered_probit"
    p[grepl("^disc", names(p))] <- NULL
  } else {
    lpdf <- paste0(bterms$family$family, "_", bterms$family$link)
  }
  if (has_thres_groups(bterms)) {
    str_add(lpdf) <- "_merged"
    n <- stan_nn(threads)
    p$Jthres <- paste0("Jthres", resp, n)
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
    str_add(p$thres) <- paste0(" - transpose(mucs", prefix, "[n])")
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
  # Stan uses shape-rate parameterization for gamma with rate = shape / mean
  sdist(lpdf, p$shape, paste0(p$shape, " / ", p$mu), p$hu)
}

stan_log_lik_hurdle_lognormal <- function(bterms, resp = "", mix = "", ...) {
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix)
  usc_logit <- stan_log_lik_dpar_usc_logit("hu", bterms)
  lpdf <- paste0("hurdle_lognormal", usc_logit)
  sdist(lpdf, p$mu, p$sigma, p$hu)
}

stan_log_lik_hurdle_cumulative <- function(bterms, resp = "", mix = "",
                                           threads = NULL, ...) {
  prefix <- paste0(str_if(nzchar(mix), paste0("_mu", mix)), resp)
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix)
  if (use_ordered_builtin(bterms, "logit")) {
    lpdf <- "hurdle_cumulative_ordered_logistic"
  } else if (use_ordered_builtin(bterms, "probit")) {
    lpdf <- "hurdle_cumulative_ordered_probit"
  } else {
    lpdf <- paste0(bterms$family$family, "_", bterms$family$link)
  }
  if (has_thres_groups(bterms)) {
    str_add(lpdf) <- "_merged"
    n <- stan_nn(threads)
    p$Jthres <- paste0("Jthres", resp, n)
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
    str_add(p$thres) <- paste0(" - transpose(mucs", prefix, "[n])")
  }
  sdist(lpdf, p$mu, p$hu, p$disc, p$thres, p$Jthres)
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
                                                threads = NULL, ...) {
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix)
  n <- stan_nn(threads)
  p$trials <- paste0("trials", resp, n)
  lpdf <- "zero_inflated_binomial"
  lpdf <- stan_log_lik_simple_lpdf(lpdf, "logit", bterms, sep = "_b")
  lpdf <- paste0(lpdf, stan_log_lik_dpar_usc_logit("zi", bterms))
  sdist(lpdf, p$trials, p$mu, p$zi)
}

stan_log_lik_zero_inflated_beta_binomial <- function(bterms, resp = "",
                                                     mix = "", threads = NULL,
                                                     ...) {
  p <- stan_log_lik_dpars(bterms, TRUE, resp, mix)
  n <- stan_nn(threads)
  p$trials <- paste0("trials", resp, n)
  lpdf <- "zero_inflated_beta_binomial"
  lpdf <- paste0(lpdf, stan_log_lik_dpar_usc_logit("zi", bterms))
  sdist(lpdf, p$trials, p$mu, p$phi, p$zi)
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

stan_log_lik_custom <- function(bterms, resp = "", mix = "", threads = NULL, ...) {
  family <- bterms$family
  no_loop <- isFALSE(family$loop)
  if (no_loop && (stan_log_lik_adj(bterms) || nzchar(mix))) {
    stop2("This model requires evaluating the custom ",
          "likelihood as a loop over observations.")
  }
  reqn <- !no_loop
  p <- stan_log_lik_dpars(bterms, reqn, resp, mix)
  dpars <- paste0(family$dpars, mix)
  if (is_ordinal(family)) {
    prefix <- paste0(resp, if (nzchar(mix)) paste0("_mu", mix))
    p$thres <- paste0("Intercept", prefix)
  }
  # insert the response name into the 'vars' strings
  # addition terms contain the response in their variable name
  n <- stan_nn(threads)
  var_names <- sub("\\[.+$", "", family$vars)
  var_indices <- get_matches("\\[.+$", family$vars, first = TRUE)
  has_n_index <- var_indices %in% "[n]"
  if (no_loop && any(has_n_index)) {
    stop2("Invalid use of index '[n]' in an unlooped custom likelihood.")
  }
  var_indices <- ifelse(has_n_index, n, var_indices)
  is_var_adterms <- var_names %in% c("se", "trials", "dec") |
    grepl("^((vint)|(vreal))[[:digit:]]+$", var_names)
  var_resps <- ifelse(is_var_adterms, resp, "")
  vars <- paste0(var_names, var_resps, var_indices)
  sdist(family$name, p[dpars], p$thres, vars)
}

# use Stan GLM primitive functions?
# @param bterms a brmsterms object
# @param allow_special_terms still use glm primitives if
#   random effects, splines, etc. are present?
# @return TRUE or FALSE
use_glm_primitive <- function(bterms, allow_special_terms = TRUE) {
  stopifnot(is.brmsterms(bterms))
  # the model can only have a single predicted parameter
  # and no additional residual or autocorrelation structure
  mu <- bterms$dpars[["mu"]]
  non_glm_adterms <- c("se", "weights", "thres", "cens", "trunc", "rate")
  if (!is.btl(mu) || length(bterms$dpars) > 1L ||
      isTRUE(bterms$rescor) || is.formula(mu$ac) ||
      any(names(bterms$adforms) %in% non_glm_adterms)) {
    return(FALSE)
  }
  # some primitives do not support special terms in the way
  # required by brms' Stan code generation
  if (!allow_special_terms && has_special_terms(mu)) {
    return(FALSE)
  }
  # supported families and link functions
  glm_links <- list(
    gaussian = "identity", bernoulli = "logit",
    poisson = "log", negbinomial = "log", negbinomial2 = "log",
    cumulative = "logit", categorical = "logit"
  )
  if (!isTRUE(glm_links[[mu$family$family]] == mu$family$link)) {
    return(FALSE)
  }
  length(all_terms(mu$fe)) > 0 && !is_sparse(mu$fe)
}

# use Stan categorical GLM primitive function?
# @param bterms a brmsterms object
# @param ... passed to use_glm_primitive
# @return TRUE or FALSE
use_glm_primitive_categorical <- function(bterms, ...) {
  # NOTE: this function is not yet in use; see stan_log_lik_categorical
  stopifnot(is.brmsterms(bterms))
  stopifnot(is_categorical(bterms))
  bterms_tmp <- bterms
  bterms_tmp$dpars <- list()
  # we know that all dpars in categorical models are mu parameters
  out <- rep(FALSE, length(bterms$dpars))
  for (i in seq_along(bterms$dpars)) {
    bterms_tmp$dpars$mu <- bterms$dpars[[i]]
    bterms_tmp$dpars$mu$family <- bterms$family
    out[i] <- use_glm_primitive(bterms_tmp, ...) &&
      # the design matrix of all mu parameters must match
      all.equal(bterms_tmp$dpars$mu$fe, bterms$dpars[[1]]$fe)
  }
  all(out)
}

# standard arguments for primitive Stan GLM functions
# @param bterms a btl object
# @param resp optional name of the response variable
# @return a named list of Stan code snippets
args_glm_primitive <- function(bterms, resp = "", threads = NULL) {
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
use_ordered_builtin <- function(bterms, link) {
  stopifnot(is.brmsterms(bterms))
  isTRUE(bterms$family$family %in% c("cumulative", "hurdle_cumulative")) &&
    isTRUE(bterms$family$link == link) &&
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
