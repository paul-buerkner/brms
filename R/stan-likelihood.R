# unless otherwise specified, functions return a single character
# string defining the likelihood of the model in Stan language

# Stan code for the log likelihood
stan_log_lik <- function(x, ...) {
  UseMethod("stan_log_lik")
}

#' @export
stan_log_lik.brmsterms <- function(x, ...) {
  if (is.mixfamily(x$family)) {
    out <- stan_log_lik_mixfamily(x, ...)
  } else {
    out <- stan_log_lik_family(x, ...)
  }
  out
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

# Stan code for the log likelihood of a regular family
stan_log_lik_family <- function(bterms, threads, ...) {
  stopifnot(is.brmsterms(bterms))
  # prepare family part of the likelihood
  log_lik_args <- nlist(bterms, threads, ...)
  log_lik_fun <- prepare_family(bterms)$fun
  log_lik_fun <- paste0("stan_log_lik_", log_lik_fun)
  ll <- do_call(log_lik_fun, log_lik_args)
  # incorporate other parts into the likelihood
  args <- nlist(ll, bterms, threads, ...)
  mix <- get_mix_id(bterms)
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
    resp <- usc(bterms$resp)
    out <- paste0(
      "  for (n in 1:N", resp, ") {\n",
      stan_nn_def(threads),
      "  ", out,
      "  }\n"
    )
  }
  out
}

# Stan code for the log likelihood of a mixture family
stan_log_lik_mixfamily <- function(bterms, threads, ...) {
  stopifnot(is.brmsterms(bterms), is.mixfamily(bterms$family))
  dp_ids <- dpar_id(names(bterms$dpars))
  fdp_ids <- dpar_id(names(bterms$fdpars))
  pred_mix_prob <- any(dpar_class(names(bterms$dpars)) %in% "theta")
  ll <- rep(NA, length(bterms$family$mix))
  for (i in seq_along(ll)) {
    sbterms <- bterms
    sbterms$family <- sbterms$family$mix[[i]]
    sbterms$dpars <- sbterms$dpars[dp_ids == i]
    sbterms$fdpars <- sbterms$fdpars[fdp_ids == i]
    ll[i] <- stan_log_lik_family(
      sbterms, pred_mix_prob = pred_mix_prob, threads = threads, ...
    )
  }
  resp <- usc(bterms$resp)
  n <- stan_nn(threads)
  has_weights <- has_ad_terms(bterms, "weights")
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

# default likelihood in Stan language
stan_log_lik_general <- function(ll, bterms, threads, normalize, ...) {
  stopifnot(is.sdist(ll))
  require_n <- grepl(stan_nn_regex(), ll$args)
  n <- str_if(require_n, stan_nn(threads), stan_slice(threads))
  lpdf <- stan_log_lik_lpdf_name(bterms, normalize, dist = ll$dist)
  Y <- stan_log_lik_Y_name(bterms)
  resp <- usc(bterms$resp)
  tr <- stan_log_lik_trunc(ll, bterms, threads = threads, ...)
  glue("{tp()}{ll$dist}_{lpdf}({Y}{resp}{n}{ll$shift} | {ll$args}){tr};\n")
}

# censored likelihood in Stan language
stan_log_lik_cens <- function(ll, bterms, threads, normalize, ...) {
  stopifnot(is.sdist(ll))
  cens <- eval_rhs(bterms$adforms$cens)
  lpdf <- stan_log_lik_lpdf_name(bterms, normalize, dist = ll$dist)
  Y <- stan_log_lik_Y_name(bterms)
  resp <- usc(bterms$resp)
  tp <- tp()
  has_weights <- has_ad_terms(bterms, "weights")
  has_trunc <- has_ad_terms(bterms, "trunc")
  has_interval_cens <- cens$vars$y2 != "NA"
  if (ll$vec && !(has_weights || has_trunc)) {
    # vectorized log-likelihood contributions
    types <- c("event", "rcens", "lcens", "icens")
    J <- args <- named_list(types)
    for (t in types) {
      Jt <- glue("J{t}{resp}[1:N{t}{resp}]")
      J[[t]] <- glue("[{Jt}]")
      if (use_threading(threads)) {
        Jtms <- glue("[add_int({Jt}, start - 1)]")
        args[[t]] <- rename(ll$args, c("[n]", "[nn]"), c(J[[t]], Jtms))
        J[[t]] <- Jtms
      } else {
        args[[t]] <- rename(ll$args, "[n]", J[[t]])
      }
    }
    out <- glue(
      "  // vectorized log-likelihood contributions of censored data\n",
      "{tp}{ll$dist}_{lpdf}(Y{resp}{J$event}{ll$shift} | {args$event});\n",
      "{tp}{ll$dist}_lccdf(Y{resp}{J$rcens}{ll$shift} | {args$rcens});\n",
      "{tp}{ll$dist}_lcdf(Y{resp}{J$lcens}{ll$shift} | {args$lcens});\n"
    )
    if (has_interval_cens) {
      rcens <- glue("rcens{resp}")
      str_add(out) <- glue(
        "{tp}log_diff_exp(\n",
        "    {ll$dist}_lcdf(rcens{resp}{J$icens}{ll$shift} | {args$icens}),\n",
        "    {ll$dist}_lcdf(Y{resp}{J$icens}{ll$shift} | {args$icens})\n",
        "  );\n"
      )
    }
  } else {
    # non-vectorized likelihood contributions
    n <- stan_nn(threads)
    w <- str_if(has_weights, glue("weights{resp}{n} * "))
    tr <- stan_log_lik_trunc(ll, bterms, threads = threads)
    out <- glue(
      "  // special treatment of censored data\n",
      "    if (cens{resp}{n} == 0) {{\n",
      "    {tp}{w}{ll$dist}_{lpdf}(Y{resp}{n}{ll$shift} | {ll$args}){tr};\n",
      "    }} else if (cens{resp}{n} == 1) {{\n",
      "    {tp}{w}{ll$dist}_lccdf(Y{resp}{n}{ll$shift} | {ll$args}){tr};\n",
      "    }} else if (cens{resp}{n} == -1) {{\n",
      "    {tp}{w}{ll$dist}_lcdf(Y{resp}{n}{ll$shift} | {ll$args}){tr};\n"
    )
    if (has_interval_cens) {
      str_add(out) <- glue(
        "    }} else if (cens{resp}{n} == 2) {{\n",
        "    {tp}{w}log_diff_exp(\n",
        "        {ll$dist}_lcdf(rcens{resp}{n}{ll$shift} | {ll$args}),\n",
        "        {ll$dist}_lcdf(Y{resp}{n}{ll$shift} | {ll$args})\n",
        "      ){tr};\n"
      )
    }
    str_add(out) <- glue("    }}\n")
  }
  out
}

# weighted likelihood in Stan language
stan_log_lik_weights <- function(ll, bterms, threads, normalize, ...) {
  stopifnot(is.sdist(ll))
  tr <- stan_log_lik_trunc(ll, bterms, threads = threads)
  lpdf <- stan_log_lik_lpdf_name(bterms, normalize, dist = ll$dist)
  Y <- stan_log_lik_Y_name(bterms)
  resp <- usc(bterms$resp)
  n <- stan_nn(threads)
  glue(
    "{tp()}weights{resp}{n} * ({ll$dist}_{lpdf}",
    "({Y}{resp}{n}{ll$shift} | {ll$args}){tr});\n"
  )
}

# likelihood of a single mixture component
# @param pred_mix_prob are mixing proportions predicted?
stan_log_lik_mix <- function(ll, bterms, pred_mix_prob, threads,
                             normalize, ...) {
  stopifnot(is.sdist(ll))
  resp <- usc(bterms$resp)
  mix <- get_mix_id(bterms)
  theta <- str_if(pred_mix_prob,
    glue("theta{mix}{resp}[n]"),
    glue("log(theta{mix}{resp})")
  )
  tr <- stan_log_lik_trunc(ll, bterms, threads = threads)
  lpdf <- stan_log_lik_lpdf_name(bterms, normalize, dist = ll$dist)
  Y <- stan_log_lik_Y_name(bterms)
  n <- stan_nn(threads)
  if (is.formula(bterms$adforms$cens)) {
    # mostly copied over from stan_log_lik_cens
    # no vectorized version available for mixture models
    cens <- eval_rhs(bterms$adforms$cens)
    out <- glue(
      "  // special treatment of censored data\n",
      "    if (cens{resp}{n} == 0) {{\n",
      "      ps[{mix}] = {theta} + ",
      "{ll$dist}_{lpdf}({Y}{resp}{n}{ll$shift} | {ll$args}){tr};\n",
      "    }} else if (cens{resp}{n} == 1) {{\n",
      "      ps[{mix}] = {theta} + ",
      "{ll$dist}_lccdf({Y}{resp}{n}{ll$shift} | {ll$args}){tr};\n",
      "    }} else if (cens{resp}{n} == -1) {{\n",
      "      ps[{mix}] = {theta} + ",
      "{ll$dist}_lcdf({Y}{resp}{n}{ll$shift} | {ll$args}){tr};\n"
    )
    has_interval_cens <- cens$vars$y2 != "NA"
    if (has_interval_cens) {
      str_add(out) <- glue(
        "    }} else if (cens{resp}{n} == 2) {{\n",
        "      ps[{mix}] = {theta} + log_diff_exp(\n",
        "        {ll$dist}_lcdf(rcens{resp}{n}{ll$shift} | {ll$args}),\n",
        "        {ll$dist}_lcdf({Y}{resp}{n}{ll$shift} | {ll$args})\n",
        "      ){tr};\n"
      )
    }
    str_add(out) <- glue("    }}\n")
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
stan_log_lik_trunc <- function(ll, bterms, threads, short = FALSE, ...) {
  stopifnot(is.sdist(ll))
  bounds <- bterms$frame$resp$bounds
  if (!any(bounds$lb > -Inf | bounds$ub < Inf)) {
    return("")
  }
  resp <- usc(bterms$resp)
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

# prepare Stan code for distributional parameters
# @param reqn will the likelihood be wrapped in a loop over n?
# @param dpars optional names of distributional parameters to be prepared
#   if not specified will prepare all distributional parameters
# @param type optional type of distribution parameters to be extract
#   see valid_dpars() for details
# @return a named list with elements containing the Stan code per parameter
stan_log_lik_dpars <- function(bterms, reqn = stan_log_lik_adj(bterms),
                               dpars = NULL, type = NULL, ...) {
  resp <- usc(bterms$resp)
  mix <- get_mix_id(bterms)
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

# stan code for log likelihood variables originating from addition terms
stan_log_lik_advars <- function(bterms, advars,
                                reqn = stan_log_lik_adj(bterms),
                                threads = NULL, ...) {
  slice <- str_if(reqn, stan_nn(threads), stan_slice(threads))
  out <- paste0(advars, usc(bterms$resp), slice)
  named_list(advars, out)
}

# adjust lpdf name if a more efficient version is available
# for a specific link. For instance 'poisson_log'
stan_log_lik_simple_lpdf <- function(lpdf, bterms, sep = "_") {
  stopifnot(is.brmsterms(bterms))
  if (stan_has_built_in_fun(bterms)) {
    lpdf <- paste0(lpdf, sep, bterms$family$link)
  }
  lpdf
}

# prepare _logit suffix for distributional parameters
# used in zero-inflated and hurdle models
stan_log_lik_dpar_usc_logit <- function(bterms, dpar) {
  stopifnot(is.brmsterms(bterms))
  stopifnot(dpar %in% c("zi", "hu"))
  has_cens_or_trunc <- has_ad_terms(bterms, c("cens", "trunc"))
  usc_logit <- isTRUE(bterms$dpars[[dpar]]$family$link == "logit")
  str_if(usc_logit && !has_cens_or_trunc, "_logit")
}

# add 'se' to 'sigma' within the Stan likelihood
stan_log_lik_add_se <- function(sigma, bterms, reqn = stan_log_lik_adj(bterms),
                                threads = NULL, ...) {
  if (!has_ad_terms(bterms, "se")) {
    return(sigma)
  }
  nse <- str_if(reqn, stan_nn(threads), stan_slice(threads))
  resp <- usc(bterms$resp)
  if (no_sigma(bterms)) {
    sigma <- glue("se{resp}{nse}")
  } else {
    sigma <- glue("sqrt(square({sigma}) + se2{resp}{nse})")
  }
  sigma
}

# multiply 'dpar' by the 'rate' denominator within the Stan likelihood
# @param log add the rate denominator on the log scale if sensible?
# @param req_dot_multiply Censoring may turn non-vectorized into vectorized
#   statements later on (see stan_log_lik_cens) which then makes the * operator
#   invalid and requires .* instead. Accordingly, req_dot_multiply should be
#   FALSE if [n] is required only because of censoring.
stan_log_lik_multiply_rate_denom <- function(
    dpar, bterms, reqn = stan_log_lik_adj(bterms),
    req_dot_multiply = stan_log_lik_adj(bterms, c("trunc", "weights")),
    log = FALSE, transform = NULL, threads = NULL, ...) {

  dpar_transform <- dpar
  if (!is.null(transform)) {
    dpar_transform <- glue("{transform}({dpar})")
  }
  if (!is.formula(bterms$adforms$rate)) {
    return(dpar_transform)
  }
  resp <- usc(bterms$resp)
  ndenom <- str_if(reqn, stan_nn(threads), stan_slice(threads))
  denom <- glue("denom{resp}{ndenom}")
  has_cens_or_trunc <- has_ad_terms(bterms, c("cens", "trunc"))
  if (log && bterms$family$link == "log" && !has_cens_or_trunc) {
    denom <- glue("log_{denom}")
    operator <- "+"
  } else {
    # dpar without resp name or index
    dpar_clean <- sub("(_|\\[).*", "", dpar)
    is_pred <- dpar_clean %in% c("mu", names(bterms$dpars))
    operator <- str_if(req_dot_multiply || !is_pred, "*", ".*")
  }
  glue("{dpar_transform} {operator} {denom}")
}

# check if the log-likelihood needs to be adjusted to a non-vectorized form
# either because of addition terms or mixture modeling
# @param terms vector of addition term names
# @return a single logical value
stan_log_lik_adj <- function(bterms, terms = c("weights", "cens", "trunc")) {
  stopifnot(is.brmsterms(bterms))
  terms <- match.arg(terms, several.ok = TRUE)
  mix <- get_mix_id(bterms)
  has_ad_terms(bterms, terms) || any(nzchar(mix))
}

# one function per family
stan_log_lik_gaussian <- function(bterms, ...) {
  if (use_glm_primitive(bterms)) {
    p <- args_glm_primitive(bterms$dpars$mu, ...)
    p$sigma <- paste0("sigma", usc(bterms$resp))
    out <- sdist("normal_id_glm", p$x, p$alpha, p$beta, p$sigma)
  } else {
    p <- stan_log_lik_dpars(bterms)
    p$sigma <- stan_log_lik_add_se(p$sigma, bterms, ...)
    out <- sdist("normal", p$mu, p$sigma)
  }
  out
}

stan_log_lik_gaussian_mv <- function(bterms,...) {
  reqn <- stan_log_lik_adj(bterms) || bterms$sigma_pred
  p <- list(Mu = paste0("Mu", str_if(reqn, "[n]")))
  p$LSigma <- paste0("LSigma", str_if(bterms$sigma_pred, "[n]"))
  sdist("multi_normal_cholesky", p$Mu, p$LSigma)
}

stan_log_lik_gaussian_time <- function(bterms, ...) {
  if (stan_log_lik_adj(bterms)) {
    stop2("Invalid addition arguments for this model.")
  }
  has_se <- is.formula(bterms$adforms$se)
  flex <- has_ac_class(bterms$frame$ac, "unstr")
  p <- stan_log_lik_dpars(bterms, reqn = FALSE)
  v <- c("Lcortime", "nobs_tg", "begin_tg", "end_tg")
  if (has_se) {
    c(v) <- "se2"
  }
  if (flex) {
    c(v) <- "Jtime_tg"
  }
  p[v] <- as.list(paste0(v, usc(bterms$resp)))
  sfx <- str_if("sigma" %in% names(bterms$dpars), "het", "hom")
  sfx <- str_if(has_se, paste0(sfx, "_se"), sfx)
  sfx <- str_if(flex, paste0(sfx, "_flex"), sfx)
  sdist(glue("normal_time_{sfx}"),
    p$mu, p$sigma, p$se2, p$Lcortime,
    p$nobs_tg, p$begin_tg, p$end_tg, p$Jtime_tg
  )
}

stan_log_lik_gaussian_fcor <- function(bterms, ...) {
  if (stan_log_lik_adj(bterms) || has_ad_terms(bterms, "se")) {
    stop2("Invalid addition arguments for this model.")
  }
  p <- stan_log_lik_dpars(bterms, reqn = FALSE)
  p$Lfcor <- paste0("Lfcor", usc(bterms$resp))
  sfx <- str_if("sigma" %in% names(bterms$dpars), "het", "hom")
  sdist(glue("normal_fcor_{sfx}"), p$mu, p$sigma, p$Lfcor)
}

stan_log_lik_gaussian_lagsar <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms, reqn = FALSE)
  p$sigma <- stan_log_lik_add_se(p$sigma, bterms, reqn = FALSE, ...)
  v <- c("lagsar", "Msar", "eigenMsar")
  p[v] <- as.list(paste0(v, usc(bterms$resp)))
  sdist("normal_lagsar", p$mu, p$sigma, p$lagsar, p$Msar, p$eigenMsar)
}

stan_log_lik_gaussian_errorsar <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms, reqn = FALSE)
  p$sigma <- stan_log_lik_add_se(p$sigma, bterms, reqn = FALSE, ...)
  v <- c("errorsar", "Msar", "eigenMsar")
  p[v] <- as.list(paste0(v, usc(bterms$resp)))
  sdist("normal_errorsar", p$mu, p$sigma, p$errorsar, p$Msar, p$eigenMsar)
}

stan_log_lik_student <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms)
  p$sigma <- stan_log_lik_add_se(p$sigma, bterms, ...)
  sdist("student_t", p$nu, p$mu, p$sigma)
}

stan_log_lik_student_mv <- function(bterms, ...) {
  reqn <- stan_log_lik_adj(bterms) || bterms$sigma_pred
  p <- stan_log_lik_dpars(bterms, reqn = reqn, dpars = "nu")
  p$Mu <- paste0("Mu", str_if(reqn, "[n]"))
  p$Sigma <- paste0("Sigma", str_if(bterms$sigma_pred, "[n]"))
  sdist("multi_student_t", p$nu, p$Mu, p$Sigma)
}

stan_log_lik_student_time <- function(bterms, ...) {
  if (stan_log_lik_adj(bterms)) {
    stop2("Invalid addition arguments for this model.")
  }
  has_se <- is.formula(bterms$adforms$se)
  flex <- has_ac_class(bterms$frame$ac, "unstr")
  p <- stan_log_lik_dpars(bterms, reqn = FALSE)
  v <- c("Lcortime", "nobs_tg", "begin_tg", "end_tg")
  if (has_se) {
    c(v) <- "se2"
  }
  if (flex) {
    c(v) <- "Jtime_tg"
  }
  p[v] <- as.list(paste0(v, usc(bterms$resp)))
  sfx <- str_if("sigma" %in% names(bterms$dpars), "het", "hom")
  sfx <- str_if(has_se, paste0(sfx, "_se"), sfx)
  sfx <- str_if(flex, paste0(sfx, "_flex"), sfx)
  sdist(glue("student_t_time_{sfx}"),
    p$nu, p$mu, p$sigma, p$se2, p$Lcortime,
    p$nobs_tg, p$begin_tg, p$end_tg, p$Jtime_tg
  )
}

stan_log_lik_student_fcor <- function(bterms, ...) {
  if (stan_log_lik_adj(bterms) || has_ad_terms(bterms, "se")) {
    stop2("Invalid addition arguments for this model.")
  }
  p <- stan_log_lik_dpars(bterms, reqn = FALSE)
  p$Lfcor <- paste0("Lfcor", usc(bterms$resp))
  sfx <- str_if("sigma" %in% names(bterms$dpars), "het", "hom")
  sdist(glue("student_t_fcor_{sfx}"), p$nu, p$mu, p$sigma, p$Lfcor)
}

stan_log_lik_student_lagsar <- function(bterms,...) {
  p <- stan_log_lik_dpars(bterms, reqn = FALSE)
  p$sigma <- stan_log_lik_add_se(p$sigma, bterms, reqn = FALSE, ...)
  v <- c("lagsar", "Msar", "eigenMsar")
  p[v] <- as.list(paste0(v, usc(bterms$resp)))
  sdist("student_t_lagsar", p$nu, p$mu, p$sigma,
        p$lagsar, p$Msar, p$eigenMsar)
}

stan_log_lik_student_errorsar <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms, reqn = FALSE)
  p$sigma <- stan_log_lik_add_se(p$sigma, bterms, reqn = FALSE, ...)
  v <- c("errorsar", "Msar", "eigenMsar")
  p[v] <- as.list(paste0(v, usc(bterms$resp)))
  sdist("student_t_errorsar", p$nu, p$mu, p$sigma,
        p$errorsar, p$Msar, p$eigenMsar)
}

stan_log_lik_lognormal <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms)
  sdist("lognormal", p$mu, p$sigma)
}

stan_log_lik_shifted_lognormal <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms)
  sdist("lognormal", p$mu, p$sigma, shift = paste0(" - ", p$ndt))
}

stan_log_lik_asym_laplace <- function(bterms,...) {
  p <- stan_log_lik_dpars(bterms, reqn = TRUE)
  sdist("asym_laplace", p$mu, p$sigma, p$quantile, vec = FALSE)
}

stan_log_lik_skew_normal <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms)
  p$sigma <- stan_log_lik_add_se(p$sigma, bterms, ...)
  # required because of CP parameterization of mu and sigma
  mix <- get_mix_id(bterms)
  resp <- usc(bterms$resp)
  reqn <- any(grepl(stan_nn_regex(), c(p$sigma, p$alpha)))
  p$omega <- paste0("omega", mix, resp, str_if(reqn, "[n]"))
  sdist("skew_normal", p$mu, p$omega, p$alpha)
}

stan_log_lik_poisson <- function(bterms, ...) {
  if (use_glm_primitive(bterms)) {
    p <- args_glm_primitive(bterms$dpars$mu, ...)
    out <- sdist("poisson_log_glm", p$x, p$alpha, p$beta)
  } else {
    p <- stan_log_lik_dpars(bterms)
    p$mu <- stan_log_lik_multiply_rate_denom(p$mu, bterms, log = TRUE, ...)
    lpdf <- stan_log_lik_simple_lpdf("poisson", bterms)
    out <- sdist(lpdf, p$mu)
  }
  out
}

stan_log_lik_negbinomial <- function(bterms, ...) {
  if (use_glm_primitive(bterms)) {
    p <- args_glm_primitive(bterms$dpars$mu, ...)
    p$shape <- paste0("shape", usc(bterms$resp))
    out <- sdist("neg_binomial_2_log_glm", p$x, p$alpha, p$beta, p$shape)
  } else {
    p <- stan_log_lik_dpars(bterms)
    p$mu <- stan_log_lik_multiply_rate_denom(p$mu, bterms, log = TRUE, ...)
    p$shape <- stan_log_lik_multiply_rate_denom(p$shape, bterms, ...)
    lpdf <- stan_log_lik_simple_lpdf("neg_binomial_2", bterms)
    out <- sdist(lpdf, p$mu, p$shape)
  }
  out
}

stan_log_lik_negbinomial2 <- function(bterms, ...) {
  if (use_glm_primitive(bterms)) {
    p <- args_glm_primitive(bterms$dpars$mu, ...)
    p$sigma <- paste0("sigma", usc(bterms$resp))
    p$shape <- paste0("inv(", p$sigma, ")")
    out <- sdist("neg_binomial_2_log_glm", p$x, p$alpha, p$beta, p$shape)
  } else {
    p <- stan_log_lik_dpars(bterms)
    p$mu <- stan_log_lik_multiply_rate_denom(p$mu, bterms, log = TRUE, ...)
    p$shape <- stan_log_lik_multiply_rate_denom(
      p$sigma, bterms, transform = "inv", ...
    )
    lpdf <- stan_log_lik_simple_lpdf("neg_binomial_2", bterms)
    out <- sdist(lpdf, p$mu, p$shape)
  }
  out
}

stan_log_lik_geometric <- function(bterms, ...) {
  if (use_glm_primitive(bterms)) {
    p <- args_glm_primitive(bterms$dpars$mu, ...)
    p$shape <- "1"
    out <- sdist("neg_binomial_2_log_glm", p$x, p$alpha, p$beta, p$shape)
  } else {
    p <- stan_log_lik_dpars(bterms)
    p$shape <- "1"
    p$mu <- stan_log_lik_multiply_rate_denom(p$mu, bterms, log = TRUE, ...)
    p$shape <- stan_log_lik_multiply_rate_denom(p$shape, bterms, ...)
    lpdf <- stan_log_lik_simple_lpdf("neg_binomial_2", bterms)
    out <- sdist(lpdf, p$mu, p$shape)
  }
}

stan_log_lik_binomial <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms)
  p$trials <- stan_log_lik_advars(bterms, "trials", ...)$trials
  lpdf <- stan_log_lik_simple_lpdf("binomial", bterms)
  sdist(lpdf, p$trials, p$mu)
}

stan_log_lik_beta_binomial <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms, reqn = TRUE)
  p$trials <- stan_log_lik_advars(bterms, "trials", ...)$trials
  sdist(
    "beta_binomial",
    p$trials,
    paste0(p$mu, " * ", p$phi),
    paste0("(1 - ", p$mu, ") * ", p$phi),
    vec = FALSE
  )
}

stan_log_lik_bernoulli <- function(bterms, ...) {
  if (use_glm_primitive(bterms)) {
    p <- args_glm_primitive(bterms$dpars$mu, ...)
    out <- sdist("bernoulli_logit_glm", p$x, p$alpha, p$beta)
  } else {
    p <- stan_log_lik_dpars(bterms)
    lpdf <- stan_log_lik_simple_lpdf("bernoulli", bterms)
    out <- sdist(lpdf, p$mu)
  }
  out
}

stan_log_lik_discrete_weibull <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms, reqn = TRUE)
  sdist("discrete_weibull", p$mu, p$shape, vec = FALSE)
}

stan_log_lik_com_poisson <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms, reqn = TRUE)
  lpdf <- stan_log_lik_simple_lpdf("com_poisson", bterms)
  sdist(lpdf, p$mu, p$shape, vec = FALSE)
}

stan_log_lik_gamma <- function(bterms, ...) {
  reqn <- stan_log_lik_adj(bterms) || is_pred_dpar(bterms, "shape")
  p <- stan_log_lik_dpars(bterms, reqn = reqn)
  # Stan uses shape-rate parameterization with rate = shape / mean
  div_op <- str_if(reqn, " / ", " ./ ")
  sdist("gamma", p$shape, paste0(p$shape, div_op, p$mu))
}

stan_log_lik_exponential <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms)
  # Stan uses rate parameterization with rate = 1 / mean
  sdist("exponential", paste0("inv(", p$mu, ")"))
}

stan_log_lik_weibull <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms)
  # Stan uses shape-scale parameterization for weibull
  need_dot_div <- !stan_log_lik_adj(bterms) && is_pred_dpar(bterms, "shape")
  div_op <- str_if(need_dot_div, " ./ ", " / ")
  p$scale <- paste0(p$mu, div_op, "tgamma(1 + 1", div_op, p$shape, ")")
  sdist("weibull", p$shape, p$scale)
}

stan_log_lik_frechet <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms)
  # Stan uses shape-scale parameterization for frechet
  need_dot_div <- !stan_log_lik_adj(bterms) && is_pred_dpar(bterms, "nu")
  div_op <- str_if(need_dot_div, " ./ ", " / ")
  p$scale <- paste0(p$mu, div_op, "tgamma(1 - 1", div_op, p$nu, ")")
  sdist("frechet", p$nu, p$scale)
}

stan_log_lik_gen_extreme_value <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms, reqn = TRUE)
  sdist("gen_extreme_value", p$mu, p$sigma, p$xi, vec = FALSE)
}

stan_log_lik_exgaussian <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms)
  sdist(
    "exp_mod_normal", paste0(p$mu, " - ", p$beta),
    p$sigma, paste0("inv(", p$beta, ")")
  )
}

stan_log_lik_inverse.gaussian <- function(bterms, ...) {
  is_pred_shape <- is_pred_dpar(bterms, "shape")
  reqn <- stan_log_lik_adj(bterms) || is_pred_shape
  p <- stan_log_lik_dpars(bterms, reqn = reqn)
  sdist("inv_gaussian", p$mu, p$shape, vec = !is_pred_shape)
}

stan_log_lik_wiener <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms, reqn = TRUE)
  p$dec <- stan_log_lik_advars(bterms, "dec", reqn = TRUE, ...)$dec
  sdist("wiener_diffusion", p$dec, p$bs, p$ndt, p$bias, p$mu, vec = FALSE)
}

stan_log_lik_beta <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms)
  req_dot_multiply <- !stan_log_lik_adj(bterms) && is_pred_dpar(bterms, "phi")
  multiply <- str_if(req_dot_multiply, " .* ", " * ")
  sdist("beta",
    paste0(p$mu, multiply, p$phi),
    paste0("(1 - ", p$mu, ")", multiply, p$phi)
  )
}

stan_log_lik_von_mises <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms)
  sdist("von_mises", p$mu, p$kappa)
}

stan_log_lik_cox <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms)
  c(p) <- stan_log_lik_advars(bterms, c("bhaz", "cbhaz"))
  lpdf <- stan_log_lik_simple_lpdf("cox", bterms)
  sdist(lpdf, p$mu, p$bhaz, p$cbhaz, vec = TRUE)
}

stan_log_lik_cumulative <- function(bterms, ...) {
  if (use_glm_primitive(bterms)) {
    p <- args_glm_primitive(bterms$dpars$mu, ...)
    out <- sdist("ordered_logistic_glm", p$x, p$beta, p$alpha)
  } else {
    out <- stan_log_lik_ordinal(bterms, ...)
  }
  out
}

stan_log_lik_sratio <- function(bterms, ...) {
  stan_log_lik_ordinal(bterms, ...)
}

stan_log_lik_cratio <- function(bterms, ...) {
  stan_log_lik_ordinal(bterms, ...)
}

stan_log_lik_acat <- function(bterms, ...) {
  stan_log_lik_ordinal(bterms, ...)
}

stan_log_lik_categorical <- function(bterms, ...) {
  stopifnot(bterms$family$link == "logit")
  if (use_glm_primitive_categorical(bterms)) {
    bterms1 <- bterms$dpars[[1]]
    bterms1$family <- bterms$family
    p <- args_glm_primitive(bterms1, ...)
    out <- sdist("categorical_logit_glm", p$x, p$alpha, p$beta)
  } else {
    p <- stan_log_lik_dpars(bterms, reqn = TRUE, dpars = "mu", type = "multi")
    out <- sdist("categorical_logit", p$mu, vec = FALSE)
  }
  out
}

stan_log_lik_multinomial <- function(bterms, ...) {
  stopifnot(bterms$family$link == "logit")
  p <- stan_log_lik_dpars(bterms, reqn = TRUE, dpars = "mu", type = "multi")
  sdist("multinomial_logit2", p$mu, vec = FALSE)
}

stan_log_lik_dirichlet <- function(bterms, ...) {
  stopifnot(bterms$family$link == "logit")
  mu <- stan_log_lik_dpars(bterms, reqn = TRUE, dpars = "mu", type = "multi")$mu
  reqn_phi <- is_pred_dpar(bterms, "phi")
  phi <- stan_log_lik_dpars(bterms, reqn = reqn_phi, dpars = "phi")$phi
  sdist("dirichlet_logit", mu, phi, vec = FALSE)
}

stan_log_lik_dirichlet2 <- function(bterms,...) {
  mu <- stan_log_lik_dpars(bterms, reqn = TRUE, dpars = "mu", type = "multi")$mu
  sdist("dirichlet", mu, vec = FALSE)
}

stan_log_lik_logistic_normal <- function(bterms, ...) {
  stopifnot(bterms$family$link == "identity")
  resp <- usc(bterms$resp)
  mix <- get_mix_id(bterms)
  p <- stan_log_lik_dpars(bterms, reqn = TRUE, type = "multi")
  p$Llncor <- glue("Llncor{mix}{resp}")
  p$refcat <- get_refcat(bterms$family, int = TRUE)
  sdist(
    "logistic_normal_cholesky_cor",
    p$mu, p$sigma, p$Llncor, p$refcat,
    vec = FALSE
  )
}

stan_log_lik_ordinal <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms, reqn = TRUE)
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
    p$Jthres <- stan_log_lik_advars(bterms, "Jthres", reqn = TRUE, ...)$Jthres
    p$thres <- "merged_Intercept"
  } else {
    p$thres <- "Intercept"
  }
  resp <- usc(bterms$resp)
  mix <- get_mix_id(bterms)
  prefix <- paste0(str_if(nzchar(mix), paste0("_mu", mix)), resp)
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
  sdist(lpdf, p$mu, p$disc, p$thres, p$Jthres, vec = FALSE)
}

stan_log_lik_hurdle_poisson <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms, reqn = TRUE)
  lpdf <- stan_log_lik_simple_lpdf("hurdle_poisson", bterms)
  lpdf <- paste0(lpdf, stan_log_lik_dpar_usc_logit(bterms, "hu"))
  sdist(lpdf, p$mu, p$hu, vec = FALSE)
}

stan_log_lik_hurdle_negbinomial <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms, reqn = TRUE)
  lpdf <- stan_log_lik_simple_lpdf("hurdle_neg_binomial", bterms)
  lpdf <- paste0(lpdf, stan_log_lik_dpar_usc_logit(bterms, "hu"))
  sdist(lpdf, p$mu, p$shape, p$hu, vec = FALSE)
}

stan_log_lik_hurdle_gamma <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms, reqn = TRUE)
  usc_logit <- stan_log_lik_dpar_usc_logit(bterms, "hu")
  lpdf <- paste0("hurdle_gamma", usc_logit)
  # Stan uses shape-rate parameterization for gamma with rate = shape / mean
  sdist(lpdf, p$shape, paste0(p$shape, " / ", p$mu), p$hu, vec = FALSE)
}

stan_log_lik_hurdle_lognormal <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms, reqn = TRUE)
  usc_logit <- stan_log_lik_dpar_usc_logit(bterms, "hu")
  lpdf <- paste0("hurdle_lognormal", usc_logit)
  sdist(lpdf, p$mu, p$sigma, p$hu, vec = FALSE)
}

stan_log_lik_hurdle_cumulative <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms, reqn = TRUE)
  if (use_ordered_builtin(bterms, "logit")) {
    lpdf <- "hurdle_cumulative_ordered_logistic"
  } else if (use_ordered_builtin(bterms, "probit")) {
    lpdf <- "hurdle_cumulative_ordered_probit"
  } else {
    lpdf <- paste0(bterms$family$family, "_", bterms$family$link)
  }
  if (has_thres_groups(bterms)) {
    str_add(lpdf) <- "_merged"
    p$Jthres <- stan_log_lik_advars(bterms, "Jthres", reqn = TRUE, ...)$Jthres
    p$thres <- "merged_Intercept"
  } else {
    p$thres <- "Intercept"
  }
  resp <- usc(bterms$resp)
  mix <- get_mix_id(bterms)
  prefix <- paste0(str_if(nzchar(mix), paste0("_mu", mix)), resp)
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
  sdist(lpdf, p$mu, p$hu, p$disc, p$thres, p$Jthres, vec = FALSE)
}

stan_log_lik_zero_inflated_poisson <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms, reqn = TRUE)
  lpdf <- stan_log_lik_simple_lpdf("zero_inflated_poisson", bterms)
  lpdf <- paste0(lpdf, stan_log_lik_dpar_usc_logit(bterms, "zi"))
  sdist(lpdf, p$mu, p$zi, vec = FALSE)
}

stan_log_lik_zero_inflated_negbinomial <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms, reqn = TRUE)
  lpdf <- stan_log_lik_simple_lpdf("zero_inflated_neg_binomial", bterms)
  lpdf <- paste0(lpdf, stan_log_lik_dpar_usc_logit(bterms, "zi"))
  sdist(lpdf, p$mu, p$shape, p$zi, vec = FALSE)
}

stan_log_lik_zero_inflated_binomial <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms, reqn = TRUE)
  p$trials <- stan_log_lik_advars(bterms, "trials", reqn = TRUE, ...)$trials
  lpdf <- "zero_inflated_binomial"
  lpdf <- stan_log_lik_simple_lpdf(lpdf, bterms, sep = "_b")
  lpdf <- paste0(lpdf, stan_log_lik_dpar_usc_logit(bterms, "zi"))
  sdist(lpdf, p$trials, p$mu, p$zi, vec = FALSE)
}

stan_log_lik_zero_inflated_beta_binomial <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms, reqn = TRUE)
  p$trials <- stan_log_lik_advars(bterms, "trials", reqn = TRUE, ...)$trials
  lpdf <- "zero_inflated_beta_binomial"
  lpdf <- paste0(lpdf, stan_log_lik_dpar_usc_logit(bterms, "zi"))
  sdist(lpdf, p$trials, p$mu, p$phi, p$zi, vec = FALSE)
}

stan_log_lik_zero_inflated_beta <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms, reqn = TRUE)
  usc_logit <- stan_log_lik_dpar_usc_logit(bterms, "zi")
  lpdf <- paste0("zero_inflated_beta", usc_logit)
  sdist(lpdf, p$mu, p$phi, p$zi, vec = FALSE)
}

stan_log_lik_zero_one_inflated_beta <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms, reqn = TRUE)
  sdist("zero_one_inflated_beta", p$mu, p$phi, p$zoi, p$coi, vec = FALSE)
}

stan_log_lik_zero_inflated_asym_laplace <- function(bterms, ...) {
  p <- stan_log_lik_dpars(bterms, reqn = TRUE)
  usc_logit <- stan_log_lik_dpar_usc_logit(bterms, "zi")
  lpdf <- paste0("zero_inflated_asym_laplace", usc_logit)
  sdist(lpdf, p$mu, p$sigma, p$quantile, p$zi, vec = FALSE)
}

stan_log_lik_custom <- function(bterms, threads = NULL, ...) {
  family <- bterms$family
  no_loop <- isFALSE(family$loop)
  if (no_loop && (stan_log_lik_adj(bterms))) {
    stop2("This model requires evaluating the custom ",
          "likelihood as a loop over observations.")
  }
  resp <- usc(bterms$resp)
  p <- stan_log_lik_dpars(bterms, reqn = !no_loop)
  mix <- get_mix_id(bterms)
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
  sdist(family$name, p[dpars], p$thres, vars, vec = no_loop)
}

# ordinal log-probability density functions in Stan language
# @return a character string
stan_ordinal_lpmf <- function(family, link) {
  family <- as_one_character(family)
  link <- as_one_character(link)
  inv_link <- stan_inv_link(link)
  th <- function(k) {
    # helper function generating stan code inside inv_link(.)
    if (family %in% c("cumulative", "sratio")) {
      out <- glue("thres[{k}] - mu")
    } else if (family %in% c("cratio", "acat")) {
      out <- glue("mu - thres[{k}]")
    }
    glue("disc * ({out})")
  }
  out <- glue(
    "  /* {family}-{link} log-PDF for a single response\n",
    "   * Args:\n",
    "   *   y: response category\n",
    "   *   mu: latent mean parameter\n",
    "   *   disc: discrimination parameter\n",
    "   *   thres: ordinal thresholds\n",
    "   * Returns:\n",
    "   *   a scalar to be added to the log posterior\n",
    "   */\n",
    "   real {family}_{link}_lpmf(int y, real mu, real disc, vector thres) {{\n"
  )
  # define the function body
  if (family == "cumulative") {
    if (inv_link == "inv_logit") {
      str_add(out) <- glue(
        "     int nthres = num_elements(thres);\n",
        "     if (y == 1) {{\n",
        "       return log_inv_logit({th(1)});\n",
        "     }} else if (y == nthres + 1) {{\n",
        "       return log1m_inv_logit({th('nthres')});\n",
        "     }} else {{\n",
        "       return log_inv_logit_diff({th('y')}, {th('y - 1')});\n",
        "     }}\n",
        "   }}\n"
      )
    } else {
      str_add(out) <- glue(
        "     int nthres = num_elements(thres);\n",
        "     real p;\n",
        "     if (y == 1) {{\n",
        "       p = {inv_link}({th(1)});\n",
        "     }} else if (y == nthres + 1) {{\n",
        "       p = 1 - {inv_link}({th('nthres')});\n",
        "     }} else {{\n",
        "       p = {inv_link}({th('y')}) -\n",
        "           {inv_link}({th('y - 1')});\n",
        "     }}\n",
        "     return log(p);\n",
        "   }}\n"
      )
    }
  } else if (family %in% c("sratio", "cratio")) {
    # TODO: support 'softit' link as well
    if (inv_link == "inv_cloglog") {
      qk <- str_if(
        family == "sratio",
        "-exp({th('k')})",
        "log1m_exp(-exp({th('k')}))"
      )
    } else if (inv_link == "inv_logit") {
      qk <- str_if(
        family == "sratio",
        "log1m_inv_logit({th('k')})",
        "log_inv_logit({th('k')})"
      )
    } else if (inv_link == "Phi") {
      qk <- str_if(
        family == "sratio",
        "std_normal_lccdf({th('k')})",
        "std_normal_lcdf({th('k')})"
      )
    } else if (inv_link == "Phi_approx") {
      qk <- str_if(
        family == "sratio",
        "log1m_inv_logit(0.07056 * pow({th('k')}, 3.0) + 1.5976 * {th('k')})",
        "log_inv_logit(0.07056 * pow({th('k')}, 3.0) + 1.5976 * {th('k')})"
      )
    } else if (inv_link == "inv_cauchit") {
      qk <- str_if(
        family == "sratio",
        "cauchy_lccdf({th('k')}|0,1)",
        "cauchy_lcdf({th('k')}|0,1)"
      )
    }
    qk <- glue(qk)
    str_add(out) <- glue(
      "     int nthres = num_elements(thres);\n",
      "     vector[nthres + 1] p;\n",
      "     vector[nthres] q;\n",
      "     int k = 1;\n",
      "     while (k <= min(y, nthres)) {{\n",
      "       q[k] = {qk};\n",
      "       p[k] = log1m_exp(q[k]);\n",
      "       for (kk in 1:(k - 1)) p[k] = p[k] + q[kk];\n",
      "       k += 1;\n",
      "     }}\n",
      "     if (y == nthres + 1) {{\n",
      "       p[nthres + 1] = sum(q);\n",
      "     }}\n",
      "     return p[y];\n",
      "   }}\n"
    )
  } else if (family == "acat") {
    if (inv_link == "inv_logit") {
      str_add(out) <- glue(
        "     int nthres = num_elements(thres);\n",
        "     vector[nthres + 1] p = append_row(0, cumulative_sum(disc * (mu - thres)));\n",
        "     return p[y] - log_sum_exp(p);\n",
        "   }}\n"
      )
    } else {
      str_add(out) <- glue(
        "     int nthres = num_elements(thres);\n",
        "     vector[nthres + 1] p;\n",
        "     vector[nthres] q;\n",
        "     for (k in 1:(nthres))\n",
        "       q[k] = {inv_link}({th('k')});\n",
        "     for (k in 1:(nthres + 1)) {{\n",
        "       p[k] = 1.0;\n",
        "       for (kk in 1:(k - 1)) p[k] = p[k] * q[kk];\n",
        "       for (kk in k:(nthres)) p[k] = p[k] * (1 - q[kk]);\n",
        "     }}\n",
        "     return log(p[y]) - log(sum(p));\n",
        "   }}\n"
      )
    }
  }
  # lpmf function for multiple merged thresholds
  str_add(out) <- glue(
    "  /* {family}-{link} log-PDF for a single response and merged thresholds\n",
    "   * Args:\n",
    "   *   y: response category\n",
    "   *   mu: latent mean parameter\n",
    "   *   disc: discrimination parameter\n",
    "   *   thres: vector of merged ordinal thresholds\n",
    "   *   j: start and end index for the applid threshold within 'thres'\n",
    "   * Returns:\n",
    "   *   a scalar to be added to the log posterior\n",
    "   */\n",
    "   real {family}_{link}_merged_lpmf(",
    "int y, real mu, real disc, vector thres, array[] int j) {{\n",
    "     return {family}_{link}_lpmf(y | mu, disc, thres[j[1]:j[2]]);\n",
    "   }}\n"
  )
  if (family == "cumulative" && link %in% c("logit", "probit")) {
    # use the more efficient ordered_link functions when disc == 1
    sfx <- str_if(link == "logit", "logistic", link)
    str_add(out) <- glue(
      "  /* ordered-{sfx} log-PDF for a single response and merged thresholds\n",
      "   * Args:\n",
      "   *   y: response category\n",
      "   *   mu: latent mean parameter\n",
      "   *   thres: vector of merged ordinal thresholds\n",
      "   *   j: start and end index for the applid threshold within 'thres'\n",
      "   * Returns:\n",
      "   *   a scalar to be added to the log posterior\n",
      "   */\n",
      "   real ordered_{sfx}_merged_lpmf(",
      "int y, real mu, vector thres, array[] int j) {{\n",
      "     return ordered_{sfx}_lpmf(y | mu, thres[j[1]:j[2]]);\n",
      "   }}\n"
    )
  }
  out
}

# log probability density for hurdle ordinal models
# TODO: generalize to non-cumulative families?
# @return a character string
stan_hurdle_ordinal_lpmf <- function(family, link) {
  family <- as_one_character(family)
  link <- as_one_character(link)
  stopifnot(family == "hurdle_cumulative")
  inv_link <- stan_inv_link(link)
  th <- function(k) {
    out <- glue("thres[{k}] - mu")
    glue("disc * ({out})")
  }
  out <- glue(
    "  /* {family}-{link} log-PDF for a single response\n",
    "   * Args:\n",
    "   *   y: response category\n",
    "   *   mu: latent mean parameter\n",
    "   *   hu: hurdle probability\n",
    "   *   disc: discrimination parameter\n",
    "   *   thres: ordinal thresholds\n",
    "   * Returns:\n",
    "   *   a scalar to be added to the log posterior\n",
    "   */\n",
    "   real {family}_{link}_lpmf(int y, real mu, real hu, real disc, vector thres) {{\n",
    "\n"
  )
  # define the function body
  if (inv_link == "inv_logit") {
    str_add(out) <- glue(
      "     int nthres = num_elements(thres);\n",
      "     if (y == 0) {{\n",
      "       return bernoulli_lpmf(1 | hu);\n",
      "     }} else if (y == 1) {{\n",
      "       return log_inv_logit({th(1)}) +\n",
      "                bernoulli_lpmf(0 | hu);\n",
      "     }} else if (y == nthres + 2) {{\n",
      "       return log1m_inv_logit({th('nthres')}) +\n",
      "                bernoulli_lpmf(0 | hu);\n",
      "     }} else {{\n",
      "       return log_inv_logit_diff({th('y')}, {th('y - 1')}) +\n",
      "                bernoulli_lpmf(0 | hu) ;\n",
      "     }}\n",
      "   }}\n"
    )
  } else {
    str_add(out) <- glue(
      "     int nthres = num_elements(thres);\n",
      "     real p;\n",
      "     if (y == 0){{\n",
      "       p = hu;\n",
      "     }} else if (y == 1) {{\n",
      "       p = {inv_link}({th(1)}) * (1 - hu);\n",
      "     }} else if (y == nthres + 1) {{\n",
      "       p = (1 - {inv_link}({th('nthres')})) * (1 - hu);\n",
      "     }} else {{\n",
      "       p = ({inv_link}({th('y')}) -\n",
      "           {inv_link}({th('y - 1')})) * (1 - hu);\n",
      "     }}\n",
      "     return log(p);\n",
      "   }}\n"
    )
  }

  # lpmf function for multiple merged thresholds
  str_add(out) <- glue(
    "  /* {family}-{link} log-PDF for a single response and merged thresholds\n",
    "   * Args:\n",
    "   *   y: response category\n",
    "   *   mu: latent mean parameter\n",
    "   *   hu: hurdle probability\n",
    "   *   disc: discrimination parameter\n",
    "   *   thres: vector of merged ordinal thresholds\n",
    "   *   j: start and end index for the applid threshold within 'thres'\n",
    "   * Returns:\n",
    "   *   a scalar to be added to the log posterior\n",
    "   */\n",
    "   real {family}_{link}_merged_lpmf(",
    "int y, real mu, real hu, real disc, vector thres, array[] int j) {{\n",
    "     return {family}_{link}_lpmf(y | mu, hu, disc, thres[j[1]:j[2]]);\n",
    "   }}\n"
  )

  if (link %in% c("logit", "probit")) {
    # use the more efficient ordered_link functions when disc == 1
    sfx <- str_if(link == "logit", "logistic", link)
    str_add(out) <- glue(
      "\n",
      "   // Use more efficient ordered_{sfx} function with disc == 1\n",
      "   real hurdle_cumulative_ordered_{sfx}_lpmf(int y, real mu, real hu, real disc, vector thres) {{\n",
      "     if (y == 0) {{\n",
      "       return bernoulli_lpmf(1 | hu);\n",
      "     }} else {{\n",
      "       return ordered_{sfx}_lpmf(y | mu, thres) +\n",
      "                bernoulli_lpmf(0 | hu);\n",
      "     }}\n",
      "   }}\n"
    )
    str_add(out) <- glue(
      "  /* use ordered-{sfx} log-PDF for a single response and merged thresholds\n",
      "   * Args:\n",
      "   *   y: response category\n",
      "   *   mu: latent mean parameter\n",
      "   *   hu: hurdle probability\n",
      "   *   thres: vector of merged ordinal thresholds\n",
      "   *   j: start and end index for the applid threshold within 'thres'\n",
      "   * Returns:\n",
      "   *   a scalar to be added to the log posterior\n",
      "   */\n",
      "   real hurdle_cumulative_ordered_{sfx}_merged_lpmf(",
      "int y, real mu, real hu, real disc, vector thres, array[] int j) {{\n",
      "     return hurdle_cumulative_ordered_{sfx}_lpmf(y | mu, hu, disc, thres[j[1]:j[2]]);\n",
      "   }}\n"
    )
  }
  out
}

# use a Stan GLM primitive function?
use_glm_primitive <- function(bterms) {
  stopifnot(is.brmsterms(bterms))
  # the model can only have a single predicted parameter
  # and no additional residual or autocorrelation structure
  mu <- bterms$dpars[["mu"]]
  non_glm_adterms <- c("se", "weights", "thres", "cens", "trunc", "rate")
  if (!is.btl(mu) || length(bterms$dpars) > 1L ||
      isTRUE(bterms$rescor) || is.formula(mu$ac) ||
      has_ad_terms(bterms, non_glm_adterms)) {
    return(FALSE)
  }
  # some primitives do not support special terms in the way
  # required by brms' Stan code generation
  allow_special_terms <- !mu$family$family %in% c("cumulative", "categorical")
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

# use Stan's categorical GLM primitive function?
use_glm_primitive_categorical <- function(bterms) {
  stopifnot(is.brmsterms(bterms))
  if (!is_categorical(bterms)) {
    return(FALSE)
  }
  tmp <- bterms
  tmp$dpars <- list()
  # we know that all dpars in categorical models are mu parameters
  out <- rep(FALSE, length(bterms$dpars))
  for (i in seq_along(bterms$dpars)) {
    tmp$dpars$mu <- bterms$dpars[[i]]
    tmp$dpars$mu$family <- bterms$family
    out[i] <- use_glm_primitive(tmp) &&
      # the design matrix of all mu parameters must match
      all.equal(tmp$dpars$mu$fe, bterms$dpars[[1]]$fe)
  }
  all(out)
}

# standard arguments for primitive Stan GLM functions
# @param bterms a btl object
# @return a named list of Stan code snippets
args_glm_primitive <- function(bterms, threads = NULL, ...) {
  stopifnot(is.btl(bterms))
  resp <- usc(bterms$resp)
  decomp <- get_decomp(bterms$fe)
  center_X <- stan_center_X(bterms)
  slice <- stan_slice(threads)
  sfx_X <- sfx_b <- ""
  if (decomp == "QR") {
    sfx_X <- sfx_b <- "Q"
  } else if (center_X) {
    sfx_X <- "c"
  }
  is_categorical <- is_categorical(bterms)
  if (is_categorical) {
    sfx_X <- glue("{sfx_X}_{bterms$dpar}")
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
      if (is_categorical) {
        alpha <- glue("rep_vector(0, ncat{resp})")
      } else {
        alpha <- "0"
      }
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
# @param dist name of the distribution in Stan language
# @param vec does the distribution have a vectorized version?
# @param shift Stan code for shifting the likelihood in shifted_* families
sdist <- function(dist, ..., vec = TRUE, shift = "") {
  args <- sargs(...)
  structure(nlist(dist, args, vec, shift), class = "sdist")
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
