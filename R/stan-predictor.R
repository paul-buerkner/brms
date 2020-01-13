# unless otherwise specified, functions return a named list 
# of Stan code snippets to be pasted together later on

# generate stan code for predictor terms
stan_predictor <- function(x, ...) {
  UseMethod("stan_predictor")
}

# combine effects for the predictors of a single (non-linear) parameter
# @param primitive set to TRUE only if primitive Stan GLM functions should 
#   be used which compute the predictor term internally
# @param ... arguments passed to the underlying effect-specific functions
#' @export
stan_predictor.btl <- function(x, primitive = FALSE, ...) {
  out <- collapse_lists(
    stan_fe(x, ...),
    stan_thres(x, ...),
    stan_sp(x, ...),
    stan_cs(x, ...),
    stan_sm(x, ...),
    stan_gp(x, ...),
    stan_ac(x, ...),
    stan_offset(x, ...),
    stan_rate(x, ...),
    stan_bhaz(x, ...),
    stan_special_prior_global(x, ...)
  )
  if (!primitive) {
    out <- stan_eta_combine(out, bterms = x, ...) 
  }
  out
}

# prepare Stan code for non-linear models
# @param names of the non-linear parameters
# @param ilink character vector of length 2 defining the link to be applied
#' @export
stan_predictor.btnl <- function(x, nlpars, ilink = c("", ""), ...) {
  stopifnot(length(ilink) == 2L)
  out <- list()
  resp <- usc(x$resp)
  par <- combine_prefix(x, keep_mu = TRUE, nlp = TRUE)
  # prepare non-linear model
  n <- str_if(x$loop, "[n] ", " ")
  new_nlpars <- glue(" nlp{usc(x$resp)}_{nlpars}{n}")
  # covariates in the non-linear model
  covars <- wsp(all.vars(rhs(x$covars)))
  new_covars <- NULL
  if (length(covars)) {
    p <- usc(combine_prefix(x))
    str_add(out$data) <- glue(
      "  int<lower=1> KC{p};  // number of covariates\n",
      "  matrix[N{resp}, KC{p}] C{p};  // matrix of covariates\n"
    )
    kc <- seq_along(covars)
    str_add(out$tdata_def) <- glue( 
      "  // extract covariate vectors for faster indexing\n",
      cglue("  vector[N{resp}] C{p}_{kc} = C{p}[, {kc}];\n")
    )
    new_covars <- glue(" C{p}_{kc}{n}")
  }
  # add whitespaces to be able to replace parameters and covariates
  syms <- c(
    "+", "-", "*", "/", "%", "^", ".*", "./", "'", ")", "(", 
    ",", "==", "!=", "<=", ">=", "<", ">", "!", "&&", "||" 
  )
  regex <- glue("(?<!\\.){escape_all(syms)}(?!=)")
  out$eta <- rm_wsp(collapse(deparse(x$formula[[2]])))
  out$eta <- wsp(rename(out$eta, regex, wsp(syms), fixed = FALSE, perl = TRUE)) 
  out$eta <- rename(out$eta, 
    c(wsp(nlpars), covars, " ( ", " ) "), 
    c(new_nlpars, new_covars, "(", ")")
  )
  str_add_list(out) <- stan_rate(x, loop = x$loop)
  # possibly transform eta in the transformed params block
  str_add(out$model_def) <- glue(
    "  // initialize non-linear predictor term\n",
    "  vector[N{resp}] {par};\n"
  )
  # make sure mu comes last as it might depend on other parameters
  is_mu <- isTRUE("mu" %in% dpar_class(x[["dpar"]]))
  position <- str_if(is_mu, "model_comp_mu_link", "model_comp_dpar_link")
  if (x$loop) {
    str_add(out[[position]]) <- glue(
      "  for (n in 1:N{resp}) {{\n",
      "    // compute non-linear predictor values\n",
      "    {par}[n] = {ilink[1]}{trimws(out$eta)}{ilink[2]};\n",
      "  }}\n"
    )
  } else {
    str_add(out[[position]]) <- glue(
      "  // compute non-linear predictor values\n",
      "  {par} = {ilink[1]}{trimws(out$eta)}{ilink[2]};\n"
    )
  }
  out$eta <- NULL
  str_add_list(out) <- stan_thres(x, ...)
  str_add_list(out) <- stan_bhaz(x, ...)
  out
}

# Stan code for distributional parameters
# @param rescor is this predictor part of a MV model estimating rescor?
#' @export
stan_predictor.brmsterms <- function(x, data, prior, rescor = FALSE, ...) {
  px <- check_prefix(x)
  resp <- usc(combine_prefix(px))
  data <- subset_data(data, x)
  out <- list(stan_response(x, data = data))
  valid_dpars <- valid_dpars(x)
  args <- nlist(data, prior, nlpars = names(x$nlpars), ...)
  args$primitive <- use_glm_primitive(x)
  for (nlp in names(x$nlpars)) {
    nlp_args <- list(x$nlpars[[nlp]])
    out[[nlp]] <- do_call(stan_predictor, c(nlp_args, args))
  }
  for (dp in valid_dpars) {
    dp_terms <- x$dpars[[dp]]
    if (is.btl(dp_terms) || is.btnl(dp_terms)) {
      ilink <- stan_eta_ilink(dp, bterms = x, resp = resp)
      dp_args <- list(dp_terms, ilink = ilink)
      out[[dp]] <- do_call(stan_predictor, c(dp_args, args))
    } else if (is.numeric(x$fdpars[[dp]]$value)) {
      dp_def <- stan_dpar_defs(dp, resp, family = x$family, fixed = TRUE)
      out[[dp]] <- list(data = dp_def)
    } else if (is.character(x$fdpars[[dp]]$value)) {
      if (!x$fdpars[[dp]]$value %in% valid_dpars) {
        stop2("Parameter '", x$fdpars[[dp]]$value, "' cannot be found.")
      }
      dp_def <- stan_dpar_defs(dp, resp, family = x$family)
      dp_co <- glue("  {dp}{resp} = {x$fdpars[[dp]]$value}{resp};\n")
      out[[dp]] <- list(tpar_def = dp_def, tpar_comp = dp_co)
    } else {
      dp_def <- stan_dpar_defs(dp, resp, family = x$family)
      dp_def_tmp <- stan_dpar_defs_tmp(dp, resp, family = x$family)
      if (nzchar(dp_def_tmp)) {
        dp_prior <- stan_prior(
          prior, dp, prefix = "tmp_", suffix = resp, px = px
        )
        out[[dp]] <- list(par = dp_def_tmp, prior = dp_prior)
      } else if (nzchar(dp_def)) {
        dp_prior <- stan_prior(prior, dp, suffix = resp, px = px)
        out[[dp]] <- list(par = dp_def, prior = dp_prior)
      }
    }
  }
  out <- collapse_lists(ls = out)
  # str_add_list(out) <- stan_autocor(x, prior = prior)
  str_add_list(out) <- stan_mixture(x, data = data, prior = prior)
  str_add_list(out) <- stan_dpar_transform(x)
  out
}

#' @export
stan_predictor.mvbrmsterms <- function(x, prior, ...) {
  out <- lapply(x$terms, stan_predictor, prior = prior, rescor = x$rescor, ...)
  out <- collapse_lists(ls = out)
  str_add(out$data, start = TRUE) <- 
    "  int<lower=1> N;  // number of observations\n"
  if (x$rescor) {
    # we already know at this point that all families are identical
    adforms <- lapply(x$terms, "[[", "adforms")
    adnames <- unique(ulapply(adforms, names))
    adallowed <- c("se", "weights", "mi")
    if (!all(adnames %in% adallowed))  {
      stop2("Only ", collapse_comma(adallowed), " are supported ", 
            "addition arguments when 'rescor' is estimated.")
    }
    family <- family_names(x)[1]
    stopifnot(family %in% c("gaussian", "student"))
    resp <- x$responses
    nresp <- length(resp)
    str_add(out$model_def) <- glue( 
      "  // multivariate predictor array\n",
      "  vector[nresp] Mu[N];\n"
    )
    str_add(out$model_comp_mvjoin) <- glue(
      "    Mu[n] = {stan_vector(glue('mu_{resp}[n]'))};\n"
    )
    str_add(out$data) <- glue(
      "  int<lower=1> nresp;  // number of responses\n",   
      "  int nrescor;  // number of residual correlations\n"
    )
    str_add(out$tdata_def) <- glue(
      "  vector[nresp] Y[N];  // response array\n"
    )
    str_add(out$tdata_comp) <- glue(
      "  for (n in 1:N) {{\n",
      "    Y[n] = {stan_vector(glue('Y_{resp}[n]'))};\n",
      "  }}\n"
    )
    if (any(adnames %in% "weights")) {
      str_add(out$tdata_def) <- glue(
        "  // weights of the pointwise log-likelihood\n",
        "  vector<lower=0>[N] weights = weights_{resp[1]};\n" 
      )
    }
    miforms <- rmNULL(lapply(adforms, "[[", "mi"))
    if (length(miforms)) {
      str_add(out$model_def) <- "  vector[nresp] Yl[N] = Y;\n"
      for (i in seq_along(miforms)) {
        j <- match(names(miforms)[i], resp)
        str_add(out$model_comp_mvjoin) <- glue(
          "    Yl[n][{j}] = Yl_{resp[j]}[n];\n"
        )
      }
    }
    str_add(out$par) <- glue(
      "  // parameters for multivariate linear models\n",
      "  cholesky_factor_corr[nresp] Lrescor;\n"
    )
    str_add(out$prior) <- stan_prior(prior, class = "Lrescor")
    if (family == "student") {
      str_add(out$par) <- stan_dpar_defs("nu")
      str_add(out$prior) <- stan_prior(prior, class = "nu")
    } 
    sigma <- ulapply(x$terms, stan_sigma_transform)
    if (any(grepl("\\[n\\]", sigma))) {
      str_add(out$model_def) <- "  vector[nresp] sigma[N];\n"
      str_add(out$model_comp_mvjoin) <- glue(
        "    sigma[n] = {stan_vector(sigma)};\n"
      )
      if (family == "gaussian") {
        str_add(out$model_def) <- glue(
          "  // cholesky factor of residual covariance matrix\n",
          "  matrix[nresp, nresp] LSigma[N];\n"
        )
        str_add(out$model_comp_mvjoin) <- glue(
          "    LSigma[n] = diag_pre_multiply(sigma[n], Lrescor);\n"
        )
      } else if (family == "student") {
        str_add(out$model_def) <- glue(
          "  // residual covariance matrix\n",
          "  matrix[nresp, nresp] Sigma[N];\n"
        )
        str_add(out$model_comp_mvjoin) <- glue(
          "    Sigma[n] = multiply_lower_tri_self_transpose(", 
          "diag_pre_multiply(sigma[n], Lrescor));\n" 
        )
      }
    } else {
      str_add(out$model_def) <- glue(
        "  vector[nresp] sigma = {stan_vector(sigma)};\n"
      )
      if (family == "gaussian") {
        str_add(out$model_def) <- glue(
          "  // cholesky factor of residual covariance matrix\n",
          "  matrix[nresp, nresp] LSigma = ",
          "diag_pre_multiply(sigma, Lrescor);\n"
        )
      } else if (family == "student") {
        str_add(out$model_def) <- glue(
          "  // residual covariance matrix\n",
          "  matrix[nresp, nresp] Sigma = ",
          "multiply_lower_tri_self_transpose(", 
          "diag_pre_multiply(sigma, Lrescor));\n"
        )
      }
    }
    str_add(out$gen_def) <- glue(
      "  // residual correlations\n",
      "  corr_matrix[nresp] Rescor",
      " = multiply_lower_tri_self_transpose(Lrescor);\n",
      "  vector<lower=-1,upper=1>[nrescor] rescor;\n"
    )
    str_add(out$gen_comp) <- stan_cor_gen_comp("rescor", "nresp")
    out$model_comp_mvjoin <- paste0(
      "  // combine univariate parameters\n",
      "  for (n in 1:N) {\n", out$model_comp_mvjoin, "  }\n"
    )
  }
  out
}

# Stan code for population-level effects
stan_fe <- function(bterms, data, prior, stanvars, ...) {
  out <- list()
  family <- bterms$family
  fixef <- colnames(data_fe(bterms, data)$X)
  sparse <- is_sparse(bterms$fe)
  decomp <- get_decomp(bterms$fe)
  if (length(fixef) < 2L) {
    # decompositions require at least two predictors
    decomp <- "none"
  }
  center_X <- stan_center_X(bterms)
  ct <- str_if(center_X, "c")
  # remove the intercept from the design matrix?
  if (center_X) {
    fixef <- setdiff(fixef, "Intercept")
  }
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  resp <- usc(px$resp)
  if (length(fixef)) {
    str_add(out$data) <- glue( 
      "  int<lower=1> K{p};",
      "  // number of population-level effects\n", 
      "  matrix[N{resp}, K{p}] X{p};",
      "  // population-level design matrix\n"
    )
    if (sparse) {
      if (decomp != "none") {
        stop2("Cannot use ", decomp, " decomposition for sparse matrices.")
      }
      str_add(out$tdata_def) <- glue(
        "  // sparse matrix representation of X{p}\n",
        "  vector[rows(csr_extract_w(X{p}))] wX{p}", 
        " = csr_extract_w(X{p});\n",
        "  int vX{p}[size(csr_extract_v(X{p}))]",
        " = csr_extract_v(X{p});\n",
        "  int uX{p}[size(csr_extract_u(X{p}))]",
        " = csr_extract_u(X{p});\n"
      )
    }
    # prepare population-level coefficients
    bound <- get_bound(prior, class = "b", px = px)
    use_horseshoe <- stan_use_horseshoe(bterms, prior)
    if (decomp == "none") {
      bsuffix <- ""
      comment_b <- "  // population-level effects"
      stan_def_b <- glue("  vector{bound}[K{ct}{p}] b{p};{comment_b}\n")
      if (use_horseshoe) {
        str_add(out$tpar_def) <- stan_def_b
      } else if (!glue("b{p}") %in% names(stanvars)) {
        str_add(out$par) <- stan_def_b
      }
    } else {
      stopifnot(decomp == "QR")
      if (nzchar(bound)) {
        stop2("Cannot impose bounds on decomposed coefficients.")
      }
      bsuffix <- "Q"
      comment_bQ <- "  // regression coefficients at QR scale"
      stan_def_bQ <- glue("  vector[K{ct}{p}] bQ{p};{comment_bQ}\n")
      if (use_horseshoe) {
        str_add(out$tpar_def) <- stan_def_bQ
      } else {
        str_add(out$par) <- stan_def_bQ
      }
      str_add(out$gen_def) <- glue(
        "  // obtain the actual coefficients\n",
        "  vector[K{ct}{p}] b{p} = XR{p}_inv * bQ{p};\n"
      )
    }
    str_add(out$prior) <- stan_prior(
      prior, class = "b", coef = fixef, px = px, 
      suffix = glue("{bsuffix}{p}")
    )
    stan_special_priors <- stan_special_prior_local(
      prior, class = "b", ncoef = length(fixef), 
      px = px, center_X = center_X, suffix = bsuffix
    )
    out <- collapse_lists(out, stan_special_priors)
  }
  
  order_intercepts <- order_intercepts(bterms)
  if (order_intercepts && !center_X) {
    stop2(
      "Identifying mixture components via ordering requires ",
      "population-level intercepts to be present.\n",
      "Try setting order = 'none' in function 'mixture'."
    )
  }
  if (center_X) {
    # centering the design matrix improves convergence
    sub_X_means <- ""
    if (length(fixef)) {
      sub_X_means <- glue(" - dot_product(means_X{p}, b{p})")
      if (is_ordinal(family)) {
        # the intercept was already removed during the data preparation
        str_add(out$tdata_def) <- glue(
          "  int Kc{p} = K{p};\n",
          "  matrix[N{resp}, Kc{p}] Xc{p};", 
          "  // centered version of X{p}\n",
          "  vector[Kc{p}] means_X{p};",
          "  // column means of X{p} before centering\n"
        )
        str_add(out$tdata_comp) <- glue(
          "  for (i in 1:K{p}) {{\n",
          "    means_X{p}[i] = mean(X{p}[, i]);\n",
          "    Xc{p}[, i] = X{p}[, i] - means_X{p}[i];\n",
          "  }}\n"
        )
      } else {
        str_add(out$tdata_def) <- glue(
          "  int Kc{p} = K{p} - 1;\n",
          "  matrix[N{resp}, Kc{p}] Xc{p};", 
          "  // centered version of X{p} without an intercept\n",
          "  vector[Kc{p}] means_X{p};",
          "  // column means of X{p} before centering\n"
        )
        str_add(out$tdata_comp) <- glue(
          "  for (i in 2:K{p}) {{\n",
          "    means_X{p}[i - 1] = mean(X{p}[, i]);\n",
          "    Xc{p}[, i - 1] = X{p}[, i] - means_X{p}[i - 1];\n",
          "  }}\n"
        )
      }
    }
    if (!is_ordinal(family)) {
      # intercepts of ordinal models are handled in 'stan_thres'
      if (order_intercepts) {
        # identify mixtures via ordering of the intercepts
        dp_id <- dpar_id(px$dpar)
        str_add(out$tpar_def) <- glue(
          "  // identify mixtures via ordering of the intercepts\n",                   
          "  real Intercept{p} = ordered_Intercept{resp}[{dp_id}];\n"
        )
      } else {
        str_add(out$par) <- glue(
          "  // temporary intercept for centered predictors\n",
          "  real Intercept{p};\n"
        )
      }
      str_add(out$eta) <- glue(" + Intercept{p}")
      str_add(out$gen_def) <- glue(
        "  // actual population-level intercept\n",
        "  real b{p}_Intercept = Intercept{p}{sub_X_means};\n"
      )
      str_add(out$prior) <- stan_prior(
        prior, class = "Intercept", px = px, suffix = p
      )
    }
  }
  if (decomp == "QR") {
    str_add(out$tdata_def) <- glue(
      "  // matrices for QR decomposition\n",
      "  matrix[N{resp}, K{ct}{p}] XQ{p};\n",
      "  matrix[K{ct}{p}, K{ct}{p}] XR{p};\n",
      "  matrix[K{ct}{p}, K{ct}{p}] XR{p}_inv;\n"
    )
    str_add(out$tdata_comp) <- glue(
      "  // compute and scale QR decomposition\n",
      "  XQ{p} = qr_thin_Q(X{ct}{p}) * sqrt(N{resp} - 1);\n",
      "  XR{p} = qr_thin_R(X{ct}{p}) / sqrt(N{resp} - 1);\n",
      "  XR{p}_inv = inverse(XR{p});\n"
    )
  }
  str_add(out$eta) <- stan_eta_fe(fixef, bterms)
  out
}

# Stan code for group-level effects
stan_re <- function(ranef, prior, ...) {
  IDs <- unique(ranef$id)
  out <- list()
  # special handling of student-t group effects
  tranef <- get_dist_groups(ranef, "student")
  if (has_rows(tranef)) {
    str_add(out$par) <- 
      "  // parameters for student-t distributed group-level effects\n"
    for (i in seq_rows(tranef)) {
      g <- usc(tranef$ggn[i])
      id <- tranef$id[i]
      str_add(out$par) <- glue(
        "  real<lower=1> df{g};\n",
        "  vector<lower=0>[N_{id}] udf{g};\n"
      )
      str_add(out$tpar_def) <- glue(
        "  vector[N_{id}] dfm{g} = sqrt(df{g} * udf{g});\n"
      )
      str_add(out$prior) <- stan_prior(
        prior, class = "df", group = tranef$group[i], suffix = g
      )
      str_add(out$prior) <- glue(
        "  target += inv_chi_square_lpdf(udf{g} | df{g});\n"
      )
    }
  }
  # the ID syntax requires group-level effects to be evaluated separately
  tmp <- lapply(IDs, .stan_re, ranef = ranef, prior = prior, ...)
  out <- collapse_lists(ls = c(list(out), tmp))
  out
}

# Stan code for group-level effects per ID
# @param id the ID of the grouping factor
# @param ranef output of tidy_ranef
# @param prior object of class brmsprior
# @param cov_ranef optional list of custom covariance matrices 
.stan_re <- function(id, ranef, prior, cov_ranef = NULL) {
  out <- list()
  r <- subset2(ranef, id = id)
  has_ccov <- r$group[1] %in% names(cov_ranef)
  has_by <- nzchar(r$by[[1]])
  Nby <- seq_along(r$bylevels[[1]]) 
  ng <- seq_along(r$gcall[[1]]$groups)
  px <- check_prefix(r)
  uresp <- usc(unique(px$resp))
  idp <- paste0(r$id, usc(combine_prefix(px)))
  # define data needed for group-level effects
  str_add(out$data) <- glue(
    "  // data for group-level effects of ID {id}\n",
    "  int<lower=1> N_{id};  // number of grouping levels\n",
    "  int<lower=1> M_{id};  // number of coefficients per level\n"
  )
  if (r$gtype[1] == "mm") {
    for (res in uresp) {
      str_add(out$data) <- cglue(
        "  int<lower=1> J_{id}{res}_{ng}[N{res}];",
        "  // grouping indicator per observation\n",
        "  real W_{id}{res}_{ng}[N{res}];",
        "  // multi-membership weights\n"
      )
    }
  } else {
    str_add(out$data) <- cglue(
      "  int<lower=1> J_{id}{uresp}[N{uresp}];",
      "  // grouping indicator per observation\n"
    )
  }
  if (has_by) {
    str_add(out$data) <- glue(
      "  int<lower=1> Nby_{id};  // number of by-factor levels\n",
      "  int<lower=1> Jby_{id}[N_{id}];", 
      "  // by-factor indicator per observation\n" 
    )
  }
  if (has_ccov) {
    str_add(out$data) <- glue(
      "  // cholesky factor of known covariance matrix\n",
      "  matrix[N_{id}, N_{id}] Lcov_{id};\n"
    )
  }
  J <- seq_rows(r)
  reqZ <- !r$type %in% "sp"
  if (any(reqZ)) {
    str_add(out$data) <- "  // group-level predictor values\n"
    if (r$gtype[1] == "mm") {
      for (i in which(reqZ)) {
        str_add(out$data) <- cglue(
          "  vector[N{usc(r$resp[i])}] Z_{idp[i]}_{r$cn[i]}_{ng};\n"
        )
      }
    } else {
      str_add(out$data) <- cglue(
        "  vector[N{usc(r$resp[reqZ])}] Z_{idp[reqZ]}_{r$cn[reqZ]};\n"
      )
    }
  }
  # define standard deviation parameters
  if (has_by) {
    str_add(out$par) <- glue(
      "  matrix<lower=0>[M_{id}, Nby_{id}] sd_{id};",
      "  // group-level standard deviations\n"
    )
    str_add(out$prior) <- stan_prior(
      prior, class = "sd", group = r$group[1], coef = r$coef,
      px = px, prefix = "to_vector(", suffix = glue("_{id})")
    )
  } else {
    str_add(out$par) <- glue(
      "  vector<lower=0>[M_{id}] sd_{id};",
      "  // group-level standard deviations\n"
    )
    str_add(out$prior) <- stan_prior(
      prior, class = "sd", group = r$group[1], coef = r$coef,
      px = px, suffix = glue("_{id}")
    )
  }
  dfm <- ""
  tr <- get_dist_groups(r, "student")
  if (nrow(r) > 1L && r$cor[1]) {
    # multiple correlated group-level effects
    str_add(out$data) <- glue( 
      "  int<lower=1> NC_{id};  // number of group-level correlations\n"
    )
    str_add(out$par) <- glue(
      "  matrix[M_{id}, N_{id}] z_{id};",
      "  // standardized group-level effects\n"
    )
    str_add(out$prior) <- glue(
      "  target += normal_lpdf(to_vector(z_{id}) | 0, 1);\n"
    )
    if (has_rows(tr)) {
      dfm <- glue("rep_matrix(dfm_{tr$ggn[1]}, M_{id}) .* ")
    }
    if (has_by) {
      if (has_ccov) {
        stop2(
          "Cannot combine 'by' variables with customized covariance ",
          "matrices when fitting multiple group-level effects."
        )
      }
      str_add(out$par) <- glue(
        "  // cholesky factor of correlation matrix\n",
        "  cholesky_factor_corr[M_{id}] L_{id}[Nby_{id}];\n"
      )
      str_add(out$tpar_def) <- glue(
        "  // actual group-level effects\n",
        "  matrix[N_{id}, M_{id}] r_{id}", 
        " = {dfm}scale_r_cor_by(z_{id}, sd_{id}, L_{id}, Jby_{id});\n"
      )
      str_add(out$prior) <- stan_prior(
        prior, class = "L", group = r$group[1],
        suffix = glue("_{id}[{Nby}]")
      )
      str_add(out$gen_def) <- cglue(
        "  // group-level correlations\n",
        "  corr_matrix[M_{id}] Cor_{id}_{Nby}",
        " = multiply_lower_tri_self_transpose(L_{id}[{Nby}]);\n",
        "  vector<lower=-1,upper=1>[NC_{id}] cor_{id}_{Nby};\n"
      )
      str_add(out$gen_comp) <- stan_cor_gen_comp(
        glue("cor_{id}_{Nby}"), glue("M_{id}")
      )
    } else {
      str_add(out$par) <- glue(
        "  // cholesky factor of correlation matrix\n",
        "  cholesky_factor_corr[M_{id}] L_{id};\n"
      )
      if (has_ccov) {
        rdef <- glue(
          "as_matrix(kronecker(Lcov_{id},", 
          " diag_pre_multiply(sd_{id}, L_{id})) *",
          " to_vector(z_{id}), N_{id}, M_{id});\n"
        )
      } else {
        rdef <- glue(
          "(diag_pre_multiply(sd_{id}, L_{id}) * z_{id})';\n"
        )
      }
      str_add(out$tpar_def) <- glue(
        "  // actual group-level effects\n",
        "  matrix[N_{id}, M_{id}] r_{id} = {dfm}{rdef}"
      )
      str_add(out$prior) <- stan_prior(
        prior, class = "L", group = r$group[1], suffix = usc(id)
      )
      str_add(out$gen_def) <- glue(
        "  // group-level correlations\n",
        "  corr_matrix[M_{id}] Cor_{id}",
        " = multiply_lower_tri_self_transpose(L_{id});\n",
        "  vector<lower=-1,upper=1>[NC_{id}] cor_{id};\n"
      )
      str_add(out$gen_comp) <- stan_cor_gen_comp(
        cor = glue("cor_{id}"), ncol = glue("M_{id}")
      )
    }
    str_add(out$tpar_def) <- 
      "  // using vectors speeds up indexing in loops\n"
    str_add(out$tpar_def) <- cglue(
      "  vector[N_{id}] r_{idp}_{r$cn} = r_{id}[, {J}];\n"
    )
  } else {
    # single or uncorrelated group-level effects
    str_add(out$par) <- glue(
      "  // standardized group-level effects\n",
      "  vector[N_{id}] z_{id}[M_{id}];\n"
    )
    str_add(out$prior) <- cglue(
      "  target += normal_lpdf(z_{id}[{seq_rows(r)}] | 0, 1);\n"
    )
    str_add(out$tpar_def) <- "  // actual group-level effects\n" 
    Lcov <- str_if(has_ccov, glue("Lcov_{id} * "))
    if (has_rows(tr)) {
      dfm <- glue("dfm_{tr$ggn[1]} .* ")
    }
    if (has_by) {
      str_add(out$tpar_def) <- cglue(
        "  vector[N_{id}] r_{idp}_{r$cn}",
        " = {dfm}(sd_{id}[{J}, Jby_{id}]' .* ({Lcov}z_{id}[{J}]));\n"
      )
    } else {
      str_add(out$tpar_def) <- cglue(
        "  vector[N_{id}] r_{idp}_{r$cn}",
        " = {dfm}(sd_{id}[{J}] * ({Lcov}z_{id}[{J}]));\n"
      )
    }
  }
  out
}

# Stan code of smooth terms
stan_sm <- function(bterms, data, prior, ...) {
  out <- list()
  smef <- tidy_smef(bterms, data)
  if (!NROW(smef)) {
    return(out)
  }
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  resp <- usc(px$resp)
  Xs_names <- attr(smef, "Xs_names")
  if (length(Xs_names)) {
    str_add(out$data) <- glue(
      "  // data for splines\n",
      "  int Ks{p};  // number of linear effects\n",
      "  matrix[N{resp}, Ks{p}] Xs{p};",
      "  // design matrix for the linear effects\n"
    )
    str_add(out$pars) <- glue(
      "  // spline coefficients\n",
      "  vector[Ks{p}] bs{p};\n"
    )
    str_add(out$eta) <- glue(" + Xs{p} * bs{p}")
    str_add(out$prior) <- stan_prior(
      prior, class = "b", coef = Xs_names,
      px = px, suffix = glue("s{p}")
    )
  }
  for (i in seq_rows(smef)) {
    pi <- glue("{p}_{i}")
    nb <- seq_len(smef$nbases[[i]])
    str_add(out$data) <- glue(
      "  // data for spline {smef$byterm[i]}\n",  
      "  int nb{pi};  // number of bases\n",
      "  int knots{pi}[nb{pi}];  // number of knots\n"
    )
    str_add(out$data) <- "  // basis function matrices\n"
    str_add(out$data) <- cglue(
      "  matrix[N{resp}, knots{pi}[{nb}]] Zs{pi}_{nb};\n"
    )
    str_add(out$par) <- glue(
      "  // parameters for spline {smef$byterm[i]}\n"
    )
    str_add(out$par) <- cglue(
      "  // standarized spline coefficients\n",
      "  vector[knots{pi}[{nb}]] zs{pi}_{nb};\n",
      "  // standard deviations of the coefficients\n",
      "  real<lower=0> sds{pi}_{nb};\n"
    )
    str_add(out$tpar_def) <- cglue(
      "  // actual spline coefficients\n",
      "  vector[knots{pi}[{nb}]] s{pi}_{nb}", 
      " = sds{pi}_{nb} * zs{pi}_{nb};\n"
    )
    str_add(out$prior) <- cglue(
      "  target += normal_lpdf(zs{pi}_{nb} | 0, 1);\n"
    )
    str_add(out$prior) <- stan_prior(
      prior, class = "sds", coef = smef$term[i], 
      px = px, suffix = glue("{pi}_{nb}")
    )
    str_add(out$eta) <- cglue(
      " + Zs{pi}_{nb} * s{pi}_{nb}"
    )
  }
  out
}

# Stan code for category specific effects
# @note not implemented for non-linear models
stan_cs <- function(bterms, data, prior, ranef, ...) {
  out <- list()
  csef <- colnames(get_model_matrix(bterms$cs, data))
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  resp <- usc(bterms$resp)
  ranef <- subset2(ranef, type = "cs", ls = px)
  if (length(csef)) {
    str_add(out$data) <- glue(
      "  int<lower=1> Kcs{p};  // number of category specific effects\n",
      "  matrix[N{resp}, Kcs{p}] Xcs{p};  // category specific design matrix\n"
    )
    bound <- get_bound(prior, class = "b", px = px)
    str_add(out$par) <- glue(
      "  matrix{bound}[Kcs{p}, nthres{resp}] bcs{p};",
      "  // category specific effects\n"
    )
    str_add(out$model_def) <- glue(
      "  // linear predictor for category specific effects\n",
      "  matrix[N{resp}, nthres{resp}] mucs{p} = Xcs{p} * bcs{p};\n"
    ) 
    str_add(out$prior) <- stan_prior(
      prior, class = "b", coef = csef,
      suffix = "cs", px = px, matrix = TRUE
    )
  }
  if (nrow(ranef)) {
    if (!length(csef)) {
      # only group-level category specific effects present
      str_add(out$model_def) <- glue(
        "  // linear predictor for category specific effects\n",               
        "  matrix[N{resp}, nthres{resp}] mucs{p}", 
        " = rep_matrix(0, N{resp}, nthres{resp});\n"
      )
    }
    thres_regex <- "(?<=\\[)[[:digit:]]+(?=\\]$)"
    thres <- get_matches(thres_regex, ranef$coef, perl = TRUE)
    nthres <- max(as.numeric(thres))
    mucs_loop <- ""
    for (i in seq_len(nthres)) {
      r_cat <- ranef[grepl(glue("\\[{i}\\]$"), ranef$coef), ]
      str_add(mucs_loop) <- glue(
        "    mucs{p}[n, {i}] = mucs{p}[n, {i}]"
      )
      for (id in unique(r_cat$id)) {
        r <- r_cat[r_cat$id == id, ]
        rpx <- check_prefix(r)
        idp <- paste0(r$id, usc(combine_prefix(rpx)))
        idresp <- paste0(r$id, usc(rpx$resp))
        str_add(mucs_loop) <- cglue(
          " + r_{idp}_{r$cn}[J_{idresp}[n]] * Z_{idp}_{r$cn}[n]"
        )
      }
      str_add(mucs_loop) <- ";\n"
    }
    str_add(out$model_comp_eta_loop) <- glue(
      "  for (n in 1:N{resp}) {{\n    ", mucs_loop, "  }\n"
    )
   
  }
  out
}

# Stan code for special effects
stan_sp <- function(bterms, data, prior, stanvars, ranef, meef, ...) {
  out <- list()
  spef <- tidy_spef(bterms, data)
  if (!nrow(spef)) return(out)
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  resp <- usc(px$resp)
  ranef <- subset2(ranef, type = "sp", ls = px)
  spef_coef <- rename(spef$term)
  invalid_coef <- setdiff(ranef$coef, spef_coef)
  if (length(invalid_coef)) {
    stop2(
      "Special group-level terms require corresponding ", 
      "population-level terms:\nOccured for ", 
      collapse_comma(invalid_coef)
    )
  }
  # prepare Stan code of the linear predictor component
  for (i in seq_rows(spef)) {
    eta <- spef$joint_call[[i]]
    if (!is.null(spef$calls_mo[[i]])) {
      new_mo <- glue("mo(simo{p}_{spef$Imo[[i]]}, Xmo{p}_{spef$Imo[[i]]}[n])")
      eta <- rename(eta, spef$calls_mo[[i]], new_mo)
    }
    if (!is.null(spef$calls_me[[i]])) {
      Kme <- seq_along(meef$term)
      Ime <- match(meef$grname, unique(meef$grname))
      nme <- ifelse(nzchar(meef$grname), glue("Jme_{Ime}[n]"), "n")
      new_me <- glue("Xme_{Kme}[{nme}]")
      eta <- rename(eta, meef$term, new_me)
    }
    if (!is.null(spef$calls_mi[[i]])) {
      new_mi <- glue("Yl_{spef$vars_mi[[i]]}[n]")
      eta <- rename(eta, spef$calls_mi[[i]], new_mi)
    }
    if (spef$Ic[i] > 0) {
      str_add(eta) <- glue(" * Csp{p}_{spef$Ic[i]}[n]")
    }
    r <- subset2(ranef, coef = spef_coef[i])
    rpars <- str_if(nrow(r), glue(" + {stan_eta_rsp(r)}"))
    str_add(out$loopeta) <- glue(" + (bsp{p}[{i}]{rpars}) * {eta}")
  }
  # prepare general Stan code
  ncovars <- max(spef$Ic)
  str_add(out$data) <- glue(
    "  int<lower=1> Ksp{p};  // number of special effects terms\n"
  )
  if (ncovars > 0L) {
    str_add(out$data) <- glue(
      "  // covariates of special effects terms\n",
      cglue("  vector[N{resp}] Csp{p}_{seq_len(ncovars)};\n")
    )
  }
  # prepare special effects coefficients
  if (!glue("bsp{p}") %in% names(stanvars)) {
    bound <- get_bound(prior, class = "b", px = px)
    if (stan_use_horseshoe(bterms, prior)) {
      str_add(out$tpar_def) <- glue(
        "  // special effects coefficients\n", 
        "  vector{bound}[Ksp{p}] bsp{p};\n"
      )
    } else {
      str_add(out$par) <- glue(
        "  // special effects coefficients\n", 
        "  vector{bound}[Ksp{p}] bsp{p};\n"
      )
    }
  }
  str_add(out$prior) <- stan_prior(
    prior, class = "b", coef = spef$coef, 
    px = px, suffix = glue("sp{p}")
  )
  # include special Stan code for monotonic effects
  I <- unlist(spef$Imo)
  if (length(I)) {
    I <- seq_len(max(I))
    str_add(out$data) <- glue(
      "  int<lower=1> Imo{p};  // number of monotonic variables\n",
      "  int<lower=1> Jmo{p}[Imo{p}];  // length of simplexes\n",
      "  // monotonic variables\n",
      cglue("  int Xmo{p}_{I}[N{resp}];\n"),
      "  // prior concentration of monotonic simplexes\n",
      cglue("  vector[Jmo{p}[{I}]] con_simo{p}_{I};\n")
    )
    str_add(out$par) <- glue(
      "  // simplexes of monotonic effects\n",
      cglue("  simplex[Jmo{p}[{I}]] simo{p}_{I};\n")
    ) 
    str_add(out$prior) <- cglue(
      "  target += dirichlet_lpdf(simo{p}_{I} | con_simo{p}_{I});\n"
    )
  }
  stan_special_priors <- stan_special_prior_local(
    prior, class = "bsp", ncoef = nrow(spef), 
    px = px, center_X = FALSE
  )  
  out <- collapse_lists(out, stan_special_priors)
  out
}

# Stan code for latent gaussian processes
stan_gp <- function(bterms, data, prior, ...) {
  out <- list()
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  resp <- usc(px$resp)
  gpef <- tidy_gpef(bterms, data)
  for (i in seq_rows(gpef)) {
    pi <- glue("{p}_{i}")
    byvar <- gpef$byvars[[i]] 
    cons <- gpef$cons[[i]]
    byfac <- length(cons) > 0L
    bynum <- !is.null(byvar) && !byfac 
    k <- gpef$k[i]
    iso <- gpef$iso[i]
    gr <- gpef$gr[i]
    sfx1 <- gpef$sfx1[[i]]
    sfx2 <- gpef$sfx2[[i]]
    str_add(out$data) <- glue(
      "  // data related to GPs\n",
      "  // number of sub-GPs (equal to 1 unless 'by' was used)\n",
      "  int<lower=1> Kgp{pi};\n",
      "  int<lower=1> Dgp{pi};  // GP dimension\n"
    )
    if (!isNA(k)) {
      # !isNA(k) indicates the use of approximate GPs
      str_add(out$data) <- glue(
        "  // number of basis functions of an approximate GP\n",
        "  int<lower=1> NBgp{pi};\n"
      )
    } 
    str_add(out$par) <- glue(
      "  // GP standard deviation parameters\n",
      "  vector<lower=0>[Kgp{pi}] sdgp{pi};\n"
    )
    if (gpef$iso[i]) {
      str_add(out$par) <- glue(
        "  // GP length-scale parameters\n",
        "  vector<lower=0>[1] lscale{pi}[Kgp{pi}];\n" 
      )
    } else {
      str_add(out$par) <- glue(
        "  // GP length-scale parameters\n",
        "  vector<lower=0>[Dgp{pi}] lscale{pi}[Kgp{pi}];\n"
      )
    }
    str_add(out$prior) <- stan_prior(
      prior, class = "sdgp", coef = sfx1, 
      px = px, suffix = pi
    )
    if (byfac) {
      J <- seq_along(cons)
      Ngp <- glue("Ngp{pi}")
      Nsubgp <- glue("N", str_if(gr, "sub"), glue("gp{pi}"))
      Igp <- glue("Igp{pi}_{J}")
      str_add(out$data) <- glue(
        "  // number of observations relevant for a certain sub-GP\n",
        "  int<lower=1> {Ngp}[Kgp{pi}];\n"
      )
      str_add(out$data) <- 
        "  // indices and contrasts of sub-GPs per observation\n"
      str_add(out$data) <- cglue(
        "  int<lower=1> {Igp}[{Ngp}[{J}]];\n",
        "  vector[{Ngp}[{J}]] Cgp{pi}_{J};\n"
      )
      if (gr) {
        str_add(out$data) <- glue(
          "  // number of latent GP groups\n",
          "  int<lower=1> Nsubgp{pi}[Kgp{pi}];\n"
        )
        str_add(out$data) <- cglue(
          "  // indices of latent GP groups per observation\n",
          "  int<lower=1> Jgp{pi}_{J}[{Ngp}[{J}]];\n"
        )
      }
      gp_call <- glue("Cgp{pi}_{J} .* ")
      if (!isNA(k)) {
        str_add(out$data) <- 
          "  // approximate GP basis matrices and eigenvalues\n"
        str_add(out$data) <- cglue(
          "  matrix[{Nsubgp}[{J}], NBgp{pi}] Xgp{pi}_{J};\n",
          "  vector[Dgp{pi}] slambda{pi}_{J}[NBgp{pi}];\n"
        )
        str_add(out$par) <- "  // latent variables of the GP\n"
        str_add(out$par) <- cglue(
          "  vector[NBgp{pi}] zgp{pi}_{J};\n"
        )
        str_add(gp_call) <- glue(
          "gpa(Xgp{pi}_{J}, sdgp{pi}[{J}], ", 
          "lscale{pi}[{J}], zgp{pi}_{J}, slambda{pi}_{J})"
        )
      } else {
        str_add(out$data) <- "  // covariates of the GP\n"
        str_add(out$data) <- cglue(
          "  vector[Dgp{pi}] Xgp{pi}_{J}[{Nsubgp}[{J}]];\n"
        )
        str_add(out$par) <- "  // latent variables of the GP\n"
        str_add(out$par) <- cglue(
          "  vector[{Nsubgp}[{J}]] zgp{pi}_{J};\n"
        )
        str_add(gp_call) <- glue(
          "gp(Xgp{pi}_{J}, sdgp{pi}[{J}], ", 
          "lscale{pi}[{J}], zgp{pi}_{J})"
        )
      }
      Jgp <- str_if(gr, glue("[Jgp{pi}_{J}]"))
      eta <- combine_prefix(px, keep_mu = TRUE, nlp = TRUE)
      eta <- glue("{eta}[{Igp}]")
      # compound '+=' statement currently causes a parser failure
      str_add(out$model_comp_basic) <- glue(
        "  {eta} = {eta} + {gp_call}{Jgp};\n"
      )
      str_add(out$prior) <- cglue(
        "{tp()}normal_lpdf(zgp{pi}_{J} | 0, 1);\n"
      )
      for (j in seq_along(cons)) {
        str_add(out$prior) <- stan_prior(
          prior, class = "lscale", coef = sfx2[j, ], 
          px = px, suffix = glue("{pi}[{j}]")
        )
      }
    } else {
      # no by-factor variable
      Nsubgp <- str_if(gr, glue("Nsubgp{pi}"), glue("N{resp}"))
      if (gr) {
        str_add(out$data) <- glue(
          "  // number of latent GP groups\n",
          "  int<lower=1> {Nsubgp};\n",
          "  // indices of latent GP groups per observation\n",
          "  int<lower=1> Jgp{pi}[N{resp}];\n"
        )
      }
      if (!isNA(k)) {
        str_add(out$data) <- glue(
          "  // approximate GP basis matrices\n",
          "  matrix[{Nsubgp}, NBgp{pi}] Xgp{pi};\n",
          "  // approximate GP eigenvalues\n",
          "  vector[Dgp{pi}] slambda{pi}[NBgp{pi}];\n"
        )
        str_add(out$par) <- glue(
          "  // latent variables of the GP\n",
          "  vector[NBgp{pi}] zgp{pi};\n"
        )
        gp_call <- glue(
          "gpa(Xgp{pi}, sdgp{pi}[1], lscale{pi}[1], zgp{pi}, slambda{pi})"
        )
      } else {
        str_add(out$data) <- glue(
          "  // covariates of the GP\n",
          "  vector[Dgp{pi}] Xgp{pi}[{Nsubgp}];\n"
        ) 
        str_add(out$par) <- glue(
          "  // latent variables of the GP\n",
          "  vector[{Nsubgp}] zgp{pi};\n"
        )
        gp_call <- glue(
          "gp(Xgp{pi}, sdgp{pi}[1], lscale{pi}[1], zgp{pi})"
        )
      }
      if (bynum) {
        str_add(out$data) <- glue(
          "  // numeric by-variable of the GP\n",
          "  vector[N{resp}] Cgp{pi};\n"
        )
      }
      Cgp <- str_if(bynum, glue("Cgp{pi} .* "))
      Jgp <- str_if(gr, glue("[Jgp{pi}]"))
      str_add(out$eta) <- glue(" + {Cgp}{gp_call}{Jgp}")
      str_add(out$prior) <- glue(
        "{tp()}normal_lpdf(zgp{pi} | 0, 1);\n"
      )
      str_add(out$prior) <- stan_prior(
        prior, class = "lscale", coef = sfx2, 
        px = px, suffix = glue("{pi}[1]")
      )
    }
  }
  out
}

# Stan code for the linear predictor of autocorrelation terms 
stan_ac <- function(bterms, data, prior, ...) {
  # stopifnot(is.btl(bterms))
  out <- list()
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  has_natural_residuals <- has_natural_residuals(bterms)
  has_latent_residuals <- has_latent_residuals(bterms)
  acef <- tidy_acef(bterms, data)
  
  acef_arma <- subset2(acef, class = "arma")
  if (NROW(acef_arma)) {
    # TODO: check if autocor is specified for mu only
    err_msg <- "ARMA models are not implemented"
    if (is.mixfamily(bterms$family)) {
      stop2(err_msg, " for mixture models.") 
    }
    str_add(out$data) <- glue( 
      "  // data needed for ARMA correlations\n",
      "  int<lower=0> Kar{p};  // AR order\n",
      "  int<lower=0> Kma{p};  // MA order\n"
    )
    str_add(out$tdata_def) <- glue( 
      "  int max_lag{p} = max(Kar{p}, Kma{p});\n"
    )
    ar_bound <- ma_bound <- "<lower=-1,upper=1>"
    if (!acef_arma$cov) {
      err_msg <- "Please set cov = TRUE in ARMA correlation structures"
      if (!has_natural_residuals) {
        stop2(err_msg, " for family '", bterms$family$family, "'.")
      }
      if (is.formula(bterms$adforms$se)) {
        stop2(err_msg, " when including known standard errors.")
      }
      str_add(out$data) <- glue( 
        "  // number of lags per observation\n",
        "  int<lower=0> J_lag{p}[N{p}];\n"                
      )
      str_add(out$model_def) <- glue(
        "  // objects storing residuals\n",
        "  matrix[N{p}, max_lag{p}] Err{p}",
        " = rep_matrix(0, N{p}, max_lag{p});\n",
        "  vector[N{p}] err{p};\n"
      )
      Y <- str_if(is.formula(bterms$adforms$mi), "Yl", "Y")
      add_ar <- str_if(acef_arma$p > 0,
        glue("    mu{p}[n] += Err{p}[n, 1:Kar{p}] * ar{p};\n")             
      )
      add_ma <- str_if(acef_arma$q > 0,
        glue("    mu{p}[n] += Err{p}[n, 1:Kma{p}] * ma{p};\n")             
      )
      str_add(out$model_comp_arma) <- glue(
        "  // include ARMA terms\n",
        "  for (n in 1:N{p}) {{\n",
        add_ma,
        "    err{p}[n] = {Y}{p}[n] - mu{p}[n];\n",
        "    for (i in 1:J_lag{p}[n]) {{\n",
        "      Err{p}[n + 1, i] = err{p}[n + 1 - i];\n",
        "    }}\n",
        add_ar,
        "  }}\n"
      )
      # in the conditional formulation no boundaries are required
      ar_bound <- ma_bound <- ""
    }
    if (acef_arma$p > 0) {
      str_add(out$par) <- glue( 
        "  vector{ar_bound}[Kar{p}] ar{p};  // autoregressive effects\n"
      )
      str_add(out$prior) <- stan_prior(
        prior, class = "ar", px = px, suffix = p
      )
    }
    if (acef_arma$q > 0) {
      str_add(out$par) <- glue( 
        "  vector{ma_bound}[Kma{p}] ma{p};  // moving-average effects\n"
      )
      str_add(out$prior) <- stan_prior(
        prior, class = "ma", px = px, suffix = p
      )
    }
  }
  
  acef_cosy <- subset2(acef, class = "cosy")
  if (NROW(acef_cosy)) {
    # compound symmetry correlation structure
    err_msg <- "Compound symmetry models are not implemented"
    if (is.mixfamily(bterms$family)) {
      stop2(err_msg, " for mixture models.")
    }
    # most code is shared with ARMA covariance models
    str_add(out$par) <- glue(
      "  real<lower=0,upper=1> cosy{p};  // compound symmetry correlation\n"
    )
    str_add(out$prior) <- stan_prior(
      prior, class = "cosy", px = px, suffix = p
    )
  }
  
  acef_time_cov <- subset2(acef, dim = "time", cov = TRUE)
  if (NROW(acef_time_cov)) {
    # use correlation structures in covariance matrix parameterization
    # optional for ARMA models and obligatory for compound symmetry models
    # TODO: move these checks to a place where all information is available
    # err_msg <- "Cannot model residual covariance matrices via 'autocor'"
    # if (isTRUE(bterms$rescor)) {
    #   stop2(err_msg, " when estimating 'rescor'.")
    # }
    # pred_other_pars <- any(!names(bterms$dpars) %in% c("mu", "sigma"))
    # if (has_natural_residuals && pred_other_pars) {
    #   stop2(err_msg, " when predicting parameters other than 'mu' and 'sigma'.")
    # }
    # can only model one covariance structure at a time
    stopifnot(NROW(acef_time_cov) == 1)
    str_add(out$data) <- glue( 
      "  // see the functions block for details\n",
      "  int<lower=1> N_tg{p};\n",
      "  int<lower=1> begin_tg{p}[N_tg{p}];\n", 
      "  int<lower=1> end_tg{p}[N_tg{p}];\n", 
      "  int<lower=1> nobs_tg{p}[N_tg{p}];\n"
    )
    if (!is.formula(bterms$adforms$se)) {
      str_add(out$tdata_def) <- glue(
        "  // no known standard errors specified by the user\n",
        "  vector[N{p}] se2{p} = rep_vector(0, N{p});\n"
      )
    }
    str_add(out$tpar_def) <- glue(
      "  // cholesky factor of the autocorrelation matrix\n",
      "  matrix[max(nobs_tg{p}), max(nobs_tg{p})] chol_cor{p};\n"               
    )
    if (acef_time_cov$class == "arma") {
      if (acef_time_cov$p > 0 && acef_time_cov$q == 0) {
        cor_fun <- "ar1"
        cor_args <- glue("ar{p}[1]")
      } else if (acef_time_cov$p == 0 && acef_time_cov$q > 0) {
        cor_fun <- "ma1"
        cor_args <- glue("ma{p}[1]")
      } else {
        cor_fun <- "arma1"
        cor_args <- glue("ar{p}[1], ma{p}[1]")
      }
    } else if (acef_time_cov$class == "cosy") {
      cor_fun <- "cosy"
      cor_args <- glue("cosy{p}")
    }
    str_add(out$tpar_comp) <- glue(
      "  // compute residual covariance matrix\n",
      "  chol_cor{p} = cholesky_cor_{cor_fun}({cor_args}, max(nobs_tg{p}));\n"
    )
    if (has_latent_residuals) {
      # TODO: support autocor terms for non-linear models via 'autocor'
      # err_msg <- "Latent residuals are not implemented"
      # if (is.btnl(bterms$dpars[["mu"]])) {
      #   stop2(err_msg, " for non-linear models.")
      # }
      # if (conv_cats_dpars(bterms)) {
      #   stop2(err_msg, " for this family.")
      # }
      str_add(out$par) <- glue(
        "  vector[N{p}] zerr{p};  // unscaled residuals\n",
        "  real<lower=0> sderr{p};  // SD of residuals\n"
      )
      str_add(out$tpar_def) <- glue(
        "  vector[N{p}] err{p};  // actual residuals\n"
      )
      str_add(out$tpar_comp) <- glue(
        "  // compute correlated residuals\n",
        "  err{p} = scale_cov_err(",
        "zerr{p}, sderr{p}, chol_cor{p}, nobs_tg{p}, begin_tg{p}, end_tg{p});\n"
      )
      str_add(out$prior) <- glue(
        "  target += normal_lpdf(zerr{p} | 0, 1);\n"
      )
      str_add(out$prior) <- stan_prior(
        prior, class = "sderr", px = px, suffix = p
      )
      str_add(out$eta) <- glue(" + err{p}")
    }
  }
  
  acef_sar <- subset2(acef, class = "sar")
  if (NROW(acef_sar)) {
    err_msg <- "SAR models are not implemented"
    if (is.mixfamily(bterms$family)) {
      stop2(err_msg, " for mixture models.") 
    }
    if (!has_natural_residuals) {
      stop2(err_msg, " for family '", bterms$family$family, "'.")
    }
    # if (isTRUE(bterms$rescor)) {
    #   stop2(err_msg, " when 'rescor' is estimated.")
    # }
    # if (any(c("sigma", "nu") %in% names(bterms$dpars))) {
    #   stop2(err_msg, " when predicting 'sigma' or 'nu'.")
    # }
    str_add(out$data) <- glue(
      "  matrix[N{p}, N{p}] W{p};  // spatial weight matrix\n",
      "  vector[N{p}] eigenW{p};  // eigenvalues of W{p}\n"
    )
    str_add(out$tdata_def) <- glue(
      "  // the eigenvalues define the boundaries of the SAR correlation\n",
      "  real min_eigenW{p} = min(eigenW{p});\n",
      "  real max_eigenW{p} = max(eigenW{p});\n"
    )
    if (identical(acef_sar$type, "lag")) {
      str_add(out$par) <- glue( 
        "  // lag-SAR correlation parameter\n",
        "  real<lower=min_eigenW{p},upper=max_eigenW{p}> lagsar{p};\n"
      )
      str_add(out$prior) <- stan_prior(
        prior, class = "lagsar", px = px, suffix = p
      )
    } else if (identical(acef_sar$type, "error")) {
      str_add(out$par) <- glue( 
        "  // error-SAR correlation parameter\n",
        "  real<lower=min_eigenW{p},upper=max_eigenW{p}> errorsar{p};\n"
      )
      str_add(out$prior) <- stan_prior(
        prior, class = "errorsar", px = px, suffix = p
      )
    }
  }
  
  acef_car <- subset2(acef, class = "car")
  if (NROW(acef_car)) {
    err_msg <- "CAR models are not implemented"
    if (is.mixfamily(bterms$family)) {
      stop2(err_msg, " for mixture models.") 
    }
    if (is.btnl(bterms)) {
      stop2(err_msg, " for non-linear models.")
    }
    str_add(out$data) <- glue(
      "  // data for the CAR structure\n",
      "  int<lower=1> Nloc{p};\n",
      "  int<lower=1> Jloc{p}[N{p}];\n",
      "  int<lower=0> Nedges{p};\n",
      "  int<lower=1> edges1{p}[Nedges{p}];\n",
      "  int<lower=1> edges2{p}[Nedges{p}];\n"
    )
    str_add(out$par) <- glue(
      "  real<lower=0> sdcar{p};  // SD of the CAR structure\n"
    )
    str_add(out$prior) <- stan_prior(
      prior, class = "sdcar", px = px, suffix = p
    )
    str_add(out$loopeta) <- glue(" + rcar{p}[Jloc{p}[n]]")
    if (acef_car$type %in% c("escar", "esicar")) {
      str_add(out$data) <- glue(
        "  vector[Nloc{p}] Nneigh{p};\n",
        "  vector[Nloc{p}] eigenW{p};\n"
      )
    }
    if (acef_car$type %in% "escar") {
      str_add(out$par) <- glue(
        "  real<lower=0, upper=1> car{p};\n",
        "  vector[Nloc{p}] rcar{p};\n"
      )
      car_args <- c(
        "car", "sdcar", "Nloc", "Nedges", 
        "Nneigh", "eigenW", "edges1", "edges2"
      )
      car_args <- paste0(car_args, p, collapse = ", ")
      str_add(out$prior) <- stan_prior(
        prior, class = "car", px = px, suffix = p
      )
      str_add(out$prior) <- glue(
        "  target += sparse_car_lpdf(\n", 
        "    rcar{p} | {car_args}\n",
        "  );\n"
      )
    } else if (acef_car$type %in% "esicar") {
      str_add(out$par) <- glue(
        "  vector[Nloc{p} - 1] zcar{p};\n"
      )
      str_add(out$tpar_def) <- glue(
        "  vector[Nloc{p}] rcar{p};\n"                
      )
      str_add(out$tpar_comp) <- glue(
        "  // sum-to-zero constraint\n",
        "  rcar[1:(Nloc{p} - 1)] = zcar{p};\n",
        "  rcar[Nloc{p}] = - sum(zcar{p});\n"
      )
      car_args <- c(
        "sdcar", "Nloc", "Nedges", 
        "Nneigh", "eigenW", "edges1", "edges2"
      )
      car_args <- paste0(car_args, p, collapse = ", ")
      str_add(out$prior) <- glue(
        "  target += sparse_icar_lpdf(\n", 
        "    rcar{p} | {car_args}\n",
        "  );\n"
      )
    } else if (acef_car$type %in% "icar") {
      # intrinsic car based on the case study of Mitzi Morris
      # http://mc-stan.org/users/documentation/case-studies/icar_stan.html
      str_add(out$par) <- glue(
        "  // parameters for the ICAR structure\n",
        "  vector[Nloc{p}] zcar{p};\n"
      )
      str_add(out$tpar_def) <- glue(
        "  // scaled parameters for the ICAR structure\n",
        "  vector[Nloc{p}] rcar{p} = zcar{p} * sdcar{p};\n"
      )
      str_add(out$prior) <- glue(
        "  // improper prior on the spatial CAR component\n",
        "  target += -0.5 * dot_self(zcar{p}[edges1{p}] - zcar{p}[edges2{p}]);\n",
        "  // soft sum-to-zero constraint\n",
        "  target += normal_lpdf(sum(zcar{p}) | 0, 0.001 * Nloc{p});\n"
      )
    } else if (acef_car$type %in% "bym2") {
      # BYM2 car based on the case study of Mitzi Morris
      # http://mc-stan.org/users/documentation/case-studies/icar_stan.html
      str_add(out$data) <- glue(
        "  // scaling factor of the spatial CAR component\n",
        "  real<lower=0> car_scale{p};\n"
      )
      str_add(out$par) <- glue(
        "  // parameters for the BYM2 structure\n",
        "  vector[Nloc{p}] zcar{p};  // spatial part\n",
        "  vector[Nloc{p}] nszcar{p};  // non-spatial part\n",
        "  // proportion of variance in the spatial part\n",
        "  real<lower=0,upper=1> rhocar{p};\n"
      )
      str_add(out$tpar_def) <- glue(
        "  // join the spatial and the non-spatial CAR component\n",
        "  vector[Nloc{p}] rcar{p} = (sqrt(1 - rhocar{p}) * nszcar{p}", 
        " + sqrt(rhocar{p} * inv(car_scale{p})) * zcar{p}) * sdcar{p};\n"
      )
      str_add(out$prior) <- stan_prior(
        prior, class = "rhocar", px = px, suffix = p
      )
      str_add(out$prior) <- glue(
        "  // improper prior on the spatial BYM2 component\n",
        "  target += -0.5 * dot_self(zcar{p}[edges1{p}] - zcar{p}[edges2{p}]);\n",
        "  // soft sum-to-zero constraint\n",
        "  target += normal_lpdf(sum(zcar{p}) | 0, 0.001 * Nloc{p});\n",
        "  // proper prior on the non-spatial BYM2 component\n",
        "  target += normal_lpdf(nszcar | 0, 1);\n"
      )
    }
  }
  
  acef_fcor <- subset2(acef, class = "fcor")
  if (NROW(acef_fcor)) {
    err_msg <- "Correlation structure 'fcor' is not implemented"
    if (is.mixfamily(bterms$family)) {
      stop2(err_msg, " for mixture models.") 
    }
    if (!has_natural_residuals) {
      stop2(err_msg, " for family '", bterms$family$family, "'.")
    }
    if (isTRUE(bterms$rescor)) {
      stop2(err_msg, " when 'rescor' is estimated.")
    }
    str_add(out$data) <- glue( 
      "  matrix[N{p}, N{p}] V{p};  // known residual covariance matrix\n"
    )
    if (bterms$family$family %in% "gaussian") {
      str_add(out$tdata_def) <- glue(
        "  matrix[N{p}, N{p}] LV{p} = cholesky_decompose(V{p});\n"
      )
    }
  }
  out
}

# stan code for offsets
stan_offset <- function(bterms, ...) {
  out <- list()
  if (is.formula(bterms$offset)) {
    p <- usc(combine_prefix(bterms))
    resp <- usc(bterms$resp)
    str_add(out$data) <- glue( "  vector[N{resp}] offset{p};\n")
    str_add(out$eta) <- glue(" + offset{p}")
  }
  out
}

# add the denominator of a rate response to the Stan predictor term
# @param loop is the denominator added within a loop over observations?
stan_rate <- function(bterms, loop = FALSE, ...) {
  loop <- as_one_logical(loop)
  out <- list()
  if (is.formula(bterms$adforms$rate)) {
    # TODO: support other link functions as well?
    if (bterms$family$link != "log") {
      stop2("The 'rate' addition term requires a log-link.")
    }
    resp <- usc(bterms$resp)
    n <- str_if(loop, "[n]")
    str_add(out$eta) <- glue(" + log_denom{resp}{n}")
  }
  out
}

# Stan code related to autocorrelation structures
stan_autocor <- function(bterms, prior) {
  stopifnot(is.brmsterms(bterms))
  family <- bterms$family
  autocor <- bterms$autocor
  resp <- bterms$response
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  has_natural_residuals <- has_natural_residuals(bterms)
  has_latent_residuals <- has_latent_residuals(bterms)
  out <- list()
  if (is.cor_arma(autocor)) {
    err_msg <- "ARMA models are not implemented"
    if (is.mixfamily(family)) {
      stop2(err_msg, " for mixture models.") 
    }
    str_add(out$data) <- glue( 
      "  // data needed for ARMA correlations\n",
      "  int<lower=0> Kar{p};  // AR order\n",
      "  int<lower=0> Kma{p};  // MA order\n"
    )
    str_add(out$tdata_def) <- glue( 
      "  int max_lag{p} = max(Kar{p}, Kma{p});\n"
    )
    ar_bound <- ma_bound <- "<lower=-1,upper=1>"
    if (!use_cov(autocor)) {
      err_msg <- "Please set cov = TRUE in ARMA correlation structures"
      if (!has_natural_residuals) {
        stop2(err_msg, " for family '", family$family, "'.")
      }
      if (is.formula(bterms$adforms$se)) {
        stop2(err_msg, " when including known standard errors.")
      }
      str_add(out$data) <- glue( 
        "  // number of lags per observation\n",
        "  int<lower=0> J_lag{p}[N{p}];\n"                
      )
      str_add(out$model_def) <- glue(
        "  // objects storing residuals\n",
        "  matrix[N{p}, max_lag{p}] Err{p}",
        " = rep_matrix(0, N{p}, max_lag{p});\n",
        "  vector[N{p}] err{p};\n"
      )
      Y <- str_if(is.formula(bterms$adforms$mi), "Yl", "Y")
      add_ar <- str_if(get_ar(autocor),
                       glue("    mu{p}[n] += Err{p}[n, 1:Kar{p}] * ar{p};\n")             
      )
      add_ma <- str_if(get_ma(autocor),
                       glue("    mu{p}[n] += Err{p}[n, 1:Kma{p}] * ma{p};\n")             
      )
      str_add(out$model_comp_arma) <- glue(
        "  // include ARMA terms\n",
        "  for (n in 1:N{p}) {{\n",
        add_ma,
        "    err{p}[n] = {Y}{p}[n] - mu{p}[n];\n",
        "    for (i in 1:J_lag{p}[n]) {{\n",
        "      Err{p}[n + 1, i] = err{p}[n + 1 - i];\n",
        "    }}\n",
        add_ar,
        "  }}\n"
      )
      # in the conditional formulation no boundaries are required
      ar_bound <- ma_bound <- ""
    }
    if (get_ar(autocor)) {
      str_add(out$par) <- glue( 
        "  vector{ar_bound}[Kar{p}] ar{p};  // autoregressive effects\n"
      )
      str_add(out$prior) <- stan_prior(
        prior, class = "ar", px = px, suffix = p
      )
    }
    if (get_ma(autocor)) {
      str_add(out$par) <- glue( 
        "  vector{ma_bound}[Kma{p}] ma{p};  // moving-average effects\n"
      )
      str_add(out$prior) <- stan_prior(
        prior, class = "ma", px = px, suffix = p
      )
    }
  }
  if (is.cor_cosy(autocor)) {
    # compound symmetry correlation structure
    err_msg <- "Compound symmetry models are not implemented"
    if (is.mixfamily(family)) {
      stop2(err_msg, " for mixture models.") 
    }
    # most code is shared with ARMA covariance models
    str_add(out$par) <- glue(
      "  real<lower=0,upper=1> cosy{p};  // compound symmetry correlation\n"
    )
    str_add(out$prior) <- stan_prior(
      prior, class = "cosy", px = px, suffix = p
    )
  }
  if (use_cov(autocor)) {
    # use correlation structures in covariance matrix parameterization
    # optional for ARMA models and obligatory for compound symmetry models
    err_msg <- "Cannot model residual covariance matrices via 'autocor'"
    if (isTRUE(bterms$rescor)) {
      stop2(err_msg, " when estimating 'rescor'.")
    }
    pred_other_pars <- any(!names(bterms$dpars) %in% c("mu", "sigma"))
    if (has_natural_residuals && pred_other_pars) {
      stop2(err_msg, " when predicting parameters other than 'mu' and 'sigma'.")
    }
    str_add(out$data) <- glue( 
      "  // see the functions block for details\n",
      "  int<lower=1> N_tg{p};\n",
      "  int<lower=1> begin_tg{p}[N_tg{p}];\n", 
      "  int<lower=1> end_tg{p}[N_tg{p}];\n", 
      "  int<lower=1> nobs_tg{p}[N_tg{p}];\n"
    )
    if (!is.formula(bterms$adforms$se)) {
      str_add(out$tdata_def) <- glue(
        "  // no known standard errors specified by the user\n",
        "  vector[N{p}] se2{p} = rep_vector(0, N{p});\n"
      )
    }
    str_add(out$tpar_def) <- glue(
      "  // cholesky factor of the autocorrelation matrix\n",
      "  matrix[max(nobs_tg{p}), max(nobs_tg{p})] chol_cor{p};\n"               
    )
    if (is.cor_arma(autocor)) {
      if (has_ar_only(autocor)) {
        cor_fun <- "ar1"
        cor_args <- glue("ar{p}[1]")
      } else if (has_ma_only(autocor)) {
        cor_fun <- "ma1"
        cor_args <- glue("ma{p}[1]")
      } else {
        cor_fun <- "arma1"
        cor_args <- glue("ar{p}[1], ma{p}[1]")
      }
    }
    if (is.cor_cosy(autocor)) {
      cor_fun <- "cosy"
      cor_args <- glue("cosy{p}")
    }
    str_add(out$tpar_comp) <- glue(
      "  // compute residual covariance matrix\n",
      "  chol_cor{p} = cholesky_cor_{cor_fun}({cor_args}, max(nobs_tg{p}));\n"
    )
    if (has_latent_residuals) {
      err_msg <- "Latent residuals are not implemented"
      if (is.btnl(bterms$dpars[["mu"]])) {
        stop2(err_msg, " for non-linear models.")
      }
      if (conv_cats_dpars(bterms)) {
        stop2(err_msg, " for this family.")
      }
      str_add(out$par) <- glue(
        "  vector[N{p}] zerr{p};  // unscaled residuals\n",
        "  real<lower=0> sderr{p};  // SD of residuals\n"
      )
      str_add(out$tpar_def) <- glue(
        "  vector[N{p}] err{p};  // actual residuals\n"
      )
      str_add(out$tpar_comp) <- glue(
        "  // compute correlated residuals\n",
        "  err{p} = scale_cov_err(",
        "zerr{p}, sderr{p}, chol_cor{p}, nobs_tg{p}, begin_tg{p}, end_tg{p});\n"
      )
      str_add(out$prior) <- glue(
        "  target += normal_lpdf(zerr{p} | 0, 1);\n"
      )
      str_add(out$prior) <- stan_prior(
        prior, class = "sderr", px = px, suffix = p
      )
    }
  }
  if (is.cor_sar(autocor)) {
    err_msg <- "SAR models are not implemented"
    if (is.mixfamily(family)) {
      stop2(err_msg, " for mixture models.") 
    }
    if (!has_natural_residuals) {
      stop2(err_msg, " for family '", family$family, "'.")
    }
    if (isTRUE(bterms$rescor)) {
      stop2(err_msg, " when 'rescor' is estimated.")
    }
    if (any(c("sigma", "nu") %in% names(bterms$dpars))) {
      stop2(err_msg, " when predicting 'sigma' or 'nu'.")
    }
    str_add(out$data) <- glue(
      "  matrix[N{p}, N{p}] W{p};  // spatial weight matrix\n",
      "  vector[N{p}] eigenW{p};  // eigenvalues of W{p}\n"
    )
    str_add(out$tdata_def) <- glue(
      "  // the eigenvalues define the boundaries of the SAR correlation\n",
      "  real min_eigenW{p} = min(eigenW{p});\n",
      "  real max_eigenW{p} = max(eigenW{p});\n"
    )
    if (identical(autocor$type, "lag")) {
      str_add(out$par) <- glue( 
        "  // lag-SAR correlation parameter\n",
        "  real<lower=min_eigenW{p},upper=max_eigenW{p}> lagsar{p};\n"
      )
      str_add(out$prior) <- stan_prior(
        prior, class = "lagsar", px = px, suffix = p
      )
    } else if (identical(autocor$type, "error")) {
      str_add(out$par) <- glue( 
        "  // error-SAR correlation parameter\n",
        "  real<lower=min_eigenW{p},upper=max_eigenW{p}> errorsar{p};\n"
      )
      str_add(out$prior) <- stan_prior(
        prior, class = "errorsar", px = px, suffix = p
      )
    }
  }
  if (is.cor_car(autocor)) {
    err_msg <- "CAR models are not implemented"
    if (is.mixfamily(family)) {
      stop2(err_msg, " for mixture models.") 
    }
    if (is.btnl(bterms$dpars[["mu"]])) {
      stop2(err_msg, " for non-linear models.")
    }
    str_add(out$data) <- glue(
      "  // data for the CAR structure\n",
      "  int<lower=1> Nloc{p};\n",
      "  int<lower=1> Jloc{p}[N{p}];\n",
      "  int<lower=0> Nedges{p};\n",
      "  int<lower=1> edges1{p}[Nedges{p}];\n",
      "  int<lower=1> edges2{p}[Nedges{p}];\n"
    )
    str_add(out$par) <- glue(
      "  real<lower=0> sdcar{p};  // SD of the CAR structure\n"
    )
    str_add(out$prior) <- stan_prior(
      prior, class = "sdcar", px = px, suffix = p
    )
    if (autocor$type %in% c("escar", "esicar")) {
      str_add(out$data) <- glue(
        "  vector[Nloc{p}] Nneigh{p};\n",
        "  vector[Nloc{p}] eigenW{p};\n"
      )
    }
    if (autocor$type %in% "escar") {
      str_add(out$par) <- glue(
        "  real<lower=0, upper=1> car{p};\n",
        "  vector[Nloc{p}] rcar{p};\n"
      )
      car_args <- c(
        "car", "sdcar", "Nloc", "Nedges", 
        "Nneigh", "eigenW", "edges1", "edges2"
      )
      car_args <- paste0(car_args, p, collapse = ", ")
      str_add(out$prior) <- stan_prior(
        prior, class = "car", px = px, suffix = p
      )
      str_add(out$prior) <- glue(
        "  target += sparse_car_lpdf(\n", 
        "    rcar{p} | {car_args}\n",
        "  );\n"
      )
    } else if (autocor$type %in% "esicar") {
      str_add(out$par) <- glue(
        "  vector[Nloc{p} - 1] zcar{p};\n"
      )
      str_add(out$tpar_def) <- glue(
        "  vector[Nloc{p}] rcar{p};\n"                
      )
      str_add(out$tpar_comp) <- glue(
        "  // sum-to-zero constraint\n",
        "  rcar[1:(Nloc{p} - 1)] = zcar{p};\n",
        "  rcar[Nloc{p}] = - sum(zcar{p});\n"
      )
      car_args <- c(
        "sdcar", "Nloc", "Nedges", 
        "Nneigh", "eigenW", "edges1", "edges2"
      )
      car_args <- paste0(car_args, p, collapse = ", ")
      str_add(out$prior) <- glue(
        "  target += sparse_icar_lpdf(\n", 
        "    rcar{p} | {car_args}\n",
        "  );\n"
      )
    } else if (autocor$type %in% "icar") {
      # intrinsic car based on the case study of Mitzi Morris
      # http://mc-stan.org/users/documentation/case-studies/icar_stan.html
      str_add(out$par) <- glue(
        "  // parameters for the ICAR structure\n",
        "  vector[Nloc{p}] zcar{p};\n"
      )
      str_add(out$tpar_def) <- glue(
        "  // scaled parameters for the ICAR structure\n",
        "  vector[Nloc{p}] rcar{p} = zcar{p} * sdcar{p};\n"
      )
      str_add(out$prior) <- glue(
        "  // improper prior on the spatial CAR component\n",
        "  target += -0.5 * dot_self(zcar{p}[edges1{p}] - zcar{p}[edges2{p}]);\n",
        "  // soft sum-to-zero constraint\n",
        "  target += normal_lpdf(sum(zcar{p}) | 0, 0.001 * Nloc{p});\n"
      )
    } else if (autocor$type %in% "bym2") {
      # BYM2 car based on the case study of Mitzi Morris
      # http://mc-stan.org/users/documentation/case-studies/icar_stan.html
      str_add(out$data) <- glue(
        "  // scaling factor of the spatial CAR component\n",
        "  real<lower=0> car_scale{p};\n"
      )
      str_add(out$par) <- glue(
        "  // parameters for the BYM2 structure\n",
        "  vector[Nloc{p}] zcar{p};  // spatial part\n",
        "  vector[Nloc{p}] nszcar{p};  // non-spatial part\n",
        "  // proportion of variance in the spatial part\n",
        "  real<lower=0,upper=1> rhocar{p};\n"
      )
      str_add(out$tpar_def) <- glue(
        "  // join the spatial and the non-spatial CAR component\n",
        "  vector[Nloc{p}] rcar{p} = (sqrt(1 - rhocar{p}) * nszcar{p}", 
        " + sqrt(rhocar{p} * inv(car_scale{p})) * zcar{p}) * sdcar{p};\n"
      )
      str_add(out$prior) <- stan_prior(
        prior, class = "rhocar", px = px, suffix = p
      )
      str_add(out$prior) <- glue(
        "  // improper prior on the spatial BYM2 component\n",
        "  target += -0.5 * dot_self(zcar{p}[edges1{p}] - zcar{p}[edges2{p}]);\n",
        "  // soft sum-to-zero constraint\n",
        "  target += normal_lpdf(sum(zcar{p}) | 0, 0.001 * Nloc{p});\n",
        "  // proper prior on the non-spatial BYM2 component\n",
        "  target += normal_lpdf(nszcar | 0, 1);\n"
      )
    }
  }
  if (is.cor_fixed(autocor)) {
    err_msg <- "Fixed residual covariance matrices are not implemented"
    if (is.mixfamily(family)) {
      stop2(err_msg, " for mixture models.") 
    }
    if (!has_natural_residuals) {
      stop2(err_msg, " for family '", family$family, "'.")
    }
    if (isTRUE(bterms$rescor)) {
      stop2(err_msg, " when 'rescor' is estimated.")
    }
    str_add(out$data) <- glue( 
      "  matrix[N{p}, N{p}] V{p};  // known residual covariance matrix\n"
    )
    if (family$family %in% "gaussian") {
      str_add(out$tdata_def) <- glue(
        "  matrix[N{p}, N{p}] LV{p} = cholesky_decompose(V{p});\n"
      )
    }
  }
  out
}

# global Stan definitions for noise-free variables
# @param meef output of tidy_meef
stan_Xme <- function(meef, prior) {
  stopifnot(is.meef_frame(meef))
  if (!nrow(meef)) {
    return(list())
  }
  out <- list()
  coefs <- rename(paste0("me", meef$xname))
  str_add(out$data) <- "  // data for noise-free variables\n"
  str_add(out$par) <- "  // parameters for noise free variables\n"
  groups <- unique(meef$grname)
  for (i in seq_along(groups)) {
    g <- groups[i]
    K <- which(meef$grname %in% g)
    if (nzchar(g)) {
      Nme <- glue("Nme_{i}")
      str_add(out$data) <- glue(
        "  int<lower=0> Nme_{i};  // number of latent values\n",
        "  int<lower=1> Jme_{i}[N];  // group index per observation\n"
      )
    } else {
      Nme <- "N"
    }
    str_add(out$data) <- glue(
      "  int<lower=1> Mme_{i};  // number of groups\n"
    )
    str_add(out$data) <- cglue(
      "  vector[{Nme}] Xn_{K};  // noisy values\n",
      "  vector<lower=0>[{Nme}] noise_{K};  // measurement noise\n"
    )
    str_add(out$par) <- cglue(
      "  vector[Mme_{i}] meanme_{i};  // latent means\n",
      "  vector<lower=0>[Mme_{i}] sdme_{i};  // latent SDs\n"
    )
    str_add(out$prior) <- glue(
      stan_prior(prior, "meanme", coef = coefs[K], suffix = usc(i)),
      stan_prior(prior, "sdme", coef = coefs[K], suffix = usc(i))
    )
    str_add(out$prior) <- cglue(
      "  target += normal_lpdf(Xn_{K} | Xme_{K}, noise_{K});\n"
    )
    if (meef$cor[K[1]] && length(K) > 1L) {
      str_add(out$data) <- glue(
        "  int<lower=1> NCme_{i};  // number of latent correlations\n"
      )
      str_add(out$par) <- glue(
        "  matrix[Mme_{i}, {Nme}] zme_{i};  // standardized latent values\n",
        "  // cholesky factor of the latent correlation matrix\n",
        "  cholesky_factor_corr[Mme_{i}] Lme_{i};\n"
      )
      str_add(out$tpar_def) <- glue(
        "  // obtain the actual latent values\n",
        "  matrix[{Nme}, Mme_{i}] Xme{i}", 
        " = rep_matrix(meanme_{i}', {Nme}) ", 
        " + (diag_pre_multiply(sdme_{i}, Lme_{i}) * zme_{i})';\n"
      )
      str_add(out$tpar_def) <- cglue(
        "  // using separate vectors increases efficiency\n",
        "  vector[{Nme}] Xme_{K} = Xme{i}[, {K}];\n"
      )
      str_add(out$prior) <- glue(
        "  target += normal_lpdf(to_vector(zme_{i}) | 0, 1);\n",
        stan_prior(prior, "Lme", group = g, suffix = usc(i))
      )
      str_add(out$gen_def) <- cglue(
        "  // obtain latent correlation matrix\n",
        "  corr_matrix[Mme_{i}] Corme_{i}", 
        " = multiply_lower_tri_self_transpose(Lme_{i});\n",
        "  vector<lower=-1,upper=1>[NCme_{i}] corme_{i};\n"
      )
      str_add(out$gen_comp) <- stan_cor_gen_comp(
        cor = glue("corme_{i}"), ncol = glue("Mme_{i}")
      )
    } else {
      str_add(out$par) <- cglue(
        "  vector[{Nme}] zme_{K};  // standardized latent values\n"
      )
      str_add(out$tpar_def) <- cglue(
        "  // obtain the actual latent values\n",
        "  vector[{Nme}] Xme_{K} = ",
        "meanme_{i}[{K}] + sdme_{i}[{K}] * zme_{K};\n"
      )
      str_add(out$prior) <- cglue(
        "  target += normal_lpdf(zme_{K} | 0, 1);\n"
      )
    }
  }
  out
}

# initialize and compute a linear predictor term in Stan language
# @param out list of character strings containing Stan code
# @param bterms btl object
# @param ranef output of tidy_ranef
# @param ilink character vector of length 2 defining the link to be applied
# @param ... currently unused
# @return list of character strings containing Stan code
stan_eta_combine <- function(out, bterms, ranef, ilink = c("", ""), ...) {
  stopifnot(is.list(out), is.btl(bterms), length(ilink) == 2L)
  px <- check_prefix(bterms)
  resp <- usc(bterms$resp)
  eta <- combine_prefix(px, keep_mu = TRUE, nlp = TRUE)
  out$eta <- sub("^[ \t\r\n]+\\+", "", out$eta, perl = TRUE)
  str_add(out$model_def) <- glue(
    "  // initialize linear predictor term\n",
    "  vector[N{resp}] {eta} ={out$eta};\n"
  )
  out$eta <- NULL
  str_add(out$loopeta) <- stan_eta_re(ranef, px = px)
  if (nzchar(out$loopeta)) {
    # parts of eta are computed in a loop over observations
    out$loopeta <- sub("^[ \t\r\n]+\\+", "", out$loopeta, perl = TRUE)
    str_add(out$model_comp_eta_loop) <- glue(
      "  for (n in 1:N{resp}) {{\n",
      "    // add more terms to the linear predictor\n",
      "    {eta}[n] +={out$loopeta};\n",
      "  }}\n"
    )
  }
  out$loopeta <- NULL
  # possibly transform eta before it is passed to the likelihood
  if (sum(nzchar(ilink))) {
    # make sure mu comes last as it might depend on other parameters
    is_mu <- isTRUE("mu" %in% dpar_class(bterms[["dpar"]]))
    position <- str_if(is_mu, "model_comp_mu_link", "model_comp_dpar_link")
    str_add(out[[position]]) <- glue(
      "  for (n in 1:N{resp}) {{\n",
      "    // apply the inverse link function\n",
      "    {eta}[n] = {ilink[1]}{eta}[n]{ilink[2]};\n",
      "  }}\n"
    )
  }
  out
}

# define Stan code to compute the fixef part of eta
# @param fixef names of the population-level effects
# @param bterms object of class 'btl'
# @return a single character string
stan_eta_fe <- function(fixef, bterms) {
  if (length(fixef)) {
    p <- usc(combine_prefix(bterms))
    center_X <- stan_center_X(bterms)
    decomp <- get_decomp(bterms$fe)
    sparse <- is_sparse(bterms$fe)
    if (sparse) {
      stopifnot(!center_X && decomp == "none")
      csr_args <- sargs(
        paste0(c("rows", "cols"), "(X", p, ")"),
        paste0(c("wX", "vX", "uX", "b"), p)
      )
      eta_fe <- glue("csr_matrix_times_vector({csr_args})")
    } else {
      sfx_X <- sfx_b <- ""
      if (decomp == "QR") {
        sfx_X <- sfx_b <- "Q"
      } else if (center_X) {
        sfx_X <- "c"
      }
      eta_fe <- glue("X{sfx_X}{p} * b{sfx_b}{p}")
    }
  } else { 
    resp <- usc(bterms$resp)
    eta_fe <- glue("rep_vector(0, N{resp})")
  }
  glue(" + {eta_fe}")
}

# write the group-level part of the linear predictor
# @return a single character string
stan_eta_re <- function(ranef, px = list()) {
  eta_re <- ""
  ranef <- subset2(ranef, type = c("", "mmc"), ls = px)
  for (id in unique(ranef$id)) {
    r <- subset2(ranef, id = id)
    rpx <- check_prefix(r)
    idp <- paste0(r$id, usc(combine_prefix(rpx)))
    idresp <- paste0(r$id, usc(rpx$resp))
    if (r$gtype[1] == "mm") {
      ng <- seq_along(r$gcall[[1]]$groups)
      for (i in seq_rows(r)) {
        str_add(eta_re) <- cglue(
          " + W_{idresp[i]}_{ng}[n]", 
          " * r_{idp[i]}_{r$cn[i]}[J_{idresp[i]}_{ng}[n]]",
          " * Z_{idp[i]}_{r$cn[i]}_{ng}[n]"
        ) 
      }
    } else {
      str_add(eta_re) <- cglue(
        " + r_{idp}_{r$cn}[J_{idresp}[n]] * Z_{idp}_{r$cn}[n]"
      )
    }
  }
  eta_re
}

# Stan code for group-level parameters in special predictor terms
# @param r data.frame created by tidy_ranef
# @return a character vector: one element per row of 'r' 
stan_eta_rsp <- function(r) {
  stopifnot(nrow(r) > 0L, length(unique(r$gtype)) == 1L)
  rpx <- check_prefix(r)
  idp <- paste0(r$id, usc(combine_prefix(rpx)))
  idresp <- paste0(r$id, usc(rpx$resp))
  if (r$gtype[1] == "mm") {
    ng <- seq_along(r$gcall[[1]]$groups)
    out <- rep("", nrow(r))
    for (i in seq_along(out)) {
      out[i] <- glue(
        "W_{idresp[i]}_{ng}[n] * r_{idp[i]}_{r$cn[i]}[J_{idresp[i]}_{ng}[n]]",
        collapse = " + "
      ) 
    }
  } else {
    out <- glue("r_{idp}_{r$cn}[J_{idresp}[n]]")
  }
  out
}

# does eta need to be transformed manually using the link functions
# @param family the model family
# @param cens_or_trunc is the model censored or truncated?
stan_eta_transform <- function(family, cens_or_trunc = FALSE) {
  transeta <- "transeta" %in% family_info(family, "specials")
  no_transform <- family$link == "identity" && !transeta || 
    (has_cat(family) || has_thres(family)) && !is.customfamily(family)
  !no_transform && !stan_has_built_in_fun(family, cens_or_trunc)
}

# correctly apply inverse link to eta
# @param dpar name of the parameter for which to define the link
# @param bterms object of class 'brmsterms'
# @param resp name of the response variable
# @return a single character string
stan_eta_ilink <- function(dpar, bterms, resp = "") {
  stopifnot(is.brmsterms(bterms))
  out <- rep("", 2)
  family <- bterms$dpars[[dpar]]$family
  cens_or_trunc <- stan_llh_adj(bterms$adforms, c("cens", "trunc"))
  if (stan_eta_transform(family, cens_or_trunc = cens_or_trunc)) {
    dpar_id <- dpar_id(dpar)
    pred_dpars <- names(bterms$dpars)
    shape <- glue("shape{dpar_id}{resp}")
    shape <- str_if(shape %in% pred_dpars, paste0(shape, "[n]"), shape)
    nu <- glue("nu{dpar_id}{resp}")
    nu <- str_if(nu %in% pred_dpars, paste0(nu, "[n]"), nu)
    family_link <- str_if(
      family$family %in% c("gamma", "hurdle_gamma", "exponential"),
      paste0(family$family, "_", family$link), family$family
    )
    ilink <- stan_ilink(family$link)
    out <- switch(family_link,
      c(glue("{ilink}("), ")"),
      gamma_log = c(glue("{shape} * exp(-("), "))"),
      gamma_inverse = c(glue("{shape} * ("), ")"),
      gamma_identity = c(glue("{shape} / ("), ")"),
      hurdle_gamma_log = c(glue("{shape} * exp(-("), "))"),
      hurdle_gamma_inverse = c(glue("{shape} * ("), ")"),
      hurdle_gamma_identity = c(glue("{shape} / ("), ")"),
      exponential_log = c("exp(-(", "))"),
      exponential_inverse = c("(", ")"),
      exponential_identity = c("inv(", ")"),
      weibull = c(glue("{ilink}("), glue(") / tgamma(1 + 1 / {shape})")),
      frechet = c(glue("{ilink}("), glue(") / tgamma(1 - 1 / {nu})"))
    )
  }
  out
}

# indicate if the population-level design matrix should be centered
# implies a temporary shift in the intercept of the model
stan_center_X <- function(x) {
  is.btl(x) && !no_center(x$fe) && has_intercept(x$fe) && 
    !fix_intercepts(x) && !is_sparse(x$fe)
}

# default Stan definitions for distributional parameters
# @param dpar name of a distributional parameter
# @param suffix optional suffix of the parameter name
# @param family optional brmsfamily object
# @param fixed should the parameter be fixed to a certain value?
stan_dpar_defs <- function(dpar, suffix = "", family = NULL, fixed = FALSE) {
  dpar <- as_one_character(dpar)
  suffix <- as_one_character(suffix)
  fixed <- as_one_logical(fixed)
  if (is.mixfamily(family)) {
    if (dpar_class(dpar) == "theta") {
      return("")  # theta is handled in stan_mixture
    }
    family <- family$mix[[as.numeric(dpar_id(dpar))]]
  }
  if (is.customfamily(family)) {
    dpar_class <- dpar_class(dpar)
    lb <- family$lb[[dpar_class]]
    ub <- family$ub[[dpar_class]]
    lb <- if (!is.na(lb)) glue("lower={lb}")
    ub <- if (!is.na(ub)) glue("upper={ub}")
    bounds <- paste0(c(lb, ub), collapse = ",")
    if (nzchar(bounds)) bounds <- glue("<{bounds}>")
    def <- glue("  real{bounds} {dpar}{suffix};\n")
    return(def)
  }
  if (fixed) {
    min_Y <- glue("min(Y{suffix})")
  } else {
    min_Y <- glue("min_Y{suffix}")
  }
  default_defs <- list(
    sigma = c(
      "  real<lower=0> ", 
      ";  // residual SD\n"
    ),
    shape = c(
      "  real<lower=0> ", 
      ";  // shape parameter\n"
    ),
    nu = c(
      "  real<lower=1> ", 
      ";  // degrees of freedom or shape\n"
    ),
    phi = c(
      "  real<lower=0> ", 
      ";  // precision parameter\n"
    ),
    kappa = c(
      "  real<lower=0> ", 
      ";  // precision parameter\n"
    ),
    beta = c(
      "  real<lower=0> ", 
      ";  // scale parameter\n"
    ),
    zi = c(
      "  real<lower=0,upper=1> ", 
      ";  // zero-inflation probability\n"
    ), 
    hu = c(
      "  real<lower=0,upper=1> ", 
      ";  // hurdle probability\n"
    ),
    zoi = c(
      "  real<lower=0,upper=1> ", 
      ";  // zero-one-inflation probability\n"
    ), 
    coi = c(
      "  real<lower=0,upper=1> ", 
      ";  // conditional one-inflation probability\n"
    ),
    bs = c(
      "  real<lower=0> ", 
      ";  // boundary separation parameter\n"
    ),
    ndt = c(
      glue("  real<lower=0,upper={min_Y}> "), 
      ";  // non-decision time parameter\n"
    ),
    bias = c(
      "  real<lower=0,upper=1> ", 
      ";  // initial bias parameter\n"
    ),
    disc = c(
      "  real<lower=0> ", 
      ";  // discrimination parameters\n"
    ),
    quantile = c(
      "  real<lower=0,upper=1> ", 
      ";  // quantile parameter\n"
    ),
    xi = c(
      "  real ", 
      ";  // shape parameter\n"
    ),
    alpha = c(
      "  real ",
      ";  // skewness parameter\n"
    )
  )
  def <- default_defs[[dpar_class(dpar)]]
  if (!is.null(def)) {
    def <- paste0(def[1], dpar, suffix, def[2])
  } else {
    def <- ""
  }
  def
}

# default Stan definitions for temporary distributional parameters
stan_dpar_defs_tmp <- function(dpar, suffix = "", family = NULL) {
  dpar <- as_one_character(dpar)
  suffix <- as_one_character(suffix)
  if (is.mixfamily(family)) {
    family <- family$mix[[as.numeric(dpar_id(dpar))]]
  }
  if (is.customfamily(family)) {
    return("")  # no temporary parameters in custom families
  }
  default_defs <- list(
    xi = c(
      "  real tmp_", 
      ";  // unscaled shape parameter\n"
    )
  )
  def <- default_defs[[dpar_class(dpar)]]
  if (!is.null(def)) {
    def <- paste0(def[1], dpar, suffix, def[2])
  } else {
    def <- ""
  }
  def
}

# Stan code for transformations of distributional parameters
stan_dpar_transform <- function(bterms) {
  stopifnot(is.brmsterms(bterms))
  out <- list()
  families <- family_names(bterms)
  p <- usc(combine_prefix(bterms))
  resp <- usc(bterms$resp)
  if (any(conv_cats_dpars(families))) {
    str_add(out$model_def) <- glue( 
      "  // linear predictor matrix\n",
      "  vector[ncat{p}] mu{p}[N{resp}];\n"
    )
    mu_dpars <- make_stan_names(glue("mu{bterms$family$cats}"))
    mu_dpars <- glue("{mu_dpars}{p}[n]")
    iref <- match(bterms$family$refcat, bterms$family$cats)
    mu_dpars[iref] <- "0" 
    str_add(out$model_comp_catjoin) <- glue(
      "  for (n in 1:N{p}) {{\n",
      "    mu{p}[n] = {stan_vector(mu_dpars)};\n",
      "  }}\n"
    )
  }
  if (any(families %in% "skew_normal")) {
    # as suggested by Stephen Martin use sigma and mu of CP 
    # but the skewness parameter alpha of DP
    dp_names <- names(bterms$dpars)
    for (i in which(families %in% "skew_normal")) {
      id <- str_if(length(families) == 1L, "", i)
      sigma <- stan_sigma_transform(bterms, id = id)
      ns <- str_if(grepl("\\[n\\]", sigma), "[n]")
      na <- str_if(glue("alpha{id}") %in% dp_names, "[n]")
      type_delta <- str_if(nzchar(na), glue("vector[N{resp}]"), "real")
      no <- str_if(any(nzchar(c(ns, na))), "[n]", "")
      type_omega <- str_if(nzchar(no), glue("vector[N{resp}]"), "real")
      str_add(out$model_def) <- glue(
        "  // parameters used to transform the skew-normal distribution\n",
        "  {type_delta} delta{id}{p};  // transformed alpha parameter\n",
        "  {type_omega} omega{id}{p};  // scale parameter\n"
      )
      alpha <- glue("alpha{id}{p}{na}")
      delta <- glue("delta{id}{p}{na}")
      omega <- glue("omega{id}{p}{no}")
      comp_delta <- glue(
        "  {delta} = {alpha} / sqrt(1 + {alpha}^2);\n"
      )
      comp_omega <- glue(
        "  {omega} = {sigma} / sqrt(1 - sqrt_2_div_pi^2 * {delta}^2);\n"
      )
      str_add(out$model_comp_dpar_trans) <- glue(
        "   // use efficient skew-normal parameterization\n",
        str_if(!nzchar(na), comp_delta),
        str_if(!nzchar(no), comp_omega),
        "  for (n in 1:N{resp}) {{\n",
        str_if(nzchar(na), glue("  ", comp_delta)),
        str_if(nzchar(no), glue("  ", comp_omega)),
        "    mu{id}{p}[n] = mu{id}{p}[n]", 
        " - {omega} * {delta} * sqrt_2_div_pi;\n",
        "  }}\n"
      )
    }
  }
  if (any(families %in% "gen_extreme_value")) {
    dp_names <- c(names(bterms$dpars), names(bterms$fdpars))
    for (i in which(families %in% "gen_extreme_value")) {
      id <- str_if(length(families) == 1L, "", i)
      xi <- glue("xi{id}")
      if (!xi %in% dp_names) {
        str_add(out$model_def) <- glue(
          "  real {xi};  // scaled shape parameter\n"
        )
        sigma <- glue("sigma{id}")
        v <- str_if(sigma %in% names(bterms$dpars), "_vector")
        args <- sargs(
          glue("tmp_{xi}"), glue("Y{p}"), 
          glue("mu{id}{p}"), glue("{sigma}{p}")
        )
        str_add(out$model_comp_dpar_trans) <- glue(
          "  {xi}{p} = scale_xi{v}({args});\n"
        )
      }
    }
  }
  out
}

# Stan code for sigma to incorporate addition argument 'se'
stan_sigma_transform <- function(bterms, id = "") {
  if (nzchar(id)) {
    family <- family_names(bterms)[as.integer(id)]
  } else {
    family <- bterms$family$family
  }
  p <- usc(combine_prefix(bterms))
  ns <- str_if(glue("sigma{id}") %in% names(bterms$dpars), "[n]")
  has_sigma <- has_sigma(family) && !no_sigma(bterms)
  sigma <- str_if(has_sigma, glue("sigma{id}{p}{ns}"))
  if (is.formula(bterms$adforms$se)) {
    sigma <- str_if(nzchar(sigma), 
      glue("sqrt({sigma}^2 + se2{p}[n])"), 
      glue("se{p}[n]")
    )
  }
  sigma
}
