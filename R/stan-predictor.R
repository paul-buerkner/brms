stan_effects.btl <- function(x, data, ranef, prior, center_X = TRUE, 
                             sparse = FALSE, threshold = "flexible",
                             nlpar = "", eta = "mu", ilink = rep("", 2),
                             ...) {
  # combine effects for the predictors of a single (non-linear) parameter
  # Args:
  #   center_X: center population-level design matrix if possible?
  #   eta: prefix of the linear predictor variable
  nlpar <- check_nlpar(nlpar)
  if (nzchar(eta) && nzchar(nlpar)) {
    eta <- usc(eta, "suffix") 
  }
  eta <- paste0(eta, nlpar)
  stopifnot(nzchar(eta))
  stopifnot(length(ilink) == 2L)
  
  ranef <- ranef[ranef$nlpar == nlpar, ]
  out <- list()
  out$modelD <- paste0("  vector[N] ", eta, "; \n")
  # include population-level effects
  center_X <- center_X && has_intercept(x$fe) && 
    !is(x$autocor, "cor_bsts") && !sparse
  rm_int <- center_X || is(x$autocor, "cor_bsts") || is_ordinal(x$family)
  cols2remove <- if (rm_int) "Intercept"
  fixef <- setdiff(colnames(data_fe(x, data)$X), cols2remove)
  text_fe <- stan_fe(
    fixef, center_X = center_X, family = x$family, 
    prior = prior, nlpar = nlpar, sparse = sparse, 
    threshold = threshold
  )
  # include smooth terms
  smooths <- get_sm_labels(x, data = data)
  text_sm <- stan_sm(smooths, prior = prior, nlpar = nlpar)
  # include category specific effects
  csef <- colnames(get_model_matrix(x$cs, data))
  text_cs <- stan_cs(csef, ranef, prior = prior)
  # include monotonic effects
  monef <- all_terms(x$mo)
  text_mo <- stan_mo(monef, ranef, prior = prior, nlpar = nlpar)
  # include measurement error variables
  meef <- get_me_labels(x, data = data)
  text_me <- stan_me(meef, ranef, prior = prior, nlpar = nlpar)
  out <- collapse_lists(list(
    out, text_fe, text_cs, text_mo, text_me, text_sm
  ))
  
  p <- usc(nlpar, "prefix")
  has_offset <- !is.null(get_offset(x$fe))
  if (has_offset) {
    out$data <- paste0(out$data, "  vector[N] offset", p, "; \n")
  }
  
  # initialize and compute eta_<nlpar>
  out$modelC1 <- paste0(
    out$modelC1, "  ", eta, " = ", 
    text_fe$eta, text_sm$eta,
    if (center_X && !is_ordinal(x$family)) 
      paste0(" + temp", p, "_Intercept"),
    if (has_offset) paste0(" + offset", p),
    if (get_arr(x$autocor)) " + Yarr * arr", 
    "; \n"
  )
  
  # repare loop over eta
  eta_ma <- ifelse(get_ma(x$autocor) && !use_cov(x$autocor),
                   paste0(" + head(E", p, "[n], Kma) * ma"), "")
  eta_loop <- paste0(
    stan_eta_re(ranef, nlpar = nlpar),
    text_mo$eta, text_me$eta,
    eta_ma, stan_eta_bsts(x$autocor)
  )
  if (nzchar(eta_loop)) {
    out$modelC2 <- paste0(out$modelC2,
      "    ", eta, "[n] = ", eta, "[n]", eta_loop, "; \n"
    )
  }
  # include autoregressive effects
  if (get_ar(x$autocor) && !use_cov(x$autocor)) {
    eta_ar <- paste0(eta, "[n] + head(E", p, "[n], Kar) * ar")
    out$modelC3 <- paste0(out$modelC3, 
      "    ", eta, "[n] = ", eta_ar, "; \n"
    )
  }
  # possibly transform eta before it is passed to the likelihood
  if (sum(nzchar(ilink))) {
    # make sure mu comes last as it might depend on other parameters
    position <- ifelse(nzchar(nlpar), "modelC3", "modelC4")
    out[[position]] <- paste0(out[[position]],
      "    ", eta, "[n] = ", ilink[1], eta, "[n]", ilink[2], "; \n"
    )
  }
  out
}

stan_effects.btnl <- function(x, data, ranef, prior, eta = "mu", 
                              ilink = rep("", 2), ...) {
  # prepare Stan code for non-linear models
  # Args:
  #   data: data.frame supplied by the user
  #   ranef: data.frame returned by tidy_ranef
  #   prior: a brmsprior object
  stopifnot(length(ilink) == 2L)
  out <- list()
  if (length(x$nlpars)) {
    nlpars <- names(x$nlpars)
    for (nlp in nlpars) {
      nl_text <- stan_effects(
        x = x$nlpars[[nlp]], data = data, 
        ranef = ranef, prior = prior, 
        nlpar = nlp, center_X = FALSE
      )
      out <- collapse_lists(list(out, nl_text))
    }
    # prepare non-linear model
    new_nlpars <- paste0(" ", eta, "_", nlpars, "[n] ")
    # covariates in the non-linear model
    covars <- wsp(setdiff(all.vars(rhs(x$formula)), nlpars))
    if (length(covars)) {
      out$data <- paste0(out$data, 
        "  int<lower=1> KC;  // number of covariates \n",
        " matrix[N, KC] C;  // covariate matrix \n"
      )
      new_covars <- paste0(" C[n, ", seq_along(covars), "] ")
    } else {
      new_covars <- NULL
    }
    # add whitespaces to be able to replace parameters and covariates
    meta_sym <- c("+", "-", "*", "/", "^", ")", "(", ",")
    nlmodel <- gsub(" ", "", collapse(deparse(x$formula[[2]])))
    nlmodel <- wsp(rename(nlmodel, meta_sym, wsp(meta_sym))) 
    nlmodel <- rename(nlmodel, 
      c(wsp(nlpars), covars, " ( ", " ) "), 
      c(new_nlpars, new_covars, "(", ")")
    )
    # possibly transform eta in the transformed params block
    out$modelD <- paste0(out$modelD, "  vector[N] ", eta, "; \n")
    out$modelC4 <- paste0(out$modelC4, 
      "    // compute non-linear predictor \n",
      "    ", eta, "[n] = ", ilink[1], trimws(nlmodel), ilink[2], "; \n"
    )
  }
  out
}

stan_effects.brmsterms <- function(x, data, ranef, prior, sparse = FALSE, 
                                   threshold = "flexible", ...) {
  # Stan code for auxiliary parameters
  # Args:
  #   bterms: object of class brmsterms
  #   other arguments: same as make_stancode
  out <- list()
  default_defs <- c(
    mu = "",  # mu is always predicted
    sigma = "  real<lower=0> sigma;  // residual SD \n",
    shape = "  real<lower=0> shape;  // shape parameter \n",
    nu = "  real<lower=1> nu;  // degrees of freedom or shape \n",
    phi = "  real<lower=0> phi;  // precision parameter \n",
    kappa = "  real<lower=0> kappa;  // precision parameter \n",
    beta = "  real<lower=0> beta;  // scale parameter \n",
    zi = "  real<lower=0,upper=1> zi;  // zero-inflation probability \n", 
    hu = "  real<lower=0,upper=1> hu;  // hurdle probability \n",
    bs = "  real<lower=0> bs;  // boundary separation parameter \n",
    ndt = "  real<lower=0,upper=min_Y> ndt;  // non-decision time parameter \n",
    bias = "  real<lower=0,upper=1> bias;  // initial bias parameter \n",
    disc = "  real<lower=0> disc;  // discrimination parameters \n",
    quantile = "  real<lower=0,upper=1> quantile;  // quantile parameter \n",
    xi = "  real xi;  // shape parameter \n"
  )
  default_defs_temp <- c(
    xi = "  real temp_xi;  // unscaled shape parameter \n"
  )
  valid_auxpars <- valid_auxpars(bterms = x)
  args <- nlist(data, ranef, prior)
  for (ap in valid_auxpars) {
    ap_terms <- x$auxpars[[ap]]
    if (is.btl(ap_terms) || is.btnl(ap_terms)) {
      ilink <- stan_eta_ilink(
        ap_terms$family, auxpars = names(x$auxpars), adforms = x$adforms
      )
      eta <- ifelse(ap == "mu", "mu", "")
      ap_args <- list(ap_terms, nlpar = ap, eta = eta, ilink = ilink)
      out[[ap]] <- do.call(stan_effects, c(ap_args, args))
    } else if (is.numeric(x$fauxpars[[ap]])) {
      out[[ap]] <- list(data = default_defs[ap]) 
    } else {
      if (ap %in% names(default_defs_temp)) {
        out[[ap]] <- list(
          par = default_defs_temp[ap],
          prior = stan_prior(prior, class = ap, prefix = "temp_")
        )
      } else {
        out[[ap]] <- list(
          par = default_defs[ap],
          prior = stan_prior(prior, class = ap)
        )
      }
    }
  }
  collapse_lists(out)
}

stan_effects_mv <- function(bterms, data, ranef, prior, sparse = FALSE) {
  # Stan code for multivariate models
  # Args:
  #   see stan_effects
  stopifnot(is.brmsterms(bterms))
  out <- list()
  resp <- bterms$response
  if (length(resp) > 1L) {
    args <- nlist(
      x = bterms$auxpars[["mu"]], data, ranef, prior, sparse,
      ilink = stan_eta_ilink(bterms$family, adforms = bterms$adforms)
    )
    resp <- bterms$response
    tmp_list <- named_list(resp)
    for (r in resp) {
      tmp_list[[r]] <- do.call(stan_effects, c(args, nlpar = r))
    }
    out <- collapse_lists(tmp_list)
    if (is_linear(bterms$family)) {
      len_Eta_n <- "nresp" 
    } else if (is_categorical(bterms$family)) {
      len_Eta_n <- "ncat - 1"
    } else {
      stop2("Multivariate models are not yet implemented ", 
            "for family '", bterms$family$family, "'.")
    }
    out$modelD <- paste0(out$modelD, 
      "  // multivariate linear predictor matrix \n",
      "  vector[", len_Eta_n, "] Mu[N]; \n"
    )
    out$modelC3 <- paste0(out$modelC3, 
      collapse("    Mu[n, ", seq_along(resp), "] = mu_", resp, "[n]; \n")
    )
  }
  out
}

stan_fe <- function(fixef, prior, family = gaussian(),
                    center_X = TRUE, nlpar = "", sparse = FALSE,
                    threshold = "flexible") {
  # Stan code for population-level effects
  # Args:
  #   fixef: names of the population-level effects
  #   center_X: center the design matrix?
  #   family: the model family
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior 
  #   threshold: either "flexible" or "equidistant" 
  # Returns:
  #   a list containing Stan code related to population-level effects
  p <- usc(nlpar, "prefix")
  ct <- ifelse(center_X, "c", "")
  out <- list()
  if (length(fixef)) {
    out$data <- paste0(out$data, 
      "  int<lower=1> K", p, ";",
      "  // number of population-level effects \n", 
      "  matrix[N, K", p, "] X", p, ";",
      "  // population-level design matrix \n"
    )
    if (sparse) {
      stopifnot(!center_X)
      if (nchar(nlpar)) {
        stop2("Sparse matrices are not yet implemented for this model.")
      }
      out$tdataD <- "  #include tdataD_sparse_X.stan \n"
      out$tdataC <- "  #include tdataC_sparse_X.stan \n"
    }
    # prepare population-level coefficients
    if (!is.null(attr(prior, "hs_df")) && !nzchar(nlpar)) {
      out$par <- paste0(out$par,
        "  // horseshoe shrinkage parameters \n",
        "  vector[K", ct, "] zb; \n",
        "  vector<lower=0>[K", ct, "] hs_local[2]; \n",
        "  real<lower=0> hs_global[2]; \n"
      )
      out$transD <- paste0(out$transD, 
        "  vector[K", ct, "] b;",
        "  // population-level effects \n"
      )
      hs_scale_global <- attr(prior, "hs_scale_global")
      hs_args <- sargs("zb", "hs_local", "hs_global", hs_scale_global)
      out$transD <- paste0(out$transD, 
        "  b = horseshoe(", hs_args, "); \n"
      )
    } else {
      bound <- get_bound(prior, class = "b", nlpar = nlpar)
      out$par <- paste0(out$par,
        "  vector", bound, "[K", ct, p, "] b", p, ";",
        "  // population-level effects \n"
      )
    }
    if (!is.null(attr(prior, "lasso_df")) && !nzchar(nlpar)) {
      out$par <- paste0(out$par,
        "  // lasso shrinkage parameter \n",
        "  real<lower=0> lasso_inv_lambda; \n"
      )
    }
    fixef_prior <- stan_prior(prior, class = "b", coef = fixef,
                              nlpar = nlpar, suffix = p)
    out$prior <- paste0(out$prior, fixef_prior)
  }
  if (center_X) {
    if (length(fixef)) {
      out$tdataD <- paste0(out$tdataD, 
        "  int Kc", p, "; \n",
        "  matrix[N, K", p, " - 1] Xc", p, ";", 
        "  // centered version of X", p, " \n",
        "  vector[K", p, " - 1] means_X", p, ";",
        "  // column means of X", p, " before centering \n"
      )
      out$tdataC <- paste0(out$tdataC, 
        "  Kc", p, " = K", p, " - 1;",
        "  // the intercept is removed from the design matrix \n",
        "  for (i in 2:K", p, ") { \n",
        "    means_X", p, "[i - 1] = mean(X", p, "[, i]); \n",
        "    Xc", p, "[, i - 1] = X", p, "[, i] - means_X", p, "[i - 1]; \n",
        "  } \n"
      )
      # cumulative and sratio models are parameterized as thres - eta
      use_plus <- family$family %in% c("cumulative", "sratio")
      sub_X_means <- paste0(
        ifelse(use_plus, " + ", " - "), 
        "dot_product(means_X", p, ", b", p, ")"
      )
    } else {
      sub_X_means <- ""
    }
    if (is_ordinal(family)) {
      # temp intercepts for ordinal models are defined in stan_ordinal
      out$genD <- "  vector[ncat - 1] b_Intercept;  // thresholds \n" 
      out$genC <- paste0(
        "  b_Intercept = temp_Intercept", sub_X_means, "; \n"
      )
    } else {
      out$par <- paste0(out$par, 
        "  real temp", p, "_Intercept;  // temporary intercept \n"
      )
      out$genD <- paste0(
        "  real b", p, "_Intercept;  // population-level intercept \n"
      )
      out$genC <- paste0(
        "  b", p, "_Intercept = temp", p, "_Intercept", sub_X_means, "; \n"
      )
    }
    # for equidistant thresholds only temp_Intercept1 is a parameter
    prefix <- paste0("temp", p, "_")
    suffix <- ifelse(threshold == "equidistant", "1", "")
    int_prior <- stan_prior(prior, class = "Intercept", nlpar = nlpar,
                            prefix = prefix, suffix = suffix)
    out$prior <- paste0(out$prior, int_prior)
  }
  out$eta <- stan_eta_fe(fixef, center_X, sparse, nlpar)
  out
}

stan_re <- function(id, ranef, prior, cov_ranef = NULL) {
  # group-level effects in Stan 
  # Args:
  #   id: the ID of the grouping factor
  #   ranef: a data.frame returned by tidy_ranef
  #   prior: object of class brmsprior
  #   cov_ranef: a list of custom covariance matrices 
  # Returns:
  #   A list of strings containing Stan code
  r <- ranef[ranef$id == id, ]
  ccov <- r$group[1] %in% names(cov_ranef)
  ng <- seq_along(r$gcall[[1]]$groups)
  idp <- paste0(r$id, usc(r$nlpar, "prefix"))
  out <- list()
  out$data <- paste0(
    "  // data for group-level effects of ID ", id, " \n",
    if (r$gtype[1] == "mm") {
      collapse(
        "  int<lower=1> J_", id, "_", ng, "[N]; \n",
        "  real<lower=0> W_", id, "_", ng, "[N]; \n"
      )
    } else {
      paste0("  int<lower=1> J_", id, "[N]; \n")
    },
    "  int<lower=1> N_", id, "; \n",
    "  int<lower=1> M_", id, "; \n",
    if (ccov) paste0(
      "  // cholesky factor of known covariance matrix \n",
      "  matrix[N_", id, ", N_", id,"] Lcov_", id,"; \n"
    )
  )
  out$prior <- stan_prior(prior, class = "sd", group = r$group[1], 
                          coef = r$coef, nlpar = r$nlpar, 
                          suffix = paste0("_", id))
  J <- seq_len(nrow(r))
  has_def_type <- !r$type %in% c("mo", "me")
  if (any(has_def_type)) {
    out$data <- paste0(out$data, 
      collapse("  vector[N] Z_", idp[has_def_type], 
               "_", r$cn[has_def_type], "; \n")) 
  }
  out$par <- paste0(
    "  vector<lower=0>[M_", id, "] sd_", id, ";",
    "  // group-level standard deviations \n"
  )
  if (nrow(r) > 1L && r$cor[1]) {
    # multiple correlated group-level effects
    out$data <- paste0(out$data, "  int<lower=1> NC_", id, "; \n")
    out$par <- paste0(out$par,
      "  matrix[M_", id, ", N_", id, "] z_", id, ";",
      "  // unscaled group-level effects \n",    
      "  // cholesky factor of correlation matrix \n",
      "  cholesky_factor_corr[M_", id, "] L_", id, "; \n"
    )
    out$prior <- paste0(out$prior, 
      stan_prior(prior, class = "L", group = r$group[1],
                 suffix = paste0("_", id)),
      "  to_vector(z_", id, ") ~ normal(0, 1); \n"
    )
    out$transD <- paste0(
      "  // group-level effects \n",
      "  matrix[N_", id, ", M_", id, "] r_", id, "; \n",
      collapse("  vector[N_", id, "] r_", idp, "_", r$cn, "; \n")
    )
    if (ccov) {  
      # customized covariance matrix supplied
      out$transC1 <- paste0(
        "  r_", id," = as_matrix(kronecker(Lcov_", id, ",", 
        " diag_pre_multiply(sd_", id,", L_", id,")) *",
        " to_vector(z_", id, "), N_", id, ", M_", id, "); \n"
      )
    } else { 
      out$transC1 <- paste0("  r_", id, " = ", 
        "(diag_pre_multiply(sd_", id, ", L_", id,") * z_", id, ")'; \n"
      )
    }
    out$transC1 <- paste0(out$transC1, 
      collapse("  r_", idp, "_", r$cn, " = r_", id, "[, ", J, "];  \n")
    )
    # return correlations above the diagonal only
    cors_genC <- ulapply(2:nrow(r), function(k) 
      lapply(1:(k - 1), function(j) paste0(
        "  cor_", id, "[", (k - 1) * (k - 2) / 2 + j, 
        "] = Cor_", id, "[", j, ",", k, "]; \n")
      )
    )
    out$genD <- paste0(
      "  corr_matrix[M_", id, "] Cor_", id, "; \n",
      "  vector<lower=-1,upper=1>[NC_", id, "] cor_", id, "; \n"
    )
    out$genC <- paste0(
      "  // take only relevant parts of correlation matrix \n",
      "  Cor_", id, " = multiply_lower_tri_self_transpose(L_", id, "); \n",
      collapse(cors_genC)
    ) 
  } else {
    # single or uncorrelated group-level effects
    out$par <- paste0(out$par,
      "  vector[N_", id, "] z_", id, "[M_", id, "];",
      "  // unscaled group-level effects \n"
    )
    out$prior <- paste0(out$prior, collapse(
      "  z_", id, "[", 1:nrow(r), "] ~ normal(0, 1); \n")
    )
    out$transD <- paste0("  // group-level effects \n", 
      collapse("  vector[N_", id, "] r_", idp, "_", r$cn, "; \n")
    )
    out$transC1 <- collapse(
      "  r_", idp, "_", r$cn, " = sd_", id, "[", J, "] * (", 
      if (ccov) paste0("Lcov_", id, " * "), "z_", id, "[", J, "]); \n"
    )
  }
  out
}

stan_sm <- function(smooths, prior, nlpar = "") {
  # Stan code of smooth terms
  # Args:
  #   smooths: names of the smooth terms
  #   prior: object of class brmsprior
  #   nlpar: optional name of a non-linear parameter
  # Returns:
  #   A list of strings containing Stan code
  out <- list()
  p <- usc(nlpar)
  if (length(smooths)) {
    stopifnot(!is.null(attr(smooths, "nbases")))
    for (i in seq_along(smooths)) {
      pi <- paste0(p, "_", i)
      nb <- seq_len(attr(smooths, "nbases")[[i]])
      out$data <- paste0(out$data,
        "  // data of smooth ", smooths[i], "\n",  
        "  int nb", pi, ";  // number of bases \n",
        "  int knots", pi, "[nb", pi, "]; \n"
      )
      out$data <- paste0(out$data, collapse(
        "  matrix[N, knots", pi, "[", nb, "]]", 
        " Zs", pi, "_", nb, "; \n")
      )
      out$par <- paste0(out$par,
        "  // parameters of smooth ", smooths[i], "\n"
      )
      out$par <- paste0(out$par, collapse(
        "  vector[knots", pi, "[", nb, "]] zs", pi,"_", nb, "; \n",
        "  real<lower=0> sds", pi, "_", nb, "; \n")
      )
      out$transD <- paste0(out$transD, collapse(
        "  vector[knots", pi, "[", nb, "]] s", pi, "_", nb, "; \n")
      )
      out$transC1 <- paste0(out$transC1, collapse(
        "  s", pi, "_", nb, " = sds", pi,  "_", nb, 
        " * zs", pi, "_", nb, "; \n")
      )
      out$prior <- paste0(out$prior, collapse(
        "  zs", pi, "_", nb, " ~ normal(0, 1); \n"),
        stan_prior(prior, class = "sds", coef = smooths[i], 
                   nlpar = nlpar, suffix = paste0(pi, "_", nb))
      )
    }
    out$eta <- stan_eta_sm(smooths, nlpar = nlpar)
  }
  out
}

stan_mo <- function(monef, ranef, prior, nlpar = "") {
  # Stan code for monotonic effects
  # Args:
  #   monef: names of the monotonic effects
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior
  p <- if (nchar(nlpar)) paste0("_", nlpar) else ""
  out <- list()
  if (length(monef)) {
    I <- seq_along(monef)
    out$data <- paste0(
      "  int<lower=1> Kmo", p, ";  // number of monotonic effects \n",
      "  int Xmo", p, "[N, Kmo", p, "];  // monotonic design matrix \n",
      "  int<lower=2> Jmo", p, "[Kmo", p, "];  // length of simplexes \n",
      collapse(
        "  vector[Jmo", p, "[", I, "]]", 
        " con_simplex", p, "_", I, "; \n"
      )
    )
    bound <- get_bound(prior, class = "b", nlpar = nlpar)
    out$par <- paste0(
      "  // monotonic effects \n", 
      "  vector", bound, "[Kmo", p, "] bmo", p, "; \n",
      collapse(
        "  simplex[Jmo", p, "[", I, "]]", 
        " simplex", p, "_", I, "; \n"
      )
    ) 
    out$prior <- paste0(
      stan_prior(prior, class = "b", coef = monef,
                 nlpar = nlpar, suffix = paste0("mo", p)),
      collapse("  simplex", p, "_", I, 
               " ~ dirichlet(con_simplex", p, "_", I, "); \n")
    )
    out$eta <- stan_eta_mo(monef, ranef = ranef, nlpar = nlpar)
  }
  out
}

stan_cs <- function(csef, ranef, prior, nlpar = "") {
  # Stan code for category specific effects
  # Args:
  #   csef: names of the category specific effects
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior
  # (!) Not yet implemented for non-linear models
  stopifnot(!nzchar(nlpar))
  ranef <- ranef[ranef$nlpar == nlpar & ranef$type == "cs", ]
  out <- list()
  if (length(csef) || nrow(ranef)) {
    out$modelD <- paste0(
      "  // linear predictor for category specific effects \n",                  
      "  matrix[N, ncat - 1] mucs; \n"
    )
  }
  if (length(csef)) {
    out$data <- paste0(
      "  int<lower=1> Kcs;  // number of category specific effects \n",
      "  matrix[N, Kcs] Xcs;  // category specific design matrix \n"
    )
    bound <- get_bound(prior, class = "b")
    out$par <- paste0(
      "  matrix", bound, "[Kcs, ncat - 1] bcs;",
      "  // category specific effects \n"
    )
    out$modelC1 <- "  mucs = Xcs * bcs; \n"
    out$prior <- stan_prior(prior, class = "b", coef = csef,
                            suffix = "cs", matrix = TRUE)
  } 
  if (nrow(ranef)) {
    if (!length(csef)) {
      # only group-level category specific effects present
      out$modelC1 <- "  mucs = rep_matrix(0, N, ncat - 1); \n"
    }
    cats <- get_matches("\\[[[:digit:]]+\\]$", ranef$coef)
    ncatM1 <- max(as.numeric(substr(cats, 2, nchar(cats) - 1)))
    for (i in seq_len(ncatM1)) {
      r_cat <- ranef[grepl(paste0("\\[", i, "\\]$"), ranef$coef), ]
      out$modelC2 <- paste0(out$modelC2,
        "    mucs[n, ", i, "] = mucs[n, ", i, "]"
      )
      for (id in unique(r_cat$id)) {
        r <- r_cat[r_cat$id == id, ]
        idp <- paste0(r$id, usc(r$nlpar, "prefix"))
        out$modelC2 <- paste0(out$modelC2, collapse(
          " + r_", idp, "_", r$cn, "[J_", r$id, "[n]]",
          " * Z_", idp, "_", r$cn, "[n]")
        )
      }
      out$modelC2 <- paste0(out$modelC2, "; \n")
    }
  }
  out
} 

stan_me <- function(meef, ranef, prior, nlpar = "") {
  # stan code for measurement error effects
  # Args:
  #   meef: vector of terms containing noisy predictors
  out <- list()
  if (length(meef)) {
    not_one <- attr(meef, "not_one")
    uni_me <- attr(meef, "uni_me")
    p <- usc(nlpar)
    pK <- paste0(p, "_", seq_along(uni_me))
    
    me_sp <- strsplit(gsub("[[:space:]]", "", meef), ":")
    meef_terms <- rep(NA, length(me_sp))
    for (i in seq_along(me_sp)) {
      # remove non-me parts from the terms
      take <- grepl_expr("^me\\([^:]*\\)$", me_sp[[i]])
      me_sp[[i]] <- me_sp[[i]][take]
      # remove 'I' (identity) function calls that 
      # were used solely to separate formula terms
      I <- grepl("^I\\(", me_sp[[i]])
      me_sp[[i]][I] <- substr(me_sp[[i]][I], 3,  nchar(me_sp[[i]][I]) - 1)
      meef_terms[i] <- paste0(me_sp[[i]], collapse = ":")
    }
    new_me <- paste0("Xme", pK, "[n]")
    meef_terms <- rename(meef_terms, uni_me, new_me)
    ci <- ulapply(seq_along(not_one), function(i) sum(not_one[1:i]))
    covars <- ifelse(not_one, paste0(" .* Cme", p, "_", ci, "[n]"), "")
    ncovars <- sum(not_one)
    
    # prepare linear predictor component
    meef <- rename(meef)
    meef_terms <- gsub(":", " .* ", meef_terms)
    ranef <- ranef[ranef$nlpar == nlpar & ranef$type == "me", ]
    invalid_coef <- setdiff(ranef$coef, meef)
    if (length(invalid_coef)) {
      stop2("Noisy group-level terms require ", 
            "corresponding population-level terms.")
    }
    for (i in seq_along(meef)) {
      r <- ranef[ranef$coef == meef[i], ]
      if (nrow(r)) {
        rpars <- paste0(" + ", stan_eta_r(r))
      } else {
        rpars <- ""
      }
      out$eta <- paste0(out$eta,
        " + (bme", p, "[", i, "]", rpars, ") * ", meef_terms[i], covars[i]
      )
    }
    
    # prepare Stan code
    out$data <- paste0(
      "  int<lower=0> Kme", p, ";",
      "  // number of terms of noise free variables \n",
      "  // noisy variables \n",
      collapse("  vector[N] Xn", pK, "; \n"),
      "  // measurement noise \n",
      collapse("  vector<lower=0>[N] noise", pK, "; \n"),
      if (ncovars > 0L) paste0(
        "  // covariates of noise free variables \n",
        collapse("  vector[N] Cme", p, "_", seq_len(ncovars), "; \n")
      )
    )
    out$par <- paste0(
      "  // noise free variables \n",
      collapse("  vector[N] Xme", pK, "; \n"),  
      "  vector[Kme", p, "] bme", p, ";",
      "  // coefficients of noise-free terms \n"
    )
    out$prior <- paste0(
      stan_prior(prior, class = "b", coef = meef, 
                 nlpar = nlpar, suffix = paste0("me", p)),
      collapse("  Xme", pK, " ~ normal(Xn", pK, ", noise", pK,"); \n")
    )
  }
  out
}

stan_eta_fe <- function(fixef, center_X = TRUE, 
                        sparse = FALSE, nlpar = "") {
  # define Stan code to compute the fixef part of eta
  # Args:
  #   fixef: names of the population-level effects
  #   center_X: use the centered design matrix?
  #   sparse: use sparse matrix multiplication?
  #   nlpar: optional name of a non-linear parameter
  p <- usc(nlpar)
  if (length(fixef)) {
    if (sparse) {
      stopifnot(!center_X, nchar(nlpar) == 0L)
      eta_fe <- "csr_matrix_times_vector(rows(X), cols(X), wX, vX, uX, b)"
    } else {
      eta_fe <- paste0("X", if (center_X) "c", p, " * b", p)
    }
  } else { 
    eta_fe <- "rep_vector(0, N)"
  }
  eta_fe
}

stan_eta_re <- function(ranef, nlpar = "") {
  # write the group-level part of the linear predictor
  # Args:
  #   ranef: a named list returned by tidy_ranef
  #   nlpar: optional name of a non-linear parameter
  eta_re <- ""
  ranef <- ranef[ranef$nlpar == nlpar & !nzchar(ranef$type), ]
  for (id in unique(ranef$id)) {
    r <- ranef[ranef$id == id, ]
    idp <- paste0(r$id, usc(r$nlpar, "prefix"))
    eta_re <- paste0(eta_re, collapse(
      " + (", stan_eta_r(r), ") * Z_", idp, "_", r$cn, "[n]"))
  }
  eta_re
}

stan_eta_r <- function(r) {
  # Stan code for r parameters in linear predictor terms
  # Args:
  #   r: data.frame created by tidy_ranef
  # Returns:
  #   A character vector, one element per row of 'r' 
  stopifnot(nrow(r) > 0L, length(unique(r$gtype)) == 1L)
  idp <- paste0(r$id, usc(r$nlpar, "prefix"))
  if (r$gtype[1] == "mm") {
    ng <- seq_along(r$gcall[[1]]$groups)
    out <- rep("", nrow(r))
    for (i in seq_along(out)) {
      out[i] <- paste0(
        "W_", r$id[i], "_", ng, "[n] * ", 
        "r_", idp[i], "_", r$cn[i], "[J_", r$id[i], "_", ng, "[n]]",
        collapse = " + ") 
    }
  } else {
    out <- paste0("r_", idp, "_", r$cn, "[J_", r$id, "[n]]")
  }
  out
}

stan_eta_mo <- function(monef, ranef, nlpar = "") {
  # write the linear predictor for monotonic effects
  # Args:
  #   monef: names of the monotonic effects
  #   nlpar: an optional character string to add to the varnames
  #         (used for non-linear models)
  p <- usc(nlpar)
  eta_mo <- ""
  ranef <- ranef[ranef$nlpar == nlpar & ranef$type == "mo", ]
  invalid_coef <- setdiff(ranef$coef, monef)
  if (length(invalid_coef)) {
    stop2("Monotonic group-level terms require ", 
          "corresponding population-level terms.")
  }
  for (i in seq_along(monef)) {
    r <- ranef[ranef$coef == monef[i], ]
    if (nrow(r)) {
      rpars <- paste0(" + ", stan_eta_r(r))
    } else {
      rpars <- ""
    }
    eta_mo <- paste0(eta_mo,
      " + (bmo", p, "[", i, "]", rpars, ") * monotonic(",
      "simplex", p, "_", i, ", Xmo", p, "[n, ", i, "])"
    )
  }
  eta_mo
}

stan_eta_sm <- function(smooths, nlpar = "") {
  # write the linear predictor for smooth terms
  # Args:
  #   smooths: names of the smooth terms
  #   nlpar: an optional character string to add to the varnames
  #         (used for non-linear models)
  p <- usc(nlpar)
  eta_smooths <- ""
  if (length(smooths)) {
    stopifnot(!is.null(attr(smooths, "nbases")))
    for (i in seq_along(smooths)) {
      pi <- paste0(p, "_", i)
      nb <- seq_len(attr(smooths, "nbases")[[smooths[i]]])
      eta_smooths <- paste0(eta_smooths, collapse(
        " + Zs", pi, "_", nb, " * s", pi, "_", nb)
      )
    }
  }
  eta_smooths
}

stan_eta_bsts <- function(autocor) {
  # write the linear predictor for bsts terms
  # Args:
  #   autocor: object of class cor_brms
  eta_bsts <- ""
  if (is(autocor, "cor_bsts")) {
    eta_bsts <- " + loclev[n]"
  }
  eta_bsts
}

stan_eta_transform <- function(family, llh_adj = FALSE) {
  # indicate whether eta needs to be transformed
  # manually using the link functions
  # Args:
  #   family: a list with elements 'family' and 'link'
  #   llh_adj: is the model censored or truncated?
  stopifnot(all(c("family", "link") %in% names(family)))
  link <- family$link
  !(!is_skewed(family) && link == "identity" ||
    is_ordinal(family) || is_categorical(family) ||
    is_zero_inflated(family) || is_hurdle(family)) &&
  (llh_adj || !stan_has_built_in_fun(family))
}

stan_eta_ilink <- function(family, auxpars = NULL, adforms = NULL) {
  # correctly apply inverse link to eta
  # Args:
  #   family: a list with elements 'family' and 'link
  #   auxpars: names of auxiliary parameters
  #   adforms: list of formulas containing addition terms
  stopifnot(all(c("family", "link") %in% names(family)))
  llh_adj <- stan_llh_adj(adforms, c("cens", "trunc"))
  if (stan_eta_transform(family, llh_adj = llh_adj)) {
    link <- family$link
    family <- family$family
    shape <- ifelse(
      is.formula(adforms$disp), "disp_shape[n]", 
      ifelse("shape" %in% auxpars, "shape[n]", "shape")
    )
    nu <- ifelse("nu" %in% auxpars, "nu[n]", "nu")
    fl <- ifelse(family %in% c("gamma", "exponential"), 
                 paste0(family, "_", link), family)
    ilink <- stan_ilink(link)
    out <- switch(fl,
      c(paste0(ilink, "("), ")"),
      gamma_log = c(paste0(shape, " * exp(-("), "))"),
      gamma_inverse = c(paste0(shape, " * ("), ")"),
      gamma_identity = c(paste0(shape, " / ("), ")"),
      exponential_log = c("exp(-(", "))"),
      exponential_inverse = c("(", ")"),
      exponential_identity = c("inv(", ")"),
      weibull = c(
        paste0(ilink, "(("), 
        paste0(") / ", shape, ")")
      ),
      frechet = c(
        paste0(ilink, "("),
        paste0(") / tgamma(1 - 1 / ", nu, ")")
      )
    )
  } else {
    out <- rep("", 2)
  }
  out
}
