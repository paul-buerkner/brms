stan_effects <- function(effects, data, family = gaussian(),
                         center_X = TRUE, ranef = empty_ranef(), 
                         prior = brmsprior(), autocor = cor_arma(), 
                         threshold = "flexible", sparse = FALSE, 
                         nlpar = "", eta = "eta") {
  # combine effects for the predictors of a single (non-linear) parameter
  # Args:
  #   center_X: center population-level design matrix if possible?
  #   eta: prefix of the linear predictor variable
  p <- usc(nlpar, "prefix")
  if (nzchar(eta) && nzchar(nlpar)) {
    eta <- usc(eta, "suffix") 
  }
  eta <- paste0(eta, nlpar)
  ranef <- ranef[ranef$nlpar == nlpar, ]
  out <- list()
  out$modelD <- paste0("  vector[N] ", eta, "; \n")
  # include fixed effects
  center_X <- center_X && has_intercept(effects$fixed) && 
              !is(autocor, "cor_bsts") && !sparse
  rm_intercept <- center_X || is(autocor, "cor_bsts") || is_ordinal(family)
  cols2remove <- if (rm_intercept) "Intercept"
  fixef <- colnames(data_fixef(effects, data, autocor = autocor)$X)
  fixef <- setdiff(fixef, cols2remove)
  text_fixef <- stan_fixef(fixef, center_X = center_X, 
                           family = family, prior = prior, nlpar = nlpar,
                           sparse = sparse, threshold = threshold)
  # include spline terms
  splines <- get_spline_labels(effects, data = data)
  text_splines <- stan_splines(splines, prior = prior, nlpar = nlpar)
  # include category specific effects
  csef <- colnames(get_model_matrix(effects$cs, data))
  text_csef <- stan_csef(csef, ranef, prior = prior)
  # include monotonic effects
  monef <- all_terms(effects$mo)
  text_monef <- stan_monef(monef, ranef, prior = prior, nlpar = nlpar)
  # include measurement error variables
  meef <- get_me_labels(effects, data = data)
  text_meef <- stan_meef(meef, ranef, prior = prior, nlpar = nlpar)
  out <- collapse_lists(list(out, text_fixef, text_csef, 
                             text_monef, text_meef, text_splines))
  
  has_offset <- !is.null(get_offset(effects$fixed))
  if (has_offset) {
    out$data <- paste0(out$data, "  vector[N] offset", p, "; \n")
  }
  
  # initialize and compute eta_<nlpar>
  out$modelC1 <- paste0(
    out$modelC1, "  ", eta, " = ", 
    text_fixef$eta, text_splines$eta,
    if (center_X && !is_ordinal(family)) 
      paste0(" + temp", p, "_Intercept"),
    if (has_offset) paste0(" + offset", p),
    if (get_arr(autocor)) " + Yarr * arr", 
    "; \n")
  
  # repare loop over eta
  eta_ma <- ifelse(get_ma(autocor) && !use_cov(autocor),
                   paste0(" + head(E", p, "[n], Kma) * ma"), "")
  eta_loop <- paste0(
    stan_eta_ranef(ranef, nlpar = nlpar),
    text_monef$eta, text_meef$eta,
    eta_ma, stan_eta_bsts(autocor))
  if (nzchar(eta_loop)) {
    out$modelC2 <- paste0(out$modelC2,
      "    ", eta, "[n] = ", eta, "[n]", eta_loop, "; \n")
  }
  # include autoregressive effects
  if (get_ar(autocor) && !use_cov(autocor)) {
    eta_ar <- paste0(eta, "[n] + head(E", p, "[n], Kar) * ar")
    out$modelC3 <- paste0(out$modelC3, 
      "    ", eta, "[n] = ", eta_ar, "; \n")
  }
  # possibly transform eta before it is passed to the likelihood
  eta_ilink <- stan_eta_ilink(family$family, family$link, effects)
  if (sum(nzchar(eta_ilink))) {
    eta_ilink <- paste0(eta_ilink[1], eta, "[n]", eta_ilink[2])
    out$modelC3 <- paste0(out$modelC3, 
      "    ", eta, "[n] = ", eta_ilink, "; \n")
  }
  out
}

stan_effects_mv <- function(effects, data, family = gaussian(), 
                            ranef = empty_ranef(), prior = brmsprior(), 
                            autocor = cor_arma(), sparse = FALSE) {
  if (sparse) {
    stop2("Sparse design matrices are not yet implemented ", 
          "for multivariate models.")
  }
  out <- list()
  resp <- effects$response
  if (length(resp) > 1L) {
    args <- nlist(effects, data, family, ranef, prior, autocor)
    resp <- effects$response
    tmp_list <- named_list(resp)
    for (r in resp) {
      tmp_list[[r]] <- do.call(stan_effects, c(args, nlpar = r))
    }
    out <- collapse_lists(tmp_list)
    if (is_linear(family)) {
      len_Eta_n <- "nresp" 
    } else if (is_categorical(family)) {
      len_Eta_n <- "ncat - 1"
    } else {
      stop2("Multivariate models are not yet implemented ", 
            "for family '", family$family, "'.")
    }
    out$modelD <- paste0(out$modelD, 
      "  // multivariate linear predictor matrix \n",
      "  vector[", len_Eta_n, "] Eta[N]; \n")
    out$modelC3 <- paste0(out$modelC3, collapse(
      "    Eta[n, ", seq_along(resp), "] = eta_", resp, "[n]; \n")) 
  }
  out
}

stan_nonlinear <- function(effects, data, family = gaussian(), 
                           ranef = empty_ranef(), prior = brmsprior()) {
  # prepare Stan code for non-linear models
  # Args:
  #   effects: a list returned by extract_effects()
  #   data: data.frame supplied by the user
  #   family: the model family
  #   cov_ranef: a list of user-defined covariance matrices
  #   prior: a brmsprior object
  out <- list()
  if (length(effects$nlpars)) {
    for (i in seq_along(effects$nlpars)) {
      nlpar <- names(effects$nlpars)[i]
      # do not pass 'family' here to avoid inverse link transformations
      nl_text <- stan_effects(effects = effects$nlpars[[i]],
                              data = data, ranef = ranef, 
                              prior = prior, nlpar = nlpar, 
                              center_X = FALSE)
      out <- collapse_lists(list(out, nl_text))
    }
    # prepare non-linear model of eta 
    nlpars <- wsp(names(effects$nlpars))
    new_nlpars <- paste0(" eta_", names(effects$nlpars), "[n] ")
    # covariates in the non-linear model
    covars <- wsp(setdiff(all.vars(effects$fixed[[3]]), 
                          names(effects$nlpars)))
    if (length(covars)) {
      out$data <- paste0(out$data, 
        "  int<lower=1> KC;  // number of covariates \n",
        "  matrix[N, KC] C;  // covariate matrix \n")
      new_covars <- paste0(" C[n, ", seq_along(covars), "] ")
    } else {
      new_covars <- NULL
    }
    # add whitespaces to be able to replace parameters and covariates
    meta_sym <- c("+", "-", "*", "/", "^", ")", "(", ",")
    nlmodel <- gsub(" ", "", collapse(deparse(effects$fixed[[3]])))
    nlmodel <- wsp(rename(nlmodel, meta_sym, wsp(meta_sym))) 
    nlmodel <- rename(nlmodel, c(nlpars, covars, " ( ", " ) "), 
                      c(new_nlpars, new_covars, "(", ")"))
    # possibly transform eta in the transformed params block
    eta_ilink <- stan_eta_ilink(family$family, family$link, effects)
    out$modelD <- paste0(out$modelD, "  vector[N] eta; \n")
    out$modelC3 <- paste0(out$modelC3, 
      "    // compute non-linear predictor \n",
      "    eta[n] = ", eta_ilink[1], trimws(nlmodel), eta_ilink[2], "; \n")
  }
  out
}

stan_auxpars <- function(effects, data, family = gaussian(),
                         ranef = empty_ranef(), prior = brmsprior(), 
                         autocor = cor_arma()) {
  # Stan code for auxiliary parameters
  # Args:
  #   effects: output of extract_effects
  #   other arguments: same as make_stancode
  out <- list()
  default_defs <- c(
    sigma = "  real<lower=0> sigma;  // residual SD \n",
    shape = "  real<lower=0> shape;  // shape parameter \n",
    nu = "  real<lower=1> nu;  // degrees of freedom \n",
    phi = "  real<lower=0> phi;  // precision parameter \n",
    kappa = "  real<lower=0> kappa;  // precision parameter \n",
    beta = "  real<lower=0> beta;  // scale parameter \n",
    zi = "  real<lower=0,upper=1> zi;  // zero-inflation probability \n", 
    hu = "  real<lower=0,upper=1> hu;  // hurdle probability \n",
    bs = "  real<lower=0> bs;  // boundary separation parameter \n",
    ndt = "  real<lower=0,upper=min_Y> ndt;  // non-decision time parameter \n",
    bias = "  real<lower=0,upper=1> bias;  // initial bias parameter \n",
    disc = "")
  valid_auxpars <- valid_auxpars(family, effects, autocor = autocor)
  # don't supply the family argument to avoid applying link functions
  args <- nlist(data, ranef, center_X = FALSE, eta = "")
  for (ap in valid_auxpars) {
    if (!is.null(effects$auxpars[[ap]])) {
      ap_ilink <- ilink_auxpars(ap, stan = TRUE)
      ap_prior <- prior[prior$nlpar == ap, ]
      ap_args <- list(effects = effects$auxpars[[ap]], 
                      nlpar = ap, prior = ap_prior)
      if (nzchar(ap_ilink)) {
        ap_ilink <- paste0("    ", ap, "[n] = ", ap_ilink, "(", ap, "[n]); \n")
      } else {
        ap_ilink <- ""
      }
      out[[ap]] <- do.call(stan_effects, c(ap_args, args))
      out[[ap]]$modelC3 <- paste0(out[[ap]]$modelC3, ap_ilink)
    } else {
      out[[ap]] <- list(par = default_defs[ap],
        prior = stan_prior(class = ap, prior = prior))
    }
  }
  collapse_lists(out)
} 

stan_fixef <- function(fixef, center_X = TRUE, family = gaussian(),
                       prior = brmsprior(), nlpar = "", sparse = FALSE, 
                       threshold = "flexible") {
  # Stan code for fixed effects
  # Args:
  #   fixef: names of the fixed effects
  #   center_X: center the design matrix?
  #   family: the model family
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior 
  #   threshold: either "flexible" or "equidistant" 
  # Returns:
  #   a list containing Stan code related to fixed effects
  p <- usc(nlpar, "prefix")
  ct <- ifelse(center_X, "c", "")
  out <- list()
  if (length(fixef)) {
    out$data <- paste0(out$data, 
      "  int<lower=1> K", p, ";",
      "  // number of population-level effects \n", 
      "  matrix[N, K", p, "] X", p, ";",
      "  // population-level design matrix \n")
    if (sparse) {
      stopifnot(!center_X)
      if (nchar(nlpar)) {
         stop2("Sparse matrices are not yet implemented for this model.")
      }
      out$tdataD <- "  #include tdataD_sparse_X.stan \n"
      out$tdataC <- "  #include tdataC_sparse_X.stan \n"
    }
    bound <- get_bound(prior, class = "b", nlpar = nlpar)
    out$par <- paste0(out$par,
      "  vector", bound, "[K", ct, p, "] b", p, ";",
      "  // population-level effects \n",
      if (!is.null(attr(prior, "hs_df")))
        paste0("  // horseshoe shrinkage parameters \n",
               "  vector<lower=0>[K", ct, "] hs_local; \n",
               "  real<lower=0> hs_global; \n"),
      if (!is.null(attr(prior, "lasso_df")))
        paste0("  // lasso shrinkage parameter \n",
               "  real<lower=0> lasso_inv_lambda; \n")
    ) 
    fixef_prior <- stan_prior(class = "b", coef = fixef, prior = prior,
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
        "  // column means of X", p, " before centering \n")
      out$tdataC <- paste0(out$tdataC, 
        "  Kc", p, " = K", p, " - 1;",
        "  // the intercept is removed from the design matrix \n",
        "  for (i in 2:K", p, ") { \n",
        "    means_X", p, "[i - 1] = mean(X", p, "[, i]); \n",
        "    Xc", p, "[, i - 1] = X", p, "[, i] - means_X", p, "[i - 1]; \n",
        "  } \n")
      # cumulative and sratio models are parameterized as thres - eta
      use_plus <- family$family %in% c("cumulative", "sratio")
      sub_X_means <- paste0(ifelse(use_plus, " + ", " - "), 
                            "dot_product(means_X", p, ", b", p, ")")
    } else {
      sub_X_means <- ""
    }
    if (is_ordinal(family)) {
      # temp intercepts for ordinal models are defined in stan_ordinal
      out$genD <- "  vector[ncat - 1] b_Intercept;  // thresholds \n" 
      out$genC <- paste0("  b_Intercept = temp_Intercept", sub_X_means, "; \n")
    } else {
      out$par <- paste0(out$par, 
        "  real temp", p, "_Intercept;  // temporary intercept \n")
      out$genD <- paste0(
        "  real b", p, "_Intercept;  // population-level intercept \n")
      out$genC <- paste0(
        "  b", p, "_Intercept = temp", p, "_Intercept", sub_X_means, "; \n")
    }
    # for equidistant thresholds only temp_Intercept1 is a parameter
    suffix <- paste0(p, "_Intercept")
    suffix <- paste0(suffix, ifelse(threshold == "equidistant", "1", ""))
    int_prior <- stan_prior("temp", prior = prior, coef = "Intercept",
                            suffix = suffix, nlpar = nlpar)
    out$prior <- paste0(out$prior, int_prior)
  }
  out$eta <- stan_eta_fixef(fixef, center_X, sparse, nlpar)
  out
}

stan_ranef <- function(id, ranef, prior = brmsprior(), 
                       cov_ranef = NULL) {
  # group-level effects in Stan 
  # Args:
  #   i: the index of the grouping factor
  #   ranef: a data.frame returned by tidy_ranef
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior
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
    if (r$gtype[1] == "mm") collapse(
      "  int<lower=1> J_", id, "_", ng, "[N]; \n",
      "  real<lower=0> W_", id, "_", ng, "[N]; \n")
    else paste0(
      "  int<lower=1> J_", id, "[N]; \n"),
    "  int<lower=1> N_", id, "; \n",
    "  int<lower=1> M_", id, "; \n",
    if (ccov) paste0(
      "  // cholesky factor of known covariance matrix \n",
      "  matrix[N_", id, ", N_", id,"] Lcov_", id,"; \n"))
  out$prior <- stan_prior(class = "sd", group = r$group[1], 
                          coef = r$coef, nlpar = r$nlpar, 
                          suffix = paste0("_", id), prior = prior)
  J <- seq_len(nrow(r))
  has_def_type <- !r$type %in% c("mo", "me")
  if (any(has_def_type)) {
    out$data <- paste0(out$data, 
      collapse("  vector[N] Z_", idp[has_def_type], 
               "_", r$cn[has_def_type], "; \n")) 
  }
  out$par <- paste0(
    "  vector<lower=0>[M_", id, "] sd_", id, ";",
    "  // group-level standard deviations \n")
  if (nrow(r) > 1L && r$cor[1]) {
    # multiple correlated group-level effects
    out$data <- paste0(out$data, "  int<lower=1> NC_", id, "; \n")
    out$par <- paste0(out$par,
      "  matrix[M_", id, ", N_", id, "] z_", id, ";",
      "  // unscaled group-level effects \n",    
      "  // cholesky factor of correlation matrix \n",
      "  cholesky_factor_corr[M_", id, "] L_", id, "; \n")
    out$prior <- paste0(out$prior, 
      stan_prior(class = "L", group = r$group[1],
                 suffix = paste0("_", id), prior = prior),
      "  to_vector(z_", id, ") ~ normal(0, 1); \n")
    out$transD <- paste0(
      "  // group-level effects \n",
      "  matrix[N_", id, ", M_", id, "] r_", id, "; \n",
      collapse("  vector[N_", id, "] r_", idp, "_", r$cn, "; \n"))
    if (ccov) {  # customized covariance matrix supplied
      out$transC1 <- paste0(
        "  r_", id," = as_matrix(kronecker(Lcov_", id, ",", 
        " diag_pre_multiply(sd_", id,", L_", id,")) *",
        " to_vector(z_", id, "), N_", id, ", M_", id, "); \n")
    } else { 
      out$transC1 <- paste0("  r_", id, " = ", 
        "(diag_pre_multiply(sd_", id, ", L_", id,") * z_", id, ")'; \n")
    }
    out$transC1 <- paste0(out$transC1, 
      collapse("  r_", idp, "_", r$cn, " = r_", id, "[, ", J, "];  \n"))
    # return correlations above the diagonal only
    cors_genC <- ulapply(2:nrow(r), function(k) 
      lapply(1:(k - 1), function(j) paste0(
        "  cor_", id, "[", (k - 1) * (k - 2) / 2 + j, 
        "] = Cor_", id, "[", j, ",", k, "]; \n")))
    out$genD <- paste0(
      "  corr_matrix[M_", id, "] Cor_", id, "; \n",
      "  vector<lower=-1,upper=1>[NC_", id, "] cor_", id, "; \n")
    out$genC <- paste0(
      "  // take only relevant parts of correlation matrix \n",
      "  Cor_", id, " = multiply_lower_tri_self_transpose(L_", id, "); \n",
      collapse(cors_genC)) 
  } else {
    # single or uncorrelated group-level effects
    out$par <- paste0(out$par,
      "  vector[N_", id, "] z_", id, "[M_", id, "];",
      "  // unscaled group-level effects \n")
    out$prior <- paste0(out$prior, collapse(
      "  z_", id, "[", 1:nrow(r), "] ~ normal(0, 1); \n"))
    out$transD <- paste0("  // group-level effects \n", 
      collapse("  vector[N_", id, "] r_", idp, "_", r$cn, "; \n"))
    out$transC1 <- collapse(
      "  r_", idp, "_", r$cn, " = sd_", id, "[", J, "] * (", 
      if (ccov) paste0("Lcov_", id, " * "), "z_", id, "[", J, "]); \n")
  }
  out
}

stan_splines <- function(splines, prior = brmsprior(), nlpar = "") {
  # Stan code of spline terms for GAMMs
  # Args:
  #   splines: names of the spline terms
  #   prior: object of class brmsprior
  #   nlpar: optional name of a non-linear parameter
  # Returns:
  #   A list of strings containing Stan code
  out <- list()
  p <- usc(nlpar)
  if (length(splines)) {
    stopifnot(!is.null(attr(splines, "nbases")))
    for (i in seq_along(splines)) {
      pi <- paste0(p, "_", i)
      nb <- seq_len(attr(splines, "nbases")[[i]])
      out$data <- paste0(out$data,
        "  // data of spline ", splines[i], "\n",  
        "  int nb", pi, ";  // number of bases \n",
        "  int knots", pi, "[nb", pi, "]; \n")
      out$data <- paste0(out$data, collapse(
        "  matrix[N, knots", pi, "[", nb, "]]", 
        " Zs", pi, "_", nb, "; \n"))
      out$par <- paste0(out$par,
        "  // parameters of spline ", splines[i], "\n")
      out$par <- paste0(out$par, collapse(
        "  vector[knots", pi, "[", nb, "]] zs", pi,"_", nb, "; \n",
        "  real<lower=0> sds", pi, "_", nb, "; \n"))
      out$transD <- paste0(out$transD, collapse(
        "  vector[knots", pi, "[", nb, "]] s", pi, "_", nb, "; \n"))
      out$transC1 <- paste0(out$transC1, collapse(
        "  s", pi, "_", nb, " = sds", pi,  "_", nb, 
        " * zs", pi, "_", nb, "; \n"))
      out$prior <- paste0(out$prior, collapse(
        "  zs", pi, "_", nb, " ~ normal(0, 1); \n"),
        stan_prior(class = "sds", coef = splines[i], 
                   nlpar = nlpar, prior = prior,
                   suffix = paste0(pi, "_", nb)))
    }
    out$eta <- stan_eta_splines(splines, nlpar = nlpar)
  }
  out
}

stan_monef <- function(monef, ranef = empty_ranef(),
                       prior = brmsprior(), nlpar = "") {
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
      collapse("  vector[Jmo", p, "[", I, "]]", 
               " con_simplex", p, "_", I, "; \n"))
    bound <- get_bound(prior, class = "b", nlpar = nlpar)
    out$par <- paste0(
      "  // monotonic effects \n", 
      "  vector", bound, "[Kmo", p, "] bmo", p, "; \n",
      collapse("  simplex[Jmo", p, "[", I, "]]", 
               " simplex", p, "_", I, "; \n")) 
    out$prior <- paste0(
      stan_prior(class = "b", coef = monef, prior = prior, 
                 nlpar = nlpar, suffix = paste0("mo", p)),
      collapse("  simplex", p, "_", I, 
               " ~ dirichlet(con_simplex", p, "_", I, "); \n"))
    out$eta <- stan_eta_monef(monef, ranef = ranef, nlpar = nlpar)
  }
  out
}

stan_csef <- function(csef, ranef = empty_ranef(), 
                      prior = brmsprior(), nlpar = "") {
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
      "  matrix[N, ncat - 1] etacs; \n")
  }
  if (length(csef)) {
    out$data <- paste0(
      "  int<lower=1> Kcs;  // number of category specific effects \n",
      "  matrix[N, Kcs] Xcs;  // category specific design matrix \n")
    bound <- get_bound(prior, class = "b")
    out$par <- paste0(
      "  matrix", bound, "[Kcs, ncat - 1] bcs;",
      "  // category specific effects \n")
    out$modelC1 <- "  etacs = Xcs * bcs; \n"
    out$prior <- stan_prior(class = "b", coef = csef, prior = prior, 
                            suffix = "cs", matrix = TRUE)
  } 
  if (nrow(ranef)) {
    if (!length(csef)) {
      # only group-level category specific effects present
      out$modelC1 <- "  etacs = rep_matrix(0, N, ncat - 1); \n"
    }
    cats <- get_matches("\\[[[:digit:]]+\\]$", ranef$coef)
    ncatM1 <- max(as.numeric(substr(cats, 2, nchar(cats) - 1)))
    for (i in seq_len(ncatM1)) {
      r_cat <- ranef[grepl(paste0("\\[", i, "\\]$"), ranef$coef), ]
      out$modelC2 <- paste0(out$modelC2,
        "    etacs[n, ", i, "] = etacs[n, ", i, "]")
      for (id in unique(r_cat$id)) {
        r <- r_cat[r_cat$id == id, ]
        idp <- paste0(r$id, usc(r$nlpar, "prefix"))
        out$modelC2 <- paste0(out$modelC2, collapse(
          " + r_", idp, "_", r$cn, "[J_", r$id, "[n]]",
          " * Z_", idp, "_", r$cn, "[n]"))
      }
      out$modelC2 <- paste0(out$modelC2, "; \n")
    }
  }
  out
} 

stan_meef <- function(meef, ranef = empty_ranef(),
                      prior = empty_brmsprior(), nlpar = "") {
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
        " + (bme", p, "[", i, "]", rpars, ") * ", meef_terms[i], covars[i])
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
        collapse("  vector[N] Cme", p, "_", seq_len(ncovars), "; \n")))
    out$par <- paste0(
      "  // noise free variables \n",
      collapse("  vector[N] Xme", pK, "; \n"),  
      "  vector[Kme", p, "] bme", p, ";",
      "  // coefficients of noise-free terms \n")
    out$prior <- paste0(
      stan_prior(class = "b", coef = meef, prior = prior, 
                 nlpar = nlpar, suffix = paste0("me", p)),
      collapse("  Xme", pK, " ~ normal(Xn", pK, ", noise", pK,"); \n"))
  }
  out
}

stan_eta_fixef <- function(fixef, center_X = TRUE, 
                           sparse = FALSE, nlpar = "") {
  # define Stan code to compute the fixef part of eta
  # Args:
  #   fixef: names of the fixed effects
  #   center_X: use the centered design matrix?
  #   sparse: use sparse matrix multiplication?
  #   nlpar: optional name of a non-linear parameter
  p <- usc(nlpar)
  if (length(fixef)) {
    if (sparse) {
      stopifnot(!center_X, nchar(nlpar) == 0L)
      eta_fixef <- "csr_matrix_times_vector(rows(X), cols(X), wX, vX, uX, b)"
    } else {
      eta_fixef <- paste0("X", if (center_X) "c", p, " * b", p)
    }
  } else { 
    eta_fixef <- "rep_vector(0, N)"
  }
  eta_fixef
}

stan_eta_ranef <- function(ranef, nlpar = "") {
  # write the group-level part of the linear predictor
  # Args:
  #   ranef: a named list returned by tidy_ranef
  #   nlpar: optional name of a non-linear parameter
  eta_ranef <- ""
  ranef <- ranef[ranef$nlpar == nlpar & !nzchar(ranef$type), ]
  for (id in unique(ranef$id)) {
    r <- ranef[ranef$id == id, ]
    idp <- paste0(r$id, usc(r$nlpar, "prefix"))
    eta_ranef <- paste0(eta_ranef, collapse(
      " + (", stan_eta_r(r), ") * Z_", idp, "_", r$cn, "[n]"))
  }
  eta_ranef
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

stan_eta_monef <- function(monef, ranef = empty_ranef(), nlpar = "") {
  # write the linear predictor for monotonic effects
  # Args:
  #   monef: names of the monotonic effects
  #   nlpar: an optional character string to add to the varnames
  #         (used for non-linear models)
  p <- usc(nlpar)
  eta_monef <- ""
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
    eta_monef <- paste0(eta_monef,
      " + (bmo", p, "[", i, "]", rpars, ") * monotonic(",
      "simplex", p, "_", i, ", Xmo", p, "[n, ", i, "])")
  }
  eta_monef
}

stan_eta_splines <- function(splines, nlpar = "") {
  # write the linear predictor for spline terms
  # Args:
  #   splines: names of the spline terms
  #   nlpar: an optional character string to add to the varnames
  #         (used for non-linear models)
  p <- usc(nlpar)
  eta_splines <- ""
  if (length(splines)) {
    stopifnot(!is.null(attr(splines, "nbases")))
    for (i in seq_along(splines)) {
      pi <- paste0(p, "_", i)
      nb <- seq_len(attr(splines, "nbases")[[splines[i]]])
      eta_splines <- paste0(eta_splines, collapse(
        " + Zs", pi, "_", nb, " * s", pi, "_", nb))
    }
  }
  eta_splines
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

stan_eta_transform <- function(family, link, llh_adj = FALSE) {
  # indicate whether eta needs to be transformed
  # in the transformed parameters block
  # Args:
  #   llh_adj: is the model censored or truncated?
  !(!is_skewed(family) && link == "identity" ||
    is_ordinal(family) || is_categorical(family) ||
    is_zero_inflated(family) || is_hurdle(family)) &&
  (llh_adj || !stan_has_built_in_fun(family, link))
}

stan_eta_ilink <- function(family, link, effects) {
  # correctly apply inverse link to eta
  # Args:
  #   family: string naming the family
  #   link: string naming the link function
  #   effects: output of extract_effects
  llh_adj <- stan_llh_adj(effects, c("cens", "trunc"))
  if (stan_eta_transform(family, link, llh_adj = llh_adj)) {
    ilink <- stan_ilink(link)
    shape <- ifelse(is.formula(effects$disp), "disp_shape[n]", 
                    ifelse("shape" %in% names(effects$auxpars), 
                           "shape[n]", "shape"))
    fl <- ifelse(family %in% c("gamma", "exponential"), 
                 paste0(family, "_", link), family)
    out <- switch(fl, 
      c(paste0(ilink, "("), ")"),
      gamma_log = c(paste0(shape, " * exp(-("), "))"),
      gamma_inverse = c(paste0(shape, " * ("), ")"),
      gamma_identity = c(paste0(shape, " / ("), ")"),
      exponential_log = c("exp(-(", "))"),
      exponential_inverse = c("(", ")"),
      exponential_identity = c("inv(", ")"),
      weibull = c(paste0(ilink, "(("), 
                  paste0(") / ", shape, ")")))
  } else {
    out <- rep("", 2)
  }
  out
}
