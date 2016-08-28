stan_effects <- function(effects, data, family = gaussian(),
                         rm_intercept = TRUE, ranef = empty_ranef(), 
                         prior = prior_frame(), autocor = cor_arma(), 
                         threshold = "flexible", sparse = FALSE, 
                         nlpar = "", eta = "eta") {
  # combine effects for the predictors of a single (non-linear) parameter
  # Args:
  #   rm_intercept: special treatment of FE intercept?
  #   eta: prefix of the linear predictor variable
  p <- usc(nlpar, "prefix")
  if (nzchar(eta) && nzchar(nlpar)) {
    eta <- usc(eta, "suffix") 
  }
  eta <- paste0(eta, nlpar)
  ranef <- ranef[ranef$nlpar == nlpar, ]
  out <- list()
  out$modelD <- paste0("  vector[N] ", eta, "; \n")
  # out$data <- paste0("  // data for effects of ", nlpar, "\n")
  # include fixed effects
  fixef <- colnames(data_fixef(effects, data = data, autocor = autocor,
                               rm_intercept = rm_intercept)$X)
  # cor_bsts parameters replace the intercept
  has_intercept <- has_intercept(effects$fixed)
  temp_intercept <- has_intercept && rm_intercept && !is(autocor, "cor_bsts")
  text_fixef <- stan_fixef(fixef, temp_intercept = temp_intercept, 
                           family = family, prior = prior, nlpar = nlpar,
                           sparse = sparse, threshold = threshold)
  # include spline terms
  splines <- get_spline_labels(effects)
  text_splines <- stan_splines(splines, prior = prior, nlpar = nlpar)
  # category specific effects
  csef <- colnames(get_model_matrix(effects$cse, data))
  text_csef <- stan_csef(csef = csef, prior = prior, nlpar = nlpar)
  if (length(csef)) {
    out$modelD <- paste0(out$modelD, 
     "  // linear predictor for category specific effects \n",                  
     "  matrix[N, ncat - 1] etap; \n")
  }
  # include monotonic effects
  monef <- colnames(data_monef(effects, data)$Xm)
  text_monef <- stan_monef(monef, prior = prior, nlpar = nlpar)
  out <- collapse_lists(list(out, text_fixef, text_csef, 
                             text_monef, text_splines))
  
  has_offset <- !is.null(get_offset(effects$fixed))
  if (has_offset) {
    out$data <- paste0(out$data, "  vector[N] offset", p, "; \n")
  }
  
  # initialize eta_<nlpar>
  out$modelC1 <- paste0(out$modelC1, "  ", eta, " = ", 
    stan_eta_fixef(fixef, sparse = sparse, nlpar = nlpar), 
    stan_eta_splines(splines, nlpar = nlpar), 
    if (temp_intercept && !is.ordinal(family)) 
      paste0(" + temp", p, "_Intercept"),
    if (has_offset) paste0(" + offset", p),
    if (get_arr(autocor)) " + Yarr * arr", 
    "; \n", 
    if (length(csef)) "  etap = Xp * bp; \n")
  
  # repare loop over eta
  eta_ma <- ifelse(get_ma(autocor) && !use_cov(autocor),
                   paste0(" + head(E", p, "[n], Kma) * ma"), "")
  eta_loop <- paste0(
    stan_eta_ranef(ranef, nlpar = nlpar),
    stan_eta_monef(monef, nlpar = nlpar),
    eta_ma, stan_eta_bsts(autocor))
  if (nzchar(eta_loop)) {
    out$modelC2 <- paste0(out$modelC2,
      "    ", eta, "[n] = ", eta, "[n]", eta_loop, "; \n")
  }
  
  # possibly transform eta before it is passed to the likelihood
  ll_adj <- any(ulapply(effects[c("cens", "trunc")], is.formula))
  transform <- stan_eta_transform(family$family, family$link, ll_adj = ll_adj)
  if (transform) {
    eta_ilink <- stan_eta_ilink(family$family, family$link, 
                                disp = is.formula(effects$disp))
  } else {
    eta_ilink <- rep("", 2)
  }
  # include autoregressive effects
  eta_ar <- ifelse(get_ar(autocor) && !use_cov(autocor), 
                   paste0(" + head(E", p, "[n], Kar) * ar"), "")
  if (sum(nzchar(c(eta_ilink, eta_ar)))) {
    out$modelC3 <- paste0(out$modelC3, "    ", eta, "[n] = ", 
      eta_ilink[1], eta, "[n]", eta_ar, eta_ilink[2], "; \n")
  }
  out
}

stan_effects_mv <- function(effects, data, family = gaussian(), 
                            ranef = empty_ranef(), prior = prior_frame(), 
                            autocor = cor_arma()) {
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
    if (is.linear(family)) {
      len_Eta_n <- "nresp" 
    } else if (is.categorical(family)) {
      len_Eta_n <- "ncat - 1"
    } else {
      stop("Invalid multivariate model", call. = FALSE)
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
                           ranef = empty_ranef(), prior = prior_frame()) {
  # prepare Stan code for non-linear models
  # Args:
  #   effects: a list returned by extract_effects()
  #   data: data.frame supplied by the user
  #   family: the model family
  #   cov_ranef: a list of user-defined covariance matrices
  #   prior: a prior_frame object
  out <- list()
  if (length(effects$nonlinear)) {
    for (i in seq_along(effects$nonlinear)) {
      nlpar <- names(effects$nonlinear)[i]
      # do not pass 'family' here to avoid inverse link transformations
      nl_text <- stan_effects(effects = effects$nonlinear[[i]],
                              data = data, ranef = ranef, 
                              prior = prior, nlpar = nlpar, 
                              rm_intercept = FALSE)
      out <- collapse_lists(list(out, nl_text))
    }
    # prepare non-linear model of eta 
    nlpars <- wsp(names(effects$nonlinear))
    new_nlpars <- paste0(" eta_", names(effects$nonlinear), "[n] ")
    # covariates in the nonlinear model
    covars <- wsp(setdiff(all.vars(effects$fixed[[3]]), 
                          names(effects$nonlinear)))
    if (length(covars)) {
      out$data <- paste0(out$data, 
        "  int<lower=1> KC;  // number of covariates \n",
        "  matrix[N, KC] C;  // covariate matrix \n")
      new_covars <- paste0(" C[n, ", seq_along(covars), "] ")
    } else new_covars <- NULL
    # add whitespaces to be able to replace parameters and covariates
    meta_sym <- c("+", "-", "*", "/", "^", ")", "(", ",")
    nlmodel <- gsub(" ", "", collapse(deparse(effects$fixed[[3]])))
    nlmodel <- wsp(rename(nlmodel, meta_sym, wsp(meta_sym))) 
    nlmodel <- rename(nlmodel, c(nlpars, covars, " ( ", " ) "), 
                      c(new_nlpars, new_covars, "(", ")"))
    # possibly transform eta in the transformed params block
    ll_adj <- any(ulapply(effects[c("cens", "trunc")], is.formula))
    transform <- stan_eta_transform(family$family, family$link, ll_adj = ll_adj)
    if (transform) {
      eta_ilink <- stan_eta_ilink(family$family, family$link, 
                                  disp = is.formula(effects$disp))
    } else eta_ilink <- rep("", 2)
    out$modelD <- paste0(out$modelD, "  vector[N] eta; \n")
    out$modelC2 <- paste0(out$modelC2, 
      "    // compute non-linear predictor \n",
      "    eta[n] = ", eta_ilink[1], trimws(nlmodel), eta_ilink[2], "; \n")
  }
  out
}

stan_auxpars <- function(effects, data, family = gaussian(),
                         ranef = empty_ranef(), prior = prior_frame(), 
                         autocor = cor_arma()) {
  # Stan code for auxiliary parameters
  # Args:
  #   effects: output of extract_effects
  #   other arguments: same as make_stancode
  out <- list()
  default_defs <- c(
    sigma = "  real<lower=0> sigma;  // residual SD \n",
    shape = "  real<lower=0> shape;  // shape parameter \n",
    nu = "  real<lower=0> nu;  // degrees of freedom \n",
    phi = "  real<lower=0> phi;  // precision parameter \n",
    kappa = "  real<lower=0> kappa;  // precision parameter \n",
    zi = "  real<lower=0,upper=1> zi;  // zero-inflation probability \n", 
    hu = "  real<lower=0,upper=1> hu;  // hurdle probability \n")
  ilinks <- c(sigma = "exp", shape = "exp", nu = "exp", 
              phi = "exp", kappa = "exp", zi = "", hu = "") 
  valid_auxpars <- valid_auxpars(family, effects, autocor = autocor)
  args <- nlist(data, family, ranef, rm_intercept = FALSE, eta = "")
  for (ap in valid_auxpars) {
    if (!is.null(effects[[ap]])) {
      ap_prior <- prior[prior$nlpar == ap, ]
      ap_args <- list(effects = effects[[ap]], nlpar = ap, prior = ap_prior)
      if (nzchar(ilinks[ap])) {
        ap_ilink <- paste0("  ", ap, " = ", ilinks[ap], "(", ap, "); \n") 
      } else ap_ilink <- ""
      out[[ap]] <- c(do.call(stan_effects, c(ap_args, args)), 
                     list(modelC4 = ap_ilink))
    } else {
      out[[ap]] <- list(par = default_defs[ap],
        prior = stan_prior(class = ap, prior = prior))
    }
  }  
  collapse_lists(out)
} 

stan_fixef <- function(fixef, temp_intercept = TRUE, family = gaussian(),
                       prior = prior_frame(), nlpar = "", sparse = FALSE, 
                       threshold = "flexible") {
  # Stan code for fixed effects
  # Args:
  #   fixef: names of the fixed effects
  #   temp_intercept: is a temporary intercept required
  #                   due to centering of the design matrix?
  #   csef: names of the category specific effects
  #   family: the model family
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior 
  #   threshold: either "flexible" or "equidistant" 
  # Returns:
  #   a list containing Stan code related to fixed effects
  p <- usc(nlpar, "prefix")
  out <- list()
  if (length(fixef)) {
    centered <- ifelse(temp_intercept, "centered", "")
    out$data <- paste0(out$data, 
      "  int<lower=1> K", p, ";",
      "  // number of population-level effects \n", 
      "  matrix[N, K", p, "] X", p, ";",
      "  // ", centered, " population-level design matrix \n")
    if (sparse) {
      if (nchar(nlpar)) {
         stop("Sparse matrices are not yet implemented for this model.",
              call. = FALSE)
      }
      out$tdataD <- "  #include tdata_def_sparse_X.stan \n"
      out$tdataC <- "  #include tdata_calc_sparse_X.stan \n"
    }
    bound <- get_bound(prior, class = "b", nlpar = nlpar)
    out$par <- paste0(out$par,
      "  vector", bound, "[K", p, "] b", p, ";",
      "  // population-level effects \n",
      if (!is.null(attr(prior, "hs_df")))
        paste0("  // horseshoe shrinkage parameters \n",
               "  vector<lower=0>[K] hs_local; \n",
               "  real<lower=0> hs_global; \n")) 
    fixef_prior <- stan_prior(class = "b", coef = fixef, prior = prior,
                              nlpar = nlpar, suffix = p)
    out$prior <- paste0(out$prior, fixef_prior)
  }
  if (temp_intercept) {
    if (length(fixef)) {
      # X_means is not defined in intercept only models
      def_X_means <- paste0("  vector[K", p, "] X_means", p, ";", 
                            "  // column means of X before centering \n")
      sub_X_means <- paste0(" - dot_product(X_means", p, ", b", p, ")")
    } else {
      def_X_means <- sub_X_means <- ""
    }
    if (is.ordinal(family)) {
      # temp intercepts for ordinal models are defined in stan_ordinal
      out$data <- paste0(out$data, def_X_means)
      out$genD <- "  vector[ncat - 1] b_Intercept;  // thresholds \n" 
      out$genC <- paste0("  b_Intercept = temp_Intercept", sub_X_means, "; \n")
    } else {
      out$data <- paste0(out$data, def_X_means)
      out$par <- paste0(out$par, 
        "  real temp", p, "_Intercept;  // temporary Intercept \n")
      out$genD <- paste0(
        "  real b", p, "_Intercept;  // population-level intercept \n")
      out$genC <- paste0(
        "  b", p, "_Intercept = temp", p, "_Intercept", sub_X_means, "; \n")
    }
    # for equidistant thresholds only temp_Intercept1 is a parameter
    suffix <- ifelse(threshold == "equidistant", "1", "")
    int_prior <- stan_prior("temp_Intercept", #coef = intercepts, 
                            prior = prior, suffix = suffix)
    out$prior <- paste0(out$prior, int_prior)
  }
  out
}

stan_ranef <- function(id, ranef, prior = prior_frame(), 
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
  cor <- r$cor[1]
  ccov <- r$group[1] %in% names(cov_ranef)
  idp <- paste0(r$id, usc(r$nlpar, "prefix"))
  out <- list()
  out$data <- paste0(
    "  // data for group-specific effects of ID ", id, " \n",
    "  int<lower=1> J_", id, "[N]; \n",
    "  int<lower=1> N_", id, "; \n",
    "  int<lower=1> M_", id, "; \n",
    if (ccov) paste0(
      "  // cholesky factor of known covariance matrix \n",
      "  matrix[N_", id, ", N_", id,"] Lcov_", id,"; \n"))
  out$prior <- stan_prior(class = "sd", group = r$group[1], 
                          coef = r$coef, nlpar = r$nlpar, 
                          suffix = paste0("_", id), prior = prior)
  if (nrow(r) == 1L) {  # only one group-level effect
    out$data <- paste0(out$data, "  vector[N] Z_", idp, "_1; \n")
    out$par <- paste0(
      "  real<lower=0> sd_", id, ";",
      "  // group-specific standard deviation \n",
      "  vector[N_", id, "] z_", id, ";",
      "  // unscaled group-specific effects \n")
    out$prior <- paste0(out$prior,"  z_", id, " ~ normal(0, 1); \n")
    out$transD <- paste0(
      "  // group-specific effects \n",
      "  vector[N_", id, "] r_", idp, "_1; \n")
    out$transC1 <- paste0("  r_", idp, "_1 = sd_", id, " * (", 
      if (ccov) paste0("Lcov_", id, " * "), "z_", id, ");\n")
  } else if (nrow(r) > 1L) {
    J <- 1:nrow(r)
    out$data <- paste0(out$data, 
      collapse("  vector[N] Z_", idp, "_", r$cn, ";  \n"))
    out$par <- paste0(
      "  vector<lower=0>[M_", id, "] sd_", id, ";",
      "  // group-specific standard deviations \n")
    if (cor) {  # multiple correlated group-level effects
      out$data <- paste0(out$data, "  int<lower=1> NC_", id, "; \n")
      out$par <- paste0(out$par,
        "  matrix[M_", id, ", N_", id, "] z_", id, ";",
        "  // unscaled group-specific effects \n",    
        "  // cholesky factor of correlation matrix \n",
        "  cholesky_factor_corr[M_", id, "] L_", id, "; \n")
      out$prior <- paste0(out$prior, 
        stan_prior(class = "L", group = r$group[1],
                   suffix = paste0("_", id), prior = prior),
        "  to_vector(z_", id, ") ~ normal(0, 1); \n")
      out$transD <- paste0(
        "  // group-specific effects \n",
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
    } else {  # multiple uncorrelated group-level effects
      out$par <- paste0(out$par,
                        "  vector[N_", id, "] z_", id, "[M_", id, "];",
                        "  // unscaled group-specific effects \n")
      out$prior <- paste0(out$prior, collapse(
        "  z_", id, "[", 1:nrow(r), "] ~ normal(0, 1); \n"))
      out$transD <- paste0("  // group-specific effects \n", 
        collapse("  vector[N_", id, "] r_", idp, "_", r$cn, "; \n"))
      out$transC1 <- collapse(
        "  r_", idp, "_", r$cn, " = sd_", id, "[", J, "] * (", 
        if (ccov) paste0("Lcov_", id, " * "), "z_", id, "[", J, "]); \n")
    }
  }
  out
}

stan_splines <- function(splines, prior = prior_frame(), nlpar = "") {
  # Stan code of spline terms for GAMMs
  # Args:
  #   splines: names of the spline terms
  #   prior: object of class prior_frame
  #   nlpar: optional name of a non-linear parameter
  # Returns:
  #   A list of strings containing Stan code
  out <- list()
  p <- if (nchar(nlpar)) paste0("_", nlpar) else ""
  if (length(splines)) {
    out$data <- paste0(
      "  int ns", p, ";  // number of splines terms \n",
      "  int knots", p, "[ns", p, "];",
      "  // number of knots per spline \n")
  }
  for (i in seq_along(splines)) {
    pi <- paste0(p, "_", i)
    out$data <- paste0(out$data, 
      "  // design matrix of spline ", splines[i], "\n",  
      "  matrix[N, knots", p, "[", i, "]] Zs", pi, "; \n")
    out$par <- paste0(out$par,
      "  // parameters of spline ", splines[i], "\n", 
      "  vector[knots", p, "[", i, "]] zs", pi, "; \n",
      "  real<lower=0> sds", pi, "; \n")
    out$transD <- paste0(out$transD,
      "  vector[knots", p, "[", i, "]] s", pi, "; \n")
    out$transC1 <- paste0(out$transC1,
      "  s", pi, " = sds", pi, " * zs", pi, "; \n")
    out$prior <- paste0(out$prior, 
      "  zs", pi, " ~ normal(0, 1); \n",
      stan_prior(class = "sds", coef = splines[i], 
                 nlpar = nlpar, suffix = pi, prior = prior))
  }
  out
}

stan_monef <- function(monef, prior = prior_frame(), nlpar = "") {
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
      "  int<lower=1> Km", p, ";  // number of monotonic effects \n",
      "  int Xm", p, "[N, Km", p, "];  // monotonic design matrix \n",
      "  int<lower=2> Jm", p, "[Km", p, "];  // length of simplexes \n",
      collapse("  vector[Jm", p, "[", I, "]]", 
               " con_simplex", p, "_", I, "; \n"))
    bound <- get_bound(prior, class = "b", nlpar = nlpar)
    out$par <- paste0(
      "  // monotonic effects \n", 
      "  vector", bound, "[Km", p, "] bm", p, "; \n",
      collapse("  simplex[Jm", p, "[", I, "]]", 
               " simplex", p, "_", I, "; \n")) 
    out$prior <- paste0(
      stan_prior(class = "b", coef = monef, prior = prior, 
                 nlpar = nlpar, suffix = paste0("m", p)),
      collapse("  simplex", p, "_", I, 
               " ~ dirichlet(con_simplex", p, "_", I, "); \n"))
  }
  out
}

stan_csef <- function(csef, prior = prior_frame(), nlpar = "") {
  # Stan code for category specific effects
  # Args:
  #   csef: names of the category specific effects
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior
  # (!) Not implemented for non-linear models
  out <- list()
  if (length(csef)) {
    stopifnot(!nzchar(nlpar))
    out$data <- paste0(
      "  int<lower=1> Kp;  // number of category specific effects \n",
      "  matrix[N, Kp] Xp;  // CSE design matrix \n")
    bound <- get_bound(prior, class = "b")
    out$par <- paste0(
      "  matrix", bound, "[Kp, ncat - 1] bp;  // category specific effects \n")
    out$prior <- stan_prior(class = "b", coef = csef, prior = prior, 
                            suffix = "p", matrix = TRUE)
  }
  out
} 

stan_eta_fixef <- function(fixef, sparse = FALSE, nlpar = "") {
  # define Stan code to compute the fixef part of eta
  # Args:
  #   fixef: names of the fixed effects
  #   nlpar: an optional character string to add to the variable names
  #         (used for non-linear models)
  #   sparse: logical; use sparse matrix multiplication?
  p <- if (nchar(nlpar)) paste0("_", nlpar) else ""
  if (length(fixef)) {
    if (sparse) {
      stopifnot(nchar(nlpar) == 0L)
      eta_fixef <- "csr_matrix_times_vector(rows(X), cols(X), wX, vX, uX, b)"
    } else {
      eta_fixef <- paste0("X", p, " * b", p)
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
  #   nlpar: currently unused
  eta_ranef <- ""
  for (id in unique(ranef$id)) {
    r <- ranef[ranef$id == id, ]
    idp <- paste0(r$id, usc(r$nlpar, "prefix"))
    eta_ranef <- paste0(eta_ranef, collapse(
      " + r_", idp, "_", r$cn, "[J_", r$id, "[n]]",
      " * Z_", idp, "_", r$cn, "[n]"))
  }
  eta_ranef
}

stan_eta_monef <- function(monef, nlpar = "") {
  # write the linear predictor for monotonic effects
  # Args:
  #   monef: names of the monotonic effects
  #   nlpar: an optional character string to add to the varnames
  #         (used for non-linear models)
  p <- if (nchar(nlpar)) paste0("_", nlpar) else ""
  eta_monef <- ""
  for (i in seq_along(monef)) {
    eta_monef <- paste0(eta_monef,
      " + bm", p, "[", i, "] * monotonic(",
      "simplex", p, "_", i, ", Xm", p, "[n, ", i, "])")
  }
  eta_monef
}

stan_eta_splines <- function(splines, nlpar = "") {
  # write the linear predictor for spline terms
  # Args:
  #   splines: names of the spline terms
  #   nlpar: an optional character string to add to the varnames
  #         (used for non-linear models)
  p <- if (nchar(nlpar)) paste0("_", nlpar) else ""
  eta_splines <- ""
  for (i in seq_along(splines)) {
    eta_splines <- paste0(eta_splines, 
      " + Zs", p, "_", i, " * s", p, "_", i)
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

stan_eta_transform <- function(family, link, ll_adj = FALSE) {
  # indicate whether eta needs to be transformed
  # in the transformed parameters block
  # Args:
  #   ll_adj: is the model censored or truncated?
  !(!is.skewed(family) && link == "identity" ||
    is.ordinal(family) || is.categorical(family) ||
    is.zero_inflated(family) || is.hurdle(family)) &&
  (ll_adj || !stan_has_built_in_fun(family, link))
}

stan_eta_ilink <- function(family, link, disp = FALSE) {
  # correctly apply inverse link to eta
  ilink <- stan_ilink(link)
  shape <- ifelse(disp, "disp_shape[n]", "shape")
  fl <- ifelse(family %in% c("gamma", "exponential"), 
               paste0(family,"_",link), family)
  switch(fl, c(paste0(ilink,"("), ")"),
         gamma_log = c(paste0(shape, " * exp(-("), "))"),
         gamma_inverse = c(paste0(shape, " * ("), ")"),
         gamma_identity = c(paste0(shape, " / ("), ")"),
         exponential_log = c("exp(-(", "))"),
         exponential_inverse = c("(", ")"),
         exponential_identity = c("inv(", ")"),
         weibull = c(paste0(ilink,"(("), 
                     paste0(") / ", shape, ")")))
}
