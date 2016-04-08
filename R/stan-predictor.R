stan_linear <- function(effects, data, family = gaussian(), 
                        intercepts = "Intercept", prior = prior_frame(), 
                        autocor = cor_arma(), threshold = "flexible",
                        cov_ranef = NULL, sparse = FALSE) {
  # prepare Stan code for (generalized) linear models
  # Args:
  #   effects: a list returned by extract_effects()
  #   data: data.frame supplied by the user
  #   family: the model family
  #   intercepts: names of the population-level intercepts
  #   prior: a prior_frame object
  #   autocor: an object of class cor_arma
  #   threshold: used in ordinal models; 
  #              either "flexible" or "equidistant" 
  #   cov_ranef: a named list of user-defined covariance matrices
  #   sparse: logical; use sparse matrix multiplication for fixed effects?
  # 
  # Return:
  #   the linear predictor in stan language
  stopifnot(is(family, "family"))
  # generate code for different types of effects
  fixef <- colnames(get_model_matrix(rhs(effects$fixed), data = data, 
                                     forked = is.forked(family),
                                     cols2remove = intercepts))
  text_fixef <- stan_fixef(fixef = fixef, intercepts = intercepts, 
                           family = family, prior = prior, 
                           sparse = sparse, threshold = threshold)
  # category specific effects
  csef <- colnames(get_model_matrix(effects$cse, data))
  text_csef <- stan_csef(csef = csef, prior = prior)
  # monotonous effects
  monef <- colnames(data_monef(effects, data)$Xm)
  text_monef <- stan_monef(monef, prior = prior)
  # group-specific effects
  ranef <- gather_ranef(effects, data = data, forked = is.forked(family))
  text_ranef <- collapse_lists(
    lapply(seq_along(ranef), stan_ranef, ranef = ranef, 
           names_cov_ranef = names(cov_ranef), prior = prior))
  out <- collapse_lists(list(text_fixef, text_csef, text_monef, text_ranef))
  
  # generate Stan code for the linear predictor 'eta'
  is_multi <- is.linear(family) && length(effects$response) > 1L
  has_offset <- !is.null(model.offset(data))
  if (has_offset) {
    out$data <- paste0(out$data, "  vector[N] offset; \n")
  }
  out$transD <- paste0(out$transD,
    "  vector[N] eta;  // linear predictor \n", 
    if (length(csef)) 
      paste0("  matrix[N, ncat - 1] etap;",
             "  // linear predictor for category specific effects \n"),
    if (is_multi) 
      paste0("  vector[K_trait] Eta[N_trait];",
             "  // multivariate linear predictor matrix \n"))
  eta_obj <- ifelse(is_multi, "Eta[m, k]", "eta[n]")
  wsp <- ifelse(is_multi, "  ", "")
  # transform eta before it is passed to the likelihood
  add <- is.formula(effects[c("weights", "cens", "trunc")])
  simplify_log <- !(is_multi || is(autocor, "cor_fixed"))
  out$transform <- stan_eta_transform(family$family, family$link, add = add,
                                      simplify_log = simplify_log)
  eta_ilink <- rep("", 2)
  if (out$transform || (get_ar(autocor) && !use_cov(autocor))) {
    eta_ilink <- stan_eta_ilink(family$family, family$link, 
                                disp = is.formula(effects$disp))
    if (get_ar(autocor)) {
      eta_ar <- ifelse(!use_cov(autocor), " + head(E[n], Kar) * ar", "")
      out$transC3 <- paste0(out$transC3,  "    ", wsp, 
        eta_obj," <- ", eta_ilink[1], eta_obj, eta_ar, eta_ilink[2], "; \n")
      eta_ilink <- rep("", 2)  # don't apply the link function twice
    }
  }
  # define fixed, random, and autocorrelation effects
  nint <- length(intercepts)
  eta_int <- if (nint > 1L) " + temp_Intercept[Jint[n]]"
  eta_ranef <- stan_eta_ranef(ranef)
  eta_monef <- stan_eta_monef(monef)
  eta_ma <- ifelse(get_ma(autocor) && !use_cov(autocor), 
                   " + head(E[n], Kma) * ma", "")
  add2eta <- any(nchar(c(eta_int, eta_monef, eta_ma, eta_ranef, eta_ilink[1])))
  if (add2eta || is_multi) {
    out$transC2 <- paste0(out$transC2,
      "    ", wsp, eta_obj," <- ", eta_ilink[1], "eta[n]", 
      eta_int, eta_monef, eta_ranef, eta_ma, eta_ilink[2],"; \n")
  }
  eta_fixef <- stan_eta_fixef(fixef, sparse = sparse)
  eta_cse <- if (length(csef)) "  etap <- Xp * bp; \n"
  out$transC1 <- paste0(out$transC1,
    "  // compute linear predictor \n",
    "  eta <- ", eta_fixef, 
    if (nint == 1L && !is.ordinal(family)) " + temp_Intercept",
    if (has_offset) " + offset",
    if (get_arr(autocor)) " + Yarr * arr", 
    "; \n", eta_cse)
  out
}

stan_nonlinear <- function(effects, data, family = gaussian(), 
                           prior = prior_frame(), autocor = cor_arma(),
                           cov_ranef = NULL) {
  # prepare Stan code for non-linear models
  # Args:
  #   effects: a list returned by extract_effects()
  #   data: data.frame supplied by the user
  #   family: the model family
  #   add: Is the model weighted, censored, or truncated?
  #   cov_ranef: a list of user-defined covariance matrices
  #   prior: a prior_frame object
  #   disp: is the 'disp' addition argument specified?
  out <- list()
  if (length(effects$nonlinear)) {
    for (i in seq_along(effects$nonlinear)) {
      nlpar <- names(effects$nonlinear)[i]
      eta <- paste0("eta_", nlpar)
      out$transD <- paste0(out$transD, "  vector[N] ", eta, "; \n")
      out$data <- paste0(out$data, 
        "  // data for non-linear effects of ", nlpar, "\n")
      out$par <- paste0(out$par, 
        "  // non-linear effects of ", nlpar, "\n")
      eta_loop <- ""
      # include fixed effects
      fixef <- colnames(get_model_matrix(effects$nonlinear[[i]]$fixed, data))
      if (length(fixef)) {
        text_fixef <- stan_fixef(fixef, intercepts = NULL, family = family,
                                 prior = prior, nlpar = nlpar)
        out <- collapse_lists(list(out, text_fixef))
        out$transC1 <- paste0(out$transC1, "  ", eta, " <- ",
                              stan_eta_fixef(fixef, nlpar = nlpar), "; \n") 
      }
      # include monotonous effects
      monef <- colnames(data_monef(effects$nonlinear[[i]], data)$Xm)
      if (length(monef)) {
        text_monef <- stan_monef(monef, prior = prior, nlpar = nlpar)
        if (length(out$fun) && grepl("#include fun_monotonous", out$fun)) {
          text_monef$fun <- NULL
        }
        out <- collapse_lists(list(out, text_monef))
        eta_loop <- paste0(eta_loop, stan_eta_monef(monef, nlpar = nlpar)) 
      }
      # include random effects
      ranef <- gather_ranef(effects$nonlinear[[i]], data = data, 
                            forked = is.forked(family))
      if (length(ranef)) {
        text_ranef <- lapply(seq_along(ranef), stan_ranef, ranef = ranef, 
                             names_cov_ranef = names(cov_ranef), 
                             prior = prior, nlpar = nlpar)
        text_ranef <- collapse_lists(text_ranef)
        out <- collapse_lists(list(out, text_ranef))
        eta_loop <- paste0(eta_loop, stan_eta_ranef(ranef, nlpar = nlpar)) 
      }
      if (nchar(eta_loop)) {
        out$transC2 <- paste0(out$transC2,
          "    ", eta, "[n] <- ", eta, "[n]", eta_loop, "; \n")
      }
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
    add <- is.formula(effects[c("weights", "cens", "trunc")])
    transform <- stan_eta_transform(family$family, family$link, add = add,
                                    simplify_log = !is(autocor, "cor_fixed"))
    if (transform) {
      eta_ilink <- stan_eta_ilink(family$family, family$link, 
                                  disp = is.formula(effects$disp))
    } else eta_ilink <- rep("", 2)
    out$transD <- paste0(out$transD, "  vector[N] eta; \n")
    out$transC2 <- paste0(out$transC2, 
      "    // compute non-linear predictor \n",
      "    eta[n] <- ", eta_ilink[1], trimws(nlmodel), eta_ilink[2], "; \n")
  }
  out
}

stan_fixef <- function(fixef, intercepts = "Intercept", 
                       family = gaussian(), prior = prior_frame(), 
                       nlpar = "", sparse = FALSE, 
                       threshold = "flexible") {
  # Stan code for fixed effects
  #
  # Args:
  #   fixef: names of the fixed effects
  #   csef: names of the category specific effects
  #   family: the model family
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior 
  #   nint: number of fixed effects intercepts
  #   threshold: either "flexible" or "equidistant" 
  #
  # Returns:
  #   a list containing Stan code related to fixed effects
  nint <- length(intercepts)
  p <- if (nchar(nlpar)) paste0("_", nlpar) else ""
  out <- list()
  if (length(fixef)) {
    centered <- ifelse(nint > 0, "centered", "")
    out$data <- paste0(out$data, 
      "  int<lower=1> K", p, ";",
      "  // number of population-level effects \n", 
      "  matrix[N, K", p, "] X", p, ";",
      "  // ", centered, " population-level design matrix \n")
    if (sparse) {
      if (nchar(nlpar)) {
        stop("Sparse matrices are not yet implemented for non-linear models.",
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
  if (nint) {
    if (nchar(nlpar)) {
      stop("no special treatment of intercepts in non-linear models")
    }  
    if (length(fixef)) {
      # X_means is not defined in intercept only models
      def_X_means <- paste0("  vector[K] X_means", 
                            if (nint > 1L) "[nint]", 
                            ";  // column means of X before centering \n")
      sub_X_means <- paste0(" - dot_product(X_means", 
                            if (nint > 1L) "[i]", ", b)")
    } else {
      def_X_means <- sub_X_means <- ""
    }
    if (is.ordinal(family)) {
      # temp intercepts for ordinal models are defined in stan_ordinal
      out$data <- paste0(out$data, def_X_means)
      out$genD <- "  vector[ncat - 1] b_Intercept;  // thresholds \n" 
      out$genC <- paste0("  b_Intercept <- temp_Intercept", sub_X_means, "; \n")
    } else {
      if (nint == 1L) {
        out$data <- paste0(out$data, def_X_means)
        out$par <- paste0(out$par, 
          "  real temp_Intercept;  // temporary Intercept \n")
        out$genD <- "  real b_Intercept;  // population-level intercept \n"
        out$genC <- paste0("  b_Intercept <- temp_Intercept", sub_X_means, "; \n")
      } else if (nint > 1L) {
        out$data <- paste0(out$data,
          "  int nint;  // number of population-level intercepts \n",
          "  int Jint[N];  // assigns intercepts to observations \n", 
          def_X_means)
        out$par <- paste0(out$par, 
          "  vector[nint] temp_Intercept;  // temporary Intercepts \n")
        out$genD <- "  vector[nint] b_Intercept;  // population-level intercepts \n"
        out$genC <- paste0(
          "  for (i in 1:nint) { \n",
          "    b_Intercept[i] <- temp_Intercept[i]", sub_X_means, "; \n",
          "  } \n")
      }
    }
    # for equidistant thresholds only temp_Intercept1 is a parameter
    suffix <- ifelse(threshold == "equidistant", "1", "")
    int_prior <- stan_prior("temp_Intercept", coef = intercepts, 
                            prior = prior, suffix = suffix)
    out$prior <- paste0(out$prior, int_prior)
  }
  out
}

stan_ranef <- function(i, ranef, prior = prior_frame(), 
                       names_cov_ranef = NULL, nlpar = "") {
  # Random effects in Stan 
  # 
  # Args:
  #   i: the index of the grouping factor
  #   ranef: a named list returned by gather_ranef
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior
  #   names_cov_ranef: names of the grouping factors 
  #                    for which custom covariance matrices are specified.
  #   par: an optional character string to add to the variable names
  #        (used for non-linear models)
  #
  # Returns:
  #   A vector of strings containing the random effects in stan language
  r <- ranef[[i]]
  g <- attr(ranef[[i]], "group")
  cor <- attr(ranef[[i]], "cor")
  ccov <- g %in% names_cov_ranef
  p <- if (nchar(nlpar)) paste0(nlpar, "_", i) else i
  out <- list()
  out$data <- paste0(
    "  // data for group-specific effects of ", g, " \n",
    "  int<lower=1> J_", p, "[N]; \n",
    "  int<lower=1> N_", p, "; \n",
    "  int<lower=1> K_", p, "; \n",
    if (ccov) paste0(
      "  // cholesky factor of known covariance matrix \n",
      "  matrix[N_", p, ", N_", p,"] Lcov_", p,"; \n"))
  out$prior <- stan_prior(class = "sd", group = g, coef = r, nlpar = nlpar,
                          suffix = paste0("_", p), prior = prior)
  if (length(r) == 1L) {  # only one random effect
    out$data <- paste0(out$data, 
                       "  vector[N] Z_", p, "; \n")
    out$par <- paste0(
      "  real<lower=0> sd_", p, ";",
      "  // group-specific standard deviation \n",
      "  vector[N_", p, "] z_", p, ";",
      "  // unscaled group-specific effects \n")
    out$prior <- paste0(out$prior,"  z_", p, " ~ normal(0, 1); \n")
    out$transD <- paste0(
      "  // group-specific effects \n",
      "  vector[N_", p, "] r_", p, "; \n")
    out$transC1 <- paste0("  r_", p,  " <- sd_", p, " * (", 
                          if (ccov) paste0("Lcov_", p, " * "), "z_", p, ");\n")
  } else if (length(r) > 1L) {
    j <- seq_along(r)
    out$data <- paste0(out$data, 
                       collapse("  vector[N] Z_", p, "_", j, ";  \n"))
    out$par <- paste0(
      "  vector<lower=0>[K_", p, "] sd_", p, ";",
      "  // group-specific standard deviations \n")
    if (cor) {  # multiple correlated random effects
      out$data <- paste0(out$data, "  int<lower=1> NC_", p, "; \n")
      out$par <- paste0(out$par,
        "  matrix[K_", p, ", N_", p, "] z_", p, ";",
        "  // unscaled group-specific effects \n",    
        "  // cholesky factor of correlation matrix \n",
        "  cholesky_factor_corr[K_", p, "] L_", p, "; \n")
      out$prior <- paste0(out$prior, 
        stan_prior(class = "L", group = g, nlpar = nlpar,
                   suffix = paste0("_", p), prior = prior),
        "  to_vector(z_", p, ") ~ normal(0, 1); \n")
      out$transD <- paste0(
        "  // group-specific effects \n",
        "  matrix[N_", p, ", K_", p, "] r_", p, "; \n",
        collapse("  vector[N_", p, "] r_", p, "_", j, "; \n"))
      if (ccov) {  # customized covariance matrix supplied
        out$transC1 <- paste0(
          "  r_", p," <- as_matrix(kronecker(Lcov_", p, ",", 
          " diag_pre_multiply(sd_", p,", L_", p,")) *",
          " to_vector(z_", p, "), N_", p, ", K_", p, "); \n")
      } else { 
        out$transC1 <- paste0("  r_", p, " <- ", 
                              "(diag_pre_multiply(sd_", p, ", L_", p,") * z_", p, ")'; \n")
      }
      out$transC1 <- paste0(out$transC1, 
                            collapse("  r_", p, "_", j, " <- r_", p, "[, ",j,"];  \n"))
      # return correlations above the diagonal only
      cors_genC <- ulapply(2:length(r), function(k) 
        lapply(1:(k - 1), function(j) paste0(
          "  cor_", p, "[", (k - 1) * (k - 2) / 2 + j, 
          "] <- Cor_", p, "[", j, ",", k, "]; \n")))
      out$genD <- paste0(
        "  corr_matrix[K_", p, "] Cor_", p, "; \n",
        "  vector<lower=-1,upper=1>[NC_", p, "] cor_", p, "; \n")
      out$genC <- paste0(
        "  // take only relevant parts of correlation matrix \n",
        "  Cor_", p, " <- multiply_lower_tri_self_transpose(L_", p, "); \n",
        collapse(cors_genC)) 
    } else {  # multiple uncorrelated random effects
      out$par <- paste0(out$par,
                        "  vector[N_", p, "] z_", p, "[K_", p, "];",
                        "  // unscaled group-specific effects \n")
      out$prior <- paste0(out$prior, collapse(
        "  z_", p, "[", j, "] ~ normal(0, 1); \n"))
      out$transD <- paste0("  // group-specific effects \n", 
                           collapse("  vector[N_", p, "] r_", p, "_", j, "; \n"))
      out$transC1 <- collapse(
        "  r_", p, "_", j, " <- sd_", p, "[", j, "] * (", 
        if (ccov) paste0("Lcov_", p, " * "), "z_", p, "[", j, "]); \n")
      out$genD <- paste0(
        "  matrix[N_", p, ", K_", p, "] r_", p, "; \n")
      out$genC <- collapse(
        "  r_", p, "[, ", j, "] <- r_", p, "_", j, "; \n")
    }
  }
  out
}

stan_monef <- function(monef, prior = prior_frame(), nlpar = "") {
  # Stan code for monotonous effects
  # Args:
  #   csef: names of the monotonous effects
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior
  p <- if (nchar(nlpar)) paste0("_", nlpar) else ""
  out <- list()
  if (length(monef)) {
    I <- seq_along(monef)
    out$fun <- "  #include fun_monotonous.stan \n"
    out$data <- paste0(
      "  int<lower=1> Km", p, ";  // number of monotonous effects \n",
      "  int Xm", p, "[N, Km", p, "];  // monotonous design matrix \n",
      "  int<lower=2> Jm", p, "[Km", p, "];  // length of simplexes \n",
      collapse("  vector[Jm", p, "[", I, "]]", 
               " con_simplex", p, "_", I, "; \n"))
    bound <- get_bound(prior, class = "b", nlpar = nlpar)
    out$par <- paste0(
      "  // monotonous effects \n", 
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

stan_csef <- function(csef, prior = prior_frame()) {
  # Stan code for category specific effects
  # Args:
  #   csef: names of the category specific effects
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior
  # (!) Not implemented for non-linear models
  out <- list()
  if (length(csef)) {
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
  # Write the random effects part of the linear predictor
  # Args:
  #   ranef: a named list returned by gather_ranef
  #   nlpar: an optional character string to add to the variable names
  #         (used for non-linear models)
  # Returns:
  #   A string containing the random effects part of the linear predictor
  eta_ranef <- ""
  for (i in seq_along(ranef)) {
    p <- if (nchar(nlpar)) paste0(nlpar, "_", i) else i
    if (length(ranef[[i]]) == 1L) {
      eta_ranef <- paste0(eta_ranef, 
        " + r_", p, "[J_", p, "[n]] * Z_", p, "[n]")
    } else {
      j <- seq_along(ranef[[i]])
      eta_ranef <- paste0(eta_ranef, collapse(
        " + r_", p, "_", j, "[J_", p, "[n]]",
        " * Z_", p, "_", j, "[n]"))
    }
  }
  eta_ranef
}

stan_eta_monef <- function(monef, nlpar = "") {
  # write the linear predictor for monotonous effects
  # Args:
  #   monef: names of the monotonous effects
  #   nlpar: an optional character string to add to the variable names
  #         (used for non-linear models)
  p <- if (nchar(nlpar)) paste0("_", nlpar) else ""
  eta_monef <- ""
  for (i in seq_along(monef)) {
    eta_monef <- paste0(eta_monef,
      " + bm", p, "[", i, "] * monotonous(",
      "simplex", p, "_", i, ", Xm", p, "[n, ", i, "])")
  }
  eta_monef
}

stan_eta_transform <- function(family, link, add = FALSE, 
                               simplify_log = TRUE) {
  # indicate whether eta needs to be transformed
  # in the transformed parameters block
  # Args:
  #   add: is the model weighted, censored, truncated?
  #   simplify_log: convert gaussian(log) to lognormal? 
  !(add || !is.skewed(family) && link == "identity" 
    || family %in% "gaussian" && link == "log" && simplify_log
    || is.count(family) && link == "log" 
    || is.binary(family) && link == "logit"
    || is.ordinal(family) || is.categorical(family) 
    || is.zero_inflated(family) || is.hurdle(family))
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
