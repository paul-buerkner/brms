stan_response <- function(bterms, data) {
  # Stan code for the response variables
  stopifnot(is.brmsterms(bterms))
  family <- bterms$family
  rtype <- ifelse(use_int(family), "int", "real")
  px <- check_prefix(bterms)
  resp <- usc(combine_prefix(px))
  out <- list()
  if (rtype == "real") {
    # don't use real Y[n]
    str_add(out$data) <- paste0(
      "  vector[N] Y", resp, ";  // response variable \n"
    )
  } else if (rtype == "int") {
    str_add(out$data) <- paste0(
      "  int Y", resp, "[N];  // response variable \n"
    )
  }
  if (has_ndt(family)) {
    str_add(out$tdataD) <- paste0(
      "  real min_Y", resp, " = min(Y", resp, "); \n"
    )
  }
  if (has_trials(family)) {
    str_add(out$data) <- paste0(
      "  int trials", resp, "[N];  // number of trials \n"
    )
  }
  if (is_ordinal(family) || is_categorical(family)) {
    str_add(out$data) <- paste0(
      "  int<lower=2> ncat", resp, ";  // number of categories \n"
    )
  }
  if (is.formula(bterms$adforms$weights)) {
    str_add(out$data) <- paste0(
      "  vector<lower=0>[N] weights", resp, ";  // model weights \n" 
    )
  }
  if (is.formula(bterms$adforms$se)) {
    str_add(out$data) <- paste0(
      "  vector<lower=0>[N] se", resp, ";  // known sampling error \n"
    )
    str_add(out$tdataD) <- paste0(
      "  vector<lower=0>[N] se2", resp, " = square(se", resp, "); \n"
    )
  }
  if (is.formula(bterms$adforms$dec)) {
    str_add(out$data) <- paste0(
      "  int<lower=0,upper=1> dec", resp, "[N];  // decisions \n"
    )
  }
  has_cens <- has_cens(bterms, data = data)
  if (has_cens) {
    str_add(out$data) <- paste0(
      "  int<lower=-1,upper=2> cens", resp, "[N];  // indicates censoring \n",
      if (isTRUE(attr(has_cens, "interval"))) {
        rcens <- ifelse(rtype == "int", 
          paste0("  int rcens", resp, "[N];"), 
          paste0("  vector[N] rcens", resp, ";")
        )
        paste0(rcens, "  // right censor points for interval censoring \n")
      }
    )
  }
  bounds <- get_bounds(bterms, data = data)
  if (any(bounds$lb > -Inf)) {
    str_add(out$data) <- paste0(
      "  ", rtype, " lb", resp, "[N];  // lower truncation bounds; \n"
    )
  }
  if (any(bounds$ub < Inf)) {
    str_add(out$data) <- paste0(
      "  ", rtype, " ub", resp, "[N];  // upper truncation bounds \n"
    )
  }
  if (is.formula(bterms$adforms$mi)) {
    Ybounds <- get_bounds(bterms, data, incl_family = TRUE, stan = TRUE)
    sdy <- get_sdy(bterms, data)
    if (is.null(sdy)) {
      # response is modeled without measurement error
      str_add(out$par) <- paste0(
        "  vector", Ybounds, "[Nmi", resp, "] Ymi", resp, ";",
        "  // estimated missings\n"
      )
      str_add(out$data) <- paste0(
        "  int<lower=0> Nmi", resp, ";  // number of missings \n",
        "  int<lower=1> Jmi", resp, "[Nmi", resp, "];",  
        "  // positions of missings \n"
      )
      str_add(out$modelD) <- paste0(
        "  vector[N] Yl", resp, " = Y", resp, ";\n" 
      )
      str_add(out$modelC1) <- paste0(
        "  Yl", resp, "[Jmi", resp, "] = Ymi", resp, ";\n"
      )
    } else {
      str_add(out$data) <- paste0(
        "  // data for measurement-error in the response\n",
        "  vector<lower=0>[N] noise", resp, ";\n",
        "  // information about non-missings\n",
        "  int<lower=0> Nme", resp, ";\n",
        "  int<lower=1> Jme", resp, "[Nme", resp, "];\n"
      )
      str_add(out$par) <- paste0(
        "  vector", Ybounds, "[N] Yl", resp, ";  // latent variable\n"
      )
      str_add(out$prior) <- paste0(
        "  target += normal_lpdf(Y", resp, "[Jme", resp, "]",
        " | Yl", resp, "[Jme", resp, "],", 
        " noise", resp, "[Jme", resp, "]);\n"
      )
    }
  }
  out
}

stan_autocor <- function(bterms, prior) {
  # Stan code related to autocorrelation structures
  # Returns: 
  #   list containing Stan code snippets
  stopifnot(is.brmsterms(bterms))
  family <- bterms$family
  autocor <- bterms$autocor
  resp <- bterms$response
  is_linear <- is_linear(family)
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  out <- list()
  Kar <- get_ar(autocor)
  Kma <- get_ma(autocor)
  if (Kar || Kma) {
    err_msg <- "ARMA models are not implemented"
    if (is.mixfamily(family)) {
      stop2(err_msg, " for mixture models.") 
    }
    if (!is_linear) {
      stop2(err_msg, " for family '", family$family, "'.") 
    }
    str_add(out$data) <- paste0( 
      "  // data needed for ARMA correlations \n",
      "  int<lower=0> Kar", p, ";  // AR order \n",
      "  int<lower=0> Kma", p, ";  // MA order \n"
    )
    str_add(out$tdataD) <- paste0( 
      "  int max_lag", p, " = max(Kar", p, ", Kma", p, "); \n"
    )
    # restrict ARMA effects to be in [-1,1] when using covariance
    # formulation as they cannot be outside this interval anyway
    if (Kar) {
      ar_bound <- subset2(prior, class = "ar", ls = px)$bound
      str_add(out$par) <- paste0( 
        "  vector", ar_bound, "[Kar", p, "] ar", p, 
        ";  // autoregressive effects \n"
      )
      str_add(out$prior) <- stan_prior(
        prior, class = "ar", px = px, suffix = p
      )
    }
    if (Kma) {
      ma_bound <- subset2(prior, class = "ma", ls = px)$bound
      str_add(out$par) <- paste0( 
        "  vector", ma_bound, "[Kma", p, "] ma", p, 
        ";  // moving-average effects \n"
      )
      str_add(out$prior) <- stan_prior(
        prior, class = "ma", px = px, suffix = p
      )
    }
    if (use_cov(autocor)) {
      # if the user wants ARMA effects to be estimated using
      # a covariance matrix for residuals
      err_msg <- "ARMA covariance matrices are not implemented"
      if (isTRUE(bterms$rescor)) {
        stop2(err_msg, " when 'rescor' is estimated.")
      }
      if (any(c("sigma", "nu") %in% names(bterms$dpars))) {
        stop2(err_msg, " when predicting 'sigma' or 'nu'.")
      }
      str_add(out$data) <- paste0( 
        "  // see the functions block for details\n",
        "  int<lower=1> N_tg", p, ";\n",
        "  int<lower=1> begin_tg", p, "[N_tg", p, "];\n", 
        "  int<lower=1> end_tg", p, "[N_tg", p, "];\n", 
        "  int<lower=1> nobs_tg", p, "[N_tg", p, "];\n"
      )
      if (!is.formula(bterms$adforms$se)) {
        str_add(out$tdataD) <- paste0(
          "  vector[N] se2", p, " = rep_vector(0, N); \n"
        )
      }
      str_add(out$tparD) <- paste0(
        "  matrix[max(nobs_tg", p, "), max(nobs_tg", p, ")] res_cov_matrix;\n"               
      )
      if (Kar && !Kma) {
        cov_mat_fun <- "ar1"
        cov_mat_args <- paste0("ar", p, "[1]")
      } else if (!Kar && Kma) {
        cov_mat_fun <- "ma1"
        cov_mat_args <- paste0("ma", p, "[1]")
      } else {
        cov_mat_fun <- "arma1"
        cov_mat_args <- paste0("ar", p, "[1], ma", p, "[1]")
      }
      str_add(out$tparC1) <- paste0(
        "  // compute residual covariance matrix \n",
        "  res_cov_matrix = cov_matrix_", cov_mat_fun, 
        "(", cov_mat_args, ", sigma", p, ", max(nobs_tg", p, ")); \n"
      )
    } else {
      err_msg <- "Please set cov = TRUE in cor_arma / cor_ar / cor_ma"
      if (is.formula(bterms$adforms$se)) {
        stop2(err_msg, " when specifying 'se'.")
      }
      if (length(bterms$dpars[["mu"]]$nlpars)) {
        stop2(err_msg, " in non-linear models.")
      }
      if (!identical(family$link, "identity")) {
        stop2(err_msg, " when using non-identity links.")
      }
      str_add(out$data) <- paste0( 
        "  int<lower=0> J_lag", p, "[N]; \n"                
      )
      str_add(out$modelD) <- paste0(
        "  // objects storing residuals \n",
        "  matrix[N, max_lag", p, "] E", p,
        " = rep_matrix(0, N, max_lag", p, "); \n",
        "  vector[N] e", p, "; \n"
      )
      str_add(out$modelC2) <- paste0(
        "    // computation of ARMA correlations \n",
        "    e", p, "[n] = Y", p, "[n] - mu", p, "[n]; \n",
        "    for (i in 1:J_lag", p, "[n]) { \n",
        "      E", p, "[n + 1, i] = e", p, "[n + 1 - i]; \n",
        "    } \n"
      )
    } 
  }
  Karr <- get_arr(autocor)
  if (Karr) {
    # autoregressive effects of the response
    err_msg <- "ARR models are not implemented"
    if (length(bterms$dpars[["mu"]]$nlpars)) {
      stop2(err_msg, " for non-linear models.")
    }
    if (is.mixfamily(family)) {
      stop2(err_msg, " for mixture models.") 
    }
    str_add(out$data) <- paste0(
      "  // data needed for ARR correlations \n",
      "  int<lower=1> Karr", p, "; \n",
      "  matrix[N, Karr", p, "] Yarr", p, ";  // ARR design matrix \n"
    )
    str_add(out$par) <- paste0(
      "  vector", subset2(prior, class = "arr", ls = px)$bound, 
      "[Karr", p, "] arr", p, ";",
      "  // autoregressive effects of the response \n"
    )
    str_add(out$prior) <- stan_prior(
      prior, class = "arr", px = px, suffix = p
    )
  }
  if (is.cor_sar(autocor)) {
    err_msg <- "SAR models are not implemented"
    if (is.mixfamily(family)) {
      stop2(err_msg, " for mixture models.") 
    }
    if (!is_linear) {
      stop2(err_msg, " for family '", family$family, "'.")
    }
    if (isTRUE(bterms$rescor)) {
      stop2(err_msg, " when 'rescor' is estimated.")
    }
    if (any(c("sigma", "nu") %in% names(bterms$dpars))) {
      stop2(err_msg, " when predicting 'sigma' or 'nu'.")
    }
    str_add(out$data) <- paste0(
      "  matrix[N, N] W", p, ";  // spatial weight matrix \n"                  
    )
    if (identical(autocor$type, "lag")) {
      str_add(out$par) <- paste0( 
        "  real<lower=0,upper=1> lagsar", p, ";  // SAR parameter \n"
      )
      str_add(out$prior) <- stan_prior(
        prior, class = "lagsar", px = px, suffix = p
      )
    } else if (identical(autocor$type, "error")) {
      str_add(out$par) <- paste0( 
        "  real<lower=0,upper=1> errorsar", p, ";  // SAR parameter \n"
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
    if (length(bterms$dpars[["mu"]]$nlpars)) {
      stop2(err_msg, " in non-linear models.")
    }
    if (isTRUE(bterms$rescor)) {
      stop2(err_msg, " when 'rescor' is estimated.")
    }
    str_add(out$data) <- paste0(
      "  // data for the CAR structure \n",
      "  int<lower=1> Nloc", p, "; \n",
      "  vector[Nloc] Nneigh", p, "; \n",
      "  vector[Nloc] eigenW", p, "; \n",
      "  int<lower=1> Jloc", p, "[N]; \n",
      "  int<lower=0> Nedges", p, "; \n",
      "  int<lower=1> edges1", p, "[Nedges", p, "]; \n",
      "  int<lower=1> edges2", p, "[Nedges", p, "]; \n"
    )
    str_add(out$par) <- paste0(
      "  // parameters for the CAR structure \n",
      "  real<lower=0> sdcar", p, "; \n"
    )
    str_add(out$prior) <- stan_prior(
      prior, class = "sdcar", px = px, suffix = p
    )
    if (identical(autocor$type, "escar")) {
      str_add(out$par) <- paste0(
        "  real<lower=0, upper=1> car", p, "; \n",
        "  vector[Nloc", p, "] rcar", p, "; \n"
      )
      car_args <- c(
        "car", "sdcar", "Nloc", "Nedges", 
        "Nneigh", "eigenW", "edges1", "edges2"
      )
      car_args <- paste0(car_args, p, collapse = ", ")
      str_add(out$prior) <- paste0(
        stan_prior(prior, class = "car", px = px, suffix = p),
        "  target += sparse_car_lpdf(\n", 
        "    rcar", p, " | ", car_args, "\n",
        "  ); \n"
      )
    } else if (identical(autocor$type, "esicar")) {
      str_add(out$par) <- paste0(
        "  vector[Nloc - 1] zcar", p, "; \n"
      )
      str_add(out$tparD) <- paste0(
        "  vector[Nloc] rcar", p, "; \n"                
      )
      str_add(out$tparC1) <- paste0(
        "  // apply sum-to-zero constraint \n",
        "  rcar[1:(Nloc", p, " - 1)] = zcar", p, "; \n",
        "  rcar[Nloc", p, "] = - sum(zcar", p, "); \n"
      )
      car_args <- c(
        "sdcar", "Nloc", "Nedges", 
        "Nneigh", "eigenW", "edges1", "edges2"
      )
      car_args <- paste0(car_args, p, collapse = ", ")
      str_add(out$prior) <- paste0(
        "  target += sparse_icar_lpdf(\n", 
        "    rcar", p, " | ", car_args, "\n",
        "  ); \n"
      )
    } 
  }
  if (is.cor_bsts(autocor)) {
    warning2(
      "The `bsts' correlation structure has been deprecated and ",
      "will be removed from the package at some point. Consider ", 
      "using splines or Gaussian processes instead."
    )
    err_msg <- "BSTS models are not implemented"
    if (is.mixfamily(family)) {
      stop2(err_msg, " for mixture models.") 
    }
    if (is_ordinal(family) || family$family %in% c("bernoulli", "categorical")) {
      stop2(err_msg, " for family '", family$family, "'.")
    }
    if (length(bterms$dpars[["mu"]]$nlpars)) {
      stop2(err_msg, " in non-linear models.")
    }
    str_add(out$data) <- paste0(
      "  vector[N] tg", p, ";  // indicates independent groups \n"
    )
    str_add(out$par) <- paste0(
      "  vector[N] loclev", p, ";  // local level terms \n",
      "  real<lower=0> sigmaLL", p, ";  // SD of local level terms \n"
    )
    if (is_linear) {
      # ensures that the means of the loclev priors start close 
      # to the response values; this often helps with convergence
      link <- stan_link(family$link)
      center <- paste0(link, "(Y", p, "[", c("1", "n"), "])")
    } else {
      center <- c("0", "0")
    }
    str_add(out$prior) <- paste0(
      stan_prior(prior, class = "sigmaLL", px = px, suffix = p),
      "  target += normal_lpdf(loclev", p, "[1]",
      " | ", center[1], ", sigmaLL", p, "); \n",
      "  for (n in 2:N) { \n",
      "    if (tg", p, "[n] == tg", p, "[n - 1]) { \n",
      "      target += normal_lpdf(loclev", p, "[n]",
      " | loclev", p, "[n - 1], sigmaLL", p, "); \n",
      "    } else { \n",
      "      target += normal_lpdf(loclev", p, "[n]",
      " | ", center[2], ", sigmaLL", p, "); \n",
      "    } \n",
      "  } \n"
    )
  }
  if (is.cor_fixed(autocor)) {
    err_msg <- "Fixed residual covariance matrices are not implemented"
    if (is.mixfamily(family)) {
      stop2(err_msg, " for mixture models.") 
    }
    if (!is_linear) {
      stop2(err_msg, " for family '", family$family, "'.")
    }
    if (isTRUE(bterms$rescor)) {
      stop2(err_msg, " when 'rescor' is estimated.")
    }
    str_add(out$data) <- paste0( 
      "  matrix[N, N] V", p, ";  // known residual covariance matrix \n"
    )
    if (family$family %in% "gaussian") {
      str_add(out$tdataD) <- paste0(
        "  matrix[N, N] LV", p, " = cholesky_decompose(V", p, "); \n"
      )
    }
  }
  out
}

stan_global_defs <- function(bterms, prior, ranef, cov_ranef) {
  # define Stan functions or globally used transformed data
  # Returns:
  #   a list of character strings
  families <- family_names(bterms)
  links <- family_names(bterms, link = TRUE)
  unique_combs <- !duplicated(paste0(families, ":", links))
  families <- families[unique_combs]
  links <- links[unique_combs]
  out <- list()
  if (any(links == "cauchit")) {
    str_add(out$fun) <- "  #include 'fun_cauchit.stan' \n"
  } else if (any(links == "cloglog")) {
    str_add(out$fun) <- "  #include 'fun_cloglog.stan' \n"
  }
  if (any(families %in% c("student", "frechet"))) {
    str_add(out$fun) <- "  #include 'fun_logm1.stan' \n"
  }
  hs_dfs <- ulapply(attr(prior, "special"), "[[", "hs_df")
  if (any(nzchar(hs_dfs))) {
    str_add(out$fun) <- "  #include 'fun_horseshoe.stan' \n"
  }
  if (stan_needs_kronecker(ranef, names(cov_ranef))) {
    str_add(out$fun) <- paste0(
      "  #include 'fun_as_matrix.stan' \n",
      "  #include 'fun_kronecker.stan' \n"
    )
  }
  if (any(families %in% "zero_inflated_poisson")) {
    str_add(out$fun) <- "  #include 'fun_zero_inflated_poisson.stan' \n"
  } 
  if (any(families %in% "zero_inflated_negbinomial")) {
    str_add(out$fun) <- "  #include 'fun_zero_inflated_negbinomial.stan' \n"
  } 
  if (any(families %in% "zero_inflated_binomial")) {
    str_add(out$fun) <- "  #include 'fun_zero_inflated_binomial.stan' \n"
  }
  if (any(families %in% "zero_inflated_beta")) {
    str_add(out$fun) <- "  #include 'fun_zero_inflated_beta.stan' \n"
  }
  if (any(families %in% "zero_one_inflated_beta")) {
    str_add(out$fun) <- "  #include 'fun_zero_one_inflated_beta.stan' \n"
  }
  if (any(families %in% "hurdle_poisson")) {
    str_add(out$fun) <- "  #include 'fun_hurdle_poisson.stan' \n"
  }
  if (any(families %in% "hurdle_negbinomial")) {
    str_add(out$fun) <- "  #include 'fun_hurdle_negbinomial.stan' \n"
  }
  if (any(families %in% "hurdle_gamma")) {
    str_add(out$fun) <- "  #include 'fun_hurdle_gamma.stan' \n"
  }
  if (any(families %in% "hurdle_lognormal")) {
    str_add(out$fun) <- "  #include 'fun_hurdle_lognormal.stan' \n"
  }
  if (any(families %in% "inverse.gaussian")) {
    str_add(out$fun) <- "  #include 'fun_inv_gaussian.stan' \n"
  } 
  if (any(families %in% "von_mises")) {
    str_add(out$fun) <- paste0(
      "  #include 'fun_tan_half.stan' \n",
      "  #include 'fun_von_mises.stan' \n"
    )
  } 
  if (any(families %in% "wiener")) {
    str_add(out$fun) <- "  #include 'fun_wiener_diffusion.stan' \n"
  }
  if (any(families %in% "asym_laplace")) {
    str_add(out$fun) <- "  #include 'fun_asym_laplace.stan' \n"
  } 
  if (any(families %in% "skew_normal")) {
    str_add(out$tdataD) <- "  real sqrt_2_div_pi = sqrt(2 / pi()); \n"
  } 
  if (any(families %in% "gen_extreme_value")) {
    str_add(out$fun) <- paste0(
      "  #include 'fun_gen_extreme_value.stan' \n",
      "  #include 'fun_scale_xi.stan' \n"
    )
  }
  is_ordinal <- ulapply(families, is_ordinal)
  if (any(is_ordinal)) {
    ord_families <- families[is_ordinal]
    ord_links <- links[is_ordinal]
    for (i in seq_along(ord_families)) {
      for (cs in c(FALSE, TRUE)) {
        str_add(out$fun) <- stan_ordinal_lpmf(ord_families[i], ord_links[i], cs)
      }
    }
  }
  uni_mo <- ulapply(get_effect(bterms, "sp"), attr, "uni_mo")
  if (length(uni_mo)) {
    str_add(out$fun) <- "  #include fun_monotonic.stan \n"
  } 
  if (length(get_effect(bterms, "gp"))) {
    str_add(out$fun) <- "  #include fun_gaussian_process.stan \n"
  }
  # functions related to autocorrelation structures
  if (is.brmsterms(bterms)) {
    autocors <- list(bterms$autocor)
  } else {
    autocors <- lapply(bterms$terms, "[[", "autocor")
  }
  if (any(ulapply(autocors, use_cov))) {
    str_add(out$fun) <- paste0(
      "  #include 'fun_normal_cov.stan'\n",
      "  #include 'fun_student_t_cov.stan'\n",
      "  #include 'fun_cov_matrix_ar1.stan'\n",
      "  #include 'fun_cov_matrix_ma1.stan'\n",
      "  #include 'fun_cov_matrix_arma1.stan'\n"
    )
  }
  if (any(ulapply(autocors, is.cor_sar))) {
    str_add(out$fun) <- paste0(
      "  #include 'fun_normal_lagsar.stan'\n",
      "  #include 'fun_student_t_lagsar.stan'\n",
      "  #include 'fun_normal_errorsar.stan'\n",
      "  #include 'fun_student_t_errorsar.stan'\n"
    )
  }
  if (any(ulapply(autocors, is.cor_car))) {
    str_add(out$fun) <- paste0(
      "  #include 'fun_sparse_car_lpdf.stan'\n",      
      "  #include 'fun_sparse_icar_lpdf.stan'\n"
    )
  }
  out
}

stan_ordinal_lpmf <- function(family, link, cs = FALSE) {
  # ordinal log-probability densitiy functions in Stan language
  # Args:
  #   cs: Logical; add category specific effects?
  stopifnot(is.character(family), is.character(link))
  cs <- as_one_logical(cs)
  ilink <- stan_ilink(link)
  th <- function(k) {
    # helper function generating stan code inside ilink(.)
    sign <- ifelse(family %in% c("cumulative", "sratio"), " - ", " + ")
    ptl <- ifelse(cs, paste0(sign, "mucs[", k, "]"), "") 
    if (sign == " - ") {
      out <- paste0("thres[", k, "]", ptl, " - mu")
    } else {
      out <- paste0("mu", ptl, " - thres[", k, "]")
    }
    paste0("disc * (", out, ")")
  }
  cs_arg <- ifelse(cs, "row_vector mucs, ", "")
  out <- paste0(
    "  /* ", family, "-", link, " log-PDF for a single response \n",
    if (cs) "   * including category specific effects \n",
    "   * Args: \n",
    "   *   y: response category \n",
    "   *   mu: linear predictor \n",
    if (cs) "   *   mucs: predictor for category specific effects \n",
    "   *   thres: ordinal thresholds \n",
    "   *   disc: discrimination parameter \n",
    "   * Returns: \n", 
    "   *   a scalar to be added to the log posterior \n",
    "   */ \n",
    "   real ", family, "_", link, ifelse(cs, "_cs", ""), 
    "_lpmf(int y, real mu, ", cs_arg, "vector thres, real disc) { \n",
    "     int ncat = num_elements(thres) + 1; \n",
    "     vector[ncat] p; \n",
    if (family != "cumulative") "     vector[ncat - 1] q; \n"
  )
  # define actual function content
  if (family == "cumulative") {
    str_add(out) <- paste0(
      "     p[1] = ", ilink, "(", th(1), "); \n",
      "     for (k in 2:(ncat - 1)) { \n", 
      "       p[k] = ", ilink, "(", th("k"), ") - \n",
      "              ", ilink, "(", th("k - 1"), "); \n", 
      "     } \n",
      "     p[ncat] = 1 - ", ilink, "(", th("ncat - 1"), "); \n"
    )
  } else if (family %in% c("sratio", "cratio")) {
    sc <- ifelse(family == "sratio", "1 - ", "")
    str_add(out) <- paste0(
      "     for (k in 1:(ncat - 1)) { \n",
      "       q[k] = ", sc, ilink, "(", th("k"), "); \n",
      "       p[k] = 1 - q[k]; \n",
      "       for (kk in 1:(k - 1)) p[k] = p[k] * q[kk]; \n", 
      "     } \n",
      "     p[ncat] = prod(q); \n"
    )
  } else if (family == "acat") {
    if (ilink == "inv_logit") {
      str_add(out) <- paste0(
        "     p[1] = 1.0; \n",
        "     for (k in 1:(ncat - 1)) { \n",
        "       q[k] = ", th("k"), "; \n",
        "       p[k + 1] = q[1]; \n",
        "       for (kk in 2:k) p[k + 1] = p[k + 1] + q[kk]; \n",
        "       p[k + 1] = exp(p[k + 1]); \n",
        "     } \n",
        "     p = p / sum(p); \n"
      )
    } else {
      str_add(out) <- paste0(   
        "     for (k in 1:(ncat - 1)) \n",
        "       q[k] = ", ilink, "(", th("k"), "); \n",
        "     for (k in 1:ncat) { \n",     
        "       p[k] = 1.0; \n",
        "       for (kk in 1:(k - 1)) p[k] = p[k] * q[kk]; \n",
        "       for (kk in k:(ncat - 1)) p[k] = p[k] * (1 - q[kk]); \n",   
        "     } \n",
        "     p = p / sum(p); \n"
      )
    }
  }
  str_add(out) <- "     return categorical_lpmf(y | p); \n   } \n" 
  out
}

stan_mixture <- function(bterms, prior) {
  # Stan code specific for mixture families
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  out <- list()
  if (is.mixfamily(bterms$family)) {
    nmix <- length(bterms$family$mix)
    theta_pred <- grepl("^theta", names(bterms$dpars))
    theta_pred <- bterms$dpars[theta_pred]
    theta_fix <- grepl("^theta", names(bterms$fdpars))
    theta_fix <- bterms$fdpars[theta_fix]
    def_thetas <- collapse(
      "  real<lower=0,upper=1> theta", 1:nmix, p, ";",
      "  // mixing proportion \n"
    )
    if (length(theta_pred)) {
      if (length(theta_pred) != nmix - 1) {
        stop2("Can only predict all but one mixing proportion.")
      }
      missing_id <- setdiff(1:nmix, dpar_id(names(theta_pred)))
      str_add(out$modelD) <- paste0(
        "  vector[N] theta", missing_id, p, " = rep_vector(0, N); \n",                   
        "  real log_sum_exp_theta; \n"      
      )
      sum_exp_theta <- paste0(
        "exp(theta", 1:nmix, p, "[n])", collapse = " + "
      )
      str_add(out$modelC3) <- paste0( 
        "    log_sum_exp_theta = log(", sum_exp_theta, "); \n",
        collapse(
          "    theta", 1:nmix, p, "[n] = theta", 1:nmix, p, "[n]", 
          " - log_sum_exp_theta; \n"
        )
      )
    } else if (length(theta_fix)) {
      if (length(theta_fix) != nmix) {
        stop2("Can only fix no or all mixing proportions.")
      }
      str_add(out$data) <- paste0(
        "  // mixing proportions \n",                
        collapse(
          "  real<lower=0,upper=1> theta", 1:nmix, p, "; \n"
        )
      )
    } else {
      str_add(out$data) <- paste0(
        "  vector[", nmix, "] con_theta", p, 
        ";  // prior concentration \n"                  
      )
      str_add(out$par) <- paste0(
        "  simplex[", nmix, "] theta", p, 
        ";  // mixing proportions \n"
      )
      str_add(out$prior) <- paste0(
        "  target += dirichlet_lpdf(theta", p, " | con_theta", p, "); \n"                
      )
      str_add(out$tparD) <- paste0(
        "  // mixing proportions \n",                
        collapse(
          "  real<lower=0,upper=1> theta", 1:nmix, p,
          " = theta", p, "[", 1:nmix, "]; \n"
        )
      )
    }
    if (bterms$family$order %in% "mu") {
      str_add(out$par) <- paste0( 
        "  ordered[", nmix, "] ordered_Intercept", p,
        ";  // to identify mixtures \n"
      )
    }
  }
  out
}

stan_Xme <- function(bterms, prior) {
  # global Stan definitions for noise-free variables
  out <- list()
  uni_me <- rename(get_uni_me(bterms))
  if (length(uni_me)) {
    K <- paste0("_", seq_along(uni_me))
    str_add(out$data) <- paste0(
      "  // noisy variables\n",
      collapse("  vector[N] Xn", K, ";\n"),
      "  // measurement noise\n",
      collapse("  vector<lower=0>[N] noise", K, ";\n")
    )
    str_add(out$par) <- paste0(
      "  // parameters for noise free variables\n",
      collapse(
        "  vector[N] zme", K, ";\n",
        "  real meanme", K, ";\n",
        "  real<lower=0> sdme", K, ";\n"
      )
    )
    str_add(out$tparD) <- collapse(
      "  vector[N] Xme", K, " = meanme", K, " + sdme", K, " * zme", K, ";\n"
    )
    for (k in seq_along(uni_me)) {
      sfx <- paste0("_", k)
      str_add(out$prior) <- stan_prior(
        prior, class = "meanme", coef = uni_me[k], suffix = sfx
      )
      str_add(out$prior) <- stan_prior(
        prior, class = "sdme", coef = uni_me[k], suffix = sfx
      )
    }
    str_add(out$prior) <- collapse(
      "  target += normal_lpdf(Xn", K, " | Xme", K, ", noise", K, ");\n",
      "  target += normal_lpdf(zme", K, " | 0, 1);\n"
    )
  }
  out
}

stan_prior <- function(prior, class, coef = "", group = "", 
                       px = list(), prefix = "", suffix = "", 
                       wsp = 2, matrix = FALSE) {
  # Define priors for parameters in Stan language
  # Args:
  #   prior: an object of class 'brmsprior'
  #   class: the parameter class
  #   coef: the coefficients of this class
  #   group: the name of a grouping factor
  #   nlpar: the name of a non-linear parameter
  #   prefix: a prefix to put at the parameter class
  #   suffix: a suffix to put at the parameter class
  #   matrix: logical; corresponds the class to a parameter matrix?
  #   wsp: an integer >= 0 defining the number of spaces 
  #        in front of the output string
  # Returns:
  #   A character strings in stan language that defines priors 
  #   for a given class of parameters. If a parameter has has 
  #   no corresponding prior in prior, an empty string is returned.
  tp <- tp(wsp)
  wsp <- wsp(nsp = wsp)
  prior_only <- identical(attr(prior, "sample_prior"), "only")
  prior <- subset2(prior, 
    class = class, coef = c(coef, ""), group = c(group, "")
  )
  if (class %in% c("sd", "cor")) {
    # only sd and cor parameters have global priors
    px_tmp <- lapply(px, function(x) c(x, ""))
    prior <- subset2(prior, ls = px_tmp)
  } else {
    prior <- subset2(prior, ls = px)
  }
  if (!nchar(class) && nrow(prior)) {
    # unchecked prior statements are directly passed to Stan
    return(collapse(wsp, prior$prior, "; \n"))
  } 
  
  px <- as.data.frame(px)
  upx <- unique(px)
  if (nrow(upx) > 1L) {
    # can only happen for SD parameters of the same ID
    base_prior <- rep(NA, nrow(upx))
    for (i in seq_len(nrow(upx))) {
      sub_upx <- lapply(upx[i, ], function(x) c(x, ""))
      sub_prior <- subset2(prior, ls = sub_upx) 
      base_prior[i] <- stan_base_prior(sub_prior)
    }
    if (length(unique(base_prior)) > 1L) {
      # define prior for single coefficients manually
      # as there is not single base_prior anymore
      prior_of_coefs <- prior[nzchar(prior$coef), vars_prefix()]
      take <- match_rows(prior_of_coefs, upx)
      prior[nzchar(prior$coef), "prior"] <- base_prior[take]
    }
    base_prior <- base_prior[1]
    bound <- ""
  } else {
    base_prior <- stan_base_prior(prior)
    bound <- prior[!nzchar(prior$coef), "bound"]
  }
  
  individual_prior <- function(i, prior, max_index) {
    # individual priors for each parameter of a class
    if (max_index > 1L || matrix) {
      index <- paste0("[", i, "]")      
    } else {
      index <- ""
    }
    if (nrow(px) > 1L) {
      prior <- subset2(prior, ls = px[i, ])
    }
    uc_prior <- prior$prior[match(coef[i], prior$coef)]
    if (!is.na(uc_prior) & nchar(uc_prior)) { 
      # user defined prior for this parameter
      coef_prior <- uc_prior
    } else { 
      # base prior for this parameter
      coef_prior <- base_prior 
    }  
    if (nzchar(coef_prior)) {  
      # implies a proper prior
      pars <- paste0(class, index)
      out <- stan_target_prior(coef_prior, pars, bound = bound)
      out <- paste0(tp, out, "; \n")
    } else {
      # implies an improper flat prior
      out <- ""
    }
    return(out)
  }
  
  # generate stan prior statements
  class <- paste0(prefix, class, suffix)
  if (any(with(prior, nchar(coef) & nchar(prior)))) {
    # generate a prior for each coefficient
    out <- sapply(
      seq_along(coef), individual_prior, 
      prior = prior, max_index = length(coef)
    )
  } else if (nchar(base_prior) > 0) {
    if (matrix) {
      class <- paste0("to_vector(", class, ")")
    }
    out <- stan_target_prior(
      base_prior, class, ncoef = length(coef), bound = bound
    )
    out <- paste0(tp, out, "; \n")
  } else {
    out <- ""
  }
  special_prior <- stan_special_prior(
    class, prior, ncoef = length(coef), px = px
  )
  out <- collapse(c(out, special_prior))
  if (prior_only && nzchar(class) && !nchar(out)) {
    stop2("Sampling from priors is not possible as ", 
          "some parameters have no proper priors. ",
          "Error occured for class '", class, "'.")
  }
  out
}

stan_base_prior <- function(prior) {
  # get base (highest level) prior of all priors
  # Args:
  #   prior: a prior.frame
  stopifnot(length(unique(prior$class)) <= 1L) 
  prior <- prior[with(prior, !nzchar(coef) & nzchar(prior)), ]
  vars <- c("group", "nlpar", "dpar", "resp", "class")
  i <- 1
  found <- FALSE
  base_prior <- ""
  take <- rep(FALSE, nrow(prior))
  while (!found && i <= length(vars)) {
    take <- nzchar(prior[[vars[i]]]) & !take
    if (any(take)) {
      base_prior <- prior[take, "prior"]
      found <- TRUE
    }
    i <- i + 1
  }
  stopifnot(length(base_prior) == 1L)
  base_prior
}

stan_target_prior <- function(prior, par, ncoef = 1, bound = "") {
  prior <- gsub("[[:space:]]+\\(", "(", prior)
  prior_name <- get_matches(
    "^[^\\(]+(?=\\()", prior, perl = TRUE, simplify = FALSE
  )
  for (i in seq_along(prior_name)) {
    if (length(prior_name[[i]]) != 1L) {
      stop2("The prior '", prior[i], "' is invalid.")
    }
  }
  prior_name <- unlist(prior_name)
  prior_args <- rep(NA, length(prior))
  for (i in seq_along(prior)) {
    prior_args[i] <- sub(
      paste0("^", prior_name[i], "\\("), "", prior[i]
    )
  }
  out <- paste0(prior_name, "_lpdf(", par, " | ", prior_args)
  par_class <- unique(get_matches("^[^_]+", par))
  par_bound <- par_bounds(par_class, bound)
  prior_bound <- prior_bounds(prior_name)
  trunc_lb <- is.character(par_bound$lb) || par_bound$lb > prior_bound$lb
  trunc_ub <- is.character(par_bound$ub) || par_bound$ub < prior_bound$ub
  if (trunc_lb || trunc_ub) {
    wsp <- wsp(nsp = 4)
    if (trunc_lb && !trunc_ub) {
      str_add(out) <- paste0(
        "\n", wsp, "- ", ncoef, " * ", prior_name, "_lccdf(", 
        par_bound$lb, " | ", prior_args
      )
    } else if (!trunc_lb && trunc_ub) {
      str_add(out) <- paste0(
        "\n", wsp, "- ", ncoef, " * ", prior_name, "_lcdf(", 
        par_bound$ub, " | ", prior_args
      )
    } else if (trunc_lb && trunc_ub) {
      str_add(out) <- paste0(
        "\n", wsp, "- ", ncoef, " * log_diff_exp(", 
        prior_name, "_lcdf(", par_bound$ub, " | ", prior_args, ", ",
        prior_name, "_lcdf(", par_bound$lb, " | ", prior_args, ")"
      )
    }
  }
  out
}

stan_special_prior <- function(class, prior, ncoef, px = list()) {
  # add special priors such as horseshoe and lasso
  out <- ""
  p <- usc(combine_prefix(px))
  if (all(class == paste0("b", p))) {
    stopifnot(length(p) == 1L)
    tp <- tp()
    # add horseshoe and lasso shrinkage priors
    prefix <- combine_prefix(px, keep_mu = TRUE)
    special <- attributes(prior)$special[[prefix]]
    if (!is.null(special$hs_df)) {
      local_args <- paste0("0.5 * hs_df", p)
      local_args <- sargs(local_args, local_args)
      global_args <- paste0("0.5 * hs_df_global", p)
      global_args <- sargs(global_args, global_args)
      c2_args <- paste0("0.5 * hs_df_slab", p)
      c2_args <- sargs(c2_args, c2_args)
      wsp <- wsp(nsp = 4)
      str_add(out) <- paste0(
        tp, "normal_lpdf(zb", p, " | 0, 1); \n",
        tp, "normal_lpdf(hs_local", p, "[1] | 0, 1)\n", 
        wsp, "- ", ncoef, " * log(0.5); \n",
        tp, "inv_gamma_lpdf(hs_local", p, "[2] | ", local_args, "); \n",
        tp, "normal_lpdf(hs_global", p, "[1] | 0, 1)\n", 
        wsp, "- 1 * log(0.5); \n",
        tp, "inv_gamma_lpdf(hs_global", p, "[2] | ", global_args, "); \n",
        tp, "inv_gamma_lpdf(hs_c2", p, " | ", c2_args, "); \n"
      )
    }
    if (!is.null(special$lasso_df)) {
      str_add(out) <- paste0(
        tp, "chi_square_lpdf(lasso_inv_lambda", p, " | lasso_df", p, "); \n"
      )
    }
  }
  out
}

stan_rngprior <- function(sample_prior, prior, par_declars,
                          gen_quantities, prior_special) {
  # stan code to sample from priors seperately
  # Args:
  #   sample_prior: take samples from priors?
  #   prior: character string taken from stan_prior
  #   par_declars: the parameters block of the Stan code
  #     required to extract boundaries
  #   gen_quantities: Stan code from the generated quantities block
  #   prior_special: a list of values pertaining to special priors
  #     such as horseshoe or lasso
  # Returns:
  #   a character string containing the priors to be sampled from in stan code
  out <- list()
  if (identical(sample_prior, "yes")) {
    prior <- strsplit(gsub(" |\\n", "", prior), ";")[[1]]
    prior <- prior[nzchar(prior)]
    pars_regex <- "(?<=_lpdf\\()[^|]+" 
    pars <- get_matches(pars_regex, prior, perl = TRUE, first = TRUE)
    pars <- gsub("to_vector\\(|\\)$", "", pars)
    excl_regex <- "^(z|zs|zb|zgp|Xme|hs)_?|^increment_log_prob\\("
    take <- !grepl(excl_regex, pars)
    prior <- prior[take]
    pars <- sub("^L_", "cor_", pars[take])
    pars <- sub("^Lrescor", "rescor", pars)
    dis_regex <- "(?<=target\\+=)[^\\(]+(?=_lpdf\\()"
    dis <- get_matches(dis_regex, prior, perl = TRUE, first = TRUE)
    dis <- sub("corr_cholesky$", "corr", dis)
    args_regex <- "(?<=\\|)[^$\\|]+(?=\\)($|-))"
    args <- get_matches(args_regex, prior, perl = TRUE, first = TRUE)
    args <- ifelse(grepl("^lkj_corr$", dis), paste0("2,", args), args)
    type <- rep("real", length(pars))
    
    # rename parameters containing indices
    has_ind <- grepl("\\[[[:digit:]]+\\]", pars)
    pars[has_ind] <- ulapply(pars[has_ind], function(par) {
      ind_regex <- "(?<=\\[)[[:digit:]]+(?=\\])"
      ind <- get_matches(ind_regex, par, perl = TRUE)
      gsub("\\[[[:digit:]]+\\]", paste0("_", ind), par)
    })
    
    # extract information from the initial parameter definition
    par_declars <- unlist(strsplit(par_declars, "\n", fixed = TRUE))
    par_declars <- gsub("^[[:blank:]]*", "", par_declars)
    par_declars <- par_declars[!grepl("^//", par_declars)]
    all_pars_regex <- "(?<= )[^[:blank:]]+(?=;)"
    all_pars <- get_matches(all_pars_regex, par_declars, perl = TRUE) 
    all_bounds <- get_matches("<.+>", par_declars, simplify = FALSE)
    all_bounds <- ulapply(all_bounds, function(x) if (length(x)) x else "")
    all_types <- get_matches("^[^[:blank:]]+", par_declars)
    
    # define parameter types and boundaries
    bounds <- rep("", length(pars))
    types <- rep("real", length(pars))
    for (i in seq_along(all_pars)) {
      k <- which(grepl(paste0("^", all_pars[i]), pars))
      bounds[k] <- all_bounds[i]
      if (grepl("^simo", all_pars[i])) {
        types[k] <- all_types[i]
      }
    }
    
    # distinguish between bounded and unbounded parameters
    has_bounds <- as.logical(nchar(bounds))
    if (any(has_bounds)) {  
      # bounded parameters have to be sampled in the model block
      str_add(out$par) <- paste0(
        "  // parameters to store prior samples\n",
        collapse(
          "  real", bounds[has_bounds], 
          " prior_", pars[has_bounds], ";\n"
        )
      )
      str_add(out$model) <- paste0(
        "  // additionally draw samples from priors\n",
        collapse(
          "  target += ", dis[has_bounds], "_lpdf(",
          "prior_", pars[has_bounds], " | ", args[has_bounds], "); \n"
        )
      )
    }
    no_bounds <- !has_bounds
    if (any(no_bounds)) {
      # use parameters sampled from priors for use in other priors
      spars <- NULL
      # cannot sample from the horseshoe prior anymore as of brms 1.5.0
      lasso_prefix <- nzchar(ulapply(prior_special, "[[", "lasso_df"))
      lasso_prefix <- names(prior_special)[lasso_prefix]
      lasso_prefix <- usc(sub("^mu(_|$)", "", lasso_prefix))
      if (length(lasso_prefix)) {
        spars <- c(spars, paste0("lasso_inv_lambda", lasso_prefix))
      }
      if (length(spars)) {
        bpars <- grepl("^b(|mo|cs|me)(_|$)", pars)
        args[bpars] <- rename(args[bpars], spars, paste0("prior_", spars))
      }
      lkj_index <- ifelse(grepl("^lkj_corr$", dis[no_bounds]), "[1, 2]", "")
      # unbounded parameters can be sampled in the generatated quantities block
      str_add(out$genD) <- paste0(
        "  // additionally draw samples from priors \n",
        collapse(
          "  ", types[no_bounds], " prior_", pars[no_bounds], 
          " = ", dis[no_bounds], "_rng(", args[no_bounds], ")",
          lkj_index, "; \n"
        )
      )
    }
    # compute priors for the actual population-level intercepts
    is_temp_intercept <- grepl("^temp.*_Intercept", pars)
    if (any(is_temp_intercept)) {
      temp_intercepts <- pars[is_temp_intercept]
      p <- gsub("^temp|_Intercept$", "", temp_intercepts)
      intercepts <- paste0("b", p, "_Intercept")
      regex <- paste0(" (\\+|-) dot_product\\(means_X", p, ", b", p, "\\)")
      sub_X_means <- lapply(regex, get_matches, gen_quantities)
      sub_X_means <- ulapply(sub_X_means, 
        function(x) if (length(x)) x[1] else ""
      )
      str_add(out$genD) <- collapse(
        "  real prior_", intercepts,
        " = prior_", temp_intercepts, sub_X_means, "; \n"
      )
    }
  }
  out
}

stan_link <- function(link) {
  # find the link in Stan language
  # Args:
  #   link: the link function
  switch(link, 
    identity = "",
    log = "log", 
    logm1 = "logm1",
    inverse = "inv",
    sqrt = "sqrt", 
    "1/mu^2" = "inv_square", 
    logit = "logit", 
    probit = "inv_Phi", 
    probit_approx = "inv_Phi", 
    cloglog = "cloglog", 
    cauchit = "cauchit",
    tan_half = "tan_half",
    log1p = "log1p"
  )
}

stan_ilink <- function(link) {
  # find the inverse link in Stan language
  # Args:
  #   link: the link function
  switch(link, 
    identity = "",
    log = "exp", 
    logm1 = "expp1",
    inverse = "inv", 
    sqrt = "square", 
    "1/mu^2" = "inv_sqrt", 
    logit = "inv_logit", 
    probit = "Phi", 
    probit_approx = "Phi_approx", 
    cloglog = "inv_cloglog",
    cauchit = "inv_cauchit",
    tan_half = "inv_tan_half",
    log1p = "expm1"
  )
}

stan_vector <- function(...) {
  # define a vector in Stan language
  paste0("[", paste0(c(...), collapse = ", "), "]'")
}

stan_has_built_in_fun <- function(family) {
  # indicates if a family-link combination has a built in 
  # function in Stan (such as binomial_logit)
  # Args:
  #   family: a list with elements 'family' and 'link'
  stopifnot(all(c("family", "link") %in% names(family)))
  link <- family$link
  dpar <- family$dpar
  family <- family$family
  log_families <- c(
    "poisson", "negbinomial", "geometric", 
    "zero_inflated_poisson", "zero_inflated_negbinomial",
    "hurdle_poisson", "hurdle_negbinomial"
  )
  logit_families <- c(
    "binomial", "bernoulli", "cumulative", "categorical",
    "zero_inflated_binomial"
  )
  isTRUE(
    family %in% log_families && link == "log" ||
    family %in% logit_families && link == "logit" ||
    isTRUE(dpar %in% c("zi", "hu")) && link == "logit"
  )
}

stan_needs_kronecker <- function(ranef, names_cov_ranef) {
  # checks if a model needs the kronecker product
  # Args: 
  #   ranef: named list returned by tidy_ranef
  #   names_cov_ranef: names of the grouping factors that
  #                    have a cov.ranef matrix 
  ids <- unique(ranef$id)
  out <- FALSE
  for (id in ids) {
    r <- ranef[ranef$id == id, ]
    out <- out || nrow(r) > 1L && r$cor[1] && 
             r$group[1] %in% names_cov_ranef
  }
  out
}
