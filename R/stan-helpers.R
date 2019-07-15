# unless otherwise specifiedm functions return a named list 
# of Stan code snippets to be pasted together later on

# Stan code for the response variables
stan_response <- function(bterms, data) {
  stopifnot(is.brmsterms(bterms))
  family <- bterms$family
  rtype <- str_if(use_int(family), "int", "real")
  multicol <- has_multicol(family)
  px <- check_prefix(bterms)
  resp <- usc(combine_prefix(px))
  out <- list()
  str_add(out$data) <- glue(
    "  int<lower=1> N{resp};  // number of observations\n"
  )
  if (has_cat(family) || is.formula(bterms$adforms$cat)) {
    str_add(out$data) <- glue(
      "  int<lower=2> ncat{resp};  // number of categories\n"
    )
  }
  if (has_multicol(family)) {
    if (rtype == "real") {
      str_add(out$data) <- glue(
        "  vector[ncat{resp}] Y{resp}[N{resp}];  // response array\n"
      )
    } else if (rtype == "int") {
      str_add(out$data) <- glue(
        "  int Y{resp}[N{resp}, ncat{resp}];  // response array\n"
      )
    }
  } else {
    if (rtype == "real") {
      # don't use real Y[n]
      str_add(out$data) <- glue(
        "  vector[N{resp}] Y{resp};  // response variable\n"
      )
    } else if (rtype == "int") {
      str_add(out$data) <- glue(
        "  int Y{resp}[N{resp}];  // response variable\n"
      )
    }
  }
  if (has_ndt(family)) {
    str_add(out$tdata_def) <- glue(
      "  real min_Y{resp} = min(Y{resp});\n"
    )
  }
  if (has_trials(family) || is.formula(bterms$adforms$trials)) {
    str_add(out$data) <- glue(
      "  int trials{resp}[N{resp}];  // number of trials\n"
    )
  }
  if (is.formula(bterms$adforms$weights)) {
    str_add(out$data) <- glue(
      "  vector<lower=0>[N{resp}] weights{resp};  // model weights\n" 
    )
  }
  if (is.formula(bterms$adforms$se)) {
    str_add(out$data) <- glue(
      "  vector<lower=0>[N{resp}] se{resp};  // known sampling error\n"
    )
    str_add(out$tdata_def) <- glue(
      "  vector<lower=0>[N{resp}] se2{resp} = square(se{resp});\n"
    )
  }
  if (is.formula(bterms$adforms$dec)) {
    str_add(out$data) <- glue(
      "  int<lower=0,upper=1> dec{resp}[N{resp}];  // decisions\n"
    )
  }
  has_cens <- has_cens(bterms, data = data)
  if (has_cens) {
    str_add(out$data) <- glue(
      "  int<lower=-1,upper=2> cens{resp}[N{resp}];  // indicates censoring\n"
    )
    if (isTRUE(attr(has_cens, "interval"))) {
      rcens <- str_if(rtype == "int", 
        glue("  int rcens{resp}[N{resp}];"), 
        glue("  vector[N{resp}] rcens{resp};")
      )
      str_add(out$data) <- glue(
        rcens, "  // right censor points for interval censoring\n"
      )
    }
  }
  bounds <- trunc_bounds(bterms, data = data)
  if (any(bounds$lb > -Inf)) {
    str_add(out$data) <- glue(
      "  {rtype} lb{resp}[N{resp}];  // lower truncation bounds;\n"
    )
  }
  if (any(bounds$ub < Inf)) {
    str_add(out$data) <- glue(
      "  {rtype} ub{resp}[N{resp}];  // upper truncation bounds\n"
    )
  }
  if (is.formula(bterms$adforms$mi)) {
    Ybounds <- trunc_bounds(bterms, data, incl_family = TRUE, stan = TRUE)
    sdy <- get_sdy(bterms, data)
    if (is.null(sdy)) {
      # response is modeled without measurement error
      str_add(out$par) <- glue(
        "  vector{Ybounds}[Nmi{resp}] Ymi{resp};",
        "  // estimated missings\n"
      )
      str_add(out$data) <- glue(
        "  int<lower=0> Nmi{resp};  // number of missings\n",
        "  int<lower=1> Jmi{resp}[Nmi{resp}];",  
        "  // positions of missings\n"
      )
      str_add(out$model_def) <- glue(
        "  // vector combining observed and missing responses\n",
        "  vector[N{resp}] Yl{resp} = Y{resp};\n" 
      )
      str_add(out$model_comp_basic) <- glue(
        "  Yl{resp}[Jmi{resp}] = Ymi{resp};\n"
      )
    } else {
      str_add(out$data) <- glue(
        "  // data for measurement-error in the response\n",
        "  vector<lower=0>[N{resp}] noise{resp};\n",
        "  // information about non-missings\n",
        "  int<lower=0> Nme{resp};\n",
        "  int<lower=1> Jme{resp}[Nme{resp}];\n"
      )
      str_add(out$par) <- glue(
        "  vector{Ybounds}[N{resp}] Yl{resp};  // latent variable\n"
      )
      str_add(out$prior) <- glue(
        "  target += normal_lpdf(Y{resp}[Jme{resp}]",
        " | Yl{resp}[Jme{resp}],", 
        " noise{resp}[Jme{resp}]);\n"
      )
    }
  }
  if (is.formula(bterms$adforms$vreal)) {
    # vectors of real values for use in custom families
    k <- length(eval_rhs(bterms$adforms$vreal, data = data))
    str_add(out$data) <- cglue(
      "  // data for custom real vectors\n",
      "  vector[N{resp}] vreal{seq_len(k)}{resp};\n"
    )
  }
  if (is.formula(bterms$adforms$vint)) {
    # vectors of integer values for use in custom families
    k <- length(eval_rhs(bterms$adforms$vint, data = data))
    str_add(out$data) <- cglue(
      "  // data for custom integer vectors\n",
      "  int vint{seq_len(k)}{resp}[N{resp}];\n"
    )
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
        glue("    mu{p}[n] += head(Err{p}[n], Kar{p}) * ar{p};\n")             
      )
      add_ma <- str_if(get_ma(autocor),
        glue("    mu{p}[n] += head(Err{p}[n], Kma{p}) * ma{p};\n")             
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
      "  matrix[max(nobs_tg{p}), max(nobs_tg{p})] chol_cor;\n"               
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
      str_add(out$par) <- glue(
        "  vector[N{p}] zerr{p};  // unscaled residuals\n",
        "  real<lower=0> sderr;  // SD of residuals\n"
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
        "  target += normal_lpdf(zerr | 0, 1);\n"
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
    if (autocor$type %in% c("escar", "esicar")) {
      str_add(out$data) <- glue(
        "  vector[Nloc{p}] Nneigh{p};\n",
        "  vector[Nloc{p}] eigenW{p};\n"
      )
      str_add(out$par) <- glue(
        "  // parameters for the CAR structure\n",
        "  real<lower=0> sdcar{p};\n"
      )
      str_add(out$prior) <- stan_prior(
        prior, class = "sdcar", px = px, suffix = p
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
        "  vector[Nloc{p}] rcar{p};\n"
      )
      str_add(out$prior) <- glue(
        "  target += -0.5 * dot_self(rcar{p}[edges1{p}]", 
        " - rcar{p}[edges2{p}]);\n",
        "  // soft sum-to-zero constraint\n",
        "  target += normal_lpdf(sum(rcar{p}) | 0, 0.001 * Nloc{p});\n"
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

# define Stan functions or globally used transformed data
stan_global_defs <- function(bterms, prior, ranef, cov_ranef) {
  families <- family_names(bterms)
  links <- family_info(bterms, "link")
  unique_combs <- !duplicated(paste0(families, ":", links))
  families <- families[unique_combs]
  links <- links[unique_combs]
  out <- list()
  if (any(links == "cauchit")) {
    str_add(out$fun) <- "  #include 'fun_cauchit.stan'\n"
  } else if (any(links == "cloglog")) {
    str_add(out$fun) <- "  #include 'fun_cloglog.stan'\n"
  } else if (any(links == "softplus")) {
    str_add(out$fun) <- "  #include 'fun_softplus.stan'\n"
  }
  hs_dfs <- ulapply(attr(prior, "special"), "[[", "hs_df")
  if (any(nzchar(hs_dfs))) {
    str_add(out$fun) <- "  #include 'fun_horseshoe.stan'\n"
  }
  if (any(nzchar(ranef$by))) {
    str_add(out$fun) <- "  #include 'fun_scale_r_cor_by.stan'\n"
  }
  if (stan_needs_kronecker(ranef, names(cov_ranef))) {
    str_add(out$fun) <- glue(
      "  #include 'fun_as_matrix.stan'\n",
      "  #include 'fun_kronecker.stan'\n"
    )
  }
  family_files <- family_info(bterms, "include")
  if (length(family_files)) {
    str_add(out$fun) <- cglue("  #include '{family_files}'\n")
  }
  const <- family_info(bterms, "const")
  if (length(const)) {
    str_add(out$tdata_def) <- cglue("  {const};\n")
  }
  is_ordinal <- ulapply(families, is_ordinal)
  if (any(is_ordinal)) {
    ord_fams <- families[is_ordinal]
    ord_links <- links[is_ordinal]
    for (i in seq_along(ord_fams)) {
      str_add(out$fun) <- stan_ordinal_lpmf(ord_fams[i], ord_links[i])
    }
  }
  uni_mo <- ulapply(get_effect(bterms, "sp"), attr, "uni_mo")
  if (length(uni_mo)) {
    str_add(out$fun) <- "  #include 'fun_monotonic.stan'\n"
  } 
  if (length(get_effect(bterms, "gp"))) {
    # TODO: include functions selectively
    str_add(out$fun) <- "  #include 'fun_gaussian_process.stan'\n"
    str_add(out$fun) <- "  #include 'fun_gaussian_process_approx.stan'\n"
  }
  # functions related to autocorrelation structures
  if (is.brmsterms(bterms)) {
    autocors <- list(bterms$autocor)
  } else {
    autocors <- lapply(bterms$terms, "[[", "autocor")
  }
  if (any(ulapply(autocors, use_cov))) {
    # TODO: include functions selectively
    str_add(out$fun) <- glue(
      "  #include 'fun_normal_cov.stan'\n",
      "  #include 'fun_student_t_cov.stan'\n",
      "  #include 'fun_scale_cov_err.stan'\n",
      "  #include 'fun_cholesky_cor_ar1.stan'\n",
      "  #include 'fun_cholesky_cor_ma1.stan'\n",
      "  #include 'fun_cholesky_cor_arma1.stan'\n",
      "  #include 'fun_cholesky_cor_cosy.stan'\n"
    )
  }
  if (any(ulapply(autocors, is.cor_sar))) {
    if ("gaussian" %in% families) {
      str_add(out$fun) <- glue(
        "  #include 'fun_normal_lagsar.stan'\n",
        "  #include 'fun_normal_errorsar.stan'\n"
      )
    }
    if ("student" %in% families) {
      str_add(out$fun) <- glue(
        "  #include 'fun_student_t_lagsar.stan'\n",
        "  #include 'fun_student_t_errorsar.stan'\n"
      )
    }
  }
  if (any(ulapply(autocors, is.cor_car))) {
    str_add(out$fun) <- glue(
      "  #include 'fun_sparse_car_lpdf.stan'\n",      
      "  #include 'fun_sparse_icar_lpdf.stan'\n"
    )
  }
  out
}

# ordinal log-probability densitiy functions in Stan language
stan_ordinal_lpmf <- function(family, link) {
  stopifnot(is.character(family), is.character(link))
  ilink <- stan_ilink(link)
  th <- function(k) {
    # helper function generating stan code inside ilink(.)
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
    "   *   mu: linear predictor\n",
    "   *   thres: ordinal thresholds\n",
    "   *   disc: discrimination parameter\n",
    "   * Returns:\n", 
    "   *   a scalar to be added to the log posterior\n",
    "   */\n",
    "   real {family}_{link}_lpmf(int y, real mu, vector thres, real disc) {{\n"
  )
  # define the function body
  if (family == "cumulative") {
    str_add(out) <- glue(
      "     int ncat = num_elements(thres) + 1;\n",
      "     real p;\n",
      "     if (y == 1) {{\n",
      "       p = {ilink}({th(1)});\n",
      "     }} else if (y == ncat) {{\n",
      "       p = 1 - {ilink}({th('ncat - 1')});\n",
      "     }} else {{\n",
      "       p = {ilink}({th('y')}) -\n",
      "           {ilink}({th('y - 1')});\n",
      "     }}\n",
      "     return log(p);\n"
    )
  } else if (family %in% c("sratio", "cratio")) {
    sc <- str_if(family == "sratio", "1 - ")
    str_add(out) <- glue(
      "     int ncat = num_elements(thres) + 1;\n",
      "     vector[ncat] p;\n",
      "     vector[ncat - 1] q;\n",
      "     int k = 1;\n",
      "     while (k <= min(y, ncat - 1)) {{\n",
      "       q[k] = {sc}{ilink}({th('k')});\n",
      "       p[k] = 1 - q[k];\n",
      "       for (kk in 1:(k - 1)) p[k] = p[k] * q[kk];\n", 
      "       k += 1;\n",
      "     }}\n",
      "     if (y == ncat) {{\n",
      "       p[ncat] = prod(q);\n",
      "     }}\n",
      "     return log(p[y]);\n"
    )
  } else if (family == "acat") {
    if (ilink == "inv_logit") {
      str_add(out) <- glue(
        "     int ncat = num_elements(thres) + 1;\n",
        "     vector[ncat] p;\n",
        "     p[1] = 0.0;\n",
        "     for (k in 1:(ncat - 1)) {{\n",
        "       p[k + 1] = p[k] + {th('k')};\n",
        "     }}\n",
        "     p = exp(p);\n",
        "     return log(p[y] / sum(p));\n"
      )
    } else {
      str_add(out) <- glue(   
        "     int ncat = num_elements(thres) + 1;\n",
        "     vector[ncat] p;\n",
        "     vector[ncat - 1] q;\n",
        "     for (k in 1:(ncat - 1))\n",
        "       q[k] = {ilink}({th('k')});\n",
        "     for (k in 1:ncat) {{\n",     
        "       p[k] = 1.0;\n",
        "       for (kk in 1:(k - 1)) p[k] = p[k] * q[kk];\n",
        "       for (kk in k:(ncat - 1)) p[k] = p[k] * (1 - q[kk]);\n",   
        "     }}\n",
        "     return log(p[y] / sum(p));\n"
      )
    }
  }
  str_add(out) <- "   }\n"
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

# link function in Stan language
# @param link name of the link function
stan_link <- function(link) {
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
    log1p = "log1p",
    softplus = "log_expm1"
  )
}

# inverse link in Stan language
# @param link name of the link function
stan_ilink <- function(link) {
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
    log1p = "expm1",
    softplus = "log1p_exp"
  )
}

# define a vector in Stan language
stan_vector <- function(...) {
  paste0("[", paste0(c(...), collapse = ", "), "]'")
}

# prepare Stan code for correlations in the generated quantities block
# @param cor name of the correlation vector
# @param ncol number of columns of the correlation matrix
stan_cor_gen_comp <- function(cor, ncol) {
  Cor <- paste0(toupper(substring(cor, 1, 1)), substring(cor, 2))
  glue(
    "  // extract upper diagonal of correlation matrix\n", 
    "  for (k in 1:{ncol}) {{\n",
    "    for (j in 1:(k - 1)) {{\n",
    "      {cor}[choose(k - 1, 2) + j] = {Cor}[j, k];\n",
    "    }}\n",
    "  }}\n"
  )
}

# indicates if a family-link combination has a built in 
# function in Stan (such as binomial_logit)
# @param family a list with elements 'family' and 'link'
# @param cens_or_trunc is the model censored or truncated?
stan_has_built_in_fun <- function(family, cens_or_trunc = FALSE) {
  stopifnot(all(c("family", "link") %in% names(family)))
  link <- family$link
  dpar <- family$dpar
  family <- family$family
  if (cens_or_trunc) {
    # only few families have special lcdf and lccdf functions
    log_families <- c("cox")
    logit_families <- character(0)
    logit_dpars <- character(0)
  } else {
    log_families <- c(
      "poisson", "negbinomial", "geometric", "com_poisson",
      "zero_inflated_poisson", "zero_inflated_negbinomial",
      "hurdle_poisson", "hurdle_negbinomial", "cox"
    )
    logit_families <- c(
      "binomial", "bernoulli", "cumulative", "categorical",
      "zero_inflated_binomial"
    )
    logit_dpars <- c("zi", "hu")
  }
  isTRUE(
    family %in% log_families && link == "log" ||
    family %in% logit_families && link == "logit" ||
    isTRUE(dpar %in% logit_dpars) && link == "logit"
  )
}

# get all variable names accepted in Stan
stan_all_vars <- function(x) {
  x <- gsub("\\.", "+", x)
  all_vars(x)
}

# transform names to be used as variable names in Stan
make_stan_names <- function(x) {
  gsub("\\.|_", "", make.names(x, unique = TRUE))
}

# checks if a model needs the kronecker product
# @param ranef output of tidy_ranef
# @param names_cov_ranef: names 'cov_ranef'
# @return a single logical value
stan_needs_kronecker <- function(ranef, names_cov_ranef) {
  ids <- unique(ranef$id)
  out <- FALSE
  for (id in ids) {
    r <- ranef[ranef$id == id, ]
    out <- out || nrow(r) > 1L && r$cor[1] && r$group[1] %in% names_cov_ranef
  }
  out
}
