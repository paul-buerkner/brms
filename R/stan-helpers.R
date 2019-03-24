stan_response <- function(bterms, data) {
  # Stan code for the response variables
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
    str_add(out$tdataD) <- glue(
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
    str_add(out$tdataD) <- glue(
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
  bounds <- get_bounds(bterms, data = data)
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
    Ybounds <- get_bounds(bterms, data, incl_family = TRUE, stan = TRUE)
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
      str_add(out$modelD) <- glue(
        "  vector[N{resp}] Yl{resp} = Y{resp};\n" 
      )
      str_add(out$modelC1) <- glue(
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
  allow_autocor <- allow_autocor(family)
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
    if (!allow_autocor) {
      stop2(err_msg, " for family '", family$family, "'.") 
    }
    str_add(out$data) <- glue( 
      "  // data needed for ARMA correlations\n",
      "  int<lower=0> Kar{p};  // AR order\n",
      "  int<lower=0> Kma{p};  // MA order\n"
    )
    str_add(out$tdataD) <- glue( 
      "  int max_lag{p} = max(Kar{p}, Kma{p});\n"
    )
    # restrict ARMA effects to be in [-1,1] when using covariance
    # formulation as they cannot be outside this interval anyway
    if (Kar) {
      ar_bound <- subset2(prior, class = "ar", ls = px)$bound
      str_add(out$par) <- glue( 
        "  vector{ar_bound}[Kar{p}] ar{p};  // autoregressive effects\n"
      )
      str_add(out$prior) <- stan_prior(
        prior, class = "ar", px = px, suffix = p
      )
    }
    if (Kma) {
      ma_bound <- subset2(prior, class = "ma", ls = px)$bound
      str_add(out$par) <- glue( 
        "  vector{ma_bound}[Kma{p}] ma{p};  // moving-average effects\n"
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
      str_add(out$data) <- glue( 
        "  // see the functions block for details\n",
        "  int<lower=1> N_tg{p};\n",
        "  int<lower=1> begin_tg{p}[N_tg{p}];\n", 
        "  int<lower=1> end_tg{p}[N_tg{p}];\n", 
        "  int<lower=1> nobs_tg{p}[N_tg{p}];\n"
      )
      if (!is.formula(bterms$adforms$se)) {
        str_add(out$tdataD) <- glue(
          "  vector[N{p}] se2{p} = rep_vector(0, N{p});\n"
        )
      }
      str_add(out$tparD) <- glue(
        "  matrix[max(nobs_tg{p}), max(nobs_tg{p})] res_cov_matrix;\n"               
      )
      if (Kar && !Kma) {
        cov_mat_fun <- "ar1"
        cov_mat_args <- glue("ar{p}[1]")
      } else if (!Kar && Kma) {
        cov_mat_fun <- "ma1"
        cov_mat_args <- glue("ma{p}[1]")
      } else {
        cov_mat_fun <- "arma1"
        cov_mat_args <- glue("ar{p}[1], ma{p}[1]")
      }
      str_add(out$tparC1) <- glue(
        "  // compute residual covariance matrix\n",
        "  res_cov_matrix = cov_matrix_{cov_mat_fun}", 
        "({cov_mat_args}, sigma{p}, max(nobs_tg{p}));\n"
      )
    } else {
      err_msg <- "Please set cov = TRUE in cor_arma / cor_ar / cor_ma"
      if (is.formula(bterms$adforms$se)) {
        stop2(err_msg, " when specifying 'se'.")
      }
      if (is.btnl(bterms$dpars[["mu"]])) {
        stop2(err_msg, " in non-linear models.")
      }
      if (!identical(family$link, "identity")) {
        stop2(err_msg, " when using non-identity links.")
      }
      str_add(out$data) <- glue( 
        "  int<lower=0> J_lag{p}[N{p}];\n"                
      )
      str_add(out$modelD) <- glue(
        "  // objects storing residuals\n",
        "  matrix[N{p}, max_lag{p}] E{p}",
        " = rep_matrix(0, N{p}, max_lag{p});\n",
        "  vector[N{p}] e{p};\n"
      )
      Y <- str_if(is.formula(bterms$adforms$mi), "Yl", "Y")
      str_add(out$modelC2) <- glue(
        "    // computation of ARMA correlations\n",
        "    e{p}[n] = {Y}{p}[n] - mu{p}[n];\n",
        "    for (i in 1:J_lag{p}[n]) {{\n",
        "      E{p}[n + 1, i] = e{p}[n + 1 - i];\n",
        "    }}\n"
      )
    } 
  }
  Karr <- get_arr(autocor)
  if (Karr) {
    # autoregressive effects of the response
    warning2(
      "The 'arr' correlation structure has been deprecated and ",
      "will be removed from the package at some point. Consider ", 
      "using lagged response values as ordinary predictors instead."
    )
    err_msg <- "ARR models are not implemented"
    if (length(bterms$dpars[["mu"]]$nlpars)) {
      stop2(err_msg, " for non-linear models.")
    }
    if (is.mixfamily(family)) {
      stop2(err_msg, " for mixture models.") 
    }
    str_add(out$data) <- glue(
      "  // data needed for ARR correlations\n",
      "  int<lower=1> Karr{p};\n",
      "  matrix[N{p}, Karr{p}] Yarr{p};  // ARR design matrix\n"
    )
    bound <- subset2(prior, class = "arr", ls = px)$bound
    str_add(out$par) <- glue(
      "  vector{bound}[Karr{p}] arr{p};  // autoregressive response effects\n"
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
    if (!allow_autocor) {
      stop2(err_msg, " for family '", family$family, "'.")
    }
    if (isTRUE(bterms$rescor)) {
      stop2(err_msg, " when 'rescor' is estimated.")
    }
    if (any(c("sigma", "nu") %in% names(bterms$dpars))) {
      stop2(err_msg, " when predicting 'sigma' or 'nu'.")
    }
    str_add(out$data) <- glue(
      "  matrix[N{p}, N{p}] W{p};  // spatial weight matrix\n"                  
    )
    if (identical(autocor$type, "lag")) {
      str_add(out$par) <- glue( 
        "  real<lower=0,upper=1> lagsar{p};  // SAR parameter\n"
      )
      str_add(out$prior) <- stan_prior(
        prior, class = "lagsar", px = px, suffix = p
      )
    } else if (identical(autocor$type, "error")) {
      str_add(out$par) <- glue( 
        "  real<lower=0,upper=1> errorsar{p};  // SAR parameter\n"
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
      str_add(out$tparD) <- glue(
        "  vector[Nloc{p}] rcar{p};\n"                
      )
      str_add(out$tparC1) <- glue(
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
  if (is.cor_bsts(autocor)) {
    warning2(
      "The 'bsts' correlation structure has been deprecated and ",
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
    str_add(out$data) <- glue(
      "  vector[N{p}] tg{p};  // indicates independent groups\n"
    )
    str_add(out$par) <- glue(
      "  vector[N{p}] loclev{p};  // local level terms\n",
      "  real<lower=0> sigmaLL{p};  // SD of local level terms\n"
    )
    if (allow_autocor) {
      # ensures that the means of the loclev priors start close 
      # to the response values; this often helps with convergence
      link <- stan_link(family$link)
      center <- c("1", "n")
      center <- glue("{link}(Y{p}[{center}])")
    } else {
      center <- c("0", "0")
    }
    str_add(out$prior) <- glue(
      stan_prior(prior, class = "sigmaLL", px = px, suffix = p),
      "  target += normal_lpdf(loclev{p}[1] | {center[1]}, sigmaLL{p});\n",
      "  for (n in 2:N{p}) {{\n",
      "    if (tg{p}[n] == tg{p}[n - 1]) {{\n",
      "    {tp()}normal_lpdf(loclev{p}[n] | loclev{p}[n - 1], sigmaLL{p});\n",
      "    }} else {{\n",
      "    {tp()}normal_lpdf(loclev{p}[n] | {center[2]}, sigmaLL{p});\n",
      "    }}\n",
      "  }}\n"
    )
  }
  if (is.cor_fixed(autocor)) {
    err_msg <- "Fixed residual covariance matrices are not implemented"
    if (is.mixfamily(family)) {
      stop2(err_msg, " for mixture models.") 
    }
    if (!allow_autocor) {
      stop2(err_msg, " for family '", family$family, "'.")
    }
    if (isTRUE(bterms$rescor)) {
      stop2(err_msg, " when 'rescor' is estimated.")
    }
    str_add(out$data) <- glue( 
      "  matrix[N{p}, N{p}] V{p};  // known residual covariance matrix\n"
    )
    if (family$family %in% "gaussian") {
      str_add(out$tdataD) <- glue(
        "  matrix[N{p}, N{p}] LV{p} = cholesky_decompose(V{p});\n"
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
    str_add(out$tdataD) <- cglue("  {const};\n")
  }
  is_ordinal <- ulapply(families, is_ordinal)
  if (any(is_ordinal)) {
    ord_fams <- families[is_ordinal]
    ord_links <- links[is_ordinal]
    for (i in seq_along(ord_fams)) {
      for (cs in 0:1) {
        str_add(out$fun) <- stan_ordinal_lpmf(ord_fams[i], ord_links[i], cs)
      }
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
    str_add(out$fun) <- glue(
      "  #include 'fun_normal_cov.stan'\n",
      "  #include 'fun_student_t_cov.stan'\n",
      "  #include 'fun_cov_matrix_ar1.stan'\n",
      "  #include 'fun_cov_matrix_ma1.stan'\n",
      "  #include 'fun_cov_matrix_arma1.stan'\n"
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

stan_ordinal_lpmf <- function(family, link, cs = FALSE) {
  # ordinal log-probability densitiy functions in Stan language
  # Args:
  #   cs: Logical; add category specific effects?
  stopifnot(is.character(family), is.character(link))
  cs <- as_one_logical(cs)
  ilink <- stan_ilink(link)
  th <- function(k) {
    # helper function generating stan code inside ilink(.)
    sign <- str_if(family %in% c("cumulative", "sratio"), " - ", " + ")
    ptl <- str_if(cs, glue("{sign}mucs[{k}]")) 
    if (sign == " - ") {
      out <- glue("thres[{k}]{ptl} - mu")
    } else {
      out <- glue("mu{ptl} - thres[{k}]")
    }
    glue("disc * ({out})")
  }
  mucs <- str_if(cs, " row_vector mucs,")
  out <- glue(
    "  /* {family}-{link} log-PDF for a single response\n",
    str_if(cs, "   * including category specific effects\n"),
    "   * Args:\n",
    "   *   y: response category\n",
    "   *   mu: linear predictor\n",
    str_if(cs, "   *   mucs: predictor for category specific effects\n"),
    "   *   thres: ordinal thresholds\n",
    "   *   disc: discrimination parameter\n",
    "   * Returns:\n", 
    "   *   a scalar to be added to the log posterior\n",
    "   */\n",
    "   real {family}_{link}", str_if(cs, "_cs"), 
    "_lpmf(int y, real mu,{mucs} vector thres, real disc) {{\n"
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
    def_thetas <- cglue(
      "  real<lower=0,upper=1> theta{1:nmix}{p};  // mixing proportion\n"
    )
    if (length(theta_pred)) {
      if (length(theta_pred) != nmix - 1) {
        stop2("Can only predict all but one mixing proportion.")
      }
      missing_id <- setdiff(1:nmix, dpar_id(names(theta_pred)))
      str_add(out$modelD) <- glue(
        "  vector[N{p}] theta{missing_id}{p} = rep_vector(0, N{p});\n",                   
        "  real log_sum_exp_theta;\n"      
      )
      sum_exp_theta <- glue("exp(theta{1:nmix}{p}[n])", collapse = " + ")
      str_add(out$modelC3) <- glue(
        "    log_sum_exp_theta = log({sum_exp_theta});\n"
      )
      str_add(out$modelC3) <- cglue(
        "    theta{1:nmix}{p}[n] = theta{1:nmix}{p}[n] - log_sum_exp_theta;\n"
      )
    } else if (length(theta_fix)) {
      if (length(theta_fix) != nmix) {
        stop2("Can only fix no or all mixing proportions.")
      }
      str_add(out$data) <- "  // mixing proportions\n"
      str_add(out$data) <- cglue(
        "  real<lower=0,upper=1> theta{1:nmix}{p};\n"
      )
    } else {
      str_add(out$data) <- glue(
        "  vector[{nmix}] con_theta{p};  // prior concentration\n"                  
      )
      str_add(out$par) <- glue(
        "  simplex[{nmix}] theta{p};  // mixing proportions\n"
      )
      str_add(out$prior) <- glue(
        "  target += dirichlet_lpdf(theta{p} | con_theta{p});\n"                
      )
      str_add(out$tparD) <- "  // mixing proportions\n"
      str_add(out$tparD) <- cglue(
        "  real<lower=0,upper=1> theta{1:nmix}{p} = theta{p}[{1:nmix}];\n"
      )
    }
    if (bterms$family$order %in% "mu") {
      str_add(out$par) <- glue( 
        "  ordered[{nmix}] ordered_Intercept{p};  // to identify mixtures\n"
      )
    }
  }
  out
}

stan_Xme <- function(meef, prior) {
  # global Stan definitions for noise-free variables
  # Args:
  #   meef: tidy data.frame as returned by tidy_meef()
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
        "  int<lower=0> Nme_{i};\n",
        "  int<lower=1> Jme_{i}[N];\n"
      )
    } else {
      Nme <- "N"
    }
    str_add(out$data) <- glue(
      "  int<lower=1> Mme_{i};\n"
    )
    str_add(out$data) <- cglue(
      "  vector[{Nme}] Xn_{K};\n",
      "  vector<lower=0>[{Nme}] noise_{K};\n"
    )
    str_add(out$par) <- cglue(
      "  vector[Mme_{i}] meanme_{i};\n",
      "  vector<lower=0>[Mme_{i}] sdme_{i};\n"
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
        "  int<lower=1> NCme_{i};\n"
      )
      str_add(out$par) <- glue(
        "  matrix[Mme_{i}, {Nme}] zme_{i};\n",
        "  cholesky_factor_corr[Mme_{i}] Lme_{i};\n"
      )
      str_add(out$tparD) <- glue(
        "  matrix[{Nme}, Mme_{i}] Xme{i}", 
        " = rep_matrix(meanme_{i}', {Nme}) ", 
        " + (diag_pre_multiply(sdme_{i}, Lme_{i}) * zme_{i})';\n"
      )
      str_add(out$tparD) <- cglue(
        "  vector[{Nme}] Xme_{K} = Xme{i}[, {K}];\n"
      )
      str_add(out$prior) <- glue(
        "  target += normal_lpdf(to_vector(zme_{i}) | 0, 1);\n",
        stan_prior(prior, "Lme", group = g, suffix = usc(i))
      )
      str_add(out$genD) <- cglue(
        "  corr_matrix[Mme_{i}] Corme_{i}", 
        " = multiply_lower_tri_self_transpose(Lme_{i});\n",
        "  vector<lower=-1,upper=1>[NCme_{i}] corme_{i};\n"
      )
      str_add(out$genC) <- stan_cor_genC(glue("corme_{i}"), glue("Mme_{i}"))
    } else {
      str_add(out$par) <- cglue(
        "  vector[{Nme}] zme_{K};\n"
      )
      str_add(out$tparD) <- cglue(
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
    log1p = "log1p",
    softplus = "log_expm1"
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
    log1p = "expm1",
    softplus = "log1p_exp"
  )
}

stan_vector <- function(...) {
  # define a vector in Stan language
  paste0("[", paste0(c(...), collapse = ", "), "]'")
}

stan_cor_genC <- function(cor, ncol) {
  # prepare Stan code for correlations in the generated quantities block
  # Args:
  #   cor: name of the correlation vector
  #   ncol: number of columns of the correlation matrix
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
    "poisson", "negbinomial", "geometric", "com_poisson",
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

stan_all_vars <- function(x) {
  # get all variable names accepted in Stan
  x <- gsub("\\.", "+", x)
  all_vars(x)
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
    out <- out || nrow(r) > 1L && r$cor[1] && r$group[1] %in% names_cov_ranef
  }
  out
}
