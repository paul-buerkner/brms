stan_predictor <- function(x, ...) {
  # generate stan code for predictor terms
  UseMethod("stan_predictor")
}

#' @export
stan_predictor.btl <- function(x, data, prior, ranef, ilink = rep("", 2), ...) {
  # combine effects for the predictors of a single (non-linear) parameter
  # Args:
  #   ranef: output of tidy_ranef()
  #   ilink: character vector of length 2 defining the link to be applied
  #   ...: arguements passed to the underlying effect-specific functions
  stopifnot(length(ilink) == 2L)
  out <- collapse_lists(
    text_fe <- stan_fe(x, data, prior, ...),
    text_sp <- stan_sp(x, data, prior, ranef = ranef, ...),
    text_cs <- stan_cs(x, data, prior, ranef = ranef, ...),
    text_sm <- stan_sm(x, data, prior, ...),
    text_gp <- stan_gp(x, data, prior, ...),
    text_ac <- stan_ac(x, ...),
    text_offset <- stan_offset(x, ...),
    stan_special_prior_global(x, data, prior, ...)
  )
  
  # initialize and compute eta
  px <- check_prefix(x)
  eta <- combine_prefix(px, keep_mu = TRUE, nlp = TRUE)
  out$eta <- sub("^[ \t\r\n]+\\+", "", out$eta, perl = TRUE)
  str_add(out$modelD) <- glue("  vector[N] {eta} ={out$eta};\n")
  str_add(out$loopeta) <- stan_eta_re(ranef, px = px)
  if (nzchar(out$loopeta)) {
    # parts of eta are computed in a loop over observations
    out$loopeta <- sub("^[ \t\r\n]+\\+", "", out$loopeta, perl = TRUE)
    str_add(out$modelC2) <- glue("    {eta}[n] +={out$loopeta};\n")
  }
  out$eta <- out$loopeta <- NULL
  
  # possibly transform eta before it is passed to the likelihood
  if (sum(nzchar(ilink))) {
    # make sure mu comes last as it might depend on other parameters
    not_mu <- nzchar(x$dpar) && dpar_class(x$dpar) != "mu"
    position <- str_if(not_mu, "modelC3", "modelC4")
    str_add(out[[position]]) <- glue(
      "    {eta}[n] = {ilink[1]}{eta}[n]{ilink[2]};\n"
    )
  }
  out
}

#' @export
stan_predictor.btnl <- function(x, data, nlpars, ilink = rep("", 2), ...) {
  # prepare Stan code for non-linear models
  # Args:
  #   data: data.frame supplied by the user
  #   ilink: character vector of length 2 containing
  #     Stan code for the link function
  #   ...: passed to stan_predictor.btl
  stopifnot(length(ilink) == 2L)
  out <- list()
  par <- combine_prefix(x, keep_mu = TRUE, nlp = TRUE)
  # prepare non-linear model
  n <- str_if(x$loop, "[n] ", " ")
  new_nlpars <- glue(" nlp{usc(x$resp)}_{nlpars}{n}")
  # covariates in the non-linear model
  covars <- wsp(all.vars(rhs(x$covars)))
  new_covars <- NULL
  if (length(covars)) {
    p <- usc(combine_prefix(x))
    covar_names <- glue("C{p}_{seq_along(covars)}")
    new_covars <- glue(" {covar_names}{n}")
    # use vectors as indexing matrices in Stan is slow
    str_add(out$data) <- glue( 
      "  // covariate vectors\n",
      cglue("  vector[N] {covar_names};\n")
    )
  }
  # add whitespaces to be able to replace parameters and covariates
  syms <- c(
    "+", "-", "*", "/", "%", "^", ".*", "./", "'", ")", "(", 
    ",", "==", "!=", "<=", ">=", "<", ">", "!", "&&", "||" 
  )
  regex <- glue("(?<!\\.){escape_all(syms)}(?!=)")
  nlmodel <- rm_wsp(collapse(deparse(x$formula[[2]])))
  nlmodel <- wsp(rename(nlmodel, regex, wsp(syms), fixed = FALSE, perl = TRUE)) 
  nlmodel <- rename(nlmodel, 
    c(wsp(nlpars), covars, " ( ", " ) "), 
    c(new_nlpars, new_covars, "(", ")")
  )
  # possibly transform eta in the transformed params block
  str_add(out$modelD) <- glue("  vector[N] {par};\n")
  if (x$loop) {
    str_add(out$modelC4) <- glue(
      "    // compute non-linear predictor\n",
      "    {par}[n] = {ilink[1]}{trimws(nlmodel)}{ilink[2]};\n"
    )
  } else {
    str_add(out$modelC5) <- glue(
      "  // compute non-linear predictor\n",
      "  {par} = {ilink[1]}{trimws(nlmodel)}{ilink[2]};\n"
    )
  }
  out
}

#' @export
stan_predictor.brmsterms <- function(x, data, prior, sparse = FALSE, 
                                     rescor = FALSE, ...) {
  # Stan code for distributional parameters
  # Args:
  #   rescor: indicate if this is part of an MV model estimating rescor
  px <- check_prefix(x)
  resp <- usc(combine_prefix(px))
  out <- list(stan_response(x, data = data))
  valid_dpars <- valid_dpars(x)
  args <- nlist(
    data, prior, sparse, nlpars = names(x$nlpars),
    order_mixture = x$family$order, ...
  )
  for (nlp in names(x$nlpars)) {
    nlp_args <- list(x$nlpars[[nlp]], center_X = FALSE)
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
      out[[dp]] <- list(tparD = dp_def, tparC1 = dp_co)
    } else {
      dp_def <- stan_dpar_defs(dp, resp, family = x$family)
      dp_def_temp <- stan_dpar_defs_temp(dp, resp, family = x$family)
      if (nzchar(dp_def_temp)) {
        dp_prior <- stan_prior(
          prior, dp, prefix = "temp_", suffix = resp, px = px
        )
        out[[dp]] <- list(par = dp_def_temp, prior = dp_prior)
      } else if (nzchar(dp_def)) {
        dp_prior <- stan_prior(prior, dp, suffix = resp, px = px)
        out[[dp]] <- list(par = dp_def, prior = dp_prior)
      }
    }
  }
  out <- lc(out,
    stan_autocor(x, prior = prior),
    stan_mixture(x, prior = prior),
    stan_dpar_transform(x)
  )
  collapse_lists(ls = out)
}

#' @export
stan_predictor.mvbrmsterms <- function(x, prior, ...) {
  out <- lapply(x$terms, stan_predictor, prior = prior, rescor = x$rescor, ...)
  out <- collapse_lists(ls = out)
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
    resp <- x$responses
    nresp <- length(resp)
    str_add(out$modelD) <- glue( 
      "  // multivariate linear predictor matrix\n",
      "  vector[nresp] Mu[N];\n"
    )
    str_add(out$modelC4) <- glue(
      "    Mu[n] = {stan_vector(glue('mu_{resp}[n]'))};\n"
    )
    str_add(out$data) <- glue(
      "  int<lower=1> nresp;  // number of responses\n",   
      "  int nrescor;  // number of residual correlations\n"
    )
    str_add(out$tdataD) <- glue(
      "  vector[nresp] Y[N];  // response matrix\n"
    )
    str_add(out$tdataC) <- glue(
      "  for (n in 1:N) {{\n",
      "    Y[n] = {stan_vector(glue('Y_{resp}[n]'))};\n",
      "  }}\n"
    )
    if (any(adnames %in% "weights")) {
      str_add(out$tdataD) <- glue(
        "  vector<lower=0>[N] weights = weights_{resp[1]};\n" 
      )
    }
    miforms <- rmNULL(lapply(adforms, "[[", "mi"))
    if (length(miforms)) {
      str_add(out$modelD) <- "  vector[nresp] Yl[N] = Y;\n"
      for (i in seq_along(miforms)) {
        j <- match(names(miforms)[i], resp)
        str_add(out$modelC2) <- glue(
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
      str_add(out$modelD) <- "  vector[nresp] sigma[N];\n"
      str_add(out$modelC4) <- glue(
        "    sigma[n] = {stan_vector(sigma)};\n"
      )
      if (family == "gaussian") {
        str_add(out$modelD) <- glue(
          "  // cholesky factor of residual covariance matrix\n",
          "  matrix[nresp, nresp] LSigma[N];\n"
        )
        str_add(out$modelC4) <- glue(
          "    LSigma[n] = diag_pre_multiply(sigma[n], Lrescor);\n"
        )
      } else if (family == "student") {
        str_add(out$modelD) <- glue(
          "  // residual covariance matrix\n",
          "  matrix[nresp, nresp] Sigma[N];\n"
        )
        str_add(out$modelC4) <- glue(
          "    Sigma[n] = multiply_lower_tri_self_transpose(", 
          "diag_pre_multiply(sigma[n], Lrescor));\n" 
        )
      }
    } else {
      str_add(out$modelD) <- glue(
        "  vector[nresp] sigma = {stan_vector(sigma)};\n"
      )
      if (family == "gaussian") {
        str_add(out$modelD) <- glue(
          "  // cholesky factor of residual covariance matrix\n",
          "  matrix[nresp, nresp] LSigma = ",
          "diag_pre_multiply(sigma, Lrescor);\n"
        )
      } else if (family == "student") {
        str_add(out$modelD) <- glue(
          "  // residual covariance matrix\n",
          "  matrix[nresp, nresp] Sigma = ",
          "multiply_lower_tri_self_transpose(", 
          "diag_pre_multiply(sigma, Lrescor));\n"
        )
      }
    }
    str_add(out$genD) <- glue(
      "  matrix[nresp, nresp] Rescor",
      " = multiply_lower_tri_self_transpose(Lrescor);\n",
      "  vector<lower=-1,upper=1>[nrescor] rescor;\n"
    )
    rescor_genC <- ulapply(seq_len(nresp), function(i) 
      lapply(seq_len(i - 1), function(j) glue(
        "  rescor[{(i - 1) * (i - 2) / 2 + j}] = Rescor[{j}, {i}];\n"
      ))
    )
    str_add(out$genC) <- glue(
      "  // take only relevant parts of residual correlation matrix\n",
      collapse(rescor_genC)
    )
  }
  out
}

stan_fe <- function(bterms, data, prior, stanvars,
                    center_X = TRUE, sparse = FALSE, 
                    order_mixture = 'none', ...) {
  # Stan code for population-level effects
  # Args:
  #   center_X: center the design matrix?
  #   sparse: should the design matrix be treated as sparse?
  #   order_mixture: order intercepts to identify mixture models?
  # Returns:
  #   a list containing Stan code related to population-level effects
  out <- list()
  family <- bterms$family
  fixef <- colnames(data_fe(bterms, data)$X)
  center_X <- center_X && has_intercept(bterms$fe) && 
    !is.cor_bsts(bterms$autocor) && !sparse
  rm_intercept <- center_X || is_ordinal(family) ||
    is.cor_bsts(bterms$autocor)
  if (rm_intercept) {
    fixef <- setdiff(fixef, "Intercept")
  }
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  resp <- usc(bterms$resp)
  if (length(fixef)) {
    str_add(out$data) <- glue( 
      "  int<lower=1> K{p};",
      "  // number of population-level effects\n", 
      "  matrix[N, K{p}] X{p};",
      "  // population-level design matrix\n"
    )
    if (sparse) {
      stopifnot(!center_X)
      str_add(out$tdataD) <- glue(
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
    prefix <- combine_prefix(px, keep_mu = TRUE)
    special <- attr(prior, "special")[[prefix]]
    define_b_in_pars <- is.null(special[["hs_df"]]) &&
      !glue("b{p}") %in% names(stanvars)
    if (define_b_in_pars) {
      ct <- str_if(center_X, "c")
      bound <- get_bound(prior, class = "b", px = px)
      str_add(out$par) <- glue(
        "  vector{bound}[K{ct}{p}] b{p};  // population-level effects\n"
      )
    }
    str_add(out$prior) <- stan_prior(
      prior, class = "b", coef = fixef, px = px, suffix = p
    )
    out <- collapse_lists(out,
      stan_special_prior_local(
        "b", prior, ncoef = length(fixef), 
        px = px, center_X = center_X  
      )                      
    )
  }
  
  if (center_X) {
    # centering of the fixed effects design matrix improves convergence
    if (length(fixef)) {
      str_add(out$tdataD) <- glue(
        "  int Kc{p} = K{p} - 1;\n",
        "  matrix[N, K{p} - 1] Xc{p};", 
        "  // centered version of X{p}\n",
        "  vector[K{p} - 1] means_X{p};",
        "  // column means of X{p} before centering\n"
      )
      str_add(out$tdataC) <- glue(
        "  for (i in 2:K{p}) {{\n",
        "    means_X{p}[i - 1] = mean(X{p}[, i]);\n",
        "    Xc{p}[, i - 1] = X{p}[, i] - means_X{p}[i - 1];\n",
        "  }}\n"
      )
      # ordinal families either use thres - mu or mu - thres
      # both implies adding <mean_X, b> to the temporary intercept
      sign <- str_if(is_ordinal(family), " + ", " - ")
      sub_X_means <- glue("{sign}dot_product(means_X{p}, b{p})")
    } else {
      sub_X_means <- ""
    }
    if (is_ordinal(family)) {
      # intercepts in ordinal models require special treatment
      type <- str_if(family$family == "cumulative", "ordered", "vector")
      intercept <- glue(
        "  {type}[ncat{resp}-1] temp{p}_Intercept;",
        "  // temporary thresholds\n"
      )
      if (family$threshold == "equidistant") {
        bound <- subset2(prior, class = "delta", ls = px)$bound
        str_add(out$par) <- glue(
          "  real temp{p}_Intercept1;  // first threshold\n",
          "  real{bound} delta{p};  // distance between thresholds\n"
        )
        str_add(out$tparD) <- intercept
        str_add(out$tparC1) <- glue(
          "  // compute equidistant thresholds\n",
          "  for (k in 1:(ncat{resp} - 1)) {{\n",
          "    temp{p}_Intercept[k] = temp{p}_Intercept1", 
          " + (k - 1.0) * delta{p};\n",
          "  }}\n"
        )
        str_add(out$prior) <- stan_prior(prior, class = "delta", px = px)
      } else {
        str_add(out$par) <- intercept
      }
      str_add(out$genD) <- glue(
        "  // compute actual thresholds\n",
        "  vector[ncat{resp} - 1] b{p}_Intercept",  
        " = temp{p}_Intercept{sub_X_means};\n" 
      )
    } else {
       if (identical(dpar_class(px$dpar), order_mixture)) {
         # identify mixtures via ordering of the intercepts
         dp_id <- dpar_id(px$dpar)
         str_add(out$tparD) <- glue(
           "  // identify mixtures via ordering of the intercepts\n",                   
           "  real temp{p}_Intercept = ordered_Intercept{resp}[{dp_id}];\n"
         )
      } else {
        str_add(out$par) <- glue(
          "  real temp{p}_Intercept;  // temporary intercept\n"
        )
      }
      str_add(out$eta) <- glue(" + temp{p}_Intercept")
      str_add(out$genD) <- glue(
        "  // actual population-level intercept\n",
        "  real b{p}_Intercept = temp{p}_Intercept{sub_X_means};\n"
      )
    }
    # for equidistant thresholds only temp_Intercept1 is a parameter
    prefix <- glue("temp{p}_")
    Icoefs <- suffix <- ""
    if (is_ordinal(family)) {
      Icoefs <- subset2(prior, class = "Intercept", ls = px)$coef
      Icoefs <- Icoefs[nzchar(Icoefs)]
      if (family$threshold == "equidistant") suffix <- "1"
    }
    str_add(out$prior) <- stan_prior(
      prior, class = "Intercept", coef = Icoefs, 
      px = px, prefix = prefix, suffix = suffix
    )
  } else {
    if (identical(dpar_class(px$dpar), order_mixture)) {
      stop2(
        "Identifying mixture components via ordering requires ",
        "population-level intercepts to be present.\n",
        "Try setting order = 'none' in function 'mixture'."
      )
    }
  }
  str_add(out$eta) <- stan_eta_fe(fixef, center_X, sparse, px = px)
  out
}

stan_re <- function(ranef, prior, ...) {
  # Stan code for group-level effects
  # the ID syntax requires group-level effects to be evaluated separately
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
      str_add(out$tparD) <- glue(
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
  tmp <- lapply(IDs, .stan_re, ranef = ranef, prior = prior, ...)
  out <- collapse_lists(ls = c(list(out), tmp))
  out
}

.stan_re <- function(id, ranef, prior, cov_ranef = NULL) {
  # Stan code for group-level effects per ID
  # Args:
  #   id: the ID of the grouping factor
  #   ranef: a data.frame returned by tidy_ranef
  #   prior: object of class brmsprior
  #   cov_ranef: a list of custom covariance matrices 
  # Returns:
  #   A list of strings containing Stan code
  out <- list()
  r <- subset2(ranef, id = id)
  has_ccov <- r$group[1] %in% names(cov_ranef)
  has_by <- nzchar(r$by[[1]])
  Nby <- seq_along(r$bylevels[[1]]) 
  ng <- seq_along(r$gcall[[1]]$groups)
  px <- check_prefix(r)
  idp <- paste0(r$id, usc(combine_prefix(px)))
  # define data needed for group-level effects
  str_add(out$data) <- glue(
    "  // data for group-level effects of ID {id}\n",
    "  int<lower=1> N_{id};\n",
    "  int<lower=1> M_{id};\n"
  )
  if (r$gtype[1] == "mm") {
    str_add(out$data) <- cglue(
      "  int<lower=1> J_{id}_{ng}[N];\n",
      "  real W_{id}_{ng}[N];\n"
    )
  } else {
    str_add(out$data) <- glue(
      "  int<lower=1> J_{id}[N];\n"
    )
  }
  if (has_by) {
    str_add(out$data) <- glue(
      "  int<lower=1> Nby_{id};\n",
      "  int<lower=1> Jby_{id}[N_{id}];\n"
    )
  }
  if (has_ccov) {
    str_add(out$data) <- glue(
      "  // cholesky factor of known covariance matrix\n",
      "  matrix[N_{id}, N_{id}] Lcov_{id};\n"
    )
  }
  J <- seq_rows(r)
  needs_Z <- !r$type %in% "sp"
  if (any(needs_Z)) {
    if (r$gtype[1] == "mm") {
      for (i in which(needs_Z)) {
        str_add(out$data) <- cglue(
          "  vector[N] Z_{idp[i]}_{r$cn[i]}_{ng};\n"
        )
      }
    } else {
      str_add(out$data) <- cglue(
        "  vector[N] Z_{idp[needs_Z]}_{r$cn[needs_Z]};\n"
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
      "  int<lower=1> NC_{id};\n"
    )
    str_add(out$par) <- glue(
      "  matrix[M_{id}, N_{id}] z_{id};",
      "  // unscaled group-level effects\n"
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
      str_add(out$tparD) <- glue(
        "  // group-level effects\n",
        "  matrix[N_{id}, M_{id}] r_{id}", 
        " = {dfm}scale_r_cor_by(z_{id}, sd_{id}, L_{id}, Jby_{id});\n"
      )
      str_add(out$prior) <- stan_prior(
        prior, class = "L", group = r$group[1],
        suffix = glue("_{id}[{Nby}]")
      )
      str_add(out$genD) <- cglue(
        "  corr_matrix[M_{id}] Cor_{id}_{Nby}",
        " = multiply_lower_tri_self_transpose(L_{id}[{Nby}]);\n",
        "  vector<lower=-1,upper=1>[NC_{id}] cor_{id}_{Nby};\n"
      )
      str_add(out$genC) <- stan_cor_genC(nrow(r), glue("{id}_{Nby}"))
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
      str_add(out$tparD) <- glue(
        "  // group-level effects\n",
        "  matrix[N_{id}, M_{id}] r_{id} = {dfm}{rdef}"
      )
      str_add(out$prior) <- stan_prior(
        prior, class = "L", group = r$group[1], suffix = usc(id)
      )
      str_add(out$genD) <- glue(
        "  corr_matrix[M_{id}] Cor_{id}",
        " = multiply_lower_tri_self_transpose(L_{id});\n",
        "  vector<lower=-1,upper=1>[NC_{id}] cor_{id};\n"
      )
      str_add(out$genC) <- stan_cor_genC(nrow(r), id)
    }
    str_add(out$tparD) <- cglue(
        "  vector[N_{id}] r_{idp}_{r$cn} = r_{id}[, {J}];\n"
    )
  } else {
    # single or uncorrelated group-level effects
    str_add(out$par) <- glue(
      "  vector[N_{id}] z_{id}[M_{id}];  // unscaled group-level effects\n"
    )
    str_add(out$prior) <- cglue(
      "  target += normal_lpdf(z_{id}[{seq_rows(r)}] | 0, 1);\n"
    )
    str_add(out$tparD) <- "  // group-level effects\n" 
    Lcov <- str_if(has_ccov, glue("Lcov_{id} * "))
    if (has_rows(tr)) {
      dfm <- glue("dfm_{tr$ggn[1]} .* ")
    }
    if (has_by) {
      str_add(out$tparD) <- cglue(
        "  vector[N_{id}] r_{idp}_{r$cn}",
        " = {dfm}(sd_{id}[{J}, Jby_{id}]' .* ({Lcov}z_{id}[{J}]));\n"
      )
    } else {
      str_add(out$tparD) <- cglue(
        "  vector[N_{id}] r_{idp}_{r$cn}",
        " = {dfm}(sd_{id}[{J}] * ({Lcov}z_{id}[{J}]));\n"
      )
    }
  }
  out
}

stan_sm <- function(bterms, data, prior, ...) {
  # Stan code of smooth terms
  out <- list()
  smef <- tidy_smef(bterms, data)
  if (!NROW(smef)) {
    return(out)
  }
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  Xs_names <- attr(smef, "Xs_names")
  if (length(Xs_names)) {
    str_add(out$data) <- glue(
      "  // data for smooth terms\n",
      "  int Ks{p};\n",
      "  matrix[N, Ks{p}] Xs{p};\n"
    )
    str_add(out$pars) <- glue(
      "  // parameters for smooth terms\n",
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
      "  // data of smooth {smef$byterm[i]}\n",  
      "  int nb{pi};  // number of bases\n",
      "  int knots{pi}[nb{pi}];\n"
    )
    str_add(out$data) <- cglue(
      "  matrix[N, knots{pi}[{nb}]] Zs{pi}_{nb};\n"
    )
    str_add(out$par) <- glue(
      "  // parameters of smooth {smef$byterm[i]}\n"
    )
    str_add(out$par) <- cglue(
      "  vector[knots{pi}[{nb}]] zs{pi}_{nb};\n",
      "  real<lower=0> sds{pi}_{nb};\n"
    )
    str_add(out$tparD) <- cglue(
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

stan_cs <- function(bterms, data, prior, ranef, ...) {
  # Stan code for category specific effects
  # (!) Not implemented for non-linear models
  out <- list()
  csef <- colnames(get_model_matrix(bterms$cs, data))
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  resp <- usc(bterms$resp)
  ranef <- subset2(ranef, type = "cs", ls = px)
  if (length(csef)) {
    str_add(out$data) <- glue(
      "  int<lower=1> Kcs{p};  // number of category specific effects\n",
      "  matrix[N, Kcs{p}] Xcs{p};  // category specific design matrix\n"
    )
    bound <- get_bound(prior, class = "b", px = px)
    str_add(out$par) <- glue(
      "  matrix{bound}[Kcs{p}, ncat{resp} - 1] bcs{p};",
      "  // category specific effects\n"
    )
    str_add(out$modelD) <- glue(
      "  // linear predictor for category specific effects\n",
      "  matrix[N, ncat{resp} - 1] mucs{p} = Xcs{p} * bcs{p};\n"
    ) 
    str_add(out$prior) <- stan_prior(
      prior, class = "b", coef = csef,
      suffix = "cs", px = px, matrix = TRUE
    )
  }
  if (nrow(ranef)) {
    if (!length(csef)) {
      # only group-level category specific effects present
      str_add(out$modelD) <- glue(
        "  // linear predictor for category specific effects\n",               
        "  matrix[N, ncat{resp} - 1] mucs{p}", 
        " = rep_matrix(0, N, ncat{resp} - 1);\n"
      )
    }
    cats_regex <- "(?<=\\[)[[:digit:]]+(?=\\]$)"
    cats <- get_matches(cats_regex, ranef$coef, perl = TRUE)
    ncatM1 <- max(as.numeric(cats))
    for (i in seq_len(ncatM1)) {
      r_cat <- ranef[grepl(glue("\\[{i}\\]$"), ranef$coef), ]
      str_add(out$modelC2) <- glue(
        "    mucs{p}[n, {i}] = mucs{p}[n, {i}]"
      )
      for (id in unique(r_cat$id)) {
        r <- r_cat[r_cat$id == id, ]
        rpx <- check_prefix(r)
        idp <- paste0(r$id, usc(combine_prefix(rpx)))
        str_add(out$modelC2) <- cglue(
          " + r_{idp}_{r$cn}[J_{r$id}[n]] * Z_{idp}_{r$cn}[n]"
        )
      }
      str_add(out$modelC2) <- ";\n"
    }
  }
  out
}

stan_sp <- function(bterms, data, prior, stanvars, ranef, meef, ...) {
  # Stan code for special effects
  out <- list()
  spef <- tidy_spef(bterms, data)
  if (!nrow(spef)) return(out)
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
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
      cglue("  vector[N] Csp{p}_{seq_len(ncovars)};\n")
    )
  }
  prefix <- combine_prefix(px, keep_mu = TRUE)
  special <- attr(prior, "special")[[prefix]]
  define_bsp_in_pars <- is.null(special[["hs_df"]]) &&
    !glue("bsp{p}") %in% names(stanvars)
  if (define_bsp_in_pars) {
    bound <- get_bound(prior, class = "b", px = px)
    str_add(out$par) <- glue(
      "  // special effects coefficients\n", 
      "  vector{bound}[Ksp{p}] bsp{p};\n"
    )
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
      "  int<lower=2> Jmo{p}[Imo{p}];  // length of simplexes\n",
      "  // monotonic variables\n",
      cglue("  int Xmo{p}_{I}[N];\n"),
      "  // prior concentration of monotonic simplexes\n",
      cglue("  vector[Jmo{p}[{I}]] con_simo{p}_{I};\n")
    )
    str_add(out$par) <- glue(
      "  // simplexes of monotonic effects\n",
      cglue("  simplex[Jmo{p}[{I}]] simo{p}_{I};\n")
    ) 
    str_add(out$prior) <- cglue(
      "  target += dirichlet_lpdf(",
      "simo{p}_{I} | con_simo{p}_{I});\n"
    )
  }
  out <- collapse_lists(out,
    stan_special_prior_local(
      "bsp", prior, ncoef = nrow(spef), 
      px = px, center_X = FALSE
    )                      
  )
  out
}

stan_gp <- function(bterms, data, prior, ...) {
  # Stan code for latent gaussian processes
  out <- list()
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
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
      "  int<lower=1> Kgp{pi};\n",
      "  int<lower=1> Dgp{pi};\n"
    )
    if (!isNA(k)) {
      # !isNA(k) indicates the use of approximate GPs
      str_add(out$data) <- glue(
        "  int<lower=1> NBgp{pi};\n"
      )
    } 
    str_add(out$par) <- glue(
      "  // GP hyperparameters\n", 
      "  vector<lower=0>[Kgp{pi}] sdgp{pi};\n"
    )
    if (gpef$iso[i]) {
      str_add(out$par) <- glue(
        "  vector<lower=0>[1] lscale{pi}[Kgp{pi}];\n" 
      )
    } else {
      str_add(out$par) <- glue(
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
        "  int<lower=1> {Ngp}[Kgp{pi}];\n"
      )
      str_add(out$data) <- cglue(
        "  int<lower=1> {Igp} [{Ngp}[{J}]];\n",
        "  vector[{Ngp}[{J}]] Cgp{pi}_{J};\n"
      )
      if (gr) {
        str_add(out$data) <- glue(
          "  int<lower=1> Nsubgp{pi}[Kgp{pi}];\n"
        )
        str_add(out$data) <- cglue(
          "  int<lower=1> Jgp{pi}_{J}[{Ngp}[{J}]];\n"
        )
      }
      gp_call <- glue("Cgp{pi}_{J} .* ")
      if (!isNA(k)) {
        str_add(out$data) <- cglue(
          "  matrix[{Nsubgp}[{J}], NBgp{pi}] Xgp{pi}_{J};\n",
          "  vector[Dgp{pi}] slambda{pi}_{J}[NBgp{pi}];\n"
        )
        str_add(out$par) <- cglue(
          "  vector[NBgp{pi}] zgp{pi}_{J};\n"
        )
        str_add(gp_call) <- glue(
          "gpa(Xgp{pi}_{J}, sdgp{pi}[{J}], ", 
          "lscale{pi}[{J}], zgp{pi}_{J}, slambda{pi}_{J})"
        )
      } else {
        str_add(out$data) <- cglue(
          "  vector[Dgp{pi}] Xgp{pi}_{J}[{Nsubgp}[{J}]];\n"
        )
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
      str_add(out$modelCgp1) <- glue(
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
      Nsubgp <- str_if(gr, glue("Nsubgp{pi}"), "N")
      if (gr) {
        str_add(out$data) <- glue(
          "  int<lower=1> {Nsubgp};\n",
          "  int<lower=1> Jgp{pi}[N];\n"
        )
      }
      if (!isNA(k)) {
        str_add(out$data) <- glue(
          "  matrix[{Nsubgp}, NBgp{pi}] Xgp{pi};\n",
          "  vector[Dgp{pi}] slambda{pi}[NBgp{pi}];\n"
        )
        str_add(out$par) <- glue(
          "  vector[NBgp{pi}] zgp{pi};\n"
        )
        gp_call <- glue(
          "gpa(Xgp{pi}, sdgp{pi}[1], lscale{pi}[1], zgp{pi}, slambda{pi})"
        )
      } else {
        str_add(out$data) <- glue(
          "  vector[Dgp{pi}] Xgp{pi}[{Nsubgp}];\n"
        ) 
        str_add(out$par) <- glue(
          "  vector[{Nsubgp}] zgp{pi};\n"
        )
        gp_call <- glue(
          "gp(Xgp{pi}, sdgp{pi}[1], lscale{pi}[1], zgp{pi})"
        )
      }
      if (bynum) {
        str_add(out$data) <- glue(
          "  vector[N] Cgp{pi};\n"
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

stan_ac <- function(bterms, ...) {
  # Stan code for the linear predictor of certain autocorrelation terms 
  out <- list()
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  autocor <- bterms$autocor
  if (get_ar(autocor) && !use_cov(autocor)) {
    eta <- combine_prefix(px, keep_mu = TRUE)
    eta_ar <- glue("head(E{p}[n], Kar{p}) * ar{p}")
    str_add(out$modelC3) <- glue("    {eta}[n] += {eta_ar};\n")
  }
  if (get_ma(autocor) && !use_cov(autocor)) {
    str_add(out$loopeta) <- glue(" + head(E{p}[n], Kma{p}) * ma{p}")
  }
  if (get_arr(autocor)) {
    str_add(out$eta) <- glue(" + Yarr{p} * arr{p}")
  }
  if (is.cor_car(autocor)) {
    str_add(out$loopeta) <- glue(" + rcar{p}[Jloc{p}[n]]")
  }
  if (is.cor_bsts(autocor)) {
    str_add(out$loopeta) <- glue(" + loclev{p}[n]")
  }
  out
}

stan_offset <- function(bterms, ...) {
  # stan code for offsets
  out <- list()
  if (is.formula(bterms$offset)) {
    p <- usc(combine_prefix(bterms))
    str_add(out$data) <- glue( "  vector[N] offset{p};\n")
    str_add(out$eta) <- glue(" + offset{p}")
  }
  out
}

stan_eta_fe <- function(fixef, center_X = TRUE, sparse = FALSE, px = list()) {
  # define Stan code to compute the fixef part of eta
  # Args:
  #   fixef: names of the population-level effects
  #   center_X: use the centered design matrix?
  #   sparse: use sparse matrix multiplication?
  p <- usc(combine_prefix(px))
  if (length(fixef)) {
    if (sparse) {
      stopifnot(!center_X)
      csr_args <- sargs(
        paste0(c("rows", "cols"), "(X", p, ")"),
        paste0(c("wX", "vX", "uX", "b"), p)
      )
      eta_fe <- glue("csr_matrix_times_vector({csr_args})")
    } else {
      eta_fe <- glue("X", str_if(center_X, "c"), "{p} * b{p}")
    }
  } else { 
    eta_fe <- "rep_vector(0, N)"
  }
  glue(" + {eta_fe}")
}

stan_eta_re <- function(ranef, px = list()) {
  # write the group-level part of the linear predictor
  # Args:
  #   ranef: a named list returned by tidy_ranef
  #   nlpar: optional name of a non-linear parameter
  eta_re <- ""
  ranef <- subset2(ranef, type = c("", "mmc"), ls = px)
  for (id in unique(ranef$id)) {
    r <- subset2(ranef, id = id)
    rpx <- check_prefix(r)
    idp <- paste0(r$id, usc(combine_prefix(rpx)))
    if (r$gtype[1] == "mm") {
      ng <- seq_along(r$gcall[[1]]$groups)
      for (i in seq_rows(r)) {
        str_add(eta_re) <- cglue(
          " + W_{r$id[i]}_{ng}[n]", 
          " * r_{idp[i]}_{r$cn[i]}[J_{r$id[i]}_{ng}[n]]",
          " * Z_{idp[i]}_{r$cn[i]}_{ng}[n]"
        ) 
      }
    } else {
      str_add(eta_re) <- cglue(
        " + r_{idp}_{r$cn}[J_{r$id}[n]] * Z_{idp}_{r$cn}[n]"
      )
    }
  }
  eta_re
}

stan_eta_rsp <- function(r) {
  # Stan code for r parameters in special predictor terms
  # Args:
  #   r: data.frame created by tidy_ranef
  # Returns:
  #   A character vector, one element per row of 'r' 
  stopifnot(nrow(r) > 0L, length(unique(r$gtype)) == 1L)
  rpx <- check_prefix(r)
  idp <- paste0(r$id, usc(combine_prefix(rpx)))
  if (r$gtype[1] == "mm") {
    ng <- seq_along(r$gcall[[1]]$groups)
    out <- rep("", nrow(r))
    for (i in seq_along(out)) {
      out[i] <- glue(
        "W_{r$id[i]}_{ng}[n] * r_{idp[i]}_{r$cn[i]}[J_{r$id[i]}_{ng}[n]]",
        collapse = " + "
      ) 
    }
  } else {
    out <- glue("r_{idp}_{r$cn}[J_{r$id}[n]]")
  }
  out
}

stan_eta_transform <- function(family, llh_adj = FALSE) {
  # indicate whether eta needs to be transformed
  # manually using the link functions
  # Args:
  #   llh_adj: is the model censored or truncated?
  transeta <- "transeta" %in% family_info(family, "specials")
  !(family$link == "identity" && !transeta || has_cat(family)) &&
    (llh_adj || !stan_has_built_in_fun(family))
}

stan_eta_ilink <- function(dpar, bterms, resp = "") {
  # correctly apply inverse link to eta
  # Args:
  #   dpar: name of the parameter for which to define the link
  #   bterms: object of class brmsterms
  #   resp: name of the response variable
  stopifnot(is.brmsterms(bterms))
  out <- rep("", 2)
  family <- bterms$dpars[[dpar]]$family
  llh_adj <- stan_llh_adj(bterms$adforms, c("cens", "trunc"))
  if (stan_eta_transform(family, llh_adj = llh_adj)) {
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

stan_dpar_defs <- function(dpar, suffix = "", family = NULL, fixed = FALSE) {
  # default Stan definitions for distributional parameters
  # Args:
  #   dpar: name of a distributional parameter
  #   suffix: optional suffix of the parameter name
  #   family: optional brmsfamily object
  #   fixed: should the parameter be fixed to a certain value?
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

stan_dpar_defs_temp <- function(dpar, suffix = "", family = NULL) {
  # default Stan definitions for temporary distributional parameters
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
      "  real temp_", 
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

stan_dpar_transform <- function(bterms) {
  # Stan code for transformations of distributional parameters
  stopifnot(is.brmsterms(bterms))
  out <- list()
  families <- family_names(bterms)
  p <- usc(combine_prefix(bterms))
  if (any(conv_cats_dpars(families))) {
    str_add(out$modelD) <- glue( 
      "  // linear predictor matrix\n",
      "  vector[ncat{p}] mu{p}[N];\n"
    )
    mu_dpars <- str_subset(names(bterms$dpars), "^mu")
    mu_vector <- stan_vector(c("0", glue("{mu_dpars}{p}[n]")))
    str_add(out$modelC4) <- glue(
      "    mu{p}[n] = {mu_vector};\n"
    )
  }
  if (any(families %in% "skew_normal")) {
    # as suggested by Stephen Martin use sigma and mu of CP 
    # but the skewness parameter alpha of DP
    ap_names <- names(bterms$dpars)
    for (i in which(families %in% "skew_normal")) {
      id <- str_if(length(families) == 1L, "", i)
      sigma <- stan_sigma_transform(bterms, id = id)
      ns <- str_if(grepl("\\[n\\]", sigma), "[n]")
      na <- str_if(glue("alpha{id}") %in% ap_names, "[n]")
      type_delta <- str_if(nzchar(na), "vector[N]", "real")
      no <- str_if(any(nzchar(c(ns, na))), "[n]", "")
      type_omega <- str_if(nzchar(no), "vector[N]", "real")
      str_add(out$modelD) <- glue(
        "  {type_delta} delta{id}{p};\n",
        "  {type_omega} omega{id}{p};\n"
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
      str_add(out$modelC5) <- glue(
        str_if(!nzchar(na), comp_delta),
        str_if(!nzchar(no), comp_omega),
        "  for (n in 1:N) {{\n",
        str_if(nzchar(na), glue("  ", comp_delta)),
        str_if(nzchar(no), glue("  ", comp_omega)),
        "    mu{id}{p}[n] = mu{id}{p}[n]", 
        " - {omega} * {delta} * sqrt_2_div_pi;\n",
        "  }}\n"
      )
    }
  }
  if (any(families %in% "gen_extreme_value")) {
    ap_names <- c(names(bterms$dpars), names(bterms$fdpars))
    for (i in which(families %in% "gen_extreme_value")) {
      id <- str_if(length(families) == 1L, "", i)
      xi <- glue("xi{id}")
      if (!xi %in% ap_names) {
        str_add(out$modelD) <- glue(
          "  real {xi};  // scaled shape parameter\n"
        )
        sigma <- glue("sigma{id}")
        v <- str_if(sigma %in% names(bterms$dpars), "_vector")
        args <- sargs(
          glue("temp_{xi}"), glue("Y{p}"), 
          glue("mu{id}{p}"), glue("{sigma}{p}")
        )
        str_add(out$modelC5) <- glue(
          "  {xi}{p} = scale_xi{v}({args});\n"
        )
      }
    }
  }
  out
}

stan_sigma_transform <- function(bterms, id = "") {
  # Stan code for sigma to incorporate addition argument 'se'
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
