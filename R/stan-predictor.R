stan_effects.btl <- function(x, data, ranef, prior, center_X = TRUE, 
                             sparse = FALSE, ilink = rep("", 2), 
                             order_mixture = 'none',  ...) {
  # combine effects for the predictors of a single (non-linear) parameter
  # Args:
  #   center_X: center population-level design matrix if possible?
  #   eta: prefix of the linear predictor variable
  stopifnot(length(ilink) == 2L)
  px <- check_prefix(x)
  eta <- combine_prefix(px, keep_mu = TRUE)
  stopifnot(nzchar(eta))
  ranef <- subset2(ranef, ls = px)
  
  out <- list()
  # include population-level effects
  center_X <- center_X && has_intercept(x$fe) && 
    !is(x$autocor, "cor_bsts") && !sparse
  rm_int <- center_X || is(x$autocor, "cor_bsts") || is_ordinal(x$family)
  cols2remove <- if (rm_int) "Intercept"
  fixef <- setdiff(colnames(data_fe(x, data)$X), cols2remove)
  text_fe <- stan_fe(
    fixef, center_X = center_X, family = x$family, 
    prior = prior, px = px, sparse = sparse, 
    order_mixture = order_mixture
  )
  # include smooth terms
  smooths <- get_sm_labels(x, data = data)
  text_sm <- stan_sm(smooths, prior = prior, px = px)
  # include category specific effects
  csef <- colnames(get_model_matrix(x$cs, data))
  text_cs <- stan_cs(csef, ranef, prior = prior)
  # include monotonic effects
  monef <- all_terms(x$mo)
  text_mo <- stan_mo(monef, ranef, prior = prior, px = px)
  # include measurement error variables
  meef <- get_me_labels(x, data = data)
  text_me <- stan_me(meef, ranef, prior = prior, px = px)
  # include gaussian processes
  gpef <- get_gp_labels(x, data = data)
  text_gp <- stan_gp(gpef, prior = prior, px = px)
  
  out <- collapse_lists(
    out, text_fe, text_cs, text_mo, text_me, text_sm, text_gp
  )
  p <- usc(combine_prefix(px))
  if (is.formula(x$offset)) {
    str_add(out$data) <- paste0( 
      "  vector[N] offset", p, "; \n"
    )
  }
  # initialize and compute eta_<nlpar>
  str_add(out$modelD) <- paste0(
    "  vector[N] ", eta, " = ", 
    text_fe$eta, text_sm$eta, text_gp$eta,
    if (center_X && !is_ordinal(x$family))
      paste0(" + temp", p, "_Intercept"),
    if (is.formula(x$offset))
      paste0(" + offset", p),
    if (get_arr(x$autocor))
      " + Yarr * arr",
    "; \n"
  )
  
  # repare loop over eta
  eta_loop <- paste0(
    stan_eta_re(ranef, px = px),
    text_mo$eta, text_me$eta,
    stan_eta_ma(x$autocor, px = px), 
    stan_eta_car(x$autocor),
    stan_eta_bsts(x$autocor)
  )
  if (nzchar(eta_loop)) {
    str_add(out$modelC2) <- paste0(
      "    ", eta, "[n] = ", eta, "[n]", eta_loop, "; \n"
    )
  }
  # include autoregressive effects
  if (get_ar(x$autocor) && !use_cov(x$autocor)) {
    eta_ar <- paste0(eta, "[n] + head(E", p, "[n], Kar) * ar")
    str_add(out$modelC3) <- paste0(
      "    ", eta, "[n] = ", eta_ar, "; \n"
    )
  }
  # possibly transform eta before it is passed to the likelihood
  if (sum(nzchar(ilink))) {
    # make sure mu comes last as it might depend on other parameters
    not_mu <- nzchar(x$dpar) && dpar_class(x$dpar) != "mu"
    position <- ifelse(not_mu, "modelC3", "modelC4")
    str_add(out[[position]]) <- paste0(
      "    ", eta, "[n] = ", ilink[1], eta, "[n]", ilink[2], "; \n"
    )
  }
  out
}

stan_effects.btnl <- function(x, data, ranef, prior,
                              ilink = rep("", 2), ...) {
  # prepare Stan code for non-linear models
  # Args:
  #   data: data.frame supplied by the user
  #   ranef: data.frame returned by tidy_ranef
  #   prior: a brmsprior object
  #   nlpar: currently unused but should not be part of ...
  #   ilink: character vector of length 2 containing
  #          Stan code for the link function
  #   ...: passed to stan_effects.btl
  stopifnot(length(ilink) == 2L)
  out <- list()
  if (!length(x$nlpars)) {
    return(out)
  }
  nlpars <- names(x$nlpars)
  for (nlp in nlpars) {
    nl_text <- stan_effects(
      x = x$nlpars[[nlp]], data = data, 
      ranef = ranef, prior = prior, 
      center_X = FALSE, ...
    )
    out <- collapse_lists(out, nl_text)
  }
  # prepare non-linear model
  par <- combine_prefix(x, keep_mu = TRUE)
  new_nlpars <- paste0(" ", par, "_", nlpars, "[n] ")
  # covariates in the non-linear model
  covars <- wsp(setdiff(all.vars(rhs(x$formula)), nlpars))
  if (length(covars)) {
    # use vectors as indexing matrices in Stan is slow
    p <- usc(combine_prefix(x), "suffix")
    str_add(out$data) <- paste0( 
      "  // covariate vectors \n",
      collapse("  vector[N] C_", p, seq_along(covars), ";\n")
    )
    new_covars <- paste0(" C_", p, seq_along(covars), "[n] ")
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
  str_add(out$modelD) <- paste0("  vector[N] ", par, "; \n")
  str_add(out$modelC4) <- paste0(
    "    // compute non-linear predictor \n",
    "    ", par, "[n] = ", ilink[1], trimws(nlmodel), ilink[2], "; \n"
  )
  out
}

stan_effects.brmsterms <- function(x, data, ranef, prior, 
                                   sparse = FALSE, ...) {
  # Stan code for auxiliary parameters
  # Args:
  #   bterms: object of class brmsterms
  #   other arguments: same as make_stancode
  out <- list()
  valid_dpars <- valid_dpars(x$family, bterms = x)
  args <- nlist(data, ranef, prior)
  for (dp in valid_dpars) {
    ap_terms <- x$dpars[[dp]]
    if (is.btl(ap_terms) || is.btnl(ap_terms)) {
      ilink <- stan_eta_ilink(
        ap_terms$family, dpars = names(x$dpars), 
        adforms = x$adforms, mix = dpar_id(dp)
      )
      eta <- ifelse(dp == "mu", "mu", "")
      ap_args <- list(
        ap_terms, eta = eta, ilink = ilink,
        sparse = sparse, order_mixture = x$family$order
      )
      out[[dp]] <- do.call(stan_effects, c(ap_args, args))
    } else if (is.numeric(x$fdpars[[dp]]$value)) {
      out[[dp]] <- list(data = stan_dpar_defs(dp))
    } else if (is.character(x$fdpars[[dp]]$value)) {
      if (!x$fdpars[[dp]]$value %in% valid_dpars) {
        stop2("Parameter '", x$fdpars[[dp]]$value, "' cannot be found.")
      }
      out[[dp]] <- list(
        tparD = stan_dpar_defs(dp),
        tparC1 = paste0("  ", dp, " = ", x$fdpars[[dp]]$value, "; \n")
      )
    } else {
      def_temp <- stan_dpar_defs_temp(dp)
      def <- stan_dpar_defs(dp)
      if (nzchar(def_temp)) {
        out[[dp]] <- list(par = def_temp,
          prior = stan_prior(prior, class = dp, prefix = "temp_")
        )
      } else if (nzchar(def)) {
        out[[dp]] <- list(par = def,
          prior = stan_prior(prior, class = dp)
        )
      }
    }
  }
  collapse_lists(ls = out)
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
      x = bterms$dpars[["mu"]], data, ranef, prior, sparse,
      ilink = stan_eta_ilink(bterms$family, adforms = bterms$adforms)
    )
    resp <- bterms$response
    tmp_list <- named_list(resp)
    for (r in resp) {
      args$x$resp <- r
      tmp_list[[r]] <- do.call(stan_effects, args)
    }
    out <- collapse_lists(ls = tmp_list)
    if (is_linear(bterms$family)) {
      len_Eta_n <- "nresp" 
    } else if (is_categorical(bterms$family)) {
      len_Eta_n <- "ncat - 1"
    } else {
      stop2("Multivariate models are not yet implemented ", 
            "for family '", bterms$family$family, "'.")
    }
    str_add(out$modelD) <- paste0( 
      "  // multivariate linear predictor matrix \n",
      "  vector[", len_Eta_n, "] Mu[N]; \n"
    )
    str_add(out$modelC3) <- paste0(
      collapse("    Mu[n, ", seq_along(resp), "] = mu_", resp, "[n]; \n")
    )
  }
  out
}

stan_fe <- function(fixef, prior, family = gaussian(),
                    center_X = TRUE, px = list(),
                    sparse = FALSE, order_mixture = 'none') {
  # Stan code for population-level effects
  # Args:
  #   fixef: names of the population-level effects
  #   center_X: center the design matrix?
  #   family: the model family
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior 
  #   order_mixture: order intercepts to identify mixture models?
  # Returns:
  #   a list containing Stan code related to population-level effects
  p <- usc(combine_prefix(px))
  ct <- ifelse(center_X, "c", "")
  out <- list()
  if (length(fixef)) {
    str_add(out$data) <- paste0( 
      "  int<lower=1> K", p, ";",
      "  // number of population-level effects \n", 
      "  matrix[N, K", p, "] X", p, ";",
      "  // population-level design matrix \n"
    )
    if (sparse) {
      stopifnot(!center_X)
      str_add(out$tdataD) <- paste0(
        "  // sparse matrix representation of X", p, "\n",
        "  vector[rows(csr_extract_w(X", p, "))] wX", p, 
        " = csr_extract_w(X", p, ");\n",
        "  int vX", p, "[size(csr_extract_v(X", p, "))]",
        " = csr_extract_v(X", p, ");\n",
        "  int uX", p, "[size(csr_extract_u(X", p, "))]",
        " = csr_extract_u(X", p, ");\n"
      )
    }
    # prepare population-level coefficients
    prefix <- combine_prefix(px, keep_mu = TRUE)
    special <- attr(prior, "special")[[prefix]]
    if (!is.null(special[["hs_df"]])) {
      str_add(out$data) <- paste0(
        "  real<lower=0> hs_df", p, "; \n",
        "  real<lower=0> hs_df_global", p, "; \n",
        "  real<lower=0> hs_df_slab", p, "; \n",
        "  real<lower=0> hs_scale_global", p, "; \n",
        "  real<lower=0> hs_scale_slab", p, "; \n"           
      )
      str_add(out$par) <- paste0(
        "  // horseshoe shrinkage parameters \n",
        "  vector[K", ct, p, "] zb", p, "; \n",
        "  vector<lower=0>[K", ct, p, "] hs_local", p, "[2]; \n",
        "  real<lower=0> hs_global", p, "[2]; \n",
        "  real<lower=0> hs_c2", p, "; \n"
      )
      hs_scale_global <- paste0("hs_scale_global", p)
      if (isTRUE(special[["hs_autoscale"]])) {
        hs_scale_global <- paste0(hs_scale_global, " * sigma")
      }
      hs_args <- sargs(
        paste0(c("zb", "hs_local", "hs_global"), p), 
        hs_scale_global, 
        paste0("hs_scale_slab", p, "^2 * hs_c2", p)
      )
      str_add(out$tparD) <- paste0(
        "  // population-level effects \n",
        "  vector[K", ct, p, "] b", p,
        " = horseshoe(", hs_args, "); \n"
      )
    } else {
      bound <- get_bound(prior, class = "b", px = px)
      str_add(out$par) <- paste0(
        "  vector", bound, "[K", ct, p, "] b", p, ";",
        "  // population-level effects \n"
      )
    }
    if (!is.null(special[["lasso_df"]])) {
      str_add(out$data) <- paste0(
        "  real<lower=0> lasso_df", p, "; \n",
        "  real<lower=0> lasso_scale", p, "; \n"
      )
      str_add(out$par) <- paste0(
        "  // lasso shrinkage parameter \n",
        "  real<lower=0> lasso_inv_lambda", p, "; \n"
      )
    }
    str_add(out$prior) <- stan_prior(
      prior, class = "b", coef = fixef, px = px, suffix = p
    )
  }
  if (center_X) {
    if (length(fixef)) {
      str_add(out$tdataD) <- paste0(
        "  int Kc", p, " = K", p, " - 1; \n",
        "  matrix[N, K", p, " - 1] Xc", p, ";", 
        "  // centered version of X", p, " \n",
        "  vector[K", p, " - 1] means_X", p, ";",
        "  // column means of X", p, " before centering \n"
      )
      str_add(out$tdataC) <- paste0(
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
      str_add(out$genD) <- paste0(
        "  // compute actual thresholds \n",
        "  vector[ncat - 1] b_Intercept",  
        " = temp_Intercept", sub_X_means, "; \n" 
      )
    } else {
       if (identical(dpar_class(px$dpar), order_mixture)) {
         # identify mixtures via ordering of the intercepts
         ap_id <- dpar_id(px$dpar)
         str_add(out$tparD) <- paste0(
           "  // identify mixtures via ordering of the intercepts \n",                   
           "  real temp", p, "_Intercept",
           " = ordered_Intercept[", ap_id, "]; \n"
         )
      } else {
        str_add(out$par) <- paste0(
          "  real temp", p, "_Intercept;  // temporary intercept \n"
        )
      }
      str_add(out$genD) <- paste0(
        "  // actual population-level intercept \n",
        "  real b", p, "_Intercept",
        " = temp", p, "_Intercept", sub_X_means, "; \n"
      )
    }
    # for equidistant thresholds only temp_Intercept1 is a parameter
    prefix <- paste0("temp", p, "_")
    suffix <- ifelse(is_equal(family$threshold, "equidistant"), "1", "")
    str_add(out$prior) <- stan_prior(
      prior, class = "Intercept", px = px,
      prefix = prefix, suffix = suffix
    )
  } else {
    if (identical(dpar_class(px$dpar), order_mixture)) {
      stop2("Identifying mixture components via ordering requires ",
            "population-level intercepts to be present.\n",
            "Try setting order = 'none' in function 'mixture'.")
    }
  }
  out$eta <- stan_eta_fe(fixef, center_X, sparse, px = px)
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
  r <- subset2(ranef, id = id)
  ccov <- r$group[1] %in% names(cov_ranef)
  ng <- seq_along(r$gcall[[1]]$groups)
  px <- check_prefix(r)
  idp <- paste0(r$id, usc(combine_prefix(px)))
  out <- list()
  str_add(out$data) <- paste0(
    "  // data for group-level effects of ID ", id, " \n",
    if (r$gtype[1] == "mm") {
      collapse(
        "  int<lower=1> J_", id, "_", ng, "[N]; \n",
        "  real W_", id, "_", ng, "[N]; \n"
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
  str_add(out$prior) <- stan_prior(
    prior, class = "sd", group = r$group[1], coef = r$coef,
    px = px, suffix = paste0("_", id)
  )
  J <- seq_len(nrow(r))
  has_def_type <- !r$type %in% c("mo", "me")
  if (any(has_def_type)) {
    str_add(out$data) <- collapse(
        "  vector[N] Z_", idp[has_def_type], 
        "_", r$cn[has_def_type], "; \n"
    ) 
  }
  str_add(out$par) <- paste0(
    "  vector<lower=0>[M_", id, "] sd_", id, ";",
    "  // group-level standard deviations \n"
  )
  if (nrow(r) > 1L && r$cor[1]) {
    # multiple correlated group-level effects
    str_add(out$data) <- paste0( 
      "  int<lower=1> NC_", id, "; \n"
    )
    str_add(out$par) <- paste0(
      "  matrix[M_", id, ", N_", id, "] z_", id, ";",
      "  // unscaled group-level effects \n",    
      "  // cholesky factor of correlation matrix \n",
      "  cholesky_factor_corr[M_", id, "] L_", id, "; \n"
    )
    str_add(out$prior) <- paste0( 
      stan_prior(prior, class = "L", group = r$group[1],
                 suffix = paste0("_", id)),
      "  target += normal_lpdf(to_vector(z_", id, ") | 0, 1); \n"
    )
    str_add(out$tparD) <- paste0(
      "  // group-level effects \n",
      "  matrix[N_", id, ", M_", id, "] r_", id, 
      if (ccov) {
        # customized covariance matrix supplied
        paste0(
          " = as_matrix(kronecker(Lcov_", id, ",", 
          " diag_pre_multiply(sd_", id,", L_", id,")) *",
          " to_vector(z_", id, "), N_", id, ", M_", id, "); \n"
        )
      } else {
        paste0(
          " = (diag_pre_multiply(sd_", id, ", L_", id,") * z_", id, ")'; \n"
        )
      },
      collapse(
        "  vector[N_", id, "] r_", idp, "_", r$cn, 
        " = r_", id, "[, ", J, "]; \n"
      )
    )
    # return correlations above the diagonal only
    cors_genC <- ulapply(2:nrow(r), function(k) 
      lapply(1:(k - 1), function(j) paste0(
        "  cor_", id, "[", (k - 1) * (k - 2) / 2 + j, 
        "] = Cor_", id, "[", j, ",", k, "]; \n"
      ))
    )
    str_add(out$genD) <- paste0(
      "  corr_matrix[M_", id, "] Cor_", id, 
      " = multiply_lower_tri_self_transpose(L_", id, "); \n",
      "  vector<lower=-1,upper=1>[NC_", id, "] cor_", id, "; \n"
    )
    str_add(out$genC) <- paste0(
      "  // take only relevant parts of correlation matrix \n",
      collapse(cors_genC)
    ) 
  } else {
    # single or uncorrelated group-level effects
    str_add(out$par) <- paste0(
      "  vector[N_", id, "] z_", id, "[M_", id, "];",
      "  // unscaled group-level effects \n"
    )
    str_add(out$prior) <- collapse(
      "  target += normal_lpdf(z_", id, "[", 1:nrow(r), "] | 0, 1); \n"
    )
    str_add(out$tparD) <- paste0(
      "  // group-level effects \n", 
      collapse(
        "  vector[N_", id, "] r_", idp, "_", r$cn,
        " = sd_", id, "[", J, "] * (", 
        if (ccov) paste0("Lcov_", id, " * "), 
        "z_", id, "[", J, "]); \n"
      )
    )
  }
  out
}

stan_sm <- function(smooths, prior, px = list()) {
  # Stan code of smooth terms
  # Args:
  #   smooths: names of the smooth terms
  #   prior: object of class brmsprior
  #   nlpar: optional name of a non-linear parameter
  # Returns:
  #   A list of strings containing Stan code
  out <- list()
  p <- usc(combine_prefix(px))
  if (length(smooths)) {
    stopifnot(!is.null(attr(smooths, "nbases")))
    for (i in seq_along(smooths)) {
      pi <- paste0(p, "_", i)
      nb <- seq_len(attr(smooths, "nbases")[[i]])
      str_add(out$data) <- paste0(
        "  // data of smooth ", smooths[i], "\n",  
        "  int nb", pi, ";  // number of bases \n",
        "  int knots", pi, "[nb", pi, "]; \n"
      )
      str_add(out$data) <- collapse(
        "  matrix[N, knots", pi, "[", nb, "]]", 
        " Zs", pi, "_", nb, "; \n"
      )
      str_add(out$par) <- paste0(
        "  // parameters of smooth ", smooths[i], "\n",
        collapse(
          "  vector[knots", pi, "[", nb, "]] zs", pi,"_", nb, "; \n",
          "  real<lower=0> sds", pi, "_", nb, "; \n"
        )
      )
      str_add(out$tparD) <- collapse(
        "  vector[knots", pi, "[", nb, "]] s", pi, "_", nb, 
        " = sds", pi,  "_", nb, " * zs", pi, "_", nb, "; \n"
      )
      str_add(out$prior) <- paste0(
        collapse(
          "  target += normal_lpdf(zs", pi, "_", nb, " | 0, 1); \n"
        ),
        stan_prior(prior, class = "sds", coef = smooths[i], 
                   px = px, suffix = paste0(pi, "_", nb))
      )
    }
    out$eta <- stan_eta_sm(smooths, px = px)
  }
  out
}

stan_mo <- function(monef, ranef, prior, px = list()) {
  # Stan code for monotonic effects
  # Args:
  #   monef: names of the monotonic effects
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior
  p <- usc(combine_prefix(px))
  out <- list()
  if (length(monef)) {
    I <- seq_along(monef)
    str_add(out$data) <- paste0(
      "  int<lower=1> Kmo", p, ";  // number of monotonic effects \n",
      "  int Xmo", p, "[N, Kmo", p, "];  // monotonic design matrix \n",
      "  int<lower=2> Jmo", p, "[Kmo", p, "];  // length of simplexes \n",
      collapse(
        "  vector[Jmo", p, "[", I, "]]", 
        " con_simplex", p, "_", I, "; \n"
      )
    )
    # FIXME: pass px
    bound <- get_bound(prior, class = "b", px = px)
    str_add(out$par) <- paste0(
      "  // monotonic effects \n", 
      "  vector", bound, "[Kmo", p, "] bmo", p, "; \n",
      collapse(
        "  simplex[Jmo", p, "[", I, "]]", 
        " simplex", p, "_", I, "; \n"
      )
    ) 
    str_add(out$prior) <- paste0(
      stan_prior(prior, class = "b", coef = monef, 
                 px = px, suffix = paste0("mo", p)),
      collapse(
        "  target += dirichlet_lpdf(",
        "simplex", p, "_", I, " | con_simplex", p, "_", I, "); \n"
      )
    )
    str_add(out$eta) <- stan_eta_mo(
      monef, ranef = ranef, px = px
    )
  }
  out
}

stan_cs <- function(csef, ranef, prior, px = list()) {
  # Stan code for category specific effects
  # Args:
  #   csef: names of the category specific effects
  #   prior: a data.frame containing user defined priors 
  #          as returned by check_prior
  # (!) Not yet implemented for non-linear models
  ranef <- subset2(ranef, type = "cs", ls = px)
  out <- list()
  if (length(csef)) {
    str_add(out$data) <- paste0(
      "  int<lower=1> Kcs;  // number of category specific effects \n",
      "  matrix[N, Kcs] Xcs;  // category specific design matrix \n"
    )
    bound <- get_bound(prior, class = "b")
    str_add(out$par) <- paste0(
      "  matrix", bound, "[Kcs, ncat - 1] bcs;",
      "  // category specific effects \n"
    )
    str_add(out$modelD) <- paste0(
      "  // linear predictor for category specific effects \n",
      "  matrix[N, ncat - 1] mucs = Xcs * bcs; \n"
    ) 
    str_add(out$prior) <- stan_prior(
      prior, class = "b", coef = csef,
      suffix = "cs", matrix = TRUE
    )
  }
  if (nrow(ranef)) {
    if (!length(csef)) {
      # only group-level category specific effects present
      str_add(out$modelD) <- paste0(
        "  // linear predictor for category specific effects \n",               
        "  matrix[N, ncat - 1] mucs = rep_matrix(0, N, ncat - 1); \n"
      )
    }
    cats <- get_matches("\\[[[:digit:]]+\\]$", ranef$coef)
    ncatM1 <- max(as.numeric(substr(cats, 2, nchar(cats) - 1)))
    for (i in seq_len(ncatM1)) {
      r_cat <- ranef[grepl(paste0("\\[", i, "\\]$"), ranef$coef), ]
      str_add(out$modelC2) <- paste0(
        "    mucs[n, ", i, "] = mucs[n, ", i, "]"
      )
      for (id in unique(r_cat$id)) {
        r <- r_cat[r_cat$id == id, ]
        rpx <- check_prefix(r)
        idp <- paste0(r$id, usc(combine_prefix(rpx)))
        str_add(out$modelC2) <- collapse(
          " + r_", idp, "_", r$cn, "[J_", r$id, "[n]]",
          " * Z_", idp, "_", r$cn, "[n]"
        )
      }
      str_add(out$modelC2) <- "; \n"
    }
  }
  out
} 

stan_me <- function(meef, ranef, prior, px = list()) {
  # stan code for measurement error effects
  # Args:
  #   meef: vector of terms containing noisy predictors
  out <- list()
  if (length(meef)) {
    not_one <- attr(meef, "not_one")
    uni_me <- attr(meef, "uni_me")
    p <- usc(combine_prefix(px))
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
      me_sp[[i]][I] <- substr(me_sp[[i]][I], 3, nchar(me_sp[[i]][I]) - 1)
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
    ranef <- subset2(ranef, type = "me", ls = px)
    invalid_coef <- setdiff(ranef$coef, meef)
    if (length(invalid_coef)) {
      stop2("Noisy group-level terms require ", 
            "corresponding population-level terms.")
    }
    for (i in seq_along(meef)) {
      r <- subset2(ranef, coef = meef[i])
      if (nrow(r)) {
        rpars <- paste0(" + ", stan_eta_r(r))
      } else {
        rpars <- ""
      }
      str_add(out$eta) <- paste0(
        " + (bme", p, "[", i, "]", rpars, ") * ", 
        meef_terms[i], covars[i]
      )
    }
    
    # prepare Stan code
    str_add(out$data) <- paste0(
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
    str_add(out$par) <- paste0(
      "  // noise free variables \n",
      collapse("  vector[N] Xme", pK, "; \n"),  
      "  vector[Kme", p, "] bme", p, ";",
      "  // coefficients of noise-free terms \n"
    )
    str_add(out$prior) <- paste0(
      stan_prior(prior, class = "b", coef = meef, 
                 px = px, suffix = paste0("me", p)),
      collapse(
        "  target += normal_lpdf(Xme", pK, " | Xn", pK, ", noise", pK,"); \n"
      )
    )
  }
  out
}

stan_gp <- function(gpef, prior, px = list()) {
  # Stan code for latent gaussian processes
  # Args:
  #   gpef: names of the gaussian process terms
  p <- usc(combine_prefix(px))
  out <- list()
  for (i in seq_along(gpef)) {
    pi <- paste0(p, "_", i)
    byvar <- attr(gpef, "byvars")[[i]]
    by_levels <- attr(gpef, "by_levels")[[i]]
    byfac <- length(by_levels) > 0L
    bynum <- !byfac && !identical(byvar, "NA")
    J <- seq_along(by_levels)
    str_add(out$data) <- paste0(
      "  int<lower=1> Kgp", pi, "; \n",
      "  int<lower=1> Mgp", pi, "; \n",
      "  vector[Mgp", pi, "] Xgp", p, "_", i, "[N]; \n",
      if (bynum) {
        paste0("  vector[N] Cgp", p, "_", i, "; \n")
      },
      if (byfac) {
        paste0(
          "  int<lower=1> Igp", pi, "[Kgp", p, "_", i, "]; \n",
          collapse(
            "  int<lower=1> Jgp", pi, "_", J, "[Igp", pi, "[", J, "]]; \n"
          )
        )
      }
    )
    str_add(out$par) <- paste0(
      "  // GP hyperparameters \n", 
      "  vector<lower=0>[Kgp", pi, "] sdgp", pi, "; \n",
      "  vector<lower=0>[Kgp", pi, "] lscale", pi, "; \n",
      "  vector[N] zgp", pi, "; \n"
    ) 
    rgpef <- rename(gpef[i])
    str_add(out$prior) <- paste0(
      stan_prior(prior, class = "sdgp", coef = rgpef, 
                 px = px, suffix = pi),
      stan_prior(prior, class = "lscale", coef = rgpef, 
                 px = px, suffix = pi),
      collapse(tp(), "normal_lpdf(zgp", pi, " | 0, 1); \n")
    )
    if (byfac) {
      Jgp <- paste0("Jgp", pi, "_", J)
      eta <- paste0(combine_prefix(px, keep_mu = TRUE), "[", Jgp, "]")
      gp_args <- paste0(
        "Xgp", pi, "[", Jgp, "], sdgp", pi, "[", J, "], ", 
        "lscale", pi, "[", J, "], zgp", pi, "[", Jgp, "]"
      )
      str_add(out$modelCgp1) <- paste0(
        collapse("  ", eta, " = ", eta, " + gp(", gp_args, "); \n")
      )
    } else {
      gp_args <- paste0(
        "Xgp", pi, ", sdgp", pi, "[1], lscale", pi, "[1], zgp", pi
      )
      Cgp <- ifelse(bynum, paste0("Cgp", pi, " .* "), "")
      str_add(out$eta) <- paste0(" + ", Cgp, "gp(", gp_args, ")")   
    }
  }
  out
}

stan_eta_fe <- function(fixef, center_X = TRUE, sparse = FALSE, 
                        px = list()) {
  # define Stan code to compute the fixef part of eta
  # Args:
  #   fixef: names of the population-level effects
  #   center_X: use the centered design matrix?
  #   sparse: use sparse matrix multiplication?
  #   nlpar: optional name of a non-linear parameter
  p <- usc(combine_prefix(px))
  if (length(fixef)) {
    if (sparse) {
      stopifnot(!center_X)
      csr_args <- sargs(
        paste0(c("rows", "cols"), "(X", p, ")"),
        paste0(c("wX", "vX", "uX", "b"), p)
      )
      eta_fe <- paste0("csr_matrix_times_vector(", csr_args, ")")
    } else {
      eta_fe <- paste0("X", if (center_X) "c", p, " * b", p)
    }
  } else { 
    eta_fe <- "rep_vector(0, N)"
  }
  eta_fe
}

stan_eta_re <- function(ranef, px = list()) {
  # write the group-level part of the linear predictor
  # Args:
  #   ranef: a named list returned by tidy_ranef
  #   nlpar: optional name of a non-linear parameter
  eta_re <- ""
  ranef <- subset2(ranef, type = "", ls = px)
  for (id in unique(ranef$id)) {
    r <- subset2(ranef, id = id)
    rpx <- check_prefix(r)
    idp <- paste0(r$id, usc(combine_prefix(rpx)))
    str_add(eta_re) <- collapse(
      " + (", stan_eta_r(r), ") * Z_", idp, "_", r$cn, "[n]"
    )
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
  rpx <- check_prefix(r)
  idp <- paste0(r$id, usc(combine_prefix(rpx)))
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

stan_eta_mo <- function(monef, ranef, px = list()) {
  # write the linear predictor for monotonic effects
  # Args:
  #   monef: names of the monotonic effects
  #   nlpar: an optional character string to add to the varnames
  #         (used for non-linear models)
  p <- usc(combine_prefix(px))
  eta_mo <- ""
  ranef <- subset2(ranef, type = "mo", ls = px)
  invalid_coef <- setdiff(ranef$coef, monef)
  if (length(invalid_coef)) {
    stop2("Monotonic group-level terms require ", 
          "corresponding population-level terms.")
  }
  for (i in seq_along(monef)) {
    r <- subset2(ranef, coef = monef[i])
    if (nrow(r)) {
      rpars <- paste0(" + ", stan_eta_r(r))
    } else {
      rpars <- ""
    }
    str_add(eta_mo) <- paste0(
      " + (bmo", p, "[", i, "]", rpars, ") * mo(",
      "simplex", p, "_", i, ", Xmo", p, "[n, ", i, "])"
    )
  }
  eta_mo
}

stan_eta_sm <- function(smooths, px = list()) {
  # write the linear predictor for smooth terms
  # Args:
  #   smooths: names of the smooth terms
  #   nlpar: optional character string to add to the varnames
  p <- usc(combine_prefix(px))
  eta_smooths <- ""
  if (length(smooths)) {
    stopifnot(!is.null(attr(smooths, "nbases")))
    for (i in seq_along(smooths)) {
      pi <- paste0(p, "_", i)
      nb <- seq_len(attr(smooths, "nbases")[[smooths[i]]])
      str_add(eta_smooths) <- collapse(
        " + Zs", pi, "_", nb, " * s", pi, "_", nb
      )
    }
  }
  eta_smooths
}

stan_eta_ma <- function(autocor, px = list()) {
  # write the linear predictor for MA terms
  # Args:
  #   autocor: object of class cor_brms
  #   nlpar: optional character string to add to the varnames
  out <- ""
  if (get_ma(autocor) && !use_cov(autocor)) {
    px <- combine_prefix(px)
    str_add(out) <- paste0(" + head(E", px, "[n], Kma) * ma")
  }
  out
}

stan_eta_car <- function(autocor) {
  # write the linear predictor for CAR terms
  # Args:
  #   autocor: object of class cor_brms
  out <- ""
  if (is.cor_car(autocor)) {
    str_add(out) <- " + rcar[Jloc[n]]"
  }
  out
}

stan_eta_bsts <- function(autocor) {
  # write the linear predictor for BSTS terms
  # Args:
  #   autocor: object of class cor_brms
  out <- ""
  if (is.cor_bsts(autocor)) {
    str_add(out) <- " + loclev[n]"
  }
  out
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
    is_ordinal(family) || is_categorical(family)) &&
  (llh_adj || !stan_has_built_in_fun(family))
}

stan_eta_ilink <- function(family, dpars = NULL, 
                           adforms = NULL, mix = "") {
  # correctly apply inverse link to eta
  # Args:
  #   family: a list with elements 'family' and 'link
  #   dpars: names of auxiliary parameters
  #   adforms: list of formulas containing addition terms
  stopifnot(all(c("family", "link") %in% names(family)))
  llh_adj <- stan_llh_adj(adforms, c("cens", "trunc"))
  if (stan_eta_transform(family, llh_adj = llh_adj)) {
    link <- family$link
    family <- family$family
    shape <- paste0("shape", mix)
    shape <- ifelse(
      is.formula(adforms$disp), paste0("disp_", shape, "[n]"), 
      ifelse(shape %in% dpars, paste0(shape, "[n]"), shape)
    )
    nu <- paste0("nu", mix)
    nu <- ifelse(nu %in% dpars, paste0(nu, "[n]"), nu)
    fl <- ifelse(
      family %in% c("gamma", "hurdle_gamma", "exponential"), 
      paste0(family, "_", link), family
    )
    ilink <- stan_ilink(link)
    out <- switch(fl,
      c(paste0(ilink, "("), ")"),
      gamma_log = c(paste0(shape, " * exp(-("), "))"),
      gamma_inverse = c(paste0(shape, " * ("), ")"),
      gamma_identity = c(paste0(shape, " / ("), ")"),
      hurdle_gamma_log = c(paste0(shape, " * exp(-("), "))"),
      hurdle_gamma_inverse = c(paste0(shape, " * ("), ")"),
      hurdle_gamma_identity = c(paste0(shape, " / ("), ")"),
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

stan_dpar_defs <- function(dpar) {
  # default Stan definitions for auxiliary parameters
  default_defs <- list(
    sigma = c(
      "  real<lower=0> ", 
      ";  // residual SD \n"
    ),
    shape = c(
      "  real<lower=0> ", 
      ";  // shape parameter \n"
    ),
    nu = c(
      "  real<lower=1> ", 
      ";  // degrees of freedom or shape \n"
    ),
    phi = c(
      "  real<lower=0> ", 
      ";  // precision parameter \n"
    ),
    kappa = c(
      "  real<lower=0> ", 
      ";  // precision parameter \n"
    ),
    beta = c(
      "  real<lower=0> ", 
      ";  // scale parameter \n"
    ),
    zi = c(
      "  real<lower=0,upper=1> ", 
      ";  // zero-inflation probability \n"
    ), 
    hu = c(
      "  real<lower=0,upper=1> ", 
      ";  // hurdle probability \n"
    ),
    zoi = c(
      "  real<lower=0,upper=1> ", 
      ";  // zero-one-inflation probability \n"
    ), 
    coi = c(
      "  real<lower=0,upper=1> ", 
      ";  // conditional one-inflation probability \n"
    ),
    bs = c(
      "  real<lower=0> ", 
      ";  // boundary separation parameter \n"
    ),
    ndt = c(
      "  real<lower=0,upper=min_Y> ", 
      ";  // non-decision time parameter \n"
    ),
    bias = c(
      "  real<lower=0,upper=1> ", 
      ";  // initial bias parameter \n"
    ),
    disc = c(
      "  real<lower=0> ", 
      ";  // discrimination parameters \n"
    ),
    quantile = c(
      "  real<lower=0,upper=1> ", 
      ";  // quantile parameter \n"
    ),
    xi = c(
      "  real ", 
      ";  // shape parameter \n"
    ),
    alpha = c(
      "  real ",
      ";  // skewness parameter \n"
    )
    # theta is handled in stan_mixture
  )
  def <- default_defs[[dpar_class(dpar)]]
  if (!is.null(def)) {
    def <- paste0(def[1], dpar, def[2])
  } else {
    def <- ""
  }
  def
}

stan_dpar_defs_temp <- function(dpar) {
  # default Stan definitions for temporary auxiliary parameters
  default_defs <- list(
    xi = c(
      "  real temp_", 
      ";  // unscaled shape parameter \n"
    )
  )
  def <- default_defs[[dpar_class(dpar)]]
  if (!is.null(def)) {
    def <- paste0(def[1], dpar, def[2])
  } else {
    def <- ""
  }
  def
}
