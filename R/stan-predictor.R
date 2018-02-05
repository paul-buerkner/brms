stan_effects <- function(x, ...) {
  # generate stan code various kind of effects 
  UseMethod("stan_effects")
}

#' @export
stan_effects.btl <- function(x, data, ranef, prior, center_X = TRUE, 
                             sparse = FALSE, ilink = rep("", 2), 
                             order_mixture = 'none',  ...) {
  # combine effects for the predictors of a single (non-linear) parameter
  # Args:
  #   center_X: center population-level design matrix if possible?
  #   sparse: should the population-level design matrix be treated as sparse?
  #   ilink: character vector of lenght 2 defining the link to be applied
  #   order_mixture: indicates how to identify mixture models via ordering
  stopifnot(length(ilink) == 2L)
  px <- check_prefix(x)
  eta <- combine_prefix(px, keep_mu = TRUE)
  stopifnot(nzchar(eta))
  ranef <- subset2(ranef, ls = px)
  center_X <- center_X && has_intercept(x$fe) && 
    !is(x$autocor, "cor_bsts") && !sparse
  out <- collapse_lists(
    text_fe <- stan_fe(
      x, data, prior = prior, center_X = center_X,
      sparse = sparse, order_mixture = order_mixture
    ),
    text_sp <- stan_sp(x, data, ranef = ranef, prior = prior),
    text_cs <- stan_cs(x, data, ranef = ranef, prior = prior),
    text_sm <- stan_sm(x, data, prior = prior),
    text_gp <- stan_gp(x, data, prior = prior)
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
      paste0(" + Yarr", p, " * arr", p),
    "; \n"
  )
  
  # repare loop over eta
  eta_loop <- paste0(
    stan_eta_re(ranef, px = px), text_sp$eta,
    stan_eta_autocor(x$autocor, px = px)
  )
  if (nzchar(eta_loop)) {
    str_add(out$modelC2) <- paste0(
      "    ", eta, "[n] = ", eta, "[n]", eta_loop, "; \n"
    )
  }
  # include autoregressive effects
  if (get_ar(x$autocor) && !use_cov(x$autocor)) {
    eta_ar <- paste0(eta, "[n] + head(E", p, "[n], Kar", p, ") * ar", p)
    str_add(out$modelC3) <- paste0("    ", eta, "[n] = ", eta_ar, "; \n")
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

#' @export
stan_effects.btnl <- function(x, data, ranef, prior,
                              ilink = rep("", 2), ...) {
  # prepare Stan code for non-linear models
  # Args:
  #   data: data.frame supplied by the user
  #   ranef: data.frame returned by tidy_ranef
  #   prior: a brmsprior object
  #   nlpar: currently unused but should not be part of ...
  #   ilink: character vector of length 2 containing
  #     Stan code for the link function
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
  nlmodel <- rm_wsp(collapse(deparse(x$formula[[2]])))
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

#' @export
stan_effects.brmsterms <- function(x, data, ranef, prior, 
                                   sparse = FALSE, rescor = FALSE, 
                                   ...) {
  # Stan code for distributional parameters
  # Args:
  #   rescor: indicate if this is part of an MV model estimating rescor
  px <- check_prefix(x)
  resp <- usc(combine_prefix(px))
  out <- list(stan_response(x, data = data))
  valid_dpars <- valid_dpars(x)
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
      out[[dp]] <- list(data = stan_dpar_defs(dp, resp))
    } else if (is.character(x$fdpars[[dp]]$value)) {
      if (!x$fdpars[[dp]]$value %in% valid_dpars) {
        stop2("Parameter '", x$fdpars[[dp]]$value, "' cannot be found.")
      }
      out[[dp]] <- list(
        tparD = stan_dpar_defs(dp, resp),
        tparC1 = paste0(
          "  ", dp, resp, " = ", x$fdpars[[dp]]$value, resp, "; \n"
        )
      )
    } else {
      def_temp <- stan_dpar_defs_temp(dp, resp)
      def <- stan_dpar_defs(dp, resp)
      if (nzchar(def_temp)) {
        out[[dp]] <- list(par = def_temp,
          prior = stan_prior(
            prior, class = dp, prefix = "temp_", suffix = resp, px = px)
        )
      } else if (nzchar(def)) {
        out[[dp]] <- list(par = def,
          prior = stan_prior(prior, class = dp, suffix = resp, px = px)
        )
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
stan_effects.mvbrmsterms <- function(x, prior, ...) {
  out <- collapse_lists(
    ls = lapply(x$terms, stan_effects, prior = prior, rescor = x$rescor, ...)
  )
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
    str_add(out$modelD) <- paste0( 
      "  // multivariate linear predictor matrix \n",
      "  vector[nresp] Mu[N]; \n"
    )
    str_add(out$modelC4) <- paste0(
      "    Mu[n] = ", stan_vector(paste0("mu_", resp, "[n]")), ";\n"
    )
    str_add(out$data) <- paste0(
      "  int<lower=1> nresp;  // number of responses\n",   
      "  int nrescor;  // number of residual correlations\n"
    )
    str_add(out$tdataD) <- "  vector[nresp] Y[N];  // response matrix\n"
    str_add(out$tdataC) <- paste0(
      "  for (n in 1:N) {\n",
      "    Y[n] = ", stan_vector(paste0("Y_", resp, "[n]")), ";\n",
      "  }\n"
    )
    if (any(adnames %in% "weights")) {
      str_add(out$tdataD) <- paste0(
        "  vector<lower=0>[N] weights = weights_", resp[1], ";\n" 
      )
    }
    miforms <- rmNULL(lapply(adforms, "[[", "mi"))
    if (length(miforms)) {
      str_add(out$modelD) <- "  vector[nresp] Yf[N] = Y;\n"
      for (i in seq_along(miforms)) {
        j <- match(names(miforms)[i], resp)
        str_add(out$modelC2) <- paste0(
          "    Yf[n][", j, "] = Yf_", resp[j], "[n];\n"
        )
      }
    }
    str_add(out$par) <- paste0(
      "  // parameters for multivariate linear models \n",
      "  cholesky_factor_corr[nresp] Lrescor; \n"
    )
    str_add(out$prior) <- stan_prior(prior, class = "Lrescor")
    if (family == "student") {
      str_add(out$par) <- stan_dpar_defs("nu")
      str_add(out$prior) <- stan_prior(prior, class = "nu")
    } 
    sigma <- ulapply(x$terms, stan_sigma_transform)
    if (any(grepl("\\[n\\]", sigma))) {
      str_add(out$modelD) <- "  vector[nresp] sigma[N];\n"
      str_add(out$modelC4) <- paste0(
        "    sigma[n] = ", stan_vector(sigma), ";\n"
      )
      if (family == "gaussian") {
        str_add(out$modelD) <- paste0(
          "  // cholesky factor of residual covariance matrix \n",
          "  matrix[nresp, nresp] LSigma[N];\n"
        )
        str_add(out$modelC4) <- 
          "    LSigma[n] = diag_pre_multiply(sigma[n], Lrescor);\n"
      } else if (family == "student") {
        str_add(out$modelD) <- paste0(
          "  // residual covariance matrix \n",
          "  matrix[nresp, nresp] Sigma[N];\n"
        )
        str_add(out$modelC4) <- paste0(
          "    Sigma[n] = multiply_lower_tri_self_transpose(", 
          "diag_pre_multiply(sigma[n], Lrescor)); \n" 
        )
      }
    } else {
      str_add(out$modelD) <- paste0(
        "  vector[nresp] sigma = ", stan_vector(sigma), ";\n"
      )
      if (family == "gaussian") {
        str_add(out$modelD) <- paste0(
          "  // cholesky factor of residual covariance matrix \n",
          "  matrix[nresp, nresp] LSigma = ",
          "diag_pre_multiply(sigma, Lrescor); \n"
        )
      } else if (family == "student") {
        str_add(out$modelD) <- paste0(
          "  // residual covariance matrix \n",
          "  matrix[nresp, nresp] Sigma = ",
          "multiply_lower_tri_self_transpose(", 
          "diag_pre_multiply(sigma, Lrescor)); \n"
        )
      }
    }
    str_add(out$genD) <- paste0(
      "  matrix[nresp, nresp] Rescor",
      " = multiply_lower_tri_self_transpose(Lrescor); \n",
      "  vector<lower=-1,upper=1>[nrescor] rescor; \n"
    )
    rescor_genC <- ulapply(2:nresp, function(i) 
      lapply(1:(i - 1), function(j) paste0(
        "  rescor[", (i - 1) * (i - 2) / 2 + j, 
        "] = Rescor[", j, ", ", i, "]; \n"
      ))
    )
    str_add(out$genC) <- paste0(
      "  // take only relevant parts of residual correlation matrix \n",
      collapse(rescor_genC)
    )
  }
  out
}

stan_fe <- function(bterms, data, prior, center_X = TRUE,
                    sparse = FALSE, order_mixture = 'none') {
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
  rm_intercept <- center_X || is.cor_bsts(bterms$autocor) || is_ordinal(family)
  if (rm_intercept) {
    fixef <- setdiff(fixef, "Intercept")
  }
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  resp <- usc(bterms$resp)
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
    ct <- ifelse(center_X, "c", "")
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
        hs_scale_global <- paste0(hs_scale_global, " * sigma", resp)
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
    # centering of the fixed effects design matrix improves convergence
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
      # ordinal families either use thres - mu or mu - thres
      # both implies adding <mean_X, b> to the temporary intercept
      sign <- if (is_ordinal(family)) " + " else " - "
      sub_X_means <- paste0(sign, "dot_product(means_X", p, ", b", p, ")")
    } else {
      sub_X_means <- ""
    }
    if (is_ordinal(family)) {
      # intercepts in ordinal models require special treatment
      type <- ifelse(family$family == "cumulative", "ordered", "vector")
      intercept <- paste0(
        "  ", type, "[ncat", p, "-1] temp", p, "_Intercept;",
        "  // temporary thresholds \n"
      )
      if (family$threshold == "flexible") {
        str_add(out$par) <- intercept
      } else if (family$threshold == "equidistant") {
        str_add(out$par) <- paste0(
          "  real temp", p, "_Intercept1;  // threshold 1 \n",
          "  real", subset2(prior, class = "delta", ls = px)$bound,
          " delta", p, ";  // distance between thresholds \n"
        )
        str_add(out$tparD) <- intercept
        str_add(out$tparC1) <- paste0(
          "  // compute equidistant thresholds \n",
          "  for (k in 1:(ncat", p, " - 1)) { \n",
          "    temp", p, "_Intercept[k] = temp", p, "_Intercept1", 
          " + (k - 1.0) * delta", p, "; \n",
          "  } \n"
        )
        str_add(out$prior) <- stan_prior(prior, class = "delta", px = px)
      }
      str_add(out$genD) <- paste0(
        "  // compute actual thresholds \n",
        "  vector[ncat", p, " - 1] b", p, "_Intercept",  
        " = temp", p, "_Intercept", sub_X_means, "; \n" 
      )
    } else {
       if (identical(dpar_class(px$dpar), order_mixture)) {
         # identify mixtures via ordering of the intercepts
         ap_id <- dpar_id(px$dpar)
         resp <- usc(px$resp)
         str_add(out$tparD) <- paste0(
           "  // identify mixtures via ordering of the intercepts \n",                   
           "  real temp", p, "_Intercept",
           " = ordered_Intercept", resp, "[", ap_id, "]; \n"
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
      stop2(
        "Identifying mixture components via ordering requires ",
        "population-level intercepts to be present.\n",
        "Try setting order = 'none' in function 'mixture'."
      )
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
  out <- list()
  r <- subset2(ranef, id = id)
  ccov <- r$group[1] %in% names(cov_ranef)
  ng <- seq_along(r$gcall[[1]]$groups)
  px <- check_prefix(r)
  idp <- paste0(r$id, usc(combine_prefix(px)))
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
  has_def_type <- !r$type %in% "sp"
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

stan_sm <- function(bterms, data, prior) {
  # Stan code of smooth terms
  out <- list()
  smooths <- get_sm_labels(bterms, data = data)
  px <- check_prefix(bterms)
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

stan_cs <- function(bterms, data, ranef, prior) {
  # Stan code for category specific effects
  # (!) Not implemented for non-linear models
  out <- list()
  csef <- colnames(get_model_matrix(bterms$cs, data))
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  ranef <- subset2(ranef, type = "cs", ls = px)
  if (length(csef)) {
    str_add(out$data) <- paste0(
      "  int<lower=1> Kcs", p, ";  // number of category specific effects\n",
      "  matrix[N, Kcs", p, "] Xcs", p, ";  // category specific design matrix\n"
    )
    bound <- get_bound(prior, class = "b", px = px)
    str_add(out$par) <- paste0(
      "  matrix", bound, "[Kcs", p, ", ncat", p, " - 1] bcs", p, ";",
      "  // category specific effects\n"
    )
    str_add(out$modelD) <- paste0(
      "  // linear predictor for category specific effects\n",
      "  matrix[N, ncat", p, " - 1] mucs", p, " = Xcs", p, " * bcs", p, ";\n"
    ) 
    str_add(out$prior) <- stan_prior(
      prior, class = "b", coef = csef,
      suffix = "cs", px = px, matrix = TRUE
    )
  }
  if (nrow(ranef)) {
    if (!length(csef)) {
      # only group-level category specific effects present
      str_add(out$modelD) <- paste0(
        "  // linear predictor for category specific effects \n",               
        "  matrix[N, ncat", p, " - 1] mucs", p, 
        " = rep_matrix(0, N, ncat", p, " - 1);\n"
      )
    }
    cats_regex <- "(?<=\\[)[[:digit:]]+(?=\\]$)"
    cats <- get_matches(cats_regex, ranef$coef, perl = TRUE)
    ncatM1 <- max(as.numeric(cats))
    for (i in seq_len(ncatM1)) {
      r_cat <- ranef[grepl(paste0("\\[", i, "\\]$"), ranef$coef), ]
      str_add(out$modelC2) <- paste0(
        "    mucs", p, "[n, ", i, "] = mucs", p, "[n, ", i, "]"
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
      str_add(out$modelC2) <- ";\n"
    }
  }
  out
}

stan_sp <- function(bterms, data, ranef, prior) {
  # Stan code for special effects
  out <- list()
  spef <- tidy_spef(bterms, data)
  if (is.null(spef)) {
    return(out)
  }
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
  for (i in seq_len(nrow(spef))) {
    eta <- spef$call_prod[[i]]
    if (!is.null(spef$call_mo[[i]])) {
      new_mo <- paste0(
        "mo(simo", p, "_", spef$Imo[[i]], 
        ", Xmo", p, "_", spef$Imo[[i]], "[n])"
      )
      eta <- rename(eta, spef$call_mo[[i]], new_mo)
    }
    if (!is.null(spef$call_me[[i]])) {
      new_me <- paste0("Xme_", seq_along(spef$uni_me[[i]]), "[n]")
      eta <- rename(eta, spef$uni_me[[i]], new_me)
    }
    if (!is.null(spef$call_mi[[i]])) {
      new_mi <- paste0("Yf_", spef$vars_mi[[i]], "[n]")
      eta <- rename(eta, spef$call_mi[[i]], new_mi)
    }
    if (spef$Ic[i] > 0) {
      str_add(eta) <- paste0(" * Csp", p, "_", spef$Ic[i], "[n]")
    }
    r <- subset2(ranef, coef = spef_coef[i])
    rpars <- if (nrow(r)) paste0(" + ", stan_eta_r(r))
    str_add(out$eta) <- paste0(" + (bsp", p, "[", i, "]", rpars, ") * ", eta)
  }
  # prepare general Stan code
  ncovars <- max(spef$Ic)
  str_add(out$data) <- paste0(
    "  int<lower=1> Ksp", p, ";  // number of special effects terms\n",
    if (ncovars > 0L) paste0(
      "  // covariates of special effects terms\n",
      collapse("  vector[N] Csp", p, "_", seq_len(ncovars), ";\n")
    )
  )
  bound <- get_bound(prior, class = "b", px = px)
  str_add(out$par) <- paste0(
    "  // special effects coefficients \n", 
    "  vector", bound, "[Ksp", p, "] bsp", p, "; \n"
  )
  str_add(out$prior) <- stan_prior(
    prior, class = "b", coef = spef$coef, 
    px = px, suffix = paste0("sp", p)
  )
  # include special Stan code for monotonic effects
  I <- unlist(spef$Imo)
  if (length(I)) {
    I <- seq_len(max(I))
    str_add(out$data) <- paste0(
      "  int<lower=1> Imo", p, ";  // number of monotonic variables\n",
      "  int<lower=2> Jmo", p, "[Imo", p, "];  // length of simplexes\n",
      "  // monotonic variables \n",
      collapse("  int Xmo", p, "_", I, "[N];\n"),
      "  // prior concentration of monotonic simplexes\n",
      collapse("  vector[Jmo", p, "[", I, "]] con_simo", p, "_", I, ";\n")
    )
    str_add(out$par) <- paste0(
      "  // simplexes of monotonic effects \n",
      collapse("  simplex[Jmo", p, "[", I, "]] simo", p, "_", I, "; \n")
    ) 
    str_add(out$prior) <- collapse(
      "  target += dirichlet_lpdf(",
      "simo", p, "_", I, " | con_simo", p, "_", I, "); \n"
    )
  }
  out
}

stan_gp <- function(bterms, data, prior) {
  # Stan code for latent gaussian processes
  out <- list()
  gpef <- get_gp_labels(bterms, data = data)
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
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
    str_add(out$prior) <- paste0(
      stan_prior(prior, class = "sdgp", coef = gpef[i], 
                 px = px, suffix = pi),
      stan_prior(prior, class = "lscale", coef = gpef[i], 
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

stan_eta_autocor <- function(autocor, px = list()) {
  # Stan code for the linear predictor of certain autocorrelation terms 
  out <- ""
  p <- usc(combine_prefix(px))
  if (get_ma(autocor) && !use_cov(autocor)) {
    str_add(out) <- paste0(" + head(E", p, "[n], Kma", p, ") * ma", p)
  }
  if (is.cor_car(autocor)) {
    str_add(out) <- paste0(" + rcar", p, "[Jloc", p, "[n]]")
  }
  if (is.cor_bsts(autocor)) {
    str_add(out) <- paste0(" + loclev", p, "[n]")
  }
  out
}

stan_eta_transform <- function(family, llh_adj = FALSE) {
  # indicate whether eta needs to be transformed
  # manually using the link functions
  # Args:
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
  #   dpars: names of distributional parameters
  #   adforms: list of formulas containing addition terms
  stopifnot(all(c("family", "link") %in% names(family)))
  llh_adj <- stan_llh_adj(adforms, c("cens", "trunc"))
  if (stan_eta_transform(family, llh_adj = llh_adj)) {
    link <- family$link
    family <- family$family
    shape <- paste0("shape", mix)
    shape <- ifelse(shape %in% dpars, paste0(shape, "[n]"), shape)
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

stan_dpar_defs <- function(dpar, suffix = "") {
  # default Stan definitions for distributional parameters
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
    def <- paste0(def[1], dpar, suffix, def[2])
  } else {
    def <- ""
  }
  def
}

stan_dpar_defs_temp <- function(dpar, suffix = "") {
  # default Stan definitions for temporary distributional parameters
  default_defs <- list(
    xi = c(
      "  real temp_", 
      ";  // unscaled shape parameter \n"
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
  if (any(families %in% "categorical")) {
    str_add(out$modelD) <- paste0( 
      "  // linear predictor matrix \n",
      "  vector[ncat", p, "] mu", p, "[N]; \n"
    )
    dpars <- names(bterms$dpars)
    mu_vector <- stan_vector(c("0", paste0(dpars, p, "[n]")))
    str_add(out$modelC4) <- paste0(
      "    mu", p, "[n] = ", mu_vector, ";\n"
    )
  }
  if (any(families %in% "skew_normal")) {
    # as suggested by Stephen Martin use sigma and mu of CP 
    # but the skewness parameter alpha of DP
    ap_names <- names(bterms$dpars)
    for (i in which(families %in% "skew_normal")) {
      id <- ifelse(length(families) == 1L, "", i)
      sigma <- stan_sigma_transform(bterms, id = id)
      ns <- ifelse(grepl("\\[n\\]", sigma), "[n]", "")
      na <- ifelse(paste0("alpha", id) %in% ap_names, "[n]", "")
      type_delta <- ifelse(nzchar(na), "vector[N]", "real")
      no <- ifelse(any(nzchar(c(ns, na))), "[n]", "")
      type_omega <- ifelse(nzchar(no), "vector[N]", "real")
      str_add(out$modelD) <- paste0(
        "  ", type_delta, " delta", id, p, "; \n",
        "  ", type_omega, " omega", id, p, "; \n"
      )
      alpha <- paste0("alpha", id, p, na)
      delta <- paste0("delta", id, p, na)
      omega <- paste0("omega", id, p, no)
      comp_delta <- paste0(
        "  ", delta, " = ", alpha, " / sqrt(1 + ", alpha, "^2); \n"
      )
      comp_omega <- paste0(
        "  ", omega, " = ", sigma, 
        " / sqrt(1 - sqrt_2_div_pi^2 * ", delta, "^2); \n"
      )
      str_add(out$modelC5) <- paste0(
        if (!nzchar(na)) comp_delta,
        if (!nzchar(no)) comp_omega,
        "  for (n in 1:N) { \n",
        if (nzchar(na)) paste0("  ", comp_delta),
        if (nzchar(no)) paste0("  ", comp_omega),
        "    mu", id, p, "[n] = mu", id, p, "[n] - ", 
        omega, " * ", delta, " * sqrt_2_div_pi; \n",
        "  } \n"
      )
    }
  }
  if (any(families %in% "gen_extreme_value")) {
    ap_names <- c(names(bterms$dpars), names(bterms$fdpars))
    for (i in which(families %in% "gen_extreme_value")) {
      id <- ifelse(length(families) == 1L, "", i)
      xi <- paste0("xi", id)
      if (!xi %in% ap_names) {
        str_add(out$modelD) <- paste0(
          "  real ", xi, ";  // scaled shape parameter \n"
        )
        sigma <- paste0("sigma", id)
        v <- ifelse(sigma %in% names(bterms$dpars), "_vector", "")
        args <- sargs(
          paste0("temp_", xi), paste0("Y", p), 
          paste0("mu", id, p), paste0(sigma, p)
        )
        str_add(out$modelC5) <- paste0(
          "  ", xi, p, " = scale_xi", v, "(", args, "); \n"
        )
      }
    }
  }
  out
}

stan_sigma_transform <- function(bterms, id = "") {
  # Stan code for sigma to incorporate addition argument 'se'
  if (nzchar(id)) {
    family <- family_names(bterms)[id]
  } else {
    family <- bterms$family$family
  }
  p <- usc(combine_prefix(bterms))
  ns <- ifelse(paste0("sigma", id) %in% names(bterms$dpars), "[n]", "")
  has_sigma <- has_sigma(family, bterms)
  sigma <- ifelse(has_sigma, paste0("sigma", id, p, ns), "")
  if (is.formula(bterms$adforms$se)) {
    sigma <- ifelse(
      nzchar(sigma), 
      paste0("sqrt(", sigma, "^2 + se2", p, "[n])"), 
      paste0("se", p, "[n]")
    )
  }
  sigma
}
