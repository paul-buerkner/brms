# unless otherwise specified, functions return a named list
# of Stan code snippets to be pasted together later on

# generate stan code for predictor terms
stan_predictor <- function(x, ...) {
  UseMethod("stan_predictor")
}

# combine effects for the predictors of a single (non-linear) parameter
# @param ... arguments passed to the underlying effect-specific functions
#' @export
stan_predictor.bframel <- function(x, ...) {
  out <- collapse_lists(
    stan_fe(x, ...),
    stan_thres(x, ...),
    stan_sp(x, ...),
    stan_cs(x, ...),
    stan_sm(x, ...),
    stan_gp(x, ...),
    stan_ac(x, ...),
    stan_offset(x, ...),
    stan_bhaz(x, ...)
  )
  out <- stan_special_prior(x, out = out, ...)
  out <- stan_eta_combine(x, out = out, ...)
  out
}

# prepare Stan code for non-linear terms
#' @export
stan_predictor.bframenl <- function(x, ...) {
  collapse_lists(
    stan_nl(x, ...),
    stan_thres(x, ...),
    stan_bhaz(x, ...),
    stan_ac(x, ...)
  )
}

#' @export
stan_predictor.brmsframe <- function(x, prior, normalize, ...,
                                     stanvars = NULL) {
  px <- check_prefix(x)
  resp <- usc(combine_prefix(px))
  out <- list()
  str_add_list(out) <- stan_response(x, normalize = normalize, ...)
  valid_dpars <- valid_dpars(x)
  family_files <- family_info(x, "include")
  if (length(family_files)) {
    str_add(out$fun) <- cglue("  #include '{family_files}'\n")
  }
  s2z_mean_noncenter <- re_s2z_cross_mean_noncenter_prefixes(
    x, prior, stanvars = stanvars
  )
  args <- nlist(
    prior, normalize, stanvars, nlpars = names(x$nlpars),
    s2z_mean_noncenter, ...
  )
  args$primitive <- use_glm_primitive(x) || use_glm_primitive_categorical(x)
  for (nlp in names(x$nlpars)) {
    nlp_args <- list(x$nlpars[[nlp]])
    str_add_list(out) <- do_call(stan_predictor, c(nlp_args, args))
  }
  for (dp in valid_dpars) {
    dp_terms <- x$dpars[[dp]]
    dp_comment <- stan_dpar_comments(dp, family = x$family)
    if (is.btl(dp_terms) || is.btnl(dp_terms)) {
      # distributional parameter is predicted
      str_add_list(out) <- do_call(stan_predictor, c(list(dp_terms), args))
    } else if (is.numeric(x$fdpars[[dp]]$value)) {
      # distributional parameter is fixed to constant
      if (is_mix_proportion(dp, family = x$family)) {
        # mixture proportions are handled in 'stan_mixture'
        next
      }
      dp_value <- x$fdpars[[dp]]$value
      dp_comment <- stan_comment(dp_comment)
      str_add(out$tpar_def) <- glue(
        "  real {dp}{resp} = {dp_value};{dp_comment}\n"
      )
      str_add(out$pll_args) <- glue(", real {dp}{resp}")
    } else if (is.character(x$fdpars[[dp]]$value)) {
      # distributional parameter is fixed to another distributional parameter
      if (!x$fdpars[[dp]]$value %in% valid_dpars) {
        stop2("Parameter '", x$fdpars[[dp]]$value, "' cannot be found.")
      }
      if (is_mix_proportion(dp, family = x$family)) {
        stop2("Cannot set mixture proportions to be equal.")
      }
      dp_value <- x$fdpars[[dp]]$value
      dp_comment <- stan_comment(dp_comment)
      str_add(out$tpar_def) <- glue(
        "  real {dp}{resp};{dp_comment}\n"
      )
      str_add(out$tpar_comp) <- glue(
        "  {dp}{resp} = {dp_value}{resp};\n"
      )
      str_add(out$pll_args) <- glue(", real {dp}{resp}")
    } else {
      # distributional parameter is estimated as a scalar
      if (is_mix_proportion(dp, family = x$family)) {
        # mixture proportions are handled in 'stan_mixture'
        next
      }
      prefix <- ""
      if (dp %in% valid_dpars(x, type = "tmp")) {
        # some parameters are fully computed only after the model is run
        prefix <- "tmp_"
        dp_comment <- paste0(dp_comment, " (temporary)")
      }
      str_add_list(out) <- stan_prior(
        prior, dp, prefix = prefix, suffix = resp,
        header_type = "real", px = px,
        comment = dp_comment, normalize = normalize
      )
    }
  }
  str_add_list(out) <- stan_dpar_transform(
    x, prior = prior, normalize = normalize, ...
  )
  str_add_list(out) <- stan_mixture(
    x, prior = prior, normalize = normalize, ...
  )
  out$model_log_lik <- stan_log_lik(
    x, normalize = normalize, ...
  )
  list(out)
}

#' @export
stan_predictor.mvbrmsframe <- function(x, prior, threads, normalize, ...) {
  out <- lapply(x$terms, stan_predictor, prior = prior, threads = threads,
                normalize = normalize, ...)
  out <- unlist(out, recursive = FALSE)
  if (!x$rescor) {
    return(out)
  }
  resp_type <- out[[1]]$resp_type
  out <- collapse_lists(ls = out)
  out$resp_type <- "vector"
  adforms <- from_list(x$terms, "adforms")
  adnames <- unique(ulapply(adforms, names))
  adallowed <- c("se", "weights", "mi")
  if (!all(adnames %in% adallowed))  {
    stop2("Only ", collapse_comma(adallowed), " are supported ",
          "addition arguments when 'rescor' is estimated.")
  }
  # we already know at this point that all families are identical
  family <- family_names(x)[1]
  stopifnot(family %in% c("gaussian", "student"))
  resp <- x$responses
  nresp <- length(resp)
  str_add(out$model_def) <- glue(
    "  // multivariate predictor array\n",
    "  array[N] vector[nresp] Mu;\n"
  )
  str_add(out$model_comp_mvjoin) <- glue(
    "    Mu[n] = {stan_vector(glue('mu_{resp}[n]'))};\n"
  )
  str_add(out$data) <- glue(
    "  int<lower=1> nresp;  // number of responses\n",
    "  int nrescor;  // number of residual correlations\n"
  )
  str_add(out$pll_args) <- glue(", data int nresp")
  str_add(out$tdata_def) <- glue(
    "  array[N] vector[nresp] Y;  // response array\n"
  )
  str_add(out$tdata_comp) <- glue(
    "  for (n in 1:N) {{\n",
    "    Y[n] = {stan_vector(glue('Y_{resp}[n]'))};\n",
    "  }}\n"
  )
  str_add(out$pll_args) <- ", data array[] vector Y"
  if (any(adnames %in% "weights")) {
    str_add(out$tdata_def) <- glue(
      "  // weights of the pointwise log-likelihood\n",
      "  vector<lower=0>[N] weights = weights_{resp[1]};\n"
    )
    str_add(out$pll_args) <- glue(", data vector weights")
  }
  miforms <- rmNULL(from_list(adforms, "mi"))
  if (length(miforms)) {
    str_add(out$model_no_pll_def) <- " array[N] vector[nresp] Yl = Y;\n"
    str_add(out$pll_args) <- ", array[] vector Yl"
    for (i in seq_along(miforms)) {
      j <- match(names(miforms)[i], resp)
      # needs to happen outside of reduce_sum
      # to maintain consistency of indexing Yl
      str_add(out$model_no_pll_comp_mvjoin) <- glue(
        "    Yl[n][{j}] = Yl_{resp[j]}[n];\n"
      )
    }
  }
  str_add_list(out) <- stan_prior(
    prior, class = "Lrescor",
    type = "cholesky_factor_corr[nresp]", header_type = "matrix",
    comment = "parameters for multivariate linear models",
    normalize = normalize
  )
  if (family == "student") {
    str_add_list(out) <- stan_prior(
      prior, class = "nu", header_type = "real",
      normalize = normalize
    )
  }
  sigma <- ulapply(x$terms, stan_sigma_transform, threads = threads)
  if (any(grepl(stan_nn_regex(), sigma))) {
    str_add(out$model_def) <- "  array[N] vector[nresp] sigma;\n"
    str_add(out$model_comp_mvjoin) <- glue(
      "    sigma[n] = {stan_vector(sigma)};\n"
    )
    if (family == "gaussian") {
      str_add(out$model_def) <- glue(
        "  // cholesky factor of residual covariance matrix\n",
        "  array[N] matrix[nresp, nresp] LSigma;\n"
      )
      str_add(out$model_comp_mvjoin) <- glue(
        "    LSigma[n] = diag_pre_multiply(sigma[n], Lrescor);\n"
      )
    } else if (family == "student") {
      str_add(out$model_def) <- glue(
        "  // residual covariance matrix\n",
        "  array[N] matrix[nresp, nresp] Sigma;\n"
      )
      str_add(out$model_comp_mvjoin) <- glue(
        "    Sigma[n] = multiply_lower_tri_self_transpose(",
        "diag_pre_multiply(sigma[n], Lrescor));\n"
      )
    }
  } else {
    str_add(out$model_def) <- glue(
      "  vector[nresp] sigma = {stan_vector(sigma)};\n"
    )
    if (family == "gaussian") {
      str_add(out$model_def) <- glue(
        "  // cholesky factor of residual covariance matrix\n",
        "  matrix[nresp, nresp] LSigma = ",
        "diag_pre_multiply(sigma, Lrescor);\n"
      )
    } else if (family == "student") {
      str_add(out$model_def) <- glue(
        "  // residual covariance matrix\n",
        "  matrix[nresp, nresp] Sigma = ",
        "multiply_lower_tri_self_transpose(",
        "diag_pre_multiply(sigma, Lrescor));\n"
      )
    }
  }
  str_add(out$gen_def) <- glue(
    "  // residual correlations\n",
    "  corr_matrix[nresp] Rescor",
    " = multiply_lower_tri_self_transpose(Lrescor);\n",
    "  vector<lower=-1,upper=1>[nrescor] rescor;\n"
  )
  str_add(out$gen_comp) <- stan_cor_gen_comp("rescor", "nresp")
  out$model_comp_mvjoin <- paste0(
    "  // combine univariate parameters\n",
    "  for (n in 1:N) {\n",
    stan_nn_def(threads),
    out$model_comp_mvjoin,
    "  }\n"
  )
  if (isTRUE(nzchar(out$model_no_pll_comp_mvjoin))) {
    out$model_no_pll_comp_mvjoin <- paste0(
      "  // combine univariate parameters\n",
      "  for (n in 1:N) {\n",
      out$model_no_pll_comp_mvjoin,
      "  }\n"
    )
  }
  out$model_log_lik <- stan_log_lik(
    x, threads = threads, normalize = normalize, ...
  )
  list(out)
}

# Stan code for population-level effects
stan_fe <- function(bframe, prior, stanvars, threads, primitive,
                    normalize, s2z_mean_noncenter = character(), ...) {
  stopifnot(is.bframel(bframe))
  out <- list()
  family <- bframe$family
  fixef <- bframe$frame$fe$vars_stan
  sparse <- bframe$frame$fe$sparse
  decomp <- bframe$frame$fe$decomp
  center <- bframe$frame$fe$center
  ct <- str_if(center, "c")
  px <- check_prefix(bframe)
  p <- usc(combine_prefix(px))
  resp <- usc(px$resp)
  lpdf <- stan_lpdf_name(normalize)
  # Strict latent-score blocks have no omitted mean and therefore do not
  # replace this predictor's ordinary population coefficients by theta_s2z.
  s2z <- has_re_s2z_conventional(bframe)
  mean_noncenter <- p %in% s2z_mean_noncenter
  info_s2z <- NULL
  ordinal_s2z <- FALSE
  if (s2z) {
    infos_s2z <- re_s2z_infos(
      bframe, prior = prior, stanvars = stanvars
    )
    stopifnot(length(infos_s2z) > 0L)
    info_s2z <- infos_s2z[[1L]]
    ordinal_s2z <- isTRUE(info_s2z$ordinal)
  }

  if (s2z) {
    K_s2z <- length(info_s2z$qnames)
    active_s2z <- info_s2z$active_q
    if (is.null(active_s2z)) {
      active_s2z <- sort(unique(unlist(lapply(infos_s2z, function(info) {
        if (ordinal_s2z) {
          which(rowSums(abs(info$H)) > 0)
        } else {
          rows <- info$match_q
          if (info$center && any(info$r$coef != "Intercept")) {
            rows <- c(rows, 1L)
          }
          rows
        }
      }), use.names = FALSE)))
    }
    coef_s2z <- if (ordinal_s2z) info_s2z$slope_q else seq_len(K_s2z)
    active_coef_s2z <- intersect(coef_s2z, active_s2z)
    inactive_coef_s2z <- setdiff(coef_s2z, active_coef_s2z)
    split_fe_s2z <- ordinal_s2z || length(inactive_coef_s2z) > 0L
    if (split_fe_s2z) {
      # Threshold primitives are owned by stan_thres(), while fixed-only
      # coefficients retain their ordinary brms priors under every physical
      # omitted-mean kernel. Assemble just the FE-owned coordinates into one
      # vector before the omitted-mean solve.
      stopifnot(!mean_noncenter)
      special_b_s2z <- has_special_prior(prior, bframe, class = "b")
      if (length(active_coef_s2z)) {
        str_add(out$par) <- glue(
          "  vector[{length(active_coef_s2z)}] theta_s2z_active{p};",
          "  // S2Z-active finite-population coefficients\n"
        )
      }
      if (length(inactive_coef_s2z)) {
        inactive_names <- info_s2z$qnames[inactive_coef_s2z]
        stopifnot(!"Intercept" %in% inactive_names)
        if (special_b_s2z) {
          stopifnot(length(inactive_names) == length(fixef))
          inactive_prior <- stan_prior_non_centered(
            class = "fixed_s2z", suffix = p, suffix_K = ct,
            normalize = normalize
          )
          inactive_prior$pll_args <- NULL
          str_add_list(out) <- inactive_prior
        } else {
          inactive_specs <- lapply(inactive_names, function(coef) {
            tryCatch(
              re_s2z_prior(
                prior, bframe, class = "b", coef = coef,
                r = info_s2z$r, stanvars = stanvars
              ),
              error = function(e) NULL
            )
          })
          broadcasted_inactive <- any(vapply(
            inactive_specs,
            function(spec) !is.null(spec) && isTRUE(spec$broadcasted),
            logical(1)
          ))
          if (broadcasted_inactive) {
            # A vector-valued global prior belongs to the original full b
            # vector. Score fixed-only coordinates one at a time so each uses
            # its original vector index rather than a shortened inactive
            # vector, which would change the ordinary brms target.
            if (any(vapply(inactive_specs, is.null, logical(1)))) {
              stop2(
                "Vector-valued global priors cannot currently be mixed with ",
                "unsupported coefficient-specific priors on S2Z-inactive ",
                "coordinates. Use scalar or supported coefficient-specific ",
                "priors for this predictor."
              )
            }
            str_add(out$par) <- glue(
              "  vector[{length(inactive_names)}] fixed_s2z{p};",
              "  // S2Z-inactive regression coefficients\n"
            )
            for (a in seq_along(inactive_specs)) {
              str_add(out$tpar_prior) <- stan_re_s2z_prior_target(
                inactive_specs[[a]],
                par = glue("fixed_s2z{p}[{a}]"),
                normalize = normalize
              )
            }
          } else {
            inactive_prior <- stan_prior(
              prior, class = "b", coef = inactive_names,
              type = glue("vector[{length(inactive_coef_s2z)}]"),
              suffix = glue("_s2z_inactive{p}"), px = px,
              comment = "S2Z-inactive regression coefficients",
              normalize = normalize
            )
            inactive_prior <- lapply(inactive_prior, function(code) {
              gsub(
                glue("b_s2z_inactive{p}"), glue("fixed_s2z{p}"), code,
                fixed = TRUE
              )
            })
            str_add_list(out) <- inactive_prior
          }
        }
      }
      str_add(out$tpar_def) <- glue(
        "  vector[{K_s2z}] theta_s2z{p};",
        "  // assembled finite-population coefficients\n"
      )
      theta_comp <- if (special_b_s2z) "tpar_special_prior" else "tpar_comp"
      for (a in seq_along(active_coef_s2z)) {
        str_add(out[[theta_comp]]) <- glue(
          "  theta_s2z{p}[{active_coef_s2z[a]}] = ",
          "theta_s2z_active{p}[{a}];\n"
        )
      }
      for (a in seq_along(inactive_coef_s2z)) {
        str_add(out[[theta_comp]]) <- glue(
          "  theta_s2z{p}[{inactive_coef_s2z[a]}] = ",
          "fixed_s2z{p}[{a}];\n"
        )
      }
    } else if (mean_noncenter) {
      str_add(out$par) <- glue(
        "  vector[{K_s2z}] z_theta_s2z{p};",
        "  // standardized finite-population coefficients for S2Z effects\n"
      )
      str_add(out$tpar_def) <- glue(
        "  vector[{K_s2z}] theta_s2z{p};",
        "  // finite-population coefficients for physical S2Z effects\n"
      )
    } else {
      str_add(out$par) <- glue(
        "  vector[{K_s2z}] theta_s2z{p};",
        "  // finite-population coefficients for physical S2Z effects\n"
      )
    }
    str_add(out$pll_args) <- glue(", vector theta_s2z{p}")
  }

  if (length(fixef)) {
    str_add(out$data) <- glue(
      "  int<lower=1> K{p};",
      "  // number of population-level effects\n",
      "  matrix[N{resp}, K{p}] X{p};",
      "  // population-level design matrix\n"
    )
    if (decomp == "none") {
      str_add(out$pll_args) <- glue(", data matrix X{ct}{p}")
    }
    if (sparse) {
      if (decomp != "none") {
        stop2("Cannot use ", decomp, " decomposition for sparse matrices.")
      }
      if (use_threading(threads)) {
        stop2("Cannot use threading and sparse matrices at the same time.")
      }
      str_add(out$tdata_def) <- glue(
        "  // sparse matrix representation of X{p}\n",
        "  vector[rows(csr_extract_w(X{p}))] wX{p}",
        " = csr_extract_w(X{p});\n",
        "  int vX{p}[size(csr_extract_v(X{p}))]",
        " = csr_extract_v(X{p});\n",
        "  int uX{p}[size(csr_extract_u(X{p}))]",
        " = csr_extract_u(X{p});\n"
      )
    }
    if (!s2z) {
      # prepare population-level coefficients
      b_type <- glue("vector[K{ct}{p}]")
      has_special_prior <- has_special_prior(prior, bframe, class = "b")
      if (decomp == "none") {
        if (has_special_prior) {
          str_add_list(out) <- stan_prior_non_centered(
            suffix = p, suffix_K = ct, normalize = normalize
          )
        } else {
          str_add_list(out) <- stan_prior(
            prior, class = "b", coef = fixef, type = b_type,
            px = px, suffix = p, header_type = "vector",
            comment = "regression coefficients", normalize = normalize
          )
        }
      } else {
        stopifnot(decomp == "QR")
        stopif_prior_bound(prior, class = "b", ls = px)
        if (has_special_prior) {
          str_add_list(out) <- stan_prior_non_centered(
            suffix = p, suffix_class = "Q", suffix_K = ct,
            normalize = normalize
          )
        } else {
          str_add_list(out) <- stan_prior(
            prior, class = "b", coef = fixef, type = b_type,
            px = px, suffix = glue("Q{p}"), header_type = "vector",
            comment = "regression coefficients on QR scale",
            normalize = normalize
          )
        }
        str_add(out$gen_def) <- glue(
          "  // obtain the actual coefficients\n",
          "  vector[K{ct}{p}] b{p} = XR{p}_inv * bQ{p};\n"
        )
      }
    }
  }

  order_intercepts <- order_intercepts(bframe)
  if (order_intercepts && !center) {
    stop2(
      "Identifying mixture components via ordering requires ",
      "population-level intercepts to be present.\n",
      "Try setting order = 'none' in function 'mixture'."
    )
  }
  if (center) {
    # centering the design matrix improves convergence
    sub_X_means <- ""
    if (length(fixef)) {
      str_add(out$data) <- glue(
        "  int<lower=1> Kc{p};",
        "  // number of population-level effects after centering\n"
      )
      sub_X_means <- glue(" - dot_product(means_X{p}, b{p})")
      if (is_ordinal(family)) {
        str_add(out$tdata_def) <- glue(
          "  matrix[N{resp}, Kc{p}] Xc{p};",
          "  // centered version of X{p}\n",
          "  vector[Kc{p}] means_X{p};",
          "  // column means of X{p} before centering\n"
        )
        str_add(out$tdata_comp) <- glue(
          "  for (i in 1:K{p}) {{\n",
          "    means_X{p}[i] = mean(X{p}[, i]);\n",
          "    Xc{p}[, i] = X{p}[, i] - means_X{p}[i];\n",
          "  }}\n"
        )
      } else {
        str_add(out$tdata_def) <- glue(
          "  matrix[N{resp}, Kc{p}] Xc{p};",
          "  // centered version of X{p} without an intercept\n",
          "  vector[Kc{p}] means_X{p};",
          "  // column means of X{p} before centering\n"
        )
        str_add(out$tdata_comp) <- glue(
          "  for (i in 2:K{p}) {{\n",
          "    means_X{p}[i - 1] = mean(X{p}[, i]);\n",
          "    Xc{p}[, i - 1] = X{p}[, i] - means_X{p}[i - 1];\n",
          "  }}\n"
        )
      }
    }
    if (!is_ordinal(family) && s2z) {
      str_add(out$eta) <- glue(" + theta_s2z{p}[1]")
    } else if (!is_ordinal(family)) {
      # intercepts of ordinal models are handled in 'stan_thres'
      intercept_type <- "real"
      if (order_intercepts) {
        # identify mixtures via ordering of the intercepts
        dp_id <- dpar_id(px$dpar)
        str_add(out$tpar_def) <- glue(
          "  // identify mixtures via ordering of the intercepts\n",
          "  real Intercept{p} = ordered_Intercept{resp}[{dp_id}];\n"
        )
        str_add(out$pll_args) <- glue(", real Intercept{p}")
        # intercept parameter needs to be defined outside of 'stan_prior'
        intercept_type <- ""
      }
      str_add(out$eta) <- glue(" + Intercept{p}")
      str_add(out$gen_def) <- glue(
        "  // actual population-level intercept\n",
        "  real b{p}_Intercept = Intercept{p}{sub_X_means};\n"
      )
      str_add_list(out) <- stan_prior(
        prior, class = "Intercept", type = intercept_type,
        suffix = p, px = px, header_type = "real",
        comment = "temporary intercept for centered predictors",
        normalize = normalize
      )
    }
  }
  if (decomp == "QR") {
    if (!length(fixef)) {
      stop2("QR decomposition requires non-intercept predictors.")
    }
    str_add(out$tdata_def) <- glue(
      "  // matrices for QR decomposition\n",
      "  matrix[N{resp}, K{ct}{p}] XQ{p};\n",
      "  matrix[K{ct}{p}, K{ct}{p}] XR{p};\n",
      "  matrix[K{ct}{p}, K{ct}{p}] XR{p}_inv;\n"
    )
    str_add(out$tdata_comp) <- glue(
      "  // compute and scale QR decomposition\n",
      "  XQ{p} = qr_thin_Q(X{ct}{p}) * sqrt(N{resp} - 1);\n",
      "  XR{p} = qr_thin_R(X{ct}{p}) / sqrt(N{resp} - 1);\n",
      "  XR{p}_inv = inverse(XR{p});\n"
    )
    str_add(out$pll_args) <- glue(", data matrix XQ{p}")
  }
  if (length(fixef) && !primitive) {
    # added in the end such that the intercept comes first in out$eta
    if (s2z) {
      slice <- stan_slice(threads)
      if (ordinal_s2z) {
        X_s2z <- str_if(center, "Xc", "X")
        eta_fe <- glue(
          " + {X_s2z}{p}{slice} * tail(theta_s2z{p}, {length(fixef)})"
        )
      } else if (center) {
        eta_fe <- glue(
          " + Xc{p}{slice} * tail(theta_s2z{p}, {length(fixef)})"
        )
      } else {
        eta_fe <- glue(" + X{p}{slice} * theta_s2z{p}")
      }
    } else if (sparse) {
      stopifnot(!center && decomp == "none")
      csr_args <- sargs(
        paste0(c("rows", "cols"), "(X", p, ")"),
        paste0(c("wX", "vX", "uX", "b"), p)
      )
      eta_fe <- glue(" + csr_matrix_times_vector({csr_args})")
    } else {
      sfx_X <- sfx_b <- ""
      if (decomp == "QR") {
        sfx_X <- sfx_b <- "Q"
      } else if (center) {
        sfx_X <- "c"
      }
      slice <- stan_slice(threads)
      eta_fe <- glue(" + X{sfx_X}{p}{slice} * b{sfx_b}{p}")
    }
    str_add(out$eta) <- eta_fe
  }
  out
}

# Stan code for group-level effects
stan_re <- function(bframe, prior, normalize, ..., stanvars = NULL) {
  lpdf <- ifelse(normalize, "lpdf", "lupdf")
  reframe <- bframe$frame$re
  stopifnot(is.reframe(reframe))
  IDs <- unique(reframe$id)
  out <- list()
  # Multiple S2Z covariance blocks in one linear predictor share the same
  # finite-population coefficients. Their omitted means must consequently be
  # integrated in one joint Gaussian system. Keep the established one-block
  # generators entirely unchanged and activate the joint path for sets with
  # more than one local S2Z ID. Ordinal thresholds also require this general
  # path for a singleton because their H map is affine and may repeat rows.
  joint_s2z_sets <- list()
  for (bfl in Filter(has_re_s2z, all_bframel(bframe))) {
    infos <- re_s2z_infos(
      bfl, prior = prior, stanvars = stanvars
    )
    ordinal <- length(infos) && isTRUE(infos[[1L]]$ordinal)
    cross <- length(infos) && is_re_s2z_cross_id(bframe, infos[[1L]]$id)
    if ((length(infos) > 1L || ordinal) && !cross) {
      set_id <- infos[[1L]]$id
      joint_s2z_sets[[as.character(set_id)]] <- nlist(
        set_id, ids = vapply(infos, `[[`, numeric(1), "id"), bfl, infos
      )
    }
  }
  joint_s2z_by_id <- list()
  for (set in joint_s2z_sets) {
    for (id in set$ids) {
      joint_s2z_by_id[[as.character(id)]] <- set
    }
  }
  # special handling of student-t group effects as their 'df' parameters
  # are defined on a per-group basis instead of a per-ID basis
  reframe_t <- subset_reframe_dist(reframe, "student")
  if (has_rows(reframe_t)) {
    str_add(out$par) <-
      "  // parameters for student-t distributed group-level effects\n"
    for (i in seq_rows(reframe_t)) {
      g <- usc(reframe_t$ggn[i])
      id <- reframe_t$id[i]
      str_add_list(out) <- stan_prior(
        prior, class = "df", group = reframe_t$group[i],
        suffix = g, normalize = normalize
      )
      str_add(out$par) <- glue(
        "  vector<lower=0>[N_{id}] udf{g};\n"
      )
      str_add(out$model_prior) <- glue(
        "  target += inv_chi_square_{lpdf}(udf{g} | df{g});\n"
      )
      # separate definition from computation to support fixed parameters
      str_add(out$tpar_def) <- glue(
        "  vector[N_{id}] dfm{g};\n"
      )
      str_add(out$tpar_comp) <- glue(
        "  dfm{g} = sqrt(df{g} * udf{g});\n"
      )
    }
  }
  # the ID syntax requires group-level effects to be evaluated separately
  tmp <- lapply(IDs, function(id) {
    .stan_re(
      id, bframe = bframe, prior = prior, normalize = normalize,
      joint_s2z = joint_s2z_by_id[[as.character(id)]],
      stanvars = stanvars, ...
    )
  })
  joint <- lapply(joint_s2z_sets, function(set) {
    .stan_re_s2z_joint(
      set, prior = prior, normalize = normalize,
      stanvars = stanvars, ...
    )
  })
  out <- collapse_lists(ls = c(list(out), tmp, joint))
  # Several independent S2Z IDs can require the same Stan helper. Includes
  # define functions in place, so retaining more than one copy would produce
  # duplicate function definitions. Keep the first occurrence after all
  # per-ID snippets have been combined.
  unique_includes <- c(
    "  #include 'fun_sum_to_zero.stan'\n"
  )
  if (!is.null(out$fun)) {
    for (include in unique_includes) {
      include_at <- regexpr(include, out$fun, fixed = TRUE)
      if (include_at[1] > 0L) {
        include_end <-
          include_at[1] + attr(include_at, "match.length") - 1L
        include_tail <- substring(out$fun, include_end + 1L)
        include_tail <- gsub(include, "", include_tail, fixed = TRUE)
        out$fun <- paste0(
          substr(out$fun, 1L, include_end), include_tail
        )
      }
    }
  }
  out
}

# Stan code for group-level effects per ID
# @param id the ID of the grouping factor
.stan_re <- function(id, bframe, prior, threads, normalize,
                     joint_s2z = NULL, ..., stanvars = NULL) {
  lpdf <- ifelse(normalize, "lpdf", "lupdf")
  out <- list()
  r <- subset2(bframe$frame$re, id = id)
  stopifnot(is.reframe(r))
  has_cov <- nzchar(r$cov[1])
  has_by <- nzchar(r$by[[1]])
  has_pw <- isTRUE(nzchar(r$gcall[[1]]$pw))
  Nby <- seq_along(r$bylevels[[1]])
  ng <- seq_along(r$gcall[[1]]$groups)
  px <- check_prefix(r)
  uresp <- usc(unique(px$resp))
  idp <- paste0(r$id, usc(combine_prefix(px)))
  s2z_latent <- isTRUE(r$s2z[1]) && all(re_s2z_latent(r))
  r_cov <- if (s2z_latent) re_s2z_latent_dimensions(r) else r
  px_cov <- check_prefix(r_cov)
  s2z_cross <- isTRUE(r$s2z[1]) && !s2z_latent &&
    is_re_s2z_cross_id(bframe, id)
  center_mode <- re_s2z_center_mode(r)
  center_fisher <- identical(center_mode, "auto")
  s2z_fisher <- isTRUE(r$s2z[1]) && center_fisher
  s2z_ordinal <- isTRUE(r$s2z[1]) && !s2z_latent && any(vapply(
    all_bframel(bframe), function(bfl) {
      r_local <- bfl$frame$re
      if (!has_rows(r_local) || !id %in% r_local$id[r_local$s2z]) {
        return(FALSE)
      }
      infos <- re_s2z_infos(bfl)
      any(vapply(infos, function(info) {
        identical(info$id, id) && isTRUE(info$ordinal)
      }, logical(1)))
    }, logical(1)
  ))
  if (s2z_ordinal && s2z_cross) {
    stop2("An S2Z ID touching ordinal thresholds cannot yet span predictors.")
  }
  if (s2z_ordinal && s2z_fisher) {
    stop2("The automatic-centering precursor is not yet supported for ",
          "ordinal S2Z effects.")
  }
  fisher_info <- NULL
  if (center_fisher) {
    fisher_info <- if (!isTRUE(r$s2z[1])) {
      stan_re_s2z_fisher_info(
        id, r = r, bframe = bframe, threads = threads
      )
    } else if (s2z_latent) {
      stan_re_s2z_latent_fisher_info(
        id, r = r, bframe = bframe, threads = threads
      )
    } else if (s2z_cross) {
      stan_re_s2z_cross_fisher_info(
        id, r = r, bframe = bframe, threads = threads
      )
    } else {
      stan_re_s2z_fisher_info(
        id, r = r, bframe = bframe, threads = threads
      )
    }
    str_add(out$fun) <- fisher_info$fun %||% ""
  }
  # define data needed for group-level effects
  str_add(out$data) <- glue(
    "  // data for group-level effects of ID {id}\n",
    "  int<lower=1> N_{id};  // number of grouping levels\n",
    "  int<lower=1> M_{id};  // number of coefficients per level\n"
  )
  if (r$gtype[1] == "mm") {
    for (res in uresp) {
      str_add(out$data) <- cglue(
        "  array[N{res}] int<lower=1> J_{id}{res}_{ng};",
        "  // grouping indicator per observation\n",
        "  array[N{res}] real W_{id}{res}_{ng};",
        "  // multi-membership weights\n"
      )
      str_add(out$pll_args) <- cglue(
        ", data array[] int J_{id}{res}_{ng}, data array[] real W_{id}{res}_{ng}"
      )
    }
  } else {
    str_add(out$data) <- cglue(
      "  array[N{uresp}] int<lower=1> J_{id}{uresp};",
      "  // grouping indicator per observation\n"
    )
    str_add(out$pll_args) <- cglue(
      ", data array[] int J_{id}{uresp}"
    )
  }
  if (has_by) {
    str_add(out$data) <- glue(
      "  int<lower=1> Nby_{id};  // number of by-factor levels\n",
      "  array[N_{id}] int<lower=1> Jby_{id};",
      "  // by-factor indicator per observation\n"
    )
  }
  if (has_cov) {
    str_add(out$data) <- glue(
      "  matrix[N_{id}, N_{id}] Lcov_{id};",
      "  // cholesky factor of known covariance matrix\n"
    )
  }
  if (has_pw) {
    str_add(out$data) <- glue(
      "  vector[N_{id}] PW_{id};",
      "  // weights for group contribution to the prior\n"
    )
  }
  J <- seq_rows(r)
  reqZ <- !r$type %in% "sp"
  if (any(reqZ)) {
    str_add(out$data) <- "  // group-level predictor values\n"
    if (r$gtype[1] == "mm") {
      for (i in which(reqZ)) {
        str_add(out$data) <- cglue(
          "  vector[N{usc(r$resp[i])}] Z_{idp[i]}_{r$cn[i]}_{ng};\n"
        )
        str_add(out$pll_args) <- cglue(
          ", data vector Z_{idp[i]}_{r$cn[i]}_{ng}"
        )
      }
    } else {
      str_add(out$data) <- cglue(
        "  vector[N{usc(r$resp[reqZ])}] Z_{idp[reqZ]}_{r$cn[reqZ]};\n"
      )
      str_add(out$pll_args) <- cglue(
        ", data vector Z_{idp[reqZ]}_{r$cn[reqZ]}"
      )
    }
  }

  # define standard deviation parameters
  has_special_prior <- has_special_prior(prior, px, class = "sd")
  if (has_by) {
    if (has_special_prior) {
      stop2("Special priors on class 'sd' are not yet compatible ",
            "with the 'by' argument.")
    }
    str_add_list(out) <- stan_prior(
      prior, class = "sd", group = r$group[1], coef = r_cov$coef,
      type = glue("matrix[M_{id}, Nby_{id}]"),
      coef_type = glue("row_vector[Nby_{id}]"),
      suffix = glue("_{id}"), px = px_cov, broadcast = "matrix",
      comment = "group-level standard deviations",
      normalize = normalize
    )
  } else {
    if (has_special_prior) {
      if (stan_has_multiple_base_priors(px)) {
        stop2("Special priors on class 'sd' are not yet compatible with ",
              "group-level coefficients correlated across formulas.")
      }
      str_add(out$tpar_def) <- glue(
        "  vector<lower=0>[M_{id}] sd_{id};  // group-level standard deviations\n"
      )
    } else {
      str_add_list(out) <- stan_prior(
        prior, class = "sd", group = r$group[1], coef = r_cov$coef,
        type = glue("vector[M_{id}]"), suffix = glue("_{id}"), px = px_cov,
        comment = "group-level standard deviations",
        normalize = normalize
      )
    }
  }

  if (s2z_latent) {
    return(.stan_re_s2z_latent(
      id, r = r, bframe = bframe, prior = prior,
      normalize = normalize, out = out, fisher_info = fisher_info
    ))
  }

  if (s2z_cross) {
    return(.stan_re_s2z_cross(
      id, bframe = bframe, prior = prior,
      normalize = normalize, out = out,
      fisher_info = fisher_info, stanvars = stanvars, ...
    ))
  }

  if (isTRUE(r$s2z[1]) && !is.null(joint_s2z)) {
    return(.stan_re_s2z_joint_block(
      id, set = joint_s2z, bframe = bframe, prior = prior,
      threads = threads, normalize = normalize, out = out,
      fisher_info = fisher_info, stanvars = stanvars
    ))
  }

  if (isTRUE(r$s2z[1]) && identical(r$scale[1], "varying")) {
    return(.stan_re_s2z_varying_scale(
      id, bframe = bframe, prior = prior, threads = threads,
      normalize = normalize, out = out, fisher_info = fisher_info,
      stanvars = stanvars
    ))
  }

  if (isTRUE(r$s2z[1])) {
    return(.stan_re_s2z(
      id, bframe = bframe, prior = prior, threads = threads,
      normalize = normalize, out = out, fisher_info = fisher_info,
      stanvars = stanvars
    ))
  }

  if (!identical(center_mode, "noncentered")) {
    return(.stan_re_centered(
      id, r = r, bframe = bframe, prior = prior,
      normalize = normalize, out = out, fisher_info = fisher_info
    ))
  }

  # define group-level coefficients
  dfm <- ""
  tr <- subset_reframe_dist(r, "student")
  if (nrow(r) > 1L && r$cor[1]) {
    # multiple correlated group-level effects
    str_add(out$data) <- glue(
      "  int<lower=1> NC_{id};  // number of group-level correlations\n"
    )
    str_add(out$par) <- glue(
      "  matrix[M_{id}, N_{id}] z_{id};",
      "  // standardized group-level effects\n"
    )
    if (has_pw) {
      str_add(out$model_prior) <- glue(
        "  for (j in 1:N_{id}) {{\n",
        "    target += PW_{id}[j] * std_normal_{lpdf}(z_{id}[, j]);\n",
        "  }\n"
      )
    } else {
      str_add(out$model_prior) <- glue(
        "  target += std_normal_{lpdf}(to_vector(z_{id}));\n"
      )
    }

    if (has_rows(tr)) {
      dfm <- glue("rep_matrix(dfm_{tr$ggn[1]}, M_{id}) .* ")
    }
    if (has_by) {
      str_add_list(out) <- stan_prior(
        prior, class = "L", group = r$group[1], coef = Nby,
        type = glue("cholesky_factor_corr[M_{id}]"),
        coef_type = glue("cholesky_factor_corr[M_{id}]"),
        suffix = glue("_{id}"), dim = glue("[Nby_{id}]"),
        comment = "cholesky factor of correlation matrix",
        normalize = normalize
      )
      # separate definition from computation to support fixed parameters
      str_add(out$tpar_def) <- glue(
        "  matrix[N_{id}, M_{id}] r_{id};  // actual group-level effects\n"
      )
      if (has_cov) {
        str_add(out$fun) <- "  #include 'fun_scale_r_cor_by_cov.stan'\n"
        rdef <- glue(
          "scale_r_cor_by_cov(z_{id}, sd_{id}, L_{id}, Jby_{id}, Lcov_{id})"
        )
      } else {
        str_add(out$fun) <- "  #include 'fun_scale_r_cor_by.stan'\n"
        rdef  <- glue("scale_r_cor_by(z_{id}, sd_{id}, L_{id}, Jby_{id})")
      }
      str_add(out$tpar_comp) <- glue(
        "  // compute actual group-level effects\n",
        "  r_{id} = {dfm}{rdef};\n"
      )
      str_add(out$gen_def) <- cglue(
        "  // compute group-level correlations\n",
        "  corr_matrix[M_{id}] Cor_{id}_{Nby}",
        " = multiply_lower_tri_self_transpose(L_{id}[{Nby}]);\n",
        "  vector<lower=-1,upper=1>[NC_{id}] cor_{id}_{Nby};\n"
      )
      str_add(out$gen_comp) <- stan_cor_gen_comp(
        glue("cor_{id}_{Nby}"), glue("M_{id}")
      )
    } else {
      str_add_list(out) <- stan_prior(
        prior, class = "L", group = r$group[1], suffix = usc(id),
        type = glue("cholesky_factor_corr[M_{id}]"),
        comment = "cholesky factor of correlation matrix",
        normalize = normalize
      )
      if (has_cov) {
        str_add(out$fun) <- "  #include 'fun_scale_r_cor_cov.stan'\n"
        rdef <- glue("scale_r_cor_cov(z_{id}, sd_{id}, L_{id}, Lcov_{id})")
      } else {
        str_add(out$fun) <- "  #include 'fun_scale_r_cor.stan'\n"
        rdef <- glue("scale_r_cor(z_{id}, sd_{id}, L_{id})")
      }
      # separate definition from computation to support fixed parameters
      str_add(out$tpar_def) <- glue(
        "  matrix[N_{id}, M_{id}] r_{id};  // actual group-level effects\n"
      )
      str_add(out$tpar_comp) <- glue(
        "  // compute actual group-level effects\n",
        "  r_{id} = {dfm}{rdef};\n"
      )
      str_add(out$gen_def) <- glue(
        "  // compute group-level correlations\n",
        "  corr_matrix[M_{id}] Cor_{id}",
        " = multiply_lower_tri_self_transpose(L_{id});\n",
        "  vector<lower=-1,upper=1>[NC_{id}] cor_{id};\n"
      )
      str_add(out$gen_comp) <- stan_cor_gen_comp(
        cor = glue("cor_{id}"), ncol = glue("M_{id}")
      )
    }
    # separate definition from computation to support fixed parameters
    str_add(out$tpar_def) <-
      "  // using vectors speeds up indexing in loops\n"
    str_add(out$tpar_def) <- cglue(
      "  vector[N_{id}] r_{idp}_{r$cn};\n"
    )
    str_add(out$tpar_comp) <- cglue(
      "  r_{idp}_{r$cn} = r_{id}[, {J}];\n"
    )
    str_add(out$pll_args) <- cglue(
      ", vector r_{idp}_{r$cn}"
    )
  } else {
    # single or uncorrelated group-level effects
    str_add(out$par) <- glue(
      "  array[M_{id}] vector[N_{id}] z_{id};",
      "  // standardized group-level effects\n"
    )
    if (has_pw) {
      str_add(out$model_prior) <- glue(
        "  for (j in 1:N_{id}) {{\n",
        cglue("    target += PW_{id}[j] * std_normal_{lpdf}(z_{id}[{seq_rows(r)}, j]);\n"),
        "  }\n"
      )
    } else {
      str_add(out$model_prior) <- cglue(
        "  target += std_normal_{lpdf}(z_{id}[{seq_rows(r)}]);\n"
      )
    }

    Lcov <- str_if(has_cov, glue("Lcov_{id} * "))
    if (has_rows(tr)) {
      dfm <- glue("dfm_{tr$ggn[1]} .* ")
    }
    if (has_by) {
      # separate definition from computation to support fixed parameters
      str_add(out$tpar_def) <- cglue(
        "  vector[N_{id}] r_{idp}_{r$cn};  // actual group-level effects\n"
      )
      str_add(out$tpar_comp) <- cglue(
        "  r_{idp}_{r$cn} = {dfm}(transpose(sd_{id}[{J}, Jby_{id}])",
        " .* ({Lcov}z_{id}[{J}]));\n"
      )
    } else {
      # separate definition from computation to support fixed parameters
      str_add(out$tpar_def) <- cglue(
        "  vector[N_{id}] r_{idp}_{r$cn};  // actual group-level effects\n"
      )
      str_add(out$tpar_comp) <- cglue(
        "  r_{idp}_{r$cn} = {dfm}(sd_{id}[{J}] * ({Lcov}z_{id}[{J}]));\n"
      )
    }
    str_add(out$pll_args) <- cglue(
      ", vector r_{idp}_{r$cn}"
    )
  }
  out
}

# Add the population-level location corresponding to each ordinary group
# coefficient. A fitted design-basis map expresses every population column in
# the span of the varying-coefficient design, covering equivalent but
# differently coded bases (for example, treatment contrasts versus cell
# means). With brms's default centered fixed-effect design, Intercept is a
# centered-design parameter, so its physical value must first be recovered.
# For QR designs, recover the physical coefficient vector in transformed
# parameters; stan_fe() otherwise does so only in generated quantities.
.stan_re_center_mean <- function(out, id, r, bframe) {
  stopifnot(is.reframe(r), has_rows(r), is.anybrmsframe(bframe))
  frames <- Filter(function(x) {
    rx <- x$frame$re
    has_rows(rx) && id %in% rx$id
  }, all_bframel(bframe))
  stopifnot(length(frames) == 1L)
  bfl <- frames[[1L]]
  fe <- bfl$frame$fe
  px <- check_prefix(bfl)
  p <- usc(combine_prefix(px))
  M <- nrow(r)
  mean_name <- glue("mean_center_re_{id}")
  coefficient_name <- glue("b{p}")

  if (has_re_s2z_conventional(bfl)) {
    stop2("Slide-style centering for ordinary group effects cannot yet be ",
          "combined with a conventional sum-to-zero group effect in the ",
          "same linear predictor.")
  }
  if (identical(fe$decomp, "QR") && length(fe$vars_stan)) {
    ct <- str_if(fe$center, "c")
    coefficient_name <- glue("b_center_re_{id}")
    str_add(out$tpar_def) <- glue(
      "  vector[K{ct}{p}] {coefficient_name};",
      "  // physical population coefficients for group centering\n"
    )
    str_add(out$tpar_comp) <- glue(
      "  {coefficient_name} = XR{p}_inv * bQ{p};\n"
    )
  }

  beta_expressions <- rep("0.0", length(fe$vars))
  names(beta_expressions) <- fe$vars
  if (isTRUE(fe$center) && !is_ordinal(bfl$family) &&
      "Intercept" %in% fe$vars) {
    intercept <- glue("Intercept{p}")
    if (length(fe$vars_stan)) {
      intercept <- glue(
        "{intercept} - dot_product(means_X{p}, {coefficient_name})"
      )
    }
    beta_expressions["Intercept"] <- intercept
  }
  remaining <- which(beta_expressions == "0.0")
  matched <- match(fe$vars[remaining], fe$vars_stan)
  present <- !is.na(matched)
  if (any(present)) {
    beta_expressions[remaining[present]] <- glue(
      "{coefficient_name}[{matched[present]}]"
    )
  }

  C <- bfl$frame$re_center_mean[[as.character(id)]]
  if (is.null(C)) {
    # Compatibility fallback for frames produced before fitted design maps
    # were stored. Current frames always take the exact branch above.
    C <- matrix(
      0, nrow = M, ncol = length(fe$vars),
      dimnames = list(r$coef, fe$vars)
    )
    matched <- match(r$coef, fe$vars)
    present <- !is.na(matched)
    C[cbind(which(present), matched[present])] <- 1
  }
  stopifnot(
    is.matrix(C), identical(dim(C), c(M, length(fe$vars))),
    identical(rownames(C), r$coef), identical(colnames(C), fe$vars)
  )
  expressions <- vapply(seq_len(M), function(j) {
    take <- which(C[j, ] != 0 & beta_expressions != "0.0")
    if (!length(take)) {
      return("0.0")
    }
    terms <- vapply(take, function(k) {
      value <- C[j, k]
      expression <- beta_expressions[k]
      if (value == 1) {
        expression
      } else if (value == -1) {
        glue("-({expression})")
      } else {
        glue("{stan_s2z_number(value)} * ({expression})")
      }
    }, character(1))
    paste(terms, collapse = " + ")
  }, character(1))

  str_add(out$tpar_def) <- glue(
    "  vector[M_{id}] {mean_name};",
    "  // matching population locations for group centering\n"
  )
  str_add(out$tpar_comp) <- glue(
    "  {mean_name} = zeros_vector(M_{id});\n"
  )
  nonzero <- which(expressions != "0.0")
  if (length(nonzero)) {
    str_add(out$tpar_comp) <- cglue(
      "  {mean_name}[{nonzero}] = {expressions[nonzero]};\n"
    )
  }
  list(out = out, name = mean_name)
}

# Exact unrestricted centering chart for ordinary group effects. For one
# grouping level with conditional Cholesky factor L, population location mu,
# and centering fractions rho, define R = diag(rho) and
# A = R L + diag(1 - rho), then map the sampled slide coordinate q to the
# conventional deviation r = L A^-1 (q - R mu). The Jacobian is |L| / |A|.
# At rho = 1, q is the total group coefficient mu + r; at rho = 0 it is the
# legacy standardized coordinate. Unlike
# the S2Z chart, this full-space map needs neither a zero-sum projection nor a
# restricted-determinant correction. Student-t effects use their existing
# conditional scale-mixture factor in L, which preserves the legacy chart at
# rho = 0 exactly.
.stan_re_centered <- function(id, r, bframe, prior, normalize,
                              out = list(), fisher_info = NULL) {
  stopifnot(is.reframe(r), has_rows(r), !any(r$s2z))
  lpdf <- stan_lpdf_name(normalize)
  mode <- re_s2z_center_mode(r)
  stopifnot(mode %in% c("centered", "partial", "auto"))
  is_cor <- nrow(r) > 1L && isTRUE(r$cor[1L])
  is_student <- identical(r$dist[1L], "student")
  has_by <- nzchar(r$by[1L])
  has_cov <- nzchar(r$cov[1L])
  has_pw <- isTRUE(nzchar(r$gcall[[1L]]$pw))
  stopifnot(!has_by, !has_cov, !has_pw)
  M <- nrow(r)
  J <- seq_rows(r)
  px <- check_prefix(r)
  idp <- paste0(r$id, usc(combine_prefix(px)))
  rho <- re_s2z_center_values(r)
  fixed_center_data <- identical(mode, "partial") && re_center_has_data(r)
  tr <- subset_reframe_dist(r, "student")
  g <- if (is_student) usc(tr$ggn[1L]) else ""

  if (fixed_center_data) {
    str_add(out$data) <- glue(
      "  matrix<lower=0,upper=1>[N_{id}, M_{id}] rho_s2z_{id};",
      "  // fixed level-by-coefficient centering fractions\n"
    )
  }

  if (!is.null(fisher_info)) {
    stopifnot(identical(mode, "auto"))
    out <- stan_re_s2z_fisher_def(out, id)
    if (isTRUE(fisher_info$fixed_design)) {
      out <- stan_re_s2z_fisher_tdata(out, id, r, fisher_info)
    }
  }
  center_mean <- .stan_re_center_mean(out, id, r = r, bframe = bframe)
  out <- center_mean$out
  mean_name <- center_mean$name

  # Keep the established internal z_ID names so default save/exclude behavior
  # and downstream public r_ID naming remain unchanged.
  chart_comment <- if (identical(mode, "centered")) {
    "centered group-level coordinates"
  } else if (identical(mode, "auto")) {
    "precursor-selected fixed group-level coordinates"
  } else {
    "partially centered group-level coordinates"
  }

  if (is_cor) {
    str_add(out$data) <- glue(
      "  int<lower=1> NC_{id};  // number of group-level correlations\n"
    )
    str_add_list(out) <- stan_prior(
      prior, class = "L", group = r$group[1L], suffix = usc(id),
      type = glue("cholesky_factor_corr[M_{id}]"),
      comment = "cholesky factor of correlation matrix",
      normalize = normalize
    )
    str_add(out$par) <- glue(
      "  matrix[M_{id}, N_{id}] z_{id};  // {chart_comment}\n"
    )
    str_add(out$tpar_def) <- glue(
      "  matrix[M_{id}, M_{id}] L_center_re_{id};",
      "  // Gaussian-reference group Cholesky factor\n",
      "  matrix[N_{id}, M_{id}] r_{id};",
      "  // actual group-level effects\n",
      "  real log_jacobian_re_{id};",
      "  // unrestricted centering log-Jacobian\n"
    )
    str_add(out$tpar_comp) <- glue(
      "  L_center_re_{id} = diag_pre_multiply(sd_{id}, L_{id});\n"
    )
    if (!is.null(fisher_info)) {
      str_add(out$gen_comp) <- stan_re_s2z_fisher_gq_comp(
        id, r = r, fisher_info = fisher_info,
        L = glue("L_center_re_{id}")
      )
    }
    if (identical(mode, "centered")) {
      str_add(out$tpar_comp) <- glue(
        "  for (j in 1:N_{id}) {{\n",
        "    r_{id}[j] = (z_{id}[, j] - {mean_name})';\n",
        "  }}\n",
        "  log_jacobian_re_{id} = 0.0;\n"
      )
    } else {
      rho_def <- if (identical(mode, "auto") || fixed_center_data) {
        ""
      } else {
        rho_code <- stan_vector(vapply(rho, stan_s2z_number, character(1)))
        glue("    vector[M_{id}] rho_center_re = {rho_code};\n")
      }
      rho_j <- if (identical(mode, "auto") || fixed_center_data) {
        "rho_level_center_re"
      } else {
        "rho_center_re"
      }
      rho_level_def <- if (identical(mode, "auto") || fixed_center_data) {
        glue(
          "      vector[M_{id}] rho_level_center_re = ",
          "rho_s2z_{id}[j]';\n"
        )
      } else {
        ""
      }
      str_add(out$tpar_comp) <- glue(
        "  {{\n",
        "{rho_def}",
        "    log_jacobian_re_{id} = 0.0;\n",
        "    for (j in 1:N_{id}) {{\n",
        "      matrix[M_{id}, M_{id}] L_group_center_re = ",
        "L_center_re_{id};\n",
        "{rho_level_def}",
        str_if(
          is_student,
          glue("      L_group_center_re *= dfm{g}[j];\n")
        ),
        "      matrix[M_{id}, M_{id}] L_partial_center_re = ",
        "diag_pre_multiply({rho_j}, L_group_center_re);\n",
        "      vector[M_{id}] white_center_re;\n",
        "      for (k in 1:M_{id}) {{\n",
        "        L_partial_center_re[k, k] += 1.0 - {rho_j}[k];\n",
        "      }}\n",
        "      white_center_re = mdivide_left_tri_low(\n",
        "        L_partial_center_re, z_{id}[, j] - ",
        "{rho_j} .* {mean_name}\n",
        "      );\n",
        "      r_{id}[j] = (L_group_center_re * white_center_re)';\n",
        "      log_jacobian_re_{id} += ",
        "sum(log(diagonal(L_group_center_re))) - ",
        "sum(log(diagonal(L_partial_center_re)));\n",
        "    }}\n",
        "  }}\n"
      )
    }
    str_add(out$model_prior) <- glue(
      "  for (j in 1:N_{id}) {{\n",
      "    matrix[M_{id}, M_{id}] L_group_center_re = ",
      "L_center_re_{id};\n",
      str_if(
        is_student,
        glue("    L_group_center_re *= dfm{g}[j];\n")
      ),
      "    target += multi_normal_cholesky_{lpdf}(r_{id}[j]' | ",
      "zeros_vector(M_{id}), L_group_center_re);\n",
      "  }}\n",
      "  target += log_jacobian_re_{id};\n"
    )
    str_add(out$gen_def) <- glue(
      "  // compute group-level correlations\n",
      "  corr_matrix[M_{id}] Cor_{id}",
      " = multiply_lower_tri_self_transpose(L_{id});\n",
      "  vector<lower=-1,upper=1>[NC_{id}] cor_{id};\n"
    )
    str_add(out$gen_comp) <- stan_cor_gen_comp(
      cor = glue("cor_{id}"), ncol = glue("M_{id}")
    )
    str_add(out$tpar_def) <-
      "  // using vectors speeds up indexing in loops\n"
    str_add(out$tpar_def) <- cglue(
      "  vector[N_{id}] r_{idp}_{r$cn};\n"
    )
    str_add(out$tpar_comp) <- cglue(
      "  r_{idp}_{r$cn} = r_{id}[, {J}];\n"
    )
  } else {
    str_add(out$par) <- glue(
      "  array[M_{id}] vector[N_{id}] z_{id};",
      "  // {chart_comment}\n"
    )
    str_add(out$tpar_def) <- glue(
      "  real log_jacobian_re_{id};",
      "  // unrestricted centering log-Jacobian\n"
    )
    str_add(out$tpar_def) <- cglue(
      "  vector[N_{id}] r_{idp}_{r$cn};",
      "  // actual group-level effects\n"
    )
    if (!is.null(fisher_info)) {
      str_add(out$gen_comp) <- stan_re_s2z_fisher_gq_comp(
        id, r = r, fisher_info = fisher_info,
        scale = if (M == 1L) glue("sd_{id}[1]") else NULL,
        diag_scale = if (M > 1L) glue("sd_{id}") else NULL
      )
    }
    if (identical(mode, "centered")) {
      str_add(out$tpar_comp) <- cglue(
        "  r_{idp}_{r$cn} = z_{id}[{J}] - ",
        "rep_vector({mean_name}[{J}], N_{id});\n"
      )
      str_add(out$tpar_comp) <- glue(
        "  log_jacobian_re_{id} = 0.0;\n"
      )
    } else {
      str_add(out$tpar_comp) <- glue(
        "  log_jacobian_re_{id} = 0.0;\n"
      )
      for (k in seq_len(M)) {
        rho_k <- if (identical(mode, "auto") || fixed_center_data) {
          glue("rho_s2z_{id}[, {k}]")
        } else {
          glue("rep_vector({stan_s2z_number(rho[k])}, N_{id})")
        }
        scale_k <- if (is_student) {
          glue("sd_{id}[{k}] * dfm{g}")
        } else {
          glue("rep_vector(sd_{id}[{k}], N_{id})")
        }
        str_add(out$tpar_comp) <- glue(
          "  {{\n",
          "    vector[N_{id}] rho_group_center_re = {rho_k};\n",
          "    vector[N_{id}] scale_group_center_re = {scale_k};\n",
          "    vector[N_{id}] denominator_center_re = 1.0 - ",
          "rho_group_center_re + rho_group_center_re .* ",
          "scale_group_center_re;\n",
          "    r_{idp[k]}_{r$cn[k]} = scale_group_center_re .* ",
          "(z_{id}[{k}] - rho_group_center_re * {mean_name}[{k}]) ./ ",
          "denominator_center_re;\n",
          "    log_jacobian_re_{id} += sum(log(scale_group_center_re)) - ",
          "sum(log(denominator_center_re));\n",
          "  }}\n"
        )
      }
    }
    for (k in seq_len(M)) {
      scale_k <- if (is_student) {
        glue("sd_{id}[{k}] * dfm{g}")
      } else {
        glue("sd_{id}[{k}]")
      }
      str_add(out$model_prior) <- glue(
        "  target += normal_{lpdf}(r_{idp[k]}_{r$cn[k]} | 0, ",
        "{scale_k});\n"
      )
    }
    str_add(out$model_prior) <- glue(
      "  target += log_jacobian_re_{id};\n"
    )
  }
  str_add(out$pll_args) <- cglue(
    ", vector r_{idp}_{r$cn}"
  )
  out
}

# Construct a link-scale reference predictor that is independent of every S2Z
# contrast in its linear predictor. Precursor proposals may depend on this
# global finite-population coordinate without changing the final fixed chart.
# Keeping this first implementation to dense population terms and fixed
# offsets also makes that independence auditable in generated Stan code.
stan_re_s2z_fisher_reference_eta <- function(bfl, n = "n") {
  stopifnot(is.bframel(bfl))
  unsupported <- function(detail) {
    stop2("Automatic centering for group-level effects ", detail, ".")
  }
  special <- c("cs", "sm", "sp", "gp")
  if (any(vapply(special, function(x) has_rows(bfl$frame[[x]]), logical(1)))) {
    unsupported(paste0(
      "currently requires a population-only reference predictor without ",
      "category-specific, smooth, spatial, or Gaussian-process terms"
    ))
  }
  fe <- bfl$frame$fe
  if (isTRUE(fe$sparse) || !identical(fe$decomp, "none")) {
    unsupported("currently requires a dense, non-QR population design")
  }
  p <- usc(combine_prefix(check_prefix(bfl)))
  K <- length(fe$vars_stan)
  s2z <- has_re_s2z_conventional(bfl)
  eta <- if (isTRUE(fe$center)) {
    intercept <- if (s2z) glue("theta_s2z{p}[1]") else glue("Intercept{p}")
    if (K) {
      coefficient <- if (s2z) {
        glue("tail(theta_s2z{p}, {K})")
      } else {
        glue("b{p}")
      }
      glue("{intercept} + dot_product(Xc{p}[{n}], {coefficient})")
    } else {
      intercept
    }
  } else {
    if (K) {
      coefficient <- if (s2z) glue("theta_s2z{p}") else glue("b{p}")
      glue("dot_product(X{p}[{n}], {coefficient})")
    } else {
      "0.0"
    }
  }
  if (is.formula(bfl$offset)) {
    eta <- glue("({eta}) + offsets{p}[{n}]")
  }
  eta
}

# Add the realized group-level contribution to the population reference at a
# precursor draw. This is safe because the resulting information is evaluated
# only in generated quantities and is frozen before the final fit. In
# particular, no derivative of the sampled chart is taken through this value.
stan_re_s2z_fisher_fitted_eta <- function(bfl, n = "n") {
  eta <- stan_re_s2z_fisher_reference_eta(bfl, n = n)
  group_eta <- stan_eta_re(bfl, threads = NULL)
  if (nzchar(group_eta)) {
    eta <- glue("({eta}{group_eta})")
  }
  eta
}

# Construct drawwise expansion points for every predicted distributional
# parameter. Likelihood curvature for one parameter may depend on the fitted
# values of the others, so all linear distributional predictors are supplied.
stan_re_s2z_fisher_fitted_references <- function(bframe, n = "n") {
  stopifnot(is.brmsframe(bframe))
  predicted <- vapply(bframe$dpars, is.bframel, logical(1))
  out <- lapply(
    bframe$dpars[predicted], stan_re_s2z_fisher_fitted_eta, n = n
  )
  names(out) <- names(bframe$dpars)[predicted]
  out
}

# Natural-scale reference and derivative for one distributional parameter.
# By default a reference excludes group-level contrasts. The automatic
# precursor overrides linear predictors with their fitted values; this is
# valid because the candidate is generated after sampling and then frozen.
stan_re_s2z_fisher_dpar_reference <- function(
    bframe, dpar, n = "n", eta = NULL
) {
  stopifnot(is.brmsframe(bframe), length(dpar) == 1L)
  bfl <- bframe$dpars[[dpar]]
  if (!is.null(eta)) {
    eta <- as_one_character(eta)
    if (!is.bframel(bfl) && !is.bframenl(bfl)) {
      stop2("An automatic-centering reference override requires a predicted ",
            "distributional parameter '", dpar, "'.")
    }
  } else if (is.bframenl(bfl)) {
    stop2("Automatic centering requires an explicit score-free reference for ",
          "nonlinear distributional parameter '", dpar, "'.")
  } else if (!is.bframel(bfl)) {
    # Unpredicted auxiliary parameters are scalar Stan parameters (or fixed
    # scalar transformed parameters) already on their natural scale.
    value <- paste0(dpar, usc(bframe$resp))
    return(nlist(
      dpar, bfl = NULL, eta = value, link = "identity", value,
      derivative = "1.0"
    ))
  } else {
    eta <- stan_re_s2z_fisher_reference_eta(bfl, n = n)
  }
  link <- bfl$family$link
  inv_link <- stan_inv_link(link, transform = TRUE)
  value <- if (nzchar(inv_link)) glue("{inv_link}({eta})") else glue("({eta})")
  derivative <- switch(
    link,
    identity = "1.0",
    log = value,
    logm1 = glue("exp({eta})"),
    log1p = glue("exp({eta})"),
    inverse = glue("-inv_square({eta})"),
    sqrt = glue("2.0 * ({eta})"),
    "1/mu^2" = glue("-0.5 * inv({eta}) * inv_sqrt({eta})"),
    logit = glue("{value} * (1.0 - {value})"),
    probit = glue("exp(std_normal_lpdf({eta}))"),
    probit_approx = glue(
      "{value} * (1.0 - {value}) * ",
      "(1.5976 + 0.21168 * square({eta}))"
    ),
    cloglog = glue("exp(({eta}) - exp({eta}))"),
    cauchit = glue("inv(pi() * (1.0 + square({eta})))"),
    softplus = glue("inv_logit({eta})"),
    squareplus = glue(
      "0.5 * (1.0 + ({eta}) / sqrt(square({eta}) + 4.0))"
    ),
    softit = glue(
      "inv_logit({eta}) / square(1.0 + log1p_exp({eta}))"
    ),
    tan_half = glue("2.0 / (1.0 + square({eta}))"),
    stop2("Automatic centering does not yet implement the derivative of link '",
          link, "'.")
  )
  nlist(dpar, bfl, eta, link, value, derivative)
}

# Response-free conditional expected information for a scalar distributional
# predictor. Exact identities are preferred; a small number of families use a
# positive analytic coarsening or moment approximation when their exact
# expectation would require a sum or integral. The result is already on the
# predictor (eta) scale for evaluation in precursor generated quantities.
stan_re_s2z_fisher_closed_form <- function(
    bframe, dpar, n = "n", reference_eta = NULL
) {
  stopifnot(is.brmsframe(bframe), length(dpar) == 1L)
  if (is.null(reference_eta)) {
    reference_eta <- list()
  }
  if (is.character(reference_eta)) {
    reference_eta <- as.list(reference_eta)
  }
  if (!is.list(reference_eta) ||
      length(reference_eta) && is.null(names(reference_eta)) ||
      any(!nzchar(names(reference_eta)))) {
    stop2("Argument 'reference_eta' must be a named list or character vector.")
  }
  family <- bframe$family$family
  references <- new.env(parent = emptyenv())
  definitions <- character()
  ref <- function(dp) {
    if (exists(dp, envir = references, inherits = FALSE)) {
      return(get(dp, envir = references, inherits = FALSE))
    }
    raw <- stan_re_s2z_fisher_dpar_reference(
      bframe, dp, n = n, eta = reference_eta[[dp]]
    )
    key <- make_stan_names(dp)
    eta_name <- paste0("eta_fisher_s2z_", key)
    value_name <- paste0("value_fisher_s2z_", key)
    derivative_name <- paste0("derivative_fisher_s2z_", key)
    value_expression <- gsub(
      raw$eta, eta_name, raw$value, fixed = TRUE
    )
    derivative_expression <- gsub(
      raw$value, value_name, raw$derivative, fixed = TRUE
    )
    derivative_expression <- gsub(
      raw$eta, eta_name, derivative_expression, fixed = TRUE
    )
    definitions <<- c(definitions, glue(
      "        real {eta_name} = {raw$eta};\n",
      "        real {value_name} = {value_expression};\n",
      "        real {derivative_name} = {derivative_expression};\n"
    ))
    out <- raw
    out$eta <- eta_name
    out$value <- value_name
    out$derivative <- derivative_name
    assign(dp, out, envir = references)
    out
  }
  target <- ref(dpar)
  sq_chain <- glue("square({target$derivative})")
  resp <- usc(bframe$resp)
  rate <- is.formula(bframe$adforms$rate)
  denom <- if (rate) glue("denom{resp}[{n}]") else "1.0"
  scaled <- function(x) if (rate) glue("({x}) * {denom}") else x
  scaled_derivative <- function(x) {
    if (rate) glue("({x}) * {denom}") else x
  }
  # Stable Bernoulli information on the predictor scale. In particular, avoid
  # the 0 / 0 form produced when an inverse-logit rounds to exactly zero or one
  # in a far-tail proposal.
  binary_information <- function(x) {
    if (identical(x$link, "logit")) {
      return(glue("({x$value}) * (1.0 - ({x$value}))"))
    }
    if (identical(x$link, "probit")) {
      return(glue(
        "exp(2.0 * std_normal_lpdf({x$eta}) - ",
        "std_normal_lcdf({x$eta}) - std_normal_lccdf({x$eta}))"
      ))
    }
    if (identical(x$link, "probit_approx")) {
      return(glue(
        "({x$value}) * (1.0 - ({x$value})) * square(",
        "1.5976 + 0.21168 * square({x$eta}))"
      ))
    }
    glue(
      "square({x$derivative}) / fmax(({x$value}) * ",
      "(1.0 - ({x$value})), 1e-12)"
    )
  }
  quantile_information <- function(x) {
    if (identical(x$link, "logit")) {
      return(glue(
        "square({x$value}) + square(1.0 - ({x$value}))"
      ))
    }
    glue(
      "square({x$derivative}) * (square({x$value}) + ",
      "square(1.0 - ({x$value}))) / fmax(square({x$value}) * ",
      "square(1.0 - ({x$value})), 1e-24)"
    )
  }
  # x^2 trigamma(x) tends to one as x tends to zero. This scaled form avoids
  # the indeterminate products arising in beta and Dirichlet information when
  # a probability predictor saturates numerically.
  scaled_trigamma <- function(x) {
    glue("(({x}) < 1e-6 ? 1.0 : square({x}) * trigamma({x}))")
  }

  value <- target$value
  derivative <- target$derivative
  out <- NULL
  if (family == "gaussian") {
    sigma <- ref("sigma")$value
    out <- switch(dpar,
      mu = glue("{sq_chain} * inv_square({sigma})"),
      sigma = glue("2.0 * {sq_chain} * inv_square({value})"),
      NULL
    )
  } else if (family == "student") {
    sigma <- ref("sigma")$value
    nu <- ref("nu")$value
    out <- switch(dpar,
      mu = glue(
        "{sq_chain} * ({nu} + 1.0) / ({nu} + 3.0) * ",
        "inv_square({sigma})"
      ),
      sigma = glue(
        "{sq_chain} * 2.0 * {nu} / ({nu} + 3.0) * ",
        "inv_square({value})"
      ),
      nu = glue(
        "{sq_chain} * ({nu} > 1e4 ? 3.5 / pow({nu}, 4) : ",
        "fmax(0.0, 0.25 * (trigamma(0.5 * {nu}) - ",
        "trigamma(0.5 * ({nu} + 1.0))) - ({nu} + 5.0) / ",
        "(2.0 * {nu} * ({nu} + 1.0) * ({nu} + 3.0))))"
      ),
      NULL
    )
  } else if (family == "skew_normal") {
    sigma <- ref("sigma")$value
    alpha <- ref("alpha")$value
    if (dpar == "alpha") {
      definitions <- c(definitions, glue(
        "        real delta_fisher_s2z_skew = ({alpha}) / ",
        "sqrt(1.0 + square({alpha}));\n",
        "        real b_fisher_s2z_skew = sqrt(2.0 / pi()) * ",
        "delta_fisher_s2z_skew;\n",
        "        real variance_factor_fisher_s2z_skew = 1.0 - ",
        "square(b_fisher_s2z_skew);\n",
        "        real db_fisher_s2z_skew = sqrt(2.0 / pi()) * ",
        "pow(1.0 + square({alpha}), -1.5);\n",
        "        real dskewness_fisher_s2z_skew = 1.5 * ",
        "(4.0 - pi()) * square(b_fisher_s2z_skew) * ",
        "db_fisher_s2z_skew * ",
        "pow(variance_factor_fisher_s2z_skew, -2.5);\n",
        "        real dkurtosis_fisher_s2z_skew = 8.0 * ",
        "(pi() - 3.0) * pow(b_fisher_s2z_skew, 3) * ",
        "db_fisher_s2z_skew * ",
        "pow(variance_factor_fisher_s2z_skew, -3.0);\n"
      ))
    }
    out <- switch(dpar,
      mu = glue("{sq_chain} * inv_square({sigma})"),
      sigma = glue("2.0 * {sq_chain} * inv_square({value})"),
      alpha = glue(
        "{sq_chain} * (square(dskewness_fisher_s2z_skew) / 6.0 + ",
        "square(dkurtosis_fisher_s2z_skew) / 24.0)"
      ),
      NULL
    )
  } else if (family == "logistic_normal") {
    all_dpars <- valid_dpars(bframe)
    mu_dpars <- all_dpars[vapply(
      all_dpars,
      function(dp) identical(dpar_class(dp, bframe$family), "mu"),
      logical(1)
    )]
    sigma_dpars <- all_dpars[vapply(
      all_dpars,
      function(dp) identical(dpar_class(dp, bframe$family), "sigma"),
      logical(1)
    )]
    Klogit <- length(mu_dpars)
    if (Klogit && length(sigma_dpars) == Klogit) {
      definitions <- c(definitions, glue(
        "        matrix[{Klogit}, {Klogit}] corr_precision_fisher_s2z_ln = ",
        "chol2inv(Llncor{resp});\n"
      ))
      if (dpar %in% mu_dpars) {
        k <- match(dpar, mu_dpars)
        sigma_k <- ref(sigma_dpars[k])$value
        out <- glue(
          "corr_precision_fisher_s2z_ln[{k}, {k}] * ",
          "inv_square({sigma_k})"
        )
      } else if (dpar %in% sigma_dpars) {
        k <- match(dpar, sigma_dpars)
        out <- glue(
          "square(({derivative}) / ({value})) * ",
          "(1.0 + corr_precision_fisher_s2z_ln[{k}, {k}])"
        )
      }
    }
  } else if (family %in% c("lognormal", "shifted_lognormal")) {
    sigma <- ref("sigma")$value
    out <- switch(dpar,
      mu = glue("{sq_chain} * inv_square({sigma})"),
      sigma = glue("2.0 * {sq_chain} * inv_square({value})"),
      ndt = str_if(
        family == "shifted_lognormal",
        glue(
          "{sq_chain} * exp(-2.0 * {ref('mu')$value} + ",
          "2.0 * square({sigma})) * (1.0 + inv_square({sigma}))"
        )
      ),
      NULL
    )
  } else if (family %in% c("bernoulli", "binomial")) {
    if (dpar == "mu") {
      trials <- if (family == "binomial") {
        glue("trials{resp}[{n}]")
      } else {
        "1.0"
      }
      out <- glue("{trials} * ({binary_information(target)})")
    }
  } else if (family %in% c("categorical", "multinomial")) {
    mu_dpars <- names(bframe$dpars)[vapply(
      names(bframe$dpars),
      function(dp) identical(dpar_class(dp, bframe$family), "mu"),
      logical(1)
    )]
    refcat <- get_refcat(bframe$family, int = TRUE)
    if (isNA(refcat)) {
      Kcat <- length(mu_dpars)
      positions <- seq_len(Kcat)
    } else {
      Kcat <- length(mu_dpars) + 1L
      positions <- setdiff(seq_len(Kcat), refcat)
    }
    target_position <- positions[match(dpar, mu_dpars)]
    if (!is.na(target_position)) {
      raw_eta <- vapply(mu_dpars, function(dp) {
        ref(dp)$eta
      }, character(1))
      assignments <- paste0(vapply(seq_along(mu_dpars), function(k) {
        glue(
          "        eta_fisher_s2z_cat[{positions[k]}] = {raw_eta[k]};\n"
        )
      }, character(1)), collapse = "")
      definitions <- c(definitions, glue(
        "        vector[{Kcat}] eta_fisher_s2z_cat = zeros_vector({Kcat});\n",
        "{assignments}",
        "        vector[{Kcat}] prob_fisher_s2z_cat = ",
        "softmax(eta_fisher_s2z_cat);\n"
      ))
      trials <- if (family == "multinomial") {
        glue("trials{resp}[{n}]")
      } else {
        "1.0"
      }
      out <- glue(
        "{trials} * prob_fisher_s2z_cat[{target_position}] * ",
        "(1.0 - prob_fisher_s2z_cat[{target_position}])"
      )
    }
  } else if (family %in% c("dirichlet", "dirichlet_multinomial")) {
    mu_dpars <- names(bframe$dpars)[vapply(
      names(bframe$dpars),
      function(dp) identical(dpar_class(dp, bframe$family), "mu"),
      logical(1)
    )]
    phi <- ref("phi")$value
    refcat <- get_refcat(bframe$family, int = TRUE)
    if (isNA(refcat)) {
      Kcat <- length(mu_dpars)
      positions <- seq_len(Kcat)
    } else {
      Kcat <- length(mu_dpars) + 1L
      positions <- setdiff(seq_len(Kcat), refcat)
    }
    raw_eta <- vapply(mu_dpars, function(dp) {
      ref(dp)$eta
    }, character(1))
    assignments <- paste0(vapply(seq_along(mu_dpars), function(k) {
      glue(
        "        eta_fisher_s2z_dir[{positions[k]}] = {raw_eta[k]};\n"
      )
    }, character(1)), collapse = "")
    definitions <- c(definitions, glue(
      "        vector[{Kcat}] eta_fisher_s2z_dir = zeros_vector({Kcat});\n",
      "{assignments}",
      "        vector[{Kcat}] prob_fisher_s2z_dir = ",
      "softmax(eta_fisher_s2z_dir);\n"
    ))
    target_position <- positions[match(dpar, mu_dpars)]
    if (!is.na(target_position) && family == "dirichlet") {
      definitions <- c(definitions, glue(
        "        real info_fisher_s2z_dir = 0.0;\n",
        "        for (c in 1:{Kcat}) {{\n",
        "          real alpha_fisher_s2z_dir = ({phi}) * ",
        "prob_fisher_s2z_dir[c];\n",
        "          real score_factor_fisher_s2z_dir = ",
        "(c == {target_position}) - ",
        "prob_fisher_s2z_dir[{target_position}];\n",
        "          info_fisher_s2z_dir += square(",
        "score_factor_fisher_s2z_dir) * ",
        "(alpha_fisher_s2z_dir < 1e-6 ? 1.0 : ",
        "square(alpha_fisher_s2z_dir) * ",
        "trigamma(alpha_fisher_s2z_dir));\n",
        "        }}\n"
      ))
      out <- "info_fisher_s2z_dir"
    } else if (dpar == "phi" && family == "dirichlet") {
      definitions <- c(definitions, glue(
        "        real info_fisher_s2z_dir = -trigamma({phi});\n",
        "        for (c in 1:{Kcat}) {{\n",
        "          real alpha_fisher_s2z_dir = ({phi}) * ",
        "prob_fisher_s2z_dir[c];\n",
        "          info_fisher_s2z_dir += ",
        "(alpha_fisher_s2z_dir < 1e-6 ? 1.0 : ",
        "square(alpha_fisher_s2z_dir) * ",
        "trigamma(alpha_fisher_s2z_dir)) * inv_square({phi});\n",
        "        }}\n"
      ))
      out <- glue("{sq_chain} * fmax(0.0, info_fisher_s2z_dir)")
    } else if (family == "dirichlet_multinomial") {
      trials <- glue("trials{resp}[{n}]")
      effective_trials <- glue(
        "(({trials}) * (1.0 + ({phi})) / (({trials}) + ({phi})))"
      )
      if (!is.na(target_position)) {
        out <- glue(
          "{effective_trials} * prob_fisher_s2z_dir[{target_position}] * ",
          "(1.0 - prob_fisher_s2z_dir[{target_position}])"
        )
      } else if (dpar == "phi") {
        out <- glue(
          "({trials} == 0 ? 0.0 : 0.5 * ({Kcat} - 1.0) * ",
          "{sq_chain} * square(inv(({trials}) + ({phi})) - ",
          "inv(1.0 + ({phi}))))"
        )
      }
    }
  } else if (family == "dirichlet2") {
    mu_dpars <- names(bframe$dpars)[vapply(
      names(bframe$dpars),
      function(dp) identical(dpar_class(dp, bframe$family), "mu"),
      logical(1)
    )]
    alpha_refs <- lapply(mu_dpars, function(dp) ref(dp))
    alpha_values <- vapply(alpha_refs, `[[`, character(1), "value")
    alpha_sum <- paste(alpha_values, collapse = " + ")
    if (dpar %in% mu_dpars) {
      own_info <- if (identical(target$link, "log")) {
        scaled_trigamma(value)
      } else {
        glue(
          "square(({derivative}) / fmax({value}, 1e-12)) * ",
          "{scaled_trigamma(value)}"
        )
      }
      out <- glue(
        "fmax(0.0, ({own_info}) - {sq_chain} * trigamma({alpha_sum}))"
      )
    }
  } else if (family == "beta_binomial") {
    mu <- ref("mu")$value
    dmu <- ref("mu")$derivative
    phi <- ref("phi")$value
    dphi <- ref("phi")$derivative
    trials <- glue("trials{resp}[{n}]")
    definitions <- c(definitions, glue(
      "        real prob_fisher_s2z_bb = fmin(1.0 - 1e-12, ",
      "fmax(1e-12, {mu}));\n",
      "        real alpha_fisher_s2z_bb = prob_fisher_s2z_bb * ({phi});\n",
      "        real beta_fisher_s2z_bb = ",
      "(1.0 - prob_fisher_s2z_bb) * ({phi});\n",
      "        real log_p0_fisher_s2z_bb = ",
      "lbeta(alpha_fisher_s2z_bb, beta_fisher_s2z_bb + ({trials})) - ",
      "lbeta(alpha_fisher_s2z_bb, beta_fisher_s2z_bb);\n",
      "        real log_pn_fisher_s2z_bb = ",
      "lbeta(alpha_fisher_s2z_bb + ({trials}), beta_fisher_s2z_bb) - ",
      "lbeta(alpha_fisher_s2z_bb, beta_fisher_s2z_bb);\n",
      "        real p0_fisher_s2z_bb = exp(log_p0_fisher_s2z_bb);\n",
      "        real pn_fisher_s2z_bb = exp(log_pn_fisher_s2z_bb);\n",
      "        real pmid_fisher_s2z_bb = fmax(0.0, ",
      "1.0 - p0_fisher_s2z_bb - pn_fisher_s2z_bb);\n"
    ))
    if (dpar == "mu") {
      definitions <- c(definitions, glue(
        "        real score0_fisher_s2z_bb = ({dmu}) * ({phi}) * (",
        "digamma(beta_fisher_s2z_bb) - ",
        "digamma(beta_fisher_s2z_bb + ({trials})));\n",
        "        real scoren_fisher_s2z_bb = ({dmu}) * ({phi}) * (",
        "digamma(alpha_fisher_s2z_bb + ({trials})) - ",
        "digamma(alpha_fisher_s2z_bb));\n"
      ))
    } else if (dpar == "phi") {
      definitions <- c(definitions, glue(
        "        real score0_fisher_s2z_bb = ({dphi}) * (",
        "digamma({phi}) - digamma(({phi}) + ({trials})) + ",
        "(1.0 - prob_fisher_s2z_bb) * (",
        "digamma(beta_fisher_s2z_bb + ({trials})) - ",
        "digamma(beta_fisher_s2z_bb)));\n",
        "        real scoren_fisher_s2z_bb = ({dphi}) * (",
        "digamma({phi}) - digamma(({phi}) + ({trials})) + ",
        "prob_fisher_s2z_bb * (",
        "digamma(alpha_fisher_s2z_bb + ({trials})) - ",
        "digamma(alpha_fisher_s2z_bb)));\n"
      ))
    }
    if (dpar %in% c("mu", "phi")) {
      definitions <- c(definitions, glue(
        "        real dpmid_fisher_s2z_bb = -p0_fisher_s2z_bb * ",
        "score0_fisher_s2z_bb - pn_fisher_s2z_bb * ",
        "scoren_fisher_s2z_bb;\n"
      ))
      coarse_info <- glue(
        "p0_fisher_s2z_bb * square(score0_fisher_s2z_bb) + ",
        "pn_fisher_s2z_bb * square(scoren_fisher_s2z_bb) + ",
        "square(dpmid_fisher_s2z_bb) / ",
        "fmax(pmid_fisher_s2z_bb, 1e-12)"
      )
      out <- if (dpar == "mu") {
        glue(
          "({trials} == 0 ? 0.0 : ({trials} == 1 ? ",
          "({binary_information(ref('mu'))}) : ({coarse_info})))"
        )
      } else {
        glue("({trials} <= 1 ? 0.0 : ({coarse_info}))")
      }
    }
  } else if (family == "poisson") {
    if (dpar == "mu") {
      lambda <- scaled(value)
      dlambda <- scaled_derivative(derivative)
      out <- if (identical(target$link, "log")) {
        lambda
      } else {
        glue("square({dlambda}) / fmax({lambda}, 1e-12)")
      }
    }
  } else if (family %in% c("negbinomial", "negbinomial2", "geometric")) {
    if (dpar == "mu") {
      mu <- scaled(value)
      dmu <- scaled_derivative(derivative)
      shape <- if (family == "geometric") {
        scaled("1.0")
      } else if (family == "negbinomial") {
        scaled(ref("shape")$value)
      } else {
        scaled(glue("inv({ref('sigma')$value})"))
      }
      out <- if (identical(target$link, "log")) {
        glue("({mu}) * {shape} / (({mu}) + {shape})")
      } else {
        glue(
          "square({dmu}) * {shape} / ",
          "fmax(({mu}) * (({mu}) + {shape}), 1e-12)"
        )
      }
    } else if (family %in% c("negbinomial", "negbinomial2")) {
      # The exact expected information for log(shape) contains
      # E[trigamma(Y + shape)].  Avoid evaluating that expectation in Stan by
      # using a positive closed-form moment approximation.  It is exact to
      # leading order as mu -> 0, includes the exact information in the event
      # Y == 0, and approaches the Gaussian variance-information limit as the
      # zero probability vanishes.
      mu <- scaled(ref("mu")$value)
      if (family == "negbinomial") {
        shape <- scaled(value)
        log_shape_chain <- glue("({derivative}) / ({value})")
      } else {
        shape <- scaled(glue("inv({value})"))
        log_shape_chain <- glue("-({derivative}) / ({value})")
      }
      definitions <- c(definitions, glue(
        "        real mean_fraction_fisher_s2z_nb = ",
        "{mu} / ({mu} + {shape});\n",
        "        real log_p0_fisher_s2z_nb = -({shape}) * ",
        "log1p(({mu}) / ({shape}));\n",
        "        real p0_fisher_s2z_nb = exp(log_p0_fisher_s2z_nb);\n",
        "        real ppos_fisher_s2z_nb = ",
        "-expm1(log_p0_fisher_s2z_nb);\n",
        "        real dlog_p0_fisher_s2z_nb = ({shape}) * ",
        "(-log1p(({mu}) / ({shape})) + ",
        "mean_fraction_fisher_s2z_nb);\n",
        "        real sparse_info_fisher_s2z_nb = 0.5 * ({shape}) / ",
        "(({shape}) + 1.0) * square(mean_fraction_fisher_s2z_nb);\n",
        "        real dense_info_fisher_s2z_nb = 0.5 * ",
        "square(mean_fraction_fisher_s2z_nb);\n",
        "        real zero_info_fisher_s2z_nb = p0_fisher_s2z_nb * ",
        "square(dlog_p0_fisher_s2z_nb) / ",
        "fmax(ppos_fisher_s2z_nb, 1e-12);\n",
        "        real log_shape_info_fisher_s2z_nb = ",
        "sparse_info_fisher_s2z_nb + zero_info_fisher_s2z_nb + ",
        "ppos_fisher_s2z_nb * (dense_info_fisher_s2z_nb - ",
        "sparse_info_fisher_s2z_nb);\n"
      ))
      out <- glue(
        "square({log_shape_chain}) * log_shape_info_fisher_s2z_nb"
      )
    }
  } else if (family == "discrete_weibull") {
    mu <- ref("mu")$value
    shape <- ref("shape")$value
    definitions <- c(definitions, glue(
      "        real prob_fisher_s2z_dw = fmin(1.0 - 1e-12, ",
      "fmax(1e-12, {mu}));\n",
      "        real tail_power_fisher_s2z_dw = pow(2.0, {shape});\n",
      "        real log_tail_prob_fisher_s2z_dw = ",
      "tail_power_fisher_s2z_dw * log(prob_fisher_s2z_dw);\n",
      "        real p2_fisher_s2z_dw = exp(log_tail_prob_fisher_s2z_dw);\n",
      "        real p1_fisher_s2z_dw = ",
      "fmax(0.0, prob_fisher_s2z_dw - p2_fisher_s2z_dw);\n",
      "        real p0_fisher_s2z_dw = 1.0 - prob_fisher_s2z_dw;\n"
    ))
    out <- switch(dpar,
      mu = glue(
        "{sq_chain} * (inv(fmax(p0_fisher_s2z_dw, 1e-12)) + ",
        "square(1.0 - tail_power_fisher_s2z_dw * p2_fisher_s2z_dw / ",
        "prob_fisher_s2z_dw) / fmax(p1_fisher_s2z_dw, 1e-12) + ",
        "square(tail_power_fisher_s2z_dw * p2_fisher_s2z_dw / ",
        "prob_fisher_s2z_dw) / fmax(p2_fisher_s2z_dw, 1e-12))"
      ),
      shape = glue(
        "{sq_chain} * square(p2_fisher_s2z_dw * ",
        "tail_power_fisher_s2z_dw * log(prob_fisher_s2z_dw) * ",
        "log(2.0)) * ",
        "(inv(fmax(p1_fisher_s2z_dw, 1e-12)) + ",
        "inv(fmax(p2_fisher_s2z_dw, 1e-12)))"
      ),
      NULL
    )
  } else if (family == "com_poisson") {
    # The generated Stan likelihood is the standard COM-Poisson exponential
    # family p(y) proportional to mu^y / (y!)^shape.  Its exact information is
    # the covariance of (Y, log(Y!)) and would require differentiating the
    # log-normalizer.  Use its standard asymptotic variance together with a
    # positive second-order delta approximation for Var(log(Y!)).  This is
    # response-free, reduces exactly to Poisson location information at
    # shape = 1, and adds no normalizer evaluation to the proposal.
    mu <- ref("mu")$value
    dmu <- ref("mu")$derivative
    shape <- ref("shape")$value
    dshape <- ref("shape")$derivative
    definitions <- c(definitions, glue(
      "        real log_mode_fisher_s2z_cmp = log(fmax({mu}, 1e-12)) / ",
      "fmax({shape}, 1e-12);\n",
      "        real mode_fisher_s2z_cmp = exp(fmin(27.6310211159, ",
      "fmax(-27.6310211159, log_mode_fisher_s2z_cmp)));\n",
      "        real variance_fisher_s2z_cmp = fmin(1e12, ",
      "mode_fisher_s2z_cmp / fmax({shape}, 1e-12));\n",
      "        real log_factorial_slope_fisher_s2z_cmp = ",
      "digamma(mode_fisher_s2z_cmp + 1.0);\n",
      "        real log_factorial_curve_fisher_s2z_cmp = ",
      "trigamma(mode_fisher_s2z_cmp + 1.0);\n",
      "        real log_factorial_variance_fisher_s2z_cmp = ",
      "square(log_factorial_slope_fisher_s2z_cmp) * ",
      "variance_fisher_s2z_cmp + 0.5 * ",
      "square(log_factorial_curve_fisher_s2z_cmp) * ",
      "square(variance_fisher_s2z_cmp);\n"
    ))
    out <- switch(dpar,
      mu = if (identical(ref("mu")$link, "log")) {
        "variance_fisher_s2z_cmp"
      } else {
        glue(
          "square(({dmu}) / fmax({mu}, 1e-12)) * ",
          "variance_fisher_s2z_cmp"
        )
      },
      shape = glue(
        "square({dshape}) * log_factorial_variance_fisher_s2z_cmp"
      ),
      NULL
    )
  } else if (family == "gamma") {
    mu <- ref("mu")$value
    shape <- ref("shape")$value
    out <- switch(dpar,
      mu = glue("{sq_chain} * {shape} * inv_square({mu})"),
      shape = glue(
        "{sq_chain} * fmax(0.0, trigamma({shape}) - inv({shape}))"
      ),
      NULL
    )
  } else if (family == "exponential") {
    if (dpar == "mu") {
      out <- glue("{sq_chain} * inv_square({value})")
    }
  } else if (family == "inverse.gaussian") {
    mu <- ref("mu")$value
    shape <- ref("shape")$value
    out <- switch(dpar,
      mu = glue("{sq_chain} * {shape} / pow({mu}, 3)"),
      shape = glue("0.5 * {sq_chain} * inv_square({shape})"),
      NULL
    )
  } else if (family == "wiener") {
    if (dpar %in% c("mu", "bs", "bias")) {
      drift <- ref("mu")$value
      boundary <- ref("bs")$value
      bias <- ref("bias")$value
      start <- glue("(({boundary}) * ({bias}))")
      definitions <- c(definitions, glue(
        "        real choice_scale_fisher_s2z_wiener = ",
        "2.0 * ({drift}) * ({boundary});\n",
        "        real p_upper_fisher_s2z_wiener;\n",
        "        real dp_dscale_fisher_s2z_wiener;\n",
        "        real dp_dbias_fisher_s2z_wiener;\n",
        "        if (abs(choice_scale_fisher_s2z_wiener) < 1e-5) {{\n",
        "          p_upper_fisher_s2z_wiener = ({bias}) + 0.5 * ",
        "({bias}) * (1.0 - ({bias})) * ",
        "choice_scale_fisher_s2z_wiener;\n",
        "          dp_dscale_fisher_s2z_wiener = ",
        "0.5 * ({bias}) * (1.0 - ({bias}));\n",
        "          dp_dbias_fisher_s2z_wiener = 1.0 + 0.5 * ",
        "(1.0 - 2.0 * ({bias})) * choice_scale_fisher_s2z_wiener;\n",
        "        }} else {{\n",
        "          if (choice_scale_fisher_s2z_wiener > 0.0) {{\n",
        "            p_upper_fisher_s2z_wiener = ",
        "-expm1(-({bias}) * choice_scale_fisher_s2z_wiener) / ",
        "-expm1(-choice_scale_fisher_s2z_wiener);\n",
        "            dp_dbias_fisher_s2z_wiener = ",
        "choice_scale_fisher_s2z_wiener * exp(-({bias}) * ",
        "choice_scale_fisher_s2z_wiener) / ",
        "-expm1(-choice_scale_fisher_s2z_wiener);\n",
        "          }} else {{\n",
        "            real abs_scale_fisher_s2z_wiener = ",
        "-choice_scale_fisher_s2z_wiener;\n",
        "            p_upper_fisher_s2z_wiener = exp(-(1.0 - ({bias})) * ",
        "abs_scale_fisher_s2z_wiener) * ",
        "-expm1(-({bias}) * abs_scale_fisher_s2z_wiener) / ",
        "-expm1(-abs_scale_fisher_s2z_wiener);\n",
        "            dp_dbias_fisher_s2z_wiener = ",
        "abs_scale_fisher_s2z_wiener * exp(-(1.0 - ({bias})) * ",
        "abs_scale_fisher_s2z_wiener) / ",
        "-expm1(-abs_scale_fisher_s2z_wiener);\n",
        "          }}\n",
        "          dp_dscale_fisher_s2z_wiener = ",
        "p_upper_fisher_s2z_wiener * (",
        "({bias}) / expm1(({bias}) * choice_scale_fisher_s2z_wiener) - ",
        "inv(expm1(choice_scale_fisher_s2z_wiener)));\n",
        "        }}\n",
        "        real prob_safe_fisher_s2z_wiener = ",
        "fmin(1.0 - 1e-12, fmax(1e-12, ",
        "p_upper_fisher_s2z_wiener));\n"
      ))
      if (dpar == "mu") {
        definitions <- c(definitions, glue(
          "        real mean_time_fisher_s2z_wiener;\n",
          "        if (abs(({drift}) * ({boundary})) < 1e-5) {{\n",
          "          mean_time_fisher_s2z_wiener = ({start}) * ",
          "(({boundary}) - ({start}));\n",
          "        }} else {{\n",
          "          mean_time_fisher_s2z_wiener = (({boundary}) * ",
          "p_upper_fisher_s2z_wiener - ({start})) / ({drift});\n",
          "        }}\n"
        ))
        out <- glue(
          "{sq_chain} * fmax(0.0, mean_time_fisher_s2z_wiener)"
        )
      } else if (dpar == "bs") {
        out <- glue(
          "square(({derivative}) * 2.0 * ({drift}) * ",
          "dp_dscale_fisher_s2z_wiener) / ",
          "fmax(prob_safe_fisher_s2z_wiener * ",
          "(1.0 - prob_safe_fisher_s2z_wiener), 1e-12)"
        )
      } else {
        out <- glue(
          "square(({derivative}) * dp_dbias_fisher_s2z_wiener) / ",
          "fmax(prob_safe_fisher_s2z_wiener * ",
          "(1.0 - prob_safe_fisher_s2z_wiener), 1e-12)"
        )
      }
    }
  } else if (family == "beta") {
    mu <- ref("mu")$value
    dmu <- ref("mu")$derivative
    phi <- ref("phi")$value
    a <- glue("({mu}) * ({phi})")
    b <- glue("(1.0 - ({mu})) * ({phi})")
    mu_info <- if (identical(ref("mu")$link, "logit")) {
      glue(
        "square(1.0 - ({mu})) * {scaled_trigamma(a)} + ",
        "square({mu}) * {scaled_trigamma(b)}"
      )
    } else {
      glue(
        "square(({dmu}) / fmax({mu}, 1e-12)) * ",
        "{scaled_trigamma(a)} + square(({dmu}) / ",
        "fmax(1.0 - ({mu}), 1e-12)) * {scaled_trigamma(b)}"
      )
    }
    out <- switch(dpar,
      mu = mu_info,
      phi = glue(
        "{sq_chain} * fmax(0.0, ({scaled_trigamma(a)} + ",
        "{scaled_trigamma(b)}) * inv_square({phi}) - trigamma({phi}))"
      ),
      NULL
    )
  } else if (family == "weibull") {
    shape <- ref("shape")$value
    out <- switch(dpar,
      mu = glue("{sq_chain} * square({shape}) * inv_square({value})"),
      shape = glue(
        "{sq_chain} * (square(pi()) / 6.0 + square(",
        "digamma(1.0 + inv({shape})) - (1.0 + digamma(1.0)))) * ",
        "inv_square({shape})"
      ),
      NULL
    )
  } else if (family == "frechet") {
    nu <- ref("nu")$value
    definitions <- c(definitions, glue(
      "        real boundary_fraction_fisher_s2z_frechet = ",
      "(({nu}) - 1.0) / ({nu});\n",
      "        real scaled_shape_score_fisher_s2z_frechet = ",
      "boundary_fraction_fisher_s2z_frechet < 1e-6 ? ",
      "-1.0 - boundary_fraction_fisher_s2z_frechet + ",
      "square(pi()) / 6.0 * ",
      "square(boundary_fraction_fisher_s2z_frechet) : ",
      "boundary_fraction_fisher_s2z_frechet * (digamma(",
      "boundary_fraction_fisher_s2z_frechet) - ",
      "(1.0 + digamma(1.0)));\n",
      "        real shape_info_fisher_s2z_frechet = ",
      "square(scaled_shape_score_fisher_s2z_frechet) + ",
      "square(pi()) / 6.0 * ",
      "square(boundary_fraction_fisher_s2z_frechet);\n"
    ))
    nu_info <- if (identical(ref("nu")$link, "logm1")) {
      "shape_info_fisher_s2z_frechet"
    } else {
      glue(
        "square(({derivative}) / fmax(({nu}) - 1.0, 1e-12)) * ",
        "shape_info_fisher_s2z_frechet"
      )
    }
    out <- switch(dpar,
      mu = glue("{sq_chain} * square({nu}) * inv_square({value})"),
      nu = nu_info,
      NULL
    )
  } else if (family == "exgaussian") {
    sigma <- ref("sigma")$value
    beta <- ref("beta")$value
    variance <- glue("(square({sigma}) + square({beta}))")
    definitions <- c(definitions, glue(
      "        real dskew_sigma_fisher_s2z_exg = -6.0 * ",
      "pow({beta}, 3) * ({sigma}) * pow({variance}, -2.5);\n",
      "        real dskew_beta_fisher_s2z_exg = 6.0 * ",
      "square({beta}) * square({sigma}) * pow({variance}, -2.5);\n",
      "        real dkurt_sigma_fisher_s2z_exg = -24.0 * ",
      "pow({beta}, 4) * ({sigma}) * pow({variance}, -3.0);\n",
      "        real dkurt_beta_fisher_s2z_exg = 24.0 * ",
      "pow({beta}, 3) * square({sigma}) * pow({variance}, -3.0);\n"
    ))
    out <- switch(dpar,
      mu = glue("{sq_chain} / {variance}"),
      sigma = glue(
        "{sq_chain} * (2.0 * square({sigma}) / square({variance}) + ",
        "square(dskew_sigma_fisher_s2z_exg) / 6.0 + ",
        "square(dkurt_sigma_fisher_s2z_exg) / 24.0)"
      ),
      beta = glue(
        "{sq_chain} * (2.0 * square({beta}) / square({variance}) + ",
        "square(dskew_beta_fisher_s2z_exg) / 6.0 + ",
        "square(dkurt_beta_fisher_s2z_exg) / 24.0)"
      ),
      NULL
    )
  } else if (family == "xbeta") {
    mu <- ref("mu")$value
    phi <- ref("phi")$value
    kappa <- ref("kappa")$value
    extension <- glue("(1.0 + 2.0 * ({kappa}))")
    definitions <- c(definitions, glue(
      "        real prob_fisher_s2z_xbeta = fmin(1.0 - 1e-12, ",
      "fmax(1e-12, {mu}));\n"
    ))
    variance <- glue(
      "(square({extension}) * prob_fisher_s2z_xbeta * ",
      "(1.0 - prob_fisher_s2z_xbeta) / ",
      "(1.0 + ({phi})))"
    )
    mu_info <- if (identical(ref("mu")$link, "logit")) {
      glue(
        "(1.0 + ({phi})) * ({mu}) * (1.0 - ({mu})) + ",
        "0.5 * square(1.0 - 2.0 * ({mu}))"
      )
    } else {
      glue(
        "{sq_chain} * (square({extension}) / fmax({variance}, 1e-12) + ",
        "0.5 * square((1.0 - 2.0 * prob_fisher_s2z_xbeta) / ",
        "(prob_fisher_s2z_xbeta * ",
        "(1.0 - prob_fisher_s2z_xbeta))))"
      )
    }
    out <- switch(dpar,
      mu = mu_info,
      phi = glue(
        "0.5 * {sq_chain} * inv_square(1.0 + ({phi}))"
      ),
      kappa = glue(
        "{sq_chain} * (square(2.0 * prob_fisher_s2z_xbeta - 1.0) / ",
        "fmax({variance}, 1e-12) + ",
        "8.0 * inv_square({extension}))"
      ),
      NULL
    )
  } else if (family == "von_mises") {
    kappa <- ref("kappa")$value
    # Differentiate the rational approximation to I1(kappa) / I0(kappa)
    # directly. This preserves the exact 1/2 endpoint at zero and the leading
    # 1 / (2 kappa^2) concentration-information limit without Bessel calls.
    definitions <- c(definitions, glue(
      "        real denominator_fisher_s2z_vm = ",
      "6.0 + 4.0 * ({kappa}) + 2.0 * square({kappa});\n",
      "        real mean_cosine_over_kappa_fisher_s2z_vm = ",
      "(3.0 + 2.0 * ({kappa})) / ",
      "denominator_fisher_s2z_vm;\n",
      "        real mean_cosine_fisher_s2z_vm = ({kappa}) * ",
      "mean_cosine_over_kappa_fisher_s2z_vm;\n",
      "        real mean_cosine_derivative_fisher_s2z_vm = ",
      "(18.0 + 24.0 * ({kappa}) + 2.0 * square({kappa})) / ",
      "square(denominator_fisher_s2z_vm);\n"
    ))
    out <- switch(dpar,
      mu = glue(
        "{sq_chain} * ({kappa}) * mean_cosine_fisher_s2z_vm"
      ),
      kappa = glue(
        "{sq_chain} * mean_cosine_derivative_fisher_s2z_vm"
      ),
      NULL
    )
  } else if (family == "asym_laplace") {
    sigma <- ref("sigma")$value
    quantile <- ref("quantile")$value
    out <- switch(dpar,
      mu = glue(
        "{sq_chain} * {quantile} * (1.0 - {quantile}) * ",
        "inv_square({sigma})"
      ),
      sigma = glue("{sq_chain} * inv_square({value})"),
      quantile = quantile_information(target),
      NULL
    )
  } else if (family == "cox") {
    if (dpar == "mu") {
      # Unit information on the ideal log-rate exponential clock. brms uses a
      # bounded spline baseline, so this is a response-free working target,
      # not the exact finite-horizon Fisher information.
      out <- if (identical(target$link, "log")) {
        "1.0"
      } else {
        glue("square(({derivative}) / fmax({value}, 1e-12))")
      }
    }
  } else if (family %in% c("hurdle_gamma", "hurdle_lognormal")) {
    hu <- ref("hu")$value
    if (dpar == "hu") {
      out <- binary_information(ref("hu"))
    } else if (family == "hurdle_gamma") {
      mu <- ref("mu")$value
      shape <- ref("shape")$value
      out <- switch(dpar,
        mu = glue(
          "(1.0 - {hu}) * {sq_chain} * {shape} * inv_square({mu})"
        ),
        shape = glue(
          "(1.0 - {hu}) * {sq_chain} * fmax(0.0, ",
          "trigamma({shape}) - inv({shape}))"
        ),
        NULL
      )
    } else {
      sigma <- ref("sigma")$value
      out <- switch(dpar,
        mu = glue("(1.0 - {hu}) * {sq_chain} * inv_square({sigma})"),
        sigma = glue("(1.0 - {hu}) * 2.0 * {sq_chain} * inv_square({value})"),
        NULL
      )
    }
  } else if (family %in% c("hurdle_poisson", "hurdle_negbinomial")) {
    hu <- ref("hu")$value
    if (dpar == "hu") {
      out <- binary_information(ref("hu"))
    } else if (family == "hurdle_poisson" && dpar == "mu") {
      lambda <- ref("mu")$value
      dlambda <- ref("mu")$derivative
      poisson_info <- if (identical(ref("mu")$link, "log")) {
        lambda
      } else {
        glue("square({dlambda}) / fmax({lambda}, 1e-12)")
      }
      definitions <- c(definitions, glue(
        "        real p0_fisher_s2z_hp = exp(-({lambda}));\n",
        "        real ppos_fisher_s2z_hp = -expm1(-({lambda}));\n",
        "        real base_info_fisher_s2z_hp = ",
        "{poisson_info};\n",
        "        real zero_score_fisher_s2z_hp = -({dlambda});\n",
        "        real positive_info_fisher_s2z_hp = fmax(0.0, ",
        "base_info_fisher_s2z_hp - p0_fisher_s2z_hp * ",
        "square(zero_score_fisher_s2z_hp) / ",
        "fmax(ppos_fisher_s2z_hp, 1e-12));\n"
      ))
      out <- glue(
        "(1.0 - {hu}) * positive_info_fisher_s2z_hp / ",
        "fmax(ppos_fisher_s2z_hp, 1e-12)"
      )
    } else if (family == "hurdle_negbinomial" &&
               dpar %in% c("mu", "shape")) {
      mu <- ref("mu")$value
      shape <- ref("shape")$value
      dmu <- ref("mu")$derivative
      dshape <- ref("shape")$derivative
      definitions <- c(definitions, glue(
        "        real mean_fraction_fisher_s2z_hnb = ",
        "({mu}) / (({mu}) + ({shape}));\n",
        "        real log_p0_fisher_s2z_hnb = -({shape}) * ",
        "log1p(({mu}) / ({shape}));\n",
        "        real p0_fisher_s2z_hnb = exp(log_p0_fisher_s2z_hnb);\n",
        "        real ppos_fisher_s2z_hnb = ",
        "-expm1(log_p0_fisher_s2z_hnb);\n"
      ))
      if (dpar == "mu") {
        nb_mean_info <- if (identical(ref("mu")$link, "log")) {
          glue("({mu}) * ({shape}) / (({mu}) + ({shape}))")
        } else {
          glue(
            "square({dmu}) * ({shape}) / fmax(({mu}) * ",
            "(({mu}) + ({shape})), 1e-12)"
          )
        }
        definitions <- c(definitions, glue(
          "        real base_info_fisher_s2z_hnb = {nb_mean_info};\n",
          "        real zero_score_fisher_s2z_hnb = -({dmu}) * ",
          "({shape}) / (({mu}) + ({shape}));\n"
        ))
      } else {
        definitions <- c(definitions, glue(
          "        real dlog_p0_log_shape_fisher_s2z_hnb = ({shape}) * ",
          "(-log1p(({mu}) / ({shape})) + ",
          "mean_fraction_fisher_s2z_hnb);\n",
          "        real sparse_info_fisher_s2z_hnb = 0.5 * ({shape}) / ",
          "(({shape}) + 1.0) * square(mean_fraction_fisher_s2z_hnb);\n",
          "        real dense_info_fisher_s2z_hnb = 0.5 * ",
          "square(mean_fraction_fisher_s2z_hnb);\n",
          "        real zero_info_fisher_s2z_hnb = p0_fisher_s2z_hnb * ",
          "square(dlog_p0_log_shape_fisher_s2z_hnb) / ",
          "fmax(ppos_fisher_s2z_hnb, 1e-12);\n",
          "        real log_shape_info_fisher_s2z_hnb = ",
          "sparse_info_fisher_s2z_hnb + zero_info_fisher_s2z_hnb + ",
          "ppos_fisher_s2z_hnb * (dense_info_fisher_s2z_hnb - ",
          "sparse_info_fisher_s2z_hnb);\n",
          "        real log_shape_chain_fisher_s2z_hnb = ",
          "({dshape}) / ({shape});\n",
          "        real base_info_fisher_s2z_hnb = ",
          "square(log_shape_chain_fisher_s2z_hnb) * ",
          "log_shape_info_fisher_s2z_hnb;\n",
          "        real zero_score_fisher_s2z_hnb = ",
          "log_shape_chain_fisher_s2z_hnb * ",
          "dlog_p0_log_shape_fisher_s2z_hnb;\n"
        ))
      }
      definitions <- c(definitions, glue(
        "        real positive_info_fisher_s2z_hnb = fmax(0.0, ",
        "base_info_fisher_s2z_hnb - p0_fisher_s2z_hnb * ",
        "square(zero_score_fisher_s2z_hnb) / ",
        "fmax(ppos_fisher_s2z_hnb, 1e-12));\n"
      ))
      out <- glue(
        "(1.0 - {hu}) * positive_info_fisher_s2z_hnb / ",
        "fmax(ppos_fisher_s2z_hnb, 1e-12)"
      )
    }
  } else if (family %in% c(
    "zero_inflated_poisson", "zero_inflated_negbinomial",
    "zero_inflated_binomial", "zero_inflated_beta_binomial"
  )) {
    zi <- ref("zi")$value
    dzi <- ref("zi")$derivative
    if (family == "zero_inflated_poisson") {
      base_value <- ref("mu")$value
      base_derivative <- ref("mu")$derivative
      definitions <- c(definitions, glue(
        "        real p0_fisher_s2z_zi = exp(-({base_value}));\n",
        "        real ppos_fisher_s2z_zi = -expm1(-({base_value}));\n",
        "        real q0_fisher_s2z_zi = ({zi}) + ",
        "(1.0 - ({zi})) * p0_fisher_s2z_zi;\n",
        "        real qpos_fisher_s2z_zi = ",
        "(1.0 - ({zi})) * ppos_fisher_s2z_zi;\n"
      ))
      if (dpar == "mu") {
        poisson_info <- if (identical(ref("mu")$link, "log")) {
          base_value
        } else {
          glue(
            "square({base_derivative}) / fmax({base_value}, 1e-12)"
          )
        }
        definitions <- c(definitions, glue(
          "        real base_info_fisher_s2z_zi = ",
          "{poisson_info};\n",
          "        real zero_score_fisher_s2z_zi = -({base_derivative});\n"
        ))
      }
    } else if (family == "zero_inflated_binomial") {
      prob <- ref("mu")$value
      dprob <- ref("mu")$derivative
      trials <- glue("trials{resp}[{n}]")
      definitions <- c(definitions, glue(
        "        real log_p0_fisher_s2z_zi = ({trials}) == 0 ? ",
        "0.0 : ({trials}) * log1m({prob});\n",
        "        real p0_fisher_s2z_zi = exp(log_p0_fisher_s2z_zi);\n",
        "        real ppos_fisher_s2z_zi = ",
        "-expm1(log_p0_fisher_s2z_zi);\n",
        "        real q0_fisher_s2z_zi = ({zi}) + ",
        "(1.0 - ({zi})) * p0_fisher_s2z_zi;\n",
        "        real qpos_fisher_s2z_zi = ",
        "(1.0 - ({zi})) * ppos_fisher_s2z_zi;\n"
      ))
      if (dpar == "mu") {
        binomial_info <- binary_information(ref("mu"))
        zero_score <- if (identical(ref("mu")$link, "logit")) {
          glue("-({trials}) * ({prob})")
        } else {
          glue(
            "-({trials}) * ({dprob}) / fmax(1.0 - ({prob}), 1e-12)"
          )
        }
        definitions <- c(definitions, glue(
          "        real base_info_fisher_s2z_zi = ({trials}) * ",
          "({binomial_info});\n",
          "        real zero_score_fisher_s2z_zi = {zero_score};\n"
        ))
      }
    } else if (family == "zero_inflated_beta_binomial") {
      prob <- ref("mu")$value
      dprob <- ref("mu")$derivative
      phi <- ref("phi")$value
      dphi <- ref("phi")$derivative
      trials <- glue("trials{resp}[{n}]")
      definitions <- c(definitions, glue(
        "        real prob_safe_fisher_s2z_zibb = fmin(1.0 - 1e-12, ",
        "fmax(1e-12, {prob}));\n",
        "        real alpha_fisher_s2z_zibb = ",
        "prob_safe_fisher_s2z_zibb * ({phi});\n",
        "        real beta_fisher_s2z_zibb = ",
        "(1.0 - prob_safe_fisher_s2z_zibb) * ({phi});\n",
        "        real log_p0_fisher_s2z_zi = ",
        "lbeta(alpha_fisher_s2z_zibb, beta_fisher_s2z_zibb + ",
        "({trials})) - lbeta(alpha_fisher_s2z_zibb, ",
        "beta_fisher_s2z_zibb);\n",
        "        real log_pn_fisher_s2z_zibb = ",
        "lbeta(alpha_fisher_s2z_zibb + ({trials}), ",
        "beta_fisher_s2z_zibb) - lbeta(alpha_fisher_s2z_zibb, ",
        "beta_fisher_s2z_zibb);\n",
        "        real p0_fisher_s2z_zi = exp(log_p0_fisher_s2z_zi);\n",
        "        real pn_fisher_s2z_zibb = exp(log_pn_fisher_s2z_zibb);\n",
        "        real pmid_fisher_s2z_zibb = fmax(0.0, ",
        "1.0 - p0_fisher_s2z_zi - pn_fisher_s2z_zibb);\n",
        "        real ppos_fisher_s2z_zi = ",
        "-expm1(log_p0_fisher_s2z_zi);\n",
        "        real q0_fisher_s2z_zi = ({zi}) + ",
        "(1.0 - ({zi})) * p0_fisher_s2z_zi;\n",
        "        real qpos_fisher_s2z_zi = ",
        "(1.0 - ({zi})) * ppos_fisher_s2z_zi;\n"
      ))
      if (dpar == "mu") {
        definitions <- c(definitions, glue(
          "        real zero_score_fisher_s2z_zi = ({dprob}) * ",
          "({phi}) * (digamma(beta_fisher_s2z_zibb) - ",
          "digamma(beta_fisher_s2z_zibb + ({trials})));\n",
          "        real end_score_fisher_s2z_zibb = ({dprob}) * ",
          "({phi}) * (digamma(alpha_fisher_s2z_zibb + ({trials})) - ",
          "digamma(alpha_fisher_s2z_zibb));\n"
        ))
      } else if (dpar == "phi") {
        definitions <- c(definitions, glue(
          "        real zero_score_fisher_s2z_zi = ({dphi}) * (",
          "digamma({phi}) - digamma(({phi}) + ({trials})) + ",
          "(1.0 - prob_safe_fisher_s2z_zibb) * (",
          "digamma(beta_fisher_s2z_zibb + ",
          "({trials})) - digamma(beta_fisher_s2z_zibb)));\n",
          "        real end_score_fisher_s2z_zibb = ({dphi}) * (",
          "digamma({phi}) - digamma(({phi}) + ({trials})) + ",
          "prob_safe_fisher_s2z_zibb * (",
          "digamma(alpha_fisher_s2z_zibb + ({trials})) - ",
          "digamma(alpha_fisher_s2z_zibb)));\n"
        ))
      }
      if (dpar %in% c("mu", "phi")) {
        definitions <- c(definitions, glue(
          "        real dpmid_fisher_s2z_zibb = -p0_fisher_s2z_zi * ",
          "zero_score_fisher_s2z_zi - pn_fisher_s2z_zibb * ",
          "end_score_fisher_s2z_zibb;\n",
          "        real coarse_info_fisher_s2z_zibb = ",
          "p0_fisher_s2z_zi * square(zero_score_fisher_s2z_zi) + ",
          "pn_fisher_s2z_zibb * square(end_score_fisher_s2z_zibb) + ",
          "square(dpmid_fisher_s2z_zibb) / ",
          "fmax(pmid_fisher_s2z_zibb, 1e-12);\n",
          "        real base_info_fisher_s2z_zi = ",
          str_if(
            dpar == "mu",
            glue(
              "({trials} == 0 ? 0.0 : ({trials} == 1 ? ",
              "({binary_information(ref('mu'))}) : ",
              "coarse_info_fisher_s2z_zibb));\n"
            ),
            glue(
              "({trials} <= 1 ? 0.0 : ",
              "coarse_info_fisher_s2z_zibb);\n"
            )
          )
        ))
      }
    } else {
      mu <- ref("mu")$value
      dmu <- ref("mu")$derivative
      shape <- ref("shape")$value
      dshape <- ref("shape")$derivative
      definitions <- c(definitions, glue(
        "        real mean_fraction_fisher_s2z_zinb = ",
        "({mu}) / (({mu}) + ({shape}));\n",
        "        real log_p0_fisher_s2z_zi = -({shape}) * ",
        "log1p(({mu}) / ({shape}));\n",
        "        real p0_fisher_s2z_zi = exp(log_p0_fisher_s2z_zi);\n",
        "        real ppos_fisher_s2z_zi = ",
        "-expm1(log_p0_fisher_s2z_zi);\n",
        "        real q0_fisher_s2z_zi = ({zi}) + ",
        "(1.0 - ({zi})) * p0_fisher_s2z_zi;\n",
        "        real qpos_fisher_s2z_zi = ",
        "(1.0 - ({zi})) * ppos_fisher_s2z_zi;\n"
      ))
      if (dpar == "mu") {
        nb_mean_info <- if (identical(ref("mu")$link, "log")) {
          glue("({mu}) * ({shape}) / (({mu}) + ({shape}))")
        } else {
          glue(
            "square({dmu}) * ({shape}) / fmax(({mu}) * ",
            "(({mu}) + ({shape})), 1e-12)"
          )
        }
        definitions <- c(definitions, glue(
          "        real base_info_fisher_s2z_zi = {nb_mean_info};\n",
          "        real zero_score_fisher_s2z_zi = -({dmu}) * ",
          "({shape}) / (({mu}) + ({shape}));\n"
        ))
      } else if (dpar == "shape") {
        definitions <- c(definitions, glue(
          "        real dlog_p0_log_shape_fisher_s2z_zinb = ({shape}) * ",
          "(-log1p(({mu}) / ({shape})) + ",
          "mean_fraction_fisher_s2z_zinb);\n",
          "        real sparse_info_fisher_s2z_zinb = 0.5 * ({shape}) / ",
          "(({shape}) + 1.0) * square(mean_fraction_fisher_s2z_zinb);\n",
          "        real dense_info_fisher_s2z_zinb = 0.5 * ",
          "square(mean_fraction_fisher_s2z_zinb);\n",
          "        real zero_info_fisher_s2z_zinb = p0_fisher_s2z_zi * ",
          "square(dlog_p0_log_shape_fisher_s2z_zinb) / ",
          "fmax(ppos_fisher_s2z_zi, 1e-12);\n",
          "        real log_shape_info_fisher_s2z_zinb = ",
          "sparse_info_fisher_s2z_zinb + zero_info_fisher_s2z_zinb + ",
          "ppos_fisher_s2z_zi * (dense_info_fisher_s2z_zinb - ",
          "sparse_info_fisher_s2z_zinb);\n",
          "        real log_shape_chain_fisher_s2z_zinb = ",
          "({dshape}) / ({shape});\n",
          "        real base_info_fisher_s2z_zi = ",
          "square(log_shape_chain_fisher_s2z_zinb) * ",
          "log_shape_info_fisher_s2z_zinb;\n",
          "        real zero_score_fisher_s2z_zi = ",
          "log_shape_chain_fisher_s2z_zinb * ",
          "dlog_p0_log_shape_fisher_s2z_zinb;\n"
        ))
      }
    }
    if (dpar == "zi") {
      out <- glue(
        "square(({dzi}) * ppos_fisher_s2z_zi) / ",
        "fmax(q0_fisher_s2z_zi * qpos_fisher_s2z_zi, 1e-12)"
      )
    } else if (dpar %in% setdiff(valid_dpars(bframe), "zi")) {
      definitions <- c(definitions, glue(
        "        real positive_info_fisher_s2z_zi = fmax(0.0, ",
        "base_info_fisher_s2z_zi - p0_fisher_s2z_zi * ",
        "square(zero_score_fisher_s2z_zi) / ",
        "fmax(ppos_fisher_s2z_zi, 1e-12));\n",
        "        real atom_derivative_fisher_s2z_zi = ",
        "(1.0 - ({zi})) * p0_fisher_s2z_zi * ",
        "zero_score_fisher_s2z_zi;\n"
      ))
      out <- glue(
        "(1.0 - ({zi})) * positive_info_fisher_s2z_zi + ",
        "square(atom_derivative_fisher_s2z_zi) / ",
        "fmax(q0_fisher_s2z_zi * qpos_fisher_s2z_zi, 1e-12)"
      )
    }
  } else if (family %in% c(
    "zero_inflated_beta", "zero_one_inflated_beta",
    "zero_inflated_asym_laplace"
  )) {
    inflation <- if (family == "zero_one_inflated_beta") "zoi" else "zi"
    atom <- ref(inflation)$value
    if (dpar == inflation) {
      out <- binary_information(ref(inflation))
    } else if (family == "zero_one_inflated_beta" && dpar == "coi") {
      coi <- ref("coi")$value
      out <- glue("{atom} * ({binary_information(ref('coi'))})")
    } else if (family %in% c("zero_inflated_beta", "zero_one_inflated_beta") &&
               dpar %in% c("mu", "phi")) {
      mu <- ref("mu")$value
      dmu <- ref("mu")$derivative
      phi <- ref("phi")$value
      a <- glue("({mu}) * ({phi})")
      b <- glue("(1.0 - ({mu})) * ({phi})")
      base <- if (dpar == "mu") {
        if (identical(ref("mu")$link, "logit")) {
          glue(
            "square(1.0 - ({mu})) * {scaled_trigamma(a)} + ",
            "square({mu}) * {scaled_trigamma(b)}"
          )
        } else {
          glue(
            "square(({dmu}) / fmax({mu}, 1e-12)) * ",
            "{scaled_trigamma(a)} + square(({dmu}) / ",
            "fmax(1.0 - ({mu}), 1e-12)) * {scaled_trigamma(b)}"
          )
        }
      } else {
        glue(
          "{sq_chain} * fmax(0.0, ({scaled_trigamma(a)} + ",
          "{scaled_trigamma(b)}) * inv_square({phi}) - trigamma({phi}))"
        )
      }
      out <- if (dpar == "mu") {
        glue("(1.0 - {atom}) * ({base})")
      } else {
        glue("(1.0 - {atom}) * ({base})")
      }
    } else if (family == "zero_inflated_asym_laplace" &&
               dpar %in% c("mu", "sigma", "quantile")) {
      sigma <- ref("sigma")$value
      quantile <- ref("quantile")$value
      base <- if (dpar == "mu") {
        glue("{quantile} * (1.0 - {quantile}) * inv_square({sigma})")
      } else if (dpar == "sigma") {
        glue("inv_square({sigma})")
      } else {
        NULL
      }
      out <- if (dpar == "quantile") {
        glue("(1.0 - {atom}) * ({quantile_information(target)})")
      } else {
        glue("(1.0 - {atom}) * {sq_chain} * ({base})")
      }
    }
  }
  if (is.null(out)) {
    return(NULL)
  }
  nlist(
    obs_prec_at_n = out,
    definitions_at_n = paste0(definitions, collapse = ""),
    target, family, dpar
  )
}

# Validate the expected-information proposal and return the observation-level
# quantities needed to construct it in precursor generated quantities. A
# transformed-data Gram fast path remains for scalar Gaussian/Student location
# information. Unsupported likelihood coordinates fail during code generation.
stan_re_s2z_fisher_info <- function(id, r, bframe, threads) {
  stopifnot(is.reframe(r), has_rows(r))
  effect_label <- if (all(r$s2z)) {
    "S2Z group-level effects"
  } else {
    "ordinary group-level effects"
  }
  unsupported <- function(detail) {
    stop2("Automatic centering for ", effect_label, " ", detail, ".")
  }
  if (any(nzchar(r$nlpar))) {
    unsupported(paste0(
      "does not yet support nonlinear predictors; the derivative of the ",
      "response mean with respect to the latent score must be represented ",
      "explicitly"
    ))
  }
  if (any(nzchar(r$gtype)) || any(nzchar(r$type)) ||
      any(nzchar(r$by)) || isTRUE(nzchar(r$gcall[[1]]$pw))) {
    unsupported(paste0(
      "currently requires an ordinary gr() block without by, pw, ",
      "multi-membership, or special-effect structure"
    ))
  }
  model_frame <- bframe
  if (is.mvbrmsframe(bframe)) {
    if (isTRUE(bframe$rescor)) {
      unsupported(paste0(
        "for predictor-local IDs in multivariate models currently requires ",
        "set_rescor(FALSE)"
      ))
    }
    responses <- unique(r$resp)
    if (length(responses) != 1L || !responses %in% names(bframe$terms)) {
      unsupported("requires one predictor-local response block")
    }
    model_frame <- bframe$terms[[responses]]
  }
  if (!is.brmsframe(model_frame)) {
    unsupported("could not identify one response model")
  }
  bfl <- re_s2z_bframel(model_frame, r)
  response_family <- model_frame$family$family
  dpar <- bfl$dpar
  if (has_rows(bfl$frame$ac)) {
    unsupported("does not yet support residual autocorrelation structures")
  }
  used_adterms <- names(Filter(is.formula, bfl$adforms))
  allowed_adterms <- c("trials", "rate", "subset", "index")
  if (identical(response_family, "cox")) {
    allowed_adterms <- c(allowed_adterms, "bhaz")
  }
  if (identical(response_family, "wiener")) {
    allowed_adterms <- c(allowed_adterms, "dec")
  }
  bad_adterms <- setdiff(used_adterms, allowed_adterms)
  if (length(bad_adterms)) {
    unsupported(paste0(
      "does not yet support response addition term",
      ifelse(length(bad_adterms) == 1L, " '", "s '"),
      paste0(bad_adterms, collapse = "', '"), "'"
    ))
  }
  resp <- usc(unique(r$resp))
  idp <- paste0(r$id, usc(combine_prefix(check_prefix(r))))
  design <- glue("Z_{idp}_{r$cn}")
  M <- nrow(r)
  N <- glue("N{resp}")
  group <- glue("J_{id}{resp}")
  design_at_n <- paste0(design, "[n]")

  fast_location_candidate <- response_family %in% c("gaussian", "student") &&
    identical(bfl$family$link, "identity") && identical(dpar, "mu")
  if (fast_location_candidate) {
    sigma <- stan_sigma_transform(model_frame, threads = threads)
    expected_sigma <- glue("sigma{resp}")
    fast_location <-
      identical(as.character(sigma), as.character(expected_sigma)) &&
      !(identical(response_family, "student") &&
        "nu" %in% names(bframe$dpars))
  } else {
    fast_location <- FALSE
  }
  if (fast_location) {
    nu <- glue("nu{resp}")
    obs_prec <- if (identical(response_family, "gaussian")) {
      glue("inv_square({sigma})")
    } else {
      # Expected location information for Student-t(nu, mu, sigma).
      glue("({nu} + 1.0) / ({nu} + 3.0) * inv_square({sigma})")
    }
    fixed_design <- TRUE
    return(nlist(
      M, fixed_design, N, group, design, design_at_n,
      sigma, obs_prec, response_family, dpar
    ))
  }

  closed <- stan_re_s2z_fisher_closed_form(
    model_frame, dpar = dpar,
    reference_eta = stan_re_s2z_fisher_fitted_references(model_frame)
  )
  if (is.null(closed)) {
    unsupported(paste0(
      "has no response-free expected-information rule for family '",
      response_family, "' and distributional parameter '", dpar,
      "'"
    ))
  }
  source <- list(
    columns = seq_len(M), N = N, group = group,
    design_at_n = design_at_n,
    obs_prec_at_n = closed$obs_prec_at_n,
    definitions_at_n = closed$definitions_at_n
  )
  fixed_design <- FALSE
  nlist(
    M, fixed_design, sources = list(source), response_family, dpar, fun = ""
  )
}

# Describe nonlinear expected information for a strict latent-score block. The
# effective design is Z_nk * d zeta_n / d xi_k, where zeta is the response
# link-scale predictor and xi is the strict score predictor. Population-only
# dependencies of that derivative (sampled loadings in particular) are
# evaluated with each precursor draw in generated quantities. Dependence on
# the latent scores themselves is rejected by the symbolic analyzer.
stan_re_s2z_latent_fisher_info <- function(id, r, bframe, threads) {
  stopifnot(is.reframe(r), has_rows(r), all(r$s2z), all(re_s2z_latent(r)))
  unsupported <- function(detail) {
    stop2("Automatic centering for strict latent-score S2Z blocks ",
          detail, ".")
  }
  if (is.mvbrmsframe(bframe) && isTRUE(bframe$rescor)) {
    unsupported(paste0(
      "currently requires set_rescor(FALSE); residual cross-response ",
      "precision is not yet included in the Fisher target"
    ))
  }
  responses <- unique(r$resp)
  dimension <- re_s2z_latent_dimension(r)
  M <- max(dimension)
  term_for_response <- function(resp) {
    if (is.brmsframe(bframe)) {
      if (length(responses) != 1L) {
        unsupported("could not match response-local nonlinear predictors")
      }
      return(bframe)
    }
    if (!is.mvbrmsframe(bframe) || !resp %in% names(bframe$terms)) {
      unsupported("could not match response-local nonlinear predictors")
    }
    bframe$terms[[resp]]
  }

  sources <- vector("list", length(responses))
  dependency_records <- list()
  for (s in seq_along(responses)) {
    response <- responses[s]
    rows <- which(r$resp == response)
    r_source <- r[rows, , drop = FALSE]
    term <- term_for_response(response)
    target_nlpars <- unique(r_source$nlpar)
    if (!length(target_nlpars) || any(!nzchar(target_nlpars)) ||
        any(!target_nlpars %in% names(term$nlpars))) {
      unsupported("must belong to explicit nonlinear score predictors")
    }
    if (has_rows(term$frame$ac)) {
      unsupported("does not yet support residual autocorrelation structures")
    }
    analyses <- lapply(target_nlpars, function(nlpar) {
      re_s2z_fisher_nl_info(
        term, term$nlpars[[nlpar]], id = id, strict = TRUE
      )
    })
    names(analyses) <- target_nlpars
    used_adterms <- unique(unlist(lapply(
      analyses, `[[`, "response_addition_terms"
    ), use.names = FALSE))
    response_family <- unique(vapply(
      analyses, function(x) x$response_family %||% "gaussian", character(1)
    ))
    if (length(response_family) != 1L) {
      unsupported("could not identify one response family")
    }
    allowed_adterms <- c("trials", "rate", "subset", "index")
    if (identical(response_family, "cox")) {
      allowed_adterms <- c(allowed_adterms, "bhaz")
    }
    if (identical(response_family, "wiener")) {
      allowed_adterms <- c(allowed_adterms, "dec")
    }
    bad_adterms <- setdiff(used_adterms, allowed_adterms)
    if (length(bad_adterms)) {
      unsupported(paste0(
        "does not yet support response addition term",
        ifelse(length(bad_adterms) == 1L, " '", "s '"),
        paste0(bad_adterms, collapse = "', '"), "'"
      ))
    }
    outer_references <- unique(vapply(
      analyses, `[[`, character(1), "outer_reference_stan"
    ))
    if (length(outer_references) != 1L) {
      unsupported(paste0(
        "could not construct one common score-free nonlinear response ",
        "reference"
      ))
    }
    closed <- stan_re_s2z_fisher_closed_form(
      term, dpar = "mu",
      reference_eta = list(mu = outer_references[[1L]])
    )
    if (is.null(closed)) {
      unsupported(paste0(
        "has no response-free expected-information rule for family '",
        response_family, "' and nonlinear location parameter 'mu'"
      ))
    }
    resp <- usc(response)
    derivatives <- vapply(r_source$nlpar, function(nlpar) {
      analyses[[nlpar]]$obs_derivative_stan
    }, character(1))
    idp <- paste0(
      r_source$id, usc(combine_prefix(check_prefix(r_source)))
    )
    group_design <- glue("Z_{idp}_{r_source$cn}[n]")
    columns <- dimension[rows]
    if (anyDuplicated(columns)) {
      stop2("A strict latent-score dimension may occur only once in each ",
            "response predictor.")
    }
    sources[[s]] <- nlist(
      response, N = glue("N{resp}"), group = glue("J_{id}{resp}"),
      columns,
      design_at_n = paste0("(", derivatives, ") * ", group_design),
      obs_prec_at_n = closed$obs_prec_at_n,
      definitions_at_n = closed$definitions_at_n,
      response_family
    )
    for (analysis in analyses) {
      for (dependency in unname(analysis$dependency_info)) {
        dependency_records[[length(dependency_records) + 1L]] <- list(
          dependency = dependency, N = glue("N{resp}")
        )
      }
    }
  }

  dependency_info <- lapply(dependency_records, `[[`, "dependency")
  dependency_tpar_def <- dependency_tpar_comp <- ""
  if (length(dependency_info)) {
    dependency_names <- vapply(
      dependency_info, function(x) x$vector_name, character(1)
    )
    keep <- !duplicated(dependency_names)
    duplicates <- unique(dependency_names[duplicated(dependency_names)])
    for (name in duplicates) {
      expressions <- unique(vapply(
        dependency_info[dependency_names == name],
        function(x) x$vector_expression, character(1)
      ))
      if (length(expressions) != 1L) {
        stop2("Internal mismatch while hoisting nonlinear Fisher dependency '",
              name, "'.")
      }
    }
    dependency_records <- dependency_records[keep]
    dependency_info <- dependency_info[keep]
    dependency_tpar_def <- paste0(vapply(
      dependency_records,
      function(x) glue(
        "  vector[{x$N}] {x$dependency$vector_name};\n"
      ),
      character(1)
    ), collapse = "")
    dependency_tpar_comp <- paste0(vapply(
      dependency_info,
      function(x) glue(
        "  {x$vector_name} = {x$vector_expression};\n"
      ),
      character(1)
    ), collapse = "")
  }

  fixed_design <- FALSE
  nlist(
    M, fixed_design, sources, dependency_info,
    dependency_tpar_def, dependency_tpar_comp
  )
}

# Stan declarations for precursor-only marginal Fisher reliabilities.  The
# model itself always consumes rho_s2z_ID as fixed data.  Candidate values are
# evaluated in generated quantities only, so neither Pathfinder's optimization
# nor the final HMC gradients differentiate through the reliability heuristic.
stan_re_s2z_fisher_def <- function(out, id) {
  str_add(out$data) <- glue(
    "  matrix<lower=0,upper=1>[N_{id}, M_{id}] rho_s2z_{id};",
    "  // fixed precursor/final centering fractions\n",
    "  int<lower=0,upper=1> compute_rho_center_candidate_{id};",
    "  // evaluate the precursor proposal in generated quantities?\n"
  )
  out <- stan_re_s2z_partial_tdata(out, id)
  str_add(out$gen_def) <- glue(
    "  matrix<lower=0,upper=1>[N_{id}, M_{id}] ",
    "rho_center_candidate_{id};\n",
    "  vector<lower=0,upper=1>[M_{id}] ",
    "mean_rho_center_candidate_{id};\n"
  )
  out
}

# Evaluate a candidate map only for the precursor run.  In the final fit the
# generated quantity is a cheap copy of the frozen input, which keeps one Stan
# executable valid for both stages without repeating Fisher work on saved HMC
# draws.
stan_re_s2z_fisher_gq_comp <- function(id, r, fisher_info, L = NULL,
                                       scale = NULL, row_var = NULL,
                                       diag_scale = NULL,
                                       precompute = "") {
  rho <- glue("rho_center_candidate_{id}")
  mean_rho <- glue("mean_rho_center_candidate_{id}")
  proposal <- stan_re_s2z_fisher_comp(
    id, r = r, fisher_info = fisher_info, L = L, scale = scale,
    row_var = row_var, diag_scale = diag_scale,
    rho = rho, mean_rho = mean_rho
  )
  glue(
    "  if (compute_rho_center_candidate_{id}) {{\n",
    "{precompute}",
    "{proposal}",
    "  }} else {{\n",
    "    {rho} = rho_s2z_{id};\n",
    "    {mean_rho} = mean_rho_s2z_{id};\n",
    "  }}\n"
  )
}

# Shape of the group contribution to the omitted-mean precision in reference-
# whitened coefficient coordinates. Shared scales are isotropic even with a
# known grouping covariance or Student mixing; independent varying scales are
# diagonal; only the remaining varying-scale cases require a dense matrix.
stan_re_s2z_group_precision_kind <- function(varying, is_cor) {
  if (!varying) {
    return("scalar")
  }
  if (!is_cor) {
    return("diagonal")
  }
  "dense"
}

stan_re_s2z_group_precision_def <- function(id, kind) {
  switch(
    kind,
    scalar = glue(
      "  real<lower=0> group_info_s2z_{id};",
      "  // isotropic omitted-mean precision\n"
    ),
    diagonal = glue(
      "  vector<lower=0>[M_{id}] group_info_s2z_{id};",
      "  // diagonal omitted-mean precision\n"
    ),
    dense = glue(
      "  matrix[M_{id}, M_{id}] P_group_s2z_{id};\n"
    ),
    stop2("Invalid S2Z group precision kind: ", kind)
  )
}

stan_re_s2z_add_group_precision <- function(base, id, kind) {
  if (kind %in% c("scalar", "diagonal")) {
    glue("add_diag({base}, group_info_s2z_{id})")
  } else {
    glue("{base} + P_group_s2z_{id}")
  }
}

stan_re_s2z_add_group_precision_block <- function(set_id, id, kind, start,
                                                   end) {
  index <- if (start == 1L) "k" else glue("{start - 1L} + k")
  if (kind == "scalar") {
    return(glue(
      "    for (k in 1:M_{id}) {{\n",
      "      P_s2z_{set_id}[{index}, {index}] += ",
      "group_info_s2z_{id};\n",
      "    }}\n"
    ))
  }
  if (kind == "diagonal") {
    return(glue(
      "    for (k in 1:M_{id}) {{\n",
      "      P_s2z_{set_id}[{index}, {index}] += ",
      "group_info_s2z_{id}[k];\n",
      "    }}\n"
    ))
  }
  glue(
    "    P_s2z_{set_id}[{start}:{end}, {start}:{end}] += ",
    "P_group_s2z_{id};\n"
  )
}


# Ordinary group-level designs and indices are fixed data. Precompute their
# unweighted within-level Gram matrices once so each generated-quantity
# proposal only combines the expected observation precision and covariance
# parameters. Strict latent-score designs may contain sampled loadings and
# deliberately bypass this path.
stan_re_s2z_fisher_tdata <- function(out, id, r, fisher_info) {
  stopifnot(
    is.reframe(r), has_rows(r), isTRUE(fisher_info$fixed_design),
    identical(fisher_info$M, nrow(r))
  )
  M <- fisher_info$M
  design_at_n <- fisher_info$design_at_n
  if (is.null(design_at_n)) {
    design_at_n <- paste0(fisher_info$design, "[n]")
  }
  stopifnot(length(design_at_n) == M)
  if (M == 1L) {
    str_add(out$tdata_def) <- glue(
      "  vector<lower=0>[N_{id}] exposure_fisher_s2z_{id};\n"
    )
    str_add(out$tdata_comp) <- glue(
      "  exposure_fisher_s2z_{id} = zeros_vector(N_{id});\n",
      "  for (n in 1:{fisher_info$N}) {{\n",
      "    exposure_fisher_s2z_{id}[{fisher_info$group}[n]] += ",
      "square({design_at_n[1]});\n",
      "  }}\n"
    )
    return(out)
  }

  design_assign <- cglue(
    "    design_fisher_s2z[{seq_len(M)}] = {design_at_n};\n"
  )
  str_add(out$tdata_def) <- glue(
    "  array[N_{id}] matrix[M_{id}, M_{id}] gram_fisher_s2z_{id};\n"
  )
  str_add(out$tdata_comp) <- glue(
    "  for (j in 1:N_{id}) {{\n",
    "    gram_fisher_s2z_{id}[j] = rep_matrix(0.0, M_{id}, M_{id});\n",
    "  }}\n",
    "  for (n in 1:{fisher_info$N}) {{\n",
    "    vector[M_{id}] design_fisher_s2z;\n",
    "{design_assign}",
    "    gram_fisher_s2z_{id}[{fisher_info$group}[n]] += ",
    "design_fisher_s2z * design_fisher_s2z';\n",
    "  }}\n"
  )
  out
}

# Compute local Gaussian variance contractions without an eigendecomposition
# or explicit inverse. For J_j = sum_n w_F z_n z_n' and Sigma = L L', the
# posterior covariance is V_j = L (I + L' J_j L)^-1 L'. For iid Gaussian S2Z
# blocks, condition these covariances on sum_j u_j = 0 before taking marginal
# contractions; other supported structures retain their local contraction.
stan_re_s2z_fisher_comp <- function(id, r, fisher_info, L = NULL,
                                    scale = NULL, row_var = NULL,
                                    diag_scale = NULL,
                                    rho = glue("rho_s2z_{id}"),
                                    mean_rho = glue("mean_rho_s2z_{id}")) {
  fixed_design <- isTRUE(fisher_info$fixed_design)
  sources <- fisher_info$sources
  M <- fisher_info$M
  stopifnot(length(M) == 1L, M >= 1L)
  row_var_j <- row_var %||% "1.0"
  constrained_s2z <- all(r$s2z) && all(r$dist == "gaussian") &&
    all(r$scale == "shared") && is.null(row_var)

  # For one coefficient, posterior variance contraction is the scalar
  # reliability I * tau^2 / (1 + I * tau^2). Ordinary designs use the
  # transformed-data exposure; sampled-loading latent designs retain only
  # their necessarily draw-specific scalar information accumulation here.
  if (M == 1L) {
    stopifnot(length(scale) == 1L)
    if (fixed_design) {
      info_def <- glue(
        "    real obs_prec_fisher_s2z = {fisher_info$obs_prec};\n"
      )
      accumulation <- ""
      information_j <- glue(
        "obs_prec_fisher_s2z * exposure_fisher_s2z_{id}[j]"
      )
    } else {
      stopifnot(length(sources) >= 1L)
      info_def <- glue(
        "    vector[N_{id}] info_fisher_s2z = zeros_vector(N_{id});\n"
      )
      accumulation <- paste0(vapply(sources, function(source) {
        stopifnot(
          identical(source$columns, 1L),
          length(source$design_at_n) == 1L
        )
        dynamic_prec <- !is.null(source$obs_prec_at_n)
        obs_prec <- source$obs_prec_at_n %||% source$obs_prec
        definitions <- source$definitions_at_n %||% ""
        glue(
          "    {{\n",
          str_if(
            !dynamic_prec,
            glue("      real obs_prec_fisher_s2z = {obs_prec};\n")
          ),
          "      for (n in 1:{source$N}) {{\n",
          "{definitions}",
          str_if(
            dynamic_prec,
            glue("        real obs_prec_fisher_s2z = {obs_prec};\n")
          ),
          "        info_fisher_s2z[{source$group}[n]] += ",
          "obs_prec_fisher_s2z * square({source$design_at_n});\n",
          "      }}\n",
          "    }}\n"
        )
      }, character(1)), collapse = "")
      information_j <- "info_fisher_s2z[j]"
    }
    if (constrained_s2z) {
      return(glue(
        "  {{\n",
        "{info_def}",
        "    vector[N_{id}] relative_post_var_fisher_s2z;\n",
        "    real restricted_prior_fraction_fisher_s2z = ",
        "1.0 - inv(N_{id});\n",
        "    real sum_relative_post_var_fisher_s2z = 0.0;\n",
        "{accumulation}",
        "    for (j in 1:N_{id}) {{\n",
        "      real scaled_info_fisher_s2z = square({scale}) * ",
        "{information_j};\n",
        "      relative_post_var_fisher_s2z[j] = ",
        "inv(1.0 + scaled_info_fisher_s2z);\n",
        "      sum_relative_post_var_fisher_s2z += ",
        "relative_post_var_fisher_s2z[j];\n",
        "    }}\n",
        "    for (j in 1:N_{id}) {{\n",
        "      real restricted_relative_post_var_fisher_s2z = ",
        "relative_post_var_fisher_s2z[j] - ",
        "square(relative_post_var_fisher_s2z[j]) / ",
        "sum_relative_post_var_fisher_s2z;\n",
        "      // Compare posterior and prior variances on the S2Z subspace.\n",
        "      {rho}[j, 1] = fmin(1.0, fmax(0.0, 1.0 - ",
        "restricted_relative_post_var_fisher_s2z / ",
        "restricted_prior_fraction_fisher_s2z));\n",
        "    }}\n",
        "    {mean_rho}[1] = mean({rho}[, 1]);\n",
        "  }}\n"
      ))
    } else {
      return(glue(
        "  {{\n",
        "{info_def}",
        "{accumulation}",
        "    for (j in 1:N_{id}) {{\n",
        "      real scaled_info_fisher_s2z = {row_var_j} * square({scale}) * ",
        "{information_j};\n",
        "      {rho}[j, 1] = 1.0 - ",
        "inv(1.0 + scaled_info_fisher_s2z);\n",
        "    }}\n",
        "    {mean_rho}[1] = mean({rho}[, 1]);\n",
        "  }}\n"
      ))
    }
  }

  diagonal <- !is.null(diag_scale)
  stopifnot(xor(!is.null(L), diagonal))

  if (fixed_design) {
    info_def <- glue(
      "    real obs_prec_fisher_s2z = {fisher_info$obs_prec};\n"
    )
    initialization <- accumulation <- ""
    K_j <- if (diagonal) {
      glue(
        "{row_var_j} * obs_prec_fisher_s2z * quad_form_diag(",
        "gram_fisher_s2z_{id}[j], {diag_scale})"
      )
    } else {
      glue(
        "{row_var_j} * obs_prec_fisher_s2z * quad_form(",
        "gram_fisher_s2z_{id}[j], {L})"
      )
    }
  } else {
    info_def <- glue(
      "    array[N_{id}] matrix[M_{id}, M_{id}] info_fisher_s2z;\n"
    )
    initialization <- glue(
      "    for (j in 1:N_{id}) {{\n",
      "      info_fisher_s2z[j] = rep_matrix(0.0, M_{id}, M_{id});\n",
      "    }}\n"
    )
    if (!is.null(fisher_info$joint_info_comp)) {
      accumulation <- fisher_info$joint_info_comp
    } else {
      stopifnot(length(sources) >= 1L)
      accumulation <- paste0(vapply(sources, function(source) {
        stopifnot(
          length(source$columns) == length(source$design_at_n),
          all(source$columns >= 1L & source$columns <= M)
        )
        dynamic_prec <- !is.null(source$obs_prec_at_n)
        obs_prec <- source$obs_prec_at_n %||% source$obs_prec
        definitions <- source$definitions_at_n %||% ""
        design_assign <- cglue(
          "        design_fisher_s2z[{source$columns}] = ",
          "{source$design_at_n};\n"
        )
        glue(
          "    {{\n",
          str_if(
            !dynamic_prec,
            glue("      real obs_prec_fisher_s2z = {obs_prec};\n")
          ),
          "      for (n in 1:{source$N}) {{\n",
          "{definitions}",
          str_if(
            dynamic_prec,
            glue("        real obs_prec_fisher_s2z = {obs_prec};\n")
          ),
          "        vector[M_{id}] design_fisher_s2z = ",
          "zeros_vector(M_{id});\n",
          "{design_assign}",
          "        info_fisher_s2z[{source$group}[n]] += ",
          "obs_prec_fisher_s2z * design_fisher_s2z * ",
          "design_fisher_s2z';\n",
          "      }}\n",
          "    }}\n"
        )
      }, character(1)), collapse = "")
    }
    K_j <- if (diagonal) {
      glue(
        "{row_var_j} * quad_form_diag(info_fisher_s2z[j], ",
        "{diag_scale})"
      )
    } else {
      glue("{row_var_j} * quad_form(info_fisher_s2z[j], {L})")
    }
  }

  if (constrained_s2z) {
    restricted_cov_def <- if (diagonal) {
      ""
    } else {
      glue(
        "      matrix[M_{id}, M_{id}] restricted_post_cov_fisher_s2z = ",
        "({L}) * restricted_white_post_cov_fisher_s2z * ({L})';\n"
      )
    }
    reliability_comp <- if (diagonal) {
      glue(
        "        {rho}[j, k] = fmin(1.0, fmax(0.0, 1.0 -\n",
        "          restricted_white_post_cov_fisher_s2z[k, k] / ",
        "restricted_prior_fraction_fisher_s2z));\n"
      )
    } else {
      glue(
        "        real restricted_prior_var_fisher_s2z = ",
        "restricted_prior_fraction_fisher_s2z * ",
        "prior_var_fisher_s2z[k];\n",
        "        {rho}[j, k] = fmin(1.0, fmax(0.0, 1.0 -\n",
        "          restricted_post_cov_fisher_s2z[k, k] / ",
        "restricted_prior_var_fisher_s2z));\n"
      )
    }
    return(glue(
      "  {{\n",
      "{info_def}",
      "    array[N_{id}] matrix[M_{id}, M_{id}] ",
      "white_post_cov_fisher_s2z;\n",
      "    matrix[M_{id}, M_{id}] sum_white_post_cov_fisher_s2z = ",
      "rep_matrix(0.0, M_{id}, M_{id});\n",
      "    matrix[M_{id}, M_{id}] L_sum_white_post_cov_fisher_s2z;\n",
      "    real restricted_prior_fraction_fisher_s2z = ",
      "1.0 - inv(N_{id});\n",
      str_if(
        !diagonal,
        glue(
          "    vector[M_{id}] prior_var_fisher_s2z = ",
          "rows_dot_self({L});\n"
        )
      ),
      "{initialization}",
      "{accumulation}",
      "    for (j in 1:N_{id}) {{\n",
      "      matrix[M_{id}, M_{id}] K_fisher_s2z = {K_j};\n",
      "      matrix[M_{id}, M_{id}] L_post_precision_fisher_s2z;\n",
      "      matrix[M_{id}, M_{id}] white_factor_fisher_s2z;\n",
      "      K_fisher_s2z = 0.5 * (K_fisher_s2z + K_fisher_s2z');\n",
      "      L_post_precision_fisher_s2z = cholesky_decompose(\n",
      "        add_diag(K_fisher_s2z, 1.0)\n",
      "      );\n",
      "      white_factor_fisher_s2z = mdivide_left_tri_low(\n",
      "        L_post_precision_fisher_s2z, ",
      "diag_matrix(rep_vector(1.0, M_{id}))\n",
      "      );\n",
      "      white_post_cov_fisher_s2z[j] = ",
      "crossprod(white_factor_fisher_s2z);\n",
      "      sum_white_post_cov_fisher_s2z += ",
      "white_post_cov_fisher_s2z[j];\n",
      "    }}\n",
      "    sum_white_post_cov_fisher_s2z = 0.5 * ",
      "(sum_white_post_cov_fisher_s2z + ",
      "sum_white_post_cov_fisher_s2z');\n",
      "    L_sum_white_post_cov_fisher_s2z = ",
      "cholesky_decompose(sum_white_post_cov_fisher_s2z);\n",
      "    for (j in 1:N_{id}) {{\n",
      "      matrix[M_{id}, M_{id}] constraint_factor_fisher_s2z = ",
      "mdivide_left_tri_low(\n",
      "        L_sum_white_post_cov_fisher_s2z, ",
      "white_post_cov_fisher_s2z[j]\n",
      "      );\n",
      "      matrix[M_{id}, M_{id}] ",
      "restricted_white_post_cov_fisher_s2z = ",
      "white_post_cov_fisher_s2z[j] - ",
      "crossprod(constraint_factor_fisher_s2z);\n",
      "{restricted_cov_def}",
      "      for (k in 1:M_{id}) {{\n",
      "        // Compare posterior and prior variances on the S2Z subspace.\n",
      "{reliability_comp}",
      "      }}\n",
      "    }}\n",
      "    for (k in 1:M_{id}) {{\n",
      "      {mean_rho}[k] = mean({rho}[, k]);\n",
      "    }}\n",
      "  }}\n"
    ))
  }

  variance_def <- if (diagonal) {
    glue(
      "      vector[M_{id}] unit_rhs_fisher_s2z = zeros_vector(M_{id});\n"
    )
  } else {
    glue(
      "      matrix[M_{id}, M_{id}] white_factor_fisher_s2z;\n",
      "      row_vector[M_{id}] post_var_fisher_s2z;\n"
    )
  }
  variance_comp <- if (diagonal) {
    ""
  } else {
    glue(
      "      white_factor_fisher_s2z = mdivide_left_tri_low(\n",
      "        L_post_precision_fisher_s2z, sqrt({row_var_j}) * ({L})'\n",
      "      );\n",
      "      post_var_fisher_s2z = columns_dot_self(",
      "white_factor_fisher_s2z);\n"
    )
  }
  reliability_comp <- if (diagonal) {
    glue(
      "      for (k in 1:M_{id}) {{\n",
      "        vector[M_{id}] unit_column_fisher_s2z;\n",
      "        unit_rhs_fisher_s2z[k] = 1.0;\n",
      "        unit_column_fisher_s2z = mdivide_left_tri_low(\n",
      "          L_post_precision_fisher_s2z, unit_rhs_fisher_s2z\n",
      "        );\n",
      "        // The exact ratio is in [0, 1]; clamp roundoff.\n",
      "        {rho}[j, k] = fmin(1.0, fmax(0.0, 1.0 -\n",
      "          dot_self(unit_column_fisher_s2z)));\n",
      "        unit_rhs_fisher_s2z[k] = 0.0;\n",
      "      }}\n"
    )
  } else {
    glue(
      "      for (k in 1:M_{id}) {{\n",
      "        // The exact ratio is in [0, 1]; clamp roundoff.\n",
      "        {rho}[j, k] = fmin(1.0, fmax(0.0, 1.0 -\n",
      "          post_var_fisher_s2z[k] / ",
      "({row_var_j} * prior_var_fisher_s2z[k])));\n",
      "      }}\n"
    )
  }

  glue(
    "  {{\n",
    "{info_def}",
    str_if(
      !diagonal,
      glue("    vector[M_{id}] prior_var_fisher_s2z = rows_dot_self({L});\n")
    ),
    "{initialization}",
    "{accumulation}",
    "    for (j in 1:N_{id}) {{\n",
    "      matrix[M_{id}, M_{id}] K_fisher_s2z = {K_j};\n",
    "      matrix[M_{id}, M_{id}] L_post_precision_fisher_s2z;\n",
    "{variance_def}",
    "      K_fisher_s2z = 0.5 * (K_fisher_s2z + K_fisher_s2z');\n",
    "      L_post_precision_fisher_s2z = cholesky_decompose(\n",
    "        add_diag(K_fisher_s2z, 1.0)\n",
    "      );\n",
    "{variance_comp}",
    "{reliability_comp}",
    "    }}\n",
    "    for (k in 1:M_{id}) {{\n",
    "      {mean_rho}[k] = mean({rho}[, k]);\n",
    "    }}\n",
    "  }}\n"
  )
}

# Precompute the column means of fixed centering fractions once. They enter
# the restricted-transform determinant but do not depend on model parameters.
stan_re_s2z_partial_tdata <- function(out, id) {
  str_add(out$tdata_def) <- glue(
    "  vector[M_{id}] mean_rho_s2z_{id};\n"
  )
  str_add(out$tdata_comp) <- glue(
    "  for (k in 1:M_{id}) {{\n",
    "    mean_rho_s2z_{id}[k] = mean(rho_s2z_{id}[, k]);\n",
    "  }}\n"
  )
  out
}

# Marginal variance of each arithmetic zero-sum deviation under a known
# grouping-level covariance. If P = I - 11' / N and Omega = Lcov Lcov', then
# diag(P Omega P) = rows_dot_self(P Lcov). Relative to the exchangeable
# no-covariance reference, whose restricted marginal variance is 1 - 1 / N,
# the level metric is diag(P Omega P) / (1 - 1 / N). Thus cov = I and an
# omitted cov argument select the same Fisher chart. Lcov is never applied to
# the physical S2Z effects because that would generally destroy their zero sum.
stan_re_s2z_cov_fisher_tdata <- function(out, id) {
  str_add(out$tdata_def) <- glue(
    "  vector<lower=0>[N_{id}] row_var_fisher_s2z_{id};\n"
  )
  str_add(out$tdata_comp) <- glue(
    "  {{\n",
    "    matrix[N_{id}, N_{id}] centered_Lcov_s2z = Lcov_{id};\n",
    "    for (k in 1:N_{id}) centered_Lcov_s2z[, k] -= ",
    "mean(centered_Lcov_s2z[, k]);\n",
    "    row_var_fisher_s2z_{id} = rows_dot_self(centered_Lcov_s2z) * ",
    "N_{id} / (N_{id} - 1.0);\n",
    "  }}\n"
  )
  out
}

# Exact Gaussian normal equations contributed by one conventional group block.
# With a known covariance across grouping levels, the public model is
#
#   Cov(u_g, u_h | scales) = Omega[g, h] A_g A_h',
#
# where A_g = diag(sd_g) L_coef and an optional Student-t mixing scale
# multiplies all coefficients of level g. The likelihood uses the arithmetic
# zero-sum deviations delta. We whiten delta and the M columns mapping the
# omitted conventional mean m = L_ref x through the same structured factor.
# This retains the nonzero mean/contrast score induced by a general Omega and
# avoids constructing or factoring an (N M) by (N M) covariance matrix.
stan_re_s2z_cov_group_comp <- function(id, varying, is_cor, is_student,
                                       delta = glue("r_s2z_{id}"),
                                       has_cov = TRUE) {
  scale_level <- if (varying) {
    glue("sd_level_s2z_{id}")
  } else {
    glue("rep_matrix(sd_{id}', N_{id})")
  }
  L_coef_def <- str_if(
    is_cor,
    glue("    matrix[M_{id}, M_{id}] L_coef_cov_s2z = L_{id};\n")
  )
  student_scale <- str_if(
    is_student,
    glue(
      "    scale_cov_s2z = diag_pre_multiply(",
      "group_scale_s2z_{id}, scale_cov_s2z);\n"
    )
  )
  if (!varying) {
    mean_rhs <- if (is_student) {
      glue("inv(group_scale_s2z_{id})")
    } else {
      glue("rep_vector(1.0, N_{id})")
    }
    white_delta <- if (has_cov) {
      glue(
        "mdivide_left_tri_low(Lcov_{id}, ",
        "coef_white_cov_s2z')"
      )
    } else {
      "coef_white_cov_s2z'"
    }
    one_white <- if (has_cov) {
      glue("mdivide_left_tri_low(Lcov_{id}, {mean_rhs})")
    } else {
      mean_rhs
    }
    coef_white <- if (is_cor) {
      "mdivide_left_tri_low(L_coef_cov_s2z, scaled_delta_cov_s2z')"
    } else {
      "scaled_delta_cov_s2z'"
    }
    return(glue(
      "  {{\n",
      "    matrix[N_{id}, M_{id}] scale_cov_s2z = {scale_level};\n",
      "{L_coef_def}",
      "    matrix[N_{id}, M_{id}] scaled_delta_cov_s2z;\n",
      "    matrix[M_{id}, N_{id}] coef_white_cov_s2z;\n",
      "    matrix[N_{id}, M_{id}] white_delta_cov_s2z;\n",
      "    vector[N_{id}] one_white_cov_s2z;\n",
      "{student_scale}",
      "    scaled_delta_cov_s2z = {delta} ./ scale_cov_s2z;\n",
      "    coef_white_cov_s2z = {coef_white};\n",
      "    white_delta_cov_s2z = {white_delta};\n",
      "    one_white_cov_s2z = {one_white};\n",
      "    group_info_s2z_{id} = dot_self(one_white_cov_s2z);\n",
      "    h_group_s2z_{id} = -white_delta_cov_s2z' * ",
      "one_white_cov_s2z;\n",
      "    group_quad_s2z_{id} = dot_self(",
      "to_vector(white_delta_cov_s2z));\n",
      "  }}\n"
    ))
  }
  if (!is_cor) {
    white_delta <- if (has_cov) {
      glue(
        "mdivide_left_tri_low(Lcov_{id}, scaled_delta_cov_s2z[, k])"
      )
    } else {
      "scaled_delta_cov_s2z[, k]"
    }
    white_basis <- if (has_cov) {
      glue("mdivide_left_tri_low(Lcov_{id}, mean_basis_cov_s2z)")
    } else {
      "mean_basis_cov_s2z"
    }
    return(glue(
      "  group_info_s2z_{id} = zeros_vector(M_{id});\n",
      "  {{\n",
      "    matrix[N_{id}, M_{id}] scale_cov_s2z = {scale_level};\n",
      "    matrix[N_{id}, M_{id}] scaled_delta_cov_s2z;\n",
      "{student_scale}",
      "    scaled_delta_cov_s2z = {delta} ./ scale_cov_s2z;\n",
      "    h_group_s2z_{id} = zeros_vector(M_{id});\n",
      "    group_quad_s2z_{id} = 0.0;\n",
      "    for (k in 1:M_{id}) {{\n",
      "      vector[N_{id}] white_delta_cov_s2z = {white_delta};\n",
      "      vector[N_{id}] mean_basis_cov_s2z = ",
      "rep_vector(reference_sd_s2z_{id}[k], N_{id}) ./ ",
      "scale_cov_s2z[, k];\n",
      "      vector[N_{id}] white_basis_cov_s2z = {white_basis};\n",
      "      group_info_s2z_{id}[k] = dot_self(white_basis_cov_s2z);\n",
      "      h_group_s2z_{id}[k] = -dot_product(",
      "white_basis_cov_s2z, white_delta_cov_s2z);\n",
      "      group_quad_s2z_{id} += dot_self(white_delta_cov_s2z);\n",
      "    }}\n",
      "  }}\n"
    ))
  }
  white_delta <- if (has_cov) {
    glue(
      "mdivide_left_tri_low(Lcov_{id}, ",
      "coef_white_cov_s2z')"
    )
  } else {
    "coef_white_cov_s2z'"
  }
  white_basis <- if (has_cov) {
    glue(
      "mdivide_left_tri_low(Lcov_{id}, ",
      "coef_basis_cov_s2z')"
    )
  } else {
    "coef_basis_cov_s2z'"
  }
  glue(
    "  {{\n",
    "    matrix[N_{id}, M_{id}] scale_cov_s2z = {scale_level};\n",
    "{L_coef_def}",
    "    matrix[N_{id}, M_{id}] scaled_delta_cov_s2z;\n",
    "    matrix[M_{id}, N_{id}] coef_white_cov_s2z;\n",
    "    matrix[N_{id}, M_{id}] white_delta_cov_s2z;\n",
    "    matrix[N_{id} * M_{id}, M_{id}] mean_factor_cov_s2z;\n",
    "{student_scale}",
    "    scaled_delta_cov_s2z = {delta} ./ scale_cov_s2z;\n",
    "    coef_white_cov_s2z = mdivide_left_tri_low(",
    "L_coef_cov_s2z, scaled_delta_cov_s2z');\n",
    "    white_delta_cov_s2z = {white_delta};\n",
    "    for (k in 1:M_{id}) {{\n",
    "      matrix[N_{id}, M_{id}] mean_basis_cov_s2z = ",
    "rep_matrix(L_Sigma_s2z_{id}[, k]', N_{id}) ./ scale_cov_s2z;\n",
    "      matrix[M_{id}, N_{id}] coef_basis_cov_s2z = ",
    "mdivide_left_tri_low(L_coef_cov_s2z, mean_basis_cov_s2z');\n",
    "      matrix[N_{id}, M_{id}] white_basis_cov_s2z = ",
    "{white_basis};\n",
    "      mean_factor_cov_s2z[, k] = to_vector(white_basis_cov_s2z);\n",
    "    }}\n",
    "    P_group_s2z_{id} = crossprod(mean_factor_cov_s2z);\n",
    "    h_group_s2z_{id} = -mean_factor_cov_s2z' * ",
    "to_vector(white_delta_cov_s2z);\n",
    "    group_quad_s2z_{id} = dot_self(to_vector(white_delta_cov_s2z));\n",
    "  }}\n"
  )
}

# Transform orthonormal S2Z coordinates with level- and coefficient-specific
# centering fractions.  Each local Cholesky interpolation is lower triangular
# with a positive diagonal for rho in [0, 1].  Projecting the transformed rows
# back onto the zero-sum subspace preserves the identifying constraint. After
# cancellation with the physical-coordinate scale determinant, the net
# determinant contribution is accumulated in log_det_partial_s2z_ID.
stan_re_s2z_partial_cor_transform <- function(id) {
  glue(
    "  {{\n",
    "    row_vector[M_{id}] mean_partial_s2z = ",
    "zeros_row_vector(M_{id});\n",
    "    log_det_partial_s2z_{id} = 0.0;\n",
    "    for (j in 1:N_{id}) {{\n",
    "      matrix[M_{id}, M_{id}] L_partial_s2z = ",
    "diag_pre_multiply(rho_s2z_{id}[j]', L_Sigma_s2z_{id});\n",
    "      for (k in 1:M_{id}) {{\n",
    "        L_partial_s2z[k, k] += 1.0 - rho_s2z_{id}[j, k];\n",
    "      }}\n",
    "      r_s2z_{id}[j] = (L_Sigma_s2z_{id} * ",
    "mdivide_left_tri_low(L_partial_s2z, r_s2z_{id}[j]'))';\n",
    "      mean_partial_s2z += r_s2z_{id}[j];\n",
    "      log_det_partial_s2z_{id} -= ",
    "sum(log(diagonal(L_partial_s2z)));\n",
    "    }}\n",
    "    mean_partial_s2z /= N_{id};\n",
    "    for (j in 1:N_{id}) {{\n",
    "      r_s2z_{id}[j] -= mean_partial_s2z;\n",
    "    }}\n",
    "    for (k in 1:M_{id}) {{\n",
    "      log_det_partial_s2z_{id} += log(\n",
    "        (1.0 - mean_rho_s2z_{id}[k]) + ",
    "mean_rho_s2z_{id}[k] * L_Sigma_s2z_{id}[k, k]\n",
    "      );\n",
    "    }}\n",
    "  }}\n"
  )
}

# Diagonal specialization of the same restricted transform.  This is used by
# scalar and independent blocks and avoids all per-level matrix operations.
stan_re_s2z_partial_independent_transform <- function(
    id, k, r_s2z, scale,
    z_s2z = glue(
      "segment(z_s2z_{id}, ({k} - 1) * (N_{id} - 1) + 1, N_{id} - 1)"
    )) {
  glue(
    "  {{\n",
    "    vector[N_{id}] centered_partial_s2z = ",
    "sum_to_zero_constrain_brms({z_s2z});\n",
    "    vector[N_{id}] scale_partial_s2z = 1.0 - ",
    "rho_s2z_{id}[, {k}] + rho_s2z_{id}[, {k}] * {scale};\n",
    "    centered_partial_s2z = {scale} * centered_partial_s2z ./ ",
    "scale_partial_s2z;\n",
    "    {r_s2z} = centered_partial_s2z - mean(centered_partial_s2z);\n",
    "    log_det_partial_s2z_{id} += -sum(log(scale_partial_s2z));\n",
    "    log_det_partial_s2z_{id} += log(\n",
    "      (1.0 - mean_rho_s2z_{id}[{k}]) + ",
    "mean_rho_s2z_{id}[{k}] * {scale}\n",
    "    );\n",
    "  }}\n"
  )
}

# The separated Gaussian mean/contrast law permits a Matheron update in the
# population-coordinate dimension instead of a complete square in the sum of
# all block dimensions. Conditional Student-t population priors remain
# eligible because their scale-mixture variables make them proper Gaussians
# here. Student-t group effects and group-varying scales do not separate their
# omitted means from their contrasts and must use the general joint system.
.stan_re_s2z_joint_matheron_info <- function(set) {
  infos <- set$infos
  if (length(infos) <= 1L || isTRUE(infos[[1L]]$ordinal)) {
    return(NULL)
  }
  separated <- all(vapply(infos, function(info) {
    all(info$r$dist == "gaussian") && all(info$r$scale == "shared") &&
      all(!nzchar(info$r$cov))
  }, logical(1)))
  if (!separated) {
    return(NULL)
  }
  # Only rows reached by an omitted group mean can couple blocks. A flat
  # population prior removes its row from that conditioning system, while a
  # proper row outside the group design remains an ordinary independent prior.
  active <- unique(unlist(lapply(infos, function(info) {
    rows <- info$match_q
    if (info$center && any(info$r$coef != "Intercept")) {
      rows <- c(rows, 1L)
    }
    rows
  }), use.names = FALSE))
  proper <- which(vapply(infos[[1L]]$prior, function(spec) {
    spec$dist %in% c("normal", "student")
  }, logical(1)))
  P <- sort(intersect(active, proper))
  total_M <- sum(vapply(infos, function(info) nrow(info$r), integer(1)))
  # The general complete square is already preferable at equal or lower
  # dimension. Strict inequality makes this dispatch a genuine fast path.
  if (length(P) >= total_M) {
    return(NULL)
  }
  list(P = P, inactive = setdiff(proper, P), total_M = total_M)
}

.stan_re_s2z_joint_uses_matheron <- function(set) {
  !is.null(.stan_re_s2z_joint_matheron_info(set))
}

# Add optional prior factors to selected realized group-level scales. These
# scales remain deterministic functions of the baseline scale, log-scale
# heterogeneity, and standardized level innovations, so evaluating a density
# at sd_level does not require a change-of-variables Jacobian. The factors are
# additive to (rather than replacements for) the exchangeable scale hierarchy.
.stan_re_s2z_sd_level_prior <- function(id, r, prior, normalize) {
  if (!"level" %in% names(prior)) {
    # Backward compatibility for serialized brmsprior objects created before
    # level-addressable priors were introduced.
    return("")
  }
  px <- check_prefix(r)
  level_prior <- subset2(
    prior, class = "sd_level", group = r$group[1],
    coef = r$coef, level = get_levels(r)[[r$group[1]]], ls = px
  )
  level_prior <- level_prior[nzchar(level_prior$prior), , drop = FALSE]
  if (!nrow(level_prior)) {
    return("")
  }
  levels <- get_levels(r)[[r$group[1]]]
  level_index <- match(level_prior$level, levels)
  r_key <- paste(
    combine_prefix(check_prefix(r)), r$coef, sep = "\r"
  )
  prior_key <- paste(
    combine_prefix(check_prefix(level_prior)), level_prior$coef, sep = "\r"
  )
  coef_index <- match(prior_key, r_key)
  if (anyNA(coef_index)) {
    # A response-free selector remains unambiguous for ordinary local blocks.
    # Cross-response IDs can repeat coefficient names, in which case the user
    # must qualify the realized-scale prior by response.
    for (i in which(is.na(coef_index))) {
      candidates <- which(r$coef == level_prior$coef[i])
      if (length(candidates) == 1L) {
        coef_index[i] <- candidates
      }
    }
  }
  if (anyNA(level_index) || anyNA(coef_index)) {
    stop2("A realized group-level scale prior in a cross-response block ",
          "must identify an unambiguous response and coefficient.")
  }
  selector <- paste(level_index, coef_index, sep = "\r")
  if (anyDuplicated(selector)) {
    stop2("Duplicated priors on a realized group-level scale are not allowed.")
  }
  out <- paste0(
    "  // additional priors on realized group-level standard deviations\n"
  )
  for (i in seq_rows(level_prior)) {
    par <- glue(
      "sd_level_s2z_{id}[{level_index[i]}, {coef_index[i]}]"
    )
    target <- stan_target_prior(
      level_prior$prior[i], par = par, bound = "<lower=0>",
      resp = level_prior$resp[i], normalize = normalize
    )
    str_add(out) <- paste0(lpp(), target, ";\n")
  }
  out
}

# A logistic population prior cannot be integrated by the conditional
# Gaussian omitted-mean kernels. Such systems retain the omitted means
# explicitly in standardized coordinates and score the exact requested prior.
.stan_re_s2z_uses_explicit_mean <- function(set) {
  any(vapply(set$infos, function(info) {
    !is.null(info$prior) && any(vapply(
      info$prior, function(spec) identical(spec$dist, "logistic"), logical(1)
    ))
  }, logical(1)))
}

# Exact scalar population-prior statement for the explicit-mean path.
stan_re_s2z_prior_target <- function(spec, par, normalize) {
  stopifnot(is.list(spec), is.character(par), length(par) == 1L)
  if (identical(spec$dist, "flat")) {
    return("")
  }
  lpdf <- stan_lpdf_name(normalize)
  location <- stan_s2z_arg_code(spec$location)
  scale <- stan_s2z_arg_code(spec$scale)
  if (identical(spec$dist, "normal")) {
    glue("  lprior += normal_{lpdf}({par} | {location}, {scale});\n")
  } else if (identical(spec$dist, "student")) {
    df <- stan_s2z_arg_code(spec$df)
    glue(
      "  lprior += student_t_{lpdf}({par} | {df}, {location}, {scale});\n"
    )
  } else if (identical(spec$dist, "logistic")) {
    glue("  lprior += logistic_{lpdf}({par} | {location}, {scale});\n")
  } else {
    stop2("Internal error: unsupported explicit S2Z population prior.")
  }
}

# One block in an exact explicit-mean system. The sampled mean coordinate is
# v = sqrt(N) L^-1 m, so m = L v / sqrt(N). With orthonormal S2Z contrasts,
# the Jacobian back to conventional group effects is |L|.
.stan_re_s2z_explicit_block <- function(id, set, prior, normalize,
                                        out = list(), fisher_info = NULL) {
  if (is.null(out[["tpar_prior"]])) {
    out[["tpar_prior"]] <- ""
  }
  take_info <- match(id, vapply(set$infos, `[[`, numeric(1), "id"))
  stopifnot(!is.na(take_info))
  info <- set$infos[[take_info]]
  r <- info$r
  q <- length(info$qnames)
  M <- nrow(r)
  J <- seq_rows(r)
  p <- info$p
  px <- check_prefix(r)
  idp <- paste0(r$id, usc(combine_prefix(px)))
  r_s2z <- glue("r_s2z_{idp}_{r$cn}")
  is_cor <- M > 1L && isTRUE(r$cor[1])
  is_student <- identical(r$dist[1], "student")
  s2z_mode <- re_s2z_center_mode(r)
  s2z_center <- identical(s2z_mode, "centered")
  s2z_fisher <- !is.null(fisher_info)
  s2z_partial <- s2z_mode %in% c("partial", "auto")
  stopifnot(!s2z_fisher || identical(s2z_mode, "auto"))

  if (s2z_fisher) {
    out <- stan_re_s2z_fisher_def(out, id)
    if (isTRUE(fisher_info$fixed_design)) {
      out <- stan_re_s2z_fisher_tdata(out, id, r, fisher_info)
    }
  } else if (s2z_partial) {
    str_add(out$data) <- glue(
      "  matrix<lower=0,upper=1>[N_{id}, M_{id}] rho_s2z_{id};",
      "  // fixed numeric centering fractions\n"
    )
    out <- stan_re_s2z_partial_tdata(out, id)
  }

  if (is_cor) {
    str_add(out$data) <- glue(
      "  int<lower=1> NC_{id};  // number of group-level correlations\n"
    )
    str_add_list(out) <- stan_prior(
      prior, class = "L", group = r$group[1], suffix = usc(id),
      type = glue("cholesky_factor_corr[M_{id}]"),
      comment = "cholesky factor of correlation matrix",
      normalize = normalize
    )
  }
  if (identical(id, set$set_id)) {
    str_add(out$fun) <- "  #include 'fun_sum_to_zero.stan'\n"
  }

  str_add(out$par) <- glue(
    "  vector[M_{id} * (N_{id} - 1)] z_s2z_{id};",
    if (s2z_partial) {
      "  // partially centered orthonormal S2Z coordinates\n"
    } else if (s2z_center) {
      "  // physical orthonormal S2Z coordinates\n"
    } else {
      "  // standardized orthonormal S2Z coordinates\n"
    },
    "  vector[M_{id}] z_mean_s2z_{id};",
    "  // standardized omitted group mean\n"
  )
  str_add(out$tpar_def) <- glue(
    "  // exact explicit-mean S2Z block {id}\n",
    "  matrix[N_{id}, M_{id}] r_s2z_{id};\n",
    "  matrix[M_{id}, M_{id}] L_Sigma_s2z_{id};\n",
    "  vector[M_{id}] mean_r_s2z_{id};\n",
    "  real<lower=0> group_quad_s2z_{id};\n",
    str_if(
      is_student,
      glue(
        "  vector<lower=0>[N_{id}] group_scale_s2z_{id};\n",
        "  vector<lower=0>[N_{id}] group_prec_s2z_{id};\n"
      )
    ),
    str_if(
      s2z_partial,
      glue("  real log_det_partial_s2z_{id};\n")
    ),
    cglue("  vector[N_{id}] {r_s2z};\n")
  )

  # H is fixed by the formula and centered predictor means. Build it once in
  # transformed data so it never enters the reverse-mode expression graph.
  str_add_list(out) <- stan_re_s2z_H_code(info)
  if (is_cor) {
    str_add(out$tpar_comp) <- glue(
      "  L_Sigma_s2z_{id} = diag_pre_multiply(sd_{id}, L_{id});\n"
    )
  } else {
    str_add(out$tpar_comp) <- glue(
      "  L_Sigma_s2z_{id} = diag_matrix(sd_{id});\n"
    )
  }
  if (s2z_fisher) {
    str_add(out$gen_comp) <- stan_re_s2z_fisher_gq_comp(
      id, r = r, fisher_info = fisher_info,
      L = if (is_cor) glue("L_Sigma_s2z_{id}") else NULL,
      scale = if (M == 1L) glue("sd_{id}[1]") else NULL,
      diag_scale = if (!is_cor && M > 1L) glue("sd_{id}") else NULL
    )
  }
  if (is_cor || !s2z_partial) {
    str_add(out$tpar_comp) <- glue(
      "  for (k in 1:M_{id}) {{\n",
      "    r_s2z_{id}[, k] = sum_to_zero_constrain_brms(",
      "segment(z_s2z_{id}, (k - 1) * (N_{id} - 1) + 1, ",
      "N_{id} - 1));\n",
      "  }}\n"
    )
  }
  if (s2z_partial) {
    if (is_cor) {
      str_add(out$tpar_comp) <- stan_re_s2z_partial_cor_transform(id)
    } else {
      str_add(out$tpar_comp) <- glue(
        "  log_det_partial_s2z_{id} = 0.0;\n"
      )
      for (k in seq_len(M)) {
        str_add(out$tpar_comp) <- stan_re_s2z_partial_independent_transform(
          id, k, glue("r_s2z_{id}[, {k}]"), glue("sd_{id}[{k}]")
        )
      }
    }
  } else if (!s2z_center) {
    if (is_cor) {
      str_add(out$tpar_comp) <- glue(
        "  r_s2z_{id} = r_s2z_{id} * L_Sigma_s2z_{id}';\n"
      )
    } else {
      for (k in seq_len(M)) {
        str_add(out$tpar_comp) <- glue(
          "  r_s2z_{id}[, {k}] *= sd_{id}[{k}];\n"
        )
      }
    }
  }
  str_add(out$tpar_comp) <- glue(
    "  mean_r_s2z_{id} = L_Sigma_s2z_{id} * z_mean_s2z_{id} / ",
    "sqrt(1.0 * N_{id});\n"
  )
  if (is_student) {
    tr <- subset_reframe_dist(r, "student")
    g <- usc(tr$ggn[1])
    str_add(out$tpar_comp) <- glue(
      "  group_scale_s2z_{id} = dfm{g};\n",
      "  group_prec_s2z_{id} = inv_square(group_scale_s2z_{id});\n"
    )
  }
  explicit_group_quad_code <- if (is_student) {
    glue(
      "    vector[M_{id}] mean_white_s2z = z_mean_s2z_{id} / ",
      "sqrt(1.0 * N_{id});\n",
      "    group_quad_s2z_{id} = 0.0;\n",
      "    for (j in 1:N_{id}) {{\n",
      "      vector[M_{id}] white_level_s2z = ",
      "white_group_s2z[, j] + mean_white_s2z;\n",
      "      group_quad_s2z_{id} += group_prec_s2z_{id}[j] * ",
      "dot_self(white_level_s2z);\n",
      "    }}\n"
    )
  } else {
    glue(
      "    group_quad_s2z_{id} = dot_self(to_vector(white_group_s2z)) ",
      "+ dot_self(z_mean_s2z_{id});\n"
    )
  }
  str_add(out$tpar_comp) <- glue(
    "  {{\n",
    "    matrix[M_{id}, N_{id}] white_group_s2z = ",
    "mdivide_left_tri_low(L_Sigma_s2z_{id}, r_s2z_{id}');\n",
    "{explicit_group_quad_code}",
    "  }}\n",
    cglue("  {r_s2z} = r_s2z_{id}[, {J}];\n")
  )
  str_add(out$pll_args) <- cglue(", vector {r_s2z}")

  str_add(out$tpar_prior) <- glue(
    "  lprior += -0.5 * group_quad_s2z_{id}\n",
    str_if(
      s2z_center,
      glue(
        "    - (N_{id} - 1) * ",
        "sum(log(diagonal(L_Sigma_s2z_{id})))\n"
      )
    ),
    str_if(
      s2z_partial,
      glue("    + log_det_partial_s2z_{id}\n")
    ),
    str_if(
      is_student,
      glue("    - M_{id} * sum(log(group_scale_s2z_{id}))\n")
    ),
    str_if(
      normalize,
      glue("    - 0.5 * N_{id} * M_{id} * log(2 * pi())\n")
    ),
    "  ;\n"
  )
  out
}

# Score and reconstruct one complete explicit-mean system after its blocks
# have been generated. This path is shared by singleton and multiblock models.
.stan_re_s2z_explicit_system <- function(set, normalize, out = list()) {
  infos <- set$infos
  info <- infos[[1L]]
  q <- length(info$qnames)
  set_id <- set$set_id
  p <- info$p
  active_q <- info$active_q
  if (is.null(active_q)) {
    active_q <- which(vapply(
      info$prior, function(spec) !identical(spec$dist, "flat"), logical(1)
    ))
  }

  str_add(out$tpar_def) <- glue(
    "  vector[{q}] q_explicit_s2z_{set_id};",
    "  // conventional coefficients for exact explicit S2Z means\n"
  )
  str_add(out$tpar_comp) <- glue(
    "  q_explicit_s2z_{set_id} = theta_s2z{p};\n"
  )
  for (block in infos) {
    str_add(out$tpar_comp) <- glue(
      "  q_explicit_s2z_{set_id} -= H_s2z_{block$id} * ",
      "mean_r_s2z_{block$id};\n"
    )
  }
  for (k in active_q) {
    str_add(out$tpar_prior) <- stan_re_s2z_prior_target(
      info$prior[[k]], glue("q_explicit_s2z_{set_id}[{k}]"), normalize
    )
  }

  str_add(out$gen_def) <- glue(
    "  vector[{q}] q_recovered_s2z_{set_id};\n"
  )
  if (info$center) {
    str_add(out$gen_def) <- glue(
      "  real Intercept{p};\n",
      str_if(length(info$fixef), glue("  vector[Kc{p}] b{p};\n")),
      "  real b{p}_Intercept;\n"
    )
  } else {
    str_add(out$gen_def) <- glue("  vector[{q}] b{p};\n")
  }
  for (block in infos) {
    r <- block$r
    id <- block$id
    J <- seq_rows(r)
    idp <- paste0(r$id, usc(combine_prefix(check_prefix(r))))
    r_s2z <- glue("r_s2z_{idp}_{r$cn}")
    r_public <- glue("r_{idp}_{r$cn}")
    is_cor <- nrow(r) > 1L && isTRUE(r$cor[1])
    if (is_cor) {
      str_add(out$gen_def) <- glue(
        "  matrix[N_{id}, M_{id}] r_{id};\n",
        "  corr_matrix[M_{id}] Cor_{id}",
        " = multiply_lower_tri_self_transpose(L_{id});\n",
        "  vector<lower=-1,upper=1>[NC_{id}] cor_{id};\n"
      )
    }
    str_add(out$gen_def) <- cglue(
      "  vector[N_{id}] {r_public};\n"
    )
    if (is_cor) {
      str_add(out$gen_comp) <- glue(
        "  r_{id} = r_s2z_{id};\n",
        "  for (j in 1:N_{id}) r_{id}[j] += mean_r_s2z_{id}';\n",
        cglue("  {r_public} = r_{id}[, {J}];\n")
      )
    } else {
      str_add(out$gen_comp) <- cglue(
        "  {r_public} = {r_s2z} + mean_r_s2z_{id}[{J}];\n"
      )
    }
  }
  str_add(out$gen_comp) <- glue(
    "  q_recovered_s2z_{set_id} = q_explicit_s2z_{set_id};\n"
  )
  for (block in infos) {
    r <- block$r
    id <- block$id
    if (nrow(r) > 1L && isTRUE(r$cor[1])) {
      str_add(out$gen_comp) <- stan_cor_gen_comp(
        cor = glue("cor_{id}"), ncol = glue("M_{id}")
      )
    }
  }
  if (info$center) {
    str_add(out$gen_comp) <- glue(
      "  Intercept{p} = q_recovered_s2z_{set_id}[1];\n",
      str_if(
        length(info$fixef),
        glue(
          "  b{p} = tail(q_recovered_s2z_{set_id}, Kc{p});\n",
          "  b{p}_Intercept = Intercept{p} - dot_product(",
          "means_X{p}, b{p});\n"
        )
      ),
      str_if(
        !length(info$fixef), glue("  b{p}_Intercept = Intercept{p};\n")
      )
    )
  } else {
    str_add(out$gen_comp) <- glue(
      "  b{p} = q_recovered_s2z_{set_id};\n"
    )
  }
  out
}

.stan_re_s2z_explicit <- function(id, bframe, prior, normalize,
                                  out = list(), fisher_info = NULL,
                                  stanvars = NULL) {
  r <- subset2(bframe$frame$re, id = id)
  bfl <- re_s2z_bframel(bframe, r)
  infos <- re_s2z_infos(
    bfl, prior = prior, stanvars = stanvars
  )
  stopifnot(length(infos) == 1L)
  set <- nlist(set_id = id, ids = id, bfl, infos)
  out <- .stan_re_s2z_explicit_block(
    id, set = set, prior = prior, normalize = normalize, out = out,
    fisher_info = fisher_info
  )
  .stan_re_s2z_explicit_system(set, normalize = normalize, out = out)
}

# Emit the deterministic omitted-mean map for one covariance block. Ordinal
# descriptors carry the fully validated affine map because their threshold
# rows are not ordinary population-intercept rows. Legacy descriptors retain
# the established name/centering construction.
stan_re_s2z_H_code <- function(info) {
  stopifnot(is.list(info), is.reframe(info$r), has_rows(info$r))
  out <- list()
  id <- info$id
  q <- length(info$qnames)
  M <- nrow(info$r)
  str_add(out$tdata_def) <- glue(
    "  matrix[{q}, M_{id}] H_s2z_{id};\n"
  )
  str_add(out$tdata_comp) <- glue(
    "  H_s2z_{id} = rep_matrix(0.0, {q}, M_{id});\n"
  )
  if (!is.null(info$H)) {
    H <- as.matrix(info$H)
    stopifnot(identical(dim(H), c(q, M)))
    nonzero <- which(H != 0, arr.ind = TRUE)
    if (nrow(nonzero)) {
      for (k in seq_len(nrow(nonzero))) {
        i <- nonzero[k, 1L]
        j <- nonzero[k, 2L]
        str_add(out$tdata_comp) <- glue(
          "  H_s2z_{id}[{i}, {j}] = {stan_s2z_number(H[i, j])};\n"
        )
      }
    }
    return(out)
  }
  for (j in seq_len(M)) {
    qi <- info$match_q[j]
    str_add(out$tdata_comp) <- glue(
      "  H_s2z_{id}[{qi}, {j}] = 1.0;\n"
    )
    if (info$center && info$r$coef[j] != "Intercept") {
      str_add(out$tdata_comp) <- glue(
        "  H_s2z_{id}[1, {j}] = means_X{info$p}[{qi - 1L}];\n"
      )
    }
  }
  out
}

# Stan code local to one covariance block participating in a joint S2Z
# omitted-mean system. The block contributes its physical zero-sum effects and
# the Gaussian normal equations for its omitted mean in reference-whitened
# coordinates. Population-level priors and mean recovery are handled once by
# .stan_re_s2z_joint().
.stan_re_s2z_joint_block <- function(id, set, bframe, prior, threads,
                                     normalize, out = list(),
                                     fisher_info = NULL, stanvars = NULL) {
  lpdf <- stan_lpdf_name(normalize)
  if (is.null(out[["tpar_prior"]])) {
    out[["tpar_prior"]] <- ""
  }
  if (.stan_re_s2z_uses_explicit_mean(set)) {
    return(.stan_re_s2z_explicit_block(
      id, set = set, prior = prior, normalize = normalize, out = out,
      fisher_info = fisher_info
    ))
  }
  r <- subset2(bframe$frame$re, id = id)
  stopifnot(is.reframe(r), has_rows(r), all(r$s2z), id %in% set$ids)
  bfl <- set$bfl
  info <- re_s2z_info(
    bfl, prior = prior, id = id, stanvars = stanvars
  )
  q <- length(info$qnames)
  M <- nrow(r)
  J <- seq_rows(r)
  px <- check_prefix(r)
  idp <- paste0(r$id, usc(combine_prefix(px)))
  r_s2z <- glue("r_s2z_{idp}_{r$cn}")
  is_cor <- M > 1L && isTRUE(r$cor[1])
  is_student <- identical(r$dist[1], "student")
  has_cov <- nzchar(r$cov[1])
  varying <- identical(r$scale[1], "varying")
  group_precision_kind <- stan_re_s2z_group_precision_kind(
    varying = varying, is_cor = is_cor
  )
  use_matheron <- .stan_re_s2z_joint_uses_matheron(set)
  s2z_mode <- re_s2z_center_mode(r)
  s2z_center <- identical(s2z_mode, "centered")
  s2z_fisher <- !is.null(fisher_info)
  s2z_partial <- s2z_mode %in% c("partial", "auto")
  stopifnot(!s2z_fisher || identical(s2z_mode, "auto"))

  if (s2z_fisher) {
    out <- stan_re_s2z_fisher_def(out, id)
    if (isTRUE(fisher_info$fixed_design)) {
      out <- stan_re_s2z_fisher_tdata(out, id, r, fisher_info)
    }
    if (has_cov) {
      out <- stan_re_s2z_cov_fisher_tdata(out, id)
    }
  } else if (s2z_partial) {
    str_add(out$data) <- glue(
      "  matrix<lower=0,upper=1>[N_{id}, M_{id}] rho_s2z_{id};",
      "  // fixed numeric centering fractions\n"
    )
    out <- stan_re_s2z_partial_tdata(out, id)
  }
  if (varying) {
    str_add_list(out) <- stan_prior(
      prior, class = "sdlog", group = r$group[1], coef = r$coef,
      type = glue("vector[M_{id}]"), suffix = glue("_{id}"), px = px,
      comment = "SDs of group-level log-scale deviations",
      normalize = normalize
    )
  }
  if (is_cor) {
    str_add(out$data) <- glue(
      "  int<lower=1> NC_{id};  // number of group-level correlations\n"
    )
    str_add_list(out) <- stan_prior(
      prior, class = "L", group = r$group[1], suffix = usc(id),
      type = glue("cholesky_factor_corr[M_{id}]"),
      comment = "cholesky factor of correlation matrix",
      normalize = normalize
    )
  }
  if (identical(id, set$set_id)) {
    str_add(out$fun) <- "  #include 'fun_sum_to_zero.stan'\n"
  }

  str_add(out$par) <- glue(
    "  vector[M_{id} * (N_{id} - 1)] z_s2z_{id};",
    if (s2z_partial) {
      "  // partially centered orthonormal S2Z coordinates\n"
    } else if (s2z_center) {
      "  // physical orthonormal S2Z coordinates\n"
    } else {
      "  // standardized orthonormal S2Z coordinates\n"
    }
  )
  if (varying) {
    str_add(out$par) <- glue(
      "  vector[M_{id} * N_{id}] z_sd_s2z_{id};",
      "  // flattened orthonormal log-scale coordinates: contrasts, then means\n"
    )
    str_add(out$model_prior) <- glue(
      "  target += std_normal_{lpdf}(z_sd_s2z_{id});\n"
    )
  }

  str_add(out$tpar_def) <- glue(
    "  // S2Z block {id} in a joint omitted-mean system\n",
    "  matrix[N_{id}, M_{id}] r_s2z_{id};\n",
    "  matrix[M_{id}, M_{id}] L_Sigma_s2z_{id};\n",
    str_if(
      !use_matheron,
      glue(
        stan_re_s2z_group_precision_def(id, group_precision_kind),
        "  vector[M_{id}] h_group_s2z_{id};\n"
      )
    ),
    "  real<lower=0> group_quad_s2z_{id};\n",
    str_if(s2z_partial, glue("  real log_det_partial_s2z_{id};\n")),
    str_if(
      is_student && !varying,
      glue(
        "  vector<lower=0>[N_{id}] group_scale_s2z_{id};\n",
        "  vector<lower=0>[N_{id}] group_prec_s2z_{id};\n"
      )
    ),
    str_if(
      varying,
      glue(
        "  matrix<lower=0>[N_{id}, M_{id}] sd_level_s2z_{id};\n",
        "  vector<lower=0>[M_{id}] reference_sd_s2z_{id};\n",
        str_if(
          is_student,
          glue(
            "  vector<lower=0>[N_{id}] group_scale_s2z_{id};\n",
            "  vector<lower=0>[N_{id}] group_prec_s2z_{id};\n"
          )
        )
      )
    ),
    "  // vectors used by the observation-level likelihood\n",
    cglue("  vector[N_{id}] {r_s2z};\n")
  )

  # H depends only on the checked formula/design mapping. Keeping it in
  # transformed data avoids rebuilding an autodiff matrix at every gradient.
  str_add_list(out) <- stan_re_s2z_H_code(info)

  if (varying) {
    str_add(out$tpar_comp) <- glue(
      "  reference_sd_s2z_{id} = sd_{id} .* exp(sdlog_{id} .* ",
      "tail(z_sd_s2z_{id}, M_{id}) / sqrt(1.0 * N_{id}));\n",
      "  for (k in 1:M_{id}) {{\n",
      "    vector[N_{id}] z_sd_centered_s2z = ",
      "sum_to_zero_constrain_brms(segment(z_sd_s2z_{id}, ",
      "(k - 1) * (N_{id} - 1) + 1, N_{id} - 1));\n",
      "    sd_level_s2z_{id}[, k] = reference_sd_s2z_{id}[k] * ",
      "exp(sdlog_{id}[k] * z_sd_centered_s2z);\n",
      "  }}\n"
    )
    if (is_student) {
      tr <- subset_reframe_dist(r, "student")
      g <- usc(tr$ggn[1])
      str_add(out$tpar_comp) <- glue(
        "  group_scale_s2z_{id} = dfm{g};\n",
        "  group_prec_s2z_{id} = inv_square(group_scale_s2z_{id});\n"
      )
    }
  } else if (is_student) {
    tr <- subset_reframe_dist(r, "student")
    g <- usc(tr$ggn[1])
    str_add(out$tpar_comp) <- glue(
      "  group_scale_s2z_{id} = dfm{g};\n",
      "  group_prec_s2z_{id} = inv_square(group_scale_s2z_{id});\n"
    )
  }
  if (varying) {
    str_add(out$tpar_prior) <- .stan_re_s2z_sd_level_prior(
      id, r = r, prior = prior, normalize = normalize
    )
  }

  scale <- if (varying) glue("reference_sd_s2z_{id}") else glue("sd_{id}")
  if (is_cor) {
    str_add(out$tpar_comp) <- glue(
      "  L_Sigma_s2z_{id} = diag_pre_multiply({scale}, L_{id});\n"
    )
  } else {
    str_add(out$tpar_comp) <- glue(
      "  L_Sigma_s2z_{id} = diag_matrix({scale});\n"
    )
  }

  if (s2z_fisher) {
    str_add(out$gen_comp) <- stan_re_s2z_fisher_gq_comp(
      id, r = r, fisher_info = fisher_info,
      L = if (is_cor) glue("L_Sigma_s2z_{id}") else NULL,
      scale = if (M == 1L) glue("{scale}[1]") else NULL,
      row_var = if (has_cov) glue("row_var_fisher_s2z_{id}[j]") else NULL,
      diag_scale = if (!is_cor && M > 1L) scale else NULL
    )
  }

  if (s2z_partial && is_cor) {
    partial_transform <- stan_re_s2z_partial_cor_transform(id)
    str_add(out$tpar_comp) <- glue(
      "  for (k in 1:M_{id}) {{\n",
      "    r_s2z_{id}[, k] = sum_to_zero_constrain_brms(",
      "segment(z_s2z_{id}, (k - 1) * (N_{id} - 1) + 1, ",
      "N_{id} - 1));\n",
      "  }}\n",
      "{partial_transform}"
    )
  } else if (s2z_partial) {
    str_add(out$tpar_comp) <- glue(
      "  log_det_partial_s2z_{id} = 0.0;\n"
    )
    for (j in seq_len(M)) {
      scale_j <- if (varying) {
        glue("reference_sd_s2z_{id}[{j}]")
      } else {
        glue("sd_{id}[{j}]")
      }
      str_add(out$tpar_comp) <- stan_re_s2z_partial_independent_transform(
        id, j, glue("r_s2z_{id}[, {j}]"), scale_j
      )
    }
  } else {
    str_add(out$tpar_comp) <- glue(
      "  for (k in 1:M_{id}) {{\n",
      "    r_s2z_{id}[, k] = sum_to_zero_constrain_brms(",
      "segment(z_s2z_{id}, (k - 1) * (N_{id} - 1) + 1, ",
      "N_{id} - 1));\n",
      "  }}\n"
    )
    if (!s2z_center) {
      if (is_cor) {
        str_add(out$tpar_comp) <- glue(
          "  r_s2z_{id} = r_s2z_{id} * L_Sigma_s2z_{id}';\n"
        )
      } else {
        str_add(out$tpar_comp) <- glue(
          "  for (k in 1:M_{id}) {{\n",
          "    r_s2z_{id}[, k] *= {scale}[k];\n",
          "  }}\n"
        )
      }
    }
  }

  if (has_cov) {
    str_add(out$tpar_comp) <- stan_re_s2z_cov_group_comp(
      id, varying = varying, is_cor = is_cor, is_student = is_student
    )
  } else if (use_matheron && !s2z_center && !s2z_partial) {
    # In non-centered coordinates the S2Z basis and coefficient transform are
    # already orthonormal/whitened, so no N by M solve or division is needed.
    str_add(out$tpar_comp) <- glue(
      "  group_quad_s2z_{id} = dot_self(z_s2z_{id});\n"
    )
  } else if (use_matheron && is_cor) {
    str_add(out$tpar_comp) <- glue(
      "  {{\n",
      "    matrix[M_{id}, N_{id}] white_group_s2z = ",
      "mdivide_left_tri_low(L_Sigma_s2z_{id}, r_s2z_{id}');\n",
      "    group_quad_s2z_{id} = dot_self(to_vector(white_group_s2z));\n",
      "  }}\n"
    )
  } else if (use_matheron) {
    str_add(out$tpar_comp) <- glue(
      "  {{\n",
      "    group_quad_s2z_{id} = 0.0;\n",
      "    for (k in 1:M_{id}) {{\n",
      "      vector[N_{id}] white_group_s2z = ",
      "r_s2z_{id}[, k] / sd_{id}[k];\n",
      "      group_quad_s2z_{id} += dot_self(white_group_s2z);\n",
      "    }}\n",
      "  }}\n"
    )
  } else if (varying && is_cor) {
    level_weight <- str_if(
      is_student, glue("group_prec_s2z_{id}[j] * ")
    )
    str_add(out$tpar_comp) <- glue(
      "  P_group_s2z_{id} = rep_matrix(0.0, M_{id}, M_{id});\n",
      "  h_group_s2z_{id} = zeros_vector(M_{id});\n",
      "  group_quad_s2z_{id} = 0.0;\n",
      "  for (j in 1:N_{id}) {{\n",
      "    matrix[M_{id}, M_{id}] L_level_s2z = ",
      "diag_pre_multiply(sd_level_s2z_{id}[j]', L_{id});\n",
      "    matrix[M_{id}, M_{id}] relative_precision_s2z = ",
      "mdivide_left_tri_low(L_level_s2z, L_Sigma_s2z_{id});\n",
      "    vector[M_{id}] white_level_s2z = ",
      "mdivide_left_tri_low(L_level_s2z, r_s2z_{id}[j]');\n",
      "    P_group_s2z_{id} += {level_weight}",
      "crossprod(relative_precision_s2z);\n",
      "    h_group_s2z_{id} -= {level_weight}",
      "relative_precision_s2z' * white_level_s2z;\n",
      "    group_quad_s2z_{id} += {level_weight}",
      "dot_self(white_level_s2z);\n",
      "  }}\n"
    )
  } else if (varying) {
    # A diagonal block requires only elementwise work at each level. This is
    # important when several high-dimensional independent blocks are combined.
    level_weight <- str_if(
      is_student, glue("group_prec_s2z_{id}[j] * ")
    )
    str_add(out$tpar_comp) <- glue(
      "  group_info_s2z_{id} = zeros_vector(M_{id});\n",
      "  {{\n",
      "    h_group_s2z_{id} = zeros_vector(M_{id});\n",
      "    group_quad_s2z_{id} = 0.0;\n",
      "    for (j in 1:N_{id}) {{\n",
      "      vector[M_{id}] relative_precision_s2z = ",
      "reference_sd_s2z_{id} ./ sd_level_s2z_{id}[j]';\n",
      "      vector[M_{id}] white_level_s2z = ",
      "r_s2z_{id}[j]' ./ sd_level_s2z_{id}[j]';\n",
      "      group_info_s2z_{id} += {level_weight}",
      "square(relative_precision_s2z);\n",
      "      h_group_s2z_{id} -= {level_weight}",
      "relative_precision_s2z .* white_level_s2z;\n",
      "      group_quad_s2z_{id} += {level_weight}",
      "dot_self(white_level_s2z);\n",
      "    }}\n",
      "  }}\n"
    )
  } else if (is_cor) {
    group_info <- if (is_student) {
      glue("sum(group_prec_s2z_{id})")
    } else {
      glue("1.0 * N_{id}")
    }
    group_score_code <- if (is_student) {
      glue(
        "    h_group_s2z_{id} = -white_group_s2z * ",
        "group_prec_s2z_{id};\n",
        "    group_quad_s2z_{id} = columns_dot_self(white_group_s2z) * ",
        "group_prec_s2z_{id};\n"
      )
    } else {
      glue(
        "    h_group_s2z_{id} = zeros_vector(M_{id});\n",
        "    group_quad_s2z_{id} = dot_self(to_vector(white_group_s2z));\n"
      )
    }
    str_add(out$tpar_comp) <- glue(
      "  {{\n",
      "    matrix[M_{id}, N_{id}] white_group_s2z = ",
      "mdivide_left_tri_low(L_Sigma_s2z_{id}, r_s2z_{id}');\n",
      "    group_info_s2z_{id} = {group_info};\n",
      "{group_score_code}",
      "  }}\n"
    )
  } else {
    group_info <- if (is_student) {
      glue("sum(group_prec_s2z_{id})")
    } else {
      glue("1.0 * N_{id}")
    }
    group_score_code <- if (is_student) {
      glue(
        "    h_group_s2z_{id} = zeros_vector(M_{id});\n",
        "    group_quad_s2z_{id} = 0.0;\n",
        "    for (k in 1:M_{id}) {{\n",
        "      vector[N_{id}] white_group_s2z = r_s2z_{id}[, k] / ",
        "{scale}[k];\n",
        "      h_group_s2z_{id}[k] = -dot_product(",
        "white_group_s2z, group_prec_s2z_{id});\n",
        "      group_quad_s2z_{id} += dot_product(",
        "group_prec_s2z_{id}, square(white_group_s2z));\n",
        "    }}\n"
      )
    } else {
      glue(
        "    h_group_s2z_{id} = zeros_vector(M_{id});\n",
        "    group_quad_s2z_{id} = 0.0;\n",
        "    for (k in 1:M_{id}) {{\n",
        "      vector[N_{id}] white_group_s2z = r_s2z_{id}[, k] / ",
        "{scale}[k];\n",
        "      group_quad_s2z_{id} += dot_self(white_group_s2z);\n",
        "    }}\n"
      )
    }
    str_add(out$tpar_comp) <- glue(
      "  {{\n",
      "    group_info_s2z_{id} = {group_info};\n",
      "{group_score_code}",
      "  }}\n"
    )
  }
  str_add(out$tpar_comp) <- cglue(
    "  {r_s2z} = r_s2z_{id}[, {J}];\n"
  )
  str_add(out$pll_args) <- cglue(", vector {r_s2z}")
  out
}

# Fast separated Gaussian path for multiple S2Z blocks. Rather than factoring
# the covariance of every omitted block mean, it factors only the induced
# covariance of theta and uses one Matheron update to recover all conventional
# means and population coefficients jointly.
.stan_re_s2z_joint_matheron <- function(set, prior, threads, normalize, ...) {
  out <- list(tpar_prior = "")
  infos <- set$infos
  matheron <- .stan_re_s2z_joint_matheron_info(set)
  stopifnot(!is.null(matheron))
  info <- infos[[1L]]
  q <- length(info$qnames)
  P <- matheron$P
  inactive <- matheron$inactive
  rdim <- length(P)
  P_index <- paste0("{", paste(P, collapse = ", "), "}")
  set_id <- set$set_id
  p <- info$p
  lpdf <- stan_lpdf_name(normalize)

  # A Student-t population prior is conditionally Gaussian, so the same
  # update applies after drawing its usual inverse-chi-square mixture scale.
  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    if (spec$dist == "student") {
      str_add(out$par) <- glue(
        "  real<lower=0> udf_b_s2z{p}_{k};",
        "  // mixing variable for population coefficient {k}\n"
      )
      str_add(out$tpar_prior) <- glue(
        "  lprior += inv_chi_square_{lpdf}(udf_b_s2z{p}_{k} | ",
        "{stan_s2z_arg_code(spec$df)});\n"
      )
    }
  }

  str_add(out$tpar_def) <- glue(
    "  // fast Gaussian Matheron system for S2Z blocks ",
    "{paste(set$ids, collapse = ', ')}\n",
    "  vector[{q}] prior_mean_s2z_{set_id};\n",
    "  vector<lower=0>[{q}] prior_scale_s2z_{set_id};\n",
    str_if(
      rdim == 1L,
      glue(
        "  real<lower=0> W_matheron_s2z_{set_id};\n",
        "  real<lower=0> sqrt_W_matheron_s2z_{set_id};\n",
        "  real theta_white_matheron_s2z_{set_id};\n"
      )
    ),
    str_if(
      rdim > 1L,
      glue(
        "  matrix[{rdim}, {rdim}] W_matheron_s2z_{set_id};\n",
        "  matrix[{rdim}, {rdim}] L_W_matheron_s2z_{set_id};\n",
        "  vector[{rdim}] theta_white_matheron_s2z_{set_id};\n"
      )
    ),
    "  real joint_quad_s2z_{set_id};\n"
  )
  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    str_add(out$tpar_comp) <- glue(
      "  prior_mean_s2z_{set_id}[{k}] = ",
      "{stan_s2z_arg_code(spec$location)};\n"
    )
    cond_scale <- if (spec$dist == "flat") {
      "1.0"
    } else if (spec$dist == "normal") {
      stan_s2z_arg_code(spec$scale)
    } else {
      glue(
        "{stan_s2z_arg_code(spec$scale)} * sqrt(",
        "{stan_s2z_arg_code(spec$df)} * udf_b_s2z{p}_{k})"
      )
    }
    str_add(out$tpar_comp) <- glue(
      "  prior_scale_s2z_{set_id}[{k}] = {cond_scale};\n"
    )
  }
  if (rdim == 1L) {
    str_add(out$tpar_comp) <- glue(
      "  W_matheron_s2z_{set_id} = ",
      "square(prior_scale_s2z_{set_id}[{P}]);\n"
    )
  }
  str_add(out$tpar_comp) <- glue(
    "  joint_quad_s2z_{set_id} = 0.0;\n"
  )
  for (b in seq_along(infos)) {
    id <- infos[[b]]$id
    if (rdim == 1L) {
      str_add(out$tpar_comp) <- glue(
        "  W_matheron_s2z_{set_id} += dot_self(",
        "H_s2z_{id}[{P}, ] * L_Sigma_s2z_{id}) / (1.0 * N_{id});\n"
      )
    } else if (rdim > 1L) {
      W_update <- if (b == 1L) {
        glue(
          "add_diag(\n",
          "      tcrossprod(H_active_s2z * L_Sigma_s2z_{id}) / ",
          "(1.0 * N_{id}),\n",
          "      square(prior_scale_s2z_{set_id}[{P_index}])\n",
          "    )"
        )
      } else {
        glue(
          "tcrossprod(",
          "H_active_s2z * L_Sigma_s2z_{id}) / (1.0 * N_{id})"
        )
      }
      W_operator <- if (b == 1L) "=" else "+="
      str_add(out$tpar_comp) <- glue(
        "  {{\n",
        "    matrix[{rdim}, M_{id}] H_active_s2z = ",
        "H_s2z_{id}[{P_index}, ];\n",
        "    W_matheron_s2z_{set_id} {W_operator} {W_update};\n",
        "  }}\n"
      )
    }
    str_add(out$tpar_comp) <- glue(
      "  joint_quad_s2z_{set_id} += group_quad_s2z_{id};\n"
    )
  }
  if (rdim == 1L) {
    str_add(out$tpar_comp) <- glue(
      "  sqrt_W_matheron_s2z_{set_id} = sqrt(W_matheron_s2z_{set_id});\n",
      "  theta_white_matheron_s2z_{set_id} = ",
      "(theta_s2z{p}[{P}] - prior_mean_s2z_{set_id}[{P}]) / ",
      "sqrt_W_matheron_s2z_{set_id};\n"
    )
  } else if (rdim > 1L) {
    str_add(out$tpar_comp) <- glue(
      "  L_W_matheron_s2z_{set_id} = ",
      "cholesky_decompose(W_matheron_s2z_{set_id});\n",
      "  {{\n",
      "    vector[{rdim}] theta_difference_s2z = ",
      "theta_s2z{p}[{P_index}] - prior_mean_s2z_{set_id}[{P_index}];\n",
      "    theta_white_matheron_s2z_{set_id} = mdivide_left_tri_low(",
      "L_W_matheron_s2z_{set_id}, theta_difference_s2z);\n",
      "  }}\n"
    )
  }

  for (k in inactive) {
    str_add(out$tpar_prior) <- glue(
      "  lprior += normal_{lpdf}(theta_s2z{p}[{k}] | ",
      "prior_mean_s2z_{set_id}[{k}], prior_scale_s2z_{set_id}[{k}]);\n"
    )
  }
  str_add(out$tpar_prior) <- glue(
    "  lprior += -0.5 * joint_quad_s2z_{set_id}\n"
  )
  if (rdim == 1L) {
    str_add(out$tpar_prior) <- glue(
      "    - 0.5 * square(theta_white_matheron_s2z_{set_id})\n",
      "    - log(sqrt_W_matheron_s2z_{set_id})\n"
    )
  } else if (rdim > 1L) {
    str_add(out$tpar_prior) <- glue(
      "    - 0.5 * dot_self(theta_white_matheron_s2z_{set_id})\n",
      "    - sum(log(diagonal(L_W_matheron_s2z_{set_id})))\n"
    )
  }
  for (b in seq_along(infos)) {
    r <- infos[[b]]$r
    id <- infos[[b]]$id
    mode <- re_s2z_center_mode(r)
    if (identical(mode, "centered")) {
      str_add(out$tpar_prior) <- glue(
        "    - (N_{id} - 1) * ",
        "sum(log(diagonal(L_Sigma_s2z_{id})))\n"
      )
    } else if (mode %in% c("partial", "auto")) {
      str_add(out$tpar_prior) <- glue(
        "    + log_det_partial_s2z_{id}\n"
      )
    }
    if (normalize) {
      str_add(out$tpar_prior) <- glue(
        "    - 0.5 * (N_{id} - 1) * M_{id} * log(2 * pi())\n"
      )
    }
  }
  if (normalize && rdim > 0L) {
    str_add(out$tpar_prior) <- glue(
      "    - 0.5 * {rdim} * log(2 * pi())\n"
    )
  }
  str_add(out$tpar_prior) <- "  ;\n"

  str_add(out$gen_def) <- glue(
    "  vector[{q}] q_recovered_s2z_{set_id};\n"
  )
  if (info$center) {
    str_add(out$gen_def) <- glue(
      "  real Intercept{p};\n",
      str_if(length(info$fixef), glue("  vector[Kc{p}] b{p};\n")),
      "  real b{p}_Intercept;\n"
    )
  } else {
    str_add(out$gen_def) <- glue("  vector[{q}] b{p};\n")
  }
  for (b in seq_along(infos)) {
    r <- infos[[b]]$r
    id <- infos[[b]]$id
    idp <- paste0(r$id, usc(combine_prefix(check_prefix(r))))
    r_public <- glue("r_{idp}_{r$cn}")
    is_cor <- nrow(r) > 1L && isTRUE(r$cor[1])
    str_add(out$gen_def) <- glue(
      "  vector[M_{id}] mean_r_s2z_{id};\n"
    )
    if (is_cor) {
      str_add(out$gen_def) <- glue(
        "  matrix[N_{id}, M_{id}] r_{id};\n"
      )
    }
    str_add(out$gen_def) <- cglue(
      "  vector[N_{id}] {r_public};\n"
    )
    if (is_cor) {
      str_add(out$gen_def) <- glue(
        "  // compute group-level correlations\n",
        "  corr_matrix[M_{id}] Cor_{id}",
        " = multiply_lower_tri_self_transpose(L_{id});\n",
        "  vector<lower=-1,upper=1>[NC_{id}] cor_{id};\n"
      )
    }
  }

  # Draw every m_f* and, when needed, the proper active part of beta*. One
  # shared innovation conditions them on theta. Flat active coordinates do not
  # enter W; beta is recovered as theta - sum(H_f m_f), which also enforces the
  # predictor identity to floating-point precision.
  str_add(out$gen_comp) <- glue(
    "  {{\n"
  )
  if (rdim > 0L) {
    str_add(out$gen_comp) <- glue(
      "    vector[{rdim}] theta_star_s2z;\n",
      "    vector[{rdim}] delta_s2z;\n"
    )
    for (a in seq_len(rdim)) {
      str_add(out$gen_comp) <- glue(
        "    theta_star_s2z[{a}] = prior_mean_s2z_{set_id}[{P[a]}] + ",
        "prior_scale_s2z_{set_id}[{P[a]}] * std_normal_rng();\n"
      )
    }
  }
  for (b in seq_along(infos)) {
    id <- infos[[b]]$id
    str_add(out$gen_comp) <- glue(
      "    {{\n",
      "      vector[M_{id}] z_mean_s2z;\n",
      "      for (k in 1:M_{id}) z_mean_s2z[k] = std_normal_rng();\n",
      "      mean_r_s2z_{id} = L_Sigma_s2z_{id} * z_mean_s2z / ",
      "sqrt(1.0 * N_{id});\n",
      "    }}\n"
    )
    if (rdim > 0L) {
      for (a in seq_len(rdim)) {
        str_add(out$gen_comp) <- glue(
          "    theta_star_s2z[{a}] += dot_product(",
          "H_s2z_{id}[{P[a]}, ], mean_r_s2z_{id}');\n"
        )
      }
    }
  }
  if (rdim == 1L) {
    str_add(out$gen_comp) <- glue(
      "    delta_s2z[1] = (theta_s2z{p}[{P}] - theta_star_s2z[1]) / ",
      "W_matheron_s2z_{set_id};\n"
    )
  } else if (rdim > 1L) {
    str_add(out$gen_comp) <- glue(
      "    {{\n",
      "      vector[{rdim}] theta_innovation_s2z;\n",
      "      vector[{rdim}] forward_solve_s2z;\n"
    )
    for (a in seq_len(rdim)) {
      str_add(out$gen_comp) <- glue(
        "      theta_innovation_s2z[{a}] = theta_s2z{p}[{P[a]}] - ",
        "theta_star_s2z[{a}];\n"
      )
    }
    str_add(out$gen_comp) <- glue(
      "      forward_solve_s2z = mdivide_left_tri_low(",
      "L_W_matheron_s2z_{set_id}, theta_innovation_s2z);\n",
      "      delta_s2z = (mdivide_right_tri_low(",
      "forward_solve_s2z', L_W_matheron_s2z_{set_id}))';\n",
      "    }}\n"
    )
  }
  str_add(out$gen_comp) <- glue(
    "    q_recovered_s2z_{set_id} = theta_s2z{p};\n"
  )
  for (b in seq_along(infos)) {
    r <- infos[[b]]$r
    id <- infos[[b]]$id
    J <- seq_rows(r)
    idp <- paste0(r$id, usc(combine_prefix(check_prefix(r))))
    r_s2z <- glue("r_s2z_{idp}_{r$cn}")
    r_public <- glue("r_{idp}_{r$cn}")
    is_cor <- nrow(r) > 1L && isTRUE(r$cor[1])
    if (rdim > 0L) {
      str_add(out$gen_comp) <- glue(
        "    {{\n",
      "      vector[M_{id}] active_score_s2z = ",
        "zeros_vector(M_{id});\n"
      )
      for (a in seq_len(rdim)) {
        str_add(out$gen_comp) <- glue(
          "      active_score_s2z += H_s2z_{id}[{P[a]}, ]' * ",
          "delta_s2z[{a}];\n"
        )
      }
      str_add(out$gen_comp) <- glue(
        "      mean_r_s2z_{id} += L_Sigma_s2z_{id} * ",
        "(L_Sigma_s2z_{id}' * active_score_s2z) / (1.0 * N_{id});\n",
        "    }}\n"
      )
    }
    str_add(out$gen_comp) <- glue(
      "    q_recovered_s2z_{set_id} -= H_s2z_{id} * mean_r_s2z_{id};\n"
    )
    if (is_cor) {
      str_add(out$gen_comp) <- glue(
        "    r_{id} = r_s2z_{id};\n",
        "    for (j in 1:N_{id}) r_{id}[j] += mean_r_s2z_{id}';\n",
        cglue("    {r_public} = r_{id}[, {J}];\n")
      )
    } else {
      str_add(out$gen_comp) <- cglue(
        "    {r_public} = {r_s2z} + mean_r_s2z_{id}[{J}];\n"
      )
    }
  }
  str_add(out$gen_comp) <- "  }\n"
  for (b in seq_along(infos)) {
    r <- infos[[b]]$r
    id <- infos[[b]]$id
    if (nrow(r) > 1L && isTRUE(r$cor[1])) {
      str_add(out$gen_comp) <- stan_cor_gen_comp(
        cor = glue("cor_{id}"), ncol = glue("M_{id}")
      )
    }
  }
  if (info$center) {
    str_add(out$gen_comp) <- glue(
      "  Intercept{p} = q_recovered_s2z_{set_id}[1];\n",
      str_if(
        length(info$fixef),
        glue(
          "  b{p} = tail(q_recovered_s2z_{set_id}, Kc{p});\n",
          "  b{p}_Intercept = Intercept{p} - dot_product(",
          "means_X{p}, b{p});\n"
        )
      ),
      str_if(
        !length(info$fixef),
        glue("  b{p}_Intercept = Intercept{p};\n")
      )
    )
  } else {
    str_add(out$gen_comp) <- glue(
      "  b{p} = q_recovered_s2z_{set_id};\n"
    )
  }
  out
}

# Jointly integrate the omitted group-effect means of all physical S2Z blocks
# in one linear predictor. Block means are represented in their own reference-
# whitened coordinates, so the independent group hierarchies contribute a
# block-diagonal precision while population-level priors provide the only
# cross-block coupling.
.stan_re_s2z_joint <- function(set, prior, threads, normalize, ...) {
  out <- list(tpar_prior = "")
  infos <- set$infos
  stopifnot(length(infos) > 1L || isTRUE(infos[[1L]]$ordinal))
  if (.stan_re_s2z_uses_explicit_mean(set)) {
    return(.stan_re_s2z_explicit_system(
      set, normalize = normalize, out = out
    ))
  }
  if (.stan_re_s2z_joint_uses_matheron(set)) {
    return(.stan_re_s2z_joint_matheron(
      set, prior = prior, threads = threads, normalize = normalize, ...
    ))
  }
  info <- infos[[1L]]
  q <- length(info$qnames)
  Ms <- vapply(infos, function(x) nrow(x$r), integer(1))
  total_M <- sum(Ms)
  starts <- cumsum(c(1L, head(Ms, -1L)))
  ends <- starts + Ms - 1L
  set_id <- set$set_id
  p <- info$p
  lpdf <- stan_lpdf_name(normalize)

  # Population Student-t and Cauchy priors use one mixing variable per fixed
  # coefficient for the entire joint system, rather than one copy per block.
  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    if (spec$dist == "student") {
      str_add(out$par) <- glue(
        "  real<lower=0> udf_b_s2z{p}_{k};",
        "  // mixing variable for population coefficient {k}\n"
      )
      str_add(out$tpar_prior) <- glue(
        "  lprior += inv_chi_square_{lpdf}(udf_b_s2z{p}_{k} | ",
        "{stan_s2z_arg_code(spec$df)});\n"
      )
    }
  }

  str_add(out$tpar_def) <- glue(
    "  // joint omitted-mean system for S2Z blocks ",
    "{paste(set$ids, collapse = ', ')}\n",
    "  vector[{q}] prior_mean_s2z_{set_id};\n",
    "  vector<lower=0>[{q}] prior_prec_s2z_{set_id};\n",
    "  matrix[{total_M}, {total_M}] P_s2z_{set_id};\n",
    "  matrix[{total_M}, {total_M}] L_P_s2z_{set_id};\n",
    "  vector[{total_M}] h_joint_s2z_{set_id};\n",
    "  vector[{total_M}] mhat_s2z_{set_id};\n",
    "  real joint_quad_s2z_{set_id};\n"
  )
  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    loc <- stan_s2z_arg_code(spec$location)
    str_add(out$tpar_comp) <- glue(
      "  prior_mean_s2z_{set_id}[{k}] = {loc};\n"
    )
    if (spec$dist == "flat") {
      prec <- "0.0"
    } else if (spec$dist == "normal") {
      prec <- glue("inv_square({stan_s2z_arg_code(spec$scale)})")
    } else {
      prec <- glue(
        "inv_square({stan_s2z_arg_code(spec$scale)} * sqrt(",
        "{stan_s2z_arg_code(spec$df)} * udf_b_s2z{p}_{k}))"
      )
    }
    str_add(out$tpar_comp) <- glue(
      "  prior_prec_s2z_{set_id}[{k}] = {prec};\n"
    )
  }

  str_add(out$tpar_comp) <- glue(
    "  joint_quad_s2z_{set_id} = 0.0;\n"
  )
  str_add(out$tpar_comp) <- glue(
    "  {{\n",
    "    matrix[{q}, {total_M}] prior_factor_s2z;\n",
    "    vector[{q}] prior_difference_s2z = ",
    "sqrt(prior_prec_s2z_{set_id}) .* (theta_s2z{p} - ",
    "prior_mean_s2z_{set_id});\n",
    "    vector[{total_M}] forward_solve_s2z;\n"
  )
  for (b in seq_along(infos)) {
    id <- infos[[b]]$id
    take <- glue("{starts[b]}:{ends[b]}")
    str_add(out$tpar_comp) <- glue(
      "    prior_factor_s2z[, {take}] = diag_pre_multiply(",
      "sqrt(prior_prec_s2z_{set_id}), H_s2z_{id} * ",
      "L_Sigma_s2z_{id});\n"
    )
  }
  str_add(out$tpar_comp) <- glue(
    "    P_s2z_{set_id} = crossprod(prior_factor_s2z);\n"
  )
  for (b in seq_along(infos)) {
    id <- infos[[b]]$id
    take <- glue("{starts[b]}:{ends[b]}")
    r <- infos[[b]]$r
    group_precision_kind <- stan_re_s2z_group_precision_kind(
      varying = identical(r$scale[1], "varying"),
      is_cor = nrow(r) > 1L && isTRUE(r$cor[1])
    )
    group_precision_code <- stan_re_s2z_add_group_precision_block(
      set_id, id, group_precision_kind, starts[b], ends[b]
    )
    str_add(out$tpar_comp) <- glue(
      "{group_precision_code}",
      "    h_joint_s2z_{set_id}[{take}] = h_group_s2z_{id};\n",
      "    joint_quad_s2z_{set_id} += group_quad_s2z_{id};\n"
    )
  }
  str_add(out$tpar_comp) <- glue(
    "    h_joint_s2z_{set_id} += prior_factor_s2z' * ",
    "prior_difference_s2z;\n",
    "    L_P_s2z_{set_id} = cholesky_decompose(P_s2z_{set_id});\n",
    "    forward_solve_s2z = mdivide_left_tri_low(",
    "L_P_s2z_{set_id}, h_joint_s2z_{set_id});\n",
    "    mhat_s2z_{set_id} = (mdivide_right_tri_low(",
    "forward_solve_s2z', L_P_s2z_{set_id}))';\n",
    "    joint_quad_s2z_{set_id} -= dot_self(forward_solve_s2z);\n",
    "  }}\n"
  )

  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    if (spec$dist == "flat") {
      next
    }
    if (spec$dist == "normal") {
      cond_scale <- stan_s2z_arg_code(spec$scale)
    } else {
      cond_scale <- glue(
        "{stan_s2z_arg_code(spec$scale)} * sqrt(",
        "{stan_s2z_arg_code(spec$df)} * udf_b_s2z{p}_{k})"
      )
    }
    str_add(out$tpar_prior) <- glue(
      "  lprior += normal_{lpdf}(theta_s2z{p}[{k}] | ",
      "{stan_s2z_arg_code(spec$location)}, {cond_scale});\n"
    )
  }
  str_add(out$tpar_prior) <- glue(
    "  lprior += -0.5 * joint_quad_s2z_{set_id}\n",
    "    - sum(log(diagonal(L_P_s2z_{set_id})))\n"
  )
  for (b in seq_along(infos)) {
    r <- infos[[b]]$r
    id <- infos[[b]]$id
    mode <- re_s2z_center_mode(r)
    if (identical(mode, "centered")) {
      str_add(out$tpar_prior) <- glue(
        "    - (N_{id} - 1) * ",
        "sum(log(diagonal(L_Sigma_s2z_{id})))\n"
      )
    } else if (mode %in% c("partial", "auto")) {
      str_add(out$tpar_prior) <- glue(
        "    + log_det_partial_s2z_{id}\n"
      )
    }
    if (identical(r$dist[1], "student")) {
      str_add(out$tpar_prior) <- glue(
        "    - M_{id} * sum(log(group_scale_s2z_{id}))\n"
      )
    }
    if (normalize && nzchar(r$cov[1])) {
      str_add(out$tpar_prior) <- glue(
        "    - M_{id} * sum(log(diagonal(Lcov_{id})))\n"
      )
    }
    if (normalize) {
      str_add(out$tpar_prior) <- glue(
        "    - 0.5 * N_{id} * M_{id} * log(2 * pi())\n",
        "    + 0.5 * M_{id} * log(1.0 * N_{id})\n"
      )
    }
  }
  if (normalize) {
    str_add(out$tpar_prior) <- glue(
      "    + 0.5 * {total_M} * log(2 * pi())\n"
    )
  }
  str_add(out$tpar_prior) <- "  ;\n"

  # One conditional draw preserves cross-block posterior covariance in the
  # conventional population/group parameterization.
  str_add(out$gen_def) <- glue(
    "  vector[{q}] q_recovered_s2z_{set_id};\n"
  )
  if (isTRUE(info$ordinal) && length(info$slope_q)) {
    str_add(out$gen_def) <- glue(
      "  vector[{length(info$slope_q)}] b{p};\n"
    )
  } else if (info$center && !isTRUE(info$ordinal)) {
    str_add(out$gen_def) <- glue(
      "  real Intercept{p};\n",
      str_if(length(info$fixef), glue("  vector[Kc{p}] b{p};\n")),
      "  real b{p}_Intercept;\n"
    )
  } else if (!isTRUE(info$ordinal)) {
    str_add(out$gen_def) <- glue("  vector[{q}] b{p};\n")
  }
  for (b in seq_along(infos)) {
    r <- infos[[b]]$r
    id <- infos[[b]]$id
    idp <- paste0(r$id, usc(combine_prefix(check_prefix(r))))
    r_public <- glue("r_{idp}_{r$cn}")
    is_cor <- nrow(r) > 1L && isTRUE(r$cor[1])
    str_add(out$gen_def) <- glue(
      "  vector[M_{id}] mean_r_s2z_{id};\n"
    )
    if (is_cor) {
      str_add(out$gen_def) <- glue(
        "  matrix[N_{id}, M_{id}] r_{id};\n"
      )
    }
    str_add(out$gen_def) <- cglue(
      "  vector[N_{id}] {r_public};\n"
    )
    if (is_cor) {
      str_add(out$gen_def) <- glue(
        "  // compute group-level correlations\n",
        "  corr_matrix[M_{id}] Cor_{id}",
        " = multiply_lower_tri_self_transpose(L_{id});\n",
        "  vector<lower=-1,upper=1>[NC_{id}] cor_{id};\n"
      )
    }
    if (identical(r$scale[1], "varying")) {
      str_add(out$gen_def) <- glue(
        "  matrix<lower=0>[N_{id}, M_{id}] sd_level_{id};\n"
      )
    }
  }

  str_add(out$gen_comp) <- glue(
    "  {{\n",
    "    vector[{total_M}] z_mean_s2z;\n",
    "    vector[{total_M}] mean_white_s2z;\n",
    "    for (k in 1:{total_M}) z_mean_s2z[k] = std_normal_rng();\n",
    "    mean_white_s2z = mhat_s2z_{set_id} + ",
    "(mdivide_right_tri_low(z_mean_s2z', L_P_s2z_{set_id}))';\n",
    "    q_recovered_s2z_{set_id} = theta_s2z{p};\n"
  )
  for (b in seq_along(infos)) {
    r <- infos[[b]]$r
    id <- infos[[b]]$id
    M <- nrow(r)
    J <- seq_rows(r)
    take <- glue("{starts[b]}:{ends[b]}")
    idp <- paste0(r$id, usc(combine_prefix(check_prefix(r))))
    r_s2z <- glue("r_s2z_{idp}_{r$cn}")
    r_public <- glue("r_{idp}_{r$cn}")
    is_cor <- M > 1L && isTRUE(r$cor[1])
    mean_transform <- if (is_cor) {
      glue("L_Sigma_s2z_{id} * mean_white_s2z[{take}]")
    } else {
      scale <- if (identical(r$scale[1], "varying")) {
        glue("reference_sd_s2z_{id}")
      } else {
        glue("sd_{id}")
      }
      glue("{scale} .* mean_white_s2z[{take}]")
    }
    str_add(out$gen_comp) <- glue(
      "    mean_r_s2z_{id} = {mean_transform};\n",
      "    q_recovered_s2z_{set_id} -= H_s2z_{id} * ",
      "mean_r_s2z_{id};\n"
    )
    if (is_cor) {
      str_add(out$gen_comp) <- glue(
        "    r_{id} = r_s2z_{id};\n",
        "    for (j in 1:N_{id}) r_{id}[j] += mean_r_s2z_{id}';\n",
        cglue("    {r_public} = r_{id}[, {J}];\n")
      )
    } else {
      str_add(out$gen_comp) <- cglue(
        "    {r_public} = {r_s2z} + mean_r_s2z_{id}[{J}];\n"
      )
    }
  }
  str_add(out$gen_comp) <- "  }\n"
  for (b in seq_along(infos)) {
    r <- infos[[b]]$r
    id <- infos[[b]]$id
    if (nrow(r) > 1L && isTRUE(r$cor[1])) {
      str_add(out$gen_comp) <- stan_cor_gen_comp(
        cor = glue("cor_{id}"), ncol = glue("M_{id}")
      )
    }
    if (identical(r$scale[1], "varying")) {
      str_add(out$gen_comp) <- glue(
        "  sd_level_{id} = sd_level_s2z_{id};\n"
      )
    }
  }
  if (isTRUE(info$ordinal) && length(info$slope_q)) {
    str_add(out$gen_comp) <- glue(
      "  b{p} = tail(q_recovered_s2z_{set_id}, ",
      "{length(info$slope_q)});\n"
    )
  } else if (info$center && !isTRUE(info$ordinal)) {
    str_add(out$gen_comp) <- glue(
      "  Intercept{p} = q_recovered_s2z_{set_id}[1];\n",
      str_if(
        length(info$fixef),
        glue(
          "  b{p} = tail(q_recovered_s2z_{set_id}, Kc{p});\n",
          "  b{p}_Intercept = Intercept{p} - dot_product(",
          "means_X{p}, b{p});\n"
        )
      ),
      str_if(
        !length(info$fixef),
        glue("  b{p}_Intercept = Intercept{p};\n")
      )
    )
  } else if (!isTRUE(info$ordinal)) {
    str_add(out$gen_comp) <- glue(
      "  b{p} = q_recovered_s2z_{set_id};\n"
    )
  }
  out
}

# Stan code for physical S2Z effects with coefficient- and group-specific
# scales. The standardized log-scale coordinates are only rotated into an
# orthonormal zero-sum part and a standardized mean; they are not marginalized.
# Thus all existing priors on the baseline sd parameters retain their exact
# meaning and the rotation has unit Jacobian. Only the common effect mean is
# integrated out.
.stan_re_s2z_varying_scale <- function(id, bframe, prior, threads, normalize,
                                       out = list(), fisher_info = NULL,
                                       stanvars = NULL) {
  lpdf <- stan_lpdf_name(normalize)
  # Avoid partial matching of $tpar_prior to $tpar_prior_const when a scale
  # parameter is fixed.
  if (is.null(out[["tpar_prior"]])) {
    out[["tpar_prior"]] <- ""
  }
  r <- subset2(bframe$frame$re, id = id)
  stopifnot(
    is.reframe(r), has_rows(r), all(r$s2z),
    all(r$scale == "varying")
  )
  bfl <- re_s2z_bframel(bframe, r)
  info <- re_s2z_info(
    bfl, prior = prior, stanvars = stanvars
  )
  stopifnot(!isTRUE(info$ordinal))
  q <- length(info$qnames)
  M <- nrow(r)
  J <- seq_rows(r)
  p <- info$p
  px <- check_prefix(r)
  idp <- paste0(r$id, usc(combine_prefix(px)))
  is_cor <- M > 1L && isTRUE(r$cor[1])
  is_student <- identical(r$dist[1], "student")
  has_cov <- nzchar(r$cov[1])
  group_precision_kind <- stan_re_s2z_group_precision_kind(
    varying = TRUE, is_cor = is_cor
  )
  s2z_mode <- re_s2z_center_mode(r)
  s2z_center <- identical(s2z_mode, "centered")
  s2z_fisher <- !is.null(fisher_info)
  s2z_partial <- s2z_mode %in% c("partial", "auto")
  stopifnot(!s2z_fisher || identical(s2z_mode, "auto"))
  r_s2z <- glue("r_s2z_{idp}_{r$cn}")
  r_public <- glue("r_{idp}_{r$cn}")

  if (s2z_fisher) {
    out <- stan_re_s2z_fisher_def(out, id)
    if (isTRUE(fisher_info$fixed_design)) {
      out <- stan_re_s2z_fisher_tdata(out, id, r, fisher_info)
    }
    if (has_cov) {
      out <- stan_re_s2z_cov_fisher_tdata(out, id)
    }
  } else if (s2z_partial) {
    str_add(out$data) <- glue(
      "  matrix<lower=0,upper=1>[N_{id}, M_{id}] rho_s2z_{id};",
      "  // fixed numeric centering fractions\n"
    )
    out <- stan_re_s2z_partial_tdata(out, id)
  }

  # sd remains brms's ordinary baseline group-level scale. sdlog controls
  # coefficient-specific variation of the log scale across observed groups.
  str_add_list(out) <- stan_prior(
    prior, class = "sdlog", group = r$group[1], coef = r$coef,
    type = glue("vector[M_{id}]"), suffix = glue("_{id}"), px = px,
    comment = "SDs of group-level log-scale deviations",
    normalize = normalize
  )

  if (is_cor) {
    str_add(out$data) <- glue(
      "  int<lower=1> NC_{id};  // number of group-level correlations\n"
    )
    str_add_list(out) <- stan_prior(
      prior, class = "L", group = r$group[1], suffix = usc(id),
      type = glue("cholesky_factor_corr[M_{id}]"),
      comment = "cholesky factor of correlation matrix",
      normalize = normalize
    )
  }

  str_add(out$fun) <- "  #include 'fun_sum_to_zero.stan'\n"
  str_add(out$par) <- glue(
    "  vector[M_{id} * (N_{id} - 1)] z_s2z_{id};",
    if (s2z_partial) {
      "  // partially centered orthonormal S2Z effect coordinates\n"
    } else if (s2z_center) {
      "  // physical orthonormal S2Z effect coordinates\n"
    } else {
      "  // standardized orthonormal S2Z effect coordinates\n"
    },
    "  vector[M_{id} * N_{id}] z_sd_s2z_{id};",
    "  // flattened orthonormal log-scale coordinates: contrasts, then means\n"
  )
  str_add(out$model_prior) <- glue(
    "  target += std_normal_{lpdf}(z_sd_s2z_{id});\n"
  )

  # Independent Student-t and Cauchy population priors are represented as
  # Gaussian scale mixtures, preserving the analytic effect-mean integral.
  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    if (spec$dist == "student") {
      str_add(out$par) <- glue(
        "  real<lower=0> udf_b_s2z{p}_{k};",
        "  // mixing variable for population coefficient {k}\n"
      )
      str_add(out$tpar_prior) <- glue(
        "  lprior += inv_chi_square_{lpdf}(udf_b_s2z{p}_{k} | ",
        "{stan_s2z_arg_code(spec$df)});\n"
      )
    }
  }

  # The scale transform is an orthogonal rotation of iid standard normals:
  # z_level = B z_centered + z_mean / sqrt(N). Relative scales have geometric
  # mean one, while sd_ref is the observed-group geometric mean scale.
  str_add(out$tpar_def) <- glue(
    "  matrix<lower=0>[N_{id}, M_{id}] sd_level_s2z_{id};\n",
    "  vector<lower=0>[M_{id}] reference_sd_s2z_{id};\n",
    str_if(
      is_student,
      glue(
        "  vector<lower=0>[N_{id}] group_scale_s2z_{id};\n",
        "  vector<lower=0>[N_{id}] group_prec_s2z_{id};\n"
      )
    )
  )
  str_add(out$tpar_comp) <- glue(
    "  reference_sd_s2z_{id} = sd_{id} .* exp(sdlog_{id} .* ",
    "tail(z_sd_s2z_{id}, M_{id}) / sqrt(1.0 * N_{id}));\n",
    "  for (k in 1:M_{id}) {{\n",
    "    vector[N_{id}] z_sd_centered_s2z = ",
    "sum_to_zero_constrain_brms(segment(z_sd_s2z_{id}, ",
    "(k - 1) * (N_{id} - 1) + 1, N_{id} - 1));\n",
    "    sd_level_s2z_{id}[, k] = reference_sd_s2z_{id}[k] * ",
    "exp(sdlog_{id}[k] * z_sd_centered_s2z);\n",
    "  }}\n"
  )
  if (is_student) {
    tr <- subset_reframe_dist(r, "student")
    g <- usc(tr$ggn[1])
    str_add(out$tpar_comp) <- glue(
      "  group_scale_s2z_{id} = dfm{g};\n",
      "  group_prec_s2z_{id} = inv_square(group_scale_s2z_{id});\n"
    )
  }
  str_add(out$tpar_prior) <- .stan_re_s2z_sd_level_prior(
    id, r = r, prior = prior, normalize = normalize
  )

  if (s2z_fisher) {
    L_fisher <- if (is_cor) {
      glue("diag_pre_multiply(reference_sd_s2z_{id}, L_{id})")
    } else {
      NULL
    }
    str_add(out$gen_comp) <- stan_re_s2z_fisher_gq_comp(
      id, r = r, fisher_info = fisher_info, L = L_fisher,
      scale = if (M == 1L) glue("reference_sd_s2z_{id}[1]") else NULL,
      row_var = if (has_cov) glue("row_var_fisher_s2z_{id}[j]") else NULL,
      diag_scale = if (!is_cor && M > 1L) {
        glue("reference_sd_s2z_{id}")
      }
    )
  }

  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    loc <- stan_s2z_arg_code(spec$location)
    str_add(out$tpar_comp) <- glue(
      "  prior_mean_s2z_{id}[{k}] = {loc};\n"
    )
    if (spec$dist == "flat") {
      prec <- "0.0"
    } else if (spec$dist == "normal") {
      prec <- glue("inv_square({stan_s2z_arg_code(spec$scale)})")
    } else {
      prec <- glue(
        "inv_square({stan_s2z_arg_code(spec$scale)} * sqrt(",
        "{stan_s2z_arg_code(spec$df)} * udf_b_s2z{p}_{k}))"
      )
    }
    str_add(out$tpar_comp) <- glue(
      "  prior_prec_s2z_{id}[{k}] = {prec};\n"
    )
  }

  if (is_cor || has_cov) {
    str_add(out$tdata_def) <- glue(
      "  matrix[{q}, M_{id}] H_s2z_{id};\n"
    )
    str_add(out$tdata_comp) <- glue(
      "  H_s2z_{id} = rep_matrix(0.0, {q}, M_{id});\n"
    )
    for (j in seq_len(M)) {
      qi <- info$match_q[j]
      str_add(out$tdata_comp) <- glue(
        "  H_s2z_{id}[{qi}, {j}] = 1.0;\n"
      )
      if (info$center && info$r$coef[j] != "Intercept") {
        str_add(out$tdata_comp) <- glue(
          "  H_s2z_{id}[1, {j}] = means_X{p}[{qi - 1L}];\n"
        )
      }
    }
    str_add(out$tpar_def) <- glue(
      "  // correlated physical S2Z effects with heterogeneous scales\n",
      "  matrix[N_{id}, M_{id}] r_s2z_{id};\n",
      "  vector[{q}] prior_mean_s2z_{id};\n",
      "  vector<lower=0>[{q}] prior_prec_s2z_{id};\n",
      "  matrix[M_{id}, M_{id}] L_Sigma_s2z_{id};\n",
      "  matrix[M_{id}, M_{id}] P_s2z_{id};\n",
      "  matrix[M_{id}, M_{id}] L_P_s2z_{id};\n",
      str_if(
        has_cov,
        glue(
          stan_re_s2z_group_precision_def(id, group_precision_kind),
          "  vector[M_{id}] h_group_s2z_{id};\n"
        )
      ),
      "  vector[M_{id}] mhat_s2z_{id};\n",
      "  real group_quad_s2z_{id};\n",
      str_if(
        s2z_partial,
        glue("  real log_det_partial_s2z_{id};\n")
      ),
      "  // using vectors speeds up likelihood indexing\n",
      cglue("  vector[N_{id}] {r_s2z};\n")
    )
    partial_transform <- if (s2z_partial) {
      stan_re_s2z_partial_cor_transform(id)
    } else {
      str_if(
        !s2z_center,
        glue("  r_s2z_{id} = r_s2z_{id} * L_Sigma_s2z_{id}';\n")
      )
    }
    L_ref_varying <- if (is_cor) {
      glue("diag_pre_multiply(reference_sd_s2z_{id}, L_{id})")
    } else {
      glue("diag_matrix(reference_sd_s2z_{id})")
    }
    str_add(out$tpar_comp) <- glue(
      "  for (k in 1:M_{id}) {{\n",
      "    r_s2z_{id}[, k] = sum_to_zero_constrain_brms(",
      "segment(z_s2z_{id}, (k - 1) * (N_{id} - 1) + 1, ",
      "N_{id} - 1));\n",
      "  }}\n",
      "  L_Sigma_s2z_{id} = {L_ref_varying};\n",
      "{partial_transform}"
    )
    if (has_cov) {
      str_add(out$tpar_comp) <- stan_re_s2z_cov_group_comp(
        id, varying = TRUE, is_cor = is_cor, is_student = is_student
      )
      str_add(out$tpar_comp) <- glue(
        "  {{\n",
        "    matrix[{q}, M_{id}] prior_factor_s2z = diag_pre_multiply(",
        "sqrt(prior_prec_s2z_{id}), H_s2z_{id}) * L_Sigma_s2z_{id};\n",
        "    vector[{q}] prior_difference_s2z = ",
        "sqrt(prior_prec_s2z_{id}) .* (theta_s2z{p} - ",
        "prior_mean_s2z_{id});\n",
        "    vector[M_{id}] h_s2z = prior_factor_s2z' * ",
        "prior_difference_s2z + h_group_s2z_{id};\n",
        "    vector[M_{id}] forward_solve_s2z;\n",
        "    P_s2z_{id} = ",
        stan_re_s2z_add_group_precision(
          "crossprod(prior_factor_s2z)", id, group_precision_kind
        ),
        ";\n",
        "    L_P_s2z_{id} = cholesky_decompose(P_s2z_{id});\n",
        "    forward_solve_s2z = mdivide_left_tri_low(",
        "L_P_s2z_{id}, h_s2z);\n",
        "    group_quad_s2z_{id} -= dot_self(forward_solve_s2z);\n",
        "    mhat_s2z_{id} = L_Sigma_s2z_{id} * ",
        "(mdivide_right_tri_low(forward_solve_s2z', L_P_s2z_{id}))';\n",
        "  }}\n"
      )
    } else {
      level_weight <- str_if(
        is_student, glue("group_prec_s2z_{id}[j] * ")
      )
      str_add(out$tpar_comp) <- glue(
        "  {{\n",
        "    matrix[{q}, M_{id}] prior_factor_s2z = diag_pre_multiply(",
        "sqrt(prior_prec_s2z_{id}), H_s2z_{id}) * L_Sigma_s2z_{id};\n",
        "    vector[{q}] prior_difference_s2z = ",
        "sqrt(prior_prec_s2z_{id}) .* (theta_s2z{p} - ",
        "prior_mean_s2z_{id});\n",
        "    vector[M_{id}] h_s2z;\n",
        "    vector[M_{id}] forward_solve_s2z;\n",
        "    P_s2z_{id} = crossprod(prior_factor_s2z);\n",
        "    h_s2z = prior_factor_s2z' * prior_difference_s2z;\n",
        "    group_quad_s2z_{id} = 0.0;\n",
        "    for (j in 1:N_{id}) {{\n",
        "      matrix[M_{id}, M_{id}] L_level_s2z = ",
        "diag_pre_multiply(sd_level_s2z_{id}[j]', L_{id});\n",
        "      matrix[M_{id}, M_{id}] relative_precision_s2z = ",
        "mdivide_left_tri_low(L_level_s2z, L_Sigma_s2z_{id});\n",
        "      vector[M_{id}] white_level_s2z = ",
        "mdivide_left_tri_low(L_level_s2z, r_s2z_{id}[j]');\n",
        "      P_s2z_{id} += {level_weight}",
        "crossprod(relative_precision_s2z);\n",
        "      h_s2z -= {level_weight}",
        "relative_precision_s2z' * white_level_s2z;\n",
        "      group_quad_s2z_{id} += {level_weight}",
        "dot_self(white_level_s2z);\n",
        "    }}\n",
        "    L_P_s2z_{id} = cholesky_decompose(P_s2z_{id});\n",
        "    forward_solve_s2z = mdivide_left_tri_low(",
        "L_P_s2z_{id}, h_s2z);\n",
        "    group_quad_s2z_{id} -= dot_self(forward_solve_s2z);\n",
        "    mhat_s2z_{id} = L_Sigma_s2z_{id} * ",
        "(mdivide_right_tri_low(forward_solve_s2z', L_P_s2z_{id}))';\n",
        "  }}\n"
      )
    }
    str_add(out$tpar_comp) <- cglue(
      "  {r_s2z} = r_s2z_{id}[, {J}];\n"
    )
    str_add(out$pll_args) <- cglue(", vector {r_s2z}")

    for (k in seq_len(q)) {
      spec <- info$prior[[k]]
      if (spec$dist == "flat") {
        next
      }
      if (spec$dist == "normal") {
        cond_scale <- stan_s2z_arg_code(spec$scale)
      } else {
        cond_scale <- glue(
          "{stan_s2z_arg_code(spec$scale)} * sqrt(",
          "{stan_s2z_arg_code(spec$df)} * udf_b_s2z{p}_{k})"
        )
      }
      str_add(out$tpar_prior) <- glue(
        "  lprior += normal_{lpdf}(theta_s2z{p}[{k}] | ",
        "{stan_s2z_arg_code(spec$location)}, {cond_scale});\n"
      )
    }
    str_add(out$tpar_prior) <- glue(
      "  lprior += -0.5 * group_quad_s2z_{id}\n",
      str_if(
        s2z_center,
        glue(
          "    - (N_{id} - 1) * ",
          "sum(log(diagonal(L_Sigma_s2z_{id})))\n"
        )
      ),
      str_if(
        s2z_partial,
        glue("    + log_det_partial_s2z_{id}\n")
      ),
      str_if(
        is_student,
        glue("    - M_{id} * sum(log(group_scale_s2z_{id}))\n")
      ),
      str_if(
        normalize && has_cov,
        glue("    - M_{id} * sum(log(diagonal(Lcov_{id})))\n")
      ),
      "    - sum(log(diagonal(L_P_s2z_{id})))",
      str_if(
        normalize,
        glue(" - 0.5 * N_{id} * M_{id} * log(2 * pi())",
             " + 0.5 * M_{id} * log(2 * pi())",
             " + 0.5 * M_{id} * log(1.0 * N_{id})")
      ),
      ";\n"
    )
  } else {
    # For diagonal scale matrices, the conditional precision is diagonal plus
    # the centered-intercept rank-one term. Scaling each equation by sd_ref^2
    # keeps the solve O(NM + M) and avoids inverse baseline scales.
    str_add(out$tpar_def) <- glue(
      "  // component-wise physical S2Z effects with heterogeneous scales\n",
      cglue("  vector[N_{id}] {r_s2z};\n"),
      "  vector[{q}] prior_mean_s2z_{id};\n",
      "  vector<lower=0>[{q}] prior_prec_s2z_{id};\n",
      "  vector<lower=0>[M_{id}] D_diag_s2z_{id};\n",
      "  real<lower=0> rank1_info_s2z_{id};\n",
      "  vector[M_{id}] mhat_s2z_{id};\n",
      "  vector[{q}] qhat_s2z_{id};\n",
      "  real<lower=0> group_quad_s2z_{id};\n",
      str_if(
        s2z_partial,
        glue("  real log_det_partial_s2z_{id};\n")
      )
    )
    str_add(out$tdata_def) <- glue(
      "  vector[M_{id}] intercept_map_s2z_{id};\n"
    )
    str_add(out$tdata_comp) <- glue(
      "  intercept_map_s2z_{id} = zeros_vector(M_{id});\n"
    )
    if (info$center) {
      for (j in seq_len(M)) {
        qi <- info$match_q[j]
        value <- if (info$r$coef[j] == "Intercept") {
          "1.0"
        } else {
          glue("means_X{p}[{qi - 1L}]")
        }
        str_add(out$tdata_comp) <- glue(
          "  intercept_map_s2z_{id}[{j}] = {value};\n"
        )
      }
    }
    if (s2z_partial) {
      str_add(out$tpar_comp) <- glue(
        "  log_det_partial_s2z_{id} = 0.0;\n"
      )
    }
    for (j in seq_len(M)) {
      if (s2z_partial) {
        str_add(out$tpar_comp) <- stan_re_s2z_partial_independent_transform(
          id, j, r_s2z[j], glue("reference_sd_s2z_{id}[{j}]")
        )
      } else {
        str_add(out$tpar_comp) <- glue(
          "  {r_s2z[j]} = sum_to_zero_constrain_brms(",
          str_if(
            !s2z_center,
            glue("reference_sd_s2z_{id}[{j}] * ")
          ),
          "segment(z_s2z_{id}, ({j} - 1) * (N_{id} - 1) + 1, ",
          "N_{id} - 1));\n"
        )
      }
    }
    str_add(out$tpar_comp) <- glue(
      "  {{\n",
      "    vector[M_{id}] base_info_s2z = zeros_vector(M_{id});\n",
      "    vector[M_{id}] base_score_s2z = zeros_vector(M_{id});\n",
      "    vector[M_{id}] group_info_s2z;\n",
      "    vector[M_{id}] group_score_s2z;\n",
      "    vector[M_{id}] scaled_score_s2z;\n",
      "    vector[M_{id}] independent_mode_s2z;\n"
    )
    for (j in seq_len(M)) {
      qi <- info$match_q[j]
      if (!info$center || qi != 1L) {
        str_add(out$tpar_comp) <- glue(
          "    base_info_s2z[{j}] = prior_prec_s2z_{id}[{qi}];\n",
          "    base_score_s2z[{j}] = prior_prec_s2z_{id}[{qi}] * ",
          "(theta_s2z{p}[{qi}] - prior_mean_s2z_{id}[{qi}]);\n"
        )
      }
      level_weight <- str_if(
        is_student, glue("group_prec_s2z_{id}[n] * ")
      )
      str_add(out$tpar_comp) <- glue(
        "    group_info_s2z[{j}] = 0.0;\n",
        "    group_score_s2z[{j}] = 0.0;\n",
        "    for (n in 1:N_{id}) {{\n",
        "      real relative_precision_s2z = reference_sd_s2z_{id}[{j}] / ",
        "sd_level_s2z_{id}[n, {j}];\n",
        "      real weighted_precision_s2z = {level_weight}",
        "square(relative_precision_s2z);\n",
        "      group_info_s2z[{j}] += weighted_precision_s2z;\n",
        "      group_score_s2z[{j}] += {r_s2z[j]}[n] * ",
        "weighted_precision_s2z;\n",
        "    }}\n"
      )
    }
    str_add(out$tpar_comp) <- glue(
      "    D_diag_s2z_{id} = group_info_s2z + ",
      "square(reference_sd_s2z_{id}) .* base_info_s2z;\n",
      "    scaled_score_s2z = square(reference_sd_s2z_{id}) .* ",
      "base_score_s2z - group_score_s2z;\n",
      "    independent_mode_s2z = scaled_score_s2z ./ D_diag_s2z_{id};\n"
    )
    coupling_prec <- str_if(
      info$center, glue("prior_prec_s2z_{id}[1]"), "0.0"
    )
    coupling_resid <- str_if(
      info$center,
      glue("theta_s2z{p}[1] - prior_mean_s2z_{id}[1]"),
      "0.0"
    )
    str_add(out$tpar_comp) <- glue(
      "    rank1_info_s2z_{id} = {coupling_prec} * dot_product(\n",
      "      square(reference_sd_s2z_{id}) .* ",
      "square(intercept_map_s2z_{id}),\n",
      "      1.0 ./ D_diag_s2z_{id}\n",
      "    );\n",
      "    mhat_s2z_{id} = independent_mode_s2z +\n",
      "      {coupling_prec} * square(reference_sd_s2z_{id}) .* ",
      "intercept_map_s2z_{id} ./\n",
      "      D_diag_s2z_{id} * ({coupling_resid} -\n",
      "      dot_product(intercept_map_s2z_{id}, ",
      "independent_mode_s2z)) /\n",
      "      (1.0 + rank1_info_s2z_{id});\n",
      "  }}\n",
      "  qhat_s2z_{id} = theta_s2z{p};\n"
    )
    if (info$center) {
      str_add(out$tpar_comp) <- glue(
        "  qhat_s2z_{id}[1] -= dot_product(",
        "intercept_map_s2z_{id}, mhat_s2z_{id});\n"
      )
      for (j in seq_len(M)) {
        qi <- info$match_q[j]
        if (qi != 1L) {
          str_add(out$tpar_comp) <- glue(
            "  qhat_s2z_{id}[{qi}] -= mhat_s2z_{id}[{j}];\n"
          )
        }
      }
    } else {
      for (j in seq_len(M)) {
        str_add(out$tpar_comp) <- glue(
          "  qhat_s2z_{id}[{info$match_q[j]}] -= ",
          "mhat_s2z_{id}[{j}];\n"
        )
      }
    }
    str_add(out$tpar_comp) <- glue(
      "  group_quad_s2z_{id} = 0.0;\n"
    )
    level_weight <- str_if(
      is_student, glue("group_prec_s2z_{id}[n] * ")
    )
    for (j in seq_len(M)) {
      str_add(out$tpar_comp) <- glue(
        "  for (n in 1:N_{id}) {{\n",
        "    real white_level_s2z = ({r_s2z[j]}[n] + ",
        "mhat_s2z_{id}[{j}]) / sd_level_s2z_{id}[n, {j}];\n",
        "    group_quad_s2z_{id} += {level_weight}",
        "square(white_level_s2z);\n",
        "  }}\n"
      )
    }
    str_add(out$pll_args) <- cglue(", vector {r_s2z}")

    for (k in seq_len(q)) {
      spec <- info$prior[[k]]
      if (spec$dist == "flat") {
        next
      }
      if (spec$dist == "normal") {
        cond_scale <- stan_s2z_arg_code(spec$scale)
      } else {
        cond_scale <- glue(
          "{stan_s2z_arg_code(spec$scale)} * sqrt(",
          "{stan_s2z_arg_code(spec$df)} * udf_b_s2z{p}_{k})"
        )
      }
      str_add(out$tpar_prior) <- glue(
        "  lprior += normal_{lpdf}(qhat_s2z_{id}[{k}] | ",
        "{stan_s2z_arg_code(spec$location)}, {cond_scale});\n"
      )
    }
    str_add(out$tpar_prior) <- glue(
      "  lprior += -0.5 * group_quad_s2z_{id}\n",
      str_if(
        s2z_center,
        glue(
          "    - (N_{id} - 1) * ",
          "sum(log(reference_sd_s2z_{id}))\n"
        )
      ),
      str_if(
        s2z_partial,
        glue("    + log_det_partial_s2z_{id}\n")
      ),
      str_if(
        is_student,
        glue("    - M_{id} * sum(log(group_scale_s2z_{id}))\n")
      ),
      "    - 0.5 * sum(log(D_diag_s2z_{id}))\n",
      "    - 0.5 * log1p(rank1_info_s2z_{id})",
      str_if(
        normalize,
        glue(" - 0.5 * N_{id} * M_{id} * log(2 * pi())",
             " + 0.5 * M_{id} * log(2 * pi())",
             " + 0.5 * M_{id} * log(1.0 * N_{id})")
      ),
      ";\n"
    )
  }

  # Reconstruct conventional public coefficients and deviations. The
  # observed-level scales are copied to generated quantities so they are
  # available under default save settings without exposing scale internals.
  str_add(out$gen_def) <- glue(
    "  vector[M_{id}] mean_r_s2z_{id};\n",
    "  vector[{q}] q_recovered_s2z_{id};\n",
    "  matrix<lower=0>[N_{id}, M_{id}] sd_level_{id};\n"
  )
  if (info$center) {
    str_add(out$gen_def) <- glue(
      "  real Intercept{p};\n",
      str_if(length(info$fixef), glue("  vector[Kc{p}] b{p};\n")),
      "  real b{p}_Intercept;\n"
    )
  } else {
    str_add(out$gen_def) <- glue("  vector[{q}] b{p};\n")
  }

  if (is_cor || has_cov) {
    str_add(out$gen_def) <- glue(
      "  matrix[N_{id}, M_{id}] r_{id};\n",
      cglue("  vector[N_{id}] {r_public};\n")
    )
    if (is_cor) {
      str_add(out$gen_def) <- glue(
        "  // compute group-level correlations\n",
        "  corr_matrix[M_{id}] Cor_{id}",
        " = multiply_lower_tri_self_transpose(L_{id});\n",
        "  vector<lower=-1,upper=1>[NC_{id}] cor_{id};\n"
      )
    }
    str_add(out$gen_comp) <- glue(
      "  {{\n",
      "    vector[M_{id}] z_mean_s2z;\n",
      "    for (k in 1:M_{id}) z_mean_s2z[k] = std_normal_rng();\n",
      "    mean_r_s2z_{id} = mhat_s2z_{id} + L_Sigma_s2z_{id} * ",
      "(mdivide_right_tri_low(z_mean_s2z', L_P_s2z_{id}))';\n",
      "  }}\n",
      "  q_recovered_s2z_{id} = theta_s2z{p} - ",
      "H_s2z_{id} * mean_r_s2z_{id};\n",
      "  r_{id} = r_s2z_{id};\n",
      "  for (j in 1:N_{id}) r_{id}[j] += mean_r_s2z_{id}';\n",
      cglue("  {r_public} = r_{id}[, {J}];\n")
    )
  } else {
    str_add(out$gen_def) <- cglue("  vector[N_{id}] {r_public};\n")
    str_add(out$gen_comp) <- glue(
      "  {{\n",
      "    vector[M_{id}] independent_noise_s2z;\n",
      "    real sqrt_rank1_s2z = sqrt(1.0 + rank1_info_s2z_{id});\n",
      "    real rank1_adjust_s2z;\n",
      "    for (k in 1:M_{id}) {{\n",
      "      independent_noise_s2z[k] = reference_sd_s2z_{id}[k] * ",
      "std_normal_rng() / sqrt(D_diag_s2z_{id}[k]);\n",
      "    }}\n",
      "    rank1_adjust_s2z = {coupling_prec} / ",
      "(sqrt_rank1_s2z * (1.0 + sqrt_rank1_s2z)) *\n",
      "      dot_product(intercept_map_s2z_{id}, ",
      "independent_noise_s2z);\n",
      "    mean_r_s2z_{id} = mhat_s2z_{id} + ",
      "independent_noise_s2z -\n",
      "      rank1_adjust_s2z * square(reference_sd_s2z_{id}) .* ",
      "intercept_map_s2z_{id} ./ D_diag_s2z_{id};\n",
      "  }}\n",
      "  q_recovered_s2z_{id} = theta_s2z{p};\n"
    )
    if (info$center) {
      str_add(out$gen_comp) <- glue(
        "  q_recovered_s2z_{id}[1] -= dot_product(",
        "intercept_map_s2z_{id}, mean_r_s2z_{id});\n"
      )
      for (j in seq_len(M)) {
        qi <- info$match_q[j]
        if (qi != 1L) {
          str_add(out$gen_comp) <- glue(
            "  q_recovered_s2z_{id}[{qi}] -= ",
            "mean_r_s2z_{id}[{j}];\n"
          )
        }
      }
    } else {
      for (j in seq_len(M)) {
        str_add(out$gen_comp) <- glue(
          "  q_recovered_s2z_{id}[{info$match_q[j]}] -= ",
          "mean_r_s2z_{id}[{j}];\n"
        )
      }
    }
    for (j in seq_len(M)) {
      str_add(out$gen_comp) <- glue(
        "  {r_public[j]} = {r_s2z[j]} + mean_r_s2z_{id}[{j}];\n"
      )
    }
  }

  str_add(out$gen_comp) <- glue(
    "  sd_level_{id} = sd_level_s2z_{id};\n"
  )
  if (info$center) {
    str_add(out$gen_comp) <- glue(
      "  Intercept{p} = q_recovered_s2z_{id}[1];\n",
      str_if(
        length(info$fixef),
        glue("  b{p} = tail(q_recovered_s2z_{id}, Kc{p});\n",
             "  b{p}_Intercept = Intercept{p} - dot_product(",
             "means_X{p}, b{p});\n")
      ),
      str_if(
        !length(info$fixef),
        glue("  b{p}_Intercept = Intercept{p};\n")
      )
    )
  } else {
    str_add(out$gen_comp) <- glue("  b{p} = q_recovered_s2z_{id};\n")
  }
  if (is_cor) {
    str_add(out$gen_comp) <- stan_cor_gen_comp(
      cor = glue("cor_{id}"), ncol = glue("M_{id}")
    )
  }
  out
}

# Stan code for a strict Gaussian latent-score block. Unlike conventional S2Z
# group effects, these scores are defined directly on the componentwise
# zero-sum subspace: there is no omitted group mean, no theta_s2z replacement,
# and no generated-quantities mean recovery. A single global ID can therefore
# provide correlated score columns to several nonlinear predictors.
.stan_re_s2z_latent <- function(id, r, bframe, prior, normalize,
                                out = list(), fisher_info = NULL) {
  stopifnot(
    is.anybrmsframe(bframe), is.reframe(r), has_rows(r), all(r$s2z),
    all(re_s2z_latent(r)), all(r$dist == "gaussian"),
    all(r$scale == "shared")
  )
  if (is.null(out[["tpar_prior"]])) {
    out[["tpar_prior"]] <- ""
  }
  lpdf <- stan_lpdf_name(normalize)
  r_dim <- re_s2z_latent_dimensions(r)
  M <- nrow(r_dim)
  J <- re_s2z_latent_dimension(r)
  idp <- paste0(r$id, usc(combine_prefix(check_prefix(r))))
  r_s2z <- glue("r_s2z_{idp}_{r$cn}")
  r_public <- glue("r_{idp}_{r$cn}")
  is_cor <- M > 1L && isTRUE(r$cor[1])
  mode <- re_s2z_center_mode(r)
  s2z_center <- identical(mode, "centered")
  s2z_fisher <- !is.null(fisher_info)
  s2z_partial <- mode %in% c("partial", "auto")
  stopifnot(
    mode %in% c("centered", "noncentered", "partial", "auto"),
    !s2z_fisher || identical(mode, "auto")
  )

  if (s2z_fisher) {
    out <- stan_re_s2z_fisher_def(out, id)
    # Nonlinear dependency values are needed only by the precursor proposal.
    # Keeping both their declarations and assignments in generated quantities
    # prevents them from entering Pathfinder or HMC gradients.
    str_add(out$gen_def) <- fisher_info$dependency_tpar_def %||% ""
  } else if (s2z_partial) {
    str_add(out$data) <- glue(
      "  matrix<lower=0,upper=1>[N_{id}, M_{id}] rho_s2z_{id};",
      "  // fixed numeric centering fractions\n"
    )
    out <- stan_re_s2z_partial_tdata(out, id)
  }
  if (is_cor) {
    str_add(out$data) <- glue(
      "  int<lower=1> NC_{id};  // number of group-level correlations\n"
    )
    str_add_list(out) <- stan_prior(
      prior, class = "L", group = r$group[1], suffix = usc(id),
      type = glue("cholesky_factor_corr[M_{id}]"),
      comment = "cholesky factor of latent-score correlation matrix",
      normalize = normalize
    )
  }

  str_add(out$fun) <- "  #include 'fun_sum_to_zero.stan'\n"
  str_add(out$par) <- glue(
    "  vector[M_{id} * (N_{id} - 1)] z_s2z_{id};",
    if (s2z_partial) {
      "  // partially centered strict latent-score coordinates\n"
    } else if (s2z_center) {
      "  // physical strict latent-score coordinates\n"
    } else {
      "  // standardized strict latent-score coordinates\n"
    }
  )
  str_add(out$tpar_def) <- glue(
    "  // strict componentwise sum-to-zero latent-score block {id}\n",
    "  matrix[N_{id}, M_{id}] r_s2z_{id};\n",
    "  matrix[M_{id}, M_{id}] L_Sigma_s2z_{id};\n",
    str_if(
      s2z_center || s2z_partial,
      glue("  real<lower=0> group_quad_s2z_{id};\n")
    ),
    str_if(s2z_partial, glue("  real log_det_partial_s2z_{id};\n")),
    "  // vectors used by the observation-level nonlinear predictors\n",
    cglue("  vector[N_{id}] {r_s2z};\n")
  )

  if (is_cor) {
    str_add(out$tpar_comp) <- glue(
      "  L_Sigma_s2z_{id} = diag_pre_multiply(sd_{id}, L_{id});\n"
    )
  } else {
    str_add(out$tpar_comp) <- glue(
      "  L_Sigma_s2z_{id} = diag_matrix(sd_{id});\n"
    )
  }
  str_add(out$tpar_comp) <- glue(
    "  for (k in 1:M_{id}) {{\n",
    "    r_s2z_{id}[, k] = sum_to_zero_constrain_brms(",
    "segment(z_s2z_{id}, (k - 1) * (N_{id} - 1) + 1, ",
    "N_{id} - 1));\n",
    "  }}\n"
  )
  if (s2z_partial) {
    if (s2z_fisher) {
      # Candidate evaluation must follow its generated-quantity dependency
      # assignments. The sampled chart below always uses fixed rho data.
      str_add(out$gen_comp) <- stan_re_s2z_fisher_gq_comp(
        id, r = r, fisher_info = fisher_info,
        L = if (is_cor) glue("L_Sigma_s2z_{id}") else NULL,
        scale = if (M == 1L) glue("sd_{id}[1]") else NULL,
        diag_scale = if (!is_cor && M > 1L) glue("sd_{id}") else NULL,
        precompute = fisher_info$dependency_tpar_comp %||% ""
      )
    }
    str_add(out$tpar_comp) <- stan_re_s2z_partial_cor_transform(id)
  } else if (!s2z_center) {
    str_add(out$tpar_comp) <- glue(
      "  r_s2z_{id} = r_s2z_{id} * L_Sigma_s2z_{id}';\n"
    )
  }
  if (s2z_center || s2z_partial) {
    str_add(out$tpar_comp) <- glue(
      "  {{\n",
      "    matrix[M_{id}, N_{id}] white_latent_s2z = ",
      "mdivide_left_tri_low(L_Sigma_s2z_{id}, r_s2z_{id}');\n",
      "    group_quad_s2z_{id} = dot_self(to_vector(white_latent_s2z));\n",
      "  }}\n"
    )
  }
  str_add(out$tpar_comp) <- cglue(
    "  {r_s2z} = r_s2z_{id}[, {J}];\n"
  )
  str_add(out$pll_args) <- cglue(", vector {r_s2z}")

  if (identical(mode, "noncentered")) {
    str_add(out$tpar_prior) <- glue(
      "  lprior += std_normal_{lpdf}(z_s2z_{id});\n"
    )
  } else {
    str_add(out$tpar_prior) <- glue(
      "  lprior += -0.5 * group_quad_s2z_{id}\n",
      str_if(
        s2z_center,
        glue(
          "    - (N_{id} - 1) * ",
          "sum(log(diagonal(L_Sigma_s2z_{id})))\n"
        )
      ),
      str_if(s2z_partial, glue("    + log_det_partial_s2z_{id}\n")),
      str_if(
        normalize,
        glue("    - 0.5 * (N_{id} - 1) * M_{id} * log(2 * pi())\n")
      ),
      "  ;\n"
    )
  }

  if (is_cor) {
    str_add(out$gen_def) <- glue(
      "  matrix[N_{id}, M_{id}] r_{id};\n",
      cglue("  vector[N_{id}] {r_public};\n"),
      "  // compute latent-score correlations\n",
      "  corr_matrix[M_{id}] Cor_{id}",
      " = multiply_lower_tri_self_transpose(L_{id});\n",
      "  vector<lower=-1,upper=1>[NC_{id}] cor_{id};\n"
    )
    str_add(out$gen_comp) <- glue(
      "  r_{id} = r_s2z_{id};\n",
      cglue("  {r_public} = r_{id}[, {J}];\n")
    )
    str_add(out$gen_comp) <- stan_cor_gen_comp(
      cor = glue("cor_{id}"), ncol = glue("M_{id}")
    )
  } else {
    str_add(out$gen_def) <- cglue("  vector[N_{id}] {r_public};\n")
    str_add(out$gen_comp) <- cglue("  {r_public} = {r_s2z};\n")
  }
  out
}

# Stan code for one physical sum-to-zero group-level block. The omitted common
# group-effect mean vector is integrated out analytically. Conditional
# Gaussian scale mixtures are handled by a group-specific scale, which makes
# this exact for both Gaussian and Student-t group effects.
.stan_re_s2z <- function(id, bframe, prior, threads, normalize, out = list(),
                         fisher_info = NULL, stanvars = NULL) {
  lpdf <- stan_lpdf_name(normalize)
  # Avoid partial matching of $tpar_prior to $tpar_prior_const when a group
  # scale is fixed. Otherwise the constant assignment is appended to the
  # model block when normalize = FALSE.
  if (is.null(out[["tpar_prior"]])) {
    out[["tpar_prior"]] <- ""
  }
  r <- subset2(bframe$frame$re, id = id)
  stopifnot(is.reframe(r), has_rows(r), all(r$s2z))
  bfl <- re_s2z_bframel(bframe, r)
  info <- re_s2z_info(
    bfl, prior = prior, stanvars = stanvars
  )
  stopifnot(!isTRUE(info$ordinal))
  q <- length(info$qnames)
  M <- nrow(r)
  J <- seq_rows(r)
  p <- info$p
  idp <- paste0(r$id, usc(combine_prefix(check_prefix(r))))
  is_cor <- M > 1L && isTRUE(r$cor[1])
  is_student <- identical(r$dist[1], "student")
  has_cov <- nzchar(r$cov[1])
  s2z_mode <- re_s2z_center_mode(r)
  s2z_center <- identical(s2z_mode, "centered")
  s2z_fisher <- !is.null(fisher_info)
  s2z_partial <- s2z_mode %in% c("partial", "auto")
  stopifnot(!s2z_fisher || identical(s2z_mode, "auto"))

  single_set <- nlist(set_id = id, ids = id, bfl, infos = list(info))
  if (.stan_re_s2z_uses_explicit_mean(single_set)) {
    return(.stan_re_s2z_explicit(
      id, bframe = bframe, prior = prior,
      normalize = normalize, out = out, fisher_info = fisher_info,
      stanvars = stanvars
    ))
  }

  if (s2z_fisher) {
    out <- stan_re_s2z_fisher_def(out, id)
    if (isTRUE(fisher_info$fixed_design)) {
      out <- stan_re_s2z_fisher_tdata(out, id, r, fisher_info)
    }
    if (has_cov) {
      out <- stan_re_s2z_cov_fisher_tdata(out, id)
    }
  }

  if (M == 1L && !has_cov) {
    return(.stan_re_s2z_scalar(
      id, r = r, info = info, normalize = normalize, out = out,
      fisher_info = fisher_info
    ))
  }

  if (!is_cor && !has_cov) {
    return(.stan_re_s2z_independent(
      id, r = r, info = info, normalize = normalize, out = out,
      fisher_info = fisher_info
    ))
  }

  if (s2z_partial && !s2z_fisher) {
    str_add(out$data) <- glue(
      "  matrix<lower=0,upper=1>[N_{id}, M_{id}] rho_s2z_{id};",
      "  // fixed numeric centering fractions\n"
    )
    out <- stan_re_s2z_partial_tdata(out, id)
  }

  # Keep brms's existing covariance parameters and priors. Only the
  # group-effect coordinates and their joint prior are replaced.
  if (is_cor) {
    str_add(out$data) <- glue(
      "  int<lower=1> NC_{id};  // number of group-level correlations\n"
    )
    str_add_list(out) <- stan_prior(
      prior, class = "L", group = r$group[1], suffix = usc(id),
      type = glue("cholesky_factor_corr[M_{id}]"),
      comment = "cholesky factor of correlation matrix",
      normalize = normalize
    )
  }

  str_add(out$fun) <- "  #include 'fun_sum_to_zero.stan'\n"
  str_add(out$par) <- glue(
    "  vector[M_{id} * (N_{id} - 1)] z_s2z_{id};",
    if (s2z_partial) {
      "  // partially centered orthonormal S2Z coordinates\n"
    } else if (s2z_center) {
      "  // physical orthonormal S2Z coordinates\n"
    } else {
      "  // standardized orthonormal S2Z coordinates\n"
    }
  )

  # Independent Student-t and Cauchy population priors are represented as
  # Gaussian scale mixtures, so the omitted group mean remains analytic.
  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    if (spec$dist == "student") {
      str_add(out$par) <- glue(
        "  real<lower=0> udf_b_s2z{p}_{k};",
        "  // mixing variable for population coefficient {k}\n"
      )
      str_add(out$tpar_prior) <- glue(
        "  lprior += inv_chi_square_{lpdf}(udf_b_s2z{p}_{k} | ",
        "{stan_s2z_arg_code(spec$df)});\n"
      )
    }
  }

  str_add(out$tpar_def) <- glue(
    "  // physical sum-to-zero group-level effects of ID {id}\n",
    "  matrix[N_{id}, M_{id}] r_s2z_{id};\n",
    "  vector[{q}] prior_mean_s2z_{id};\n",
    "  vector<lower=0>[{q}] prior_prec_s2z_{id};\n",
    str_if(
      is_student,
      glue(
        "  vector<lower=0>[N_{id}] group_scale_s2z_{id};\n",
        "  vector<lower=0>[N_{id}] group_prec_s2z_{id};\n"
      )
    ),
    "  matrix[M_{id}, M_{id}] L_Sigma_s2z_{id};\n",
    "  // precision of the omitted mean in L_Sigma-whitened coordinates\n",
    "  matrix[M_{id}, M_{id}] P_s2z_{id};\n",
    "  matrix[M_{id}, M_{id}] L_P_s2z_{id};\n",
    str_if(
      has_cov,
      glue(
        stan_re_s2z_group_precision_def(id, "scalar"),
        "  vector[M_{id}] h_group_s2z_{id};\n"
      )
    ),
    str_if(
      is_student,
      glue("  vector[M_{id}] contrast_score_s2z_{id};\n")
    ),
    "  real group_quad_s2z_{id};\n",
    str_if(
      s2z_partial,
      glue("  real log_det_partial_s2z_{id};\n")
    ),
    "  // using vectors speeds up indexing in loops\n",
    cglue("  vector[N_{id}] r_s2z_{idp}_{r$cn};\n")
  )

  # The mapping H absorbs the omitted raw group mean into brms's population
  # coordinates. For centered predictors its first row contains the actual
  # means of every matching raw design column, including interactions. It is
  # fixed by the data and belongs outside the autodiff graph.
  str_add(out$tdata_def) <- glue(
    "  matrix[{q}, M_{id}] H_s2z_{id};\n"
  )
  str_add(out$tdata_comp) <- glue(
    "  H_s2z_{id} = rep_matrix(0.0, {q}, M_{id});\n"
  )
  for (j in seq_len(M)) {
    qi <- info$match_q[j]
    str_add(out$tdata_comp) <- glue(
      "  H_s2z_{id}[{qi}, {j}] = 1.0;\n"
    )
    if (info$center && info$r$coef[j] != "Intercept") {
      str_add(out$tdata_comp) <- glue(
        "  H_s2z_{id}[1, {j}] = means_X{p}[{qi - 1L}];\n"
      )
    }
  }
  if (is_cor) {
    str_add(out$tpar_comp) <- glue(
      "  L_Sigma_s2z_{id} = diag_pre_multiply(sd_{id}, L_{id});\n"
    )
  } else {
    str_add(out$tpar_comp) <- glue(
      "  L_Sigma_s2z_{id} = diag_matrix(sd_{id});\n"
    )
  }
  if (s2z_fisher) {
    str_add(out$gen_comp) <- stan_re_s2z_fisher_gq_comp(
      id, r = r, fisher_info = fisher_info,
      L = if (is_cor) glue("L_Sigma_s2z_{id}") else NULL,
      scale = if (M == 1L) glue("sd_{id}[1]") else NULL,
      row_var = if (has_cov) glue("row_var_fisher_s2z_{id}[j]") else NULL,
      diag_scale = if (!is_cor && M > 1L) glue("sd_{id}") else NULL
    )
  }
  partial_transform <- if (s2z_partial) {
    stan_re_s2z_partial_cor_transform(id)
  } else {
    str_if(
      !s2z_center,
      glue("  r_s2z_{id} = r_s2z_{id} * L_Sigma_s2z_{id}';\n")
    )
  }
  str_add(out$tpar_comp) <- glue(
    "  for (k in 1:M_{id}) {{\n",
    "    r_s2z_{id}[, k] = sum_to_zero_constrain_brms(",
    "segment(z_s2z_{id}, (k - 1) * (N_{id} - 1) + 1, ",
    "N_{id} - 1));\n",
    "  }}\n",
    "{partial_transform}"
  )

  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    loc <- stan_s2z_arg_code(spec$location)
    str_add(out$tpar_comp) <- glue(
      "  prior_mean_s2z_{id}[{k}] = {loc};\n"
    )
    if (spec$dist == "flat") {
      prec <- "0.0"
    } else if (spec$dist == "normal") {
      prec <- glue("inv_square({stan_s2z_arg_code(spec$scale)})")
    } else {
      prec <- glue(
        "inv_square({stan_s2z_arg_code(spec$scale)} * sqrt(",
        "{stan_s2z_arg_code(spec$df)} * udf_b_s2z{p}_{k}))"
      )
    }
    str_add(out$tpar_comp) <- glue(
      "  prior_prec_s2z_{id}[{k}] = {prec};\n"
    )
  }

  if (is_student) {
    tr <- subset_reframe_dist(r, "student")
    g <- usc(tr$ggn[1])
    str_add(out$tpar_comp) <- glue(
      "  group_scale_s2z_{id} = dfm{g};\n",
      "  group_prec_s2z_{id} = inv_square(group_scale_s2z_{id});\n"
    )
  }
  if (has_cov) {
    str_add(out$tpar_comp) <- stan_re_s2z_cov_group_comp(
      id, varying = FALSE, is_cor = is_cor, is_student = is_student
    )
    contrast_score <- glue(" + h_group_s2z_{id}")
    str_add(out$tpar_comp) <- glue(
      "  {{\n",
      "    matrix[{q}, M_{id}] prior_factor_s2z = diag_pre_multiply(",
      "sqrt(prior_prec_s2z_{id}), H_s2z_{id}) * L_Sigma_s2z_{id};\n",
      "    vector[{q}] prior_difference_s2z = sqrt(prior_prec_s2z_{id}) .* ",
      "(theta_s2z{p} - prior_mean_s2z_{id});\n",
      "    vector[M_{id}] h_s2z;\n",
      "    vector[M_{id}] whitened_h_s2z;\n",
      "    P_s2z_{id} = add_diag(crossprod(prior_factor_s2z), ",
      "group_info_s2z_{id});\n",
      "    h_s2z = prior_factor_s2z' * prior_difference_s2z + ",
      "h_group_s2z_{id};\n",
      "    L_P_s2z_{id} = cholesky_decompose(P_s2z_{id});\n",
      "    whitened_h_s2z = mdivide_left_tri_low(L_P_s2z_{id}, h_s2z);\n",
      "    group_quad_s2z_{id} -= dot_self(whitened_h_s2z);\n",
      "  }}\n",
      cglue("  r_s2z_{idp}_{r$cn} = r_s2z_{id}[, {J}];\n")
    )
  } else {
    group_info <- str_if(
      is_student, glue("sum(group_prec_s2z_{id})"), glue("1.0 * N_{id}")
    )
    contrast_score <- str_if(
      is_student, glue(" - contrast_score_s2z_{id}")
    )
    contrast_score_code <- str_if(
      is_student,
      glue(
        "    contrast_score_s2z_{id} = white_s2z * ",
        "group_prec_s2z_{id};\n"
      )
    )
    group_quad_code <- if (is_student) {
      glue(
        "    group_quad_s2z_{id} += columns_dot_self(white_s2z) * ",
        "group_prec_s2z_{id};\n"
      )
    } else {
      glue(
        "    group_quad_s2z_{id} += dot_self(to_vector(white_s2z));\n"
      )
    }
    str_add(out$tpar_comp) <- glue(
      "  {{\n",
      "    matrix[{q}, M_{id}] prior_factor_s2z = diag_pre_multiply(",
      "sqrt(prior_prec_s2z_{id}), H_s2z_{id}) * L_Sigma_s2z_{id};\n",
      "    vector[{q}] prior_difference_s2z = ",
      "sqrt(prior_prec_s2z_{id}) .* ",
      "(theta_s2z{p} - prior_mean_s2z_{id});\n",
      "    matrix[M_{id}, N_{id}] white_s2z = ",
      "mdivide_left_tri_low(L_Sigma_s2z_{id}, r_s2z_{id}');\n",
      "    vector[M_{id}] h_s2z;\n",
      "    vector[M_{id}] whitened_h_s2z;\n",
      "{contrast_score_code}",
      "    P_s2z_{id} = add_diag(crossprod(prior_factor_s2z), ",
      "{group_info});\n",
      "    h_s2z = prior_factor_s2z' * prior_difference_s2z",
      "{contrast_score};\n",
      "    L_P_s2z_{id} = cholesky_decompose(P_s2z_{id});\n",
      "    whitened_h_s2z = mdivide_left_tri_low(L_P_s2z_{id}, h_s2z);\n",
      "    group_quad_s2z_{id} = -dot_self(whitened_h_s2z);\n",
      "{group_quad_code}",
      "  }}\n",
      cglue("  r_s2z_{idp}_{r$cn} = r_s2z_{id}[, {J}];\n")
    )
  }
  str_add(out$pll_args) <- cglue(
    ", vector r_s2z_{idp}_{r$cn}"
  )

  # Score the original conditional Gaussian hierarchy at an omitted mean of
  # zero. group_quad stores the base group quadratic minus the completed-square
  # correction in L_Sigma-whitened coordinates. This needs one triangular solve
  # and never forms the group covariance or precision. The final log(N) is the
  # measure correction for the orthonormal contrast coordinates.
  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    if (spec$dist == "flat") {
      next
    }
    if (spec$dist == "normal") {
      cond_scale <- stan_s2z_arg_code(spec$scale)
    } else {
      cond_scale <- glue(
        "{stan_s2z_arg_code(spec$scale)} * sqrt(",
        "{stan_s2z_arg_code(spec$df)} * udf_b_s2z{p}_{k})"
      )
    }
    str_add(out$tpar_prior) <- glue(
      "  lprior += normal_{lpdf}(theta_s2z{p}[{k}] | ",
      "{stan_s2z_arg_code(spec$location)}, {cond_scale});\n"
    )
  }
  str_add(out$tpar_prior) <- glue(
    "  lprior += -0.5 * group_quad_s2z_{id}\n",
    str_if(
      s2z_center,
      glue(
        "    - (N_{id} - 1) * ",
        "sum(log(diagonal(L_Sigma_s2z_{id})))\n"
      )
    ),
    str_if(
      s2z_partial,
      glue("    + log_det_partial_s2z_{id}\n")
    ),
    str_if(
      is_student,
      glue("    - M_{id} * sum(log(group_scale_s2z_{id}))\n")
    ),
    str_if(
      normalize && has_cov,
      glue("    - M_{id} * sum(log(diagonal(Lcov_{id})))\n")
    ),
    "    - sum(log(diagonal(L_P_s2z_{id})))",
    str_if(
      normalize,
      glue(" - 0.5 * N_{id} * M_{id} * log(2 * pi())",
           " + 0.5 * M_{id} * log(2 * pi())",
           " + 0.5 * M_{id} * log(1.0 * N_{id})")
    ),
    ";\n"
  )

  # Recover conventional super-population coefficients and deviations using
  # one shared conditional draw of the analytically omitted group mean.
  str_add(out$gen_def) <- glue(
    "  vector[M_{id}] mean_r_s2z_{id};\n",
    "  vector[{q}] q_recovered_s2z_{id};\n",
    "  matrix[N_{id}, M_{id}] r_{id};\n"
  )
  if (info$center) {
    str_add(out$gen_def) <- glue(
      "  real Intercept{p};\n",
      str_if(length(info$fixef), glue("  vector[Kc{p}] b{p};\n")),
      "  real b{p}_Intercept;\n"
    )
  } else {
    str_add(out$gen_def) <- glue("  vector[{q}] b{p};\n")
  }
  str_add(out$gen_def) <- cglue(
    "  vector[N_{id}] r_{idp}_{r$cn};\n"
  )
  if (is_cor) {
    str_add(out$gen_def) <- glue(
      "  // compute group-level correlations\n",
      "  corr_matrix[M_{id}] Cor_{id}",
      " = multiply_lower_tri_self_transpose(L_{id});\n",
      "  vector<lower=-1,upper=1>[NC_{id}] cor_{id};\n"
    )
  }

  str_add(out$gen_comp) <- glue(
    "  {{\n",
    "    matrix[{q}, M_{id}] prior_factor_s2z = diag_pre_multiply(",
    "sqrt(prior_prec_s2z_{id}), H_s2z_{id}) * L_Sigma_s2z_{id};\n",
    "    vector[{q}] prior_difference_s2z = sqrt(prior_prec_s2z_{id}) .* ",
    "(theta_s2z{p} - prior_mean_s2z_{id});\n",
    "    vector[M_{id}] h_s2z = prior_factor_s2z' * ",
    "prior_difference_s2z{contrast_score};\n",
    "    vector[M_{id}] forward_solve_s2z = mdivide_left_tri_low(",
    "L_P_s2z_{id}, h_s2z);\n",
    "    vector[M_{id}] r_mean_s2z = (mdivide_right_tri_low(",
    "forward_solve_s2z', L_P_s2z_{id}))';\n",
    "    vector[M_{id}] z_mean_s2z;\n",
    "    for (k in 1:M_{id}) z_mean_s2z[k] = std_normal_rng();\n",
    "    mean_r_s2z_{id} = L_Sigma_s2z_{id} * (r_mean_s2z + ",
    "(mdivide_right_tri_low(z_mean_s2z', L_P_s2z_{id}))');\n",
    "  }}\n",
    "  q_recovered_s2z_{id} = theta_s2z{p} - H_s2z_{id} * mean_r_s2z_{id};\n",
    "  r_{id} = r_s2z_{id};\n",
    "  for (j in 1:N_{id}) r_{id}[j] += mean_r_s2z_{id}';\n"
  )
  if (info$center) {
    str_add(out$gen_comp) <- glue(
      "  Intercept{p} = q_recovered_s2z_{id}[1];\n",
      str_if(
        length(info$fixef),
        glue("  b{p} = tail(q_recovered_s2z_{id}, Kc{p});\n",
             "  b{p}_Intercept = Intercept{p} - dot_product(",
             "means_X{p}, b{p});\n")
      ),
      str_if(
        !length(info$fixef),
        glue("  b{p}_Intercept = Intercept{p};\n")
      )
    )
  } else {
    str_add(out$gen_comp) <- glue("  b{p} = q_recovered_s2z_{id};\n")
  }
  str_add(out$gen_comp) <- cglue(
    "  r_{idp}_{r$cn} = r_{id}[, {J}];\n"
  )
  if (is_cor) {
    str_add(out$gen_comp) <- stan_cor_gen_comp(
      cor = glue("cor_{id}"), ncol = glue("M_{id}")
    )
  }
  out
}

# Component-wise specialization for multiple independent group-level effects.
# The omitted means are conditionally independent except for the single
# rank-one coupling induced by brms's centered population intercept.  A
# diagonal-plus-rank-one solve therefore replaces all M x M factorizations.
.stan_re_s2z_independent <- function(id, r, info, normalize, out = list(),
                                     fisher_info = NULL) {
  lpdf <- stan_lpdf_name(normalize)
  stopifnot(
    is.reframe(r), nrow(r) > 1L, !isTRUE(r$cor[1]), all(r$s2z),
    !isTRUE(info$ordinal)
  )
  q <- length(info$qnames)
  M <- nrow(r)
  J <- seq_rows(r)
  p <- info$p
  idp <- paste0(r$id, usc(combine_prefix(check_prefix(r))))
  is_student <- identical(r$dist[1], "student")
  s2z_mode <- re_s2z_center_mode(r)
  s2z_center <- identical(s2z_mode, "centered")
  s2z_fisher <- !is.null(fisher_info)
  s2z_partial <- s2z_mode %in% c("partial", "auto")
  stopifnot(!s2z_fisher || identical(s2z_mode, "auto"))
  r_s2z <- glue("r_s2z_{idp}_{r$cn}")
  r_public <- glue("r_{idp}_{r$cn}")

  if (s2z_partial && !s2z_fisher) {
    str_add(out$data) <- glue(
      "  matrix<lower=0,upper=1>[N_{id}, M_{id}] rho_s2z_{id};",
      "  // fixed numeric centering fractions\n"
    )
    out <- stan_re_s2z_partial_tdata(out, id)
  }

  str_add(out$fun) <- "  #include 'fun_sum_to_zero.stan'\n"
  str_add(out$par) <- glue(
    "  vector[M_{id} * (N_{id} - 1)] z_s2z_{id};",
    if (s2z_partial) {
      "  // partially centered orthonormal independent S2Z coordinates\n"
    } else if (s2z_center) {
      "  // physical orthonormal independent S2Z coordinates\n"
    } else {
      "  // standardized orthonormal independent S2Z coordinates\n"
    }
  )

  # Independent Student-t and Cauchy population priors remain conditionally
  # Gaussian after introducing one scalar mixing variable per coefficient.
  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    if (spec$dist == "student") {
      str_add(out$par) <- glue(
        "  real<lower=0> udf_b_s2z{p}_{k};",
        "  // mixing variable for population coefficient {k}\n"
      )
      str_add(out$tpar_prior) <- glue(
        "  lprior += inv_chi_square_{lpdf}(udf_b_s2z{p}_{k} | ",
        "{stan_s2z_arg_code(spec$df)});\n"
      )
    }
  }

  str_add(out$tpar_def) <- glue(
    "  // component-wise physical S2Z effects of ID {id}\n",
    cglue("  vector[N_{id}] {r_s2z};\n"),
    "  vector[{q}] prior_mean_s2z_{id};\n",
    "  vector<lower=0>[{q}] prior_prec_s2z_{id};\n",
    str_if(
      is_student,
      glue(
        "  vector<lower=0>[N_{id}] group_scale_s2z_{id};\n",
        "  vector<lower=0>[N_{id}] group_prec_s2z_{id};\n"
      )
    ),
    "  vector<lower=0>[M_{id}] D_diag_s2z_{id};\n",
    "  real<lower=0> rank1_info_s2z_{id};\n",
    "  vector[M_{id}] mhat_s2z_{id};\n",
    "  vector[{q}] qhat_s2z_{id};\n",
    "  real<lower=0> group_quad_s2z_{id};\n",
    str_if(
      s2z_partial,
      glue("  real log_det_partial_s2z_{id};\n")
    )
  )

  str_add(out$tdata_def) <- glue(
    "  vector[M_{id}] intercept_map_s2z_{id};\n"
  )
  str_add(out$tdata_comp) <- glue(
    "  intercept_map_s2z_{id} = zeros_vector(M_{id});\n"
  )
  if (info$center) {
    for (j in seq_len(M)) {
      qi <- info$match_q[j]
      value <- if (info$r$coef[j] == "Intercept") {
        "1.0"
      } else {
        glue("means_X{p}[{qi - 1L}]")
      }
      str_add(out$tdata_comp) <- glue(
        "  intercept_map_s2z_{id}[{j}] = {value};\n"
      )
    }
  }

  if (s2z_fisher) {
    str_add(out$gen_comp) <- stan_re_s2z_fisher_gq_comp(
      id, r = r, fisher_info = fisher_info,
      diag_scale = glue("sd_{id}")
    )
  }
  if (s2z_partial) {
    str_add(out$tpar_comp) <- glue(
      "  log_det_partial_s2z_{id} = 0.0;\n"
    )
  }
  for (j in seq_len(M)) {
    if (s2z_partial) {
      str_add(out$tpar_comp) <- stan_re_s2z_partial_independent_transform(
        id, j, r_s2z[j], glue("sd_{id}[{j}]")
      )
    } else {
      str_add(out$tpar_comp) <- glue(
        "  {r_s2z[j]} = sum_to_zero_constrain_brms(",
        str_if(!s2z_center, glue("sd_{id}[{j}] * ")),
        "segment(z_s2z_{id}, ({j} - 1) * (N_{id} - 1) + 1, ",
        "N_{id} - 1));\n"
      )
    }
  }
  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    loc <- stan_s2z_arg_code(spec$location)
    str_add(out$tpar_comp) <- glue(
      "  prior_mean_s2z_{id}[{k}] = {loc};\n"
    )
    if (spec$dist == "flat") {
      prec <- "0.0"
    } else if (spec$dist == "normal") {
      prec <- glue("inv_square({stan_s2z_arg_code(spec$scale)})")
    } else {
      prec <- glue(
        "inv_square({stan_s2z_arg_code(spec$scale)} * sqrt(",
        "{stan_s2z_arg_code(spec$df)} * udf_b_s2z{p}_{k}))"
      )
    }
    str_add(out$tpar_comp) <- glue(
      "  prior_prec_s2z_{id}[{k}] = {prec};\n"
    )
  }

  if (is_student) {
    tr <- subset_reframe_dist(r, "student")
    g <- usc(tr$ggn[1])
    str_add(out$tpar_comp) <- glue(
      "  group_scale_s2z_{id} = dfm{g};\n",
      "  group_prec_s2z_{id} = inv_square(group_scale_s2z_{id});\n"
    )
  }

  # After multiplying each normal equation by its group SD squared, its
  # diagonal is D.  The only remaining coupling is lambda_Intercept * v v',
  # where v maps omitted means into the centered intercept.  Sherman-Morrison
  # gives both the mode and determinant in linear time and avoids unstable
  # inverse-square group SDs during warmup.
  str_add(out$tpar_comp) <- glue(
    "  {{\n",
    "    vector[M_{id}] base_info_s2z = zeros_vector(M_{id});\n",
    "    vector[M_{id}] base_score_s2z = zeros_vector(M_{id});\n",
    "    vector[M_{id}] scaled_score_s2z;\n",
    "    vector[M_{id}] independent_mode_s2z;\n",
    "    real group_info_s2z = ",
    str_if(is_student, glue("sum(group_prec_s2z_{id})"), glue("N_{id}")),
    ";\n"
  )
  for (j in seq_len(M)) {
    qi <- info$match_q[j]
    if (!info$center || qi != 1L) {
      str_add(out$tpar_comp) <- glue(
        "    base_info_s2z[{j}] = prior_prec_s2z_{id}[{qi}];\n",
        "    base_score_s2z[{j}] = prior_prec_s2z_{id}[{qi}] * ",
        "(theta_s2z{p}[{qi}] - prior_mean_s2z_{id}[{qi}]);\n"
      )
    }
  }
  str_add(out$tpar_comp) <- glue(
    "    D_diag_s2z_{id} = group_info_s2z + ",
    "square(sd_{id}) .* base_info_s2z;\n"
  )
  for (j in seq_len(M)) {
    contrast_score <- if (is_student) {
      glue(" - dot_product({r_s2z[j]}, group_prec_s2z_{id})")
    } else {
      ""
    }
    str_add(out$tpar_comp) <- glue(
      "    scaled_score_s2z[{j}] = square(sd_{id}[{j}]) * ",
      "base_score_s2z[{j}]{contrast_score};\n"
    )
  }
  coupling_prec <- str_if(info$center, glue("prior_prec_s2z_{id}[1]"), "0.0")
  coupling_resid <- str_if(
    info$center,
    glue("theta_s2z{p}[1] - prior_mean_s2z_{id}[1]"),
    "0.0"
  )
  str_add(out$tpar_comp) <- glue(
    "    independent_mode_s2z = scaled_score_s2z ./ D_diag_s2z_{id};\n",
    "    rank1_info_s2z_{id} = {coupling_prec} * dot_product(\n",
    "      square(sd_{id}) .* square(intercept_map_s2z_{id}),\n",
    "      1.0 ./ D_diag_s2z_{id}\n",
    "    );\n",
    "    mhat_s2z_{id} = independent_mode_s2z +\n",
    "      {coupling_prec} * square(sd_{id}) .* intercept_map_s2z_{id} ./\n",
    "      D_diag_s2z_{id} * ({coupling_resid} -\n",
    "      dot_product(intercept_map_s2z_{id}, independent_mode_s2z)) /\n",
    "      (1.0 + rank1_info_s2z_{id});\n",
    "  }}\n",
    "  qhat_s2z_{id} = theta_s2z{p};\n"
  )
  if (info$center) {
    str_add(out$tpar_comp) <- glue(
      "  qhat_s2z_{id}[1] -= dot_product(",
      "intercept_map_s2z_{id}, mhat_s2z_{id});\n"
    )
    for (j in seq_len(M)) {
      qi <- info$match_q[j]
      if (qi != 1L) {
        str_add(out$tpar_comp) <- glue(
          "  qhat_s2z_{id}[{qi}] -= mhat_s2z_{id}[{j}];\n"
        )
      }
    }
  } else {
    for (j in seq_len(M)) {
      str_add(out$tpar_comp) <- glue(
        "  qhat_s2z_{id}[{info$match_q[j]}] -= mhat_s2z_{id}[{j}];\n"
      )
    }
  }
  str_add(out$tpar_comp) <- glue(
    "  group_quad_s2z_{id} = 0.0;\n"
  )
  for (j in seq_len(M)) {
    standardized <- glue(
      "({r_s2z[j]} + mhat_s2z_{id}[{j}]) / sd_{id}[{j}]",
      str_if(is_student, glue(" ./ group_scale_s2z_{id}"))
    )
    str_add(out$tpar_comp) <- glue(
      "  group_quad_s2z_{id} += dot_self({standardized});\n"
    )
  }
  str_add(out$pll_args) <- cglue(", vector {r_s2z}")

  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    if (spec$dist == "flat") {
      next
    }
    if (spec$dist == "normal") {
      cond_scale <- stan_s2z_arg_code(spec$scale)
    } else {
      cond_scale <- glue(
        "{stan_s2z_arg_code(spec$scale)} * sqrt(",
        "{stan_s2z_arg_code(spec$df)} * udf_b_s2z{p}_{k})"
      )
    }
    str_add(out$tpar_prior) <- glue(
      "  lprior += normal_{lpdf}(qhat_s2z_{id}[{k}] | ",
      "{stan_s2z_arg_code(spec$location)}, {cond_scale});\n"
    )
  }
  str_add(out$tpar_prior) <- glue(
    "  lprior += -0.5 * group_quad_s2z_{id}\n",
    str_if(
      s2z_center,
      glue("    - (N_{id} - 1) * sum(log(sd_{id}))\n")
    ),
    str_if(
      s2z_partial,
      glue("    + log_det_partial_s2z_{id}\n")
    ),
    str_if(
      is_student,
      glue("    - M_{id} * sum(log(group_scale_s2z_{id}))\n")
    ),
    "    - 0.5 * sum(log(D_diag_s2z_{id}))\n",
    "    - 0.5 * log1p(rank1_info_s2z_{id})",
    str_if(
      normalize,
      glue(" - 0.5 * N_{id} * M_{id} * log(2 * pi())",
           " + 0.5 * M_{id} * log(2 * pi())",
           " + 0.5 * M_{id} * log(1.0 * N_{id})")
    ),
    ";\n"
  )

  # Draw from the diagonal-minus-rank-one conditional covariance without a
  # dense Cholesky, then reconstruct conventional public b/r coordinates.
  str_add(out$gen_def) <- glue(
    "  vector[M_{id}] mean_r_s2z_{id};\n",
    "  vector[{q}] q_recovered_s2z_{id};\n"
  )
  if (info$center) {
    str_add(out$gen_def) <- glue(
      "  real Intercept{p};\n",
      str_if(length(info$fixef), glue("  vector[Kc{p}] b{p};\n")),
      "  real b{p}_Intercept;\n"
    )
  } else {
    str_add(out$gen_def) <- glue("  vector[{q}] b{p};\n")
  }
  str_add(out$gen_def) <- cglue("  vector[N_{id}] {r_public};\n")

  str_add(out$gen_comp) <- glue(
    "  {{\n",
    "    vector[M_{id}] independent_noise_s2z;\n",
    "    real sqrt_rank1_s2z = sqrt(1.0 + rank1_info_s2z_{id});\n",
    "    real rank1_adjust_s2z;\n",
    "    for (k in 1:M_{id}) {{\n",
    "      independent_noise_s2z[k] = sd_{id}[k] * ",
    "std_normal_rng() / sqrt(D_diag_s2z_{id}[k]);\n",
    "    }}\n",
    "    rank1_adjust_s2z = {coupling_prec} / ",
    "(sqrt_rank1_s2z * (1.0 + sqrt_rank1_s2z)) *\n",
    "      dot_product(intercept_map_s2z_{id}, independent_noise_s2z);\n",
    "    mean_r_s2z_{id} = mhat_s2z_{id} + independent_noise_s2z -\n",
    "      rank1_adjust_s2z * square(sd_{id}) .* ",
    "intercept_map_s2z_{id} ./ D_diag_s2z_{id};\n",
    "  }}\n",
    "  q_recovered_s2z_{id} = theta_s2z{p};\n"
  )
  if (info$center) {
    str_add(out$gen_comp) <- glue(
      "  q_recovered_s2z_{id}[1] -= dot_product(",
      "intercept_map_s2z_{id}, mean_r_s2z_{id});\n"
    )
    for (j in seq_len(M)) {
      qi <- info$match_q[j]
      if (qi != 1L) {
        str_add(out$gen_comp) <- glue(
          "  q_recovered_s2z_{id}[{qi}] -= mean_r_s2z_{id}[{j}];\n"
        )
      }
    }
  } else {
    for (j in seq_len(M)) {
      str_add(out$gen_comp) <- glue(
        "  q_recovered_s2z_{id}[{info$match_q[j]}] -= ",
        "mean_r_s2z_{id}[{j}];\n"
      )
    }
  }
  for (j in seq_len(M)) {
    str_add(out$gen_comp) <- glue(
      "  {r_public[j]} = {r_s2z[j]} + mean_r_s2z_{id}[{j}];\n"
    )
  }
  if (info$center) {
    str_add(out$gen_comp) <- glue(
      "  Intercept{p} = q_recovered_s2z_{id}[1];\n",
      str_if(
        length(info$fixef),
        glue("  b{p} = tail(q_recovered_s2z_{id}, Kc{p});\n",
             "  b{p}_Intercept = Intercept{p} - dot_product(",
             "means_X{p}, b{p});\n")
      ),
      str_if(
        !length(info$fixef),
        glue("  b{p}_Intercept = Intercept{p};\n")
      )
    )
  } else {
    str_add(out$gen_comp) <- glue("  b{p} = q_recovered_s2z_{id};\n")
  }
  out
}

# Scalar specialization of the physical sum-to-zero block.
# In addition to avoiding all 1 x 1 matrix algebra, the Gaussian branch uses
# the exact zero sum of the orthonormal contrasts to remove a weighted dot
# product and the group-scale vectors altogether.
.stan_re_s2z_scalar <- function(id, r, info, normalize, out = list(),
                                fisher_info = NULL) {
  lpdf <- stan_lpdf_name(normalize)
  stopifnot(
    is.reframe(r), nrow(r) == 1L, all(r$s2z), !isTRUE(info$ordinal)
  )
  q <- length(info$qnames)
  p <- info$p
  idp <- paste0(r$id, usc(combine_prefix(check_prefix(r))))
  cn <- r$cn[1]
  is_student <- identical(r$dist[1], "student")
  s2z_mode <- re_s2z_center_mode(r)
  s2z_center <- identical(s2z_mode, "centered")
  s2z_fisher <- !is.null(fisher_info)
  s2z_partial <- s2z_mode %in% c("partial", "auto")
  stopifnot(!s2z_fisher || identical(s2z_mode, "auto"))
  r_s2z <- glue("r_s2z_{idp}_{cn}")
  r_public <- glue("r_{idp}_{cn}")

  if (s2z_partial && !s2z_fisher) {
    str_add(out$data) <- glue(
      "  matrix<lower=0,upper=1>[N_{id}, M_{id}] rho_s2z_{id};",
      "  // fixed numeric centering fractions\n"
    )
    out <- stan_re_s2z_partial_tdata(out, id)
  }

  str_add(out$fun) <- "  #include 'fun_sum_to_zero.stan'\n"
  str_add(out$par) <- glue(
    "  vector[N_{id} - 1] z_s2z_{id};",
    if (s2z_partial) {
      "  // partially centered orthonormal scalar S2Z coordinates\n"
    } else if (s2z_center) {
      "  // physical orthonormal scalar S2Z coordinates\n"
    } else {
      "  // standardized orthonormal scalar S2Z coordinates\n"
    }
  )

  # Independent Student-t and Cauchy population priors remain conditionally
  # Gaussian after introducing one scalar mixing variable per coefficient.
  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    if (spec$dist == "student") {
      str_add(out$par) <- glue(
        "  real<lower=0> udf_b_s2z{p}_{k};",
        "  // mixing variable for population coefficient {k}\n"
      )
      str_add(out$tpar_prior) <- glue(
        "  lprior += inv_chi_square_{lpdf}(udf_b_s2z{p}_{k} | ",
        "{stan_s2z_arg_code(spec$df)});\n"
      )
    }
  }

  str_add(out$tpar_def) <- glue(
    "  // specialized scalar physical S2Z effects of ID {id}\n",
    "  vector[N_{id}] {r_s2z};\n",
    "  vector[{q}] prior_mean_s2z_{id};\n",
    "  vector<lower=0>[{q}] prior_prec_s2z_{id};\n",
    str_if(
      is_student,
      glue(
        "  vector<lower=0>[N_{id}] group_scale_s2z_{id};\n",
        "  vector<lower=0>[N_{id}] group_prec_s2z_{id};\n"
      )
    ),
    "  real<lower=0> D_s2z_{id};\n",
    "  real<lower=0> sqrt_D_s2z_{id};\n",
    "  real mhat_s2z_{id};\n",
    "  vector[{q}] qhat_s2z_{id};\n",
    "  real<lower=0> group_quad_s2z_{id};\n",
    str_if(
      s2z_partial,
      glue("  real log_det_partial_s2z_{id};\n")
    )
  )

  # H maps the omitted scalar group mean into all matching population
  # coordinates. A centered varying slope also shifts brms's temporary
  # centered intercept by the mean of its raw design column.
  qi <- info$match_q[1]
  str_add(out$tdata_def) <- glue(
    "  vector[{q}] H_s2z_{id};\n"
  )
  str_add(out$tdata_comp) <- glue(
    "  H_s2z_{id} = zeros_vector({q});\n",
    "  H_s2z_{id}[{qi}] = 1.0;\n",
    str_if(
      info$center && info$r$coef[1] != "Intercept",
      glue("  H_s2z_{id}[1] = means_X{p}[{qi - 1L}];\n")
    )
  )
  if (s2z_fisher) {
    str_add(out$gen_comp) <- stan_re_s2z_fisher_gq_comp(
      id, r = r, fisher_info = fisher_info,
      scale = glue("sd_{id}[1]")
    )
  }
  if (s2z_partial) {
    str_add(out$tpar_comp) <- glue(
      "  log_det_partial_s2z_{id} = 0.0;\n"
    )
    str_add(out$tpar_comp) <- stan_re_s2z_partial_independent_transform(
      id, 1, r_s2z, glue("sd_{id}[1]"), glue("z_s2z_{id}")
    )
  } else {
    str_add(out$tpar_comp) <- glue(
      "  {r_s2z} = sum_to_zero_constrain_brms(",
      str_if(!s2z_center, glue("sd_{id}[1] * ")),
      "z_s2z_{id});\n"
    )
  }
  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    loc <- stan_s2z_arg_code(spec$location)
    str_add(out$tpar_comp) <- glue(
      "  prior_mean_s2z_{id}[{k}] = {loc};\n"
    )
    if (spec$dist == "flat") {
      prec <- "0.0"
    } else if (spec$dist == "normal") {
      prec <- glue("inv_square({stan_s2z_arg_code(spec$scale)})")
    } else {
      prec <- glue(
        "inv_square({stan_s2z_arg_code(spec$scale)} * sqrt(",
        "{stan_s2z_arg_code(spec$df)} * udf_b_s2z{p}_{k}))"
      )
    }
    str_add(out$tpar_comp) <- glue(
      "  prior_prec_s2z_{id}[{k}] = {prec};\n"
    )
  }

  # D = square(sd) * P is the scaled conditional precision of the omitted
  # mean. Working with D avoids forming inv_square(sd), which can overflow
  # near zero during warmup. For Gaussian effects, sum(r_s2z) is exactly zero
  # in the physical coordinates, so its contribution to the mode vanishes.
  if (is_student) {
    tr <- subset_reframe_dist(r, "student")
    g <- usc(tr$ggn[1])
    str_add(out$tpar_comp) <- glue(
      "  group_scale_s2z_{id} = dfm{g};\n",
      "  group_prec_s2z_{id} = inv_square(group_scale_s2z_{id});\n",
      "  {{\n",
      "    real tau_sq_s2z = square(sd_{id}[1]);\n",
      "    real prior_info_s2z = dot_product(prior_prec_s2z_{id}, ",
      "square(H_s2z_{id}));\n",
      "    real prior_score_s2z = dot_product(H_s2z_{id}, ",
      "prior_prec_s2z_{id} .* (theta_s2z{p} - ",
      "prior_mean_s2z_{id}));\n",
      "    D_s2z_{id} = tau_sq_s2z * prior_info_s2z + ",
      "sum(group_prec_s2z_{id});\n",
      "    mhat_s2z_{id} = (tau_sq_s2z * prior_score_s2z - ",
      "dot_product({r_s2z}, group_prec_s2z_{id})) / D_s2z_{id};\n",
      "  }}\n"
    )
  } else {
    str_add(out$tpar_comp) <- glue(
      "  {{\n",
      "    real tau_sq_s2z = square(sd_{id}[1]);\n",
      "    real prior_info_s2z = dot_product(prior_prec_s2z_{id}, ",
      "square(H_s2z_{id}));\n",
      "    real prior_score_s2z = dot_product(H_s2z_{id}, ",
      "prior_prec_s2z_{id} .* (theta_s2z{p} - ",
      "prior_mean_s2z_{id}));\n",
      "    D_s2z_{id} = tau_sq_s2z * prior_info_s2z + N_{id};\n",
      "    mhat_s2z_{id} = tau_sq_s2z * prior_score_s2z / D_s2z_{id};\n",
      "  }}\n"
    )
  }
  str_add(out$tpar_comp) <- glue(
    "  sqrt_D_s2z_{id} = sqrt(D_s2z_{id});\n",
    "  qhat_s2z_{id} = theta_s2z{p} - H_s2z_{id} * ",
    "mhat_s2z_{id};\n",
    "  {{\n",
    "    vector[N_{id}] white_s2z = ({r_s2z} + mhat_s2z_{id}) / ",
    "sd_{id}[1]",
    str_if(is_student, glue(" ./ group_scale_s2z_{id}")),
    ";\n",
    "    group_quad_s2z_{id} = dot_self(white_s2z);\n",
    "  }}\n"
  )
  str_add(out$pll_args) <- glue(", vector {r_s2z}")

  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    if (spec$dist == "flat") {
      next
    }
    if (spec$dist == "normal") {
      cond_scale <- stan_s2z_arg_code(spec$scale)
    } else {
      cond_scale <- glue(
        "{stan_s2z_arg_code(spec$scale)} * sqrt(",
        "{stan_s2z_arg_code(spec$df)} * udf_b_s2z{p}_{k})"
      )
    }
    str_add(out$tpar_prior) <- glue(
      "  lprior += normal_{lpdf}(qhat_s2z_{id}[{k}] | ",
      "{stan_s2z_arg_code(spec$location)}, {cond_scale});\n"
    )
  }
  str_add(out$tpar_prior) <- glue(
    "  lprior += -0.5 * group_quad_s2z_{id}\n",
    str_if(
      s2z_center,
      glue("    - (N_{id} - 1) * log(sd_{id}[1])\n")
    ),
    str_if(
      s2z_partial,
      glue("    + log_det_partial_s2z_{id}\n")
    ),
    str_if(
      is_student,
      glue("    - sum(log(group_scale_s2z_{id}))\n")
    ),
    "    - 0.5 * log(D_s2z_{id})",
    str_if(
      normalize,
      glue(" - 0.5 * N_{id} * log(2 * pi())",
           " + 0.5 * log(2 * pi())",
           " + 0.5 * log(1.0 * N_{id})")
    ),
    ";\n"
  )

  # Reconstruct the conventional coefficient and deviations from one draw
  # of the analytically omitted scalar group mean.
  str_add(out$gen_def) <- glue(
    "  real mean_r_s2z_{id};\n",
    "  vector[{q}] q_recovered_s2z_{id};\n"
  )
  if (info$center) {
    str_add(out$gen_def) <- glue(
      "  real Intercept{p};\n",
      str_if(length(info$fixef), glue("  vector[Kc{p}] b{p};\n")),
      "  real b{p}_Intercept;\n"
    )
  } else {
    str_add(out$gen_def) <- glue("  vector[{q}] b{p};\n")
  }
  str_add(out$gen_def) <- glue("  vector[N_{id}] {r_public};\n")

  str_add(out$gen_comp) <- glue(
    "  mean_r_s2z_{id} = mhat_s2z_{id} + ",
    "sd_{id}[1] * std_normal_rng() / sqrt_D_s2z_{id};\n",
    "  q_recovered_s2z_{id} = theta_s2z{p} - H_s2z_{id} * ",
    "mean_r_s2z_{id};\n",
    "  {r_public} = {r_s2z} + mean_r_s2z_{id};\n"
  )
  if (info$center) {
    str_add(out$gen_comp) <- glue(
      "  Intercept{p} = q_recovered_s2z_{id}[1];\n",
      str_if(
        length(info$fixef),
        glue("  b{p} = tail(q_recovered_s2z_{id}, Kc{p});\n",
             "  b{p}_Intercept = Intercept{p} - dot_product(",
             "means_X{p}, b{p});\n")
      ),
      str_if(
        !length(info$fixef),
        glue("  b{p}_Intercept = Intercept{p};\n")
      )
    )
  } else {
    str_add(out$gen_comp) <- glue("  b{p} = q_recovered_s2z_{id};\n")
  }
  out
}

# Render one validated scalar population-prior argument. Symbolic arguments are
# checked at runtime because transformed-data and parameter expressions do not
# necessarily have a value that R can verify before sampling.
stan_s2z_arg_code <- function(x) {
  stopifnot(is.re_s2z_prior_arg(x))
  if (!nzchar(x$code)) {
    stopifnot(x$known)
    return(stan_s2z_number(x$value))
  }
  if (identical(x$role, "location")) {
    as.character(glue("s2z_require_finite_brms({x$code})"))
  } else {
    as.character(glue("s2z_require_positive_brms({x$code})"))
  }
}

# Stable formatting for numeric constants inserted into generated Stan code.
stan_s2z_number <- function(x) {
  stopifnot(length(x) == 1L, is.numeric(x), is.finite(x))
  trimws(formatC(x, digits = 17, format = "g"))
}

# Stan code of smooth terms
stan_sm <- function(bframe, prior, threads, normalize, ...) {
  stopifnot(is.bframel(bframe))
  lpdf <- ifelse(normalize, "lpdf", "lupdf")
  out <- list()
  smframe <- bframe$frame$sm
  if (!has_rows(smframe)) {
    return(out)
  }
  px <- check_prefix(bframe)
  p <- usc(combine_prefix(px))
  resp <- usc(px$resp)
  slice <- stan_slice(threads)
  Xs_names <- attr(smframe, "Xs_names")
  if (length(Xs_names)) {
    str_add(out$data) <- glue(
      "  // data for splines\n",
      "  int Ks{p};  // number of linear effects\n",
      "  matrix[N{resp}, Ks{p}] Xs{p};",
      "  // design matrix for the linear effects\n"
    )
    str_add(out$pll_args) <- glue(", data matrix Xs{p}")
    if (has_special_prior(prior, px, class = "b")) {
      str_add_list(out) <- stan_prior_non_centered(
        suffix = glue("s{p}"), normalize = normalize
      )
    } else {
      str_add_list(out) <- stan_prior(
        prior, class = "b", coef = Xs_names,
        type = glue("vector[Ks{p}]"), suffix = glue("s{p}"),
        header_type = "vector", px = px,
        comment = "unpenalized spline coefficients",
        normalize = normalize
      )
    }
    str_add(out$eta) <- glue(" + Xs{p}{slice} * bs{p}")
  }
  for (i in seq_rows(smframe)) {
    if (smframe$nbases[[i]] == 0) {
      next  # no penalized spline components present
    }
    pi <- glue("{p}_{i}")
    nb <- seq_len(smframe$nbases[[i]])
    str_add(out$data) <- glue(
      "  // data for spline {i}\n",
      "  int nb{pi};  // number of bases\n",
      "  array[nb{pi}] int knots{pi};  // number of knots\n"
    )
    str_add(out$data) <- "  // basis function matrices\n"
    str_add(out$data) <- cglue(
      "  matrix[N{resp}, knots{pi}[{nb}]] Zs{pi}_{nb};\n"
    )
    str_add(out$pll_args) <- cglue(", data matrix Zs{pi}_{nb}")
    str_add(out$par) <- glue(
      "  // parameters for spline {i}\n"
    )
    str_add(out$par) <- cglue(
      "  // standardized penalized spline coefficients\n",
      "  vector[knots{pi}[{nb}]] zs{pi}_{nb};\n"
    )
    if (has_special_prior(prior, px, class = "sds")) {
      str_add(out$tpar_def) <- glue(
        "  // SDs of penalized spline coefficients\n",
        "  vector<lower=0>[nb{pi}] sds{pi};\n"
      )
      str_add(out$prior_global_scales) <- glue(" sds{pi}")
      str_add(out$prior_global_lengths) <- glue(" nb{pi}")
    } else {
      str_add_list(out) <- stan_prior(
        prior, class = "sds", coef = smframe$term[i], suffix = pi, px = px,
        type = glue("vector[nb{pi}]"), coef_type = glue("vector[nb{pi}]"),
        comment = "SDs of penalized spline coefficients",
        normalize = normalize
      )
    }
    # separate definition from computation to support fixed parameters
    str_add(out$tpar_def) <- cglue(
      "  // penalized spline coefficients\n",
      "  vector[knots{pi}[{nb}]] s{pi}_{nb};\n"
    )
    str_add(out$tpar_special_prior) <- cglue(
      "  // compute penalized spline coefficients\n",
      "  s{pi}_{nb} = sds{pi}[{nb}] * zs{pi}_{nb};\n"
    )
    str_add(out$pll_args) <- cglue(", vector s{pi}_{nb}")
    str_add(out$model_prior) <- cglue(
      "  target += std_normal_{lpdf}(zs{pi}_{nb});\n"
    )
    str_add(out$eta) <- cglue(
      " + Zs{pi}_{nb}{slice} * s{pi}_{nb}"
    )
  }
  out
}

# Stan code for category specific effects
# @note not implemented for non-linear models
stan_cs <- function(bframe, prior, threads, normalize, ...) {
  stopifnot(is.bframel(bframe))
  out <- list()
  csef <- bframe$frame$cs$vars
  px <- check_prefix(bframe)
  p <- usc(combine_prefix(px))
  resp <- usc(bframe$resp)
  slice <- stan_slice(threads)
  reframe <- subset2(bframe$frame$re, type = "cs")
  if (length(csef)) {
    str_add(out$data) <- glue(
      "  int<lower=1> Kcs{p};  // number of category specific effects\n",
      "  matrix[N{resp}, Kcs{p}] Xcs{p};  // category specific design matrix\n"
    )
    str_add(out$pll_args) <- glue(", data matrix Xcs{p}")
    str_add_list(out) <- stan_prior(
      prior, class = "b", coef = csef,
      type = glue("matrix[Kcs{p}, nthres{resp}]"),
      coef_type = glue("row_vector[nthres{resp}]"),
      suffix = glue("cs{p}"), px = px, broadcast = "matrix",
      header_type = "matrix", comment = "category specific effects",
      normalize = normalize
    )
    str_add(out$model_def) <- glue(
      "  // linear predictor for category specific effects\n",
      "  matrix[N{resp}, nthres{resp}] mucs{p} = Xcs{p}{slice} * bcs{p};\n"
    )
  }
  if (has_rows(reframe)) {
    if (!length(csef)) {
      # only group-level category specific effects present
      str_add(out$model_def) <- glue(
        "  // linear predictor for category specific effects\n",
        "  matrix[N{resp}, nthres{resp}] mucs{p}",
        " = rep_matrix(0, N{resp}, nthres{resp});\n"
      )
    }
    n <- stan_nn(threads)
    thres_regex <- "(?<=\\[)[[:digit:]]+(?=\\]$)"
    thres <- get_matches(thres_regex, reframe$coef, perl = TRUE)
    nthres <- max(as.numeric(thres))
    mucs_loop <- ""
    for (i in seq_len(nthres)) {
      r_cat <- reframe[grepl(glue("\\[{i}\\]$"), reframe$coef), ]
      str_add(mucs_loop) <- glue(
        "    mucs{p}[n, {i}] = mucs{p}[n, {i}]"
      )
      for (id in unique(r_cat$id)) {
        r <- r_cat[r_cat$id == id, ]
        rpx <- check_prefix(r)
        idp <- paste0(r$id, usc(combine_prefix(rpx)))
        idresp <- paste0(r$id, usc(rpx$resp))
        str_add(mucs_loop) <- cglue(
          " + r_{idp}_{r$cn}[J_{idresp}{n}] * Z_{idp}_{r$cn}{n}"
        )
      }
      str_add(mucs_loop) <- ";\n"
    }
    str_add(out$model_comp_eta_loop) <- glue(
      "  for (n in 1:N{resp}) {{\n",
      stan_nn_def(threads), mucs_loop,
      "  }\n"
    )
  }
  out
}

# Stan code for special effects
stan_sp <- function(bframe, prior, stanvars, threads, normalize, ...) {
  stopifnot(is.bframel(bframe))
  out <- list()
  spframe <- bframe$frame$sp
  reframe <- bframe$frame$re
  meframe <- bframe$frame$me
  if (!has_rows(spframe)) {
    return(out)
  }
  px <- check_prefix(bframe)
  p <- usc(combine_prefix(px))
  resp <- usc(px$resp)
  lpdf <- stan_lpdf_name(normalize)
  n <- stan_nn(threads)
  reframe <- subset2(reframe, type = "sp")
  spframe_coef <- rename(spframe$term)
  invalid_coef <- setdiff(reframe$coef, spframe_coef)
  if (length(invalid_coef)) {
    stop2(
      "Special group-level terms require corresponding ",
      "population-level terms:\nOccured for ",
      collapse_comma(invalid_coef)
    )
  }
  # prepare Stan code of the linear predictor component
  for (i in seq_rows(spframe)) {
    eta <- spframe$joint_call[[i]]
    if (!is.null(spframe$calls_mo[[i]])) {
      new_mo <- glue("mo(simo{p}_{spframe$Imo[[i]]}, Xmo{p}_{spframe$Imo[[i]]}{n})")
      eta <- rename(eta, spframe$calls_mo[[i]], new_mo)
    }
    if (!is.null(spframe$calls_me[[i]])) {
      Kme <- seq_along(meframe$term)
      Ime <- match(meframe$grname, unique(meframe$grname))
      nme <- ifelse(nzchar(meframe$grname), glue("[Jme_{Ime}{n}]"), n)
      new_me <- glue("Xme_{Kme}{nme}")
      eta <- rename(eta, meframe$term, new_me)
    }
    if (!is.null(spframe$calls_mi[[i]])) {
      is_na_idx <- is.na(spframe$idx2_mi[[i]])
      idx_mi <- glue("[idxl{p}_{spframe$vars_mi[[i]]}_{spframe$idx2_mi[[i]]}{n}]")
      idx_mi <- ifelse(is_na_idx, n, idx_mi)
      new_mi <- glue("Yl_{spframe$vars_mi[[i]]}{idx_mi}")
      eta <- rename(eta, spframe$calls_mi[[i]], new_mi)
      str_add(out$pll_args) <- glue(", vector Yl_{spframe$vars_mi[[i]]}")
    }
    if (spframe$Ic[i] > 0) {
      str_add(eta) <- glue(" * Csp{p}_{spframe$Ic[i]}{n}")
    }
    r <- subset2(reframe, coef = spframe_coef[i])
    rpars <- str_if(nrow(r), cglue(" + {stan_eta_rsp(r, threads)}"))
    str_add(out$loopeta) <- glue(" + (bsp{p}[{i}]{rpars}) * {eta}")
  }

  # prepare general Stan code
  ncovars <- max(spframe$Ic)
  str_add(out$data) <- glue(
    "  int<lower=1> Ksp{p};  // number of special effects terms\n"
  )
  if (ncovars > 0L) {
    str_add(out$data) <- "  // covariates of special effects terms\n"
    str_add(out$data) <- cglue(
      "  vector[N{resp}] Csp{p}_{seq_len(ncovars)};\n"
    )
    str_add(out$pll_args) <- cglue(", data vector Csp{p}_{seq_len(ncovars)}")
  }

  # include special Stan code for monotonic effects
  which_Imo <- which(lengths(spframe$Imo) > 0)
  if (any(which_Imo)) {
    str_add(out$fun) <- "  #include 'fun_monotonic.stan'\n"
    str_add(out$data) <- glue(
      "  int<lower=1> Imo{p};  // number of monotonic variables\n",
      "  array[Imo{p}] int<lower=1> Jmo{p};  // length of simplexes\n"
    )
    ids <- unlist(spframe$ids_mo)
    lpdf <- stan_lpdf_name(normalize)
    for (i in which_Imo) {
      for (k in seq_along(spframe$Imo[[i]])) {
        j <- spframe$Imo[[i]][[k]]
        id <- spframe$ids_mo[[i]][[k]]
        # index of first ID appearance
        j_id <- match(id, ids)
        str_add(out$data) <- glue(
          "  array[N{resp}] int Xmo{p}_{j};  // monotonic variable\n"
        )
        str_add(out$pll_args) <- glue(
          ", array[] int Xmo{p}_{j}, vector simo{p}_{j}"
        )
        if (is.na(id) || j_id == j) {
          # no ID or first appearance of the ID
          str_add(out$data) <- glue(
            "  vector[Jmo{p}[{j}]] con_simo{p}_{j};",
            "  // prior concentration of monotonic simplex\n"
          )
          str_add(out$par) <- glue(
            "  simplex[Jmo{p}[{j}]] simo{p}_{j};  // monotonic simplex\n"
          )
          str_add(out$tpar_prior) <- glue(
            "  lprior += dirichlet_{lpdf}(simo{p}_{j} | con_simo{p}_{j});\n"
          )
        } else {
          # use the simplex shared across all terms of the same ID
          str_add(out$tpar_def) <- glue(
            "  simplex[Jmo{p}[{j}]] simo{p}_{j} = simo{p}_{j_id};\n"
          )
        }
      }
    }
  }

  # include special Stan code for missing value terms
  uni_mi <- na.omit(attr(spframe, "uni_mi"))
  for (j in seq_rows(uni_mi)) {
    idxl <- glue("idxl{p}_{uni_mi$var[j]}_{uni_mi$idx2[j]}")
    str_add(out$data) <- glue(
      "  array[N{resp}] int {idxl};  // matching indices\n"
    )
    str_add(out$pll_args) <- glue(", data array[] int {idxl}")
  }

  # prepare special effects coefficients
  if (has_special_prior(prior, bframe, class = "b")) {
    stopif_prior_bound(prior, class = "b", ls = px)
    str_add_list(out) <- stan_prior_non_centered(
      suffix = glue("sp{p}"), normalize = normalize
    )
  } else {
    str_add_list(out) <- stan_prior(
      prior, class = "b", coef = spframe$coef,
      type = glue("vector[Ksp{p}]"), px = px,
      suffix = glue("sp{p}"), header_type = "vector",
      comment = "special effects coefficients",
      normalize = normalize
    )
  }
  out
}

# Stan code for latent gaussian processes
stan_gp <- function(bframe, prior, threads, normalize, ...) {
  stopifnot(is.bframel(bframe))
  lpdf <- stan_lpdf_name(normalize)
  out <- list()
  px <- check_prefix(bframe)
  p <- usc(combine_prefix(px))
  resp <- usc(px$resp)
  slice <- stan_slice(threads)
  gpframe <- bframe$frame$gp
  # kernel methods cannot simply be split up into partial sums
  for (i in seq_rows(gpframe)) {
    pi <- glue("{p}_{i}")
    byvar <- gpframe$byvars[[i]]
    cons <- gpframe$cons[[i]]
    byfac <- length(cons) > 0L
    bynum <- !is.null(byvar) && !byfac
    k <- gpframe$k[i]
    is_approx <- !isNA(k)
    iso <- gpframe$iso[i]
    gr <- gpframe$gr[i]
    cov <- gpframe$cov[i]
    gp <- glue("gp_{cov}")
    sfx1 <- gpframe$sfx1[[i]]
    sfx2 <- gpframe$sfx2[[i]]
    str_add(out$data) <- glue(
      "  // data related to GPs\n",
      "  int<lower=1> Kgp{pi};",
      "  // number of sub-GPs (equal to 1 unless 'by' was used)\n",
      "  int<lower=1> Dgp{pi};  // GP dimension\n"
    )
    if (is_approx) {
      str_add(out$fun) <- glue("  #include 'fun_spd_gp_{cov}.stan'\n")
      str_add(out$data) <- glue(
        "  // number of basis functions of an approximate GP\n",
        "  int<lower=1> NBgp{pi};\n"
      )
    } else {
      str_add(out$fun) <- glue("  #include 'fun_gp_{cov}.stan'\n")
    }
    if (has_special_prior(prior, px, class = "sdgp")) {
      str_add(out$tpar_def) <- glue(
        "  vector<lower=0>[Kgp{pi}] sdgp{pi};  // GP standard deviation parameters\n"
      )
      str_add(out$prior_global_scales) <- glue(" sdgp{pi}")
      str_add(out$prior_global_lengths) <- glue(" Kgp{pi}")
    } else {
      str_add_list(out) <- stan_prior(
        prior, class = "sdgp", coef = sfx1, px = px, suffix = pi,
        type = glue("vector[Kgp{pi}]"), coef_type = glue("vector[Kgp{pi}]"),
        comment = "GP standard deviation parameters",
        normalize = normalize
      )
    }
    if (gpframe$iso[i]) {
      lscale_type <- "vector[1]"
      lscale_dim <- glue("[Kgp{pi}]")
      lscale_comment <- "GP length-scale parameters"
    } else {
      lscale_type <- glue("vector[Dgp{pi}]")
      lscale_dim <- glue("[Kgp{pi}]")
      lscale_comment <- "GP length-scale parameters"
    }
    if (byfac) {
      J <- seq_along(cons)
      Ngp <- glue("Ngp{pi}")
      Nsubgp <- glue("N", str_if(gr, "sub"), glue("gp{pi}"))
      Igp <- glue("Igp{pi}_{J}")
      str_add(out$data) <- glue(
        "  // number of observations relevant for a certain sub-GP\n",
        "  array[Kgp{pi}] int<lower=1> {Ngp};\n"
      )
      str_add(out$data) <-
        "  // indices and contrasts of sub-GPs per observation\n"
      str_add(out$data) <- cglue(
        "  array[{Ngp}[{J}]] int<lower=1> {Igp};\n",
        "  vector[{Ngp}[{J}]] Cgp{pi}_{J};\n"
      )
      str_add(out$pll_args) <- cglue(
        ", data array[] int {Igp}, data vector Cgp{pi}_{J}"
      )
      str_add_list(out) <- stan_prior(
        prior, class = "lscale", coef = sfx2,
        type = lscale_type, dim = lscale_dim, suffix = glue("{pi}"),
        px = px, comment = lscale_comment, normalize = normalize
      )
      if (gr) {
        str_add(out$data) <- glue(
          "  // number of latent GP groups\n",
          "  array[Kgp{pi}] int<lower=1> Nsubgp{pi};\n"
        )
        str_add(out$data) <- cglue(
          "  // indices of latent GP groups per observation\n",
          "  array[{Ngp}[{J}]] int<lower=1> Jgp{pi}_{J};\n"
        )
        str_add(out$pll_args) <- cglue(", data array[] int Jgp{pi}_{J}")
      }
      if (is_approx) {
        str_add(out$data) <-
          "  // approximate GP basis matrices and eigenvalues\n"
        str_add(out$data) <- cglue(
          "  matrix[{Nsubgp}[{J}], NBgp{pi}] Xgp{pi}_{J};\n",
          "  array[NBgp{pi}] vector[Dgp{pi}] slambda{pi}_{J};\n"
        )
        str_add(out$par) <- "  // latent variables of the GP\n"
        str_add(out$par) <- cglue(
          "  vector[NBgp{pi}] zgp{pi}_{J};\n"
        )
        str_add(out$model_no_pll_def) <- "  // scale latent variables of the GP\n"
        str_add(out$model_no_pll_def) <- cglue(
          "  vector[NBgp{pi}] rgp{pi}_{J} = sqrt(spd_gp_{cov}(",
          "slambda{pi}_{J}, sdgp{pi}[{J}], lscale{pi}[{J}])) .* zgp{pi}_{J};\n"
        )
        gp_call <- glue("Xgp{pi}_{J} * rgp{pi}_{J}")
      } else {
        # exact GPs
        str_add(out$data) <- "  // covariates of the GP\n"
        str_add(out$data) <- cglue(
          "  array[{Nsubgp}[{J}]] vector[Dgp{pi}] Xgp{pi}_{J};\n"
        )
        str_add(out$par) <- "  // latent variables of the GP\n"
        str_add(out$par) <- cglue(
          "  vector[{Nsubgp}[{J}]] zgp{pi}_{J};\n"
        )
        gp_call <- glue(
          "gp_{cov}(Xgp{pi}_{J}, sdgp{pi}[{J}], lscale{pi}[{J}], zgp{pi}_{J})"
        )
      }
      slice2 <- ""
      Igp_sub <- Igp
      if (use_threading(threads)) {
        str_add(out$fun) <- "  #include 'fun_which_range.stan'\n"
        str_add(out$model_comp_basic) <- cglue(
          "  array[size_range({Igp}, start, end)] int which_gp{pi}_{J} =",
          " which_range({Igp}, start, end);\n"
        )
        slice2 <- glue("[which_gp{pi}_{J}]")
        Igp_sub <- glue("start_at_one({Igp}{slice2}, start)")
      }
      # TODO: add all GP elements to 'eta' at the same time?
      eta <- combine_prefix(px, keep_mu = TRUE, nlp = TRUE)
      eta <- glue("{eta}[{Igp_sub}]")
      str_add(out$model_no_pll_def) <- cglue(
        "  vector[{Nsubgp}[{J}]] gp_pred{pi}_{J} = {gp_call};\n"
      )
      str_add(out$pll_args) <- cglue(", vector gp_pred{pi}_{J}")
      Cgp <- glue("Cgp{pi}_{J}{slice2} .* ")
      Jgp <- str_if(gr, glue("[Jgp{pi}_{J}{slice2}]"), slice)
      str_add(out$model_comp_basic) <- cglue(
        "  {eta} += {Cgp}gp_pred{pi}_{J}{Jgp};\n"
      )
      str_add(out$model_prior) <- cglue(
        "{tp()}std_normal_{lpdf}(zgp{pi}_{J});\n"
      )
    } else {
      # no by-factor variable
      str_add_list(out) <- stan_prior(
        prior, class = "lscale", coef = sfx2,
        type = lscale_type, dim = lscale_dim, suffix = glue("{pi}"),
        px = px, comment = lscale_comment, normalize = normalize
      )
      Nsubgp <- glue("N{resp}")
      if (gr) {
        Nsubgp <- glue("Nsubgp{pi}")
        str_add(out$data) <- glue(
          "  // number of latent GP groups\n",
          "  int<lower=1> {Nsubgp};\n",
          "  // indices of latent GP groups per observation\n",
          "  array[N{resp}] int<lower=1> Jgp{pi};\n"
        )
        str_add(out$pll_args) <- glue(", data array[] int Jgp{pi}")
      }
      Cgp <- ""
      if (bynum) {
        str_add(out$data) <- glue(
          "  // numeric by-variable of the GP\n",
          "  vector[N{resp}] Cgp{pi};\n"
        )
        str_add(out$pll_args) <- glue(", data vector Cgp{pi}")
        Cgp <- glue("Cgp{pi}{slice} .* ")
      }
      if (is_approx) {
        str_add(out$data) <- glue(
          "  // approximate GP basis matrices\n",
          "  matrix[{Nsubgp}, NBgp{pi}] Xgp{pi};\n",
          "  // approximate GP eigenvalues\n",
          "  array[NBgp{pi}] vector[Dgp{pi}] slambda{pi};\n"
        )
        str_add(out$par) <- glue(
          "  vector[NBgp{pi}] zgp{pi};  // latent variables of the GP\n"
        )
        str_add(out$model_no_pll_def) <- glue(
          "  // scale latent variables of the GP\n",
          "  vector[NBgp{pi}] rgp{pi} = sqrt(spd_gp_{cov}(",
          "slambda{pi}, sdgp{pi}[1], lscale{pi}[1])) .* zgp{pi};\n"
        )
        if (gr) {
          # grouping prevents GPs to be computed efficiently inside reduce_sum
          str_add(out$model_no_pll_def) <- glue(
            "  vector[{Nsubgp}] gp_pred{pi} = Xgp{pi} * rgp{pi};\n"
          )
          str_add(out$eta) <- glue(" + {Cgp}gp_pred{pi}[Jgp{pi}{slice}]")
          str_add(out$pll_args) <- glue(", vector gp_pred{pi}")
        } else {
          # efficient computation of approx GPs inside reduce_sum is possible
          str_add(out$model_def) <- glue(
            "  vector[N{resp}] gp_pred{pi} = Xgp{pi}{slice} * rgp{pi};\n"
          )
          str_add(out$eta) <- glue(" + {Cgp}gp_pred{pi}")
          str_add(out$pll_args) <- glue(", data matrix Xgp{pi}, vector rgp{pi}")
        }
      } else {
        # exact GPs
        str_add(out$data) <- glue(
          "  array[{Nsubgp}] vector[Dgp{pi}] Xgp{pi};  // covariates of the GP\n"
        )
        str_add(out$par) <- glue(
          "  vector[{Nsubgp}] zgp{pi};  // latent variables of the GP\n"
        )
        gp_call <- glue("gp_{cov}(Xgp{pi}, sdgp{pi}[1], lscale{pi}[1], zgp{pi})")
        # exact GPs are kernel based methods which
        # need to be computed outside of reduce_sum
        str_add(out$model_no_pll_def) <- glue(
          "  vector[{Nsubgp}] gp_pred{pi} = {gp_call};\n"
        )
        Jgp <- str_if(gr, glue("[Jgp{pi}{slice}]"), slice)
        str_add(out$eta) <- glue(" + {Cgp}gp_pred{pi}{Jgp}")
        str_add(out$pll_args) <- glue(", vector gp_pred{pi}")
      }
      str_add(out$model_prior) <- glue(
        "{tp()}std_normal_{lpdf}(zgp{pi});\n"
      )
    }
  }
  out
}

# Stan code for the linear predictor of autocorrelation terms
stan_ac <- function(bframe, prior, threads, normalize, ...) {
  lpdf <- stan_lpdf_name(normalize)
  out <- list()
  px <- check_prefix(bframe)
  p <- usc(combine_prefix(px))
  resp <- usc(px$resp)
  n <- stan_nn(threads)
  slice <- stan_slice(threads)
  acframe <- bframe$frame$ac
  stopifnot(is.acframe(acframe))
  has_natural_residuals <- has_ac_natural_residuals(acframe)
  has_latent_residuals <- has_ac_latent_residuals(acframe)
  families <- family_names(bframe)
  # TODO: include family-specific functions inside the corresponding
  # stan_log_lik functions once they return lists of character vectors

  if (has_latent_residuals) {
    # families that do not have natural residuals require latent
    # residuals for residual-based autocor structures
    err_msg <- "Latent residuals are not implemented"
    if (is.btnl(bframe)) {
      stop2(err_msg, " for non-linear models.")
    }
    str_add(out$par) <- glue(
      "  vector[N{resp}] zerr{p};  // unscaled residuals\n"
    )
    if (has_special_prior(prior, px, class = "sderr")) {
      str_add(out$tpar_def) <- glue(
        "  real<lower=0> sderr{p};  // SD of residuals\n"
      )
      str_add(out$prior_global_scales) <- glue(" sderr{p}")
      str_add(out$prior_global_lengths) <- glue(" 1")
    } else {
      str_add_list(out) <- stan_prior(
        prior, class = "sderr", px = px, suffix = p,
        comment = "SD of residuals", normalize = normalize
      )
    }
    str_add(out$tpar_def) <- glue(
      "  vector[N{resp}] err{p};  // actual residuals\n"
    )
    str_add(out$pll_args) <- glue(", vector err{p}")
    str_add(out$model_prior) <- glue(
      "  target += std_normal_{lpdf}(zerr{p});\n"
    )
    str_add(out$eta) <- glue(" + err{p}{slice}")
  }

  # validity of the autocor terms has already been checked before
  acframe_arma <- subset2(acframe, class = "arma")
  if (has_rows(acframe_arma)) {
    if (use_threading(threads) && (!acframe_arma$cov || has_natural_residuals)) {
      stop2("Threading is not supported for this ARMA model.")
    }
    str_add(out$data) <- glue(
      "  // data needed for ARMA correlations\n",
      "  int<lower=0> Kar{p};  // AR order\n",
      "  int<lower=0> Kma{p};  // MA order\n"
    )
    str_add(out$tdata_def) <- glue(
      "  int max_lag{p} = max(Kar{p}, Kma{p});\n"
    )
    if (!acframe_arma$cov) {
      err_msg <- "Please set cov = TRUE in ARMA structures"
      if (is.formula(bframe$adforms$se)) {
        stop2(err_msg, " when including known standard errors.")
      }
      str_add(out$data) <- glue(
        "  // number of lags per observation\n",
        "  array[N{resp}] int<lower=0> J_lag{p};\n"
      )
      str_add(out$model_def) <- glue(
        "  // matrix storing lagged residuals\n",
        "  matrix[N{resp}, max_lag{p}] Err{p}",
        " = rep_matrix(0, N{resp}, max_lag{p});\n"
      )
      if (has_natural_residuals) {
        str_add(out$model_def) <- glue(
          "  vector[N{resp}] err{p};  // actual residuals\n"
        )
        Y <- str_if(is.formula(bframe$adforms$mi), "Yl", "Y")
        comp_err <- glue("    err{p}[n] = {Y}{p}[n] - mu{p}[n];\n")
      } else {
        if (acframe_arma$q > 0) {
          # AR and MA structures cannot be distinguished when
          # using a single vector of latent residuals
          stop2("Please set cov = TRUE when modeling MA structures ",
                "for this family.")
        }
        str_add(out$tpar_comp) <- glue(
          "  // compute ctime-series residuals\n",
          "  err{p} = sderr{p} * zerr{p};\n"
        )
        comp_err <- ""
      }
      add_ar <- str_if(acframe_arma$p > 0,
        glue("    mu{p}[n] += Err{p}[n, 1:Kar{p}] * ar{p};\n")
      )
      add_ma <- str_if(acframe_arma$q > 0,
        glue("    mu{p}[n] += Err{p}[n, 1:Kma{p}] * ma{p};\n")
      )
      str_add(out$model_comp_arma) <- glue(
        "  // include ARMA terms\n",
        "  for (n in 1:N{resp}) {{\n",
        add_ma,
        comp_err,
        "    for (i in 1:J_lag{p}[n]) {{\n",
        "      Err{p}[n + 1, i] = err{p}[n + 1 - i];\n",
        "    }}\n",
        add_ar,
        "  }}\n"
      )
    }
    if (acframe_arma$p > 0) {
      if (has_special_prior(prior, px, class = "ar")) {
        if (acframe_arma$cov) {
          stop2("Cannot use shrinkage priors on 'ar' if cov = TRUE.")
        }
        str_add_list(out) <- stan_prior_non_centered(
          class = "ar", suffix = p, suffix_K = "ar"
        )
      } else {
        str_add_list(out) <- stan_prior(
          prior, class = "ar", px = px, suffix = p,
          coef = seq_along(acframe_arma$p),
          type = glue("vector[Kar{p}]"),
          header_type = "vector",
          comment = "autoregressive coefficients",
          normalize = normalize
        )
      }
    }
    if (acframe_arma$q > 0) {
      if (has_special_prior(prior, px, class = "ma")) {
        if (acframe_arma$cov) {
          stop2("Cannot use shrinkage priors on 'ma' if cov = TRUE.")
        }
        str_add_list(out) <- stan_prior_non_centered(
          class = "ma", suffix = p, suffix_K = "ma"
        )
      } else {
        str_add_list(out) <- stan_prior(
          prior, class = "ma", px = px, suffix = p,
          coef = seq_along(acframe_arma$q),
          type = glue("vector[Kma{p}]"),
          header_type = "vector",
          comment = "moving-average coefficients",
          normalize = normalize
        )
      }
    }
  }

  acframe_cosy <- subset2(acframe, class = "cosy")
  if (has_rows(acframe_cosy)) {
    # compound symmetry correlation structure
    # most code is shared with ARMA covariance models
    str_add_list(out) <- stan_prior(
      prior, class = "cosy", px = px, suffix = p,
      comment = "compound symmetry correlation",
      normalize = normalize
    )
  }

  acframe_unstr <- subset2(acframe, class = "unstr")
  if (has_rows(acframe_unstr)) {
    # unstructured correlation matrix
    # most code is shared with ARMA covariance models
    # define prior on the Cholesky scale to consistency across
    # autocorrelation structures
    str_add_list(out) <- stan_prior(
      prior, class = "Lcortime", px = px, suffix = p,
      type = glue("cholesky_factor_corr[n_unique_t{p}]"),
      header_type = "matrix",
      comment = "cholesky factor of unstructured autocorrelation matrix",
      normalize = normalize
    )
  }

  acframe_time_cov <- subset2(acframe, dim = "time", cov = TRUE)
  if (has_rows(acframe_time_cov)) {
    # use correlation structures in covariance matrix parameterization
    # optional for ARMA models and obligatory for COSY and UNSTR models
    # can only model one covariance structure at a time
    stopifnot(nrow(acframe_time_cov) == 1)
    if (use_threading(threads)) {
      stop2("Threading is not supported for covariance-based autocorrelation models.")
    }
    str_add(out$fun) <- glue(
      "  #include 'fun_sequence.stan'\n",
      "  #include 'fun_is_equal.stan'\n",
      "  #include 'fun_stack_vectors.stan'\n"
    )
    if ("gaussian" %in% families) {
      str_add(out$fun) <- glue(
        "  #include 'fun_normal_time.stan'\n",
        "  #include 'fun_normal_time_se.stan'\n"
      )
    }
    if ("student" %in% families) {
      str_add(out$fun) <- glue(
        "  #include 'fun_student_t_time.stan'\n",
        "  #include 'fun_student_t_time_se.stan'\n"
      )
    }
    str_add(out$data) <- glue(
      "  // see the functions block for details\n",
      "  int<lower=1> N_tg{p};\n",
      "  array[N_tg{p}] int<lower=1> begin_tg{p};\n",
      "  array[N_tg{p}] int<lower=1> end_tg{p};\n",
      "  array[N_tg{p}] int<lower=1> nobs_tg{p};\n"
    )
    str_add(out$pll_args) <- glue(
      ", array[] int begin_tg{p}, array[] int end_tg{p}, array[] int nobs_tg{p}"
    )
    str_add(out$tdata_def) <- glue(
      "  int max_nobs_tg{p} = max(nobs_tg{p});",
      "  // maximum dimension of the autocorrelation matrix\n"
    )
    if (acframe_time_cov$class == "unstr") {
      # unstructured time-covariances require additional data and cannot
      # be represented directly via Cholesky factors due to potentially
      # different time subsets
      str_add(out$data) <- glue(
        "  array[N_tg{p}, max(nobs_tg{p})] int<lower=0> Jtime_tg{p};\n",
        "  int n_unique_t{p};  // total number of unique time points\n",
        "  int n_unique_cortime{p};  // number of unique correlations\n"
      )
      str_add(out$pll_args) <- glue(", array[,] int Jtime_tg{p}")
      if (has_latent_residuals) {
        str_add(out$fun) <- "  #include 'fun_scale_time_err_flex.stan'\n"
        str_add(out$tpar_comp) <- glue(
          "  // compute correlated time-series residuals\n",
          "  err{p} = scale_time_err_flex(",
          "zerr{p}, sderr{p}, Lcortime{p}, nobs_tg{p}, begin_tg{p}, end_tg{p}, Jtime_tg{p});\n"
        )
      }
      str_add(out$gen_def) <- glue(
        "  // compute group-level correlations\n",
        "  corr_matrix[n_unique_t{p}] Cortime{p}",
        " = multiply_lower_tri_self_transpose(Lcortime{p});\n",
        "  vector<lower=-1,upper=1>[n_unique_cortime{p}] cortime{p};\n"
      )
      str_add(out$gen_comp) <- stan_cor_gen_comp(
        glue("cortime{p}"), glue("n_unique_t{p}")
      )
    } else {
      # all other time-covariance structures can be represented directly
      # through Cholesky factors of the correlation matrix
      if (acframe_time_cov$class == "arma") {
        if (acframe_time_cov$p > 0 && acframe_time_cov$q == 0) {
          cor_fun <- "ar1"
          cor_args <- glue("ar{p}[1]")
        } else if (acframe_time_cov$p == 0 && acframe_time_cov$q > 0) {
          cor_fun <- "ma1"
          cor_args <- glue("ma{p}[1]")
        } else {
          cor_fun <- "arma1"
          cor_args <- glue("ar{p}[1], ma{p}[1]")
        }
      } else if (acframe_time_cov$class == "cosy") {
        cor_fun <- "cosy"
        cor_args <- glue("cosy{p}")
      }
      str_add(out$fun) <- glue(
        "  #include 'fun_cholesky_cor_{cor_fun}.stan'\n"
      )
      str_add(out$tpar_def) <- glue(
        "  // cholesky factor of the autocorrelation matrix\n",
        "  matrix[max_nobs_tg{p}, max_nobs_tg{p}] Lcortime{p};\n"
      )
      str_add(out$pll_args) <- glue(", matrix Lcortime{p}")
      str_add(out$tpar_comp) <- glue(
        "  // compute residual covariance matrix\n",
        "  Lcortime{p} = cholesky_cor_{cor_fun}({cor_args}, max_nobs_tg{p});\n"
      )
      if (has_latent_residuals) {
        str_add(out$fun) <- "  #include 'fun_scale_time_err.stan'\n"
        str_add(out$tpar_comp) <- glue(
          "  // compute correlated time-series residuals\n",
          "  err{p} = scale_time_err(",
          "zerr{p}, sderr{p}, Lcortime{p}, nobs_tg{p}, begin_tg{p}, end_tg{p});\n"
        )
      }
    }

  }

  acframe_sar <- subset2(acframe, class = "sar")
  if (has_rows(acframe_sar)) {
    if (!has_natural_residuals) {
      stop2("SAR terms are not implemented for this family.")
    }
    if (use_threading(threads)) {
      stop2("Threading is not supported for SAR models.")
    }
    str_add(out$data) <- glue(
      "  matrix[N{resp}, N{resp}] Msar{p};  // spatial weight matrix\n",
      "  vector[N{resp}] eigenMsar{p};  // eigenvalues of Msar{p}\n"
    )
    str_add(out$tdata_def) <- glue(
      "  // the eigenvalues define the boundaries of the SAR correlation\n",
      "  real min_eigenMsar{p} = min(eigenMsar{p});\n",
      "  real max_eigenMsar{p} = max(eigenMsar{p});\n"
    )
    if (acframe_sar$type == "lag") {
      if ("gaussian" %in% families) {
        str_add(out$fun) <- "  #include 'fun_normal_lagsar.stan'\n"
      }
      if ("student" %in% families) {
        str_add(out$fun) <- "  #include 'fun_student_t_lagsar.stan'\n"
      }
      str_add_list(out) <- stan_prior(
        prior, class = "lagsar", px = px, suffix = p,
        comment = "lag-SAR correlation parameter",
        normalize = normalize
      )
    } else if (acframe_sar$type == "error") {
      if ("gaussian" %in% families) {
        str_add(out$fun) <- "  #include 'fun_normal_errorsar.stan'\n"
      }
      if ("student" %in% families) {
        str_add(out$fun) <- "  #include 'fun_student_t_errorsar.stan'\n"
      }
      str_add_list(out) <- stan_prior(
        prior, class = "errorsar", px = px, suffix = p,
        comment = "error-SAR correlation parameter",
        normalize = normalize
      )
    }
  }

  acframe_car <- subset2(acframe, class = "car")
  if (has_rows(acframe_car)) {
    if (is.btnl(bframe)) {
      stop2("CAR terms are not implemented for non-linear models.")
    }
    str_add(out$data) <- glue(
      "  // data for the CAR structure\n",
      "  int<lower=1> Nloc{p};\n",
      "  array[N{resp}] int<lower=1> Jloc{p};\n",
      "  int<lower=0> Nedges{p};\n",
      "  array[Nedges{p}] int<lower=1> edges1{p};\n",
      "  array[Nedges{p}] int<lower=1> edges2{p};\n"
    )
    if (has_special_prior(prior, px, class = "sdcar")) {
      str_add(out$tpar_def) <- glue(
        "  real<lower=0> sdcar{p};  // SD of the CAR structure\n"
      )
      str_add(out$prior_global_scales) <- glue(" sdcar{p}")
      str_add(out$prior_global_lengths) <- glue(" 1")
    } else {
      str_add_list(out) <- stan_prior(
        prior, class = "sdcar", px = px, suffix = p,
        comment = "SD of the CAR structure", normalize = normalize
      )
    }
    str_add(out$pll_args) <- glue(", vector rcar{p}, data array[] int Jloc{p}")
    str_add(out$loopeta) <- glue(" + rcar{p}[Jloc{p}{n}]")
    if (acframe_car$type %in% c("escar", "esicar")) {
      str_add(out$data) <- glue(
        "  vector[Nloc{p}] Nneigh{p};\n",
        "  vector[Nloc{p}] eigenMcar{p};\n"
      )
    }
    if (acframe_car$type == "escar") {
      str_add(out$fun) <- "  #include 'fun_sparse_car_lpdf.stan'\n"
      str_add(out$par) <- glue(
        "  vector[Nloc{p}] rcar{p};\n"
      )
      str_add_list(out) <- stan_prior(
        prior, class = "car", px = px, suffix = p,
        normalize = normalize
      )
      car_args <- c(
        "car", "sdcar", "Nloc", "Nedges",
        "Nneigh", "eigenMcar", "edges1", "edges2"
      )
      car_args <- paste0(car_args, p, collapse = ", ")
      str_add(out$model_prior) <- glue(
        "  target += sparse_car_lpdf(\n",
        "    rcar{p} | {car_args}\n",
        "  );\n"
      )
    } else if (acframe_car$type == "esicar") {
      str_add(out$fun) <- "  #include 'fun_sparse_icar_lpdf.stan'\n"
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
        "sdcar", "Nloc", "Nedges", "Nneigh",
        "eigenMcar", "edges1", "edges2"
      )
      car_args <- paste0(car_args, p, collapse = ", ")
      str_add(out$model_prior) <- glue(
        "  target += sparse_icar_lpdf(\n",
        "    rcar{p} | {car_args}\n",
        "  );\n"
      )
    } else if (acframe_car$type %in% "icar") {
      # intrinsic car based on the case study of Mitzi Morris
      # http://mc-stan.org/users/documentation/case-studies/icar_stan.html
      str_add(out$par) <- glue(
        "  // parameters for the ICAR structure\n",
        "  vector[Nloc{p}] zcar{p};\n"
      )
      # separate definition from computation to support fixed parameters
      str_add(out$tpar_def) <- glue(
        "  // scaled parameters for the ICAR structure\n",
        "  vector[Nloc{p}] rcar{p};\n"
      )
      str_add(out$tpar_comp) <- glue(
        "  // compute scaled parameters for the ICAR structure\n",
        "  rcar{p} = zcar{p} * sdcar{p};\n"
      )
      str_add(out$model_prior) <- glue(
        "  // improper prior on the spatial CAR component\n",
        "  target += -0.5 * dot_self(zcar{p}[edges1{p}] - zcar{p}[edges2{p}]);\n",
        "  // soft sum-to-zero constraint\n",
        "  target += normal_{lpdf}(sum(zcar{p}) | 0, 0.001 * Nloc{p});\n"
      )
    } else if (acframe_car$type == "bym2") {
      # BYM2 car based on the case study of Mitzi Morris
      # http://mc-stan.org/users/documentation/case-studies/icar_stan.html
      str_add(out$data) <- glue(
        "  // scaling factor of the spatial CAR component\n",
        "  real<lower=0> car_scale{p};\n"
      )
      str_add(out$par) <- glue(
        "  // parameters for the BYM2 structure\n",
        "  vector[Nloc{p}] zcar{p};  // spatial part\n",
        "  vector[Nloc{p}] nszcar{p};  // non-spatial part\n",
        "  // proportion of variance in the spatial part\n"
      )
      str_add_list(out) <- stan_prior(
        prior, class = "rhocar", px = px, suffix = p,
        normalize = normalize
      )
      # separate definition from computation to support fixed parameters
      str_add(out$tpar_def) <- glue(
        "  // scaled parameters for the BYM2 structure\n",
        "  vector[Nloc{p}] rcar{p};\n"
      )
      str_add(out$tpar_comp) <- glue(
        "  // join the spatial and the non-spatial CAR component\n",
        "  rcar{p} = (sqrt(1 - rhocar{p}) * nszcar{p}",
        " + sqrt(rhocar{p} * inv(car_scale{p})) * zcar{p}) * sdcar{p};\n"
      )
      str_add(out$model_prior) <- glue(
        "  // improper prior on the spatial BYM2 component\n",
        "  target += -0.5 * dot_self(zcar{p}[edges1{p}] - zcar{p}[edges2{p}]);\n",
        "  // soft sum-to-zero constraint\n",
        "  target += normal_{lpdf}(sum(zcar{p}) | 0, 0.001 * Nloc{p});\n",
        "  // proper prior on the non-spatial BYM2 component\n",
        "  target += std_normal_{lpdf}(nszcar{p});\n"
      )
    }
  }

  acframe_fcor <- subset2(acframe, class = "fcor")
  if (has_rows(acframe_fcor)) {
    if (!has_natural_residuals) {
      stop2("FCOR terms are not implemented for this family.")
    }
    if (use_threading(threads)) {
      stop2("Threading is not supported for FCOR models.")
    }
    if ("gaussian" %in% families) {
      str_add(out$fun) <- "  #include 'fun_normal_fcor.stan'\n"
    }
    if ("student" %in% families) {
      str_add(out$fun) <- "  #include 'fun_student_t_fcor.stan'\n"
    }
    str_add(out$data) <- glue(
      "  matrix[N{resp}, N{resp}] Mfcor{p};  // known residual covariance matrix\n"
    )
    str_add(out$tdata_def) <- glue(
      "  matrix[N{resp}, N{resp}] Lfcor{p} = cholesky_decompose(Mfcor{p});\n"
    )
  }
  out
}

# stan code for offsets
stan_offset <- function(bframe, threads, ...) {
  stopifnot(is.bframel(bframe))
  out <- list()
  if (is.formula(bframe$offset)) {
    p <- usc(combine_prefix(bframe))
    resp <- usc(bframe$resp)
    slice <- stan_slice(threads)
    # use 'offsets' as 'offset' is reserved in stanc3
    str_add(out$data) <- glue( "  vector[N{resp}] offsets{p};\n")
    str_add(out$pll_args) <- glue(", data vector offsets{p}")
    str_add(out$eta) <- glue(" + offsets{p}{slice}")
  }
  out
}

# Stan code for non-linear predictor terms
# @param nlpars names of the non-linear parameters
stan_nl <- function(bframe, nlpars, threads, ...) {
  stopifnot(is.bframenl(bframe))
  out <- list()
  resp <- usc(bframe$resp)
  par <- combine_prefix(bframe, keep_mu = TRUE, nlp = TRUE)
  # prepare non-linear model
  n <- paste0(str_if(bframe$loop, "[n]"), " ")
  new_nlpars <- glue(" nlp{resp}_{nlpars}{n}")
  # covariates in the non-linear model
  covars <- all.vars(bframe$covars)
  new_covars <- NULL
  if (length(covars)) {
    p <- usc(combine_prefix(bframe))
    new_covars <- rep(NA, length(covars))
    frame <- bframe$frame$cnl
    # data_cnl <- data_cnl(bframe, data)
    if (bframe$loop) {
      slice <- stan_nn(threads)
    } else {
      slice <- stan_slice(threads)
    }
    slice <- paste0(slice, " ")
    str_add(out$data) <- "  // covariates for non-linear functions\n"
    for (i in seq_along(covars)) {
      if (frame$integer[i]) {
        if (frame$matrix[i]) {
          str_add(out$data) <- glue(
            "  array[N{resp}, {frame$dim2[i]}] int C{p}_{i};\n"
          )
          str_add(out$pll_args) <- glue(", data array[,] int C{p}_{i}")
        } else {
          str_add(out$data) <- glue(
            "  array[N{resp}] int C{p}_{i};\n"
          )
          str_add(out$pll_args) <- glue(", data array[] int C{p}_{i}")
        }
      } else {
        if (frame$matrix[i]) {
          str_add(out$data) <- glue(
            "  matrix[N{resp}, {frame$dim2[i]}] C{p}_{i};\n"
          )
          str_add(out$pll_args) <- glue(", data matrix C{p}_{i}")
        } else {
          str_add(out$data) <- glue(
            "  vector[N{resp}] C{p}_{i};\n"
          )
          str_add(out$pll_args) <- glue(", data vector C{p}_{i}")
        }
      }
      new_covars[i] <- glue(" C{p}_{i}{slice}")
    }
  }
  # add white spaces to be able to replace parameters and covariates
  syms <- c(
    "+", "-", "*", "/", "%", "^", ".*", "./", "'", ")", "(",
    ",", "==", "!=", "<=", ">=", "<", ">", "!", "&&", "||"
  )
  regex <- glue("(?<!\\.){escape_all(syms)}(?!=)")
  eta <- rm_wsp(deparse0(bframe$formula[[2]]))
  eta <- wsp(rename(eta, regex, wsp(syms), fixed = FALSE, perl = TRUE))
  vars <- c(wsp(nlpars), wsp(covars), " ( ", " ) ")
  new_vars <- c(new_nlpars, new_covars, "(", ")")
  eta <- trimws(rename(eta, vars, new_vars))
  # possibly transform eta in the transformed params block
  str_add(out$model_def) <- glue(
    "  // initialize non-linear predictor term\n",
    "  vector[N{resp}] {par};\n"
  )
  if (bframe$loop) {
    inv_link <- stan_inv_link(bframe$family$link, transform = bframe$transform)
    str_add(out$model_comp_dpar_link) <- glue(
      "  for (n in 1:N{resp}) {{\n",
      stan_nn_def(threads),
      "    // compute non-linear predictor values\n",
      "    {par}[n] = {inv_link}({eta});\n",
      "  }}\n"
    )
  } else {
    inv_link <- stan_inv_link(bframe$family$link, transform = bframe$transform)
    str_add(out$model_comp_dpar_link) <- glue(
      "  // compute non-linear predictor values\n",
      "  {par} = {inv_link}({eta});\n"
    )
  }
  out
}

# global Stan definitions for noise-free variables
stan_Xme <- function(bframe, prior, threads, normalize) {
  meframe <- bframe$frame$me
  stopifnot(is.meframe(meframe))
  if (!has_rows(meframe)) {
    return(list())
  }
  lpdf <- stan_lpdf_name(normalize)
  out <- list()
  coefs <- rename(paste0("me", meframe$xname))
  str_add(out$data) <- "  // data for noise-free variables\n"
  str_add(out$par) <- "  // parameters for noise free variables\n"
  groups <- unique(meframe$grname)
  for (i in seq_along(groups)) {
    g <- groups[i]
    # K are the global and J the local (within group) indices
    K <- which(meframe$grname %in% g)
    J <- seq_along(K)
    if (nzchar(g)) {
      Nme <- glue("Nme_{i}")
      str_add(out$data) <- glue(
        "  int<lower=0> Nme_{i};  // number of latent values\n",
        "  array[N] int<lower=1> Jme_{i};  // group index per observation\n"
      )
      str_add(out$pll_args) <- glue(", data array[] int Jme_{i}")
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
    str_add_list(out) <- stan_prior(
      prior, "meanme", coef = coefs[K], suffix = usc(i),
      type = glue("vector[Mme_{i}]"), comment = "latent means",
      normalize = normalize
    )
    str_add_list(out) <- stan_prior(
      prior, "sdme", coef = coefs[K], suffix = usc(i),
      type = glue("vector[Mme_{i}]"), comment = "latent SDs",
      normalize = normalize
    )
    str_add(out$model_prior) <- cglue(
      "  target += normal_{lpdf}(Xn_{K} | Xme_{K}, noise_{K});\n"
    )
    if (meframe$cor[K[1]] && length(K) > 1L) {
      str_add(out$data) <- glue(
        "  int<lower=1> NCme_{i};  // number of latent correlations\n"
      )
      str_add(out$par) <- glue(
        "  matrix[Mme_{i}, {Nme}] zme_{i};  // standardized latent values\n"
      )
      str_add_list(out) <- stan_prior(
        prior, "Lme", group = g, suffix = usc(i),
        type = glue("cholesky_factor_corr[Mme_{i}]"),
        comment = "cholesky factor of the latent correlation matrix",
        normalize = normalize
      )
      # separate definition from computation to support fixed parameters
      str_add(out$tpar_def) <- glue(
        "  matrix[{Nme}, Mme_{i}] Xme{i};  // actual latent values\n"
      )
      str_add(out$tpar_comp) <- glue(
        "  // compute actual latent values\n",
        "  Xme{i} = rep_matrix(transpose(meanme_{i}), {Nme})",
        " + transpose(diag_pre_multiply(sdme_{i}, Lme_{i}) * zme_{i});\n"
      )
      str_add(out$tpar_def) <- cglue(
        "  // using separate vectors increases efficiency\n",
        "  vector[{Nme}] Xme_{K};\n"
      )
      str_add(out$tpar_comp) <- cglue(
        "  Xme_{K} = Xme{i}[, {J}];\n"
      )
      str_add(out$pll_args) <- cglue(", vector Xme_{K}")
      str_add(out$model_prior) <- glue(
        "  target += std_normal_{lpdf}(to_vector(zme_{i}));\n"
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
      # separate definition from computation to support fixed parameters
      str_add(out$tpar_def) <- cglue(
        "  vector[{Nme}] Xme_{K};  // actual latent values\n"
      )
      str_add(out$tpar_comp) <- cglue(
        "  // compute actual latent values\n",
        "  Xme_{K} = meanme_{i}[{J}] + sdme_{i}[{J}] * zme_{K};\n"
      )
      str_add(out$pll_args) <- cglue(", vector Xme_{K}")
      str_add(out$model_prior) <- cglue(
        "  target += std_normal_{lpdf}(zme_{K});\n"
      )
    }
  }
  out
}

# initialize and compute a linear predictor term in Stan language
# @param out list of character strings containing Stan code
# @param bframe btl object
# @param primitive use Stan's GLM likelihood primitives?
# @param ... currently unused
# @return list of character strings containing Stan code
stan_eta_combine <- function(bframe, out, threads, primitive, ...) {
  stopifnot(is.btl(bframe), is.list(out))
  if (primitive && !has_special_terms(bframe)) {
    # only overall effects and perhaps an intercept are present
    # which will be evaluated directly in the GLM primitive likelihood
    return(out)
  }
  px <- check_prefix(bframe)
  resp <- usc(bframe$resp)
  eta <- combine_prefix(px, keep_mu = TRUE, nlp = TRUE)
  out$eta <- sub("^[ \t\r\n]+\\+", "", out$eta, perl = TRUE)
  str_add(out$model_def) <- glue(
    "  // initialize linear predictor term\n",
    "  vector[N{resp}] {eta} = rep_vector(0.0, N{resp});\n"
  )
  if (isTRUE(nzchar(out$eta))) {
    str_add(out$model_comp_eta_basic) <- glue("  {eta} +={out$eta};\n")
  }
  out$eta <- NULL
  str_add(out$loopeta) <- stan_eta_re(bframe, threads = threads)
  if (isTRUE(nzchar(out$loopeta))) {
    # parts of eta are computed in a loop over observations
    out$loopeta <- sub("^[ \t\r\n]+\\+", "", out$loopeta, perl = TRUE)
    str_add(out$model_comp_eta_loop) <- glue(
      "  for (n in 1:N{resp}) {{\n",
      "    // add more terms to the linear predictor\n",
      stan_nn_def(threads),
      "    {eta}[n] +={out$loopeta};\n",
      "  }}\n"
    )
  }
  out$loopeta <- NULL
  # some links need custom Stan functions
  link <- bframe$family$link
  link_names <- c("cauchit", "cloglog", "softplus", "squareplus", "softit", "tan_half")
  needs_link_fun <- isTRUE(link %in% link_names)
  if (needs_link_fun) {
    str_add(out$fun) <- glue("  #include 'fun_{link}.stan'\n")
  }
  # possibly transform eta before it is passed to the likelihood
  inv_link <- stan_inv_link(bframe$family$link, transform = bframe$transform)
  if (nzchar(inv_link)) {
    str_add(out$model_comp_dpar_link) <- glue(
      "  {eta} = {inv_link}({eta});\n"
    )
  }
  out
}

# write the group-level part of the linear predictor
# @return a single character string
stan_eta_re <- function(bframe, threads) {
  eta_re <- ""
  n <- stan_nn(threads)
  reframe <- subset2(bframe$frame$re, type = c("", "mmc"))
  for (id in unique(reframe$id)) {
    r <- subset2(reframe, id = id)
    rpx <- check_prefix(r)
    idp <- paste0(r$id, usc(combine_prefix(rpx)))
    idresp <- paste0(r$id, usc(rpx$resp))
    rprefix <- str_if(isTRUE(r$s2z[1]), "r_s2z_", "r_")
    if (r$gtype[1] == "mm") {
      ng <- seq_along(r$gcall[[1]]$groups)
      for (i in seq_rows(r)) {
        str_add(eta_re) <- cglue(
          " + W_{idresp[i]}_{ng}{n}",
          " * {rprefix}{idp[i]}_{r$cn[i]}[J_{idresp[i]}_{ng}{n}]",
          " * Z_{idp[i]}_{r$cn[i]}_{ng}{n}"
        )
      }
    } else {
      str_add(eta_re) <- cglue(
        " + {rprefix}{idp}_{r$cn}[J_{idresp}{n}] * Z_{idp}_{r$cn}{n}"
      )
    }
  }
  eta_re
}

# Stan code for group-level parameters in special predictor terms
# @param r data.frame created by frame_re
# @return a character vector: one element per row of 'r'
stan_eta_rsp <- function(r, threads) {
  stopifnot(is.reframe(r))
  stopifnot(nrow(r) > 0L, length(unique(r$gtype)) == 1L)
  n <- stan_nn(threads)
  rpx <- check_prefix(r)
  idp <- paste0(r$id, usc(combine_prefix(rpx)))
  idresp <- paste0(r$id, usc(rpx$resp))
  if (r$gtype[1] == "mm") {
    ng <- seq_along(r$gcall[[1]]$groups)
    out <- rep("", nrow(r))
    for (i in seq_along(out)) {
      out[i] <- glue(
        "W_{idresp[i]}_{ng}{n} * r_{idp[i]}_{r$cn[i]}[J_{idresp[i]}_{ng}{n}]",
        collapse = " + "
      )
    }
  } else {
    out <- glue("r_{idp}_{r$cn}[J_{idresp}{n}]")
  }
  out
}

# does eta need to be transformed manually using the inv_link function
stan_eta_transform <- function(bframe, family) {
  no_transform <- family$link == "identity" ||
    has_joint_link(family) && !is.customfamily(family)
  !no_transform && !stan_has_built_in_fun(bframe, family)
}

# indicate if the population-level design matrix should be centered
# implies a temporary shift in the intercept of the model
stan_center_X <- function(x) {
  is.btl(x) && !no_center(x$fe) && has_intercept(x$fe) &&
    !fix_intercepts(x) && !is_sparse(x$fe) && !has_sum_to_zero_thres(x)
}

stan_dpar_comments <- function(dpar, family) {
  dpar_class <- dpar_class(dpar, family)
  out <- switch(dpar_class, "",
    sigma = "dispersion parameter",
    shape = "shape parameter",
    nu = "degrees of freedom or shape",
    phi = "precision parameter",
    kappa = "precision parameter",
    beta = "scale parameter",
    zi = "zero-inflation probability",
    hu = "hurdle probability",
    zoi = "zero-one-inflation probability",
    coi = "conditional one-inflation probability",
    bs = "boundary separation parameter",
    ndt = "non-decision time parameter",
    bias = "initial bias parameter",
    disc = "discrimination parameters",
    quantile = "quantile parameter",
    xi = "shape parameter",
    alpha = "skewness parameter"
  )
  out
}

# Stan code for transformations of distributional parameters
# TODO: refactor into family-specific functions
# TODO: add gamma and related families here to compute rate = shape / mean
stan_dpar_transform <- function(bframe, prior, threads, normalize, ...) {
  stopifnot(is.brmsterms(bframe))
  out <- list()
  families <- family_names(bframe)
  px <- check_prefix(bframe)
  p <- usc(combine_prefix(px))
  resp <- usc(bframe$resp)
  if (any(conv_cats_dpars(families))) {
    stopifnot(length(families) == 1L)
    iref <- get_refcat(bframe$family, int = TRUE)
    mus <- make_stan_names(glue("mu{bframe$family$cats}"))
    mus <- glue("{mus}{p}")
    if (use_glm_primitive_categorical(bframe)) {
      bterms1 <- bframe$dpars[[1]]
      center <- stan_center_X(bterms1)
      ct <- str_if(center, "c")
      K <- glue("K{ct}_{bterms1$dpar}{p}")
      str_add(out$pll_args) <- glue(", int {K}")
      str_add(out$model_def) <- glue(
        "  // joint regression coefficients over categories\n",
        "  matrix[{K}, ncat{p}] b{p};\n"
      )
      bnames <- glue("b_{mus}")
      bnames[iref] <- glue("rep_vector(0, {K})")
      str_add(out$model_comp_catjoin) <- cglue(
        "  b{p}[, {seq_along(bnames)}] = {bnames};\n"
      )
      if (center) {
        Inames <- glue("Intercept_{mus}")
        Inames[iref] <- "0"
        str_add(out$model_def) <- glue(
          "  // joint intercepts over categories\n",
          "  vector[ncat{p}] Intercept{p};\n"
        )
        str_add(out$model_comp_catjoin) <- glue(
          "  Intercept{p} = {stan_vector(Inames)};\n"
        )
      }
    } else {
      is_logistic_normal <- any(is_logistic_normal(families))
      len_mu <- glue("ncat{p}{str_if(is_logistic_normal, '-1')}")
      str_add(out$model_def) <- glue(
        "  // linear predictor matrix\n",
        "  array[N{resp}] vector[{len_mu}] mu{p};\n"
      )
      mus <- glue("{mus}[n]")
      if (is_logistic_normal) {
        mus <- mus[-iref]
      } else {
        mus[iref] <- "0"
      }
      str_add(out$model_comp_catjoin) <- glue(
        "  for (n in 1:N{resp}) {{\n",
        "    mu{p}[n] = {stan_vector(mus)};\n",
        "  }}\n"
      )
    }
  }
  if (any(families %in% "skew_normal")) {
    # as suggested by Stephen Martin use sigma and mu of CP
    # but the skewness parameter alpha of DP
    dp_names <- names(bframe$dpars)
    for (i in which(families %in% "skew_normal")) {
      id <- str_if(length(families) == 1L, "", i)
      sigma <- stan_sigma_transform(bframe, id = id, threads = threads)
      ns <- str_if(grepl(stan_nn_regex(), sigma), "[n]")
      na <- str_if(glue("alpha{id}") %in% dp_names, "[n]")
      type_delta <- str_if(nzchar(na), glue("vector[N{resp}]"), "real")
      no <- str_if(any(nzchar(c(ns, na))), "[n]", "")
      type_omega <- str_if(nzchar(no), glue("vector[N{resp}]"), "real")
      str_add(out$model_def) <- glue(
        "  // parameters used to transform the skew-normal distribution\n",
        "  {type_delta} delta{id}{p};  // transformed alpha parameter\n",
        "  {type_omega} omega{id}{p};  // scale parameter\n"
      )
      alpha <- glue("alpha{id}{p}{na}")
      delta <- glue("delta{id}{p}{na}")
      omega <- glue("omega{id}{p}{no}")
      comp_delta <- glue(
        "  {delta} = {alpha} / sqrt(1 + {alpha}^2);\n"
      )
      comp_omega <- glue(
        "  {omega} = {sigma} / sqrt(1 - sqrt(2 / pi())^2 * {delta}^2);\n"
      )
      str_add(out$model_comp_dpar_trans) <- glue(
        "  // use efficient skew-normal parameterization\n",
        str_if(!nzchar(na), comp_delta),
        str_if(!nzchar(no), comp_omega),
        "  for (n in 1:N{resp}) {{\n",
        stan_nn_def(threads),
        str_if(nzchar(na), glue("  ", comp_delta)),
        str_if(nzchar(no), glue("  ", comp_omega)),
        "    mu{id}{p}[n] = mu{id}{p}[n]",
        " - {omega} * {delta} * sqrt(2 / pi());\n",
        "  }}\n"
      )
    }
  }
  if (any(families %in% "gen_extreme_value")) {
    # TODO: remove the gen_extreme_value family in brms 3.0
    warning2("The 'gen_extreme_value' family is deprecated ",
             "and will be removed in the future.")
    dp_names <- c(names(bframe$dpars), names(bframe$fdpars))
    for (i in which(families %in% "gen_extreme_value")) {
      id <- str_if(length(families) == 1L, "", i)
      xi <- glue("xi{id}")
      if (!xi %in% dp_names) {
        str_add(out$model_def) <- glue(
          "  real {xi}{p};  // scaled shape parameter\n"
        )
        args <- sargs(
          glue("tmp_{xi}{p}"), glue("Y{p}"),
          glue("mu{id}{p}"), glue("sigma{id}{p}")
        )
        str_add(out$model_comp_dpar_trans) <- glue(
          "  {xi}{p} = scale_xi({args});\n"
        )
      }
    }
  }
  if (any(families %in% "logistic_normal")) {
    stopifnot(length(families) == 1L)
    predcats <- make_stan_names(get_predcats(bframe$family))
    sigma_dpars <- glue("sigma{predcats}")
    reqn <- sigma_dpars %in% names(bframe$dpars)
    n <- ifelse(reqn, "[n]", "")
    sigma_dpars <- glue("{sigma_dpars}{p}{n}")
    ncatm1 <- glue("ncat{p}-1")
    if (any(reqn)) {
      # some of the sigmas are predicted
      str_add(out$model_def) <- glue(
        "  // sigma parameter matrix\n",
        "  array[N{resp}] vector[{ncatm1}] sigma{p};\n"
      )
      str_add(out$model_comp_catjoin) <- glue(
        "  for (n in 1:N{resp}) {{\n",
        "    sigma{p}[n] = {stan_vector(sigma_dpars)};\n",
        "  }}\n"
      )
    } else {
      # none of the sigmas is predicted
      str_add(out$model_def) <- glue(
        "  // sigma parameter vector\n",
        "  vector[{ncatm1}] sigma{p} = {stan_vector(sigma_dpars)};\n"
      )
    }
    # handle the latent correlation matrix 'lncor'
    str_add(out$tdata_def) <- glue(
      "  // number of logistic normal correlations\n",
      "  int nlncor{p} = choose({ncatm1}, 2);\n"
    )
    str_add_list(out) <- stan_prior(
      prior, "Llncor", suffix = p, px = px,
      type = glue("cholesky_factor_corr[{ncatm1}]"),
      header_type = "matrix",
      comment = "logistic-normal Cholesky correlation matrix",
      normalize = normalize
    )
    str_add(out$gen_def) <- glue(
      "  // logistic normal correlations\n",
      "  corr_matrix[{ncatm1}] Lncor",
      " = multiply_lower_tri_self_transpose(Llncor);\n",
      "  vector<lower=-1,upper=1>[nlncor] lncor;\n"
    )
    str_add(out$gen_comp) <- stan_cor_gen_comp("lncor", ncatm1)
  }
  out
}

# Stan code for sigma to incorporate addition argument 'se'
stan_sigma_transform <- function(bframe, id = "", threads = NULL) {
  if (nzchar(id)) {
    # find the right family in mixture models
    family <- family_names(bframe)[as.integer(id)]
  } else {
    family <- bframe$family$family
    stopifnot(!isTRUE(family == "mixture"))
  }
  p <- usc(combine_prefix(bframe))
  ns <- str_if(glue("sigma{id}") %in% names(bframe$dpars), "[n]")
  has_sigma <- has_sigma(family) && !no_sigma(bframe)
  sigma <- str_if(has_sigma, glue("sigma{id}{p}{ns}"))
  if (is.formula(bframe$adforms$se)) {
    nse <- stan_nn(threads)
    sigma <- str_if(nzchar(sigma),
      glue("sqrt(square({sigma}) + se2{p}{nse})"),
      glue("se{p}{nse}")
    )
  }
  sigma
}
