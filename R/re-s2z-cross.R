# Helpers for conventional S2Z IDs spanning response predictors

# Does one conventional S2Z covariance ID occur in more than one linear
# predictor? The ID remains one ordinary brms covariance block; only its
# sampling coordinates are assembled globally.
is_re_s2z_cross_id <- function(bframe, id) {
  stopifnot(is.anybrmsframe(bframe), length(id) == 1L)
  r <- subset2(bframe$frame$re, id = id)
  if (!has_rows(r) || !all(r$s2z) || any(re_s2z_latent(r))) {
    return(FALSE)
  }
  length(unique(combine_prefix(check_prefix(r)))) > 1L
}

# Data and Stan generation must agree on whether a cross-predictor block can
# use the restricted-covariance eigenmode chart. The remaining Fisher checks
# (Gaussian identity responses and scalar residual scales) issue their usual
# contextual diagnostics while building the response information.
is_re_s2z_cross_spectral_candidate <- function(bframe, r, id = unique(r$id)) {
  stopifnot(is.anybrmsframe(bframe), is.reframe(r), length(id) == 1L)
  is.mvbrmsframe(bframe) && isTRUE(bframe$rescor) &&
    is_re_s2z_cross_id(bframe, id) &&
    nrow(r) == length(unique(r$resp)) && all(r$coef == "Intercept")
}

# Return the response-local frames touched by one cross-predictor ID, ordered
# by the first occurrence of their prefix in the global covariance block.
re_s2z_cross_frames <- function(bframe, r) {
  stopifnot(is.anybrmsframe(bframe), is.reframe(r), has_rows(r))
  keys <- unique(combine_prefix(check_prefix(r)))
  frames <- all_bframel(bframe)
  out <- lapply(keys, function(key) {
    found <- Filter(
      function(x) identical(combine_prefix(check_prefix(x)), key),
      frames
    )
    if (length(found) != 1L) {
      stop2("Internal error while matching a cross-predictor sum-to-zero ",
            "block to its linear predictors.")
    }
    found[[1L]]
  })
  names(out) <- keys
  out
}

# Describe the global omitted-mean map for one conventional covariance block.
# Population coordinates remain predictor-local, but are stacked internally
# so the correlated conventional group mean can be integrated in one system.
re_s2z_cross_info <- function(bframe, prior, id) {
  stopifnot(
    is.anybrmsframe(bframe), is.brmsprior(prior), length(id) == 1L
  )
  r <- subset2(bframe$frame$re, id = id)
  if (!is_re_s2z_cross_id(bframe, id)) {
    stop2("Internal error: expected a cross-predictor sum-to-zero ID.")
  }
  frames <- re_s2z_cross_frames(bframe, r)
  infos <- lapply(frames, re_s2z_info, prior = prior, id = id)
  q <- vapply(infos, function(x) length(x$qnames), integer(1))
  starts <- cumsum(c(1L, head(q, -1L)))
  ends <- starts + q - 1L
  nlist(
    id, r, frames, infos, q, starts, ends,
    Q = sum(q), M = nrow(r)
  )
}

# Assemble expected Fisher information for one conventional S2Z covariance ID
# spanning response predictors. Conditionally independent responses contribute
# response-local outer products. For an ordinary multivariate Gaussian model
# with residual correlation, the response designs are whitened jointly by the
# current residual Cholesky factor, retaining its off-diagonal information.
stan_re_s2z_cross_fisher_info <- function(id, r, bframe, threads) {
  stopifnot(
    is.mvbrmsframe(bframe), is.reframe(r), has_rows(r), all(r$s2z),
    is_re_s2z_cross_id(bframe, id)
  )
  unsupported <- function(detail) {
    stop2("Fisher centering for a cross-predictor S2Z ID ", detail, ".")
  }
  responses <- unique(r$resp)
  local <- lapply(responses, function(response) {
    term <- bframe$terms[[response]]
    if (!is.brmsframe(term)) {
      unsupported("could not match a response-local model frame")
    }
    r_source <- r[r$resp == response, , drop = FALSE]
    info <- stan_re_s2z_fisher_info(
      id, r = r_source, bframe = term, threads = threads
    )
    nlist(response, term, r_source, info)
  })
  M <- nrow(r)

  if (!isTRUE(bframe$rescor)) {
    sources <- unlist(lapply(local, function(x) {
      info <- x$info
      if (isTRUE(info$fixed_design)) {
        return(list(nlist(
          columns = as.integer(x$r_source$cn), N = info$N,
          group = info$group, design_at_n = info$design_at_n,
          obs_prec = info$obs_prec, definitions_at_n = ""
        )))
      }
      lapply(info$sources, function(source) {
        source$columns <- as.integer(x$r_source$cn[source$columns])
        source
      })
    }), recursive = FALSE)
    return(nlist(
      M, fixed_design = FALSE, sources,
      fun = paste0(vapply(
        local, function(x) x$info$fun %||% "", character(1)
      ), collapse = "")
    ))
  }

  # Residual correlation is currently modeled only for multivariate Gaussian
  # or Student responses. Restrict Fisher centering to the Gaussian identity-
  # link case with observation-invariant residual scales; this admits an exact
  # joint location-information matrix in transformed parameters.
  response_order <- bframe$responses
  terms <- bframe$terms[response_order]
  families <- vapply(terms, function(x) x$family$family, character(1))
  links <- vapply(terms, function(x) x$family$link, character(1))
  if (any(families != "gaussian") || any(links != "identity")) {
    unsupported(paste0(
      "with residual response correlation currently requires Gaussian ",
      "identity-link responses"
    ))
  }
  if (any(!vapply(local, function(x) isTRUE(x$info$fixed_design), logical(1)))) {
    unsupported(paste0(
      "with residual response correlation requires observation-invariant ",
      "residual scales"
    ))
  }
  sigma <- vapply(terms, function(term) {
    value <- stan_sigma_transform(term, threads = threads)
    if (length(value) != 1L || grepl(stan_nn_regex(), value)) {
      unsupported(paste0(
        "with residual response correlation requires one scalar residual ",
        "scale per response"
      ))
    }
    as.character(value)
  }, character(1))
  first <- local[[1L]]$info
  assignments <- paste0(vapply(local, function(x) {
    response_index <- match(x$response, response_order)
    cglue(
      "        design_fisher_s2z[{response_index}, ",
      "{as.integer(x$r_source$cn)}] = {x$info$design_at_n};\n"
    )
  }, character(1)), collapse = "")
  # A residual-correlated intercept-only block has one observation-invariant
  # embedded design matrix.  Its residual Fisher matrix can be formed once and
  # then carried into the restricted grouping-covariance eigenmodes.  Unequal
  # group counts remain valid: their diagonal modal information is represented
  # by a transformed-data exposure.  With equal counts this is the exact modal
  # Fisher decomposition.
  spectral <- NULL
  spectral_intercept <-
    is_re_s2z_cross_spectral_candidate(bframe, r, id) &&
    length(local) == M && all(vapply(
    local,
    function(x) {
      nrow(x$r_source) == 1L &&
        identical(x$r_source$coef, "Intercept") &&
        isTRUE(x$info$fixed_design)
    },
    logical(1)
  ))
  if (spectral_intercept) {
    spectral_assignments <- paste0(vapply(local, function(x) {
      response_index <- match(x$response, response_order)
      cglue(
        "      design_modal_fisher_s2z[{response_index}, ",
        "{as.integer(x$r_source$cn)}] = 1.0;\n"
      )
    }, character(1)), collapse = "")
    spectral <- nlist(
      N = first$N, group = first$group,
      assignments = spectral_assignments, sigma
    )
  }
  joint_info_comp <- glue(
    "    {{\n",
    "      vector[nresp] sigma_fisher_s2z = {stan_vector(sigma)};\n",
    "      matrix[nresp, nresp] L_residual_fisher_s2z = ",
    "diag_pre_multiply(sigma_fisher_s2z, Lrescor);\n",
    "      for (n in 1:{first$N}) {{\n",
    "        matrix[nresp, M_{id}] design_fisher_s2z = ",
    "rep_matrix(0.0, nresp, M_{id});\n",
    "        matrix[nresp, M_{id}] white_design_fisher_s2z;\n",
    "{assignments}",
    "        white_design_fisher_s2z = mdivide_left_tri_low(\n",
    "          L_residual_fisher_s2z, design_fisher_s2z\n",
    "        );\n",
    "        info_fisher_s2z[{first$group}[n]] += ",
    "crossprod(white_design_fisher_s2z);\n",
    "      }}\n",
    "    }}\n"
  )
  nlist(
    M, fixed_design = FALSE, sources = list(), joint_info_comp,
    spectral, fun = ""
  )
}

# Validate one conventional covariance block spanning response predictors.
# The coefficient covariance, group distribution, scales, and known grouping
# covariance retain their ordinary meanings; only the sampling chart changes.
validate_re_s2z_cross_id <- function(bframe, prior, id) {
  stopifnot(
    is.anybrmsframe(bframe), is.brmsprior(prior), length(id) == 1L
  )
  r <- subset2(bframe$frame$re, id = id)
  stopifnot(is.reframe(r), has_rows(r))
  if (!is.mvbrmsframe(bframe) || any(!nzchar(r$resp)) ||
      any(nzchar(r$dpar)) || any(nzchar(r$nlpar))) {
    stop2("A cross-predictor sum-to-zero ID currently supports only ",
          "multivariate response-location predictors.")
  }
  if (length(unique(r$group)) != 1L) {
    stop2("All coefficients sharing a cross-predictor sum-to-zero ID must ",
          "use the same grouping factor.")
  }
  if (length(unique(r$cov)) != 1L) {
    stop2("All coefficients sharing a cross-predictor sum-to-zero ID must ",
          "use the same 'cov' setting.")
  }
  if (length(unique(r$scale)) != 1L) {
    stop2("All coefficients sharing a cross-predictor sum-to-zero ID must ",
          "use the same 'scale' setting.")
  }
  if (length(unique(r$dist)) != 1L) {
    stop2("All coefficients sharing a cross-predictor sum-to-zero ID must ",
          "use the same group-level distribution.")
  }
  if (!all(r$cor)) {
    stop2("A cross-predictor sum-to-zero ID currently requires correlated ",
          "group effects.")
  }
  frames <- re_s2z_cross_frames(bframe, r)
  for (x in frames) {
    rx <- x$frame$re
    other <- unique(rx$id[
      rx$s2z & !re_s2z_latent(rx) & rx$id != id
    ])
    if (length(other)) {
      stop2("A predictor participating in a cross-predictor sum-to-zero ID ",
            "cannot yet contain another conventional sum-to-zero block.")
    }
  }
  # Constructing the descriptors here also validates every response-local
  # population-column match and effective population prior.
  re_s2z_cross_info(bframe, prior = prior, id = id)
  invisible(NULL)
}

# Stan code for a conventional correlated block spanning response predictors.
# Let delta be the componentwise zero-sum group deviations, m the omitted
# conventional group mean, beta the conventional population coefficients, and
# theta = beta + H m the finite-population coordinates used by the likelihood.
# The code integrates m exactly and reconstructs beta and b = delta + m in
# generated quantities while retaining sd and L on their conventional scale.
.stan_re_s2z_cross <- function(id, bframe, prior, normalize,
                               out = list(), fisher_info = NULL, ...) {
  if (is.null(out[["tpar_prior"]])) {
    out[["tpar_prior"]] <- ""
  }
  lpdf <- stan_lpdf_name(normalize)
  cross <- re_s2z_cross_info(bframe, prior = prior, id = id)
  r <- cross$r
  Q <- cross$Q
  M <- cross$M
  J <- seq_rows(r)
  px <- check_prefix(r)
  idp <- paste0(r$id, usc(combine_prefix(px)))
  has_cov <- nzchar(r$cov[1])
  varying <- identical(r$scale[1], "varying")
  is_student <- identical(r$dist[1], "student")
  s2z_mode <- re_s2z_center_mode(r)
  s2z_center <- identical(s2z_mode, "centered")
  s2z_fisher <- !is.null(fisher_info)
  s2z_partial <- s2z_mode %in% c("partial", "auto")
  spectral_fisher <- s2z_fisher && has_cov && !varying && !is_student &&
    !is.null(fisher_info$spectral)
  scale <- if (varying) {
    glue("reference_sd_s2z_{id}")
  } else {
    glue("sd_{id}")
  }

  stopifnot(
    M > 1L, all(r$s2z), !any(re_s2z_latent(r)), all(r$cor),
    s2z_mode %in% c("centered", "noncentered", "partial", "auto"),
    !s2z_fisher || identical(s2z_mode, "auto")
  )

  if (s2z_fisher) {
    out <- stan_re_s2z_fisher_def(out, id)
    if (spectral_fisher) {
      str_add(out$data) <- glue(
        "  matrix[N_{id}, N_{id} - 1] Ecov_s2z_{id};",
        "  // eigenvectors of the covariance restricted to the S2Z space\n",
        "  vector<lower=0>[N_{id} - 1] lambda_cov_s2z_{id};",
        "  // restricted covariance eigenvalues\n",
        "  vector[N_{id} - 1] coupling_cov_s2z_{id};",
        "  // mean-contrast coupling E' Omega^-1 1\n",
        "  real<lower=0> mean_prec_cov_s2z_{id};",
        "  // precision 1' Omega^-1 1 of the omitted arithmetic mean\n"
      )
      out <- stan_re_s2z_modal_fisher_tdata(
        out, id, spectral = fisher_info$spectral
      )
    } else if (has_cov) {
      out <- stan_re_s2z_cov_fisher_tdata(out, id)
    }
  } else if (s2z_partial) {
    str_add(out$data) <- glue(
      "  matrix<lower=0,upper=1>[N_{id}, M_{id}] rho_s2z_{id};",
      "  // fixed numeric centering fractions\n"
    )
    out <- stan_re_s2z_partial_tdata(out, id)
  }

  str_add(out$data) <- glue(
    "  int<lower=1> NC_{id};  // number of group-level correlations\n"
  )
  str_add_list(out) <- stan_prior(
    prior, class = "L", group = r$group[1], suffix = usc(id),
    type = glue("cholesky_factor_corr[M_{id}]"),
    comment = "cholesky factor of correlation matrix",
    normalize = normalize
  )
  if (varying) {
    str_add_list(out) <- stan_prior(
      prior, class = "sdlog", group = r$group[1], coef = r$coef,
      type = glue("vector[M_{id}]"), suffix = glue("_{id}"), px = px,
      comment = "SDs of group-level log-scale deviations",
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
    },
    str_if(
      varying,
      glue(
        "  vector[M_{id} * N_{id}] z_sd_s2z_{id};",
        "  // flattened orthonormal log-scale coordinates: ",
        "contrasts, then means\n"
      )
    )
  )
  if (varying) {
    str_add(out$model_prior) <- glue(
      "  target += std_normal_{lpdf}(z_sd_s2z_{id});\n"
    )
  }

  # Conditional Gaussian mixtures preserve exact Student-t and Cauchy priors
  # on the conventional population coefficients. Names remain predictor-local
  # so existing exclusion and public-draw machinery continues to apply.
  for (a in seq_along(cross$infos)) {
    info <- cross$infos[[a]]
    p <- info$p
    for (k in seq_along(info$prior)) {
      spec <- info$prior[[k]]
      if (spec$dist == "student") {
        str_add(out$par) <- glue(
          "  real<lower=0> udf_b_s2z{p}_{k};",
          "  // mixing variable for population coefficient {k}\n"
        )
        str_add(out$tpar_prior) <- glue(
          "  lprior += inv_chi_square_{lpdf}(udf_b_s2z{p}_{k} | ",
          "{stan_s2z_number(spec$df)});\n"
        )
      }
    }
  }

  str_add(out$tpar_def) <- glue(
    "  // conventional cross-response block in S2Z sampling coordinates\n",
    "  matrix[N_{id}, M_{id}] r_s2z_{id};\n",
    "  vector[{Q}] prior_mean_s2z_{id};\n",
    "  vector<lower=0>[{Q}] prior_prec_s2z_{id};\n",
    "  matrix[M_{id}, M_{id}] L_Sigma_s2z_{id};\n",
    "  matrix[M_{id}, M_{id}] P_s2z_{id};\n",
    "  matrix[M_{id}, M_{id}] L_P_s2z_{id};\n",
    str_if(
      varying,
      glue(
        "  matrix<lower=0>[N_{id}, M_{id}] sd_level_s2z_{id};\n",
        "  vector<lower=0>[M_{id}] reference_sd_s2z_{id};\n"
      )
    ),
    str_if(
      is_student,
      glue(
        "  vector<lower=0>[N_{id}] group_scale_s2z_{id};\n",
        "  vector<lower=0>[N_{id}] group_prec_s2z_{id};\n"
      )
    ),
    "  matrix[M_{id}, M_{id}] P_group_s2z_{id};\n",
    "  vector[M_{id}] h_group_s2z_{id};\n",
    "  real group_quad_s2z_{id};\n",
    str_if(s2z_partial, glue("  real log_det_partial_s2z_{id};\n")),
    "  // zero-sum deviations used by the observation-level likelihood\n",
    cglue("  vector[N_{id}] r_s2z_{idp}_{r$cn};\n")
  )
  str_add(out$tdata_def) <- glue(
    "  matrix[{Q}, M_{id}] H_s2z_{id};\n"
  )
  str_add(out$tdata_comp) <- glue(
    "  H_s2z_{id} = rep_matrix(0.0, {Q}, M_{id});\n"
  )
  for (a in seq_along(cross$infos)) {
    info <- cross$infos[[a]]
    offset <- cross$starts[a] - 1L
    for (j in seq_rows(info$r)) {
      qi <- offset + info$match_q[j]
      column <- info$r$cn[j]
      str_add(out$tdata_comp) <- glue(
        "  H_s2z_{id}[{qi}, {column}] = 1.0;\n"
      )
      if (info$center && info$r$coef[j] != "Intercept") {
        str_add(out$tdata_comp) <- glue(
          "  H_s2z_{id}[{offset + 1L}, {column}] = ",
          "means_X{info$p}[{info$match_q[j] - 1L}];\n"
        )
      }
    }
  }

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
    str_add(out$tpar_prior) <- .stan_re_s2z_sd_level_prior(
      id, r = r, prior = prior, normalize = normalize
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
  partial_transform <- if (s2z_partial && !spectral_fisher) {
    stan_re_s2z_partial_cor_transform(id)
  } else {
    str_if(
      !s2z_center,
      glue("  r_s2z_{id} = r_s2z_{id} * L_Sigma_s2z_{id}';\n")
    )
  }
  str_add(out$tpar_comp) <- glue(
    "  L_Sigma_s2z_{id} = diag_pre_multiply({scale}, L_{id});\n"
  )
  if (s2z_fisher) {
    if (spectral_fisher) {
      str_add(out$tpar_comp) <- stan_re_s2z_modal_fisher_comp(
        id, spectral = fisher_info$spectral,
        L = glue("L_Sigma_s2z_{id}"), M = M
      )
    } else {
      str_add(out$tpar_comp) <- stan_re_s2z_fisher_comp(
        id, r = r, fisher_info = fisher_info,
        L = glue("L_Sigma_s2z_{id}"),
        row_var = if (has_cov) glue("row_var_fisher_s2z_{id}[j]") else NULL
      )
    }
  }
  if (spectral_fisher) {
    str_add(out$tpar_comp) <- stan_re_s2z_modal_transform_group_comp(
      id, L = glue("L_Sigma_s2z_{id}"), M = M
    )
  } else {
    str_add(out$tpar_comp) <- glue(
      "  for (k in 1:M_{id}) {{\n",
      "    r_s2z_{id}[, k] = sum_to_zero_constrain_brms(",
      "segment(z_s2z_{id}, (k - 1) * (N_{id} - 1) + 1, ",
      "N_{id} - 1));\n",
      "  }}\n",
      "{partial_transform}"
    )
  }

  for (a in seq_along(cross$infos)) {
    info <- cross$infos[[a]]
    p <- info$p
    for (k in seq_along(info$prior)) {
      spec <- info$prior[[k]]
      index <- cross$starts[a] + k - 1L
      loc <- stan_s2z_number(spec$location)
      if (spec$dist == "flat") {
        prec <- "0.0"
      } else if (spec$dist == "normal") {
        prec <- glue("inv_square({stan_s2z_number(spec$scale)})")
      } else {
        prec <- glue(
          "inv_square({stan_s2z_number(spec$scale)} * sqrt(",
          "{stan_s2z_number(spec$df)} * udf_b_s2z{p}_{k}))"
        )
      }
      str_add(out$tpar_comp) <- glue(
        "  prior_mean_s2z_{id}[{index}] = {loc};\n",
        "  prior_prec_s2z_{id}[{index}] = {prec};\n"
      )
    }
  }

  if (!spectral_fisher) {
    str_add(out$tpar_comp) <- stan_re_s2z_cov_group_comp(
      id, varying = varying, is_cor = TRUE, is_student = is_student,
      has_cov = has_cov
    )
  }
  str_add(out$tpar_comp) <- glue(
    "  {{\n",
    "    matrix[{Q}, M_{id}] prior_factor_s2z = diag_pre_multiply(",
    "sqrt(prior_prec_s2z_{id}), H_s2z_{id}) * L_Sigma_s2z_{id};\n",
    "    vector[{Q}] prior_difference_s2z;\n",
    "    vector[M_{id}] h_s2z;\n",
    "    vector[M_{id}] whitened_h_s2z;\n"
  )
  for (a in seq_along(cross$infos)) {
    info <- cross$infos[[a]]
    str_add(out$tpar_comp) <- glue(
      "    prior_difference_s2z[{cross$starts[a]}:{cross$ends[a]}] = ",
      "sqrt(segment(prior_prec_s2z_{id}, ",
      "{cross$starts[a]}, {cross$q[a]})) .* (theta_s2z{info$p} - ",
      "segment(prior_mean_s2z_{id}, {cross$starts[a]}, ",
      "{cross$q[a]}));\n"
    )
  }
  str_add(out$tpar_comp) <- glue(
    "    P_s2z_{id} = crossprod(prior_factor_s2z) + ",
    "P_group_s2z_{id};\n",
    "    h_s2z = prior_factor_s2z' * prior_difference_s2z + ",
    "h_group_s2z_{id};\n",
    "    L_P_s2z_{id} = cholesky_decompose(P_s2z_{id});\n",
    "    whitened_h_s2z = mdivide_left_tri_low(L_P_s2z_{id}, h_s2z);\n",
    "    group_quad_s2z_{id} -= dot_self(whitened_h_s2z);\n",
    "  }}\n",
    cglue("  r_s2z_{idp}_{r$cn} = r_s2z_{id}[, {J}];\n")
  )
  str_add(out$pll_args) <- cglue(
    ", vector r_s2z_{idp}_{r$cn}"
  )

  for (a in seq_along(cross$infos)) {
    info <- cross$infos[[a]]
    p <- info$p
    for (k in seq_along(info$prior)) {
      spec <- info$prior[[k]]
      if (spec$dist == "flat") {
        next
      }
      if (spec$dist == "normal") {
        cond_scale <- stan_s2z_number(spec$scale)
      } else {
        cond_scale <- glue(
          "{stan_s2z_number(spec$scale)} * sqrt(",
          "{stan_s2z_number(spec$df)} * udf_b_s2z{p}_{k})"
        )
      }
      str_add(out$tpar_prior) <- glue(
        "  lprior += normal_{lpdf}(theta_s2z{p}[{k}] | ",
        "{stan_s2z_number(spec$location)}, {cond_scale});\n"
      )
    }
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

  str_add(out$gen_def) <- glue(
    "  vector[M_{id}] mean_r_s2z_{id};\n",
    "  vector[{Q}] q_recovered_s2z_{id};\n",
    "  matrix[N_{id}, M_{id}] r_{id};\n",
    str_if(
      varying,
      glue("  matrix<lower=0>[N_{id}, M_{id}] sd_level_{id};\n")
    )
  )
  for (a in seq_along(cross$infos)) {
    info <- cross$infos[[a]]
    p <- info$p
    if (info$center) {
      str_add(out$gen_def) <- glue(
        "  real Intercept{p};\n",
        str_if(length(info$fixef), glue("  vector[Kc{p}] b{p};\n")),
        "  real b{p}_Intercept;\n"
      )
    } else {
      str_add(out$gen_def) <- glue("  vector[{cross$q[a]}] b{p};\n")
    }
  }
  str_add(out$gen_def) <- cglue(
    "  vector[N_{id}] r_{idp}_{r$cn};\n"
  )
  str_add(out$gen_def) <- glue(
    "  // conventional group-level correlations\n",
    "  corr_matrix[M_{id}] Cor_{id}",
    " = multiply_lower_tri_self_transpose(L_{id});\n",
    "  vector<lower=-1,upper=1>[NC_{id}] cor_{id};\n"
  )

  str_add(out$gen_comp) <- glue(
    "  {{\n",
    "    matrix[{Q}, M_{id}] prior_factor_s2z = diag_pre_multiply(",
    "sqrt(prior_prec_s2z_{id}), H_s2z_{id}) * L_Sigma_s2z_{id};\n",
    "    vector[{Q}] prior_difference_s2z;\n",
    "    vector[M_{id}] h_s2z;\n",
    "    vector[M_{id}] forward_solve_s2z;\n",
    "    vector[M_{id}] r_mean_s2z;\n",
    "    vector[M_{id}] z_mean_s2z;\n"
  )
  for (a in seq_along(cross$infos)) {
    info <- cross$infos[[a]]
    str_add(out$gen_comp) <- glue(
      "    prior_difference_s2z[{cross$starts[a]}:{cross$ends[a]}] = ",
      "sqrt(segment(prior_prec_s2z_{id}, ",
      "{cross$starts[a]}, {cross$q[a]})) .* (theta_s2z{info$p} - ",
      "segment(prior_mean_s2z_{id}, {cross$starts[a]}, ",
      "{cross$q[a]}));\n"
    )
  }
  str_add(out$gen_comp) <- glue(
    "    h_s2z = prior_factor_s2z' * prior_difference_s2z",
    " + h_group_s2z_{id}",
    ";\n",
    "    forward_solve_s2z = mdivide_left_tri_low(L_P_s2z_{id}, h_s2z);\n",
    "    r_mean_s2z = (mdivide_right_tri_low(",
    "forward_solve_s2z', L_P_s2z_{id}))';\n",
    "    for (k in 1:M_{id}) z_mean_s2z[k] = std_normal_rng();\n",
    "    mean_r_s2z_{id} = L_Sigma_s2z_{id} * (r_mean_s2z + ",
    "(mdivide_right_tri_low(z_mean_s2z', L_P_s2z_{id}))');\n",
    "  }}\n"
  )
  for (a in seq_along(cross$infos)) {
    info <- cross$infos[[a]]
    str_add(out$gen_comp) <- glue(
      "  q_recovered_s2z_{id}[{cross$starts[a]}:{cross$ends[a]}] = ",
      "theta_s2z{info$p};\n"
    )
  }
  str_add(out$gen_comp) <- glue(
    "  q_recovered_s2z_{id} -= H_s2z_{id} * mean_r_s2z_{id};\n",
    "  r_{id} = r_s2z_{id};\n",
    "  for (j in 1:N_{id}) r_{id}[j] += mean_r_s2z_{id}';\n",
    str_if(
      varying,
      glue("  sd_level_{id} = sd_level_s2z_{id};\n")
    )
  )
  for (a in seq_along(cross$infos)) {
    info <- cross$infos[[a]]
    p <- info$p
    start <- cross$starts[a]
    if (info$center) {
      str_add(out$gen_comp) <- glue(
        "  Intercept{p} = q_recovered_s2z_{id}[{start}];\n",
        str_if(
          length(info$fixef),
          glue(
            "  b{p} = segment(q_recovered_s2z_{id}, {start + 1L}, ",
            "Kc{p});\n",
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
        "  b{p} = segment(q_recovered_s2z_{id}, {start}, ",
        "{cross$q[a]});\n"
      )
    }
  }
  str_add(out$gen_comp) <- cglue(
    "  r_{idp}_{r$cn} = r_{id}[, {J}];\n"
  )
  str_add(out$gen_comp) <- stan_cor_gen_comp(
    cor = glue("cor_{id}"), ncol = glue("M_{id}")
  )
  out
}
