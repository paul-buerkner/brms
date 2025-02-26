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
stan_predictor.brmsframe <- function(x, prior, normalize, ...) {
  px <- check_prefix(x)
  resp <- usc(combine_prefix(px))
  out <- list()
  str_add_list(out) <- stan_response(x, normalize = normalize, ...)
  valid_dpars <- valid_dpars(x)
  family_files <- family_info(x, "include")
  if (length(family_files)) {
    str_add(out$fun) <- cglue("  #include '{family_files}'\n")
  }
  args <- nlist(prior, normalize, nlpars = names(x$nlpars), ...)
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
                    normalize, ...) {
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
    if (!is_ordinal(family)) {
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
    if (sparse) {
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
stan_re <- function(bframe, prior, normalize, ...) {
  lpdf <- ifelse(normalize, "lpdf", "lupdf")
  reframe <- bframe$frame$re
  stopifnot(is.reframe(reframe))
  IDs <- unique(reframe$id)
  out <- list()
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
  tmp <- lapply(
    IDs, .stan_re, bframe = bframe, prior = prior,
    normalize = normalize, ...
  )
  out <- collapse_lists(ls = c(list(out), tmp))
  out
}

# Stan code for group-level effects per ID
# @param id the ID of the grouping factor
.stan_re <- function(id, bframe, prior, threads, normalize, ...) {
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
      prior, class = "sd", group = r$group[1], coef = r$coef,
      type = glue("matrix[M_{id}, Nby_{id}]"),
      coef_type = glue("row_vector[Nby_{id}]"),
      suffix = glue("_{id}"), px = px, broadcast = "matrix",
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
        prior, class = "sd", group = r$group[1], coef = r$coef,
        type = glue("vector[M_{id}]"), suffix = glue("_{id}"), px = px,
        comment = "group-level standard deviations",
        normalize = normalize
      )
    }
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
    if (r$gtype[1] == "mm") {
      ng <- seq_along(r$gcall[[1]]$groups)
      for (i in seq_rows(r)) {
        str_add(eta_re) <- cglue(
          " + W_{idresp[i]}_{ng}{n}",
          " * r_{idp[i]}_{r$cn[i]}[J_{idresp[i]}_{ng}{n}]",
          " * Z_{idp[i]}_{r$cn[i]}_{ng}{n}"
        )
      }
    } else {
      str_add(eta_re) <- cglue(
        " + r_{idp}_{r$cn}[J_{idresp}{n}] * Z_{idp}_{r$cn}{n}"
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
