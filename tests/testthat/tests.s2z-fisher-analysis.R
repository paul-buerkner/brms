context("Tests for S2Z Fisher analysis")

s2z_fisher_factor_data <- local({
  out <- expand.grid(
    person = factor(seq_len(4L)), item = factor(letters[1:3]),
    KEEP.OUT.ATTRS = FALSE
  )
  out$batch <- factor(rep(seq_len(3L), length.out = nrow(out)))
  out$x <- seq_len(nrow(out)) / nrow(out)
  out$y <- sin(seq_len(nrow(out)))
  out
})

s2z_fisher_factor_frame <- function(outer = alpha + loading * eta,
                                    loading = loading ~ 0 + item,
                                    family = gaussian()) {
  outer <- substitute(outer)
  loading <- substitute(loading)
  outer <- as.formula(call("~", quote(y), outer), env = parent.frame())
  loading <- as.formula(loading, env = parent.frame())
  form <- bf(
    outer,
    alpha ~ 0 + item,
    loading,
    eta ~ 1 + (1 | gr(
      person, s2z = TRUE, center = "auto", latent = TRUE
    )),
    nl = TRUE,
    family = family
  )
  brmsframe(brmsterms(form), data = s2z_fisher_factor_data)
}

test_that("closed-form Fisher rules are response-free", {
  dat <- data.frame(
    y = rep(0:1, 3),
    trials = rep(6L, 6),
    g = factor(rep(seq_len(2), each = 3))
  )

  make_frame <- function(form) {
    brmsframe(brmsterms(form), data = dat)
  }
  com_poisson_frame <- make_frame(bf(
    y ~ 1 + (1 | gr(g, s2z = TRUE, center = "auto")),
    shape ~ 1 + (1 | gr(g, s2z = TRUE, center = "auto")),
    family = brmsfamily("com_poisson")
  ))
  cases <- list(
    saturated_logit = list(
      frame = make_frame(bf(
        y ~ 1 + (1 | gr(g, s2z = TRUE, center = "auto")),
        family = bernoulli()
      )),
      dpar = "mu",
      marker = paste0(
        "(value_fisher_s2z_mu) * ",
        "(1.0 - (value_fisher_s2z_mu))"
      ),
      absent = "square(derivative_fisher_s2z_mu) /"
    ),
    exact = list(
      frame = make_frame(bf(
        y ~ 1 + (1 | gr(g, s2z = TRUE, center = "auto")),
        family = poisson()
      )),
      dpar = "mu", marker = "value_fisher_s2z_mu"
    ),
    coarsened = list(
      frame = make_frame(bf(
        y | trials(trials) ~ 1 +
          (1 | gr(g, s2z = TRUE, center = "auto")),
        phi ~ 1 + (1 | gr(g, s2z = TRUE, center = "auto")),
        family = beta_binomial()
      )),
      dpar = "phi",
      marker = c("trials[n] <= 1 ? 0.0", "pmid_fisher_s2z_bb")
    ),
    moment = list(
      frame = make_frame(bf(
        y ~ 1 + (1 | gr(g, s2z = TRUE, center = "auto")),
        shape ~ 1 + (1 | gr(g, s2z = TRUE, center = "auto")),
        family = negbinomial()
      )),
      dpar = "shape", marker = "log_shape_info_fisher_s2z_nb"
    ),
    atom = list(
      frame = make_frame(bf(
        y ~ 1 + (1 | gr(g, s2z = TRUE, center = "auto")),
        zi ~ 1 + (1 | gr(g, s2z = TRUE, center = "auto")),
        family = zero_inflated_poisson()
      )),
      dpar = "mu", marker = "atom_derivative_fisher_s2z_zi"
    ),
    com_poisson_location = list(
      frame = com_poisson_frame,
      dpar = "mu",
      marker = c(
        paste0(
          "log(fmax(value_fisher_s2z_mu, 1e-12)) / ",
          "fmax(value_fisher_s2z_shape, 1e-12)"
        ),
        paste0(
          "mode_fisher_s2z_cmp / ",
          "fmax(value_fisher_s2z_shape, 1e-12)"
        )
      ),
      obs_prec = "variance_fisher_s2z_cmp"
    ),
    com_poisson_shape = list(
      frame = com_poisson_frame,
      dpar = "shape",
      marker = c(
        "log_factorial_slope_fisher_s2z_cmp",
        "log_factorial_curve_fisher_s2z_cmp",
        "log_factorial_variance_fisher_s2z_cmp"
      ),
      obs_prec = paste0(
        "square(derivative_fisher_s2z_shape) * ",
        "log_factorial_variance_fisher_s2z_cmp"
      )
    )
  )

  for (name in names(cases)) {
    case <- cases[[name]]
    got <- stan_re_s2z_fisher_closed_form(case$frame, case$dpar)
    expect_false(is.null(got), info = name)
    emitted <- paste0(got$definitions_at_n, got$obs_prec_at_n)
    expect_false(grepl("Y[n]", emitted, fixed = TRUE), info = name)
    expect_true(all(vapply(
      case$marker, grepl, logical(1), x = emitted, fixed = TRUE
    )), info = name)
    if (!is.null(case$absent)) {
      expect_false(grepl(case$absent, emitted, fixed = TRUE), info = name)
    }
    if (!is.null(case$obs_prec)) {
      expect_identical(
        as.character(got$obs_prec_at_n), case$obs_prec, info = name
      )
    }
  }
})

test_that("the scalar native Fisher catalog is exhaustive", {
  N <- 12L
  dat <- data.frame(
    g = factor(rep(letters[1:4], each = 3L)),
    y_real = seq(-1.2, 1.4, length.out = N),
    y_pos = seq(0.2, 3, length.out = N),
    y_unit = seq(0.05, 0.95, length.out = N),
    y_count = rep(0:5, length.out = N),
    y_bin = rep(0:1, length.out = N),
    trials = rep(6L, N),
    dec = rep(0:1, length.out = N)
  )
  term <- "(1 | gr(g, s2z = TRUE, center = 'auto'))"
  make_frame <- function(response, family, target) {
    rhs <- if (target == "mu") paste("1 +", term) else "1"
    main <- as.formula(paste(response, "~", rhs))
    extra <- if (target == "mu") list() else {
      list(as.formula(paste(target, "~ 1 +", term)))
    }
    form <- do.call(bf, c(list(main), extra, list(family = family)))
    brmsframe(brmsterms(form), data = dat)
  }
  cases <- list(
    gaussian = list(gaussian(), c("mu", "sigma"), "y_real"),
    student = list(student(), c("mu", "sigma", "nu"), "y_real"),
    skew_normal = list(
      skew_normal(), c("mu", "sigma", "alpha"), "y_real"
    ),
    bernoulli = list(bernoulli(), "mu", "y_bin"),
    binomial = list(binomial(), "mu", "y_bin | trials(trials)"),
    xbeta = list(xbeta(), c("mu", "phi", "kappa"), "y_unit"),
    beta_binomial = list(
      beta_binomial(), c("mu", "phi"), "y_count | trials(trials)"
    ),
    beta = list(Beta(), c("mu", "phi"), "y_unit"),
    poisson = list(poisson(), "mu", "y_count"),
    negbinomial = list(
      negbinomial(), c("mu", "shape"), "y_count"
    ),
    negbinomial2 = list(
      negbinomial2(), c("mu", "sigma"), "y_count"
    ),
    geometric = list(geometric(), "mu", "y_count"),
    discrete_weibull = list(
      discrete_weibull(), c("mu", "shape"), "y_count"
    ),
    com_poisson = list(
      com_poisson(), c("mu", "shape"), "y_count"
    ),
    gamma = list(Gamma(), c("mu", "shape"), "y_pos"),
    weibull = list(weibull(), c("mu", "shape"), "y_pos"),
    exponential = list(exponential(), "mu", "y_pos"),
    frechet = list(frechet(), c("mu", "nu"), "y_pos"),
    inverse_gaussian = list(
      brmsfamily("inverse.gaussian"), c("mu", "shape"), "y_pos"
    ),
    lognormal = list(lognormal(), c("mu", "sigma"), "y_pos"),
    shifted_lognormal = list(
      shifted_lognormal(), c("mu", "sigma", "ndt"), "y_pos"
    ),
    exgaussian = list(
      exgaussian(), c("mu", "sigma", "beta"), "y_real"
    ),
    von_mises = list(von_mises(), c("mu", "kappa"), "y_real"),
    asym_laplace = list(
      asym_laplace(), c("mu", "sigma", "quantile"), "y_real"
    ),
    hurdle_poisson = list(
      hurdle_poisson(), c("mu", "hu"), "y_count"
    ),
    hurdle_negbinomial = list(
      hurdle_negbinomial(), c("mu", "shape", "hu"), "y_count"
    ),
    hurdle_gamma = list(
      hurdle_gamma(), c("mu", "shape", "hu"), "y_pos"
    ),
    hurdle_lognormal = list(
      hurdle_lognormal(), c("mu", "sigma", "hu"), "y_pos"
    ),
    zero_inflated_poisson = list(
      zero_inflated_poisson(), c("mu", "zi"), "y_count"
    ),
    zero_inflated_negbinomial = list(
      zero_inflated_negbinomial(), c("mu", "shape", "zi"), "y_count"
    ),
    zero_inflated_binomial = list(
      zero_inflated_binomial(), c("mu", "zi"),
      "y_count | trials(trials)"
    ),
    zero_inflated_beta_binomial = list(
      zero_inflated_beta_binomial(), c("mu", "phi", "zi"),
      "y_count | trials(trials)"
    ),
    zero_inflated_beta = list(
      zero_inflated_beta(), c("mu", "phi", "zi"), "y_unit"
    ),
    zero_one_inflated_beta = list(
      zero_one_inflated_beta(), c("mu", "phi", "zoi", "coi"),
      "y_unit"
    ),
    zero_inflated_asym_laplace = list(
      zero_inflated_asym_laplace(),
      c("mu", "sigma", "quantile", "zi"), "y_real"
    ),
    cox = list(cox(), "mu", "y_pos"),
    wiener = list(wiener(), c("mu", "bs", "bias"),
                  "y_pos | dec(dec)")
  )

  checked <- character()
  for (family_name in names(cases)) {
    case <- cases[[family_name]]
    for (target in case[[2L]]) {
      label <- paste(family_name, target, sep = "/")
      bframe <- make_frame(case[[3L]], case[[1L]], target)
      got <- stan_re_s2z_fisher_closed_form(bframe, target)
      expect_false(is.null(got), info = label)
      emitted <- paste0(got$definitions_at_n, got$obs_prec_at_n)
      expect_false(grepl("Y[n]", emitted, fixed = TRUE), info = label)
      checked <- c(checked, label)
    }
  }
  expect_length(checked, 85L)
})

test_that("the multicategory native Fisher catalog is exhaustive", {
  N <- 12L
  g <- factor(rep(letters[1:4], each = 3L))
  categorical_data <- data.frame(
    y = factor(rep(c("a", "b", "c"), 4L)), g = g
  )
  count_data <- data.frame(
    y1 = rep(c(1L, 2L, 3L), 4L),
    y2 = rep(c(2L, 1L, 2L), 4L),
    y3 = rep(c(3L, 2L, 1L), 4L),
    g = g
  )
  count_data$size <- rowSums(count_data[c("y1", "y2", "y3")])
  count_data$y <- I(as.matrix(count_data[c("y1", "y2", "y3")]))
  simplex_data <- data.frame(
    y1 = seq(0.15, 0.25, length.out = N),
    y2 = seq(0.25, 0.35, length.out = N),
    g = g
  )
  simplex_data$y3 <- 1 - simplex_data$y1 - simplex_data$y2
  simplex_data$y <- I(as.matrix(
    simplex_data[c("y1", "y2", "y3")]
  ))
  term <- "(1 | gr(g, s2z = TRUE, center = 'auto'))"
  cases <- list(
    categorical = list(
      bf(y ~ 1, as.formula(paste("mub ~ 1 +", term)),
         family = categorical()),
      categorical_data, "mub"
    ),
    multinomial = list(
      bf(
        y | trials(size) ~ 1,
        as.formula(paste("muy2 ~ 1 +", term)),
        family = multinomial()
      ),
      count_data, "muy2"
    ),
    dirichlet_mu = list(
      bf(y ~ 1, as.formula(paste("muy2 ~ 1 +", term)),
         family = dirichlet()),
      simplex_data, "muy2"
    ),
    dirichlet_phi = list(
      bf(y ~ 1, as.formula(paste("phi ~ 1 +", term)),
         family = dirichlet()),
      simplex_data, "phi"
    ),
    dirichlet2 = list(
      bf(y ~ 1, as.formula(paste("muy1 ~ 1 +", term)),
         family = brmsfamily("dirichlet2")),
      simplex_data, "muy1"
    ),
    dirichlet_multinomial_mu = list(
      bf(
        y | trials(size) ~ 1,
        as.formula(paste("muy2 ~ 1 +", term)),
        family = dirichlet_multinomial()
      ),
      count_data, "muy2"
    ),
    dirichlet_multinomial_phi = list(
      bf(
        y | trials(size) ~ 1,
        as.formula(paste("phi ~ 1 +", term)),
        family = dirichlet_multinomial()
      ),
      count_data, "phi"
    ),
    logistic_normal_mu = list(
      bf(y ~ 1, as.formula(paste("muy2 ~ 1 +", term)),
         family = logistic_normal()),
      simplex_data, "muy2"
    ),
    logistic_normal_sigma = list(
      bf(y ~ 1, as.formula(paste("sigmay2 ~ 1 +", term)),
         family = logistic_normal()),
      simplex_data, "sigmay2"
    )
  )

  for (name in names(cases)) {
    case <- cases[[name]]
    formula <- validate_formula(case[[1L]], data = case[[2L]])
    bframe <- brmsframe(brmsterms(formula), data = case[[2L]])
    got <- stan_re_s2z_fisher_closed_form(bframe, case[[3L]])
    expect_false(is.null(got), info = name)
    emitted <- paste0(got$definitions_at_n, got$obs_prec_at_n)
    expect_false(grepl("Y[n]", emitted, fixed = TRUE), info = name)
  }
})

test_that("Gaussian nonlinear loading derivatives are analyzed for codegen", {
  bframe <- s2z_fisher_factor_frame()
  got <- re_s2z_fisher_nl_info(bframe, bframe$nlpars$eta)

  expect_true(got$supported)
  expect_equal(got$id, 1)
  expect_identical(got$target_nlpar, "eta")
  expect_identical(got$derivative, quote(loading))
  expect_identical(got$dependencies, "loading")
  expect_identical(got$nlpar_dependencies, "loading")
  expect_identical(
    got$outer_nlpar_dependencies, c("alpha", "loading", "eta")
  )
  expect_length(got$covariate_dependencies, 0L)
  expect_true(got$is_gaussian_identity)
  expect_false(got$is_student_identity)
  expect_true(got$is_location_identity)
  expect_identical(got$response_family, "gaussian")
  expect_true(got$is_observation_local)
  expect_true(got$score_independent)
  expect_false(got$has_response_addition_terms)
  expect_length(got$response_addition_terms, 0L)
  expect_identical(got$sigma_kind, "scalar_parameter")
  expect_true(got$sigma_is_scalar)
  expect_null(got$sigma_fixed_value)
  expect_identical(got$sigma_parameter, "sigma")
  expect_identical(got$obs_derivative_stan, "fisher_s2z_1_nlp_loading[n]")
  expect_match(
    got$outer_reference_stan, "fisher_s2z_1_nlp_alpha[n]", fixed = TRUE
  )
  expect_match(
    got$outer_reference_stan, "fisher_s2z_1_nlp_loading[n]", fixed = TRUE
  )
  expect_match(got$outer_reference_stan, "X_eta[n]", fixed = TRUE)
  expect_identical(got$dzeta_dxi_stan, "fisher_s2z_1_nlp_loading[n]")
  expect_identical(got$cn, 1L)
  expect_identical(got$frame_rows, "1")
  expect_identical(got$row_keys, "1:eta:1")
  expect_identical(
    got$row_metadata[c("row_key", "nlpar", "coef", "cn")],
    data.frame(
      row_key = "1:eta:1", nlpar = "eta", coef = "Intercept", cn = 1L,
      stringsAsFactors = FALSE
    )
  )

  loading <- got$dependency_info$loading
  expect_setequal(names(got$dependency_info), c("alpha", "loading"))
  expect_identical(loading$x_name, "X_loading")
  expect_identical(loading$coefficient_name, "b_loading")
  expect_identical(loading$vector_name, "fisher_s2z_1_nlp_loading")
  expect_identical(
    loading$vector_expression, "X_loading * b_loading"
  )
})

test_that("nonlinear Fisher dependency references retain fixed offsets", {
  bframe <- s2z_fisher_factor_frame(
    loading = loading ~ 0 + item + offset(x)
  )
  got <- re_s2z_fisher_nl_info(bframe, bframe$nlpars$eta)

  expect_identical(
    got$dependency_info$loading$vector_expression,
    "(X_loading * b_loading) + offsets_loading"
  )
  expect_match(
    got$outer_reference_stan,
    "fisher_s2z_1_nlp_loading[n]", fixed = TRUE
  )

  offset_only <- s2z_fisher_factor_frame(
    loading = loading ~ 0 + offset(x)
  )
  got <- re_s2z_fisher_nl_info(
    offset_only, offset_only$nlpars$eta
  )
  expect_identical(
    got$dependency_info$loading$vector_expression,
    "offsets_loading"
  )
})

test_that("dependency aliases are isolated across Fisher factor IDs", {
  form <- bf(
    y ~ alpha + loading * eta,
    alpha ~ 0 + item,
    loading ~ 0 + item,
    eta ~ 1 +
      (1 | gr(
        person, id = "person_factor", s2z = TRUE,
        center = "auto", latent = TRUE
      )) +
      (1 | gr(
        batch, id = "batch_factor", s2z = TRUE,
        center = "auto", latent = TRUE
      )),
    nl = TRUE
  )
  bframe <- brmsframe(brmsterms(form), data = s2z_fisher_factor_data)
  ids <- unique(bframe$nlpars$eta$frame$re$id)
  expect_length(ids, 2L)
  expect_error(
    re_s2z_fisher_nl_info(bframe, bframe$nlpars$eta),
    "requires an explicit group-level ID"
  )

  got <- lapply(ids, function(id) {
    re_s2z_fisher_nl_info(bframe, bframe$nlpars$eta, id = id)
  })
  aliases <- vapply(
    got, function(x) x$dependency_info$loading$vector_name, character(1)
  )
  expect_identical(
    aliases, paste0("fisher_s2z_", ids, "_nlp_loading")
  )
  expect_length(unique(aliases), 2L)
})

test_that("row metadata aligns dimensions shared across nonlinear targets", {
  form <- bf(
    y ~ alpha + loading1 * eta1 + loading2 * eta2,
    alpha ~ 0 + item,
    loading1 ~ 0 + item,
    loading2 ~ 0 + item,
    eta1 ~ 1 + (1 | gr(
      person, id = "factor", s2z = TRUE,
      center = "auto", latent = TRUE
    )),
    eta2 ~ 1 + (1 | gr(
      person, id = "factor", s2z = TRUE,
      center = "auto", latent = TRUE
    )),
    nl = TRUE
  )
  bframe <- brmsframe(brmsterms(form), data = s2z_fisher_factor_data)
  r <- bframe$frame$re
  id <- unique(r$id)
  expect_length(id, 1L)

  got <- lapply(c("eta1", "eta2"), function(nlpar) {
    re_s2z_fisher_nl_info(
      bframe, bframe$nlpars[[nlpar]], id = id
    )
  })
  metadata <- do.call(rbind, lapply(got, `[[`, "row_metadata"))
  global_keys <- paste(r$id, r$nlpar, r$cn, sep = ":")
  aligned <- metadata[match(global_keys, metadata$row_key), , drop = FALSE]

  expect_false(anyNA(match(global_keys, metadata$row_key)))
  expect_identical(aligned$row_key, global_keys)
  expect_identical(aligned$cn, unname(r$cn))
  expect_identical(aligned$nlpar, unname(r$nlpar))
  expect_identical(vapply(got, `[[`, character(1), "target_nlpar"),
                   c("eta1", "eta2"))
  expect_identical(got[[1]]$dzeta_dxi_stan,
                   "fisher_s2z_1_nlp_loading1[n]")
  expect_identical(got[[2]]$dzeta_dxi_stan,
                   "fisher_s2z_1_nlp_loading2[n]")
})

test_that("nonlinear Fisher analysis maps safe observation covariates", {
  bframe <- s2z_fisher_factor_frame(
    outer = alpha + x * loading * eta
  )
  got <- re_s2z_fisher_nl_info(bframe, bframe$nlpars$eta)

  expect_setequal(got$dependencies, c("x", "loading"))
  expect_identical(got$covariate_dependencies, "x")
  expect_match(got$obs_derivative_stan, "C_1[n]", fixed = TRUE)
  expect_match(
    got$obs_derivative_stan, "fisher_s2z_1_nlp_loading[n]", fixed = TRUE
  )
})

test_that("nonlinear Fisher analysis rejects score-dependent derivatives", {
  bframe <- s2z_fisher_factor_frame(
    outer = alpha + loading * eta^2
  )
  expect_error(
    re_s2z_fisher_nl_info(bframe, bframe$nlpars$eta),
    "depends on target nonlinear parameter 'eta'"
  )

  soft <- re_s2z_fisher_nl_info(
    bframe, bframe$nlpars$eta, strict = FALSE
  )
  expect_false(soft$supported)
  expect_match(soft$reason_unsupported, "target nonlinear parameter 'eta'")
})

test_that("nonlinear Fisher analysis rejects latent loading predictors", {
  bframe <- s2z_fisher_factor_frame(
    loading = loading ~ 0 + item + (1 | batch)
  )
  expect_error(
    re_s2z_fisher_nl_info(bframe, bframe$nlpars$eta),
    "latent or group-varying nonlinear parameter 'loading'"
  )
})

test_that("nonlinear Fisher analysis rejects centered loading predictors", {
  form <- bf(
    y ~ alpha + loading * eta,
    alpha ~ 0 + item,
    lf(loading ~ item, center = TRUE),
    eta ~ 1 + (1 | gr(
      person, s2z = TRUE, center = "auto", latent = TRUE
    )),
    nl = TRUE
  )
  bframe <- brmsframe(
    brmsterms(form), data = s2z_fisher_factor_data
  )
  expect_error(
    re_s2z_fisher_nl_info(bframe, bframe$nlpars$eta),
    "centered population-level predictors.*'loading'"
  )
})

test_that("nonlinear Fisher analysis reports likelihood capability metadata", {
  weighted <- bf(
    y | weights(x) ~ alpha + loading * eta,
    alpha ~ 0 + item,
    loading ~ 0 + item,
    eta ~ 1 + (1 | gr(
      person, s2z = TRUE, center = "auto", latent = TRUE
    )),
    nl = TRUE
  )
  bframe <- brmsframe(
    brmsterms(weighted), data = s2z_fisher_factor_data
  )
  got <- re_s2z_fisher_nl_info(bframe, bframe$nlpars$eta)
  expect_true(got$has_response_addition_terms)
  expect_identical(got$response_addition_terms, "weights")

  predicted_sigma <- bf(
    y ~ alpha + loading * eta,
    alpha ~ 0 + item,
    loading ~ 0 + item,
    eta ~ 1 + (1 | gr(
      person, s2z = TRUE, center = "auto", latent = TRUE
    )),
    sigma ~ 0 + item,
    nl = TRUE
  )
  bframe <- brmsframe(
    brmsterms(predicted_sigma), data = s2z_fisher_factor_data
  )
  got <- re_s2z_fisher_nl_info(bframe, bframe$nlpars$eta)
  expect_identical(got$sigma_kind, "predictor")
  expect_false(got$sigma_is_scalar)

  student_frame <- s2z_fisher_factor_frame(family = student())
  got <- re_s2z_fisher_nl_info(
    student_frame, student_frame$nlpars$eta
  )
  expect_identical(got$response_family, "student")
  expect_false(got$is_gaussian_identity)
  expect_true(got$is_student_identity)
  expect_true(got$is_location_identity)
})

test_that("constant sampled loadings retain the ordinary b symbol", {
  form <- bf(
    y ~ alpha + loading * eta,
    alpha ~ 0 + item,
    loading ~ 0 + item,
    eta ~ 1 + (1 | gr(
      person, s2z = TRUE, center = FALSE, latent = TRUE
    )),
    nl = TRUE
  )
  scode <- stancode(
    form, data = s2z_fisher_factor_data,
    prior = prior(constant(1), class = "b", nlpar = "loading")
  )
  expect_match2(scode, "vector[K_loading] b_loading;")
  expect_match2(scode, "b_loading = rep_vector(1, rows(b_loading));")
  expect_match2(scode, "nlp_loading += X_loading * b_loading;")
})

test_that("nonlinear Fisher analysis supports scalar non-Gaussian locations", {
  bernoulli_data <- s2z_fisher_factor_data
  bernoulli_data$y <- rep(0:1, length.out = nrow(bernoulli_data))
  form <- bf(
    y ~ alpha + loading * eta,
    alpha ~ 0 + item,
    loading ~ 0 + item,
    eta ~ 1 + (1 | gr(
      person, s2z = TRUE, center = "auto", latent = TRUE
    )),
    nl = TRUE,
    family = bernoulli()
  )
  bframe <- brmsframe(brmsterms(form), data = bernoulli_data)
  got <- re_s2z_fisher_nl_info(bframe, bframe$nlpars$eta)
  expect_true(got$supported)
  expect_identical(got$response_family, "bernoulli")
  expect_false(got$is_location_identity)
  expect_match(
    got$outer_reference_stan, "fisher_s2z_1_nlp_alpha[n]", fixed = TRUE
  )

  unknown <- bf(
    y ~ alpha + mystery(eta),
    alpha ~ 0 + item,
    eta ~ 1 + (1 | gr(
      person, s2z = TRUE, center = "auto", latent = TRUE
    )),
    nl = TRUE
  )
  bframe <- brmsframe(
    brmsterms(unknown), data = s2z_fisher_factor_data
  )
  expect_error(
    re_s2z_fisher_nl_info(bframe, bframe$nlpars$eta),
    "unsupported derivative function 'mystery'"
  )
})

test_that("strict latent references retain population score coordinates", {
  with_intercept <- s2z_fisher_factor_frame(family = bernoulli())
  got <- re_s2z_fisher_nl_info(
    with_intercept, with_intercept$nlpars$eta
  )
  expect_match(got$outer_reference_stan, "X_eta[n]", fixed = TRUE)
  expect_match(got$outer_reference_stan, "b_eta", fixed = TRUE)

  no_intercept <- bf(
    y ~ alpha + loading * eta,
    alpha ~ 0 + item,
    loading ~ 0 + item,
    eta ~ 0 + (1 | gr(
      person, s2z = TRUE, center = "auto", latent = TRUE
    )),
    nl = TRUE,
    family = bernoulli()
  )
  dat <- s2z_fisher_factor_data
  dat$y <- rep(0:1, length.out = nrow(dat))
  bframe <- brmsframe(brmsterms(no_intercept), data = dat)
  got <- re_s2z_fisher_nl_info(bframe, bframe$nlpars$eta)
  expect_match(
    got$outer_reference_stan,
    "fisher_s2z_1_nlp_loading[n] * 0", fixed = TRUE
  )
  expect_false(grepl("X_eta", got$outer_reference_stan, fixed = TRUE))
})
