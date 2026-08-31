context("Tests for physical sum-to-zero group-level effects")

expect_match2 <- brms:::expect_match2

s2z_count_fixed <- function(x, pattern) {
  unname(lengths(regmatches(
    x, gregexpr(pattern, x, fixed = TRUE)
  )))
}

s2z_stan_between <- function(x, start, end) {
  start_at <- regexpr(start, x, fixed = TRUE)[1L]
  expect_gt(start_at, 0L)
  end_at <- regexpr(
    end, substring(x, start_at + nchar(start)), fixed = TRUE
  )[1L]
  expect_gt(end_at, 0L)
  substring(
    x, start_at + nchar(start),
    start_at + nchar(start) + end_at - 2L
  )
}

s2z_dat <- local({
  set.seed(1916)
  n <- 72
  data.frame(
    y = rnorm(n),
    x = seq(1.25, 4.75, length.out = n),
    z = rep(c(-0.5, 2), length.out = n),
    f = factor(rep(c("a", "b", "c"), length.out = n)),
    g = factor(rep(seq_len(6), each = 12)),
    h = factor(rep(seq_len(8), each = 9)),
    w = rep(seq(0.8, 1.2, length.out = 6), each = 12)
  )
})

s2z_ten_dat <- local({
  data.frame(
    y = seq(-1, 1, length.out = 80),
    ten = factor(rep(letters[1:10], 8)),
    g = factor(rep(seq_len(8), each = 10))
  )
})

test_that("S2Z centering API defaults compatibly and reaches the reframe", {
  default_form <- y ~ x + (1 + x | gr(g, s2z = TRUE))
  centered_form <- y ~ x +
    (1 + x | gr(g, s2z = TRUE, center = TRUE))
  noncentered_form <- y ~ x +
    (1 + x | gr(g, s2z = TRUE, center = FALSE))
  partial_form <- y ~ x +
    (1 + x | gr(g, s2z = TRUE, center = 0.35))
  auto_form <- y ~ x +
    (1 + x | gr(g, s2z = TRUE, center = "auto"))
  fisher_form <- y ~ x +
    (1 + x | gr(g, s2z = TRUE, center = "fisher"))

  default_terms <- brmsterms(default_form)
  centered_terms <- brmsterms(centered_form)
  noncentered_terms <- brmsterms(noncentered_form)
  partial_terms <- brmsterms(partial_form)
  auto_terms <- brmsterms(auto_form)
  fisher_terms <- brmsterms(fisher_form)
  expect_equal(default_terms$dpars$mu$re$gcall[[1]]$s2z_center, 1)
  expect_equal(centered_terms$dpars$mu$re$gcall[[1]]$s2z_center, 1)
  expect_equal(noncentered_terms$dpars$mu$re$gcall[[1]]$s2z_center, 0)
  expect_equal(partial_terms$dpars$mu$re$gcall[[1]]$s2z_center, 0.35)
  expect_equal(auto_terms$dpars$mu$re$gcall[[1]]$s2z_center, 0.5)
  expect_false(default_terms$dpars$mu$re$gcall[[1]]$s2z_center_auto)
  expect_true(auto_terms$dpars$mu$re$gcall[[1]]$s2z_center_auto)
  expect_true(fisher_terms$dpars$mu$re$gcall[[1]]$s2z_center_auto)
  expect_equal(brms:::frame_re(default_terms, s2z_dat)$s2z_center, c(1, 1))
  expect_equal(brms:::frame_re(centered_terms, s2z_dat)$s2z_center, c(1, 1))
  expect_equal(
    brms:::frame_re(noncentered_terms, s2z_dat)$s2z_center, c(0, 0)
  )
  expect_equal(
    brms:::frame_re(partial_terms, s2z_dat)$s2z_center, c(0.35, 0.35)
  )
  expect_true(all(
    brms:::frame_re(auto_terms, s2z_dat)$s2z_center_auto
  ))
  auto_reframe <- brms:::frame_re(auto_terms, s2z_dat)
  expect_identical(brms:::re_s2z_center_mode(auto_reframe), "auto")
  expect_identical(
    stancode(fisher_form, data = s2z_dat),
    stancode(auto_form, data = s2z_dat)
  )

  default_code <- stancode(default_form, data = s2z_dat)
  centered_code <- stancode(centered_form, data = s2z_dat)
  noncentered_code <- stancode(noncentered_form, data = s2z_dat)
  expect_identical(default_code, centered_code)
  expect_false(identical(default_code, noncentered_code))
  expect_identical(
    stancode(
      y ~ x + (1 + x | gr(g, s2z = TRUE, center = 1)),
      data = s2z_dat
    ),
    centered_code
  )
  expect_identical(
    stancode(
      y ~ x + (1 + x | gr(g, s2z = TRUE, center = 0)),
      data = s2z_dat
    ),
    noncentered_code
  )
  expect_null(standata(centered_form, data = s2z_dat)$rho_s2z_1)
  expect_null(standata(noncentered_form, data = s2z_dat)$rho_s2z_1)
  expect_null(standata(auto_form, data = s2z_dat)$rho_s2z_1)

  # Reframes made from the original S2Z formula representation did not carry
  # s2z_center. They must retain the original centered behavior.
  legacy_terms <- default_terms
  legacy_terms$dpars$mu$re$gcall[[1]]$s2z_center <- NULL
  legacy_terms$dpars$mu$re$gcall[[1]]$s2z_center_auto <- NULL
  expect_equal(
    brms:::frame_re(legacy_terms, s2z_dat)$s2z_center, c(1, 1)
  )
  expect_false(any(
    brms:::frame_re(legacy_terms, s2z_dat)$s2z_center_auto
  ))

  conventional <- stancode(y ~ x + (1 + x | gr(g)), data = s2z_dat)
  conventional_false <- stancode(
    y ~ x + (1 + x | gr(g, center = FALSE)), data = s2z_dat
  )
  conventional_zero <- stancode(
    y ~ x + (1 + x | gr(g, center = 0)), data = s2z_dat
  )
  expect_identical(conventional_false, conventional)
  expect_identical(conventional_zero, conventional)
  conventional_terms <- brmsterms(y ~ x + (1 + x | gr(g)))
  expect_equal(
    brms:::frame_re(conventional_terms, s2z_dat)$s2z_center, c(0, 0)
  )
  conventional_centered <- stancode(
    y ~ x + (1 + x | gr(g, center = TRUE)), data = s2z_dat
  )
  expect_identical(
    conventional_centered,
    stancode(y ~ x + (1 + x | gr(g, center = 1)), data = s2z_dat)
  )
  expect_false(identical(conventional_centered, conventional))
  conventional_partial <- stancode(
    y ~ x + (1 + x | gr(g, center = 0.35)), data = s2z_dat
  )
  expect_false(identical(conventional_partial, conventional_centered))
  expect_identical(
    stancode(
      y ~ x + (1 + x | gr(g, center = "fisher")), data = s2z_dat
    ),
    stancode(
      y ~ x + (1 + x | gr(g, center = "auto")), data = s2z_dat
    )
  )
  conventional_partial_terms <- brmsterms(
    y ~ x + (1 + x | gr(g, center = 0.35))
  )
  expect_equal(
    brms:::frame_re(
      conventional_partial_terms, s2z_dat
    )$s2z_center,
    c(0.35, 0.35)
  )
  # Old ordinary reframes did not carry centering metadata and must retain
  # their historical non-centered default.
  legacy_conventional <- conventional_terms
  legacy_conventional$dpars$mu$re$gcall[[1]]$s2z_center <- NULL
  legacy_conventional$dpars$mu$re$gcall[[1]]$s2z_center_auto <- NULL
  expect_equal(
    brms:::frame_re(legacy_conventional, s2z_dat)$s2z_center, c(0, 0)
  )
  expect_error(gr(g, s2z = TRUE, center = NA), "center")
  expect_error(
    gr(g, s2z = TRUE, center = c(TRUE, FALSE)), "center"
  )
  for (value in list(-0.01, 1.01, NaN, Inf, "AUTO", "partial")) {
    expect_error(gr(g, s2z = TRUE, center = value), "center")
  }

  mixed_data <- standata(
    y ~ x +
      (1 | gr(g, id = "mixed", s2z = TRUE, center = 0.2)) +
      (0 + x | gr(g, id = "mixed", s2z = TRUE, center = 0.8)),
    data = s2z_dat
  )
  expect_equal(unname(mixed_data$rho_s2z_1[, 1]), rep(0.2, 6))
  expect_equal(unname(mixed_data$rho_s2z_1[, 2]), rep(0.8, 6))

  expect_error(
    standata(
      y ~ x +
        (1 | gr(g, id = "mixed-auto", s2z = TRUE,
                center = "auto")) +
        (0 + x | gr(g, id = "mixed-auto", s2z = TRUE,
                    center = 0.8)),
      data = s2z_dat
    ),
    "must use Fisher centering if any coefficient does",
    fixed = TRUE
  )
})

test_that("scalar Fisher S2Z hoists exposure and uses closed-form rho", {
  form <- y ~ 1 +
    (1 | gr(g, s2z = TRUE, center = "fisher"))
  scode <- stancode(form, data = s2z_dat)
  tdata <- s2z_stan_between(
    scode, "transformed data {", "\nparameters {"
  )
  tpar <- s2z_stan_between(
    scode, "transformed parameters {", "\nmodel {"
  )

  expect_match2(
    scode, "matrix<lower=0,upper=1>[N_1, M_1] rho_s2z_1;"
  )
  expect_match2(
    scode, "vector<lower=0,upper=1>[M_1] mean_rho_s2z_1;"
  )
  expect_match2(
    tdata, "vector<lower=0>[N_1] exposure_fisher_s2z_1;"
  )
  expect_match2(tdata, "exposure_fisher_s2z_1 = zeros_vector(N_1);")
  expect_match2(
    tdata,
    "exposure_fisher_s2z_1[J_1[n]] += square(Z_1_1[n]);"
  )
  expect_match2(tpar, "inv_square(sigma)")
  expect_match2(tpar, "rho_s2z_1[j, 1]")
  expect_match2(tpar, "exposure_fisher_s2z_1[j]")
  expect_match2(
    tpar,
    paste0(
      "real scaled_info_fisher_s2z = 1.0 * square(sd_1[1]) * ",
      "obs_prec_fisher_s2z * exposure_fisher_s2z_1[j];"
    )
  )
  expect_match2(
    tpar,
    "rho_s2z_1[j, 1] = 1.0 - inv(1.0 + scaled_info_fisher_s2z);"
  )
  expect_false(grepl("J_1[n]", tpar, fixed = TRUE))
  expect_false(grepl("mdivide_left_spd", scode, fixed = TRUE))
  expect_false(grepl("quad_form(", scode, fixed = TRUE))
  expect_match2(scode, "scale_partial_s2z = 1.0 - rho_s2z_1[, 1]")
  expect_match2(scode, "+ log_det_partial_s2z_1")
  expect_false(grepl("eigenvectors_sym", scode, fixed = TRUE))
  expect_false(grepl("inverse_spd", scode, fixed = TRUE))
  expect_false(grepl("if (mean_rho_s2z_1", scode, fixed = TRUE))

  student_code <- stancode(form, data = s2z_dat, family = student())
  student_tpar <- s2z_stan_between(
    student_code, "transformed parameters {", "\nmodel {"
  )
  expect_match2(
    student_tpar,
    "(nu + 1.0) / (nu + 3.0) * inv_square(sigma)"
  )
  expect_match2(student_code, "+ log_det_partial_s2z_1")
})

test_that("multivariate Fisher S2Z hoists Gram matrices and solves diagonals", {
  form <- y ~ x +
    (1 + x | gr(g, s2z = TRUE, center = "fisher"))
  scode <- stancode(form, data = s2z_dat)
  tdata <- s2z_stan_between(
    scode, "transformed data {", "\nparameters {"
  )
  tpar <- s2z_stan_between(
    scode, "transformed parameters {", "\nmodel {"
  )

  expect_match2(tdata, "array[N_1] matrix[M_1, M_1]")
  expect_match2(tdata, "J_1[n]")
  expect_match2(tdata, "Z_1_1[n]")
  expect_match2(tdata, "Z_1_2[n]")
  expect_false(grepl("J_1[n]", tpar, fixed = TRUE))
  expect_false(grepl("design_fisher_s2z", tpar, fixed = TRUE))
  expect_match2(tpar, "inv_square(sigma)")
  expect_match2(tpar, "cholesky_decompose(")
  expect_match2(tpar, "mdivide_left_tri_low(")
  expect_match2(tpar, "L_post_precision_fisher_s2z")
  expect_match2(tpar, "white_factor_fisher_s2z")
  expect_match2(tpar, "post_var_fisher_s2z = columns_dot_self(")
  expect_false(grepl("mdivide_left_spd", scode, fixed = TRUE))
  expect_equal(s2z_count_fixed(tpar, "quad_form("), 1L)
  expect_false(grepl("identity_matrix(M_1)", scode, fixed = TRUE))
  expect_false(grepl("white_cov_fisher_s2z", scode, fixed = TRUE))
  expect_false(grepl("post_cov_fisher_s2z", scode, fixed = TRUE))
  expect_match2(
    scode, "diag_pre_multiply(rho_s2z_1[j]', L_Sigma_s2z_1);"
  )
  expect_match2(scode, "+ log_det_partial_s2z_1")
  expect_null(standata(form, data = s2z_dat)$rho_s2z_1)

  independent_form <- y ~ x +
    (1 + x || gr(g, s2z = TRUE, center = "fisher"))
  independent_code <- stancode(independent_form, data = s2z_dat)
  independent_tdata <- s2z_stan_between(
    independent_code, "transformed data {", "\nparameters {"
  )
  independent_tpar <- s2z_stan_between(
    independent_code, "transformed parameters {", "\nmodel {"
  )
  expect_match2(
    independent_tdata,
    "array[N_1] matrix[M_1, M_1] gram_fisher_s2z_1;"
  )
  expect_false(grepl("J_1[n]", independent_tpar, fixed = TRUE))
  expect_match2(
    independent_tpar,
    "quad_form_diag(gram_fisher_s2z_1[j], sd_1)"
  )
  expect_match2(independent_tpar, "unit_rhs_fisher_s2z")
  expect_match2(independent_tpar, "dot_self(unit_column_fisher_s2z)")
  expect_equal(
    s2z_count_fixed(independent_tpar, "quad_form("), 0L
  )
  expect_false(grepl("diag_matrix(sd_1)", independent_tpar, fixed = TRUE))
  expect_false(grepl("identity_matrix(M_1)", independent_tpar, fixed = TRUE))
  expect_false(grepl(
    "mdivide_left_spd", independent_code, fixed = TRUE
  ))
  expect_false(grepl(
    "post_cov_fisher_s2z", independent_code, fixed = TRUE
  ))

  varying_form <- y ~ x + (1 + x || gr(
    g, s2z = TRUE, center = "fisher", scale = "varying"
  ))
  varying_code <- stancode(varying_form, data = s2z_dat)
  varying_tpar <- s2z_stan_between(
    varying_code, "transformed parameters {", "\nmodel {"
  )
  expect_match2(
    varying_tpar,
    paste0(
      "quad_form_diag(gram_fisher_s2z_1[j], ",
      "reference_sd_s2z_1)"
    )
  )
  expect_false(grepl(
    "diag_matrix(reference_sd_s2z_1)", varying_tpar, fixed = TRUE
  ))
  expect_false(grepl(
    "identity_matrix(M_1)", varying_tpar, fixed = TRUE
  ))
})

test_that("Fisher S2Z supports Bernoulli joint blocks and NB2 varying scales", {
  bern_dat <- data.frame(
    y = rep(c(0L, 1L), 8),
    x = seq(-1, 1, length.out = 16),
    g = factor(rep(seq_len(4), each = 4)),
    h = factor(rep(seq_len(8), each = 2))
  )
  bern_form <- y ~ x +
    (1 | gr(g, s2z = TRUE, center = "fisher")) +
    (1 | gr(h, s2z = TRUE, center = "fisher"))
  bern_prior <-
    prior(normal(0, 2), class = "b") +
    prior(normal(0, 2), class = "Intercept") +
    prior(exponential(1), class = "sd")
  bern_code <- stancode(
    bern_form, data = bern_dat, family = bernoulli(), prior = bern_prior
  )
  bern_tpar <- s2z_stan_between(
    bern_code, "transformed parameters {", "\nmodel {"
  )

  expect_match2(
    bern_tpar,
    paste0(
      "real eta_fisher_s2z_mu = theta_s2z[1] + ",
      "dot_product(Xc[n], tail(theta_s2z, 1));"
    )
  )
  expect_match2(
    bern_tpar,
    "real value_fisher_s2z_mu = inv_logit(eta_fisher_s2z_mu);"
  )
  expect_match2(bern_tpar, "info_fisher_s2z[J_1[n]] +=")
  expect_match2(bern_tpar, "info_fisher_s2z[J_2[n]] +=")
  expect_match2(
    bern_tpar,
    paste0(
      "real obs_prec_fisher_s2z = 1.0 * ",
      "((value_fisher_s2z_mu) * (1.0 - (value_fisher_s2z_mu)));"
    )
  )
  expect_false(grepl(
    "square(derivative_fisher_s2z_mu) /",
    bern_tpar, fixed = TRUE
  ))
  expect_match2(bern_code, "+ log_det_partial_s2z_1")
  expect_match2(bern_code, "+ log_det_partial_s2z_2")
  expect_match2(bern_code, "fast Gaussian Matheron system")
  expect_match2(bern_code, "real<lower=0> W_matheron_s2z_1;")
  expect_false(grepl("exposure_fisher_s2z", bern_code, fixed = TRUE))
  expect_null(standata(
    bern_form, data = bern_dat, family = bernoulli(), prior = bern_prior
  )$rho_s2z_1)

  # This mirrors the nested factor-slope structure in test_model2.R. Locally
  # absent surgeon columns remain exact structural zeros in the Fisher Gram.
  nested_dat <- data.frame(
    y = rep(c(0L, 1L), 8),
    x = seq(-1, 1, length.out = 16),
    hospital = factor(rep(c("a", "b"), each = 8)),
    surgeon = factor(rep(letters[1:8], each = 2))
  )
  nested_form <-
    y ~ surgeon + x +
      (1 + surgeon | gr(
        hospital, s2z = TRUE, center = "fisher"
      )) +
      (1 | gr(surgeon, s2z = TRUE, center = "fisher"))
  nested_prior <-
    prior(normal(0, 10), class = "b") +
    prior(exponential(1), class = "sd")
  nested_code <- stancode(
    nested_form, data = nested_dat, family = bernoulli(),
    prior = nested_prior
  )
  expect_match2(nested_code, "array[N_1] matrix[M_1, M_1] info_fisher_s2z")
  expect_match2(nested_code, "design_fisher_s2z[8] = Z_1_8[n];")
  expect_match2(nested_code, "+ log_det_partial_s2z_1")
  expect_match2(nested_code, "+ log_det_partial_s2z_2")
  expect_match2(nested_code, "matrix[8, 8] W_matheron_s2z_1;")
  expect_match2(nested_code, "matrix[8, 8] L_W_matheron_s2z_1;")

  nested_sdata <- standata(
    nested_form, data = nested_dat, family = bernoulli(),
    prior = nested_prior
  )
  nested_Z <- do.call(cbind, nested_sdata[paste0("Z_1_", seq_len(8))])
  nested_grams <- lapply(seq_len(nested_sdata$N_1), function(j) {
    crossprod(nested_Z[nested_sdata$J_1 == j, , drop = FALSE])
  })
  expect_true(all(vapply(
    nested_grams, function(x) any(diag(x) == 0), logical(1)
  )))

  nb_dat <- data.frame(
    y = c(2L, 8L, 3L, 7L, 12L, 6L, 10L, 5L),
    x = seq(-1, 1, length.out = 8),
    g = factor(rep(seq_len(2), each = 4))
  )
  nb_form <- bf(
    y ~ x + (1 | gr(
      g, s2z = TRUE, center = "fisher", scale = "varying"
    )),
    shape ~ (1 | gr(
      g, s2z = TRUE, center = "fisher", scale = "varying"
    ))
  )
  nb_code <- stancode(nb_form, data = nb_dat, family = negbinomial())
  nb_tpar <- s2z_stan_between(
    nb_code, "transformed parameters {", "\nmodel {"
  )
  expect_match2(nb_tpar, "real mean_fraction_fisher_s2z_nb =")
  expect_match2(nb_tpar, "real log_p0_fisher_s2z_nb =")
  expect_match2(nb_tpar, "real sparse_info_fisher_s2z_nb =")
  expect_match2(nb_tpar, "real dense_info_fisher_s2z_nb =")
  expect_match2(nb_tpar, "real log_shape_info_fisher_s2z_nb =")
  expect_false(grepl("Y[n]", nb_tpar, fixed = TRUE))
  expect_lt(
    regexpr("reference_sd_s2z_1 =", nb_tpar, fixed = TRUE)[1],
    regexpr("rho_s2z_1[j, 1] =", nb_tpar, fixed = TRUE)[1]
  )
  expect_lt(
    regexpr("rho_s2z_1[j, 1] =", nb_tpar, fixed = TRUE)[1],
    regexpr("scale_partial_s2z =", nb_tpar, fixed = TRUE)[1]
  )
  nb_data <- standata(nb_form, data = nb_dat, family = negbinomial())
  expect_null(nb_data$rho_s2z_1)
  expect_null(nb_data$rho_s2z_2)

  # This is the mixed varying/shared scale configuration of test_model.R mod3.
  nb_mixed_form <- bf(
    y ~ x + (1 | gr(
      g, s2z = TRUE, center = "fisher", scale = "varying"
    )),
    shape ~ (1 | gr(g, s2z = TRUE, center = "fisher"))
  )
  nb_mixed_code <- stancode(
    nb_mixed_form, data = nb_dat, family = negbinomial()
  )
  expect_match2(nb_mixed_code, "reference_sd_s2z_1")
  expect_false(grepl("reference_sd_s2z_2", nb_mixed_code, fixed = TRUE))
  expect_match2(nb_mixed_code, "rho_s2z_1[j, 1]")
  expect_match2(nb_mixed_code, "rho_s2z_2[j, 1]")
})

test_that("response-free Fisher rules cover representative likelihoods", {
  dat <- data.frame(
    y = c(0L, 1L, 2L, 4L, 3L, 5L, 1L, 2L, 0L, 4L, 2L, 3L),
    trials = rep(6L, 12),
    g = factor(rep(seq_len(3), each = 4)),
    category = factor(rep(c("a", "b", "c"), 4)),
    p1 = seq(0.15, 0.25, length.out = 12),
    p2 = seq(0.25, 0.35, length.out = 12)
  )
  dat$p3 <- 1 - dat$p1 - dat$p2
  dat$simplex <- I(as.matrix(dat[c("p1", "p2", "p3")]))

  check_code <- function(form, family, markers, prior = empty_prior()) {
    scode <- stancode(form, data = dat, family = family, prior = prior)
    tpar <- s2z_stan_between(
      scode, "transformed parameters {", "\nmodel {"
    )
    sdata <- standata(form, data = dat, family = family, prior = prior)
    expect_false(any(startsWith(names(sdata), "rho_s2z_")))
    expect_false(grepl("Y[n]", tpar, fixed = TRUE))
    for (marker in markers) {
      expect_match2(tpar, marker)
    }
    invisible(tpar)
  }

  # Exact conditional expected information.
  check_code(
    y ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher")),
    poisson(), "real obs_prec_fisher_s2z = value_fisher_s2z_mu;"
  )

  # Three-bin coarsening: {0}, {trials}, and the interior counts.
  check_code(
    bf(
      y | trials(trials) ~ 1 +
        (1 | gr(g, s2z = TRUE, center = "fisher")),
      phi ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher"))
    ),
    beta_binomial(),
    c("pmid_fisher_s2z_bb", "dpmid_fisher_s2z_bb")
  )

  # Structural-atom decomposition for both the count and atom predictors.
  zip_prior <-
    prior(normal(0, 2), class = "Intercept") +
    prior(normal(0, 2), class = "Intercept", dpar = "zi")
  check_code(
    bf(
      y ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher")),
      zi ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher"))
    ),
    zero_inflated_poisson(),
    c("q0_fisher_s2z_zi", "atom_derivative_fisher_s2z_zi"),
    prior = zip_prior
  )

  # Marginal category information from the population-only softmax.
  check_code(
    category ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher")),
    categorical(), "prob_fisher_s2z_cat"
  )

  # Simplex-coordinate information from the Dirichlet trigamma identity.
  check_code(
    simplex ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher")),
    dirichlet(),
    c(
      "prob_fisher_s2z_dir",
      "alpha_fisher_s2z_dir < 1e-6 ? 1.0 :",
      "square(alpha_fisher_s2z_dir) * trigamma(alpha_fisher_s2z_dir)"
    )
  )

  # COM-Poisson moment information; at shape one the location variance
  # marker reduces to the Poisson mean information.
  check_code(
    bf(
      y ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher")),
      shape ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher"))
    ),
    brmsfamily("com_poisson"),
    c(
      paste0(
        "log(fmax(value_fisher_s2z_mu, 1e-12)) / ",
        "fmax(value_fisher_s2z_shape, 1e-12);"
      ),
      paste0(
        "mode_fisher_s2z_cmp / ",
        "fmax(value_fisher_s2z_shape, 1e-12));"
      ),
      "real log_factorial_variance_fisher_s2z_cmp =",
      "real obs_prec_fisher_s2z = variance_fisher_s2z_cmp;",
      paste0(
        "square(derivative_fisher_s2z_shape) * ",
        "log_factorial_variance_fisher_s2z_cmp"
      )
    )
  )
})

test_that("response-free Fisher rules guard zero-trial observations", {
  dat <- data.frame(
    y = c(0L, 0L, 1L, 2L, 0L, 4L),
    trials = c(0L, 1L, 2L, 4L, 0L, 5L),
    g = factor(rep(seq_len(2), each = 3))
  )
  beta_binomial_form <- bf(
    y | trials(trials) ~ 1 +
      (1 | gr(g, s2z = TRUE, center = "fisher")),
    phi ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher"))
  )
  beta_binomial_code <- stancode(
    beta_binomial_form, data = dat, family = beta_binomial()
  )
  beta_binomial_tpar <- s2z_stan_between(
    beta_binomial_code, "transformed parameters {", "\nmodel {"
  )
  expect_match2(
    beta_binomial_tpar,
    "(trials[n] == 0 ? 0.0 : (trials[n] == 1 ?"
  )
  expect_match2(
    beta_binomial_tpar,
    "(trials[n] <= 1 ? 0.0 :"
  )
  expect_match2(
    beta_binomial_tpar,
    paste0(
      "((value_fisher_s2z_mu) * ",
      "(1.0 - (value_fisher_s2z_mu)))"
    )
  )
  expect_match2(
    beta_binomial_tpar,
    paste0(
      "prob_fisher_s2z_bb = fmin(1.0 - 1e-12, ",
      "fmax(1e-12, value_fisher_s2z_mu));"
    )
  )
  beta_binomial_data <- standata(
    beta_binomial_form, data = dat, family = beta_binomial()
  )
  expect_equal(as.vector(beta_binomial_data$trials), dat$trials)
  expect_false(any(startsWith(
    names(beta_binomial_data), "rho_s2z_"
  )))
  expect_false(grepl("Y[n]", beta_binomial_tpar, fixed = TRUE))

  dm_dat <- data.frame(
    trials = c(0L, 3L, 4L, 0L, 5L, 2L),
    g = factor(rep(seq_len(2), each = 3))
  )
  dm_dat$counts <- I(rbind(
    c(0L, 0L, 0L), c(1L, 1L, 1L), c(1L, 2L, 1L),
    c(0L, 0L, 0L), c(2L, 1L, 2L), c(1L, 0L, 1L)
  ))
  dm_form <- bf(
    counts | trials(trials) ~ 1 +
      (1 | gr(g, s2z = TRUE, center = "fisher")),
    phi ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher"))
  )
  dm_code <- stancode(
    dm_form, data = dm_dat, family = dirichlet_multinomial()
  )
  dm_tpar <- s2z_stan_between(
    dm_code, "transformed parameters {", "\nmodel {"
  )
  expect_match2(
    dm_tpar,
    paste0(
      "((trials[n]) * (1.0 + (value_fisher_s2z_phi)) / ",
      "((trials[n]) + (value_fisher_s2z_phi)))"
    )
  )
  expect_match2(dm_tpar, "(trials[n] == 0 ? 0.0 : 0.5 *")
  dm_data <- standata(
    dm_form, data = dm_dat, family = dirichlet_multinomial()
  )
  expect_equal(as.vector(dm_data$trials), dm_dat$trials)
  expect_false(any(startsWith(names(dm_data), "rho_s2z_")))
  expect_false(grepl("Y[n]", dm_tpar, fixed = TRUE))
})

test_that("Fisher link algebra is stable at parameter boundaries", {
  dat <- data.frame(
    y_unit = seq(0.1, 0.9, length.out = 8),
    y_pos = seq(0.3, 2, length.out = 8),
    y_bin = rep(0:1, 4),
    g = factor(rep(seq_len(2), each = 4))
  )

  xbeta_code <- stancode(
    y_unit ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher")),
    data = dat, family = xbeta()
  )
  for (marker in c(
    paste0(
      "prob_fisher_s2z_xbeta = fmin(1.0 - 1e-12, ",
      "fmax(1e-12, value_fisher_s2z_mu));"
    ),
    paste0(
      "(1.0 + (value_fisher_s2z_phi)) * (value_fisher_s2z_mu) * ",
      "(1.0 - (value_fisher_s2z_mu))"
    )
  )) {
    expect_match2(xbeta_code, marker)
  }

  frechet_code <- stancode(
    bf(
      y_pos ~ 1,
      nu ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher"))
    ),
    data = dat, family = frechet()
  )
  for (marker in c(
    "real boundary_fraction_fisher_s2z_frechet =",
    "boundary_fraction_fisher_s2z_frechet < 1e-6 ?",
    "real shape_info_fisher_s2z_frechet =",
    "real obs_prec_fisher_s2z = shape_info_fisher_s2z_frechet;"
  )) {
    expect_match2(frechet_code, marker)
  }

  softit_code <- stancode(
    y_bin ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher")),
    data = dat, family = bernoulli("softit")
  )
  expect_match2(softit_code, "return log(expm1(-p ./ (p - 1)));")
  expect_match2(softit_code, "return log1p_exp(y) ./ (1 + log1p_exp(y));")
})

test_that("Wiener Fisher rules use exact decision coarsening", {
  dat <- data.frame(
    q = seq(0.5, 1.6, length.out = 12),
    decision = rep(0:1, 6),
    g = factor(rep(seq_len(3), each = 4))
  )
  form <- bf(
    q | dec(decision) ~ 1,
    bs ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher")),
    bias ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher")),
    ndt = 0.1
  )
  s2z_prior <-
    prior(normal(0, 2), class = "Intercept", dpar = "bs") +
    prior(normal(0, 2), class = "Intercept", dpar = "bias")
  scode <- stancode(
    form, data = dat, family = wiener(), prior = s2z_prior
  )
  tpar <- s2z_stan_between(
    scode, "transformed parameters {", "\nmodel {"
  )
  for (marker in c(
    "real choice_scale_fisher_s2z_wiener =",
    "if (abs(choice_scale_fisher_s2z_wiener) < 1e-5)",
    "real dp_dscale_fisher_s2z_wiener;",
    "real dp_dbias_fisher_s2z_wiener;",
    paste0(
      "prob_safe_fisher_s2z_wiener = fmin(1.0 - 1e-12, ",
      "fmax(1e-12, p_upper_fisher_s2z_wiener));"
    ),
    paste0(
      "square((derivative_fisher_s2z_bs) * 2.0 * ",
      "(value_fisher_s2z_mu) * dp_dscale_fisher_s2z_wiener)"
    ),
    paste0(
      "square((derivative_fisher_s2z_bias) * ",
      "dp_dbias_fisher_s2z_wiener)"
    )
  )) {
    expect_match2(tpar, marker)
  }
  expect_false(grepl("Y[n]", tpar, fixed = TRUE))
  expect_false(grepl("dec[n]", tpar, fixed = TRUE))
  sdata <- standata(
    form, data = dat, family = wiener(), prior = s2z_prior
  )
  expect_false(any(startsWith(names(sdata), "rho_s2z_")))

  expect_error(
    stancode(
      bf(
        q | dec(decision) ~ 1,
        ndt ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher"))
      ),
      data = dat, family = wiener()
    ),
    paste0(
      "has no response-free expected-information rule for family ",
      "'wiener' and distributional parameter 'ndt'"
    ),
    fixed = TRUE
  )
})

test_that("sampled-loading Fisher S2Z information remains dynamic", {
  dat <- expand.grid(
    person = factor(seq_len(8)),
    item = factor(letters[1:4])
  )
  dat$y <- seq_len(nrow(dat)) / 10
  form <- bf(
    y ~ alpha + loading * eta,
    alpha ~ 0 + item,
    loading ~ 0 + item,
    eta ~ 0 + (1 | gr(
      person, s2z = TRUE, latent = TRUE, center = "fisher"
    )),
    nl = TRUE
  )
  scode <- stancode(form, data = dat)
  tdata <- s2z_stan_between(
    scode, "transformed data {", "\nparameters {"
  )
  tpar <- s2z_stan_between(
    scode, "transformed parameters {", "\nmodel {"
  )

  expect_false(grepl("fisher_s2z_1_nlp_loading", tdata, fixed = TRUE))
  expect_false(grepl("info_fisher_s2z", tdata, fixed = TRUE))
  expect_false(grepl("gram_fisher_s2z_1", scode, fixed = TRUE))
  expect_false(grepl("exposure_fisher_s2z_1", scode, fixed = TRUE))
  expect_match2(tpar, "vector[N] fisher_s2z_1_nlp_loading;")
  expect_match2(
    tpar, "fisher_s2z_1_nlp_loading = X_loading * b_loading;"
  )
  expect_match2(tpar, "for (n in 1:N)")
  expect_match2(tpar, "fisher_s2z_1_nlp_loading[n]")
  expect_match2(tpar, "J_1[n]")
  expect_match2(tpar, "vector[N_1] info_fisher_s2z = zeros_vector(N_1);")
  expect_match2(tpar, "info_fisher_s2z[J_1[n]] +=")
  expect_match2(tpar, "rho_s2z_1[j, 1]")
  expect_match2(scode, "+ log_det_partial_s2z_1")
  expect_null(standata(form, data = dat)$rho_s2z_1)
})

test_that("strict latent Fisher uses every scalar native location rule", {
  N <- 12L
  dat <- data.frame(
    person = factor(rep(letters[1:4], each = 3L)),
    item = factor(rep(letters[1:3], 4L)),
    y_real = seq(-1.2, 1.4, length.out = N),
    y_pos = seq(0.5, 3, length.out = N),
    y_unit = seq(0.05, 0.95, length.out = N),
    y_count = rep(0:5, length.out = N),
    y_bin = rep(0:1, length.out = N),
    trials = rep(6L, N),
    dec = rep(0:1, length.out = N)
  )
  cases <- list(
    gaussian = list(gaussian(), "y_real"),
    student = list(student(), "y_real"),
    skew_normal = list(skew_normal(), "y_real"),
    bernoulli = list(bernoulli(), "y_bin"),
    binomial = list(binomial(), "y_bin | trials(trials)"),
    xbeta = list(xbeta(), "y_unit"),
    beta_binomial = list(
      beta_binomial(), "y_count | trials(trials)"
    ),
    beta = list(Beta(), "y_unit"),
    poisson = list(poisson(), "y_count"),
    negbinomial = list(negbinomial(), "y_count"),
    negbinomial2 = list(negbinomial2(), "y_count"),
    geometric = list(geometric(), "y_count"),
    discrete_weibull = list(discrete_weibull(), "y_count"),
    com_poisson = list(com_poisson(), "y_count"),
    gamma = list(Gamma(), "y_pos"),
    weibull = list(weibull(), "y_pos"),
    exponential = list(exponential(), "y_pos"),
    frechet = list(frechet(), "y_pos"),
    inverse_gaussian = list(brmsfamily("inverse.gaussian"), "y_pos"),
    lognormal = list(lognormal(), "y_pos"),
    shifted_lognormal = list(shifted_lognormal(), "y_pos"),
    exgaussian = list(exgaussian(), "y_real"),
    von_mises = list(von_mises(), "y_real"),
    asym_laplace = list(asym_laplace(), "y_real"),
    hurdle_poisson = list(hurdle_poisson(), "y_count"),
    hurdle_negbinomial = list(hurdle_negbinomial(), "y_count"),
    hurdle_gamma = list(hurdle_gamma(), "y_pos"),
    hurdle_lognormal = list(hurdle_lognormal(), "y_pos"),
    zero_inflated_poisson = list(zero_inflated_poisson(), "y_count"),
    zero_inflated_negbinomial = list(
      zero_inflated_negbinomial(), "y_count"
    ),
    zero_inflated_binomial = list(
      zero_inflated_binomial(), "y_count | trials(trials)"
    ),
    zero_inflated_beta_binomial = list(
      zero_inflated_beta_binomial(), "y_count | trials(trials)"
    ),
    zero_inflated_beta = list(zero_inflated_beta(), "y_unit"),
    zero_one_inflated_beta = list(zero_one_inflated_beta(), "y_unit"),
    zero_inflated_asym_laplace = list(
      zero_inflated_asym_laplace(), "y_real"
    ),
    cox = list(cox(), "y_pos"),
    wiener = list(wiener(), "y_pos | dec(dec)")
  )

  for (name in names(cases)) {
    case <- cases[[name]]
    response <- as.formula(
      paste(case[[2L]], "~ baseline + loading * eta")
    )
    form <- bf(
      response,
      baseline ~ 0 + item,
      loading ~ 0 + item,
      eta ~ 0 + (1 | gr(
        person, id = "score", s2z = TRUE, latent = TRUE,
        center = "fisher"
      )),
      nl = TRUE,
      family = case[[1L]]
    )
    scode <- stancode(form, data = dat)
    tpar <- s2z_stan_between(
      scode, "transformed parameters {", "\nmodel {"
    )
    expect_match2(tpar, "info_fisher_s2z")
    expect_match2(tpar, "fisher_s2z_1_nlp_loading")
    expect_false(grepl("Y[n]", tpar, fixed = TRUE), info = name)
  }
  expect_length(cases, 37L)
})

test_that("Fisher S2Z rejects nonlocal likelihood structures", {
  expect_error(
    stancode(
      y | weights(w) ~ 1 +
        (1 | gr(g, s2z = TRUE, center = "fisher")),
      data = s2z_dat
    ),
    "does not yet support response addition term 'weights'",
    fixed = TRUE
  )

  latent_data <- expand.grid(
    person = factor(seq_len(5)), item = factor(seq_len(3))
  )
  latent_data$y <- seq_len(nrow(latent_data)) / 10
  latent_formula <- bf(
    y ~ alpha + loading * eta,
    alpha ~ 0 + item,
    loading ~ 0 + item,
    eta ~ 1 +
      (1 | gr(person, s2z = TRUE, center = "fisher")),
    nl = TRUE
  )
  latent_prior <- c(
    prior(normal(0, 1), nlpar = "alpha"),
    prior(normal(0, 1), nlpar = "loading"),
    prior(normal(0, 1), nlpar = "eta")
  )
  expect_error(
    stancode(latent_formula, data = latent_data, prior = latent_prior),
    paste0(
      "does not yet support nonlinear predictors; the derivative of the ",
      "response mean with respect to the latent score must be represented ",
      "explicitly"
    ),
    fixed = TRUE
  )

  mv_data <- transform(s2z_dat, y2 = y + 0.2 * x)
  mv_form <- bf(
    y ~ 1 + (1 | gr(g, id = "y_block", s2z = TRUE,
                    center = "fisher"))
  ) +
    bf(
      y2 ~ 1 + (1 | gr(g, id = "y2_block", s2z = TRUE,
                       center = "fisher"))
    ) +
    set_rescor(FALSE)
  mv_code <- stancode(mv_form, data = mv_data)
  expect_match2(mv_code, "rho_s2z_1[j, 1]")
  expect_match2(mv_code, "rho_s2z_2[j, 1]")

  expect_error(
    stancode(mv_form + set_rescor(TRUE), data = mv_data),
    "currently requires set_rescor(FALSE)",
    fixed = TRUE
  )
})

test_that("explicit centered S2Z preserves default code for every kernel", {
  cases <- list(
    scalar = list(
      data = s2z_dat,
      default = y ~ 1 + (1 | gr(g, s2z = TRUE)),
      explicit = y ~ 1 +
        (1 | gr(g, s2z = TRUE, center = TRUE))
    ),
    independent_four = list(
      data = s2z_dat,
      default = y ~ x * z + (1 + x * z || gr(g, s2z = TRUE)),
      explicit = y ~ x * z +
        (1 + x * z || gr(g, s2z = TRUE, center = TRUE))
    ),
    independent_ten = list(
      data = s2z_ten_dat,
      default = y ~ 0 + ten + (0 + ten || gr(g, s2z = TRUE)),
      explicit = y ~ 0 + ten +
        (0 + ten || gr(g, s2z = TRUE, center = TRUE))
    ),
    correlated = list(
      data = s2z_dat,
      default = y ~ x * z + (1 + x * z | gr(g, s2z = TRUE)),
      explicit = y ~ x * z +
        (1 + x * z | gr(g, s2z = TRUE, center = TRUE))
    ),
    student = list(
      data = s2z_dat,
      default = y ~ x +
        (1 + x | gr(g, s2z = TRUE, dist = "student")),
      explicit = y ~ x +
        (1 + x | gr(
          g, s2z = TRUE, center = TRUE, dist = "student"
        ))
    ),
    varying_correlated = list(
      data = s2z_dat,
      default = y ~ x * z +
        (1 + x * z | gr(g, s2z = TRUE, scale = "varying")),
      explicit = y ~ x * z +
        (1 + x * z | gr(
          g, s2z = TRUE, center = TRUE, scale = "varying"
        ))
    ),
    varying_independent = list(
      data = s2z_ten_dat,
      default = y ~ 0 + ten +
        (0 + ten || gr(g, s2z = TRUE, scale = "varying")),
      explicit = y ~ 0 + ten +
        (0 + ten || gr(
          g, s2z = TRUE, center = TRUE, scale = "varying"
        ))
    )
  )

  for (case in cases) {
    expect_identical(
      stancode(case$explicit, data = case$data),
      stancode(case$default, data = case$data)
    )
  }
  expect_identical(
    standata(cases$correlated$default, data = s2z_dat),
    standata(
      y ~ x * z +
        (1 + x * z | gr(g, s2z = TRUE, center = FALSE)),
      data = s2z_dat
    )
  )
})

test_that("partial S2Z code covers every specialized covariance path", {
  scalar <- stancode(
    y ~ 1 + (1 | gr(g, s2z = TRUE, center = 0.35)),
    data = s2z_dat,
    prior = prior(normal(0, 2), class = Intercept)
  )
  scalar_data <- standata(
    y ~ 1 + (1 | gr(g, s2z = TRUE, center = 0.35)),
    data = s2z_dat
  )
  expect_equal(
    unname(scalar_data$rho_s2z_1),
    matrix(0.35, nrow = 6L, ncol = 1L)
  )
  for (term in c(
    "matrix<lower=0,upper=1>[N_1, M_1] rho_s2z_1;",
    "vector[M_1] mean_rho_s2z_1;",
    "mean_rho_s2z_1[k] = mean(rho_s2z_1[, k]);",
    "real log_det_partial_s2z_1;",
    paste0(
      "vector[N_1] scale_partial_s2z = 1.0 - rho_s2z_1[, 1] + ",
      "rho_s2z_1[, 1] * sd_1[1];"
    ),
    paste0(
      "centered_partial_s2z = sd_1[1] * centered_partial_s2z ./ ",
      "scale_partial_s2z;"
    ),
    "r_s2z_1_1 = centered_partial_s2z - mean(centered_partial_s2z);",
    "log_det_partial_s2z_1 += -sum(log(scale_partial_s2z));",
    "log_det_partial_s2z_1 += log(",
    paste0(
      "(1.0 - mean_rho_s2z_1[1]) + ",
      "mean_rho_s2z_1[1] * sd_1[1]"
    ),
    "+ log_det_partial_s2z_1"
  )) {
    expect_true(grepl(term, scalar, fixed = TRUE), info = term)
  }
  expect_false(grepl("matrix[M_1, M_1]", scalar, fixed = TRUE))

  independent_ten <- stancode(
    y ~ 0 + ten +
      (0 + ten || gr(g, s2z = TRUE, center = 0.62)),
    data = s2z_ten_dat,
    prior = prior(normal(0, 2), class = b)
  )
  independent_ten_data <- standata(
    y ~ 0 + ten +
      (0 + ten || gr(g, s2z = TRUE, center = 0.62)),
    data = s2z_ten_dat
  )
  expect_equal(
    unname(independent_ten_data$rho_s2z_1),
    matrix(0.62, nrow = 8L, ncol = 10L)
  )
  for (k in seq_len(10L)) {
    expect_match2(
      independent_ten,
      sprintf(
        paste0(
          "rho_s2z_1[, %1$s] + rho_s2z_1[, %1$s] * sd_1[%1$s]"
        ),
        k
      )
    )
    expect_match2(
      independent_ten,
      sprintf(
        "r_s2z_1_%s = centered_partial_s2z - mean(centered_partial_s2z)",
        k
      )
    )
  }
  expect_match2(independent_ten, "+ log_det_partial_s2z_1")
  expect_false(grepl(
    "matrix[M_1, M_1] L_partial_s2z", independent_ten, fixed = TRUE
  ))
  expect_false(grepl("cholesky_decompose(", independent_ten, fixed = TRUE))

  correlated <- stancode(
    y ~ x * z +
      (1 + x * z | gr(g, s2z = TRUE, center = 0.44)),
    data = s2z_dat
  )
  for (term in c(
    paste0(
      "matrix[M_1, M_1] L_partial_s2z = ",
      "diag_pre_multiply(rho_s2z_1[j]', L_Sigma_s2z_1);"
    ),
    "L_partial_s2z[k, k] += 1.0 - rho_s2z_1[j, k];",
    paste0(
      "r_s2z_1[j] = (L_Sigma_s2z_1 * ",
      "mdivide_left_tri_low(L_partial_s2z, r_s2z_1[j]'))';"
    ),
    "mean_partial_s2z += r_s2z_1[j];",
    "mean_partial_s2z /= N_1;",
    "r_s2z_1[j] -= mean_partial_s2z;",
    "log_det_partial_s2z_1 -= sum(log(diagonal(L_partial_s2z)));",
    "log_det_partial_s2z_1 += log(",
    paste0(
      "(1.0 - mean_rho_s2z_1[k]) + ",
      "mean_rho_s2z_1[k] * L_Sigma_s2z_1[k, k]"
    ),
    "+ log_det_partial_s2z_1"
  )) {
    expect_true(grepl(term, correlated, fixed = TRUE), info = term)
  }

  student <- stancode(
    y ~ x + (1 + x | gr(
      g, s2z = TRUE, center = 0.44, dist = "student"
    )),
    data = s2z_dat
  )
  for (term in c(
    "group_prec_s2z_1 = inv_square(group_scale_s2z_1);",
    "- M_1 * sum(log(group_scale_s2z_1))",
    "+ log_det_partial_s2z_1"
  )) {
    expect_true(grepl(term, student, fixed = TRUE), info = term)
  }

  varying <- stancode(
    y ~ 0 + ten + (0 + ten || gr(
      g, s2z = TRUE, center = 0.62, scale = "varying"
    )),
    data = s2z_ten_dat,
    prior = prior(normal(0, 2), class = b)
  )
  for (term in c(
    "rho_s2z_1[, 10] * reference_sd_s2z_1[10]",
    paste0(
      "sd_level_s2z_1[, k] = reference_sd_s2z_1[k] * ",
      "exp(sdlog_1[k] * z_sd_centered_s2z);"
    ),
    "+ log_det_partial_s2z_1"
  )) {
    expect_true(grepl(term, varying, fixed = TRUE), info = term)
  }
  expect_false(grepl(
    "matrix[M_1, M_1] L_partial_s2z", varying, fixed = TRUE
  ))

  unnormalized <- stancode(
    y ~ 1 + (1 | gr(g, s2z = TRUE, center = 0.35)),
    data = s2z_dat, normalize = FALSE,
    prior = prior(normal(0, 2), class = Intercept)
  )
  expect_match2(unnormalized, "+ log_det_partial_s2z_1")
  expect_false(grepl("log(2 * pi())", unnormalized, fixed = TRUE))
})

test_that("partial S2Z log determinants are stable at centering endpoints", {
  rho <- c(0, 1 - 1e-8, 1)
  scale <- rep(1e-20, length(rho))
  stable_log_term <- vapply(seq_along(rho), function(i) {
    if (rho[i] == 1) {
      log(scale[i])
    } else {
      log1p(-rho[i] * (1 - scale[i]))
    }
  }, numeric(1))
  reference <- log((1 - rho) + rho * scale)

  expect_equal(stable_log_term, reference, tolerance = 1e-10)
  expect_true(all(is.finite(stable_log_term)))
  expect_identical(log1p(-(1 - 1e-20)), -Inf)
})

test_that("partial varying-scale Student kernels retain every measure term", {
  correlated_form <- y ~ x + (1 + x | gr(
    g, s2z = TRUE, center = 0.43, scale = "varying",
    dist = "student"
  ))
  scalar_form <- y ~ 1 + (1 | gr(
    g, s2z = TRUE, center = 0.43, scale = "varying",
    dist = "student"
  ))

  for (normalize in c(TRUE, FALSE)) {
    correlated <- stancode(
      correlated_form, data = s2z_dat, normalize = normalize
    )
    scalar <- stancode(
      scalar_form, data = s2z_dat, normalize = normalize
    )

    for (term in c(
      "matrix<lower=0,upper=1>[N_1, M_1] rho_s2z_1;",
      paste0(
        "reference_sd_s2z_1 = sd_1 .* exp(sdlog_1 .* ",
        "tail(z_sd_s2z_1, M_1) / sqrt(1.0 * N_1));"
      ),
      paste0(
        "sd_level_s2z_1[, k] = reference_sd_s2z_1[k] * ",
        "exp(sdlog_1[k] * z_sd_centered_s2z);"
      ),
      "group_prec_s2z_1 = inv_square(group_scale_s2z_1);",
      "+ log_det_partial_s2z_1",
      "- M_1 * sum(log(group_scale_s2z_1))"
    )) {
      expect_true(grepl(term, correlated, fixed = TRUE), info = term)
      expect_true(grepl(term, scalar, fixed = TRUE), info = term)
    }

    for (term in c(
      paste0(
        "L_Sigma_s2z_1 = diag_pre_multiply(",
        "reference_sd_s2z_1, L_1);"
      ),
      paste0(
        "matrix[M_1, M_1] L_partial_s2z = ",
        "diag_pre_multiply(rho_s2z_1[j]', L_Sigma_s2z_1);"
      ),
      paste0(
        "matrix[M_1, M_1] relative_precision_s2z = ",
        "mdivide_left_tri_low(L_level_s2z, L_Sigma_s2z_1);"
      ),
      paste0(
        "P_s2z_1 += group_prec_s2z_1[j] * ",
        "crossprod(relative_precision_s2z);"
      ),
      paste0(
        "group_quad_s2z_1 += group_prec_s2z_1[j] * ",
        "dot_self(white_level_s2z);"
      )
    )) {
      expect_true(grepl(term, correlated, fixed = TRUE), info = term)
    }
    expect_false(grepl(
      paste0(
        "- (N_1 - 1) * ",
        "sum(log(diagonal(L_Sigma_s2z_1)))"
      ),
      correlated, fixed = TRUE
    ))

    for (term in c(
      paste0(
        "scale_partial_s2z = 1.0 - ",
        "rho_s2z_1[, 1] + rho_s2z_1[, 1] * ",
        "reference_sd_s2z_1[1];"
      ),
      paste0(
        "centered_partial_s2z = reference_sd_s2z_1[1] * ",
        "centered_partial_s2z ./ scale_partial_s2z;"
      ),
      paste0(
        "real weighted_precision_s2z = group_prec_s2z_1[n] * ",
        "square(relative_precision_s2z);"
      ),
      paste0(
        "group_quad_s2z_1 += group_prec_s2z_1[n] * ",
        "square(white_level_s2z);"
      )
    )) {
      expect_true(grepl(term, scalar, fixed = TRUE), info = term)
    }
    expect_false(grepl(
      "- (N_1 - 1) * sum(log(reference_sd_s2z_1))",
      scalar, fixed = TRUE
    ))
    expect_false(grepl(
      "matrix[M_1, M_1] L_partial_s2z", scalar, fixed = TRUE
    ))
    expect_false(grepl("cholesky_decompose(", scalar, fixed = TRUE))

    if (normalize) {
      expect_match2(correlated, "+ 0.5 * M_1 * log(1.0 * N_1)")
      expect_match2(scalar, "+ 0.5 * M_1 * log(1.0 * N_1)")
    } else {
      expect_false(grepl("log(1.0 * N_1)", correlated, fixed = TRUE))
      expect_false(grepl("log(1.0 * N_1)", scalar, fixed = TRUE))
      expect_false(grepl("log(2 * pi())", correlated, fixed = TRUE))
      expect_false(grepl("log(2 * pi())", scalar, fixed = TRUE))
    }
  }
})

test_that("direct non-centered scalar S2Z scales contrasts exactly", {
  centered <- stancode(
    y ~ 1 + (1 | gr(g, s2z = TRUE)), data = s2z_dat,
    prior = prior(normal(0, 2), class = Intercept)
  )
  direct <- stancode(
    y ~ 1 + (1 | gr(g, s2z = TRUE, center = FALSE)),
    data = s2z_dat,
    prior = prior(normal(0, 2), class = Intercept)
  )

  expect_match2(
    centered,
    "r_s2z_1_1 = sum_to_zero_constrain_brms(z_s2z_1);"
  )
  expect_match2(
    direct,
    paste0(
      "r_s2z_1_1 = sum_to_zero_constrain_brms(",
      "sd_1[1] * z_s2z_1);"
    )
  )
  expect_match2(centered, "- (N_1 - 1) * log(sd_1[1])")
  expect_false(grepl(
    "- (N_1 - 1) * log(sd_1[1])", direct, fixed = TRUE
  ))
  for (term in c(
    "-0.5 * group_quad_s2z_1",
    "- 0.5 * log(D_s2z_1)",
    "+ 0.5 * log(1.0 * N_1)",
    "mean_r_s2z_1 = mhat_s2z_1 + sd_1[1] * std_normal_rng()",
    "q_recovered_s2z_1 = theta_s2z - H_s2z_1 * mean_r_s2z_1",
    "r_1_1 = r_s2z_1_1 + mean_r_s2z_1"
  )) {
    expect_true(grepl(term, direct, fixed = TRUE), info = term)
  }

  student <- stancode(
    y ~ 1 + (1 | gr(
      g, s2z = TRUE, center = FALSE, dist = "student"
    )),
    data = s2z_dat,
    prior = prior(normal(0, 2), class = Intercept)
  )
  expect_match2(
    student,
    paste0(
      "r_s2z_1_1 = sum_to_zero_constrain_brms(",
      "sd_1[1] * z_s2z_1);"
    )
  )
  expect_false(grepl(
    "- (N_1 - 1) * log(sd_1[1])", student, fixed = TRUE
  ))
  for (term in c(
    "dot_product(r_s2z_1_1, group_prec_s2z_1)",
    "- sum(log(group_scale_s2z_1))",
    "- 0.5 * log(D_s2z_1)",
    "+ 0.5 * log(1.0 * N_1)"
  )) {
    expect_true(grepl(term, student, fixed = TRUE), info = term)
  }
})

test_that("direct independent S2Z scales K4 and K10 component-wise", {
  four_centered <- stancode(
    y ~ x * z + (1 + x * z || gr(g, s2z = TRUE)),
    data = s2z_dat
  )
  four_direct <- stancode(
    y ~ x * z +
      (1 + x * z || gr(g, s2z = TRUE, center = FALSE)),
    data = s2z_dat
  )
  ten_direct <- stancode(
    y ~ 0 + ten +
      (0 + ten || gr(g, s2z = TRUE, center = FALSE)),
    data = s2z_ten_dat,
    prior = prior(normal(0, 2), class = b)
  )

  expect_match2(four_centered, "- (N_1 - 1) * sum(log(sd_1))")
  expect_false(grepl(
    "- (N_1 - 1) * sum(log(sd_1))", four_direct, fixed = TRUE
  ))
  expect_false(grepl(
    "- (N_1 - 1) * sum(log(sd_1))", ten_direct, fixed = TRUE
  ))
  for (code_and_k in list(c(four_direct, 4L), c(ten_direct, 10L))) {
    code <- code_and_k[[1]]
    k <- as.integer(code_and_k[[2]])
    for (j in seq_len(k)) {
      expect_match2(
        code,
        sprintf(
          paste0(
            "r_s2z_1_%1$s = sum_to_zero_constrain_brms(",
            "sd_1[%1$s] * segment(z_s2z_1, ",
            "(%1$s - 1) * (N_1 - 1) + 1, N_1 - 1));"
          ),
          j
        )
      )
    }
    for (term in c(
      "- 0.5 * sum(log(D_diag_s2z_1))",
      "- 0.5 * log1p(rank1_info_s2z_1)",
      "q_recovered_s2z_1", "mean_r_s2z_1"
    )) {
      expect_true(grepl(term, code, fixed = TRUE), info = term)
    }
    expect_false(grepl("matrix[M_1, M_1]", code, fixed = TRUE))
    expect_false(grepl("cholesky_decompose(", code, fixed = TRUE))
  }

  student_direct <- stancode(
    y ~ x * z + (1 + x * z || gr(
      g, s2z = TRUE, center = FALSE, dist = "student"
    )),
    data = s2z_dat
  )
  expect_false(grepl(
    "- (N_1 - 1) * sum(log(sd_1))", student_direct, fixed = TRUE
  ))
  for (term in c(
    paste0(
      "sd_1[4] * segment(z_s2z_1, ",
      "(4 - 1) * (N_1 - 1) + 1, N_1 - 1)"
    ),
    "- M_1 * sum(log(group_scale_s2z_1))",
    "- 0.5 * sum(log(D_diag_s2z_1))",
    "- 0.5 * log1p(rank1_info_s2z_1)"
  )) {
    expect_true(grepl(term, student_direct, fixed = TRUE), info = term)
  }
})

test_that("direct correlated S2Z applies the reference Cholesky", {
  centered <- stancode(
    y ~ x * z + (1 + x * z | gr(g, s2z = TRUE)),
    data = s2z_dat
  )
  direct <- stancode(
    y ~ x * z +
      (1 + x * z | gr(g, s2z = TRUE, center = FALSE)),
    data = s2z_dat
  )

  expect_match2(
    centered,
    paste0(
      "r_s2z_1[, k] = sum_to_zero_constrain_brms(segment(z_s2z_1, ",
      "(k - 1) * (N_1 - 1) + 1, N_1 - 1));"
    )
  )
  expect_false(grepl(
    "r_s2z_1 = r_s2z_1 * L_Sigma_s2z_1';",
    centered, fixed = TRUE
  ))
  expect_match2(
    direct, "r_s2z_1 = r_s2z_1 * L_Sigma_s2z_1';"
  )
  expect_match2(
    centered,
    "- (N_1 - 1) * sum(log(diagonal(L_Sigma_s2z_1)))"
  )
  expect_false(grepl(
    "- (N_1 - 1) * sum(log(diagonal(L_Sigma_s2z_1)))",
    direct, fixed = TRUE
  ))
  for (term in c(
    "-0.5 * group_quad_s2z_1",
    "- sum(log(diagonal(L_P_s2z_1)))",
    "+ 0.5 * M_1 * log(1.0 * N_1)",
    "q_recovered_s2z_1 = theta_s2z - H_s2z_1 * mean_r_s2z_1",
    "r_1 = r_s2z_1;",
    "for (j in 1:N_1) r_1[j] += mean_r_s2z_1';"
  )) {
    expect_true(grepl(term, direct, fixed = TRUE), info = term)
  }

  student <- stancode(
    y ~ x + (1 + x | gr(
      g, s2z = TRUE, center = FALSE, dist = "student"
    )),
    data = s2z_dat
  )
  expect_match2(student, "r_s2z_1 = r_s2z_1 * L_Sigma_s2z_1';")
  expect_false(grepl(
    "- (N_1 - 1) * sum(log(diagonal(L_Sigma_s2z_1)))",
    student, fixed = TRUE
  ))
  for (term in c(
    "contrast_score_s2z_1 = white_s2z * group_prec_s2z_1",
    "- M_1 * sum(log(group_scale_s2z_1))",
    "- sum(log(diagonal(L_P_s2z_1)))",
    "+ 0.5 * M_1 * log(1.0 * N_1)"
  )) {
    expect_true(grepl(term, student, fixed = TRUE), info = term)
  }
})

test_that("direct varying-scale S2Z cancels only its reference determinant", {
  scalar_centered <- stancode(
    y ~ 1 +
      (1 | gr(g, s2z = TRUE, scale = "varying")),
    data = s2z_dat
  )
  scalar_direct <- stancode(
    y ~ 1 +
      (1 | gr(
        g, s2z = TRUE, scale = "varying", center = FALSE
      )),
    data = s2z_dat
  )
  expect_match2(
    scalar_direct,
    paste0(
      "r_s2z_1_1 = sum_to_zero_constrain_brms(",
      "reference_sd_s2z_1[1] * segment(z_s2z_1, ",
      "(1 - 1) * (N_1 - 1) + 1, N_1 - 1));"
    )
  )
  expect_match2(
    scalar_centered,
    "- (N_1 - 1) * sum(log(reference_sd_s2z_1))"
  )
  expect_false(grepl(
    "- (N_1 - 1) * sum(log(reference_sd_s2z_1))",
    scalar_direct, fixed = TRUE
  ))
  for (term in c(
    paste0(
      "sd_level_s2z_1[, k] = reference_sd_s2z_1[k] * ",
      "exp(sdlog_1[k] * z_sd_centered_s2z);"
    ),
    "- 0.5 * sum(log(D_diag_s2z_1))",
    "- 0.5 * log1p(rank1_info_s2z_1)",
    "+ 0.5 * M_1 * log(1.0 * N_1)"
  )) {
    expect_true(grepl(term, scalar_direct, fixed = TRUE), info = term)
  }

  ten_centered <- stancode(
    y ~ 0 + ten +
      (0 + ten || gr(g, s2z = TRUE, scale = "varying")),
    data = s2z_ten_dat,
    prior = prior(normal(0, 2), class = b)
  )
  ten_direct <- stancode(
    y ~ 0 + ten +
      (0 + ten || gr(
        g, s2z = TRUE, scale = "varying", center = FALSE
      )),
    data = s2z_ten_dat,
    prior = prior(normal(0, 2), class = b)
  )
  for (j in seq_len(10L)) {
    expect_match2(
      ten_direct,
      sprintf(
        paste0(
          "r_s2z_1_%1$s = sum_to_zero_constrain_brms(",
          "reference_sd_s2z_1[%1$s] * segment(z_s2z_1, ",
          "(%1$s - 1) * (N_1 - 1) + 1, N_1 - 1));"
        ),
        j
      )
    )
  }
  expect_match2(
    ten_centered,
    "- (N_1 - 1) * sum(log(reference_sd_s2z_1))"
  )
  expect_false(grepl(
    "- (N_1 - 1) * sum(log(reference_sd_s2z_1))",
    ten_direct, fixed = TRUE
  ))
  for (term in c(
    paste0(
      "real relative_precision_s2z = reference_sd_s2z_1[1] / ",
      "sd_level_s2z_1[n, 1];"
    ),
    "- 0.5 * sum(log(D_diag_s2z_1))",
    "- 0.5 * log1p(rank1_info_s2z_1)"
  )) {
    expect_true(grepl(term, ten_direct, fixed = TRUE), info = term)
  }
  expect_false(grepl("matrix[M_1, M_1]", ten_direct, fixed = TRUE))

  correlated_centered <- stancode(
    y ~ x * z +
      (1 + x * z | gr(g, s2z = TRUE, scale = "varying")),
    data = s2z_dat
  )
  correlated_direct <- stancode(
    y ~ x * z +
      (1 + x * z | gr(
        g, s2z = TRUE, scale = "varying", center = FALSE
      )),
    data = s2z_dat
  )
  expect_match2(
    correlated_direct, "r_s2z_1 = r_s2z_1 * L_Sigma_s2z_1';"
  )
  expect_match2(
    correlated_centered,
    "- (N_1 - 1) * sum(log(diagonal(L_Sigma_s2z_1)))"
  )
  expect_false(grepl(
    "- (N_1 - 1) * sum(log(diagonal(L_Sigma_s2z_1)))",
    correlated_direct, fixed = TRUE
  ))
  for (term in c(
    "relative_precision_s2z = mdivide_left_tri_low(",
    "P_s2z_1 += crossprod(relative_precision_s2z);",
    "group_quad_s2z_1 -= dot_self(forward_solve_s2z)",
    "- sum(log(diagonal(L_P_s2z_1)))",
    "+ 0.5 * M_1 * log(1.0 * N_1)"
  )) {
    expect_true(grepl(term, correlated_direct, fixed = TRUE), info = term)
  }

  student_direct <- stancode(
    y ~ x * z + (1 + x * z | gr(
      g, s2z = TRUE, scale = "varying", center = FALSE,
      dist = "student"
    )),
    data = s2z_dat
  )
  expect_match2(
    student_direct, "r_s2z_1 = r_s2z_1 * L_Sigma_s2z_1';"
  )
  expect_false(grepl(
    "- (N_1 - 1) * sum(log(diagonal(L_Sigma_s2z_1)))",
    student_direct, fixed = TRUE
  ))
  for (term in c(
    "group_prec_s2z_1 = inv_square(group_scale_s2z_1)",
    "relative_precision_s2z = mdivide_left_tri_low(",
    "- M_1 * sum(log(group_scale_s2z_1))",
    "- sum(log(diagonal(L_P_s2z_1)))"
  )) {
    expect_true(grepl(term, student_direct, fixed = TRUE), info = term)
  }
})

test_that("endpoint and fixed-partial S2Z keep recovery and public names", {
  centered_form <- y ~ x * z +
    (1 + x * z | gr(g, s2z = TRUE, scale = "varying"))
  direct_form <- y ~ x * z +
    (1 + x * z | gr(
      g, s2z = TRUE, scale = "varying", center = FALSE
    ))
  partial_form <- y ~ x * z +
    (1 + x * z | gr(
      g, s2z = TRUE, scale = "varying", center = 0.4
    ))
  second_partial_form <- y ~ x * z +
    (1 + x * z | gr(
      g, s2z = TRUE, scale = "varying", center = 0.65
    ))
  centered_fit <- brm(centered_form, data = s2z_dat, empty = TRUE)
  centered_excluded <- unlist(
    brms:::exclude_pars(centered_fit), use.names = FALSE
  )

  prior_columns <- c(
    "prior", "class", "coef", "group", "resp", "dpar", "nlpar",
    "lb", "ub"
  )
  centered_prior <- as.data.frame(get_prior(
    centered_form, data = s2z_dat
  ))[, prior_columns]

  for (form in list(direct_form, partial_form, second_partial_form)) {
    fit <- brm(form, data = s2z_dat, empty = TRUE)
    expect_identical(
      unlist(brms:::exclude_pars(fit), use.names = FALSE),
      centered_excluded
    )
    expect_identical(
      as.data.frame(get_prior(form, data = s2z_dat))[, prior_columns],
      centered_prior
    )
    code <- stancode(form, data = s2z_dat)
    for (term in c(
      "vector[M_1] mean_r_s2z_1;",
      "vector[4] q_recovered_s2z_1;",
      "matrix<lower=0>[N_1, M_1] sd_level_1;",
      "real Intercept;", "vector[Kc] b;", "real b_Intercept;",
      "matrix[N_1, M_1] r_1;", "corr_matrix[M_1] Cor_1"
    )) {
      expect_true(grepl(term, code, fixed = TRUE), info = term)
    }
  }
})

test_that("partial S2Z Jacobian state follows save_pars", {
  form <- y ~ x +
    (1 + x | gr(g, s2z = TRUE, center = 0.4))
  default_fit <- brm(form, data = s2z_dat, empty = TRUE)
  saved_fit <- brm(
    form, data = s2z_dat, empty = TRUE,
    save_pars = save_pars(all = TRUE)
  )
  default_excluded <- unlist(
    brms:::exclude_pars(default_fit), use.names = FALSE
  )
  saved_excluded <- unlist(
    brms:::exclude_pars(saved_fit), use.names = FALSE
  )

  expect_true("log_det_partial_s2z_1" %in% default_excluded)
  expect_false("log_det_partial_s2z_1" %in% saved_excluded)
})

test_that("S2Z uses the dedicated Gaussian scalar kernel", {
  form <- y ~ 1 + (1 | gr(g, s2z = TRUE))
  scode <- stancode(
    form, data = s2z_dat,
    prior = prior(normal(0, 2), class = Intercept)
  )
  sdata <- standata(form, data = s2z_dat)

  expect_equal(sdata$M_1, 1)
  expect_match2(scode, "vector[N_1 - 1] z_s2z_1;")
  expect_match2(scode, "vector[N_1] r_s2z_1_1;")
  expect_match2(scode, "vector[1] H_s2z_1;")
  expect_match2(scode, "real<lower=0> D_s2z_1;")
  expect_match2(scode, "real<lower=0> sqrt_D_s2z_1;")
  expect_match2(scode, "real mhat_s2z_1;")
  expect_match2(scode, "real<lower=0> group_quad_s2z_1;")
  expect_match2(scode, "vector[N_1] white_s2z =")
  expect_match2(
    scode,
    "r_s2z_1_1 = sum_to_zero_constrain_brms(z_s2z_1);"
  )
  expect_match2(
    scode,
    "D_s2z_1 = tau_sq_s2z * prior_info_s2z + N_1;"
  )
  expect_false(grepl(
    "dot_product(r_s2z_1_1, group_prec_s2z_1)", scode, fixed = TRUE
  ))
  expect_false(grepl("group_scale_s2z_1", scode, fixed = TRUE))
  expect_false(grepl("group_prec_s2z_1", scode, fixed = TRUE))
  expect_match2(scode, "- (N_1 - 1) * log(sd_1[1])")
  expect_match2(scode, "- 0.5 * log(D_s2z_1)")
  expect_match2(scode, "+ 0.5 * log(1.0 * N_1)")

  # The scalar branch must not instantiate one-by-one matrix factorizations.
  for (term in c(
    "array[M_1] vector[N_1 - 1] z_s2z_1;",
    "vector[M_1 * (N_1 - 1)] z_s2z_1;",
    "matrix[N_1, M_1] r_s2z_1;",
    "matrix[M_1, M_1] Q_Sigma_s2z_1;",
    "cholesky_decompose(P_s2z_1)",
    "mdivide_left_spd(P_s2z_1",
    "mdivide_left_tri_low(L_Sigma_s2z_1"
  )) {
    expect_false(grepl(term, scode, fixed = TRUE))
  }

  expect_match2(scode, "real mean_r_s2z_1;")
  expect_match2(
    scode,
    paste0(
      "mean_r_s2z_1 = mhat_s2z_1 + sd_1[1] * std_normal_rng() / ",
      "sqrt_D_s2z_1;"
    )
  )
  expect_match2(
    scode,
    "r_1_1 = r_s2z_1_1 + mean_r_s2z_1;"
  )
  expect_match2(scode, "Intercept = q_recovered_s2z_1[1];")
  expect_match2(scode, "b_Intercept = Intercept;")
})

test_that("S2Z scalar kernel handles Student group scales exactly", {
  scode <- stancode(
    y ~ 1 + (1 | gr(g, dist = "student", s2z = TRUE)),
    data = s2z_dat,
    prior = prior(normal(0, 2), class = Intercept)
  )

  expect_match2(scode, "vector<lower=0>[N_1] udf_1;")
  expect_match2(scode, "group_scale_s2z_1 = dfm_1;")
  expect_match2(
    scode,
    paste0(
      "D_s2z_1 = tau_sq_s2z * prior_info_s2z + ",
      "sum(group_prec_s2z_1);"
    )
  )
  expect_match2(
    scode,
    "dot_product(r_s2z_1_1, group_prec_s2z_1)"
  )
  expect_match2(scode, "- sum(log(group_scale_s2z_1))")
  expect_match2(scode, "inv_chi_square_lpdf(udf_1 | df_1)")
  expect_false(grepl("matrix[M_1, M_1] P_s2z_1;", scode, fixed = TRUE))
})

test_that("S2Z scalar kernel handles a varying slope without an intercept", {
  form <- y ~ 0 + x + (0 + x | gr(g, s2z = TRUE))
  scode <- stancode(
    form, data = s2z_dat,
    prior = prior(normal(0, 2), class = b, coef = x)
  )
  sdata <- standata(form, data = s2z_dat)

  expect_equal(sdata$M_1, 1)
  expect_match2(scode, "vector[N_1 - 1] z_s2z_1;")
  expect_match2(scode, "H_s2z_1[1] = 1.0;")
  expect_match2(
    scode,
    "normal_id_glm_lpdf(Y | X, mu, theta_s2z, sigma)"
  )
  expect_match2(scode, "vector[1] b;")
  expect_match2(scode, "b = q_recovered_s2z_1;")
  expect_match2(scode, "vector[N_1] r_1_1;")
  expect_false(grepl("real b_Intercept;", scode, fixed = TRUE))
})

test_that("S2Z scalar kernel maps a centered varying interaction", {
  form <- y ~ x * z + (0 + x:z | gr(g, s2z = TRUE))
  scode <- stancode(
    form, data = s2z_dat,
    prior = prior(normal(0, 2), class = b)
  )
  sdata <- standata(form, data = s2z_dat)

  expect_equal(sdata$M_1, 1)
  expect_equal(sdata$Kc, 3)
  expect_match2(scode, "vector[4] H_s2z_1;")
  expect_match2(scode, "H_s2z_1[4] = 1.0;")
  expect_match2(scode, "H_s2z_1[1] = means_X[3];")
  expect_match2(
    scode,
    "normal_id_glm_lpdf(Y | Xc, mu, tail(theta_s2z, 3), sigma)"
  )
})

test_that("S2Z scalar fixed SD and unnormalized code retain the kernel", {
  form <- y ~ 1 + (1 | gr(g, s2z = TRUE))
  bprior <- prior(normal(0, 2), class = Intercept) +
    prior(constant(1.25), class = sd)
  scode <- stancode(
    form, data = s2z_dat, prior = bprior, normalize = FALSE
  )

  expect_match2(scode, "sd_1 = rep_vector(1.25, rows(sd_1));")
  expect_equal(
    unname(lengths(regmatches(
      scode,
      gregexpr(
        "sd_1 = rep_vector(1.25, rows(sd_1));", scode, fixed = TRUE
      )
    ))),
    1L
  )
  expect_match2(scode, "real tau_sq_s2z = square(sd_1[1]);")
  expect_match2(scode, "- (N_1 - 1) * log(sd_1[1])")
  expect_match2(scode, "- 0.5 * log(D_s2z_1)")
  expect_false(grepl("log(1.0 * N_1)", scode, fixed = TRUE))
  expect_false(grepl(
    "+ 0.5 * log(2 * pi())", scode, fixed = TRUE
  ))
  expect_error(
    stancode(
      form, data = s2z_dat,
      prior = prior(constant(0), class = sd)
    ),
    "standard deviations fixed with 'constant' must be positive"
  )
})

test_that("threaded S2Z scalar likelihood receives the internal vector", {
  scode <- stancode(
    y ~ 1 + (1 | gr(g, s2z = TRUE)), data = s2z_dat,
    threads = threading(2), parse = FALSE
  )

  expect_match2(
    scode,
    paste0(
      "real partial_log_lik_lpmf(array[] int seq, int start, int end, ",
      "data vector Y, vector theta_s2z, real sigma, data array[] int J_1, ",
      "data vector Z_1_1, vector r_s2z_1_1)"
    )
  )
  expect_match2(
    scode,
    paste0(
      "target += reduce_sum(partial_log_lik_lpmf, seq, grainsize, Y, ",
      "theta_s2z, sigma, J_1, Z_1_1, r_s2z_1_1);"
    )
  )
})

test_that("S2Z Stan code covers intercepts, slopes, and interactions", {
  scode <- stancode(
    y ~ x * z + (1 + x * z | gr(g, s2z = TRUE)), data = s2z_dat
  )
  sdata <- standata(
    y ~ x * z + (1 + x * z | gr(g, s2z = TRUE)), data = s2z_dat
  )

  expect_equal(sdata$M_1, 4)
  expect_equal(sdata$NC_1, 6)
  expect_match2(scode, "vector[M_1 * (N_1 - 1)] z_s2z_1;")
  expect_match2(scode, "matrix[N_1, M_1] r_s2z_1;")
  expect_match2(scode, "H_s2z_1[1, 2] = means_X[1];")
  expect_match2(scode, "H_s2z_1[1, 3] = means_X[2];")
  expect_match2(scode, "H_s2z_1[1, 4] = means_X[3];")
  expect_match2(scode, "r_s2z_1[, k] = sum_to_zero_constrain_brms")
  expect_match2(scode, "mu += theta_s2z[1]")
  expect_match2(
    scode,
    "normal_id_glm_lpdf(Y | Xc, mu, tail(theta_s2z, 3), sigma)"
  )
  expect_match2(scode, "q_recovered_s2z_1 = theta_s2z")
  expect_match2(scode, "r_1 = r_s2z_1;")
  expect_match2(scode, "for (j in 1:N_1) r_1[j] += mean_r_s2z_1';")
  expect_match2(scode, "b_Intercept = Intercept - dot_product(means_X, b)")
  expect_true(grepl("normal_id_glm", scode, fixed = TRUE))
})

test_that("S2Z uses actual mixed-interaction design columns", {
  form <- y ~ f * x + (1 + f * x | gr(g, s2z = TRUE))
  sdata <- standata(form, data = s2z_dat)
  Z <- model.matrix(~ f * x, s2z_dat)

  expect_equal(sdata$M_1, ncol(Z))
  expect_equal(sdata$NC_1, choose(ncol(Z), 2))
  for (k in seq_len(ncol(Z))) {
    expect_equal(
      unname(sdata[[sprintf("Z_1_%s", k)]]),
      as.array(unname(Z[, k]))
    )
  }

  scode <- stancode(form, data = s2z_dat)
  for (k in seq_len(ncol(Z) - 1L)) {
    expect_match2(
      scode,
      sprintf("H_s2z_1[1, %s] = means_X[%s];", k + 1L, k)
    )
  }
})

test_that("S2Z supports correlated and diagonal Gaussian effects", {
  sc_cor <- stancode(
    y ~ x + (1 + x | gr(g, s2z = TRUE)), data = s2z_dat
  )
  sc_diag <- stancode(
    y ~ x + (1 + x || gr(g, s2z = TRUE)), data = s2z_dat
  )

  expect_match2(sc_cor, "cholesky_factor_corr[M_1] L_1;")
  expect_match2(sc_cor, "diag_pre_multiply(sd_1, L_1)")
  expect_match2(sc_cor, "corr_matrix[M_1] Cor_1")
  expect_match2(
    sc_cor,
    paste0(
      "prior_factor_s2z = diag_pre_multiply(",
      "sqrt(prior_prec_s2z_1), H_s2z_1) * L_Sigma_s2z_1;"
    )
  )
  expect_match2(
    sc_cor,
    paste0(
      "matrix[M_1, N_1] white_s2z = mdivide_left_tri_low(",
      "L_Sigma_s2z_1, r_s2z_1');"
    )
  )
  expect_match2(
    sc_cor,
    paste0(
      "P_s2z_1 = add_diag(crossprod(prior_factor_s2z), ",
      "1.0 * N_1);"
    )
  )
  expect_match2(sc_cor, "L_P_s2z_1 = cholesky_decompose(P_s2z_1);")
  expect_match2(
    sc_cor,
    paste0(
      "whitened_h_s2z = mdivide_left_tri_low(",
      "L_P_s2z_1, h_s2z);"
    )
  )
  expect_match2(
    sc_cor, "group_quad_s2z_1 = -dot_self(whitened_h_s2z);"
  )
  expect_match2(sc_cor, "lprior += -0.5 * group_quad_s2z_1")
  expect_match2(
    sc_cor,
    paste0(
      "mean_r_s2z_1 = L_Sigma_s2z_1 * (r_mean_s2z + ",
      "(mdivide_right_tri_low(z_mean_s2z', L_P_s2z_1))');"
    )
  )
  for (term in c(
    "Q_Sigma_s2z_1",
    "L_inv_s2z",
    "mdivide_left_spd(",
    "mdivide_left_tri_low(L_Sigma_s2z_1, diag_matrix",
    "qhat_s2z_1"
  )) {
    expect_false(grepl(term, sc_cor, fixed = TRUE), info = term)
  }
  expect_false(grepl("cholesky_factor_corr[M_1] L_1;", sc_diag, fixed = TRUE))
  expect_match2(sc_diag, "D_diag_s2z_1")
  expect_match2(sc_diag, "rank1_info_s2z_1")
  expect_false(grepl("L_Sigma_s2z_1", sc_diag, fixed = TRUE))
})

test_that("independent S2Z specializes centered slopes and interactions", {
  form <- y ~ x * z + (1 + x * z || gr(g, s2z = TRUE))
  scode <- stancode(
    form, data = s2z_dat,
    prior = prior(normal(0, 2), class = Intercept) +
      prior(normal(0, 1.5), class = b)
  )
  sdata <- standata(form, data = s2z_dat)

  expect_equal(sdata$M_1, 4)
  expect_match2(scode, "vector[M_1 * (N_1 - 1)] z_s2z_1;")
  for (k in seq_len(sdata$M_1)) {
    expect_match2(scode, sprintf("vector[N_1] r_s2z_1_%s;", k))
    expect_match2(
      scode,
      sprintf(
        paste0(
          "r_s2z_1_%1$s = sum_to_zero_constrain_brms(",
          "segment(z_s2z_1, (%1$s - 1) * (N_1 - 1) + 1, N_1 - 1));"
        ),
        k
      )
    )
  }
  expect_match2(scode, "intercept_map_s2z_1[1] = 1.0;")
  expect_match2(scode, "intercept_map_s2z_1[2] = means_X[1];")
  expect_match2(scode, "intercept_map_s2z_1[3] = means_X[2];")
  expect_match2(scode, "intercept_map_s2z_1[4] = means_X[3];")
  expect_match2(
    scode,
    "D_diag_s2z_1 = group_info_s2z + square(sd_1) .* base_info_s2z;"
  )
  expect_match2(scode, "real group_info_s2z = N_1;")
  expect_false(grepl("group_scale_s2z_1", scode, fixed = TRUE))
  expect_false(grepl("group_prec_s2z_1", scode, fixed = TRUE))
  expect_match2(
    scode,
    "rank1_info_s2z_1 = prior_prec_s2z_1[1] * dot_product("
  )
  expect_match2(scode, "- 0.5 * sum(log(D_diag_s2z_1))")
  expect_match2(scode, "- 0.5 * log1p(rank1_info_s2z_1)")
  expect_match2(scode, "+ 0.5 * M_1 * log(1.0 * N_1)")
  expect_match2(scode, "real sqrt_rank1_s2z")
  expect_match2(scode, "real rank1_adjust_s2z;")
  expect_match2(
    scode,
    "q_recovered_s2z_1[1] -= dot_product(intercept_map_s2z_1"
  )
  expect_match2(scode, "r_1_4 = r_s2z_1_4 + mean_r_s2z_1[4];")
  expect_false(grepl("corr_matrix[M_1] Cor_1", scode, fixed = TRUE))

  # The independent specialization must be O(K): no dense K by K storage,
  # factorization, or positive-definite solve is emitted anywhere in the block.
  for (term in c(
    "matrix[N_1, M_1] r_s2z_1;",
    "matrix[M_1, M_1]",
    "L_Sigma_s2z_1",
    "Q_Sigma_s2z_1",
    "P_s2z_1",
    "L_P_s2z_1",
    "cholesky_decompose(",
    "mdivide_left_spd(",
    "mdivide_left_tri_low(",
    "mdivide_right_tri_low("
  )) {
    expect_false(grepl(term, scode, fixed = TRUE), info = term)
  }
})

test_that("independent S2Z uses H identity for ten no-intercept effects", {
  ten_dat <- data.frame(
    y = seq(-1, 1, length.out = 80),
    ten = factor(rep(letters[1:10], 8)),
    g = factor(rep(seq_len(8), each = 10))
  )
  form <- y ~ 0 + ten + (0 + ten || gr(g, s2z = TRUE))
  scode <- stancode(
    form, data = ten_dat,
    prior = prior(normal(0, 2), class = b), normalize = FALSE
  )
  sdata <- standata(form, data = ten_dat)

  expect_equal(sdata$M_1, 10)
  expect_equal(sdata$K, 10)
  expect_match2(scode, "vector[10] theta_s2z;")
  expect_match2(scode, "intercept_map_s2z_1 = zeros_vector(M_1);")
  expect_false(grepl("intercept_map_s2z_1[1] =", scode, fixed = TRUE))
  expect_match2(scode, "rank1_info_s2z_1 = 0.0 * dot_product(")
  for (k in seq_len(10L)) {
    expect_match2(
      scode,
      sprintf("base_info_s2z[%1$s] = prior_prec_s2z_1[%1$s];", k)
    )
    expect_match2(
      scode,
      sprintf("qhat_s2z_1[%1$s] -= mhat_s2z_1[%1$s];", k)
    )
    expect_match2(
      scode,
      sprintf(
        "q_recovered_s2z_1[%1$s] -= mean_r_s2z_1[%1$s];", k
      )
    )
  }
  expect_match2(scode, "b = q_recovered_s2z_1;")
  expect_false(grepl("real b_Intercept;", scode, fixed = TRUE))
  expect_false(grepl("matrix[M_1, M_1]", scode, fixed = TRUE))
  expect_false(grepl("cholesky_decompose(", scode, fixed = TRUE))
  expect_false(grepl("mdivide_left_spd(", scode, fixed = TRUE))
  expect_false(grepl("log(1.0 * N_1)", scode, fixed = TRUE))
})

test_that("independent Student S2Z retains weighted scales unnormalized", {
  scode <- stancode(
    y ~ x * z +
      (1 + x * z || gr(g, dist = "student", s2z = TRUE)),
    data = s2z_dat, normalize = FALSE
  )

  expect_match2(scode, "group_scale_s2z_1 = dfm_1;")
  expect_match2(
    scode,
    "group_prec_s2z_1 = inv_square(group_scale_s2z_1);"
  )
  expect_match2(scode, "real group_info_s2z = sum(group_prec_s2z_1);")
  for (k in seq_len(4L)) {
    expect_match2(
      scode,
      sprintf(
        "dot_product(r_s2z_1_%s, group_prec_s2z_1)", k
      )
    )
  }
  expect_match2(scode, "- M_1 * sum(log(group_scale_s2z_1))")
  expect_match2(scode, "- (N_1 - 1) * sum(log(sd_1))")
  expect_match2(scode, "- 0.5 * sum(log(D_diag_s2z_1))")
  expect_match2(scode, "- 0.5 * log1p(rank1_info_s2z_1)")
  expect_false(grepl("log(1.0 * N_1)", scode, fixed = TRUE))
  expect_false(grepl("matrix[M_1, M_1]", scode, fixed = TRUE))
  expect_false(grepl("cholesky_decompose(", scode, fixed = TRUE))
  expect_false(grepl("mdivide_left_spd(", scode, fixed = TRUE))
  expect_false(grepl("mdivide_left_tri_low(", scode, fixed = TRUE))
  expect_false(grepl("mdivide_right_tri_low(", scode, fixed = TRUE))
})

test_that("independent S2Z handles slope subsets and heavy-tailed priors", {
  form <- y ~ x * z + (0 + x + x:z || gr(g, s2z = TRUE))
  bprior <- c(
    prior(cauchy(0, 1), class = Intercept),
    prior(student_t(5, 0, 1.5), class = b, coef = x),
    prior(normal(0, 2), class = b, coef = "x:z"),
    prior(constant(0.4), class = sd, group = g, coef = x),
    prior(constant(0.7), class = sd, group = g, coef = "x:z")
  )
  scode <- stancode(form, data = s2z_dat, prior = bprior)
  sdata <- standata(form, data = s2z_dat, prior = bprior)

  expect_equal(sdata$M_1, 2)
  expect_equal(sdata$Kc, 3)
  expect_match2(scode, "intercept_map_s2z_1[1] = means_X[1];")
  expect_match2(scode, "intercept_map_s2z_1[2] = means_X[3];")
  expect_match2(scode, "base_info_s2z[1] = prior_prec_s2z_1[2];")
  expect_match2(scode, "base_info_s2z[2] = prior_prec_s2z_1[4];")
  expect_match2(scode, "qhat_s2z_1[2] -= mhat_s2z_1[1];")
  expect_match2(scode, "qhat_s2z_1[4] -= mhat_s2z_1[2];")
  expect_false(grepl("qhat_s2z_1[3] -=", scode, fixed = TRUE))
  expect_match2(scode, "inv_chi_square_lpdf(udf_b_s2z_1 | 1)")
  expect_match2(scode, "inv_chi_square_lpdf(udf_b_s2z_2 | 5)")
  expect_false(grepl("matrix[M_1, M_1]", scode, fixed = TRUE))
})

test_that("S2Z preserves ordinary priors on every varying effect", {
  form <- y ~ x * z + (1 + x * z | gr(g, s2z = TRUE))
  varying_priors <- c(
    prior(exponential(2), class = "sd", group = "g", coef = "Intercept"),
    prior(normal(0, 0.4), class = "sd", group = "g", coef = "x"),
    prior(student_t(4, 0, 0.3), class = "sd", group = "g", coef = "z"),
    prior(exponential(3), class = "sd", group = "g", coef = "x:z"),
    prior(lkj(4), class = "cor", group = "g")
  )
  available <- as.data.frame(get_prior(form, data = s2z_dat))
  group_sd <- subset(
    available, class == "sd" & group == "g" & nzchar(coef)
  )
  expect_setequal(
    group_sd$coef, c("Intercept", "x", "z", "x:z")
  )

  scode <- stancode(form, data = s2z_dat, prior = varying_priors)
  expect_match2(scode, "exponential_lpdf(sd_1[1] | 2)")
  expect_match2(scode, "normal_lpdf(sd_1[2] | 0, 0.4)")
  expect_match2(scode, "student_t_lpdf(sd_1[3] | 4, 0, 0.3)")
  expect_match2(scode, "exponential_lpdf(sd_1[4] | 3)")
  expect_match2(scode, "lkj_corr_cholesky_lpdf(L_1 | 4)")

  scalar_code <- stancode(
    y ~ 0 + x + (0 + x | gr(g, s2z = TRUE)), data = s2z_dat,
    prior = prior(exponential(7), class = "sd", group = "g", coef = "x")
  )
  expect_match2(scalar_code, "exponential_lpdf(sd_1[1] | 7)")

  student_code <- stancode(
    y ~ x + (1 + x | gr(g, dist = "student", s2z = TRUE)),
    data = s2z_dat,
    prior = prior(gamma(3, 0.25), class = "df", group = "g")
  )
  expect_match2(student_code, "gamma_lpdf(df_1 | 3, 0.25)")
  expect_match2(student_code, "inv_chi_square_lpdf(udf_1 | df_1)")
})

test_that("group-varying scales preserve baseline priors and add sdlog", {
  form <- y ~ x * z +
    (1 + x * z | gr(g, s2z = TRUE, scale = "varying"))
  available <- as.data.frame(get_prior(form, data = s2z_dat))
  group_sd <- subset(
    available, class == "sd" & group == "g" & nzchar(coef)
  )
  group_sdlog <- subset(
    available, class == "sdlog" & group == "g" & nzchar(coef)
  )
  default_sdlog <- subset(
    available, class == "sdlog" & !nzchar(group) & !nzchar(coef)
  )

  expect_setequal(group_sd$coef, c("Intercept", "x", "z", "x:z"))
  expect_setequal(group_sdlog$coef, group_sd$coef)
  expect_identical(default_sdlog$prior, "normal(0, 0.25)")
  expect_identical(default_sdlog$lb, "0")

  varying_priors <- c(
    prior(exponential(2), class = "sd", group = "g",
          coef = "Intercept"),
    prior(lognormal(-1, 0.4), class = "sd", group = "g", coef = "x"),
    prior(normal(0, 0.1), class = "sdlog", group = "g",
          coef = "Intercept"),
    prior(exponential(5), class = "sdlog", group = "g", coef = "x"),
    prior(normal(0, 0.3), class = "sdlog", group = "g", coef = "z"),
    prior(student_t(4, 0, 0.2), class = "sdlog", group = "g",
          coef = "x:z"),
    prior(lkj(3), class = "cor", group = "g")
  )
  scode <- stancode(form, data = s2z_dat, prior = varying_priors)

  expect_match2(scode, "vector<lower=0>[M_1] sd_1;")
  expect_match2(scode, "vector<lower=0>[M_1] sdlog_1;")
  expect_match2(scode, "exponential_lpdf(sd_1[1] | 2)")
  expect_match2(scode, "lognormal_lpdf(sd_1[2] | -1, 0.4)")
  expect_match2(scode, "normal_lpdf(sdlog_1[1] | 0, 0.1)")
  expect_match2(scode, "exponential_lpdf(sdlog_1[2] | 5)")
  expect_match2(scode, "normal_lpdf(sdlog_1[3] | 0, 0.3)")
  expect_match2(scode, "student_t_lpdf(sdlog_1[4] | 4, 0, 0.2)")
  expect_match2(scode, "lkj_corr_cholesky_lpdf(L_1 | 3)")

  implicit_shared <- stancode(y ~ x + (1 + x | g), data = s2z_dat)
  explicit_shared <- stancode(
    y ~ x + (1 + x | gr(g, scale = "shared")), data = s2z_dat
  )
  expect_identical(explicit_shared, implicit_shared)
  expect_error(
    stancode(
      y ~ x + (1 + x | gr(g, scale = "varying")), data = s2z_dat
    ),
    "require 's2z = TRUE'"
  )
  expect_error(
    stancode(
      form, data = s2z_dat,
      prior = prior(constant(-0.1), class = "sdlog")
    ),
    "must be non-negative"
  )
  expect_error(
    stancode(
      form, data = s2z_dat,
      prior = prior(normal(0, 0.2), class = "sdlog", lb = -1)
    ),
    "finite non-negative lower bound"
  )
  expect_error(
    stancode(
      form, data = s2z_dat,
      prior = prior(normal(0, 0.2), class = "sdlog", lb = 0, ub = 0)
    ),
    "finite positive upper bound"
  )
  expect_error(
    stancode(
      form, data = s2z_dat,
      prior = prior(normal(1, 0.2), class = "sd", lb = -1)
    ),
    "finite non-negative lower bound"
  )

  zero_code <- stancode(
    form, data = s2z_dat,
    prior = prior(constant(0), class = "sdlog")
  )
  expect_match2(zero_code, "sdlog_1 = rep_vector(0, rows(sdlog_1));")
  expect_match2(
    zero_code,
    paste0(
      "reference_sd_s2z_1 = sd_1 .* exp(sdlog_1 .* ",
      "tail(z_sd_s2z_1, M_1) / sqrt(1.0 * N_1));"
    )
  )

  mixed_code <- stancode(
    form, data = s2z_dat, normalize = FALSE,
    prior = c(
      prior(constant(0), class = "sdlog", group = "g",
            coef = "Intercept"),
      prior(normal(0, 0.3), class = "sdlog", group = "g", coef = "x")
    )
  )
  mixed_lines <- strsplit(mixed_code, "\n", fixed = TRUE)[[1]]
  expect_equal(sum(grepl("sdlog_1[1] = 0;", mixed_lines, fixed = TRUE)), 1)
  expect_match2(mixed_code, "normal_lupdf(sdlog_1[2] | 0, 0.3)")
})

test_that("realized scale priors are compactly discoverable and additive", {
  form <- y ~ x +
    (1 + x | gr(g, s2z = TRUE, scale = "varying"))
  available <- as.data.frame(get_prior(form, data = s2z_dat))
  level_rows <- subset(
    available, class == "sd_level" & group == "g"
  )

  # Availability stays O(M), not O(N * M). Users select fitted-data levels
  # explicitly when constructing an additional prior.
  expect_equal(nrow(level_rows), 2L)
  expect_setequal(level_rows$coef, c("Intercept", "x"))
  expect_identical(level_rows$level, rep("", 2L))
  expect_identical(level_rows$prior, rep("", 2L))
  expect_false(any(get_prior(
    y ~ x + (1 + x | gr(g, s2z = TRUE, scale = "shared")),
    data = s2z_dat
  )$class == "sd_level"))

  one_prior <- prior(
    normal(0.4, 0.2), class = "sd_level", group = "g",
    coef = "Intercept", level = "2"
  )
  expect_identical(one_prior$level, "2")
  validated <- validate_prior(one_prior, form, data = s2z_dat)
  selected <- subset(validated, class == "sd_level" & nzchar(prior))
  expect_equal(nrow(selected), 1L)
  expect_identical(selected$level, "2")

  base_code <- stancode(form, data = s2z_dat)
  one_code <- stancode(form, data = s2z_dat, prior = one_prior)
  target <- "normal_lpdf(sd_level_s2z_1[2, 1] | 0.4, 0.2)"
  truncation <- "- 1 * normal_lccdf(0 | 0.4, 0.2)"
  expect_equal(s2z_count_fixed(one_code, target), 1L)
  expect_equal(s2z_count_fixed(one_code, truncation), 1L)
  expect_false(grepl(target, base_code, fixed = TRUE))
  expect_false(grepl(
    "_lpdf(sd_level_s2z_", base_code, fixed = TRUE
  ))

  # Adding the factor must not alter parameter declarations, hierarchy terms,
  # or add a change-of-variables Jacobian. Apart from the explanatory comment,
  # the lpdf and its positive-support normalizer are the entire code delta.
  base_lines <- strsplit(base_code, "\n", fixed = TRUE)[[1]]
  one_lines <- strsplit(one_code, "\n", fixed = TRUE)[[1]]
  added <- setdiff(one_lines, base_lines)
  expect_length(setdiff(base_lines, one_lines), 0L)
  expect_equal(length(added), 3L)
  expect_true(any(grepl(
    "additional priors on realized group-level standard deviations",
    added, fixed = TRUE
  )))
  expect_true(any(grepl(target, added, fixed = TRUE)))
  expect_true(any(grepl(truncation, added, fixed = TRUE)))
  expect_false(any(grepl("jacobian", added, ignore.case = TRUE)))

  many_prior <- set_prior(
    c("lognormal(-1, 0.3)", "exponential(4)"),
    class = "sd_level", group = "g", coef = "x",
    level = c("1", "6")
  )
  many_code <- stancode(
    form, data = s2z_dat, prior = many_prior, normalize = FALSE
  )
  expect_match2(
    many_code,
    "lognormal_lupdf(sd_level_s2z_1[1, 2] | -1, 0.3)"
  )
  expect_match2(
    many_code,
    "exponential_lupdf(sd_level_s2z_1[6, 2] | 4)"
  )
  expect_false(grepl("_lccdf(", many_code, fixed = TRUE))
  expect_equal(
    s2z_count_fixed(
      many_code,
      "additional priors on realized group-level standard deviations"
    ),
    1L
  )

  named_dat <- transform(
    s2z_dat,
    surgeon = factor(g, labels = paste("Mr", LETTERS[seq_len(6L)]))
  )
  named_code <- stancode(
    y ~ 1 +
      (1 | gr(surgeon, s2z = TRUE, scale = "varying")),
    data = named_dat,
    prior = prior(
      lognormal(-1, 0.2), class = "sd_level", group = "surgeon",
      coef = "Intercept", level = "Mr B"
    )
  )
  expect_match2(
    named_code,
    "lognormal_lpdf(sd_level_s2z_1[2, 1] | -1, 0.2)"
  )
})

test_that("realized scale priors target coefficients, blocks, and dpars", {
  interaction_form <- y ~ x * z +
    (1 + x * z | gr(g, s2z = TRUE, scale = "varying"))
  interaction_prior <- c(
    prior(
      lognormal(-0.5, 0.25), class = "sd_level", group = "g",
      coef = "Intercept", level = "1"
    ),
    prior(
      gamma(3, 4), class = "sd_level", group = "g",
      coef = "x", level = "3"
    ),
    prior(
      exponential(5), class = "sd_level", group = "g",
      coef = "x:z", level = "6"
    )
  )
  interaction_code <- stancode(
    interaction_form, data = s2z_dat, prior = interaction_prior
  )
  for (term in c(
    "lognormal_lpdf(sd_level_s2z_1[1, 1] | -0.5, 0.25)",
    "gamma_lpdf(sd_level_s2z_1[3, 2] | 3, 4)",
    "exponential_lpdf(sd_level_s2z_1[6, 4] | 5)",
    "cholesky_factor_corr[M_1] L_1;"
  )) {
    expect_true(grepl(term, interaction_code, fixed = TRUE), info = term)
  }

  independent_form <- y ~ x +
    (1 + x || gr(g, s2z = TRUE, scale = "varying"))
  independent_code <- stancode(
    independent_form, data = s2z_dat,
    prior = prior(
      lognormal(-1, 0.2), class = "sd_level", group = "g",
      coef = "x", level = "4"
    )
  )
  expect_match2(
    independent_code,
    "lognormal_lpdf(sd_level_s2z_1[4, 2] | -1, 0.2)"
  )
  expect_false(grepl(
    "cholesky_factor_corr[M_1] L_1;", independent_code, fixed = TRUE
  ))

  multiblock_form <- y ~ x +
    (1 + x | gr(g, s2z = TRUE, scale = "varying")) +
    (1 + x || gr(h, s2z = TRUE, scale = "varying"))
  multiblock_prior <- c(
    prior(
      lognormal(-1, 0.2), class = "sd_level", group = "g",
      coef = "Intercept", level = "2"
    ),
    prior(
      exponential(3), class = "sd_level", group = "h",
      coef = "x", level = "8"
    )
  )
  multiblock_code <- stancode(
    multiblock_form, data = s2z_dat, prior = multiblock_prior
  )
  expect_match2(
    multiblock_code,
    "lognormal_lpdf(sd_level_s2z_1[2, 1] | -1, 0.2)"
  )
  expect_match2(
    multiblock_code,
    "exponential_lpdf(sd_level_s2z_2[8, 2] | 3)"
  )

  dpar_form <- bf(
    y ~ x +
      (1 + x | gr(g, s2z = TRUE, scale = "varying")),
    sigma ~ z +
      (1 + z | gr(h, s2z = TRUE, scale = "varying"))
  )
  dpar_available <- as.data.frame(get_prior(dpar_form, data = s2z_dat))
  sigma_level_rows <- subset(
    dpar_available,
    class == "sd_level" & group == "h" & dpar == "sigma"
  )
  expect_equal(nrow(sigma_level_rows), 2L)
  expect_setequal(sigma_level_rows$coef, c("Intercept", "z"))
  dpar_code <- stancode(
    dpar_form, data = s2z_dat,
    prior = prior(
      lognormal(-0.7, 0.15), class = "sd_level", group = "h",
      coef = "z", level = "8", dpar = "sigma"
    )
  )
  expect_match2(
    dpar_code,
    "lognormal_lpdf(sd_level_s2z_2[8, 2] | -0.7, 0.15)"
  )
})

test_that("realized scale prior selectors and densities validate strictly", {
  form <- y ~ x +
    (1 + x | gr(g, s2z = TRUE, scale = "varying"))
  selected <- function(..., prior = "normal(0, 1)") {
    set_prior(
      prior, class = "sd_level", group = "g", coef = "Intercept",
      level = "1", ...
    )
  }

  missing_selectors <- list(
    set_prior(
      "normal(0, 1)", class = "sd_level", coef = "Intercept", level = "1"
    ),
    set_prior(
      "normal(0, 1)", class = "sd_level", group = "g", level = "1"
    ),
    set_prior(
      "normal(0, 1)", class = "sd_level", group = "g",
      coef = "Intercept"
    )
  )
  for (bad_prior in missing_selectors) {
    expect_error(
      stancode(form, data = s2z_dat, prior = bad_prior),
      "require nonempty 'group', 'coef', and 'level'"
    )
  }
  for (bad_prior in list(
    set_prior(
      "normal(0, 1)", class = "sd_level", group = "unknown",
      coef = "Intercept", level = "1"
    ),
    set_prior(
      "normal(0, 1)", class = "sd_level", group = "g",
      coef = "unknown", level = "1"
    )
  )) {
    expect_error(
      stancode(form, data = s2z_dat, prior = bad_prior),
      "does not correspond to a coefficient"
    )
  }
  expect_error(
    stancode(
      form, data = s2z_dat,
      prior = set_prior(
        "normal(0, 1)", class = "sd_level", group = "g",
        coef = "Intercept", level = "unknown"
      )
    ),
    "not found in grouping factor"
  )
  expect_error(
    stancode(
      form, data = s2z_dat,
      prior = set_prior(
        "", class = "sd_level", group = "g",
        coef = "Intercept", level = "1"
      )
    ),
    "must specify a nonempty distribution"
  )

  expect_error(
    stancode(
      y ~ x + (1 + x | gr(g, s2z = TRUE, scale = "shared")),
      data = s2z_dat, prior = selected()
    ),
    "does not correspond to a coefficient"
  )
  expect_error(
    stancode(
      form, data = s2z_dat,
      prior = set_prior(
        "normal(0, 1)", class = "sd", group = "g",
        coef = "Intercept", level = "1"
      )
    ),
    "Argument 'level' is only supported for class 'sd_level'"
  )
  expect_error(
    stancode(
      form, data = s2z_dat,
      prior = selected(prior = "constant(0.7)")
    ),
    "Constant priors are not supported for class 'sd_level'"
  )
  expect_error(
    stancode(
      form, data = s2z_dat,
      prior = selected(prior = "horseshoe(1)")
    ),
    "Special shrinkage priors are not supported for class 'sd_level'"
  )
  for (discrete in c(
    "bernoulli(0.5)", "poisson(2)", "yule_simon(2)",
    "dirichlet_multinomial(rep_vector(1, 2))"
  )) {
    expect_error(
      stancode(
        form, data = s2z_dat, prior = selected(prior = discrete)
      ),
      "Discrete distributions are not supported for class 'sd_level'"
    )
  }
  expect_error(
    stancode(
      form, data = s2z_dat,
      prior = selected(tag = "level_scale")
    ),
    "Prior argument 'tag' is not supported for class 'sd_level'"
  )
  expect_error(
    stancode(
      form, data = s2z_dat,
      prior = selected(lb = 0)
    ),
    "Prior bounds are not supported for class 'sd_level'"
  )
  expect_error(
    stancode(
      form, data = s2z_dat,
      prior = selected(ub = 3)
    ),
    "Prior bounds are not supported for class 'sd_level'"
  )

  dpar_form <- bf(
    y ~ x,
    sigma ~ z +
      (1 + z | gr(h, s2z = TRUE, scale = "varying"))
  )
  expect_error(
    stancode(
      dpar_form, data = s2z_dat,
      prior = prior(
        normal(0, 1), class = "sd_level", group = "h",
        coef = "z", level = "1"
      )
    ),
    "does not correspond to a coefficient"
  )
})

test_that("legacy prior objects need no level column", {
  form <- y ~ x +
    (1 + x | gr(g, s2z = TRUE, scale = "varying"))
  current <- prior(
    exponential(2), class = "sd", group = "g", coef = "Intercept"
  )
  legacy <- current
  legacy$level <- NULL
  expect_false("level" %in% names(legacy))

  current_code <- stancode(form, data = s2z_dat, prior = current)
  legacy_code <- stancode(form, data = s2z_dat, prior = legacy)
  expect_identical(legacy_code, current_code)
  validated <- validate_prior(legacy, form, data = s2z_dat)
  expect_true("level" %in% names(validated))
  expect_true(all(!nzchar(validated$level)))

  combined <- c(
    legacy,
    prior(
      normal(0.5, 0.2), class = "sd_level", group = "g",
      coef = "Intercept", level = "1"
    )
  )
  expect_true("level" %in% names(combined))
  expect_equal(nrow(combined), 2L)

  obsolete <- prior(
    normal(0.5, 0.2), class = "sd_level", group = "g",
    coef = "Intercept", level = "removed level"
  )
  attr(obsolete, "allow_invalid_prior") <- TRUE
  updated <- validate_prior(obsolete, form, data = s2z_dat)
  expect_false(any(
    updated$class == "sd_level" & nzchar(updated$level)
  ))

  malformed <- set_prior(
    "normal(0.5, 0.2)", class = "sd_level", group = "g",
    coef = "Intercept"
  )
  attr(malformed, "allow_invalid_prior") <- TRUE
  updated <- validate_prior(malformed, form, data = s2z_dat)
  expect_false(any(
    updated$class == "sd_level" & updated$source == "user"
  ))
})

test_that("correlated group-varying scales use the heterogeneous kernel", {
  scode <- stancode(
    y ~ x * z +
      (1 + x * z | gr(g, s2z = TRUE, scale = "varying")),
    data = s2z_dat,
    prior = prior(normal(0, 1), class = b) +
      prior(normal(0, 1.5), class = Intercept)
  )

  for (term in c(
    "vector[M_1 * N_1] z_sd_s2z_1;",
    "matrix<lower=0>[N_1, M_1] sd_level_s2z_1;",
    "vector<lower=0>[M_1] reference_sd_s2z_1;"
  )) {
    expect_match2(scode, term)
  }
  expect_match2(
    scode,
    paste0(
      "reference_sd_s2z_1 = sd_1 .* exp(sdlog_1 .* ",
      "tail(z_sd_s2z_1, M_1) / sqrt(1.0 * N_1));"
    )
  )
  expect_match2(
    scode,
    paste0(
      "sd_level_s2z_1[, k] = reference_sd_s2z_1[k] * ",
      "exp(sdlog_1[k] * z_sd_centered_s2z);"
    )
  )
  expect_match2(
    scode,
    paste0(
      "L_level_s2z = diag_pre_multiply(",
      "sd_level_s2z_1[j]', L_1);"
    )
  )
  expect_match2(
    scode,
    paste0(
      "relative_precision_s2z = mdivide_left_tri_low(",
      "L_level_s2z, L_Sigma_s2z_1);"
    )
  )
  expect_match2(
    scode,
    paste0(
      "P_s2z_1 += crossprod(relative_precision_s2z);"
    )
  )
  expect_match2(
    scode,
    paste0(
      "h_s2z -= relative_precision_s2z' * white_level_s2z;"
    )
  )
  expect_match2(scode, "group_quad_s2z_1 -= dot_self(forward_solve_s2z)")
  expect_match2(
    scode,
    "- (N_1 - 1) * sum(log(diagonal(L_Sigma_s2z_1)))"
  )
  expect_match2(scode, "- sum(log(diagonal(L_P_s2z_1)))")
  expect_match2(scode, "+ 0.5 * M_1 * log(1.0 * N_1)")
  expect_match2(scode, "sd_level_1 = sd_level_s2z_1;")
  expect_match2(scode, "target += std_normal_lpdf(z_sd_s2z_1);")
  expect_equal(
    s2z_count_fixed(
      scode, "target += std_normal_lpdf(z_sd_s2z_1);"
    ),
    1L
  )
  expect_match2(
    scode,
    paste0(
      "sum_to_zero_constrain_brms(segment(z_sd_s2z_1, ",
      "(k - 1) * (N_1 - 1) + 1, N_1 - 1));"
    )
  )
  expect_false(grepl(
    "std_normal_lpdf(z_sd_s2z_1[k])", scode, fixed = TRUE
  ))
  expect_false(grepl("z_sd_mean_s2z_1", scode, fixed = TRUE))
  expect_false(grepl("group_prec_s2z_1", scode, fixed = TRUE))
  expect_false(grepl("relative_sd_s2z_1", scode, fixed = TRUE))
  expect_false(grepl("jacobian", scode, ignore.case = TRUE))
  expect_false(grepl("mdivide_left_spd(", scode, fixed = TRUE))
})

test_that("independent varying scales retain the O(K) specialization", {
  ten_dat <- data.frame(
    y = seq(-1, 1, length.out = 80),
    ten = factor(rep(letters[1:10], 8)),
    g = factor(rep(seq_len(8), each = 10))
  )
  scode <- stancode(
    y ~ 0 + ten +
      (0 + ten || gr(g, s2z = TRUE, scale = "varying")),
    data = ten_dat, prior = prior(normal(0, 2), class = b),
    normalize = FALSE
  )

  expect_match2(
    scode,
    paste0(
      "real relative_precision_s2z = reference_sd_s2z_1[1] / ",
      "sd_level_s2z_1[n, 1];"
    )
  )
  expect_match2(
    scode,
    paste0(
      "D_diag_s2z_1 = group_info_s2z + ",
      "square(reference_sd_s2z_1) .* base_info_s2z;"
    )
  )
  expect_match2(
    scode,
    "- (N_1 - 1) * sum(log(reference_sd_s2z_1))"
  )
  expect_match2(scode, "- 0.5 * sum(log(D_diag_s2z_1))")
  expect_match2(scode, "- 0.5 * log1p(rank1_info_s2z_1)")
  expect_false(grepl("log(1.0 * N_1)", scode, fixed = TRUE))
  expect_false(grepl("group_prec_s2z_1", scode, fixed = TRUE))
  expect_false(grepl("relative_sd_s2z_1", scode, fixed = TRUE))
  for (term in c(
    "matrix[M_1, M_1]", "cholesky_decompose(",
    "mdivide_left_spd(", "mdivide_left_tri_low(",
    "mdivide_right_tri_low("
  )) {
    expect_false(grepl(term, scode, fixed = TRUE), info = term)
  }

  student_code <- stancode(
    y ~ x * z +
      (1 + x * z || gr(
        g, dist = "student", s2z = TRUE, scale = "varying"
      )),
    data = s2z_dat, normalize = FALSE
  )
  expect_match2(student_code, "group_scale_s2z_1 = dfm_1;")
  expect_match2(
    student_code,
    "group_prec_s2z_1 = inv_square(group_scale_s2z_1);"
  )
  expect_match2(
    student_code,
    paste0(
      "real weighted_precision_s2z = group_prec_s2z_1[n] * ",
      "square(relative_precision_s2z);"
    )
  )
  expect_match2(student_code, "- M_1 * sum(log(group_scale_s2z_1))")
})

test_that("varying-scale public names stay separate from kernel internals", {
  form <- y ~ x +
    (1 + x | gr(g, s2z = TRUE, scale = "varying"))
  fit <- brm(form, data = s2z_dat, empty = TRUE)
  excluded <- unlist(brms:::exclude_pars(fit), use.names = FALSE)

  for (name in c(
    "z_sd_s2z_1", "reference_sd_s2z_1", "sd_level_s2z_1"
  )) {
    expect_true(name %in% excluded, info = name)
  }
  # Retain obsolete internal names for old serialized varying-scale fits.
  for (name in c("z_sd_mean_s2z_1", "relative_sd_s2z_1")) {
    expect_true(name %in% excluded, info = name)
  }
  expect_false("sdlog_1" %in% excluded)
  expect_false("sd_level_1" %in% excluded)

  fit_no_group <- brm(
    form, data = s2z_dat, empty = TRUE,
    save_pars = save_pars(group = FALSE)
  )
  excluded_no_group <- unlist(
    brms:::exclude_pars(fit_no_group), use.names = FALSE
  )
  expect_true("sd_level_1" %in% excluded_no_group)
  expect_false("sdlog_1" %in% excluded_no_group)

  bframe <- brms:::brmsframe(brmsterms(form), s2z_dat)
  raw_names <- c(
    "sdlog_1[1]", "sdlog_1[2]",
    sprintf("sd_level_1[%d,1]", seq_len(6)),
    sprintf("sd_level_1[%d,2]", seq_len(6))
  )
  rename_map <- brms:::rename_re(bframe, pars = raw_names)
  renamed <- unlist(lapply(rename_map, `[[`, "fnames"), use.names = FALSE)
  expect_true(all(c(
    "sdlog_g__Intercept", "sdlog_g__x",
    "sd_level_g[1,Intercept]", "sd_level_g[6,Intercept]",
    "sd_level_g[1,x]", "sd_level_g[6,x]"
  ) %in% renamed))

  candidates <- c(
    "b_Intercept", "sd_g__Intercept", "sdlog_g__Intercept",
    "sd_level_g[1,Intercept]", "sd_level_s2z_1[1,1]"
  )
  plot_regex <- brms:::default_plot_variables(gaussian())
  selected <- candidates[vapply(candidates, function(x) {
    any(vapply(plot_regex, grepl, logical(1), x = x))
  }, logical(1))]
  expect_setequal(
    selected, c("b_Intercept", "sd_g__Intercept", "sdlog_g__Intercept")
  )
})

test_that("new Gaussian levels draw fresh group scales", {
  form <- y ~ 1 +
    (1 | gr(g, s2z = TRUE, scale = "varying"))
  bframe <- brms:::brmsframe(brmsterms(form), s2z_dat)
  reframe <- bframe$frame$re
  old_levels <- levels(s2z_dat$g)
  used_levels <- c(old_levels, "new")
  ndraws <- 8L
  baseline_sd <- seq(0.7, 1.4, length.out = ndraws)
  sdlog <- seq(0.1, 0.55, length.out = ndraws)
  draws <- posterior::as_draws_matrix(cbind(
    sd_g__Intercept = baseline_sd,
    sdlog_g__Intercept = sdlog
  ))
  rdraws <- matrix(0, nrow = ndraws, ncol = length(old_levels))

  set.seed(1916)
  scale_z <- rnorm(ndraws)
  effect_z <- rnorm(ndraws)
  expected <- baseline_sd * exp(sdlog * scale_z) * effect_z
  set.seed(1916)
  actual <- brms:::get_new_rdraws(
    reframe = reframe, gf = list(length(old_levels) + 1L),
    rdraws = rdraws, used_levels = used_levels,
    old_levels = old_levels, sample_new_levels = "gaussian",
    draws = draws
  )
  expect_equal(dim(actual), c(ndraws, 1L))
  expect_equal(as.numeric(actual[, 1]), expected, tolerance = 1e-14)

  student_form <- y ~ 1 + (1 | gr(
    g, dist = "student", s2z = TRUE, scale = "varying"
  ))
  student_frame <- brms:::brmsframe(brmsterms(student_form), s2z_dat)
  expect_error(
    brms:::get_new_rdraws(
      reframe = student_frame$frame$re,
      gf = list(length(old_levels) + 1L), rdraws = rdraws,
      used_levels = used_levels, old_levels = old_levels,
      sample_new_levels = "gaussian", draws = draws
    ),
    "not available for non-gaussian"
  )
})

test_that("S2Z handles Student-t effects by conditional Gaussian integration", {
  scode <- stancode(
    y ~ x + (1 + x | gr(g, dist = "student", s2z = TRUE)),
    data = s2z_dat
  )
  expect_match2(scode, "vector<lower=0>[N_1] udf_1;")
  expect_match2(scode, "dfm_1 = sqrt(df_1 * udf_1);")
  expect_match2(scode, "group_scale_s2z_1 = dfm_1;")
  expect_match2(scode, "group_prec_s2z_1 = inv_square(group_scale_s2z_1)")
  expect_match2(
    scode,
    paste0(
      "matrix[M_1, N_1] white_s2z = mdivide_left_tri_low(",
      "L_Sigma_s2z_1, r_s2z_1');"
    )
  )
  expect_match2(
    scode,
    paste0(
      "P_s2z_1 = add_diag(crossprod(prior_factor_s2z), ",
      "sum(group_prec_s2z_1));"
    )
  )
  expect_match2(
    scode,
    paste0(
      "h_s2z = prior_factor_s2z' * prior_difference_s2z - ",
      "contrast_score_s2z_1;"
    )
  )
  expect_match2(
    scode,
    paste0(
      "group_quad_s2z_1 += columns_dot_self(white_s2z) * ",
      "group_prec_s2z_1;"
    )
  )
  expect_match2(scode, "lprior += -0.5 * group_quad_s2z_1")
  expect_match2(scode, "M_1 * sum(log(group_scale_s2z_1))")
  expect_match2(scode, "inv_chi_square_lpdf(udf_1 | df_1)")
  for (term in c(
    "Q_Sigma_s2z_1",
    "L_inv_s2z",
    "mdivide_left_spd(",
    "mdivide_left_tri_low(L_Sigma_s2z_1, diag_matrix",
    "qhat_s2z_1"
  )) {
    expect_false(grepl(term, scode, fixed = TRUE), info = term)
  }
})

test_that("S2Z supports Gaussian scale-mixture population priors", {
  priors <- c(
    prior(student_t(7, 1, 2), class = Intercept),
    prior(cauchy(0, 1), class = b, coef = x)
  )
  scode <- stancode(
    y ~ x + (1 + x | gr(g, s2z = TRUE)),
    data = s2z_dat, prior = priors
  )
  expect_match2(scode, "udf_b_s2z_1;")
  expect_match2(scode, "udf_b_s2z_2;")
  expect_match2(scode, "inv_chi_square_lpdf(udf_b_s2z_1 | 7)")
  expect_match2(scode, "inv_chi_square_lpdf(udf_b_s2z_2 | 1)")
  expect_match2(scode, "normal_lpdf(theta_s2z[1]")
  expect_match2(scode, "normal_lpdf(theta_s2z[2]")
  expect_false(grepl("qhat_s2z_1", scode, fixed = TRUE))

  scode_normal <- stancode(
    y ~ x + (1 + x | gr(g, s2z = TRUE)), data = s2z_dat,
    prior = c(
      prior(normal(0, 2), class = Intercept),
      prior(normal(0, 1), class = b)
    )
  )
  expect_false(grepl("udf_b_s2z", scode_normal, fixed = TRUE))
  expect_match2(scode_normal, "prior_prec_s2z_1[1] = inv_square")
  expect_match2(scode_normal, "prior_prec_s2z_1[2] = inv_square")
})

test_that("S2Z keeps all parameter-dependent normalizers", {
  sc_norm <- stancode(
    y ~ x + (1 + x | gr(g, s2z = TRUE)),
    data = s2z_dat, normalize = TRUE
  )
  sc_kernel <- stancode(
    y ~ x + (1 + x | gr(g, s2z = TRUE)),
    data = s2z_dat, normalize = FALSE
  )

  for (term in c(
    "(N_1 - 1) * sum(log(diagonal(L_Sigma_s2z_1)))",
    "sum(log(diagonal(L_P_s2z_1)))"
  )) {
    expect_true(grepl(term, sc_norm, fixed = TRUE))
    expect_true(grepl(term, sc_kernel, fixed = TRUE))
  }
  expect_match2(sc_norm, "0.5 * M_1 * log(1.0 * N_1)")
  expect_false(grepl("log(1.0 * N_1)", sc_kernel, fixed = TRUE))
})

test_that("s2z FALSE leaves conventional code unchanged", {
  implicit <- stancode(y ~ x + (1 + x | g), data = s2z_dat)
  explicit <- stancode(
    y ~ x + (1 + x | gr(g, s2z = FALSE)), data = s2z_dat
  )
  expect_identical(implicit, explicit)
})

test_that("S2Z blocks can belong to distinct distributional predictors", {
  form <- bf(
    y ~ x + (1 + x | gr(g, s2z = TRUE)),
    sigma ~ z + (1 + z | gr(h, s2z = TRUE))
  )
  scode <- stancode(form, data = s2z_dat)

  expect_match2(scode, "vector[2] theta_s2z;")
  expect_match2(scode, "vector[2] theta_s2z_sigma;")
  expect_match2(scode, "matrix[N_1, M_1] r_s2z_1;")
  expect_match2(scode, "matrix[N_2, M_2] r_s2z_2;")
  expect_equal(
    lengths(regmatches(
      scode, gregexpr("vector sum_to_zero_constrain_brms", scode, fixed = TRUE)
    )),
    1L
  )
})

test_that("mixed S2Z IDs cannot span linear predictors in either term order", {
  mixed_ids <- list(
    bf(
      y ~ 1 + (1 | gr(g, id = "across", s2z = TRUE)),
      sigma ~ 1 + (1 | gr(g, id = "across", s2z = FALSE))
    ),
    bf(
      y ~ 1 + (1 | gr(g, id = "across", s2z = FALSE)),
      sigma ~ 1 + (1 | gr(g, id = "across", s2z = TRUE))
    )
  )
  for (form in mixed_ids) {
    expect_error(
      stancode(form, data = s2z_dat),
      "cannot span multiple linear predictors"
    )
  }
})

test_that("a conventional correlated S2Z ID spans response predictors", {
  mv_data <- transform(
    s2z_dat,
    phen = y,
    cofactor = 0.4 * y + sin(seq_len(nrow(s2z_dat)) * 0.31),
    phylo = g
  )
  form <- bf(
    mvbind(phen, cofactor) ~ 1 + (1 | q | gr(phylo, s2z = TRUE))
  )
  cross_prior <-
    prior(normal(0, 1.25), class = "Intercept", resp = "phen") +
    prior(normal(0.5, 2), class = "Intercept", resp = "cofactor") +
    prior(exponential(2), class = "sd", group = "phylo", resp = "phen") +
    prior(exponential(3), class = "sd", group = "phylo", resp = "cofactor") +
    prior(lkj(3), class = "cor", group = "phylo")
  scode <- stancode(form, data = mv_data, prior = cross_prior)

  expect_match2(scode, "vector[1] theta_s2z_phen;")
  expect_match2(scode, "vector[1] theta_s2z_cofactor;")
  expect_match2(scode, "cholesky_factor_corr[M_1] L_1;")
  expect_match2(
    scode,
    "L_Sigma_s2z_1 = diag_pre_multiply(sd_1, L_1);"
  )
  expect_match2(
    scode,
    "vector[M_1 * (N_1 - 1)] z_s2z_1;"
  )
  expect_match2(
    scode,
    "real<lower=0> group_info_s2z_1;"
  )
  expect_match2(
    scode, "group_info_s2z_1 = dot_self(one_white_cov_s2z);"
  )
  expect_match2(
    scode,
    paste0(
      "P_s2z_1 = add_diag(crossprod(prior_factor_s2z), ",
      "group_info_s2z_1);"
    )
  )
  expect_false(grepl("P_group_s2z_1", scode, fixed = TRUE))
  expect_match2(
    scode, "group_quad_s2z_1 -= dot_self(whitened_h_s2z);"
  )
  expect_match2(
    scode,
    "- (N_1 - 1) * sum(log(diagonal(L_Sigma_s2z_1)))"
  )
  expect_match2(scode, "- sum(log(diagonal(L_P_s2z_1)))")
  expect_match2(scode, "+ 0.5 * M_1 * log(1.0 * N_1)")
  expect_match2(scode, "matrix[2, M_1] H_s2z_1;")
  expect_match2(scode, "H_s2z_1[1, 1] = 1.0;")
  expect_match2(scode, "H_s2z_1[2, 2] = 1.0;")
  expect_match2(scode, "r_s2z_1_phen_1 = r_s2z_1[, 1];")
  expect_match2(scode, "r_s2z_1_cofactor_2 = r_s2z_1[, 2];")
  expect_match2(scode, "mu_phen += theta_s2z_phen[1];")
  expect_match2(scode, "mu_cofactor += theta_s2z_cofactor[1];")
  expect_match2(
    scode,
    "q_recovered_s2z_1 -= H_s2z_1 * mean_r_s2z_1;"
  )
  expect_match2(scode, "b_phen_Intercept = Intercept_phen;")
  expect_match2(scode, "b_cofactor_Intercept = Intercept_cofactor;")
  expect_match2(scode, "r_1_phen_1 = r_1[, 1];")
  expect_match2(scode, "r_1_cofactor_2 = r_1[, 2];")
  expect_match2(
    scode,
    "corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);"
  )
  # Priors retain their conventional response-specific coefficient, scale,
  # and cross-response correlation meanings under the coordinate map.
  expect_match2(
    scode, "lprior += normal_lpdf(theta_s2z_phen[1] | 0, 1.25);"
  )
  expect_match2(
    scode, "lprior += normal_lpdf(theta_s2z_cofactor[1] | 0.5, 2);"
  )
  expect_match2(scode, "lprior += exponential_lpdf(sd_1[1] | 2);")
  expect_match2(scode, "lprior += exponential_lpdf(sd_1[2] | 3);")
  expect_match2(scode, "lprior += lkj_corr_cholesky_lpdf(L_1 | 3);")

  mixture_code <- stancode(
    form, data = mv_data,
    prior =
      prior(
        student_t(7, 0.2, 1.3), class = "Intercept", resp = "phen"
      ) +
      prior(
        cauchy(-0.1, 0.8), class = "Intercept", resp = "cofactor"
      )
  )
  for (term in c(
    "real<lower=0> udf_b_s2z_phen_1;",
    "real<lower=0> udf_b_s2z_cofactor_1;",
    "inv_chi_square_lpdf(udf_b_s2z_phen_1 | 7)",
    "inv_chi_square_lpdf(udf_b_s2z_cofactor_1 | 1)",
    "normal_lpdf(theta_s2z_phen[1]",
    "normal_lpdf(theta_s2z_cofactor[1]"
  )) {
    expect_match2(mixture_code, term)
  }

  kernel_code <- stancode(
    form, data = mv_data, prior = cross_prior, normalize = FALSE
  )
  for (term in c(
    "- (N_1 - 1) * sum(log(diagonal(L_Sigma_s2z_1)))",
    "- sum(log(diagonal(L_P_s2z_1)))",
    "lkj_corr_cholesky_lupdf(L_1 | 3)"
  )) {
    expect_match2(kernel_code, term)
  }
  expect_false(grepl(
    "0.5 * M_1 * log(1.0 * N_1)", kernel_code, fixed = TRUE
  ))

  sdata <- standata(form, data = mv_data, prior = cross_prior)
  expect_equal(sdata$M_1, 2L)
  expect_equal(sdata$NC_1, 1L)

  noncentered_code <- stancode(
    bf(
      mvbind(phen, cofactor) ~ 1 +
        (1 | q | gr(phylo, s2z = TRUE, center = FALSE))
    ),
    data = mv_data
  )
  expect_match2(
    noncentered_code,
    "r_s2z_1 = r_s2z_1 * L_Sigma_s2z_1';"
  )
  expect_false(grepl(
    "- (N_1 - 1) * sum(log(diagonal(L_Sigma_s2z_1)))",
    noncentered_code, fixed = TRUE
  ))
  expect_error(
    stancode(
      bf(
        y ~ 1 + (1 | gr(g, id = "across", s2z = TRUE)),
        sigma ~ 1 + (1 | gr(g, id = "across", s2z = TRUE))
      ),
      data = s2z_dat
    ),
    "multivariate response-location predictors"
  )
})

test_that("cross-response S2Z composes charts, scales, Student, and cov", {
  mv_data <- transform(
    s2z_dat,
    phen = y,
    cofactor = 0.4 * y + sin(seq_len(nrow(s2z_dat)) * 0.31),
    phylo = g
  )
  levels_g <- levels(mv_data$phylo)
  Omega <- outer(seq_along(levels_g), seq_along(levels_g), function(i, j) {
    (-0.3)^abs(i - j)
  })
  dimnames(Omega) <- list(rev(levels_g), rev(levels_g))
  form <-
    bf(
      phen ~ 1 + (1 | q | gr(
        phylo, s2z = TRUE, center = 0.2, scale = "varying",
        dist = "student", cov = Omega
      ))
    ) +
    bf(
      cofactor ~ 1 + (1 | q | gr(
        phylo, s2z = TRUE, center = 0.8, scale = "varying",
        dist = "student", cov = Omega
      ))
    ) +
    set_rescor(FALSE)
  cross_prior <-
    prior(
      normal(0.4, 0.2), class = "sd_level", group = "phylo",
      coef = "Intercept", level = "1", resp = "phen"
    ) +
    prior(
      normal(0.7, 0.2), class = "sd_level", group = "phylo",
      coef = "Intercept", level = "1", resp = "cofactor"
    )
  scode <- stancode(
    form, data = mv_data, data2 = list(Omega = Omega), prior = cross_prior
  )

  for (term in c(
    "matrix<lower=0,upper=1>[N_1, M_1] rho_s2z_1;",
    "vector[M_1 * N_1] z_sd_s2z_1;",
    "reference_sd_s2z_1 = sd_1 .* exp(sdlog_1 .*",
    "group_scale_s2z_1 = dfm_1;",
    "L_Sigma_s2z_1 = diag_pre_multiply(reference_sd_s2z_1, L_1);",
    "matrix[N_1 * M_1, M_1] mean_factor_cov_s2z;",
    "P_group_s2z_1 = crossprod(mean_factor_cov_s2z);",
    "+ log_det_partial_s2z_1",
    "- M_1 * sum(log(group_scale_s2z_1))",
    "- M_1 * sum(log(diagonal(Lcov_1)))",
    "sd_level_1 = sd_level_s2z_1;",
    "sd_level_s2z_1[1, 1] | 0.4, 0.2",
    "sd_level_s2z_1[1, 2] | 0.7, 0.2"
  )) {
    expect_match2(scode, term)
  }

  sdata <- standata(
    form, data = mv_data, data2 = list(Omega = Omega), prior = cross_prior
  )
  expect_equal(dim(sdata$rho_s2z_1), c(length(levels_g), 2L))
  expect_equal(
    unname(sdata$rho_s2z_1),
    cbind(rep(0.2, length(levels_g)), rep(0.8, length(levels_g)))
  )
  expect_equal(
    unname(sdata$Lcov_1),
    unname(t(chol(Omega[levels_g, levels_g])))
  )

  noncentered_code <- stancode(
    bf(mvbind(phen, cofactor) ~ 1 + (1 | q | gr(
      phylo, s2z = TRUE, center = FALSE, dist = "student", cov = Omega
    ))) + set_rescor(FALSE),
    data = mv_data, data2 = list(Omega = Omega)
  )
  expect_match2(
    noncentered_code,
    "r_s2z_1 = r_s2z_1 * L_Sigma_s2z_1';"
  )
  expect_match2(noncentered_code, "one_white_cov_s2z")
  expect_false(grepl(
    "- (N_1 - 1) * sum(log(diagonal(L_Sigma_s2z_1)))",
    noncentered_code, fixed = TRUE
  ))

  fisher_independent <- stancode(
    bf(mvbind(phen, cofactor) ~ 1 + (1 | q | gr(
      phylo, s2z = TRUE, center = "fisher", scale = "varying",
      dist = "student", cov = Omega
    ))) + set_rescor(FALSE),
    data = mv_data, data2 = list(Omega = Omega)
  )
  for (term in c(
    "matrix<lower=0,upper=1>[N_1, M_1] rho_s2z_1;",
    "info_fisher_s2z[J_1_phen[n]] += obs_prec_fisher_s2z",
    "info_fisher_s2z[J_1_cofactor[n]] += obs_prec_fisher_s2z",
    "row_var_fisher_s2z_1[j] * quad_form(info_fisher_s2z[j]",
    "group_scale_s2z_1 = dfm_1;",
    "+ log_det_partial_s2z_1",
    "- M_1 * sum(log(diagonal(Lcov_1)))"
  )) {
    expect_match2(fisher_independent, term)
  }

  fisher_rescor <- stancode(
    bf(mvbind(phen, cofactor) ~ 1 + (1 | q | gr(
      phylo, s2z = TRUE, center = "fisher", cov = Omega
    ))) + set_rescor(TRUE),
    data = mv_data, data2 = list(Omega = Omega)
  )
  for (term in c(
    "matrix[N_1, N_1 - 1] Ecov_s2z_1;",
    "vector<lower=0>[N_1 - 1] lambda_cov_s2z_1;",
    "vector[N_1 - 1] coupling_cov_s2z_1;",
    "real<lower=0> mean_prec_cov_s2z_1;",
    "exposure_modal_fisher_s2z_1[ell] = dot_product(",
    "L_residual_modal_fisher_s2z = diag_pre_multiply(",
    "white_design_modal_fisher_s2z = mdivide_left_tri_low(",
    "info_modal_fisher_s2z = crossprod(",
    "unit_info_modal_fisher_s2z = quad_form(",
    "conditional_precision_modal_fisher_s2z =",
    "rho_s2z_1[ell, 1] = fmin(1.0, fmax(0.0,",
    "denominator11_modal_s2z =",
    "r_s2z_1 = Ecov_s2z_1 * effect_modal_s2z;",
    "coupling_score_modal_s2z = coef_white_modal_s2z' *",
    "group_info_s2z_1 = mean_prec_cov_s2z_1;",
    paste0(
      "P_s2z_1 = add_diag(crossprod(prior_factor_s2z), ",
      "group_info_s2z_1);"
    ),
    "group_quad_s2z_1 = dot_product(",
    "+ log_det_partial_s2z_1"
  )) {
    expect_match2(fisher_rescor, term)
  }
  fisher_rescor_default <- suppressWarnings(stancode(
    bf(mvbind(phen, cofactor) ~ 1 + (1 | q | gr(
      phylo, s2z = TRUE, center = "fisher", cov = Omega
    ))),
    data = mv_data, data2 = list(Omega = Omega)
  ))
  expect_identical(fisher_rescor_default, fisher_rescor)
  expect_match2(
    fisher_rescor_default,
    "real conditional_unit_info_modal_fisher_s2z = fmax(0.0,"
  )
  for (term in c(
    "info_fisher_s2z[J_1_phen[n]]",
    "row_var_fisher_s2z_1",
    "one_white_cov_s2z",
    "mean_factor_cov_s2z",
    "P_group_s2z_1",
    "L_post_modal_fisher_s2z"
  )) {
    expect_false(grepl(term, fisher_rescor, fixed = TRUE))
  }
  fisher_data <- standata(
    bf(mvbind(phen, cofactor) ~ 1 + (1 | q | gr(
      phylo, s2z = TRUE, center = "fisher", cov = Omega
    ))) + set_rescor(TRUE),
    data = mv_data, data2 = list(Omega = Omega)
  )
  expect_false("rho_s2z_1" %in% names(fisher_data))

  G <- length(levels_g)
  one <- rep(1, G)
  P_s2z <- diag(G) - matrix(1 / G, G, G)
  Omega_ordered <- Omega[levels_g, levels_g, drop = FALSE]
  E <- unname(fisher_data$Ecov_s2z_1)
  lambda <- unname(fisher_data$lambda_cov_s2z_1)
  coupling <- unname(fisher_data$coupling_cov_s2z_1)
  kappa <- unname(fisher_data$mean_prec_cov_s2z_1)
  Omega_inv_one <- drop(solve(Omega_ordered, one))

  expect_equal(dim(E), c(G, G - 1L))
  expect_equal(length(lambda), G - 1L)
  expect_equal(crossprod(E), diag(G - 1L), tolerance = 1e-13)
  expect_equal(drop(crossprod(E, one)), numeric(G - 1L), tolerance = 1e-13)
  expect_true(all(lambda > 0))
  expect_equal(
    E %*% diag(lambda) %*% t(E),
    P_s2z %*% Omega_ordered %*% P_s2z,
    tolerance = 2e-13
  )
  expect_equal(
    coupling, drop(crossprod(E, Omega_inv_one)), tolerance = 2e-13
  )
  expect_equal(kappa, drop(crossprod(one, Omega_inv_one)), tolerance = 1e-13)
  expect_equal(
    crossprod(E, solve(Omega_ordered, E)),
    diag(1 / lambda) + tcrossprod(coupling) / kappa,
    tolerance = 3e-13
  )

  # Each condition excluded by the spectral dispatch retains the established
  # levelwise Fisher path: independent residuals, varying scales, or Student-t
  # group effects.
  fisher_fallback <- list(
    independent = stancode(
      bf(mvbind(phen, cofactor) ~ 1 + (1 | q | gr(
        phylo, s2z = TRUE, center = "fisher", cov = Omega
      ))) + set_rescor(FALSE),
      data = mv_data, data2 = list(Omega = Omega)
    ),
    varying = stancode(
      bf(mvbind(phen, cofactor) ~ 1 + (1 | q | gr(
        phylo, s2z = TRUE, center = "fisher", scale = "varying",
        cov = Omega
      ))) + set_rescor(TRUE),
      data = mv_data, data2 = list(Omega = Omega)
    ),
    student = stancode(
      bf(mvbind(phen, cofactor) ~ 1 + (1 | q | gr(
        phylo, s2z = TRUE, center = "fisher", dist = "student",
        cov = Omega
      ))) + set_rescor(TRUE),
      data = mv_data, data2 = list(Omega = Omega)
    )
  )
  for (fallback_code in fisher_fallback) {
    expect_false(grepl("modal_fisher_s2z", fallback_code, fixed = TRUE))
    expect_false(grepl("Ecov_s2z_1", fallback_code, fixed = TRUE))
    expect_match2(fallback_code, "row_var_fisher_s2z_1[j] * quad_form(")
    expect_match2(fallback_code, "+ log_det_partial_s2z_1")
  }
  independent_data <- standata(
    bf(mvbind(phen, cofactor) ~ 1 + (1 | q | gr(
      phylo, s2z = TRUE, center = "fisher", cov = Omega
    ))) + set_rescor(FALSE),
    data = mv_data, data2 = list(Omega = Omega)
  )
  expect_false(any(c(
    "Ecov_s2z_1", "lambda_cov_s2z_1", "coupling_cov_s2z_1",
    "mean_prec_cov_s2z_1"
  ) %in% names(independent_data)))

  non_gaussian_data <- transform(
    mv_data,
    yb = as.integer(phen > median(phen)),
    yp = as.integer(rank(cofactor) %% 5L)
  )
  fisher_non_gaussian <- stancode(
    bf(yb ~ 1 + (1 | q | gr(
      phylo, s2z = TRUE, center = "fisher", scale = "varying",
      dist = "student", cov = Omega
    )), family = bernoulli()) +
      bf(yp ~ 1 + (1 | q | gr(
        phylo, s2z = TRUE, center = "fisher", scale = "varying",
        dist = "student", cov = Omega
      )), family = poisson()) +
      set_rescor(FALSE),
    data = non_gaussian_data, data2 = list(Omega = Omega)
  )
  for (term in c(
    "inv_logit(eta_fisher_s2z_mu)",
    "exp(eta_fisher_s2z_mu)",
    "info_fisher_s2z[J_1_yb[n]]",
    "info_fisher_s2z[J_1_yp[n]]",
    "group_scale_s2z_1 = dfm_1;",
    "row_var_fisher_s2z_1"
  )) {
    expect_match2(fisher_non_gaussian, term)
  }
})

test_that("known group covariance composes with conventional S2Z charts", {
  levels_g <- levels(s2z_dat$g)
  G <- length(levels_g)
  P_s2z <- diag(G) - matrix(1 / G, G, G)
  expect_equal(
    diag(P_s2z %*% diag(G) %*% P_s2z) / (1 - 1 / G),
    rep(1, G), tolerance = 1e-14
  )
  rho <- -0.35
  Omega <- outer(seq_along(levels_g), seq_along(levels_g), function(i, j) {
    rho^abs(i - j)
  })
  dimnames(Omega) <- list(rev(levels_g), rev(levels_g))
  data2 <- list(Omega = Omega)

  forms <- list(
    centered = y ~ x + (1 | gr(g, s2z = TRUE, cov = Omega)),
    noncentered = y ~ x + (1 | gr(
      g, s2z = TRUE, center = FALSE, cov = Omega
    )),
    partial = y ~ x + (1 | gr(
      g, s2z = TRUE, center = 0.4, cov = Omega
    )),
    fisher = y ~ x + (1 | gr(
      g, s2z = TRUE, center = "fisher", cov = Omega
    )),
    student_partial = y ~ x + (1 + x | gr(
      g, s2z = TRUE, center = 0.35, dist = "student", cov = Omega
    )),
    student_fisher = y ~ x + (1 + x | gr(
      g, s2z = TRUE, center = "fisher", dist = "student", cov = Omega
    )),
    varying_fisher = y ~ x + (1 + x || gr(
      g, s2z = TRUE, center = "fisher", scale = "varying", cov = Omega
    )),
    varying_student_fisher = y ~ x + (1 + x | gr(
      g, s2z = TRUE, center = "fisher", scale = "varying",
      dist = "student", cov = Omega
    ))
  )
  for (i in seq_along(forms)) {
    form <- forms[[i]]
    scode <- stancode(form, data = s2z_dat, data2 = data2)
    for (term in c(
      "matrix[N_1, N_1] Lcov_1;",
      "h_group_s2z_1",
      "- M_1 * sum(log(diagonal(Lcov_1)))"
    )) {
      expect_match2(scode, term)
    }
    if (!startsWith(names(forms)[i], "varying")) {
      expect_match2(scode, "one_white_cov_s2z")
      expect_match2(scode, paste0(
        "group_info_s2z_1 = dot_self(one_white_cov_s2z);"
      ))
      expect_match2(scode, paste0(
        "P_s2z_1 = add_diag(crossprod(prior_factor_s2z), ",
        "group_info_s2z_1);"
      ))
      expect_false(grepl("P_group_s2z_1", scode, fixed = TRUE))
    } else if (names(forms)[i] == "varying_fisher") {
      expect_match2(
        scode, "vector<lower=0>[M_1] group_info_s2z_1;"
      )
      expect_match2(
        scode, "group_info_s2z_1[k] = dot_self(white_basis_cov_s2z);"
      )
      expect_match2(scode, paste0(
        "P_s2z_1 = add_diag(crossprod(prior_factor_s2z), ",
        "group_info_s2z_1);"
      ))
      expect_false(grepl("P_group_s2z_1", scode, fixed = TRUE))
      expect_false(grepl("mean_factor_cov_s2z", scode, fixed = TRUE))
    } else {
      expect_match2(
        scode, "matrix[N_1 * M_1, M_1] mean_factor_cov_s2z;"
      )
      expect_match2(
        scode, "P_group_s2z_1 = crossprod(mean_factor_cov_s2z);"
      )
      expect_match2(scode, "h_group_s2z_1 = -mean_factor_cov_s2z'")
    }
  }

  fisher_code <- stancode(forms[["fisher"]], data = s2z_dat, data2 = data2)
  expect_match2(fisher_code, "row_var_fisher_s2z_1")
  expect_match2(fisher_code, "rows_dot_self(centered_Lcov_s2z)")
  expect_match2(fisher_code, "N_1 / (N_1 - 1.0)")
  kernel_code <- stancode(
    forms[["centered"]], data = s2z_dat, data2 = data2, normalize = FALSE
  )
  expect_false(grepl(
    "sum(log(diagonal(Lcov_1)))", kernel_code, fixed = TRUE
  ))

  sdata <- standata(forms[["centered"]], data = s2z_dat, data2 = data2)
  expect_equal(
    unname(sdata$Lcov_1),
    unname(t(chol(Omega[levels_g, levels_g])))
  )

  joint_form <- y ~ x +
    (1 | gr(g, id = "a", s2z = TRUE, cov = Omega)) +
    (0 + x | gr(g, id = "b", s2z = TRUE, cov = Omega))
  joint_code <- stancode(joint_form, data = s2z_dat, data2 = data2)
  expect_false(grepl("Matheron system", joint_code, fixed = TRUE))
  expect_match2(joint_code, "group_info_s2z_1")
  expect_match2(joint_code, "group_info_s2z_2")
  expect_false(grepl("P_group_s2z_", joint_code, fixed = TRUE))

  mv_data <- transform(
    s2z_dat, phen = y, cofactor = 0.4 * y + sin(seq_len(nrow(s2z_dat)))
  )
  mv_code <- stancode(
    bf(mvbind(phen, cofactor) ~ 1 + (1 | q | gr(
      g, s2z = TRUE, cov = Omega
    ))),
    data = mv_data, data2 = data2
  )
  expect_match2(mv_code, "one_white_cov_s2z")
  expect_match2(mv_code, "h_group_s2z_1")
  expect_match2(mv_code, "q_recovered_s2z_1 -= H_s2z_1")

  for (chunk in c(
    "fun_scale_r_cor_cov.stan", "fun_scale_r_cor_by_cov.stan"
  )) {
    source <- readLines(system.file("chunks", chunk, package = "brms"))
    expect_true(any(grepl("Lcov[icov, jcov] != 0", source, fixed = TRUE)))
    expect_false(any(grepl("Lcov[icov, jcov] > 1e-10", source, fixed = TRUE)))
  }
})

test_that("large crossed scalar Gaussian systems use one-dimensional Matheron work", {
  n <- 120L
  many_data <- data.frame(y = sin(seq_len(n) * 0.07))
  n_factor <- 10L
  for (k in seq_len(n_factor)) {
    many_data[[paste0("g", k)]] <- factor(
      (seq_len(n) - 1L) %% (k + 1L) + 1L
    )
  }
  terms <- sprintf(
    "(1 | gr(g%s, s2z = TRUE, center = FALSE))",
    seq_len(n_factor)
  )
  form <- stats::as.formula(paste("y ~ 1 +", paste(terms, collapse = "+")))
  scode <- stancode(
    form, data = many_data,
    prior = prior(normal(0, 1.4), class = Intercept)
  )

  expect_match2(scode, "// fast Gaussian Matheron system for S2Z blocks")
  expect_match2(scode, "real<lower=0> W_matheron_s2z_1;")
  expect_match2(scode, "real<lower=0> sqrt_W_matheron_s2z_1;")
  expect_equal(
    s2z_count_fixed(scode, "W_matheron_s2z_1 += dot_self("),
    n_factor
  )
  expect_equal(
    s2z_count_fixed(scode, "theta_star_s2z[1] += dot_product("),
    n_factor
  )
  expect_equal(s2z_count_fixed(scode, "cholesky_decompose("), 0L)
  expect_false(grepl("matrix[10, 10] P_s2z_1", scode, fixed = TRUE))
  expect_false(grepl("L_P_s2z_1", scode, fixed = TRUE))
  expect_false(grepl("P_group_s2z_", scode, fixed = TRUE))
  expect_match2(
    scode,
    "delta_s2z[1] = (theta_s2z[1] - theta_star_s2z[1]) / "
  )
  expect_match2(scode, "q_recovered_s2z_1 = theta_s2z;")
})

test_that("Matheron supports overlapping blocks and all centering modes", {
  form <- y ~ x * z +
    (1 + x | gr(g, s2z = TRUE, center = TRUE)) +
    (1 + x * z || gr(h, s2z = TRUE, center = FALSE)) +
    (1 | gr(f, s2z = TRUE, center = 0.4))
  bprior <- prior(normal(0, 2), class = Intercept) +
    prior(normal(0, 1), class = b)
  scode <- stancode(form, data = s2z_dat, prior = bprior)

  for (term in c(
    "// fast Gaussian Matheron system for S2Z blocks 1, 2, 3",
    "matrix[4, 4] W_matheron_s2z_1;",
    "matrix[4, 4] L_W_matheron_s2z_1;",
    "W_matheron_s2z_1 = add_diag(",
    "square(prior_scale_s2z_1[{1, 2, 3, 4}])",
    "H_active_s2z = H_s2z_1[{1, 2, 3, 4}, ];",
    paste0(
      "theta_difference_s2z = theta_s2z[{1, 2, 3, 4}] - ",
      "prior_mean_s2z_1[{1, 2, 3, 4}];"
    ),
    "- (N_2 - 1) * sum(log(diagonal(L_Sigma_s2z_2)))",
    "+ log_det_partial_s2z_1",
    "- 0.5 * (N_1 - 1) * M_1 * log(2 * pi())",
    "- 0.5 * (N_2 - 1) * M_2 * log(2 * pi())",
    "- 0.5 * (N_3 - 1) * M_3 * log(2 * pi())",
    "- 0.5 * 4 * log(2 * pi())",
    "mean_r_s2z_1 += L_Sigma_s2z_1 *",
    "mean_r_s2z_2 += L_Sigma_s2z_2 *",
    "mean_r_s2z_3 += L_Sigma_s2z_3 *",
    "group_quad_s2z_3 = dot_self(z_s2z_3);",
    "q_recovered_s2z_1 = theta_s2z;"
  )) {
    expect_true(grepl(term, scode, fixed = TRUE), info = term)
  }
  expect_false(grepl(
    "- (N_3 - 1) * sum(log(diagonal(L_Sigma_s2z_3)))",
    scode, fixed = TRUE
  ))
  expect_false(grepl("log_det_partial_s2z_3", scode, fixed = TRUE))
  expect_false(grepl(
    "r_s2z_3 ./ rep_matrix", scode, fixed = TRUE
  ))
  expect_false(grepl(
    "L_Sigma_s2z_3, r_s2z_3", scode, fixed = TRUE
  ))
  expect_false(grepl("matrix[7, 7] P_s2z_1", scode, fixed = TRUE))
  expect_false(grepl("L_P_s2z_1", scode, fixed = TRUE))
  expect_false(grepl("H_joint_s2z_1", scode, fixed = TRUE))
  expect_false(grepl(
    "W_matheron_s2z_1 = rep_matrix", scode, fixed = TRUE
  ))
  expect_false(grepl(
    "diag_matrix(square(prior_scale_s2z_1", scode, fixed = TRUE
  ))
  expect_equal(
    s2z_count_fixed(
      scode,
      "W_matheron_s2z_1 += tcrossprod(H_active_s2z * "
    ),
    2L
  )
  expect_equal(s2z_count_fixed(scode, "cholesky_decompose("), 1L)

  selective_prior <- prior(normal(0, 2), class = Intercept) +
    prior(normal(0, 1), class = b, coef = "x:z")
  selective_code <- stancode(form, data = s2z_dat, prior = selective_prior)
  for (term in c(
    "W_matheron_s2z_1 = add_diag(",
    "square(prior_scale_s2z_1[{1, 4}])",
    "H_active_s2z = H_s2z_1[{1, 4}, ];",
    paste0(
      "theta_difference_s2z = theta_s2z[{1, 4}] - ",
      "prior_mean_s2z_1[{1, 4}];"
    )
  )) {
    expect_true(grepl(term, selective_code, fixed = TRUE), info = term)
  }
})

test_that("conditional Student and Cauchy population priors use Matheron", {
  form <- y ~ x +
    (1 + x || gr(g, s2z = TRUE)) +
    (1 + x | gr(h, s2z = TRUE, center = FALSE))
  bprior <- prior(student_t(5, 0.2, 1.7), class = Intercept) +
    prior(cauchy(-0.1, 0.8), class = b)
  scode <- stancode(form, data = s2z_dat, prior = bprior)

  expect_match2(scode, "// fast Gaussian Matheron system")
  expect_match2(scode, "real<lower=0> udf_b_s2z_1;")
  expect_match2(scode, "real<lower=0> udf_b_s2z_2;")
  expect_match2(scode, "inv_chi_square_lpdf(udf_b_s2z_1 | 5)")
  expect_match2(scode, "inv_chi_square_lpdf(udf_b_s2z_2 | 1)")
  expect_match2(scode, "prior_scale_s2z_1[1] = 1.7 * sqrt(5 *")
  expect_match2(scode, "prior_scale_s2z_1[2] = ")
  expect_match2(scode, "sqrt(1 * udf_b_s2z_2)")
  expect_match2(scode, "matrix[2, 2] W_matheron_s2z_1;")
  expect_match2(scode, "matrix[2, 2] L_W_matheron_s2z_1;")
  expect_false(grepl("matrix[4, 4] P_s2z_1", scode, fixed = TRUE))
})

test_that("proper population coordinates outside group maps score independently", {
  form <- y ~ x + z +
    (1 | gr(g, s2z = TRUE)) +
    (1 | gr(h, s2z = TRUE, center = FALSE))
  bprior <- prior(normal(0.2, 1.5), class = Intercept) +
    prior(normal(-0.1, 0.7), class = b)
  scode <- stancode(form, data = s2z_dat, prior = bprior)

  expect_match2(scode, "// fast Gaussian Matheron system")
  expect_match2(scode, "real<lower=0> W_matheron_s2z_1;")
  expect_false(grepl("matrix[2, 2] P_s2z_1", scode, fixed = TRUE))
  expect_equal(
    s2z_count_fixed(
      scode,
      "normal_lpdf(fixed_s2z | -0.1, 0.7)"
    ),
    1L
  )
  expect_match2(scode, "theta_s2z[2] = fixed_s2z[1];")
  expect_match2(scode, "theta_s2z[3] = fixed_s2z[2];")
  expect_false(grepl("normal_lpdf(theta_s2z[1]", scode, fixed = TRUE))
  expect_match2(scode, "theta_s2z[1] - prior_mean_s2z_1[1]")
})

test_that("flat active population coordinates use the zero-rank fast path", {
  form <- y ~ 0 + x +
    (0 + x | gr(g, s2z = TRUE)) +
    (0 + x | gr(h, s2z = TRUE))
  scode <- stancode(form, data = s2z_dat)

  expect_match2(scode, "// fast Gaussian Matheron system")
  expect_false(grepl("W_matheron_s2z_1", scode, fixed = TRUE))
  expect_false(grepl("theta_star_s2z", scode, fixed = TRUE))
  expect_false(grepl("delta_s2z", scode, fixed = TRUE))
  expect_false(grepl("P_s2z_1", scode, fixed = TRUE))
  expect_match2(scode, "q_recovered_s2z_1 = theta_s2z;")
  expect_match2(
    scode,
    "q_recovered_s2z_1 -= H_s2z_1 * mean_r_s2z_1;"
  )
  expect_match2(
    scode,
    "q_recovered_s2z_1 -= H_s2z_2 * mean_r_s2z_2;"
  )
})

test_that("nonseparable or nonbeneficial multiblock models keep dense fallback", {
  cases <- list(
    list(
      form = y ~ 0 + x + z +
        (0 + x | gr(g, s2z = TRUE)) +
        (0 + z | gr(h, s2z = TRUE)),
      prior = prior(normal(0, 1), class = b), total = 2L
    ),
    list(
      form = y ~ x +
        (1 + x | gr(g, s2z = TRUE, dist = "student")) +
        (1 + x | gr(h, s2z = TRUE)),
      prior = prior(normal(0, 2), class = Intercept) +
        prior(normal(0, 1), class = b), total = 4L
    ),
    list(
      form = y ~ x +
        (1 + x || gr(g, s2z = TRUE, scale = "varying")) +
        (1 + x | gr(h, s2z = TRUE)),
      prior = prior(normal(0, 2), class = Intercept) +
        prior(normal(0, 1), class = b), total = 4L
    )
  )
  for (case in cases) {
    scode <- stancode(case$form, data = s2z_dat, prior = case$prior)
    expect_false(
      grepl("fast Gaussian Matheron system", scode, fixed = TRUE),
      info = deparse0(case$form)
    )
    expect_match2(
      scode, sprintf("matrix[%s, %s] P_s2z_1;", case$total, case$total)
    )
    expect_match2(scode, "L_P_s2z_1 = cholesky_decompose(P_s2z_1);")
  }
})

test_that("joint Student varying scales retain precision weights", {
  form <- y ~ x +
    (1 + x || gr(
      g, s2z = TRUE, scale = "varying", dist = "student"
    )) +
    (1 + x | gr(h, s2z = TRUE))
  bprior <- prior(normal(0, 2), class = Intercept) +
    prior(normal(0, 1), class = b)
  scode <- stancode(form, data = s2z_dat, prior = bprior)

  for (term in c(
    "group_prec_s2z_1 = inv_square(group_scale_s2z_1);",
    paste0(
      "group_info_s2z_1 += group_prec_s2z_1[j] * ",
      "square(relative_precision_s2z);"
    ),
    paste0(
      "h_group_s2z_1 -= group_prec_s2z_1[j] * ",
      "relative_precision_s2z .* white_level_s2z;"
    ),
    paste0(
      "group_quad_s2z_1 += group_prec_s2z_1[j] * ",
      "dot_self(white_level_s2z);"
    ),
    "L_P_s2z_1 = cholesky_decompose(P_s2z_1);"
  )) {
    expect_true(grepl(term, scode, fixed = TRUE), info = term)
  }
  expect_false(grepl("fast Gaussian Matheron system", scode, fixed = TRUE))
  expect_false(grepl("H_joint_s2z_1", scode, fixed = TRUE))
  expect_match2(
    scode,
    paste0(
      "prior_factor_s2z[, 1:2] = diag_pre_multiply(",
      "sqrt(prior_prec_s2z_1), H_s2z_1 * L_Sigma_s2z_1);"
    )
  )
  expect_match2(
    scode,
    paste0(
      "prior_factor_s2z[, 3:4] = diag_pre_multiply(",
      "sqrt(prior_prec_s2z_1), H_s2z_2 * L_Sigma_s2z_2);"
    )
  )
  expect_match2(scode, "P_s2z_1 = crossprod(prior_factor_s2z);")
  expect_false(grepl(
    "prior_factor_s2z = rep_matrix", scode, fixed = TRUE
  ))
  expect_match2(
    scode,
    paste0(
      "mean_r_s2z_1 = reference_sd_s2z_1 .* ",
      "mean_white_s2z[1:2];"
    )
  )
  expect_match2(
    scode,
    "mean_r_s2z_2 = L_Sigma_s2z_2 * mean_white_s2z[3:4];"
  )
  expect_match2(
    scode, "q_recovered_s2z_1 -= H_s2z_1 * mean_r_s2z_1;"
  )
  expect_match2(
    scode, "q_recovered_s2z_1 -= H_s2z_2 * mean_r_s2z_2;"
  )

  correlated_code <- stancode(
    y ~ x +
      (1 + x | gr(
        g, s2z = TRUE, scale = "varying", dist = "student"
      )) +
      (1 + x | gr(h, s2z = TRUE)),
    data = s2z_dat, prior = bprior
  )
  for (term in c(
    paste0(
      "P_group_s2z_1 += group_prec_s2z_1[j] * ",
      "crossprod(relative_precision_s2z);"
    ),
    paste0(
      "h_group_s2z_1 -= group_prec_s2z_1[j] * ",
      "relative_precision_s2z' * white_level_s2z;"
    ),
    paste0(
      "group_quad_s2z_1 += group_prec_s2z_1[j] * ",
      "dot_self(white_level_s2z);"
    )
  )) {
    expect_true(grepl(term, correlated_code, fixed = TRUE), info = term)
  }
})

test_that("crossed scalar S2Z factors share one omitted-mean system", {
  form <- y ~ 1 +
    (1 | gr(g, s2z = TRUE, center = 0.25)) +
    (1 | gr(h, s2z = TRUE, center = 9 / 34))
  bprior <- prior(normal(0, 2), class = Intercept)
  scode <- stancode(form, data = s2z_dat, prior = bprior)
  sdata <- standata(form, data = s2z_dat, prior = bprior)

  expect_equal(
    unname(sdata$rho_s2z_1), matrix(0.25, nrow = 6L, ncol = 1L)
  )
  expect_equal(
    unname(sdata$rho_s2z_2), matrix(9 / 34, nrow = 8L, ncol = 1L),
    tolerance = 2e-14
  )
  for (term in c(
    "// fast Gaussian Matheron system for S2Z blocks 1, 2",
    "real<lower=0> W_matheron_s2z_1;",
    "real<lower=0> sqrt_W_matheron_s2z_1;",
    "W_matheron_s2z_1 += dot_self(H_s2z_1[1, ] * L_Sigma_s2z_1)",
    "W_matheron_s2z_1 += dot_self(H_s2z_2[1, ] * L_Sigma_s2z_2)",
    "+ log_det_partial_s2z_1",
    "+ log_det_partial_s2z_2",
    "mean_r_s2z_1 += L_Sigma_s2z_1 *",
    "mean_r_s2z_2 += L_Sigma_s2z_2 *",
    "q_recovered_s2z_1 = theta_s2z;",
    "r_1_1 = r_s2z_1_1 + mean_r_s2z_1[1];",
    "r_2_1 = r_s2z_2_1 + mean_r_s2z_2[1];"
  )) {
    expect_true(grepl(term, scode, fixed = TRUE), info = term)
  }
  expect_equal(
    s2z_count_fixed(scode, "cholesky_decompose("), 0L
  )
  expect_equal(
    s2z_count_fixed(
      scode,
      "q_recovered_s2z_1 -= H_s2z_"
    ),
    2L
  )
  expect_match2(scode, "prior_mean_s2z_1[1] = 0;")
  expect_match2(scode, "prior_scale_s2z_1[1] = 2;")
  expect_false(grepl("P_s2z_1", scode, fixed = TRUE))
  expect_false(grepl("H_joint_s2z_1", scode, fixed = TRUE))
  expect_equal(s2z_count_fixed(scode, "real Intercept;"), 1L)
  expect_equal(s2z_count_fixed(scode, "real b_Intercept;"), 1L)

  # Every fitted fixed-partial coordinate map must remain available. Silently
  # recomputing just one block would reinterpret its saved latent variables.
  partial_form <- y ~ 1 +
    (1 | gr(g, s2z = TRUE, center = 0.3)) +
    (1 | gr(h, s2z = TRUE, center = 0.6))
  fit <- brm(partial_form, data = s2z_dat, empty = TRUE)
  expect_setequal(
    names(fit$basis$dpars$mu$re_s2z_center),
    c("rho_s2z_1", "rho_s2z_2")
  )
  corrupt_fit <- fit
  corrupt_fit$basis$dpars$mu$re_s2z_center$rho_s2z_2 <- NULL
  expect_error(
    standata(corrupt_fit),
    "do not contain the fitted coordinate map for group-level ID 2"
  )
})

test_that("independent and correlated interaction blocks stay specialized", {
  form <- y ~ x * z +
    (1 + x * z || gr(g, s2z = TRUE, center = 0.35)) +
    (1 + x * z | gr(h, s2z = TRUE, center = 0.65))
  bprior <- prior(normal(0, 2), class = Intercept) +
    prior(normal(0, 1), class = b)
  scode <- stancode(form, data = s2z_dat, prior = bprior)
  sdata <- standata(form, data = s2z_dat, prior = bprior)

  expect_equal(
    unname(sdata$rho_s2z_1), matrix(0.35, nrow = 6L, ncol = 4L)
  )
  expect_equal(
    unname(sdata$rho_s2z_2), matrix(0.65, nrow = 8L, ncol = 4L)
  )
  for (term in c(
    "// fast Gaussian Matheron system for S2Z blocks 1, 2",
    "matrix[4, 4] W_matheron_s2z_1;",
    "matrix[4, 4] L_W_matheron_s2z_1;",
    "W_matheron_s2z_1 += tcrossprod(H_active_s2z *",
    "L_Sigma_s2z_1 = diag_matrix(sd_1);",
    "L_Sigma_s2z_2 = diag_pre_multiply(sd_2, L_2);",
    "matrix[M_2, M_2] L_partial_s2z",
    "+ log_det_partial_s2z_1",
    "+ log_det_partial_s2z_2",
    "r_1_1 = r_s2z_1_1 + mean_r_s2z_1[1];",
    "r_1_4 = r_s2z_1_4 + mean_r_s2z_1[4];",
    "r_2 = r_s2z_2;",
    "for (j in 1:N_2) r_2[j] += mean_r_s2z_2';",
    "vector[Kc] b;",
    "b = tail(q_recovered_s2z_1, Kc);"
  )) {
    expect_true(grepl(term, scode, fixed = TRUE), info = term)
  }
  # The K=4 independent block remains elementwise. Only the correlated block
  # may perform coefficient-space triangular work before the Matheron update.
  expect_false(grepl(
    "matrix[M_1, M_1] L_partial_s2z", scode, fixed = TRUE
  ))
  expect_false(grepl(
    "mdivide_left_tri_low(L_Sigma_s2z_1", scode, fixed = TRUE
  ))
  expect_match2(
    scode,
    "vector[N_1] white_group_s2z = r_s2z_1[, k] / sd_1[k];"
  )
  expect_match2(
    scode, "group_quad_s2z_1 += dot_self(white_group_s2z);"
  )
  expect_false(grepl("corr_matrix[M_1] Cor_1", scode, fixed = TRUE))
  expect_match2(scode, "corr_matrix[M_2] Cor_2")
  expect_equal(s2z_count_fixed(scode, "cholesky_decompose("), 1L)
  expect_equal(
    s2z_count_fixed(
      scode,
      "q_recovered_s2z_1 -= H_s2z_"
    ),
    2L
  )
  expect_match2(scode, "q_recovered_s2z_1 = theta_s2z;")
  expect_false(grepl("matrix[8, 8] P_s2z_1", scode, fixed = TRUE))
  expect_false(grepl("H_joint_s2z_1", scode, fixed = TRUE))
  expect_equal(s2z_count_fixed(scode, "vector[Kc] b;"), 1L)
  expect_equal(
    s2z_count_fixed(scode, "b = tail(q_recovered_s2z_1, Kc);"), 1L
  )
  for (k in seq_len(4L)) {
    expect_equal(
      s2z_count_fixed(
        scode, sprintf("prior_scale_s2z_1[%s] =", k)
      ),
      1L
    )
  }
})

test_that("Gaussian and Student blocks contribute separately to one solve", {
  form <- y ~ x +
    (1 + x | gr(g, s2z = TRUE, dist = "gaussian")) +
    (1 + x | gr(
      h, s2z = TRUE, center = FALSE, dist = "student"
    ))
  bprior <- prior(normal(0, 2), class = Intercept) +
    prior(normal(0, 1), class = b)
  scode <- stancode(form, data = s2z_dat, prior = bprior)

  for (term in c(
    "group_info_s2z_1 = 1.0 * N_1;",
    "group_scale_s2z_2 = dfm_2;",
    "group_prec_s2z_2 = inv_square(group_scale_s2z_2);",
    "group_info_s2z_2 = sum(group_prec_s2z_2);",
    "h_group_s2z_2 = -white_group_s2z * group_prec_s2z_2;",
    paste0(
      "group_quad_s2z_2 = columns_dot_self(white_group_s2z) * ",
      "group_prec_s2z_2;"
    ),
    "- M_2 * sum(log(group_scale_s2z_2))",
    paste0(
      "P_s2z_1[k, k] += ",
      "group_info_s2z_1;"
    ),
    paste0(
      "P_s2z_1[2 + k, 2 + k] += ",
      "group_info_s2z_2;"
    )
  )) {
    expect_true(grepl(term, scode, fixed = TRUE), info = term)
  }
  expect_false(grepl("group_scale_s2z_1", scode, fixed = TRUE))
  expect_false(grepl("group_prec_s2z_1", scode, fixed = TRUE))
  expect_equal(s2z_count_fixed(scode, "cholesky_decompose("), 1L)
  expect_equal(
    s2z_count_fixed(scode, "normal_lpdf(theta_s2z[1]"), 1L
  )
  expect_equal(
    s2z_count_fixed(scode, "normal_lpdf(theta_s2z[2]"), 1L
  )
})

test_that("shared and varying scales compose in a joint S2Z model", {
  form <- y ~ x +
    (1 + x || gr(
      g, s2z = TRUE, center = 0.4, scale = "varying"
    )) +
    (1 + x | gr(h, s2z = TRUE, center = FALSE))
  bprior <- prior(normal(0, 2), class = Intercept) +
    prior(normal(0, 1), class = b)
  scode <- stancode(form, data = s2z_dat, prior = bprior)
  sdata <- standata(form, data = s2z_dat, prior = bprior)

  expect_equal(dim(sdata$rho_s2z_1), c(6L, 2L))
  expect_null(sdata$rho_s2z_2)
  for (term in c(
    "vector<lower=0>[M_1] sdlog_1;",
    "vector[M_1 * N_1] z_sd_s2z_1;",
    "target += std_normal_lpdf(z_sd_s2z_1);",
    paste0(
      "sum_to_zero_constrain_brms(segment(z_sd_s2z_1, ",
      "(k - 1) * (N_1 - 1) + 1, N_1 - 1));"
    ),
    paste0(
      "reference_sd_s2z_1 = sd_1 .* exp(sdlog_1 .* ",
      "tail(z_sd_s2z_1, M_1) / sqrt(1.0 * N_1));"
    ),
    "matrix<lower=0>[N_1, M_1] sd_level_s2z_1;",
    "L_Sigma_s2z_1 = diag_matrix(reference_sd_s2z_1);",
    paste0(
      "vector[M_1] relative_precision_s2z = ",
      "reference_sd_s2z_1 ./ sd_level_s2z_1[j]';"
    ),
    "group_info_s2z_1 = zeros_vector(M_1);",
    "L_Sigma_s2z_2 = diag_pre_multiply(sd_2, L_2);",
    paste0(
      "P_s2z_1[k, k] += ",
      "group_info_s2z_1[k];"
    ),
    paste0(
      "P_s2z_1[2 + k, 2 + k] += ",
      "group_info_s2z_2;"
    ),
    "matrix<lower=0>[N_1, M_1] sd_level_1;",
    "sd_level_1 = sd_level_s2z_1;"
  )) {
    expect_true(grepl(term, scode, fixed = TRUE), info = term)
  }
  expect_false(grepl(
    "matrix[M_1, M_1] L_level_s2z", scode, fixed = TRUE
  ))
  expect_false(grepl(
    "mdivide_left_tri_low(L_level_s2z", scode, fixed = TRUE
  ))
  expect_false(grepl(
    "std_normal_lpdf(z_sd_s2z_1[k])", scode, fixed = TRUE
  ))
  expect_false(grepl("z_sd_mean_s2z_1", scode, fixed = TRUE))
  expect_equal(
    s2z_count_fixed(
      scode, "target += std_normal_lpdf(z_sd_s2z_1);"
    ),
    1L
  )
  expect_equal(s2z_count_fixed(scode, "cholesky_decompose("), 1L)
})

test_that("Matheron S2Z supports threading without normalizing constants", {
  form <- y ~ 1 +
    (1 | gr(g, s2z = TRUE, center = 0.20)) +
    (1 | gr(h, s2z = TRUE, center = 0.80))
  scode <- stancode(
    form, data = s2z_dat,
    prior = prior(normal(0, 2), class = Intercept),
    threads = threading(2), normalize = FALSE, parse = FALSE
  )

  expect_match2(
    scode,
    paste0(
      "data vector Z_1_1, vector r_s2z_1_1, data array[] int J_2, ",
      "data vector Z_2_1, vector r_s2z_2_1"
    )
  )
  expect_match2(
    scode,
    paste0(
      "target += reduce_sum(partial_log_lik_lpmf, seq, grainsize, Y, ",
      "theta_s2z, sigma, J_1, Z_1_1, r_s2z_1_1, J_2, Z_2_1, ",
      "r_s2z_2_1);"
    )
  )
  expect_match2(scode, "// fast Gaussian Matheron system")
  expect_match2(scode, "+ log_det_partial_s2z_1")
  expect_match2(scode, "+ log_det_partial_s2z_2")
  expect_false(grepl("log(1.0 * N_1)", scode, fixed = TRUE))
  expect_false(grepl("log(1.0 * N_2)", scode, fixed = TRUE))
  expect_false(grepl("log(2 * pi())", scode, fixed = TRUE))
  expect_equal(
    s2z_count_fixed(scode, "cholesky_decompose("), 0L
  )
  expect_equal(
    s2z_count_fixed(
      scode,
      "q_recovered_s2z_1 -= H_s2z_"
    ),
    2L
  )
})

test_that("joint S2Z implementation details follow save_pars", {
  form <- y ~ x * z +
    (1 + x * z || gr(g, s2z = TRUE, center = 0.35)) +
    (1 + x * z | gr(h, s2z = TRUE, center = 0.65))
  default_fit <- brm(form, data = s2z_dat, empty = TRUE)
  saved_fit <- brm(
    form, data = s2z_dat, empty = TRUE,
    save_pars = save_pars(all = TRUE)
  )
  default_excluded <- unlist(
    brms:::exclude_pars(default_fit), use.names = FALSE
  )
  saved_excluded <- unlist(
    brms:::exclude_pars(saved_fit), use.names = FALSE
  )

  internal <- c(
    "H_s2z_1", "L_Sigma_s2z_1", "group_quad_s2z_1",
    "H_s2z_2", "L_Sigma_s2z_2", "group_quad_s2z_2",
    "prior_mean_s2z_1", "prior_scale_s2z_1",
    "W_matheron_s2z_1", "L_W_matheron_s2z_1",
    "theta_white_matheron_s2z_1", "joint_quad_s2z_1",
    "mean_r_s2z_1", "mean_r_s2z_2", "q_recovered_s2z_1"
  )
  for (name in internal) {
    expect_true(name %in% default_excluded, info = name)
    expect_false(name %in% saved_excluded, info = name)
  }
})

test_that("split same-ID S2Z terms validate every constituent", {
  split_dist <- list(
    y ~ x +
      (1 | gr(g, id = "split", s2z = TRUE, dist = "gaussian")) +
      (0 + x | gr(g, id = "split", s2z = TRUE, dist = "student")),
    y ~ x +
      (0 + x | gr(g, id = "split", s2z = TRUE, dist = "student")) +
      (1 | gr(g, id = "split", s2z = TRUE, dist = "gaussian"))
  )
  for (form in split_dist) {
    expect_error(
      stancode(form, data = s2z_dat),
      "same group-level distribution"
    )
  }

  split_cor <- list(
    y ~ x +
      (1 | gr(g, id = "split", s2z = TRUE)) +
      (0 + x || gr(g, id = "split", s2z = TRUE)),
    y ~ x +
      (0 + x || gr(g, id = "split", s2z = TRUE)) +
      (1 | gr(g, id = "split", s2z = TRUE))
  )
  for (form in split_cor) {
    expect_error(stancode(form, data = s2z_dat), "same 'cor' setting")
  }

  by_dat <- transform(
    s2z_dat, f_by = factor(rep(c("a", "b"), each = nrow(s2z_dat) / 2))
  )
  split_by <- list(
    y ~ x +
      (1 | gr(g, id = "split", s2z = TRUE, by = f_by)) +
      (0 + x | gr(g, id = "split", s2z = TRUE)),
    y ~ x +
      (1 | gr(g, id = "split", s2z = TRUE)) +
      (0 + x | gr(g, id = "split", s2z = TRUE, by = f_by))
  )
  for (form in split_by) {
    expect_error(stancode(form, data = by_dat), "Argument 'by'")
  }

  split_pw <- list(
    y ~ x +
      (1 | gr(g, id = "split", s2z = TRUE, pw = w)) +
      (0 + x | gr(g, id = "split", s2z = TRUE)),
    y ~ x +
      (1 | gr(g, id = "split", s2z = TRUE)) +
      (0 + x | gr(g, id = "split", s2z = TRUE, pw = w))
  )
  for (form in split_pw) {
    expect_error(stancode(form, data = s2z_dat), "Argument 'pw'")
  }

  covariance <- diag(nlevels(s2z_dat$g))
  dimnames(covariance) <- list(levels(s2z_dat$g), levels(s2z_dat$g))
  split_cov <- list(
    y ~ x +
      (1 | gr(g, id = "split", s2z = TRUE, cov = covariance)) +
      (0 + x | gr(g, id = "split", s2z = TRUE)),
    y ~ x +
      (1 | gr(g, id = "split", s2z = TRUE)) +
      (0 + x | gr(g, id = "split", s2z = TRUE, cov = covariance))
  )
  for (form in split_cov) {
    expect_error(
      stancode(
        form, data = s2z_dat, data2 = list(covariance = covariance)
      ),
      "same 'cov' setting"
    )
  }

  split_group <- list(
    y ~ x +
      (1 | gr(g, id = "split", s2z = TRUE)) +
      (0 + x | gr(h, id = "split", s2z = TRUE)),
    y ~ x +
      (0 + x | gr(h, id = "split", s2z = TRUE)) +
      (1 | gr(g, id = "split", s2z = TRUE))
  )
  for (form in split_group) {
    expect_error(
      stancode(form, data = s2z_dat),
      "same grouping factor"
    )
  }
})

test_that("unsupported S2Z structures fail clearly", {
  expect_error(
    stancode(y ~ x + (1 + x + z | gr(g, s2z = TRUE)), s2z_dat),
    "matching population-level design column"
  )
  by_dat <- transform(
    s2z_dat, f_by = factor(rep(c("a", "b"), each = nrow(s2z_dat) / 2))
  )
  expect_error(
    stancode(y ~ x + (1 + x | gr(g, by = f_by, s2z = TRUE)), by_dat),
    "not yet supported"
  )
  expect_error(
    stancode(y ~ x + (1 + x | gr(g, pw = w, s2z = TRUE)), s2z_dat),
    "not yet supported"
  )
  expect_error(
    stancode(
      y ~ x + (1 + x | gr(g, s2z = TRUE)), s2z_dat,
      sample_prior = "yes"
    ),
    "sample_prior"
  )
  expect_error(
    stancode(
      y ~ x + (1 + x | gr(g, s2z = TRUE)), s2z_dat,
      prior = prior(laplace(0, 1), class = b, coef = x)
    ),
    "not supported by the sum-to-zero"
  )
  expect_error(
    stancode(
      y ~ x + (1 + x | gr(g, s2z = TRUE)), s2z_dat,
      prior = prior(constant(0), class = sd)
    ),
    "standard deviations fixed with 'constant' must be positive"
  )
  expect_match2(
    stancode(
      y ~ x + (1 + x | gr(g, s2z = TRUE)), s2z_dat,
      prior = prior(constant(1), class = sd)
    ),
    "vector<lower=0>[M_1] sd_1;"
  )
  ordered_mix <- bf(
    y ~ x + (1 + x | gr(g, s2z = TRUE)),
    family = mixture(gaussian, gaussian, order = "mu")
  )
  expect_error(
    stancode(ordered_mix, s2z_dat),
    "Ordered mixture intercepts are not yet supported"
  )
  unordered_mix <- bf(
    y ~ x + (1 + x | gr(g, s2z = TRUE, center = "fisher")),
    family = mixture(gaussian, student, order = "none")
  )
  expect_error(
    stancode(unordered_mix, s2z_dat),
    "has no response-free expected-information rule for family 'mixture'"
  )
  one_group <- transform(s2z_dat, g = factor("only"))
  expect_error(
    stancode(y ~ x + (1 + x | gr(g, s2z = TRUE)), one_group),
    "at least two observed grouping levels"
  )
})

test_that("strict latent S2Z blocks span nonlinear score predictors", {
  expect_false(gr(g, s2z = TRUE)$latent)
  expect_true(gr(g, s2z = TRUE, latent = TRUE)$latent)
  expect_error(
    gr(g, latent = TRUE),
    "latent = TRUE.*requires.*s2z = TRUE"
  )

  dat <- expand.grid(
    person = factor(seq_len(8)),
    item = factor(letters[1:4])
  )
  dat$y <- seq_len(nrow(dat)) / 10
  centered_form <- bf(
    y ~ alpha + loading1 * eta1 + loading2 * eta2,
    alpha ~ 0 + item,
    loading1 ~ 0 + item,
    loading2 ~ 0 + item,
    eta1 ~ 0 +
      (1 | gr(
        person, id = "score", s2z = TRUE, latent = TRUE, center = TRUE
      )),
    eta2 ~ 0 +
      (1 | gr(
        person, id = "score", s2z = TRUE, latent = TRUE, center = TRUE
      )),
    nl = TRUE
  )
  noncentered_form <- bf(
    y ~ alpha + loading1 * eta1 + loading2 * eta2,
    alpha ~ 0 + item,
    loading1 ~ 0 + item,
    loading2 ~ 0 + item,
    eta1 ~ 0 +
      (1 | gr(
        person, id = "score", s2z = TRUE, latent = TRUE, center = FALSE
      )),
    eta2 ~ 0 +
      (1 | gr(
        person, id = "score", s2z = TRUE, latent = TRUE, center = FALSE
      )),
    nl = TRUE
  )

  centered_code <- stancode(
    centered_form, data = dat, normalize = TRUE
  )
  expect_match2(
    centered_code,
    "vector[M_1 * (N_1 - 1)] z_s2z_1;"
  )
  expect_match2(centered_code, "matrix[N_1, M_1] r_s2z_1;")
  expect_match2(centered_code, "vector[N_1] r_s2z_1_eta1_1;")
  expect_match2(centered_code, "vector[N_1] r_s2z_1_eta2_2;")
  expect_match2(
    centered_code,
    "white_latent_s2z = mdivide_left_tri_low("
  )
  expect_match2(
    centered_code,
    "- (N_1 - 1) * sum(log(diagonal(L_Sigma_s2z_1)))"
  )
  expect_match2(
    centered_code,
    "- 0.5 * (N_1 - 1) * M_1 * log(2 * pi())"
  )
  for (term in c("theta_s2z_eta", "H_s2z_1", "mean_r_s2z_1")) {
    expect_false(grepl(term, centered_code, fixed = TRUE), info = term)
  }

  noncentered_code <- stancode(
    noncentered_form, data = dat, normalize = TRUE
  )
  expect_match2(
    noncentered_code,
    "r_s2z_1 = r_s2z_1 * L_Sigma_s2z_1';"
  )
  expect_match2(noncentered_code, "std_normal_lpdf(z_s2z_1)")
  expect_false(grepl(
    "group_quad_s2z_1", noncentered_code, fixed = TRUE
  ))
  expect_false(grepl(
    "sum(log(diagonal(L_Sigma_s2z_1)))",
    noncentered_code, fixed = TRUE
  ))
  expect_equal(
    s2z_count_fixed(
      centered_code, "vector sum_to_zero_constrain_brms(vector y)"
    ),
    1L
  )
  expect_equal(
    s2z_count_fixed(
      noncentered_code, "vector sum_to_zero_constrain_brms(vector y)"
    ),
    1L
  )

  fisher_form <- bf(
    y ~ alpha + loading1 * eta1 + loading2 * eta2,
    alpha ~ 0 + item,
    loading1 ~ 0 + item,
    loading2 ~ 0 + item,
    eta1 ~ 0 +
      (1 | gr(
        person, id = "score", s2z = TRUE, latent = TRUE,
        center = "fisher"
      )),
    eta2 ~ 0 +
      (1 | gr(
        person, id = "score", s2z = TRUE, latent = TRUE,
        center = "fisher"
      )),
    nl = TRUE
  )
  fisher_code <- stancode(fisher_form, data = dat)
  for (term in c(
    "cholesky_factor_corr[M_1] L_1;",
    "fisher_s2z_1_nlp_loading1 = X_loading1 * b_loading1;",
    "fisher_s2z_1_nlp_loading2 = X_loading2 * b_loading2;",
    paste0(
      "design_fisher_s2z[1] = ",
      "(fisher_s2z_1_nlp_loading1[n]) * Z_1_eta1_1[n];"
    ),
    paste0(
      "design_fisher_s2z[2] = ",
      "(fisher_s2z_1_nlp_loading2[n]) * Z_1_eta2_2[n];"
    ),
    paste0(
      "K_fisher_s2z = 1.0 * quad_form(",
      "info_fisher_s2z[j], L_Sigma_s2z_1);"
    )
  )) {
    expect_true(grepl(term, fisher_code, fixed = TRUE), info = term)
  }
})

test_that("strict latent S2Z validation stays narrow and explicit", {
  dat <- expand.grid(
    person = factor(seq_len(5)),
    item = factor(letters[1:3])
  )
  dat$y <- seq_len(nrow(dat)) / 10
  mismatch_form <- bf(
    y ~ loading1 * eta1 + loading2 * eta2,
    loading1 ~ 0 + item,
    loading2 ~ 0 + item,
    eta1 ~ 0 +
      (1 | gr(person, id = "score", s2z = TRUE, latent = TRUE)),
    eta2 ~ 0 + (1 | gr(person, id = "score", s2z = TRUE)),
    nl = TRUE
  )
  expect_error(
    stancode(mismatch_form, data = dat),
    "same 'latent' setting"
  )

  partial_form <- bf(
    y ~ loading1 * eta1 + loading2 * eta2,
    loading1 ~ 0 + item,
    loading2 ~ 0 + item,
    eta1 ~ 0 + (1 | gr(
      person, id = "score", s2z = TRUE, latent = TRUE, center = 0.3
    )),
    eta2 ~ 0 + (1 | gr(
      person, id = "score", s2z = TRUE, latent = TRUE, center = 0.7
    )),
    nl = TRUE
  )
  expect_error(
    stancode(partial_form, data = dat),
    "support only centered, noncentered, and Fisher centering modes"
  )

  fisher_form <- bf(
    y ~ loading * eta,
    loading ~ 0 + item,
    eta ~ 0 + (1 | gr(
      person, id = "score", s2z = TRUE, latent = TRUE, center = "fisher"
    )),
    nl = TRUE
  )
  fisher_code <- stancode(fisher_form, data = dat)
  for (term in c(
    "matrix<lower=0,upper=1>[N_1, M_1] rho_s2z_1;",
    "vector[N] fisher_s2z_1_nlp_loading;",
    "fisher_s2z_1_nlp_loading = X_loading * b_loading;",
    paste0(
      "info_fisher_s2z[J_1[n]] += obs_prec_fisher_s2z * square(",
      "(fisher_s2z_1_nlp_loading[n]) * Z_1_eta_1[n]);"
    ),
    "rho_s2z_1[j, 1] = 1.0 - inv(1.0 + scaled_info_fisher_s2z);",
    "+ log_det_partial_s2z_1"
  )) {
    expect_true(grepl(term, fisher_code, fixed = TRUE), info = term)
  }
  expect_null(
    standata(fisher_form, data = dat)$rho_s2z_1
  )
  fisher_student_code <- stancode(
    fisher_form, data = dat, family = student()
  )
  expect_match2(
    fisher_student_code,
    paste0(
      "(value_fisher_s2z_nu + 1.0) / ",
      "(value_fisher_s2z_nu + 3.0) * ",
      "inv_square(value_fisher_s2z_sigma)"
    )
  )
  bernoulli_dat <- dat
  bernoulli_dat$y <- rep(0:1, length.out = nrow(bernoulli_dat))
  fisher_bernoulli_code <- stancode(
    fisher_form, data = bernoulli_dat, family = bernoulli()
  )
  for (term in c(
    "real value_fisher_s2z_mu = inv_logit(eta_fisher_s2z_mu);",
    paste0(
      "real obs_prec_fisher_s2z = 1.0 * ",
      "((value_fisher_s2z_mu) * (1.0 - (value_fisher_s2z_mu)));"
    ),
    paste0(
      "info_fisher_s2z[J_1[n]] += obs_prec_fisher_s2z * square(",
      "(fisher_s2z_1_nlp_loading[n]) * Z_1_eta_1[n]);"
    )
  )) {
    expect_match2(fisher_bernoulli_code, term)
  }
  expect_false(grepl(
    "Y[n]", strsplit(fisher_bernoulli_code, "model {", fixed = TRUE)[[1]][1],
    fixed = TRUE
  ))
  expect_error(
    suppressWarnings(stancode(
      fisher_form + cor_ar(~ 1 | person), data = dat
    )),
    "does not yet support residual autocorrelation structures"
  )

  student_form <- bf(
    y ~ loading * eta,
    loading ~ 0 + item,
    eta ~ 0 + (1 | gr(
      person, id = "score", s2z = TRUE, latent = TRUE, dist = "student"
    )),
    nl = TRUE
  )
  expect_error(
    stancode(student_form, data = dat),
    "currently require Gaussian group effects"
  )
})

test_that("strict Fisher scores are shared across wide response predictors", {
  dat <- data.frame(
    person = factor(seq_len(7)),
    y1 = sin(seq_len(7) / 3),
    y2 = cos(seq_len(7) / 4)
  )
  response_factor <- function(response) {
    bf(
      as.formula(paste0(response, " ~ loading * eta")),
      lf(loading ~ 1, center = FALSE),
      eta ~ 0 + (1 | gr(
        person, id = "score", s2z = TRUE, latent = TRUE,
        center = "fisher"
      )),
      nl = TRUE
    )
  }
  form <- response_factor("y1") + response_factor("y2") +
    set_rescor(FALSE)
  sdata <- standata(form, data = dat)
  scode <- stancode(form, data = dat)

  expect_equal(sdata$M_1, 1L)
  expect_equal(sdata$NC_1, 0L)
  available <- as.data.frame(get_prior(form, data = dat))
  dimension_sd <- subset(
    available, class == "sd" & group == "person" & nzchar(coef)
  )
  # Every response alias is an available prior scope even though the fitted
  # covariance contains only one dimension. Validation below requires aliases
  # to imply the same prior.
  expect_equal(nrow(dimension_sd), 2L)
  expect_setequal(dimension_sd$resp, c("y1", "y2"))
  expect_false(any(available$class == "cor" & available$group == "person"))
  expect_error(
    stancode(
      form, data = dat,
      prior = prior(lkj(2), class = cor, group = person)
    ),
    "do not correspond to any model parameter"
  )
  for (term in c(
    "vector[N_y1] fisher_s2z_1_nlp_y1_loading;",
    "vector[N_y2] fisher_s2z_1_nlp_y2_loading;",
    "info_fisher_s2z[J_1_y1[n]] += obs_prec_fisher_s2z",
    "info_fisher_s2z[J_1_y2[n]] += obs_prec_fisher_s2z",
    "r_s2z_1_y1_eta_1 = r_s2z_1[, 1];",
    "r_s2z_1_y2_eta_2 = r_s2z_1[, 1];"
  )) {
    expect_true(grepl(term, scode, fixed = TRUE), info = term)
  }

  bframe <- brms:::brmsframe(brmsterms(form), dat)
  covariance_reframe <- brms:::re_s2z_covariance_dimensions(
    bframe$frame$re
  )
  expect_equal(nrow(covariance_reframe), 1L)
  expect_identical(
    brms:::re_s2z_covariance_dimension(
      bframe$frame$re, covariance_r = covariance_reframe
    ),
    c(1L, 1L)
  )
  raw_names <- c(
    "sd_1[1]",
    paste0("r_1_y1_eta_1[", seq_len(7), "]"),
    paste0("r_1_y2_eta_2[", seq_len(7), "]")
  )
  rename_map <- brms:::rename_re(bframe, pars = raw_names)
  sd_map <- Filter(function(x) {
    brms:::is.rlist(x) && any(startsWith(x$fnames, "sd_person__"))
  }, rename_map)
  expect_length(sd_map, 1L)
  expect_equal(sum(sd_map[[1]]$pos), 1L)
  expect_identical(sd_map[[1]]$fnames, "sd_person__y1_eta_Intercept")
  expect_false(any(vapply(rename_map, function(x) {
    brms:::is.rlist(x) && any(startsWith(x$fnames, "cor_person__"))
  }, logical(1))))

  # Strict score defaults are response-scale independent so a wide factor is
  # usable without first overriding response-specific automatic sd priors.
  scaled_dat <- dat
  scaled_dat$y1 <- scaled_dat$y1 / 1000
  scaled_dat$y2 <- scaled_dat$y2 * 1000
  scaled_available <- as.data.frame(get_prior(form, data = scaled_dat))
  scaled_defaults <- subset(
    scaled_available,
    class == "sd" & !nzchar(group) & nlpar == "eta"
  )
  expect_equal(nrow(scaled_defaults), 2L)
  expect_identical(
    unique(scaled_defaults$prior), "student_t(3, 0, 2.5)"
  )
  expect_silent(stancode(form, data = scaled_dat))

  conflicting_sd <- c(
    prior(
      normal(0, 0.7), class = sd, group = person,
      resp = y1, nlpar = eta
    ),
    prior(
      normal(0, 1.4), class = sd, group = person,
      resp = y2, nlpar = eta
    )
  )
  expect_error(
    stancode(form, data = dat, prior = conflicting_sd),
    "shared across responses.*requires identical 'sd' priors"
  )
  expect_error(
    suppressWarnings(stancode(
      response_factor("y1") + response_factor("y2"), data = dat
    )),
    "requires set_rescor\\(FALSE\\)"
  )

  response_two_factor <- function(response) {
    bf(
      as.formula(paste0(
        response, " ~ loading1 * eta1 + loading2 * eta2"
      )),
      lf(loading1 ~ 1, center = FALSE),
      lf(loading2 ~ 1, center = FALSE),
      eta1 ~ 0 + (1 | gr(
        person, id = "score", s2z = TRUE, latent = TRUE,
        center = "fisher"
      )),
      eta2 ~ 0 + (1 | gr(
        person, id = "score", s2z = TRUE, latent = TRUE,
        center = "fisher"
      )),
      nl = TRUE
    )
  }
  two_factor <- response_two_factor("y1") +
    response_two_factor("y2") + set_rescor(FALSE)
  two_data <- standata(two_factor, data = dat)
  two_code <- stancode(two_factor, data = dat)
  expect_equal(two_data$M_1, 2L)
  expect_equal(two_data$NC_1, 1L)
  two_available <- as.data.frame(get_prior(two_factor, data = dat))
  two_dimension_sd <- subset(
    two_available, class == "sd" & group == "person" & nzchar(coef)
  )
  expect_equal(nrow(two_dimension_sd), 4L)
  expect_equal(sum(
    two_available$class == "cor" & two_available$group == "person"
  ), 1L)
  for (term in c(
    "cholesky_factor_corr[M_1] L_1;",
    "design_fisher_s2z[1] = (fisher_s2z_1_nlp_y2_loading1[n]",
    "design_fisher_s2z[2] = (fisher_s2z_1_nlp_y2_loading2[n]",
    "r_s2z_1_y2_eta1_3 = r_s2z_1[, 1];",
    "r_s2z_1_y2_eta2_4 = r_s2z_1[, 2];"
  )) {
    expect_true(grepl(term, two_code, fixed = TRUE), info = term)
  }
  two_frame <- brms:::brmsframe(brmsterms(two_factor), dat)
  expect_identical(
    brms:::re_s2z_covariance_dimension(two_frame$frame$re),
    c(1L, 2L, 1L, 2L)
  )

  # Coefficient-specific prior aliases remain expressible when one nonlinear
  # score predictor contributes more than one covariance dimension.
  dat$x <- seq(-1, 1, length.out = nrow(dat))
  response_slope_factor <- function(response) {
    bf(
      as.formula(paste0(response, " ~ loading * eta")),
      lf(loading ~ 1, center = FALSE),
      eta ~ 0 + (1 + x | gr(
        person, id = "score", s2z = TRUE, latent = TRUE,
        center = "fisher"
      )),
      nl = TRUE
    )
  }
  slope_factor <- response_slope_factor("y1") +
    response_slope_factor("y2") + set_rescor(FALSE)
  shared_slope_prior <- c(
    prior(
      normal(0, 0.5), class = sd, coef = "Intercept", group = person,
      resp = y1, nlpar = eta
    ),
    prior(
      normal(0, 1.5), class = sd, coef = "x", group = person,
      resp = y1, nlpar = eta
    ),
    prior(
      normal(0, 0.5), class = sd, coef = "Intercept", group = person,
      resp = y2, nlpar = eta
    ),
    prior(
      normal(0, 1.5), class = sd, coef = "x", group = person,
      resp = y2, nlpar = eta
    )
  )
  expect_silent(stancode(
    slope_factor, data = dat, prior = shared_slope_prior
  ))

  # Gaussian new-level simulation draws one covariance coordinate and expands
  # it back to both public response aliases. A response-only prediction still
  # finds the representative parameter named from the complete fitted frame.
  ndraws <- 6L
  old_levels <- levels(dat$person)
  used_levels <- c(old_levels, "new")
  covariance_draws <- posterior::as_draws_matrix(matrix(
    1, nrow = ndraws, ncol = 1L,
    dimnames = list(NULL, "sd_person__y1_eta_Intercept")
  ))
  full_reframe <- bframe$frame$re
  full_rdraws <- matrix(
    0, nrow = ndraws, ncol = length(old_levels) * nrow(full_reframe)
  )
  set.seed(90210)
  new_full <- brms:::get_new_rdraws(
    reframe = full_reframe, gf = list(length(old_levels) + 1L),
    rdraws = full_rdraws, used_levels = used_levels,
    old_levels = old_levels, sample_new_levels = "gaussian",
    draws = covariance_draws, covariance_reframe = full_reframe
  )
  expect_equal(dim(new_full), c(ndraws, 2L))
  expect_equal(new_full[, 1], new_full[, 2], tolerance = 0)

  y2_reframe <- subset(full_reframe, resp == "y2")
  y2_rdraws <- matrix(0, nrow = ndraws, ncol = length(old_levels))
  expect_silent(brms:::get_new_rdraws(
    reframe = y2_reframe, gf = list(length(old_levels) + 1L),
    rdraws = y2_rdraws, used_levels = used_levels,
    old_levels = old_levels, sample_new_levels = "gaussian",
    draws = covariance_draws, covariance_reframe = full_reframe
  ))

  empty_fit <- brm(form, data = dat, empty = TRUE)
  excluded <- unlist(
    brms:::exclude_pars(empty_fit, bframe = bframe), use.names = FALSE
  )
  expect_true(all(c(
    "rho_s2z_1", "mean_rho_s2z_1",
    "fisher_s2z_1_nlp_y1_loading",
    "fisher_s2z_1_nlp_y2_loading"
  ) %in% excluded))
})

test_that("strict Fisher scores combine heterogeneous wide families", {
  dat <- data.frame(
    person = factor(seq_len(8)),
    yb = rep(0:1, 4),
    yp = c(0L, 1L, 3L, 2L, 5L, 1L, 4L, 2L)
  )
  response_factor <- function(response, family) {
    bf(
      as.formula(paste0(response, " ~ loading * eta")),
      lf(loading ~ 1, center = FALSE),
      eta ~ 0 + (1 | gr(
        person, id = "score", s2z = TRUE, latent = TRUE,
        center = "fisher"
      )),
      nl = TRUE, family = family
    )
  }
  form <- response_factor("yb", bernoulli()) +
    response_factor("yp", poisson()) + set_rescor(FALSE)
  scode <- stancode(form, data = dat)

  for (term in c(
    "real value_fisher_s2z_mu = inv_logit(eta_fisher_s2z_mu);",
    "real value_fisher_s2z_mu = exp(eta_fisher_s2z_mu);",
    "info_fisher_s2z[J_1_yb[n]] += obs_prec_fisher_s2z",
    "info_fisher_s2z[J_1_yp[n]] += obs_prec_fisher_s2z"
  )) {
    expect_true(grepl(term, scode, fixed = TRUE), info = term)
  }
  transformed_code <- strsplit(scode, "model {", fixed = TRUE)[[1]][1]
  expect_false(grepl("Y_yb[n]", transformed_code, fixed = TRUE))
  expect_false(grepl("Y_yp[n]", transformed_code, fixed = TRUE))
})

test_that("strict Fisher loading references include nonlinear offsets", {
  dat <- expand.grid(
    person = factor(seq_len(4)), item = factor(letters[1:3]),
    KEEP.OUT.ATTRS = FALSE
  )
  dat$x <- seq(0.1, 0.8, length.out = nrow(dat))
  dat$y <- sin(seq_len(nrow(dat)))
  form <- bf(
    y ~ loading * eta,
    loading ~ 0 + item + offset(x),
    eta ~ 0 + (1 | gr(
      person, id = "score", s2z = TRUE, latent = TRUE,
      center = "fisher"
    )),
    nl = TRUE
  )
  scode <- stancode(form, data = dat)

  expect_match2(
    scode,
    paste0(
      "fisher_s2z_1_nlp_loading = ",
      "(X_loading * b_loading) + offsets_loading;"
    )
  )
})
