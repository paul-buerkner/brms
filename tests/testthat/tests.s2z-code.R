context("Tests for marginalized sum-to-zero group-level effects")

expect_match2 <- brms:::expect_match2

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

test_that("S2Z Stan code covers intercepts, slopes, and interactions", {
  scode <- stancode(
    y ~ x * z + (1 + x * z | gr(g, s2z = TRUE)), data = s2z_dat
  )
  sdata <- standata(
    y ~ x * z + (1 + x * z | gr(g, s2z = TRUE)), data = s2z_dat
  )

  expect_equal(sdata$M_1, 4)
  expect_equal(sdata$NC_1, 6)
  expect_match2(scode, "array[M_1] vector[N_1 - 1] z_s2z_1;")
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
  expect_match2(scode, "r_1 = r_s2z_1 + rep_matrix(mean_r_s2z_1'")
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
  expect_false(grepl("cholesky_factor_corr[M_1] L_1;", sc_diag, fixed = TRUE))
  expect_match2(sc_diag, "L_Sigma_s2z_1 = diag_matrix(sd_1);")
})

test_that("S2Z handles Student-t effects by conditional marginalization", {
  scode <- stancode(
    y ~ x + (1 + x | gr(g, dist = "student", s2z = TRUE)),
    data = s2z_dat
  )
  expect_match2(scode, "vector<lower=0>[N_1] udf_1;")
  expect_match2(scode, "dfm_1 = sqrt(df_1 * udf_1);")
  expect_match2(scode, "group_scale_s2z_1 = dfm_1;")
  expect_match2(scode, "group_prec_s2z_1 = inv_square(group_scale_s2z_1)")
  expect_match2(scode, "r_s2z_1' * group_prec_s2z_1")
  expect_match2(scode, "M_1 * sum(log(group_scale_s2z_1))")
  expect_match2(scode, "inv_chi_square_lpdf(udf_1 | df_1)")
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
  expect_match2(scode, "normal_lpdf(qhat_s2z_1[1]")
  expect_match2(scode, "normal_lpdf(qhat_s2z_1[2]")

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
    "N_1 * sum(log(diagonal(L_Sigma_s2z_1)))",
    "M_1 * sum(log(group_scale_s2z_1))",
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

test_that("unsupported S2Z structures fail clearly", {
  expect_error(
    stancode(y ~ x + (1 + x + z | gr(g, s2z = TRUE)), s2z_dat),
    "matching population-level design column"
  )
  expect_error(
    stancode(
      y ~ x + (1 + x | gr(g, s2z = TRUE)) +
        (1 + x | gr(h, s2z = TRUE)),
      s2z_dat
    ),
    "Only one marginalized sum-to-zero"
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
    "not supported by the marginalized"
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
  one_group <- transform(s2z_dat, g = factor("only"))
  expect_error(
    stancode(y ~ x + (1 + x | gr(g, s2z = TRUE)), one_group),
    "at least two observed grouping levels"
  )
})
