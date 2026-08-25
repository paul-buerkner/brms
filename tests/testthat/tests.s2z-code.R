context("Tests for physical sum-to-zero group-level effects")

expect_match2 <- brms:::expect_match2

s2z_count_fixed <- function(x, pattern) {
  unname(lengths(regmatches(
    x, gregexpr(pattern, x, fixed = TRUE)
  )))
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
          "r_s2z_1_%1$s = sum_to_zero_constrain_brms(segment(z_s2z_1, ",
          "(%1$s - 1) * (N_1 - 1) + 1, N_1 - 1));"
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

test_that("S2Z handles Student-t effects by conditional Gaussian integration", {
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

test_that("an S2Z ID cannot span linear predictors in either term order", {
  mixed_ids <- list(
    bf(
      y ~ 1 + (1 | gr(g, id = "across", s2z = TRUE)),
      sigma ~ 1 + (1 | gr(g, id = "across", s2z = FALSE))
    ),
    bf(
      y ~ 1 + (1 | gr(g, id = "across", s2z = FALSE)),
      sigma ~ 1 + (1 | gr(g, id = "across", s2z = TRUE))
    ),
    bf(
      y ~ 1 + (1 | gr(g, id = "across", s2z = TRUE)),
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
    "(1 | gr(g%s, s2z = TRUE))",
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

test_that("Matheron supports overlapping physical S2Z blocks", {
  form <- y ~ x * z +
    (1 + x | gr(g, s2z = TRUE)) +
    (1 + x * z || gr(h, s2z = TRUE)) +
    (1 | gr(f, s2z = TRUE))
  bprior <- prior(normal(0, 2), class = Intercept) +
    prior(normal(0, 1), class = b)
  scode <- stancode(form, data = s2z_dat, prior = bprior)

  for (term in c(
    "// fast Gaussian Matheron system for S2Z blocks 1, 2, 3",
    "matrix[4, 4] W_matheron_s2z_1;",
    "matrix[4, 4] L_W_matheron_s2z_1;",
    "- (N_1 - 1) * sum(log(diagonal(L_Sigma_s2z_1)))",
    "- (N_2 - 1) * sum(log(diagonal(L_Sigma_s2z_2)))",
    "- (N_3 - 1) * sum(log(diagonal(L_Sigma_s2z_3)))",
    "- 0.5 * (N_1 - 1) * M_1 * log(2 * pi())",
    "- 0.5 * (N_2 - 1) * M_2 * log(2 * pi())",
    "- 0.5 * (N_3 - 1) * M_3 * log(2 * pi())",
    "- 0.5 * 4 * log(2 * pi())",
    "mean_r_s2z_1 += L_Sigma_s2z_1 *",
    "mean_r_s2z_2 += L_Sigma_s2z_2 *",
    "mean_r_s2z_3 += L_Sigma_s2z_3 *",
    "r_s2z_3[, k] / sd_3[k]",
    "q_recovered_s2z_1 = theta_s2z;"
  )) {
    expect_true(grepl(term, scode, fixed = TRUE), info = term)
  }
  expect_false(grepl(
    "L_Sigma_s2z_3, r_s2z_3", scode, fixed = TRUE
  ))
  expect_false(grepl("matrix[7, 7] P_s2z_1", scode, fixed = TRUE))
  expect_false(grepl("L_P_s2z_1", scode, fixed = TRUE))
  expect_false(grepl("H_joint_s2z_1", scode, fixed = TRUE))
  expect_equal(
    s2z_count_fixed(
      scode,
      "W_matheron_s2z_1 += tcrossprod(H_active_s2z * "
    ),
    3L
  )
  expect_equal(s2z_count_fixed(scode, "cholesky_decompose("), 1L)
})

test_that("conditional Student and Cauchy population priors use Matheron", {
  form <- y ~ x +
    (1 + x || gr(g, s2z = TRUE)) +
    (1 + x | gr(h, s2z = TRUE))
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
    (1 | gr(h, s2z = TRUE))
  bprior <- prior(normal(0.2, 1.5), class = Intercept) +
    prior(normal(-0.1, 0.7), class = b)
  scode <- stancode(form, data = s2z_dat, prior = bprior)

  expect_match2(scode, "// fast Gaussian Matheron system")
  expect_match2(scode, "real<lower=0> W_matheron_s2z_1;")
  expect_false(grepl("matrix[2, 2] P_s2z_1", scode, fixed = TRUE))
  expect_match2(
    scode,
    "vector[1] theta_s2z_active;  // S2Z-active finite-population coefficients"
  )
  expect_match2(
    scode,
    "vector[2] fixed_s2z;  // S2Z-inactive regression coefficients"
  )
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

test_that("crossed scalar S2Z factors share one omitted-mean system", {
  form <- y ~ 1 +
    (1 | gr(g, s2z = TRUE)) +
    (1 | gr(h, s2z = TRUE))
  bprior <- prior(normal(0, 2), class = Intercept)
  scode <- stancode(form, data = s2z_dat, prior = bprior)
  for (term in c(
    "// fast Gaussian Matheron system for S2Z blocks 1, 2",
    "real<lower=0> W_matheron_s2z_1;",
    "real<lower=0> sqrt_W_matheron_s2z_1;",
    "W_matheron_s2z_1 += dot_self(H_s2z_1[1, ] * L_Sigma_s2z_1)",
    "W_matheron_s2z_1 += dot_self(H_s2z_2[1, ] * L_Sigma_s2z_2)",
    "- (N_1 - 1) * sum(log(diagonal(L_Sigma_s2z_1)))",
    "- (N_2 - 1) * sum(log(diagonal(L_Sigma_s2z_2)))",
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

})

test_that("independent and correlated interaction blocks stay specialized", {
  form <- y ~ x * z +
    (1 + x * z || gr(g, s2z = TRUE)) +
    (1 + x * z | gr(h, s2z = TRUE))
  bprior <- prior(normal(0, 2), class = Intercept) +
    prior(normal(0, 1), class = b)
  scode <- stancode(form, data = s2z_dat, prior = bprior)
  for (term in c(
    "// fast Gaussian Matheron system for S2Z blocks 1, 2",
    "matrix[4, 4] W_matheron_s2z_1;",
    "matrix[4, 4] L_W_matheron_s2z_1;",
    "W_matheron_s2z_1 += tcrossprod(H_active_s2z *",
    "L_Sigma_s2z_1 = diag_matrix(sd_1);",
    "L_Sigma_s2z_2 = diag_pre_multiply(sd_2, L_2);",
    "- (N_1 - 1) * sum(log(diagonal(L_Sigma_s2z_1)))",
    "- (N_2 - 1) * sum(log(diagonal(L_Sigma_s2z_2)))",
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
    "mdivide_left_tri_low(L_Sigma_s2z_1", scode, fixed = TRUE
  ))
  expect_match2(
    scode,
    "vector[N_1] white_group_s2z = r_s2z_1[, k] / sd_1[k];"
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
    (1 + x | gr(h, s2z = TRUE, dist = "student"))
  bprior <- prior(normal(0, 2), class = Intercept) +
    prior(normal(0, 1), class = b)
  scode <- stancode(form, data = s2z_dat, prior = bprior)

  for (term in c(
    "P_group_s2z_1 = diag_matrix(rep_vector(1.0 * N_1, M_1));",
    "group_scale_s2z_2 = dfm_2;",
    "group_prec_s2z_2 = inv_square(group_scale_s2z_2);",
    "P_group_s2z_2 = diag_matrix(rep_vector(",
    "sum(group_prec_s2z_2), M_2));",
    "h_group_s2z_2 = -white_group_s2z * group_prec_s2z_2;",
    "- M_2 * sum(log(group_scale_s2z_2))",
    "P_s2z_1[1:2, 1:2] += P_group_s2z_1;",
    "P_s2z_1[3:4, 3:4] += P_group_s2z_2;"
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

test_that("Matheron S2Z supports threading without normalizing constants", {
  form <- y ~ 1 +
    (1 | gr(g, s2z = TRUE)) +
    (1 | gr(h, s2z = TRUE))
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
  expect_match2(
    scode, "- (N_1 - 1) * sum(log(diagonal(L_Sigma_s2z_1)))"
  )
  expect_match2(
    scode, "- (N_2 - 1) * sum(log(diagonal(L_Sigma_s2z_2)))"
  )
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
    (1 + x * z || gr(g, s2z = TRUE)) +
    (1 + x * z | gr(h, s2z = TRUE))
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

test_that("fixed-only S2Z internals stay out of public coefficient discovery", {
  form <- y ~ x + z + (1 + x | gr(g, s2z = TRUE))
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
    "theta_s2z", "theta_s2z_active", "fixed_s2z",
    "zfixed_s2z", "sdfixed_s2z", "par_fixed_s2z_1"
  )
  expect_true(all(internal %in% default_excluded))
  expect_false(any(internal %in% saved_excluded))

  draw_names <- c(
    "b_Intercept", "b_x", "b_z", "fixed_s2z[1]",
    "zfixed_s2z[1]", "sdfixed_s2z[1]", "par_fixed_s2z_1"
  )
  expect_identical(
    draw_names[grepl(brms:::fixef_pars(), draw_names)],
    c("b_Intercept", "b_x", "b_z")
  )
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
      "Argument 'cov'"
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
      prior = prior(double_exponential(0, 1), class = b, coef = x)
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
  one_group <- transform(s2z_dat, g = factor("only"))
  expect_error(
    stancode(y ~ x + (1 + x | gr(g, s2z = TRUE)), one_group),
    "at least two observed grouping levels"
  )
})
