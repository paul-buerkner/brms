context("Tests for nonlinear cross-predictor S2Z blocks")

expect_match2 <- brms:::expect_match2

nl_s2z_loss_data <- local({
  group_size <- 6:2
  AY <- factor(rep(1991:1995, group_size))
  dev <- unlist(lapply(group_size, seq_len), use.names = FALSE) * 6
  mean_cum <- 5000 * (1 - exp(-(dev / 45)))
  data.frame(cum = mean_cum, dev, AY)
})

nl_s2z_loss_formula <- bf(
  cum ~ exp(ult) *
    (1 - exp(-(dev / exp(theta))^exp(omega))),
  ult ~ 1 +
    (1 | z | gr(AY, s2z = TRUE, center = FALSE)),
  omega ~ 1 +
    (1 | z | gr(AY, s2z = TRUE, center = FALSE)),
  theta ~ 1 +
    (1 | z | gr(AY, s2z = TRUE, center = FALSE)),
  nl = TRUE
)

nl_s2z_loss_prior <- c(
  prior(normal(8.517193, 0.25), nlpar = "ult"),
  prior(normal(0, 0.5), nlpar = "omega"),
  prior(normal(3.806662, 0.3), nlpar = "theta"),
  prior(exponential(2), class = "sd", group = "AY", nlpar = "ult"),
  prior(exponential(3), class = "sd", group = "AY", nlpar = "omega"),
  prior(exponential(4), class = "sd", group = "AY", nlpar = "theta"),
  prior(lkj(2), class = "cor", group = "AY")
)

test_that("a correlated S2Z ID spans nonlinear predictors", {
  scode <- stancode(
    nl_s2z_loss_formula,
    data = nl_s2z_loss_data,
    family = gaussian(),
    prior = nl_s2z_loss_prior
  )

  # The three nonlinear parameters form one ordinary correlated covariance
  # block. The ID `z` labels that block; it is not a data variable.
  for (term in c(
    "vector[1] z_theta_s2z_ult;",
    "vector[1] z_theta_s2z_omega;",
    "vector[1] z_theta_s2z_theta;",
    "cholesky_factor_corr[M_1] L_1;",
    "vector[M_1 * (N_1 - 1)] z_s2z_1;",
    "matrix[N_1, M_1] r_s2z_1;",
    "vector[1] theta_s2z_ult;",
    "vector[1] theta_s2z_omega;",
    "vector[1] theta_s2z_theta;",
    "matrix[3, M_1] H_s2z_1;",
    "matrix[M_1, M_1] L_mean_s2z_1;",
    "H_s2z_1[1, 1] = 1.0;",
    "H_s2z_1[2, 2] = 1.0;",
    "H_s2z_1[3, 3] = 1.0;",
    "r_s2z_1 = r_s2z_1 * L_Sigma_s2z_1';",
    paste0(
      "L_mean_s2z_1 = cholesky_decompose(diag_matrix(1.0 ./ ",
      "prior_prec_s2z_1) + tcrossprod(L_Sigma_s2z_1) / N_1);"
    ),
    "z_mean_s2z[1] = z_theta_s2z_ult[1];",
    "z_mean_s2z[2] = z_theta_s2z_omega[1];",
    "z_mean_s2z[3] = z_theta_s2z_theta[1];",
    "theta_mean_s2z = prior_mean_s2z_1 + L_mean_s2z_1 * z_mean_s2z;",
    "theta_s2z_ult[1] = theta_mean_s2z[1];",
    "theta_s2z_omega[1] = theta_mean_s2z[2];",
    "theta_s2z_theta[1] = theta_mean_s2z[3];",
    "r_s2z_1_ult_1 = r_s2z_1[, 1];",
    "r_s2z_1_omega_2 = r_s2z_1[, 2];",
    "r_s2z_1_theta_3 = r_s2z_1[, 3];",
    "nlp_ult += X_ult * theta_s2z_ult;",
    "nlp_omega += X_omega * theta_s2z_omega;",
    "nlp_theta += X_theta * theta_s2z_theta;",
    "nlp_ult[n] += r_s2z_1_ult_1[J_1[n]]",
    "nlp_omega[n] += r_s2z_1_omega_2[J_1[n]]",
    "nlp_theta[n] += r_s2z_1_theta_3[J_1[n]]"
  )) {
    expect_match2(scode, term)
  }

  # Population, scale, and correlation priors retain their conventional
  # meanings while the likelihood is evaluated in restricted coordinates.
  for (term in c(
    "normal_lpdf(theta_s2z_ult[1] | 8.517193",
    "normal_lpdf(theta_s2z_omega[1] | 0, 0.5)",
    "normal_lpdf(theta_s2z_theta[1] | 3.806662",
    "exponential_lpdf(sd_1[1] | 2)",
    "exponential_lpdf(sd_1[2] | 3)",
    "exponential_lpdf(sd_1[3] | 4)",
    "lkj_corr_cholesky_lpdf(L_1 | 2)"
  )) {
    expect_match2(scode, term)
  }
  expect_match2(
    scode, "+ sum(log(diagonal(L_mean_s2z_1)))"
  )

  kernel_code <- stancode(
    nl_s2z_loss_formula,
    data = nl_s2z_loss_data,
    family = gaussian(),
    prior = nl_s2z_loss_prior,
    normalize = FALSE
  )
  expect_match2(
    kernel_code, "+ sum(log(diagonal(L_mean_s2z_1)))"
  )
  expect_false(grepl(
    "std_normal_lpdf(z_theta_s2z", scode, fixed = TRUE
  ))

  # Conventional public coefficients and correlated level effects are
  # reconstructed jointly from the omitted group mean.
  for (term in c(
    "vector[M_1] mean_r_s2z_1;",
    "q_recovered_s2z_1 -= H_s2z_1 * mean_r_s2z_1;",
    "b_ult = segment(q_recovered_s2z_1, 1, 1);",
    "b_omega = segment(q_recovered_s2z_1, 2, 1);",
    "b_theta = segment(q_recovered_s2z_1, 3, 1);",
    "r_1_ult_1 = r_1[, 1];",
    "r_1_omega_2 = r_1[, 2];",
    "r_1_theta_3 = r_1[, 3];",
    "corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);"
  )) {
    expect_match2(scode, term)
  }

  expect_false(grepl(
    "- (N_1 - 1) * sum(log(diagonal(L_Sigma_s2z_1)))",
    scode, fixed = TRUE
  ))

  sdata <- standata(
    nl_s2z_loss_formula,
    data = nl_s2z_loss_data,
    family = gaussian(),
    prior = nl_s2z_loss_prior
  )
  expect_identical(sdata$N_1, nlevels(nl_s2z_loss_data$AY))
  expect_identical(sdata$M_1, 3L)
  expect_identical(sdata$NC_1, 3L)

  default_fit <- brm(
    nl_s2z_loss_formula,
    data = nl_s2z_loss_data,
    family = gaussian(),
    prior = nl_s2z_loss_prior,
    empty = TRUE
  )
  saved_fit <- brm(
    nl_s2z_loss_formula,
    data = nl_s2z_loss_data,
    family = gaussian(),
    prior = nl_s2z_loss_prior,
    empty = TRUE,
    save_pars = save_pars(all = TRUE)
  )
  internal <- c(
    "z_theta_s2z_ult", "z_theta_s2z_omega", "z_theta_s2z_theta",
    "L_mean_s2z_1"
  )
  default_excluded <- unlist(
    brms:::exclude_pars(default_fit), use.names = FALSE
  )
  saved_excluded <- unlist(
    brms:::exclude_pars(saved_fit), use.names = FALSE
  )
  expect_true(all(internal %in% default_excluded))
  expect_false(any(internal %in% saved_excluded))
})

test_that("nonlinear cross-predictor S2Z reaches the centered endpoint", {
  form <- bf(
    cum ~ exp(ult) *
      (1 - exp(-(dev / exp(theta))^exp(omega))),
    ult ~ 1 + (1 | z | gr(AY, s2z = TRUE, center = TRUE)),
    omega ~ 1 + (1 | z | gr(AY, s2z = TRUE, center = TRUE)),
    theta ~ 1 + (1 | z | gr(AY, s2z = TRUE, center = TRUE)),
    nl = TRUE
  )
  scode <- stancode(form, data = nl_s2z_loss_data, family = gaussian())

  expect_match2(
    scode,
    "vector[M_1 * (N_1 - 1)] z_s2z_1;  // physical orthonormal"
  )
  expect_match2(
    scode,
    "- (N_1 - 1) * sum(log(diagonal(L_Sigma_s2z_1)))"
  )
  expect_false(grepl(
    "r_s2z_1 = r_s2z_1 * L_Sigma_s2z_1';",
    scode, fixed = TRUE
  ))
  expect_false(grepl("z_theta_s2z_", scode, fixed = TRUE))
  expect_false(grepl("L_mean_s2z_1", scode, fixed = TRUE))
})

test_that("nonlinear cross IDs retain the general endpoint group model", {
  Omega <- diag(nlevels(nl_s2z_loss_data$AY))
  dimnames(Omega) <- rep(list(levels(nl_s2z_loss_data$AY)), 2L)
  form <- bf(
    cum ~ exp(ult) *
      (1 - exp(-(dev / exp(theta))^exp(omega))),
    ult ~ 1 + (1 | gr(
      AY, id = "z", s2z = TRUE, center = FALSE,
      dist = "student", scale = "varying", cov = Omega
    )),
    omega ~ 1 + (1 | gr(
      AY, id = "z", s2z = TRUE, center = FALSE,
      dist = "student", scale = "varying", cov = Omega
    )),
    theta ~ 1 + (1 | gr(
      AY, id = "z", s2z = TRUE, center = FALSE,
      dist = "student", scale = "varying", cov = Omega
    )),
    nl = TRUE
  )
  scode <- stancode(
    form,
    data = nl_s2z_loss_data,
    data2 = list(Omega = Omega),
    family = gaussian(),
    prior = nl_s2z_loss_prior
  )

  for (term in c(
    "group_scale_s2z_1 = dfm_1;",
    "matrix<lower=0>[N_1, M_1] sd_level_s2z_1;",
    "L_Sigma_s2z_1 = diag_pre_multiply(reference_sd_s2z_1, L_1);",
    "- M_1 * sum(log(diagonal(Lcov_1)))"
  )) {
    expect_match2(scode, term)
  }
  expect_false(grepl("z_theta_s2z_", scode, fixed = TRUE))
  expect_false(grepl("L_mean_s2z_1", scode, fixed = TRUE))
})

test_that("nonlinear mean standardization selects only its exact case", {
  flat_code <- stancode(
    nl_s2z_loss_formula,
    data = nl_s2z_loss_data,
    family = gaussian()
  )
  student_prior <- c(
    prior(student_t(4, 8.517193, 0.25), nlpar = "ult"),
    prior(normal(0, 0.5), nlpar = "omega"),
    prior(normal(3.806662, 0.3), nlpar = "theta")
  )
  student_code <- stancode(
    nl_s2z_loss_formula,
    data = nl_s2z_loss_data,
    family = gaussian(),
    prior = student_prior
  )
  slope_form <- bf(
    cum ~ exp(ult) *
      (1 - exp(-(dev / exp(theta))^exp(omega))),
    ult ~ 1 + dev +
      (1 + dev | z | gr(AY, s2z = TRUE, center = FALSE)),
    omega ~ 1 + dev +
      (1 + dev | z | gr(AY, s2z = TRUE, center = FALSE)),
    theta ~ 1 + dev +
      (1 + dev | z | gr(AY, s2z = TRUE, center = FALSE)),
    nl = TRUE
  )
  slope_prior <- c(
    prior(normal(8.517193, 0.25), nlpar = "ult"),
    prior(normal(0, 0.5), nlpar = "omega"),
    prior(normal(3.806662, 0.3), nlpar = "theta")
  )
  slope_code <- stancode(
    slope_form,
    data = nl_s2z_loss_data,
    family = gaussian(),
    prior = slope_prior
  )

  for (scode in list(flat_code, student_code, slope_code)) {
    expect_false(grepl("z_theta_s2z_", scode, fixed = TRUE))
    expect_false(grepl("L_mean_s2z_1", scode, fixed = TRUE))
  }
  expect_match2(flat_code, "vector[1] theta_s2z_ult;")
  expect_match2(student_code, "real<lower=0> udf_b_s2z_ult_1;")
  expect_match2(slope_code, "vector[2] theta_s2z_ult;")
  expect_match2(slope_code, "matrix[6, M_1] H_s2z_1;")
})

test_that("reordered nonlinear predictors keep chart coordinates aligned", {
  form <- bf(
    cum ~ exp(ult) *
      (1 - exp(-(dev / exp(theta))^exp(omega))),
    theta ~ 1 + (1 | z | gr(AY, s2z = TRUE, center = FALSE)),
    ult ~ 1 + (1 | z | gr(AY, s2z = TRUE, center = FALSE)),
    omega ~ 1 + (1 | z | gr(AY, s2z = TRUE, center = FALSE)),
    nl = TRUE
  )
  scode <- stancode(
    form,
    data = nl_s2z_loss_data,
    family = gaussian(),
    prior = nl_s2z_loss_prior
  )

  for (term in c(
    "z_mean_s2z[1] = z_theta_s2z_theta[1];",
    "z_mean_s2z[2] = z_theta_s2z_ult[1];",
    "z_mean_s2z[3] = z_theta_s2z_omega[1];",
    "theta_s2z_theta[1] = theta_mean_s2z[1];",
    "theta_s2z_ult[1] = theta_mean_s2z[2];",
    "theta_s2z_omega[1] = theta_mean_s2z[3];"
  )) {
    expect_match2(scode, term)
  }
})

test_that("the nonlinear finite-mean chart has the exact full target", {
  set.seed(8124)
  M <- 3L
  G <- 10L
  A <- matrix(rnorm(M * M), M)
  Sigma <- tcrossprod(A) + diag(M)
  L_Sigma <- t(chol(Sigma))
  prior_mean <- c(8.517193, 0, 3.806662)
  population_var <- c(0.25^2, 0.5^2, 0.3^2)
  prior_prec <- 1 / population_var

  L_mean <- t(chol(diag(population_var) + Sigma / G))
  z_theta <- rnorm(M)
  theta <- drop(prior_mean + L_mean %*% z_theta)

  z_contrast <- rnorm(M * (G - 1L))
  B <- brms:::re_s2z_basis(G)
  standard_delta <- matrix(NA_real_, G, M)
  for (k in seq_len(M)) {
    at <- (k - 1L) * (G - 1L) + seq_len(G - 1L)
    standard_delta[, k] <- B %*% z_contrast[at]
  }
  delta <- standard_delta %*% t(L_Sigma)
  white_delta <- t(forwardsolve(L_Sigma, t(delta)))

  prior_factor <- diag(sqrt(prior_prec)) %*% L_Sigma
  prior_difference <- sqrt(prior_prec) * (theta - prior_mean)
  P_group <- diag(G, M)
  h_group <- -drop(crossprod(white_delta, rep(1, G)))
  group_quad <- sum(white_delta^2)
  P <- crossprod(prior_factor) + P_group
  h <- drop(crossprod(prior_factor, prior_difference)) + h_group
  L_P <- t(chol(P))
  group_quad <- group_quad - drop(crossprod(h, solve(P, h)))

  log_population <- sum(dnorm(
    theta, mean = prior_mean, sd = sqrt(population_var), log = TRUE
  ))
  log_integrated_group <-
    -0.5 * group_quad -
    sum(log(diag(L_P))) -
    0.5 * G * M * log(2 * pi) +
    0.5 * M * log(2 * pi) +
    0.5 * M * log(G)
  log_jacobian <- sum(log(diag(L_mean)))
  log_generated <- log_population + log_integrated_group + log_jacobian
  log_standard <-
    sum(dnorm(z_theta, log = TRUE)) +
    sum(dnorm(z_contrast, log = TRUE))

  expect_equal(log_generated, log_standard, tolerance = 1e-10)
  expect_equal(h_group, numeric(M), tolerance = 1e-12)
})

test_that("nonlinear cross-predictor S2Z diagnoses non-endpoint charts", {
  partial <- bf(
    cum ~ exp(ult) * (1 - exp(-(dev / exp(theta))^exp(omega))),
    ult ~ 1 + (1 | z | gr(AY, s2z = TRUE, center = 0.25)),
    omega ~ 1 + (1 | z | gr(AY, s2z = TRUE, center = 0.25)),
    theta ~ 1 + (1 | z | gr(AY, s2z = TRUE, center = 0.25)),
    nl = TRUE
  )
  automatic <- bf(
    cum ~ exp(ult) * (1 - exp(-(dev / exp(theta))^exp(omega))),
    ult ~ 1 + (1 | z | gr(AY, s2z = TRUE, center = "auto")),
    omega ~ 1 + (1 | z | gr(AY, s2z = TRUE, center = "auto")),
    theta ~ 1 + (1 | z | gr(AY, s2z = TRUE, center = "auto")),
    nl = TRUE
  )
  diagnostic <- paste0(
    "A cross-predictor sum-to-zero ID spanning nonlinear parameters ",
    "currently supports only center = TRUE or FALSE."
  )

  for (form in list(partial, automatic)) {
    expect_error(
      stancode(form, data = nl_s2z_loss_data, family = gaussian()),
      diagnostic,
      fixed = TRUE
    )
  }
})
