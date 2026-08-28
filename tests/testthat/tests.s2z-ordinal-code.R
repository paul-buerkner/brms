context("Ordinal threshold code for physical sum-to-zero group effects")

expect_match2 <- brms:::expect_match2

s2z_ordinal_count_fixed <- function(x, pattern) {
  unname(lengths(regmatches(
    x, gregexpr(pattern, x, fixed = TRUE)
  )))
}

s2z_ordinal_error_message <- function(expr) {
  error <- tryCatch(expr, error = identity)
  expect_s3_class(error, "error")
  conditionMessage(error)
}

s2z_ordinal_dat <- local({
  n <- 24L
  data.frame(
    y = ordered(rep(1:4, 6L)),
    x = seq(-1, 2, length.out = n),
    g = factor(rep(seq_len(6L), each = 4L))
  )
})

test_that("flexible ordinal S2Z thresholds use finite likelihood coordinates", {
  bprior <- prior(normal(0.5, 1.5), class = Intercept) +
    prior(normal(0, 1), class = b)
  code <- stancode(
    y ~ x + (1 + x | gr(g, id = "ordinal", s2z = TRUE)),
    data = s2z_ordinal_dat, family = cumulative(), prior = bprior
  )

  expect_match2(
    code,
    "ordered[nthres] theta_Intercept;  // S2Z finite-population thresholds"
  )
  expect_match2(
    code,
    "ordered[nthres] finite_Intercept;  // finite thresholds used by the ordinal likelihood"
  )
  expect_match2(code, "finite_Intercept = theta_Intercept;")
  expect_match2(code, "theta_s2z[1] = theta_Intercept[1];")
  expect_match2(code, "theta_s2z[2] = theta_Intercept[2];")
  expect_match2(code, "theta_s2z[3] = theta_Intercept[3];")
  expect_match2(code, "theta_s2z[4] = theta_s2z_active[1];")
  expect_match2(code, "H_s2z_1[1, 1] = -1;")
  expect_match2(code, "H_s2z_1[4, 2] = 1;")
  expect_match2(
    code,
    "ordered_logistic_lpmf(Y[n] | mu[n], finite_Intercept)"
  )

  # The joint omitted-mean system is the sole owner of conventional
  # threshold priors. The constrained finite parameter is not scored again.
  expect_false(grepl(
    "normal_lpdf(theta_Intercept", code, fixed = TRUE
  ))
  for (k in 1:3) {
    expect_match2(
      code, sprintf("normal_lpdf(theta_s2z[%s] | 0.5, 1.5)", k)
    )
  }

  expect_match2(code, "vector[1] b;")
  expect_match2(code, "b = tail(q_recovered_s2z_1, 1);")
  expect_match2(
    code,
    paste0(
      "Intercept[1] = q_recovered_s2z_1[1];\n",
      "  b_Intercept[1] = Intercept[1] + ",
      "dot_product(means_X, b);"
    )
  )
  recovered_at <- regexpr(
    "q_recovered_s2z_1 -= H_s2z_1 * mean_r_s2z_1;",
    code, fixed = TRUE
  )[[1L]]
  slopes_at <- regexpr(
    "b = tail(q_recovered_s2z_1, 1);", code, fixed = TRUE
  )[[1L]]
  thresholds_at <- regexpr(
    "Intercept[1] = q_recovered_s2z_1[1]", code, fixed = TRUE
  )[[1L]]
  expect_true(0 < recovered_at && recovered_at < slopes_at)
  expect_true(slopes_at < thresholds_at)
})

test_that("equidistant ordinal S2Z translates only the first threshold", {
  bprior <- prior(normal(0.5, 1.5), class = Intercept) +
    prior(normal(1, 0.5), class = delta) +
    prior(normal(0, 1), class = b)
  code <- stancode(
    y ~ x + (1 + x | gr(g, id = "ordinal", s2z = TRUE)),
    data = s2z_ordinal_dat,
    family = cumulative(threshold = "equidistant"), prior = bprior
  )

  expect_match2(
    code,
    "real first_theta_Intercept;  // S2Z finite-population first threshold"
  )
  expect_match2(code, "vector[2] theta_s2z;")
  expect_match2(code, "theta_s2z[1] = first_theta_Intercept;")
  expect_match2(code, "theta_s2z[2] = theta_s2z_active[1];")
  expect_match2(
    code,
    "finite_Intercept[k] = first_theta_Intercept + (k - 1.0) * delta;"
  )
  expect_match2(code, "normal_lpdf(delta | 1, 0.5)")
  expect_equal(
    s2z_ordinal_count_fixed(code, "normal_lpdf(delta | 1, 0.5)"), 1L
  )
  expect_false(grepl("theta_s2z[1] = delta", code, fixed = TRUE))
  expect_match2(
    code,
    paste0(
      "Intercept[k] = q_recovered_s2z_1[1] + (k - 1.0) * delta;\n",
      "    b_Intercept[k] = Intercept[k] + dot_product(means_X, b);"
    )
  )
})

test_that("grouped ordinal S2Z keeps unequal threshold vectors separate", {
  grouped_dat <- data.frame(
    y = c(rep(1:5, 4L), rep(1:6, 4L)),
    x = seq(-1, 2, length.out = 44L),
    threshold_group = rep(c("a", "b"), c(20L, 24L)),
    ncat_group = rep(c(5L, 6L), c(20L, 24L)),
    subject = factor(rep(seq_len(11L), each = 4L))
  )
  form <- y | thres(ncat_group, threshold_group) ~ x +
    (1 + x | gr(subject, id = "ordinal", s2z = TRUE))
  code <- stancode(
    form, data = grouped_dat, family = cumulative(),
    prior = prior(normal(0, 1), class = Intercept) +
      prior(normal(0, 1), class = b), parse = TRUE
  )

  expect_match2(code, "ordered[nthres[1]] theta_Intercept_1;")
  expect_match2(code, "ordered[nthres[2]] theta_Intercept_2;")
  expect_match2(code, "vector[12] theta_s2z;")
  expect_match2(code, "theta_s2z[5] = theta_Intercept_1[5];")
  expect_match2(code, "theta_s2z[6] = theta_Intercept_2[1];")
  expect_match2(code, "theta_s2z[11] = theta_Intercept_2[6];")
  expect_match2(code, "theta_s2z[12] = theta_s2z_active[1];")
  expect_match2(
    code,
    paste0(
      "merged_Intercept[Kthres_start[1]:Kthres_end[1]] = ",
      "finite_Intercept_1;"
    )
  )
  expect_match2(
    code,
    paste0(
      "merged_Intercept[Kthres_start[2]:Kthres_end[2]] = ",
      "finite_Intercept_2;"
    )
  )
  expect_match2(
    code, "Intercept_1[5] = q_recovered_s2z_1[5];"
  )
  expect_match2(
    code, "Intercept_2[1] = q_recovered_s2z_1[6];"
  )
  expect_match2(
    code, "Intercept_2[6] = q_recovered_s2z_1[11];"
  )
})

test_that("ordinal S2Z parses across native family likelihoods", {
  ordinary <- list(cumulative(), cratio(), sratio(), acat())
  for (family in ordinary) {
    code <- stancode(
      y ~ x + (1 + x | gr(g, id = "ordinal", s2z = TRUE)),
      data = s2z_ordinal_dat, family = family,
      prior = prior(normal(0, 1), class = Intercept) +
        prior(normal(0, 1), class = b),
      parse = TRUE
    )
    expect_match2(code, "q_recovered_s2z_1")
    expect_match2(code, "Intercept[1] = q_recovered_s2z_1[1]")
  }

  hurdle_dat <- transform(s2z_ordinal_dat, y = rep(0:4, length.out = 24L))
  hurdle_code <- stancode(
    bf(
      y ~ x + (1 + x | gr(g, id = "ordinal", s2z = TRUE)),
      hu ~ x
    ),
    data = hurdle_dat, family = hurdle_cumulative(),
    prior = prior(normal(0, 1), class = Intercept) +
      prior(normal(0, 1), class = b),
    parse = TRUE
  )
  expect_match2(hurdle_code, "hurdle_cumulative_ordered_logistic_lpmf")
  expect_match2(hurdle_code, "Intercept[1] = q_recovered_s2z_1[1]")
})

test_that("ordinal families parse across supported threshold structures", {
  skip_if_not_installed("rstan")

  ordinary_dat <- data.frame(
    y = ordered(rep(1:5, 8L)),
    x = seq(-1.7, 2.3, length.out = 40L),
    g = factor(rep(seq_len(8L), each = 5L))
  )
  hurdle_dat <- ordinary_dat
  hurdle_dat$y <- rep(0:4, 8L)

  grouped_dat <- data.frame(
    y = c(rep(1:4, 5L), rep(1:5, 5L)),
    x = seq(-1.7, 2.3, length.out = 45L),
    threshold_group = rep(c("a", "b"), c(20L, 25L)),
    nthres_group = rep(c(3L, 4L), c(20L, 25L)),
    g = factor(rep(seq_len(9L), each = 5L))
  )
  grouped_hurdle_dat <- grouped_dat
  grouped_hurdle_dat$y[c(
    1L, 6L, 11L, 16L, 21L, 27L, 33L, 39L, 45L
  )] <- 0L

  families <- list(
    cumulative = cumulative,
    cratio = cratio,
    sratio = sratio,
    acat = acat,
    hurdle_cumulative = hurdle_cumulative
  )
  structures <- c("flexible", "equidistant", "grouped")
  structure_markers <- c(
    flexible = "theta_Intercept",
    equidistant = "first_theta_Intercept",
    grouped = "merged_Intercept"
  )

  for (family_name in names(families)) {
    for (threshold_structure in structures) {
      case <- paste(family_name, threshold_structure, sep = " / ")
      family <- if (threshold_structure == "equidistant") {
        families[[family_name]](threshold = "equidistant")
      } else {
        families[[family_name]]()
      }
      is_hurdle <- family_name == "hurdle_cumulative"
      if (threshold_structure == "grouped") {
        formula <- y | thres(nthres_group, threshold_group) ~ x +
          (1 + x | gr(g, id = "ordinal", s2z = TRUE))
        data <- if (is_hurdle) grouped_hurdle_dat else grouped_dat
      } else {
        formula <- y ~ x +
          (1 + x | gr(g, id = "ordinal", s2z = TRUE))
        data <- if (is_hurdle) hurdle_dat else ordinary_dat
      }

      code <- try(stancode(
        formula, data = data, family = family, parse = TRUE
      ), silent = TRUE)
      expect_false(inherits(code, "try-error"), info = case)
      if (inherits(code, "try-error")) {
        next
      }
      expect_s3_class(code, "brmsmodel")
      expect_match2(code, "q_recovered_s2z_1", info = case)
      expect_match2(
        code, structure_markers[[threshold_structure]], info = case
      )
    }
  }
})

test_that("hurdle ordinal mu and hu S2Z coexist with varying disc", {
  skip_if_not_installed("rstan")

  n <- 40L
  dat <- data.frame(
    y = rep(0:4, 8L),
    x = seq(-1.7, 2.3, length.out = n),
    z = cos(seq_len(n) / 4),
    g = factor(rep(seq_len(8L), each = 5L)),
    h = factor(rep(seq_len(5L), each = 8L))
  )
  form <- bf(
    y ~ x + (1 + x | gr(g, id = "ordinal-mu", s2z = TRUE)),
    hu ~ z + (1 + z | gr(h, id = "hurdle-hu", s2z = TRUE)),
    disc ~ x
  )
  code <- stancode(
    form, data = dat, family = hurdle_cumulative(), parse = TRUE
  )

  # Ordinal mu uses the threshold-aware dense system, while the default
  # logistic hu intercept uses PR1's exact explicit-mean system.
  expect_match2(code, "vector[4] theta_s2z;")
  expect_match2(code, "vector[2] theta_s2z_hu;")
  expect_match2(code, "q_recovered_s2z_1 -= H_s2z_1 * mean_r_s2z_1;")
  expect_match2(code, "q_explicit_s2z_2 = theta_s2z_hu;")
  expect_match2(
    code, "lprior += logistic_lpdf(q_explicit_s2z_2[1] | 0, 1);"
  )
  expect_match2(code, "vector[N] disc = rep_vector(0.0, N);")
  expect_match2(code, "disc += Intercept_disc + Xc_disc * b_disc;")
  expect_match2(
    code,
    paste0(
      "hurdle_cumulative_logit_lpmf(Y[n] | mu[n], hu[n], ",
      "disc[n], finite_Intercept)"
    )
  )
  expect_match2(code, "q_recovered_s2z_2 = q_explicit_s2z_2;")
})

test_that("multiple ordinal S2Z blocks share one weighted dense solve", {
  dat <- transform(
    s2z_ordinal_dat,
    h = factor(rep(seq_len(4L), each = 6L))
  )
  code <- stancode(
    y ~ x +
      (1 + x | gr(g, id = "gaussian", s2z = TRUE)) +
      (1 | gr(
        h, id = "student", dist = "student", s2z = TRUE
      )),
    data = dat, family = cumulative(),
    prior = prior(normal(0, 1), class = Intercept) +
      prior(normal(0, 1), class = b),
    parse = TRUE
  )

  expect_match2(code, "// joint omitted-mean system for S2Z blocks 1, 2")
  expect_match2(code, "group_prec_s2z_2 = inv_square(group_scale_s2z_2);")
  expect_match2(
    code,
    "h_group_s2z_2[k] = -dot_product(white_group_s2z, group_prec_s2z_2);"
  )
  expect_match2(code, "prior_factor_s2z[, 1:2]")
  expect_match2(code, "prior_factor_s2z[, 3:3]")
  expect_equal(
    s2z_ordinal_count_fixed(
      code, "q_recovered_s2z_1 -= H_s2z_"
    ),
    2L
  )
})

test_that("fixed ordinal centering charts retain the local dense system", {
  bprior <- prior(normal(0, 1), class = Intercept) +
    prior(normal(0, 1), class = b)
  forms <- list(
    noncentered = y ~ x + (1 + x | gr(
      g, id = "ordinal-noncentered", s2z = TRUE, center = FALSE
    )),
    partial = y ~ x + (1 + x | gr(
      g, id = "ordinal-partial", s2z = TRUE, center = 0.35
    ))
  )
  code <- lapply(forms, function(form) {
    stancode(
      form, data = s2z_ordinal_dat, family = cumulative(),
      prior = bprior, parse = TRUE
    )
  })

  for (name in names(code)) {
    for (term in c(
      "// S2Z block 1 in a joint omitted-mean system",
      "// joint omitted-mean system for S2Z blocks 1",
      "matrix[2, 2] P_s2z_1;",
      "L_P_s2z_1 = cholesky_decompose(P_s2z_1);",
      "H_s2z_1[1, 1] = -1;",
      "H_s2z_1[4, 2] = 1;",
      "q_recovered_s2z_1 -= H_s2z_1 * mean_r_s2z_1;",
      "ordered_logistic_lpmf(Y[n] | mu[n], finite_Intercept)"
    )) {
      expect_match2(code[[name]], term, info = name)
    }
    expect_false(
      grepl("fast Gaussian Matheron system", code[[name]], fixed = TRUE),
      info = name
    )
  }

  expect_match2(
    code$noncentered, "// standardized orthonormal S2Z coordinates"
  )
  expect_match2(
    code$noncentered, "r_s2z_1 = r_s2z_1 * L_Sigma_s2z_1';"
  )
  expect_false(grepl("rho_s2z_1", code$noncentered, fixed = TRUE))
  expect_false(grepl(
    "log_det_partial_s2z_1", code$noncentered, fixed = TRUE
  ))

  for (term in c(
    "matrix<lower=0,upper=1>[N_1, M_1] rho_s2z_1;",
    "vector[M_1] mean_rho_s2z_1;",
    "real log_det_partial_s2z_1;",
    "+ log_det_partial_s2z_1"
  )) {
    expect_match2(code$partial, term)
  }
  partial_data <- standata(
    forms$partial, data = s2z_ordinal_dat, family = cumulative()
  )
  expect_equal(dim(partial_data$rho_s2z_1), c(6L, 2L))
  expect_equal(unname(partial_data$rho_s2z_1), matrix(0.35, 6L, 2L))
})

test_that("ordinal omitted-mean precision uses its smallest exact shape", {
  bprior <- prior(normal(0, 1), class = Intercept) +
    prior(normal(0, 1), class = b)
  varying_form <- y ~ x + (1 + x || gr(
    g, id = "ordinal-varying", s2z = TRUE,
    center = 0.4, scale = "varying"
  ))
  varying_code <- stancode(
    varying_form, data = s2z_ordinal_dat, family = cumulative(),
    prior = bprior, parse = TRUE
  )
  for (term in c(
    "// S2Z block 1 in a joint omitted-mean system",
    "// joint omitted-mean system for S2Z blocks 1",
    "vector[M_1 * N_1] z_sd_s2z_1;",
    "matrix<lower=0>[N_1, M_1] sd_level_s2z_1;",
    "reference_sd_s2z_1 = sd_1 .* exp(sdlog_1 .*",
    "group_info_s2z_1 = zeros_vector(M_1);",
    paste0(
      "P_s2z_1[k, k] += ",
      "group_info_s2z_1[k];"
    ),
    "sd_level_1 = sd_level_s2z_1;",
    "+ log_det_partial_s2z_1",
    "finite_Intercept"
  )) {
    expect_match2(varying_code, term)
  }
  expect_false(grepl(
    "fast Gaussian Matheron system", varying_code, fixed = TRUE
  ))

  levels_g <- levels(s2z_ordinal_dat$g)
  Omega <- outer(seq_along(levels_g), seq_along(levels_g), function(i, j) {
    0.3^abs(i - j)
  })
  dimnames(Omega) <- list(rev(levels_g), rev(levels_g))
  covariance_form <- y ~ x + (1 + x | gr(
    g, id = "ordinal-covariance", s2z = TRUE,
    center = FALSE, cov = Omega
  ))
  covariance_code <- stancode(
    covariance_form, data = s2z_ordinal_dat, family = cumulative(),
    data2 = list(Omega = Omega), prior = bprior, parse = TRUE
  )
  for (term in c(
    "// S2Z block 1 in a joint omitted-mean system",
    "// joint omitted-mean system for S2Z blocks 1",
    "matrix[N_1, N_1] Lcov_1;",
    "vector[N_1] one_white_cov_s2z;",
    "group_info_s2z_1 = dot_self(one_white_cov_s2z);",
    paste0(
      "P_s2z_1[k, k] += ",
      "group_info_s2z_1;"
    ),
    "h_group_s2z_1 = -white_delta_cov_s2z' * one_white_cov_s2z;",
    "- M_1 * sum(log(diagonal(Lcov_1)))",
    "q_recovered_s2z_1 -= H_s2z_1 * mean_r_s2z_1;"
  )) {
    expect_match2(covariance_code, term)
  }
  expect_false(grepl("P_group_s2z_1", covariance_code, fixed = TRUE))
  covariance_data <- standata(
    covariance_form, data = s2z_ordinal_dat, family = cumulative(),
    data2 = list(Omega = Omega)
  )
  expect_equal(
    unname(covariance_data$Lcov_1),
    unname(t(chol(Omega[levels_g, levels_g])))
  )
})

test_that("ordinal Fisher charts reject with local predictor context", {
  forms <- list(
    fisher = y ~ x + (1 + x | gr(
      g, id = "ordinal-fisher", s2z = TRUE, center = "fisher"
    )),
    auto = y ~ x + (1 + x | gr(
      g, id = "ordinal-auto", s2z = TRUE, center = "auto"
    ))
  )
  for (name in names(forms)) {
    msg <- s2z_ordinal_error_message(stancode(
      forms[[name]], data = s2z_ordinal_dat, family = cumulative()
    ))
    expect_match(
      msg, "S2Z capability 'ordinal_fisher_centering'", fixed = TRUE
    )
    expect_match(msg, "response 'y'", fixed = TRUE)
    expect_match(msg, "family 'cumulative'", fixed = TRUE)
    expect_match(msg, "dpar 'mu'", fixed = TRUE)
    expect_match(msg, "group 'g'", fixed = TRUE)
    expect_match(msg, paste0("ID 'ordinal-", name, "'"), fixed = TRUE)
    expect_match(msg, "coefficient(s) 'Intercept, x'", fixed = TRUE)
    expect_match(msg, "use a fixed center value in [0, 1]", fixed = TRUE)
  }
})

test_that("cross-predictor IDs touching ordinal location reject contextually", {
  hurdle_data <- transform(s2z_ordinal_dat, y = rep(0:4, length.out = 24L))
  hurdle_form <- bf(
    y ~ x + (1 | gr(g, id = "ordinal-shared", s2z = TRUE)),
    hu ~ 1 + (1 | gr(g, id = "ordinal-shared", s2z = TRUE))
  )
  mv_data <- transform(
    s2z_ordinal_dat,
    y_continuous = seq(-1, 1, length.out = nrow(s2z_ordinal_dat))
  )
  mv_form <- bf(
    y ~ x + (1 | gr(g, id = "ordinal-shared", s2z = TRUE)),
    family = cumulative()
  ) + bf(
    y_continuous ~ x +
      (1 | gr(g, id = "ordinal-shared", s2z = TRUE)),
    family = gaussian()
  ) + set_rescor(FALSE)

  messages <- list(
    hurdle = s2z_ordinal_error_message(stancode(
      hurdle_form, data = hurdle_data, family = hurdle_cumulative()
    )),
    multivariate = s2z_ordinal_error_message(stancode(
      mv_form, data = mv_data
    ))
  )
  for (name in names(messages)) {
    msg <- messages[[name]]
    expect_match(
      msg, "S2Z capability 'ordinal_cross_predictor_id'", fixed = TRUE,
      info = name
    )
    expect_match(msg, "response 'y'", fixed = TRUE, info = name)
    expect_match(msg, "family '", fixed = TRUE, info = name)
    expect_match(msg, "dpar 'mu'", fixed = TRUE, info = name)
    expect_match(msg, "group 'g'", fixed = TRUE, info = name)
    expect_match(msg, "ID 'ordinal-shared'", fixed = TRUE, info = name)
    expect_match(
      msg, "use predictor-local S2Z IDs for ordinal location predictors",
      fixed = TRUE, info = name
    )
  }
})

test_that("ordinal-only gates leave nonordinal advanced charts available", {
  dat <- transform(
    s2z_ordinal_dat,
    y_continuous = seq(-1, 1, length.out = nrow(s2z_ordinal_dat))
  )
  levels_g <- levels(dat$g)
  Omega <- outer(seq_along(levels_g), seq_along(levels_g), function(i, j) {
    0.25^abs(i - j)
  })
  dimnames(Omega) <- list(levels_g, levels_g)
  form <- y_continuous ~ x + (1 + x | gr(
    g, id = "nonordinal-advanced", s2z = TRUE,
    center = "fisher", scale = "varying", cov = Omega
  ))
  code <- stancode(
    form, data = dat, data2 = list(Omega = Omega), parse = TRUE
  )

  for (term in c(
    "matrix<lower=0,upper=1>[N_1, M_1] rho_s2z_1;",
    "row_var_fisher_s2z_1",
    "matrix<lower=0>[N_1, M_1] sd_level_s2z_1;",
    "matrix[N_1, N_1] Lcov_1;",
    "+ log_det_partial_s2z_1"
  )) {
    expect_match2(code, term)
  }
  expect_false(grepl(
    "ordinal_fisher_centering", code, fixed = TRUE
  ))
})

test_that("centered slope-only ordinal blocks allow a zero threshold map", {
  dat <- s2z_ordinal_dat
  dat$x <- seq(-1, 1, length.out = nrow(dat))
  code <- stancode(
    y ~ x + (0 + x | gr(g, id = "slope", s2z = TRUE)),
    data = dat, family = cumulative(),
    prior = prior(normal(0, 1), class = Intercept) +
      prior(normal(0, 1), class = b),
    parse = TRUE
  )

  expect_match2(code, "H_s2z_1[4, 1] = 1;")
  expect_false(grepl("H_s2z_1[1, 1] =", code, fixed = TRUE))
  expect_match2(code, "Intercept[1] = q_recovered_s2z_1[1]")
})

test_that("ordinal S2Z solves the general affine design map", {
  scaled <- validate_formula(
    y ~ I(2 * x) + (0 + x | gr(g, id = "scaled", s2z = TRUE)),
    data = s2z_ordinal_dat, family = cumulative()
  )
  scaled_frame <- brmsframe(
    brmsterms(scaled), data = s2z_ordinal_dat
  )$dpars$mu
  scaled_info <- re_s2z_infos(scaled_frame)[[1L]]
  expect_equal(unname(scaled_info$C), matrix(0.5, 1L, 1L))
  expect_equal(unname(scaled_info$a), mean(s2z_ordinal_dat$x))
  expect_true(scaled_info$affine_ok)
  scaled_code <- stancode(
    scaled, data = s2z_ordinal_dat,
    prior = prior(normal(0, 1), class = Intercept) +
      prior(normal(0, 1), class = b),
    parse = TRUE
  )
  expect_match2(scaled_code, "H_s2z_1[4, 1] = 0.5;")

  factor_dat <- transform(
    s2z_ordinal_dat,
    f = factor(rep(c("a", "b", "c"), length.out = nrow(s2z_ordinal_dat)))
  )
  factor_formula <- validate_formula(
    y ~ f + (0 + f | gr(g, id = "factor", s2z = TRUE)),
    data = factor_dat, family = cumulative()
  )
  factor_frame <- brmsframe(
    brmsterms(factor_formula), data = factor_dat
  )$dpars$mu
  factor_info <- re_s2z_infos(factor_frame)[[1L]]
  expect_equal(
    unname(factor_info$C),
    matrix(c(-1, 1, 0, -1, 0, 1), nrow = 2L, byrow = TRUE)
  )
  expect_equal(unname(factor_info$a), rep(1 / 3, 3L))
  expect_true(factor_info$affine_ok)
  expect_silent(stancode(
    factor_formula, data = factor_dat,
    prior = prior(normal(0, 1), class = Intercept) +
      prior(normal(0, 1), class = b),
    parse = TRUE
  ))
})

test_that("multivariate equidistant thresholds use predictor-local deltas", {
  dat <- transform(
    s2z_ordinal_dat,
    y2 = ordered(rep(1:4, 6L)),
    h = factor(rep(seq_len(4L), each = 6L))
  )
  form <- bf(
    y ~ x + (1 + x | gr(g, id = "first", s2z = TRUE)),
    family = cumulative(threshold = "equidistant")
  ) + bf(
    y2 ~ x + (1 + x | gr(h, id = "second", s2z = TRUE)),
    family = cratio(threshold = "equidistant")
  ) + set_rescor(FALSE)
  code <- stancode(
    form, data = dat,
    prior = prior(normal(0, 1), class = Intercept, resp = y) +
      prior(normal(0, 1), class = b, resp = y) +
      prior(normal(0, 1), class = Intercept, resp = y2) +
      prior(normal(0, 1), class = b, resp = y2),
    parse = TRUE
  )

  expect_match2(code, "real<lower=0> delta_y;")
  expect_match2(code, "real delta_y2;")
  expect_match2(
    code,
    "finite_Intercept_y[k] = first_theta_Intercept_y + (k - 1.0) * delta_y;"
  )
  expect_match2(
    code,
    "finite_Intercept_y2[k] = first_theta_Intercept_y2 + (k - 1.0) * delta_y2;"
  )
})

test_that("ordinal finite coordinates stay internal by default", {
  formula <- validate_formula(
    y ~ x + (1 + x | gr(g, id = "ordinal", s2z = TRUE)),
    data = s2z_ordinal_dat, family = cumulative()
  )
  bframe <- brmsframe(
    brmsterms(formula), data = s2z_ordinal_dat
  )$dpars$mu
  excluded <- exclude_pars(bframe, save_pars = save_pars())

  expect_true("theta_Intercept" %in% excluded)
  expect_true("finite_Intercept" %in% excluded)
  expect_true(all(paste0("udf_b_s2z_", 1:4) %in% excluded))
  expect_false("Intercept" %in% excluded)
  expect_false("b_Intercept" %in% excluded)
  expect_false("b" %in% excluded)
})
