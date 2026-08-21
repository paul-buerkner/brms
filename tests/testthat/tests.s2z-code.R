context("Tests for physical sum-to-zero group-level effects")

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

s2z_ten_dat <- local({
  data.frame(
    y = seq(-1, 1, length.out = 80),
    ten = factor(rep(letters[1:10], 8)),
    g = factor(rep(seq_len(8), each = 10))
  )
})

s2z_auto_dat <- local({
  counts <- c(5L, 7L, 10L, 15L, 25L)
  g_index <- rep(seq_along(counts), counts)
  within <- unlist(lapply(counts, seq_len), use.names = FALSE)
  centered_within <- unlist(lapply(counts, function(n) {
    seq_len(n) - mean(seq_len(n))
  }), use.names = FALSE)
  x <- centered_within + 0.17 * g_index
  z <- sin(0.61 * within + 0.23 * g_index) +
    0.08 * centered_within^2
  z[g_index == 1L] <- 2 * x[g_index == 1L]
  w <- cos(0.43 * within - 0.19 * g_index) + 0.6
  w[g_index == 2L] <- 0
  data.frame(
    y = sin(seq_along(g_index) * 0.17), x = x, z = z, w = w,
    g = factor(g_index, levels = seq_along(counts))
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

  default_terms <- brmsterms(default_form)
  centered_terms <- brmsterms(centered_form)
  noncentered_terms <- brmsterms(noncentered_form)
  partial_terms <- brmsterms(partial_form)
  auto_terms <- brmsterms(auto_form)
  expect_equal(default_terms$dpars$mu$re$gcall[[1]]$s2z_center, 1)
  expect_equal(centered_terms$dpars$mu$re$gcall[[1]]$s2z_center, 1)
  expect_equal(noncentered_terms$dpars$mu$re$gcall[[1]]$s2z_center, 0)
  expect_equal(partial_terms$dpars$mu$re$gcall[[1]]$s2z_center, 0.35)
  expect_equal(auto_terms$dpars$mu$re$gcall[[1]]$s2z_center, 0.5)
  expect_false(default_terms$dpars$mu$re$gcall[[1]]$s2z_center_auto)
  expect_true(auto_terms$dpars$mu$re$gcall[[1]]$s2z_center_auto)
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
  expect_error(
    gr(g, center = TRUE),
    "not yet supported for ordinary group-level effects"
  )
  expect_error(gr(g, center = 0.01), "ordinary group-level effects")
  expect_error(gr(g, center = "auto"), "ordinary group-level effects")
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
})

test_that("automatic S2Z weights reflect groupwise design exposure", {
  intercept_form <- y ~ 1 +
    (1 | gr(g, s2z = TRUE, center = "auto"))
  intercept_data <- standata(intercept_form, data = s2z_auto_dat)
  counts <- as.numeric(table(s2z_auto_dat$g))
  expect_equal(dim(intercept_data$rho_s2z_1), c(5L, 1L))
  expect_equal(
    unname(drop(intercept_data$rho_s2z_1)), counts / (counts + 25),
    tolerance = 2e-14
  )
  expect_true(all(diff(drop(intercept_data$rho_s2z_1)) > 0))

  design_form <- y ~ 0 + x + z + x:w +
    (0 + x + z + x:w || gr(g, s2z = TRUE, center = "auto"))
  design_data <- standata(design_form, data = s2z_auto_dat)
  expect_equal(dim(design_data$rho_s2z_1), c(5L, 3L))
  expect_true(all(design_data$rho_s2z_1 >= 0))
  expect_true(all(design_data$rho_s2z_1 <= 1))
  # In group 1, z is exactly twice x. Residualizing either column against
  # the other must therefore report no independent local exposure.
  expect_lt(max(design_data$rho_s2z_1[1L, 1:2]), 1e-20)
  # The interaction is structurally zero throughout group 2.
  expect_equal(design_data$rho_s2z_1[2L, 3L], 0)
  expect_gt(max(design_data$rho_s2z_1[-2L, 3L]), 0)

  rescaled <- s2z_auto_dat
  rescaled$x <- -37 * rescaled$x
  rescaled$z <- 0.021 * rescaled$z
  rescaled$w <- 4.7 * rescaled$w
  rescaled_data <- standata(design_form, data = rescaled)
  expect_equal(
    rescaled_data$rho_s2z_1, design_data$rho_s2z_1,
    tolerance = 2e-13
  )

  permutation <- c(seq(2L, nrow(s2z_auto_dat), by = 2L),
                   seq(1L, nrow(s2z_auto_dat), by = 2L))
  permuted_data <- standata(
    design_form, data = s2z_auto_dat[permutation, , drop = FALSE]
  )
  expect_equal(
    permuted_data$rho_s2z_1, design_data$rho_s2z_1,
    tolerance = 2e-13
  )
  releveled <- s2z_auto_dat
  releveled$g <- factor(
    releveled$g, levels = rev(levels(s2z_auto_dat$g))
  )
  releveled_data <- standata(design_form, data = releveled)
  expect_equal(
    unname(releveled_data$rho_s2z_1),
    unname(design_data$rho_s2z_1[5:1, , drop = FALSE]),
    tolerance = 2e-13
  )

  both_auto <- standata(
    y ~ x +
      (1 | gr(g, id = "mixed-auto", s2z = TRUE, center = "auto")) +
      (0 + x | gr(
        g, id = "mixed-auto", s2z = TRUE, center = "auto"
      )),
    data = s2z_auto_dat
  )
  mixed_auto <- standata(
    y ~ x +
      (1 | gr(g, id = "mixed-auto", s2z = TRUE, center = "auto")) +
      (0 + x | gr(g, id = "mixed-auto", s2z = TRUE, center = 0.8)),
    data = s2z_auto_dat
  )
  expect_equal(
    mixed_auto$rho_s2z_1[, 1L], both_auto$rho_s2z_1[, 1L],
    tolerance = 2e-14
  )
  expect_equal(unname(mixed_auto$rho_s2z_1[, 2L]), rep(0.8, 5L))

  # Once a model is fitted, newdata must not reinterpret its latent
  # coordinates by recomputing the automatic fractions from new predictors.
  fit <- brm(design_form, data = s2z_auto_dat, empty = TRUE)
  shifted <- s2z_auto_dat
  shifted$x <- shifted$x^2 + 3
  shifted$z <- cos(shifted$x)
  shifted$w <- seq_len(nrow(shifted)) / 7
  fitted_data <- standata(fit)
  prediction_data <- standata(fit, newdata = shifted)
  expect_equal(
    prediction_data$rho_s2z_1, fitted_data$rho_s2z_1,
    tolerance = 0
  )

  releveled_prediction <- standata(fit, newdata = releveled)
  expect_equal(
    releveled_prediction$rho_s2z_1, fitted_data$rho_s2z_1,
    tolerance = 0
  )

  new_level_data <- s2z_auto_dat[seq_len(6L), , drop = FALSE]
  new_level_data$g <- factor(
    c(as.character(new_level_data$g[-6L]), "6"),
    levels = as.character(seq_len(6L))
  )
  new_level_prediction <- standata(
    fit, newdata = new_level_data, allow_new_levels = TRUE
  )
  expect_equal(new_level_prediction$N_1, fitted_data$N_1)
  expect_equal(
    new_level_prediction$rho_s2z_1, fitted_data$rho_s2z_1,
    tolerance = 0
  )
  expect_equal(dim(new_level_prediction$rho_s2z_1), c(5L, 3L))
  expect_true(any(new_level_prediction$J_1 > fitted_data$N_1))

  slope_only <- standata(
    fit,
    re_formula = ~ (0 + x | gr(g, s2z = TRUE, center = "auto"))
  )
  expect_equal(slope_only$M_1, 1L)
  expect_equal(
    slope_only$rho_s2z_1[, 1L], fitted_data$rho_s2z_1[, 1L],
    tolerance = 0
  )

  corrupt_fit <- fit
  cached_rho <- corrupt_fit$basis$dpars$mu$re_s2z_center$rho_s2z_1
  corrupt_fit$basis$dpars$mu$re_s2z_center$rho_s2z_1 <-
    cached_rho[-1L, , drop = FALSE]
  expect_error(
    standata(corrupt_fit),
    "do not match the fitted grouping levels and coefficients"
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
    "real log_det_partial_s2z_1;",
    paste0(
      "vector[N_1] scale_partial_s2z = rep_vector(1.0, N_1) - ",
      "rho_s2z_1[, 1] + rho_s2z_1[, 1] * sd_1[1];"
    ),
    paste0(
      "vector[N_1] raw_partial_s2z = sd_1[1] * ",
      "centered_partial_s2z ./ scale_partial_s2z;"
    ),
    paste0(
      "r_s2z_1_1 = raw_partial_s2z - ",
      "rep_vector(mean(raw_partial_s2z), N_1);"
    ),
    paste0(
      "log_det_partial_s2z_1 += -sum(log(scale_partial_s2z)) + ",
      "log(mean(scale_partial_s2z));"
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
      sprintf("r_s2z_1_%s = raw_partial_s2z", k)
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
      "raw_partial_s2z[j] = (L_Sigma_s2z_1 * ",
      "mdivide_left_tri_low(L_partial_s2z, r_s2z_1[j]'))';"
    ),
    "r_s2z_1 = raw_partial_s2z - rep_matrix(mean_partial_s2z, N_1);",
    "log_det_partial_s2z_1 -= sum(log(diagonal(L_partial_s2z)));",
    paste0(
      "log_det_partial_s2z_1 += ",
      "sum(log(diagonal(L_partial_mean_s2z)));"
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
    "relative_sd_s2z_1[, k] = exp(sdlog_1[k]",
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
        "reference_sd_s2z_1[k] = sd_1[k] * exp(sdlog_1[k] * ",
        "z_sd_mean_s2z_1[k] / sqrt(1.0 * N_1));"
      ),
      paste0(
        "relative_sd_s2z_1[, k] = exp(sdlog_1[k] * ",
        "z_sd_centered_s2z);"
      ),
      paste0(
        "sd_level_s2z_1[, k] = reference_sd_s2z_1[k] * ",
        "relative_sd_s2z_1[, k];"
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
        "scale_partial_s2z = rep_vector(1.0, N_1) - ",
        "rho_s2z_1[, 1] + rho_s2z_1[, 1] * ",
        "reference_sd_s2z_1[1];"
      ),
      paste0(
        "raw_partial_s2z = reference_sd_s2z_1[1] * ",
        "centered_partial_s2z ./ scale_partial_s2z;"
      ),
      paste0(
        "group_info_s2z[1] = dot_product(group_prec_s2z_1, ",
        "inv_square(relative_sd_s2z_1[, 1]));"
      ),
      paste0(
        "group_quad_s2z_1 += dot_product(group_prec_s2z_1, ",
        "square((r_s2z_1_1 + mhat_s2z_1[1]) ./ ",
        "sd_level_s2z_1[, 1]));"
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
    "-0.5 * dot_self(white_s2z_1)",
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
            "sd_1[%1$s] * z_s2z_1[%1$s]);"
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
    "sd_1[4] * z_s2z_1[4]",
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
    "r_s2z_1[, k] = sum_to_zero_constrain_brms(z_s2z_1[k]);"
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
    "r_1 = r_s2z_1 + rep_matrix(mean_r_s2z_1', N_1)"
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
    "white_s2z_1 * group_prec_s2z_1",
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
      "reference_sd_s2z_1[1] * z_s2z_1[1]);"
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
    "relative_sd_s2z_1[, k] = exp(sdlog_1[k]",
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
          "reference_sd_s2z_1[%1$s] * z_s2z_1[%1$s]);"
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
    "group_info_s2z[1] = dot_product(group_prec_s2z_1",
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
    "P_s2z_1 += group_prec_s2z_1[j] *",
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

test_that("all S2Z centering modes keep recovery and public names", {
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
  auto_form <- y ~ x * z +
    (1 + x * z | gr(
      g, s2z = TRUE, scale = "varying", center = "auto"
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

  for (form in list(direct_form, partial_form, auto_form)) {
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
    (1 + x | gr(g, s2z = TRUE, center = "auto"))
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
  expect_match2(scode, "vector[N_1] white_s2z_1;")
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
      "white_s2z_1 = mdivide_left_tri_low(",
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
  expect_match2(scode, "array[M_1] vector[N_1 - 1] z_s2z_1;")
  for (k in seq_len(sdata$M_1)) {
    expect_match2(scode, sprintf("vector[N_1] r_s2z_1_%s;", k))
    expect_match2(
      scode,
      sprintf(
        "r_s2z_1_%1$s = sum_to_zero_constrain_brms(z_s2z_1[%1$s]);",
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
  expect_match2(scode, "intercept_map_s2z_1 = rep_vector(0.0, M_1);")
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
      "reference_sd_s2z_1[k] = sd_1[k] * exp(sdlog_1[k] * ",
      "z_sd_mean_s2z_1[k] / sqrt(1.0 * N_1));"
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

test_that("correlated group-varying scales use the heterogeneous kernel", {
  scode <- stancode(
    y ~ x * z +
      (1 + x * z | gr(g, s2z = TRUE, scale = "varying")),
    data = s2z_dat,
    prior = prior(normal(0, 1), class = b) +
      prior(normal(0, 1.5), class = Intercept)
  )

  for (term in c(
    "array[M_1] vector[N_1 - 1] z_sd_s2z_1;",
    "vector[M_1] z_sd_mean_s2z_1;",
    "matrix<lower=0>[N_1, M_1] relative_sd_s2z_1;",
    "matrix<lower=0>[N_1, M_1] sd_level_s2z_1;",
    "vector<lower=0>[M_1] reference_sd_s2z_1;"
  )) {
    expect_match2(scode, term)
  }
  expect_match2(
    scode,
    paste0(
      "relative_sd_s2z_1[, k] = exp(sdlog_1[k] * ",
      "z_sd_centered_s2z);"
    )
  )
  expect_match2(
    scode,
    paste0(
      "reference_sd_s2z_1[k] = sd_1[k] * exp(sdlog_1[k] * ",
      "z_sd_mean_s2z_1[k] / sqrt(1.0 * N_1));"
    )
  )
  expect_match2(
    scode,
    paste0(
      "sd_level_s2z_1[, k] = reference_sd_s2z_1[k] * ",
      "relative_sd_s2z_1[, k];"
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
      "P_s2z_1 += group_prec_s2z_1[j] * ",
      "crossprod(relative_precision_s2z);"
    )
  )
  expect_match2(
    scode,
    paste0(
      "h_s2z -= group_prec_s2z_1[j] * ",
      "relative_precision_s2z' * white_level_s2z;"
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
  expect_match2(scode, "target += std_normal_lpdf(z_sd_mean_s2z_1);")
  expect_match2(scode, "target += std_normal_lpdf(z_sd_s2z_1[k]);")
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
      "group_info_s2z[1] = dot_product(group_prec_s2z_1, ",
      "inv_square(relative_sd_s2z_1[, 1]));"
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
      "group_score_s2z[1] = dot_product(r_s2z_1_1, ",
      "group_prec_s2z_1 .* inv_square(relative_sd_s2z_1[, 1]));"
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
    "z_sd_s2z_1", "z_sd_mean_s2z_1", "relative_sd_s2z_1",
    "reference_sd_s2z_1", "sd_level_s2z_1"
  )) {
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
      "white_s2z_1 = mdivide_left_tri_low(",
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
      "white_s2z_1 * group_prec_s2z_1;"
    )
  )
  expect_match2(
    scode,
    paste0(
      "group_quad_s2z_1 += group_prec_s2z_1[j] * ",
      "dot_self(white_s2z_1[, j]);"
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
    expect_error(stancode(form, data = by_dat), "'by', 'cov', and 'pw'")
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
    expect_error(stancode(form, data = s2z_dat), "'by', 'cov', and 'pw'")
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
      "'by', 'cov', and 'pw'"
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
  expect_error(
    stancode(
      y ~ x + (1 + x | gr(g, s2z = TRUE)) +
        (1 + x | gr(h, s2z = TRUE)),
      s2z_dat
    ),
    "Only one sum-to-zero"
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
  one_group <- transform(s2z_dat, g = factor("only"))
  expect_error(
    stancode(y ~ x + (1 + x | gr(g, s2z = TRUE)), one_group),
    "at least two observed grouping levels"
  )
})
