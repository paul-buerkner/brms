context("Allocation-safe physical S2Z Stan code")

expect_match2 <- brms:::expect_match2

s2z_opt_dat <- local({
  set.seed(6129)
  n <- 72L
  data.frame(
    y = rnorm(n),
    x = seq(-1.5, 2.5, length.out = n),
    z = rep(c(-0.75, 1.25), length.out = n),
    g = factor(rep(seq_len(6L), each = 12L)),
    h = factor(rep(seq_len(8L), each = 9L))
  )
})

s2z_fixed_position <- function(code, pattern) {
  out <- regexpr(pattern, code, fixed = TRUE)[1]
  expect_gt(out, 0L)
  out
}

test_that("fixed S2Z maps live in transformed data", {
  code <- stancode(
    y ~ x + z + (1 + x + z | gr(g, s2z = TRUE)),
    data = s2z_opt_dat
  )

  means_pos <- s2z_fixed_position(
    code, "means_X[i - 1] = mean(X[, i]);"
  )
  map_pos <- s2z_fixed_position(
    code, "H_s2z_1[1, 2] = means_X[1];"
  )
  par_pos <- s2z_fixed_position(code, "\nparameters {\n")
  expect_lt(means_pos, map_pos)
  expect_lt(map_pos, par_pos)
  expect_match2(code, "matrix[3, M_1] H_s2z_1;")

  expect_match2(
    code, "vector[M_1 * (N_1 - 1)] z_s2z_1;"
  )
  expect_match2(
    code,
    paste0(
      "segment(z_s2z_1, (k - 1) * (N_1 - 1) + 1, ",
      "N_1 - 1)"
    )
  )
  expect_false(grepl(
    "array[M_1] vector[N_1 - 1] z_s2z_1;", code, fixed = TRUE
  ))
  expect_match2(code, "vector[M_1] mean_white_s2z =")
  expect_false(grepl(
    "r_s2z_1 + rep_matrix(mhat_s2z_1'", code, fixed = TRUE
  ))
})

test_that("scalar and independent kernels localize fixed and white work", {
  scalar <- stancode(
    y ~ 1 + (1 | gr(g, s2z = TRUE)), data = s2z_opt_dat,
    prior = prior(normal(0, 2), class = Intercept)
  )
  scalar_par <- s2z_fixed_position(scalar, "\nparameters {\n")
  expect_lt(
    s2z_fixed_position(scalar, "vector[1] H_s2z_1;"), scalar_par
  )
  expect_match2(scalar, "real<lower=0> group_quad_s2z_1;")
  expect_match2(scalar, "vector[N_1] white_s2z =")
  expect_false(grepl(
    "vector[N_1] white_s2z_1;", scalar, fixed = TRUE
  ))

  independent <- stancode(
    y ~ x + z + (1 + x + z || gr(g, s2z = TRUE)),
    data = s2z_opt_dat
  )
  independent_par <- s2z_fixed_position(independent, "\nparameters {\n")
  expect_lt(
    s2z_fixed_position(
      independent, "vector[M_1] intercept_map_s2z_1;"
    ),
    independent_par
  )
  expect_match2(
    independent, "intercept_map_s2z_1 = zeros_vector(M_1);"
  )
  expect_match2(
    independent, "vector[M_1 * (N_1 - 1)] z_s2z_1;"
  )
  expect_match2(
    independent,
    paste0(
      "segment(z_s2z_1, (1 - 1) * (N_1 - 1) + 1, ",
      "N_1 - 1)"
    )
  )
  expect_false(grepl(
    "matrix[N_1, M_1] r_s2z_1;", independent, fixed = TRUE
  ))
})

test_that("dense multiblock systems build the prior factor blockwise", {
  code <- stancode(
    y ~ x +
      (1 + x | gr(g, dist = "student", s2z = TRUE)) +
      (1 + x | gr(h, s2z = TRUE)),
    data = s2z_opt_dat,
    prior = prior(normal(0, 1.5), class = Intercept) +
      prior(normal(0, 0.8), class = b)
  )

  expect_match2(code, "matrix[2, 4] prior_factor_s2z;")
  expect_match2(
    code,
    paste0(
      "prior_factor_s2z[, 1:2] = diag_pre_multiply(",
      "sqrt(prior_prec_s2z_1), H_s2z_1 * L_Sigma_s2z_1);"
    )
  )
  expect_match2(code, "P_s2z_1 = crossprod(prior_factor_s2z);")
  expect_match2(code, "P_s2z_1[k, k] += group_info_s2z_1;")
  expect_false(grepl("P_group_s2z_", code, fixed = TRUE))
  expect_false(grepl("H_joint_s2z_1", code, fixed = TRUE))
  expect_match2(code, "q_recovered_s2z_1 = theta_s2z;")
  expect_match2(
    code, "q_recovered_s2z_1 -= H_s2z_1 * mean_r_s2z_1;"
  )
  expect_match2(code, "for (j in 1:N_1) r_1[j] += mean_r_s2z_1';")
})

test_that("Matheron kernels use indexed active rows and row recovery", {
  code <- stancode(
    y ~ x +
      (1 + x | gr(g, s2z = TRUE)) +
      (1 + x | gr(h, s2z = TRUE)),
    data = s2z_opt_dat,
    prior = prior(normal(0, 1.5), class = Intercept) +
      prior(normal(0, 0.8), class = b)
  )

  expect_match2(
    code,
    "W_matheron_s2z_1 = add_diag("
  )
  expect_match2(code, "square(prior_scale_s2z_1[{1, 2}])")
  expect_false(grepl(
    "diag_matrix(square(prior_scale_s2z_1", code, fixed = TRUE
  ))
  expect_match2(code, "H_s2z_1[{1, 2}, ];")
  expect_match2(
    code,
    "theta_s2z[{1, 2}] - prior_mean_s2z_1[{1, 2}]"
  )
  expect_match2(code, "active_score_s2z = zeros_vector(M_1);")
  expect_match2(code, "for (j in 1:N_1) r_1[j] += mean_r_s2z_1';")
  expect_false(grepl(
    "r_1 = r_s2z_1 + rep_matrix(mean_r_s2z_1'", code, fixed = TRUE
  ))
})

test_that("Matheron indexing preserves noncontiguous active rows", {
  code <- stancode(
    y ~ x * z +
      (1 + x * z | gr(g, s2z = TRUE)) +
      (1 + x * z | gr(h, s2z = TRUE)),
    data = s2z_opt_dat,
    prior = prior(normal(0, 1.5), class = Intercept) +
      prior(normal(0, 0.8), class = b, coef = "x:z")
  )

  expect_match2(
    code,
    "W_matheron_s2z_1 = add_diag("
  )
  expect_match2(code, "square(prior_scale_s2z_1[{1, 4}])")
  expect_false(grepl(
    "diag_matrix(square(prior_scale_s2z_1", code, fixed = TRUE
  ))
  expect_match2(code, "H_s2z_1[{1, 4}, ];")
  expect_match2(
    code,
    "theta_s2z[{1, 4}] - prior_mean_s2z_1[{1, 4}]"
  )
  expect_false(grepl("{1, 2}", code, fixed = TRUE))
})

test_that("exact logistic blocks use the same allocation-safe foundation", {
  code <- stancode(
    y ~ x + (1 + x | gr(g, s2z = TRUE)),
    data = s2z_opt_dat,
    prior = prior(normal(0, 2), class = Intercept) +
      prior(logistic(0, 1), class = b)
  )

  expect_match2(
    code, "vector[M_1 * (N_1 - 1)] z_s2z_1;"
  )
  expect_match2(code, "vector[M_1] z_mean_s2z_1;")
  expect_lt(
    s2z_fixed_position(code, "matrix[2, M_1] H_s2z_1;"),
    s2z_fixed_position(code, "\nparameters {\n")
  )
  expect_match2(
    code,
    "group_quad_s2z_1 = dot_self(to_vector(white_group_s2z)) + "
  )
  expect_match2(code, "logistic_lpdf(q_explicit_s2z_1[2] | 0, 1);")
  expect_false(grepl("H_joint_s2z_1", code, fixed = TRUE))

  student_code <- stancode(
    y ~ x + (1 + x | gr(g, dist = "student", s2z = TRUE)),
    data = s2z_opt_dat,
    prior = prior(normal(0, 2), class = Intercept) +
      prior(logistic(0, 1), class = b)
  )
  expect_match2(
    student_code,
    "group_quad_s2z_1 += group_prec_s2z_1[j] * "
  )
  expect_match2(student_code, "vector[M_1] white_level_s2z =")
})

test_that("fixed-only special priors retain ordinary brms scaling", {
  form <- y ~ x + z + (1 | gr(g, s2z = TRUE))
  priors <- list(
    prior(horseshoe(), class = b) +
      prior(normal(0, 2), class = Intercept),
    prior(R2D2(), class = b) +
      prior(normal(0, 2), class = Intercept)
  )

  for (bprior in priors) {
    code <- stancode(form, data = s2z_opt_dat, prior = bprior)
    sdata <- standata(form, data = s2z_opt_dat, prior = bprior)
    expect_equal(sdata$Kscales, sdata$Kc)
    expect_match2(code, "vector[Kc] zfixed_s2z;")
    expect_match2(code, "vector[Kc] fixed_s2z;")
    expect_match2(code, "vector<lower=0>[Kc] sdfixed_s2z;")
    expect_match2(code, "fixed_s2z = zfixed_s2z .* sdfixed_s2z;")
    expect_match2(code, "theta_s2z[2] = fixed_s2z[1];")
    expect_false(grepl("b_s2z_inactive", code, fixed = TRUE))

    scale_pos <- s2z_fixed_position(
      code, "fixed_s2z = zfixed_s2z .* sdfixed_s2z;"
    )
    theta_pos <- s2z_fixed_position(
      code, "theta_s2z[2] = fixed_s2z[1];"
    )
    expect_lt(scale_pos, theta_pos)
  }
})

test_that("blockwise dense algebra equals the concatenated-map algebra", {
  set.seed(902)
  q <- 4L
  H1 <- matrix(rnorm(q * 2L), q)
  H2 <- matrix(rnorm(q * 3L), q)
  L1 <- matrix(c(1.2, 0.1, 0, 0.8), 2L, byrow = TRUE)
  L2 <- matrix(c(1.1, 0, 0, -0.2, 0.9, 0, 0.1, 0.3, 1.3), 3L,
               byrow = TRUE)
  prior_prec <- rexp(q)
  prior_sqrt <- diag(sqrt(prior_prec), q)

  H_joint <- cbind(H1 %*% L1, H2 %*% L2)
  old_factor <- prior_sqrt %*% H_joint
  new_factor <- cbind(
    prior_sqrt %*% H1 %*% L1,
    prior_sqrt %*% H2 %*% L2
  )
  expect_equal(crossprod(new_factor), crossprod(old_factor), tolerance = 1e-12)

  theta <- rnorm(q)
  mean_white <- rnorm(5L)
  old_q <- theta - H_joint %*% mean_white
  new_q <- theta - H1 %*% (L1 %*% mean_white[1:2]) -
    H2 %*% (L2 %*% mean_white[3:5])
  expect_equal(drop(new_q), drop(old_q), tolerance = 1e-12)
})

test_that("explicit Gaussian quadratic shortcut preserves the target", {
  set.seed(1803)
  N <- 9L
  M <- 3L
  raw <- matrix(rnorm(N * M), N, M)
  contrasts <- sweep(raw, 2L, colMeans(raw))
  L <- matrix(c(1.2, 0, 0, 0.1, 0.8, 0, -0.2, 0.3, 1.1), M,
              byrow = TRUE)
  z_mean <- rnorm(M)
  white <- solve(L, t(contrasts))
  mean_white <- z_mean / sqrt(N)

  loop_quad <- sum(vapply(seq_len(N), function(j) {
    sum((white[, j] + mean_white)^2)
  }, numeric(1)))
  shortcut_quad <- sum(white^2) + sum(z_mean^2)
  expect_equal(shortcut_quad, loop_quad, tolerance = 1e-12)
})
