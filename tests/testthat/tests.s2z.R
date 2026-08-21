context("Tests for physical sum-to-zero group-level effects")

# R mirror of inst/chunks/fun_sum_to_zero.stan. Keeping this small reference
# implementation in the tests makes the geometry checks independent of Stan.
.s2z_constrain_r <- function(y) {
  n <- length(y)
  z <- numeric(n + 1L)
  sum_w <- 0
  for (ii in seq_len(n)) {
    i <- n - ii + 1L
    w <- y[i] / sqrt(i * (i + 1))
    sum_w <- sum_w + w
    z[i] <- z[i] + sum_w
    z[i + 1L] <- z[i + 1L] - i * w
  }
  z
}

.s2z_basis <- function(n) {
  stopifnot(n >= 2L)
  vapply(
    seq_len(n - 1L),
    function(j) .s2z_constrain_r(diag(n - 1L)[, j]),
    numeric(n)
  )
}

.s2z_log_mvn <- function(x, mean, sigma) {
  R <- chol(sigma)
  z <- forwardsolve(t(R), x - mean)
  -0.5 * length(x) * log(2 * pi) - sum(log(diag(R))) -
    0.5 * sum(z^2)
}

# Compare the component-wise diagonal-plus-rank-one implementation with the
# dense complete square it replaces.  This deliberately mirrors the scaled
# equations in .stan_re_s2z_independent without calling code under test.
.s2z_independent_case <- function(k, student = FALSE, center = TRUE) {
  n <- 11L
  basis <- .s2z_basis(n)
  contrast <- matrix(
    sin(seq_len((n - 1L) * k) * 0.71) +
      cos(seq_len((n - 1L) * k) * 0.19),
    nrow = n - 1L, ncol = k
  ) / 3
  r_s2z <- basis %*% contrast
  theta <- seq(-0.8, 1.1, length.out = k)
  prior_mean <- seq(0.45, -0.35, length.out = k)
  prior_sd <- seq(0.65, 1.55, length.out = k)
  prior_prec <- 1 / prior_sd^2
  group_sd <- seq(0.4, 1.25, length.out = k)
  intercept_map <- if (center) {
    c(1, seq(-0.55, 0.65, length.out = k - 1L))
  } else {
    numeric(k)
  }
  H <- diag(k)
  if (center) {
    H[1L, ] <- intercept_map
  }
  group_scale <- if (student) {
    sqrt(4.25 * seq(0.09, 1.37, length.out = n)^1.35)
  } else {
    rep(1, n)
  }
  group_prec <- 1 / group_scale^2
  weighted_contrast <- drop(crossprod(r_s2z, group_prec))

  # Dense K-dimensional complete square.
  Q_group <- diag(1 / group_sd^2, k)
  P <- crossprod(H, prior_prec * H) + sum(group_prec) * Q_group
  h <- drop(crossprod(H, prior_prec * (theta - prior_mean))) -
    drop(Q_group %*% weighted_contrast)
  dense_mode <- drop(solve(P, h))
  dense_qhat <- drop(theta - H %*% dense_mode)

  # Scaled diagonal-plus-rank-one equations used by the specialization.
  base_info <- prior_prec
  base_score <- prior_prec * (theta - prior_mean)
  coupling_prec <- 0
  coupling_resid <- 0
  if (center) {
    base_info[1L] <- 0
    base_score[1L] <- 0
    coupling_prec <- prior_prec[1L]
    coupling_resid <- theta[1L] - prior_mean[1L]
  }
  D <- sum(group_prec) + group_sd^2 * base_info
  scaled_score <- group_sd^2 * base_score
  if (student) {
    scaled_score <- scaled_score - weighted_contrast
  }
  independent_mode <- scaled_score / D
  rank1_info <- coupling_prec * sum(
    group_sd^2 * intercept_map^2 / D
  )
  fast_mode <- independent_mode +
    coupling_prec * group_sd^2 * intercept_map / D *
    (coupling_resid - sum(intercept_map * independent_mode)) /
    (1 + rank1_info)
  fast_qhat <- theta - H %*% fast_mode

  log_joint_measure <- function(group_mean) {
    q <- drop(theta - H %*% group_mean)
    group_log <- sum(vapply(
      seq_len(n),
      function(i) {
        sum(stats::dnorm(
          r_s2z[i, ] + group_mean, 0,
          group_sd * group_scale[i], log = TRUE
        ))
      },
      numeric(1)
    ))
    sum(stats::dnorm(q, prior_mean, prior_sd, log = TRUE)) +
      group_log + 0.5 * k * log(n)
  }
  dense_log <- log_joint_measure(dense_mode) +
    0.5 * k * log(2 * pi) -
    0.5 * as.numeric(determinant(P, logarithm = TRUE)$modulus)
  group_quad <- sum(vapply(
    seq_len(k),
    function(j) {
      sum(((r_s2z[, j] + fast_mode[j]) /
             (group_sd[j] * group_scale))^2)
    },
    numeric(1)
  ))
  fast_log <- sum(stats::dnorm(
    fast_qhat, prior_mean, prior_sd, log = TRUE
  )) - 0.5 * group_quad - (n - 1L) * sum(log(group_sd)) -
    ifelse(student, k * sum(log(group_scale)), 0) -
    0.5 * sum(log(D)) - 0.5 * log1p(rank1_info) -
    0.5 * n * k * log(2 * pi) + 0.5 * k * log(2 * pi) +
    0.5 * k * log(n)

  # The generated-quantities draw is a square-root rank-one downdate of the
  # independent covariance diag(group_sd^2 / D).
  independent_cov <- diag(group_sd^2 / D, k)
  sqrt_rank1 <- sqrt(1 + rank1_info)
  rank1_coef <- coupling_prec /
    (sqrt_rank1 * (1 + sqrt_rank1))
  adjustment <- diag(k) - rank1_coef * outer(
    group_sd^2 * intercept_map / D, intercept_map
  )
  gq_cov <- adjustment %*% independent_cov %*% t(adjustment)

  list(
    r_s2z = r_s2z,
    weighted_contrast = weighted_contrast,
    dense_mode = dense_mode,
    fast_mode = fast_mode,
    dense_qhat = dense_qhat,
    fast_qhat = drop(fast_qhat),
    dense_log = dense_log,
    fast_log = fast_log,
    conditional_cov = solve(P),
    gq_cov = gq_cov,
    rank1_info = rank1_info
  )
}

test_that("the physical sum-to-zero transform is orthonormal", {
  for (n in c(2L, 3L, 7L)) {
    basis <- .s2z_basis(n)
    y <- (seq_len(n - 1L) - n / 3) / 2
    z <- .s2z_constrain_r(y)

    expect_equal(crossprod(basis), diag(n - 1L), tolerance = 1e-14)
    expect_equal(
      drop(crossprod(rep(1, n), basis)), numeric(n - 1L),
      tolerance = 1e-14
    )
    expect_equal(z, drop(basis %*% y), tolerance = 1e-14)
    expect_equal(sum(z), 0, tolerance = 1e-14)
    expect_equal(sum(z^2), sum(y^2), tolerance = 1e-14)
  }
})

test_that("Gaussian completed-square density is normalized in S2Z coordinates", {
  n <- 5L
  m <- 2L
  q <- 3L
  basis <- .s2z_basis(n)

  # H is deliberately rectangular and non-identity. This exercises both a
  # population coefficient shared by several group effects and a coefficient
  # that is only partly shifted by the omitted group mean.
  H <- matrix(c(1, 0.35, -0.2, -0.45, 0.8, 0.3), nrow = q)
  prior_prec <- diag(c(0.8, 1.7, 2.4))
  prior_cov <- solve(prior_prec)
  prior_mean <- c(0.4, -0.7, 0.2)
  Sigma <- matrix(c(1.3, 0.28, 0.28, 0.75), nrow = m)
  Q_Sigma <- solve(Sigma)
  theta_s2z <- c(1.2, -0.3, 0.9)
  contrast <- matrix(
    c(0.2, -0.4, 0.8, -0.6, 0.1, 0.35, -0.2, 0.5),
    nrow = n - 1L
  )
  r_s2z <- basis %*% contrast

  P <- t(H) %*% prior_prec %*% H + n * Q_Sigma
  h <- t(H) %*% prior_prec %*% (theta_s2z - prior_mean) -
    Q_Sigma %*% (t(r_s2z) %*% rep(1, n))
  mhat <- drop(solve(P, h))
  qhat <- drop(theta_s2z - H %*% mhat)

  log_complete <- .s2z_log_mvn(qhat, prior_mean, prior_cov) +
    sum(vapply(
      seq_len(n),
      function(j) {
        .s2z_log_mvn(r_s2z[j, ] + mhat, numeric(m), Sigma)
      },
      numeric(1)
    )) +
    0.5 * m * log(2 * pi) - sum(log(diag(chol(P)))) +
    0.5 * m * log(n)

  # Independently derive the normalized Gaussian density of
  # (theta_s2z, z_s2z)
  # as a linear transformation of conventional coefficients and effects.
  mean_map <- kronecker(
    diag(m), matrix(rep(1 / n, n), nrow = 1L)
  )
  contrast_map <- kronecker(diag(m), t(basis))
  transform <- rbind(
    cbind(diag(q), H %*% mean_map),
    cbind(matrix(0, m * (n - 1L), q), contrast_map)
  )
  conventional_cov <- matrix(0, q + n * m, q + n * m)
  iq <- seq_len(q)
  iu <- q + seq_len(n * m)
  conventional_cov[iq, iq] <- prior_cov
  conventional_cov[iu, iu] <- kronecker(Sigma, diag(n))
  s2z_cov <- transform %*% conventional_cov %*% t(transform)
  s2z_value <- c(theta_s2z, as.vector(t(basis) %*% r_s2z))
  s2z_mean <- c(prior_mean, rep(0, m * (n - 1L)))
  log_direct <- .s2z_log_mvn(s2z_value, s2z_mean, s2z_cov)

  expect_equal(log_complete, log_direct, tolerance = 1e-10)
  expect_equal(
    log_direct - (log_complete - 0.5 * m * log(n)),
    0.5 * m * log(n), tolerance = 1e-12
  )
})

test_that("independent Gaussian effects match the dense complete square", {
  for (k in c(2L, 5L, 10L)) {
    ans <- .s2z_independent_case(k, student = FALSE, center = TRUE)

    expect_equal(colSums(ans$r_s2z), numeric(k), tolerance = 2e-14)
    expect_equal(ans$dense_mode, ans$fast_mode, tolerance = 2e-13)
    expect_equal(ans$dense_qhat, ans$fast_qhat, tolerance = 2e-13)
    expect_equal(ans$dense_log, ans$fast_log, tolerance = 2e-10)
    expect_gt(ans$rank1_info, 0)
  }
})

test_that("independent Student effects match the weighted dense square", {
  for (k in c(2L, 5L, 10L)) {
    ans <- .s2z_independent_case(k, student = TRUE, center = TRUE)

    # Unequal group scales make a physical zero-sum contrast have a nonzero
    # precision-weighted sum; this is the term the Gaussian shortcut omits.
    expect_gt(max(abs(ans$weighted_contrast)), 0.05)
    expect_equal(ans$dense_mode, ans$fast_mode, tolerance = 2e-13)
    expect_equal(ans$dense_qhat, ans$fast_qhat, tolerance = 2e-13)
    expect_equal(ans$dense_log, ans$fast_log, tolerance = 2e-10)
  }
})

test_that("rank-one generated-quantities covariance is exact", {
  for (k in c(2L, 5L, 10L)) {
    for (student in c(FALSE, TRUE)) {
      ans <- .s2z_independent_case(k, student = student, center = TRUE)
      expect_equal(
        ans$gq_cov, ans$conditional_cov,
        tolerance = 3e-14, scale = 1
      )
      expect_gt(min(eigen(ans$gq_cov, symmetric = TRUE)$values), 0)
    }
  }
})

test_that("independent no-intercept effects reduce to diagonal solves", {
  for (k in c(2L, 5L, 10L)) {
    for (student in c(FALSE, TRUE)) {
      ans <- .s2z_independent_case(k, student = student, center = FALSE)
      expect_equal(ans$rank1_info, 0)
      expect_equal(ans$dense_mode, ans$fast_mode, tolerance = 2e-13)
      expect_equal(ans$dense_log, ans$fast_log, tolerance = 2e-10)
      expect_equal(
        ans$gq_cov, ans$conditional_cov,
        tolerance = 3e-14, scale = 1
      )
    }
  }
})

test_that("scalar completed-square density is normalized exactly", {
  n <- 6L
  basis <- .s2z_basis(n)
  prior_mean <- -0.35
  prior_sd <- 1.4
  group_sd <- 0.65
  H <- 0.8
  theta <- 1.15
  contrast <- c(-0.6, 0.2, 0.75, -0.1, 0.45)
  r_s2z <- drop(basis %*% contrast)

  prior_prec <- 1 / prior_sd^2
  group_prec <- rep(1, n)
  Q_sigma <- 1 / group_sd^2
  P <- H^2 * prior_prec + sum(group_prec) * Q_sigma
  h <- H * prior_prec * (theta - prior_mean) -
    Q_sigma * sum(r_s2z * group_prec)
  mhat <- h / P
  A <- H^2 * prior_prec
  a <- H * prior_prec * (theta - prior_mean)
  D <- group_sd^2 * A + n
  mhat_scaled <- group_sd^2 * a / D
  qhat <- theta - H * mhat

  log_complete <- stats::dnorm(
    qhat, prior_mean, prior_sd, log = TRUE
  ) + sum(stats::dnorm(
    r_s2z + mhat, 0, group_sd, log = TRUE
  )) + 0.5 * log(2 * pi) - 0.5 * log(P) + 0.5 * log(n)

  # Independently transform (q, u_1, ..., u_n) to the finite-population
  # coefficient theta = q + H * mean(u) and orthonormal contrasts B' u.
  transform <- rbind(
    c(1, rep(H / n, n)),
    cbind(0, t(basis))
  )
  conventional_cov <- diag(c(prior_sd^2, rep(group_sd^2, n)))
  transformed_cov <- transform %*% conventional_cov %*% t(transform)
  transformed_value <- c(theta, drop(crossprod(basis, r_s2z)))
  transformed_mean <- c(prior_mean, numeric(n - 1L))
  log_direct <- .s2z_log_mvn(
    transformed_value, transformed_mean, transformed_cov
  )
  log_scaled <- stats::dnorm(
    qhat, prior_mean, prior_sd, log = TRUE
  ) - 0.5 * sum(((r_s2z + mhat_scaled) / group_sd)^2) -
    (n - 1) * log(group_sd) - 0.5 * log(D) -
    0.5 * (n - 1) * log(2 * pi) + 0.5 * log(n)

  expect_equal(mhat_scaled, mhat, tolerance = 1e-14)
  expect_equal(log_complete, log_direct, tolerance = 1e-11)
  expect_equal(log_scaled, log_direct, tolerance = 1e-11)
  expect_equal(sum(r_s2z), 0, tolerance = 1e-14)
})

test_that("scalar Student mixture uses the exact heterogeneous complete square", {
  n <- 6L
  basis <- .s2z_basis(n)
  prior_mean <- 0.25
  prior_sd <- 1.1
  group_sd <- 0.7
  H <- -0.45
  theta <- -0.8
  r_s2z <- drop(basis %*% c(0.4, -0.7, 0.15, 0.9, -0.2))

  # Conditional on these scales this is the scalar Student-t group model.
  df <- 4.25
  udf <- c(0.11, 0.35, 0.8, 1.6, 0.27, 0.52)
  group_scale <- sqrt(df * udf)
  group_prec <- 1 / group_scale^2
  prior_prec <- 1 / prior_sd^2
  Q_sigma <- 1 / group_sd^2
  P <- H^2 * prior_prec + sum(group_prec) * Q_sigma
  h <- H * prior_prec * (theta - prior_mean) -
    Q_sigma * sum(r_s2z * group_prec)
  mhat <- h / P
  A <- H^2 * prior_prec
  a <- H * prior_prec * (theta - prior_mean)
  D <- group_sd^2 * A + sum(group_prec)
  mhat_scaled <- (
    group_sd^2 * a - sum(r_s2z * group_prec)
  ) / D

  log_joint_measure <- function(group_mean) {
    q <- theta - H * group_mean
    stats::dnorm(q, prior_mean, prior_sd, log = TRUE) +
      sum(stats::dnorm(
        r_s2z + group_mean, 0, group_sd * group_scale, log = TRUE
      )) + 0.5 * log(n)
  }
  log_marginal <- log_joint_measure(mhat) +
    0.5 * log(2 * pi) - 0.5 * log(P)
  qhat <- theta - H * mhat_scaled
  log_marginal_scaled <- stats::dnorm(
    qhat, prior_mean, prior_sd, log = TRUE
  ) - 0.5 * sum(
    ((r_s2z + mhat_scaled) / (group_sd * group_scale))^2
  ) - (n - 1) * log(group_sd) - sum(log(group_scale)) -
    0.5 * log(D) - 0.5 * (n - 1) * log(2 * pi) + 0.5 * log(n)

  expect_equal(mhat_scaled, mhat, tolerance = 1e-14)
  expect_equal(log_marginal_scaled, log_marginal, tolerance = 1e-11)
  expect_gt(diff(range(group_prec)), 1)
  expect_gt(abs(sum(r_s2z * group_prec)), 0.1)
  for (group_mean in c(mhat, mhat + 0.55, -1.2)) {
    log_conditional <- stats::dnorm(
      group_mean, mhat, sqrt(1 / P), log = TRUE
    )
    expect_equal(
      log_joint_measure(group_mean) - log_conditional,
      log_marginal,
      tolerance = 1e-11
    )
  }
})

test_that("heterogeneous Student mixture scales obey the complete square", {
  n <- 5L
  m <- 2L
  q <- 3L
  basis <- .s2z_basis(n)
  H <- matrix(c(1, 0.2, -0.35, -0.4, 0.75, 0.5), nrow = q)
  prior_mean <- c(-0.2, 0.6, 0.1)
  prior_sd <- c(0.7, 1.1, 1.6)
  prior_prec <- diag(1 / prior_sd^2)
  prior_cov <- diag(prior_sd^2)
  Sigma <- matrix(c(0.9, -0.24, -0.24, 1.4), nrow = m)
  Q_Sigma <- solve(Sigma)
  theta_s2z <- c(0.8, -0.5, 1.1)
  contrast <- matrix(
    c(-0.4, 0.3, 0.9, -0.2, 0.5, -0.7, 0.1, 0.45),
    nrow = n - 1L
  )
  r_s2z <- basis %*% contrast

  # Conditional on the inverse-chi-square mixing variables, Student group
  # effects are Gaussian with level-specific scales sqrt(df * udf).
  df <- 4.5
  udf <- c(0.13, 0.31, 0.8, 1.4, 0.22)
  group_scale <- sqrt(df * udf)
  group_prec <- 1 / group_scale^2
  expect_gt(diff(range(group_prec)), 1)
  expect_gt(sum(abs(drop(t(r_s2z) %*% group_prec))), 0.1)

  P <- t(H) %*% prior_prec %*% H +
    sum(group_prec) * Q_Sigma
  h <- t(H) %*% prior_prec %*% (theta_s2z - prior_mean) -
    Q_Sigma %*% (t(r_s2z) %*% group_prec)
  mhat <- drop(solve(P, h))

  log_joint_measure <- function(group_mean) {
    q_value <- drop(theta_s2z - H %*% group_mean)
    group_log <- sum(vapply(
      seq_len(n),
      function(j) {
        .s2z_log_mvn(
          r_s2z[j, ] + group_mean, numeric(m),
          group_scale[j]^2 * Sigma
        )
      },
      numeric(1)
    ))
    .s2z_log_mvn(q_value, prior_mean, prior_cov) + group_log +
      0.5 * m * log(n)
  }
  log_marginal <- log_joint_measure(mhat) +
    0.5 * m * log(2 * pi) - sum(log(diag(chol(P))))
  conditional_cov <- solve(P)

  # Removing the normalized conditional density of the omitted mean from the
  # joint density must give the same marginal at every possible mean.
  candidate_means <- list(
    mhat,
    mhat + c(0.3, -0.7),
    c(-1.1, 0.4)
  )
  for (group_mean in candidate_means) {
    by_identity <- log_joint_measure(group_mean) -
      .s2z_log_mvn(group_mean, mhat, conditional_cov)
    expect_equal(by_identity, log_marginal, tolerance = 1e-10)
  }
})

test_that("centering preserves numeric and mixed-interaction predictors", {
  dat <- data.frame(
    x = c(-1.7, 0.2, 2.3, 1.1, -0.4, 3.2, 0.7, -2.1, 1.8),
    f = factor(c("a", "b", "c", "a", "c", "b", "a", "b", "c")),
    g = factor(c("g1", "g1", "g2", "g3", "g3", "g3", "g4", "g4", "g4"))
  )
  g_index <- match(dat$g, levels(dat$g))

  for (form in list(~ x, ~ f * x)) {
    X <- model.matrix(form, dat)
    p <- ncol(X)
    means_X <- colMeans(X[, -1L, drop = FALSE])
    Xc <- sweep(X[, -1L, drop = FALSE], 2L, means_X)
    q_value <- seq(-0.8, 0.7, length.out = p)
    r <- matrix(seq_len(nlevels(dat$g) * p), nrow = nlevels(dat$g)) / 7
    r <- sweep(r, 2L, seq_len(p) / 5)
    group_mean <- colMeans(r)
    r_s2z <- sweep(r, 2L, group_mean)

    H <- diag(p)
    H[1L, -1L] <- means_X
    theta_s2z <- drop(q_value + H %*% group_mean)

    eta_conventional <- q_value[1L] + drop(Xc %*% q_value[-1L]) +
      rowSums(X * r[g_index, , drop = FALSE])
    eta_s2z <- theta_s2z[1L] + drop(Xc %*% theta_s2z[-1L]) +
      rowSums(X * r_s2z[g_index, , drop = FALSE])
    raw_q <- c(
      q_value[1L] - sum(means_X * q_value[-1L]), q_value[-1L]
    )

    expect_equal(colSums(r_s2z), numeric(p), tolerance = 1e-14)
    expect_equal(drop(theta_s2z - H %*% group_mean), q_value,
                 tolerance = 1e-14)
    expect_equal(eta_s2z, eta_conventional, tolerance = 1e-13)
    expect_equal(drop(X %*% raw_q),
                 q_value[1L] + drop(Xc %*% q_value[-1L]),
                 tolerance = 1e-14)
  }

  mixed_X <- model.matrix(~ f * x, dat)
  expect_true(any(grepl(":", colnames(mixed_X), fixed = TRUE)))
  expect_gt(max(abs(colMeans(mixed_X[, -1L, drop = FALSE]))), 0.1)
})

test_that("generated Stan code contains the S2Z algebra and measure", {
  expect_match2 <- brms:::expect_match2
  dat <- data.frame(
    y = seq(-1, 1, length.out = 12),
    x = c(-2, -1, 0.5, 1, 2, 3, -0.5, 1.5, 2.5, 0.25, 1.25, 3.25),
    f = factor(rep(c("a", "b"), 6)),
    g = factor(rep(letters[1:4], each = 3))
  )
  bprior <- prior(normal(0, 1), class = "Intercept") +
    prior(normal(0, 1), class = "b")

  gaussian_code <- stancode(
    y ~ f * x + (1 + f * x | gr(g, s2z = TRUE)),
    data = dat, prior = bprior, normalize = TRUE
  )
  expect_match2(
    gaussian_code,
    "array[M_1] vector[N_1 - 1] z_s2z_1;"
  )
  expect_match2(
    gaussian_code,
    "r_s2z_1[, k] = sum_to_zero_constrain_brms(z_s2z_1[k]);"
  )
  expect_match2(gaussian_code, "H_s2z_1[1, 2] = means_X[1];")
  expect_match2(gaussian_code, "H_s2z_1[1, 4] = means_X[3];")
  expect_match2(
    gaussian_code,
    "P_s2z_1 = crossprod(diag_pre_multiply(sqrt(prior_prec_s2z_1), H_s2z_1)) + sum(group_prec_s2z_1) * Q_Sigma_s2z_1;"
  )
  expect_match2(
    gaussian_code,
    "+ 0.5 * M_1 * log(1.0 * N_1)"
  )
  expect_match2(
    gaussian_code,
    "normal_id_glm_lpdf(Y | Xc, mu, tail(theta_s2z, 3), sigma)"
  )

  student_code <- stancode(
    y ~ f * x +
      (1 + f * x | gr(g, s2z = TRUE, dist = "student")),
    data = dat, prior = bprior, normalize = TRUE
  )
  expect_match2(student_code, "group_scale_s2z_1 = dfm_1;")
  expect_match2(
    student_code,
    "group_prec_s2z_1 = inv_square(group_scale_s2z_1);"
  )
  expect_match2(
    student_code,
    "r_s2z_1' * group_prec_s2z_1"
  )
  expect_match2(
    student_code,
    "- M_1 * sum(log(group_scale_s2z_1))"
  )
})

test_that("save_pars(all = TRUE) keeps S2Z internals extractor-safe", {
  dat <- data.frame(
    y = seq_len(8),
    x = seq(-1, 1, length.out = 8),
    g = factor(rep(seq_len(4), each = 2))
  )
  fit <- brm(
    y ~ x + (1 + x | gr(g, s2z = TRUE)), data = dat,
    empty = TRUE, save_pars = save_pars(all = TRUE)
  )
  excluded <- unlist(brms:::exclude_pars(fit), use.names = FALSE)

  # save_pars(all = TRUE) retains the internal parameter for operations such
  # as bridge sampling, but its name must not look like a recovered b_* draw.
  expect_false("theta_s2z" %in% excluded)
  draw_names <- c(
    "b_Intercept", "b_x", "theta_s2z", "theta_s2z[1]",
    "theta_s2z[2]", "q_recovered_s2z_1[1]"
  )
  expect_identical(
    draw_names[grepl(brms:::fixef_pars(), draw_names)],
    c("b_Intercept", "b_x")
  )
})
