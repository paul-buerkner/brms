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

# Small dense-matrix helpers for multiblock reference calculations. These
# deliberately do not call the production S2Z code.
.s2z_block_diag <- function(blocks) {
  stopifnot(length(blocks) >= 1L)
  nr <- vapply(blocks, nrow, integer(1))
  nc <- vapply(blocks, ncol, integer(1))
  out <- matrix(0, nrow = sum(nr), ncol = sum(nc))
  row_offset <- 0L
  col_offset <- 0L
  for (i in seq_along(blocks)) {
    rows <- row_offset + seq_len(nr[i])
    cols <- col_offset + seq_len(nc[i])
    out[rows, cols] <- blocks[[i]]
    row_offset <- row_offset + nr[i]
    col_offset <- col_offset + nc[i]
  }
  out
}

# Dense transformation from conventional population and group effects to one
# finite-population coefficient vector and orthonormal contrasts for every
# independent grouping factor.
.s2z_multiblock_dense_map <- function(prior_cov, H, groups, effect_cov) {
  stopifnot(
    is.matrix(prior_cov), nrow(prior_cov) == ncol(prior_cov),
    length(H) == length(groups), length(H) == length(effect_cov)
  )
  q <- nrow(prior_cov)
  m <- vapply(H, ncol, integer(1))
  stopifnot(
    all(vapply(H, nrow, integer(1)) == q),
    all(groups >= 2L),
    all(vapply(seq_along(H), function(i) {
      all(dim(effect_cov[[i]]) == groups[i] * m[i])
    }, logical(1)))
  )
  n_conventional <- q + sum(groups * m)
  n_transformed <- q + sum((groups - 1L) * m)
  transform <- matrix(
    0, nrow = n_transformed, ncol = n_conventional
  )
  transform[seq_len(q), seq_len(q)] <- diag(q)
  bases <- lapply(groups, .s2z_basis)
  col_offset <- q
  row_offset <- q
  for (i in seq_along(H)) {
    effect_columns <- col_offset + seq_len(groups[i] * m[i])
    contrast_rows <- row_offset + seq_len((groups[i] - 1L) * m[i])
    mean_map <- kronecker(
      diag(m[i]), matrix(rep(1 / groups[i], groups[i]), nrow = 1L)
    )
    contrast_map <- kronecker(diag(m[i]), t(bases[[i]]))
    transform[seq_len(q), effect_columns] <- H[[i]] %*% mean_map
    transform[contrast_rows, effect_columns] <- contrast_map
    col_offset <- col_offset + groups[i] * m[i]
    row_offset <- row_offset + (groups[i] - 1L) * m[i]
  }
  conventional_cov <- .s2z_block_diag(c(list(prior_cov), effect_cov))
  list(
    transform = transform,
    covariance = transform %*% conventional_cov %*% t(transform),
    bases = bases
  )
}

# Independent reference for the Gaussian Matheron shortcut with arbitrary
# rectangular loading maps.  The implementation below deliberately forms the
# full omitted-mean precision as a check on the lower-dimensional update; the
# production fast path must not need this matrix.
.s2z_matheron_case <- function(
    theta, prior_mean, prior_cov, H, groups, covariance,
    prior_draw = NULL) {
  stopifnot(
    length(H) == length(groups), length(H) == length(covariance),
    length(theta) == length(prior_mean),
    all(dim(prior_cov) == length(theta))
  )
  q <- length(theta)
  dimensions <- vapply(H, ncol, integer(1))
  stopifnot(
    all(vapply(H, nrow, integer(1)) == q),
    all(vapply(seq_along(H), function(i) {
      all(dim(covariance[[i]]) == dimensions[i])
    }, logical(1)))
  )
  A <- do.call(cbind, H)
  S <- lapply(seq_along(H), function(i) covariance[[i]] / groups[i])
  S_stack <- .s2z_block_diag(S)
  W <- prior_cov + A %*% S_stack %*% t(A)
  W_inv <- solve(W)
  difference <- theta - prior_mean

  # Dense completed-square reference in the stacked omitted means.
  prior_prec <- solve(prior_cov)
  P <- solve(S_stack) + crossprod(A, prior_prec %*% A)
  h <- drop(crossprod(A, prior_prec %*% difference))
  dense_mean <- drop(solve(P, h))
  dense_cov <- solve(P)

  # Matheron conditional moments use only the q by q induced covariance W.
  fast_mean <- drop(S_stack %*% t(A) %*% W_inv %*% difference)
  fast_cov <- S_stack -
    S_stack %*% t(A) %*% W_inv %*% A %*% S_stack

  out <- list(
    A = A, S = S, S_stack = S_stack, W = W,
    dense_mean = dense_mean, dense_cov = dense_cov,
    fast_mean = fast_mean, fast_cov = fast_cov
  )
  if (!is.null(prior_draw)) {
    stopifnot(length(prior_draw) == q + sum(dimensions))
    beta_star <- prior_draw[seq_len(q)]
    m_star <- prior_draw[q + seq_len(sum(dimensions))]
    theta_star <- drop(beta_star + A %*% m_star)
    delta <- drop(W_inv %*% (theta - theta_star))
    beta_recovered <- drop(beta_star + prior_cov %*% delta)
    m_recovered <- drop(m_star + S_stack %*% t(A) %*% delta)
    out$beta_recovered <- beta_recovered
    out$m_recovered <- m_recovered
    out$theta_recovered <- drop(beta_recovered + A %*% m_recovered)
  }
  out
}

# Compare the component-wise diagonal-plus-rank-one implementation with the
# dense complete square it replaces.  This deliberately mirrors the scaled
# equations in .stan_re_s2z_independent without calling code under test.
.s2z_independent_case <- function(k, student = FALSE,
                                  has_intercept = TRUE) {
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
  intercept_map <- if (has_intercept) {
    c(1, seq(-0.55, 0.65, length.out = k - 1L))
  } else {
    numeric(k)
  }
  H <- diag(k)
  if (has_intercept) {
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
  if (has_intercept) {
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
    ans <- .s2z_independent_case(
      k, student = FALSE, has_intercept = TRUE
    )

    expect_equal(colSums(ans$r_s2z), numeric(k), tolerance = 2e-14)
    expect_equal(ans$dense_mode, ans$fast_mode, tolerance = 2e-13)
    expect_equal(ans$dense_qhat, ans$fast_qhat, tolerance = 2e-13)
    expect_equal(ans$dense_log, ans$fast_log, tolerance = 2e-10)
    expect_gt(ans$rank1_info, 0)
  }
})

test_that("independent Student effects match the weighted dense square", {
  for (k in c(2L, 5L, 10L)) {
    ans <- .s2z_independent_case(
      k, student = TRUE, has_intercept = TRUE
    )

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
      ans <- .s2z_independent_case(
        k, student = student, has_intercept = TRUE
      )
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
      ans <- .s2z_independent_case(
        k, student = student, has_intercept = FALSE
      )
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

test_that("multiple scalar factor means match the closed-form recovery law", {
  groups <- c(3L, 5L, 8L)
  tau <- c(0.45, 1.10, 0.72)
  n_block <- length(groups)
  v <- tau^2 / groups
  total_v <- sum(v)
  alpha <- v / (1 + total_v)
  expected_conditional_cov <- diag(v) - tcrossprod(v) / (1 + total_v)

  # This is exactly the scalar formula in the multifactored derivation:
  # theta = mu_pop + sum(m), mu_pop ~ N(0, 1), and m_k ~ N(0, v_k).
  P <- diag(1 / v) + matrix(1, nrow = n_block, ncol = n_block)
  conditional_cov <- solve(P)
  theta <- 0.83
  conditional_mean <- drop(conditional_cov %*% rep(theta, n_block))
  expect_equal(conditional_mean, alpha * theta, tolerance = 2e-14)
  expect_equal(
    conditional_cov, expected_conditional_cov, tolerance = 2e-14
  )
  expect_true(all(conditional_cov[row(conditional_cov) !=
                                   col(conditional_cov)] < 0))

  # Jointly drawing theta and then the conditional means must recover the
  # original mutually independent population mean and factor means.
  joint_cov <- rbind(
    c(1 + total_v, v),
    cbind(v, diag(v))
  )
  recovery_map <- rbind(
    c(1, rep(-1, n_block)),
    cbind(rep(0, n_block), diag(n_block))
  )
  recovered_cov <- recovery_map %*% joint_cov %*% t(recovery_map)
  expect_equal(recovered_cov, diag(c(1, v)), tolerance = 2e-14)

  # Independently compare the integrated density with a dense transformation
  # of all conventional raw effects. This also checks every sqrt(G_k) measure
  # correction rather than only the conditional moments above.
  H <- replicate(n_block, matrix(1, nrow = 1L), simplify = FALSE)
  effect_cov <- Map(
    function(n, scale) diag(scale^2, nrow = n), groups, tau
  )
  dense <- .s2z_multiblock_dense_map(
    matrix(1, nrow = 1L), H, groups, effect_cov
  )
  r_s2z <- lapply(seq_along(groups), function(i) {
    z <- matrix(
      sin(seq_len(groups[i] - 1L) * (0.31 + 0.07 * i)),
      ncol = 1L
    )
    dense$bases[[i]] %*% z
  })
  transformed_value <- c(
    theta,
    unlist(Map(
      function(basis, r) as.vector(crossprod(basis, r)),
      dense$bases, r_s2z
    ))
  )
  transformed_mean <- drop(dense$transform %*% numeric(1L + sum(groups)))
  log_direct <- .s2z_log_mvn(
    transformed_value, transformed_mean, dense$covariance
  )
  mhat <- conditional_mean
  log_joint_at_mode <- stats::dnorm(
    theta - sum(mhat), 0, 1, log = TRUE
  ) + sum(vapply(seq_along(groups), function(i) {
    sum(stats::dnorm(
      r_s2z[[i]][, 1L] + mhat[i], 0, tau[i], log = TRUE
    ))
  }, numeric(1))) + 0.5 * sum(log(groups))
  log_integrated <- log_joint_at_mode +
    0.5 * n_block * log(2 * pi) - sum(log(diag(chol(P))))
  expect_equal(log_integrated, log_direct, tolerance = 2e-11)
})

test_that("large scalar Gaussian systems have exact Matheron recovery", {
  n_factor <- 31L
  groups <- 2L + (seq_len(n_factor) %% 13L)
  tau <- 0.25 + (seq_len(n_factor) %% 9L) / 8
  H <- replicate(n_factor, matrix(1, nrow = 1L), simplify = FALSE)
  covariance <- lapply(tau, function(x) matrix(x^2, nrow = 1L))
  theta <- 0.87
  prior_mean <- -0.23
  prior_cov <- matrix(1.7^2, nrow = 1L)
  v <- tau^2 / groups

  beta_star <- 0.41
  m_star <- sin(seq_len(n_factor) * 0.37) * sqrt(v)
  ans <- .s2z_matheron_case(
    theta, prior_mean, prior_cov, H, groups, covariance,
    prior_draw = c(beta_star, m_star)
  )
  W <- drop(prior_cov) + sum(v)
  expected_mean <- v / W * (theta - prior_mean)
  expected_cov <- diag(v) - tcrossprod(v) / W

  expect_equal(ans$W, matrix(W, nrow = 1L), tolerance = 2e-14)
  expect_equal(ans$fast_mean, expected_mean, tolerance = 2e-13)
  expect_equal(ans$dense_mean, expected_mean, tolerance = 2e-13)
  expect_equal(ans$fast_cov, expected_cov, tolerance = 3e-13)
  expect_equal(ans$dense_cov, expected_cov, tolerance = 3e-13)
  expect_true(all(ans$fast_cov[row(ans$fast_cov) !=
                                 col(ans$fast_cov)] < 0))
  expect_equal(ans$theta_recovered, theta, tolerance = 3e-14)

  # This explicit update is the scalar Matheron algorithm.  It needs only W,
  # even though the omitted means have dimension 31.
  theta_star <- beta_star + sum(m_star)
  delta <- (theta - theta_star) / W
  expect_equal(
    ans$beta_recovered, beta_star + drop(prior_cov) * delta,
    tolerance = 2e-14
  )
  expect_equal(ans$m_recovered, m_star + v * delta, tolerance = 2e-14)
})

test_that("rectangular overlapping blocks use exact Matheron moments and density", {
  q <- 5L
  groups <- c(3L, 8L, 5L, 11L)
  dimensions <- c(2L, 1L, 3L, 2L)
  H <- lapply(dimensions, function(m) matrix(0, nrow = q, ncol = m))
  H[[1L]][cbind(c(1L, 2L, 1L), c(1L, 1L, 2L))] <- c(1, 1, 0.31)
  H[[1L]][3L, 2L] <- 1
  H[[2L]][c(1L, 2L), 1L] <- c(-0.27, 1)
  H[[3L]][1L, ] <- c(1, 0.18, -0.22)
  H[[3L]][cbind(c(2L, 4L, 5L), seq_len(3L))] <- 1
  H[[4L]][1L, ] <- c(1, 0.36)
  H[[4L]][cbind(c(3L, 5L), seq_len(2L))] <- 1
  A <- do.call(cbind, H)
  expect_gt(sum(rowSums(A != 0) > 1L), 2L)

  prior_mean <- c(-0.4, 0.25, 0.15, -0.6, 0.75)
  prior_cov <- diag(c(0.8, 1.1, 0.65, 1.35, 0.9)^2)
  L <- list(
    matrix(c(0.62, 0, 0.17, 0.91), nrow = 2L, byrow = TRUE),
    matrix(0.73, nrow = 1L),
    matrix(c(
      0.55, 0, 0, -0.12, 0.84, 0, 0.20, 0.09, 0.68
    ), nrow = 3L, byrow = TRUE),
    matrix(c(1.02, 0, -0.23, 0.59), nrow = 2L, byrow = TRUE)
  )
  covariance <- lapply(L, tcrossprod)
  theta <- c(0.95, -0.35, 0.57, 0.12, -0.81)
  prior_draw <- c(
    prior_mean + c(0.14, -0.31, 0.22, 0.17, -0.08),
    sin(seq_len(sum(dimensions)) * 0.43) / 3
  )
  ans <- .s2z_matheron_case(
    theta, prior_mean, prior_cov, H, groups, covariance,
    prior_draw = prior_draw
  )

  expect_equal(ans$fast_mean, ans$dense_mean, tolerance = 3e-13)
  expect_equal(ans$fast_cov, ans$dense_cov, tolerance = 4e-13)
  expect_equal(ans$theta_recovered, theta, tolerance = 3e-14)

  # Independently construct the covariance of the shared Matheron update.
  C <- .s2z_block_diag(list(prior_cov, ans$S_stack))
  B <- cbind(diag(q), A)
  update <- diag(nrow(C)) - C %*% t(B) %*% solve(ans$W, B)
  recovered_cov <- update %*% C %*% t(update)
  expected_cov <- C - C %*% t(B) %*% solve(ans$W, B %*% C)
  expect_equal(recovered_cov, expected_cov, tolerance = 5e-13)
  expect_equal(
    recovered_cov[q + seq_len(sum(dimensions)),
                  q + seq_len(sum(dimensions))],
    ans$dense_cov, tolerance = 5e-13
  )
  expect_equal(B %*% update, matrix(0, nrow = q, ncol = nrow(C)),
               tolerance = 3e-14)

  # The separated Matheron density equals the dense pushforward of ordinary
  # Gaussian population and level effects, including unequal group sizes.
  effect_cov <- Map(
    function(n, sigma) kronecker(sigma, diag(n)), groups, covariance
  )
  dense <- .s2z_multiblock_dense_map(prior_cov, H, groups, effect_cov)
  contrast <- lapply(seq_along(groups), function(i) {
    matrix(
      sin(seq_len((groups[i] - 1L) * dimensions[i]) *
            (0.19 + 0.03 * i)),
      nrow = groups[i] - 1L, ncol = dimensions[i]
    ) / 2
  })
  transformed_value <- c(theta, unlist(contrast))
  transformed_mean <- drop(dense$transform %*%
    c(prior_mean, numeric(sum(groups * dimensions))))
  log_dense <- .s2z_log_mvn(
    transformed_value, transformed_mean, dense$covariance
  )
  log_fast <- .s2z_log_mvn(theta, prior_mean, ans$W) +
    sum(vapply(seq_along(groups), function(i) {
      .s2z_log_mvn(
        as.vector(contrast[[i]]), numeric((groups[i] - 1L) * dimensions[i]),
        kronecker(covariance[[i]], diag(groups[i] - 1L))
      )
    }, numeric(1)))
  expect_equal(log_fast, log_dense, tolerance = 5e-10)
})

test_that("Matheron conditions only proper active population coordinates", {
  groups <- c(4L, 7L, 5L)
  dimensions <- c(2L, 1L, 2L)
  H <- list(
    matrix(c(1, 0.25, 1, 0, 0, 0, 0, -0.3, 0, 1), nrow = 5L),
    matrix(c(-0.2, 1, 0, 0, 0), nrow = 5L),
    matrix(c(0.4, 0, 0, 0, 1, -0.1, 0, 0, 0, 0.7), nrow = 5L)
  )
  A <- do.call(cbind, H)
  covariance <- list(
    matrix(c(0.64, 0.12, 0.12, 0.49), nrow = 2L),
    matrix(0.81, nrow = 1L),
    matrix(c(0.36, -0.08, -0.08, 0.72), nrow = 2L)
  )
  S <- .s2z_block_diag(Map(
    function(sigma, n) sigma / n, covariance, groups
  ))
  # Rows 1 and 4 have proper priors, but row 4 is inactive. Rows 2 and 5
  # are active and flat. Thus the shared conditioning system has rank one.
  proper <- c(TRUE, FALSE, FALSE, TRUE, FALSE)
  active_proper <- which(proper & rowSums(A != 0) > 0)
  inactive_proper <- which(proper & rowSums(A != 0) == 0)
  expect_identical(active_proper, 1L)
  expect_identical(inactive_proper, 4L)
  prior_mean <- c(-0.35, NA, NA, 0.42, NA)
  prior_var <- c(1.3^2, NA, NA, 0.75^2, NA)
  theta <- c(0.8, -0.2, 0.55, 0.1, -0.65)

  A_active <- A[active_proper, , drop = FALSE]
  W <- prior_var[active_proper] +
    drop(A_active %*% S %*% t(A_active))
  expected_mean <- drop(
    S %*% t(A_active) / W *
      (theta[active_proper] - prior_mean[active_proper])
  )
  expected_cov <- S -
    S %*% t(A_active) %*% A_active %*% S / W

  # Dense conditional reference excludes flat and proper-inactive rows.
  P <- solve(S) + crossprod(A_active) / prior_var[active_proper]
  h <- drop(t(A_active) *
    (theta[active_proper] - prior_mean[active_proper]) /
    prior_var[active_proper])
  expect_equal(drop(solve(P, h)), expected_mean, tolerance = 2e-13)
  expect_equal(solve(P), expected_cov, tolerance = 3e-13)

  # The normalized density on the proper theta coordinates and every
  # contrast agrees with direct integration over the omitted means. Flat
  # active theta rows correctly contribute only their Lebesgue constant.
  bases <- lapply(groups, .s2z_basis)
  contrast <- lapply(seq_along(groups), function(i) {
    matrix(
      cos(seq_len((groups[i] - 1L) * dimensions[i]) *
            (0.22 + 0.04 * i)),
      nrow = groups[i] - 1L, ncol = dimensions[i]
    ) / 2.4
  })
  r_s2z <- Map(`%*%`, bases, contrast)
  mhat <- expected_mean
  offsets <- cumsum(c(0L, dimensions))
  log_joint_mode <- sum(stats::dnorm(
    (theta - A %*% mhat)[proper], prior_mean[proper],
    sqrt(prior_var[proper]), log = TRUE
  ))
  for (i in seq_along(groups)) {
    mi <- mhat[offsets[i] + seq_len(dimensions[i])]
    log_joint_mode <- log_joint_mode + sum(vapply(
      seq_len(groups[i]),
      function(j) .s2z_log_mvn(
        r_s2z[[i]][j, ] + mi, numeric(dimensions[i]), covariance[[i]]
      ),
      numeric(1)
    )) + 0.5 * dimensions[i] * log(groups[i])
  }
  log_integrated <- log_joint_mode +
    0.5 * sum(dimensions) * log(2 * pi) -
    0.5 * as.numeric(determinant(P, logarithm = TRUE)$modulus)
  log_fast <- stats::dnorm(
    theta[active_proper], prior_mean[active_proper], sqrt(W), log = TRUE
  ) + stats::dnorm(
    theta[inactive_proper], prior_mean[inactive_proper],
    sqrt(prior_var[inactive_proper]), log = TRUE
  ) + sum(vapply(seq_along(groups), function(i) {
    .s2z_log_mvn(
      as.vector(contrast[[i]]),
      numeric((groups[i] - 1L) * dimensions[i]),
      kronecker(covariance[[i]], diag(groups[i] - 1L))
    )
  }, numeric(1)))
  expect_equal(log_fast, log_integrated, tolerance = 4e-10)

  beta_star <- prior_mean[active_proper] + 0.27
  m_star <- cos(seq_len(sum(dimensions)) * 0.31) / 4
  theta_star <- beta_star + drop(A_active %*% m_star)
  delta <- (theta[active_proper] - theta_star) / W
  m_recovered <- drop(m_star + S %*% t(A_active) * delta)
  beta_recovered <- theta - drop(A %*% m_recovered)
  beta_recovered[active_proper] <- beta_star +
    prior_var[active_proper] * delta
  beta_recovered[inactive_proper] <- theta[inactive_proper]
  expect_equal(
    beta_recovered + drop(A %*% m_recovered), theta,
    tolerance = 3e-14
  )

  # With no proper active coordinate, theta supplies no information about the
  # omitted means: draw them directly and recover every active flat beta.
  proper_r0 <- c(FALSE, FALSE, FALSE, TRUE, FALSE)
  expect_length(which(proper_r0 & rowSums(A != 0) > 0), 0L)
  no_active_proper <- which(proper_r0 & rowSums(A != 0) == 0)
  expect_identical(no_active_proper, 4L)
  P_r0 <- solve(S)
  expect_equal(
    drop(solve(P_r0, numeric(sum(dimensions)))),
    numeric(sum(dimensions)), tolerance = 2e-14
  )
  expect_equal(solve(P_r0), S, tolerance = 2e-14)
  beta_r0 <- theta - drop(A %*% m_star)
  beta_r0[no_active_proper] <- theta[no_active_proper]
  expect_equal(beta_r0 + drop(A %*% m_star), theta, tolerance = 3e-14)
})

test_that("multiple vector blocks use one overlapping completed square", {
  q <- 4L
  groups <- c(3L, 5L, 4L)
  dimensions <- c(2L, 1L, 3L)
  H <- list(
    matrix(0, nrow = q, ncol = dimensions[1L]),
    matrix(0, nrow = q, ncol = dimensions[2L]),
    matrix(0, nrow = q, ncol = dimensions[3L])
  )
  H[[1L]][1L, 1L] <- 1
  H[[1L]][c(1L, 2L), 2L] <- c(0.35, 1)
  H[[2L]][c(1L, 2L), 1L] <- c(-0.20, 1)
  H[[3L]][1L, 1L] <- 1
  H[[3L]][c(1L, 2L), 2L] <- c(0.10, 1)
  H[[3L]][c(1L, 4L), 3L] <- c(-0.30, 1)
  A <- do.call(cbind, H)
  expect_gte(sum(A[2L, ] != 0), 3L)

  prior_mean <- c(0.30, -0.45, 0.15, 0.70)
  prior_cov <- diag(c(0.75, 1.10, 0.90, 1.35)^2)
  prior_prec <- solve(prior_cov)
  L3 <- rbind(
    c(0.80, 0, 0),
    c(0.15, 0.65, 0),
    c(-0.10, 0.12, 0.55)
  )
  Sigma <- list(
    matrix(c(0.81, 0.18, 0.18, 0.49), nrow = 2L),
    matrix(0.64, nrow = 1L),
    tcrossprod(L3)
  )
  effect_cov <- Map(
    function(n, sigma) kronecker(sigma, diag(n)), groups, Sigma
  )
  dense <- .s2z_multiblock_dense_map(
    prior_cov, H, groups, effect_cov
  )
  r_s2z <- lapply(seq_along(groups), function(i) {
    z <- matrix(
      sin(seq_len((groups[i] - 1L) * dimensions[i]) *
            (0.21 + 0.08 * i)) +
        cos(seq_len((groups[i] - 1L) * dimensions[i]) * 0.17),
      nrow = groups[i] - 1L, ncol = dimensions[i]
    ) / 2
    dense$bases[[i]] %*% z
  })
  expect_equal(
    unlist(lapply(r_s2z, colSums)), numeric(sum(dimensions)),
    tolerance = 2e-14
  )
  theta <- c(0.95, -0.20, 0.55, -0.65)

  group_information <- Map(
    function(n, sigma) n * solve(sigma), groups, Sigma
  )
  group_score <- Map(
    function(r, sigma) solve(sigma, colSums(r)), r_s2z, Sigma
  )
  P <- crossprod(A, prior_prec %*% A) +
    .s2z_block_diag(group_information)
  h <- drop(crossprod(A, prior_prec %*% (theta - prior_mean))) -
    unlist(group_score)
  mhat <- drop(solve(P, h))

  block_id <- rep(seq_along(dimensions), dimensions)
  cross_block <- outer(block_id, block_id, `!=`)
  expect_gt(max(abs(P[cross_block])), 0.1)

  offsets <- cumsum(c(0L, dimensions))
  log_joint_at_mode <- .s2z_log_mvn(
    theta - A %*% mhat, prior_mean, prior_cov
  )
  for (i in seq_along(groups)) {
    mi <- mhat[offsets[i] + seq_len(dimensions[i])]
    log_joint_at_mode <- log_joint_at_mode + sum(vapply(
      seq_len(groups[i]),
      function(j) .s2z_log_mvn(
        r_s2z[[i]][j, ] + mi, numeric(dimensions[i]), Sigma[[i]]
      ),
      numeric(1)
    ))
  }
  log_integrated <- log_joint_at_mode +
    0.5 * sum(dimensions * log(groups)) +
    0.5 * sum(dimensions) * log(2 * pi) -
    sum(log(diag(chol(P))))

  transformed_value <- c(
    theta,
    unlist(Map(
      function(basis, r) as.vector(crossprod(basis, r)),
      dense$bases, r_s2z
    ))
  )
  conventional_mean <- c(
    prior_mean, numeric(sum(groups * dimensions))
  )
  transformed_mean <- drop(dense$transform %*% conventional_mean)
  log_direct <- .s2z_log_mvn(
    transformed_value, transformed_mean, dense$covariance
  )
  expect_equal(log_integrated, log_direct, tolerance = 3e-10)
})

test_that("shared-scale multiblock covariances share one exact integral", {
  q <- 3L
  groups <- c(4L, 3L, 5L)
  dimensions <- c(2L, 1L, 2L)
  H <- list(
    matrix(c(1, 0, 0, 0.25, 1, 0), nrow = q),
    matrix(c(-0.15, 1, 0), nrow = q),
    matrix(c(1, 0, 0, -0.30, 0, 1), nrow = q)
  )
  A <- do.call(cbind, H)
  prior_mean <- c(-0.20, 0.55, 0.10)
  prior_cov <- diag(c(0.65, 1.20, 0.85)^2)
  prior_prec <- solve(prior_cov)

  cor_1 <- matrix(c(1, 0.35, 0.35, 1), nrow = 2L)
  cor_3 <- matrix(c(1, -0.28, -0.28, 1), nrow = 2L)
  Sigma <- list(
    diag(c(0.64, 0.78)) %*% cor_1 %*% diag(c(0.64, 0.78)),
    matrix(0.71^2, nrow = 1L),
    diag(c(0.72, 0.59)) %*% cor_3 %*% diag(c(0.72, 0.59))
  )
  effect_cov <- Map(
    function(n, sigma) kronecker(sigma, diag(n)), groups, Sigma
  )
  dense <- .s2z_multiblock_dense_map(
    prior_cov, H, groups, effect_cov
  )
  r_s2z <- lapply(seq_along(groups), function(i) {
    z <- matrix(
      sin(seq_len((groups[i] - 1L) * dimensions[i]) *
            (0.34 + 0.05 * i)) -
        cos(seq_len((groups[i] - 1L) * dimensions[i]) * 0.29),
      nrow = groups[i] - 1L, ncol = dimensions[i]
    ) / 1.7
    dense$bases[[i]] %*% z
  })
  expect_equal(
    unlist(lapply(r_s2z, colSums)), numeric(sum(dimensions)),
    tolerance = 2e-14
  )

  group_information <- vector("list", length(groups))
  group_score <- vector("list", length(groups))
  for (i in seq_along(groups)) {
    precision <- solve(Sigma[[i]])
    group_information[[i]] <- groups[i] * precision
    group_score[[i]] <- precision %*% colSums(r_s2z[[i]])
  }
  expect_equal(
    unlist(group_score), numeric(sum(dimensions)), tolerance = 2e-13
  )

  theta <- c(0.85, -0.35, 0.62)
  P <- crossprod(A, prior_prec %*% A) +
    .s2z_block_diag(group_information)
  h <- drop(crossprod(A, prior_prec %*% (theta - prior_mean))) -
    unlist(group_score)
  mhat <- drop(solve(P, h))
  conditional_cov <- solve(P)
  offsets <- cumsum(c(0L, dimensions))
  log_joint_measure <- function(group_mean) {
    out <- .s2z_log_mvn(
      theta - A %*% group_mean, prior_mean, prior_cov
    )
    for (i in seq_along(groups)) {
      mi <- group_mean[offsets[i] + seq_len(dimensions[i])]
      out <- out + sum(vapply(seq_len(groups[i]), function(j) {
        .s2z_log_mvn(
          r_s2z[[i]][j, ] + mi, numeric(dimensions[i]),
          Sigma[[i]]
        )
      }, numeric(1)))
    }
    out + 0.5 * sum(dimensions * log(groups))
  }
  log_marginal <- log_joint_measure(mhat) +
    0.5 * sum(dimensions) * log(2 * pi) -
    sum(log(diag(chol(P))))

  candidate_means <- list(
    mhat,
    mhat + seq(-0.30, 0.25, length.out = sum(dimensions)),
    c(-0.80, 0.35, 0.20, -0.45, 0.70)
  )
  for (group_mean in candidate_means) {
    by_identity <- log_joint_measure(group_mean) -
      .s2z_log_mvn(group_mean, mhat, conditional_cov)
    expect_equal(by_identity, log_marginal, tolerance = 3e-10)
  }

  transformed_value <- c(
    theta,
    unlist(Map(
      function(basis, r) as.vector(crossprod(basis, r)),
      dense$bases, r_s2z
    ))
  )
  conventional_mean <- c(
    prior_mean, numeric(sum(groups * dimensions))
  )
  transformed_mean <- drop(dense$transform %*% conventional_mean)
  log_direct <- .s2z_log_mvn(
    transformed_value, transformed_mean, dense$covariance
  )
  expect_equal(log_marginal, log_direct, tolerance = 4e-10)
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
