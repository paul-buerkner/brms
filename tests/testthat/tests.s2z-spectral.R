context("Tests for spectral S2Z covariance geometry")

# Independent R references for the proposed spectral known-covariance path.
# These deliberately do not call the production S2Z helpers: their purpose is
# to check the matrix identities used by generated Stan code.
.s2z_spectral_constrain <- function(y) {
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

.s2z_spectral_basis <- function(n) {
  stopifnot(n >= 2L)
  vapply(
    seq_len(n - 1L),
    function(j) .s2z_spectral_constrain(diag(n - 1L)[, j]),
    numeric(n)
  )
}

.s2z_spectral_case <- function(G = 6L, M = 3L) {
  stopifnot(G >= 3L, M >= 1L)
  B <- .s2z_spectral_basis(G)

  # The unequal diagonal is intentional: Omega * 1 is not proportional to 1,
  # so the rank-one correction in the modal quadratic is nonzero.
  X <- matrix(
    sin(seq_len(G * G) * 0.37) + cos(seq_len(G * G) * 0.19),
    G, G
  )
  Omega <- tcrossprod(X) / G + diag(seq(0.7, 1.4, length.out = G))
  restricted <- crossprod(B, Omega %*% B)
  eig <- eigen(restricted, symmetric = TRUE)
  E <- B %*% eig$vectors

  Y <- matrix(
    sin(seq_len(M * M) * 0.29) - cos(seq_len(M * M) * 0.43),
    M, M
  )
  Sigma <- tcrossprod(Y) / M + diag(seq(0.8, 1.3, length.out = M))
  L <- t(chol(Sigma))
  T <- matrix(
    sin(seq_len((G - 1L) * M) * 0.31) +
      0.4 * cos(seq_len((G - 1L) * M) * 0.53),
    G - 1L, M
  )
  one <- rep(1, G)
  Omega_inv_one <- solve(Omega, one)

  list(
    G = G, M = M, B = B, Omega = Omega,
    restricted = restricted, E = E, lambda = eig$values,
    Sigma = Sigma, L = L, T = T, delta = E %*% T,
    kappa = drop(crossprod(one, Omega_inv_one)),
    c = drop(crossprod(E, Omega_inv_one))
  )
}

.s2z_spectral_group_terms <- function(case) {
  one <- rep(1, case$G)
  delta_white <- case$delta %*% solve(t(case$L))
  W <- case$T %*% solve(t(case$L))

  q_dense <- sum(delta_white * solve(case$Omega, delta_white))
  h_dense <- -drop(crossprod(
    one, solve(case$Omega, delta_white)
  ))
  P_dense <- case$kappa * diag(case$M)

  modal_cross <- drop(crossprod(W, case$c))
  q_diagonal <- sum(rowSums(W^2) / case$lambda)
  list(
    delta_white = delta_white,
    W = W,
    q_dense = q_dense,
    h_dense = h_dense,
    P_dense = P_dense,
    q_diagonal = q_diagonal,
    q_modal = q_diagonal + sum(modal_cross^2) / case$kappa,
    h_modal = -modal_cross,
    P_modal = case$kappa * diag(case$M)
  )
}

.s2z_spectral_modal_transform <- function(z, lambda, L, rho) {
  stopifnot(
    is.matrix(z), all(dim(rho) == dim(z)),
    length(lambda) == nrow(z), ncol(z) == nrow(L),
    all(dim(L) == ncol(z))
  )
  out <- matrix(NA_real_, nrow(z), ncol(z))
  log_jacobian <- 0
  for (ell in seq_len(nrow(z))) {
    S <- sqrt(lambda[ell]) * L
    R <- diag(rho[ell, ], ncol(z))
    D <- diag(ncol(z)) - R + R %*% S
    out[ell, ] <- drop(S %*% solve(D, z[ell, ]))
    log_jacobian <- log_jacobian +
      as.numeric(determinant(S, logarithm = TRUE)$modulus) -
      as.numeric(determinant(D, logarithm = TRUE)$modulus)
  }
  list(value = out, log_jacobian = log_jacobian)
}

.s2z_spectral_fd_jacobian <- function(fun, x, epsilon = 1e-6) {
  out <- matrix(NA_real_, length(fun(x)), length(x))
  for (j in seq_along(x)) {
    step <- numeric(length(x))
    step[j] <- epsilon
    out[, j] <- (fun(x + step) - fun(x - step)) / (2 * epsilon)
  }
  out
}

test_that("restricted eigensystem diagonalizes projected covariance", {
  case <- .s2z_spectral_case()
  P <- diag(case$G) - matrix(1 / case$G, case$G, case$G)

  expect_equal(crossprod(case$B), diag(case$G - 1L), tolerance = 1e-14)
  expect_equal(
    drop(crossprod(case$B, rep(1, case$G))),
    numeric(case$G - 1L), tolerance = 1e-14
  )
  expect_equal(crossprod(case$E), diag(case$G - 1L), tolerance = 1e-13)
  expect_equal(
    drop(crossprod(case$E, rep(1, case$G))),
    numeric(case$G - 1L), tolerance = 1e-13
  )
  expect_true(all(case$lambda > 0))
  expect_equal(
    case$E %*% diag(case$lambda) %*% t(case$E),
    P %*% case$Omega %*% P, tolerance = 2e-13
  )

  # Restricting Omega^{-1} is not the inverse of restricting Omega. The
  # difference is exactly the rank-one omitted-mean correction.
  expect_equal(
    crossprod(case$E, solve(case$Omega, case$E)),
    diag(1 / case$lambda) + tcrossprod(case$c) / case$kappa,
    tolerance = 2e-13
  )
  expect_gt(sum(abs(case$c)), 1e-3)
})

test_that("modal q h and P equal dense known-covariance calculations", {
  case <- .s2z_spectral_case()
  terms <- .s2z_spectral_group_terms(case)

  expect_equal(terms$q_modal, terms$q_dense, tolerance = 2e-12)
  expect_equal(terms$h_modal, terms$h_dense, tolerance = 2e-12)
  expect_equal(terms$P_modal, terms$P_dense, tolerance = 1e-14)
  expect_gt(abs(terms$q_dense - terms$q_diagonal), 1e-4)

  # Check the entire omitted-mean quadratic, not only its coefficients.
  x <- seq(-0.35, 0.45, length.out = case$M)
  shifted <- terms$delta_white + outer(rep(1, case$G), x)
  dense_quadratic <- sum(shifted * solve(case$Omega, shifted))
  modal_quadratic <- terms$q_modal -
    2 * drop(crossprod(terms$h_modal, x)) +
    drop(crossprod(x, terms$P_modal %*% x))
  expect_equal(modal_quadratic, dense_quadratic, tolerance = 2e-12)
})

test_that("two-coefficient modal fast paths equal dense matrix operations", {
  case <- .s2z_spectral_case(G = 7L, M = 2L)
  J_seed <- matrix(c(0.9, -0.35, 0.2, 1.1), 2L, 2L)
  J <- tcrossprod(J_seed) + diag(c(0.4, 0.7))
  unit_info <- crossprod(case$L, J %*% case$L)
  prior_var <- diag(case$Sigma)
  exposure <- seq(0.6, 1.4, length.out = case$G - 1L)
  rho <- matrix(NA_real_, case$G - 1L, 2L)

  for (ell in seq_len(case$G - 1L)) {
    scaled_lambda <- exposure[ell] * case$lambda[ell]
    a <- 1 + scaled_lambda * unit_info[1, 1]
    inverse_a <- 1 / a
    scaled_over_a <- scaled_lambda * inverse_a
    regression <- scaled_over_a * unit_info[1, 2]
    conditional_unit_info <- max(
      0, unit_info[2, 2] - unit_info[1, 2]^2 / unit_info[1, 1]
    )
    conditional_precision <- 1 +
      scaled_over_a * unit_info[2, 2] +
      scaled_over_a * unit_info[1, 1] * scaled_lambda *
        conditional_unit_info
    post_ratio1 <- inverse_a + regression^2 / conditional_precision
    post_var2 <- case$L[2, 1]^2 / a +
      (case$L[2, 2] - regression * case$L[2, 1])^2 /
        conditional_precision
    rho[ell, ] <- c(
      1 - post_ratio1,
      1 - post_var2 / prior_var[2]
    )

    precision <- diag(2L) + scaled_lambda * unit_info
    posterior <- case$lambda[ell] *
      case$L %*% solve(precision, t(case$L))
    expect_equal(
      rho[ell, ], 1 - diag(posterior) /
        (case$lambda[ell] * prior_var),
      tolerance = 2e-14
    )
  }

  # Preserve the identity contribution when scaled information is huge and
  # rank one. Directly subtracting the two O(s^2) determinant terms yields
  # zero in double precision for this fixture, although the exact Schur
  # complement tends to two.
  extreme_scale <- 1e16
  extreme_a <- 1 + extreme_scale
  extreme_scaled_over_a <- extreme_scale / extreme_a
  extreme_conditional <- max(0, 1 - 1^2 / 1)
  extreme_q <- 1 + extreme_scaled_over_a +
    extreme_scaled_over_a * extreme_scale * extreme_conditional
  expect_equal(extreme_q, 2, tolerance = 1e-15)

  z <- matrix(
    cos(seq_len((case$G - 1L) * 2L) * 0.37), case$G - 1L, 2L
  )
  generic <- .s2z_spectral_modal_transform(z, case$lambda, case$L, rho)
  specialized <- matrix(NA_real_, nrow(z), 2L)
  specialized_log_jacobian <- 0
  for (ell in seq_len(nrow(z))) {
    sqrt_lambda <- sqrt(case$lambda[ell])
    s11 <- sqrt_lambda * case$L[1, 1]
    s21 <- sqrt_lambda * case$L[2, 1]
    s22 <- sqrt_lambda * case$L[2, 2]
    d11 <- 1 - rho[ell, 1] + rho[ell, 1] * s11
    d21 <- rho[ell, 2] * s21
    d22 <- 1 - rho[ell, 2] + rho[ell, 2] * s22
    w1 <- z[ell, 1] / d11
    w2 <- (z[ell, 2] - d21 * w1) / d22
    specialized[ell, ] <- c(s11 * w1, s21 * w1 + s22 * w2)
    specialized_log_jacobian <- specialized_log_jacobian +
      0.5 * 2L * log(case$lambda[ell]) + sum(log(diag(case$L))) -
      log(d11) - log(d22)
  }
  expect_equal(specialized, generic$value, tolerance = 2e-14)
  expect_equal(
    specialized_log_jacobian, generic$log_jacobian, tolerance = 2e-14
  )
})

test_that("modal partial-centering Jacobian has exact endpoints", {
  case <- .s2z_spectral_case(G = 5L, M = 3L)
  modes <- case$G - 1L
  z <- matrix(
    sin(seq_len(modes * case$M) * 0.41), modes, case$M
  )

  noncentered <- .s2z_spectral_modal_transform(
    z, case$lambda, case$L, matrix(0, modes, case$M)
  )
  expected_noncentered <- matrix(NA_real_, modes, case$M)
  for (ell in seq_len(modes)) {
    expected_noncentered[ell, ] <-
      drop(sqrt(case$lambda[ell]) * case$L %*% z[ell, ])
  }
  expected_noncentered_jacobian <-
    case$M / 2 * sum(log(case$lambda)) +
    modes * sum(log(diag(case$L)))
  expect_equal(noncentered$value, expected_noncentered, tolerance = 1e-13)
  expect_equal(
    noncentered$log_jacobian, expected_noncentered_jacobian,
    tolerance = 1e-13
  )

  centered <- .s2z_spectral_modal_transform(
    z, case$lambda, case$L, matrix(1, modes, case$M)
  )
  expect_equal(centered$value, z, tolerance = 1e-13)
  expect_equal(centered$log_jacobian, 0, tolerance = 1e-13)

  rho <- matrix(
    seq(0.08, 0.92, length.out = modes * case$M),
    modes, case$M
  )
  partial <- .s2z_spectral_modal_transform(
    z, case$lambda, case$L, rho
  )
  map <- function(z_vector) {
    z_matrix <- matrix(z_vector, modes, case$M)
    transformed <- .s2z_spectral_modal_transform(
      z_matrix, case$lambda, case$L, rho
    )$value
    # Include both synthesis into the physical S2Z space and projection back
    # to square restricted coordinates in the finite-difference check.
    as.vector(crossprod(case$E, case$E %*% transformed))
  }
  numeric_jacobian <- .s2z_spectral_fd_jacobian(map, as.vector(z))
  numeric_log_jacobian <- as.numeric(
    determinant(numeric_jacobian, logarithm = TRUE)$modulus
  )
  expect_equal(
    numeric_log_jacobian, partial$log_jacobian,
    tolerance = 2e-8
  )
})

test_that("nonflat population priors complete the omitted-mean square", {
  case <- .s2z_spectral_case(G = 6L, M = 3L)
  terms <- .s2z_spectral_group_terms(case)
  Q <- 4L
  H <- matrix(
    sin(seq_len(Q * case$M) * 0.23) +
      cos(seq_len(Q * case$M) * 0.47),
    Q, case$M
  )
  V_seed <- matrix(
    sin(seq_len(Q * Q) * 0.17) - cos(seq_len(Q * Q) * 0.31),
    Q, Q
  )
  V_beta <- tcrossprod(V_seed) / Q + diag(seq(0.9, 1.5, length.out = Q))
  Lambda_beta <- solve(V_beta)
  theta <- seq(-0.6, 0.75, length.out = Q)
  prior_mean <- seq(0.25, -0.35, length.out = Q)
  difference <- theta - prior_mean

  # x is the coefficient-whitened omitted mean: m = L x.
  A <- H %*% case$L
  P_prior <- crossprod(A, Lambda_beta %*% A)
  h_prior <- drop(crossprod(A, Lambda_beta %*% difference))
  P <- terms$P_modal + P_prior
  h <- terms$h_modal + h_prior
  base <- -0.5 * (
    terms$q_modal + drop(crossprod(difference, Lambda_beta %*% difference))
  )

  direct_log_kernel <- function(x) {
    beta_difference <- difference - drop(A %*% x)
    shifted <- terms$delta_white + outer(rep(1, case$G), x)
    -0.5 * (
      drop(crossprod(beta_difference, Lambda_beta %*% beta_difference)) +
      sum(shifted * solve(case$Omega, shifted))
    )
  }
  completed_log_kernel <- function(x) {
    base + drop(crossprod(h, x)) -
      0.5 * drop(crossprod(x, P %*% x))
  }

  points <- rbind(
    rep(0, case$M),
    seq(-0.3, 0.4, length.out = case$M),
    seq(0.5, -0.2, length.out = case$M)
  )
  for (i in seq_len(nrow(points))) {
    expect_equal(
      completed_log_kernel(points[i, ]), direct_log_kernel(points[i, ]),
      tolerance = 3e-12
    )
  }
  expect_gt(sum(abs(P_prior)), 1e-3)
  expect_gt(sum(abs(h_prior)), 1e-3)
  expect_true(all(eigen(P, symmetric = TRUE)$values > 0))

  conditional_mean <- drop(solve(P, h))
  epsilon <- 1e-6
  gradient <- vapply(seq_len(case$M), function(j) {
    step <- numeric(case$M)
    step[j] <- epsilon
    (direct_log_kernel(conditional_mean + step) -
      direct_log_kernel(conditional_mean - step)) / (2 * epsilon)
  }, numeric(1))
  expect_equal(gradient, numeric(case$M), tolerance = 2e-8)

  integrated_log_kernel <- base +
    0.5 * drop(crossprod(h, solve(P, h))) +
    0.5 * case$M * log(2 * pi) -
    0.5 * as.numeric(determinant(P, logarithm = TRUE)$modulus)
  log_at_mode <- direct_log_kernel(conditional_mean)
  expect_equal(
    integrated_log_kernel,
    log_at_mode + 0.5 * case$M * log(2 * pi) -
      0.5 * as.numeric(determinant(P, logarithm = TRUE)$modulus),
    tolerance = 3e-12
  )

  physical_mean <- drop(case$L %*% conditional_mean)
  beta <- theta - drop(H %*% physical_mean)
  group_effects <- case$delta + outer(rep(1, case$G), physical_mean)
  expect_equal(
    beta + drop(H %*% colMeans(group_effects)), theta,
    tolerance = 2e-13
  )
})
