context("Reference density checks for ordinal physical S2Z effects")

# These helpers intentionally use only base R and stats.  In particular, the
# expected values below are not obtained from the production S2Z generators or
# their algebra helpers.
.ordinal_ref_block_diag <- function(...) {
  blocks <- list(...)
  nrows <- vapply(blocks, nrow, integer(1))
  ncols <- vapply(blocks, ncol, integer(1))
  out <- matrix(0, sum(nrows), sum(ncols))
  row_start <- 1L
  col_start <- 1L
  for (i in seq_along(blocks)) {
    rows <- row_start - 1L + seq_len(nrows[i])
    cols <- col_start - 1L + seq_len(ncols[i])
    out[rows, cols] <- blocks[[i]]
    row_start <- row_start + nrows[i]
    col_start <- col_start + ncols[i]
  }
  out
}

.ordinal_ref_log_inv_chisq <- function(x, df) {
  ifelse(
    x > 0,
    -0.5 * df * log(2) - lgamma(0.5 * df) -
      (0.5 * df + 1) * log(x) - 0.5 / x,
    -Inf
  )
}

test_that("negative threshold rows preserve the ordinal likelihood", {
  x_raw <- c(-1.6, -0.4, 0.2, 0.9, 1.7, 2.4)
  x_mean <- mean(x_raw)
  X <- matrix(x_raw - x_mean, ncol = 1L)
  Z <- cbind(1, x_raw)
  a <- c(1, x_mean)
  C <- matrix(c(0, 1), nrow = 1L)

  # Check the affine design identity rather than assuming a particular
  # intercept-plus-slope layout.
  affine_Z <- matrix(1, nrow(X), 1L) %*% t(a) + X %*% C
  expect_equal(unname(affine_Z), unname(Z), tolerance = 2e-15)

  group <- c(1L, 2L, 3L, 1L, 2L, 3L)
  delta <- matrix(
    c(-0.55, 0.35, 0.20, 0.45, -0.70, 0.25),
    nrow = 3L
  )
  expect_equal(colSums(delta), c(0, 0), tolerance = 1e-15)
  omitted_mean <- c(0.42, -0.31)
  conventional_r <- sweep(delta, 2L, omitted_mean, "+")
  conventional_beta <- 0.73
  conventional_threshold <- c(-1.25, 0.15, 1.60)

  finite_beta <- drop(conventional_beta + C %*% omitted_mean)
  common_shift <- sum(a * omitted_mean)
  finite_threshold <- conventional_threshold - common_shift
  eta_conventional <- drop(X %*% conventional_beta) +
    rowSums(Z * conventional_r[group, , drop = FALSE])
  eta_finite <- drop(X %*% finite_beta) +
    rowSums(Z * delta[group, , drop = FALSE])

  conventional_difference <-
    outer(rep(1, nrow(X)), conventional_threshold) - eta_conventional
  finite_difference <-
    outer(rep(1, nrow(X)), finite_threshold) - eta_finite
  expect_equal(
    conventional_difference, finite_difference, tolerance = 2e-15
  )

  # Observation-specific discrimination is harmless because it multiplies
  # the already preserved threshold-minus-predictor difference.
  discrimination <- seq(0.55, 1.45, length.out = nrow(X))
  conventional_cdf <- plogis(sweep(
    conventional_difference, 1L, discrimination, "*"
  ))
  finite_cdf <- plogis(sweep(
    finite_difference, 1L, discrimination, "*"
  ))
  expect_equal(conventional_cdf, finite_cdf, tolerance = 2e-15)

  # q = theta - H m has one negative row per threshold and a positive
  # matching-slope row.  Reversing only the threshold sign breaks the
  # per-draw likelihood identity.
  H <- rbind(
    matrix(
      rep(-a, times = length(conventional_threshold)),
      nrow = length(conventional_threshold), byrow = TRUE
    ),
    C
  )
  theta <- c(finite_threshold, finite_beta)
  recovered_q <- theta - drop(H %*% omitted_mean)
  expect_equal(
    recovered_q,
    c(conventional_threshold, conventional_beta),
    tolerance = 2e-15
  )

  H_wrong <- H
  H_wrong[seq_along(conventional_threshold), ] <- -
    H_wrong[seq_along(conventional_threshold), ]
  wrong_q <- theta - drop(H_wrong %*% omitted_mean)
  expect_gt(
    max(abs(wrong_q[seq_along(conventional_threshold)] -
      conventional_threshold)),
    0.1
  )
})

test_that("equidistant recovery translates only the first coordinate", {
  x_raw <- c(-1.8, -0.9, -0.1, 0.4, 0.8, 1.3, 1.9, 2.6)
  x_mean <- mean(x_raw)
  X <- matrix(x_raw - x_mean, ncol = 1L)
  Z <- cbind(1, x_raw)
  a <- c(1, x_mean)
  C <- matrix(c(0, 1), nrow = 1L)
  expect_equal(
    unname(matrix(1, nrow(X), 1L) %*% t(a) + X %*% C),
    unname(Z), tolerance = 2e-15
  )

  group <- rep(seq_len(4L), 2L)
  physical_r <- rbind(
    c(-0.50, 0.20), c(0.30, -0.60),
    c(0.40, 0.10), c(-0.20, 0.30)
  )
  expect_equal(colSums(physical_r), c(0, 0), tolerance = 1e-15)
  omitted_mean <- c(0.37, -0.24)
  conventional_r <- sweep(physical_r, 2L, omitted_mean, "+")
  conventional_beta <- 0.68
  conventional_first <- -1.42
  spacing <- 0.83
  threshold_index <- 0:3
  conventional_threshold <- conventional_first + threshold_index * spacing

  finite_beta <- drop(conventional_beta + C %*% omitted_mean)
  finite_first <- conventional_first - sum(a * omitted_mean)
  finite_threshold <- finite_first + threshold_index * spacing
  eta_conventional <- drop(X %*% conventional_beta) +
    rowSums(Z * conventional_r[group, , drop = FALSE])
  eta_finite <- drop(X %*% finite_beta) +
    rowSums(Z * physical_r[group, , drop = FALSE])

  conventional_difference <-
    outer(rep(1, nrow(X)), conventional_threshold) - eta_conventional
  finite_difference <-
    outer(rep(1, nrow(X)), finite_threshold) - eta_finite
  expect_equal(
    conventional_difference, finite_difference, tolerance = 2e-15
  )
  conventional_probability <- t(vapply(
    eta_conventional,
    function(eta) diff(c(0, plogis(conventional_threshold - eta), 1)),
    numeric(length(conventional_threshold) + 1L)
  ))
  finite_probability <- t(vapply(
    eta_finite,
    function(eta) diff(c(0, plogis(finite_threshold - eta), 1)),
    numeric(length(finite_threshold) + 1L)
  ))
  expect_equal(
    conventional_probability, finite_probability, tolerance = 2e-15
  )

  # The equidistant q vector contains first_Intercept and slopes only. The
  # spacing is unchanged and is deliberately absent from H.
  H <- rbind(-a, C)
  theta <- c(finite_first, finite_beta)
  recovered_q <- theta - drop(H %*% omitted_mean)
  expect_equal(
    recovered_q, c(conventional_first, conventional_beta),
    tolerance = 2e-15
  )
  recovered_threshold <- recovered_q[1L] + threshold_index * spacing
  expect_equal(
    recovered_threshold, conventional_threshold, tolerance = 2e-15
  )
  expect_equal(diff(finite_threshold), rep(spacing, 3L))
  expect_equal(nrow(H), 2L)
})

test_that("grouped thresholds share one shift and recover in flattened order", {
  x_raw <- c(-1.5, -0.8, -0.2, 0.1, 0.5, 0.9, 1.4, 1.8, 2.2, 2.7)
  x_mean <- mean(x_raw)
  X <- matrix(x_raw - x_mean, ncol = 1L)
  Z <- cbind(1, x_raw)
  a <- c(1, x_mean)
  C <- matrix(c(0, 1), nrow = 1L)
  expect_equal(
    unname(matrix(1, nrow(X), 1L) %*% t(a) + X %*% C),
    unname(Z), tolerance = 2e-15
  )

  group <- rep(seq_len(5L), 2L)
  threshold_group <- rep(c("short", "long"), each = 5L)
  physical_r <- rbind(
    c(-0.40, 0.25), c(0.15, -0.50), c(0.30, 0.10),
    c(0.20, 0.35), c(-0.25, -0.20)
  )
  expect_equal(colSums(physical_r), c(0, 0), tolerance = 1e-15)
  omitted_mean <- c(0.31, -0.27)
  conventional_r <- sweep(physical_r, 2L, omitted_mean, "+")
  conventional_beta <- -0.56
  conventional_threshold <- list(
    short = c(-0.92, 0.44),
    long = c(-1.74, -0.58, 0.31, 1.49)
  )
  common_shift <- sum(a * omitted_mean)
  finite_threshold <- lapply(
    conventional_threshold, function(x) x - common_shift
  )
  finite_beta <- drop(conventional_beta + C %*% omitted_mean)

  eta_conventional <- drop(X %*% conventional_beta) +
    rowSums(Z * conventional_r[group, , drop = FALSE])
  eta_finite <- drop(X %*% finite_beta) +
    rowSums(Z * physical_r[group, , drop = FALSE])
  for (label in names(conventional_threshold)) {
    take <- threshold_group == label
    conventional_difference <- outer(
      rep(1, sum(take)), conventional_threshold[[label]]
    ) - eta_conventional[take]
    finite_difference <- outer(
      rep(1, sum(take)), finite_threshold[[label]]
    ) - eta_finite[take]
    expect_equal(
      conventional_difference, finite_difference,
      tolerance = 2e-15, info = label
    )

    conventional_probability <- t(vapply(
      eta_conventional[take],
      function(eta) diff(c(
        0, plogis(conventional_threshold[[label]] - eta), 1
      )),
      numeric(length(conventional_threshold[[label]]) + 1L)
    ))
    finite_probability <- t(vapply(
      eta_finite[take],
      function(eta) diff(c(
        0, plogis(finite_threshold[[label]] - eta), 1
      )),
      numeric(length(finite_threshold[[label]]) + 1L)
    ))
    expect_equal(
      conventional_probability, finite_probability,
      tolerance = 2e-15, info = label
    )
  }

  counts <- lengths(conventional_threshold)
  n_threshold <- sum(counts)
  H <- rbind(
    matrix(rep(-a, times = n_threshold), nrow = n_threshold, byrow = TRUE),
    C
  )
  theta <- c(
    unlist(finite_threshold, use.names = FALSE), finite_beta
  )
  recovered_q <- theta - drop(H %*% omitted_mean)
  expect_equal(
    recovered_q,
    c(unlist(conventional_threshold, use.names = FALSE), conventional_beta),
    tolerance = 2e-15
  )

  starts <- cumsum(c(1L, head(counts, -1L)))
  recovered_threshold <- Map(
    function(start, n) recovered_q[start - 1L + seq_len(n)],
    starts, counts
  )
  names(recovered_threshold) <- names(conventional_threshold)
  expect_equal(recovered_threshold, conventional_threshold, tolerance = 2e-15)
  expect_equal(
    unlist(conventional_threshold, use.names = FALSE) -
      unlist(finite_threshold, use.names = FALSE),
    rep(common_shift, n_threshold), tolerance = 2e-15
  )
})

test_that("the Gaussian omitted-mean density has every normalizer", {
  groups <- 4L
  L <- 0.82
  H <- matrix(c(-1.15, -1.15, -1.15, 0.70), ncol = 1L)
  theta <- c(-1.30, 0.20, 1.55, 0.65)
  prior_mean <- c(-1.05, -0.15, 1.10, -0.35)
  prior_sd <- c(0.58, 1.05, 0.76, 0.91)
  prior_precision <- 1 / prior_sd^2
  delta <- c(-0.72, 0.13, 0.44, 0.15)
  expect_equal(sum(delta), 0, tolerance = 1e-15)

  d <- sqrt(prior_precision) * (theta - prior_mean)
  F <- sweep(H * L, 1L, sqrt(prior_precision), "*")
  P <- groups + drop(crossprod(F))
  h <- drop(crossprod(F, d))
  log_formula <- sum(stats::dnorm(
    theta, prior_mean, prior_sd, log = TRUE
  )) - 0.5 * sum((delta / L)^2) + 0.5 * h^2 / P -
    0.5 * log(P) - (groups - 1L) * log(L) -
    0.5 * (groups - 1L) * log(2 * pi) + 0.5 * log(groups)

  # Numerically integrate the original conventional density over its omitted
  # arithmetic mean.  sqrt(G) is the Jacobian from group effects to an
  # orthonormal contrast plus the arithmetic mean.
  log_joint_measure <- function(mean_effect) {
    conventional_q <- theta - drop(H * mean_effect)
    sum(stats::dnorm(
      conventional_q, prior_mean, prior_sd, log = TRUE
    )) + sum(stats::dnorm(
      delta + mean_effect, 0, L, log = TRUE
    )) + 0.5 * log(groups)
  }
  optimum <- optimize(
    log_joint_measure, interval = c(-8, 8), maximum = TRUE,
    tol = 1e-13
  )
  peak <- optimum$objective
  integral <- stats::integrate(
    function(mean_effect) {
      vapply(
        mean_effect,
        function(value) exp(log_joint_measure(value) - peak),
        numeric(1)
      )
    },
    lower = -Inf, upper = Inf, rel.tol = 2e-12,
    subdivisions = 500L, stop.on.error = TRUE
  )$value
  log_numerical <- peak + log(integral)

  expect_equal(log_formula, log_numerical, tolerance = 2e-11)

  # Omitting the contrast/common-mean measure factor changes the normalized
  # log density by exactly log(sqrt(G)).
  expect_equal(
    log_formula - (log_numerical - 0.5 * log(groups)),
    0.5 * log(groups),
    tolerance = 2e-11
  )
})

test_that("conditional recovery has the direct Gaussian moments", {
  groups <- 5L
  dimensions <- 2L
  L <- matrix(c(0.78, 0, -0.21, 1.08), dimensions, byrow = TRUE)
  Sigma <- tcrossprod(L)
  a <- c(1, 0.37)
  C <- diag(dimensions)
  H <- rbind(
    matrix(rep(-a, times = 3L), nrow = 3L, byrow = TRUE),
    C
  )
  theta <- c(-1.35, -0.10, 1.25, 0.62, -0.48)
  prior_mean <- c(-1.10, 0.25, 0.95, -0.20, 0.35)
  prior_precision <- c(1.8, 0.55, 2.4, 0.82, 1.35)
  delta <- rbind(
    c(-0.42, 0.25), c(0.31, -0.63), c(0.54, 0.18),
    c(-0.16, 0.47), c(-0.27, -0.27)
  )
  expect_equal(colSums(delta), numeric(dimensions), tolerance = 1e-15)

  d <- sqrt(prior_precision) * (theta - prior_mean)
  F <- sweep(H %*% L, 1L, sqrt(prior_precision), "*")
  P <- groups * diag(dimensions) + crossprod(F)
  h <- drop(crossprod(F, d))
  mean_white <- drop(solve(P, h))
  covariance_white <- solve(P)
  recovery_mean <- drop(L %*% mean_white)
  recovery_covariance <- L %*% covariance_white %*% t(L)

  # Expand the conventional conditional density directly in m.  This uses
  # q = theta - H m and r_g = delta_g + m, without whitening by L.
  raw_precision <- crossprod(
    H, sweep(H, 1L, prior_precision, "*")
  ) + groups * solve(Sigma)
  raw_score <- drop(crossprod(
    H, prior_precision * (theta - prior_mean)
  )) - drop(solve(Sigma, colSums(delta)))
  direct_mean <- drop(solve(raw_precision, raw_score))
  direct_covariance <- solve(raw_precision)

  expect_equal(recovery_mean, direct_mean, tolerance = 2e-13)
  expect_equal(
    recovery_covariance, direct_covariance, tolerance = 2e-13
  )

  conventional_q <- theta - drop(H %*% recovery_mean)
  conventional_r <- sweep(delta, 2L, recovery_mean, "+")
  direct_gradient <- drop(crossprod(
    H, prior_precision * (conventional_q - prior_mean)
  )) - drop(solve(Sigma, colSums(conventional_r)))
  expect_equal(direct_gradient, numeric(dimensions), tolerance = 3e-13)

  # R's upper Cholesky R has P = R'R, hence R^-1 is the exact white-noise
  # recovery root corresponding to the lower-Cholesky formula in the plan.
  recovery_root <- L %*% solve(chol(P))
  expect_equal(
    tcrossprod(recovery_root), recovery_covariance,
    tolerance = 2e-13
  )
})

test_that("inverse-chi-square threshold mixtures are exactly Student t", {
  cases <- list(
    list(value = -0.35, location = 0.40, scale = 1.25, df = 1),
    list(value = 1.70, location = -0.25, scale = 0.72, df = 2.75),
    list(value = -1.10, location = 0.65, scale = 1.60, df = 8)
  )

  for (case in cases) {
    log_integrand <- function(log_u) {
      u <- exp(log_u)
      stats::dnorm(
        case$value, case$location,
        case$scale * sqrt(case$df * u), log = TRUE
      ) + .ordinal_ref_log_inv_chisq(u, case$df) + log_u
    }
    optimum <- optimize(
      log_integrand, interval = c(-35, 35), maximum = TRUE,
      tol = 1e-13
    )
    peak <- optimum$objective
    mixture <- stats::integrate(
      function(log_u) {
        vapply(
          log_u,
          function(value) exp(log_integrand(value) - peak),
          numeric(1)
        )
      },
      lower = -40, upper = 40, rel.tol = 2e-12,
      subdivisions = 500L, stop.on.error = TRUE
    )$value
    log_mixture <- peak + log(mixture)
    log_student <- stats::dt(
      (case$value - case$location) / case$scale,
      df = case$df, log = TRUE
    ) - log(case$scale)

    expect_equal(log_mixture, log_student, tolerance = 2e-10)
  }
})

test_that("Student group recovery retains the weighted contrast score", {
  groups <- 5L
  dimensions <- 2L
  L <- matrix(c(0.66, 0, 0.24, 1.12), dimensions, byrow = TRUE)
  Sigma <- tcrossprod(L)
  H <- matrix(
    c(-1, -0.35, -1, -0.35, -1, -0.35, 0, 1, 1, 0.20),
    ncol = dimensions, byrow = TRUE
  )
  theta <- c(-1.15, 0.05, 1.30, 0.58, -0.44)
  prior_mean <- c(-0.80, -0.20, 0.95, 0.15, 0.31)
  prior_precision <- c(1.6, 0.48, 2.1, 0.75, 1.35)
  delta <- rbind(
    c(-0.48, 0.30), c(0.37, -0.62), c(0.55, 0.11),
    c(-0.19, 0.46), c(-0.25, -0.25)
  )
  expect_equal(colSums(delta), numeric(dimensions), tolerance = 1e-15)

  # Conditional inverse-chi-square draws produce unequal level precisions.
  df <- 4.5
  inverse_chisq <- c(0.10, 0.28, 0.77, 1.45, 0.21)
  level_scale <- sqrt(df * inverse_chisq)
  weights <- 1 / level_scale^2
  white_delta <- t(forwardsolve(L, t(delta)))
  weighted_contrast <- colSums(sweep(
    white_delta, 1L, weights, "*"
  ))
  expect_gt(max(abs(weighted_contrast)), 0.1)

  d <- sqrt(prior_precision) * (theta - prior_mean)
  F <- sweep(H %*% L, 1L, sqrt(prior_precision), "*")
  P <- crossprod(F) + sum(weights) * diag(dimensions)
  h <- drop(crossprod(F, d)) - weighted_contrast
  mean_white <- drop(solve(P, h))
  recovery_mean <- drop(L %*% mean_white)

  # The same mode derived directly from the heteroscedastic conventional
  # group densities.  The arithmetic contrast sums to zero, but its weighted
  # sum generally does not.
  raw_group_score <- drop(solve(
    Sigma, colSums(sweep(delta, 1L, weights, "*"))
  ))
  raw_precision <- crossprod(
    H, sweep(H, 1L, prior_precision, "*")
  ) + sum(weights) * solve(Sigma)
  raw_score <- drop(crossprod(
    H, prior_precision * (theta - prior_mean)
  )) - raw_group_score
  direct_mean <- drop(solve(raw_precision, raw_score))
  expect_equal(recovery_mean, direct_mean, tolerance = 3e-13)

  # Omitting the contrast term leaves a visibly nonzero score in the original
  # conditional density, making this assertion sensitive to both its presence
  # and its sign.
  incomplete_mean <- drop(L %*% solve(P, crossprod(F, d)))
  incomplete_q <- theta - drop(H %*% incomplete_mean)
  incomplete_r <- sweep(delta, 2L, incomplete_mean, "+")
  incomplete_gradient <- drop(crossprod(
    H, prior_precision * (incomplete_q - prior_mean)
  )) - drop(solve(
    Sigma,
    colSums(sweep(incomplete_r, 1L, weights, "*"))
  ))
  expect_gt(max(abs(incomplete_gradient)), 0.1)
})

test_that("flexible thresholds retain unequal precisions and locations", {
  threshold <- c(-1.75, -0.30, 0.62, 1.85)
  location <- c(-1.20, 0.18, -0.05, 1.10)
  precision <- c(0.38, 1.75, 0.64, 3.20)
  a <- c(1, -0.46)
  L <- matrix(c(0.74, 0, 0.22, 1.05), 2L, byrow = TRUE)
  H_threshold <- matrix(
    rep(-a, times = length(threshold)),
    nrow = length(threshold), byrow = TRUE
  )

  d_threshold <- sqrt(precision) * (threshold - location)
  F_threshold <- sweep(
    H_threshold %*% L, 1L, sqrt(precision), "*"
  )
  information_uncollapsed <- crossprod(F_threshold)
  score_uncollapsed <- drop(crossprod(F_threshold, d_threshold))

  v_plus <- drop(crossprod(L, a))
  information_collapsed <- sum(precision) * tcrossprod(v_plus)
  score_collapsed <- -v_plus * sum(precision * (threshold - location))
  expect_equal(
    information_uncollapsed, information_collapsed,
    tolerance = 2e-14
  )
  expect_equal(score_uncollapsed, score_collapsed, tolerance = 2e-14)

  log_rows <- sum(stats::dnorm(
    threshold, location, 1 / sqrt(precision), log = TRUE
  ))
  log_rows_manual <- 0.5 * sum(log(precision)) -
    0.5 * length(threshold) * log(2 * pi) -
    0.5 * sum(precision * (threshold - location)^2)
  expect_equal(log_rows, log_rows_manual, tolerance = 2e-14)

  # Add proper matching-slope priors so reversing all threshold rows changes
  # h without changing P and therefore changes the integrated target.
  H_slope <- matrix(c(0, 1, 1, 0.25), nrow = 2L, byrow = TRUE)
  slope <- c(0.72, -0.41)
  slope_location <- c(-0.20, 0.33)
  slope_precision <- c(0.86, 1.45)
  d_slope <- sqrt(slope_precision) * (slope - slope_location)
  F_slope <- sweep(H_slope %*% L, 1L, sqrt(slope_precision), "*")
  information <- 6 * diag(2L) + information_uncollapsed +
    crossprod(F_slope)
  slope_score <- drop(crossprod(F_slope, d_slope))
  correct_score <- score_uncollapsed + slope_score
  wrong_score <- -score_uncollapsed + slope_score
  correct_adjustment <- 0.5 * sum(correct_score * solve(
    information, correct_score
  ))
  wrong_adjustment <- 0.5 * sum(wrong_score * solve(
    information, wrong_score
  ))

  expect_gt(abs(correct_adjustment - wrong_adjustment), 0.01)
  expect_gt(
    max(abs(score_collapsed +
      v_plus * sum(precision * threshold))),
    0.05
  )
  expect_gt(diff(range(precision)), 1)
})

test_that("shared threshold priors create multiple-block cross terms", {
  threshold_count <- 3L
  dimensions <- c(1L, 2L)
  groups <- c(4L, 7L)
  a1 <- 0.85
  a2 <- c(1, -0.32)
  L1 <- matrix(0.68, 1L, 1L)
  L2 <- matrix(c(0.81, 0, -0.17, 1.16), 2L, byrow = TRUE)
  L_joint <- .ordinal_ref_block_diag(L1, L2)

  H_threshold <- matrix(
    rep(c(-a1, -a2), times = threshold_count),
    nrow = threshold_count, byrow = TRUE
  )
  # All other active rows are local to one block.  They cannot create an
  # off-diagonal block in the population-prior information.
  H_local <- rbind(
    c(1, 0, 0),
    c(0, 1, 0),
    c(0, 0, 1)
  )
  H <- rbind(H_threshold, H_local)
  theta <- c(-1.30, 0.10, 1.42, 0.64, -0.35, 0.48)
  prior_mean <- c(-0.95, -0.18, 1.05, -0.22, 0.15, -0.31)
  prior_precision <- c(0.55, 1.80, 2.35, 0.72, 1.25, 0.91)

  F <- sweep(H %*% L_joint, 1L, sqrt(prior_precision), "*")
  group_information <- diag(c(groups[1], rep(groups[2], dimensions[2])))
  P <- crossprod(F) + group_information
  d <- sqrt(prior_precision) * (theta - prior_mean)
  h <- drop(crossprod(F, d))

  v1 <- drop(crossprod(L1, a1))
  v2 <- drop(crossprod(L2, a2))
  expected_cross <- sum(prior_precision[seq_len(threshold_count)]) *
    v1 * v2
  expect_equal(P[1L, 2:3], expected_cross, tolerance = 2e-14)
  expect_gt(max(abs(P[1L, 2:3])), 0.1)

  # A direct raw-scale expansion gives the same joint system.  Its Cholesky
  # transformation makes clear that these are genuine shared-threshold cross
  # terms, rather than group-prior coupling.
  Sigma1 <- tcrossprod(L1)
  Sigma2 <- tcrossprod(L2)
  raw_group_information <- .ordinal_ref_block_diag(
    groups[1] * solve(Sigma1),
    groups[2] * solve(Sigma2)
  )
  raw_precision <- crossprod(
    H, sweep(H, 1L, prior_precision, "*")
  ) + raw_group_information
  raw_score <- drop(crossprod(
    H, prior_precision * (theta - prior_mean)
  ))
  expect_equal(
    P, crossprod(L_joint, raw_precision %*% L_joint),
    tolerance = 3e-13
  )
  expect_equal(h, drop(crossprod(L_joint, raw_score)), tolerance = 3e-13)

  joint_mean <- drop(L_joint %*% solve(P, h))
  P_blockwise <- P
  P_blockwise[1L, 2:3] <- 0
  P_blockwise[2:3, 1L] <- 0
  blockwise_mean <- drop(L_joint %*% solve(P_blockwise, h))
  expect_gt(max(abs(joint_mean - blockwise_mean)), 0.01)
})
