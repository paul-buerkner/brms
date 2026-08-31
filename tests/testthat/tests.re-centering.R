context("Tests for centered ordinary group-level effects")

expect_match2 <- brms:::expect_match2

.re_center_map <- function(q, L, rho) {
  stopifnot(length(q) == nrow(L), length(rho) == length(q))
  A <- diag(1 - rho, nrow(L)) + diag(rho, nrow(L)) %*% L
  white <- forwardsolve(A, q)
  effect <- drop(L %*% white)
  log_jacobian <- sum(log(diag(L))) - sum(log(diag(A)))
  list(A = A, white = white, effect = effect,
       log_jacobian = log_jacobian)
}

.re_center_log_density <- function(q, L, rho) {
  mapped <- .re_center_map(q, L, rho)
  physical <- sum(dnorm(forwardsolve(L, mapped$effect), log = TRUE)) -
    sum(log(diag(L)))
  chart <- sum(dnorm(mapped$white, log = TRUE)) -
    sum(log(diag(mapped$A)))
  c(physical_plus_jacobian = physical + mapped$log_jacobian,
    chart = chart)
}

.re_center_numeric_jacobian <- function(fun, x, eps = 1e-6) {
  out <- matrix(NA_real_, length(fun(x)), length(x))
  for (j in seq_along(x)) {
    step <- rep(0, length(x))
    step[j] <- eps
    out[, j] <- (fun(x + step) - fun(x - step)) / (2 * eps)
  }
  out
}

.re_center_stan_between <- function(x, start, end) {
  start_at <- regexpr(start, x, fixed = TRUE)[1L]
  stopifnot(start_at > 0L)
  from <- start_at + nchar(start)
  end_at <- regexpr(end, substring(x, from), fixed = TRUE)[1L]
  stopifnot(end_at > 0L)
  substring(x, from, from + end_at - 2L)
}

test_that("unrestricted centering map has exact density and endpoints", {
  L <- matrix(c(
    1.7, -0.65, 0.45,
    0, 0.8, 0.35,
    0, 0, 1.25
  ), 3, 3)
  q <- c(-0.7, 1.1, 0.25)

  for (rho in list(rep(0, 3), c(0.15, 0.63, 0.9), rep(1, 3))) {
    mapped <- .re_center_map(q, L, rho)
    density <- .re_center_log_density(q, L, rho)
    expect_equal(density[[1]], density[[2]], tolerance = 1e-12)
    dense_map <- L %*% solve(mapped$A)
    expect_equal(
      mapped$log_jacobian,
      as.numeric(determinant(dense_map, logarithm = TRUE)$modulus),
      tolerance = 1e-12
    )
  }

  expect_equal(.re_center_map(q, L, rep(0, 3))$effect, drop(L %*% q))
  expect_equal(.re_center_map(q, L, rep(1, 3))$effect, q)
})

test_that("conditional Student scales remain inside the exact chart", {
  L_reference <- matrix(c(1.4, 0.55, 0, 0.75), 2, 2)
  q <- list(c(-0.4, 0.8), c(1.2, -0.3), c(0.1, 0.6))
  mixers <- c(0.45, 1.3, 2.2)
  rho <- c(0.27, 0.81)

  densities <- vapply(seq_along(q), function(j) {
    L_j <- mixers[j] * L_reference
    values <- .re_center_log_density(q[[j]], L_j, rho)
    expect_equal(values[[1]], values[[2]], tolerance = 1e-12)
    values[[1]]
  }, numeric(1))
  expect_true(all(is.finite(densities)))

  for (j in seq_along(q)) {
    L_j <- mixers[j] * L_reference
    expect_equal(
      .re_center_map(q[[j]], L_j, c(0, 0))$effect,
      drop(L_j %*% q[[j]])
    )
    expect_equal(
      .re_center_map(q[[j]], L_j, c(1, 1))$effect,
      q[[j]]
    )
  }
})

test_that("parameter-dependent Fisher charts have block-triangular Jacobian", {
  transform <- function(x) {
    q <- x[1:2]
    theta <- x[3:4]
    L <- matrix(c(
      exp(theta[1]), 0.6 * tanh(theta[2]),
      0, exp(-0.25 * theta[1])
    ), 2, 2)
    rho <- plogis(c(theta[1] - theta[2], theta[2] + 0.2))
    c(.re_center_map(q, L, rho)$effect, theta)
  }
  x <- c(-0.35, 0.9, 0.25, -0.45)
  jacobian <- .re_center_numeric_jacobian(transform, x)
  numeric_log_jacobian <- log(abs(det(jacobian)))

  theta <- x[3:4]
  L <- matrix(c(
    exp(theta[1]), 0.6 * tanh(theta[2]),
    0, exp(-0.25 * theta[1])
  ), 2, 2)
  rho <- plogis(c(theta[1] - theta[2], theta[2] + 0.2))
  analytic <- .re_center_map(x[1:2], L, rho)$log_jacobian
  expect_equal(numeric_log_jacobian, analytic, tolerance = 2e-8)
  expect_equal(jacobian[3:4, 1:2], matrix(0, 2, 2), tolerance = 1e-10)
})

re_center_dat <- local({
  set.seed(1923)
  data.frame(
    y = rnorm(48),
    y2 = rnorm(48),
    x = seq(-1, 1, length.out = 48),
    g = factor(rep(seq_len(6), each = 8)),
    h = factor(rep(seq_len(8), each = 6)),
    f = factor(rep(c("a", "b"), 24)),
    f_nested = factor(rep(c("a", "b"), each = 24)),
    w = rep(seq(0.8, 1.2, length.out = 6), each = 8)
  )
})

test_that("ordinary correlated charts generate exact unrestricted maps", {
  centered <- stancode(
    y ~ x + (1 + x | gr(g, center = TRUE)), data = re_center_dat
  )
  partial <- stancode(
    y ~ x +
      (1 | gr(g, id = "mixed", center = 0.2)) +
      (0 + x | gr(g, id = "mixed", center = 0.8)),
    data = re_center_dat
  )
  fisher <- stancode(
    y ~ x + (1 + x | gr(g, center = "fisher")),
    data = re_center_dat
  )

  expect_match2(centered, "matrix[M_1, N_1] z_1;")
  expect_match2(centered, "r_1 = z_1';")
  expect_match2(centered, "vector[N_1] r_1_1;")
  expect_match2(centered, "vector[N_1] r_1_2;")
  expect_match2(centered, "corr_matrix[M_1] Cor_1")
  expect_match2(centered, "vector<lower=-1,upper=1>[NC_1] cor_1;")
  expect_match2(centered, "multi_normal_cholesky_lpdf(r_1[j]'")
  expect_match2(centered, "target += log_jacobian_re_1;")

  expect_match2(partial, "vector[M_1] rho_center_re")
  expect_match2(
    partial,
    "diag_pre_multiply(rho_center_re, L_group_center_re)"
  )
  expect_match2(
    partial,
    "white_center_re = mdivide_left_tri_low(L_partial_center_re, z_1[, j]);"
  )
  expect_match2(
    partial,
    "sum(log(diagonal(L_group_center_re))) - sum(log(diagonal(L_partial_center_re)))"
  )
  expect_false(grepl("mean_partial_s2z", partial, fixed = TRUE))

  expect_match2(
    fisher, "matrix<lower=0,upper=1>[N_1, M_1] rho_s2z_1;"
  )
  expect_match2(fisher, "rho_level_center_re = rho_s2z_1[j]';")
  expect_null(standata(
    y ~ x + (1 + x | gr(g, center = 0.35)), data = re_center_dat
  )$rho_s2z_1)
  expect_null(standata(
    y ~ x + (1 + x | gr(g, center = "fisher")), data = re_center_dat
  )$rho_s2z_1)
  tpar <- .re_center_stan_between(
    fisher, "transformed parameters {", "\nmodel {"
  )
  expect_false(grepl("Y[n]", tpar, fixed = TRUE))
})

test_that("ordinary independent and Student charts retain conditional scales", {
  independent <- stancode(
    y ~ x + (1 + x || gr(g, center = 0.35)), data = re_center_dat
  )
  student <- stancode(
    y ~ x + (1 + x | gr(g, dist = "student", center = 0.35)),
    data = re_center_dat
  )
  student_independent <- stancode(
    y ~ x + (1 + x || gr(g, dist = "student", center = 0.35)),
    data = re_center_dat
  )
  student_fisher <- stancode(
    y ~ x + (1 + x | gr(
      g, dist = "student", center = "fisher"
    )),
    data = re_center_dat
  )

  expect_match2(independent, "vector[N_1] denominator_center_re")
  expect_match2(independent, "r_1_1 = scale_group_center_re .* z_1[1]")
  expect_match2(independent, "target += normal_lpdf(r_1_1 | 0, sd_1[1]);")
  expect_match2(student, "L_group_center_re *= dfm_1[j];")
  expect_match2(
    student, "multi_normal_cholesky_lpdf(r_1[j]' | zeros_vector(M_1)"
  )
  expect_match2(
    student_independent, "scale_group_center_re = sd_1[1] * dfm_1;"
  )
  expect_match2(
    student_independent, "normal_lpdf(r_1_1 | 0, sd_1[1] * dfm_1)"
  )
  expect_match2(
    student_fisher,
    "vector[M_1] prior_var_fisher_s2z = rows_dot_self(L_center_re_1);"
  )
  expect_match2(student_fisher, "L_group_center_re *= dfm_1[j];")
  fisher_comp <- .re_center_stan_between(
    student_fisher, "L_center_re_1 =", "L_group_center_re *= dfm_1[j];"
  )
  expect_false(grepl("dfm_1", fisher_comp, fixed = TRUE))
})

test_that("fixed ordinary centering composes with supported predictor types", {
  ordinal_dat <- transform(
    re_center_dat,
    yo = ordered(rep(c("low", "middle", "high"), length.out = 48))
  )
  ordinal <- stancode(
    yo ~ x + (1 + x | gr(g, center = 0.4)),
    data = ordinal_dat, family = cumulative()
  )
  expect_match2(ordinal, "ordered[nthres] Intercept")
  expect_match2(ordinal, "L_partial_center_re")

  distributional <- stancode(
    bf(
      y ~ x + (1 | gr(g, id = "mu-local", center = 0.3)),
      sigma ~ x + (1 | gr(g, id = "sigma-local", center = 0.7))
    ),
    data = re_center_dat, family = gaussian()
  )
  expect_match2(distributional, "log_jacobian_re_1")
  expect_match2(distributional, "log_jacobian_re_2")

  multivariate <- stancode(
    bf(y ~ x + (1 | gr(g, id = "y-local", center = 0.25))) +
      bf(y2 ~ x + (1 | gr(g, id = "y2-local", center = 0.75))) +
      set_rescor(FALSE),
    data = re_center_dat
  )
  expect_match2(multivariate, "log_jacobian_re_1")
  expect_match2(multivariate, "log_jacobian_re_2")

  multivariate_fisher <- stancode(
    bf(y ~ x + (1 | gr(g, id = "y-fisher", center = "fisher"))) +
      bf(y2 ~ x + (1 | gr(g, id = "y2-fisher", center = "fisher"))) +
      set_rescor(FALSE),
    data = re_center_dat
  )
  expect_match2(multivariate_fisher, "rho_s2z_1[j, 1]")
  expect_match2(multivariate_fisher, "rho_s2z_2[j, 1]")
  expect_error(
    stancode(
      bf(y ~ x + (1 | gr(g, center = "fisher"))) +
        bf(y2 ~ x) + set_rescor(TRUE),
      data = re_center_dat
    ),
    "requires set_rescor\\(FALSE\\)"
  )

  multiple <- stancode(
    y ~ x + (1 | gr(g, center = 0.25)) +
      (1 | gr(h, center = 0.75)),
    data = re_center_dat
  )
  expect_match2(multiple, "log_jacobian_re_1")
  expect_match2(multiple, "log_jacobian_re_2")

  nonlinear_form <- bf(
    y ~ eta,
    eta ~ 1 + x + (1 | gr(g, center = 0.4)),
    nl = TRUE
  )
  nonlinear_prior <- prior(normal(0, 1), class = "b", nlpar = "eta")
  nonlinear <- stancode(
    nonlinear_form, data = re_center_dat, prior = nonlinear_prior
  )
  expect_match2(nonlinear, "log_jacobian_re_1")
  expect_error(
    stancode(
      bf(
        y ~ eta,
        eta ~ 1 + x + (1 | gr(g, center = "fisher")),
        nl = TRUE
      ),
      data = re_center_dat, prior = nonlinear_prior
    ),
    "not yet supported.*nonlinear predictors"
  )
})

test_that("ordinary centering gates later structural combinations", {
  expect_error(
    stancode(
      y ~ x + (1 + x | gr(g, by = f_nested, center = 0.4)),
      data = re_center_dat
    ),
    "'by', 'cov', or 'pw'"
  )
  expect_error(
    stancode(
      y ~ x + (1 | gr(g, pw = w, center = 0.4)),
      data = re_center_dat
    ),
    "'by', 'cov', or 'pw'"
  )
  expect_error(
    stancode(
      y ~ x + (1 | gr(g, cov = A, center = 0.4)),
      data = re_center_dat,
      data2 = list(A = structure(
        diag(nlevels(re_center_dat$g)),
        dimnames = rep(list(levels(re_center_dat$g)), 2)
      ))
    ),
    "'by', 'cov', or 'pw'"
  )

  cross <- bf(
    y ~ x + (1 | gr(g, id = "cross", center = 0.4)),
    sigma ~ x + (1 | gr(g, id = "cross", center = 0.4))
  )
  expect_error(
    stancode(cross, data = re_center_dat),
    "must be predictor-local"
  )

  ordinal_dat <- transform(
    re_center_dat,
    yo = ordered(rep(c("low", "middle", "high"), length.out = 48))
  )
  expect_error(
    stancode(
      yo ~ x + (1 | gr(g, center = "fisher")),
      data = ordinal_dat, family = cumulative()
    ),
    "not yet supported.*ordinal location"
  )
  expect_error(
    stancode(
      y ~ x +
        (1 | gr(g, id = "mixed", center = "fisher")) +
        (0 + x | gr(g, id = "mixed", center = 0.5)),
      data = re_center_dat
    ),
    "must use Fisher centering if any coefficient does"
  )
  expect_error(
    stancode(
      y ~ x + (1 | gr(g, center = TRUE)),
      data = re_center_dat,
      prior = prior(constant(0), class = "sd", group = "g",
                    coef = "Intercept")
    ),
    "must be positive numeric scalars"
  )
  expect_silent(stancode(
    y ~ x + (1 | gr(g, center = TRUE)),
    data = re_center_dat,
    prior = prior(constant(1.25), class = "sd", group = "g",
                  coef = "Intercept")
  ))
})

test_that("ordinary centering state is excluded unless all state is saved", {
  fit <- brm(
    y ~ x + (1 + x | gr(g, center = "fisher")),
    data = re_center_dat, empty = TRUE, cores = 1
  )
  excluded <- brms:::exclude_pars(fit)
  expect_true(all(c(
    "z_1", "L_center_re_1", "log_jacobian_re_1",
    "rho_s2z_1", "mean_rho_s2z_1"
  ) %in% excluded))

  fit_all <- brm(
    y ~ x + (1 + x | gr(g, center = "fisher")),
    data = re_center_dat, empty = TRUE, cores = 1,
    save_pars = save_pars(all = TRUE)
  )
  saved_exclusions <- brms:::exclude_pars(fit_all)
  expect_false(any(c(
    "z_1", "L_center_re_1", "log_jacobian_re_1",
    "rho_s2z_1", "mean_rho_s2z_1"
  ) %in% saved_exclusions))
})
