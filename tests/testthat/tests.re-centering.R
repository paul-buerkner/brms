context("Tests for centered ordinary group-level effects")

expect_match2 <- brms:::expect_match2

.re_center_map <- function(q, L, rho, mu = rep(0, length(q))) {
  stopifnot(length(q) == nrow(L), length(rho) == length(q),
            length(mu) == length(q))
  A <- diag(1 - rho, nrow(L)) + diag(rho, nrow(L)) %*% L
  white <- forwardsolve(A, q - rho * mu)
  effect <- drop(L %*% white)
  log_jacobian <- sum(log(diag(L))) - sum(log(diag(A)))
  list(A = A, white = white, effect = effect,
       log_jacobian = log_jacobian)
}

.re_center_log_density <- function(q, L, rho, mu = rep(0, length(q))) {
  mapped <- .re_center_map(q, L, rho, mu = mu)
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

test_that("scalar chart is the tutorial construction with rho = 1 - w", {
  mu <- -0.45
  sigma <- 1.7
  chi <- 0.8
  w <- 0.37
  rho <- 1 - w
  denominator <- sigma * (1 - w) + w
  tutorial_alpha <- (mu * w + sigma * chi) / denominator
  mapped_u <- .re_center_map(
    chi, matrix(sigma, 1, 1), rho = rho, mu = mu
  )$effect
  expect_equal(mu + mapped_u, tutorial_alpha, tolerance = 1e-12)
})

test_that("unrestricted centering map has exact density and endpoints", {
  L <- matrix(c(
    1.7, -0.65, 0.45,
    0, 0.8, 0.35,
    0, 0, 1.25
  ), 3, 3)
  q <- c(-0.7, 1.1, 0.25)
  mu <- c(0.6, -0.35, 1.4)

  for (rho in list(rep(0, 3), c(0.15, 0.63, 0.9), rep(1, 3))) {
    mapped <- .re_center_map(q, L, rho, mu = mu)
    density <- .re_center_log_density(q, L, rho, mu = mu)
    expect_equal(density[[1]], density[[2]], tolerance = 1e-12)
    dense_map <- L %*% solve(mapped$A)
    expect_equal(
      mapped$log_jacobian,
      as.numeric(determinant(dense_map, logarithm = TRUE)$modulus),
      tolerance = 1e-12
    )
  }

  expect_equal(
    .re_center_map(q, L, rep(0, 3), mu = mu)$effect,
    drop(L %*% q)
  )
  expect_equal(.re_center_map(q, L, rep(1, 3), mu = mu)$effect, q - mu)
  expect_equal(
    mu + .re_center_map(q, L, rep(1, 3), mu = mu)$effect,
    q
  )
})

test_that("conditional Student scales remain inside the exact chart", {
  L_reference <- matrix(c(1.4, 0.55, 0, 0.75), 2, 2)
  q <- list(c(-0.4, 0.8), c(1.2, -0.3), c(0.1, 0.6))
  mixers <- c(0.45, 1.3, 2.2)
  rho <- c(0.27, 0.81)
  mu <- c(0.75, -0.2)

  densities <- vapply(seq_along(q), function(j) {
    L_j <- mixers[j] * L_reference
    values <- .re_center_log_density(q[[j]], L_j, rho, mu = mu)
    expect_equal(values[[1]], values[[2]], tolerance = 1e-12)
    values[[1]]
  }, numeric(1))
  expect_true(all(is.finite(densities)))

  for (j in seq_along(q)) {
    L_j <- mixers[j] * L_reference
    expect_equal(
      .re_center_map(q[[j]], L_j, c(0, 0), mu = mu)$effect,
      drop(L_j %*% q[[j]])
    )
    expect_equal(
      .re_center_map(q[[j]], L_j, c(1, 1), mu = mu)$effect,
      q[[j]] - mu
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
    mu <- c(theta[1] + 0.4 * theta[2], -0.3 * theta[1])
    c(.re_center_map(q, L, rho, mu = mu)$effect, theta)
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
  mu <- c(theta[1] + 0.4 * theta[2], -0.3 * theta[1])
  analytic <- .re_center_map(x[1:2], L, rho, mu = mu)$log_jacobian
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
  auto <- stancode(
    y ~ x + (1 + x | gr(g, center = "auto")),
    data = re_center_dat
  )

  expect_match2(centered, "matrix[M_1, N_1] z_1;")
  expect_match2(centered, "vector[M_1] mean_center_re_1;")
  expect_match2(
    centered,
    "mean_center_re_1[1] = Intercept - dot_product(means_X, b);"
  )
  expect_match2(centered, "mean_center_re_1[2] = b[1];")
  expect_match2(centered, "r_1[j] = (z_1[, j] - mean_center_re_1)';")
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
    "z_1[, j] - rho_center_re .* mean_center_re_1"
  )
  expect_match2(
    partial,
    "sum(log(diagonal(L_group_center_re))) - sum(log(diagonal(L_partial_center_re)))"
  )
  expect_false(grepl("mean_partial_s2z", partial, fixed = TRUE))

  expect_match2(
    auto, "matrix<lower=0,upper=1>[N_1, M_1] rho_s2z_1;"
  )
  expect_match2(auto, "rho_level_center_re = rho_s2z_1[j]';")
  expect_null(standata(
    y ~ x + (1 + x | gr(g, center = 0.35)), data = re_center_dat
  )$rho_s2z_1)
  auto_data <- standata(
    y ~ x + (1 + x | gr(g, center = "auto")), data = re_center_dat
  )
  expect_equal(unname(auto_data$rho_s2z_1), matrix(0, 6L, 2L))
  expect_identical(auto_data$compute_rho_center_candidate_1, 0L)
  resolved_rho <- brms:::.as_brmsautocenter_resolved(matrix(
    seq(0.1, 0.9, length.out = 12L), 6L, 2L,
    dimnames = list(levels(re_center_dat$g), c("Intercept", "x"))
  ))
  resolved_formula <- y ~ x + (1 + x | gr(g, center = resolved_rho))
  expect_identical(stancode(resolved_formula, data = re_center_dat), auto)
  resolved_data <- standata(resolved_formula, data = re_center_dat)
  expect_equal(
    unname(resolved_data$rho_s2z_1), unname(unclass(resolved_rho))
  )
  expect_identical(resolved_data$compute_rho_center_candidate_1, 0L)
  tpar <- .re_center_stan_between(
    auto, "transformed parameters {", "\nmodel {"
  )
  expect_false(grepl("Y[n]", tpar, fixed = TRUE))
  expect_false(grepl("rho_center_candidate_1", tpar, fixed = TRUE))
  generated <- strsplit(
    auto, "generated quantities {", fixed = TRUE
  )[[1L]][2L]
  before_generated <- strsplit(
    auto, "generated quantities {", fixed = TRUE
  )[[1L]][1L]
  expect_false(grepl(
    "matrix<lower=0,upper=1>[N_1, M_1] rho_center_candidate_1;",
    before_generated, fixed = TRUE
  ))
  expect_match2(generated, "if (compute_rho_center_candidate_1)")
  expect_match2(generated, "rho_center_candidate_1[j, k]")
  expect_match2(generated, "rho_center_candidate_1 = rho_s2z_1;")
})

test_that("ordinary independent and Student charts retain conditional scales", {
  rho <- matrix(
    seq(0.05, 0.95, length.out = 12), nrow = 6,
    dimnames = list(levels(re_center_dat$g), c("Intercept", "x"))
  )
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
      g, dist = "student", center = "auto"
    )),
    data = re_center_dat
  )
  student_fisher_independent <- stancode(
    y ~ x + (1 + x || gr(
      g, dist = "student", center = "auto"
    )),
    data = re_center_dat
  )
  student_heterogeneous <- stancode(
    y ~ x + (1 + x | gr(g, dist = "student", center = rho)),
    data = re_center_dat
  )
  student_independent_data <- standata(
    y ~ x + (1 + x || gr(g, dist = "student", center = rho)),
    data = re_center_dat
  )

  expect_match2(independent, "vector[N_1] denominator_center_re")
  expect_match2(independent, "vector[N_1] rho_group_center_re")
  expect_match2(
    independent,
    "z_1[1] - rho_group_center_re * mean_center_re_1[1]"
  )
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
  expect_match2(
    student_fisher_independent,
    "rho_group_center_re = rho_s2z_1[, 1];"
  )
  expect_match2(
    student_fisher_independent,
    "rho_group_center_re + rho_group_center_re .* scale_group_center_re"
  )
  expect_match2(
    student_fisher_independent,
    "z_1[1] - rho_group_center_re * mean_center_re_1[1]"
  )
  expect_match2(
    student_heterogeneous,
    "rho_level_center_re = rho_s2z_1[j]';"
  )
  expect_equal(unname(student_independent_data$rho_s2z_1), unname(rho))
  fisher_comp <- .re_center_stan_between(
    student_fisher, "L_center_re_1 =", "L_group_center_re *= dfm_1[j];"
  )
  expect_false(grepl("dfm_1", fisher_comp, fixed = TRUE))
  student_tpar <- .re_center_stan_between(
    student_fisher, "transformed parameters {", "\nmodel {"
  )
  expect_false(grepl(
    "rho_center_candidate_1", student_tpar, fixed = TRUE
  ))
})

test_that("level-specific fixed fractions align and reach ordinary charts", {
  levels_g <- levels(re_center_dat$g)
  rho <- matrix(
    c(seq(0.1, 0.6, length.out = 6), seq(0.85, 0.35, length.out = 6)),
    nrow = 6, dimnames = list(levels_g, c("Intercept", "x"))
  )
  reversed <- rho[rev(levels_g), c("x", "Intercept"), drop = FALSE]
  sdata <- standata(
    y ~ x + (1 + x | gr(g, center = reversed)), data = re_center_dat
  )
  expect_equal(unname(sdata$rho_s2z_1), unname(rho))
  code <- stancode(
    y ~ x + (1 + x | gr(g, center = reversed)), data = re_center_dat
  )
  expect_match2(
    code, "matrix<lower=0,upper=1>[N_1, M_1] rho_s2z_1;"
  )
  expect_match2(code, "rho_level_center_re = rho_s2z_1[j]';")
  expect_match2(
    code,
    "z_1[, j] - rho_level_center_re .* mean_center_re_1"
  )

  s2z_data <- standata(
    y ~ x + (1 + x | gr(g, s2z = TRUE, center = reversed)),
    data = re_center_dat
  )$rho_s2z_1
  expect_equal(unname(s2z_data), unname(rho))

  fit <- brm(
    y ~ x + (1 + x | gr(g, center = reversed)),
    data = re_center_dat, empty = TRUE, cores = 1
  )
  expect_equal(
    unname(fit$basis$dpars$mu$re_s2z_center$rho_s2z_1),
    unname(rho)
  )

  shared <- setNames(rev(rho[, 1]), rev(levels_g))
  shared_data <- standata(
    y ~ x + (1 + x || gr(g, center = shared)), data = re_center_dat
  )$rho_s2z_1
  expect_equal(unname(shared_data[, 1]), unname(rho[, 1]))
  expect_equal(unname(shared_data[, 2]), unname(rho[, 1]))

  programmatic_data <- (function(center) {
    standata(
      y ~ x + (1 + x | gr(g, center = center)), data = re_center_dat
    )$rho_s2z_1
  })(shared)
  expect_equal(unname(programmatic_data[, 1]), unname(rho[, 1]))
  expect_equal(unname(programmatic_data[, 2]), unname(rho[, 1]))

  expect_error(
    standata(
      y ~ x + (1 + x | gr(g, center = rep(0.5, 2))),
      data = re_center_dat
    ),
    "one value per grouping level"
  )
  malformed_constant <- matrix(
    0.5, 2, 2, dimnames = list(c("wrong-1", "wrong-2"), c("Intercept", "x"))
  )
  expect_error(
    standata(
      y ~ x + (1 + x | gr(g, center = malformed_constant)),
      data = re_center_dat
    ),
    "Row names.*must match"
  )
  bad_coef <- rho
  colnames(bad_coef) <- c("Intercept", "wrong")
  expect_error(
    standata(
      y ~ x + (1 + x | gr(g, center = bad_coef)), data = re_center_dat
    ),
    "Column names.*must match"
  )
})

test_that("population locations follow the fitted fixed-to-random basis", {
  cell_means <- stancode(
    y ~ f + (0 + f | gr(g, center = TRUE)), data = re_center_dat
  )
  expect_match2(
    cell_means,
    "mean_center_re_1[1] = Intercept - dot_product(means_X, b);"
  )
  expect_match2(
    cell_means,
    paste0(
      "mean_center_re_1[2] = Intercept - ",
      "dot_product(means_X, b) + b[1];"
    )
  )

  bframe <- brms:::brmsframe(
    brms:::brmsterms(bf(y ~ f + (0 + f | gr(g, center = TRUE)))),
    re_center_dat
  )
  expect_equal(
    bframe$dpars$mu$frame$re_center_mean[["1"]],
    matrix(
      c(1, 1, 0, 1), nrow = 2,
      dimnames = list(c("fa", "fb"), c("Intercept", "fb"))
    )
  )
})

test_that("fitted level-specific maps are self-contained", {
  levels_g <- levels(re_center_dat$g)
  rho <- matrix(
    seq(0.1, 0.9, length.out = 12), nrow = 6,
    dimnames = list(levels_g, c("Intercept", "x"))
  )
  key <- ".brms_re_center_test_weights"
  assign(key, rho, envir = .GlobalEnv)
  on.exit({
    if (exists(key, envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = key, envir = .GlobalEnv)
    }
  }, add = TRUE)
  form <- y ~ x +
    (1 + x | gr(g, center = .brms_re_center_test_weights))
  environment(form) <- .GlobalEnv
  fit <- brm(form, data = re_center_dat, empty = TRUE, cores = 1)
  rm(list = key, envir = .GlobalEnv)

  expect_equal(unname(standata(fit)$rho_s2z_1), unname(rho))
  expect_silent(standata(fit, re_formula = ~(1 + x | g)))
  path <- tempfile(fileext = ".rds")
  saveRDS(fit, path)
  restored <- readRDS(path)
  expect_equal(unname(standata(restored)$rho_s2z_1), unname(rho))
  expect_silent(standata(restored, re_formula = ~(1 + x | g)))

  rho_new <- setNames(seq(0.2, 0.7, length.out = 6), levels_g)
  updated <- update(
    fit,
    formula. = y ~ x + (1 + x | gr(g, center = rho_new)),
    testmode = TRUE, recompile = FALSE, algorithm = "meanfield", cores = 1
  )
  rm(rho_new)
  expect_equal(
    unname(standata(updated, newdata = re_center_dat)$rho_s2z_1[, 1]),
    seq(0.2, 0.7, length.out = 6)
  )
})

test_that("fixed ordinary centering composes with supported predictor types", {
  rho_two <- matrix(
    seq(0.05, 0.95, length.out = 12), nrow = 6,
    dimnames = list(levels(re_center_dat$g), c("Intercept", "x"))
  )
  rho_one <- matrix(
    seq(0.1, 0.9, length.out = 6), ncol = 1,
    dimnames = list(levels(re_center_dat$g), "shared")
  )
  ordinal_dat <- transform(
    re_center_dat,
    yo = ordered(rep(c("low", "middle", "high"), length.out = 48))
  )
  ordinal <- stancode(
    yo ~ x + (1 + x | gr(g, center = rho_two)),
    data = ordinal_dat, family = cumulative()
  )
  expect_match2(ordinal, "ordered[nthres] Intercept")
  expect_match2(ordinal, "L_partial_center_re")
  expect_match2(ordinal, "mean_center_re_1[2] = b[1];")
  expect_false(grepl("mean_center_re_1[1] = Intercept", ordinal, fixed = TRUE))

  qr <- stancode(
    bf(y ~ x + (1 + x | gr(g, center = 0.4)), decomp = "QR"),
    data = re_center_dat
  )
  expect_match2(qr, "vector[Kc] b_center_re_1;")
  expect_match2(qr, "b_center_re_1 = XR_inv * bQ;")
  expect_match2(
    qr,
    "Intercept - dot_product(means_X, b_center_re_1)"
  )

  no_population_intercept <- stancode(
    y ~ 0 + x + (1 + x | gr(g, center = 0.4)), data = re_center_dat
  )
  expect_match2(no_population_intercept, "mean_center_re_1[2] = b[1];")
  expect_false(grepl(
    "mean_center_re_1[1] =", no_population_intercept, fixed = TRUE
  ))

  distributional <- stancode(
    bf(
      y ~ x + (1 | gr(g, id = "mu-local", center = rho_one)),
      sigma ~ x + (1 | gr(g, id = "sigma-local", center = rho_one))
    ),
    data = re_center_dat, family = gaussian()
  )
  expect_match2(distributional, "log_jacobian_re_1")
  expect_match2(distributional, "log_jacobian_re_2")
  expect_match2(
    distributional,
    "mean_center_re_2[1] = Intercept_sigma - dot_product(means_X_sigma, b_sigma);"
  )

  multivariate <- stancode(
    bf(y ~ x + (1 | gr(g, id = "y-local", center = rho_one))) +
      bf(y2 ~ x + (1 | gr(g, id = "y2-local", center = rho_one))) +
      set_rescor(FALSE),
    data = re_center_dat
  )
  expect_match2(multivariate, "log_jacobian_re_1")
  expect_match2(multivariate, "log_jacobian_re_2")
  expect_match2(multivariate, "mean_center_re_1[1] = Intercept_y")
  expect_match2(multivariate, "mean_center_re_2[1] = Intercept_y2")

  multivariate_fisher <- stancode(
    bf(y ~ x + (1 | gr(g, id = "y-fisher", center = "auto"))) +
      bf(y2 ~ x + (1 | gr(g, id = "y2-fisher", center = "auto"))) +
      set_rescor(FALSE),
    data = re_center_dat
  )
  expect_match2(multivariate_fisher, "rho_s2z_1[j, 1]")
  expect_match2(multivariate_fisher, "rho_s2z_2[j, 1]")
  expect_error(
    stancode(
      bf(y ~ x + (1 | gr(g, center = "auto"))) +
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
    eta ~ 1 + x + (1 | gr(g, center = rho_one)),
    nl = TRUE
  )
  nonlinear_prior <- prior(normal(0, 1), class = "b", nlpar = "eta")
  nonlinear <- stancode(
    nonlinear_form, data = re_center_dat, prior = nonlinear_prior
  )
  expect_match2(nonlinear, "log_jacobian_re_1")
  expect_match2(nonlinear, "mean_center_re_1[1] = b_eta[1];")
  expect_error(
    stancode(
      bf(
        y ~ eta,
        eta ~ 1 + x + (1 | gr(g, center = "auto")),
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
      yo ~ x + (1 | gr(g, center = "auto")),
      data = ordinal_dat, family = cumulative()
    ),
    "not yet supported.*ordinal location"
  )
  expect_error(
    stancode(
      y ~ x +
        (1 | gr(g, id = "mixed", center = "auto")) +
        (0 + x | gr(g, id = "mixed", center = 0.5)),
      data = re_center_dat
    ),
    "must use.*centering if any coefficient does"
  )
  expect_error(
    stancode(
      y ~ x +
        (1 | gr(g, id = "ordinary", center = 0.4)) +
        (1 | gr(h, id = "s2z", s2z = TRUE)),
      data = re_center_dat
    ),
    "cannot yet be combined with a conventional sum-to-zero"
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
  resolved <- brms:::.as_brmsautocenter_resolved(matrix(
    seq(0.1, 0.9, length.out = 12L), 6L, 2L,
    dimnames = list(levels(re_center_dat$g), c("Intercept", "x"))
  ))
  fit <- brm(
    y ~ x + (1 + x | gr(g, center = resolved)),
    data = re_center_dat, empty = TRUE, cores = 1
  )
  excluded <- brms:::exclude_pars(fit)
  expect_true(all(c(
    "z_1", "L_center_re_1", "log_jacobian_re_1",
    "mean_center_re_1", "b_center_re_1",
    "rho_center_candidate_1", "mean_rho_center_candidate_1"
  ) %in% excluded))

  fit_all <- brm(
    y ~ x + (1 + x | gr(g, center = resolved)),
    data = re_center_dat, empty = TRUE, cores = 1,
    save_pars = save_pars(all = TRUE)
  )
  saved_exclusions <- brms:::exclude_pars(fit_all)
  expect_false(any(c(
    "z_1", "L_center_re_1", "log_jacobian_re_1",
    "mean_center_re_1", "b_center_re_1",
    "rho_center_candidate_1", "mean_rho_center_candidate_1"
  ) %in% saved_exclusions))
})
