context("Distribution validation: analytical / known values")

test_that("gaussian matches stats:: reference", {
  entry <- dist_registry_get("gaussian")[[1]]
  x <- c(-2, -0.5, 0, 1.2)
  expect_equal(
    call_dist(entry$d, x, entry$params),
    dnorm(x, mean = 0, sd = 1)
  )
  expect_equal(
    call_dist(entry$p, x, entry$params),
    pnorm(x, mean = 0, sd = 1)
  )
  p <- c(0.1, 0.5, 0.9)
  expect_equal(
    call_dist(entry$q, p, entry$params),
    qnorm(p, mean = 0, sd = 1)
  )
})

test_that("student_t matches scaled stats::t", {
  entry <- dist_registry_get("student_t")[[1]]
  x <- c(-1, 0, 2)
  df <- entry$params$df
  sigma <- entry$params$sigma
  mu <- entry$params$mu
  expect_equal(
    call_dist(entry$d, x, entry$params),
    dt((x - mu) / sigma, df = df) / sigma
  )
  expect_equal(
    call_dist(entry$p, x, entry$params),
    pt((x - mu) / sigma, df = df)
  )
})

test_that("com_poisson reduces to Poisson when shape = 1", {
  entry <- dist_registry_get("com_poisson")[[1]]
  mu <- 3
  p <- c(0.1, 0.5, 0.9)
  expect_equal(
    brms:::qcom_poisson(p, mu = mu, shape = 1),
    qpois(p, lambda = mu)
  )
  x <- 0:8
  expect_equal(
    brms:::dcom_poisson(x, mu = mu, shape = 1),
    dpois(x, lambda = mu),
    tolerance = 1e-6
  )
})

test_that("qexgaussian handles probabilities outside (0, 1)", {
  entry <- dist_registry_get("exgaussian")[[1]]
  res <- SW(call_dist(entry$q, c(-1, 0, 1, 1.5), entry$params))
  expect_equal(res, c(NA, -Inf, Inf, NA))
})

test_that("qvon_mises handles probabilities outside (0, 1) and stays in (-pi, pi)", {
  entry <- dist_registry_get("von_mises")[[1]]
  q <- call_dist(entry$q, c(0.1, 0.5, 0.9), entry$params)
  expect_true(all(q >= -pi & q <= pi))
  res <- SW(call_dist(entry$q, c(-1, 0, 0.5, 1, 1.5), entry$params))
  expect_equal(res[1], NA_real_)
  expect_equal(res[2], -pi)
  expect_equal(res[4], pi)
  expect_equal(res[5], NA_real_)
})


test_that("qzero_one_inflated_beta matches mixture CDF boundaries", {
  entry <- dist_registry_get("zero_one_inflated_beta")[[1]]
  zoi <- entry$params$zoi
  coi <- entry$params$coi
  p0 <- zoi * (1 - coi)
  p1 <- zoi * coi
  p <- c(0.01, 0.1, 0.2, 0.5, 0.9, 0.95, 0.99)
  q <- call_dist(entry$q, p, entry$params)
  expect_equal(q[p <= p0], rep(0, sum(p <= p0)))
  expect_equal(q[p >= 1 - p1], rep(1, sum(p >= 1 - p1)))
  expect_true(all(q[p > p0 & p < 1 - p1] > 0 & q[p > p0 & p < 1 - p1] < 1))
  expect_true(all(diff(q) >= 0))
  F_q <- call_dist(entry$p, q, entry$params)
  expect_true(all(as.numeric(F_q) + 1e-8 >= p))
})

test_that("qzero_inflated_asym_laplace matches mixture CDF regions", {
  entry <- dist_registry_get("zero_inflated_asym_laplace")[[1]]
  mu <- entry$params$mu
  sigma <- entry$params$sigma
  quantile <- entry$params$quantile
  zi <- entry$params$zi
  p <- c(0.01, 0.1, 0.3, 0.5, 0.65, 0.85, 0.99)
  q <- call_dist(entry$q, p, entry$params)
  F0_base <- pasym_laplace(0, mu = mu, sigma = sigma, quantile = quantile)
  F0_minus <- (1 - zi) * F0_base
  F0 <- zi + (1 - zi) * F0_base
  expect_true(all(q[p < F0_minus] < 0))
  expect_equal(q[p >= F0_minus & p <= F0],
               rep(0, sum(p >= F0_minus & p <= F0)))
  expect_true(all(q[p > F0] > 0))
  expect_true(all(diff(q) >= 0))
  F_q <- call_dist(entry$p, q, entry$params)
  expect_true(all(as.numeric(F_q) + 1e-8 >= p))
})

test_that("qordinal matches the ordinal CDF", {
  set.seed(11)
  ns <- 7
  nthres <- 3
  eta <- rnorm(ns)
  thres <- matrix(rep(c(-1, 0, 1), each = ns), nrow = ns)
  disc <- rep(1, ns)
  p <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  for (fam in c("cumulative", "sratio", "cratio", "acat")) {
    F_mat <- brms:::pordinal(
      seq_len(nthres + 1), eta = eta, thres = thres, disc = disc,
      family = fam, link = "logit"
    )
    for (pj in p) {
      qj <- brms:::qordinal(
        pj, eta = eta, thres = thres, disc = disc,
        family = fam, link = "logit"
      )
      expect_length(qj, ns)
      expect_true(all(qj >= 1 & qj <= nthres + 1))
      expect_true(all(F_mat[cbind(seq_len(ns), qj)] + 1e-10 >= pj))
      expect_true(all(
        ifelse(qj > 1, F_mat[cbind(seq_len(ns), qj - 1)] < pj + 1e-10, TRUE)
      ))
    }
  }
})

test_that("qcategorical matches the categorical CDF", {
  set.seed(12)
  ns <- 7
  ncat <- 4
  eta <- matrix(rnorm(ns * ncat), nrow = ns)
  p <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  F_mat <- brms:::pcategorical(seq_len(ncat), eta = eta)
  for (pj in p) {
    qj <- brms:::qcategorical(pj, eta = eta)
    expect_length(qj, ns)
    expect_true(all(qj >= 1 & qj <= ncat))
    expect_true(all(F_mat[cbind(seq_len(ns), qj)] + 1e-10 >= pj))
    expect_true(all(
      ifelse(qj > 1, F_mat[cbind(seq_len(ns), qj - 1)] < pj + 1e-10, TRUE)
    ))
  }
})

test_that("qhurdle_cumulative matches the mixture CDF", {
  set.seed(13)
  ns <- 7
  nthres <- 3
  ncat <- nthres + 1
  eta <- rnorm(ns)
  thres <- matrix(rep(c(-1, 0, 1), each = ns), nrow = ns)
  disc <- rep(1, ns)
  hu <- rep(0.25, ns)
  p <- c(0.05, 0.2, 0.25, 0.5, 0.9, 0.99)
  for (pj in p) {
    qj <- brms:::qhurdle_cumulative(
      pj, eta = eta, thres = thres, hu = hu, disc = disc, link = "logit"
    )
    expect_length(qj, ns)
    expect_true(all(qj >= 0 & qj <= ncat))
    F_q <- vapply(seq_len(ns), function(s) {
      brms:::phurdle_cumulative(
        qj[s], eta = eta[s], thres = thres[s, , drop = FALSE],
        hu = hu[s], disc = disc[s], link = "logit"
      )
    }, numeric(1))
    expect_true(all(F_q + 1e-10 >= pj))
  }
  expect_equal(
    brms:::qhurdle_cumulative(
      0.1, eta = eta, thres = thres, hu = hu, disc = disc, link = "logit"
    ),
    rep(0, ns)
  )
})

test_that("beta_binomial matches extraDistr when available", {
  skip_if_not_installed("extraDistr")
  entry <- dist_registry_get("beta_binomial")[[1]]
  size <- entry$params$size
  mu <- entry$params$mu
  phi <- entry$params$phi
  # extraDistr uses alpha, beta parameterization
  alpha <- mu * phi
  beta <- (1 - mu) * phi
  x <- 0:size
  expect_equal(
    as.numeric(call_dist(entry$d, x, entry$params)),
    extraDistr::dbbinom(x, size = size, alpha = alpha, beta = beta)
  )
  expect_equal(
    as.numeric(call_dist(entry$p, x, entry$params)),
    extraDistr::pbbinom(x, size = size, alpha = alpha, beta = beta)
  )
})

test_that("lognormal matches stats::lnorm", {
  entry <- dist_registry_get("lognormal")[[1]]
  x <- c(0.5, 1, 2, 4)
  expect_equal(
    call_dist(entry$d, x, entry$params),
    dlnorm(x, meanlog = 0, sdlog = 1)
  )
  expect_equal(
    call_dist(entry$p, x, entry$params),
    plnorm(x, meanlog = 0, sdlog = 1)
  )
  p <- c(0.1, 0.5, 0.9)
  expect_equal(
    call_dist(entry$q, p, entry$params),
    qlnorm(p, meanlog = 0, sdlog = 1)
  )
})

test_that("geometric equals nbinom with size = 1", {
  entry <- dist_registry_get("geometric")[[1]]
  mu <- entry$params$mu
  x <- 0:12
  expect_equal(
    call_dist(entry$d, x, entry$params),
    dnbinom(x, mu = mu, size = 1)
  )
  expect_equal(
    call_dist(entry$p, x, entry$params),
    pnbinom(x, mu = mu, size = 1)
  )
})

test_that("zi = 0 reduces zero_inflated_poisson to poisson", {
  entry <- dist_registry_get("zero_inflated_poisson")[[1]]
  lambda <- entry$params$lambda
  x <- 0:10
  expect_equal(
    dzero_inflated_poisson(x, lambda = lambda, zi = 0),
    dpois(x, lambda = lambda)
  )
  expect_equal(
    pzero_inflated_poisson(x, lambda = lambda, zi = 0),
    ppois(x, lambda = lambda)
  )
  p <- c(0.1, 0.5, 0.9)
  expect_equal(
    qzero_inflated_poisson(p, lambda = lambda, zi = 0),
    qpois(p, lambda = lambda)
  )
})

test_that("qhurdle_negbinomial recovers special cases", {
  mu <- 4
  shape <- 2.5
  # hu = 0 => same as truncated-at-zero-removed nbinom quantile path via
  # qhurdle; for p > 0 the positive part matches renormalized nbinom.
  # hu = 1 => all mass at 0
  expect_equal(
    qhurdle_negbinomial(c(0.1, 0.5, 0.9), mu = mu, shape = shape, hu = 1),
    c(0, 0, 0)
  )
  # hu small: values above 0 for p > hu
  q <- qhurdle_negbinomial(c(0.05, 0.3, 0.8), mu = mu, shape = shape, hu = 0.1)
  expect_equal(q[1], 0)
  expect_true(all(q[2:3] > 0))
  expect_true(all(diff(q) >= 0))
})
