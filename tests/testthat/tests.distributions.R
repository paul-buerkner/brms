context("Tests for distribution functions")

test_that("student distribution works correctly", {
  expect_equal(integrate(dstudent_t, -100, 100, df = 15, mu = 10, sigma = 5)$value, 1)
  expect_equal(dstudent_t(1, df = 10, mu = 0, sigma = 5), dt(1/5, df = 10)/5)
  expect_equal(pstudent_t(2, df = 20, mu = 2, sigma = 0.4), pt(0, df = 20))
  ps <- c(-1,0,0.7,1,1.5)
  SW(expect_equal(qstudent_t(ps, df = 5, mu = 2, sigma = 3), 2 + 3*qt(ps, df = 5)))
  expect_equal(length(rstudent_t(10, df = 10, mu = rnorm(10), sigma = 1:10)), 10)
})

test_that("multivariate normal and student distributions work correctly", {
  mu <- rnorm(3)
  Sigma <- cov(matrix(rnorm(300), ncol = 3))
  expect_equal(dmulti_normal(1:3, mu = mu, Sigma = Sigma),
               mnormt::dmnorm(1:3, mu, Sigma))
  expect_equal(dmulti_student_t(1:3, mu = mu, Sigma = Sigma, df = 10, log = TRUE),
               mnormt::dmt(1:3, df = 10, mean = mu, S = Sigma, log = TRUE))
  expect_equal(dim(rmulti_normal(7, mu = mu, Sigma = Sigma)), c(7, 3))
  expect_equal(dim(rmulti_student_t(7, mu = mu, Sigma = Sigma, df = 10)),
               c(7, 3))
  # test errors
  expect_error(dmulti_normal(1:3, mu = rnorm(2), Sigma = Sigma, check = TRUE),
               "Dimension of mu is incorrect")
  expect_error(dmulti_normal(1:3, mu = mu, Sigma = Sigma[1:2, 1:2],
                             check = TRUE),
               "Dimension of Sigma is incorrect")
  expect_error(dmulti_normal(1:3, mu = mu, Sigma = Sigma[1:3, 3:1],
                             check = TRUE),
               "Sigma must be a symmetric matrix")
  expect_error(rmulti_normal(1.5, mu = mu, Sigma = Sigma, check = TRUE),
               "n must be a positive integer")
  expect_error(rmulti_normal(10, mu = mu, Sigma = Sigma[1:3, 3:1],
                             check = TRUE),
               "Sigma must be a symmetric matrix")
  expect_error(dmulti_student_t(rnorm(3), mu = mu, Sigma = Sigma,
                                df = -1, check = TRUE),
               "df must be greater than 0")
  expect_error(dmulti_student_t(rnorm(3), mu = mu, Sigma = Sigma[1:3, 3:1],
                                df = 30, check = TRUE),
               "Sigma must be a symmetric matrix")
  expect_error(rmulti_student_t(10, mu = mu, Sigma = Sigma,
                                df = -1, check = TRUE),
               "df must be greater than 0")
})

test_that("von_mises distribution functions run without errors", {
  n <- 10
  res <- dvon_mises(runif(n, -pi, pi), mu = 1, kappa = 1:n)
  expect_true(length(res) == n)
  res <- pvon_mises(runif(n, -pi, pi), mu = rnorm(n), kappa = 0:(n-1))
  expect_true(length(res) == n)
  res <- rvon_mises(n, mu = rnorm(n), kappa = 0:(n-1))
  expect_true(length(res) == n)
})

test_that("skew_normal distribution functions run without errors", {
  n <- 10
  x <- rnorm(n, 10, 3)
  res <- dskew_normal(x, mu = 1, sigma = 2, alpha = 1)
  expect_true(length(res) == n)
  res <- pskew_normal(x, mu = rnorm(n), sigma = 1:n,
                      alpha = 3, log.p = TRUE)
  expect_true(length(res) == n)
  p <- log(runif(n, 0, 1))
  res <- qskew_normal(p, mu = rnorm(n), sigma = 1:n,
                      alpha = 3, log.p = TRUE)
  expect_true(length(res) == n)
  ps <- c(-1, 0, 0.5, 1, 1.5)
  res <- SW(qskew_normal(ps))
  expect_equal(res, c(NA, -Inf, 0, Inf, NA))
  res <- rskew_normal(n, mu = rnorm(n), sigma = 10, alpha = -4:5)
  expect_true(length(res) == n)
})

test_that("exgaussian distribution functions run without errors", {
  n <- 10
  x <- rnorm(n, 10, 3)
  res <- dexgaussian(x, mu = 1, sigma = 2, beta = 1)
  expect_true(length(res) == n)
  res <- pexgaussian(x, mu = rnorm(n), sigma = 1:n,
                     beta = 3, log.p = TRUE)
  expect_true(length(res) == n)
  res <- rexgaussian(n, mu = rnorm(n), sigma = 10, beta = 1:10)
  expect_true(length(res) == n)
})

test_that("frechet distribution functions run without errors", {
  n <- 10
  x <- 21:30
  res <- dfrechet(x, loc = 1, scale = 2, shape = 1, log = TRUE)
  expect_true(length(res) == n)
  loc <- 1:10
  res <- pfrechet(x, loc = loc, scale = 1:n, shape = 3)
  expect_true(length(res) == n)
  q <- qfrechet(res, loc = loc, scale = 1:n, shape = 3)
  expect_equal(x, q)
  ps <- c(-1, 0, 1, 1.5)
  res <- SW(qfrechet(ps))
  expect_equal(res, c(NA, 0, Inf, NA))
  res <- rfrechet(n, loc = loc, scale = 10, shape = 1:10)
  expect_true(length(res) == n)
})

test_that("inv_gaussian distribution functions run without errors", {
  n <- 10
  x <- rgamma(n, 10, 3)
  res <- dinv_gaussian(x, mu = 1, shape = 1)
  expect_true(length(res) == n)
  res <- pinv_gaussian(x, mu = abs(rnorm(n)), shape = 3)
  expect_true(length(res) == n)
  res <- rinv_gaussian(n, mu = abs(rnorm(n)), shape = 1:10)
  expect_true(length(res) == n)
})

test_that("beta_binomial distribution functions run without errors", {
  skip_if_not_installed("extraDistr")

  n <- 10
  x <- rpois(n, lambda = 1)

  res <- dbeta_binomial(x, c(2, 10), mu = 0.4, phi = 1)
  expect_true(length(res) == n)
  res <- pbeta_binomial(x, c(2, 10), mu = 0.4, phi = 1)
  expect_true(length(res) == n)
  res <- rbeta_binomial(n, c(2, 10), mu = 0.4, phi = 1)
  expect_true(length(res) == n)
})

test_that("gen_extreme_value distribution functions run without errors", {
  n <- 10
  x <- rgamma(n, 10, 3)
  res <- dgen_extreme_value(x, mu = 1, sigma = 2, xi = 1)
  expect_true(length(res) == n)
  mu <- rnorm(n)
  res <- pgen_extreme_value(x, mu = mu, sigma = 1:n, xi = 3)
  expect_true(length(res) == n)
  q <- qgen_extreme_value(res, mu = mu, sigma = 1:n, xi = 3)
  expect_equal(x, q)
  res <- rgen_extreme_value(n, mu = mu, sigma = 10, xi = 1:10)
  expect_true(length(res) == n)
})

test_that("asym_laplace distribution functions run without errors", {
  n <- 10
  x <- rnorm(n, 10, 3)
  res <- dasym_laplace(x, mu = 1, sigma = 2, quantile = 0.5)
  expect_true(length(res) == n)
  res <- pasym_laplace(x, mu = rnorm(n), sigma = 1:n, quantile = 0.3)
  expect_true(length(res) == n)
  res <- rasym_laplace(n, mu = rnorm(n), sigma = 10,
                       quantile = runif(n, 0, 1))
  expect_true(length(res) == n)
})

test_that("zero-inflated distribution functions run without errors", {
  skip_if_not_installed("extraDistr")
  n <- 10
  x <- rpois(n, lambda = 1)

  res <- dzero_inflated_poisson(x, lambda = 1, zi = 0.1)
  expect_true(length(res) == n)
  res <- pzero_inflated_poisson(x, lambda = 1, zi = 0.1)
  expect_true(length(res) == n)

  res <- dzero_inflated_negbinomial(x, mu = 2, shape = 5, zi = 0.1)
  expect_true(length(res) == n)
  res <- pzero_inflated_negbinomial(x, mu = 2, shape = 5, zi = 0.1)
  expect_true(length(res) == n)

  res <- dzero_inflated_binomial(x, size = c(2, 10), prob = 0.4, zi = 0.1)
  expect_true(length(res) == n)
  res <- pzero_inflated_binomial(x, size = c(2, 10), prob = 0.4, zi = 0.1)
  expect_true(length(res) == n)

  res <- dzero_inflated_beta_binomial(x, c(2, 10), mu = 0.4, phi = 1, zi = 0.1)
  expect_true(length(res) == n)
  res <- pzero_inflated_beta_binomial(x, c(2, 10), mu = 0.4, phi = 1, zi = 0.1)
  expect_true(length(res) == n)

  x <- c(rbeta(n - 2, shape1 = 2, shape2 = 3), 0, 0)
  res <- dzero_inflated_beta(x, shape1 = 2, shape2 = 3, zi = 0.1)
  expect_true(length(res) == n)
  res <- pzero_inflated_beta(x, shape1 = 2, shape2 = 3, zi = 0.1)
  expect_true(length(res) == n)
})

test_that("hurdle distribution functions run without errors", {
  n <- 10
  x <- rpois(n, lambda = 1)

  res <- dhurdle_poisson(x, lambda = 1, hu = 0.1)
  expect_true(length(res) == n)
  res <- phurdle_poisson(x, lambda = 1, hu = 0.1)
  expect_true(length(res) == n)

  res <- dhurdle_negbinomial(x, mu = 2, shape = 5, hu = 0.1)
  expect_true(length(res) == n)
  res <- phurdle_negbinomial(x, mu = 2, shape = 5, hu = 0.1)
  expect_true(length(res) == n)

  res <- dhurdle_gamma(x, shape = 1, scale = 3, hu = 0.1)
  expect_true(length(res) == n)
  res <- phurdle_gamma(x, shape = 1, scale = 3, hu = 0.1)
  expect_true(length(res) == n)

  res <- dhurdle_lognormal(x, mu = 2, sigma = 5, hu = 0.1)
  expect_true(length(res) == n)
  res <- phurdle_lognormal(x, mu = 2, sigma = 5, hu = 0.1)
  expect_true(length(res) == n)
})

test_that("wiener distribution functions run without errors", {
  skip_if_not_installed("RWiener")

  set.seed(1234)
  n <- 10
  x <- seq(0.1, 1, 0.1)
  alpha <- rexp(n)
  tau <- 0.05
  beta <- 0.5
  delta <- rnorm(n)
  resp <- sample(c(0, 1), n, TRUE)

  d1 <- dwiener(x, alpha, tau, beta, delta, resp, backend = "Rwiener")
  d2 <- dwiener(x, alpha, tau, beta, delta, resp, backend = "rtdists")
  expect_equal(d1, d2)

  r1 <- rwiener(n, alpha, tau, beta, delta, backend = "Rwiener")
  r2 <- rwiener(n, alpha, tau, beta, delta, backend = "rtdists")
  expect_equal(names(r1), names(r2))
  expect_equal(dim(r1), dim(r2))
})

test_that("d<ordinal_family>() works correctly", {
  source(testthat::test_path(file.path("helpers", "inv_link_ordinal_ch.R")))
  source(testthat::test_path(file.path("helpers", "simopts_catlike_oneobs.R")))
  for (ndraws in ndraws_vec) {
    for (ncat in ncat_vec) {
      thres_test <- matrix(rnorm(ndraws * (ncat - 1)), nrow = ndraws)
      # Emulate no category-specific effects (i.e., only a single vector of
      # linear predictors) as well as category-specific effects (i.e., a matrix
      # of linear predictors):
      eta_test_list <- list(
        rnorm(ndraws),
        matrix(rnorm(ndraws * (ncat - 1)), nrow = ndraws)
      )
      for (eta_test in eta_test_list) {
        thres_eta <- if (is.matrix(eta_test)) {
          stopifnot(identical(dim(eta_test), dim(thres_test)))
          thres_test - eta_test
        } else {
          # Just to try something different:
          sweep(thres_test, 1, as.array(eta_test))
        }
        eta_thres <- if (is.matrix(eta_test)) {
          stopifnot(identical(dim(eta_test), dim(thres_test)))
          eta_test - thres_test
        } else {
          # Just to try something different:
          sweep(-thres_test, 1, as.array(eta_test), FUN = "+")
        }
        for (link in c("logit", "probit", "cauchit", "cloglog")) {
          # cumulative():
          d_cumul <- dcumulative(seq_len(ncat),
                                 eta_test, thres_test, link = link)
          d_cumul_ch <- inv_link_cumulative_ch(thres_eta, link = link)
          expect_equivalent(d_cumul, d_cumul_ch)
          expect_equal(dim(d_cumul), c(ndraws, ncat))

          # sratio():
          d_sratio <- dsratio(seq_len(ncat),
                              eta_test, thres_test, link = link)
          d_sratio_ch <- inv_link_sratio_ch(thres_eta, link = link)
          expect_equivalent(d_sratio, d_sratio_ch)
          expect_equal(dim(d_sratio), c(ndraws, ncat))

          # cratio():
          d_cratio <- dcratio(seq_len(ncat),
                              eta_test, thres_test, link = link)
          d_cratio_ch <- inv_link_cratio_ch(eta_thres, link = link)
          expect_equivalent(d_cratio, d_cratio_ch)
          expect_equal(dim(d_cratio), c(ndraws, ncat))

          # acat():
          d_acat <- dacat(seq_len(ncat),
                          eta_test, thres_test, link = link)
          d_acat_ch <- inv_link_acat_ch(eta_thres, link = link)
          expect_equivalent(d_acat, d_acat_ch)
          expect_equal(dim(d_acat), c(ndraws, ncat))
        }
      }
    }
  }
})

test_that("inv_link_<ordinal_family>() works correctly for arrays", {
  source(testthat::test_path(file.path("helpers", "inv_link_ordinal_ch.R")))
  source(testthat::test_path(file.path("helpers", "simopts_catlike.R")))
  for (ndraws in ndraws_vec) {
    for (nobsv in nobsv_vec) {
      for (ncat in ncat_vec) {
        x_test <- array(rnorm(ndraws * nobsv * (ncat - 1)),
                        dim = c(ndraws, nobsv, ncat - 1))
        nx_test <- -x_test
        for (link in c("logit", "probit", "cauchit", "cloglog")) {
          # cumulative():
          il_cumul <- inv_link_cumulative(x_test, link = link)
          il_cumul_ch <- inv_link_cumulative_ch(x_test, link = link)
          expect_equivalent(il_cumul, il_cumul_ch)
          expect_equal(dim(il_cumul), c(ndraws, nobsv, ncat))

          # sratio():
          il_sratio <- inv_link_sratio(x_test, link = link)
          il_sratio_ch <- inv_link_sratio_ch(x_test, link = link)
          expect_equivalent(il_sratio, il_sratio_ch)
          expect_equal(dim(il_sratio), c(ndraws, nobsv, ncat))

          # cratio():
          il_cratio <- inv_link_cratio(nx_test, link = link)
          il_cratio_ch <- inv_link_cratio_ch(nx_test, link = link)
          expect_equivalent(il_cratio, il_cratio_ch)
          expect_equal(dim(il_cratio), c(ndraws, nobsv, ncat))

          # acat():
          il_acat <- inv_link_acat(nx_test, link = link)
          il_acat_ch <- inv_link_acat_ch(nx_test, link = link)
          expect_equivalent(il_acat, il_acat_ch)
          expect_equal(dim(il_acat), c(ndraws, nobsv, ncat))
        }
      }
    }
  }
})

test_that("link_<ordinal_family>() works correctly for arrays", {
  source(testthat::test_path(file.path("helpers", "link_ordinal_ch.R")))
  source(testthat::test_path(file.path("helpers", "simopts_catlike.R")))
  for (ndraws in ndraws_vec) {
    for (nobsv in nobsv_vec) {
      for (ncat in ncat_vec) {
        x_test <- array(rdirichlet(ndraws * nobsv, alpha = rep(1, ncat)),
                        dim = c(ndraws, nobsv, ncat))
        for (link in c("logit", "probit", "cauchit", "cloglog")) {
          # cumulative():
          l_cumul <- link_cumulative(x_test, link = link)
          l_cumul_ch <- link_cumulative_ch(x_test, link = link)
          expect_equivalent(l_cumul, l_cumul_ch)
          expect_equal(dim(l_cumul), c(ndraws, nobsv, ncat - 1))

          # sratio():
          l_sratio <- link_sratio(x_test, link = link)
          l_sratio_ch <- link_sratio_ch(x_test, link = link)
          expect_equivalent(l_sratio, l_sratio_ch)
          expect_equal(dim(l_sratio), c(ndraws, nobsv, ncat - 1))

          # cratio():
          l_cratio <- link_cratio(x_test, link = link)
          l_cratio_ch <- link_cratio_ch(x_test, link = link)
          expect_equivalent(l_cratio, l_cratio_ch)
          expect_equal(dim(l_cratio), c(ndraws, nobsv, ncat - 1))

          # acat():
          l_acat <- link_acat(x_test, link = link)
          l_acat_ch <- link_acat_ch(x_test, link = link)
          expect_equivalent(l_acat, l_acat_ch)
          expect_equal(dim(l_acat), c(ndraws, nobsv, ncat - 1))
        }
      }
    }
  }
})

test_that("inv_link_<ordinal_family>() inverts link_<ordinal_family>()", {
  source(testthat::test_path(file.path("helpers", "simopts_catlike.R")))
  for (ndraws in ndraws_vec) {
    for (nobsv in nobsv_vec) {
      for (ncat in ncat_vec) {
        x_test <- array(rdirichlet(ndraws * nobsv, alpha = rep(1, ncat)),
                        dim = c(ndraws, nobsv, ncat))
        for (link in c("logit", "probit", "cauchit", "cloglog")) {
          # cumulative():
          l_cumul <- link_cumulative(x_test, link = link)
          il_cumul <- inv_link_cumulative(l_cumul, link = link)
          expect_equivalent(il_cumul, x_test)

          # sratio():
          l_sratio <- link_sratio(x_test, link = link)
          il_sratio <- inv_link_sratio(l_sratio, link = link)
          expect_equivalent(il_sratio, x_test)

          # cratio():
          l_cratio <- link_cratio(x_test, link = link)
          il_cratio <- inv_link_cratio(l_cratio, link = link)
          expect_equivalent(il_cratio, x_test)

          # acat():
          l_acat <- link_acat(x_test, link = link)
          il_acat <- inv_link_acat(l_acat, link = link)
          expect_equivalent(il_acat, x_test)
        }
      }
    }
  }
})

test_that("link_<ordinal_family>() inverts inv_link_<ordinal_family>()", {
  source(testthat::test_path(file.path("helpers", "simopts_catlike.R")))
  for (ndraws in ndraws_vec) {
    for (nobsv in nobsv_vec) {
      for (ncat in ncat_vec) {
        x_test <- array(rnorm(ndraws * nobsv * (ncat - 1)),
                        dim = c(ndraws, nobsv, ncat - 1))
        nx_test <- -x_test
        for (link in c("logit", "probit", "cauchit", "cloglog")) {
          # cumulative():
          il_cumul <- inv_link_cumulative(x_test, link = link)
          l_cumul <- link_cumulative(il_cumul, link = link)
          expect_equivalent(l_cumul, x_test)

          # sratio():
          il_sratio <- inv_link_sratio(x_test, link = link)
          l_sratio <- link_sratio(il_sratio, link = link)
          expect_equivalent(l_sratio, x_test)

          # cratio():
          il_cratio <- inv_link_cratio(x_test, link = link)
          l_cratio <- link_cratio(il_cratio, link = link)
          expect_equivalent(l_cratio, x_test)

          # acat():
          il_acat <- inv_link_acat(x_test, link = link)
          l_acat <- link_acat(il_acat, link = link)
          expect_equivalent(l_acat, x_test)
        }
      }
    }
  }
})

test_that(paste(
  "dsratio() and dcratio() give the same results for symmetric distribution",
  "functions"
), {
  source(testthat::test_path(file.path("helpers", "simopts_catlike_oneobs.R")))
  for (ndraws in ndraws_vec) {
    for (ncat in ncat_vec) {
      thres_test <- matrix(rnorm(ndraws * (ncat - 1)), nrow = ndraws)
      # Emulate no category-specific effects (i.e., only a single vector of
      # linear predictors) as well as category-specific effects (i.e., a matrix
      # of linear predictors):
      eta_test_list <- list(
        rnorm(ndraws),
        matrix(rnorm(ndraws * (ncat - 1)), nrow = ndraws)
      )
      for (eta_test in eta_test_list) {
        for (link in c("logit", "probit", "cauchit", "cloglog")) {
          d_sratio <- dsratio(seq_len(ncat),
                              eta_test, thres_test, link = link)
          d_cratio <- dcratio(seq_len(ncat),
                              eta_test, thres_test, link = link)
          if (link != "cloglog") {
            expect_equal(d_sratio, d_cratio)
          } else {
            expect_false(isTRUE(all.equal(d_sratio, d_cratio)))
          }
        }
      }
    }
  }
})

test_that(paste(
  "inv_link_sratio() and inv_link_cratio() applied to arrays give the same",
  "results for symmetric distribution functions (when respecting the sign",
  "appropriately)."
), {
  source(testthat::test_path(file.path("helpers", "simopts_catlike.R")))
  for (ndraws in ndraws_vec) {
    for (nobsv in nobsv_vec) {
      for (ncat in ncat_vec) {
        x_test <- array(rnorm(ndraws * nobsv * (ncat - 1)),
                        dim = c(ndraws, nobsv, ncat - 1))
        nx_test <- -x_test
        for (link in c("logit", "probit", "cauchit", "cloglog")) {
          il_sratio <- inv_link_sratio(x_test, link = link)
          il_cratio <- inv_link_cratio(nx_test, link = link)
          if (link != "cloglog") {
            expect_equal(il_sratio, il_cratio)
          } else {
            expect_false(isTRUE(all.equal(il_sratio, il_cratio)))
          }
        }
      }
    }
  }
})

test_that(paste(
  "link_sratio() and link_cratio() applied to arrays give the same",
  "results for symmetric distribution functions (when respecting the sign",
  "appropriately)."
), {
  source(testthat::test_path(file.path("helpers", "simopts_catlike.R")))
  for (ndraws in ndraws_vec) {
    for (nobsv in nobsv_vec) {
      for (ncat in ncat_vec) {
        x_test <- array(rdirichlet(ndraws * nobsv, alpha = rep(1, ncat)),
                        dim = c(ndraws, nobsv, ncat))
        for (link in c("logit", "probit", "cauchit", "cloglog")) {
          l_sratio <- link_sratio(x_test, link = link)
          l_cratio <- link_cratio(x_test, link = link)
          if (link != "cloglog") {
            expect_equal(l_sratio, -l_cratio)
          } else {
            expect_false(isTRUE(all.equal(l_sratio, -l_cratio)))
          }
        }
      }
    }
  }
})

test_that("dcategorical() works correctly", {
  source(testthat::test_path(file.path("helpers", "inv_link_categorical_ch.R")))
  source(testthat::test_path(file.path("helpers", "simopts_catlike_oneobs.R")))
  for (ndraws in ndraws_vec) {
    for (ncat in ncat_vec) {
      eta_test_list <- list(cbind(
        0, matrix(rnorm(ndraws * (ncat - 1)), nrow = ndraws)
      ))
      if (ndraws == 1) {
        eta_test_list <- c(eta_test_list, list(c(0, rnorm(ncat - 1))))
      }
      for (eta_test in eta_test_list) {
        d_categorical <- dcategorical(seq_len(ncat), eta_test)
        d_categorical_ch <- inv_link_categorical_ch(eta_test,
                                                    refcat_ins = FALSE)
        expect_equivalent(d_categorical, d_categorical_ch)
        expect_equal(dim(d_categorical), c(ndraws, ncat))
      }
    }
  }
})

test_that("inv_link_categorical() works correctly for arrays", {
  source(testthat::test_path(file.path("helpers", "inv_link_categorical_ch.R")))
  source(testthat::test_path(file.path("helpers", "simopts_catlike.R")))
  for (ndraws in ndraws_vec) {
    for (nobsv in nobsv_vec) {
      for (ncat in ncat_vec) {
        x_test <- array(rnorm(ndraws * nobsv * (ncat - 1)),
                        dim = c(ndraws, nobsv, ncat - 1))
        il_categorical <- inv_link_categorical(x_test)
        il_categorical_ch <- inv_link_categorical_ch(x_test)
        expect_equivalent(il_categorical, il_categorical_ch)
        expect_equal(dim(il_categorical), c(ndraws, nobsv, ncat))
      }
    }
  }
})

test_that("link_categorical() works correctly for arrays", {
  source(testthat::test_path(file.path("helpers", "link_categorical_ch.R")))
  source(testthat::test_path(file.path("helpers", "simopts_catlike.R")))
  for (ndraws in ndraws_vec) {
    for (nobsv in nobsv_vec) {
      for (ncat in ncat_vec) {
        x_test <- array(rdirichlet(ndraws * nobsv, alpha = rep(1, ncat)),
                        dim = c(ndraws, nobsv, ncat))
        l_categorical <- link_categorical(x_test)
        l_categorical_ch <- link_categorical_ch(x_test)
        expect_equivalent(l_categorical, l_categorical_ch)
        expect_equal(dim(l_categorical), c(ndraws, nobsv, ncat - 1))
      }
    }
  }
})

test_that("inv_link_categorical() inverts link_categorical()", {
  source(testthat::test_path(file.path("helpers", "simopts_catlike.R")))
  for (ndraws in ndraws_vec) {
    for (nobsv in nobsv_vec) {
      for (ncat in ncat_vec) {
        x_test <- array(rdirichlet(ndraws * nobsv, alpha = rep(1, ncat)),
                        dim = c(ndraws, nobsv, ncat))
        l_categorical <- link_categorical(x_test)
        il_categorical <- inv_link_categorical(l_categorical)
        expect_equivalent(il_categorical, x_test)
      }
    }
  }
})

test_that("link_categorical() inverts inv_link_categorical()", {
  source(testthat::test_path(file.path("helpers", "simopts_catlike.R")))
  for (ndraws in ndraws_vec) {
    for (nobsv in nobsv_vec) {
      for (ncat in ncat_vec) {
        x_test <- array(rnorm(ndraws * nobsv * (ncat - 1)),
                        dim = c(ndraws, nobsv, ncat - 1))
        il_categorical <- inv_link_categorical(x_test)
        l_categorical <- link_categorical(il_categorical)
        expect_equivalent(l_categorical, x_test)
      }
    }
  }
})
