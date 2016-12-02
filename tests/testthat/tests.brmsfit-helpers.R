test_that("first_greater returns expected results", {
  set.seed(12345)
  A <- cbind(1:10, 11:20, 21:30)
  x <- sample(1:25, size = 10, replace = TRUE)
  expect_equal(first_greater(A, x), c(3,3,3,3,2,1,2,2,2,3))
  expect_equal(first_greater(A, x, i = 2), c(3,3,3,3,2,2,2,2,2,3))
})

test_that("array2list performs correct conversion", {
  A <- array(1:27, dim = c(3,3,3))
  B <- list(matrix(1:9,3,3), matrix(10:18,3,3), matrix(19:27,3,3))
  expect_equal(array2list(A), B)
})

test_that("ilink(x, probit) and ilink(x, probit_approx) produce similar results", {
  all.equal(ilink(-10:10, "probit"), ilink(-10:10, "probit_approx"), tolerance = 1e-3)
})

test_that("get_cornames returns desired correlation names", {
  names <- c("Intercept", "x", "y")
  expect_equal(get_cornames(names), c("cor(Intercept,x)", "cor(Intercept,y)", "cor(x,y)"))
  expect_equal(get_cornames(names, brackets = FALSE), 
               c("cor__Intercept__x", "cor__Intercept__y", "cor__x__y"))
  expect_equal(get_cornames(names, brackets = FALSE, sep = "_"), 
               c("cor_Intercept_x", "cor_Intercept_y", "cor_x_y"))
  expect_equal(get_cornames(names, type = "rescor"),
               c("rescor(Intercept,x)", "rescor(Intercept,y)", "rescor(x,y)"))
})

test_that("get_cov_matrix returns appropriate dimensions", {
  sd <- cbind(1:10,11:20); cor <- cbind(seq(-0.5, 0.4, 0.1))
  expect_equal(dim(get_cov_matrix(sd = sd, cor = cor)$cov), c(10,2,2))
  expect_equal(dim(get_cov_matrix(sd = sd)$cov), c(10,2,2))
})

test_that("ARMA covariance matrices are computed correctly", {
  ar <- 0.5
  ma <- 0.3
  sigma <- 2
  se2 <- 1:4
  # test for AR1 cov matrix
  ar_mat <- get_cov_matrix_ar1(ar = matrix(ar), sigma = matrix(sigma), 
                               se2 = se2, nrows = length(se2))
  expected_ar_mat <- sigma^2 / (1 - ar^2) * 
                     cbind(c(1, ar, ar^2, ar^3),
                           c(ar, 1, ar, ar^2),
                           c(ar^2, ar, 1, ar),
                           c(ar^3, ar^2, ar, 1))
  expected_ar_mat <- expected_ar_mat + diag(se2)
  expect_equal(ar_mat[1, , ], expected_ar_mat)
  # test for MA1 cov matrix
  ma_mat <- get_cov_matrix_ma1(ma = matrix(ma), sigma = matrix(sigma), 
                               se2 = se2, nrows = length(se2))
  expected_ma_mat <- sigma^2 * 
                     cbind(c(1+ma^2, ma, 0, 0),
                           c(ma, 1+ma^2, ma, 0),
                           c(0, ma, 1+ma^2, ma),
                           c(0, 0, ma, 1+ma^2))
  expected_ma_mat <- expected_ma_mat + diag(se2)
  expect_equal(ma_mat[1, , ], expected_ma_mat)
  # test for ARMA1 cov matrix
  arma_mat <- get_cov_matrix_arma1(ar = matrix(ar), ma = matrix(ma), 
                                 sigma = matrix(sigma), 
                                 se2 = se2, nrows = length(se2))
  g0 <- 1 + ma^2 + 2 * ar * ma
  g1 <- (1 + ar * ma) * (ar + ma)
  expected_arma_mat <- sigma^2 / (1 - ar^2) * 
                       cbind(c(g0, g1, g1 * ar, g1 * ar^2),
                             c(g1, g0, g1, g1 * ar),
                             c(g1 * ar, g1, g0, g1),
                             c(g1 * ar^2, g1 * ar, g1, g0))
  expected_arma_mat <- expected_arma_mat + diag(se2)
  expect_equal(arma_mat[1, , ], expected_arma_mat)
  # test for identity matrix
  ident_mat <- get_cov_matrix_ident(sigma = matrix(sigma), 
                                    se2 = se2, nrows = length(se2))
  expected_ident_mat <- diag(sigma^2 + se2)
  expect_equal(ident_mat[1, , ], expected_ident_mat)
})

test_that("evidence_ratio returns expected results", {
  post_samples <- c(-4:10); prior_samples <- c(-2:12)
  expect_true(evidence_ratio(x = post_samples, prior_samples = prior_samples) > 1)
  expect_equal(evidence_ratio(x = post_samples, cut = 0.5, wsign = "greater"), 10/5)
  expect_equal(evidence_ratio(x = post_samples, cut = 0.5, wsign = "less"), 5/10)
})

test_that("expand_matrix returns expected results", {
  A <- matrix(1:6, 3, 2); x <- c(1,2,1)
  expect_equivalent(as.matrix(expand_matrix(A, x)), 
                    rbind(c(1, 4, 0, 0), c(0, 0, 2, 5), c(3, 6, 0, 0)))
})

test_that("find_names finds all valid variable names in a string", {
  expect_equal(find_names("x + b.x - .5 + abc(a__3) : 1/2 - 0.2"), c("x", "b.x", "a__3"))
})

test_that("get_summary returns correct dims and names", {
  col_names <- c("Estimate", "Est.Error", "5%ile", "50%ile", "95%ile")
  samples_2dim <- matrix(rnorm(1000), nrow = 100, ncol = 10)
  summary_2dim <- get_summary(samples_2dim, probs = c(0.05, 0.5, 0.95))
  expect_equal(dim(summary_2dim), c(10, 5))
  expect_equal(dimnames(summary_2dim), 
               list(as.character(1:10), col_names))
  
  col_names <- c("Estimate", "Est.Error", "2.5%ile", "97.5%ile")
  samples_3dim <- array(sample(1:10, x = 4000, replace = TRUE), 
                        dim = c(100, 10, 4))
  summary_3dim <- get_summary(samples_3dim)
  expect_equal(dim(summary_3dim), c(10, 4, 4))
  expect_equal(dimnames(summary_3dim), 
               list(as.character(1:10), col_names, 
                    paste0("P(Y = ", 1:4, ")")))
  expect_error(get_summary(rnorm(100)))
})

test_that("get_table returns correct dims and names", {
  samples <- matrix(sample(1:5, 1000, replace = TRUE), 
                    nrow = 100, ncol = 10)
  res_table <- get_table(samples)
  expect_equal(dim(res_table), c(10, 5))
  expect_equal(dimnames(res_table), 
               list(as.character(1:10), paste0("P(Y = ", 1:5, ")")))
  expect_equal(sum(res_table[5, ]), 1)
})

test_that("evidence_ratio runs without errors", {
  post <- rnorm(100, mean = 1)
  prior <- rnorm(100, sd = 2)
  expect_true(is.na(evidence_ratio(post, wsign = "equal")))
  expect_true(is.numeric(evidence_ratio(post, wsign = "equal",
                                        prior_samples = prior)))
  expect_true(is.numeric(evidence_ratio(post, wsign = "greater")))
  expect_true(is.numeric(evidence_ratio(post, wsign = "less")))
})

test_that("get_sigma correctly extract residual SDs", {
  nsamples <- 40
  sigma <- matrix(rexp(nsamples, 1))
  expect_equal(length(get_sigma(sigma, data = list(), i = 2)), 
               nsamples)
  expect_equal(length(get_sigma(sigma, data = list(se = 2:11), i = 3)), 
               nsamples)
  expect_equal(dim(get_sigma(NULL, data = list(se = 2:11, N = 10), 
                             dim = c(5, 10))), c(5, 10))
  expect_equal(get_sigma(NULL, data = list(sigma = 2:11), i = 5), 6)
})

test_that("arma_predictor runs without errors", {
  ns <- 30
  nobs <- 18
  data <- list(Y = rnorm(nobs), tgroup = rep(1:3, each = 6))
  ar <- matrix(rnorm(ns * nobs), nrow = ns, ncol = nobs)
  ma <- matrix(rnorm(ns * nobs), nrow = ns, ncol = nobs)
  eta <- matrix(rnorm(ns * nobs), nrow = ns, ncol = nobs)
  expect_equal(arma_predictor(standata = data, eta = eta), eta)
  expect_silent(arma_predictor(standata = data, eta = eta, ar = ar))
  expect_silent(arma_predictor(standata = data, eta = eta, ma = ma))
  expect_silent(arma_predictor(standata = data, eta = eta, ar = ar, ma = ma))
})

test_that("cs_predictor runs without errors", {
  X <- matrix(rnorm(300), nrow = 100, ncol = 3)
  b <- matrix(rnorm(30 * 9), nrow = 30, ncol = 9)
  eta <- matrix(rnorm(3000), nrow = 30, ncol = 100)
  expect_equal(dim(cs_predictor(X = X, b = b, eta = eta, ncat = 4)),
               c(30, 100, 3))
  rc <- brms:::named_list(1:4, list(g = eta))
  expect_equal(dim(cs_predictor(X = X, b = b, eta = eta, 
                                 ncat = 4, r = rc)),
               c(30, 100, 3))
})

test_that("extract_pars returns correct parameter names", {
  all_pars <- c("ab", "ba", "bac")
  expect_equal(extract_pars("^b", all_pars = all_pars),
               c("ba", "bac"))
  expect_equal(extract_pars(c("ab", "ba$", "c$"), all_pars = all_pars),
               all_pars)
  expect_equal(extract_pars(c("ab", "ba", "cd"), all_pars = all_pars, 
                            exact_match = TRUE), 
               c("ab", "ba"))
  expect_equal(extract_pars(NA, all_pars = all_pars, na_value = NA),
               NA)
})