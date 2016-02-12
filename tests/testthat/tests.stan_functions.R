test_that("self-defined Stan functions work correctly", {
  skip_on_cran()
  rstan::expose_stan_functions(new_stan_functions)
  
  # ARMA matrix generating functions
  cov_ar1_R <- get_cov_matrix_ar1(ar = matrix(0.5), sigma = 2, 
                                  sq_se = 0, nrows = 3)[1, , ]
  expect_equal(cov_matrix_ar1(0.5, 2, 3), cov_ar1_R)
  cov_ma1_R <- matrix(get_cov_matrix_ma1(ma = matrix(-0.3), sigma = 3, 
                                         sq_se = 0, nrows = 1)[1, , ])
  expect_equal(cov_matrix_ma1(-0.3, 3, 1), cov_ma1_R)
  cov_arma1_R <- get_cov_matrix_arma1(ar = matrix(-0.5), ma = matrix(0.7), 
                                      sigma = 4, sq_se = 0, nrows = 5)[1, , ]
  expect_equal(cov_matrix_arma1(-0.5, 0.7, 4, 5), cov_arma1_R)
  
  # log-likelihood functions for covariance models
  y <- rnorm(9)
  eta <- rnorm(9)
  ll_stan <- normal_cov_log(y, eta = eta, squared_se = 1:9, 
                            N_tg = 2, begin = c(1, 5), end = c(4, 9),
                            nrows = c(4, 5), res_cov_matrix = cov_arma1_R)
  ll_R <- c(dmulti_normal(y[1:4], eta[1:4], cov_arma1_R[1:4, 1:4] + diag(1:4)),
            dmulti_normal(y[5:9], eta[5:9], cov_arma1_R[1:5, 1:5] + diag(5:9)))
  expect_equal(ll_stan, sum(ll_R))
  ll_stan <- student_t_cov_log(y, nu = 10, eta = eta, squared_se = 1:9, 
                               N_tg = 2, begin = c(1, 5), end = c(4, 9),
                               nrows = c(4, 5), res_cov_matrix = cov_arma1_R)
  ll_R <- c(dmulti_student(y[1:4], df = 10, mu = eta[1:4], 
                           Sigma = cov_arma1_R[1:4, 1:4] + diag(1:4)),
            dmulti_student(y[5:9], df = 10, mu = eta[5:9], 
                           Sigma = cov_arma1_R[1:5, 1:5] + diag(5:9)))
  expect_equal(ll_stan, sum(ll_R))
  
  # inverse gaussian functions
  shape <- rgamma(1, 20, 1)
  mu <- 20
  y <- statmod::rinvgauss(1, mean = mu, shape = shape)
  expect_equal(inv_gaussian_cdf_log(y, mu, shape, log(y), sqrt(y)),
               pinvgauss(y, mean = mu, shape = shape, log = TRUE))
  expect_equal(inv_gaussian_ccdf_log(y, mu, shape, log(y), sqrt(y)),
               log(1 - pinvgauss(y, mean = mu, shape = shape))) 
  expect_equal(inv_gaussian_log(y, mu, shape, log(y), sqrt(y)),
               dinvgauss(y, mean = mu, shape = shape, log = TRUE))
  mu <- 18:22
  y <- statmod::rinvgauss(5, mean = mu, shape = shape)
  expect_equal(inv_gaussian_vector_log(y, mu, shape, sum(log(y)), sqrt(y)),
               sum(dinvgauss(y, mean = mu, shape = shape, log = TRUE)))
  
  # zero-inflated and hurdle log-densities
  dat <- list(Y = c(0, 10), N_trait = 2, max_obs = 15)
  dat2 <- list(Y = c(0, 0.5), N_trait = 2)
  samp <- list(eta = matrix(rnorm(4), ncol = 4), shape = 2, phi = 2)
  for (i in seq_along(dat$Y)) {
    # zero-inflated
    args <- list(y = dat$Y[i], eta = samp$eta[i], eta_zi = samp$eta[i+2])
    expect_equal(do.call(zero_inflated_poisson_log, args),
                 loglik_zero_inflated_poisson(i, dat, samp))
    expect_equal(do.call(zero_inflated_neg_binomial_2_log, 
                         c(args, shape = samp$shape)),
                 loglik_zero_inflated_negbinomial(i, dat, samp))
    expect_equal(do.call(zero_inflated_binomial_log, 
                         c(args, trials = dat$max_obs)),
                 loglik_zero_inflated_binomial(i, dat, samp))
    # zero_inflated_beta requires Y to be in (0,1)
    args <- list(y = dat2$Y[i], eta = samp$eta[i], eta_zi = samp$eta[i+2])
    expect_equal(do.call(zero_inflated_beta_log, c(args, phi = samp$phi)),
                 loglik_zero_inflated_beta(i, dat2, samp))
    # hurdle
    args <- list(y = dat$Y[i], eta = samp$eta[i], eta_hu = samp$eta[i+2])
    expect_equal(do.call(hurdle_poisson_log, args),
                 loglik_hurdle_poisson(i, dat, samp))
    expect_equal(do.call(hurdle_neg_binomial_2_log, 
                         c(args, shape = samp$shape)),
                 loglik_hurdle_negbinomial(i, dat, samp))
    expect_equal(do.call(hurdle_gamma_log, 
                         c(args, shape = samp$shape)),
                 loglik_hurdle_gamma(i, dat, samp))
  }
  
  # ordinal log-densities
  dat <- list(Y = 2, max_obs = 4)
  eta <- rnorm(1)
  etap <- array(rnorm(6), dim = c(2, 1, 3))
  thres <- sort(rnorm(3))
  # cumulative and sratio require thres - eta
  samp <- list(eta = rep(thres, each = 2) - array(eta, dim = c(2, 1, 3)))
  expect_equal(cumulative_log(dat$Y, eta, thres),
               loglik_cumulative(1, dat, samp, link = "probit")[1])
  expect_equal(sratio_log(dat$Y, eta, thres),
               loglik_sratio(1, dat, samp, link = "logit")[1])
  # acat and cratio require eta - thres
  # also category specific effects are included here
  samp <- list(eta = eta + etap - rep(thres, each = 2))
  expect_equal(cratio_log(dat$Y, eta, etap[1, , ], thres),
               loglik_cratio(1, dat, samp, link = "cloglog")[1])
  expect_equal(acat_log(dat$Y, eta, etap[1, , ], thres),
               loglik_acat(1, dat, samp, link = "cauchit")[1])
 
  # kronecker product
  A <- matrix(c(3, 2, 1, 2, 4, 1, 1, 1, 5), nrow = 3)
  B <- matrix(c(3, 2, 2, 4), nrow = 2)
  sd <- c(2, 7)
  expect_equal(t(kronecker_cholesky(A, t(chol(B)), sd)),
               chol(kronecker(A, diag(sd) %*% B %*% diag(sd))))
  
  # cauchit link
  expect_equal(inv_cauchit(1.5), pcauchy(1.5)) 
})

