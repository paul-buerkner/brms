context("Tests for self-defined Stan functions")

test_that("self-defined Stan functions work correctly", {
  # for some reason expose_stan_functions doesn't work within R CMD CHECK
  skip_if_not(exists("new_stan_functions", asNamespace("brms")))
  rstan::expose_stan_functions(brms:::new_stan_functions)

  # ARMA matrix generating functions
  cov_ar1_R <- get_cov_matrix_ar1(ar = matrix(0.5), sigma = 2,
                                  se2 = 0, nrows = 3)[1, , ]
  expect_equal(cov_matrix_ar1(0.5, 2, 3), cov_ar1_R)
  cov_ma1_R <- matrix(get_cov_matrix_ma1(ma = matrix(-0.3), sigma = 3,
                                         se2 = 0, nrows = 1)[1, , ])
  expect_equal(cov_matrix_ma1(-0.3, 3, 1), cov_ma1_R)
  cov_arma1_R <- get_cov_matrix_arma1(ar = matrix(-0.5), ma = matrix(0.7),
                                      sigma = 4, se2 = 0, nrows = 5)[1, , ]
  expect_equal(cov_matrix_arma1(-0.5, 0.7, 4, 5), cov_arma1_R)

  # log-likelihood functions for covariance models
  y <- rnorm(9)
  eta <- rnorm(9)
  ll_stan <- normal_cov_lpdf(y, eta = eta, se2 = 1:9, I = 2,
                            begin = c(1, 5), end = c(4, 9), nobs = c(4, 5),
                            res_cov_matrix = cov_arma1_R)
  ll_R <- c(dmulti_normal(y[1:4], eta[1:4], cov_arma1_R[1:4, 1:4] + diag(1:4)),
            dmulti_normal(y[5:9], eta[5:9], cov_arma1_R[1:5, 1:5] + diag(5:9)))
  expect_equal(ll_stan, sum(ll_R))
  ll_stan <- student_t_cov_lpdf(y, nu = 10, eta = eta, se2 = 1:9, I = 2,
                               begin = c(1, 5), end = c(4, 9), nobs = c(4, 5),
                               res_cov_matrix = cov_arma1_R)
  ll_R <- c(dmulti_student(y[1:4], df = 10, mu = eta[1:4],
                           Sigma = cov_arma1_R[1:4, 1:4] + diag(1:4)),
            dmulti_student(y[5:9], df = 10, mu = eta[5:9],
                           Sigma = cov_arma1_R[1:5, 1:5] + diag(5:9)))
  expect_equal(ll_stan, sum(ll_R))

  # inverse gaussian functions
  shape <- rgamma(1, 20, 1)
  mu <- 20
  y <- statmod::rinvgauss(1, mean = mu, shape = shape)
  expect_equal(inv_gaussian_lpdf(y, mu, shape, log(y), sqrt(y)),
               dinvgauss(y, mean = mu, shape = shape, log = TRUE))
  expect_equal(inv_gaussian_lcdf(y, mu, shape, log(y), sqrt(y)),
               pinvgauss(y, mean = mu, shape = shape, log = TRUE))
  expect_equal(inv_gaussian_lccdf(y, mu, shape, log(y), sqrt(y)),
               log(1 - pinvgauss(y, mean = mu, shape = shape)))

  mu <- 18:22
  y <- statmod::rinvgauss(5, mean = mu, shape = shape)
  expect_equal(inv_gaussian_vector_lpdf(y, mu, shape, sum(log(y)), sqrt(y)),
               sum(dinvgauss(y, mean = mu, shape = shape, log = TRUE)))

  # exgaussian functions
  beta <- rgamma(1, 1, 0.1)
  sigma <- rgamma(1, 10, 0.1)
  mu <- 10
  y <- rexgaussian(1, mu = mu, sigma = sigma, beta = beta)
  expect_equal(exgaussian_lpdf(y, mu, sigma, beta),
               dexgaussian(y, mu, sigma, beta, log = TRUE))
  expect_equal(exgaussian_lcdf(y, mu, sigma, beta),
               pexgaussian(y, mu, sigma, beta, log = TRUE))
  expect_equal(exgaussian_lccdf(y, mu, sigma, beta),
               pexgaussian(y, mu, sigma, beta,
                        lower.tail = FALSE, log = TRUE))

  # asym_laplace functions
  mu <- 10
  quantile <- rbeta(1, 2, 1)
  sigma <- rgamma(1, 10, 0.1)
  y <- rasym_laplace(1, mu = mu, sigma = sigma, quantile = quantile)
  expect_equal(asym_laplace_lpdf(y, mu, sigma, quantile),
               dasym_laplace(y, mu, sigma, quantile, log = TRUE))
  expect_equal(asym_laplace_lcdf(y, mu, sigma, quantile),
               pasym_laplace(y, mu, sigma, quantile, log = TRUE))
  expect_equal(asym_laplace_lccdf(y, mu, sigma, quantile),
               pasym_laplace(y, mu, sigma, quantile,
                             lower.tail = FALSE, log = TRUE))

  # wiener diffusion model functions
  alpha = 2
  tau = 0.5
  beta = 0.5
  delta = 0.5
  y <- rWiener(1, alpha, tau, beta, delta)
  y$resp <- ifelse(y$resp == "lower", 0, 1)
  expect_equal(wiener_diffusion_lpdf(y$q, y$resp, alpha, tau, beta, delta),
               dWiener(y$q, alpha, tau, beta, delta,
                       resp = y$resp, log = TRUE))

  # zero-inflated and hurdle log-densities
  draws <- draws2 <- list(eta = matrix(rnorm(4), ncol = 4),
                          shape = 2, phi = 2, sigma = 2)
  draws$data <- list(Y = c(0, 10), N_trait = 2, max_obs = 15)
  draws2$data <- list(Y = c(0, 0.5), N_trait = 2)
  for (i in seq_along(draws$data$Y)) {
    eta_zi_args <- list(y = draws$data$Y[i], eta = draws$eta[i],
                        eta_zi = draws$eta[i+2])
    zi_args <- list(y = draws$data$Y[i], eta = draws$eta[i],
                    zi = inv_logit(eta_zi_args$eta_zi))
    eta_hu_args <- list(y = draws$data$Y[i], eta = draws$eta[i],
                    eta_hu = draws$eta[i+2])
    hu_args <- list(y = draws$data$Y[i], eta = draws$eta[i],
                    hu = inv_logit(eta_hu_args$eta_hu))
    draws$f$link <- "log"

    expect_equal(do.call(zero_inflated_poisson_lpmf, zi_args),
                 loglik_zero_inflated_poisson(i, draws))
    expect_equal(do.call(zero_inflated_poisson_logit_lpmf, eta_zi_args),
                 loglik_zero_inflated_poisson(i, draws))

    expect_equal(do.call(zero_inflated_neg_binomial_lpmf,
                         c(zi_args, shape = draws$shape)),
                 loglik_zero_inflated_negbinomial(i, draws))
    expect_equal(do.call(zero_inflated_neg_binomial_logit_lpmf,
                         c(eta_zi_args, shape = draws$shape)),
                 loglik_zero_inflated_negbinomial(i, draws))

    expect_equal(do.call(hurdle_poisson_lpmf, hu_args),
                 loglik_hurdle_poisson(i, draws))
    expect_equal(do.call(hurdle_poisson_logit_lpmf, eta_hu_args),
                 loglik_hurdle_poisson(i, draws))

    expect_equal(do.call(hurdle_neg_binomial_lpmf,
                         c(hu_args, shape = draws$shape)),
                 loglik_hurdle_negbinomial(i, draws))
    expect_equal(do.call(hurdle_neg_binomial_logit_lpmf,
                         c(eta_hu_args, shape = draws$shape)),
                 loglik_hurdle_negbinomial(i, draws))

    expect_equal(do.call(hurdle_gamma_lpdf,
                         c(hu_args, shape = draws$shape)),
                 loglik_hurdle_gamma(i, draws))
    expect_equal(do.call(hurdle_gamma_logit_lpdf,
                         c(eta_hu_args, shape = draws$shape)),
                 loglik_hurdle_gamma(i, draws))

    draws$f$link <- "identity"
    expect_equal(do.call(hurdle_lognormal_lpdf,
                         c(hu_args, sigma = draws$sigma)),
                 loglik_hurdle_lognormal(i, draws))
    expect_equal(do.call(hurdle_lognormal_logit_lpdf,
                         c(eta_hu_args, sigma = draws$sigma)),
                 loglik_hurdle_lognormal(i, draws))

    draws$f$link <- "logit"
    expect_equal(do.call(zero_inflated_binomial_lpmf,
                         c(zi_args, trials = draws$data$max_obs)),
                 loglik_zero_inflated_binomial(i, draws))
    expect_equal(do.call(zero_inflated_binomial_logit_lpmf,
                         c(eta_zi_args, trials = draws$data$max_obs)),
                 loglik_zero_inflated_binomial(i, draws))

    # zero_inflated_beta requires Y to be in (0,1)
    draws2$f$link <- "logit"
    eta_zi_args <- list(y = draws2$data$Y[i], eta = draws$eta[i],
                        eta_zi = draws$eta[i+2])
    zi_args <- list(y = draws2$data$Y[i], eta = draws$eta[i],
                    zi = inv_logit(eta_zi_args$eta_zi))
    expect_equal(do.call(zero_inflated_beta_lpdf,
                         c(zi_args, phi = draws$phi)),
                 loglik_zero_inflated_beta(i, draws2))
    expect_equal(do.call(zero_inflated_beta_logit_lpdf,
                         c(eta_zi_args, phi = draws$phi)),
                 loglik_zero_inflated_beta(i, draws2))
  }

  # ordinal log-densities
  eta <- rnorm(1)
  etap <- array(rnorm(6), dim = c(2, 1, 3))
  thres <- sort(rnorm(3))
  # cumulative and sratio require thres - eta
  draws <- list(eta = rep(thres, each = 2) - array(eta, dim = c(2, 1, 3)))
  draws$data <- list(Y = 2, max_obs = 4)
  draws$f$link <- "probit"
  expect_equal(cumulative_lpmf(draws$data$Y, eta, thres),
               loglik_cumulative(1, draws)[1])
  draws$f$link <- "logit"
  expect_equal(sratio_lpmf(draws$data$Y, eta, thres),
               loglik_sratio(1, draws)[1])
  # acat and cratio require eta - thres
  # also category specific effects are included here
  draws$eta <- eta + etap - rep(thres, each = 2)
  draws$f$link <- "cloglog"
  expect_equal(cratio_lpmf(draws$data$Y, eta, etap[1, , ], thres),
               loglik_cratio(1, draws)[1])
  draws$f$link <- "cauchit"
  expect_equal(acat_lpmf(draws$data$Y, eta, etap[1, , ], thres),
               loglik_acat(1, draws)[1])

  # kronecker product
  A <- matrix(c(3, 2, 1, 2, 4, 1, 1, 1, 5), nrow = 3)
  B <- matrix(c(3, 2, 2, 4), nrow = 2)
  sd <- c(2, 7)
  expect_equal(t(chol(base::kronecker(A, diag(sd) %*% B %*% diag(sd)))),
               kronecker(t(chol(A)), diag(sd) %*% t(chol(B))))

  # as_matrix
  expect_equal(as_matrix(1:28, 4, 7),
               rbind(1:7, 8:14, 15:21, 22:28))
  expect_equal(as_matrix(1:28, 3, 4), rbind(1:4, 5:8, 9:12))

  # cauchit and cloglog link
  expect_equal(inv_cauchit(1.5), pcauchy(1.5))
  expect_equal(cauchit(0.7), qcauchy(0.7))
  expect_equal(cloglog(0.2), link(0.2, "cloglog"))

  # monotonic
  # slightly arkward way to call this function to make sure
  # is doesn't conflict with the brms R function of the same name
  monotonic_temp <- get("monotonic", globalenv())
  expect_equal(monotonic_temp(1:10, 4), sum(1:4))
  expect_equal(monotonic_temp(rnorm(5), 0), 0)
})

