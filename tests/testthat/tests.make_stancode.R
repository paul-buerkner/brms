# simplifies manual calling of tests
expect_match2 <- brms:::expect_match2

test_that("specified priors appear in the Stan code", {
  dat <- data.frame(y = 1:10, x1 = rnorm(10), x2 = rnorm(10), 
                    g = rep(1:5, 2), h = factor(rep(1:5, each = 2)))
  prior <- c(prior(normal(0,1), coef = x1),
             prior(normal(0,2), coef = x2),
             prior(cauchy(0,1), sd, group = g),
             prior(cauchy(0,2), sd, group = g, coef = x1),
             prior(gamma(1, 1), class = sd, group = h))
  scode <- make_stancode(y ~ x1*x2 + (x1*x2|g) + (1 | h), dat,
                         prior = prior, sample_prior = TRUE)
  expect_match2(scode, "b[1] ~ normal(0, 1)")
  expect_match2(scode, "b[2] ~ normal(0, 2)")
  expect_match2(scode, "sd_1[1] ~ cauchy(0, 1)")
  expect_match2(scode, "sd_1[2] ~ cauchy(0, 2)")
  expect_match2(scode, "sigma ~ student_t(3, 0, 10)")
  expect_match2(scode, "sd_2 ~ gamma(1, 1)")
  expect_match2(scode, "prior_b_1 = normal_rng(0,1);")
  expect_match2(scode, "prior_sd_1_1 ~ cauchy(0,1)")
  expect_match2(scode, "prior_sd_2 ~ gamma(1,1)")
  
  prior <- c(prior(lkj(0.5), class = cor, group = g),
             prior(normal(0, 1), class = b),
             prior(normal(0, 5), class = Intercept),
             prior(cauchy(0, 5), class = sd))
  scode <- make_stancode(y ~ x1 + cs(x2) + (0 + x1 + x2 | g), 
                         dat, family = acat(), 
                         prior = prior, sample_prior = TRUE)
  expect_match2(scode, "b ~ normal(0, 1)")
  expect_match2(scode, "temp_Intercept ~ normal(0, 5)")
  expect_match2(scode, "sd_1 ~ cauchy(0, 5)")
  expect_match2(scode, "L_1 ~ lkj_corr_cholesky(0.5)")
  expect_match2(scode, "to_vector(bcs) ~ normal(0, 1)")
  expect_match2(scode, "prior_bcs = normal_rng(0,1)")
  expect_match2(scode, 
    "prior_b_Intercept = prior_temp_Intercept - dot_product(means_X, b)")
  
  prior <- c(prior(normal(0,5), nlpar = a),
             prior(normal(0,10), nlpar = b),
             prior(cauchy(0,1), class = sd, nlpar = a),
             prior(lkj(2), class = cor, group = g))
  scode <- make_stancode(bf(y ~ a * exp(-b * x1), 
                            a + b ~ (1|ID|g), nl = TRUE),
                         data = dat, prior = prior,
                         sample_prior = TRUE)
  expect_match2(scode, "b_a ~ normal(0, 5)")
  expect_match2(scode, "b_b ~ normal(0, 10)")
  expect_match2(scode, "sd_1[1] ~ cauchy(0, 1)")
  expect_match2(scode, "L_1 ~ lkj_corr_cholesky(2)")
  expect_match2(scode, "prior_b_a = normal_rng(0,5)")
  expect_match2(scode, "prior_sd_1_2 ~ student_t(3,0,10)")
  expect_match2(scode, "prior_cor_1 = lkj_corr_rng(2,2)[1, 2]")
  
  prior <- c(prior(lkj(2), rescor),
             prior(cauchy(0, 5), sigma),
             prior(cauchy(0, 1), sigma, coef = x1))
  scode <- make_stancode(cbind(y, x1) ~ x2, dat, prior = prior, 
                         sample_prior = TRUE)
  expect_match2(scode, "sigma[1] ~ cauchy(0, 5)")
  expect_match2(scode, "sigma[2] ~ cauchy(0, 1)")
  expect_match2(scode, "Lrescor ~ lkj_corr_cholesky(2)")
  expect_match2(scode, "prior_rescor = lkj_corr_rng(2,2)[1, 2]")
  
  prior <- c(prior(uniform(-1, 1), ar, lb = -0.7, ub = 0.5),
             prior(normal(0, 0.5), ma),
             prior(double_exponential(0, 1), arr),
             prior(normal(0, 5)))
  autocor <- cor_arma(p = 1, q = 2, r = 3)
  expect_warning(
    scode <- make_stancode(y ~ mo(g), dat, autocor = autocor,
                           prior = prior, sample_prior = TRUE),
    "Changing the boundaries of autocorrelation parameters"
  )
  expect_match2(scode, "vector<lower=-0.7,upper=0.5>[Kar] ar;")
  expect_match2(scode, "ar ~ uniform(-1, 1)")
  expect_match2(scode, "ma ~ normal(0, 0.5)")
  expect_match2(scode, "arr ~ double_exponential(0, 1)")
  expect_match2(scode, "bmo ~ normal(0, 5)")
  expect_match2(scode, "simplex_1 ~ dirichlet(con_simplex_1)")
  expect_match2(scode, "prior_simplex_1 = dirichlet_rng(con_simplex_1)")
  expect_match2(scode, "prior_ar ~ uniform(-1,1)")
  
  prior <- c(set_prior("target += normal_lpdf(b[1] | 0, 1)", check = FALSE),
             set_prior("", class = "sigma"))
  scode <- make_stancode(y ~ x1, dat, prior = prior,
                         sample_prior = TRUE)
  expect_match2(scode, "target += normal_lpdf(b[1] | 0, 1)")
  expect_true(!grepl("sigma ~", scode))
  
  prior <- prior(gamma(0, 1), coef = x1)
  expect_warning(make_stancode(y~x1, dat, prior = prior),
                 "no natural lower bound")
  prior <- prior(uniform(0,5), class = sd)
  expect_warning(make_stancode(y~x1 + (1|g), dat, prior = prior),
                  "no natural upper bound")
})

test_that("special shrinkage priors appear in the Stan code", {
  dat <- data.frame(y = 1:10, x1 = rnorm(10), x2 = rnorm(10))
  scode <- make_stancode(y ~ x1*x2, data = dat,
                         prior = prior(horseshoe(7, scale_global = 2)),
                         sample_prior = TRUE)
  expect_match2(scode,
    "  vector<lower=0>[Kc] hs_local; \n  real<lower=0> hs_global;")
  expect_match2(scode, 
    "  hs_local ~ student_t(7, 0, 1); \n  hs_global ~ cauchy(0, 2);")
  expect_match2(scode, "  b ~ normal(0, hs_local * hs_global);")
  
  scode <- make_stancode(y ~ x1*x2, data = dat,
                         prior = prior(lasso(2, scale = 10)),
                         sample_prior = TRUE)
  expect_match2(scode, "  lasso_inv_lambda ~ chi_square(2);")
  expect_match2(scode, "  b ~ double_exponential(0, 10 * lasso_inv_lambda);")
  expect_match2(scode,
    "  prior_b = double_exponential_rng(0,10*prior_lasso_inv_lambda);")
  
  expect_error(make_stancode(y ~ x1*x2, data = dat, 
                             prior = prior(horseshoe(-1))),
               "Degrees of freedom of the local priors")
  expect_error(make_stancode(y ~ x1*x2, data = dat, 
                             prior = prior(horseshoe(1, -1))),
               "Scale of the global prior")
  expect_error(make_stancode(y ~ x1*x2, data = dat, 
                             prior = prior(lasso(-1))),
               "Degrees of freedom of the shrinkage parameter prior")
})

test_that("link functions appear in the Stan code", {
  dat <- data.frame(y = 1:10, x = rnorm(10))
  expect_match2(make_stancode(y ~ x, dat, family = poisson()), 
               "Y ~ poisson_log(eta);")
  expect_match2(make_stancode(cbind(y, y + 1) ~ x, dat, family = gaussian("log")), 
               "eta_y[n] = exp(eta_y[n]);")
  expect_match2(make_stancode(y ~ x, dat, family = von_mises(tan_half)), 
               "eta[n] = inv_tan_half(eta[n]);")
  expect_match2(make_stancode(y ~ x, dat, family = weibull()),
                "eta[n] = exp((eta[n]) / shape);")
  expect_match2(make_stancode(y ~ x, dat, family = exponential("identity")),
               "eta[n] = inv(eta[n]);")
  expect_match2(make_stancode(y ~ x, dat, family = poisson("sqrt")),
               "eta[n] = square(eta[n]);")
  expect_match2(make_stancode(y ~ x, dat, family = bernoulli()),
                "Y ~ bernoulli_logit(eta);")
})

test_that("customized covariances appear in the Stan code", {
  scode <- make_stancode(rating ~ treat + period + carry + (1|subject), 
                         data = inhaler, cov_ranef = list(subject = 1))
  expect_match2(scode, "r_1_1 = sd_1[1] * (Lcov_1 * z_1[1])")
  
  scode <- make_stancode(rating ~ treat + period + carry + (1+carry|subject), 
                         data = inhaler, cov_ranef = list(subject = 1))
  expect_match2(scode,
    "kronecker(Lcov_1, diag_pre_multiply(sd_1, L_1)) * to_vector(z_1)")
  
  scode <- make_stancode(rating ~ treat + period + carry + (1+carry||subject), 
                         data = inhaler, cov_ranef = list(subject = 1))
  expect_match2(scode, 
               paste0("  r_1_1 = sd_1[1] * (Lcov_1 * z_1[1]); \n",
                      "  r_1_2 = sd_1[2] * (Lcov_1 * z_1[2]);"))
})

test_that("truncation appears in the Stan code", {
  expect_match2(make_stancode(time | trunc(0) ~ age + sex + disease,
                              data = kidney, family = "gamma"), 
               "Y[n] ~ gamma(shape, eta[n]) T[lb[n], ];")
  expect_match2(make_stancode(time | trunc(ub = 100) ~ age + sex + disease, 
                             data = kidney, family = student("log")), 
               "Y[n] ~ student_t(nu, eta[n], sigma) T[, ub[n]];")
  expect_match2(make_stancode(count | trunc(0, 150) ~ Trt_c, 
                             data = epilepsy, family = "poisson"), 
               "Y[n] ~ poisson(eta[n]) T[lb[n], ub[n]];")
})

test_that("make_stancode combines strings of multiple grouping factors", {
  expect_match2(make_stancode(count ~ (1|patient) + (1 + Trt_c | visit), 
                             data = epilepsy, family = "poisson"), 
               paste0("  vector[N] Z_1_1; \n",
                      "  // data for group-level effects of ID 2"))
  expect_match2(make_stancode(count ~ (1|visit) + (1+Trt_c|patient), 
                             data = epilepsy, family = "poisson"), 
               paste0("  int<lower=1> NC_1; \n",
                      "  // data for group-level effects of ID 2"))
})

test_that("make_stancode handles models without fixed effects", {
  expect_match2(make_stancode(count ~ 0 + (1|patient) + (1+Trt_c|visit), 
                             data = epilepsy, family = "poisson"), 
               "eta = rep_vector(0, N);")
})

test_that("make_stancode correctly restricts FE parameters", {
  data <- data.frame(y = rep(0:1, each = 5), x = rnorm(10))
  
  scode <- make_stancode(y ~ x, data, prior = set_prior("", lb = 2))
  expect_match2(scode, "vector<lower=2>[Kc] b")
  
  scode <- make_stancode(y ~ x, data, prior = set_prior("normal(0,2)", ub = "4"))
  expect_match2(scode, "vector<upper=4>[Kc] b")
  
  prior <- set_prior("normal(0,5)", lb = "-3", ub = 5)
  scode <- make_stancode(y ~ 0 + x, data, prior = prior)
  expect_match2(scode, "vector<lower=-3,upper=5>[K] b")
})

test_that("self-defined functions appear in the Stan code", {
  # cauchit link
  scode <- make_stancode(rating ~ treat, data = inhaler,
                         family = bernoulli("cauchit"))
  expect_match2(scode, "real inv_cauchit(real y)")
  
  # tan_half link
  expect_match2(make_stancode(rating ~ treat, data = inhaler,
                              family = von_mises("tan_half")),
               "real inv_tan_half(real y)")
  
  # logm1 link
  expect_match2(make_stancode(rating ~ treat, data = inhaler,
                              family = frechet()),
                "real expp1(real y)")
  
  # inverse gaussian models
  scode <- make_stancode(time | cens(censored) ~ age, data = kidney,
                                 family = inverse.gaussian)
  expect_match2(scode, "real inv_gaussian_lpdf(real y")
  expect_match2(scode, "real inv_gaussian_lcdf(real y")
  expect_match2(scode, "real inv_gaussian_lccdf(real y")
  expect_match2(scode, "real inv_gaussian_vector_lpdf(vector y")
  
  # von Mises models
  scode <- make_stancode(time ~ age, data = kidney, family = von_mises)
  expect_match2(scode, "real von_mises_real_lpdf(real y")
  expect_match2(scode, "real von_mises_vector_lpdf(vector y")
  
  # zero-inflated and hurdle models
  expect_match2(make_stancode(count ~ Trt_c, data = epilepsy, 
                             family = "zero_inflated_poisson"),
               "real zero_inflated_poisson_lpmf(int y")
  expect_match2(make_stancode(count ~ Trt_c, data = epilepsy, 
                             family = "zero_inflated_negbinomial"),
               "real zero_inflated_neg_binomial_lpmf(int y")
  expect_match2(make_stancode(count ~ Trt_c, data = epilepsy, 
                             family = "zero_inflated_binomial"),
               "real zero_inflated_binomial_lpmf(int y")
  expect_match2(make_stancode(count ~ Trt_c, data = epilepsy, 
                             family = "zero_inflated_beta"),
               "real zero_inflated_beta_lpdf(real y")
  expect_match2(make_stancode(count ~ Trt_c, data = epilepsy, 
                             family = hurdle_poisson()),
               "real hurdle_poisson_lpmf(int y")
  expect_match2(make_stancode(count ~ Trt_c, data = epilepsy, 
                             family = hurdle_negbinomial),
               "real hurdle_neg_binomial_lpmf(int y")
  expect_match2(make_stancode(count ~ Trt_c, data = epilepsy, 
                             family = hurdle_gamma("log")),
               "real hurdle_gamma_lpdf(real y")
  expect_match2(make_stancode(count ~ Trt_c, data = epilepsy, 
                             family = hurdle_lognormal("identity")),
               "real hurdle_lognormal_lpdf(real y")
  
  # linear models with special covariance structures
  expect_match2(make_stancode(rating ~ treat, data = inhaler, 
                             autocor = cor_ma(cov = TRUE)),
               "real normal_cov_lpdf(vector y")
  expect_match2(make_stancode(time ~ age, data = kidney, family = "student", 
                             autocor = cor_ar(cov = TRUE)),
               "real student_t_cov_lpdf(vector y")
  
  # ARMA covariance matrices
  expect_match2(make_stancode(rating ~ treat, data = inhaler, 
                             autocor = cor_ar(cov = TRUE)),
               "matrix cov_matrix_ar1(real ar")
  expect_match2(make_stancode(time ~ age, data = kidney, family = "student", 
                             autocor = cor_ma(cov = TRUE)),
               "matrix cov_matrix_ma1(real ma")
  expect_match2(make_stancode(time ~ age, data = kidney, family = "student", 
                             autocor = cor_arma(p = 1, q = 1, cov = TRUE)),
               "matrix cov_matrix_arma1(real ar, real ma")
  
  # kronecker matrices
  expect_match(make_stancode(rating ~ treat + period + carry + (1+carry|subject), 
                             data = inhaler, cov_ranef = list(subject = 1)), 
              "matrix as_matrix.*matrix kronecker")
})

test_that("invalid combinations of modeling options are detected", {
  data <- data.frame(y1 = rnorm(10), y2 = rnorm(10), 
                     wi = 1:10, ci = sample(-1:1, 10, TRUE))
  expect_error(make_stancode(y1 | cens(ci) ~ y2, data = data,
                             autocor = cor_ar(cov = TRUE)),
               "Invalid addition arguments for this model")
  expect_error(make_stancode(cbind(y1, y2) ~ 1, data = data,
                             autocor = cor_ar(cov = TRUE)),
               "ARMA covariance matrices are not yet working")
  expect_error(make_stancode(y1 | resp_se(wi) ~ y2, data = data,
                             autocor = cor_ma()),
               "Please set cov = TRUE")
  expect_error(make_stancode(y1 | trunc(lb = -50) | weights(wi) ~ y2,
                             data = data),
               "Truncation is not yet possible")
})

test_that("Stan code for multivariate models is correct", {
  dat <- data.frame(y1 = rnorm(10), y2 = rnorm(10), x = 1:10)
  scode <- make_stancode(cbind(y1, y2) ~ x, dat)
  expect_match2(scode, "Y ~ multi_normal_cholesky(Eta, LSigma);")
  expect_match2(scode, "LSigma = diag_pre_multiply(sigma, Lrescor);")
  
  scode <- make_stancode(cbind(y1, y2) ~ x, dat, student())
  expect_match2(scode, "Y ~ multi_student_t(nu, Eta, Sigma);")
  expect_match2(scode, "cov_matrix[nresp] Sigma;")
  
  expect_match2(make_stancode(cbind(y1, y2) | weights(x) ~ 1, dat),
                "lp_pre[n] = multi_normal_cholesky_lpdf(Y[n] | Eta[n], LSigma);")
})

test_that("Stan code for categorical models is correct", {
  dat <- data.frame(y = rep(1:4, 2), x = 1:8, g = 1:8)
  scode <- make_stancode(y ~ x + (1|ID|g), dat, categorical())
  expect_match2(scode, "Y[n] ~ categorical_logit(append_row(zero, Eta[n]));")
  expect_match2(scode, "eta_2 = Xc_2 * b_2 + temp_2_Intercept;")
  expect_match2(scode, "eta_4[n] = eta_4[n] + (r_1_4_3[J_1[n]]) * Z_1_4_3[n];")
})

test_that("Stan code for autocorrelated models is correct", {
  dat <- data.frame(y = rep(1:4, 2), x = 1:8, time = 1:8)
  scode <- make_stancode(y ~ x, dat, student(log), 
                         autocor = cor_ar(~time))
  expect_match2(scode, "e[n] = log(Y[n]) - eta[n];")
  expect_match2(scode,
    "eta[n] = eta[n] + head(E[n], Kar) * ar; \n    eta[n] = exp(eta[n]);")
  
  scode <- make_stancode(y ~ x, dat, student(log), 
                         autocor = cor_ma(~time, q = 2))
  expect_match2(scode, "eta[n] = eta[n] + head(E[n], Kma) * ma;")
  
  scode <- make_stancode(y ~ x, dat, student(log), 
                         autocor = cor_arr(~time, r = 2))
  expect_match2(scode, "eta = Xc * b + temp_Intercept + Yarr * arr;")
  
  scode <- make_stancode(cbind(y, x) ~ 1, dat, gaussian(inverse),
                         autocor = cor_ar())
  expect_match2(scode, "e_y[n] = inv(Y[n, 1]) - eta_y[n];")
})

test_that("the Stan code for intercept only models is correct", {
  expect_match2(make_stancode(rating ~ 1, data = inhaler),
               "b_Intercept = temp_Intercept;") 
  expect_match2(make_stancode(rating ~ 1, data = inhaler, family = cratio()),
               "b_Intercept = temp_Intercept;") 
  expect_match2(make_stancode(rating ~ 1, data = inhaler, family = categorical()),
               "b_3_Intercept = temp_3_Intercept;") 
})

test_that("make_stancode returns correct code for spline only models", {
  expect_match2(make_stancode(count ~ s(log_Age_c), data = epilepsy),
               "matrix[N, K - 1] Xc;")
})

test_that("Stan code of ordinal models is correct", {
  dat <- data.frame(y = c(rep(1:4, 2), 1, 1), x1 = rnorm(10), 
                    x2 = rnorm(10), g = rep(1:2, 5))
  
  scode <- make_stancode(y ~ x1, dat, family = cumulative())
  expect_match2(scode, "Y[n] ~ ordered_logistic(eta[n], temp_Intercept);")
  
  scode <- make_stancode(y ~ x1, dat, family = cumulative("probit"),
                         threshold = "equidistant")
  expect_match2(scode, "real cumulative_lpmf(int y")
  expect_match2(scode, "p[1] = Phi(disc * (thres[1] - eta));")
  expect_match2(scode, "temp_Intercept[k] = temp_Intercept1 + (k - 1.0) * delta;")
  expect_match2(scode, "b_Intercept = temp_Intercept + dot_product(means_X, b);")
  
  scode <- make_stancode(y ~ x1, dat, family = cratio("probit_approx"))
  expect_match2(scode, "real cratio_lpmf(int y")
  expect_match2(scode, "q[k] = Phi_approx(disc * (eta - thres[k]));")

  scode <- make_stancode(y ~ x1 + cs(x2), dat, family = sratio())
  expect_match2(scode, "real sratio_lpmf(int y")
  expect_match2(scode, "matrix[N, Kcs] Xcs;")
  expect_match2(scode, "matrix[Kcs, ncat - 1] bcs;")
  expect_match2(scode, "etacs = Xcs * bcs;")
  expect_match2(scode, "sratio(eta[n], etacs[n], temp_Intercept, disc);")
  
  scode <- make_stancode(y ~ x1 + cse(x2) + (cse(1)|g), dat, family = acat())
  expect_match2(scode, "real acat_lpmf(int y")
  expect_match2(scode, "q[k] = disc * (eta + etacs[k] - thres[k]);")
  expect_match2(scode, "etacs[n, 1] = etacs[n, 1] + r_1_1[J_1[n]] * Z_1_1[n];")
  expect_match2(scode, "b_Intercept = temp_Intercept - dot_product(means_X, b);")
  
  scode <- make_stancode(y ~ x1 + (cse(x2)||g), dat, family = acat("probit"))
  expect_match2(scode, "q[k] = Phi(disc * (eta + etacs[k] - thres[k]));")
  expect_match2(scode, 
    paste("etacs[n, 3] = etacs[n, 3] + r_1_3[J_1[n]] * Z_1_3[n]", 
          "+ r_1_6[J_1[n]] * Z_1_6[n];"))
  expect_match2(scode, "Y[n] ~ acat(eta[n], etacs[n], temp_Intercept, disc);")
})

test_that("ordinal disc parameters appear in the Stan code", {
  scode <- make_stancode(bf(rating ~ period + carry + treat, disc ~ 1),
                         data = inhaler, family = cumulative(), 
                         prior = c(prior(normal(0,5), nlpar = disc)))
  expect_match2(scode, "Y[n] ~ cumulative(eta[n], temp_Intercept, disc[n])")
  expect_match2(scode, "b_disc ~ normal(0, 5)")
  expect_match2(scode, "disc[n] = exp(disc[n])")
})

test_that("monotonic effects appear in the Stan code", {
  prior <- c(prior(normal(0,1), class = b, coef = x1),
             prior(dirichlet(c(1,0.5,2)), simplex, coef = x1),
             prior(dirichlet(c(1,0.5,2)), simplex, coef = x2))
  dat <- data.frame(y = rpois(120, 10), x1 = rep(1:4, 30), 
                    x2 = factor(rep(c("a", "b", "c"), 40), ordered = TRUE))
  scode <- make_stancode(y ~ mo(x1 + x2), dat, prior = prior)
  expect_match2(scode, "int Xmo[N, Kmo];")
  expect_match2(scode, "simplex[Jmo[1]] simplex_1;")
  expect_match2(scode, "(bmo[2]) * monotonic(simplex_2, Xmo[n, 2]);")
  expect_match2(scode, "bmo[1] ~ normal(0, 1)")
  expect_match2(scode, "simplex_1 ~ dirichlet(con_simplex_1);")
  expect_match2(scode, "simplex_2 ~ dirichlet(con_simplex_2);")
  scode <- make_stancode(y ~ mono(x1) + (mono(x1)|x2) + (1|x2), dat)
  expect_match2(scode, "(bmo[1] + r_1_1[J_1[n]]) * monotonic(simplex_1, Xmo[n, 1]);")
  # test that Z_1_1 is (correctly) undefined
  expect_match2(scode, paste0("  int<lower=1> M_1; \n",
    "  // data for group-level effects of ID 2"))
  expect_error(make_stancode(y ~ mono(x1) + (mono(x1+x2)|x2), dat),
               "Monotonic group-level terms require")
  expect_error(make_stancode(y ~ mo(x1), dat, 
                             prior = prior(beta(1, 1), simplex, coef = x1)),
               "'dirichlet' is the only valid prior for simplex parameters")
})

test_that("Stan code for non-linear models is correct", {
  flist <- list(a ~ x, b ~ z + (1|g))
  data <- data.frame(y = rgamma(9, 1, 1), x = rnorm(9), z = rnorm(9), 
                     g = rep(1:3, 3))
  prior <- c(set_prior("normal(0,5)", nlpar = "a"),
             set_prior("normal(0,1)", nlpar = "b"))
  # syntactic validity is already checked within make_stancode
  scode <- make_stancode(bf(y ~ a - exp(b^z), flist = flist, nl = TRUE), 
                            data = data, prior = prior)
  expect_match2(scode, "eta[n] = eta_a[n] - exp(eta_b[n] ^ C[n, 1]);")
  
  flist <- list(a1 ~ 1, a2 ~ z + (x|g))
  prior <- c(set_prior("beta(1,1)", nlpar = "a1", lb = 0, ub = 1),
             set_prior("normal(0,1)", nlpar = "a2"))
  scode <- make_stancode(bf(y ~ a1 * exp(-x/(a2 + z)), 
                               flist = flist, nl = TRUE),
                            data = data, family = Gamma("log"), prior = prior)
  expect_match2(scode,
    paste("eta[n] = shape * exp(-(eta_a1[n] *", 
          "exp( - C[n, 1] / (eta_a2[n] + C[n, 2]))));"))
})

test_that("make_stancode accepts very long non-linear formulas", {
  data <- data.frame(y = rnorm(10), this_is_a_very_long_predictor = rnorm(10))
  expect_silent(make_stancode(bf(y ~ b0 + this_is_a_very_long_predictor + 
                                 this_is_a_very_long_predictor +
                                 this_is_a_very_long_predictor,
                                 b0 ~ 1, nl = TRUE),
                data = data, prior = prior(normal(0,1), nlpar = "b0")))
})

test_that("no loop in trans-par is defined for simple 'identity' models", {
  expect_true(!grepl(make_stancode(time ~ age, data = kidney),
                     "eta[n] <- (eta[n]);", fixed = TRUE))
  expect_true(!grepl(make_stancode(time ~ age, data = kidney, 
                                   family = poisson("identity")), 
                     "eta[n] <- (eta[n]);", fixed = TRUE))
})

test_that("known standard errors appear in the Stan code", {
  scode <- make_stancode(time | se(age) ~ sex, data = kidney)
  expect_match2(scode, "Y ~ normal(eta, se)")
  scode <- make_stancode(time | se(age) | weights(age) ~ sex, data = kidney)
  expect_match2(scode, "lp_pre[n] = normal_lpdf(Y[n] | eta[n], se[n])")
  scode <- make_stancode(time | se(age, sigma = TRUE) ~ sex, data = kidney)
  expect_match2(scode, "Y[n] ~ normal(eta[n], sqrt(sigma^2 + se2[n]))")
  scode <- make_stancode(bf(time | se(age, sigma = TRUE) ~ sex, sigma ~ sex),
                         data = kidney)
  expect_match2(scode, "Y[n] ~ normal(eta[n], sqrt(sigma[n]^2 + se2[n]))")
})

test_that("Addition term 'disp' appears in the Stan code", {
  scode <- make_stancode(time | disp(sqrt(age)) ~ sex + age, data = kidney)
  expect_match2(scode, "disp_sigma = sigma * disp;")
  expect_match2(scode, "Y ~ normal(eta, disp_sigma)")
  
  scode <- make_stancode(time | disp(1/age) ~ sex + age, 
                         data = kidney, family = lognormal())
  expect_match2(scode, "Y ~ lognormal(eta, disp_sigma);")
  
  scode <- make_stancode(time | disp(1/age) ~ sex + age + (1|patient), 
                         data = kidney, family = Gamma())
  expect_match2(scode, "eta[n] = disp_shape[n] * (")
  expect_match2(scode, "Y ~ gamma(disp_shape, eta)")
  
  scode <- make_stancode(time | disp(1/age) ~ sex + age, 
                         data = kidney, family = weibull())
  expect_match2(scode, "eta[n] = exp((eta[n]) / disp_shape[n]);")
  
  scode <- make_stancode(time | disp(1/age) ~ sex + age, 
                         data = kidney, family = negbinomial())
  expect_match2(scode, "Y ~ neg_binomial_2_log(eta, disp_shape);")
  
  scode <- make_stancode(bf(y | disp(y) ~ a - b^x, a + b ~ 1, nl = TRUE),
                            family = weibull(),
                            data = data.frame(y = rpois(10, 10), x = rnorm(10)),
                            prior = c(set_prior("normal(0,1)", nlpar = "a"),
                                      set_prior("normal(0,1)", nlpar = "b")))
  expect_match2(scode,
    "eta[n] = exp((eta_a[n] - eta_b[n] ^ C[n, 1]) / disp_shape[n]);")
})

test_that("functions defined in 'stan_funs' appear in the functions block", {
  test_fun <- paste0("  real test_fun(real a, real b) { \n",
                     "    return a + b; \n",
                     "  } \n")
  expect_match2(make_stancode(time ~ age, data = kidney, stan_funs = test_fun),
               test_fun)
})

test_that("fixed residual covariance matrices appear in the Stan code", {
  data <- data.frame(y = 1:5)
  V <- diag(5)
  expect_match2(make_stancode(y~1, data = data, family = gaussian(), 
                             autocor = cor_fixed(V)),
               "Y ~ multi_normal_cholesky(eta, LV)")
  expect_match2(make_stancode(y~1, data = data, family = student(),
                             autocor = cor_fixed(V)),
               "Y ~ multi_student_t(nu, eta, V)")
  expect_match2(make_stancode(y~1, data = data, family = student(),
                             autocor = cor_fixed(V)),
               "Y ~ multi_student_t(nu, eta, V)")
})

test_that("Stan code for bsts models is correct", {
  dat <- data.frame(y = rnorm(10), x = rnorm(10))
  scode <- make_stancode(y ~ x, data = dat, autocor = cor_bsts(),
                         prior = prior(normal(0, 10), sigmaLL))
  expect_match2(scode, "+ loclev[n]")
  expect_match2(scode, "loclev[n] ~ normal(loclev[n - 1], sigmaLL)")
  expect_match2(scode, "sigmaLL ~ normal(0, 10)")
  
  dat <- data.frame(y = rexp(1), x = rnorm(10))
  scode <- make_stancode(y~x, data = dat, family = student("log"),
                         autocor = cor_bsts())
  expect_match2(scode, "loclev[n] ~ normal(log(Y[n]), sigmaLL)")
})

test_that("Stan code for GAMMs is correct", {
  dat <- data.frame(y = rnorm(10), x = rnorm(10), g = rep(1:2, 5))
  scode <- make_stancode(y ~ s(x) + (1|g), data = dat,
                            prior = set_prior("normal(0,2)", "sds"))
  expect_match2(scode, "Zs_1_1 * s_1_1")
  expect_match2(scode, "matrix[N, knots_1[1]] Zs_1_1")
  expect_match2(scode, "zs_1_1 ~ normal(0, 1)")
  expect_match2(scode, "sds_1_1 ~ normal(0,2)")
  
  prior <- c(set_prior("normal(0,5)", nlpar = "lp"),
             set_prior("normal(0,2)", "sds", nlpar = "lp"))
  scode <- make_stancode(bf(y ~ lp, lp ~ s(x) + (1|g), nl = TRUE), 
                            data = dat, prior = prior)
  expect_match2(scode, "Zs_lp_1_1 * s_lp_1_1")
  expect_match2(scode, "matrix[N, knots_lp_1[1]] Zs_lp_1_1")
  expect_match2(scode, "zs_lp_1_1 ~ normal(0, 1)")
  expect_match2(scode, "sds_lp_1_1 ~ normal(0,2)")
  
  scode <- make_stancode(y ~ s(x) + t2(x,y), data = dat,
                            prior = set_prior("normal(0,2)", "sds"))
  expect_match2(scode, "Zs_2_2 * s_2_2")
  expect_match2(scode, "matrix[N, knots_2[2]] Zs_2_2")
  expect_match2(scode, "zs_2_2 ~ normal(0, 1)")
  expect_match2(scode, "sds_2_2 ~ normal(0,2)")
})

test_that("Stan code of exgaussian models is correct", {
  dat <- epilepsy
  dat$cens <- sample(-1:1, nrow(dat), TRUE)
  scode <- make_stancode(count ~ Trt_c + (1|patient),
                      data = dat, family = exgaussian("log"),
                      prior = prior(gamma(1,1), class = beta))
  expect_match2(scode, "Y[n] ~ exgaussian(eta[n], sigma, beta)")
  expect_match2(scode, "eta[n] = exp(eta[n])")
  expect_match2(scode, "beta ~ gamma(1, 1)")
  
  scode <- make_stancode(bf(count ~ Trt_c + (1|patient),
                         sigma ~ Trt_c, beta ~ Trt_c),
                      data = dat, family = exgaussian())
  expect_match2(scode, "Y[n] ~ exgaussian(eta[n], sigma[n], beta[n])")
  expect_match2(scode, "beta[n] = exp(beta[n])")
  
  scode <- make_stancode(count | cens(cens) ~ Trt_c + (1|patient),
                      data = dat, family = exgaussian("inverse"))
  expect_match2(scode, "exgaussian_lccdf(Y[n] | eta[n], sigma, beta)")
})

test_that("Stan code of wiener diffusion models is correct", {
  dat <- RWiener::rwiener(n=100, alpha=2, tau=.3, beta=.5, delta=.5)
  dat$x <- rnorm(100)
  scode <- make_stancode(q | dec(resp) ~ x, data = dat, family = wiener())
  expect_match2(scode, "Y[n] ~ wiener_diffusion(dec[n], bs, ndt, bias, eta[n])")
  
  scode <- make_stancode(bf(q | dec(resp) ~ x, bs ~ x, ndt ~ x, bias ~ x), 
                         data = dat, family = wiener())
  expect_match2(scode, "Y[n] ~ wiener_diffusion(dec[n], bs[n], ndt[n], bias[n], eta[n])")
  expect_match2(scode, "bias[n] = inv_logit(bias[n]);")
  
  expect_error(make_stancode(q ~ x, data = dat, family = wiener()),
               "Addition argument 'dec' is required for family 'wiener'")
})

test_that("Group IDs appear in the Stan code", {
  form <- bf(count ~ Trt_c + (1+Trt_c|3|visit) + (1|patient), 
             shape ~ (1|3|visit) + (Trt_c||patient))
  scode <- make_stancode(form, data = epilepsy, family = negbinomial())
  expect_match2(scode, "r_2_1 = r_2[, 1]")
  expect_match2(scode, "r_2_shape_3 = r_2[, 3]")
  
  form <- bf(count ~ a, sigma ~ (1|3|visit) + (Trt_c||patient),
             a ~ Trt_c + (1+Trt_c|3|visit) + (1|patient), nl = TRUE)
  scode <- make_stancode(form, data = epilepsy, family = student(),
                      prior = set_prior("normal(0,5)", nlpar = "a"))
  expect_match2(scode, "r_2_a_3 = r_2[, 3];")
  expect_match2(scode, "r_1_sigma_2 = sd_1[2] * (z_1[2]);")
})

test_that("distributional gamma models are handled correctly", {
  # test fix of issue #124
  scode <- make_stancode(bf(time ~ age * sex + disease + (1|patient), 
                            shape ~ age + (1|patient)), 
                         data = kidney, family = Gamma("log"))
  expect_match2(scode, paste0(
    "    shape[n] = exp(shape[n]); \n", 
    "    eta[n] = shape[n] * exp(-(eta[n]));"))
  
  scode <- make_stancode(bf(time ~ inv_logit(a) * exp(b * age),
                            a + b ~ sex + (1|patient), nl = TRUE, 
                            shape ~ age + (1|patient)), 
                         data = kidney, family = Gamma("identity"),
                         prior = c(set_prior("normal(2,2)", nlpar = "a"),
                                   set_prior("normal(0,3)", nlpar = "b")))
  expect_match2(scode, paste0(
    "    shape[n] = exp(shape[n]); \n", 
    "    // compute non-linear predictor \n",
    "    eta[n] = shape[n] / (inv_logit(eta_a[n]) * exp(eta_b[n] * C[n, 1]));"))
})

test_that("weighted and censored responses appear in the Stan code", {
  dat <- data.frame(y = 1:9, x = rep(-1:1, 3), y2 = 10:18)
  
  scode <- make_stancode(y | weights(y2) ~ 1, dat, poisson())
  expect_match2(scode, "lp_pre[n] = poisson_log_lpmf(Y[n] | eta[n]);")
  expect_match2(scode, "target += dot_product(weights, lp_pre);")
  
  scode <- make_stancode(y | trials(y2) + weights(y2) ~ 1, dat, binomial())
  expect_match2(scode, "lp_pre[n] = binomial_logit_lpmf(Y[n] | trials[n], eta[n]);")
  
  expect_match2(make_stancode(y | cens(x, y2) ~ 1, dat, poisson()),
               "target += poisson_lpmf(Y[n] | eta[n]); \n")
  expect_match2(make_stancode(y | cens(x) ~ 1, dat, weibull()), 
               "target += weibull_lccdf(Y[n] | shape, eta[n]); \n")
  dat$x[1] <- 2
  expect_match2(make_stancode(y | cens(x, y2) ~ 1, dat, gaussian()), 
               "target += log_diff_exp(normal_lcdf(rcens[n] | eta[n], sigma),")
  dat$x <- 1
  make_stancode(y | cens(x) + weights(x) ~ 1, dat, weibull())
  expect_match2(make_stancode(y | cens(x) + weights(x) ~ 1, dat, weibull()),
   "target += weights[n] * weibull_lccdf(Y[n] | shape, eta[n]); \n")
})

test_that("priors on intercepts appear in the Stan code", {
  scode <- make_stancode(count ~ log_Age_c + log_Base4_c * Trt_c,
                      data = epilepsy, family = gaussian(), 
                      prior = c(prior(student_t(5,0,10), class = b),
                                prior(normal(0,5), class = Intercept)),
                      sample_prior = TRUE)
  expect_match2(scode, paste0("prior_b_Intercept = prior_temp_Intercept ", 
                          "- dot_product(means_X, b);"))
  scode <- make_stancode(cbind(count, count) ~ log_Age_c + log_Base4_c * Trt_c,
                      data = epilepsy, family = gaussian(), 
                      prior = c(prior(student_t(5,0,10), class = b),
                                prior(normal(0,5), class = Intercept)),
                      sample_prior = TRUE)
  expect_match2(scode, paste0("prior_b_count_Intercept = prior_temp_count_Intercept ", 
                          "- dot_product(means_X_count, b_count);"))
})

test_that("noise-free terms appear in the Stan code", {
  N <- 30
  dat <- data.frame(y = rnorm(N), x = rnorm(N), z = rnorm(N),
                    xsd = abs(rnorm(N, 1)), zsd = abs(rnorm(N, 1)),
                    ID = rep(1:5, each = N / 5))
  scode <- make_stancode(y ~ me(x, xsd)*me(z, zsd)*x, data = dat,
                      prior = prior(normal(0,5)))
  expect_match2(scode,
    "(bme[1]) * Xme_1[n] + (bme[2]) * Xme_2[n] + (bme[3]) * Xme_1[n] .* Xme_2[n]")
  expect_match2(scode, 
    "(bme[6]) * Xme_1[n] .* Xme_2[n] .* Cme_3[n]")
  expect_match2(scode, "Xme_2 ~ normal(Xn_2, noise_2)")
  expect_match2(scode, "bme ~ normal(0, 5)")
  scode <- make_stancode(y ~ me(x, xsd)*me(z, zsd) + (me(x, xsd)|ID), data = dat)
  expect_match2(scode, "(bme[1] + r_1_1[J_1[n]]) * Xme_1[n]")
  expect_match2(make_stancode(y ~ I(me(x, xsd)^2), data = dat),
               "(bme[1]) * Xme_1[n]^2")
})

test_that("Stan code of multi-membership models is correct", {
  dat <- data.frame(y = rnorm(10), g1 = sample(1:10, 10, TRUE),
                    g2 = sample(1:10, 10, TRUE), w1 = rep(1, 10),
                    w2 = rep(abs(rnorm(10))))
  expect_match2(make_stancode(y ~ (1|mm(g1, g2)), data = dat), 
               "(W_1_1[n] * r_1_1[J_1_1[n]] + W_1_2[n] * r_1_1[J_1_2[n]]) * Z_1_1[n]")
  expect_match2(make_stancode(y ~ (1+w1|mm(g1,g2)), data = dat), 
               "(W_1_1[n] * r_1_2[J_1_1[n]] + W_1_2[n] * r_1_2[J_1_2[n]]) * Z_1_2[n]")
})

test_that("Group syntax | and || is handled correctly,", {
  data <- data.frame(y = rnorm(10), x = rnorm(10),
                     g1 = rep(1:5, each = 2), g2 = rep(1:2, 5))
  scode <- make_stancode(y ~ x + (1+x||g1) + (I(x/4)|g2), data)
  expect_match2(scode, "r_1_2 = sd_1[2] * (z_1[2]);")
  expect_match2(scode, "r_2_1 = r_2[, 1];")
  expect_match2(scode, "r_2 = (diag_pre_multiply(sd_2, L_2) * z_2)';")
})

test_that("fixing auxiliary parameters is possible", {
  scode <- make_stancode(bf(y ~ 1, sigma = 0.5), data = list(y = rnorm(10)))
  expect_match(scode, "data \\{[^\\}]*real<lower=0> sigma;")
})
