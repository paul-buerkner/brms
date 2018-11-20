context("Tests for make_stancode")

# simplifies manual calling of tests
expect_match2 <- brms:::expect_match2
SW <- brms:::SW

test_that("specified priors appear in the Stan code", {
  dat <- data.frame(y = 1:10, x1 = rnorm(10), x2 = rnorm(10), 
                    g = rep(1:5, 2), h = factor(rep(1:5, each = 2)))
  prior <- c(prior(normal(0,1), coef = x1),
             prior(normal(0,2), coef = x2),
             prior(normal(0,5), Intercept),
             prior(cauchy(0,1), sd, group = g),
             prior(cauchy(0,2), sd, group = g, coef = x1),
             prior(gamma(1, 1), class = sd, group = h))
  scode <- make_stancode(y ~ x1*x2 + (x1*x2|g) + (1 | h), dat,
                         prior = prior, sample_prior = "yes")
  expect_match2(scode, "target += normal_lpdf(b[1] | 0, 1)")
  expect_match2(scode, "target += normal_lpdf(b[2] | 0, 2)")
  expect_match2(scode, "target += normal_lpdf(temp_Intercept | 0, 5)")
  expect_match2(scode, "target += cauchy_lpdf(sd_1[1] | 0, 1)")
  expect_match2(scode, "- 1 * cauchy_lccdf(0 | 0, 1)")
  expect_match2(scode, "target += cauchy_lpdf(sd_1[2] | 0, 2)")
  expect_match2(scode, "target += student_t_lpdf(sigma | 3, 0, 10)")
  expect_match2(scode, "- 1 * student_t_lccdf(0 | 3, 0, 10)")
  expect_match2(scode, "target += gamma_lpdf(sd_2 | 1, 1)")
  expect_match2(scode, "prior_b_1 = normal_rng(0,1);")
  expect_match2(scode, "prior_sd_1_1 = cauchy_rng(0,1)")
  expect_match2(scode, "while (prior_sd_1_1 < 0)")
  expect_match2(scode, "prior_sd_2 = gamma_rng(1,1)")
  expect_match2(scode, "while (prior_sd_2 < 0)")
  
  prior <- c(prior(lkj(0.5), class = cor, group = g),
             prior(normal(0, 1), class = b),
             prior(normal(0, 5), class = Intercept),
             prior(cauchy(0, 5), class = sd))
  scode <- make_stancode(y ~ x1 + cs(x2) + (0 + x1 + x2 | g), 
                         dat, family = acat(), 
                         prior = prior, sample_prior = TRUE)
  expect_match2(scode, "target += normal_lpdf(b | 0, 1)")
  expect_match2(scode, "target += normal_lpdf(temp_Intercept | 0, 5)")
  expect_match2(scode, "target += cauchy_lpdf(sd_1 | 0, 5)")
  expect_match2(scode, "target += lkj_corr_cholesky_lpdf(L_1 | 0.5)")
  expect_match2(scode, "target += normal_lpdf(to_vector(bcs) | 0, 1)")
  expect_match2(scode, "prior_bcs = normal_rng(0,1)")
  
  prior <- c(prior(normal(0,5), nlpar = a),
             prior(normal(0,10), nlpar = b),
             prior(cauchy(0,1), class = sd, nlpar = a),
             prior(lkj(2), class = cor, group = g))
  scode <- make_stancode(
    bf(y ~ a * exp(-b * x1), a + b ~ (1|ID|g), nl = TRUE),
    data = dat, prior = prior, sample_prior = TRUE
  )
  expect_match2(scode, "target += normal_lpdf(b_a | 0, 5)")
  expect_match2(scode, "target += normal_lpdf(b_b | 0, 10)")
  expect_match2(scode, "target += cauchy_lpdf(sd_1[1] | 0, 1)")
  expect_match2(scode, "target += lkj_corr_cholesky_lpdf(L_1 | 2)")
  expect_match2(scode, "prior_b_a = normal_rng(0,5)")
  expect_match2(scode, "prior_sd_1_2 = student_t_rng(3,0,10)")
  expect_match2(scode, "prior_cor_1 = lkj_corr_rng(M_1,2)[1, 2]")
  
  prior <- c(prior(lkj(2), rescor),
             prior(cauchy(0, 5), sigma, resp = y),
             prior(cauchy(0, 1), sigma, resp = x1))
  scode <- make_stancode(cbind(y, x1) ~ x2, dat, prior = prior, 
                         sample_prior = TRUE)
  expect_match2(scode, "target += lkj_corr_cholesky_lpdf(Lrescor | 2)")
  expect_match2(scode, "prior_sigma_y = cauchy_rng(0,5)")
  expect_match2(scode, "prior_rescor = lkj_corr_rng(nresp,2)[1, 2]")
  
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
  expect_match2(scode, "target += uniform_lpdf(ar | -1, 1)")
  expect_match2(scode, "target += normal_lpdf(ma | 0, 0.5)")
  expect_match2(scode, 
    "- 1 * log_diff_exp(normal_lcdf(1 | 0, 0.5), normal_lcdf(-1 | 0, 0.5))"
  )
  expect_match2(scode, "target += double_exponential_lpdf(arr | 0, 1)")
  expect_match2(scode, "target += normal_lpdf(bsp | 0, 5)")
  expect_match2(scode, "target += dirichlet_lpdf(simo_1 | con_simo_1)")
  expect_match2(scode, "prior_simo_1 = dirichlet_rng(con_simo_1)")
  expect_match2(scode, "prior_ar = uniform_rng(-1,1)")
  expect_match2(scode, "while (prior_ar < -0.7 || prior_ar > 0.5)")
  
  # test for problem described in #213
  prior <- c(prior(normal(0, 1), coef = x1),
             prior(normal(0, 2), coef = x1, dpar = sigma))
  scode <- make_stancode(bf(y ~ x1, sigma ~ x1), dat, prior = prior)
  expect_match2(scode, "target += normal_lpdf(b | 0, 1);")
  expect_match2(scode, "target += normal_lpdf(b_sigma | 0, 2);")
  
  prior <- c(set_prior("target += normal_lpdf(b[1] | 0, 1)", check = FALSE),
             set_prior("", class = "sigma"))
  scode <- make_stancode(y ~ x1, dat, prior = prior, sample_prior = TRUE)
  expect_match2(scode, "target += normal_lpdf(b[1] | 0, 1)")
  expect_true(!grepl("sigma \\|", scode))
  
  prior <- prior(gamma(0, 1), coef = x1)
  expect_warning(make_stancode(y ~ x1, dat, prior = prior),
                 "no natural lower bound")
  prior <- prior(uniform(0,5), class = sd)
  expect_warning(make_stancode(y ~ x1 + (1|g), dat, prior = prior),
                  "no natural upper bound")
  prior <- prior(uniform(-1, 1), class = cor)
  expect_error(
    make_stancode(y ~ x1 + (x1|g), dat, prior = prior),
    "prior for correlation matrices is the 'lkj' prior"
  )
})

test_that("special shrinkage priors appear in the Stan code", {
  dat <- data.frame(y = 1:10, x1 = rnorm(10), x2 = rnorm(10))
  
  # horseshoe prior
  hs <- horseshoe(7, scale_global = 2, df_global = 3,
                  df_slab = 6, scale_slab = 3)
  scode <- make_stancode(y ~ x1*x2, data = dat, 
                         prior = set_prior(hs),
                         sample_prior = TRUE)
  expect_match2(scode, "vector<lower=0>[Kc] hs_local[2];") 
  expect_match2(scode, "real<lower=0> hs_global[2];") 
  expect_match2(scode, 
    "target += inv_gamma_lpdf(hs_local[2] | 0.5 * hs_df, 0.5 * hs_df);"
  )
  expect_match2(scode, 
    "target += inv_gamma_lpdf(hs_global[2] | 0.5 * hs_df_global, 0.5 * hs_df_global);"
  )
  expect_match2(scode, 
    "target += inv_gamma_lpdf(hs_c2 | 0.5 * hs_df_slab, 0.5 * hs_df_slab);"
  )
  expect_match2(scode, 
    paste0(
      "b = horseshoe(zb, hs_local, hs_global, ", 
      "hs_scale_global * sigma, hs_scale_slab^2 * hs_c2);"
    )
  )
  
  scode <- make_stancode(y ~ x1*x2, data = dat, poisson(),
                         prior = prior(horseshoe(scale_global = 3)))
  expect_match2(scode, 
    paste0(
      "b = horseshoe(zb, hs_local, hs_global, ", 
      "hs_scale_global, hs_scale_slab^2 * hs_c2);"
    )
  )
  
  scode <- make_stancode(x1 ~ mo(y), dat, prior = prior(horseshoe()))
  expect_match2(scode, "target += normal_lpdf(zbsp | 0, 1);")
  expect_match2(scode,
    "target += normal_lpdf(hs_localsp[1] | 0, 1)\n    - 1 * log(0.5);"          
  )
  expect_match2(scode,
    "target += inv_gamma_lpdf(hs_localsp[2] | 0.5 * hs_df, 0.5 * hs_df);"          
  )
  expect_match2(scode,
    paste0(
      "bsp = horseshoe(zbsp, hs_localsp, hs_global, ", 
      "hs_scale_global * sigma, hs_scale_slab^2 * hs_c2);"
    )
  )
  
  # lasso prior
  scode <- make_stancode(y ~ x1*x2, data = dat,
                         prior = prior(lasso(2, scale = 10)),
                         sample_prior = TRUE)
  expect_match2(scode, "target += chi_square_lpdf(lasso_inv_lambda | lasso_df);")
  expect_match2(scode, 
    "target += double_exponential_lpdf(b | 0, lasso_scale * lasso_inv_lambda);"
  )
  
  scode <- make_stancode(x1 ~ mo(y), dat, prior = prior(lasso()))
  expect_match2(scode, 
    "double_exponential_lpdf(bsp | 0, lasso_scale * lasso_inv_lambda)"
  )
  
  # horseshoe and lasso prior applied in a non-linear model
  hs_a1 <- horseshoe(7, scale_global = 2, df_global = 3)
  lasso_a2 <- lasso(2, scale = 10)
  scode <- make_stancode(
    bf(y ~ a1 + a2, a1 ~ x1, a2 ~ 0 + x2, nl = TRUE),
    data = dat, sample_prior = TRUE,
    prior = c(set_prior(hs_a1, nlpar = "a1"),
              set_prior(lasso_a2, nlpar = "a2"))
  )
  expect_match2(scode, "vector<lower=0>[K_a1] hs_local_a1[2];")
  expect_match2(scode, "real<lower=0> hs_global_a1[2];")
  expect_match2(scode, 
    "target += inv_gamma_lpdf(hs_local_a1[2] | 0.5 * hs_df_a1, 0.5 * hs_df_a1);"
  )
  expect_match2(scode, 
    "target += inv_gamma_lpdf(hs_global_a1[2] | 0.5 * hs_df_global_a1, 0.5 * hs_df_global_a1);"
  )
  expect_match2(scode, 
    "target += inv_gamma_lpdf(hs_c2_a1 | 0.5 * hs_df_slab_a1, 0.5 * hs_df_slab_a1);"
  )
  expect_match2(scode, 
    paste0(
      "b_a1 = horseshoe(zb_a1, hs_local_a1, hs_global_a1, ", 
      "hs_scale_global_a1 * sigma, hs_scale_slab_a1^2 * hs_c2_a1);"
    )
  )
  expect_match2(scode,
    "target += chi_square_lpdf(lasso_inv_lambda_a2 | lasso_df_a2);"
  )
  expect_match2(scode, 
    "target += double_exponential_lpdf(b_a2 | 0, lasso_scale_a2 * lasso_inv_lambda_a2);"
  )
  
  # check error messages
  expect_error(make_stancode(y ~ x1*x2, data = dat, 
                             prior = prior(horseshoe(-1))),
               "Degrees of freedom of the local priors")
  expect_error(make_stancode(y ~ x1*x2, data = dat, 
                             prior = prior(horseshoe(1, -1))),
               "Scale of the global prior")
  expect_error(make_stancode(y ~ x1*x2, data = dat, prior = prior(lasso(-1))),
               "Degrees of freedom of the shrinkage parameter prior")
  expect_error(make_stancode(x1 ~ cs(y), dat, acat(), prior = prior(lasso())),
               "Horseshoe or lasso priors are not yet allowed")
  bprior <- prior(horseshoe()) + prior(normal(0, 1), coef = "y")
  expect_error(make_stancode(x1 ~ y, dat, prior = bprior),
               "Defining priors for single population-level parameters is not")
  expect_error(make_stancode(x1 ~ y, dat, prior = prior(lasso(), lb = 0)),
               "Boundaries for population-level effects are not allowed")
})

test_that("link functions appear in the Stan code", {
  dat <- data.frame(y = 1:10, x = rnorm(10))
  expect_match2(make_stancode(y ~ x, dat, family = poisson()), 
               "target += poisson_log_lpmf(Y | mu);")
  expect_match2(make_stancode(cbind(y, y + 1) ~ x, dat, family = gaussian("log")), 
               "mu_y[n] = exp(mu_y[n]);")
  expect_match2(make_stancode(y ~ x, dat, family = von_mises(tan_half)), 
               "mu[n] = inv_tan_half(mu[n]);")
  expect_match2(make_stancode(y ~ x, dat, family = weibull()),
                "mu[n] = exp(mu[n]) / tgamma(1 + 1 / shape);")
  expect_match2(make_stancode(y ~ x, dat, family = exponential("identity")),
               "mu[n] = inv(mu[n]);")
  expect_match2(make_stancode(y ~ x, dat, family = poisson("sqrt")),
               "mu[n] = square(mu[n]);")
  expect_match2(make_stancode(y ~ x, dat, family = bernoulli()),
                "target += bernoulli_logit_lpmf(Y | mu);")
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
  expect_match2(scode, " r_1_1 = sd_1[1] * (Lcov_1 * z_1[1]);")
  expect_match2(scode, " r_1_2 = sd_1[2] * (Lcov_1 * z_1[2]);")
})

test_that("truncation appears in the Stan code", {
  scode <- make_stancode(time | trunc(0) ~ age + sex + disease,
                         data = kidney, family = "gamma")
  expect_match2(scode, "target += gamma_lpdf(Y[n] | shape, mu[n]) -")
  expect_match2(scode, "gamma_lccdf(lb[n] | shape, mu[n]);")
  
  scode <- make_stancode(time | trunc(ub = 100) ~ age + sex + disease, 
                         data = kidney, family = student("log"))
  
  expect_match2(scode, "target += student_t_lpdf(Y[n] | nu, mu[n], sigma) -")
  expect_match2(scode, "student_t_lcdf(ub[n] | nu, mu[n], sigma);")
  
  scode <- make_stancode(count | trunc(0, 150) ~ Trt_c, 
                         data = epilepsy, family = "poisson")
  expect_match2(scode, "target += poisson_lpmf(Y[n] | mu[n]) -")
  expect_match2(scode, 
    "log_diff_exp(poisson_lcdf(ub[n] | mu[n]), poisson_lcdf(lb[n] | mu[n]));"
  )
})

test_that("make_stancode combines strings of multiple grouping factors", {
  expect_match2(make_stancode(count ~ (1|patient) + (1 + Trt_c | visit), 
                             data = epilepsy, family = "poisson"), 
               paste0("  vector[N] Z_1_1;\n",
                      "  // data for group-level effects of ID 2"))
  expect_match2(make_stancode(count ~ (1|visit) + (1+Trt_c|patient), 
                             data = epilepsy, family = skew_normal()), 
               paste0("  int<lower=1> NC_1;\n",
                      "  // data for group-level effects of ID 2"))
})

test_that("make_stancode handles models without fixed effects", {
  expect_match2(make_stancode(count ~ 0 + (1|patient) + (1+Trt_c|visit), 
                             data = epilepsy, family = "poisson"), 
               "mu = rep_vector(0, N);")
})

test_that("make_stancode correctly restricts FE parameters", {
  data <- data.frame(y = rep(0:1, each = 5), x = rnorm(10))
  
  scode <- make_stancode(y ~ x, data, prior = set_prior("", lb = 2))
  expect_match2(scode, "vector<lower=2>[Kc] b")
  
  scode <- make_stancode(
    y ~ x, data, prior = set_prior("normal (0, 2)", ub = "4")
  )
  expect_match2(scode, "vector<upper=4>[Kc] b")
  expect_match2(scode, "- 1 * normal_lcdf(4 | 0, 2)")
  
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
                              family = "zero_one_inflated_beta"),
                "real zero_one_inflated_beta_lpdf(real y")
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
  
  # kronecker products
  expect_match(make_stancode(rating ~ treat + period + carry + (1+carry|subject), 
                             data = inhaler, cov_ranef = list(subject = 1)), 
              "matrix as_matrix.*matrix kronecker")
})

test_that("invalid combinations of modeling options are detected", {
  data <- data.frame(y1 = rnorm(10), y2 = rnorm(10), 
                     wi = 1:10, ci = sample(-1:1, 10, TRUE))
  expect_error(
    make_stancode(y1 | cens(ci) ~ y2, data = data, autocor = cor_ar(cov = TRUE)),
    "Invalid addition arguments for this model"
  )
  expect_error(
    make_stancode(cbind(y1, y2) ~ 1, data = data, autocor = cor_ar(cov = TRUE)),
    "ARMA covariance matrices are not implemented when 'rescor' is estimated."
  )
  expect_error(
    make_stancode(y1 | resp_se(wi) ~ y2, data = data, autocor = cor_ma()),
    "Please set cov = TRUE"
  )
})

test_that("Stan code for multivariate models is correct", {
  dat <- data.frame(
    y1 = rnorm(10), y2 = rnorm(10), 
    x = 1:10, g = rep(1:2, each = 5),
    censi = sample(0:1, 10, TRUE)
  )
  # models with residual correlations
  scode <- make_stancode(cbind(y1, y2) ~ x, dat, prior = prior(horseshoe(2)))
  expect_match2(scode, "target += multi_normal_cholesky_lpdf(Y | Mu, LSigma);")
  expect_match2(scode, "LSigma = diag_pre_multiply(sigma, Lrescor);")
  expect_match2(scode, "target += normal_lpdf(hs_local_y1[1] | 0, 1)")
  expect_match2(scode, 
    "target += inv_gamma_lpdf(hs_local_y2[2] | 0.5 * hs_df_y2, 0.5 * hs_df_y2)"
  )
  
  scode <- make_stancode(cbind(y1, y2) ~ x, dat, student(),
                         prior = prior(lasso(2, 10)))
  expect_match2(scode, "target += multi_student_t_lpdf(Y | nu, Mu, Sigma);")
  expect_match2(scode, "matrix[nresp, nresp] Sigma = multiply_lower")
  expect_match2(scode, "target += gamma_lpdf(nu | 2, 0.1)")
  expect_match2(scode, 
    "target += chi_square_lpdf(lasso_inv_lambda_y1 | lasso_df_y1)"
  )
  expect_match2(scode, 
    "target += chi_square_lpdf(lasso_inv_lambda_y2 | lasso_df_y2)"
  )
  expect_match2(make_stancode(cbind(y1, y2) | weights(x) ~ 1, dat),
    "target += weights[n] * multi_normal_cholesky_lpdf(Y[n] | Mu[n], LSigma);"
  )
  
  # models without residual correlations
  bform <- bf(y1 | cens(censi) ~ x + y2 + (x|2|g)) + 
    gaussian() + cor_ar() +
    (bf(x ~ 1) + mixture(poisson, nmix = 2)) +
    (bf(y2 ~ s(y2) + (1|2|g)) + skew_normal())
  bprior <- prior(normal(0, 5), resp = y1) +
    prior(normal(0, 10), resp = y2)
  scode <- make_stancode(bform, dat, prior = bprior)
  expect_match2(scode, "vector[N_1] r_1_y2_3 = r_1[, 3]")
  expect_match2(scode, "e_y1[n] = Y_y1[n] - mu_y1[n]")
  expect_match2(scode, "target += normal_lccdf(Y_y1[n] | mu_y1[n], sigma_y1)")
  expect_match2(scode, "target += skew_normal_lpdf(Y_y2 | mu_y2, omega_y2, alpha_y2)")
  expect_match2(scode, "ps[1] = log(theta1_x) + poisson_log_lpmf(Y_x[n] | mu1_x[n])")
  expect_match2(scode, "target += normal_lpdf(b_y1 | 0, 5)")
  expect_match2(scode, "target += normal_lpdf(bs_y2 | 0, 10)")
  
  # multivariate binomial models
  bform <- bf(x ~ 1) + bf(g ~ 1) + binomial()
  scode <- make_stancode(bform, dat)
  expect_match2(scode, "binomial_logit_lpmf(Y_x | trials_x, mu_x)")
  expect_match2(scode, "binomial_logit_lpmf(Y_g | trials_g, mu_g)")
  
  bform <- bform + weibull()
  scode <- make_stancode(bform, dat)
  expect_match2(scode, "mu_g[n] = exp(mu_g[n]) / tgamma(1 + 1 / shape_g)")
})

test_that("Stan code for categorical models is correct", {
  dat <- data.frame(y = rep(1:4, 2), x = 1:8, g = 1:8)
  prior <- c(
    prior(normal(0, 5), "b"),
    prior(normal(0, 10), "b", dpar = mu2),
    prior(cauchy(0, 1), "Intercept"),
    prior(normal(0, 2), "Intercept", dpar = mu3)
  )
  scode <- make_stancode(y ~ x + (1|ID|g), data = dat, 
                         family = categorical(), prior = prior)
  expect_match2(scode, "target += categorical_logit_lpmf(Y[n] | mu[n]);")
  expect_match2(scode, "mu2 = temp_mu2_Intercept + Xc_mu2 * b_mu2;")
  expect_match2(scode, "mu4[n] += r_1_mu4_3[J_1[n]] * Z_1_mu4_3[n];")
  expect_match2(scode, "target += normal_lpdf(b_mu2 | 0, 10);")
  expect_match2(scode, "target += normal_lpdf(b_mu4 | 0, 5);")
  expect_match2(scode, "target += cauchy_lpdf(temp_mu2_Intercept | 0, 1);")
  expect_match2(scode, "target += normal_lpdf(temp_mu3_Intercept | 0, 2);")
})

test_that("Stan code for ARMA models is correct", {
  dat <- data.frame(y = rep(1:4, 2), x = 1:8, time = 1:8)
  scode <- make_stancode(y ~ x, dat, student(), 
                         autocor = cor_ar(~time))
  expect_match2(scode, "e[n] = Y[n] - mu[n];")
  expect_match2(scode, "mu[n] += head(E[n], Kar) * ar;")
  
  scode <- make_stancode(y ~ x, dat, student(), 
                         autocor = cor_ma(~time, q = 2))
  expect_match2(scode, "mu[n] += head(E[n], Kma) * ma;")
  
  scode <- expect_warning(
    make_stancode(y ~ x, dat, student(), autocor = cor_arr(~time, r = 2)),
    "The 'arr' correlation structure has been deprecated"
  )
  expect_match2(scode, "mu = temp_Intercept + Xc * b + Yarr * arr;")
  
  scode <- make_stancode(cbind(y, x) ~ 1, dat, gaussian(),
                         autocor = cor_ar())
  expect_match2(scode, "e_y[n] = Y_y[n] - mu_y[n];")
})

test_that("Stan code for intercept only models is correct", {
  expect_match2(make_stancode(rating ~ 1, data = inhaler),
               "b_Intercept = temp_Intercept;") 
  expect_match2(make_stancode(rating ~ 1, data = inhaler, family = cratio()),
               "b_Intercept = temp_Intercept;") 
  expect_match2(make_stancode(rating ~ 1, data = inhaler, family = categorical()),
               "b_mu3_Intercept = temp_mu3_Intercept;")
})

test_that("Stan code of ordinal models is correct", {
  dat <- data.frame(y = c(rep(1:4, 2), 1, 1), x1 = rnorm(10), 
                    x2 = rnorm(10), g = rep(1:2, 5))
  
  scode <- make_stancode(
    y ~ x1, dat, family = cumulative(),
    prior = prior(normal(0, 2), Intercept, coef = 2)
  )
  expect_match2(scode, 
    "target += ordered_logistic_lpmf(Y[n] | mu[n], temp_Intercept);"
  )
  expect_match2(scode, "target += student_t_lpdf(temp_Intercept[1] | 3, 0, 10);")
  expect_match2(scode, "target += normal_lpdf(temp_Intercept[2] | 0, 2);")
  
  scode <- make_stancode(
    y ~ x1, dat, cumulative("probit", threshold = "equidistant"),
    prior = prior(normal(0, 2), Intercept, coef = 1)
  )
  expect_match2(scode, "real cumulative_probit_lpmf(int y")
  expect_match2(scode, "p = Phi(disc * (thres[1] - mu));")
  expect_match2(scode, "real<lower=0> delta;")
  expect_match2(scode, "temp_Intercept[k] = temp_Intercept1 + (k - 1.0) * delta;")
  expect_match2(scode, "b_Intercept = temp_Intercept + dot_product(means_X, b);")
  expect_match2(scode, "target += normal_lpdf(temp_Intercept1 | 0, 2);")
  
  scode <- make_stancode(y ~ x1, dat, family = cratio("probit_approx"))
  expect_match2(scode, "real cratio_probit_approx_lpmf(int y")
  expect_match2(scode, "q[k] = Phi_approx(disc * (mu - thres[k]));")

  scode <- make_stancode(y ~ x1 + cs(x2) + cs(g), dat, family = sratio())
  expect_match2(scode, "real sratio_logit_lpmf(int y")
  expect_match2(scode, "matrix[N, Kcs] Xcs;")
  expect_match2(scode, "matrix[Kcs, ncat - 1] bcs;")
  expect_match2(scode, "mucs = Xcs * bcs;")
  expect_match2(scode, 
    "target += sratio_logit_cs_lpmf(Y[n] | mu[n], mucs[n], temp_Intercept, disc);"
  )
  
  scode <- make_stancode(y ~ x1 + cse(x2) + (cse(1)|g), dat, family = acat())
  expect_match2(scode, "real acat_logit_lpmf(int y")
  expect_match2(scode, "p[k + 1] = p[k] + disc * (mu + mucs[k] - thres[k]);")
  expect_match2(scode, "mucs[n, 1] = mucs[n, 1] + r_1_1[J_1[n]] * Z_1_1[n];")
  expect_match2(scode, "b_Intercept = temp_Intercept + dot_product(means_X, b);")
  
  scode <- make_stancode(y ~ x1 + (cse(x2)||g), dat, family = acat("probit"))
  expect_match2(scode, "q[k] = Phi(disc * (mu + mucs[k] - thres[k]));")
  expect_match2(scode, 
    paste("mucs[n, 3] = mucs[n, 3] + r_1_3[J_1[n]] * Z_1_3[n]", 
          "+ r_1_6[J_1[n]] * Z_1_6[n];"))
  expect_match2(scode, 
    "target += acat_probit_cs_lpmf(Y[n] | mu[n], mucs[n], temp_Intercept, disc);"
  )
})

test_that("ordinal disc parameters appear in the Stan code", {
  scode <- make_stancode(
    bf(rating ~ period + carry + treat, disc ~ period),
    data = inhaler, family = cumulative(), 
    prior = prior(normal(0,5), dpar = disc)
  )
  expect_match2(scode, 
    "target += cumulative_logit_lpmf(Y[n] | mu[n], temp_Intercept, disc[n])"
  )
  expect_match2(scode, "target += normal_lpdf(b_disc | 0, 5)")
  expect_match2(scode, "disc[n] = exp(disc[n])")
})

test_that("monotonic effects appear in the Stan code", {
  prior <- c(prior(normal(0,1), class = b, coef = mox1),
             prior(dirichlet(c(1,0.5,2)), simo, coef = mox11),
             prior(dirichlet(c(1,0.5,2)), simo, coef = mox21))
  dat <- data.frame(y = rpois(120, 10), x1 = rep(1:4, 30), 
                    x2 = factor(rep(c("a", "b", "c"), 40), ordered = TRUE))
  scode <- make_stancode(y ~ y*mo(x1)*mo(x2), dat, prior = prior)
  expect_match2(scode, "int Xmo_3[N];")
  expect_match2(scode, "simplex[Jmo[1]] simo_1;")
  expect_match2(scode, "(bsp[2]) * mo(simo_2, Xmo_2[n])")
  expect_match2(scode, 
    "(bsp[6]) * mo(simo_7, Xmo_7[n]) * mo(simo_8, Xmo_8[n]) * Csp_3[n]"
  )
  expect_match2(scode, "target += normal_lpdf(bsp[1] | 0, 1)")
  expect_match2(scode, "target += dirichlet_lpdf(simo_1 | con_simo_1);")
  expect_match2(scode, "target += dirichlet_lpdf(simo_8 | con_simo_8);")
  
  scode <- make_stancode(y ~ mono(x1) + (mono(x1)|x2), dat)
  expect_match2(scode, "(bsp[1] + r_1_2[J_1[n]]) * mo(simo_1, Xmo_1[n])")
  expect_true(!grepl("Z_1_w", scode))
  
  expect_error(
    make_stancode(y ~ mono(x1) + (mono(x2)|x2), dat),
    "Special group-level terms require"
  )
  
  prior <- prior(beta(1, 1), simo, coef = mox11)
  expect_error(
    make_stancode(y ~ mo(x1), dat, prior = prior),
    "'dirichlet' is the only valid prior for simplex parameters"
  )
})

test_that("Stan code for non-linear models is correct", {
  flist <- list(a ~ x, b ~ z + (1|g))
  data <- data.frame(
    y = rgamma(9, 1, 1), x = rnorm(9), 
    z = rnorm(9), g = rep(1:3, 3)
  )
  prior <- c(set_prior("normal(0,5)", nlpar = "a"),
             set_prior("normal(0,1)", nlpar = "b"))
  # syntactic validity is already checked within make_stancode
  scode <- make_stancode(
    bf(y ~ a - exp(b^z) * (z <= a), flist = flist, nl = TRUE), 
    data = data, prior = prior
  )
  expect_match2(scode, 
    "mu[n] = nlp_a[n] - exp(nlp_b[n] ^ C_1[n]) * (C_1[n] <= nlp_a[n]);"
  )
  
  # non-linear predictor can be computed outside a loop
  scode <- make_stancode(bf(y ~ a - exp(b + z), flist = flist, 
                            nl = TRUE, loop = FALSE), 
                         data = data, prior = prior)
  expect_match2(scode, "\n  mu = nlp_a - exp(nlp_b + C_1);")
  
  flist <- list(a1 ~ 1, a2 ~ z + (x|g))
  prior <- c(set_prior("beta(1,1)", nlpar = "a1", lb = 0, ub = 1),
             set_prior("normal(0,1)", nlpar = "a2"))
  scode <- make_stancode(
    bf(y ~ a1 * exp(-x/(a2 + z)), 
       flist = flist, nl = TRUE),
    data = data, family = Gamma("log"), 
    prior = prior
  )
  expect_match2(scode,
    paste("mu[n] = shape * exp(-(nlp_a1[n] *", 
          "exp( - C_1[n] / (nlp_a2[n] + C_2[n]))));"))
  
  bform <- bf(y ~ x) + 
    nlf(sigma ~ a1 * exp(-x/(a2 + z))) +
    lf(a1 ~ 1, a2 ~ z + (x|g)) +
    lf(alpha ~ x)
  scode <- make_stancode(
    bform, data, family = skew_normal(),
    prior = c(
      prior(normal(0, 1), nlpar = a1),
      prior(normal(0, 5), nlpar = a2)
    )
  )
  expect_match2(scode, "nlp_a1 = X_a1 * b_a1")
  expect_match2(scode,
    "sigma[n] = exp(nlp_a1[n] * exp( - C_sigma_1[n] / (nlp_a2[n] + C_sigma_2[n])))"
  )
  expect_match2(scode, "target += normal_lpdf(b_a2 | 0, 5)")
  
  expect_error(make_stancode(bform, data, family = skew_normal()),
               "Priors on population-level effects are required")
})

test_that("Stan code for nested non-linear parameters is correct", {
  dat <- data.frame(y = rnorm(10), x = rnorm(10), z = 1:5)
  bform <- bf(
    y ~ lb + (1 - lb) * inv_logit(b * x),
    b + a ~ 1 + (1 | z), nlf(lb ~ inv_logit(a / x)),
    nl = TRUE
  )
  bprior <- prior(normal(0, 1), nlpar = "a") +
    prior(normal(0, 1), nlpar = "b")
  scode <- make_stancode(bform, dat, prior = bprior)
  expect_match2(scode, "nlp_lb[n] = inv_logit(nlp_a[n] / C_lb_1[n]);")
  expect_match2(scode, 
    "mu[n] = nlp_lb[n] + (1 - nlp_lb[n]) * inv_logit(nlp_b[n] * C_1[n]);"
  )
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
                     "mu[n] = (mu[n]);", fixed = TRUE))
  expect_true(!grepl(make_stancode(time ~ age, data = kidney, 
                                   family = poisson("identity")), 
                     "mu[n] = (mu[n]);", fixed = TRUE))
})

test_that("known standard errors appear in the Stan code", {
  scode <- make_stancode(time | se(age) ~ sex, data = kidney)
  expect_match2(scode, "target += normal_lpdf(Y | mu, se)")
  scode <- make_stancode(time | se(age) + weights(age) ~ sex, data = kidney)
  expect_match2(scode, "target += weights[n] * normal_lpdf(Y[n] | mu[n], se[n])")
  scode <- make_stancode(time | se(age, sigma = TRUE) ~ sex, data = kidney)
  expect_match2(scode, "target += normal_lpdf(Y | mu, sqrt(square(sigma) + se2))")
  scode <- make_stancode(bf(time | se(age, sigma = TRUE) ~ sex, sigma ~ sex), 
                         data = kidney)
  expect_match2(scode, "target += normal_lpdf(Y | mu, sqrt(square(sigma) + se2))")
})

test_that("functions defined in 'stan_funs' appear in the functions block", {
  test_fun <- paste0("  real test_fun(real a, real b) { \n",
                     "    return a + b; \n",
                     "  } \n")
  scode <- SW(make_stancode(time ~ age, data = kidney, stan_funs = test_fun))
  expect_match2(scode, test_fun)
})

test_that("fixed residual covariance matrices appear in the Stan code", {
  data <- data.frame(y = 1:5)
  V <- diag(5)
  expect_match2(make_stancode(y~1, data = data, family = gaussian(), 
                             autocor = cor_fixed(V)),
               "target += multi_normal_cholesky_lpdf(Y | mu, LV)")
  expect_match2(make_stancode(y~1, data = data, family = student(),
                             autocor = cor_fixed(V)),
               "target += multi_student_t_lpdf(Y | nu, mu, V)")
  expect_match2(make_stancode(y~1, data = data, family = student(),
                             autocor = cor_fixed(V)),
               "target += multi_student_t_lpdf(Y | nu, mu, V)")
})

test_that("Stan code for BSTS models is correct", {
  dat <- data.frame(y = rnorm(10), x = rnorm(10))
  scode <- SW(make_stancode(
    y ~ x, data = dat, autocor = cor_bsts(),
    prior = prior(normal(0, 10), sigmaLL)
  ))
  expect_match2(scode, "mu[n] += loclev[n]")
  expect_match2(scode, "target += normal_lpdf(loclev[n] | loclev[n - 1], sigmaLL)")
  expect_match2(scode, "target += normal_lpdf(sigmaLL | 0, 10)")
  
  dat <- data.frame(y = rexp(1), x = rnorm(10))
  scode <- SW(make_stancode(
    y~x, data = dat, family = student("log"),
    autocor = cor_bsts()
  ))
  expect_match2(scode, "target += normal_lpdf(loclev[n] | log(Y[n]), sigmaLL)")
})

test_that("Stan code for GAMMs is correct", {
  dat <- data.frame(y = rnorm(10), x = rnorm(10), g = factor(rep(1:2, 5)))
  scode <- make_stancode(y ~ s(x) + (1|g), data = dat,
                         prior = set_prior("normal(0,2)", "sds"))
  expect_match2(scode, "Zs_1_1 * s_1_1")
  expect_match2(scode, "matrix[N, knots_1[1]] Zs_1_1")
  expect_match2(scode, "target += normal_lpdf(zs_1_1 | 0, 1)")
  expect_match2(scode, "target += normal_lpdf(sds_1_1 | 0,2)")
  
  prior <- c(set_prior("normal(0,5)", nlpar = "lp"),
             set_prior("normal(0,2)", "sds", nlpar = "lp"))
  scode <- make_stancode(bf(y ~ lp, lp ~ s(x) + (1|g), nl = TRUE), 
                         data = dat, prior = prior)
  expect_match2(scode, "Zs_lp_1_1 * s_lp_1_1")
  expect_match2(scode, "matrix[N, knots_lp_1[1]] Zs_lp_1_1")
  expect_match2(scode, "target += normal_lpdf(zs_lp_1_1 | 0, 1)")
  expect_match2(scode, "target += normal_lpdf(sds_lp_1_1 | 0,2)")
  
  scode <- make_stancode(y ~ s(x) + t2(x,y), data = dat,
                        prior = set_prior("normal(0,2)", "sds"))
  expect_match2(scode, "Zs_2_2 * s_2_2")
  expect_match2(scode, "matrix[N, knots_2[2]] Zs_2_2")
  expect_match2(scode, "target += normal_lpdf(zs_2_2 | 0, 1)")
  expect_match2(scode, "target += normal_lpdf(sds_2_2 | 0,2)")
  
  scode <- make_stancode(y ~ g + s(x, by = g), data = dat)
  expect_match2(scode, "vector[knots_2[1]] zs_2_1")
  expect_match2(scode, "s_2_1 = sds_2_1 * zs_2_1")
})

test_that("Stan code of response times models is correct", {
  dat <- epilepsy
  dat$cens <- sample(-1:1, nrow(dat), TRUE)
  scode <- make_stancode(count ~ Trt_c + (1|patient),
                      data = dat, family = exgaussian("log"),
                      prior = prior(gamma(1,1), class = beta))
  expect_match2(scode,
    "target += exp_mod_normal_lpdf(Y | mu - beta, sigma, inv(beta))"
  )
  expect_match2(scode, "mu[n] = exp(mu[n])")
  expect_match2(scode, "target += gamma_lpdf(beta | 1, 1)")
  
  scode <- make_stancode(bf(count ~ Trt_c + (1|patient),
                         sigma ~ Trt_c, beta ~ Trt_c),
                      data = dat, family = exgaussian())
  expect_match2(scode, 
    "target += exp_mod_normal_lpdf(Y | mu - beta, sigma, inv(beta))"
  )
  expect_match2(scode, "beta[n] = exp(beta[n])")
  
  scode <- make_stancode(count | cens(cens) ~ Trt_c + (1|patient),
                      data = dat, family = exgaussian("inverse"))
  expect_match2(scode, "exp_mod_normal_lccdf(Y[n] | mu[n] - beta, sigma, inv(beta))")
  
  scode <- make_stancode(count ~ Trt_c, dat, family = shifted_lognormal())
  expect_match2(scode, "target += lognormal_lpdf(Y - ndt | mu, sigma)")
  
  scode <- make_stancode(count | cens(cens) ~ Trt_c, dat, family = shifted_lognormal())
  expect_match2(scode, "target += lognormal_lcdf(Y[n] - ndt | mu[n], sigma)")
})

test_that("Stan code of wiener diffusion models is correct", {
  dat <- data.frame(q = 1:10, resp = sample(0:1, 10, TRUE), x = rnorm(10))
  scode <- make_stancode(q | dec(resp) ~ x, data = dat, family = wiener())
  expect_match2(scode, 
    "target += wiener_diffusion_lpdf(Y[n] | dec[n], bs, ndt, bias, mu[n])"
  )
  
  scode <- make_stancode(bf(q | dec(resp) ~ x, bs ~ x, ndt ~ x, bias ~ x), 
                         data = dat, family = wiener())
  expect_match2(scode,
    "target += wiener_diffusion_lpdf(Y[n] | dec[n], bs[n], ndt[n], bias[n], mu[n])"
  )
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
  expect_match2(scode, "r_2_a_2 = r_2[, 2];")
  expect_match2(scode, "r_1_sigma_2 = sd_1[2] * (z_1[2]);")
})

test_that("distributional gamma models are handled correctly", {
  # test fix of issue #124
  scode <- make_stancode(
    bf(time ~ age * sex + disease + (1|patient), 
       shape ~ age + (1|patient)), 
    data = kidney, family = Gamma("log")
  )
  expect_match2(scode, paste0(
    "    shape[n] = exp(shape[n]); \n", 
    "    mu[n] = shape[n] * exp(-(mu[n]));"))
  
  scode <- make_stancode(
    bf(time ~ inv_logit(a) * exp(b * age),
       a + b ~ sex + (1|patient), nl = TRUE, 
       shape ~ age + (1|patient)), 
    data = kidney, family = Gamma("identity"),
    prior = c(set_prior("normal(2,2)", nlpar = "a"),
              set_prior("normal(0,3)", nlpar = "b"))
  )
  expect_match2(scode, paste0(
    "    shape[n] = exp(shape[n]); \n", 
    "    // compute non-linear predictor \n",
    "    mu[n] = shape[n] / (inv_logit(nlp_a[n]) * exp(nlp_b[n] * C_1[n]));"
  ))
  
  scode <- make_stancode(
    bf(time ~ age, shape ~ age), 
    data = kidney,
    family = brmsfamily("gamma", link_shape = "identity")
  )
  # test that no link function is applied on 'shape'
  expect_match2(scode, paste0(
    "  for (n in 1:N) { \n",
    "    mu[n] = shape[n] * exp(-(mu[n])); \n",
    "  } \n"
  ))
})

test_that("weighted, censored, and truncated likelihoods are correct", {
  dat <- data.frame(y = 1:9, x = rep(-1:1, 3), y2 = 10:18)
  
  scode <- make_stancode(y | weights(y2) ~ 1, dat, poisson())
  expect_match2(scode, "target += weights[n] * poisson_log_lpmf(Y[n] | mu[n]);")
  
  scode <- make_stancode(y | trials(y2) + weights(y2) ~ 1, dat, binomial())
  expect_match2(scode, 
    "target += weights[n] * binomial_logit_lpmf(Y[n] | trials[n], mu[n]);"
  )
  
  expect_match2(make_stancode(y | cens(x, y2) ~ 1, dat, poisson()),
                "target += poisson_lpmf(Y[n] | mu[n]);")
  expect_match2(make_stancode(y | cens(x) ~ 1, dat, weibull()), 
                "target += weibull_lccdf(Y[n] | shape, mu[n]);")
  dat$x[1] <- 2
  scode <- make_stancode(y | cens(x, y2) ~ 1, dat, gaussian())
  expect_match2(scode, paste0(
    "target += log_diff_exp(\n", 
    "          normal_lcdf(rcens[n] | mu[n], sigma),"
  ))
  dat$x <- 1
  expect_match2(make_stancode(y | cens(x) + weights(x) ~ 1, dat, weibull()),
   "target += weights[n] * weibull_lccdf(Y[n] | shape, mu[n]);")
  
  scode <- make_stancode(y | cens(x) + trunc(0.1) ~ 1, dat, weibull())
  expect_match2(scode, "target += weibull_lccdf(Y[n] | shape, mu[n]) -")
  expect_match2(scode, "  weibull_lccdf(lb[n] | shape, mu[n]);")
  
  scode <- make_stancode(y | cens(x) + trunc(ub = 30) ~ 1, dat)
  expect_match2(scode, "target += normal_lccdf(Y[n] | mu[n], sigma) -")
  expect_match2(scode, "  normal_lcdf(ub[n] | mu[n], sigma);")
  
  scode <- make_stancode(y | weights(x) + trunc(0, 30) ~ 1, dat)
  expect_match2(scode, "target += weights[n] * normal_lpdf(Y[n] | mu[n], sigma) -")
  expect_match2(scode, "  log_diff_exp(normal_lcdf(ub[n] | mu[n], sigma),")
})

test_that("noise-free terms appear in the Stan code", {
  N <- 30
  dat <- data.frame(
    y = rnorm(N), x = rnorm(N), z = rnorm(N),
    xsd = abs(rnorm(N, 1)), zsd = abs(rnorm(N, 1)),
    ID = rep(1:5, each = N / 5)
  )
  me_prior <- prior(normal(0,5)) + 
    prior(normal(0, 10), "meanme") +
    prior(cauchy(0, 5), "sdme", coef = "mez") +
    prior(lkj(2), "corme")
  scode <- make_stancode(
    y ~ me(x, xsd)*me(z, zsd)*x, data = dat, prior = me_prior,
    sample_prior = "yes"
  )
  expect_match2(scode, 
    "(bsp[1]) * Xme_1[n] + (bsp[2]) * Xme_2[n] + (bsp[3]) * Xme_1[n] * Xme_2[n]"
  )
  expect_match2(scode, "(bsp[6]) * Xme_1[n] * Xme_2[n] * Csp_3[n]")
  expect_match2(scode, "target += normal_lpdf(Xn_2 | Xme_2, noise_2)")
  expect_match2(scode, "target += normal_lpdf(bsp | 0, 5)")
  expect_match2(scode, "target += normal_lpdf(to_vector(zme_1) | 0, 1)")
  expect_match2(scode, "target += normal_lpdf(meanme_1 | 0, 10)")
  expect_match2(scode, "target += cauchy_lpdf(sdme_1[2] | 0, 5)")
  expect_match2(scode, "target += lkj_corr_cholesky_lpdf(Lme_1 | 2)")
  expect_match2(scode, "+ (diag_pre_multiply(sdme_1, Lme_1) * zme_1)'")
  
  scode <- make_stancode(
    y ~ me(x, xsd)*z + (me(x, xsd)*z|ID), data = dat
  )
  expect_match2(scode, "(bsp[1] + r_1_3[J_1[n]]) * Xme_1[n]")
  expect_match2(scode, "(bsp[2] + r_1_4[J_1[n]]) * Xme_1[n] * Csp_1[n]")
  
  expect_match2(make_stancode(y ~ I(me(x, xsd)^2), data = dat),
               "(bsp[1]) * (Xme_1[n]^2)")
  
  # test that noise-free variables are unique across model parts
  scode <- make_stancode(
    bf(y ~ me(x, xsd)*me(z, zsd)*x, sigma ~ me(x, xsd)), 
    data = dat, prior = prior(normal(0,5))
  )
  expect_match2(scode, "mu[n] += (bsp[1]) * Xme_1[n]")
  expect_match2(scode, "sigma[n] += (bsp_sigma[1]) * Xme_1[n]")
  
  scode <- make_stancode(
    bf(y ~ a * b, a + b ~ me(x, xsd), nl = TRUE), 
    data = dat, 
    prior = prior(normal(0,5), nlpar = a) + 
      prior(normal(0, 5), nlpar = b)
  )
  expect_match2(scode, "nlp_a[n] += (bsp_a[1]) * Xme_1[n]")
  expect_match2(scode, "nlp_b[n] += (bsp_b[1]) * Xme_1[n]")
  
  bform <- bf(cbind(y, z) ~ me(x, xsd)) + set_mecor(FALSE)
  scode <- make_stancode(cbind(y, z) ~ me(x, xsd), dat)
  expect_match2(scode, "mu_y[n] += (bsp_y[1]) * Xme_1[n]")
  expect_match2(scode, "mu_z[n] += (bsp_z[1]) * Xme_1[n]")
  expect_match2(scode, "vector[N] Xme_1 = meanme_1[1] + sdme_1[1] * zme_1;")
  
  # noise-free terms with grouping factors
  bform <- bf(y ~ me(x, xsd, ID) + (me(x, xsd, ID) | ID))
  scode <- make_stancode(bform, dat)
  expect_match2(scode, "vector[Nme_1] Xn_1;")
  expect_match2(scode, "(bsp[1] + r_1_2[J_1[n]]) * Xme_1[Jme_1[n]]")
  
  bform <- bform + set_mecor(FALSE)
  scode <- make_stancode(bform, dat)
  expect_match2(scode, "vector[Nme_1] Xme_1 = meanme_1[1] + sdme_1[1] * zme_1;")
})

test_that("Stan code of multi-membership models is correct", {
  dat <- data.frame(y = rnorm(10), g1 = sample(1:10, 10, TRUE),
                    g2 = sample(1:10, 10, TRUE), w1 = rep(1, 10),
                    w2 = rep(abs(rnorm(10))))
  expect_match2(make_stancode(y ~ (1|mm(g1, g2)), data = dat), 
    paste0(" W_1_1[n] * r_1_1[J_1_1[n]] * Z_1_1_1[n]",
           " + W_1_2[n] * r_1_1[J_1_2[n]] * Z_1_1_2[n]")
  )
  expect_match2(make_stancode(y ~ (1+w1|mm(g1,g2)), data = dat), 
    paste0(" W_1_1[n] * r_1_2[J_1_1[n]] * Z_1_2_1[n]",
           " + W_1_2[n] * r_1_2[J_1_2[n]] * Z_1_2_2[n]")
  )
  expect_match2(make_stancode(y ~ (1+mmc(w1, w2)|mm(g1,g2)), data = dat), 
    " W_1_2[n] * r_1_2[J_1_2[n]] * Z_1_2_2[n];"
  )
})

test_that("by variables in grouping terms are handled correctly", {
  dat <- data.frame(
    y = rnorm(100), x = rnorm(100),
    g = rep(1:10, each = 10),
    z = factor(rep(c(0, 4.5, 3, 2, 5), each = 20))
  )
  scode <- make_stancode(y ~ x + (1 | gr(g, by = z)), dat)
  expect_match2(scode, "r_1_1 = sd_1[1, Jby_1]' .* (z_1[1]);")
  scode <- make_stancode(y ~ x + (x | gr(g, by = z)), dat)
  expect_match2(scode, "r_1 = scale_r_cor_by(z_1, sd_1, L_1, Jby_1);")
  expect_match2(scode, "target += student_t_lpdf(to_vector(sd_1) | 3, 0, 10);")
  expect_match2(scode, "target += lkj_corr_cholesky_lpdf(L_1[5] | 1);")
})

test_that("Group syntax | and || is handled correctly,", {
  data <- data.frame(y = rnorm(10), x = rnorm(10),
                     g1 = rep(1:5, each = 2), g2 = rep(1:2, 5))
  scode <- make_stancode(y ~ x + (1+x||g1) + (I(x/4)|g2), data)
  expect_match2(scode, "r_1_2 = sd_1[2] * (z_1[2]);")
  expect_match2(scode, "r_2_1 = r_2[, 1];")
  expect_match2(scode, "r_2 = (diag_pre_multiply(sd_2, L_2) * z_2)';")
})

test_that("predicting zi and hu works correctly", {
  scode <- make_stancode(bf(count ~ Trt_c, zi ~ Trt_c), epilepsy, 
                         family = "zero_inflated_poisson")
  expect_match2(scode, 
    "target += zero_inflated_poisson_log_logit_lpmf(Y[n] | mu[n], zi[n])"
  )
  expect_true(!grepl("inv_logit\\(", scode))
  expect_true(!grepl("exp(mu[n])", scode, fixed = TRUE))
  
  scode <- make_stancode(bf(count ~ Trt_c, zi ~ Trt_c), epilepsy, 
                         family = zero_inflated_poisson(identity))
  expect_match2(scode, 
    "target += zero_inflated_poisson_logit_lpmf(Y[n] | mu[n], zi[n])"
  )
  
  scode <- make_stancode(bf(count ~ Trt_c, zi ~ Trt_c), epilepsy, 
                         family = "zero_inflated_binomial")
  expect_match2(scode, 
    "target += zero_inflated_binomial_blogit_logit_lpmf(Y[n] | trials[n], mu[n], zi[n])"
  )
  expect_true(!grepl("inv_logit\\(", scode))
  
  fam <- zero_inflated_binomial("probit", link_zi = "identity")
  scode <- make_stancode(bf(count ~ Trt_c, zi ~ Trt_c), epilepsy, 
                         family = fam)
  expect_match2(scode, 
    "target += zero_inflated_binomial_lpmf(Y[n] | trials[n], mu[n], zi[n])"
  )
  expect_match2(scode, "mu[n] = Phi(mu[n]);")
  
  scode <- make_stancode(
    bf(count ~ Trt_c, zi ~ Trt_c), epilepsy, 
    family = zero_inflated_beta()
  )
  expect_match2(scode,
    "target += zero_inflated_beta_logit_lpdf(Y[n] | mu[n], phi, zi[n])"      
  )
  
  scode <- make_stancode(bf(count ~ Trt_c, hu ~ Trt_c), epilepsy, 
                         family = "hurdle_negbinomial")
  expect_match2(scode, 
    "target += hurdle_neg_binomial_log_logit_lpmf(Y[n] | mu[n], shape, hu[n])"
  )
  expect_true(!grepl("inv_logit\\(", scode))
  expect_true(!grepl("exp(mu[n])", scode, fixed = TRUE))
  
  scode <- make_stancode(bf(count ~ Trt_c, hu ~ Trt_c), epilepsy, 
                         family = "hurdle_gamma")
  expect_match2(scode, 
    "target += hurdle_gamma_logit_lpdf(Y[n] | shape, mu[n], hu[n])"
  )
  expect_true(!grepl("inv_logit\\(", scode))
  expect_match2(scode, "mu[n] = shape * exp(-(mu[n]));")
  
  scode <- make_stancode(
    bf(count ~ Trt_c, hu ~ Trt_c), epilepsy, 
    family = hurdle_gamma(link_hu = "identity")
  )
  expect_match2(scode, "target += hurdle_gamma_lpdf(Y[n] | shape, mu[n], hu[n])")
  expect_true(!grepl("inv_logit\\(", scode))
  expect_match2(scode, "mu[n] = shape * exp(-(mu[n]));")
})

test_that("fixing auxiliary parameters is possible", {
  scode <- make_stancode(bf(y ~ 1, sigma = 0.5), data = list(y = rnorm(10)))
  expect_match(scode, "data \\{[^\\}]*real<lower=0> sigma;")
})

test_that("Stan code of quantile regression models is correct", {
  data <- data.frame(y = rnorm(10), x = rnorm(10), c = 1)
  scode <- make_stancode(y ~ x, data, family = asym_laplace())
  expect_match2(scode, "target += asym_laplace_lpdf(Y[n] | mu[n], sigma, quantile)")
  
  scode <- make_stancode(bf(y ~ x, quantile = 0.75), data, family = asym_laplace())
  expect_match(scode, "data \\{[^\\}]*real<lower=0,upper=1> quantile;")
  
  scode <- make_stancode(y | cens(c) ~ x, data, family = asym_laplace())
  expect_match2(scode, "target += asym_laplace_lccdf(Y[n] | mu[n], sigma, quantile)")
  
  scode <- make_stancode(bf(y ~ x, sigma ~ x), data, family = asym_laplace())
  expect_match2(scode, "target += asym_laplace_lpdf(Y[n] | mu[n], sigma[n], quantile)")
})

test_that("Stan code of GEV models is correct", {
  data <- data.frame(y = rnorm(10), x = rnorm(10), c = 1)
  scode <- make_stancode(y ~ x, data, gen_extreme_value())
  expect_match2(scode, "target += gen_extreme_value_lpdf(Y[n] | mu[n], sigma, xi)")
  expect_match2(scode, "xi = scale_xi(temp_xi, Y, mu, sigma)")
  
  scode <- make_stancode(bf(y ~ x, sigma ~ x), data, gen_extreme_value())
  expect_match2(scode, "xi = scale_xi_vector(temp_xi, Y, mu, sigma)")
  
  scode <- make_stancode(bf(y ~ x, xi ~ x), data, gen_extreme_value())
  expect_match2(scode, "xi[n] = expm1(xi[n])")
  
  scode <- make_stancode(bf(y ~ x, xi = 0), data, gen_extreme_value())
  expect_match(scode, "data \\{[^\\}]*real xi;  // shape parameter")
  
  scode <- make_stancode(y | cens(c) ~ x, data, gen_extreme_value())
  expect_match2(scode, "target += gen_extreme_value_lccdf(Y[n] | mu[n], sigma, xi)")
})

test_that("offsets appear in the Stan code", {
  data <- data.frame(y = rnorm(10), x = rnorm(10), c = 1)
  scode <- make_stancode(y ~ x + offset(c), data)
  expect_match2(scode, "temp_Intercept + Xc * b + offset;")
  scode <- make_stancode(bf(y ~ a, a ~ offset(log(c + 1)), nl = TRUE),
                         data, prior = prior(normal(0,1), nlpar = a))
  expect_match2(scode, "X_a * b_a + offset_a;")
})

test_that("prior only models are correctly checked", {
  data <- data.frame(y = rnorm(10), x = rnorm(10), c = 1)
  prior <- prior(normal(0, 5), b) + prior("", Intercept)
  expect_error(make_stancode(y ~ x, data, prior = prior,
                             sample_prior = "only"),
               "Sampling from priors is not possible")
  prior <- prior(normal(0, 5), b) + prior(normal(0, 10), Intercept)
  scode <- make_stancode(y ~ x, data, prior = prior,
                         sample_prior = "only")
  expect_match2(scode, "target += normal_lpdf(temp_Intercept | 0, 10)")
})

test_that("Stan code of mixture model is correct", {
  data <- data.frame(y = 1:10, x = rnorm(10), c = 1)
  scode <- make_stancode(
    bf(y ~ x,  sigma2 ~ x), data, 
    family = mixture(gaussian, gaussian),
    sample_prior = TRUE
  )
  expect_match2(scode, "ordered[2] ordered_Intercept;")
  expect_match2(scode, "temp_mu2_Intercept = ordered_Intercept[2];")
  expect_match2(scode, "target += dirichlet_lpdf(theta | con_theta);")
  expect_match2(scode, "ps[1] = log(theta1) + normal_lpdf(Y[n] | mu1[n], sigma1);")
  expect_match2(scode, "ps[2] = log(theta2) + normal_lpdf(Y[n] | mu2[n], sigma2[n]);")
  expect_match2(scode, "target += log_sum_exp(ps);")
  expect_match2(scode, "simplex[2] prior_theta = dirichlet_rng(con_theta);")
  
  data$z <- abs(data$y)
  scode <- make_stancode(bf(z | weights(c) ~ x, shape1 ~ x, theta1 = 1, theta2 = 2), 
                         data = data, mixture(Gamma("log"), weibull))
  expect_match(scode, "data \\{[^\\}]*real<lower=0,upper=1> theta1;")
  expect_match(scode, "data \\{[^\\}]*real<lower=0,upper=1> theta2;")
  expect_match(scode, "shape1\\[n\\] = exp\\(shape1\\[n\\]\\); \\\n    mu1\\[n\\] = ")
  expect_match2(scode, "ps[1] = log(theta1) + gamma_lpdf(Y[n] | shape1[n], mu1[n]);")
  expect_match2(scode, "target += weights[n] * log_sum_exp(ps);")
  
  scode <- make_stancode(bf(abs(y) | se(c) ~ x), data = data, 
                         mixture(gaussian, student))
  expect_match2(scode, "ps[1] = log(theta1) + normal_lpdf(Y[n] | mu1[n], se[n]);")
  expect_match2(scode, "ps[2] = log(theta2) + student_t_lpdf(Y[n] | nu2, mu2[n], se[n]);")
  
  fam <- mixture(gaussian, student, exgaussian)
  scode <- make_stancode(bf(y ~ x), data = data, family = fam)
  expect_match(scode, "parameters \\{[^\\}]*real temp_mu3_Intercept;")
  expect_match2(scode, 
    "ps[2] = log(theta2) + student_t_lpdf(Y[n] | nu2, mu2[n], sigma2);"
  )
  expect_match2(scode, 
    "ps[3] = log(theta3) + exp_mod_normal_lpdf(Y[n] | mu3[n] - beta3, sigma3, inv(beta3));"
  )
  
  scode <- make_stancode(bf(y ~ x, theta1 ~ x, theta3 ~ x), 
                         data = data, family = fam)
  expect_match2(scode, "log_sum_exp_theta = log(exp(theta1[n]) + exp(theta2[n]) + exp(theta3[n]));")
  expect_match2(scode, "theta2 = rep_vector(0, N);")
  expect_match2(scode, "theta3[n] = theta3[n] - log_sum_exp_theta;")
  expect_match2(scode, "ps[1] = theta1[n] + normal_lpdf(Y[n] | mu1[n], sigma1);")
  
  fam <- mixture(cumulative, sratio)
  scode <- make_stancode(y ~ x, data, family = fam)
  expect_match2(scode, "ordered_logistic_lpmf(Y[n] | mu1[n], temp_mu1_Intercept);")
  expect_match2(scode, "sratio_logit_lpmf(Y[n] | mu2[n], temp_mu2_Intercept, disc2);")
  
  # censored mixture model
  fam <- mixture(gaussian, gaussian)
  scode <- make_stancode(y | cens(2, y2 = 2) ~ x, data, fam)
  expect_match2(scode, 
    "ps[2] = log(theta2) + normal_lccdf(Y[n] | mu2[n], sigma2);"
  )
  expect_match2(scode, paste0(
    "ps[2] = log(theta2) + log_diff_exp(\n",
    "          normal_lcdf(rcens[n] | mu2[n], sigma2),"
  ))
  
  # truncated mixture model
  scode <- make_stancode(y | trunc(3) ~ x, data, fam)
  expect_match2(scode, paste0(
    "ps[1] = log(theta1) + normal_lpdf(Y[n] | mu1[n], sigma1) -\n", 
    "        normal_lccdf(lb[n] | mu1[n], sigma1);"
  ))
  
  # non-linear mixture model
  bform <- bf(y ~ 1) + 
    nlf(mu1 ~ eta^2) +
    nlf(mu2 ~ log(eta) + a) +
    lf(eta + a ~ x) +
    mixture(gaussian, nmix = 2)
  bprior <- prior(normal(0, 1), nlpar = "eta") + 
    prior(normal(0, 1), nlpar = "a")
  scode <- make_stancode(bform, data = data, prior = bprior)
  expect_match2(scode, "mu1[n] = nlp_eta[n] ^ 2;")
  expect_match2(scode, "mu2[n] = log(nlp_eta[n]) + nlp_a[n];")
})

test_that("sparse matrix multiplication is applied correctly", {
  data <- data.frame(y = rnorm(10), x = rnorm(10))
  # linear model
  scode <- make_stancode(
    bf(y ~ x, sigma ~ x), data, sparse = TRUE,
    prior = prior(normal(0, 5), coef = "Intercept")
  )
  expect_match2(scode, "wX = csr_extract_w(X);")
  expect_match2(scode, 
    "mu = csr_matrix_times_vector(rows(X), cols(X), wX, vX, uX, b);"
  )
  expect_match2(scode, 
    "uX_sigma[size(csr_extract_u(X_sigma))] = csr_extract_u(X_sigma);"
  )
  expect_match2(scode,
    paste0(
      "sigma = csr_matrix_times_vector(rows(X_sigma), cols(X_sigma), ", 
      "wX_sigma, vX_sigma, uX_sigma, b_sigma);"
    )
  )
  expect_match2(scode, "target += normal_lpdf(b[1] | 0, 5);")
  # non-linear model
  scode <- make_stancode(
    bf(y ~ a, a ~ x, nl = TRUE), 
    data, sparse = TRUE, 
    prior = prior(normal(0, 1), nlpar = a)
  )
  expect_match2(scode, 
    "vX_a[size(csr_extract_v(X_a))] = csr_extract_v(X_a);"
  )
  expect_match2(scode, 
    "nlp_a = csr_matrix_times_vector(rows(X_a), cols(X_a), wX_a, vX_a, uX_a, b_a);"
  )
})

test_that("Stan code for Gaussian processes is correct", {
  dat <- data.frame(y = rnorm(40), x1 = rnorm(40), x2 = rnorm(40),
                    z = factor(rep(3:6, each = 10)))
  
  prior <- prior(gamma(0.1, 0.1), sdgp)
  scode <- make_stancode(y ~ gp(x1) + gp(x2, by = x1), dat, prior = prior)
  expect_match2(scode, "target += inv_gamma_lpdf(lscale_1")
  expect_match2(scode, "target += gamma_lpdf(sdgp_1 | 0.1, 0.1)")
  expect_match2(scode, "Cgp_2 .* gp(Xgp_2, sdgp_2[1], lscale_2[1], zgp_2)")
  
  prior <- prior + prior(normal(0, 1), lscale, coef = gpx1)
  scode <- make_stancode(y ~ gp(x1) + gp(x2, by = x1, gr = TRUE), 
                         data = dat, prior = prior)
  expect_match2(scode, "target += normal_lpdf(lscale_1[1] | 0, 1)")
  expect_match2(scode, "+ Cgp_2 .* gp(Xgp_2, sdgp_2[1], lscale_2[1], zgp_2)[Jgp_2]")
  
  # non-isotropic GP
  scode <- make_stancode(y ~ gp(x1, x2, by = z, iso = FALSE), dat)
  expect_match2(scode, "target += inv_gamma_lpdf(lscale_1[1][2]")
  expect_match2(scode, "target += inv_gamma_lpdf(lscale_1[4][2]")
  
  # Suppress Stan parser warnings that can currently not be avoided
  scode <- make_stancode(y ~ gp(x1, x2) + gp(x1, by = z), 
                         dat, silent = TRUE)
  expect_match2(scode, "gp(Xgp_1, sdgp_1[1], lscale_1[1], zgp_1)")
  expect_match2(scode, paste0(
    "mu[Igp_2_2] = mu[Igp_2_2] + Cgp_2_2 .* gp(Xgp_2_2, ", 
    "sdgp_2[2], lscale_2[2], zgp_2_2);"
  ))
  
  # approximate GPS
  scode <- make_stancode(y ~ gp(x1, k = 10) + gp(x2, by = x1, k = 10), dat)
  expect_match2(scode, "target += inv_gamma_lpdf(lscale_1")
  expect_match2(scode, 
    "gpa(Xgp_1, sdgp_1[1], lscale_1[1], zgp_1, slambda_1)"           
  )
  expect_match2(scode, 
    "Cgp_2 .* gpa(Xgp_2, sdgp_2[1], lscale_2[1], zgp_2, slambda_2)"
  )
  
  prior <- c(prior(normal(0, 10), lscale, coef = gpx1, nlpar = a),
             prior(gamma(0.1, 0.1), sdgp, nlpar = a),
             prior(normal(0, 1), b, nlpar = a))
  scode <- make_stancode(bf(y ~ a, a ~ gp(x1), nl = TRUE), 
                         data = dat, prior = prior)
  expect_match2(scode, "target += normal_lpdf(lscale_a_1[1] | 0, 10)")
  expect_match2(scode, "target += gamma_lpdf(sdgp_a_1 | 0.1, 0.1)")
  expect_match2(scode, "gp(Xgp_a_1, sdgp_a_1[1], lscale_a_1[1], zgp_a_1)")
  
  prior <- prior(gamma(2, 2), lscale, coef = gpx1z5, nlpar = "a")
  scode <- make_stancode(bf(y ~ a, a ~ gp(x1, by = z, gr = TRUE), nl = TRUE),
                         data = dat, prior = prior, silent = TRUE)
  expect_match2(scode, 
    "nlp_a[Igp_a_1_1] = nlp_a[Igp_a_1_1] + Cgp_a_1_1 .* gp(Xgp_a_1_1,"
  )
  expect_match2(scode,
    "gp(Xgp_a_1_3, sdgp_a_1[3], lscale_a_1[3], zgp_a_1_3)[Jgp_a_1_3]"             
  )
  expect_match2(scode, "target += gamma_lpdf(lscale_a_1[3] | 2, 2);")
  expect_match2(scode, "target += normal_lpdf(zgp_a_1_3 | 0, 1);")
})

test_that("Stan code for SAR models is correct", {
  data(oldcol, package = "spdep")
  scode <- make_stancode(CRIME ~ INC + HOVAL, data = COL.OLD, 
                         autocor = cor_lagsar(COL.nb),
                         prior = prior(normal(0.5, 1), lagsar))
  expect_match2(scode, 
    "target += normal_lagsar_lpdf(Y | mu, sigma, lagsar, W)"
  )
  expect_match2(scode, "target += normal_lpdf(lagsar | 0.5, 1)")
  
  scode <- make_stancode(CRIME ~ INC + HOVAL, data = COL.OLD, 
                         family = student(), autocor = cor_lagsar(COL.nb))
  expect_match2(scode, 
    "target += student_t_lagsar_lpdf(Y | nu, mu, sigma, lagsar, W)"
  )
  
  scode <- make_stancode(CRIME ~ INC + HOVAL, data = COL.OLD, 
                         autocor = cor_errorsar(COL.nb))
  expect_match2(scode, 
    "target += normal_errorsar_lpdf(Y | mu, sigma, errorsar, W)"
  )
  
  scode <- make_stancode(CRIME ~ INC + HOVAL, data = COL.OLD, 
                         family = student(), autocor = cor_errorsar(COL.nb),
                         prior = prior(beta(2, 3), errorsar))
  expect_match2(scode, 
    "target += student_t_errorsar_lpdf(Y | nu, mu, sigma, errorsar, W)"
  )
  expect_match2(scode, "target += beta_lpdf(errorsar | 2, 3)")
  
  expect_error(
    make_stancode(bf(CRIME ~ INC + HOVAL, sigma ~ INC),
                  data = COL.OLD, autocor = cor_lagsar(COL.nb)),
    "SAR models are not implemented when predicting 'sigma'" 
  )
})

test_that("Stan code for CAR models is correct", {
  dat = data.frame(y = rnorm(10), x = rnorm(10))
  edges <- cbind(1:10, 10:1)
  W <- matrix(0, nrow = 10, ncol = 10)
  for (i in seq_len(nrow(edges))) {
    W[edges[i, 1], edges[i, 2]] <- 1 
  }
  
  scode <- make_stancode(y ~ x, dat, autocor = cor_car(W))
  expect_match2(scode, "real sparse_car_lpdf(vector phi")
  expect_match2(scode, "target += sparse_car_lpdf(")
  expect_match2(scode, "mu[n] += rcar[Jloc[n]]")

  scode <- make_stancode(y ~ x, dat, autocor = cor_car(W, type = "esicar"))
  expect_match2(scode, "real sparse_icar_lpdf(vector phi")
  expect_match2(scode, "target += sparse_icar_lpdf(")
  expect_match2(scode, "mu[n] += rcar[Jloc[n]]")
  expect_match2(scode, "rcar[Nloc] = - sum(zcar)")
  
  scode <- make_stancode(y ~ x, dat, autocor = cor_icar(W))
  expect_match2(scode, "target += -0.5 * dot_self(rcar[edges1] - rcar[edges2])")
  expect_match2(scode, "target += normal_lpdf(sum(rcar) | 0, 0.001 * Nloc)")
  expect_match2(scode, "mu[n] += rcar[Jloc[n]]")
})

test_that("Stan code for skew_normal models is correct", {
  dat = data.frame(y = rnorm(10), x = rnorm(10))
  scode <- make_stancode(y ~ x, dat, skew_normal())
  expect_match2(scode, "delta = alpha / sqrt(1 + alpha^2);")
  expect_match2(scode, "omega = sigma / sqrt(1 - sqrt_2_div_pi^2 * delta^2);")
  expect_match2(scode, "mu[n] = mu[n] - omega * delta * sqrt_2_div_pi;")
  
  scode <- make_stancode(bf(y ~ x, sigma ~ x), dat, skew_normal())
  expect_match2(scode, "omega[n] = sigma[n] / sqrt(1 - sqrt_2_div_pi^2 * delta^2);")
  expect_match2(scode, "mu[n] = mu[n] - omega[n] * delta * sqrt_2_div_pi;")
  
  scode <- make_stancode(bf(y | se(x) ~ x, alpha ~ x), dat, skew_normal())
  expect_match2(scode, "delta[n] = alpha[n] / sqrt(1 + alpha[n]^2);")
  expect_match2(scode, "omega[n] = se[n] / sqrt(1 - sqrt_2_div_pi^2 * delta[n]^2);")
  expect_match2(scode, "mu[n] = mu[n] - omega[n] * delta[n] * sqrt_2_div_pi;")
  
  scode <- make_stancode(y ~ x, dat, mixture(skew_normal, nmix = 2))
  expect_match2(scode, "omega1 = sigma1 / sqrt(1 - sqrt_2_div_pi^2 * delta1^2);")
  expect_match2(scode, "mu2[n] = mu2[n] - omega2 * delta2 * sqrt_2_div_pi;")
})

test_that("Stan code for missing value terms works correctly", {
  dat = data.frame(y = rnorm(10), x = rnorm(10), g = 1:10, z = 1)
  dat$x[c(1, 3, 9)] <- NA
  
  bform <- bf(y ~ mi(x)*g) + bf(x | mi() ~ g) + set_rescor(FALSE)
  scode <- make_stancode(bform, dat)
  expect_match2(scode, "Yl_x[Jmi_x] = Ymi_x;")
  expect_match2(scode, "(bsp_y[1]) * Yl_x[n] + (bsp_y[2]) * Yl_x[n] * Csp_y_1[n];")
  expect_match2(scode, "target += normal_lpdf(Yl_x | mu_x, sigma_x);")

  bform <- bf(y ~ mi(x) + (mi(x) | g)) + bf(x | mi() ~ 1) + set_rescor(FALSE)
  scode <- make_stancode(bform, dat)
  expect_match2(scode, 
    "(bsp_y[1] + r_1_y_2[J_1[n]]) * Yl_x[n] + r_1_y_1[J_1[n]] * Z_1_y_1[n];"
  )
  
  bform <- bf(y ~ a, a ~ mi(x), nl = TRUE) + bf(x | mi() ~ 1) + set_rescor(FALSE)
  bprior <- prior(normal(0, 1), nlpar = "a", resp = "y")
  scode <- make_stancode(bform, dat, prior = bprior)
  expect_match2(scode, "nlp_y_a[n] += (bsp_y_a[1]) * Yl_x[n];")
  expect_match2(scode, "target += normal_lpdf(bsp_y_a | 0, 1);")
  
  bform <- bf(y ~ mi(x)*mo(g)) + bf(x | mi() ~ 1) + set_rescor(FALSE)
  scode <- make_stancode(bform, dat)
  expect_match2(scode, "(bsp_y[3]) * Yl_x[n] * mo(simo_y_2, Xmo_y_2[n]);")
  
  bform <- bf(x | mi() ~ y, family = "lognormal")
  scode <- make_stancode(bform, dat)
  expect_match2(scode, "vector<lower=0>[Nmi] Ymi;")
  
  bform <- bf(y ~ I(log(mi(x))) * g) + 
    bf(x | mi() + trunc(lb = 1) ~ y, family = "lognormal")
  scode <- make_stancode(bform, dat)
  expect_match2(scode, "vector<lower=1>[Nmi_x] Ymi_x;")
  expect_match2(scode, 
    "(bsp_y[1]) * (log(Yl_x[n])) + (bsp_y[2]) * (log(Yl_x[n])) * Csp_y_1[n]"              
  )
  
  bform <- bf(y ~ mi(x)*g) + 
    bf(x | mi() + cens(z) ~ y, family = "beta")
  scode <- make_stancode(bform, dat)
  expect_match2(scode, "vector<lower=0,upper=1>[Nmi_x] Ymi_x;")
  expect_match2(scode, 
    "target += beta_lpdf(Yl_x[n] | mu_x[n] * phi_x, (1 - mu_x[n]) * phi_x);"
  )
})

test_that("Stan code for overimputation works correctly", {
  dat = data.frame(y = rnorm(10), x_x = rnorm(10), g = 1:10, z = 1)
  dat$x[c(1, 3, 9)] <- NA
  bform <- bf(y ~ mi(x_x)*g) + bf(x_x | mi(g) ~ 1) + set_rescor(FALSE)
  scode <- make_stancode(bform, dat, sample_prior = "yes")
  expect_match2(scode, "target += normal_lpdf(Yl_xx | mu_xx, sigma_xx)")
  expect_match2(scode, 
    "target += normal_lpdf(Y_xx[Jme_xx] | Yl_xx[Jme_xx], noise_xx[Jme_xx])"
  )
  expect_match2(scode, "vector[N] Yl_xx;")
})

test_that("argument 'stanvars' is handled correctly", {
  bprior <- prior(normal(mean_intercept, 10), class = "Intercept")
  mean_intercept <- 5
  stanvars <- stanvar(mean_intercept)
  scode <- make_stancode(count ~ Trt, data = epilepsy, 
                         prior = bprior, stanvars = stanvars)
  expect_match2(scode, "real mean_intercept;")
  
  # define a multi_normal prior with known covariance matrix
  bprior <- prior(multi_normal(M, V), class = "b")
  stanvars <- stanvar(rep(0, 2), "M", scode = "  vector[K] M;") +
    stanvar(diag(2), "V", scode = "  matrix[K, K] V;") 
  scode <- make_stancode(count ~ Trt + log_Base4_c, epilepsy,
                         prior = bprior, stanvars = stanvars)
  expect_match2(scode, "vector[K] M;")
  expect_match2(scode, "matrix[K, K] V;")
  
  # define a hierachical prior on the regression coefficients
  bprior <- set_prior("normal(0, tau)", class = "b") +
    set_prior("target += normal_lpdf(tau | 0, 10)", check = FALSE)
  stanvars <- stanvar(scode = "real<lower=0> tau;", 
                      block = "parameters")
  scode <- make_stancode(count ~ Trt + log_Base4_c, epilepsy,
                         prior = bprior, stanvars = stanvars)
  expect_match2(scode, "real<lower=0> tau;")
  expect_match2(scode, "target += normal_lpdf(b | 0, tau);")
  
  # use the non-centered parameterization for 'b'
  bprior <- set_prior("target += normal_lpdf(zb | 0, 1)", check = FALSE) +
    set_prior("target += normal_lpdf(tau | 0, 10)", check = FALSE)
  stanvars <- stanvar(scode = "vector[Kc] zb;", block = "parameters") +
    stanvar(scode = "real<lower=0> tau;", block = "parameters") +
    stanvar(scode = "vector[Kc] b = zb * tau;", 
            block="tparameters", name = "b")
  scode <- make_stancode(count ~ Trt, epilepsy, 
                         prior = bprior, stanvars = stanvars)
  expect_match2(scode, "vector[Kc] b = zb * tau;")
  
  stanvars <- stanvar(scode = "vector[Ksp] zbsp;", block = "parameters") +
    stanvar(scode = "real<lower=0> tau;", block = "parameters") +
    stanvar(scode = "vector[Ksp] bsp = zbsp * tau;", 
            block="tparameters", name = "bsp")
  scode <- make_stancode(count ~ mo(Base), epilepsy, stanvars = stanvars)
  expect_match2(scode, "vector[Ksp] bsp = zbsp * tau;")
})

test_that("custom families are handled correctly", {
  dat <- data.frame(size = 10, y = sample(0:10, 20, TRUE), x = rnorm(20))
  
  # define a custom beta-binomial family
  log_lik_beta_binomial2 <- function(i, draws) {
    mu <- draws$dpars$mu[, i]
    tau <- draws$dpars$tau
    trials <- draws$data$trials[i]
    y <- draws$data$Y[i]
    beta_binomial2_lpmf(y, mu, tau, trials)
  }
  predict_beta_binomial2 <- function(i, draws, ...) {
    mu <- draws$dpars$mu[, i]
    tau <- draws$dpars$tau
    trials <- draws$data$trials[i]
    beta_binomial2_rng(mu, tau, trials)
  }
  fitted_beta_binomial2 <- function(draws) {
    mu <- draws$dpars$mu
    trials <- draws$data$trials
    trials <- matrix(
      trials, nrow = nrow(mu), ncol = ncol(mu), byrow = TRUE
    )
    mu * trials
  }
  beta_binomial2 <- custom_family(
    "beta_binomial2", dpars = c("mu", "tau"),
    links = c("logit", "log"), lb = c(NA, 0),
    type = "int", vars = "trials[n]",
    log_lik = log_lik_beta_binomial2,
    fitted = fitted_beta_binomial2,
    predict = predict_beta_binomial2
  )
  # define custom stan functions
  stan_funs <- "
    real beta_binomial2_lpmf(int y, real mu, real phi, int N) {
      return beta_binomial_lpmf(y | N, mu * phi, (1 - mu) * phi);
    }
    int beta_binomial2_rng(real mu, real phi, int N) {
      return beta_binomial_rng(N, mu * phi, (1 - mu) * phi);
    }
  "
  stanvars <- stanvar(as.integer(dat$size), "trials") +
    stanvar(scode = stan_funs, block = "functions")
  scode <- make_stancode(
    y ~ x, data = dat, family = beta_binomial2, 
    prior = prior(gamma(0.1, 0.1), class = "tau"),
    stanvars = stanvars
  )
  expect_match2(scode, "int trials[20];")
  expect_match2(scode, "real<lower=0> tau;")
  expect_match2(scode, "mu[n] = inv_logit(mu[n]);")
  expect_match2(scode, "target += gamma_lpdf(tau | 0.1, 0.1);")
  expect_match2(scode, "target += beta_binomial2_lpmf(Y[n] | mu[n], tau, trials[n]);")
  
  scode <- make_stancode(
    bf(y ~ x, tau ~ x), data = dat, family = beta_binomial2, 
    stanvars = stanvars
  )
  expect_match2(scode, "tau[n] = exp(tau[n]); ")
  expect_match2(scode, "target += beta_binomial2_lpmf(Y[n] | mu[n], tau[n], trials[n]);")
  
  # check custom families in mixture models
  scode <- make_stancode(
    y ~ x, data = dat, family = mixture(binomial, beta_binomial2), 
    stanvars = stanvars
  )
  expect_match2(scode, 
    "log(theta2) + beta_binomial2_lpmf(Y[n] | mu2[n], tau2, trials[n]);"
  )
})

test_that("likelihood of distributional beta models is correct", {
  # test issue #404
  dat <- data.frame(prop = rbeta(100, shape1 = 2, shape2 = 2))
  scode <- make_stancode(
    bf(prop ~ 1, phi ~ 1), data = dat, family = Beta()
  )
  expect_match2(scode, "beta_lpdf(Y[n] | mu[n] * phi[n], (1 - mu[n]) * phi[n])")
})

test_that("student-t group-level effects work without errors", {
  scode <- make_stancode(count ~ Trt + (1|gr(patient, dist = "st")), epilepsy)
  expect_match2(scode, "sqrt(df_1 * udf_1) * sd_1[1] * (z_1[1]);")
  expect_match2(scode, "target += gamma_lpdf(df_1 | 2, 0.1);")
  expect_match2(scode, "target += inv_chi_square_lpdf(udf_1 | df_1);")
  
  bprior <- prior(normal(20, 5), class = df, group = patient)
  scode <- make_stancode(
    count ~ Trt + (Trt|gr(patient, dist = "st")), 
    epilepsy, prior = bprior
  )
  expect_match2(scode, "sqrt(df_1 * udf_1) * (diag_pre_multiply(sd_1, L_1) * z_1)';")
  expect_match2(scode, "target += normal_lpdf(df_1 | 20, 5);")
})
