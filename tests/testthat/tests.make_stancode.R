test_that("make_stancode handles horseshoe priors correctly", {
  prior <- prior_frame(prior = "normal(0, hs_local * hs_global", class = "b")
  attr(prior, "hs_df") <- 7
  temp_stancode <- make_stancode(rating ~ treat*period*carry, data = inhaler,
                                 prior = prior)
  expect_match(temp_stancode, fixed = TRUE,
               "  vector<lower=0>[K] hs_local; \n  real<lower=0> hs_global; \n")
  expect_match(temp_stancode, fixed = TRUE,
               "  hs_local ~ student_t(7, 0, 1); \n  hs_global ~ cauchy(0, 1); \n")
})

test_that("make_stancode accepts supported links", {
  expect_match(make_stancode(rating ~ treat + period + carry, 
                             data = inhaler, family = sratio("probit_approx")), 
               "Phi_approx")
  expect_match(make_stancode(rating ~ treat + period + carry, 
                             data = inhaler, family = c("cumulative", "probit")), 
               "Phi")
  expect_match(make_stancode(rating ~ treat + period + carry, 
                             data = inhaler, family = "poisson"), 
               "log")
})

test_that(paste("make_stancode returns correct strings", 
                "for customized covariances"), {
                  expect_match(make_stancode(rating ~ treat + period + carry + (1|subject), 
                                             data = inhaler, cov_ranef = list(subject = 1)), 
                               "r_1 <- sd_1 * (cov_1 * pre_1)", fixed = TRUE)
                  expect_match(make_stancode(rating ~ treat + period + carry + (1+carry|subject), 
                                             data = inhaler, cov_ranef = list(subject = 1)),
                               paste0("r_1 <- to_array(kronecker_cholesky(cov_1, L_1, sd_1) * ",
                                      "to_vector(pre_1), N_1, K_1"),
                               fixed = TRUE)
                  expect_match(make_stancode(rating ~ treat + period + carry + (1+carry||subject), 
                                             data = inhaler, cov_ranef = list(subject = 1)), 
                               paste0("  r_1_1 <- sd_1[1] * (cov_1 * pre_1[1]);  // scale REs \n",
                                      "  r_1_2 <- sd_1[2] * (cov_1 * pre_1[2]);"),
                               fixed = TRUE)
                })

test_that("make_stancode handles addition arguments correctly", {
  expect_match(make_stancode(time | cens(censored) ~ age + sex + disease, 
                             data = kidney, family = c("weibull", "log")), 
               "vector[N] cens;", fixed = TRUE)
  expect_match(make_stancode(time | trunc(0) ~ age + sex + disease,
                             data = kidney, family = "gamma"), 
               "T[lb, ];", fixed = TRUE)
  expect_match(make_stancode(time | trunc(ub = 100) ~ age + sex + disease, 
                             data = kidney, family = cauchy("log")), 
               "T[, ub];", fixed = TRUE)
  expect_match(make_stancode(count | trunc(0, 150) ~ Trt_c, 
                             data = epilepsy, family = "poisson"), 
               "T[lb, ub];", fixed = TRUE)
})

test_that("make_stancode correctly combines strings of multiple grouping factors", {
  expect_match(make_stancode(count ~ (1|patient) + (1+Trt_c|visit), 
                             data = epilepsy, family = "poisson"), 
               paste0("  real Z_1[N];  // RE design matrix \n",
                      "  // data for random effects of visit \n"), 
               fixed = TRUE)
  expect_match(make_stancode(count ~ (1|visit) + (1+Trt_c|patient), 
                             data = epilepsy, family = "poisson"), 
               paste0("  int NC_1;  // number of correlations \n",
                      "  // data for random effects of visit \n"), 
               fixed = TRUE)
})

test_that("make_stancode handles models without fixed effects correctly", {
  expect_match(make_stancode(count ~ 0 + (1|patient) + (1+Trt_c|visit), 
                             data = epilepsy, family = "poisson"), 
               "  eta <- rep_vector(0, N); \n", fixed = TRUE)
})

test_that("make_stancode returns expected code for 2PL models", {
  data <- data.frame(y = rep(0:1, each = 5), x = rnorm(10))
  stancode <- make_stancode(y ~ x, data = data, 
                            family = bernoulli(type = "2PL"))
  expect_match(stancode, paste0("eta_2PL <- head(eta, N_trait)", 
                                " .* exp(tail(eta, N_trait))"),
               fixed = TRUE)
  expect_match(stancode, "Y ~ bernoulli_logit(eta_2PL);",
               fixed = TRUE)
})

test_that("make_stancode correctly restricts FE parameters", {
  data <- data.frame(y = rep(0:1, each = 5), x = rnorm(10))
  sc <- make_stancode(y ~ x, data, prior = set_prior("", lb = 2))
  expect_match(sc, "vector<lower=2>[K] b", fixed = TRUE)
  sc <- make_stancode(y ~ x, data, prior = set_prior("p", ub = "4"))
  expect_match(sc, "vector<upper=4>[K] b", fixed = TRUE)
  prior <- set_prior("normal(0,5)", lb = "-3", ub = 5)
  sc <- make_stancode(y ~ x, data, prior = prior)
  expect_match(sc, "vector<lower=-3,upper=5>[K] b", fixed = TRUE)
})

test_that("make_stancode returns correct selfmade functions", {
  # cauchit link
  expect_match(make_stancode(rating ~ treat, data = inhaler,
                             family = cumulative("cauchit")),
               "real inv_cauchit(real y)", fixed = TRUE)
  # inverse gaussian models
  temp_stancode <- make_stancode(time | cens(censored) ~ age, data = kidney,
                                 family = inverse.gaussian)
  expect_match(temp_stancode, "real inv_gaussian_log(real y", fixed = TRUE)
  expect_match(temp_stancode, "real inv_gaussian_cdf_log(real y", fixed = TRUE)
  expect_match(temp_stancode, "real inv_gaussian_ccdf_log(real y", fixed = TRUE)
  expect_match(make_stancode(time ~ 1, data = kidney, family = inverse.gaussian),
               "real inv_gaussian_log(vector y", fixed = TRUE)
  # zero-inflated and hurdle models
  expect_match(make_stancode(count ~ Trt_c, data = epilepsy, 
                             family = "zero_inflated_poisson"),
               "real zero_inflated_poisson_log(int y", fixed = TRUE)
  expect_match(make_stancode(count ~ Trt_c, data = epilepsy, 
                             family = "zero_inflated_negbinomial"),
               "real zero_inflated_neg_binomial_2_log(int y", fixed = TRUE)
  expect_match(make_stancode(count ~ Trt_c, data = epilepsy, 
                             family = "zero_inflated_binomial"),
               "real zero_inflated_binomial_log(int y", fixed = TRUE)
  expect_match(make_stancode(count ~ Trt_c, data = epilepsy, 
                             family = "zero_inflated_beta"),
               "real zero_inflated_beta_log(real y", fixed = TRUE)
  expect_match(make_stancode(count ~ Trt_c, data = epilepsy, 
                             family = hurdle_poisson()),
               "real hurdle_poisson_log(int y", fixed = TRUE)
  expect_match(make_stancode(count ~ Trt_c, data = epilepsy, 
                             family = hurdle_negbinomial),
               "real hurdle_neg_binomial_2_log(int y", fixed = TRUE)
  expect_match(make_stancode(count ~ Trt_c, data = epilepsy, 
                             family = hurdle_gamma("log")),
               "real hurdle_gamma_log(real y", fixed = TRUE)
  # linear models with special covariance structures
  expect_match(make_stancode(rating ~ treat, data = inhaler, 
                             autocor = cor_ma(cov = TRUE)),
               "real normal_cov_log(vector y", fixed = TRUE)
  expect_match(make_stancode(time ~ age, data = kidney, family = "student", 
                             autocor = cor_ar(cov = TRUE)),
               "real student_t_cov_log(vector y", fixed = TRUE)
  # ARMA covariance matrices
  expect_match(make_stancode(rating ~ treat, data = inhaler, 
                             autocor = cor_ar(cov = TRUE)),
               "matrix cov_matrix_ar1(real ar", fixed = TRUE)
  expect_match(make_stancode(time ~ age, data = kidney, family = "student", 
                             autocor = cor_ma(cov = TRUE)),
               "matrix cov_matrix_ma1(real ma", fixed = TRUE)
  expect_match(make_stancode(time ~ age, data = kidney, family = "cauchy", 
                             autocor = cor_arma(p = 1, q = 1, cov = TRUE)),
               "matrix cov_matrix_arma1(real ar, real ma", fixed = TRUE)
  # kronecker matrices
  expect_match(make_stancode(rating ~ treat + period + carry + (1+carry|subject), 
                             data = inhaler, cov_ranef = list(subject = 1)), 
               "matrix kronecker_cholesky.*vector\\[\\] to_array")
})

test_that("make_stancode detects invalid combinations of modeling options", {
  data <- data.frame(y1 = rnorm(10), y2 = rnorm(10), 
                     wi = 1:10, ci = sample(-1:1, 10, TRUE))
  expect_error(make_stancode(y1 | cens(ci) ~ y2, data = data,
                             autocor = cor_ar(cov = TRUE)),
               "Invalid addition arguments")
  expect_error(make_stancode(cbind(y1, y2) ~ 1, data = data,
                             autocor = cor_ar(cov = TRUE)),
               "multivariate models are not yet allowed")
  expect_error(make_stancode(y1 | se(wi) ~ y2, data = data,
                             autocor = cor_ma()),
               "Please set cov = TRUE", fixed = TRUE)
  expect_error(make_stancode(y1 | trunc(lb = -50) | weights(wi) ~ y2,
                             data = data),
               "truncation is not yet possible")
})

test_that("make_stancode is silent for multivariate models", {
  data <- data.frame(y1 = rnorm(10), y2 = rnorm(10), x = 1:10)
  expect_silent(make_stancode(cbind(y1, y2) ~ x, data = data))
})

test_that("make_stancode is silent for categorical models", {
  data <- data.frame(y = sample(1:4, 10, TRUE), x = 1:10)
  expect_silent(make_stancode(y ~ x, data = data, family = categorical()))
})

test_that("make_stancode returns correct code for intercept only models", {
  expect_match(make_stancode(rating ~ 1, data = inhaler),
               "b_Intercept <- temp_Intercept;", fixed = TRUE) 
  expect_match(make_stancode(rating ~ 1, data = inhaler, family = sratio()),
               "b_Intercept <- temp_Intercept;", fixed = TRUE) 
  expect_match(make_stancode(rating ~ 1, data = inhaler, family = categorical()),
               "b_Intercept <- temp_Intercept;", fixed = TRUE) 
})
