test_that("Test that stan_prior accepts supported prior classes", {
  prior <- prior_frame(prior = "uniform(0,10)", class = "b")
  expect_equal(stan_prior(class = "b", coef = "x1", prior = prior), 
               "  b ~ uniform(0,10); \n")
  prior <- prior_frame(prior = c("uniform(0,10)", "normal(0,1)"), 
                       class = "b", coef = c("", "x1"))
  expect_equal(stan_prior(class = "b", coef = c("x1","x2"), prior = prior),
               "  b[1] ~ normal(0,1); \n  b[2] ~ uniform(0,10); \n")
  expect_equal(stan_prior("ar", prior = prior_frame("uniform(0,1)", class = "ar")),
               "  ar ~ uniform(0,1); \n")
  expect_equal(stan_prior("ma", prior = prior_frame("normal(0,5)", class = "ma")),
               "  ma ~ normal(0,5); \n")
  prior <- prior_frame("lkj_corr_cholesky(2)", class = "rescor")
  expect_equal(stan_prior("rescor", prior = prior),
               "  rescor ~ lkj_corr_cholesky(2); \n")
})

test_that("Test that stan_prior returns the correct indices", {
  prior <- prior_frame(prior = c("cauchy(0,5)", "normal(0,1)", "normal(0,1)"), 
                       class = c("sd", "sd", "bp"), coef = c("", "x2", "z")) 
  expect_equal(stan_prior(class = "sd", coef = "Intercept", prior = prior), 
               "  sd ~ cauchy(0,5); \n")
  expect_equal(stan_prior(class = "sd", coef = c("x1", "x2"), prior = prior), 
               "  sd[1] ~ cauchy(0,5); \n  sd[2] ~ normal(0,1); \n")
  expect_equal(stan_prior("bp", coef = "z", prior = prior),
               "  bp[1] ~ normal(0,1); \n")
})

test_that("Test that stan_prior can remove default priors", {
  prior <- prior_frame(prior = "", class = c("sigma", "sd", "shape"), 
                       group = c("", "g", ""))
  expect_equal(stan_prior("sigma", prior = prior), "")
  expect_equal(stan_prior("sd", group = "g", prior = prior), "")
  expect_equal(stan_prior("shape", prior = prior), "")
})

test_that("Test that stan_prior passes increment_log_prob statements without changes", {
  prior <- prior_frame(prior = c("increment_log_prob(a)", "increment_log_prob(b)"), 
                       class = rep("", 2))
  expect_equal(stan_prior("", prior = prior),
               "  increment_log_prob(a); \n  increment_log_prob(b); \n")
})

test_that("Test that stan_eta returns correct strings for autocorrelation models", {
  expect_match(stan_eta(family = "student", link = "log", f = c("Trt_c"),
                        autocor = cor_arma(~visit|patient, p = 2))$transC3,
               "eta[n] <- exp(eta[n] + head(E[n], Kar) * ar)", fixed = TRUE)
  expect_match(stan_eta(family = "gaussian", link = "log", f = c("Trt_c"),
                        autocor = cor_arma(~visit|patient, q = 1))$transC2,
               "eta[n] <- eta[n] + head(E[n], Kma) * ma", fixed = TRUE)
  expect_match(stan_eta(family = "poisson", link = "log", f = c("Trt_c"),
                        autocor = cor_arma(~visit|patient, r = 3))$transC1,
               "eta <- X * b + b_Intercept + Yarr * arr", fixed = TRUE)
})

test_that("Test_that stan_arma returns correct strings (or errors)", {
  expect_equal(stan_arma(family = "gaussian", link = "log", 
                         autocor = cor.arma()), list())
  expect_match(stan_arma(family = "gaussian", link = "log", 
                         autocor = cor.arma(~visit|patient, q = 1))$transC2,
               "E[n + 1, i] <- e[n + 1 - i]", fixed = TRUE)
  expect_match(stan_arma(family = "multinormal", link = "inverse", 
                       autocor = cor.arma(~visit|patient, p = 1))$transC2,
               "e[n] <- inv(Y[m, k]) - eta[n]", fixed = TRUE)
  expect_error(stan_arma(family = "poisson", link = "log", 
                       autocor = cor.arma(~visit|patient, p = 1, q = 1)),
               "ARMA effects for family poisson are not yet implemented")
})  

test_that("Test that make_stancode accepts supported links", {
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

test_that("Test that make_stancode returns correct strings for customized covariances", {
  expect_match(make_stancode(rating ~ treat + period + carry + (1|subject), 
                             data = inhaler, cov.ranef = list(subject = 1)), 
              "r_1 <- sd_1 * (cov_1 * pre_1)", fixed = TRUE)
  expect_match(make_stancode(rating ~ treat + period + carry + (1+carry|subject), 
                             data = inhaler, cov.ranef = list(subject = 1)),
               paste0("r_1 <- to_array(kronecker_cholesky(cov_1, L_1, sd_1) * ",
                      "to_vector(pre_1), N_1, K_1"),
               fixed = TRUE)
  expect_match(make_stancode(rating ~ treat + period + carry + (1+carry||subject), 
                             data = inhaler, cov.ranef = list(subject = 1)), 
               paste0("r_1 <- to_array(to_vector(rep_matrix(sd_1, N_1)) .* ",
                      "(cov_1 * to_vector(pre_1)), N_1, K_1)"),
               fixed = TRUE)
})

test_that("Test that make_stancode handles addition arguments correctly", {
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

test_that("Test that make_stancode correctly combines strings of multiple grouping factors", {
  expect_match(make_stancode(count ~ (1|patient) + (1+Trt_c|visit), 
                             data = epilepsy, family = "poisson"), 
               "  real Z_1[N];  # RE design matrix \n  # data for random effects of visit \n", 
               fixed = TRUE)
  expect_match(make_stancode(count ~ (1|visit) + (1+Trt_c|patient), 
                             data = epilepsy, family = "poisson"), 
               "  int NC_1;  # number of correlations \n  # data for random effects of visit \n", 
               fixed = TRUE)
})

test_that("Test that stan_ordinal returns correct strings", {
  expect_match(stan_ordinal(family = "sratio", link = "logit")$par, "")
  
})

test_that("Test that stan_llh uses simplifications when possible", {
  expect_equal(stan_llh(family = "bernoulli", link = "logit"), "  Y ~ bernoulli_logit(eta); \n")
  expect_equal(stan_llh(family = "gaussian", link = "log"), "  Y ~ lognormal(eta, sigma); \n")
  expect_match(stan_llh(family = "gaussian", link = "log", weights = TRUE), 
               "lognormal_log(Y[n], (eta[n]), sigma); \n", fixed = TRUE)
  expect_equal(stan_llh(family = "poisson", link = "log"), "  Y ~ poisson_log(eta); \n")
  expect_match(stan_llh(family = "cumulative", link = "logit"), fixed = TRUE,
               "  Y[n] ~ ordered_logistic(eta[n], b_Intercept); \n")
})

test_that("Test that stan_llh returns correct llhs under weights and censoring", {
  expect_equal(stan_llh(family = "cauchy", link = "inverse", weights = TRUE),
               "  lp_pre[n] <- cauchy_log(Y[n], inv(eta[n]), sigma); \n")
  expect_equal(stan_llh(family = "poisson", link = "log", weights = TRUE),
               "  lp_pre[n] <- poisson_log_log(Y[n], eta[n]); \n")
  expect_match(stan_llh(family = "poisson", link = "log", cens = TRUE),
               "Y[n] ~ poisson(exp(eta[n])); \n", fixed = TRUE)
  expect_equal(stan_llh(family = "binomial", link = "logit", trials = TRUE, weights = TRUE),
               "  lp_pre[n] <- binomial_logit_log(Y[n], trials[n], eta[n]); \n")
  expect_match(stan_llh(family = "weibull", link = "log", cens = TRUE), fixed = TRUE,
               "increment_log_prob(weibull_ccdf_log(Y[n], shape, exp(eta[n] / shape))); \n")
  expect_match(stan_llh(family = "weibull", link = "inverse", cens = TRUE, weights = TRUE), fixed = TRUE,
               "increment_log_prob(weights[n] * weibull_ccdf_log(Y[n], shape, inv(eta[n] / shape))); \n")
})

test_that("Test that stan_llh returns correct llhs under truncation", {
  expect_equal(stan_llh(family = "cauchy", link = "inverse", trunc = .trunc(0)),
               "  Y[n] ~ cauchy(inv(eta[n]), sigma) T[lb, ]; \n")
  expect_equal(stan_llh(family = "poisson", link = "log", trunc = .trunc(ub = 100)),
               "  Y[n] ~ poisson(exp(eta[n])) T[, ub]; \n")
  expect_equal(stan_llh(family = "gaussian", link = "identity", 
                        se = TRUE, trunc = .trunc(0, 100)),
               "  Y[n] ~ normal((eta[n]), se[n]) T[lb, ub]; \n")
  expect_equal(stan_llh(family = "binomial", link = "logit", 
                        trials = TRUE, trunc = .trunc(0, 100)),
               "  Y[n] ~ binomial(trials[n], inv_logit(eta[n])) T[lb, ub]; \n")
})

test_that("Test that stan_llh returns correct llhs for zero-inflated an hurdle models", {
  expect_equal(stan_llh(family = "zero_inflated_poisson", link = "log"),
               "  Y[n] ~ zero_inflated_poisson(eta[n], eta[n + N_trait]); \n")
  expect_equal(stan_llh(family = "hurdle_negbinomial", link = "log"),
               "  Y[n] ~ hurdle_neg_binomial_2(eta[n], eta[n + N_trait], shape); \n")
  expect_equal(stan_llh(family = "hurdle_gamma", link = "log"),
               "  Y[n] ~ hurdle_gamma(shape, eta[n], eta[n + N_trait]); \n")
})

test_that("Test that stan_rngprior returns correct sampling statements for priors", {
  c1 <- "  # parameters to store prior samples \n"
  c2 <- "  # additionally draw samples from priors \n"
  expect_equal(stan_rngprior(TRUE, prior = "nu ~ uniform(0,100); \n"),
               list(par = paste0(c1,"  real<lower=0> prior_nu; \n"), 
                    model = paste0(c2,"  prior_nu ~ uniform(0,100); \n")))
  expect_equal(stan_rngprior(TRUE, prior = "delta ~ normal(0,1); \n", family = "cumulative"),
               list(par = paste0(c1,"  real<lower=0> prior_delta; \n"), 
                    model = paste0(c2,"  prior_delta ~ normal(0,1); \n")))
  expect_equal(stan_rngprior(TRUE, prior = "b ~ normal(0,5); \n"),
               list(genD = "  real prior_b; \n", 
                    genC = paste0(c2,"  prior_b <- normal_rng(0,5); \n")))
  expect_equal(stan_rngprior(TRUE, prior = "b[1] ~ normal(0,5); \n"),
               list(genD = "  real prior_b_1; \n", 
                    genC = paste0(c2,"  prior_b_1 <- normal_rng(0,5); \n")))
  expect_equal(stan_rngprior(TRUE, prior = "bp[1] ~ normal(0,5); \n"),
               list(genD = "  real prior_bp_1; \n", 
                    genC = paste0(c2,"  prior_bp_1 <- normal_rng(0,5); \n")))
  expect_equal(stan_rngprior(TRUE, prior = "sigma[2] ~ normal(0,5); \n"),
               list(par = paste0(c1,"  real<lower=0> prior_sigma_2; \n"), 
                    model = paste0(c2,"  prior_sigma_2 ~ normal(0,5); \n")))
  expect_equal(stan_rngprior(TRUE, prior = "sd_1[1] ~ normal(0,5); \n  sd_1[2] ~ cauchy(0,2); \n"),
               list(par = paste0(c1,"  real<lower=0> prior_sd_1_1; \n  real<lower=0> prior_sd_1_2; \n"), 
                    model = paste0(c2,"  prior_sd_1_1 ~ normal(0,5); \n  prior_sd_1_2 ~ cauchy(0,2); \n")))
})

test_that("Test that stan_functions returns relevant user defined functions", {
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
                             data = inhaler, cov.ranef = list(subject = 1)), 
               "matrix kronecker_cholesky.*vector\\[\\] to_array")
})