test_that("stan_prior accepts supported prior classes", {
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
  
  prior <- prior_frame("normal(0, 1)", class = "bp")
  expect_equal(stan_prior(class = "bp", coef = c("x1", "x2"), prior = prior),
               "  to_vector(bp) ~ normal(0, 1); \n")
})

test_that("stan_prior returns the correct indices", {
  prior <- prior_frame(prior = c("cauchy(0,5)", "normal(0,1)", "normal(0,1)"), 
                       class = c("sd", "sd", "bp"), coef = c("", "x2", "z")) 
  expect_equal(stan_prior(class = "sd", coef = "Intercept", prior = prior), 
               "  sd ~ cauchy(0,5); \n")
  expect_equal(stan_prior(class = "sd", coef = c("x1", "x2"), prior = prior), 
               "  sd[1] ~ cauchy(0,5); \n  sd[2] ~ normal(0,1); \n")
  expect_equal(stan_prior("bp", coef = "z", prior = prior),
               "  bp[1] ~ normal(0,1); \n")
})

test_that("stan_prior can remove default priors", {
  prior <- prior_frame(prior = "", class = c("sigma", "sd", "shape"), 
                       group = c("", "g", ""))
  expect_equal(stan_prior("sigma", prior = prior), "")
  expect_equal(stan_prior("sd", group = "g", prior = prior), "")
  expect_equal(stan_prior("shape", prior = prior), "")
})

test_that("stan_prior passes increment_log_prob statements without changes", {
  prior <- prior_frame(prior = c("increment_log_prob(a)", "increment_log_prob(b)"), 
                       class = rep("", 2))
  expect_equal(stan_prior("", prior = prior),
               "  increment_log_prob(a); \n  increment_log_prob(b); \n")
})

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

test_that("stan_eta returns correct strings for autocorrelation models", {
  expect_match(stan_eta(family = student(log), f = c("Trt_c"),
                        autocor = cor_arma(~visit|patient, p = 2))$transC3,
               "eta[n] <- exp(eta[n] + head(E[n], Kar) * ar)", fixed = TRUE)
  expect_match(stan_eta(family = gaussian(log), f = c("Trt_c"),
                        autocor = cor_arma(~visit|patient, q = 1))$transC2,
               "eta[n] <- eta[n] + head(E[n], Kma) * ma", fixed = TRUE)
  expect_match(stan_eta(family = poisson(), f = c("Trt_c"),
                        autocor = cor_arma(~visit|patient, r = 3))$transC1,
               "eta <- X * b + temp_Intercept + Yarr * arr", fixed = TRUE)
})

test_that("Test_that stan_arma returns correct strings (or errors)", {
  expect_equal(stan_arma(family = gaussian(log), 
                         autocor = cor.arma()), list())
  prior <- c(set_prior("normal(0,2)", class = "ar"),
             set_prior("cauchy(0,1)", class = "ma"))
  
  temp_arma <- stan_arma(family = gaussian(log), prior = prior,
                         autocor = cor_arma(~visit|patient, q = 1))
  expect_match(temp_arma$transC2, "E[n + 1, i] <- e[n + 1 - i]", 
               fixed = TRUE)
  expect_match(temp_arma$prior, "ma ~ cauchy(0,1)", fixed = TRUE)
  
  temp_arma <- stan_arma(family = gaussian(log), is_multi = TRUE, 
                         autocor = cor_arma(~visit|patient, p = 1),
                         prior = prior)
  expect_match(temp_arma$transC2, "e[n] <- log(Y[m, k]) - eta[n]", 
               fixed = TRUE)
  expect_match(temp_arma$prior, "ar ~ normal(0,2)", fixed = TRUE)
  
  temp_arma <- stan_arma(family = gaussian(log), prior = prior,
                         autocor = cor_arr(~visit|patient))
  expect_match(temp_arma$data, fixed = TRUE,
               "int<lower=1> Karr; \n  matrix[N, Karr] Yarr;")
  expect_match(temp_arma$par, "vector[Karr] arr;", fixed = TRUE)
  
  expect_error(stan_arma(family = poisson(),
                         autocor = cor.arma(~visit|patient, p = 1, q = 1)),
               "ARMA effects for family poisson are not yet implemented")
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

test_that("stan_ordinal returns correct strings", {
  expect_match(stan_ordinal(family = sratio())$par, "")
  out <- stan_ordinal(family = acat(), threshold = "equidistant")
  expect_match(out$par, "real delta;")
  expect_match(out$transC1, fixed = TRUE, 
               "temp_Intercept[k] <- temp_Intercept1 + (k - 1.0) * delta;")
  expect_match(stan_ordinal(family = cumulative("cloglog"))$fun, 
               "cumulative_log.*inv_cloglog")
  expect_match(stan_ordinal(family = sratio("logit"))$fun, 
               "sratio_log.*inv_logit")
  expect_match(stan_ordinal(family = cratio("cauchit"))$fun, 
               "cratio_log.*inv_cauchit")
  expect_match(stan_ordinal(family = acat("logit"))$fun, 
               "p[k + 1] <- exp(p[k + 1])", fixed = TRUE)
  expect_match(stan_ordinal(family = acat("probit_approx"))$fun, 
               "acat_log.*Phi_approx")
  
})

test_that("stan_llh uses simplifications when possible", {
  expect_equal(stan_llh(family = bernoulli("logit")), 
               "  Y ~ bernoulli_logit(eta); \n")
  expect_equal(stan_llh(family = gaussian("log")), 
               "  Y ~ lognormal(eta, sigma); \n")
  expect_match(stan_llh(family = gaussian("log"), weights = TRUE), 
               "lognormal_log(Y[n], (eta[n]), sigma); \n", fixed = TRUE)
  expect_equal(stan_llh(family = poisson()), 
               "  Y ~ poisson_log(eta); \n")
  expect_match(stan_llh(family = cumulative("logit")), fixed = TRUE,
               "  Y[n] ~ ordered_logistic(eta[n], temp_Intercept); \n")
})

test_that("stan_llh returns correct llhs under weights and censoring", {
  expect_equal(stan_llh(family = cauchy("inverse"), weights = TRUE),
               "  lp_pre[n] <- cauchy_log(Y[n], inv(eta[n]), sigma); \n")
  expect_equal(stan_llh(family = poisson(), weights = TRUE),
               "  lp_pre[n] <- poisson_log_log(Y[n], eta[n]); \n")
  expect_match(stan_llh(family = poisson(), cens = TRUE),
               "Y[n] ~ poisson(exp(eta[n])); \n", fixed = TRUE)
  expect_equal(stan_llh(family = binomial(logit), trials = TRUE, weights = TRUE),
               "  lp_pre[n] <- binomial_logit_log(Y[n], trials[n], eta[n]); \n")
  expect_match(stan_llh(family = weibull("log"), cens = TRUE), fixed = TRUE,
               paste("increment_log_prob(weibull_ccdf_log(Y[n],", 
                     "shape, exp(eta[n] / shape))); \n"))
  expect_match(stan_llh(family = weibull("inverse"), cens = TRUE, weights = TRUE),
               paste("increment_log_prob(weights[n] * weibull_ccdf_log(Y[n],", 
                     "shape, inv(eta[n] / shape))); \n"),
               fixed = TRUE)
})

test_that("stan_llh returns correct llhs under truncation", {
  expect_equal(stan_llh(family = cauchy(inverse), trunc = .trunc(0)),
               "  Y[n] ~ cauchy(inv(eta[n]), sigma) T[lb, ]; \n")
  expect_equal(stan_llh(family = poisson(), trunc = .trunc(ub = 100)),
               "  Y[n] ~ poisson(exp(eta[n])) T[, ub]; \n")
  expect_equal(stan_llh(family = gaussian(), 
                        se = TRUE, trunc = .trunc(0, 100)),
               "  Y[n] ~ normal((eta[n]), se[n]) T[lb, ub]; \n")
  expect_equal(stan_llh(family = binomial(), trials = TRUE, 
                        trunc = .trunc(0, 100)),
               "  Y[n] ~ binomial(trials[n], inv_logit(eta[n])) T[lb, ub]; \n")
})

test_that("stan_llh returns correct llhs for zero-inflated an hurdle models", {
  expect_equal(stan_llh(family = zero_inflated_poisson()),
               "  Y[n] ~ zero_inflated_poisson(eta[n], eta[n + N_trait]); \n")
  expect_equal(stan_llh(family = hurdle_negbinomial()),
               "  Y[n] ~ hurdle_neg_binomial_2(eta[n], eta[n + N_trait], shape); \n")
  expect_equal(stan_llh(family = hurdle_gamma()),
               "  Y[n] ~ hurdle_gamma(shape, eta[n], eta[n + N_trait]); \n")
})

test_that("stan_llh returns correct llhs for multivariate models", {
  expect_equal(stan_llh(family = gaussian(), is_multi = TRUE),
               "  Y ~ multi_normal_cholesky(Eta, LSigma); \n")
  expect_equal(stan_llh(family = student(), is_multi = TRUE),
               "  Y ~ multi_student_t(nu, Eta, Sigma); \n")
  expect_equal(stan_llh(family = cauchy(),
                        is_multi = TRUE, weights = TRUE),
               "  lp_pre[n] <- multi_student_t_log(Y[n], 1.0, Eta[n], Sigma); \n")
})

test_that("stan_rngprior returns correct sampling statements for priors", {
  c1 <- "  // parameters to store prior samples \n"
  c2 <- "  // additionally draw samples from priors \n"
  expect_equal(stan_rngprior(TRUE, prior = "nu ~ gamma(2,0.1); \n"),
               list(par = paste0(c1,"  real<lower=1> prior_nu; \n"), 
                    model = paste0(c2,"  prior_nu ~ gamma(2,0.1); \n")))
  expect_equal(stan_rngprior(TRUE, prior = "delta ~ normal(0,1); \n", 
                             family = cumulative()),
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
  expect_equal(stan_rngprior(TRUE, prior = paste0("sd_1[1] ~ normal(0,5); \n",
                                                  "  sd_1[2] ~ cauchy(0,2); \n")),
               list(par = paste0(c1, paste0("  real<lower=0> prior_sd_1_1; \n",
                                            "  real<lower=0> prior_sd_1_2; \n")), 
                    model = paste0(c2, paste0("  prior_sd_1_1 ~ normal(0,5); \n",
                                              "  prior_sd_1_2 ~ cauchy(0,2); \n"))))
  prior_code <- "b ~ normal(0, hs_local * hs_global); \n"
  expect_match(stan_rngprior(TRUE, prior = prior_code, hs_df = 3)$genC,
               "prior_b <- normal_rng(0, prior_hs_local * prior_hs_global);",
               fixed = TRUE)
  
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

test_that("stan_multi returns correct Stan code (or errors)", {
  expect_equal(stan_multi(gaussian(), "y"), list())
  expect_match(stan_multi(gaussian(), c("y1", "y2"))$transC, 
               "LSigma <- diag_pre_multiply(sigma, Lrescor); \n",
               fixed = TRUE)
  expect_equal(stan_multi(student(), c("y1", "y2"))$transD, 
               "  cov_matrix[K_trait] Sigma; \n")
  expect_equal(stan_multi(hurdle_gamma(), c("y", "huy")), list())
  expect_error(stan_multi(poisson(), c("y1", "y2")),
               "invalid multivariate model")
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
