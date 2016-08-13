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
  
  prior <- prior_frame("normal(0, 1)", class = "b")
  expect_equal(stan_prior(class = "b", coef = c("x1", "x2"), suffix = "p", 
                          matrix = TRUE, prior = prior),
               "  to_vector(bp) ~ normal(0, 1); \n")
})

test_that("stan_prior returns the correct indices", {
  prior <- prior_frame(prior = c("cauchy(0,5)", "normal(0,1)", "normal(0,1)"), 
                       class = c("sd", "sd", "bp"), coef = c("", "x2", "z")) 
  expect_equal(stan_prior(class = "sd", coef = "Intercept", prior = prior), 
               "  sd ~ cauchy(0,5); \n")
  expect_equal(stan_prior(class = "sd", coef = c("x1", "x2"), prior = prior), 
               "  sd[1] ~ cauchy(0,5); \n  sd[2] ~ normal(0,1); \n")
  expect_equal(stan_prior("bp", coef = "z", prior = prior, matrix = TRUE),
               "  bp[1] ~ normal(0,1); \n")
})

test_that("stan_prior can remove default priors", {
  prior <- prior_frame(prior = "", class = c("sigma", "sd", "shape"), 
                       group = c("", "g", ""))
  expect_equal(stan_prior("sigma", prior = prior), "")
  expect_equal(stan_prior("sd", group = "g", suffix = "_1", prior = prior), "")
  expect_equal(stan_prior("shape", prior = prior), "")
})

test_that("stan_prior passes increment_log_prob statements without changes", {
  prior <- prior_frame(prior = c("increment_log_prob(a)", "increment_log_prob(b)"), 
                       class = rep("", 2))
  expect_equal(stan_prior("", prior = prior),
               "  increment_log_prob(a); \n  increment_log_prob(b); \n")
})

test_that("stan_eta returns correct strings for autocorrelation models", {
  ee <- extract_effects(count ~ Trt_c)
  expect_match(stan_linear(ee, data = epilepsy, family = student(log),
                           autocor = cor_arma(~visit|patient, p = 2))$modelC3,
               "eta[n] = exp(eta[n] + head(E[n], Kar) * ar)", fixed = TRUE)
  expect_match(stan_linear(ee, data = epilepsy, family = gaussian(log),
                           autocor = cor_arma(~visit|patient, q = 1))$modelC2,
               "eta[n] = exp(eta[n] + head(E[n], Kma) * ma)", fixed = TRUE)
  expect_match(stan_linear(ee, data = epilepsy, family = poisson(),
                           autocor = cor_arma(~visit|patient, r = 3))$modelC1,
               "eta = X * b + temp_Intercept + Yarr * arr", fixed = TRUE)
})

test_that("stan_autocor returns correct strings (or errors)", {
  expect_equal(stan_autocor(family = gaussian(log), 
                            autocor = cor.arma()), list())
  prior <- c(set_prior("normal(0,2)", class = "ar"),
             set_prior("cauchy(0,1)", class = "ma"))
  
  temp_arma <- stan_autocor(family = gaussian(log), prior = prior,
                            autocor = cor_arma(~visit|patient, q = 1))
  expect_match(temp_arma$modelC2, "E[n + 1, i] = e[n + 1 - i]", 
               fixed = TRUE)
  expect_match(temp_arma$prior, "ma ~ cauchy(0,1)", fixed = TRUE)
  
  effects <- list(response = c("y1", "y2"))
  temp_arma <- stan_autocor(family = gaussian(log), effects = effects, 
                            autocor = cor_arma(~visit|patient, p = 1),
                            prior = prior)
  expect_match(temp_arma$modelC2, "e[n] = log(Y[m, k]) - eta[n]", 
               fixed = TRUE)
  expect_match(temp_arma$prior, "ar ~ normal(0,2)", fixed = TRUE)
  
  temp_arma <- stan_autocor(family = gaussian(log), prior = prior,
                            autocor = cor_arr(~visit|patient))
  expect_match(temp_arma$data, fixed = TRUE,
               "int<lower=1> Karr; \n  matrix[N, Karr] Yarr;")
  expect_match(temp_arma$par, "vector[Karr] arr;", fixed = TRUE)
  
  expect_error(stan_autocor(family = poisson(),
                            autocor = cor.arma(~visit|patient, p = 1, q = 1)),
               "ARMA effects for family poisson are not yet implemented")
})  

test_that("stan_ordinal returns correct strings", {
  expect_match(stan_ordinal(family = sratio())$par, "")
  out <- stan_ordinal(family = acat(), threshold = "equidistant")
  expect_match(out$par, "real delta;")
  expect_match(out$transC1, fixed = TRUE, 
               "temp_Intercept[k] = temp_Intercept1 + (k - 1.0) * delta;")
  expect_match(stan_ordinal(family = cumulative("cloglog"))$fun, 
               "cumulative_lpmf.*inv_cloglog")
  expect_match(stan_ordinal(family = sratio("logit"))$fun, 
               "sratio_lpmf.*inv_logit")
  expect_match(stan_ordinal(family = cratio("cauchit"))$fun, 
               "cratio_lpmf.*inv_cauchit")
  expect_match(stan_ordinal(family = acat("logit"))$fun, 
               "p[k + 1] = exp(p[k + 1])", fixed = TRUE)
  expect_match(stan_ordinal(family = acat("probit_approx"))$fun, 
               "acat_lpmf.*Phi_approx")
  
})

test_that("stan_llh uses simplifications when possible", {
  expect_equal(stan_llh(bernoulli("logit")), 
               "  Y ~ bernoulli_logit(eta); \n")
  expect_match(stan_llh(lognormal(), effects = list(weights = ~x)), 
               "lognormal_lpdf(Y[n] | (eta[n]), sigma); \n", fixed = TRUE)
  expect_equal(stan_llh(poisson()), "  Y ~ poisson_log(eta); \n")
  expect_match(stan_llh(cumulative("logit")), fixed = TRUE,
               "  Y[n] ~ ordered_logistic(eta[n], temp_Intercept); \n")
})

test_that("stan_llh returns correct llhs under weights and censoring", {
  expect_match(stan_llh(cauchy("inverse"), effects = list(weights = ~x)),
               "  lp_pre[n] = cauchy_lpdf(Y[n] | inv(eta[n]), sigma); \n",
               fixed = TRUE)
  expect_match(stan_llh(poisson(), effects = list(weights = ~x)),
               "  lp_pre[n] = poisson_log_lpmf(Y[n] | (eta[n])); \n",
               fixed = TRUE)
  expect_match(stan_llh(poisson(), effects = list(cens = ~x)),
               "Y[n] ~ poisson(exp(eta[n])); \n", fixed = TRUE)
  expect_match(stan_llh(binomial(logit), list(weights = ~x, trials = ~x)),
               "  lp_pre[n] = binomial_logit_lpmf(Y[n] | trials[n], (eta[n])); \n",
               fixed = TRUE)
  expect_match(stan_llh(weibull("log"), effects = list(cens = ~x)), fixed = TRUE,
               "target += weibull_lccdf(Y[n] | shape, exp(eta[n] / shape)); \n")
  expect_match(stan_llh(weibull("inverse"), list(cens = ~x, weights = ~x)),
               paste("target += weights[n] * weibull_lccdf(Y[n] |", 
                     "shape, inv(eta[n] / shape)); \n"), fixed = TRUE)
})

test_that("stan_llh returns correct llhs under truncation", {
  
  expect_match(stan_llh(cauchy(inverse), trunc_bounds = list(lb = 0)),
               "  Y[n] ~ cauchy(inv(eta[n]), sigma) T[lb[n], ];", fixed = TRUE)
  expect_match(stan_llh(poisson(), trunc_bounds = list(ub = 100)),
               "  Y[n] ~ poisson(exp(eta[n])) T[, ub[n]]; \n", fixed = TRUE)
  expect_match(stan_llh(gaussian(), effects = list(se = ~x),
                        trunc_bounds = list(lb = 0, ub = 100)),
               "  Y[n] ~ normal((eta[n]), se[n]) T[lb[n], ub[n]];", fixed = TRUE)
  expect_match(stan_llh(binomial(), effects = list(trials = ~x),
                        trunc_bounds = list(lb = 0, ub = 100)),
               "  Y[n] ~ binomial(trials[n], inv_logit(eta[n])) T[lb[n], ub[n]];",
               fixed = TRUE)
})

test_that("stan_llh returns correct llhs for zero-inflated an hurdle models", {
  expect_match(stan_llh(zero_inflated_poisson()), fixed = TRUE,
               "  Y[n] ~ zero_inflated_poisson(eta[n], eta[n + N_trait]);")
               
  expect_match(stan_llh(hurdle_negbinomial()), fixed = TRUE,
               "  Y[n] ~ hurdle_neg_binomial_2(eta[n], eta[n + N_trait], shape);")
  expect_match(stan_llh(hurdle_gamma()), fixed = TRUE,
               "  Y[n] ~ hurdle_gamma(shape, eta[n], eta[n + N_trait]);")
})

test_that("stan_llh returns correct llhs for multivariate models", {
  expect_match(stan_llh(gaussian(), list(response = c("y1", "y2"))),
               "  Y ~ multi_normal_cholesky(Eta, LSigma); \n", fixed = TRUE)
  expect_match(stan_llh(student(), list(response = c("y1", "y2", "y3"))),
               "  Y ~ multi_student_t(nu, Eta, Sigma); \n", fixed = TRUE)
  expect_match(stan_llh(cauchy(), list(response = c("y1", "y2"), weights = ~x)),
               "  lp_pre[n] = multi_student_t_lpdf(Y[n] | 1.0, (Eta[n]), Sigma);",
               fixed = TRUE)
})

test_that("stan_rngprior returns correct sampling statements for priors", {
  c1 <- "  // parameters to store prior samples \n"
  c2 <- "  // additionally draw samples from priors \n"
  expect_equal(stan_rngprior(TRUE, prior = "nu ~ gamma(2,0.1); \n",
                             par_declars = "real<lower=1> nu; \n"),
               list(par = paste0(c1,"  real<lower=1> prior_nu; \n"), 
                    model = paste0(c2,"  prior_nu ~ gamma(2,0.1); \n")))
  
  expect_equal(stan_rngprior(TRUE, prior = "delta ~ normal(0,1); \n", 
                             par_declars = "real<lower=0> delta; \n",
                             family = cumulative()),
               list(par = paste0(c1,"  real<lower=0> prior_delta; \n"), 
                    model = paste0(c2,"  prior_delta ~ normal(0,1); \n")))
  
  expect_equal(stan_rngprior(TRUE, prior = "b ~ normal(0,5); \n",
                             par_declars = paste("vector[K] b; \n",  
                                                 "real<lower=0> sigma;")),
               list(genD = "  real prior_b; \n", 
                    genC = paste0(c2,"  prior_b = normal_rng(0,5); \n")))
  
  expect_equal(stan_rngprior(TRUE, prior = "b[1] ~ normal(0,5); \n",
                             par_declars = "vector[K] b; \n"),
               list(genD = "  real prior_b_1; \n", 
                    genC = paste0(c2,"  prior_b_1 = normal_rng(0,5); \n")))
  
  expect_equal(stan_rngprior(TRUE, prior = "b_m[1] ~ normal(0,5); \n",
                             par_declars = "vector[K] b_m; \n"),
               list(genD = "  real prior_b_m_1; \n", 
                    genC = paste0(c2,"  prior_b_m_1 = normal_rng(0,5); \n")))
  
  expect_equal(stan_rngprior(TRUE, prior = "bp[1] ~ normal(0,5); \n",
                             par_declars = paste("vector[K] b; \n")),
               list(genD = "  real prior_bp_1; \n", 
                    genC = paste0(c2,"  prior_bp_1 = normal_rng(0,5); \n")))
  
  expect_equal(stan_rngprior(TRUE, prior = "sigma[2] ~ normal(0,5); \n",
                             par_declars = paste("vector[K] b; \n",  
                                                 "real<lower=0> sigma;")),
               list(par = paste0(c1,"  real<lower=0> prior_sigma_2; \n"), 
                    model = paste0(c2,"  prior_sigma_2 ~ normal(0,5); \n")))
  
  expect_equal(stan_rngprior(TRUE, 
                 prior = "simplex_1 ~ dirichlet(con_simplex_1); \n",
                 par_declars = "vector[Jm[1]] simplex_1;"),
    list(genD = "  vector[Jm[1]] prior_simplex_1; \n", 
         genC = paste0(c2, "  prior_simplex_1 = dirichlet_rng(con_simplex_1); \n")))
  
  expect_equal(stan_rngprior(TRUE, prior = paste("sd_1[1] ~ normal(0,5); \n",
                                                 "sd_1[2] ~ cauchy(0,2); \n"),
                             par_declars = "vector<lower=0> sd_1;"),
               list(par = paste0(c1, paste0("  real<lower=0> prior_sd_1_1; \n",
                                            "  real<lower=0> prior_sd_1_2; \n")), 
                    model = paste0(c2, paste0("  prior_sd_1_1 ~ normal(0,5); \n",
                                              "  prior_sd_1_2 ~ cauchy(0,2); \n"))))
  
  prior_code <- "b ~ normal(0, hs_local * hs_global); \n"
  expect_match(stan_rngprior(TRUE, prior = prior_code, hs_df = 3)$genC,
               "prior_b = normal_rng(0, prior_hs_local * prior_hs_global);",
               fixed = TRUE)
  
})

test_that("stan_multi returns correct Stan code (or errors)", {
  expect_equal(stan_multi(gaussian(), "y"), list())
  expect_match(stan_multi(gaussian(), c("y1", "y2"))$transC, 
               "LSigma = diag_pre_multiply(sigma, Lrescor); \n",
               fixed = TRUE)
  expect_equal(stan_multi(student(), c("y1", "y2"))$transD, 
               "  cov_matrix[K_trait] Sigma; \n")
  expect_equal(stan_multi(hurdle_gamma(), c("y", "huy")), list())
  expect_error(stan_multi(poisson(), c("y1", "y2")),
               "invalid multivariate model")
})
