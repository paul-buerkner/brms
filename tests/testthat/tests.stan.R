test_that("Test that stan.prior accepts supported prior families", {
  expect_equal(stan.prior("b_x1", prior = list(b = "uniform(0,10)")), 
               "  b ~ uniform(0,10); \n")
  expect_equal(stan.prior(c("b_x1","b_x2"), prior = list(b = "uniform(0,10)", 
               b_x1 = "normal(0,1)"), ind = 1:2), 
               "  b[1] ~ normal(0,1); \n  b[2] ~ uniform(0,10); \n")
  expect_equal(stan.prior("ar", prior = list(ar = "uniform(0,1)")),
               "  ar ~ uniform(0,1); \n")
  expect_equal(stan.prior("ma", prior = list(ma = "normal(0,5)")),
               "  ma ~ normal(0,5); \n")
  expect_equal(stan.prior("rescor", prior = list(rescor = "lkj_corr_cholesky(2)")),
               "  rescor ~ lkj_corr_cholesky(2); \n")
})

test_that("Test that stan.prior returns the correct indices", {
  expect_equal(stan.prior("sd_Intercept"), 
               "  sd ~ cauchy(0,5); \n")
  expect_equal(stan.prior("sd_Intercept", ind = "k"), 
               "  sd ~ cauchy(0,5); \n")
  expect_equal(stan.prior("sd_Intercept", ind = "k", prior = list(sd_Intercept = "normal(0,1)")), 
               "  sd[k] ~ normal(0,1); \n")
  expect_equal(stan.prior(c("sd_x1","sd_x2"), ind = 1:2, prior = list(sd_x1 = "normal(0,1)")),
               "  sd[1] ~ normal(0,1); \n  sd[2] ~ cauchy(0,5); \n") 
  expect_equal(stan.prior("sigma_y", prior = list(sigma_y = "cauchy(0,1)"), ind = 1),
               "  sigma[1] ~ cauchy(0,1); \n")
})

test_that("Test that stan.prior can remove default priors", {
  expect_equal(stan.prior("sigma_y", prior = list(sigma = "")), "")
  expect_equal(stan.prior("sd_y_Intercept", prior = list(sd = "")), "")
  expect_equal(stan.prior("sd_y_Intercept", prior = list(sd_y = ""), add.type = "y"), "")
  expect_equal(stan.prior("shape", prior = list(shape = "")), "")
})

test_that("Test that stan.eta returns correct strings for autocorrelation models", {
  expect_match(stan.eta(family = "poisson", link = "log", f = c("Trt_c"),
                        autocor = cor.arma(~visit|patient, p=1))$transC1,
               "eta <- X\\*b \\+ Yar\\*ar")
  expect_match(stan.eta(family = "poisson", link = "log", f = c("Trt_c"),
                        autocor = cor.arma(~visit|patient, q=1))$transC2,
               "eta\\[n\\] <- eta\\[n\\] \\+ Ema\\[n\\]\\*ma")
})

test_that("Test_that stan.ma returns correct strings (or errors) for moving average models", {
  expect_equal(stan.ma(family = "gaussian", link = "log", group = list("g1", "g2"), 
                       levels = c(120,60), N = 240, autocor = cor.arma()), list())
  expect_match(stan.ma(family = "poisson", link = "log", group = list("g1", "g2"), 
                 levels = c(240,60), N = 240, autocor = cor.arma(~visit|patient, q=1))$transC2,
               "Ema[n+1,i] <- r_g1[n+1-i]", fixed = TRUE)
  expect_match(stan.ma(family = "gaussian", link = "log", group = list("g1", "g2"), 
                       levels = c(120,60), N = 240, autocor = cor.arma(~visit|patient, q=1))$transC2,
               "Ema[n+1,i] <- e[n+1-i]", fixed = TRUE)
  expect_match(stan.ma(family = "multigaussian", link = "log", group = "g1", 
                       levels = 60, N = 240, autocor = cor.arma(~visit|patient, q=1))$transC2,
               "e[n] <- log(Y[m,k]) - eta[n]", fixed = TRUE)
  expect_error(stan.ma(family = "poisson", link = "log", group = list("g1", "g2"), 
                       levels = c(120,60), N = 240, autocor = cor.arma(~visit|patient, p=1, q=1)),
               paste0("moving-average models for family poisson require a random effect with the same number \n",
                      "of levels as observations in the data"))
})

test_that("Test that stan.genquant returns correct strings", {
  expect_equal(stan.genquant(family = "multigaussian", link = "identity"), list())
  expect_match(stan.genquant(family = "multigaussian", link = "identity", predict = TRUE)$genD, 
               "Y_pred[N_trait]", fixed = TRUE)
  expect_match(stan.genquant(family = "poisson", link = "log", logllh = TRUE)$genD, 
               "vector[N] log_llh", fixed = TRUE)
  expect_match(stan.genquant(family = "poisson", link = "log", logllh = TRUE)$genC, 
               "poisson_log(Y[n],exp(eta[n]))", fixed = TRUE)
  expect_match(stan.genquant(family = "weibull", link = "log", logllh = TRUE)$genC, 
               "weibull_log(Y[n],shape,eta[n])", fixed = TRUE)
})

test_that("Test that stan.model accepts supported links", {
  expect_match(stan.model(rating ~ treat + period + carry, data = inhaler, family = "sratio", 
                          link="probit_approx"), "Phi_approx")
  expect_match(stan.model(rating ~ treat + period + carry, data = inhaler, family = "cumulative", 
                          link="probit"), "Phi")
  expect_match(stan.model(rating ~ treat + period + carry, data = inhaler, family = "poisson", 
                          link="log"), "log")
})

test_that("Test that stan.model returns correct strings for customized covariances", {
  expect_match(stan.model(rating ~ treat + period + carry + (1|subject), data = inhaler,
                          cov.ranef = "subject"), fixed = TRUE,
             "r_subject <- sd_subject * (CFcov_subject*pre_subject)")
})

test_that("Test that stan.model handles addition arguments correctly", {
  expect_match(stan.model(time | cens(censored) ~ age + sex + disease, data = kidney,
                          family = "weibull", link = "log"), "vector[N] cens;", fixed = TRUE)
})

test_that("Test that stan.ord returns correct strings", {
  expect_match(stan.ord(family = "sratio", link = "logit")$par, "")
  
})

test_that("Test that stan.llh uses simplifications when possible", {
  expect_equal(stan.llh(family = "bernoulli", link = "logit"), "  Y ~ bernoulli_logit(eta); \n")
  expect_equal(stan.llh(family = "gaussian", link = "log"), "  Y ~ lognormal(eta,sigma); \n")
  expect_match(stan.llh(family = "gaussian", link = "log", weights = TRUE), 
               "lognormal_log(Y[n],eta[n],sigma); \n", fixed = TRUE)
  expect_equal(stan.llh(family = "poisson", link = "log"), "  Y ~ poisson_log(eta); \n")
  expect_match(stan.llh(family = "cumulative", link = "logit"), fixed = TRUE,
               "  Y[n] ~ ordered_logistic(eta[n],b_Intercept); \n")
})

test_that("Test that stan.llh returns correct llhs under weights and censoring", {
  expect_equal(stan.llh(family = "cauchy", link = "inverse", weights = TRUE),
               "  lp_pre[n] <- cauchy_log(Y[n],eta[n],sigma); \n")
  expect_equal(stan.llh(family = "poisson", link = "log", weights = TRUE),
               "  lp_pre[n] <- poisson_log_log(Y[n],eta[n]); \n")
  expect_match(stan.llh(family = "poisson", link = "log", cens = TRUE),
               "Y[n] ~ poisson(exp(eta[n])); \n", fixed = TRUE)
  expect_equal(stan.llh(family = "binomial", link = "logit", add = TRUE, weights = TRUE),
               "  lp_pre[n] <- binomial_logit_log(Y[n],max_obs[n],eta[n]); \n")
  expect_match(stan.llh(family = "weibull", link = "inverse", cens = TRUE), fixed = TRUE,
               "increment_log_prob(weibull_ccdf_log(Y[n],shape,eta[n])); \n")
  expect_match(stan.llh(family = "weibull", link = "inverse", cens = TRUE, weights = TRUE), fixed = TRUE,
               "increment_log_prob(weights[n] * weibull_ccdf_log(Y[n],shape,eta[n])); \n")
})

test_that("Test that stan.rngprior returns correct sampling statements for priors", {
  expect_equal(stan.rngprior(prior = "nu ~ uniform(0,100); \n"), list())
  expect_equal(stan.rngprior(prior = "nu ~ uniform(0,100); \n", sample.prior = TRUE, family = "student"),
               list(par = "  real<lower=0> prior_nu; \n", model = "  prior_nu ~ uniform(0,100); \n"))
  expect_equal(stan.rngprior(prior = "b ~ normal(0,5); \n", sample.prior = TRUE),
               list(genD = "  real prior_b; \n", genC = "  prior_b <- normal_rng(0,5); \n"))
})