test_that("Test that stan.ranef contains the correct strings", {
  expect_match(stan.ranef(list(c("Intercept","PROD"), "site"), f = c("Intercept","Prod"))$genC, 
               "cor_site\\[1,2\\]", all = FALSE)
  expect_equal(stan.ranef(list(c("Intercept"), "site"), f = c("Intercept","Prod"))$genC, 
               "  sd_site_Intercept <- sd_site; \n", all = FALSE)
})

test_that("Test that stan.model accepts supported links", {
  expect_match(stan.model(rating ~ treat + period + carry, data = inhaler, family = "sratio", 
                        link="probit_approx"), "Phi_approx")
  expect_match(stan.model(rating ~ treat + period + carry, data = inhaler, family = "cumulative", 
                        link="probit"), "Phi")
  expect_match(stan.model(rating ~ treat + period + carry, data = inhaler, family = "poisson", 
                        link="log"), "log")
})

test_that("Test that stan.prior accepts supported prior families", {
  expect_equal(stan.prior("b_x1", prior = list(b = "uniform(0,10)")), 
               "  b ~ uniform(0,10); \n")
  expect_equal(stan.prior(c("b_x1","b_x2"), prior = list(b = "uniform(0,10)", 
               b_x1 = "normal(0,1)"), ind = 1:2), 
               c("  b[1] ~ normal(0,1); \n", "  b[2] ~ uniform(0,10); \n"))
  expect_equal(stan.prior("ar", prior = list(ar = "uniform(0,1)")),
               "  ar ~ uniform(0,1); \n")
  expect_equal(stan.prior("ma", prior = list(ma = "normal(0,5)")),
               "  ma ~ normal(0,5); \n")
})

test_that("Test that stan.prior returns the correct indices", {
  expect_equal(stan.prior("sd_Intercept"), 
               "  sd ~ cauchy(0,5); \n")
  expect_equal(stan.prior("sd_Intercept", ind = "k"), 
               "  sd ~ cauchy(0,5); \n")
  expect_equal(stan.prior("sd_Intercept", ind = "k", prior = list(sd_Intercept = "normal(0,1)")), 
               "  sd[k] ~ normal(0,1); \n")
  expect_equal(stan.prior(c("sd_x1","sd_x2"), ind = 1:2, prior = list(sd_x1 = "normal(0,1)")),
               c("  sd[1] ~ normal(0,1); \n","  sd[2] ~ cauchy(0,5); \n"))                                                       
})

test_that("Test that stan.model returns correct strings (or errors) for autocorrelation models", {
  expect_match(stan.model(count~Trt_c, data=epilepsy, family = "poisson", link = "log",
                          autocor = cor.arma(~visit|patient, p=1)),
               "eta <- X\\*b \\+ Yar\\*ar")
  expect_match(stan.model(rating ~ treat + period + carry + (1|subject), data = inhaler,
                          autocor = cor.arma(~visit|patient, p=1, q=2)),
               "eta\\[n\\] <- eta\\[n\\] \\+ Ema\\[n\\]\\*ma")
  expect_error(stan.model(count~Trt_c, data=epilepsy, family = "poisson", link = "log",
                          autocor = cor.arma(~visit|patient, p=1, q=1)),
               paste0("moving-average models for family poisson require a random effect with the same number \n",
                      "of levels as observations in the data"))
  
})

test_that("Test that stan.model returns correct strings for customized covariances", {
  expect_match(stan.model(rating ~ treat + period + carry + (1|subject), data = inhaler,
                          cov.ranef = "subject"),
             "r_subject <- b\\[1\\] \\+ sd_subject \\* \\(CF_cov_subject\\*pre_subject\\)")
})