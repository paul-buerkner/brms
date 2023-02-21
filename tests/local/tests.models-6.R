source("setup_tests_local.R")
library(mixcure)
data(leukaemia)

test_that("The generation of Stan code for mixture cure models work correctly", {
  brms:::expect_match2(make_stancode(
    formula = bf(
      time | cens(1 - cens) ~ transplant,
      inc ~ transplant
    ),
    data = leukaemia, family = mixcure_lognormal(link_inc = "logit")
  ),
  "target += mixcure_lognormal_logit_lccdf(Y[n] | mu[n], sigma, inc[n])")
  
  brms:::expect_match2(make_stancode(
    formula = bf(
      time | cens(1 - cens) ~ transplant,
      inc ~ transplant
    ),
    data = leukaemia, family = mixcure_lognormal(link_inc = "cloglog")
  ),
  "inc = inv_cloglog(inc)")
  
  brms:::expect_match2(make_stancode(
    formula = bf(
      time | cens(1 - cens) ~ transplant,
      inc ~ transplant
    ),
    data = leukaemia, family = mixcure_lognormal(link_inc = "probit")
  ),
  "inc = Phi(inc)")
  
  brms:::expect_match2(make_stancode(
    formula = bf(
      time | cens(1 - cens) ~ transplant,
      inc ~ transplant
    ),
    data = leukaemia, family = mixcure_lognormal(link_inc = "probit_approx")
  ),
  "inc = Phi_approx(inc)")
  
  brms:::expect_match2(make_stancode(
    formula = bf(
      time | cens(1 - cens) ~ transplant,
      inc ~ transplant
    ),
    data = leukaemia, family = mixcure_weibull(link_inc = "cloglog")
  ),
  "target += mixcure_weibull_lccdf(Y[n] | mu[n] / tgamma(1 + 1 / shape), shape, inc[n])")  
})

test_that("Mixture cure models work correctly", {
  fit1 <- brm(
    formula = bf(
      time | cens(1 - cens) ~ transplant,
      inc ~ transplant
    ),
    data = leukaemia, family = mixcure_lognormal(link_inc = "cloglog"),
    prior = c(set_prior("normal(0, 5)", class = "b"),
              set_prior("normal(0, 5)", class = "b", dpar = "inc")
             ),
    save_pars = save_pars(all = TRUE),
    warmup = 500, iter = 1000, refresh = 500, 
    chains = 8, cores = 8,
    control = list(adapt_delta = 0.85),
  )
  print(prior_summary(fit1))
  print(fit1, digits = 3)
  expect_range(posterior_summary(fit1)["b_Intercept", "Estimate"], 4.9, 5.1)
  expect_range(posterior_summary(fit1)["b_inc_Intercept", "Estimate"], 0.25, 0.35)
    
  fit2 <- brm(
    formula = bf(
      time | cens(1 - cens) ~ transplant,
      inc ~ transplant
    ),
    data = leukaemia, family = mixcure_weibull(link_inc = "cloglog"),
    prior = c(set_prior("normal(0, 5)", class = "b"),
              set_prior("normal(0, 5)", class = "b", dpar = "inc")
             ),
    save_pars = save_pars(all = TRUE),
    warmup = 500, iter = 1000, refresh = 500, 
    chains = 8, cores = 8,
    control = list(adapt_delta = 0.85),
  )
  print(prior_summary(fit2))
  print(fit2, digits = 3)
  expect_range(posterior_summary(fit2)["b_Intercept", "Estimate"], 5.3, 5.9)
  expect_range(posterior_summary(fit2)["b_inc_Intercept", "Estimate"], 0.25, 0.35)
  
  fit3 <- brm(
    formula = bf(
      time | cens(1 - cens) ~ transplant,
      inc ~ transplant
    ),
    data = leukaemia, family = mixcure_weibull(link_inc = "logit"),
    prior = c(set_prior("normal(0, 5)", class = "b"),
              set_prior("normal(0, 5)", class = "b", dpar = "inc")
             ),
    save_pars = save_pars(all = TRUE),
    warmup = 500, iter = 1000, refresh = 500, 
    chains = 8, cores = 8,
    control = list(adapt_delta = 0.85),
  )
  print(prior_summary(fit3))
  print(fit3, digits = 3)
  expect_range(posterior_summary(fit3)["b_Intercept", "Estimate"], 4.8, 5.4)
  expect_range(posterior_summary(fit3)["b_inc_Intercept", "Estimate"], 0.9, 1.1)
    
  fit1 <- SW(add_criterion(fit1, criterion = "loo", moment_match = TRUE))
  fit2 <- SW(add_criterion(fit2, criterion = "loo", moment_match = TRUE))
  fit3 <- SW(add_criterion(fit2, criterion = "loo", moment_match = TRUE))
  print(loo(fit1), digits = 3)
  print(loo(fit2), digits = 3)
  print(loo(fit3), digits = 3)
  print(loo_compare(fit1, fit2, fit3), digits = 3)
})
