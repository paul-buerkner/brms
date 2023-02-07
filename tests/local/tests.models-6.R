source("setup_tests_local.R")

test_that("Mixture cure models work correctly", {
    data(leukaemia)
    fit1 <- brm(
        formula = bf(
            time | cens(1 - cens) ~ transplant,
            inc ~ transplant
        ),
        data = leukaemia, family = mixcure_lognormal(),
        prior = c(set_prior("normal(0, 10)", class = "b"),
                  set_prior("normal(0, 10)", class = "b", dpar = "inc")
                 ),
        save_pars = save_pars(all = TRUE),
        warmup = 1000, iter = 2000, refresh = 2000, 
        chains = 8, cores = 8,
        control = list(adapt_delta = 0.825),
    )
    print(fit1, digits = 3)
    print(prior_summary(fit1))
    expect_range(posterior_summary(fit1)["b_Intercept", "Estimate"], 4.9, 5.1)
    expect_range(posterior_summary(fit1)["b_inc_Intercept", "Estimate"], 0.95, 1.05)
    
    fit2 <- brm(
        formula = bf(
            time | cens(1 - cens) ~ transplant,
            inc ~ transplant
        ),
        data = leukaemia, family = mixcure_weibull(),
        prior = c(set_prior("normal(0, 10)", class = "b"),
                  set_prior("normal(0, 10)", class = "b", dpar = "inc")
                 ),
        save_pars = save_pars(all = TRUE),
        warmup = 1000, iter = 2000, refresh = 2000, 
        chains = 8, cores = 8,
        control = list(adapt_delta = 0.825),
    )
    print(fit2, digits = 3)
    print(prior_summary(fit2))
    expect_range(posterior_summary(fit2)["b_Intercept", "Estimate"], 4.9, 5.1)
    expect_range(posterior_summary(fit2)["b_inc_Intercept", "Estimate"], 0.95, 1.05)
})
