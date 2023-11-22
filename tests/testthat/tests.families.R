context("Tests for family functions")

test_that("family functions returns expected results", {
  expect_equal(student(identity)$link, "identity")
  expect_equal(student()$link, "identity")
  expect_equal(bernoulli(logit)$link, "logit")
  expect_error(bernoulli("sqrt"), "bernoulli")
  expect_equal(negbinomial(sqrt)$link, "sqrt")
  expect_error(negbinomial(inverse), "inverse")
  expect_equal(geometric(identity)$link, "identity")
  expect_error(geometric("inv"), "geometric")
  expect_equal(exponential(log)$link, "log")
  expect_error(exponential("cloglog"), "exponential")
  expect_equal(weibull()$family, "weibull")
  expect_error(weibull(sqrt), "weibull")
  expect_equal(Beta("probit")$link, "probit")
  expect_error(Beta("1/mu^2"), "beta")
  expect_equal(hurdle_poisson()$link, "log")
  expect_equal(hurdle_negbinomial(log)$link, "log")
  expect_error(hurdle_negbinomial("inverse"), "hurdle_negbinomial")
  expect_equal(hurdle_gamma()$family, "hurdle_gamma")
  expect_error(hurdle_gamma(sqrt), "hurdle_gamma")
  expect_equal(zero_inflated_poisson(log)$link, "log")
  expect_error(zero_inflated_poisson(list(1)), "zero_inflated_poisson")
  expect_equal(zero_inflated_negbinomial("log")$link, "log")
  expect_error(zero_inflated_negbinomial("logit"),
               "zero_inflated_negbinomial")
  expect_equal(zero_inflated_beta(logit)$family,
               "zero_inflated_beta")
  expect_equivalent(zero_inflated_binomial()$link_zi, "logit")
  expect_error(zero_inflated_binomial(y~x), "zero_inflated_binomial")
  expect_equal(categorical()$link, "logit")
  expect_error(categorical(probit), "probit")
  expect_equal(cumulative(cauchit)$family, "cumulative")
  expect_equal(sratio(probit_approx)$link, "probit_approx")
  expect_equal(cratio("cloglog")$family, "cratio")
  expect_equal(acat(cloglog)$link, "cloglog")
  expect_equal(brmsfamily("gaussian", inverse)$link, "inverse")
  expect_equal(brmsfamily("geometric", "identity")$family, "geometric")
  expect_equal(brmsfamily("zi_poisson")$link_zi, "logit")

  expect_error(weibull(link_shape = "logit"),
               "'logit' is not a supported link for parameter 'shape'")
  expect_error(weibull(link_shape = c("log", "logit")),
               "Cannot coerce 'alink' to a single character value")

  expect_equal(beta_binomial()$link, "logit")
  expect_equal(beta_binomial('probit')$link, "probit")
  expect_equal(beta_binomial()$link_phi, "log")
  expect_error(beta_binomial('log'))
  expect_error(beta_binomial(link_phi = 'logit'))
  expect_equal(zero_inflated_beta_binomial()$link, "logit")
  expect_equal(zero_inflated_beta_binomial('probit')$link, "probit")
  expect_equal(zero_inflated_beta_binomial()$link_phi, "log")
  expect_equal(zero_inflated_beta_binomial()$link_zi, "logit")
  expect_equal(zero_inflated_beta_binomial(link_zi = "identity")$link_zi, "identity")
  expect_error(zero_inflated_beta_binomial('sqrt'))
  expect_error(zero_inflated_beta_binomial(link_phi = 'logit'))
  expect_error(zero_inflated_beta_binomial(link_zi = 'log'))
  expect_equal(hurdle_cumulative()$link, "logit")
  expect_equal(hurdle_cumulative('probit')$link, "probit")
  expect_equal(hurdle_cumulative('cauchit')$link, "cauchit")
  expect_equal(hurdle_cumulative()$link_hu, "logit")
  expect_equal(hurdle_cumulative()$link_disc, "log")
  expect_error(hurdle_cumulative(link = "log")$link)
  expect_error(hurdle_cumulative(link_hu = "probit")$link_hu)
  expect_error(hurdle_cumulative(link_disc = "logit")$link_disc)

})

test_that("print brmsfamily works correctly", {
  expect_output(print(weibull()), "Family: weibull \nLink function: log")
})

test_that("mixture returns expected results and errors", {
  mix <- mixture(gaussian, nmix = 3)
  expect_equal(brms:::family_names(mix), rep("gaussian", 3))
  mix <- mixture(gaussian, student, weibull, nmix = 3:1)
  expect_equal(
    brms:::family_names(mix),
    c(rep("gaussian", 3), rep("student", 2), "weibull")
  )
  expect_error(mixture(gaussian, "x"),
               "x is not a supported family")
  expect_error(mixture(poisson(), categorical()),
               "Some of the families are not allowed in mixture models")
  expect_error(mixture(poisson, "cumulative"),
               "Cannot mix ordinal and non-ordinal families")
  expect_error(mixture(lognormal, exgaussian, poisson()),
               "Cannot mix families with real and integer support")
  expect_error(mixture(lognormal),
               "Expecting at least 2 mixture components")
  expect_error(mixture(poisson, binomial, order = "x"),
               "Argument 'order' is invalid")
})
