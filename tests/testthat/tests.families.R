test_that("family functions returns expected results", {
  expect_equal(student(identity)$link, "identity")
  expect_equal(student()$link, "identity")
  expect_error(student("logit"), "student")
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
  expect_error(Beta(log), "beta")
  expect_equal(hurdle_poisson()$link, "log")
  expect_error(hurdle_poisson(identity), "hurdle_poisson")
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
  zi_binom <- list(family = "zero_inflated_binomial", link = "logit",
                   link_zi = "logit")
  class(zi_binom) <- "family"
  expect_equivalent(zero_inflated_binomial(), zi_binom)
  expect_error(zero_inflated_binomial(y~x), "zero_inflated_binomial")
  expect_equal(categorical()$link, "logit")
  expect_error(categorical(probit), "probit")
  expect_equal(cumulative(cauchit)$family, "cumulative")
  expect_equal(sratio(probit_approx)$link, "probit_approx")
  expect_equal(cratio("cloglog")$family, "cratio")
  expect_equal(acat(cloglog)$link, "cloglog")
  expect_equivalent(brmsfamily("gaussian", inverse),
                    list(family = "gaussian", link = "inverse",
                         link_sigma = "log"))
  expect_equivalent(brmsfamily("geometric", "identity"),
                    list(family = "geometric", link = "identity"))
  expect_equivalent(brmsfamily("zi_poisson"),
                    list(family = "zero_inflated_poisson", link = "log",
                         link_zi = "logit"))
  
  expect_error(weibull(link_shape = "logit"), 
               "Link 'logit' is invalid for parameter 'shape'")
  expect_error(weibull(link_shape = c("log", "logit")),
               "Link functions must be of length 1")
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
  expect_error(mixture(gaussian, categorical()), 
               "Families 'categorical' are currently not allowed in mixture models")
  expect_error(mixture(poisson, "cumulative"), 
               "Cannot mix ordinal and non-ordinal families")
  expect_error(mixture(lognormal, exgaussian, poisson()), 
               "Cannot mix families with real and integer support")
  expect_error(mixture(lognormal), 
               "Expecting at least 2 mixture components")
  expect_error(mixture(gaussian(), student(), theta = c(1, 2, 3)), 
               "The length of 'theta' should be the same as the number")
  expect_error(mixture(gaussian(), student(), theta = c(1, NA)), 
               "'theta' should contain positive values only")
  expect_error(mixture(poisson, binomial, order = "x"),
               "Argument 'order' must be either TRUE or FALSE")
})
