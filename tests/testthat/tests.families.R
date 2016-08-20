test_that("family functions returns expected results", {
  expect_equal(student(identity)$link, "identity")
  expect_equal(student()$link, "identity")
  expect_error(student("logit"), "student")
  expect_equal(cauchy(log)$family, "cauchy")
  expect_error(cauchy("inv"), "cauchy")
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
  expect_error(zero_inflated_poisson(list()), "zero_inflated_poisson")
  expect_equal(zero_inflated_negbinomial("log")$link, "log")
  expect_error(zero_inflated_negbinomial("logit"), 
               "zero_inflated_negbinomial")
  expect_equal(zero_inflated_beta(logit)$family, 
               "zero_inflated_beta")
  zi_binom <- list(family = "zero_inflated_binomial", link = "logit")
  class(zi_binom) <- "family"
  expect_equivalent(zero_inflated_binomial(), zi_binom)
  expect_error(zero_inflated_binomial(y~x), "zero_inflated_binomial")
  expect_equal(categorical()$link, "logit")
  expect_error(categorical(probit), "probit")
  expect_equal(cumulative(cauchit)$family, "cumulative")
  expect_equal(sratio(probit_approx)$link, "probit_approx")
  expect_equal(cratio("cloglog")$family, "cratio")
  expect_equal(acat(cloglog)$link, "cloglog")
})

test_that("check_family returns correct links", {
  expect_equal(check_family("gaussian")$link, "identity")
  expect_equal(check_family("weibull")$link, "log")
  expect_equal(check_family(binomial)$link, "logit")
  expect_equal(check_family(binomial("probit"))$link, "probit")
  expect_equal(check_family(c("acat", "cloglog"))$link, "cloglog")
})

test_that("check_family return an error on wrong links", {
  expect_error(check_family(gaussian("logit")), 
               "logit is not a supported link for family gaussian")
  expect_error(check_family(poisson("inverse")), 
               "inverse is not a supported link for family poisson")
  expect_error(check_family(c("weibull", "sqrt")), 
               "sqrt is not a supported link for family weibull")
  expect_error(check_family(c("categorical","probit")), 
               "probit is not a supported link for family categorical")
})

test_that("check_family rejects invalid families", {
  expect_error(check_family("multigaussian"),
               "multigaussian is not a supported family")
  expect_error(check_family("ordinal"),
               "ordinal is not a supported family")
})

test_that("print brmsfamily works correctly", {
  expect_output(print(weibull()), "Family: weibull \nLink function: log")
})