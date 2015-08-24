test_that("Test that link4familys returns correct links", {
  expect_equal(link4family("gaussian"), "identity")
  expect_equal(link4family("weibull"), "log")
  expect_equal(link4family("binomial"), "logit")
  expect_equal(link4family(c("binomial", "probit")), "probit")
  expect_warning(link4family(c("poisson", "sqrt")), "poisson model with sqrt link may not be uniquely identified")
})

test_that("Test that link4familys return an error on wrong links", {
  expect_error(link4family(c("gaussian","logit")), "logit is not a valid link for family gaussian")
  expect_error(link4family(c("poisson","inverse")), "inverse is not a valid link for family poisson")
  expect_error(link4family(c("weibull","sqrt")), "sqrt is not a valid link for family weibull")
  expect_error(link4family(c("categorical","probit")), "probit is not a valid link for family categorical")
})
