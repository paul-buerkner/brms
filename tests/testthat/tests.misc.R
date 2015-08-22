test_that("Test that brm.links returns correct links", {
  expect_equal(brm.link("gaussian"), "identity")
  expect_equal(brm.link("weibull"), "log")
  expect_equal(brm.link("binomial"), "logit")
  expect_equal(brm.link(c("binomial", "probit")), "probit")
  expect_warning(brm.link(c("poisson", "sqrt")), "poisson model with sqrt link may not be uniquely identified")
})

test_that("Test that brm.links return an error on wrong links", {
  expect_error(brm.link(c("gaussian","logit")), "logit is not a valid link for family gaussian")
  expect_error(brm.link(c("poisson","inverse")), "inverse is not a valid link for family poisson")
  expect_error(brm.link(c("weibull","sqrt")), "sqrt is not a valid link for family weibull")
  expect_error(brm.link(c("categorical","probit")), "probit is not a valid link for family categorical")
})
