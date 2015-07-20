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

test_that("Test that rename return an error on duplicated names", {
  expect_error(rename(c(letters[1:4],"a()","a[")), fixed = TRUE,
    "Internal renaming of variables led to duplicated names. \nOccured for variables: a, a(), a[")
  expect_error(rename(c("aDb","a/b","b")), fixed = TRUE,
    "Internal renaming of variables led to duplicated names. \nOccured for variables: aDb, a/b")
  expect_error(rename(c("log(a,b)","logab","bac","ba")), fixed = TRUE,
    "Internal renaming of variables led to duplicated names. \nOccured for variables: log(a,b), logab")
})
