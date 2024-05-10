test_that('set_inits produces the correct format', {
  res <- set_inits('normal(0, 1)', class = "Intercept", dpar = "mu")
  expect_s3_class(res, "data.frame")
  expect_s3_class(res, "brmsinits")
  out <- as.data.frame(res)
  expect_equal(out, data.frame(distribution = "normal(0, 1)",
                               class = "Intercept",
                               coef = "",
                               group = "",
                               dpar = "",
                               nlpar = ""))

  res2 <- set_inits('normal(0, 1)', class = "sd", dpar = "sigma")
  out <- res + res2
  out <- as.data.frame(out)
  expect_equal(out, data.frame(distribution = c("normal(0, 1)", "normal(0, 1)"),
                               class = c("Intercept", "sd"),
                               coef = c("", ""),
                               group = c("", ""),
                               dpar = c("", "sigma"),
                               nlpar = c("", "")))
})


test_that('parse_dist works', {
  d <- 'normal(0, 1)'
  res <- parse_dist(d)
  expect_equal(res, list(fun = "rnorm", args = list(0, 1)))

  d <- 'uniform(-1, 1.5)'
  res <- parse_dist(d)
  expect_equal(res, list(fun = "runif", args = list(-1, 1.5)))
})


test_that('.inits_fun works', {
  data <- epilepsy
  data$cat <- factor(sample(1:5, nrow(data), replace = TRUE))
  formula <- formula <- bf(count ~ cat + Age,
                           sigma ~ cat + Age)
  bterms <- brmsterms(formula)
  sdata <- standata(formula, data = data)

  inits <- set_inits('normal(0, 1)', class = "Intercept", dpar = "mu") +
    set_inits('uniform(-1, 1)', class = "b", dpar = "sigma")

  out <- .inits_fun(inits, bterms = bterms, data = data, sdata = sdata)
  expect_type(out, "list")
  expect_equal(names(out), c("Intercept", "b_sigma"))
  expect_length(out$Intercept, 1)
  expect_equal(class(out$Intercept), "numeric")
  expect_length(out$b_sigma, 5)
  expect_equal(class(out$b_sigma), "array")
})
