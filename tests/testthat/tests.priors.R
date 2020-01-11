# most tests of prior related stuff can be found in tests.make_stancode.R
context("Tests for prior generating functions")

test_that("get_prior finds all classes for which priors can be specified", {
  expect_equal(
    sort(
      get_prior(
        count ~ zBase * Trt + (1|patient) + (1+Trt|visit),
        data = epilepsy, family = "poisson"
      )$class
    ),
    sort(c(rep("b", 4), c("cor", "cor"), "Intercept", rep("sd", 6)))
  )
  expect_equal(
    sort(
      get_prior(
        rating ~ treat + period + cse(carry), data = inhaler, 
        family = sratio(threshold = "equidistant")
      )$class
    ),
    sort(c(rep("b", 4), "delta", rep("Intercept", 1)))
  )
})

test_that("set_prior allows arguments to be vectors", {
  bprior <- set_prior("normal(0, 2)", class = c("b", "sd"))
  expect_is(bprior, "brmsprior")
  expect_equal(bprior$prior, rep("normal(0, 2)", 2))
  expect_equal(bprior$class, c("b", "sd"))
})

test_that("print for class brmsprior works correctly", {
  expect_output(print(set_prior("normal(0,1)")), fixed = TRUE,
                "b ~ normal(0,1)")
  expect_output(print(set_prior("normal(0,1)", coef = "x")), 
                "b_x ~ normal(0,1)", fixed = TRUE)
  expect_output(print(set_prior("cauchy(0,1)", class = "sd", group = "x")), 
                "sd_x ~ cauchy(0,1)", fixed = TRUE)
  expect_output(print(set_prior("increment_log_prob(normal_log(0,1))")), 
                "increment_log_prob(normal_log(0,1))", fixed = TRUE)
})

test_that("get_prior returns correct nlpar names for random effects pars", {
  # reported in issue #47
  data <- data.frame(y = rnorm(10), x = rnorm(10), g = rep(1:2, 5))
  gp <- get_prior(bf(y ~ a - b^x, a + b ~ (1+x|g), nl = TRUE), 
                  data = data)
  expect_equal(sort(unique(gp$nlpar)), c("", "a", "b"))
})

test_that("get_prior returns correct fixed effect names for GAMMs", {
  dat <- data.frame(y = rnorm(10), x = rnorm(10), 
                    z = rnorm(10), g = rep(1:2, 5))
  prior <- get_prior(y ~ z + s(x) + (1|g), data = dat)
  expect_equal(prior[prior$class == "b", ]$coef, 
               c("", "sx_1", "z"))
  prior <- get_prior(bf(y ~ lp, lp ~ z + s(x) + (1|g), nl = TRUE), 
                     data = dat)
  expect_equal(prior[prior$class == "b", ]$coef, 
               c("", "Intercept", "sx_1", "z"))
})

test_that("get_prior returns correct prior names for auxiliary parameters", {
  dat <- data.frame(y = rnorm(10), x = rnorm(10), 
                    z = rnorm(10), g = rep(1:2, 5))
  prior <- get_prior(bf(y ~ 1, phi ~ z + (1|g)), data = dat, family = Beta())
  prior <- prior[prior$dpar == "phi", ]
  pdata <- data.frame(class = c("b", "b", "Intercept", rep("sd", 3)), 
                      coef = c("", "z", "", "", "", "Intercept"),
                      group = c(rep("", 4), "g", "g"),
                      stringsAsFactors = FALSE)
  pdata <- pdata[with(pdata, order(class, group, coef)), ]
  expect_equivalent(prior[, c("class", "coef", "group")], pdata)
})

test_that("get_prior returns correct priors for multivariate models", {
  dat <- data.frame(y1 = rnorm(10), y2 = c(1, rep(1:3, 3)), 
                    x = rnorm(10), g = rep(1:2, 5))
  bform <- bf(mvbind(y1, y2) ~ x + (x|ID1|g))
  
  # check global priors
  prior <- get_prior(bform, dat, family = gaussian())
  expect_equal(prior[prior$resp == "y1" & prior$class == "b", "coef"], c("", "x"))
  expect_equal(prior[prior$class == "rescor", "prior"], "lkj(1)")
  
  # check family and autocor specific priors
  family <- list(gaussian, Beta())
  autocor <- list(cor_ar(), NULL)
  prior <- get_prior(bform, dat, family = family, autocor = autocor)
  expect_true(any(with(prior, class == "sigma" & resp == "y1")))
  expect_true(any(with(prior, class == "ar" & resp == "y1")))
  expect_true(any(with(prior, class == "phi" & resp == "y2")))
  expect_true(!any(with(prior, class == "ar" & resp == "y2")))
})

test_that("get_prior returns correct priors for categorical models", {
  # check global priors
  dat <- data.frame(y2 = c(1, rep(1:3, 3)), x = rnorm(10), g = rep(1:2, 5))
  prior <- get_prior(y2 ~ x + (x|ID1|g), data = dat, family = categorical())
  expect_equal(prior[prior$dpar == "mu2" & prior$class == "b", "coef"], c("", "x"))
})

test_that("set_prior alias functions produce equivalent results", {
  expect_equal(set_prior("normal(0, 1)", class = "sd"),
               prior(normal(0, 1), class = sd))
  expect_equal(set_prior("normal(0, 1)", class = "sd", nlpar = "a"),
               prior(normal(0, 1), class = "sd", nlpar = a))
  expect_equal(set_prior("normal(0, 1)", class = "sd", nlpar = "a"),
               prior_(~normal(0, 1), class = ~sd, nlpar = quote(a)))
  expect_equal(set_prior("normal(0, 1)", class = "sd"),
               prior_string("normal(0, 1)", class = "sd"))
})