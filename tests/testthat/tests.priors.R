# most tests of prior related stuff can be found in tests.stancode.R
context("Tests for prior generating functions")

test_that("default_prior finds all classes for which priors can be specified", {
  expect_equal(
    sort(
      default_prior(
        count ~ zBase * Trt + (1|patient) + (1+Trt|visit),
        data = epilepsy, family = "poisson"
      )$class
    ),
    sort(c(rep("b", 4), c("cor", "cor"), "Intercept", rep("sd", 6)))
  )
  expect_equal(
    sort(
      default_prior(
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
  expect_output(print(set_prior("target += normal_lpdf(x | 0,1))", check = FALSE)),
                "target += normal_lpdf(x | 0,1))", fixed = TRUE)
})

test_that("default_prior returns correct nlpar names for random effects pars", {
  # reported in issue #47
  dat <- data.frame(y = rnorm(10), x = rnorm(10), g = rep(1:2, 5))
  bform <- bf(y ~ a - b^x, a + b ~ (1+x|g), nl = TRUE)
  gp <- default_prior(bform, data = dat)
  expect_equal(sort(unique(gp$nlpar)), c("", "a", "b"))
})

test_that("default_prior returns correct fixed effect names for GAMMs", {
  dat <- data.frame(y = rnorm(10), x = rnorm(10),
                    z = rnorm(10), g = rep(1:2, 5))
  prior <- default_prior(y ~ z + s(x) + (1|g), data = dat)
  expect_equal(prior[prior$class == "b", ]$coef,
               c("", "sx_1", "z"))
  prior <- default_prior(bf(y ~ lp, lp ~ z + s(x) + (1|g), nl = TRUE),
                     data = dat)
  expect_equal(prior[prior$class == "b", ]$coef,
               c("", "Intercept", "sx_1", "z"))
})

test_that("default_prior returns correct prior names for auxiliary parameters", {
  dat <- data.frame(y = rnorm(10), x = rnorm(10),
                    z = rnorm(10), g = rep(1:2, 5))
  bform <- bf(y ~ 1, phi ~ z + (1|g), family = Beta())
  prior <- default_prior(bform, data = dat)
  prior <- prior[prior$dpar == "phi", ]
  pdata <- data.frame(class = c("b", "b", "Intercept", rep("sd", 3)),
                      coef = c("", "z", "", "", "", "Intercept"),
                      group = c(rep("", 4), "g", "g"),
                      stringsAsFactors = FALSE)
  pdata <- pdata[with(pdata, order(class, group, coef)), ]
  expect_equivalent(prior[, c("class", "coef", "group")], pdata)
})

test_that("default_prior returns correct priors for multivariate models", {
  dat <- data.frame(y1 = rnorm(10), y2 = c(1, rep(1:3, 3)),
                    x = rnorm(10), g = rep(1:2, 5))
  bform <- bf(mvbind(y1, y2) ~ x + (x|ID1|g)) + set_rescor(TRUE)

  # check global priors
  prior <- default_prior(bform, dat, family = gaussian())
  expect_equal(prior[prior$resp == "y1" & prior$class == "b", "coef"], c("", "x"))
  expect_equal(prior[prior$class == "rescor", "prior"], "lkj(1)")

  # check family and autocor specific priors
  family <- list(gaussian, Beta())
  bform <- bf(y1 ~ x + (x|ID1|g) + ar()) + bf(y2 ~ 1)
  prior <- default_prior(bform, dat, family = family)
  expect_true(any(with(prior, class == "sigma" & resp == "y1")))
  expect_true(any(with(prior, class == "ar" & resp == "y1")))
  expect_true(any(with(prior, class == "phi" & resp == "y2")))
  expect_true(!any(with(prior, class == "ar" & resp == "y2")))
})

test_that("default_prior returns correct priors for categorical models", {
  # check global priors
  dat <- data.frame(y2 = c(1, rep(1:3, 3)), x = rnorm(10), g = rep(1:2, 5))
  prior <- default_prior(y2 ~ x + (x|ID1|g), data = dat, family = categorical())
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

test_that("external interface of validate_prior works correctly", {
  prior1 <- prior(normal(0,10), class = b) +
    prior(cauchy(0,2), class = sd)
  prior1 <- validate_prior(
    prior1, count ~ zAge + zBase * Trt + (1|patient),
    data = epilepsy, family = poisson()
  )
  expect_true(all(c("b", "Intercept", "sd") %in% prior1$class))
  expect_equal(nrow(prior1), 9)
})

test_that("overall intercept priors are adjusted for the intercept", {
  dat <- data.frame(y = rep(c(1, 3), each = 5), off = 10)
  prior1 <- default_prior(y ~ 1 + offset(off), dat)
  int_prior <- prior1$prior[prior1$class == "Intercept"]
  expect_equal(int_prior, "student_t(3, -8, 2.5)")
})

test_that("as.brmsprior works correctly", {
  dat <- data.frame(prior = "normal(0,1)", x = "test", coef = c("a", "b"))
  bprior <- as.brmsprior(dat)
  expect_equal(bprior$prior, rep("normal(0,1)", 2))
  expect_equal(bprior$class, rep("b", 2))
  expect_equal(bprior$coef, c("a", "b"))
  expect_equal(bprior$x, NULL)
  expect_equal(bprior$lb, rep(NA_character_, 2))
})

test_that("prior tags are correctly applied", {

  prior1 <- prior(normal(0, 1), class = sd, tag = "prior_tag1")
  prior2 <- prior(normal(0, 5), class = b, tag = "prior_tag2")
  prior3 <- prior(normal(0, 0.5), coef = "Trt1", tag = "prior_tag3")
  prior4 <- prior(normal(0, 10), class = "Intercept", tag = "prior_tag4")
  prior5 <- prior(lkj_corr_cholesky(3), class = "L", group = "visit", tag = "prior_tag5")

  v <- validate_prior(
    c(prior1, prior2, prior3, prior4, prior5),
    formula = count ~ zBase * Trt + (1 | patient) + (1 + Trt | visit),
    data = epilepsy, family = poisson())

  expect_equal(v[which(v$class == "sd"),]$tag[[1]], "prior_tag1")
  expect_equal(v[which(v$class == "b" & v$coef != "Trt1"),]$tag[[1]], "prior_tag2")
  expect_equal(v[which(v$class == "b" & v$coef == "Trt1"),]$tag, "prior_tag3")
  expect_equal(v[which(v$class == "Intercept"),]$tag, "prior_tag4")
  expect_equal(v[which(v$class == "L"),]$tag[[2]], "prior_tag5")
})
