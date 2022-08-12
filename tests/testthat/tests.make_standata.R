context("Tests for make_standata")

test_that(paste("make_standata returns correct data names ",
                "for fixed and random effects"), {
  expect_equal(names(make_standata(rating ~ treat + period + carry
                                   + (1|subject), data = inhaler)),
               c("N", "Y",  "K", "X", "Z_1_1",
                 "J_1", "N_1", "M_1", "NC_1", "prior_only"))
  expect_equal(names(make_standata(rating ~ treat + period + carry
                                   + (1+treat|id|subject), data = inhaler,
                                   family = "categorical")),
               c("N", "Y", "ncat", "K_mu2", "X_mu2", "Z_1_mu2_1",
                 "Z_1_mu2_2", "K_mu3", "X_mu3", "Z_1_mu3_3", "Z_1_mu3_4",
                 "K_mu4", "X_mu4", "Z_1_mu4_5", "Z_1_mu4_6",
                 "J_1", "N_1", "M_1", "NC_1", "prior_only"))
  expect_equal(names(make_standata(rating ~ treat + period + carry
                                   + (1+treat|subject), data = inhaler)),
               c("N", "Y", "K", "X", "Z_1_1", "Z_1_2", "J_1", "N_1", "M_1",
                 "NC_1", "prior_only"))

  dat <- data.frame(y = 1:10, g = 1:10, h = 11:10, x = rep(0,10))
  expect_equal(names(make_standata(y ~ x + (1|g) + (1|h), dat, "poisson")),
               c("N", "Y", "K", "X", "Z_1_1", "Z_2_1",
                 "J_1", "J_2", "N_1", "M_1", "NC_1", "N_2", "M_2", "NC_2",
                 "prior_only"))
  expect_true(all(c("Z_1_1", "Z_1_2", "Z_2_1", "Z_2_2") %in%
                  names(make_standata(y ~ x + (1+x|g/h), dat))))
  expect_equal(make_standata(y ~ x + (1+x|g+h), dat),
               make_standata(y ~ x + (1+x|g) + (1+x|h), dat))
})

test_that(paste("make_standata handles variables used as fixed effects",
                "and grouping factors at the same time"), {
  data <- data.frame(y = 1:9, x = factor(rep(c("a","b","c"), 3)))
  standata <- make_standata(y ~ x + (1|x), data = data)
  expect_equal(colnames(standata$X), c("Intercept", "xb", "xc"))
  expect_equal(standata$J_1, as.array(rep(1:3, 3)))
  standata2 <- make_standata(y ~ x + (1|x), data = data,
                             control = list(not4stan = TRUE))
  expect_equal(colnames(standata2$X), c("Intercept", "xb", "xc"))
})

test_that("make_standata returns correct data names for addition terms", {
  dat <- data.frame(y = 1:10, w = 1:10, t = 1:10, x = rep(0,10),
                          c = sample(-1:1,10,TRUE))
  expect_equal(names(make_standata(y | se(w) ~ x, dat, gaussian())),
               c("N", "Y", "se", "K", "X", "sigma", "prior_only"))
  expect_equal(names(make_standata(y | weights(w) ~ x, dat, "gaussian")),
               c("N", "Y", "weights", "K", "X",  "prior_only"))
  expect_equal(names(make_standata(y | cens(c) ~ x, dat, "student")),
               c("N", "Y", "cens", "K", "X", "prior_only"))
  expect_equal(names(make_standata(y | trials(t) ~ x, dat, "binomial")),
               c("N", "Y", "trials", "K", "X", "prior_only"))
  expect_equal(names(make_standata(y | trials(10) ~ x, dat, "binomial")),
               c("N", "Y", "trials", "K", "X", "prior_only"))
  expect_equal(names(make_standata(y | thres(11) ~ x, dat, "acat")),
               c("N", "Y", "nthres", "K", "X", "disc", "prior_only"))
  expect_equal(names(make_standata(y | thres(10) ~ x, dat, cumulative())),
               c("N", "Y", "nthres", "K", "X", "disc", "prior_only"))
  sdata <- make_standata(y | trunc(0,20) ~ x, dat, "gaussian")
  expect_true(all(sdata$lb == 0) && all(sdata$ub == 20))
  sdata <- make_standata(y | trunc(ub = 21:30) ~ x, dat)
  expect_true(all(all(sdata$ub == 21:30)))
})

test_that(paste("make_standata accepts correct response variables",
                "depending on the family"), {
  expect_equal(make_standata(y ~ 1, data = data.frame(y = seq(-9.9,0,0.1)),
                             family = "student")$Y, as.array(seq(-9.9,0,0.1)))
  expect_equal(make_standata(y | trials(10) ~ 1, data = data.frame(y = 1:10),
                             family = "binomial")$Y, as.array(1:10))
  expect_equal(make_standata(y ~ 1, data = data.frame(y = 10:20),
                             family = "poisson")$Y, as.array(10:20))
  expect_equal(make_standata(y ~ 1, data = data.frame(y = rep(-c(1:2),5)),
                             family = "bernoulli")$Y, as.array(rep(1:0,5)))
  expect_equal(make_standata(y ~ 1, data = data.frame(y = rep(c(TRUE, FALSE),5)),
                             family = "bernoulli")$Y, as.array(rep(1:0,5)))
  expect_equal(make_standata(y ~ 1, data = data.frame(y = rep(1,5)),
                             family = "bernoulli")$Y, as.array(rep(1, 5)))
  expect_equal(make_standata(y ~ 1, data = data.frame(y = rep(0,5)),
                             family = "bernoulli")$Y, as.array(rep(0, 5)))
  expect_equal(make_standata(y ~ 1, data = data.frame(y = rep(1:10,5)),
                             family = "categorical")$Y, as.array(rep(1:10,5)))
  expect_equal(make_standata(y ~ 1, data = data.frame(y = rep(11:20,5)),
                             family = "categorical")$Y, as.array(rep(1:10,5)))
  expect_equal(make_standata(y ~ 1, data = data.frame(y = factor(rep(11:20,5))),
                             family = "categorical")$Y, as.array(rep(1:10,5)))
  expect_equal(make_standata(y ~ 1, data = data.frame(y = rep(1:10,5)),
                             family = "cumulative")$Y, as.array(rep(1:10,5)))
  dat <- data.frame(y = factor(rep(-4:5,5), order = TRUE))
  expect_equal(make_standata(y ~ 1, data = dat, family = "acat")$Y,
               as.array(rep(1:10,5)))
  expect_equal(make_standata(y ~ 1, data = data.frame(y = seq(1,10,0.1)),
                             family = "exponential")$Y, as.array(seq(1,10,0.1)))

  dat <- data.frame(y1 = 1:10, y2 = 11:20, x = rep(0,10))
  form <- bf(mvbind(y1, y2) ~ x) + set_rescor(TRUE)
  sdata <- make_standata(form, data = dat)
  expect_equal(sdata$Y_y1, as.array(1:10))
  expect_equal(sdata$Y_y2, as.array(11:20))
})

test_that(paste("make_standata rejects incorrect response variables",
                "depending on the family"), {
  expect_error(make_standata(y ~ 1, data = data.frame(y = factor(1:10)),
                             family = "student"),
               "Family 'student' requires numeric responses")
  expect_error(make_standata(y ~ 1, data = data.frame(y = -5:5),
                             family = "geometric"),
               "Family 'geometric' requires response greater than or equal to 0")
  expect_error(make_standata(y ~ 1, data = data.frame(y = -1:1),
                             family = "bernoulli"),
               "contain only two different values")
  expect_error(make_standata(y ~ 1, data = data.frame(y = factor(-1:1)),
                             family = "cratio"),
               "Family 'cratio' requires either positive integers or ordered factors")
  expect_error(make_standata(y ~ 1, data = data.frame(y = rep(0.5:7.5), 2),
                             family = "sratio"),
               "Family 'sratio' requires either positive integers or ordered factors")
  expect_error(make_standata(y ~ 1, data = data.frame(y = rep(-7.5:7.5), 2),
                             family = "gamma"),
               "Family 'gamma' requires response greater than 0")
  expect_error(make_standata(y ~ 1, data = data.frame(y = c(0.1, 0.5, 1)),
                             family = Beta()),
               "Family 'beta' requires response smaller than 1")
  expect_error(make_standata(y ~ 1, data = data.frame(y = c(0, 0.5, 4)),
                             family = von_mises()),
               "Family 'von_mises' requires response smaller than or equal to 3.14")
  expect_error(make_standata(y ~ 1, data = data.frame(y = c(-1, 2, 5)),
                             family = hurdle_gamma()),
               "Family 'hurdle_gamma' requires response greater than or equal to 0")
})

test_that("make_standata suggests using family bernoulli if appropriate", {
  expect_message(make_standata(y | trials(1) ~ 1, data = list(y = rep(0:1,5)),
                               family = "binomial"),
                 "family 'bernoulli' might be a more efficient choice.")
  expect_message(make_standata(y ~ 1, data = data.frame(y = rep(1:2, 5)),
                               family = "acat"),
                 "family 'bernoulli' might be a more efficient choice.")
  expect_message(make_standata(y ~ 1, data = data.frame(y = rep(0:1,5)),
                             family = "categorical"),
                "family 'bernoulli' might be a more efficient choice.")
})

test_that("make_standata returns correct values for addition terms", {
  dat <- data.frame(y = rnorm(9), s = 1:9, w = 1:9, c1 = rep(-1:1, 3),
                    c2 = rep(c("left","none","right"), 3),
                    c3 = c(rep(c(TRUE, FALSE), 4), FALSE),
                    c4 = c(sample(-1:1, 5, TRUE), rep(2, 4)),
                    t = 11:19)
  expect_equivalent(make_standata(y | se(s) ~ 1, data = dat)$se,
                    as.array(1:9))
  expect_equal(make_standata(y | weights(w) ~ 1, data = dat)$weights,
               as.array(1:9))
  expect_equal(make_standata(y | cens(c1) ~ 1, data = dat)$cens,
               as.array(rep(-1:1, 3)))
  expect_equal(make_standata(y | cens(c2) ~ 1, data = dat)$cens,
               as.array(rep(-1:1, 3)))
  expect_equal(make_standata(y | cens(c3) ~ 1, data = dat)$cens,
               as.array(c(rep(1:0, 4), 0)))
  expect_equal(make_standata(y | cens(c4, y + 2) ~ 1, data = dat)$rcens,
               as.array(c(rep(0, 5), dat$y[6:9] + 2)))
  sdata <- suppressWarnings(make_standata(s ~ 1, dat, family = "binomial"))
  expect_equal(sdata$trials, as.array(rep(9, 9)))
  expect_equal(make_standata(s | trials(10) ~ 1, dat,
                             family = "binomial")$trials,
               as.array(rep(10, 9)))
  expect_equal(make_standata(s | trials(t) ~ 1, data = dat,
                             family = "binomial")$trials,
               as.array(11:19))
  expect_equal(SW(make_standata(s | cat(19) ~ 1, data = dat,
                      family = "cumulative"))$nthres,
               18)
})

test_that("make_standata rejects incorrect addition terms", {
  dat <- data.frame(y = rnorm(9), s = -(1:9), w = -(1:9),
                    c = rep(-2:0, 3), t = 9:1, z = 1:9)
  expect_error(make_standata(y | se(s) ~ 1, data = dat),
               "Standard errors must be non-negative")
  expect_error(make_standata(y | weights(w) ~ 1, data = dat),
               "Weights must be non-negative")
  expect_error(make_standata(y | cens(c) ~ 1, data = dat))
  expect_error(make_standata(z | trials(t) ~ 1, data = dat,
                             family = "binomial"),
               "Number of trials is smaller than the number of events")
})

test_that("make_standata handles multivariate models", {
  dat <- data.frame(
    y1 = 1:10, y2 = 11:20,
    x = rep(0, 10), g = rep(1:2, 5),
    censi = sample(0:1, 10, TRUE),
    tim = 10:1, w = 1:10
  )

  form <- bf(mvbind(y1, y2) | weights(w) ~ x) + set_rescor(TRUE)
  sdata <- make_standata(form, data = dat)
  expect_equal(sdata$Y_y1, as.array(dat$y1))
  expect_equal(sdata$Y_y2, as.array(dat$y2))
  expect_equal(sdata$weights_y1, as.array(1:10))

  expect_error(make_standata(bf(mvbind(y1, y2, y2) ~ x) + set_resor(FALSE),
                             data = dat),
               "Cannot use the same response variable twice")

  form <- bf(mvbind(y1 / y2, y2, y1 * 3) ~ x) + set_rescor(FALSE)
  sdata <- make_standata(form, data = dat)
  expect_equal(sdata$Y_y1y2, as.array(dat$y1 / dat$y2))

  sdata <- suppressWarnings(
    make_standata(mvbind(y1, y2) ~ x, dat, autocor = cor_ar(~ tim | g))
  )
  target1 <- c(seq(9, 1, -2), seq(10, 2, -2))
  expect_equal(sdata$Y_y1, as.array(target1))
  target2 <- c(seq(19, 11, -2), seq(20, 12, -2))
  expect_equal(sdata$Y_y2, as.array(target2))

  # models without residual correlations
  expect_warning(
    bform <- bf(y1 | cens(censi) ~ x + y2 + (x|2|g)) +
      gaussian() + cor_ar() +
      (bf(x ~ 1) + mixture(poisson, nmix = 2)) +
      (bf(y2 ~ s(y2) + (1|2|g)) + skew_normal()),
    "Using 'cor_brms' objects for 'autocor' is deprecated"
  )
  bprior <- prior(normal(0, 5), resp = y1) +
    prior(normal(0, 10), resp = y2) +
    prior(dirichlet(2, 1), theta, resp = x)
  sdata <- make_standata(bform, dat, prior = bprior)
  sdata_names <- c(
    "N", "J_1_y1",  "cens_y1", "Kma_y1", "Z_1_y2_3",
    "Zs_y2_1_1", "Y_y2", "con_theta_x", "X_mu2_x"
  )
  expect_true(all(sdata_names %in% names(sdata)))
  expect_equal(sdata$con_theta_x, as.array(c(2, 1)))
})

test_that("make_standata removes NAs correctly", {
  dat <- data.frame(y = c(rnorm(9), NA))
  sdata <- suppressWarnings(make_standata(y ~ 1, dat))
  expect_equal(as.numeric(sdata$Y), dat$y[1:9])
})

test_that("make_standata handles the 'subset' addition argument correctly", {
  dat1 <- data.frame(
    y1 = rnorm(15), y2 = NA,
    x1 = rnorm(15), x2 = NA, x3 = rnorm(15),
    sub1 = 1, sub2 = 0
  )
  dat2 <- data.frame(
    y1 = NA, y2 = rnorm(10),
    x1 = NA, x2 = rnorm(10), x3 = NA,
    sub1 = 0, sub2 = 1
  )
  dat <- rbind(dat1, dat2)

  bform <-
    bf(y1 | subset(sub1) ~ x1*x3 + sin(x1), family = gaussian()) +
    bf(y2 | subset(sub2) ~ x2, family = gaussian()) +
    set_rescor(FALSE)

  sdata <- make_standata(bform, dat)
  nsub1 <- sum(dat$sub1)
  nsub2 <- sum(dat$sub2)
  expect_equal(sdata$N_y1, nsub1)
  expect_equal(sdata$N_y2, nsub2)
  expect_equal(length(sdata$Y_y1), nsub1)
  expect_equal(nrow(sdata$X_y2), nsub2)
})

test_that("make_standata returns correct data for ARMA terms", {
  dat <- data.frame(y = 1:10, x = rep(0, 10), tim = 10:1, g = rep(3:4, 5))

  sdata <- make_standata(y ~ x + ma(tim, g), data = dat)
  expect_equal(sdata$J_lag, as.array(c(1, 1, 1, 1, 0, 1, 1, 1, 1, 0)))

  sdata <- make_standata(y ~ x + ar(tim, g, p = 2), data = dat)
  expect_equal(sdata$J_lag, as.array(c(1, 2, 2, 2, 0, 1, 2, 2, 2, 0)))

  sdata <- make_standata(y ~ x + ar(tim, g, cov = TRUE), data = dat)
  expect_equal(sdata$begin_tg, as.array(c(1, 6)))
  expect_equal(sdata$nobs_tg, as.array(c(5, 5)))

  bform <- bf(y ~ exp(b * x), b ~ 1, nl = TRUE, autocor = ~arma())
  sdata <- make_standata(bform, dat)
})

test_that("make_standata allows to retrieve the initial data order", {
  dat <- data.frame(y1 = rnorm(100), y2 = rnorm(100),
                          id = sample(1:10, 100, TRUE),
                          time = sample(1:100, 100))
  # univariate model
  sdata1 <- make_standata(y1 ~ ar(time, id), data = dat, internal = TRUE)
  expect_equal(dat$y1, as.numeric(sdata1$Y[attr(sdata1, "old_order")]))

  # multivariate model
  form <- bf(mvbind(y1, y2) ~ ma(time, id)) + set_rescor(FALSE)
  sdata2 <- make_standata(form, data = dat, internal = TRUE)
  expect_equal(sdata2$Y_y1[attr(sdata2, "old_order")], as.array(dat$y1))
  expect_equal(sdata2$Y_y2[attr(sdata2, "old_order")], as.array(dat$y2))
})

test_that("make_standata handles covariance matrices correctly", {
  A <- structure(diag(1, 4), dimnames = list(1:4, NULL))
  sdata <- make_standata(count ~ Trt + (1|gr(visit, cov = A)),
                         data = epilepsy, data2 = list(A = A))
  expect_equivalent(sdata$Lcov_1, t(chol(A)))

  B <- structure(diag(1:5), dimnames = list(c(1,5,2,4,3), NULL))
  sdata <- make_standata(count ~ Trt + (1|gr(visit, cov = B)),
                         data = epilepsy, data2 = list(B = B))
  expect_equivalent(sdata$Lcov_1, t(chol(B[c(1,3,5,4), c(1,3,5,4)])))

  B <- diag(1, 4)
  expect_error(make_standata(count ~ Trt + (1|gr(visit, cov = B)),
                             data = epilepsy, data2 = list(B = B)),
               "Row or column names are required")

  B <- structure(diag(1, 4), dimnames = list(2:5, NULL))
  expect_error(make_standata(count ~ Trt + (1|gr(visit, cov = B)),
                             data = epilepsy, data2 = list(B = B)),
               "Levels of .* do not match")

  B <- A
  B[1,2] <- 0.5
  expect_error(make_standata(count ~ Trt + (1|gr(visit, cov = B)),
                             data = epilepsy,  data2 = list(B = B)),
               "must be symmetric")

  expect_warning(
    sdata <- make_standata(count ~ Trt + (1|visit), data = epilepsy,
                           cov_ranef = list(visit = A)),
    "Argument 'cov_ranef' is deprecated"
  )
  expect_equivalent(sdata$Lcov_1, t(chol(A)))
})

test_that("make_standata correctly prepares data for non-linear models", {
  flist <- list(a ~ x + (1|1|g), b ~ mo(z) + (1|1|g))
  dat <- data.frame(
    y = rnorm(9), x = rnorm(9), z = sample(1:9, 9), g = rep(1:3, 3)
  )
  bform <- bf(y ~ a - b^z, flist = flist, nl = TRUE)
  sdata <- make_standata(bform, data = dat)
  expect_equal(names(sdata),
    c("N", "Y", "C_1", "K_a", "X_a", "Z_1_a_1",
      "K_b", "X_b", "Ksp_b", "Imo_b", "Xmo_b_1", "Jmo_b",
      "con_simo_b_1", "Z_1_b_2", "J_1", "N_1",
      "M_1", "NC_1", "prior_only")
  )
  expect_equal(colnames(sdata$X_a), c("Intercept", "x"))
  expect_equal(sdata$J_1, as.array(dat$g))

  bform <- bf(y ~ x) +
    nlf(sigma ~ a1 * exp(-x/(a2 + z))) +
    lf(a1 ~ 1, a2 ~ z + (x|g)) +
    lf(alpha ~ x)
  sdata <- make_standata(bform, dat, family = skew_normal())
  sdata_names <- c("C_sigma_1", "X_a2", "Z_1_a2_1")
  expect_true(all(sdata_names %in% names(sdata)))
})

test_that("make_standata correctly prepares data for monotonic effects", {
  data <- data.frame(
    y = rpois(120, 10), x1 = rep(1:4, 30), z = rnorm(10),
    x2 = factor(rep(c("a", "b", "c"), 40), ordered = TRUE)
  )
  sdata <- make_standata(y ~ mo(x1)*mo(x2)*y, data = data)
  sdata_names <- c("Xmo_1", "Imo", "Jmo",  "con_simo_8", "con_simo_5")
  expect_true(all(sdata_names %in% names(sdata)))
  expect_equivalent(sdata$Xmo_1, as.array(data$x1 - 1))
  expect_equivalent(sdata$Xmo_2, as.array(as.numeric(data$x2) - 1))
  expect_equal(
    as.vector(unname(sdata$Jmo)),
    rep(c(max(data$x1) - 1, length(unique(data$x2)) - 1), 4)
  )
  expect_equal(sdata$con_simo_1, as.array(rep(1, 3)))

  prior <- set_prior("dirichlet(1:3)", coef = "mox11",
                     class = "simo", dpar = "sigma")
  sdata <- make_standata(bf(y ~ 1, sigma ~ mo(x1)),
                         data = data, prior = prior)
  expect_equal(sdata$con_simo_sigma_1, as.array(1:3))

  prior <- c(
    set_prior("normal(0,1)", class = "b", coef = "mox1"),
    set_prior("dirichlet(c(1, 0.5, 2))", class = "simo", coef = "mox11"),
    prior_(~dirichlet(c(1, 0.5, 2)), class = "simo", coef = "mox1:mox21")
  )
  sdata <- make_standata(y ~ mo(x1)*mo(x2), data = data, prior = prior)
  expect_equal(sdata$con_simo_1, as.array(c(1, 0.5, 2)))
  expect_equal(sdata$con_simo_3, as.array(c(1, 0.5, 2)))

  # test issue #924 (conditional monotonicity)
  prior <- c(prior(dirichlet(c(1, 0.5, 2)), simo, coef = "v"),
             prior(dirichlet(c(1,3)), simo, coef = "w"))
  sdata <- make_standata(y ~ y*mo(x1, id = "v")*mo(x2, id = "w"),
                         data, prior = prior)
  expect_equal(sdata$con_simo_1, as.array(c(1, 0.5, 2)))
  expect_equal(sdata$con_simo_2, as.array(c(1, 3)))
  expect_true(!"sdata$con_simo_3" %in% names(sdata))

  expect_error(
    make_standata(y ~ mo(z), data = data),
    "Monotonic predictors must be integers or ordered factors"
  )

  prior <- c(set_prior("dirichlet(c(1,0.5,2))", class = "simo", coef = "mox21"))
  expect_error(
    make_standata(y ~ mo(x2), data = data, prior = prior),
    "Invalid Dirichlet prior"
  )
})

test_that("make_standata returns FCOR covariance matrices", {
  data <- data.frame(y = 1:5)
  data2 <- list(V = diag(5))
  expect_equal(make_standata(y ~ fcor(V), data, data2 = data2)$Mfcor,
               data2$V, check.attributes = FALSE)

  expect_warning(
    expect_error(
      make_standata(y~1, data, autocor = cor_fixed(diag(2))),
      "Dimensions of 'M' for FCOR terms must be equal"
    ),
    "Using 'cor_brms' objects for 'autocor' is deprecated"
  )
})

test_that("make_standata returns data for GAMMs", {
  dat <- data.frame(y = rnorm(10), x1 = rnorm(10), x2 = rnorm(10),
                    x3 = rnorm(10), z = rnorm(10), g = factor(rep(1:2, 5)))
  sdata <- make_standata(y ~ s(x1) + z + s(x2, by = x3), data = dat)
  expect_equal(sdata$nb_1, 1)
  expect_equal(as.vector(sdata$knots_2), 8)
  expect_equal(dim(sdata$Zs_1_1), c(10, 8))
  expect_equal(dim(sdata$Zs_2_1), c(10, 8))

  bform <- bf(y ~ lp, lp ~ s(x1) + z + s(x2, by = x3), nl = TRUE)
  sdata <- make_standata(bform, dat)
  expect_equal(sdata$nb_lp_1, 1)
  expect_equal(as.vector(sdata$knots_lp_2), 8)
  expect_equal(dim(sdata$Zs_lp_1_1), c(10, 8))
  expect_equal(dim(sdata$Zs_lp_2_1), c(10, 8))

  sdata <- make_standata(y ~ g + s(x2, by = g), data = dat)
  expect_true(all(c("knots_1", "knots_2") %in% names(sdata)))

  # test issue #562
  dat$g <- as.character(dat$g)
  sdata <- make_standata(y ~ g + s(x2, by = g), data = dat)
  expect_true(all(c("knots_1", "knots_2") %in% names(sdata)))

  sdata <- make_standata(y ~ t2(x1, x2), data = dat)
  expect_equal(sdata$nb_1, 3)
  expect_equal(as.vector(sdata$knots_1), c(9, 6, 6))
  expect_equal(dim(sdata$Zs_1_1), c(10, 9))
  expect_equal(dim(sdata$Zs_1_3), c(10, 6))

  expect_error(make_standata(y ~ te(x1, x2), data = dat),
               "smooths 'te' and 'ti' are not yet implemented")
})

test_that("make_standata returns correct group ID data", {
  form <- bf(count ~ Trt + (1+Trt|3|visit) + (1|patient),
             shape ~ (1|3|visit) + (Trt||patient))
  sdata <- make_standata(form, data = epilepsy, family = negbinomial())
  expect_true(all(c("Z_1_1", "Z_2_2", "Z_3_shape_1", "Z_2_shape_3") %in%
                    names(sdata)))

  form <- bf(count ~ a, sigma ~ (1|3|visit) + (Trt||patient),
             a ~ Trt + (1+Trt|3|visit) + (1|patient), nl = TRUE)
  sdata <- make_standata(form, data = epilepsy, family = student())
  expect_true(all(c("Z_1_sigma_1", "Z_2_a_3", "Z_2_sigma_1",
                    "Z_3_a_1") %in% names(sdata)))
})

test_that("make_standata handles population-level intercepts", {
  dat <- data.frame(y = 10:1, x = 1:10)
  sdata <- make_standata(y ~ 0 + x, data = dat)
  expect_equal(unname(sdata$X[, 1]), dat$x)

  sdata <- make_standata(y ~ x, dat, cumulative(),
                         control = list(not4stan = TRUE))
  expect_equal(unname(sdata$X[, 1]), dat$x)

  sdata <- make_standata(y ~ 0 + Intercept + x, data = dat)
  expect_equal(unname(sdata$X), cbind(1, dat$x))
})

test_that("make_standata handles category specific effects", {
  sdata <- make_standata(rating ~ period + carry + cse(treat),
                         data = inhaler, family = sratio())
  expect_equivalent(sdata$Xcs, matrix(inhaler$treat))
  sdata <- make_standata(rating ~ period + carry + cs(treat) + (cs(1)|subject),
                         data = inhaler, family = acat())
  expect_equivalent(sdata$Z_1_3, as.array(rep(1, nrow(inhaler))))
  sdata <- make_standata(rating ~ period + carry + (cs(treat)|subject),
                         data = inhaler, family = cratio())
  expect_equivalent(sdata$Z_1_4, as.array(inhaler$treat))
  expect_warning(
    make_standata(rating ~ 1 + cs(treat), data = inhaler,
                  family = "cumulative"),
    "Category specific effects for this family should be considered experimental"
  )
  expect_error(make_standata(rating ~ 1 + (treat + cs(1)|subject),
                             data = inhaler, family = "cratio"),
               "category specific effects in separate group-level terms")
})

test_that("make_standata handles wiener diffusion models", {
  dat <- data.frame(q = 1:10, resp = sample(0:1, 10, TRUE), x = rnorm(10))
  dat$dec <- ifelse(dat$resp == 0, "lower", "upper")
  dat$test <- "a"
  sdata <- make_standata(q | dec(resp) ~ x, data = dat, family = wiener())
  expect_equal(sdata$dec, as.array(dat$resp))
  sdata <- make_standata(q | dec(dec) ~ x, data = dat, family = wiener())
  expect_equal(sdata$dec, as.array(dat$resp))
  expect_error(make_standata(q | dec(test) ~ x, data = dat, family = wiener()),
               "Decisions should be 'lower' or 'upper'")
})

test_that("make_standata handles noise-free terms", {
  N <- 30
  dat <- data.frame(
    y = rnorm(N), x = rnorm(N), z = rnorm(N),
    xsd = abs(rnorm(N, 1)), zsd = abs(rnorm(N, 1)),
    ID = rep(1:5, each = N / 5)
  )
  sdata <- make_standata(
    bf(y ~ me(x, xsd)*me(z, zsd)*x, sigma ~ me(x, xsd)),
    data = dat
  )
  expect_equal(sdata$Xn_1, as.array(dat$x))
  expect_equal(sdata$noise_2, as.array(dat$zsd))
  expect_equal(unname(sdata$Csp_3), as.array(dat$x))
  expect_equal(sdata$Ksp, 6)
  expect_equal(sdata$NCme_1, 1)
})

test_that("make_standata handles noise-free terms with grouping factors", {
  dat <- data.frame(
    y = rnorm(10),
    x1 = rep(1:5, each = 2),
    sdx = rep(1:5, each = 2),
    g = rep(c("b", "c", "a", "d", 1), each = 2)
  )
  sdata <- make_standata(y ~ me(x1, sdx, gr = g), dat)
  expect_equal(unname(sdata$Xn_1), as.array(c(5, 3, 1, 2, 4)))
  expect_equal(unname(sdata$noise_1), as.array(c(5, 3, 1, 2, 4)))

  dat$sdx[2] <- 10
  expect_error(
    make_standata(y ~ me(x1, sdx, gr = g), dat),
    "Measured values and measurement error should be unique"
  )
})

test_that("make_standata handles missing value terms", {
  dat = data.frame(y = rnorm(10), x = rnorm(10), g = 1:10)
  miss <- c(1, 3, 9)
  dat$x[miss] <- NA
  bform <- bf(y ~ mi(x)*g) + bf(x | mi() ~ g) + set_rescor(FALSE)
  sdata <- make_standata(bform, dat)
  expect_equal(sdata$Jmi_x, as.array(miss))
  expect_true(all(is.infinite(sdata$Y_x[miss])))

  # dots in variable names are correctly handled #452
  dat$x.2 <- dat$x
  bform <- bf(y ~ mi(x.2)*g) + bf(x.2 | mi() ~ g) + set_rescor(FALSE)
  sdata <- make_standata(bform, dat)
  expect_equal(sdata$Jmi_x, as.array(miss))

  dat$z <- rbeta(10, 1, 1)
  dat$z[miss] <- NA
  bform <- bf(exp(y) ~ mi(z)*g) + bf(z | mi() ~ g, family = Beta()) +
    set_rescor(FALSE)
  sdata <- make_standata(bform, dat)
  expect_equal(sdata$Jmi_z, as.array(miss))
})

test_that("make_standata handles overimputation", {
  dat = data.frame(y = rnorm(10), x = rnorm(10), g = 1:10, sdy = 1)
  miss <- c(1, 3, 9)
  dat$x[miss] <- dat$sdy[miss] <- NA
  bform <- bf(y ~ mi(x)*g) + bf(x | mi(sdy) ~ g) + set_rescor(FALSE)
  sdata <- make_standata(bform, dat)
  expect_equal(sdata$Jme_x, as.array(setdiff(1:10, miss)))
  expect_true(all(is.infinite(sdata$Y_x[miss])))
  expect_true(all(is.infinite(sdata$noise_x[miss])))
})

test_that("make_standata handles 'mi' terms with 'subset'", {
  dat <- data.frame(
    y = rnorm(10), x = c(rnorm(9), NA), z = rnorm(10),
    g1 = sample(1:5, 10, TRUE), g2 = 10:1, g3 = 1:10,
    s = c(FALSE, rep(TRUE, 9))
  )

  bform <- bf(y ~ mi(x, idx = g1)) +
    bf(x | mi() + index(g2) + subset(s)  ~ 1) +
    set_rescor(FALSE)
  sdata <- make_standata(bform, dat)
  expect_true(all(sdata$idxl_y_x_1 %in% 9:5))

  # test a bunch of errors
  # fails on CRAN for some reason
  # bform <- bf(y ~ mi(x, idx = g1)) +
  #   bf(x | mi() + index(g3) + subset(s) ~ 1) +
  #   set_rescor(FALSE)
  # expect_error(make_standata(bform, dat),
  #   "Could not match all indices in response 'x'"
  # )

  bform <- bf(y ~ mi(x, idx = g1)) +
    bf(x | mi() + subset(s) ~ 1) +
    set_rescor(FALSE)
  expect_error(make_standata(bform, dat),
    "Response 'x' needs to have an 'index' addition term"
  )

  bform <- bf(y ~ mi(x)) +
    bf(x | mi() + subset(s) + index(g2)  ~ 1) +
    set_rescor(FALSE)
  expect_error(make_standata(bform, dat),
    "mi() terms of subsetted variables require the 'idx' argument",
    fixed = TRUE
  )

  bform <- bf(y | mi() ~ mi(x, idx = g1)) +
    bf(x | mi() + subset(s) + index(g2)  ~ mi(y)) +
    set_rescor(FALSE)
  expect_error(make_standata(bform, dat),
    "mi() terms in subsetted formulas require the 'idx' argument",
    fixed = TRUE
  )
})

test_that("make_standata handles multi-membership models", {
  dat <- data.frame(y = rnorm(10), g1 = c(7:2, rep(10, 4)),
                    g2 = 1:10, w1 = rep(1, 10),
                    w2 = rep(abs(rnorm(10))))
  sdata <- make_standata(y ~ (1|mm(g1,g2,g1,g2)), data = dat)
  expect_true(all(paste0(c("W_1_", "J_1_"), 1:4) %in% names(sdata)))
  expect_equal(sdata$W_1_4, as.array(rep(0.25, 10)))
  expect_equal(unname(sdata$Z_1_1_1), as.array(rep(1, 10)))
  expect_equal(unname(sdata$Z_1_1_2), as.array(rep(1, 10)))
  # this checks whether combintation of factor levels works as intended
  expect_equal(sdata$J_1_1, as.array(c(6, 5, 4, 3, 2, 1, 7, 7, 7, 7)))
  expect_equal(sdata$J_1_2, as.array(c(8, 1, 2, 3, 4, 5, 6, 9, 10, 7)))

  sdata <- make_standata(y ~ (1|mm(g1,g2, weights = cbind(w1, w2))), dat)
  expect_equal(sdata$W_1_1, as.array(dat$w1 / (dat$w1 + dat$w2)))

  # tests mmc terms
  sdata <- make_standata(y ~ (1+mmc(w1, w2)|mm(g1,g2)), data = dat)
  expect_equal(unname(sdata$Z_1_2_1), as.array(dat$w1))
  expect_equal(unname(sdata$Z_1_2_2), as.array(dat$w2))

  expect_error(
    make_standata(y ~ (mmc(w1, w2, y)|mm(g1,g2)), data = dat),
    "Invalid term 'mmc(w1, w2, y)':", fixed = TRUE
  )
  expect_error(
    make_standata(y ~ (mmc(w1, w2)*y|mm(g1,g2)), data = dat),
    "The term 'mmc(w1, w2):y' is invalid", fixed = TRUE
  )

  # tests if ":" works in multi-membership models
  sdata <- make_standata(y ~ (1|mm(w1:g1,w1:g2)), dat)
  expect_true(all(c("J_1_1", "J_1_2") %in% names(sdata)))
})

test_that("by variables in grouping terms are handled correctly", {
  gvar <- c("1A", "1B", "2A", "2B", "3A", "3B", "10", "100", "2", "3")
  gvar <- rep(gvar, each = 10)
  g_order <- order(gvar)
  byvar <- c(0, 4.5, 3, 2, "x 1")
  byvar <- factor(rep(byvar, each = 20))
  dat <- data.frame(
    y = rnorm(100), x = rnorm(100),
    g = gvar,
    g2 = gvar[g_order],
    z = byvar,
    z2 = byvar[g_order],
    z3 = factor(1:2)
  )
  sdata <- make_standata(y ~ x + (x | gr(g, by = z)), dat)
  expect_equal(sdata$Nby_1, 5)
  expect_equal(sdata$Jby_1, as.array(c(2, 2, 1, 1, 5, 4, 4, 5, 3, 3)))

  sdata <- make_standata(y ~ x + (x | mm(g, g2, by = cbind(z, z2))), dat)
  expect_equal(sdata$Nby_1, 5)
  expect_equal(sdata$Jby_1, as.array(c(2, 2, 1, 1, 5, 4, 4, 5, 3, 3)))

  expect_error(make_standata(y ~ x + (1|gr(g, by = z3)), dat),
               "Some levels of 'g' correspond to multiple levels of 'z3'")
})

test_that("make_standata handles calls to the 'poly' function", {
  dat <- data.frame(y = rnorm(10), x = rnorm(10))
  expect_equal(colnames(make_standata(y ~ 1 + poly(x, 3), dat)$X),
               c("Intercept", "polyx31", "polyx32", "polyx33"))
})

test_that("make_standata allows fixed distributional parameters", {
  dat <- list(y = 1:10)
  expect_equal(make_standata(bf(y ~ 1, nu = 3), dat, student())$nu, 3)
  expect_equal(make_standata(y ~ 1, dat, acat())$disc, 1)
  expect_error(make_standata(bf(y ~ 1, bias = 0.5), dat),
               "Invalid fixed parameters: 'bias'")
})

test_that("Cell-mean coding can be disabled", {
  df <- data.frame(y = 1:10, g = rep(c("a", "b"), 5))
  bform <- bf(y ~ g) +
    lf(disc ~ 0 + g + (0 + g | y), cmc = FALSE) +
    cumulative()

  sdata <- make_standata(bform, df)
  target <- matrix(rep(0:1, 5), dimnames = list(1:10, "gb"))
  expect_equal(sdata$X_disc, target)
  expect_equal(unname(sdata$Z_1_disc_1), as.array(rep(0:1, 5)))
  expect_true(!"Z_1_disc_2" %in% names(sdata))

  bform <- bf(y ~ 0 + g + (1 | y), cmc = FALSE)
  sdata <- make_standata(bform, df)
  expect_equal(sdata$X, target)
  expect_equal(unname(sdata$Z_1_1), as.array(rep(1, 10)))
})

test_that("make_standata correctly includes offsets", {
  data <- data.frame(y = rnorm(10), x = rnorm(10), c = 1)
  sdata <- make_standata(bf(y ~ x + offset(c), sigma ~ offset(c + 1)), data)
  expect_equal(sdata$offsets, as.array(data$c))
  expect_equal(sdata$offsets_sigma, as.array(data$c + 1))
  sdata <- make_standata(y ~ x + offset(c) + offset(x), data)
  expect_equal(sdata$offsets, as.array(data$c + data$x))
})

test_that("make_standata includes data for mixture models", {
  data <- data.frame(y = rnorm(10), x = rnorm(10), c = 1)
  form <- bf(y ~ x, mu1 ~ 1, family = mixture(gaussian, gaussian))
  sdata <- make_standata(form, data)
  expect_equal(sdata$con_theta, as.array(c(1, 1)))
  expect_equal(dim(sdata$X_mu1), c(10, 1))
  expect_equal(dim(sdata$X_mu2), c(10, 2))

  form <- bf(y ~ x, family = mixture(gaussian, gaussian))
  sdata <- make_standata(form, data, prior = prior(dirichlet(10, 2), theta))
  expect_equal(sdata$con_theta, as.array(c(10, 2)))

  form <- bf(y ~ x, theta1 = 1, theta2 = 3, family = mixture(gaussian, gaussian))
  sdata <- make_standata(form, data)
  expect_equal(sdata$theta1, 1/4)
  expect_equal(sdata$theta2, 3/4)
})

test_that("make_standata includes data for Gaussian processes", {
  dat <- data.frame(y = rnorm(10), x1 = rnorm(10),
                    z = factor(c(2, 2, 2, 3, 4, rep(5, 5))))
  sdata <- make_standata(y ~ gp(x1), dat)
  expect_equal(max(sdata$Xgp_1) - min(sdata$Xgp_1), 1)
  sdata <- make_standata(y ~ gp(x1, scale = FALSE), dat)
  expect_equal(max(sdata$Xgp_1) - min(sdata$Xgp_1), max(dat$x1) - min(dat$x1))

  sdata <- SW(make_standata(y ~ gp(x1, by = z, gr = TRUE, scale = FALSE), dat))
  expect_equal(sdata$Igp_1_2, as.array(4))
  expect_equal(sdata$Jgp_1_4, as.array(1:5))
  expect_equal(sdata$Igp_1_4, as.array(6:10))

  sdata <- SW(make_standata(y ~ gp(x1, by = y, gr = TRUE), dat))
  expect_equal(sdata$Cgp_1, as.array(dat$y))
})

test_that("make_standata includes data for approximate Gaussian processes", {
  dat <- data.frame(y = rnorm(10), x1 = sample(1:10, 10),
                    z = factor(c(2, 2, 2, 3, 4, rep(5, 5))))

  sdata <- make_standata(y ~ gp(x1, k = 5, c = 5/4), dat)
  expect_equal(sdata$NBgp_1, 5)
  expect_equal(dim(sdata$Xgp_1), c(10, 5))
  expect_equal(dim(sdata$slambda_1), c(5, 1))

  sdata <- SW(make_standata(y ~ gp(x1, by = z, k = 5, c = 5/4, scale = FALSE), dat))
  expect_equal(sdata$Igp_1_2, as.array(4))
  expect_equal(sdata$Cgp_1_2, as.array(1))
  expect_equal(sdata$Igp_1_4, as.array(6:10))
})

test_that("make_standata includes data for SAR models", {
  dat <- data.frame(y = rnorm(10), x = rnorm(10))
  W <- matrix(0, nrow = 10, ncol = 10)
  dat2 <- list(W = W)

  sdata <- make_standata(y ~ x + sar(W), data = dat, data2 = dat2)
  expect_equal(dim(sdata$M), rep(nrow(W), 2))

  dat2 <- list(W = matrix(0, 2, 2))
  expect_error(
    make_standata(y ~ x + sar(W), data = dat, data2 = dat2),
    "Dimensions of 'M' for SAR terms must be equal"
  )
})

test_that("make_standata includes data for CAR models", {
  dat = data.frame(y = rnorm(10), x = rnorm(10), obs = 1:10)
  edges <- cbind(1:10, 10:1)
  W <- matrix(0, nrow = 10, ncol = 10)
  for (i in seq_len(nrow(edges))) {
    W[edges[i, 1], edges[i, 2]] <- 1
  }
  rownames(W) <- 1:nrow(W)
  dat2 <- list(W = W)

  sdata <- make_standata(y ~ x + car(W, gr = obs), dat, data2 = dat2)
  expect_equal(sdata$Nloc, 10)
  expect_equal(unname(sdata$Nneigh), rep(1, 10))
  expect_equal(unname(sdata$edges1), as.array(10:6))
  expect_equal(unname(sdata$edges2), as.array(1:5))

  sdata_old <- SW(make_standata(y ~ x, dat, autocor = cor_car(W)))
  expect_equal(sdata, sdata_old)

  rownames(dat2$W) <- c("a", 2:9, "b")
  dat$group <- rep(c("a", "b"), each = 5)
  sdata <- make_standata(y ~ x + car(W, gr = group), dat, data2 = dat2)
  expect_equal(sdata$Nloc, 2)
  expect_equal(sdata$edges1, as.array(2))
  expect_equal(sdata$edges2, as.array(1))

  sdata <- make_standata(y ~ x + car(W, group, type = "bym2"),
                         data = dat, data2 = dat2)
  expect_equal(length(sdata$car_scale), 1L)

  dat2$W[1, 10] <- 4
  dat2$W[10, 1] <- 4
  expect_message(make_standata(y ~ car(W, gr = group), dat, data2 = dat2),
               "Converting all non-zero values in 'M' to 1")

  # test error messages
  rownames(dat2$W) <- c(1:9, "a")
  expect_error(make_standata(y ~ car(W, gr = group), dat, data2 = dat2),
               "Row names of 'M' for CAR terms do not match")
  rownames(dat2$W) <- NULL
  expect_error(make_standata(y ~ car(W, gr = group), dat, data2 = dat2),
               "Row names are required for 'M'")
  dat2$W[1, 10] <- 0
  expect_error(make_standata(y ~ car(W), dat, data2 = dat2),
               "'M' for CAR terms must be symmetric")
  dat2$W[10, 1] <- 0
  expect_error(SW(make_standata(y ~ x + car(W), dat, data2 = dat2)),
               "all locations should have at least one neighbor")
})

test_that("make_standata includes data of special priors", {
  dat <- data.frame(y = 1:10, x1 = rnorm(10), x2 = rnorm(10))

  # horseshoe prior
  hs <- horseshoe(7, scale_global = 2, df_global = 3,
                  df_slab = 6, scale_slab = 3)
  sdata <- make_standata(y ~ x1*x2, data = dat,
                         prior = set_prior(hs))
  expect_equal(sdata$hs_df, 7)
  expect_equal(sdata$hs_df_global, 3)
  expect_equal(sdata$hs_df_slab, 6)
  expect_equal(sdata$hs_scale_global, 2)
  expect_equal(sdata$hs_scale_slab, 3)

  hs <- horseshoe(par_ratio = 0.1)
  sdata <- make_standata(y ~ x1*x2, data = dat, prior = set_prior(hs))
  expect_equal(sdata$hs_scale_global, 0.1 / sqrt(nrow(dat)))

  # R2D2 prior
  sdata <- make_standata(y ~ x1*x2, data = dat,
                         prior = prior(R2D2(0.5, 10)))
  expect_equal(sdata$R2D2_mean_R2, 0.5)
  expect_equal(sdata$R2D2_prec_R2, 10)
  expect_equal(sdata$R2D2_cons_D2, as.array(rep(1, 3)))

  # lasso prior
  sdata <- make_standata(y ~ x1*x2, data = dat,
                         prior = prior(lasso(2, scale = 10)))
  expect_equal(sdata$lasso_df, 2)
  expect_equal(sdata$lasso_scale, 10)

  # horseshoe and lasso prior applied in a non-linear model
  hs_a1 <- horseshoe(7, scale_global = 2, df_global = 3)
  lasso_a2 <- lasso(2, scale = 10)
  sdata <- make_standata(
    bf(y ~ a1 + a2, a1 ~ x1, a2 ~ 0 + x2, nl = TRUE),
    data = dat, sample_prior = TRUE,
    prior = c(set_prior(hs_a1, nlpar = "a1"),
              set_prior(lasso_a2, nlpar = "a2"))
  )
  expect_equal(sdata$hs_df_a1, 7)
  expect_equal(sdata$lasso_df_a2, 2)
})

test_that("dots in formula are correctly expanded", {
  dat <- data.frame(y = 1:10, x1 = 1:10, x2 = 1:10)
  sdata <- make_standata(y ~ ., dat)
  expect_equal(colnames(sdata$X), c("Intercept", "x1", "x2"))
})

test_that("argument 'stanvars' is handled correctly", {
  bprior <- prior(normal(mean_intercept, 10), class = "Intercept")
  mean_intercept <- 5
  stanvars <- stanvar(mean_intercept)
  sdata <- make_standata(count ~ Trt, data = epilepsy,
                         prior = bprior, stanvars = stanvars)
  expect_equal(sdata$mean_intercept, 5)

  # define a multi_normal prior with known covariance matrix
  bprior <- prior(multi_normal(M, V), class = "b")
  stanvars <- stanvar(rep(0, 2), "M", scode = "  vector[K] M;") +
    stanvar(diag(2), "V", scode = "  matrix[K, K] V;")
  sdata <- make_standata(count ~ Trt + zBase, epilepsy,
                         prior = bprior, stanvars = stanvars)
  expect_equal(sdata$M, rep(0, 2))
  expect_equal(sdata$V, diag(2))
})

test_that("addition arguments 'vint' and 'vreal' work correctly", {
  dat <- data.frame(size = 10, y = sample(0:10, 20, TRUE), x = rnorm(20))
  beta_binomial2 <- custom_family(
    "beta_binomial2",
    dpars = c("mu", "tau"),
    links = c("logit", "log"),
    lb = c(NA, 0),
    type = "int",
    vars = c("vint1[n]", "vreal1[n]")
  )
  sdata <- make_standata(
    y | vint(size) + vreal(x, size) ~ 1,
    data = dat, family = beta_binomial2,
  )
  expect_equal(sdata$vint1, as.array(rep(10, 20)))
  expect_equal(sdata$vreal1, as.array(dat$x))
  expect_equal(sdata$vreal2, as.array(rep(10, 20)))
})

test_that("reserved variables 'Intercept' is handled correctly", {
  dat <- data.frame(y = 1:10)
  expect_warning(
    sdata <- make_standata(y ~ 0 + intercept, dat),
    "Reserved variable name 'intercept' is deprecated."
  )
  expect_true(all(sdata$X[, "intercept"] == 1))
  sdata <- make_standata(y ~ 0 + Intercept, dat)
  expect_true(all(sdata$X[, "Intercept"] == 1))
})

test_that("data for multinomial and dirichlet models is correct", {
  N <- 15
  dat <- as.data.frame(rdirichlet(N, c(3, 2, 1)))
  names(dat) <- c("y1", "y2", "y3")
  dat$t1 <- round(dat$y1 * rpois(N, 10))
  dat$t2 <- round(dat$y2 * rpois(N, 10))
  dat$t3 <- round(dat$y3 * rpois(N, 10))
  dat$x <- rnorm(N)
  dat$y <- with(dat, cbind(y1, y2, y3))
  dat$t <- with(dat, cbind(t1, t2, t3))
  dat$size <- rowSums(dat$t)

  sdata <- make_standata(t | trials(size) ~ x, dat, multinomial())
  expect_equal(sdata$trials, as.array(dat$size))
  expect_equal(sdata$ncat, 3)
  expect_equal(sdata$Y, unname(dat$t))

  sdata <- make_standata(y ~ x, data = dat, family = dirichlet())
  expect_equal(sdata$ncat, 3)
  expect_equal(sdata$Y, unname(dat$y))

  expect_error(
    make_standata(t | trials(10) ~ x, data = dat, family = multinomial()),
    "Number of trials does not match the number of events"
  )
  expect_error(make_standata(t ~ x, data = dat, family = dirichlet()),
               "Response values in simplex models must sum to 1")
})

test_that("make_standata handles cox models correctly", {
  data <- data.frame(y = rexp(100), x = rnorm(100))
  bform <- bf(y ~ x)
  bprior <- prior(dirichlet(3), sbhaz)
  sdata <- make_standata(bform, data, brmsfamily("cox"), prior = bprior)
  expect_equal(dim(sdata$Zbhaz), c(100, 5))
  expect_equal(dim(sdata$Zcbhaz), c(100, 5))
  expect_equal(sdata$con_sbhaz, as.array(rep(3, 5)))

  sdata <- make_standata(bform, data, brmsfamily("cox", bhaz = list(df = 6)))
  expect_equal(dim(sdata$Zbhaz), c(100, 6))
  expect_equal(dim(sdata$Zcbhaz), c(100, 6))
})

test_that("make_standata handles addition term 'rate' is correctly", {
  data <- data.frame(y = rpois(10, 1), x = rnorm(10), time = 1:10)
  sdata <- make_standata(y | rate(time) ~ x, data, poisson())
  expect_equal(sdata$denom, as.array(data$time))
})

test_that("make_standata handles grouped ordinal thresholds correctly", {
  dat <- data.frame(
    y = c(1:5, 1:4, 4),
    gr = rep(c("a", "b"), each = 5),
    th = rep(5:6, each = 5),
    x = rnorm(10)
  )

  # thresholds without a grouping factor
  sdata <- make_standata(y ~ x, dat, cumulative())
  expect_equal(sdata$nthres, 4)

  sdata <- make_standata(y | thres(5) ~ x, dat, cumulative())
  expect_equal(sdata$nthres, 5)

  expect_error(
    make_standata(y | thres(th) ~ x, dat, cumulative()),
    "Number of thresholds needs to be a single value"
  )

  # thresholds with a grouping factor
  sdata <- make_standata(y | thres(th, gr) ~ x, dat, cumulative())
  expect_equal(sdata$nthres, as.array(c(5, 6)))
  expect_equal(sdata$ngrthres, 2)
  expect_equal(unname(sdata$Jthres[1, ]), c(1, 5))
  expect_equal(unname(sdata$Jthres[10, ]), c(6, 11))

  sdata <- make_standata(y | thres(gr = gr) ~ x, dat, cumulative())
  expect_equal(sdata$nthres, as.array(c(4, 3)))
  expect_equal(sdata$ngrthres, 2)

  sdata <- make_standata(y | thres(6, gr = gr) ~ x, dat, cumulative())
  expect_equal(sdata$nthres, as.array(c(6, 6)))
  expect_equal(sdata$ngrthres, 2)
})

test_that("information for threading is handled correctly", {
  dat <- data.frame(y = 1:10)
  sdata <- make_standata(y ~ 1, dat, threads = threading(2, grainsize = 3))
  expect_equal(sdata$grainsize, 3)
})

test_that("variables in data2 can be used in population-level effects", {
  dat <- data.frame(y = 1:10, x1 = rnorm(10), x2 = rnorm(10), x3 = rnorm(10))
  foo <- function(..., idx = NULL) {
    out <- cbind(...)
    if (!is.null(idx)) {
      out <- out[, idx, drop = FALSE]
    }
    out
  }
  sdata <- make_standata(y ~ foo(x1, x2, x3, idx = id), data = dat,
                         data2 = list(id = c(3, 1)))
  target <- c("Intercept", "foox1x2x3idxEQidx3", "foox1x2x3idxEQidx1")
  expect_equal(colnames(sdata$X), target)
  expect_equivalent(sdata$X[, 2], dat$x3)
  expect_equivalent(sdata$X[, 3], dat$x1)
})

test_that("NAs are allowed in unused interval censoring variables", {
  dat <- data.frame(y = rnorm(10), ce = c(1, rep(2, 9)))
  dat$y2 <- dat$y + 2
  dat$y2[1] <- NA
  sdata <- make_standata(y | cens(ce, y2 = y2) ~ 1, data = dat)
  expect_equal(sdata$N, 10L)
  expect_equal(sdata$rcens[1], 0)

  dat$ce[1] <- 2
  expect_error(
    make_standata(y | cens(ce, y2 = y2) ~ 1, data = dat),
    "'y2' should not be NA for interval censored observations"
  )
})

test_that("drop_unused_factor levels works correctly", {
  dat <- data.frame(y = rnorm(10), x = factor(c("a", "b"), levels = c("a", "b", "c")))

  # should drop level "c"
  sdata <- make_standata(y ~ x, data = dat)
  expect_equal(colnames(sdata$X), c("Intercept", "xb"))

  # should not drop level "c"
  sdata <- make_standata(y ~ x, data = dat, drop_unused_levels = FALSE)
  expect_equal(colnames(sdata$X), c("Intercept", "xb", "xc"))
})
