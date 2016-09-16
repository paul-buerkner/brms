test_that(paste("make_standata returns correct data names", 
                "for fixed and random effects"), {
  expect_equal(names(make_standata(rating ~ treat + period + carry 
                                   + (1|subject), data = inhaler)),
               c("N", "Y",  "K", "X", "Z_1_1",
                 "J_1", "N_1", "M_1", "NC_1", "prior_only"))
  expect_equal(names(make_standata(rating ~ treat + period + carry 
                                   + (1+treat|id|subject), data = inhaler,
                                   family = "categorical")),
               c("N", "Y", "K_2", "X_2", "Z_1_2_1", "Z_1_2_2", 
                 "K_3", "X_3", "Z_1_3_3", "Z_1_3_4",
                 "K_4", "X_4", "Z_1_4_5", "Z_1_4_6",
                 "J_1", "N_1", "M_1", "NC_1", "ncat", "max_obs", 
                 "prior_only"))
  expect_equal(names(make_standata(rating ~ treat + period + carry 
                                   + (1+treat|subject), data = inhaler,
                                   control = list(not4stan = TRUE))),
               c("N", "Y", "K", "X", "Z_1", "J_1", "N_1", "M_1",
                 "NC_1", "prior_only"))
  temp_data <- data.frame(y = 1:10, g = 1:10, h = 11:10, x = rep(0,10))
  expect_equal(names(make_standata(y ~ x + (1|g) + (1|h), family = "poisson",
                                   data = temp_data)),
               c("N", "Y", "K", "X", "Z_1_1", "Z_2_1",
                 "J_1", "N_1", "M_1", "NC_1", "J_2", "N_2", "M_2", "NC_2", 
                 "prior_only"))
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

test_that(paste("make_standata returns correct data names", 
                "for addition and cse variables"), {
  temp_data <- data.frame(y = 1:10, w = 1:10, t = 1:10, x = rep(0,10), 
                          c = sample(-1:1,10,TRUE))
  expect_equal(names(make_standata(y | se(w) ~ x, family = "gaussian", 
                                   data = temp_data)), 
               c("N", "Y", "K", "X", "se", "prior_only"))
  expect_equal(names(make_standata(y | weights(w) ~ x, family = "gaussian", 
                                   data = temp_data)), 
               c("N", "Y", "K", "X", "weights", "prior_only"))
  expect_equal(names(make_standata(y | cens(c) ~ x, family = "student", 
                                   data = temp_data)), 
               c("N", "Y", "K", "X", "cens", "prior_only"))
  expect_equal(names(make_standata(y | trials(t) ~ x, family = "binomial", 
                                   data = temp_data)), 
               c("N", "Y", "K", "X", "trials", "max_obs", "prior_only"))
  expect_equal(names(make_standata(y | trials(10) ~ x, family = "binomial", 
                                   data = temp_data)), 
               c("N", "Y", "K", "X", "trials", "max_obs", "prior_only"))
  expect_equal(names(make_standata(y | cat(11) ~ x, family = "acat", 
                                   data = temp_data)), 
               c("N", "Y", "K", "X", "ncat", "max_obs", "prior_only"))
  expect_equal(names(make_standata(y | cat(10) ~ x, family = "cumulative", 
                                   data = temp_data)), 
               c("N", "Y", "K", "X", "ncat", "max_obs", "prior_only"))
  standata <- make_standata(y | trunc(0,20) ~ x, family = "gaussian", 
                            data = temp_data)
  expect_true(all(standata$lb == 0) && all(standata$ub == 20))
  standata <- make_standata(y | trunc(ub = 21:30) ~ x, family = "gaussian", 
                            data = temp_data)
  expect_true(all(all(standata$ub == 21:30)))
})

test_that(paste("make_standata accepts correct response variables", 
                "depending on the family"), {
  expect_equal(make_standata(y ~ 1, data = data.frame(y = seq(-9.9,0,0.1)), 
                             family = "student")$Y, as.array(seq(-9.9,0,0.1)))
  expect_equal(make_standata(y ~ 1, data = data.frame(y = 1:10), 
                             family = "binomial")$Y, as.array(1:10))
  expect_equal(make_standata(y ~ 1, data = data.frame(y = 10:20), 
                             family = "poisson")$Y, as.array(10:20))
  expect_equal(make_standata(y ~ 1, data = data.frame(y = rep(-c(1:2),5)), 
                             family = "bernoulli")$Y, as.array(rep(1:0,5)))
  expect_equal(make_standata(y ~ 1, data = data.frame(y = rep(c(TRUE, FALSE),5)),
                             family = "bernoulli")$Y, as.array(rep(1:0,5)))
  expect_equal(make_standata(y ~ 1, data = data.frame(y = rep(1:10,5)), 
                             family = "categorical")$Y, as.array(rep(1:10,5)))
  expect_equal(make_standata(y ~ 1, data = data.frame(y = rep(-4:5,5)), 
                             family = "categorical")$Y, as.array(rep(1:10,5)))
  expect_equal(make_standata(y ~ 1, data = data.frame(y = factor(rep(-4:5,5))), 
                             family = "categorical")$Y, as.array(rep(1:10,5)))
  expect_equal(make_standata(y ~ 1, data = data.frame(y = rep(1:10,5)), 
                             family = "cumulative")$Y, as.array(rep(1:10,5)))
  temp_data <- data.frame(y = factor(rep(-4:5,5), order = TRUE))
  expect_equal(make_standata(y ~ 1, data = temp_data, family = "acat")$Y, 
               as.array(rep(1:10,5)))
  expect_equal(make_standata(y ~ 1, data = data.frame(y = seq(1,10,0.1)), 
                             family = "exponential")$Y, as.array(seq(1,10,0.1)))
  temp_data <- data.frame(y1 = 1:10, y2 = 11:20, w = 1:10, x = rep(0,10))
  expect_equal(make_standata(cbind(y1,y2) | weights(w) ~ x, family = "gaussian",
                             data = temp_data)$Y, 
               cbind(1:10, 11:20))
})

test_that(paste("make_standata rejects incorrect response variables", 
                "depending on the family"), {
  expect_error(make_standata(y ~ 1, data = data.frame(y = factor(1:10)), 
                             family = "student"),
               "family student expects numeric response variable")
  expect_error(make_standata(y ~ 1, data = data.frame(y = -5:5), 
                             family = "geometric"),
               "family geometric expects response variable of non-negative integers")
  expect_error(make_standata(y ~ 1, data = data.frame(y = -1:1), 
                             family = "bernoulli"),
               "contain only two different values")
  expect_error(make_standata(y ~ 1, data = data.frame(y = factor(-1:1)), 
                             family = "cratio"),
               "family cratio expects either integers or ordered factors")
  expect_error(make_standata(y ~ 1, data = data.frame(y = rep(0.5:7.5), 2), 
                             family = "sratio"),
               "family sratio expects either integers or ordered factors")
  expect_error(make_standata(y ~ 1, data = data.frame(y = rep(-7.5:7.5), 2), 
                             family = "gamma"),
               "family gamma requires response variable to be positive")
  expect_error(make_standata(y ~ 1, data = data.frame(y = c(0, 0.5, 1)),
                             family = Beta()),
               "requires responses between 0 and 1")
  expect_error(make_standata(y ~ 1, data = data.frame(y = c(0, 0.5, 4)),
                             family = von_mises()),
               "requires responses between -pi and pi")
  expect_error(make_standata(y ~ 1, data = data.frame(y = c(-1, 2, 5)),
                             family = hurdle_gamma()),
               "requires response variable to be non-negative")
})

test_that("make_standata suggests using family bernoulli if appropriate", {
  expect_message(make_standata(y ~ 1, data = data.frame(y = rep(0:1,5)), 
                               family = "binomial"),
                 paste("family bernoulli might be a more efficient choice."))
  expect_message(make_standata(y ~ 1, data = data.frame(y = rep(0:1,5)), 
                               family = "acat"),
                 paste("family bernoulli might be a more efficient choice."))
  expect_error(make_standata(y ~ 1, data = data.frame(y = rep(0:1,5)), 
                             family = "categorical"),
                 paste("At least 3 response categories are required"))
})

test_that("make_standata returns correct values for addition arguments", {
  temp_data <- data.frame(y = rnorm(9), s = 1:9, w = 1:9, c1 = rep(-1:1, 3), 
                          c2 = rep(c("left","none","right"), 3),
                          c3 = c(rep(c(TRUE, FALSE), 4), FALSE),
                          t = 11:19)
  expect_equal(make_standata(y | se(s) ~ 1, data = temp_data)$se, 
               1:9)
  expect_equal(make_standata(y | weights(w) ~ 1, data = temp_data)$weights, 
               1:9)
  expect_equal(make_standata(y | disp(w) ~ 1, data = temp_data)$disp, 
               1:9)
  expect_equal(make_standata(y | cens(c1) ~ 1, data = temp_data)$cens, 
               rep(-1:1, 3))
  expect_equal(make_standata(y | cens(c2) ~ 1, data = temp_data)$cens,
               rep(-1:1, 3))
  expect_equal(make_standata(y | cens(c3) ~ 1, data = temp_data)$cens, 
               c(rep(1:0, 4), 0))
  expect_equal(make_standata(s ~ 1, data = temp_data, 
                             family = "binomial")$max_obs, 9)
  expect_equal(make_standata(s | trials(10) ~ 1, data = temp_data, 
                             family = "binomial")$max_obs, 10)
  expect_equal(make_standata(s | trials(t) ~ 1, data = temp_data, 
                             family = "binomial")$max_obs, 11:19)
  expect_equal(make_standata(s | cat(19) ~ 1, data = temp_data, 
                             family = "cumulative")$ncat, 19)
})

test_that("make_standata rejects incorrect addition arguments", {
  temp_data <- data.frame(y = rnorm(9), s = -(1:9), w = -(1:9), 
                          c = rep(-2:0, 3), t = 9:1, z = 1:9)
  expect_error(make_standata(y | se(s) ~ 1, data = temp_data), 
               "standard errors must be non-negative")
  expect_error(make_standata(y | weights(w) ~ 1, data = temp_data), 
               "weights must be non-negative")
  expect_error(make_standata(y | cens(c) ~ 1, data = temp_data))
  expect_error(make_standata(z | trials(t) ~ 1, data = temp_data, 
                             family = "binomial"),
               "Number of trials is smaller than the response variable")
})

test_that(paste("make_standata handles addition arguments", 
                "and autocorrelation in multinormal models"), {
  temp_data <- data.frame(y1 = 1:10, y2 = 11:20, w = 1:10, x = rep(0,10), 
                          tim = 10:1, g = rep(1:2,5))
  expect_equal(make_standata(cbind(y1,y2) | weights(w) ~ x, 
                             family = "gaussian", data = temp_data)$weights, 
               1:10)
  expect_equal(make_standata(cbind(y1,y2) | weights(w) ~ x, 
                             family = "gaussian", data = temp_data,
                             autocor = cor_ar(~tim | g))$Y,
               cbind(c(seq(9,1,-2), seq(10,2,-2)), 
                     c(seq(19,11,-2), seq(20,12,-2))))
})

test_that(paste("make_standata returns correct data", 
                "for autocorrelations structures"), {
  temp_data <- data.frame(y=1:10, x=rep(0,10), tim=10:1, g = rep(3:4,5))
  expect_equal(make_standata(y ~ x, data = temp_data,
                             autocor = cor_arr(~tim|g))$Yarr,
               cbind(c(0,9,7,5,3,0,10,8,6,4)))
  expect_equal(make_standata(y ~ x, data = temp_data,
                             autocor = cor_arr(~tim|g, r = 2))$Yarr,
               cbind(c(0,9,7,5,3,0,10,8,6,4), c(0,0,9,7,5,0,0,10,8,6)))
  expect_equal(make_standata(y ~ x, data = temp_data,
                             autocor = cor_ma(~tim|g))$tg,
               c(rep(1,5), rep(2,5)))
  expect_equal(make_standata(y ~ x, data = temp_data,
                             autocor = cor_ar(~tim|g))$tg,
               c(rep(1,5), rep(2,5)))
  standata <- make_standata(y ~ x, data = temp_data,
                            autocor = cor_ar(~tim|g, cov = TRUE))
  expect_equal(standata$begin_tg, as.array(c(1, 6)))
  expect_equal(standata$nobs_tg, as.array(c(5, 5)))
  expect_equal(standata$se2, rep(0, 10))
})

test_that("make_standata allows to retrieve the initial data order", {
  temp_data <- data.frame(y1 = rnorm(100), y2 = rnorm(100), 
                          id = sample(1:10, 100, TRUE), 
                          time = sample(1:100, 100))
  # univariate model
  sdata1 <- make_standata(y1 ~ 1, data = temp_data, 
                          autocor = cor_ar(~time|id),
                          control = list(save_order = TRUE))
  expect_equal(temp_data$y1, as.numeric(sdata1$Y[attr(sdata1, "old_order")]))
  # multivariate model
  sdata2 <- make_standata(cbind(y1, y2) ~ 1, data = temp_data, 
                          autocor = cor_ma(~time|id),
                          control = list(save_order = TRUE))
  expect_equal(c(temp_data$y1, temp_data$y2), 
               as.numeric(sdata2$Y[attr(sdata2, "old_order"), ]))
})

test_that("make_standata rejects invalid input for cse effects", {
  expect_error(make_standata(rating ~ 1 + cse(treat), data = inhaler,
                             family = "gaussian"), "only meaningful")
  expect_error(make_standata(rating ~ 1 + cse(1), data = inhaler,
                             family = "acat"), "invalid input")
})

test_that("make_standata handles covariance matrices correctly", {
  A <- structure(diag(1, 4), dimnames = list(1:4, NULL))
  expect_equivalent(make_standata(count ~ Trt_c + (1|visit), data = epilepsy,
                                  cov_ranef = list(visit = A))$Lcov_1, A)
  B <- diag(1, 4)
  expect_error(make_standata(count ~ Trt_c + (1|visit), data = epilepsy,
                             cov_ranef = list(visit = B)),
               "rownames are required")
  B <- structure(diag(1, 4), dimnames = list(2:5, NULL))
  expect_error(make_standata(count ~ Trt_c + (1|visit), data = epilepsy,
                             cov_ranef = list(visit = B)),
               "rownames .* do not match")
  B <- structure(diag(1:5), dimnames = list(c(1,5,2,4,3), NULL))
  expect_equivalent(make_standata(count ~ Trt_c + (1|visit), data = epilepsy,
                             cov_ranef = list(visit = B))$Lcov_1,
                    t(chol(B[c(1,3,5,4), c(1,3,5,4)])))
  B <- A
  B[1,2] <- 0.5
  expect_error(make_standata(count ~ Trt_c + (1|visit), data = epilepsy,
                             cov_ranef = list(visit = B)),
               "not symmetric")
})

test_that("make_standata computes data for inverse.gaussian models", {
  temp_data <- data.frame(y = 1:10, x = rep(0,10), w = 1:10)
  standata <- make_standata(y ~ x, data = temp_data, 
                            family = inverse.gaussian)
  expect_equal(standata$log_Y, sum(log(temp_data$y)))
  expect_equal(as.numeric(standata$sqrt_Y), sqrt(temp_data$y))
  standata <- make_standata(y | weights(w) ~ x, data = temp_data,
                            family = inverse.gaussian)
  expect_equal(as.numeric(standata$log_Y), log(temp_data$y))                         
})

test_that("brmdata is backwards compatible", {
  temp_data <- data.frame(y = 1:10, x = sample(1:5, 10, TRUE))
  expect_identical(SW(brmdata(y ~ x + (1|x), data = temp_data, 
                           family = "poisson")), 
                   make_standata(y ~ x + (1|x), data = temp_data, 
                                 family = "poisson"))
  expect_identical(SW(brmdata(y ~ 1, data = temp_data, 
                              family = "acat", partial = ~ x)), 
                   SW(make_standata(y ~ 1, data = temp_data, 
                                    family = "acat", partial = ~ x)))
})

test_that("make_standata correctly prepares data for non-linear models", {
  nonlinear <- list(a ~ x + (1|1|g), b ~ mono(z) + (1|1|g))
  data <- data.frame(y = rnorm(9), x = rnorm(9), z = sample(1:4, 9, TRUE), 
                     g = rep(1:3, 3))
  standata <- make_standata(y ~ a - b^z, data = data, nonlinear = nonlinear)
  expect_equal(names(standata), c("N", "Y", "KC", "C", "K_a", "X_a", "Z_1_a_1", 
                                  "K_b", "X_b", "Km_b", "Xm_b", "Jm_b", 
                                  "con_simplex_b_1", "Z_1_b_2", "J_1", "N_1", 
                                  "M_1", "NC_1", "prior_only"))
  expect_equal(colnames(standata$X_a), c("Intercept", "x"))
  expect_equal(colnames(standata$C), "z")
  expect_equal(standata$J_1, as.array(data$g))
})

test_that("make_standata correctly prepares data for monotonic effects", {
  data <- data.frame(y = rpois(120, 10), x1 = rep(1:4, 30), 
                     x2 = factor(rep(c("a", "b", "c"), 40), ordered = TRUE))
  sdata <- make_standata(y ~ mono(x1 + x2), data = data)
  expect_true(all(c("Xm", "Jm", "con_simplex_1", "con_simplex_2") %in% names(sdata)))
  expect_equivalent(sdata$Xm, cbind(data$x1 - 1, as.numeric(data$x2) - 1))
  expect_equal(as.vector(unname(sdata$Jm)), 
               c(max(data$x1) - 1, length(unique(data$x2)) - 1))
  expect_equal(sdata$con_simplex_1, rep(1, 3))
  
  prior <- set_prior("dirichlet(1:3)", coef = "x1", 
                     class = "simplex", nlpar = "sigma")
  sdata <- make_standata(bf(y ~ 1, sigma ~ mono(x1)), 
                         data = data, prior = prior)
  expect_equal(sdata$con_simplex_sigma_1, 1:3)
  
  prior <- c(set_prior("normal(0,1)", class = "b", coef = "x"),
             set_prior("dirichlet(c(1,0.5,2))", class = "simplex", coef = "x1"))
  sdata <- make_standata(y ~ monotonic(x1 + x2), data = data, prior = prior)
  expect_equal(sdata$con_simplex_1, c(1,0.5,2))
  
  prior <- c(set_prior("dirichlet(c(1,0.5,2))", class = "simplex", coef = "x2"))
  expect_error(make_standata(y ~ monotonic(x1 + x2), data = data, prior = prior),
               "Invalid dirichlet prior for the simplex of x2", fixed = TRUE)
})

test_that("make_standata returns fixed residual covariance matrices", {
  data <- data.frame(y = 1:5)
  V <- diag(5)
  expect_equal(make_standata(y~1, data, autocor = SW(cor_fixed(V)))$V, V)
  expect_error(make_standata(y~1, data, autocor = cov_fixed(diag(2))),
               "'V' must have the same number of rows as 'data'")
})

test_that("make_standata returns data for bsts models", {
  dat <- data.frame(y = 1:5, g = c(1:3, sample(1:3, 2, TRUE)), t = 1:5)
  expect_equal(make_standata(y~1, data = dat, autocor = cor_bsts(~t|g))$tg,
               sort(dat$g))
  expect_equivalent(make_standata(bf(y~1, sigma ~ 1), data = dat, 
                                  autocor = cor_bsts(~t|g))$X_sigma[, 1],
                    rep(1, nrow(dat)))
})

test_that("make_standata returns data for GAMMs", {
  dat <- data.frame(y = rnorm(10), x1 = rnorm(10), x2 = rnorm(10),
                    z = rnorm(10), g = rep(1:2, 5))
  standata <- make_standata(y ~ s(x1) + z + s(x2) + (1|g), data = dat)
  expect_true(all(c("ns", "knots", "Zs_1", "Zs_2") %in% names(standata)))
  expect_equal(standata$ns, 2)
  expect_equal(as.vector(standata$knots), c(8, 8))
  expect_equal(dim(standata$Zs_1), c(10, 8))
  expect_equal(dim(standata$Zs_2), c(10, 8))
  
  standata <- make_standata(y ~ lp, nonlinear = lp ~ s(x1) + z + s(x2) + (1|g), 
                            data = dat)
  expect_true(all(c("ns_lp", "knots_lp", "Zs_lp_1", "Zs_lp_2") %in% 
                    names(standata)))
  expect_equal(standata$ns_lp, 2)
  expect_equal(as.vector(standata$knots_lp), c(8, 8))
  expect_equal(dim(standata$Zs_lp_1), c(10, 8))
  expect_equal(dim(standata$Zs_lp_2), c(10, 8))
})

test_that("make_standata returns correct group ID data", {
  form <- bf(count ~ Trt_c + (1+Trt_c|3|visit) + (1|patient), 
             shape ~ (1|3|visit) + (Trt_c||patient))
  sdata <- make_standata(form, data = epilepsy, family = negbinomial())
  expect_true(all(c("Z_1_1", "Z_2_2", "Z_3_shape_1", "Z_2_shape_3") %in% 
                    names(sdata)))
  
  form <- bf(count ~ a, sigma ~ (1|3|visit) + (Trt_c||patient),
             nonlinear = a ~ Trt_c + (1+Trt_c|3|visit) + (1|patient))
  sdata <- make_standata(form, data = epilepsy, family = student())
  expect_true(all(c("Z_1_sigma_1", "Z_2_a_3", "Z_2_sigma_1",  
                    "Z_3_a_1") %in% names(sdata)))
})

test_that("make_standata does not center X in models without an intercept", {
  dat <- data.frame(y = rnorm(10), x = 1:10)
  sdata <- make_standata(y~0+x, data = dat)
  expect_equal(unname(sdata$X[, 1]), dat$x)
})

test_that("make_standata handles variable 'intercept' correclty", {
  dat <- data.frame(y = rnorm(10), x = 1:10)
  sdata <- make_standata(y~0+intercept + x, data = dat)
  expect_equal(unname(sdata$X), cbind(1, dat$x))
})
