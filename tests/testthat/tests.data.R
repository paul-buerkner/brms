test_that("melt_data returns data in long format", {
  data <- data.frame(x = rep(c("a","b"), 5), y1 = 1:10, y2 = 11:20, 
                     y3 = 21:30, z = 100:91)
  effects <- extract_effects(y ~ x, family = "poisson")
  expect_equal(melt_data(data, effects = effects, 
                         family = "poisson"), data)
  
  target1 <- data.frame(x = rep(c("a","b"), 10), y1 = rep(1:10, 2), 
                        y2 = rep(11:20, 2), y3 = rep(21:30, 2),
                        z = rep(100:91, 2), 
                        trait = factor(rep(c("y3", "y1"), each = 10),
                                       levels = c("y3", "y1")), 
                        response = c(21:30, 1:10))
  effects <- extract_effects(cbind(y3,y1) ~ x, family = "gaussian")
  expect_equal(melt_data(data, effects = effects, family = "gaussian"), 
               target1)
  
  target2 <- data.frame(x = rep(c("a","b"), 15), y1 = rep(1:10, 3), 
                        y2 = rep(11:20, 3), y3 = rep(21:30, 3),
                        z = rep(100:91, 3),
                        trait = factor(rep(c("y2", "y1", "y3"), each = 10),
                                       levels = c("y2", "y1", "y3")), 
                        response = c(11:20, 1:10, 21:30))
  effects <- extract_effects(cbind(y2,y1,y3) ~ x, family = "gaussian")
  expect_equal(melt_data(data, effects = effects, family = "gaussian"), 
               target2)
})

test_that("melt_data returns expected errors", {
  ee <- extract_effects(y1 ~ x, family = hurdle_poisson())
  data <- data.frame(y1 = rnorm(10), y2 = rnorm(10), x = 1:10)
  expect_error(melt_data(data = NULL, family = hurdle_poisson(), effects = ee),
               "data must be a data.frame for multivariate models", 
               fixed = TRUE)
  data$main <- 1:10 
  expect_error(melt_data(data = data, family = hurdle_poisson(), effects = ee),
               "main is a resevered variable name", 
               fixed = TRUE)
  data$response <- 1:10 
  expect_error(melt_data(data = data, family = hurdle_poisson(), effects = ee),
               "response is a resevered variable name in multivariate models", 
               fixed = TRUE)
  data$trait <- 1:10 
  expect_error(melt_data(data = data, family = hurdle_poisson(), effects = ee),
               "trait is a resevered variable name in multivariate models", 
               fixed = TRUE)
  
  ee <- extract_effects(cbind(y1, y2) ~ x)
  data <- data.frame(y1 = rnorm(10), y2 = rnorm(10), x = 1:10)
  expect_error(melt_data(data = data, family = poisson(), effects = ee),
               "Invalid multivariate model", fixed = TRUE)
})

test_that("combine_groups does the expected", {
  data <- data.frame(x = rep(c("a","b"), 5), y1 = 1:10, 
                     y2 = 11:20, y3 = 21:30, z = 100:91)
  expected <- data 
  expected[["y1:y2"]] <- paste0(data$y1, "_", data$y2)
  expected[["y1:y2:y3"]] <- paste0(data$y1, "_", data$y2, "_", data$y3)
  expect_equal(combine_groups(data, "y1:y2", "y1:y2:y3"), expected)
})

test_that("get_model_matrix removes intercepts correctly", {
  data <- data.frame(x = factor(rep(1:2, 5)), y = 11:20)
  expect_equal(get_model_matrix(y ~ x, data, rm_intercept = TRUE),
               structure(matrix(rep(0:1, 5)), dimnames = list(1:10, "x2")))
})

test_that(paste("arr_design_matrix returns correct design", 
                "matrices for autoregressive effects"), {
  expect_equal(arr_design_matrix(1:10, 0, sort(rep(1:2, 5))), NULL)
  expect_equal(arr_design_matrix(1:10, 1, sort(rep(1:2, 5))), 
               matrix(c(0,1:4.5,0,6:9.5)))
  expect_equal(arr_design_matrix(1:10, 2, sort(rep(1:2, 5))), 
               cbind(c(0,1:4.5,0,6:9), c(0,0,1:3,0,0,6:8)))
})

test_that(paste("make_standata returns correct data names", 
                "for fixed and random effects"), {
  expect_equal(names(make_standata(rating ~ treat + period + carry 
                                   + (1|subject), data = inhaler)),
               c("N","Y","K","X","J_1","N_1","K_1","Z_1","NC_1"))
  expect_equal(names(make_standata(rating ~ treat + period + carry 
                                   + (1+treat|subject), data = inhaler,
                                   family = "categorical")),
               c("N","Y","Kp","Xp","J_1","N_1","K_1",
                 "Z_1","NC_1", "ncat", "max_obs"))
  temp_data <- data.frame(y = 1:10, g = 1:10, h = 11:10, x = rep(0,10))
  expect_equal(names(make_standata(y ~ x + (1|g) + (1|h), family = "poisson",
                                   data = temp_data)),
               c("N","Y","K","X","J_1","N_1","K_1","Z_1","NC_1",
                 "J_2","N_2","K_2","Z_2","NC_2"))
})

test_that(paste("make_standata handles variables used as fixed effects", 
                "and grouping factors at the same time"), {
  data <- data.frame(y = 1:9, x = factor(rep(c("a","b","c"), 3)))
  standata <- make_standata(y ~ x + (1|x), data = data)
  expect_equal(colnames(standata$X), c("xb", "xc"))
  expect_equal(standata$J_1, rep(1:3, 3))
  standata2 <- make_standata(y ~ x + (1|x), data = data, 
                             control = list(keep_intercept = TRUE))
  expect_equal(colnames(standata2$X), c("Intercept", "xb", "xc"))
})

test_that(paste("make_standata returns correct data names", 
                "for addition and partial variables"), {
  temp_data <- data.frame(y = 1:10, w = 1:10, t = 1:10, x = rep(0,10), 
                          c = sample(-1:1,10,TRUE))
  expect_equal(names(make_standata(y | se(w) ~ x, family = "gaussian", 
                                   data = temp_data)), 
               c("N","Y","K","X","se"))
  expect_equal(names(make_standata(y | weights(w) ~ x, family = "gaussian", 
                                   data = temp_data)), 
               c("N","Y","K","X","weights"))
  expect_equal(names(make_standata(y | cens(c) ~ x, family = "cauchy", 
                                   data = temp_data)), 
               c("N","Y","K","X","cens"))
  expect_equal(names(make_standata(y | trials(t) ~ x, family = "binomial", 
                                   data = temp_data)), 
               c("N","Y","K","X","trials","max_obs"))
  expect_equal(names(make_standata(y | trials(10) ~ x, family = "binomial", 
                                   data = temp_data)), 
               c("N","Y","K","X","trials","max_obs"))
  expect_equal(names(make_standata(y | cat(11) ~ x, family = "acat", 
                                   data = temp_data)), 
               c("N","Y","K","X","ncat","max_obs"))
  expect_equal(names(make_standata(y | cat(10) ~ x, family = "cumulative", 
                                   data = temp_data)), 
               c("N","Y","K","X","ncat","max_obs"))
  expect_warning(names(make_standata(y | cat(t) ~ x, family = "cumulative", 
                                     data = temp_data)),
                 "no longer have different numbers of categories")
  standata <- make_standata(y | trunc(0,20) ~ x, family = "gaussian", 
                            data = temp_data)
  expect_true(standata$lb == 0 && standata$ub == 20)
})

test_that(paste("make_standata accepts correct response variables", 
                "depending on the family"), {
  expect_equal(make_standata(y ~ 1, data = data.frame(y = seq(-9.9,0,0.1)), 
                             family = "student")$Y, seq(-9.9,0,0.1))
  expect_equal(make_standata(y ~ 1, data = data.frame(y = 1:10), 
                             family = "binomial")$Y, 1:10)
  expect_equal(make_standata(y ~ 1, data = data.frame(y = 10:20), 
                             family = "poisson")$Y, 10:20)
  expect_equal(make_standata(y ~ 1, data = data.frame(y = rep(-c(1:2),5)), 
                             family = "bernoulli")$Y, rep(1:0,5))
  expect_equal(make_standata(y ~ 1, data = data.frame(y = rep(c(TRUE, FALSE),5)),
                             family = "bernoulli")$Y, rep(1:0,5))
  expect_equal(make_standata(y ~ 1, data = data.frame(y = rep(1:10,5)), 
                             family = "categorical")$Y, rep(1:10,5))
  expect_equal(make_standata(y ~ 1, data = data.frame(y = rep(-4:5,5)), 
                             family = "categorical")$Y, rep(1:10,5))
  expect_equal(make_standata(y ~ 1, data = data.frame(y = factor(rep(-4:5,5))), 
                             family = "categorical")$Y, rep(1:10,5))
  expect_equal(make_standata(y ~ 1, data = data.frame(y = rep(1:10,5)), 
                             family = "cumulative")$Y, rep(1:10,5))
  temp_data <- data.frame(y = factor(rep(-4:5,5), order = TRUE))
  expect_equal(make_standata(y ~ 1, data = temp_data, family = "acat")$Y, 
               rep(1:10,5))
  expect_equal(make_standata(y ~ 1, data = data.frame(y = seq(0,10,0.1)), 
                             family = "exponential")$Y, seq(0,10,0.1))
  temp_data <- data.frame(y1 = 1:10, y2 = 11:20, w = 1:10, x = rep(0,10))
  expect_equal(make_standata(cbind(y1,y2) | weights(w) ~ x, family = "gaussian",
                             data = temp_data)$Y, 
               cbind(1:10,11:20))
})

test_that(paste("make_standata rejects incorrect response variables", 
                "depending on the family"), {
  expect_error(make_standata(y ~ 1, data = data.frame(y = factor(1:10)), 
                             family = "cauchy"),
               "family cauchy expects numeric response variable")
  expect_error(make_standata(y ~ 1, data = data.frame(y = -5:5), 
                             family = "geometric"),
               "family geometric expects response variable of non-negative integers")
  expect_error(make_standata(y ~ 1, data = data.frame(y = -1:1), 
                             family = "bernoulli"),
               "family bernoulli expects response variable to contain only two different values")
  expect_error(make_standata(y ~ 1, data = data.frame(y = factor(-1:1)), 
                             family = "cratio"),
               "family cratio requires factored response variables to be ordered")
  expect_error(make_standata(y ~ 1, data = data.frame(y = rep(0.5:7.5), 2), 
                             family = "sratio"),
               paste("family sratio expects either integers or ordered factors", 
                     "as response variables"))
  expect_error(make_standata(y ~ 1, data = data.frame(y = rep(-7.5:7.5), 2), 
                             family = "gamma"),
               "family gamma requires response variable to be non-negative")
})

test_that("make_standata suggests using family bernoulli if appropriate", {
  expect_message(make_standata(y ~ 1, data = data.frame(y = rep(0:1,5)), 
                               family = "binomial"),
                 paste("Only 2 levels detected so that family bernoulli", 
                       "might be a more efficient choice."))
  expect_message(make_standata(y ~ 1, data = data.frame(y = rep(0:1,5)), 
                               family = "categorical"),
                 paste("Only 2 levels detected so that family bernoulli", 
                       "might be a more efficient choice."))
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
                             family = "categorical")$max_obs, 19)
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
                             autocor = cor_ar(~tim | g:trait))$Y,
               cbind(c(seq(9,1,-2), seq(10,2,-2)), 
                     c(seq(19,11,-2), seq(20,12,-2))))
  expect_error(make_standata(cbind(y1,y2) | weights(w) ~ x, 
                             family = "gaussian", data = temp_data,
                             autocor = cor.ar(~tim | g)),
               paste("autocorrelation structures for multiple responses", 
                     "must contain 'trait' as grouping variable"))
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
                             autocor = cor_ma(~tim|g))$tgroup,
               c(rep(1,5), rep(2,5)))
  expect_equal(make_standata(y ~ x, data = temp_data,
                             autocor = cor_ar(~tim|g))$tgroup,
               c(rep(1,5), rep(2,5)))
  standata <- make_standata(y ~ x, data = temp_data,
                            autocor = cor_ar(~tim|g, cov = TRUE))
  expect_equal(standata$begin_tg, c(1, 6))
  expect_equal(standata$nrows_tg, c(5, 5))
  expect_equal(standata$squared_se, rep(0, 10))
})

test_that("make_standata allows to retrieve the initial data order", {
  temp_data <- data.frame(y1 = rnorm(100), y2 = rnorm(100), 
                          id = sample(1:10, 100, TRUE), 
                          time = sample(1:100, 100))
  # univariate model
  sdata1 <- make_standata(y1 ~ 1, data = temp_data, 
                          autocor = cor_ar(~time|id),
                          control = list(save_order = TRUE))
  expect_equal(temp_data$y1, sdata1$Y[attr(sdata1, "old_order")])
  # multivariate model
  sdata2 <- make_standata(cbind(y1, y2) ~ 1, data = temp_data, 
                          autocor = cor_ma(~time|id:trait),
                          control = list(save_order = TRUE))
  expect_equal(c(temp_data$y1, temp_data$y2), 
               sdata2$Y[attr(sdata2, "old_order")])
})

test_that("make_standata rejects invalid input for argument partial", {
  expect_error(make_standata(rating ~ 1, data = inhaler,
                             partial = ~treat, family = "gaussian"))
  expect_error(make_standata(rating ~ 1, data = inhaler,
                             partial = 1, family = "acat"))
})

test_that("make_standata handles covariance matrices correctly", {
  A <- structure(diag(1, 4), dimnames = list(1:4, NULL))
  expect_equivalent(make_standata(count ~ Trt_c + (1|visit), data = epilepsy,
                             cov.ranef = list(visit = A))$cov_1, A)
  B <- diag(1, 4)
  expect_error(make_standata(count ~ Trt_c + (1|visit), data = epilepsy,
                             cov.ranef = list(visit = B)),
               "rownames are required")
  B <- structure(diag(1, 4), dimnames = list(2:5, NULL))
  expect_error(make_standata(count ~ Trt_c + (1|visit), data = epilepsy,
                             cov.ranef = list(visit = B)),
               "rownames .* do not match")
  B <- structure(diag(1, 5), dimnames = list(1:5, NULL))
  expect_error(make_standata(count ~ Trt_c + (1|visit), data = epilepsy,
                             cov.ranef = list(visit = B)),
               "dimension .* is incorrect")
  B <- A
  B[1,2] <- 0.5
  expect_error(make_standata(count ~ Trt_c + (1|visit), data = epilepsy,
                             cov.ranef = list(visit = B)),
               "not symmetric")
})

test_that("make_standata computes data for inverse.gaussian models", {
  temp_data <- data.frame(y = 1:10, x = rep(0,10), w = 1:10)
  standata <- make_standata(y ~ x, data = temp_data, 
                            family = inverse.gaussian)
  expect_equal(standata$log_Y, sum(log(temp_data$y)))
  expect_equal(standata$sqrt_Y, sqrt(temp_data$y))
  standata <- make_standata(y | weights(w) ~ x, data = temp_data,
                            family = inverse.gaussian)
  expect_equal(standata$log_Y, log(temp_data$y))                         
})

test_that("make_standata computes data for 2PL models", {
  temp_data <- data.frame(y = sample(0:1, 10, TRUE), x = rnorm(10))
  standata <- make_standata(y ~ x, data = temp_data,
                            family = bernoulli(type = "2PL"))
  expect_equal(standata$Y, temp_data$y)
  expect_equal(standata$N_trait, 10)
  temp_data$y <- factor(rep(c("N", "Y"), each = 5), levels = c("N", "Y"))
  standata <- make_standata(y ~ x, data = temp_data,
                            family = bernoulli(type = "2PL"))
  expect_equal(standata$Y, rep(0:1, each = 5))
})

test_that("amend_newdata handles factors correctly", {
  fit <- rename_pars(brmsfit_example)
  fit$data$fac <- factor(sample(1:3, nrow(fit$data), replace = TRUE))
  newdata <- fit$data[1:5, ]
  expect_silent(amend_newdata(newdata, fit))
  newdata$visit <- 1:5
  expect_error(amend_newdata(newdata, fit), fixed = TRUE,
               "levels 5 of grouping factor visit not found")
  newdata$fac <- 1:5
  expect_error(amend_newdata(newdata, fit), fixed = TRUE,
               "New factor levels are not allowed")
})

test_that("brmdata and brm.data are backwards compatible", {
  temp_data <- data.frame(y = 1:10, x = sample(1:5, 10, TRUE))
  expect_identical(brmdata(y ~ x + (1|x), data = temp_data, 
                           family = "poisson"), 
                   make_standata(y ~ x + (1|x), data = temp_data, 
                                 family = "poisson"))
  expect_identical(brmdata(y ~ 1, data = temp_data, 
                            family = "acat", partial = ~ x), 
                   make_standata(y ~ 1, data = temp_data, 
                                 family = "acat", partial = ~ x))
  expect_identical(brm.data(y ~ x + (1|x), data = temp_data, 
                            family = "poisson"), 
                   make_standata(y ~ x + (1|x), data = temp_data, 
                                 family = "poisson"))
  expect_identical(brm.data(y ~ 1, data = temp_data, 
                            family = "acat", partial = ~ x), 
                   make_standata(y ~ 1, data = temp_data, 
                                 family = "acat", partial = ~ x))
})