test_that("Test that melt returns data in correct long format", {
  data <- data.frame(x = rep(c("a","b"), 5), y1 = 1:10, y2 = 11:20, y3 = 21:30, z = 100:91)
  expect_equal(melt(data, response = "y1", family = "poisson"), data)
  
  target1 <- data.frame(x = rep(c("a","b"), 10), y2 = rep(11:20, 2), z = rep(100:91, 2),
                        trait = c(rep("y3", 10), rep("y1", 10)), y3 = c(21:30,1:10))
  expect_equal(melt(data, response = c("y3", "y1"), family = "gaussian"), target1)
  
  target2 <- data.frame(x = rep(c("a","b"), 15), z = rep(100:91, 3),
                        trait = c(rep("y2", 10), rep("y1", 10), rep("y3", 10)), 
                        y2 = c(11:20, 1:10, 21:30))
  expect_equal(melt(data, response = c("y2", "y1", "y3"), family = "gaussian"), target2)
})

test_that("Test that combine_groups does the expected", {
  data <- data.frame(x = rep(c("a","b"), 5), y1 = 1:10, y2 = 11:20, y3 = 21:30, z = 100:91)
  expected <- data 
  expected[["y1:y2"]] <- paste0(data$y1, "_", data$y2)
  expected[["y1:y2:y3"]] <- paste0(data$y1, "_", data$y2, "_", data$y3)
  expect_equal(combine_groups(data, "y1:y2", "y1:y2:y3"), expected)
})

test_that("Test that get_model_matrix removes intercepts correctly", {
  data <- data.frame(x = factor(rep(1:2, 5)), y = 11:20)
  expect_equal(get_model_matrix(y ~ x, data, rm_intercept = TRUE),
               structure(matrix(rep(0:1, 5)), dimnames = list(1:10, "x2")))
})

test_that("Test that arr_design_matrix returns correct design matrices for autoregressive effects", {
  expect_equal(arr_design_matrix(1:10, 0, sort(rep(1:2, 5))), NULL)
  expect_equal(arr_design_matrix(1:10, 1, sort(rep(1:2, 5))), 
               matrix(c(0,1:4.5,0,6:9.5)))
  expect_equal(arr_design_matrix(1:10, 2, sort(rep(1:2, 5))), 
               cbind(c(0,1:4.5,0,6:9), c(0,0,1:3,0,0,6:8)))
})

test_that("Test that generate_standata returns correct data names for fixed and random effects", {
  expect_equal(names(generate_standata(rating ~ treat + period + carry 
                                       + (1|subject), data = inhaler)),
               c("N","Y","K","X","J_1","N_1","K_1","Z_1","NC_1"))
  expect_equal(names(generate_standata(rating ~ treat + period + carry 
                                       + (1+treat|subject), data = inhaler,
                                       family = "categorical")),
               c("N","Y","Kp","Xp","J_1","N_1","K_1",
                 "Z_1","NC_1", "ncat", "max_obs"))
  temp_data <- data.frame(y = 1:10, g = 1:10, h = 11:10, x = rep(0,10))
  expect_equal(names(generate_standata(y ~ x + (1|g) + (1|h), 
                                       family = "poisson",
                                       data = temp_data)),
               c("N","Y","K","X","J_1","N_1","K_1","Z_1","NC_1",
                 "J_2","N_2","K_2","Z_2","NC_2"))
})

test_that(paste0("Test that generate_standata handles variables used as fixed effects", 
                 "and grouping factors at the same time"), {
  data <- data.frame(y = 1:9, x = factor(rep(c("a","b","c"), 3)))
  standata <- generate_standata(y ~ x + (1|x), data = data)
  expect_equal(colnames(standata$X), c("xb", "xc"))
  expect_equal(standata$J_1, rep(1:3, 3))
  standata2 <- generate_standata(y ~ x + (1|x), data = data, 
                                 keep_intercept = TRUE)
  expect_equal(colnames(standata2$X), c("Intercept", "xb", "xc"))
})

test_that("Test that generate_standata returns correct data names for addition and partial variables", {
  data <- data.frame(y = 1:10, w = 1:10, t = 1:10, x = rep(0,10), c = sample(-1:1,10,TRUE))
  expect_equal(names(generate_standata(y | se(w) ~ x, family = "gaussian", 
                                       data = data)), 
               c("N","Y","K","X","se"))
  expect_equal(names(generate_standata(y | weights(w) ~ x, family = "gaussian", 
                                       data = data)), 
               c("N","Y","K","X","weights"))
  expect_equal(names(generate_standata(y | cens(c) ~ x, family = "cauchy", 
                                       data = data)), 
               c("N","Y","K","X","cens"))
  expect_equal(names(generate_standata(y | trials(t) ~ x, family = "binomial", 
                                       data = data)), 
               c("N","Y","K","X","trials","max_obs"))
  expect_equal(names(generate_standata(y | trials(10) ~ x, family = "binomial", 
                                       data = data)), 
               c("N","Y","K","X","trials","max_obs"))
  expect_equal(names(generate_standata(y | cat(11) ~ x, family = "acat", 
                                       data = data)), 
               c("N","Y","K","X","ncat","max_obs"))
  expect_equal(names(generate_standata(y | cat(10) ~ x, family = "cumulative", 
                                       data = data)), 
               c("N","Y","K","X","ncat","max_obs"))
  temp_data <- generate_standata(y | trunc(0,20) ~ x, family = "gaussian", 
                                 data = data)
  expect_true(temp_data$lb == 0 && temp_data$ub == 20)
})

test_that("Test that generate_standata accepts correct response variables depending on the family", {
  expect_equal(generate_standata(y ~ 1, data = data.frame(y = seq(-9.9,0,0.1)), 
                                 family = "student")$Y, seq(-9.9,0,0.1))
  expect_equal(generate_standata(y ~ 1, data = data.frame(y = 1:10), 
                                 family = "binomial")$Y, 1:10)
  expect_equal(generate_standata(y ~ 1, data = data.frame(y = 10:20), 
                                 family = "poisson")$Y, 10:20)
  expect_equal(generate_standata(y ~ 1, data = data.frame(y = rep(-c(1:2),5)), 
                                 family = "bernoulli")$Y, rep(1:0,5))
  expect_equal(generate_standata(y ~ 1, data = data.frame(y = rep(c(TRUE, FALSE),5)),
                                 family = "bernoulli")$Y, rep(1:0,5))
  expect_equal(generate_standata(y ~ 1, data = data.frame(y = rep(1:10,5)), 
                                 family = "categorical")$Y, rep(1:10,5))
  expect_equal(generate_standata(y ~ 1, data = data.frame(y = rep(-4:5,5)), 
                                 family = "categorical")$Y, rep(1:10,5))
  expect_equal(generate_standata(y ~ 1, data = data.frame(y = factor(rep(-4:5,5))), 
                                 family = "categorical")$Y, rep(1:10,5))
  expect_equal(generate_standata(y ~ 1, data = data.frame(y = rep(1:10,5)), 
                                 family = "cumulative")$Y, rep(1:10,5))
  temp_data <- data.frame(y = factor(rep(-4:5,5), order = TRUE))
  expect_equal(generate_standata(y ~ 1, data = temp_data, family = "acat")$Y, 
               rep(1:10,5))
  expect_equal(generate_standata(y ~ 1, data = data.frame(y = seq(0,10,0.1)), 
                                 family = "exponential")$Y, seq(0,10,0.1))
  expect_equal(generate_standata(cbind(y1,y2) | weights(w) ~ x, family = "gaussian",
               data = data.frame(y1 = 1:10, y2 = 11:20, w = 1:10, x = rep(0,10)))$Y, 
               cbind(1:10,11:20))
})

test_that(paste("Test that generate_standata rejects incorrect", 
                "response variables depending on the family"), {
  expect_error(generate_standata(y ~ 1, data = data.frame(y = factor(1:10)), 
                                 family = "cauchy"),
               "family cauchy expects numeric response variable")
  expect_error(generate_standata(y ~ 1, data = data.frame(y = -5:5), 
                                 family = "geometric"),
               "family geometric expects response variable of non-negative integers")
  expect_error(generate_standata(y ~ 1, data = data.frame(y = -1:1), 
                                 family = "bernoulli"),
               "family bernoulli expects response variable to contain only two different values")
  expect_error(generate_standata(y ~ 1, data = data.frame(y = factor(-1:1)), 
                                 family = "cratio"),
               "family cratio requires factored response variables to be ordered")
  expect_error(generate_standata(y ~ 1, data = data.frame(y = rep(0.5:7.5), 2), 
                                 family = "sratio"),
               "family sratio expects either integers or ordered factors as response variables")
  expect_error(generate_standata(y ~ 1, data = data.frame(y = rep(-7.5:7.5), 2), 
                       family = "gamma"),
               "family gamma requires response variable to be non-negative")
})

test_that("Test that generate_standata suggests using family bernoulli if appropriate", {
  expect_message(generate_standata(y ~ 1, data = data.frame(y = rep(0:1,5)), 
                                   family = "binomial"),
                 paste("Only 2 levels detected so that family 'bernoulli'", 
                       "might be a more efficient choice."))
  expect_message(generate_standata(y ~ 1, data = data.frame(y = rep(0:1,5)), 
                                   family = "categorical"),
                 paste("Only 2 levels detected so that family 'bernoulli'", 
                       "might be a more efficient choice."))
})

test_that("Test that generate_standata returns correct values for addition arguments", {
  temp_data <- data.frame(y = rnorm(9), s = 1:9, w = 1:9, c1 = rep(-1:1, 3), 
                          c2 = rep(c("left","none","right"), 3),
                          c3 = c(rep(c(TRUE, FALSE), 4), FALSE),
                          t = 11:19)
  expect_equal(generate_standata(y | se(s) ~ 1, data = temp_data)$se, 
               1:9)
  expect_equal(generate_standata(y | weights(w) ~ 1, data = temp_data)$weights, 
               1:9)
  expect_equal(generate_standata(y | cens(c1) ~ 1, data = temp_data)$cens, 
               rep(-1:1, 3))
  expect_equal(generate_standata(y | cens(c2) ~ 1, data = temp_data)$cens,
               rep(-1:1, 3))
  expect_equal(generate_standata(y | cens(c3) ~ 1, data = temp_data)$cens, 
               c(rep(1:0, 4), 0))
  expect_equal(generate_standata(s ~ 1, data = temp_data, 
                                 family = "binomial")$max_obs, 9)
  expect_equal(generate_standata(s | trials(10) ~ 1, data = temp_data, 
                                 family = "binomial")$max_obs, 10)
  expect_equal(generate_standata(s | trials(t) ~ 1, data = temp_data, 
                        family = "binomial")$max_obs, 11:19)
  expect_equal(generate_standata(s | cat(19) ~ 1, data = temp_data, 
                        family = "categorical")$max_obs, 19)
})

test_that("Test that generate_standata rejects incorrect addition arguments", {
  temp_data <- data.frame(y = rnorm(9), s = -(1:9), w = -(1:9), 
                          c = rep(-2:0, 3), t = 9:1, z = 1:9)
  expect_error(generate_standata(y | se(s) ~ 1, data = temp_data), 
               "standard errors must be non-negative")
  expect_error(generate_standata(y | weights(w) ~ 1, data = temp_data), 
               "weights must be non-negative")
  expect_error(generate_standata(y | cens(c) ~ 1, data = temp_data))
  expect_error(generate_standata(z | trials(t) ~ 1, data = temp_data, 
                                 family = "binomial"),
               "Number of trials is smaller than the response variable would suggest.")
})

test_that(paste("Test that generate_standata handles addition arguments", 
                "and autocorrelation in multinormal models"), {
  temp_data <- data.frame(y1 = 1:10, y2 = 11:20, w = 1:10, x = rep(0,10), 
                          tim = 10:1, g = rep(1:2,5))
  expect_equal(generate_standata(cbind(y1,y2) | weights(w) ~ x, 
                                 family = "gaussian", data = temp_data)$weights, 
               1:10)
  expect_equal(generate_standata(cbind(y1,y2) | weights(w) ~ x, 
                                 family = "gaussian", data = temp_data,
                                 autocor = cor_ar(~tim | g:trait))$Y,
               cbind(c(seq(9,1,-2), seq(10,2,-2)), c(seq(19,11,-2), seq(20,12,-2))))
  expect_error(generate_standata(cbind(y1,y2) | weights(w) ~ x, 
                                 family = "gaussian", data = temp_data,
                                 autocor = cor.ar(~tim | g)),
               paste("autocorrelation structures for multiple responses", 
                     "must contain 'trait' as grouping variable"))
})

test_that(paste("Test that generate_standata returns correct data", 
                "for autocorrelations structures"), {
  temp_data <- data.frame(y=1:10, x=rep(0,10), tim=10:1, g = rep(3:4,5))
  expect_equal(generate_standata(y ~ x, family = "gaussian", data = temp_data,
                                 autocor = cor_arr(~tim|g))$Yarr,
               cbind(c(0,9,7,5,3,0,10,8,6,4)))
  expect_equal(generate_standata(y ~ x, family = "gaussian", data = temp_data,
                                 autocor = cor_arr(~tim|g, r = 2))$Yarr,
               cbind(c(0,9,7,5,3,0,10,8,6,4), c(0,0,9,7,5,0,0,10,8,6)))
  expect_equal(generate_standata(y ~ x, family = "gaussian", data = temp_data,
                                 autocor = cor_ma(~tim|g))$tgroup,
               c(rep(1,5), rep(2,5)))
  expect_equal(generate_standata(y ~ x, family = "gaussian", data = temp_data,
                                 autocor = cor_ar(~tim|g))$tgroup,
               c(rep(1,5), rep(2,5)))
})

test_that("Test brmdata and brm.data for backwards compatibility", {
  temp_data <- data.frame(y = 1:10, x = sample(1:5, 10, TRUE))
  expect_identical(brmdata(y ~ x + (1|x), data = temp_data, 
                            family = "poisson"), 
                   generate_standata(y ~ x + (1|x), data = temp_data, 
                                     family = "poisson"))
  expect_identical(brmdata(y ~ 1, data = temp_data, 
                            family = "acat", partial = ~ x), 
                   generate_standata(y ~ 1, data = temp_data, 
                                     family = "acat", partial = ~ x))
  expect_identical(brm.data(y ~ x + (1|x), data = temp_data, 
                            family = "poisson"), 
                   generate_standata(y ~ x + (1|x), data = temp_data, 
                                     family = "poisson"))
  expect_identical(brm.data(y ~ 1, data = temp_data, 
                            family = "acat", partial = ~ x), 
                   generate_standata(y ~ 1, data = temp_data, 
                                     family = "acat", partial = ~ x))
})