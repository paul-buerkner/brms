test_that("melt_data returns data in long format", {
  data <- data.frame(x = rep(c("a","b"), 5), y1 = 1:10, y2 = 11:20, 
                     y3 = 21:30, z = 100:91)
  effects <- extract_effects(y ~ x, family = "poisson")
  expect_equivalent(melt_data(data, effects = effects, 
                         family = "poisson"), data)
  
  target1 <- data.frame(x = rep(c("a","b"), 10), y1 = rep(1:10, 2), 
                        y2 = rep(11:20, 2), y3 = rep(21:30, 2),
                        z = rep(100:91, 2), 
                        trait = factor(rep(c("y3", "y1"), each = 10),
                                       levels = c("y3", "y1")), 
                        response = c(21:30, 1:10))
  effects <- extract_effects(cbind(y3,y1) ~ x, family = "gaussian")
  expect_equivalent(melt_data(data, effects = effects, family = "gaussian"), 
               target1[, c("x", "y1", "y3", "trait", "response")])
  
  target2 <- data.frame(x = rep(c("a","b"), 15), y1 = rep(1:10, 3), 
                        y2 = rep(11:20, 3), y3 = rep(21:30, 3),
                        z = rep(100:91, 3),
                        trait = factor(rep(c("y2", "y1", "y3"), each = 10),
                                       levels = c("y2", "y1", "y3")), 
                        response = c(11:20, 1:10, 21:30))
  effects <- extract_effects(cbind(y2,y1,y3) ~ x, family = "gaussian")
  expect_equivalent(melt_data(data, effects = effects, family = "gaussian"), 
               target2[, c("x", "y1", "y2", "y3", "trait", "response")])
})

test_that("melt_data keeps factor contrasts", {
  data <- data.frame(y1 = rnorm(10), y2 = rnorm(10),
                     x = factor(rep(1:2, each = 5)))
  contrasts(data$x) <- contr.sum(2)
  effects <- extract_effects(cbind(y1,y2) ~ x)
  newdata <- melt_data(data, family = "gaussian", effects = effects)
  expect_equal(attr(newdata$x, "contrasts"), attr(data$x, "contrasts"))
})

test_that("melt_data (deprecated) returns expected errors", {
  ee <- extract_effects(structure(y1 ~ x:main, old_mv = TRUE), 
                        family = hurdle_poisson())
  data <- data.frame(y1 = rnorm(10), y2 = rnorm(10), x = 1:10)
  expect_error(melt_data(data = NULL, family = hurdle_poisson(), effects = ee),
               "'data' must be a data.frame", fixed = TRUE)
  data$main <- 1:10 
  expect_error(melt_data(data = data, family = hurdle_poisson(), effects = ee),
               "'main' is a reserved variable name", fixed = TRUE)
  data$response <- 1:10 
  ee <- extract_effects(structure(response ~ x:main, old_mv = TRUE), 
                        family = hurdle_poisson())
  expect_error(melt_data(data = data, family = hurdle_poisson(), effects = ee),
               "'response' is a reserved variable name", fixed = TRUE)
  data$trait <- 1:10 
  ee <- extract_effects(structure(y ~ 0 + x*trait, old_mv = TRUE), 
                        family = hurdle_poisson())
  expect_error(melt_data(data = data, family = hurdle_poisson(), effects = ee),
               "'trait', 'response' is a reserved variable name", fixed = TRUE)
  
  ee <- extract_effects(structure(cbind(y1, y2) ~ x, old_mv = TRUE))
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
  expect_equal(get_model_matrix(y ~ x, data, cols2remove = "(Intercept)"),
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


test_that("amend_newdata handles factors correctly", {
  fit <- brms:::rename_pars(brmsfit_example1)
  fit$data$fac <- factor(sample(1:3, nrow(fit$data), replace = TRUE))
  newdata <- fit$data[1:5, ]
  expect_silent(amend_newdata(newdata, fit))
  newdata$visit <- 1:5
  expect_error(amend_newdata(newdata, fit), fixed = TRUE,
               "Levels '5' of grouping factor 'visit' cannot")
  newdata$fac <- 1:5
  expect_error(amend_newdata(newdata, fit), fixed = TRUE,
               "New factor levels are not allowed")
})

test_that("update_data handles NAs correctly in old MV models", {
  data <- data.frame(y1 = c(1, NA, 3), y2 = 4:6, x = 10:12, z = NA)
  ee <- extract_effects(structure(cbind(y1, y2) ~ x, old_mv = TRUE), 
                        family = "gaussian")
  expect_warning(mf <- update_data(data, family = "gaussian", effects = ee),
                 "NAs were excluded")
  expect_equivalent(mf, data.frame(response = c(1, 3, 4, 6), y1 = c(1, 3, 1, 3), 
                                   y2 = c(4, 6, 4, 6), x = c(10, 12, 10, 12)))
  
  ee <- extract_effects(structure(y1 ~ x, old_mv = TRUE), 
                        family = "hurdle_gamma")
  expect_warning(mf <- update_data(data, family = "hurdle_gamma", 
                                   effects = ee),
                 "NAs were excluded")
  expect_equivalent(mf, data.frame(response = c(1, 3, 1, 3), 
                                   y1 = c(1, 3, 1, 3),
                                   x = c(10, 12, 10, 12)))
  
  ee <- extract_effects(structure(y1 ~ x, old_mv = TRUE), 
                        family = "zero_inflated_poisson")
  expect_warning(mf <- update_data(data, family = "zero_inflated_poisson", 
                                   effects = ee), "NAs were excluded")
  expect_equivalent(mf, data.frame(response = c(1, 3, 1, 3), y1 = c(1, 3, 1, 3),
                                   x = c(10, 12, 10, 12)))
})