test_that("(deprecated) melt_data returns data in long format", {
  data <- data.frame(x = rep(c("a","b"), 5), y1 = 1:10, y2 = 11:20, 
                     y3 = 21:30, z = 100:91)
  bterms <- parse_bf(y ~ x, family = "poisson")
  expect_equivalent(melt_data(data, bterms = bterms, 
                         family = "poisson"), data)
  
  target1 <- data.frame(x = rep(c("a","b"), 10), y1 = rep(1:10, 2), 
                        y2 = rep(11:20, 2), y3 = rep(21:30, 2),
                        z = rep(100:91, 2), 
                        trait = factor(rep(c("y3", "y1"), each = 10),
                                       levels = c("y3", "y1")), 
                        response = c(21:30, 1:10))
  bterms <- parse_bf(cbind(y3,y1) ~ x, family = "gaussian")
  expect_equivalent(melt_data(data, bterms = bterms, family = "gaussian"), 
               target1[, c("x", "y1", "y3", "trait", "response")])
  
  target2 <- data.frame(x = rep(c("a","b"), 15), y1 = rep(1:10, 3), 
                        y2 = rep(11:20, 3), y3 = rep(21:30, 3),
                        z = rep(100:91, 3),
                        trait = factor(rep(c("y2", "y1", "y3"), each = 10),
                                       levels = c("y2", "y1", "y3")), 
                        response = c(11:20, 1:10, 21:30))
  bterms <- parse_bf(cbind(y2,y1,y3) ~ x, family = "gaussian")
  expect_equivalent(melt_data(data, bterms = bterms, family = "gaussian"), 
               target2[, c("x", "y1", "y2", "y3", "trait", "response")])
})

test_that("(deprecated) melt_data keeps factor contrasts", {
  data <- data.frame(y1 = rnorm(10), y2 = rnorm(10),
                     x = factor(rep(1:2, each = 5)))
  contrasts(data$x) <- contr.sum(2)
  bterms <- parse_bf(cbind(y1,y2) ~ x)
  newdata <- melt_data(data, family = "gaussian", bterms = bterms)
  expect_equal(attr(newdata$x, "contrasts"), attr(data$x, "contrasts"))
})

test_that("(deprecated) melt_data returns expected errors", {
  data <- data.frame(y1 = rnorm(10), y2 = rnorm(10), x = 1:10)
  
  formula <- bf(y1 ~ x:main, family = hurdle_poisson())
  formula$old_mv <- TRUE
  bterms <- parse_bf(formula)
  expect_error(melt_data(data = NULL, family = hurdle_poisson(), bterms = bterms),
               "'data' must be a data.frame", fixed = TRUE)
  
  data$main <- 1:10 
  expect_error(melt_data(data = data, family = hurdle_poisson(), bterms = bterms),
               "'main' is a reserved variable name", fixed = TRUE)
  
  data$response <- 1:10
  formula <- bf(response ~ x:main, family = hurdle_poisson())
  formula$old_mv <- TRUE
  bterms <- parse_bf(formula)
  expect_error(melt_data(data = data, family = hurdle_poisson(), bterms = bterms),
               "'response' is a reserved variable name", fixed = TRUE)
  
  data$trait <- 1:10
  formula <- bf(y ~ 0 + x*trait, family = hurdle_poisson())
  formula$old_mv <- TRUE
  bterms <- parse_bf(formula)
  expect_error(melt_data(data = data, family = hurdle_poisson(), bterms = bterms),
               "'trait', 'response' is a reserved variable name", fixed = TRUE)
  
  data <- data.frame(y1 = rnorm(10), y2 = rnorm(10), x = 1:10)
  formula <- bf(cbind(y1, y2) ~ x)
  formula$old_mv <- TRUE
  bterms <- parse_bf(formula)
  expect_error(melt_data(data = data, family = poisson(), bterms = bterms),
               "Invalid multivariate model", fixed = TRUE)
})

test_that("arr_design_matrix works correctly", {
  expect_equal(arr_design_matrix(1:10, 0, sort(rep(1:2, 5))), NULL)
  expect_equal(arr_design_matrix(1:10, 1, sort(rep(1:2, 5))), 
               matrix(c(0,1:4.5,0,6:9.5)))
  expect_equal(arr_design_matrix(1:10, 2, sort(rep(1:2, 5))), 
               cbind(c(0, 1:4.5, 0, 6:9), c(0, 0, 1:3, 0 ,0, 6:8)))
})

test_that("amend_newdata handles factors correctly", {
  fit <- brms:::rename_pars(brms:::brmsfit_example1)
  fit$data$fac <- factor(sample(1:3, nrow(fit$data), TRUE))
  newdata <- fit$data[1:5, ]
  expect_silent(brms:::amend_newdata(newdata, fit))
  newdata$visit <- 1:5
  expect_error(brms:::amend_newdata(newdata, fit),
               "Levels '5' of grouping factor 'visit' cannot")
  newdata$fac <- 1:5
  expect_error(brms:::amend_newdata(newdata, fit),
               "New factor levels are not allowed")
})

test_that("update_data returns correct model.frames", {
  dat <- data.frame(y = 1:5, x = 1:5, z = 6:10, g = 5:1)
  
  bterms <- brms:::parse_bf(y ~ as.numeric(x) + (as.factor(z) | g))
  mf <- brms:::update_data(dat, family = gaussian(), bterms = bterms)
  expect_true(all(c("x", "z") %in% names(mf)))
  
  bterms <- brms:::parse_bf(y ~ 1 + (1|g/x/z))
  mf <- brms:::update_data(dat, family = gaussian(), bterms = bterms)
  expect_equal(mf[["g:x"]], paste0(dat$g, "_", dat$x))
  expect_equal(mf[["g:x:z"]], paste0(dat$g, "_", dat$x, "_", dat$z))
})

test_that("(deprecated) update_data handles NAs correctly in old MV models", {
  data <- data.frame(y1 = c(1, NA, 3), y2 = 4:6, x = 10:12, z = NA)
  formula <- bf(cbind(y1, y2) ~ x)
  formula$old_mv <- TRUE
  bterms <- parse_bf(formula, family = "gaussian")
  expect_warning(mf <- update_data(data, family = "gaussian", bterms = bterms),
                 "NAs were excluded")
  expect_equivalent(mf, data.frame(response = c(1, 3, 4, 6), y1 = c(1, 3, 1, 3), 
                                   y2 = c(4, 6, 4, 6), x = c(10, 12, 10, 12)))
  
  formula <- bf(y1 ~ x, family = "hurdle_gamma")
  formula$old_mv <- TRUE
  bterms <- parse_bf(formula)
  expect_warning(mf <- update_data(data, family = "hurdle_gamma", 
                                   bterms = bterms),
                 "NAs were excluded")
  expect_equivalent(mf, data.frame(response = c(1, 3, 1, 3), 
                                   y1 = c(1, 3, 1, 3),
                                   x = c(10, 12, 10, 12)))
  
  formula$family <- zero_inflated_poisson()
  bterms <- parse_bf(formula)
  expect_warning(mf <- update_data(data, family = "zero_inflated_poisson", 
                                   bterms = bterms), "NAs were excluded")
  expect_equivalent(mf, data.frame(response = c(1, 3, 1, 3), y1 = c(1, 3, 1, 3),
                                   x = c(10, 12, 10, 12)))
})
