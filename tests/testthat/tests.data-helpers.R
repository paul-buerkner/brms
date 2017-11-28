test_that("arr_design_matrix works correctly", {
  expect_equal(arr_design_matrix(1:10, 0, sort(rep(1:2, 5))), NULL)
  expect_equal(arr_design_matrix(1:10, 1, sort(rep(1:2, 5))), 
               matrix(c(0,1:4.5,0,6:9.5)))
  expect_equal(arr_design_matrix(1:10, 2, sort(rep(1:2, 5))), 
               cbind(c(0, 1:4.5, 0, 6:9), c(0, 0, 1:3, 0 ,0, 6:8)))
})

test_that("validate_newdata handles factors correctly", {
  fit <- brms:::rename_pars(brms:::brmsfit_example1)
  fit$data$fac <- factor(sample(1:3, nrow(fit$data), TRUE))
  newdata <- fit$data[1:5, ]
  expect_silent(brms:::validate_newdata(newdata, fit))
  newdata$visit <- 1:5
  expect_error(brms:::validate_newdata(newdata, fit),
               "Levels '5' of grouping factor 'visit' cannot")
  newdata$fac <- 1:5
  expect_error(brms:::validate_newdata(newdata, fit),
               "New factor levels are not allowed")
})

test_that("update_data returns correct model.frames", {
  dat <- data.frame(y = 1:5, x = 1:5, z = 6:10, g = 5:1)
  
  bterms <- parse_bf(y ~ as.numeric(x) + (as.factor(z) | g),
                     family = gaussian())
  mf <- brms:::update_data(dat, bterms = bterms)
  expect_true(all(c("x", "z") %in% names(mf)))
  
  bterms <- parse_bf(y ~ 1 + (1|g/x/z), family = gaussian())
  mf <- brms:::update_data(dat, bterms = bterms)
  expect_equal(mf[["g:x"]], paste0(dat$g, "_", dat$x))
  expect_equal(mf[["g:x:z"]], paste0(dat$g, "_", dat$x, "_", dat$z))
})

