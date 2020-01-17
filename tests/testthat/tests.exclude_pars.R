test_that("exclude_pars returns expected parameter names", {
  ranef <- data.frame(
    id = c(1, 1, 2), group = c("g1", "g1", "g2"),
    gn = c(1, 1, 2), coef = c("x", "z", "x"), 
    cn = c(1, 2, 1), nlpar = "", ggn = c(1, 2, 2), 
    cor = c(TRUE, TRUE, FALSE)
  )
  fit <- brms:::brmsfit()
  fit$ranef <- ranef
  fit$formula <- y ~ 1
  # empty_effects <- structure(list(), class = "brmsterms")
  ep <- brms:::exclude_pars(fit)
  expect_true(all(c("r_1", "r_2") %in% ep))
  ep <- brms:::exclude_pars(fit, ranef = ranef, save_ranef = FALSE)
  expect_true("r_1_1" %in% ep)
  
  fit$ranef$nlpar <- c("a", "a", "")
  ep <- brms:::exclude_pars(fit, save_ranef = FALSE)
  expect_true(all(c("r_1_a_1", "r_1_a_2") %in% ep))
  
  fit$formula <- y ~ x + s(z)
  fit$data <- data.frame(y = rnorm(20), x = rnorm(20), z = rnorm(20))
  expect_true("zs_1_1" %in% brms:::exclude_pars(fit))
  fit$formula <- bf(y ~ eta, eta ~ x + s(z), family = gaussian(), nl = TRUE)
  expect_true("zs_eta_1_1" %in% brms:::exclude_pars(fit))
})
