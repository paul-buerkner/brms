context("Tests for exclude_pars helper functions")

test_that("exclude_pars returns expected parameter names", {
  dat <- data.frame(
    y = 1:10, x1 = rnorm(10), x2 = rnorm(10),
    g = rep(1:5, 2), h = factor(rep(1:5, each = 2))
  )

  fit <- brm(y ~ x1 * x2 + (x1 * x2 | g) + (1 | h), dat,
    empty = TRUE
  )
  ep <- brms:::exclude_pars(fit)
  expect_true(all(c("r_1", "r_2") %in% ep))

  fit <- brm(y ~ x1 * x2 + (x1 * x2 | g) + (1 | h), dat,
    empty = TRUE, save_pars = save_pars(all = TRUE)
  )
  ep <- brms:::exclude_pars(fit)
  expect_true(!any(c("z_1", "z_2") %in% ep))

  fit <- brm(y ~ x1 * x2 + (x1 * x2 | g) + (1 | h), dat,
    empty = TRUE, save_pars = save_pars(group = FALSE)
  )
  ep <- brms:::exclude_pars(fit)
  expect_true("r_1_1" %in% ep)

  fit <- brm(y ~ x1 * x2 + (x1 | g) + (1 | h), dat,
    empty = TRUE, save_pars = save_pars(group = "h")
  )
  ep <- brms:::exclude_pars(fit)
  expect_true(!"r_1_3" %in% ep)

  fit <- brm(y ~ s(x1) + x2, dat, empty = TRUE)
  ep <- brms:::exclude_pars(fit)
  expect_true("zs_1_1" %in% ep)

  fit <- brm(bf(y ~ eta, eta ~ x1 + s(x2), nl = TRUE), dat, empty = TRUE)
  ep <- brms:::exclude_pars(fit)
  expect_true("zs_eta_1_1" %in% ep)

  fit <- brm(y ~ me(x1, g), dat, empty = TRUE)
  ep <- brms:::exclude_pars(fit)
  expect_true("Xme_1" %in% ep)

  fit <- brm(y ~ me(x1, g), dat,
    empty = TRUE,
    save_pars = save_pars(latent = "x1")
  )
  ep <- brms:::exclude_pars(fit)
  expect_true(!"Xme_1" %in% ep)

  fit <- brm(y ~ me(x1, g), dat,
    empty = TRUE,
    save_pars = save_pars(manual = "Lme_1")
  )
  ep <- brms:::exclude_pars(fit)
  expect_true(!"Lme_1" %in% ep)
})
