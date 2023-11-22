context("Tests for restructuring of old brmsfit objects")

test_that("restructure can be run without error", {
  # This test does not check if old models can really be restructured
  # since restructure is called with an already up-to-date model.
  fit2 <- brms:::rename_pars(brms:::brmsfit_example2)
  fit2$version <- NULL
  fit2$exclude <- c("L_1", "zs_1")
  expect_warning(
    fit2_up <- restructure(fit2),
    "Models fitted with brms < 1.0 are no longer offically supported"
  )
  expect_is(fit2_up, "brmsfit")
})

test_that("restructure_formula_v1 works correctly", {
  form <- structure(
    y ~ x + z, sigma = sigma ~ x,
    class = c("brmsformula", "formula")
  )
  form <- brms:::restructure_formula_v1(form)
  expect_equal(form$formula, y ~ x + z)
  expect_equal(form$pforms, list(sigma = sigma ~ x))
  expect_true(!attr(form$formula, "nl"))

  form <- structure(
    y ~ a * exp(-b * x),
    nonlinear = list(a = a ~ x, b = b ~ 1),
    class = c("brmsformula", "formula")
  )
  form <- brms:::restructure_formula_v1(form)
  expect_equal(form$formula, y ~ a * exp(-b * x))
  expect_equal(form$pforms, list(a = a ~ x, b = b ~ 1))
  expect_true(attr(form$formula, "nl"))
})

test_that("rename_prior returns expected lists", {
  pars <- c("b", "b_1", "bp", "bp_1", "prior_b", "prior_b__1",
            "prior_b__3", "sd_x[1]", "prior_bp__1")
  expect_equivalent(
    brms:::rename_prior(
      class = "b", pars = pars, names = c("x1", "x3", "x2")
    ),
    list(list(pos = 6, fnames = "prior_b_x1"),
         list(pos = 7, fnames = "prior_b_x2"))
  )
  expect_equivalent(
    brms:::rename_prior(
      class = "bp", pars = pars,
      names = c("x1", "x2"), new_class = "b"
    ),
    list(list(pos = 9, fnames = "prior_b_x1")))
})

test_that("rename_old_re and rename_old_re2 return expected lists", {
  data <- data.frame(y = rnorm(10), x = rnorm(10), g = 1:10)
  bterms <- brmsterms(bf(y ~ a, a ~ x + (1+x|g),
                        family = gaussian(), nl = TRUE))
  ranef <- brms:::tidy_ranef(bterms, data = data)
  target <- list(
    list(pos = c(rep(FALSE, 2), TRUE, rep(FALSE, 22)),
         oldname = "sd_a_g_Intercept", pnames = "sd_g_a_Intercept",
         fnames = "sd_g_a_Intercept", dims = numeric(0)),
    list(pos = c(rep(FALSE, 3), TRUE, rep(FALSE, 21)),
         oldname = "sd_a_g_x", pnames = "sd_g_a_x",
         fnames = "sd_g_a_x", dims = numeric(0)),
    list(pos = c(rep(FALSE, 4), TRUE, rep(FALSE, 20)),
         oldname = "cor_a_g_Intercept_x", pnames = "cor_g_a_Intercept_a_x",
         fnames = "cor_g_a_Intercept_a_x", dims = numeric(0)),
    list(pos = c(rep(FALSE, 5), rep(TRUE, 20)), oldname = "r_a_g",
         pnames = "r_g_a",
         fnames = c(paste0("r_g_a[", 1:10, ",Intercept]"),
                    paste0("r_g_a[", 1:10, ",x]")),
         dims = c(10, 2)))

  pars <- c("b_a_Intercept", "b_a_x", "sd_a_g_Intercept", "sd_a_g_x",
            "cor_a_g_Intercept_x", paste0("r_a_g[", 1:10, ",Intercept]"),
            paste0("r_a_g[", 1:10, ",x]"))
  dims <- list("sd_a_g_Intercept" = numeric(0), "sd_a_g_x" = numeric(0),
               "cor_a_g_Intercept_x" = numeric(0), "r_a_g" = c(10, 2))
  expect_equivalent(brms:::rename_old_re(ranef, pars = pars, dims = dims), target)

  target <- list(
    list(pos = c(rep(FALSE, 2), TRUE, rep(FALSE, 22)),
         oldname = "sd_g_a_Intercept", pnames = "sd_g__a_Intercept",
         fnames = "sd_g__a_Intercept", dims = numeric(0)),
    list(pos = c(rep(FALSE, 3), TRUE, rep(FALSE, 21)),
         oldname = "sd_g_a_x", pnames = "sd_g__a_x",
         fnames = "sd_g__a_x", dims = numeric(0)),
    list(pos = c(rep(FALSE, 4), TRUE, rep(FALSE, 20)),
         oldname = "cor_g_a_Intercept_a_x", pnames = "cor_g__a_Intercept__a_x",
         fnames = "cor_g__a_Intercept__a_x", dims = numeric(0)),
    list(pos = c(rep(FALSE, 5), rep(TRUE, 20)), oldname = "r_g_a",
         pnames = "r_g__a",
         fnames = c(paste0("r_g__a[", 1:10, ",Intercept]"),
                    paste0("r_g__a[", 1:10, ",x]")),
         dims = c(10, 2)))

  pars <- c("b_a_Intercept", "b_a_x", "sd_g_a_Intercept", "sd_g_a_x",
            "cor_g_a_Intercept_a_x", paste0("r_g_a[", 1:10, ",Intercept]"),
            paste0("r_g_a[", 1:10, ",x]"))
  dims <- list("sd_g_a_Intercept" = numeric(0), "sd_g_a_x" = numeric(0),
               "cor_g_a_Intercept_a_x" = numeric(0), "r_g_a" = c(10, 2))
  expect_equivalent(brms:::rename_old_re2(ranef, pars = pars, dims = dims), target)
})

test_that("rename_old_sm return expected lists", {
  target <- list(
    list(pos = c(FALSE, TRUE, rep(FALSE, 15)),
         oldname = "sds_sx1kEQ9",
         pnames = "sds_sx1_1",
         fnames = "sds_sx1_1",
         dims = numeric(0)),
    list(pos = c(rep(FALSE, 8), rep(TRUE, 9)),
         oldname = "s_sx1kEQ9",
         pnames = "s_sx1_1",
         fnames = paste0("s_sx1_1[", 1:9, "]"),
         dims = 9),
    list(pos = c(TRUE, rep(FALSE, 16)),
         oldname = "sds_sigma_t2x0",
         pnames = "sds_sigma_t2x0_1",
         fnames = "sds_sigma_t2x0_1",
         dims = numeric(0)),
    list(pos = c(FALSE, FALSE, rep(TRUE, 6), rep(FALSE, 9)),
         oldname = "s_sigma_t2x0",
         pnames = "s_sigma_t2x0_1",
         fnames = paste0("s_sigma_t2x0_1[", 1:6, "]"),
         dims = 6)
  )
  pars <- c("sds_sigma_t2x0", "sds_sx1kEQ9",
            paste0("s_sigma_t2x0[", 1:6, "]"),
            paste0("s_sx1kEQ9[", 1:9, "]"))
  dims <- list(sds_sigma_t2x0 = numeric(0), sds_sx1kEQ9 = numeric(0),
               s_sigma_t2x0 = 6, s_sx1kEQ9 = 9)
  bterms <- brmsterms(bf(y ~ s(x1, k = 9), sigma ~ t2(x0)), family = gaussian())
  dat <- data.frame(y = rnorm(100), x1 = rnorm(100), x0 = rnorm(100))
  expect_equivalent(brms:::rename_old_sm(bterms, dat, pars, dims), target)
})

test_that("rename_old_mo returns expected lists", {
  bterms <- brmsterms(bf(y ~ mo(x), sigma ~ mo(x)), family = gaussian())
  data <- data.frame(y = rnorm(10), x = rep(1:5, 2))
  pars <- c(
    "bmo_x", "bmo_sigma_x",
    paste0("simplex_x[", 1:5, "]"),
    paste0("simplex_sigma_x[", 1:5, "]")
  )
  target <- list(
    list(
      pos = c(TRUE, rep(FALSE, 11)),
      fnames = "bmo_mox"
    ),
    list(
      pos = c(FALSE, FALSE, rep(TRUE, 5), rep(FALSE, 5)),
      fnames = paste0("simo_mox1[", 1:5, "]")
    ),
    list(
      pos = c(FALSE, TRUE, rep(FALSE, 10)),
      fnames = "bmo_sigma_mox"
    ),
    list(
      pos = c(rep(FALSE, 7), rep(TRUE, 5)),
      fnames = paste0("simo_sigma_mox1[", 1:5, "]")
    )
  )
  expect_equivalent(brms:::rename_old_mo(bterms, data, pars), target)
})

test_that("rename_old_categorical works correctly", {
  dat <- data.frame(
    y = rep(c("cat1", "cat2", "cat3"), 3),
    x = rnorm(9)
  )
  fam <- categorical()
  fam$dpars <- c("mucat2", "mucat3")
  bterms <- brmsterms(bf(y ~ x) + fam)
  pars <- c("b_cat2_Intercept", "b_cat3_Intercept",
            "b_cat2_x", "b_cat3_x")
  res <- brms:::rename_old_categorical(bterms, dat, pars)
  target <- list(
    list(
      pos = rep(TRUE, 4),
      fnames = c(
        "b_mucat2_Intercept", "b_mucat3_Intercept",
        "b_mucat2_x", "b_mucat3_x"
      )
    )
  )
  expect_equivalent(res, target)
})

