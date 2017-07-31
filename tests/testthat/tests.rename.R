test_that("make_index_names returns correct 1 and 2 dimensional indices", {
  expect_equal(make_index_names(rownames = 1:2), c("[1]", "[2]"))
  expect_equal(make_index_names(rownames = 1:2, colnames = 1:3, dim = 1), 
               c("[1]", "[2]"))
  expect_equal(make_index_names(rownames = c("a","b"), colnames = 1:3, dim = 2), 
               c("[a,1]", "[b,1]", "[a,2]", "[b,2]", "[a,3]", "[b,3]"))
})

test_that("combine_duplicates works as expected", {
  expect_equal(combine_duplicates(list(a = c(2,2), b = c("a", "c"))),
               list(a = c(2,2), b = c("a", "c")))
  expect_equal(combine_duplicates(list(a = c(2,2), b = c("a", "c"), a = c(4,2))),
               list(a = c(2,2,4,2), b = c("a", "c")))
})

test_that("rm_int_fe works as expected", {
  dat <- data.frame(y = 1:3, x = rnorm(3))
  code <- make_stancode(y ~ 1, data = dat)
  expect_equal(rm_int_fe("Intercept", code), character(0))
  code <- make_stancode(cbind(y, y) ~ x, data = dat)
  expect_equal(rm_int_fe(c("Intercept", "x"), code, nlpar = "y"), "x")
  code <- make_stancode(y ~ x, data = dat, family = sratio())
  expect_equal(rm_int_fe(c("Intercept", "x"), code), "x")
  code <- make_stancode(y ~ 0 + intercept + x, data = dat)
  expect_equal(rm_int_fe(c("intercept", "x"), code), 
               c("intercept", "x"))
})

test_that("change_prior returns expected lists", {
  pars <- c("b", "b_1", "bp", "bp_1", "prior_b", "prior_b_1", 
            "prior_b_3", "sd_x[1]", "prior_bp_1")
  expect_equal(
    change_prior(class = "b", pars = pars, names = c("x1", "x3", "x2")),
    list(list(pos = 6, fnames = "prior_b_x1"),
         list(pos = 7, fnames = "prior_b_x2"))
  )
  expect_equal(
    change_prior(class = "bp", pars = pars, 
                 names = c("x1", "x2"), new_class = "b"),
    list(list(pos = 9, fnames = "prior_b_x1")))
})

test_that("change_old_re and change_old_re2 return expected lists", {
  data <- data.frame(y = rnorm(10), x = rnorm(10), g = 1:10)
  bterms <- parse_bf(bf(y ~ a, a ~ x + (1+x|g), 
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
  expect_equal(brms:::change_old_re(ranef, pars = pars, dims = dims), target)
  
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
  expect_equal(brms:::change_old_re2(ranef, pars = pars, dims = dims), target)
})

test_that("change_old_sm return expected lists", {
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
  bterms <- parse_bf(bf(y ~ s(x1, k = 9), sigma ~ t2(x0)), 
                     family = gaussian())
  expect_equal(brms:::change_old_sm(bterms, pars, dims), target)
})
