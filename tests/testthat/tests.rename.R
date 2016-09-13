test_that("rename returns an error on duplicated names", {
  expect_error(rename(c(letters[1:4],"a()","a["), check_dup = TRUE), fixed = TRUE,
               paste("Internal renaming of variables led to duplicated names.", 
                     "\nOccured for variables: a, a(), a["))
  expect_error(rename(c("aDb","a/b","b"), check_dup = TRUE), fixed = TRUE,
               paste("Internal renaming of variables led to duplicated names.", 
                     "\nOccured for variables: aDb, a/b"))
  expect_error(rename(c("log(a,b)","logab","bac","ba"), check_dup = TRUE), fixed = TRUE,
               paste("Internal renaming of variables led to duplicated names.", 
                     "\nOccured for variables: log(a,b), logab"))
})

test_that("rename perform correct renaming", {
  names <- c("acd", "a[23]", "b__")
  expect_equal(rename(names, symbols = c("[", "]", "__"), subs = c(".", ".", ":")),
               c("acd", "a.23.", "b:"))
  expect_equal(rename(names, symbols = c("^\\[", "\\]", "__$"), 
                      subs = c(".", ".", ":"), fixed = FALSE),
               c("acd", "a[23.", "b:"))
})

test_that("model_names works correctly", {
  expect_equal(model_name(NA), "brms-model")
  expect_equal(model_name(gaussian()), "gaussian(identity) brms-model")
})

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

test_that("change_prior returns expected lists", {
  pars <- c("b", "b_1", "bp", "bp_1", "prior_b", "prior_b_1", 
            "prior_b_3", "sd_x[1]", "prior_bp_1")
  expect_equal(change_prior(class = "b", pars = pars, 
                           names = c("x1", "x3", "x2")),
               list(list(pos = 6, oldname = "prior_b_1", 
                         pnames = "prior_b_x1", fnames = "prior_b_x1"),
                    list(pos = 7, oldname = "prior_b_3", 
                         pnames = "prior_b_x2", fnames = "prior_b_x2")))
  expect_equal(change_prior(class = "bp", pars = pars, 
                           names = c("x1", "x2"), new_class = "b"),
               list(list(pos = 9, oldname = "prior_bp_1", 
                         pnames = "prior_b_x1", fnames = "prior_b_x1")))
})

test_that("change_old_ranef returns expected lists", {
  data <- data.frame(y = rnorm(10), x = rnorm(10), g = 1:10)
  ee <- extract_effects(y ~ a, nonlinear = a ~ x + (1|g))
  ranef <- tidy_ranef(ee, data = data)
  pars <- c("b_a_Intercept", "b_a_x", "sd_a_g_Intercept", 
            paste0("r_a_g[", 1:10, ",Intercept]"))
  dims <- list("sd_a_g_Intercept" = numeric(0), "r_a_g" = c(10, 1))
  target <- list(
    list(pos = c(rep(FALSE, 2), TRUE, rep(FALSE, 10)),
         oldname = "sd_a_g_Intercept", pnames = "sd_g_a_Intercept",
         fnames = "sd_g_a_Intercept", dims = numeric(0)),
    list(pos = c(rep(FALSE, 3), rep(TRUE, 10)), oldname = "r_a_g",
         pnames = "r_g_a", fnames = paste0("r_g_a[", 1:10, ",Intercept]"),
         dims = c(10, 1)))
  expect_equal(change_old_ranef(ranef, pars = pars, dims = dims), target)
})
