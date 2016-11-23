test_that("p performs correct indexing", {
  expect_equal(p(1:10), 1:10)
  x <- rnorm(10)
  expect_equal(p(x, i = 3), x[3])
  A <- matrix(x, nrow = 5)
  expect_equal(p(A, i = 3), A[3, , drop = FALSE]) 
  expect_equal(p(A, i = 2, row = FALSE), A[, 2, drop = FALSE]) 
})

test_that("rmNULL removes all NULL entries", {
  expect_equal(rmNULL(list(a = NULL, b = 1, c = list(NULL, 1))),
               list(b = 1, c = list(1)))
  expect_equal(rmNULL(list(a = NULL, b = 1, c = NULL)),
               list(b = 1))
})

test_that("rmNum remove all numeric entries", {
  expect_equal(rmNum(list(1, "a", 2.3, "b")), list("a","b"))
  expect_equal(rmNum(list(x = 1.5, y = "abc", z = pi)), list(y = "abc"))
})

test_that("forumla2string performs correct conversion", {
  expect_equal(formula2string("y~x"), "y~x")
  expect_equal(formula2string(y ~ x + c), "y~x+c")
  expect_equal(formula2string(abc ~ x + cd, rm = c(3,2)), "~x+")
})

test_that("collapse_lists performs correct collapsing after names", {
  x <- list(a = "a <- ", b = "b <- ")
  y <- list(b = "cauchy(1,2)", c = "normal(0,1)", a = "gamma(1,1)")
  expect_equal(collapse_lists(list()), list())
  expect_equal(collapse_lists(list(x, y)), 
               list(a = "a <- gamma(1,1)", b = "b <- cauchy(1,2)", c = "normal(0,1)"))
  expect_equal(collapse_lists(list(c(x, c = "c <- "), y)),
               list(a = "a <- gamma(1,1)", b = "b <- cauchy(1,2)", c = "c <- normal(0,1)"))
})

test_that("subset_attr works correctly", {
  x <- list(a = 1, b = 2, c = 3)
  attr(x, "att") <- "x"
  res <- subset_attr(x, c("a", "c"))
  expect_equivalent(res, list(a = 1, c = 3))
  expect_equal(attr(res, "att"), "x")
  expect_equal(names(res), c("a", "c"))
})

test_that("nlist works correctly", {
  x <- 1
  y <- 2:3
  exlist <- list(x = x, y = y)
  expect_equal(nlist(x = x, y = y), exlist)
  expect_equal(nlist(x, y), exlist)
  expect_equal(nlist(x = x, y), exlist)
})

test_that("convenience functions for model families work correctly", {
  expect_true(use_real(gaussian()))
  expect_true(!use_real("poisson"))
  expect_true(use_int(binomial()))
  expect_true(has_trials("zero_inflated_binomial"))
  expect_true(has_cat("acat"))
  expect_true(has_sigma(student()))
  effects <- list(se = TRUE)
  expect_true(!has_sigma("cauchy", effects = effects))
  expect_true(has_sigma("cauchy", effects = effects, 
                        autocor = cor_ar(cov = TRUE)))
})

test_that("use_alias works correctly", {
  a <- 2
  b <- 3
  expect_warning(use_alias(a, b), fixed = TRUE,
                 "'b' is deprecated. Please use argument 'a' instead.")
  dots <- list(c = 1)
  expect_warning(use_alias(a, dots$c), fixed = TRUE,
                 "'c' is deprecated. Please use argument 'a' instead.")
  expect_equal(use_alias(a, dots$c, warn = FALSE), dots$c)
})

test_that("rhs keeps attributes", {
  form <- structure(y~x, test = TRUE)
  expect_equal(attributes(form), attributes(rhs(form)))
})

test_that("lsp works correctly", {
  expect_equal(lsp("base", pattern = "^log"),
               c("log", "log10", "log1p", "log2", "logb", "logical"))
  expect_equal(lsp("brms", pattern = "^log_"),
               c("log_diff_exp", "log_inv_logit", 
                 "log_lik.brmsfit", "log_posterior.brmsfit",
                 "log_sum_exp"))
})

test_that(".addition and .cat works correctly", {
  expect_equal(.addition(~ brms:::.cat(x), data = data.frame(x = 2:3)), 2:3)
  expect_error(.addition(~ brms:::.cat(x), data = data.frame(x = -2)),
               "number of categories must be positive integers")
})