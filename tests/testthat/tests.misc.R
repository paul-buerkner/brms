context("Tests for miscellaneous functions")

test_that("p performs correct indexing", {
  expect_equal(p(1:10), 1:10)
  x <- rnorm(10)
  expect_equal(p(x, i = 3), x[3])
  A <- matrix(x, nrow = 5)
  expect_equal(p(A, i = 3), A[3, , drop = FALSE])
  expect_equal(p(A, i = 2, row = FALSE), A[, 2, drop = FALSE])
})

test_that("rmNULL removes all NULL entries", {
  expect_equal(
    rmNULL(list(a = NULL, b = 1, c = list(NULL, 1))),
    list(b = 1, c = list(1))
  )
  expect_equal(
    rmNULL(list(a = NULL, b = 1, c = NULL)),
    list(b = 1)
  )
})

test_that("rename returns an error on duplicated names", {
  expect_error(rename(c(letters[1:4], "a()", "a["), check_dup = TRUE),
    fixed = TRUE,
    paste("Occured for: 'a', 'a()', 'a['")
  )
  expect_error(rename(c("aDb", "a/b", "b"), check_dup = TRUE),
    fixed = TRUE,
    paste("Occured for: 'aDb', 'a/b'")
  )
  expect_error(rename(c("log(a,b)", "logab", "bac", "ba"), check_dup = TRUE),
    fixed = TRUE,
    paste("Occured for: 'log(a,b)', 'logab'")
  )
})

test_that("rename perform correct renaming", {
  names <- c("acd", "a[23]", "b__")
  expect_equal(
    rename(names, c("[", "]", "__"), c(".", ".", ":")),
    c("acd", "a.23.", "b:")
  )
  expect_equal(
    rename(names, c("^\\[", "\\]", "__$"), c(".", ".", ":"), fixed = FALSE),
    c("acd", "a[23.", "b:")
  )
})

test_that("collapse_lists performs correct collapsing after names", {
  x <- list(a = "a <- ", b = "b <- ")
  y <- list(b = "cauchy(1,2)", c = "normal(0,1)", a = "gamma(1,1)")
  expect_equal(collapse_lists(list()), list())
  expect_equal(
    collapse_lists(x, y),
    list(
      a = "a <- gamma(1,1)", b = "b <- cauchy(1,2)",
      c = "normal(0,1)"
    )
  )
  expect_equal(
    collapse_lists(ls = list(c(x, c = "c <- "), y)),
    list(
      a = "a <- gamma(1,1)", b = "b <- cauchy(1,2)",
      c = "c <- normal(0,1)"
    )
  )
})

test_that("nlist works correctly", {
  x <- 1
  y <- 2:3
  exlist <- list(x = x, y = y)
  expect_equal(nlist(x = x, y = y), exlist)
  expect_equal(nlist(x, y), exlist)
  expect_equal(nlist(x = x, y), exlist)
})

test_that("use_alias works correctly", {
  a <- 2
  b <- 3
  expect_warning(use_alias(a, b),
    fixed = TRUE,
    "'b' is deprecated. Please use argument 'a' instead."
  )
  dots <- list(c = 1)
  expect_warning(use_alias(a, dots$c),
    fixed = TRUE,
    "'c' is deprecated. Please use argument 'a' instead."
  )
  expect_equal(use_alias(a, dots$c, warn = FALSE), dots$c)
})

test_that("rhs keeps attributes", {
  form <- structure(y ~ x, test = TRUE)
  expect_equal(attributes(form), attributes(rhs(form)))
})

test_that("lsp works correctly", {
  expect_equal(
    lsp("base", pattern = "^log"),
    c("log", "log10", "log1p", "log2", "logb", "logical")
  )
  expect_equal(
    lsp("brms", pattern = "^inv_logit"),
    c("inv_logit", "inv_logit_scaled")
  )
})
