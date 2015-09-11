test_that("Test that rename returns an error on duplicated names", {
  expect_error(rename(c(letters[1:4],"a()","a["), check_dup = TRUE), fixed = TRUE,
               "Internal renaming of variables led to duplicated names. \nOccured for variables: a, a(), a[")
  expect_error(rename(c("aDb","a/b","b"), check_dup = TRUE), fixed = TRUE,
               "Internal renaming of variables led to duplicated names. \nOccured for variables: aDb, a/b")
  expect_error(rename(c("log(a,b)","logab","bac","ba"), check_dup = TRUE), fixed = TRUE,
               "Internal renaming of variables led to duplicated names. \nOccured for variables: log(a,b), logab")
})

test_that("Test that rename perform correct renaming", {
  names <- c("acd", "a[23]", "b__")
  expect_equal(rename(names, symbols = c("[", "]", "__"), subs = c(".", ".", ":")),
               c("acd", "a.23.", "b:"))
  expect_equal(rename(names, symbols = c("^\\[", "\\]", "__$"), subs = c(".", ".", ":"), fixed = FALSE),
               c("acd", "a[23.", "b:"))
})

test_that("Test that make_group_frame returns correct first and last indices", {
  expect_equal(make_group_frame(list(a = c("x","Int"), b = c("x"))),
               data.frame(g = c("a", "b"), first = c(1, 1), last = c(2, 1)))
  expect_equal(make_group_frame(list(a = c("x","Int"), b = c("x"), a = list("y","z"), b = list("b"))),
               data.frame(g = c("a", "b", "a", "b"), first = c(1, 1, 3, 2), last = c(2, 1, 4, 2)))
  expect_equal(make_group_frame(list(a = c("x","Int"), b = c("x"), a = list("y","z"), a = list("b"))),
               data.frame(g = c("a", "b", "a", "a"), first = c(1, 1, 3, 5), last = c(2, 1, 4, 5)))
})

test_that("Test that make_indices returns correct 1 and 2 dimensional indices", {
  expect_equal(make_indices(rows = 1:2), c("[1]", "[2]"))
  expect_equal(make_indices(rows = 1:2, cols = 1:3, dim = 1), c("[1]", "[2]"))
  expect_equal(make_indices(rows = 1:2, cols = 1:3, dim = 2), 
               c("[1,1]", "[2,1]", "[1,2]", "[2,2]", "[1,3]", "[2,3]"))
})

test_that("Test that combine_duplicates works as expected", {
  expect_equal(combine_duplicates(list(a = c(2,2), b = c("a", "c"))),
               list(a = c(2,2), b = c("a", "c")))
  expect_equal(combine_duplicates(list(a = c(2,2), b = c("a", "c"), a = c(4,2))),
               list(a = c(2,2,4,2), b = c("a", "c")))
})

test_that("Test that prior_names returns correct lists to be understood by rename_pars", {
  pars <- c("b", "b_1", "bp", "bp_1", "prior_b", "prior_b_1", "prior_b_3", "sd_x[1]", "prior_bp_1")
  expect_equal(prior_names(class = "b", pars = pars, names = c("x1", "x3", "x2")),
               list(list(pos = 6, oldname = "prior_b_1", pnames = "prior_b_x1", fnames = "prior_b_x1"),
                    list(pos = 7, oldname = "prior_b_3", pnames = "prior_b_x2", fnames = "prior_b_x2")))
  expect_equal(prior_names(class = "bp", pars = pars, names = c("x1", "x2"), new_class = "b"),
               list(list(pos = 9, oldname = "prior_bp_1", pnames = "prior_b_x1", fnames = "prior_b_x1")))
})


