context("Tests for renaming helper functions")

test_that("make_index_names returns correct 1 and 2 dimensional indices", {
  expect_equal(make_index_names(rownames = 1:2), c("[1]", "[2]"))
  expect_equal(
    make_index_names(rownames = 1:2, colnames = 1:3, dim = 1),
    c("[1]", "[2]")
  )
  expect_equal(
    make_index_names(rownames = c("a", "b"), colnames = 1:3, dim = 2),
    c("[a,1]", "[b,1]", "[a,2]", "[b,2]", "[a,3]", "[b,3]")
  )
})
