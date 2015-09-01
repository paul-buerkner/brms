test_that("Test that rename returns an error on duplicated names", {
  expect_error(rename(c(letters[1:4],"a()","a["), check_dup = TRUE), fixed = TRUE,
               "Internal renaming of variables led to duplicated names. \nOccured for variables: a, a(), a[")
  expect_error(rename(c("aDb","a/b","b"), check_dup = TRUE), fixed = TRUE,
               "Internal renaming of variables led to duplicated names. \nOccured for variables: aDb, a/b")
  expect_error(rename(c("log(a,b)","logab","bac","ba"), check_dup = TRUE), fixed = TRUE,
               "Internal renaming of variables led to duplicated names. \nOccured for variables: log(a,b), logab")
})

