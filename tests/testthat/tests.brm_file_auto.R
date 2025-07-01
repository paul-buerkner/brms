
# file_auto
test_that("file_auto option works", {
   skip("Temporarily disabled for debugging reasons")
  # test_cache_dir <- tempdir()
  # options(brms.cache_folder = test_cache_dir)
    epilepsy2 <- epilepsy[-c(1), ]
    # same
    f1 <- brm_mock_for_file_auto(count ~ zAge + zBase * Trt + (1|patient),
                           data = epilepsy, family = poisson() ,    file_auto = TRUE)
    f2 <- brm_mock_for_file_auto(count ~ zAge + zBase * Trt + (1|patient),
                                 data = epilepsy, family = gaussian() ,    file_auto = TRUE)
    expect_equal(f1$file, f2$file)
    # different data
    f1 <- brm_mock_for_file_auto(count ~ zAge + zBase * Trt + (1|patient),
                                 data = epilepsy, family = poisson(), file_auto = TRUE)
    f2 <- brm_mock_for_file_auto(count ~ zAge + zBase * Trt + (1|patient),
                                 data = epilepsy2, family = poisson(), file_auto = TRUE)
    expect_false(f1$file == f2$file)
    # different family
    f1 <- brm_mock_for_file_auto(count ~ zAge + zBase * Trt + (1|patient),
                                 data = epilepsy, family = poisson(), file_auto = TRUE)
    f2 <- brm_mock_for_file_auto(count ~ zAge + zBase * Trt + (1|patient),
                                 data = epilepsy, family = gaussian(), file_auto = TRUE)
    expect_equal(f1$file, f2$file)
})

test_that("hash function check" , {
  aa1 =   hash_dots( count ~  zAge + zBase * Trt + (1|patient),
                         data = epilepsy, family = poisson() ,    file_auto = TRUE )
  aa2 =   hash_dots( count ~  zAge + zBase * Trt + (1|patient),
                         data = epilepsy, family = poisson() ,    file_auto = TRUE )
  expect_equal( aa1 , aa2 )
  g1 =  gaussian()
  g2 = gaussian()
  expect_true(hash_dots(g1) == hash_dots(g2))
  expect_true(hash_dots(gaussian()) == hash_dots(gaussian()))
  expect_false(hash_dots(poisson()) == hash_dots(gaussian()))
  expect_true(hash_dots(epilepsy) == hash_dots(epilepsy))
  expect_false(hash_dots(epilepsy[-c(2) , ]) == hash_dots(epilepsy))
})

test_that("identical arguments give identical hash", {
  a <- nlist(formula = f1, data = d1, family = fam1)
  h1 <- hash_dots(a, data_policy = "hash")
  h2 <- hash_dots(a, data_policy = "hash")
  expect_identical(h1, h2)
})

test_that("ordering of list elements is irrelevant", {
  a <- nlist(formula = f1, data = d1, family = fam1)
  b <- a[c("family", "data", "formula")]   # reorder
  h1 <- hash_dots(a, data_policy = "hash")
  h2 <- hash_dots(b, data_policy = "hash")
  expect_identical(h1, h2)
})

test_that("formula environment ignored", {
  env <- new.env()
  environment(f1) <- env
  a <- nlist(formula = f1, data = d1, family = fam1)
  b <- nlist(formula = mpg ~ wt, data = d1, family = fam1)
  h1 <- hash_dots(a, data_policy = "hash")
  h2 <- hash_dots(b, data_policy = "hash")
  expect_identical(h1, h2)
})

test_that("changing *one* statistical input changes the hash", {
  # toy objects reused across tests
  f1 <- mpg ~ wt
  f2 <- mpg ~ wt + cyl            # same data, different formula
  d1 <- mtcars
  d2 <- mtcars[sample(nrow(mtcars)), ]   # same rows, shuffled order
  fam1 <- gaussian()
  fam2 <- student()
  base <- nlist(formula = f1, data = d1, family = fam1)
  h_base <- hash_dots(base, data_policy = "hash")
  # 1. Different formula
  h_formula <- hash_dots(
    modifyList(base, list(formula = f2)), data_policy = "hash"
  )
  expect_false(identical(h_base, h_formula))
  # 2. Different family
  h_family <- hash_dots(
    modifyList(base, list(family = fam2)), data_policy = "hash"
  )
  expect_false(identical(h_base, h_family))
  # 3. Different data (row order → different digest with data_policy = "hash")
  h_data <- hash_dots(
    modifyList(base, list(data = d2)), data_policy = "hash"
  )
  expect_false(identical(h_base, h_data))
})



