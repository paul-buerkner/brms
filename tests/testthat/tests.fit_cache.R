context("Tests for fit cache")


test_that("normalize_stancode", {
  expect_equal(
    normalize_stancode("// a\nb;\n  b + c = 4; // kde\ndata"),
    normalize_stancode("// dasflkjldl\n   // adsfadsfa\n b;\n\n  \n  \t\rb + c = 4;\ndata")
  )
  expect_equal(
    normalize_stancode("data /* adfa */ {\nint a;\n /* asdddede \n asdfas \n asf */}\n"),
    normalize_stancode("data {\nint a;\n} /* aa \n adfasdf \n asdfadsf ddd */\n")
  )
  expect_equal(
    normalize_stancode("data \n {\nint a;\n\n }  \t\n"),
    normalize_stancode("data {\nint a;\n} \n")
  )
  expect_equal(
    normalize_stancode("/* \n\n */\na*/"),
    normalize_stancode("a*/")
  )
  expect_equal(
    normalize_stancode("//adsfadf \ra // asdfasdf\r\n"),
    normalize_stancode("a")
  )
  expect_equal(
    normalize_stancode("/* * \n * \n * fg / */hhh"),
    normalize_stancode("hhh")
  )
  expect_equal(
    normalize_stancode("a //b"),
    normalize_stancode("a")
  )
  expect_false(normalize_stancode("// a\ndata {\nint a;\n}\n") ==
                 normalize_stancode("// a\ndata {\nint b;\n}\n"))
  #Should not remove single whitespace
  expect_false(normalize_stancode("da ta") ==
                 normalize_stancode("data"))
  #Should handle wrong nested comments
  expect_false(normalize_stancode("/* \n\n */\na*/") ==
                 normalize_stancode("b*/"))
})

test_that("cache_read_write_check", {
  cache_tmp <- tempfile(fileext = ".rds")
  
  expect_null(read_brmsfit(cache_tmp))

  saveRDS(list(a = 1), file = cache_tmp)
  expect_error(read_brmsfit(cache_tmp))
  
  data_model1 <- data.frame(y = rnorm(10), x = rnorm(10))
  fake_fit <- brm(y ~ x, data = data_model1,
                      empty = TRUE)
  
  fake_fit_file <- fake_fit
  fake_fit_file$file <- cache_tmp
  
  scode_model1 <- make_stancode(y ~ x, data = data_model1)
  sdata_model1 <- make_standata(y ~ x, data = data_model1)

  data_model2 <- data_model1
  data_model2$x[1] <- data_model2$x[1] + 1
  scode_model2 <- make_stancode(y ~ 0 + x, data = data_model2)
  sdata_model2 <- make_standata(y ~ 0 + x, data = data_model2)
  
    
  write_brmsfit(fake_fit, file = cache_tmp)
  
  cache_res <- read_brmsfit(file = cache_tmp)
  
  expect_equal(cache_res, fake_fit_file)
  
  expect_false(brmsfit_needs_refit(
    cache_res, sdata = sdata_model1, scode = scode_model1, 
    algorithm = "sampling", silent = TRUE))
  expect_false(brmsfit_needs_refit(
    cache_res, sdata = sdata_model1, scode = scode_model1, algorithm = NULL,
    silent = TRUE))
  expect_false(brmsfit_needs_refit(
    cache_res, sdata = sdata_model1, scode = NULL, 
    algorithm = "sampling", silent = TRUE))
  expect_false(brmsfit_needs_refit(
    cache_res, sdata = NULL, scode = scode_model1, algorithm = "sampling",
    silent = TRUE))
  
  
  
  expect_true(brmsfit_needs_refit(
    cache_res, sdata = sdata_model2, scode = scode_model1, 
    algorithm = "sampling", silent = TRUE))
  expect_true(brmsfit_needs_refit(
    cache_res, sdata = sdata_model1, scode = scode_model2, 
    algorithm = "sampling", silent = TRUE))
  expect_true(brmsfit_needs_refit(
    cache_res, sdata = sdata_model2, scode = scode_model2, 
    algorithm = "sampling", silent = TRUE))
  expect_true(brmsfit_needs_refit(
    cache_res, sdata = sdata_model1, scode = scode_model1, 
    algorithm = "optimize", silent = TRUE))

  expect_true(brmsfit_needs_refit(
    cache_res, sdata = make_standata(y ~ x, data = data_model1, 
                                     sample_prior = "only"), 
    scode = scode_model1, algorithm = NULL, silent = TRUE))
  
})