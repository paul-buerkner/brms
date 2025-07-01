context("Tests for brmsfit helper functions")

test_that("first_greater returns expected results", {
  A <- cbind(1:10, 11:20, 21:30)
  x <- c(5, 25, 7, 15, 7, 10, 15, 19, 3, 11)
  expect_equal(first_greater(A, x), c(2, 3, 2, 3, 2, 2, 2, 3, 1, 2))
  expect_equal(first_greater(A, x, i = 2), c(2, 3, 2, 3, 2, 2, 2, 3, 2, 2))
})

test_that("array2list performs correct conversion", {
  A <- array(1:27, dim = c(3,3,3))
  B <- list(matrix(1:9,3,3), matrix(10:18,3,3), matrix(19:27,3,3))
  expect_equal(brms:::array2list(A), B)
})

test_that("probit and probit_approx produce similar results", {
  expect_equal(brms:::inv_link(-10:10, "probit"),
               brms:::inv_link(-10:10, "probit_approx"),
               tolerance = 1e-3)
})

test_that("autocorrelation matrices are computed correctly", {
  ar <- 0.5
  ma <- 0.3

  ar_mat <- brms:::get_cor_matrix_ar1(ar = matrix(ar), nobs = 4)
  expected_ar_mat <- 1 / (1 - ar^2) *
    cbind(c(1, ar, ar^2, ar^3),
          c(ar, 1, ar, ar^2),
          c(ar^2, ar, 1, ar),
          c(ar^3, ar^2, ar, 1))
  expect_equal(ar_mat[1, , ], expected_ar_mat)

  ma_mat <- brms:::get_cor_matrix_ma1(ma = matrix(ma), nobs = 4)
  expected_ma_mat <- cbind(c(1+ma^2, ma, 0, 0),
                           c(ma, 1+ma^2, ma, 0),
                           c(0, ma, 1+ma^2, ma),
                           c(0, 0, ma, 1+ma^2))
  expect_equal(ma_mat[1, , ], expected_ma_mat)

  arma_mat <- brms:::get_cor_matrix_arma1(
    ar = matrix(ar), ma = matrix(ma), nobs = 4
  )
  g0 <- 1 + ma^2 + 2 * ar * ma
  g1 <- (1 + ar * ma) * (ar + ma)
  expected_arma_mat <- 1 / (1 - ar^2) *
    cbind(c(g0, g1, g1 * ar, g1 * ar^2),
          c(g1, g0, g1, g1 * ar),
          c(g1 * ar, g1, g0, g1),
          c(g1 * ar^2, g1 * ar, g1, g0))
  expect_equal(arma_mat[1, , ], expected_arma_mat)

  cosy <- 0.6
  cosy_mat <- brms:::get_cor_matrix_cosy(cosy = as.matrix(cosy), nobs = 4)
  expected_cosy_mat <- matrix(cosy, 4, 4)
  diag(expected_cosy_mat) <- 1
  expect_equal(cosy_mat[1, , ], expected_cosy_mat)

  ident_mat <- brms:::get_cor_matrix_ident(ndraws = 10, nobs = 4)
  expected_ident_mat <- diag(1, 4)
  expect_equal(ident_mat[1, , ], expected_ident_mat)
})

test_that("evidence_ratio returns expected results", {
  ps <- -4:10
  prs <- -2:12
  expect_true(evidence_ratio(ps, prior_samples = prs) > 1)
  expect_true(is.na(evidence_ratio(ps)))
  expect_equal(evidence_ratio(ps, cut = 0.5, wsign = "greater"), 10/5)
  expect_equal(evidence_ratio(ps, cut = 0.5, wsign = "less"), 5/10)
})

test_that("find_vars finds all valid variable names in a string", {
  string <- "x + b.x - .5 + abc(a__3) : 1/2 - 0.2"
  expect_equal(find_vars(string), c("x", "b.x", "a__3"))
})

test_that(".predictor_arma runs without errors", {
  ns <- 20
  nobs <- 30
  Y = rnorm(nobs)
  J_lag = c(1:3, 3, 3, rep(c(0:3, 3), 4), 0:3, 0)
  ar <- matrix(rnorm(ns * 3), nrow = ns, ncol = 3)
  ma <- matrix(rnorm(ns * 1), nrow = ns, ncol = 1)
  eta <- matrix(rnorm(ns * nobs), nrow = ns, ncol = nobs)
  expect_equal(.predictor_arma(eta, Y = Y, J_lag = J_lag), eta)
  expect_silent(.predictor_arma(eta, Y = Y, J_lag = J_lag, ar = ar))
  expect_silent(.predictor_arma(eta, Y = Y, J_lag = J_lag, ma = ma))
  expect_silent(.predictor_arma(eta, Y = Y, J_lag = J_lag, ar = ar, ma = ma))
})

test_that("make_conditions works correctly", {
  conds <- make_conditions(epilepsy, c("zBase", "zAge"))
  expect_equal(dim(conds), c(9, 3))
  expect_equal(conds$cond__[3], "zBase = -1 & zAge = 1")
})

test_that("brmsfit_needs_refit works correctly", {
  cache_tmp <- tempfile(fileext = ".rds")

  expect_null(read_brmsfit(cache_tmp))

  saveRDS(list(a = 1), file = cache_tmp)
  expect_error(read_brmsfit(cache_tmp))

  data_model1 <- data.frame(y = rnorm(10), x = rnorm(10))
  fake_fit <- brm(y ~ x, data = data_model1, empty = TRUE)
  test_that("summary works without error", {
    expect_error(summary(fake_fit), NA)
  })
  fake_fit_file <- fake_fit
  # align windows with unix encoding of file paths
  fake_fit_file$file <- gsub("\\", "/", cache_tmp, fixed = TRUE)

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

test_that("insert_refcat() works correctly", {
  source(testthat::test_path(file.path("helpers", "insert_refcat_ch.R")))
  source(testthat::test_path(file.path("helpers", "simopts_catlike.R")))
  for (ndraws in ndraws_vec) {
    for (nobsv in nobsv_vec) {
      for (ncat in ncat_vec) {
        cats <- paste0("cat", 1:ncat)
        ref_list <- list(
          ref1 = 1,
          reflast = ncat
        )
        fam_list <- list(
          fam_ref1 = categorical(refcat = cats[1]),
          fam_reflast = categorical(refcat = cats[ncat])
        )
        if (ncat > 2) {
          ref_list <- c(ref_list, list(ref2 = 2))
          fam_list <- c(fam_list, list(fam_ref2 = categorical(refcat = cats[2])))
        }
        eta_test_list <- list(array(rnorm(ndraws * nobsv * (ncat - 1)),
                                    dim = c(ndraws, nobsv, ncat - 1)))
        if (nobsv == 1) {
          eta_test_list <- c(
            eta_test_list,
            list(matrix(rnorm(ndraws * (ncat - 1)), nrow = ndraws))
          )
        }
        for (eta_test in eta_test_list) {
          for (i in seq_along(fam_list)) {
            # Emulate content of `fam` after fit:
            fam <- fam_list[[i]]
            if (is.null(fam$refcat)) {
              fam$refcat <- cats[1]
            }
            fam$cats <- cats
            ref <- ref_list[[i]]

            # Perform the check:
            eta_ref <- insert_refcat(eta_test, ref)
            eta_ref_ch <- insert_refcat_ch(eta_test, fam)
            expect_equivalent(eta_ref, eta_ref_ch)
            if (length(dim(eta_test)) == 3) {
              expect_equal(dim(eta_ref), c(ndraws, nobsv, ncat))
            } else if (length(dim(eta_test)) == 2) {
              expect_equal(dim(eta_ref), c(ndraws, ncat))
            }
          }
        }
      }
    }
  }
})

# split_folder_and_file
test_that("split_folder_and_file returns expected results", {
  cache_folder <- getOption('brms.cache_folder' , default = '.')
  files <- c("somefile",  "./somefile" , "somepath/somefolder/somefile")
  result <- base::lapply(files, split_folder_and_file)
  exp_result <-   list( list(folder = cache_folder, file = "somefile") ,
                        list(folder = cache_folder, file = "somefile") ,
                        list(folder = "somepath/somefolder", file = "somefile")
  )
  expect_equal(result, exp_result)
})

# check_brmsfit_file
test_that("check_brmsfit_file returns expected results", {
  cache_folder <- getOption('brms.cache_folder', default = '.')
  files <- c("somefile",  "./somefile", "somefile.rds", "somepath/somefolder/somefile")
  result <- base::lapply(files, function(x) check_brmsfit_file(x, .check_folder = F))
  exp_result <-   list( file.path(cache_folder, "somefile.rds"),
                        file.path(cache_folder, "somefile.rds"),
                        file.path(cache_folder, "somefile.rds"),
                        "somepath/somefolder/somefile.rds")
  expect_equal(result, exp_result)
})

# get_cache_folder
test_that("get_cache_folder returns expected results", {
  cache_folder <- getOption('brms.cache_folder', default = '.')
  files <- c("somefile", "./somefile", "abcde/somefile.rds",
             "somepath/somefolder/somefile")
  result <- base::lapply(files, brms:::get_cache_folder)
  exp_result <- list(cache_folder, cache_folder, "abcde", "somepath/somefolder")
  expect_equal(result, exp_result)
})
