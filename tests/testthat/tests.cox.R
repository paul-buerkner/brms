context("Native Cox predictions")

test_that("Cox Stan data use response-specific names", {
  skip_if_not_installed("splines2")
  dat <- data.frame(y = 1:10, y2 = 2:11, x = seq_len(10))
  form <- bf(y ~ x, family = cox()) +
    bf(y2 ~ x, family = cox()) + set_rescor(FALSE)
  sdata <- standata(form, dat)
  expect_false(anyDuplicated(names(sdata)) > 0L)
  expect_true(all(c("Zbhaz_y", "Zcbhaz_y", "Kbhaz_y", "con_sbhaz_y",
                    "Zbhaz_y2", "Zcbhaz_y2") %in% names(sdata)))
  fit <- brm(form, dat, empty = TRUE)
  frame <- brms:::brmsframe(fit$formula, data = fit$data, basis = fit$basis)
  pars <- c("b_y_Intercept", "b_y_x", "b_y2_Intercept", "b_y2_x",
             paste0("sbhaz_y[", 1:5, "]"), paste0("sbhaz_y2[", 1:5, "]"))
  draws <- posterior::as_draws_matrix(
    matrix(0.2, 3, length(pars), dimnames = list(NULL, pars))
  )
  prep <- prepare_predictions(
    frame, draws = draws, sdata = sdata,
    cox_version = fit$version$brms, tail = "constant"
  )
  expect_equal(prep$resps$y$bhaz$basis, fit$basis$resps$y$bhaz$basis_matrix)
  expect_equal(prep$resps$y2$bhaz$basis, fit$basis$resps$y2$bhaz$basis_matrix)
})

test_that("Cox preparation preserves singleton and absent strata", {
  skip_if_not_installed("splines2")
  dat <- data.frame(y = 1:10, x = seq_len(10),
                    g = factor(rep(c("a", "b"), 5)))
  fit <- brm(y | bhaz(gr = g) ~ x, data = dat, family = cox(), empty = TRUE)
  frame <- brms:::brmsframe(fit$formula, data = fit$data, basis = fit$basis)
  pars <- c(paste0("sbhaz[a,", 1:5, "]"), paste0("sbhaz[b,", 1:5, "]"))
  draws <- posterior::as_draws_matrix(
    matrix(0.2, 3, 10, dimnames = list(NULL, pars))
  )
  for (rows in list(1:2, 1L)) {
    sdata <- standata(fit, newdata = dat[rows, ])
    prep <- brms:::prepare_predictions_bhaz(
      frame$dpars$mu, draws, sdata, cox_version = fit$version$brms
    )
    expect_equal(dim(prep$bhaz), c(3L, length(rows)))
    expect_equal(prep$basis, fit$basis$bhaz$basis_matrix)
    expect_equal(prep$bhaz, unname(tcrossprod(matrix(0.2, 3, 5), sdata$Zbhaz)))
    expect_equal(prep$cbhaz, unname(tcrossprod(matrix(0.2, 3, 5), sdata$Zcbhaz)))
  }
})

test_that("Unchecked Cox responses do not extrapolate the spline", {
  skip_if_not_installed("splines2")
  dat <- data.frame(y = 1:10, x = seq_len(10))
  fit <- brm(y ~ x, dat, family = cox(), empty = TRUE)
  expect_silent(sdata <- standata(
    fit, newdata = data.frame(y = c(0, 100), x = 1:2),
    check_response = FALSE
  ))
  expect_equal(as.vector(sdata$Zcbhaz %*% rep(0.2, 5)), c(0, 1))
  expect_silent(sdata <- standata(
    fit, newdata = data.frame(x = 1:2), check_response = FALSE
  ))
  expect_true(all(is.na(sdata$Y)))
  expect_equal(as.vector(sdata$Zcbhaz %*% rep(0.2, 5)), c(0, 0))
})
