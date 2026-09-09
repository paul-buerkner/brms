context("Centering values in updated formula environments")

re_center_update_dat <- data.frame(
  y = seq_len(12), x = seq_len(12) / 10,
  g = factor(rep(letters[1:3], each = 4)),
  h = factor(rep(letters[1:3], 4))
)

test_that("formula updates preserve old and new centering bindings", {
  old <- local({
    rho <- setNames(c(0.1, 0.2, 0.3), letters[1:3])
    sentinel <- "original"
    brms:::materialize_re_center(bf(y ~ x + (1 | gr(g, center = rho))))
  })
  incoming <- local({
    rho <- setNames(c(0.6, 0.7, 0.8), letters[1:3])
    sentinel <- "incoming"
    ~ . + (1 | gr(h, center = rho))
  })
  old_env <- environment(old$formula)
  updated <- update(old, incoming)
  rm(rho, envir = environment(incoming))
  sdata <- standata(updated, re_center_update_dat)
  expect_equal(sdata$rho_s2z_1[, 1], c(0.1, 0.2, 0.3))
  expect_equal(sdata$rho_s2z_2[, 1], c(0.6, 0.7, 0.8))
  expect_identical(environment(old$formula), old_env)
  expect_equal(get("rho", old_env),
               setNames(c(0.1, 0.2, 0.3), letters[1:3]))
  expect_equal(get("sentinel", environment(updated$formula)), "original")

  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(updated, path)
  restored <- standata(readRDS(path), re_center_update_dat)
  expect_equal(restored$rho_s2z_1, sdata$rho_s2z_1)
  expect_equal(restored$rho_s2z_2, sdata$rho_s2z_2)
})

test_that("removed centering terms do not evaluate incoming bindings", {
  old <- local({
    rho <- setNames(rep(0.2, 3), letters[1:3])
    brms:::materialize_re_center(bf(
      y ~ x + (1 | gr(g, center = rho)) + (0 + x | gr(g, center = rho))
    ))
  })
  incoming <- as.formula(
    "~ . - (1 | gr(g, center = rho))", env = new.env(parent = baseenv())
  )
  updated <- update(old, incoming)
  expect_identical(environment(updated$formula), environment(old$formula))
  sdata <- standata(updated, re_center_update_dat)
  expect_equal(sdata$M_1, 1L)
  expect_equal(sdata$rho_s2z_1[, 1], rep(0.2, 3))
})

test_that("new centering symbols still support later formula subtraction", {
  old <- bf(y ~ x)
  incoming <- local({
    rho_new <- setNames(c(0.2, 0.3, 0.4), letters[1:3])
    ~ . + (1 | gr(g, center = rho_new))
  })
  added <- update(old, incoming)
  removed <- update(added, ~ . - (1 | gr(g, center = rho_new)))
  expect_equal(removed$formula, old$formula)

  repeated <- ~ . + (1 | gr(h, center = rho_new))
  added_again <- update(added, repeated)
  expect_equal(standata(added_again, re_center_update_dat)$rho_s2z_2[, 1],
               c(0.2, 0.3, 0.4))
  removed_again <- update(added_again, ~ . - (1 | gr(h, center = rho_new)))
  expect_equal(removed_again$formula, added$formula, ignore_attr = TRUE)

  shifted <- local({
    shift <- 0.1
    ~ . + (1 | gr(h, center = rho_new + shift))
  })
  shifted <- update(added, shifted)
  expect_equal(standata(shifted, re_center_update_dat)$rho_s2z_2[, 1],
               c(0.3, 0.4, 0.5))
})

test_that("centering expressions use incoming environments in formula updates", {
  old <- bf(y ~ x)
  incoming <- local({
    weights <- c(0.2, 0.4, 0.6)
    bf(
      ~ . + (1 | gr(g, center = setNames(weights, letters[1:3]))),
      sigma ~ 1 + (1 | gr(h, center = weights))
    )
  })
  updated <- brms:::materialize_re_center(update(old, incoming))
  rm(weights, envir = environment(incoming$formula))
  sdata <- standata(updated, re_center_update_dat)
  expect_equal(sdata$rho_s2z_1[, 1], c(0.2, 0.4, 0.6))
  expect_equal(sdata$rho_s2z_2[, 1], c(0.2, 0.4, 0.6))

  ignored <- as.formula(
    "~ (1 | gr(g, center = absent_rho))", env = new.env(parent = baseenv())
  )
  expect_equal(update(old, ignored, mode = "keep"), old)
  replacement <- local({
    rho <- setNames(rep(0.7, 3), letters[1:3])
    y ~ x + (1 | gr(g, center = rho))
  })
  replaced <- update(old, replacement, mode = "replace")
  expect_equal(standata(replaced, re_center_update_dat)$rho_s2z_1[, 1],
               rep(0.7, 3))
})
