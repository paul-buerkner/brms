context("Predictor plans for physical sum-to-zero group-level effects")

s2z_plan_dat <- data.frame(
  y = sin(seq_len(24L)),
  x = rep(c(-1, 1), 12L),
  w = cos(seq_len(24L)),
  z = seq_len(24L) / 24,
  g = factor(rep(seq_len(4L), each = 6L)),
  h = factor(rep(seq_len(6L), each = 4L))
)

test_that("one plan owns the original coordinates of all local blocks", {
  form <- y ~ x + w + z +
    (1 | gr(g, s2z = TRUE)) + (0 + z | gr(h, s2z = TRUE))
  bfl <- brmsframe(brmsterms(form), data = s2z_plan_dat)$dpars$mu
  plan <- re_s2z_plan(bfl)

  expect_identical(plan, bfl$frame$re_s2z_plan)
  expect_null(bfl$frame[["re_s2z"]])
  expect_identical(plan$qnames, c("Intercept", "x", "w", "z"))
  expect_identical(plan$original_q, 1:4)
  expect_identical(plan$fixef_q, 2:4)
  expect_identical(plan$active_q, c(1L, 4L))
  expect_identical(plan$inactive_q, c(2L, 3L))
  expect_identical(plan$active_index, c(1L, NA_integer_, NA_integer_, 2L))
  expect_identical(plan$inactive_index, c(NA_integer_, 1L, 2L, NA_integer_))
  expect_length(plan$blocks, 2L)
  expect_identical(plan$blocks[[1L]]$match_active, 1L)
  expect_identical(plan$blocks[[2L]]$match_active, 2L)
  expect_identical(plan$blocks[[2L]]$block_active_q, c(1L, 4L))
  for (block in plan$blocks) {
    expect_false(any(c(
      "qnames", "active_q", "inactive_q", "prior", "center", "p"
    ) %in% names(block)))
  }

  infos <- re_s2z_plan_infos(plan)
  expect_identical(infos, re_s2z_infos(bfl))
  expect_identical(infos[[2L]]$qnames, plan$qnames)
  expect_identical(infos[[2L]]$active_q, plan$active_q)
  expect_identical(re_s2z_info(bfl, id = plan$ids[2L]), infos[[2L]])

  # Frames without the new cache still recover the same structural plan.
  bfl$frame$re_s2z_plan <- NULL
  expect_identical(re_s2z_plan(bfl), plan)
})

test_that("active priors are resolved fresh without caching fixed-only priors", {
  form <- y ~ x + w + z +
    (1 | gr(g, s2z = TRUE)) + (0 + z | gr(h, s2z = TRUE))
  bfl <- brmsframe(brmsterms(form), data = s2z_plan_dat)$dpars$mu
  structural <- re_s2z_plan(bfl)
  p1 <- validate_prior(
    prior(normal(0, 1), class = Intercept) +
      prior(normal(0, 2), class = b, coef = z) +
      prior(double_exponential(0, 1), class = b, coef = x),
    formula = form, data = s2z_plan_dat
  )
  p2 <- p1
  p2$prior[p2$class == "b" & p2$coef == "z"] <- "logistic(1, 3)"
  plan1 <- re_s2z_plan(bfl, prior = p1)
  plan2 <- re_s2z_plan(bfl, prior = p2)

  expect_named(plan1$prior, c("Intercept", "z"))
  expect_identical(plan1$prior$z$dist, "normal")
  expect_identical(plan2$prior$z$dist, "logistic")
  expect_equal(plan2$prior$z$location, 1)
  expect_equal(plan2$prior$z$scale, 3)
  expect_identical(re_s2z_plan(bfl), structural)
  expect_null(bfl$frame$re_s2z_plan$prior)
  expect_true(all(vapply(plan2$blocks, function(block) {
    is.null(block$prior)
  }, logical(1))))

  # Existing kernels see the same original-length vector and neutral
  # placeholders; unsupported fixed-only priors do not enter the solver.
  legacy <- re_s2z_plan_infos(plan2)[[1L]]$prior
  expect_identical(vapply(legacy, `[[`, character(1), "dist"),
                   c("normal", "flat", "flat", "logistic"))
})

test_that("structural centering maps do not depend on realized column means", {
  form <- y ~ x + w + (0 + x | gr(g, s2z = TRUE))
  bfl <- brmsframe(brmsterms(form), data = s2z_plan_dat)$dpars$mu
  plan <- re_s2z_plan(bfl)
  expect_equal(mean(s2z_plan_dat$x), 0)
  expect_identical(plan$active_q, 1:2)
  expect_identical(plan$blocks[[1L]]$match_q, 2L)
  expect_identical(plan$blocks[[1L]]$block_active_q, 1:2)

  uncentered <- bf(form, center = FALSE)
  bfl <- brmsframe(brmsterms(uncentered), data = s2z_plan_dat)$dpars$mu
  plan <- re_s2z_plan(bfl)
  expect_false(plan$center)
  expect_identical(plan$active_names, "x")
  expect_identical(plan$inactive_names, c("Intercept", "w"))
})

test_that("ordinary predictors do not build or cache S2Z plans", {
  for (form in list(
    y ~ x,
    y ~ x + (1 | g),
    y ~ x + (1 | gr(g, s2z = FALSE))
  )) {
    bframe <- brmsframe(brmsterms(form), data = s2z_plan_dat)
    bfl <- bframe$dpars$mu
    expect_null(bframe$frame$s2z_context)
    expect_null(bfl$frame$re_s2z_plan)
    expect_null(re_s2z_plan(bfl))
    expect_identical(re_s2z_infos(bfl), list())
  }
})
