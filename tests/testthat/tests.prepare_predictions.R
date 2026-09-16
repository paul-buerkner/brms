context("Tests for prediction preparation")

test_that("reduced prediction formulas retain fitted group distributions", {
  dat <- data.frame(y = seq_len(6), g = factor(rep(c("a", "b"), each = 3)))
  draws <- posterior::as_draws_matrix(cbind(
    `r_g[a,Intercept]` = rep(0.1, 4),
    `r_g[b,Intercept]` = rep(-0.1, 4),
    sd_g__Intercept = rep(1, 4)
  ))
  sdata <- list(J_1 = 3L)
  brms:::set_levels(sdata, "used") <- list(g = c("a", "b", "new"))
  forms <- list(
    student = bf(y ~ 1 + (1 | gr(g, dist = "student"))),
    gaussian = bf(y ~ 1 + (1 | gr(g)))
  )
  for (name in names(forms)) {
    fitted <- brms:::brmsframe(brmsterms(forms[[name]]), dat)
    reduced <- brms:::brmsframe(brmsterms(
      brms:::update_re_terms(forms[[name]], ~(1 | g))
    ), dat)
    args <- list(bframe = reduced, draws = draws, sdata = sdata,
                 old_reframe = fitted$frame$re)
    if (name == "student") {
      expect_error(do.call(brms:::prepare_predictions_re_global,
                           c(args, sample_new_levels = "gaussian")),
                   "not available for non-gaussian group-level effects")
    } else {
      out <- do.call(brms:::prepare_predictions_re_global,
                     c(args, sample_new_levels = "gaussian"))
      expect_equal(dim(out$g$rdraws), c(4L, 1L))
      expect_true(all(is.finite(out$g$rdraws)))
    }
    out <- do.call(brms:::prepare_predictions_re_global,
                   c(args, sample_new_levels = "uncertainty"))
    expect_equal(dim(out$g$rdraws), c(4L, 1L))
    expect_true(all(is.finite(out$g$rdraws)))
  }
})
