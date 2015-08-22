test_that("Test that brm.pars returns correct parameter names", {
  expect_equal(brm.pars(rating ~ treat + period + carry + (1|subject), data = inhaler),
               c("b", "sigma", "sd_subject", "r_subject"))
  expect_equal(brm.pars(rating ~ treat + period + carry + (1+treat|subject), data = inhaler),
               c("b", "sigma", "sd_subject", "cor_subject", "r_subject"))
  expect_equal(brm.pars(rating ~ treat + period + carry + (1+treat|subject), data = inhaler, autocor = cor.ma()),
               c("b", "sigma", "ma", "sd_subject", "cor_subject", "r_subject"))
  expect_equal(brm.pars(rating ~ treat + period + carry + (1+treat|subject), data = inhaler, 
                        autocor = cor.ma(), ranef = FALSE), c("b", "sigma", "ma", "sd_subject", "cor_subject"))
})