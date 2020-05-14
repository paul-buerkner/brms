context("Tests for emmeans support")

set.seed(2020.0515)
fit3 <- brm(time | cens(censored) ~ age * sex + disease + (1|patient),
            data = kidney, family = lognormal())

rg <- ref_grid(fit3)
smry <- summary(emmeans(rg, "disease"), point.est = mean)

test_that("We can get basic emmeans estimates", {
    expect_equal(rg@levels$age, 43.697, tolerance = 0.002)
    expect_equal(rg@levels$sex, c("male", "female"))
    expect_equal(smry$emmean, c(4.12, 3.68, 3.58, 4.70), tolerance = 0.02)
})

