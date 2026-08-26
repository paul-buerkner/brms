context("Tests for cross-predictor S2Z data")

test_that("cross-predictor partial centering follows covariance columns", {
  dat <- data.frame(
    phen = seq(-1, 1, length.out = 12),
    cofactor = cos(seq_len(12)),
    phylo = factor(rep(c("c", "a", "b"), 4))
  )
  form <-
    bf(
      phen ~ 1 +
        (1 | q | gr(phylo, s2z = TRUE, center = 0.2))
    ) +
    bf(
      cofactor ~ 1 +
        (1 | q | gr(phylo, s2z = TRUE, center = 0.8))
    ) +
    set_rescor(FALSE)

  sdata <- standata(form, data = dat)
  expect_identical(sum(names(sdata) == "rho_s2z_1"), 1L)
  expect_equal(dim(sdata$rho_s2z_1), c(3L, 2L))
  expect_equal(
    unname(sdata$rho_s2z_1),
    cbind(rep(0.2, 3L), rep(0.8, 3L))
  )
  expect_identical(
    colnames(sdata$rho_s2z_1),
    c("phen:Intercept[1]", "cofactor:Intercept[2]")
  )

  scode <- stancode(form, data = dat)
  expect_identical(
    lengths(regmatches(
      scode,
      gregexpr(
        "matrix<lower=0,upper=1>[N_1, M_1] rho_s2z_1;",
        scode, fixed = TRUE
      )
    )),
    1L
  )
})
