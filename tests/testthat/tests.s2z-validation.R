context("Validation of physical sum-to-zero group-level effects")

expect_match2 <- brms:::expect_match2

s2z_validation_dat <- local({
  n <- 36L
  data.frame(
    y = sin(seq_len(n) / 5),
    y_ord = ordered(rep(letters[1:3], length.out = n)),
    x = seq(-1.5, 1.5, length.out = n),
    z = cos(seq_len(n) / 7),
    w = seq(0.2, 2.1, length.out = n),
    sx = rep(0.1, n),
    g = factor(rep(seq_len(6L), each = 6L)),
    h = factor(rep(seq_len(4L), each = 9L))
  )
})

s2z_error_message <- function(expr) {
  error <- tryCatch(expr, error = identity)
  expect_s3_class(error, "error")
  conditionMessage(error)
}

test_that("S2Z logistic prior specifications are exact and validated", {
  form <- y ~ 1 + (1 | gr(g, s2z = TRUE))
  code <- stancode(
    form, data = s2z_validation_dat,
    prior = prior(logistic(0, 1), class = Intercept)
  )
  expect_match2(code, "logistic_lpdf(q_explicit_s2z_1[1] | 0, 1)")
  expect_error(
    stancode(form, data = s2z_validation_dat,
             prior = prior(logistic(0, 0), class = Intercept)),
    "Scale and degrees-of-freedom arguments must be positive"
  )
  expect_error(
    stancode(form, data = s2z_validation_dat,
             prior = prior(logistic(location, 1), class = Intercept)),
    "must currently be numeric constants"
  )
})

test_that("S2Z preserves original coordinates across blocks and prior changes", {
  form <- y ~ x + w + z +
    (1 | gr(g, s2z = TRUE)) + (0 + z | gr(h, s2z = TRUE))
  bprior <- prior(normal(0, 1), class = Intercept) +
    prior(double_exponential(0, 1), class = b, coef = x)
  normal_code <- stancode(
    form, data = s2z_validation_dat,
    prior = bprior + prior(normal(0, 2), class = b, coef = z)
  )
  logistic_code <- stancode(
    form, data = s2z_validation_dat,
    prior = bprior + prior(logistic(1, 3), class = b, coef = z)
  )

  for (code in list(normal_code, logistic_code)) {
    # The two fixed-only slopes keep their positions when the group maps
    # involve nonadjacent population coefficients.
    expect_match2(code, "theta_s2z[1] = theta_s2z_active[1];")
    expect_match2(code, "theta_s2z[4] = theta_s2z_active[2];")
    expect_match2(code, "theta_s2z[2] = fixed_s2z[1];")
    expect_match2(code, "theta_s2z[3] = fixed_s2z[2];")
    expect_match2(code, "double_exponential_lpdf(fixed_s2z[1] | 0, 1)")
  }
  expect_false(grepl("q_explicit_s2z", normal_code, fixed = TRUE))
  expect_match2(logistic_code, "logistic_lpdf(q_explicit_s2z_1[4] | 1, 3)")
})

test_that("S2Z centering maps include intercepts even for mean-zero slopes", {
  form <- y ~ x + w + (0 + x | gr(g, s2z = TRUE))
  dat <- transform(s2z_validation_dat, x = rep(c(-1, 1), 18))
  code <- stancode(form, data = dat)
  expect_match2(code, "H_s2z_1[1] = means_X[1];")
  expect_match2(code, "theta_s2z[1] = theta_s2z_active[1];")
  expect_match2(code, "theta_s2z[2] = theta_s2z_active[2];")

  code <- stancode(bf(form, center = FALSE), data = dat)
  expect_match2(code, "theta_s2z[1] = fixed_s2z[1];")
  expect_match2(code, "theta_s2z[2] = theta_s2z_active[1];")
  expect_match2(code, "theta_s2z[3] = fixed_s2z[2];")
})

test_that("S2Z coefficient priors retain global bounds", {
  form <- y ~ x + z + (1 + x | gr(g, s2z = TRUE))
  bprior <- prior(normal(0, 1), class = b, lb = 0) +
    prior(normal(0, 2), class = b, coef = x)
  # Overriding a coefficient's density does not remove its class-level bound.
  expect_error(
    validate_prior(bprior, formula = form, data = s2z_validation_dat),
    "S2Z capability 'active_prior_bounds'", fixed = TRUE
  )
  expect_error(
    stancode(form, data = s2z_validation_dat, prior = bprior),
    "S2Z capability 'active_prior_bounds'", fixed = TRUE
  )
})

test_that("multi-block prior diagnostics identify a touching S2Z block", {
  form <- y ~ x + z +
    (1 + x | gr(g, id = "gblock", s2z = TRUE)) +
    (0 + z | gr(h, id = "hblock", s2z = TRUE))
  msg <- s2z_error_message(stancode(
    form, data = s2z_validation_dat,
    prior = prior(double_exponential(0, 1), class = b, coef = z)
  ))

  expect_match(msg, "group 'h'", fixed = TRUE)
  expect_match(msg, "ID 'hblock'", fixed = TRUE)
  expect_match(msg, "coefficient(s) 'z'", fixed = TRUE)
  expect_false(grepl("group 'g'", msg, fixed = TRUE))

  special_form <- y ~ z +
    (1 | gr(g, id = "intercept-block", s2z = TRUE)) +
    (0 + z | gr(h, id = "slope-block", s2z = TRUE))
  msg <- s2z_error_message(validate_prior(
    prior(horseshoe(), class = b),
    formula = special_form, data = s2z_validation_dat
  ))
  expect_match(msg, "S2Z capability 'active_special_prior'", fixed = TRUE)
  expect_match(msg, "group 'h'", fixed = TRUE)
  expect_match(msg, "ID 'slope-block'", fixed = TRUE)
  expect_match(msg, "coefficient(s) 'z'", fixed = TRUE)
})

test_that("S2Z structure errors precede design errors and carry context", {
  ordinal_form <- y_ord ~ x +
    (1 | gr(g, id = "score", s2z = TRUE))
  msg <- s2z_error_message(stancode(
    ordinal_form, data = s2z_validation_dat, family = cumulative()
  ))
  expect_match(msg, "S2Z capability 'ordinal_location'", fixed = TRUE)
  expect_match(msg, "response 'y_ord'", fixed = TRUE)
  expect_match(msg, "family 'cumulative'", fixed = TRUE)
  expect_match(msg, "dpar 'mu'", fixed = TRUE)
  expect_match(msg, "group 'g'", fixed = TRUE)
  expect_match(msg, "ID 'score'", fixed = TRUE)
  expect_match(msg, "Remedy:", fixed = TRUE)
  expect_false(grepl("matching population-level", msg, fixed = TRUE))

  special_form <- y ~ x +
    (1 + me(x, sx) | gr(g, id = "me-score", s2z = TRUE))
  msg <- s2z_error_message(stancode(
    special_form, data = s2z_validation_dat
  ))
  expect_match(msg, "S2Z capability 'ordinary_gr_only'", fixed = TRUE)
  expect_match(msg, "ID 'me-score'", fixed = TRUE)
  expect_false(grepl("matching population-level", msg, fixed = TRUE))
})

test_that("the VerbAgg ordinal model reaches the intentional S2Z gate", {
  skip_if_not_installed("lme4")
  data_env <- new.env(parent = emptyenv())
  utils::data("VerbAgg", package = "lme4", envir = data_env)
  verbagg <- data_env$VerbAgg
  form <- resp ~ (Anger + Gender + btype + situ)^2 +
    (1 | gr(id, id = "person-s2z", s2z = TRUE)) +
    (1 | gr(item, id = "item-s2z", s2z = TRUE))

  msg <- s2z_error_message(stancode(
    form, data = verbagg, family = cumulative()
  ))
  expect_match(msg, "S2Z capability 'ordinal_location'", fixed = TRUE)
  expect_match(msg, "response 'resp'", fixed = TRUE)
  expect_match(msg, "family 'cumulative'", fixed = TRUE)
  expect_match(msg, "group 'id'", fixed = TRUE)
  expect_match(msg, "ID 'person-s2z'", fixed = TRUE)
  expect_match(msg, "ordinal S2Z support", fixed = TRUE)
  expect_false(grepl("matching population-level", msg, fixed = TRUE))
})

test_that("S2Z design and prior/global phases use capability diagnostics", {
  missing_form <- y ~ x +
    (1 + x + z | gr(g, id = "score", s2z = TRUE))
  msg <- s2z_error_message(stancode(
    missing_form, data = s2z_validation_dat
  ))
  expect_match(msg, "S2Z capability 'matching_name'", fixed = TRUE)
  expect_match(msg, "coefficient(s) 'z'", fixed = TRUE)
  expect_match(msg, "Remedy:", fixed = TRUE)

  active_prior_form <- y ~ x + w +
    (1 + x | gr(g, id = "score", s2z = TRUE))
  msg <- s2z_error_message(stancode(
    active_prior_form, data = s2z_validation_dat,
    prior = prior(double_exponential(0, 1), class = b, coef = x)
  ))
  expect_match(
    msg, "S2Z capability 'active_prior_distribution'", fixed = TRUE
  )
  expect_match(msg, "prior 'double_exponential(0,1)'", fixed = TRUE)
  expect_match(msg, "Remedy:", fixed = TRUE)

  msg <- s2z_error_message(validate_prior(
    prior(double_exponential(0, 1), class = b, coef = x),
    formula = active_prior_form, data = s2z_validation_dat
  ))
  expect_match(
    msg, "S2Z capability 'active_prior_distribution'", fixed = TRUE
  )
  expect_match(msg, "coefficient(s) 'x'", fixed = TRUE)

  cross_form <- bf(
    y ~ 1 + (1 | gr(g, id = "across", s2z = TRUE)),
    sigma ~ 1 + (1 | gr(g, id = "across", s2z = FALSE))
  )
  msg <- s2z_error_message(stancode(
    cross_form, data = s2z_validation_dat
  ))
  expect_match(msg, "S2Z capability 'cross_predictor_id'", fixed = TRUE)
  expect_match(msg, "ID 'across'", fixed = TRUE)
  expect_match(msg, "Affected linear predictors:", fixed = TRUE)
  expect_match(msg, "dpar 'mu'", fixed = TRUE)
  expect_match(msg, "dpar 'sigma'", fixed = TRUE)
  expect_match(msg, "Remedy:", fixed = TRUE)
})

test_that("cross-predictor diagnostics enumerate mvbind responses", {
  form <- bf(
    mvbind(y, z) ~ 1 + (1 | shared | gr(g, s2z = TRUE))
  ) + set_rescor(FALSE)
  msg <- s2z_error_message(stancode(
    form, data = s2z_validation_dat
  ))

  expect_match(msg, "Affected linear predictors:", fixed = TRUE)
  expect_match(msg, "response 'y'", fixed = TRUE)
  expect_match(msg, "response 'z'", fixed = TRUE)
  expect_match(msg, "ID 'shared'", fixed = TRUE)
  expect_match(msg, "s2z = FALSE", fixed = TRUE)
})

test_that("mvbind shorthand without a shared ID uses local S2Z blocks", {
  form <- bf(
    mvbind(y, z) ~ 1 + (1 | gr(g, s2z = TRUE))
  ) + set_rescor(FALSE)
  code <- stancode(form, data = s2z_validation_dat)
  sdata <- standata(form, data = s2z_validation_dat)

  expect_match2(code, "theta_s2z_y")
  expect_match2(code, "theta_s2z_z")
  expect_match2(code, "q_recovered_s2z_1")
  expect_match2(code, "q_recovered_s2z_2")
  expect_equal(sdata$N_1, nlevels(s2z_validation_dat$g))
  expect_equal(sdata$N_2, nlevels(s2z_validation_dat$g))
})

test_that("cov diagnostics identify the unsupported phylogenetic argument", {
  phylo_dat <- data.frame(
    phen = seq(-1, 1, length.out = 12),
    cofactor = rep(c(-0.5, 0.5), 6),
    phylo = factor(rep(letters[1:4], each = 3))
  )
  A <- diag(4)
  dimnames(A) <- list(levels(phylo_dat$phylo), levels(phylo_dat$phylo))
  msg <- s2z_error_message(stancode(
    phen ~ cofactor + (1 | gr(phylo, cov = A, s2z = TRUE)),
    data = phylo_dat, data2 = list(A = A)
  ))

  expect_match(msg, "Argument 'cov' is not yet supported", fixed = TRUE)
  expect_match(msg, "S2Z capability 'cov'", fixed = TRUE)
  expect_match(msg, "response 'phen'", fixed = TRUE)
  expect_match(msg, "group 'phylo'", fixed = TRUE)
  expect_match(msg, "s2z = FALSE", fixed = TRUE)
  expect_match(msg, "supplied group covariance matrix", fixed = TRUE)
  expect_false(grepl("Arguments 'by', 'cov', and 'pw'", msg, fixed = TRUE))
})
