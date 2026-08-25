context("Phased validation of physical sum-to-zero group-level effects")

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
  expect_equal(
    parse_re_s2z_prior("logistic(0, 1)"),
    list(dist = "logistic", location = 0, scale = 1, df = NA_real_)
  )
  expect_error(
    parse_re_s2z_prior("logistic(0, 0)"),
    "Scale and degrees-of-freedom arguments must be positive"
  )
  expect_error(
    parse_re_s2z_prior("logistic(location, 1)"),
    "must currently be numeric constants"
  )
})

test_that("S2Z descriptors separate active and fixed-only coordinates", {
  form <- y ~ x + z + w +
    (1 + x | gr(g, id = "first", s2z = TRUE)) +
    (0 + z | gr(h, id = "second", s2z = TRUE))
  bframe <- brmsframe(brmsterms(form), data = s2z_validation_dat)
  bfl <- bframe$dpars$mu
  infos <- re_s2z_infos(bfl)

  expect_length(infos, 2L)
  for (info in infos) {
    expect_equal(info$qnames, c("Intercept", "x", "z", "w"))
    expect_equal(info$active_q, 1:3)
    expect_equal(info$inactive_q, 4L)
    expect_equal(info$active_names, c("Intercept", "x", "z"))
    expect_equal(info$inactive_names, "w")
    expect_equal(info$active_index, c(1L, 2L, 3L, NA_integer_))
    expect_equal(info$match_active, match(info$match_q, info$active_q))
  }

  bprior <- prior(normal(0, 2), class = Intercept) +
    prior(normal(0, 1), class = b, coef = x) +
    prior(normal(0, 1.5), class = b, coef = z) +
    prior(double_exponential(0, 1), class = b, coef = w)
  effective_prior <- validate_prior(
    bprior, formula = form, data = s2z_validation_dat
  )
  expect_no_error(
    validate_re_s2z_prior_global(bframe, prior = effective_prior)
  )
  specs <- re_s2z_infos(bfl, prior = effective_prior)[[1L]]$prior
  expect_equal(vapply(specs, `[[`, character(1), "dist"),
               c("normal", "normal", "normal", "flat"))
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

test_that("S2Z design matching is numerically exact", {
  form <- y ~ x + (1 + x | gr(g, id = "exact", s2z = TRUE))
  bframe <- brmsframe(brmsterms(form), data = s2z_validation_dat)$dpars$mu
  expect_no_error(
    validate_re_s2z_design(bframe, data = s2z_validation_dat)
  )

  bframe$sdata$fe$X[, 2L] <- bframe$sdata$fe$X[, 2L] + 1e-12
  msg <- s2z_error_message(
    validate_re_s2z_design(bframe, data = s2z_validation_dat)
  )
  expect_match(msg, "S2Z capability 'matching_values'", fixed = TRUE)
  expect_match(msg, "coefficient(s) 'x'", fixed = TRUE)
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
  expect_match(msg, '`id = "y_mu_s2z"`', fixed = TRUE)
  expect_match(msg, '`id = "y_sigma_s2z"`', fixed = TRUE)
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
  expect_match(msg, '`id = "y_mu_s2z"`', fixed = TRUE)
  expect_match(msg, '`id = "z_mu_s2z"`', fixed = TRUE)
  expect_match(msg, "For `mvbind(...)` or category shorthand", fixed = TRUE)
  expect_match(msg, "omit the shared `| shared |` tag", fixed = TRUE)
  expect_match(msg, "separate `bf()` response formulas", fixed = TRUE)
  expect_match(msg, "do not retain cross-predictor group-effect", fixed = TRUE)
  expect_match(msg, "use s2z = FALSE", fixed = TRUE)
})

test_that("mvbind shorthand without a shared ID uses local S2Z blocks", {
  form <- bf(
    mvbind(y, z) ~ 1 + (1 | gr(g, s2z = TRUE))
  ) + set_rescor(FALSE)
  code <- stancode(form, data = s2z_validation_dat)
  bframe <- brmsframe(brmsterms(form), data = s2z_validation_dat)
  local_re <- lapply(bframe$terms, function(x) x$frame$re)

  expect_s3_class(code, "brmsmodel")
  expect_true(all(vapply(local_re, function(x) any(x$s2z), logical(1))))
  expect_length(unique(unlist(lapply(local_re, `[[`, "id"))), 2L)
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
  expect_match(msg, "use s2z = FALSE", fixed = TRUE)
  expect_match(msg, "supplied group covariance matrix", fixed = TRUE)
  expect_false(grepl("Arguments 'by', 'cov', and 'pw'", msg, fixed = TRUE))
})
