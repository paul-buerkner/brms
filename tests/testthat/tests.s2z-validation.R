context("Phased validation of physical sum-to-zero group-level effects")

s2z_validation_dat <- local({
  n <- 36L
  data.frame(
    y = sin(seq_len(n) / 5),
    y_ord = ordered(rep(letters[1:3], length.out = n)),
    y_ord4 = ordered(rep(letters[1:4], length.out = n)),
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

s2z_ordinal_frame <- function(formula, family, data = s2z_validation_dat) {
  formula <- validate_formula(formula, family = family, data = data)
  bterms <- brmsterms(formula)
  data <- validate_data(data, bterms = bterms)
  brmsframe(bterms, data = data)
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
  special_form <- y ~ x +
    (1 + me(x, sx) | gr(g, id = "me-score", s2z = TRUE))
  msg <- s2z_error_message(stancode(
    special_form, data = s2z_validation_dat
  ))
  expect_match(msg, "S2Z capability 'ordinary_gr_only'", fixed = TRUE)
  expect_match(msg, "ID 'me-score'", fixed = TRUE)
  expect_false(grepl("matching population-level", msg, fixed = TRUE))
})

test_that("the VerbAgg ordinal model reaches threshold-aware descriptors", {
  skip_if_not_installed("lme4")
  data_env <- new.env(parent = emptyenv())
  utils::data("VerbAgg", package = "lme4", envir = data_env)
  verbagg <- data_env$VerbAgg
  form <- resp ~ (Anger + Gender + btype + situ)^2 +
    (1 | gr(id, id = "person-s2z", s2z = TRUE)) +
    (1 | gr(item, id = "item-s2z", s2z = TRUE))

  bframe <- s2z_ordinal_frame(form, cumulative(), data = verbagg)$dpars$mu
  infos <- re_s2z_infos(bframe)
  expect_length(infos, 2L)
  expect_true(all(vapply(infos, `[[`, logical(1), "ordinal")))
  expect_true(all(vapply(infos, function(x) {
    identical(x$r$coef, "Intercept") &&
      all(x$H[x$threshold_q, 1L] == -1) &&
      all(x$H[x$slope_q, 1L] == 0)
  }, logical(1))))
})

test_that("ordinal q coordinates put threshold primitives before slopes", {
  form <- y_ord ~ w + z +
    (1 + w | gr(g, id = "ordinal-map", s2z = TRUE))
  bframe <- s2z_ordinal_frame(form, cumulative())$dpars$mu
  info <- re_s2z_info(bframe)

  expect_true(info$ordinal)
  expect_equal(info$qtype, c("threshold", "threshold", "b", "b"))
  expect_equal(info$threshold_q, 1:2)
  expect_equal(info$slope_q, 3:4)
  expect_equal(info$qnames, c("Intercept[1]", "Intercept[2]", "w", "z"))
  expect_equal(
    names(info$threshold),
    c(
      "q", "group", "group_index", "coef", "threshold_index", "kind",
      "source", "qname"
    )
  )
  expect_equal(info$threshold$coef, c("1", "2"))
  expect_equal(info$threshold$kind, rep("flexible", 2L))
  expect_equal(info$match_q, c(NA_integer_, 3L))

  expect_equal(info$C, matrix(
    c(0, 1, 0, 0), nrow = 2L, byrow = TRUE,
    dimnames = list(c("w", "z"), c("Intercept", "w"))
  ))
  expect_equal(info$a, c(1, mean(s2z_validation_dat$w)))
  expect_equal(
    unname(info$H[info$threshold_q, , drop = FALSE]),
    matrix(rep(-info$a, 2L), nrow = 2L, byrow = TRUE)
  )
  expect_equal(
    unname(info$H[info$slope_q, , drop = FALSE]), unname(info$C)
  )
  expect_equal(info$active_q, 1:3)
  expect_equal(info$inactive_q, 4L)

  Z <- re_s2z_design_matrix(bframe, s2z_validation_dat)
  Xc <- re_s2z_likelihood_design(bframe)
  expect_equal(
    unname(Z),
    unname(matrix(1, nrow(Z), 1L) %*% t(info$a) + Xc %*% info$C),
    tolerance = info$affine_tolerance
  )
})

test_that("ordinal varying intercepts do not require a fixed Intercept", {
  form <- y_ord ~ 1 + (1 | gr(g, id = "threshold-only", s2z = TRUE))
  bframe <- s2z_ordinal_frame(form, cumulative())$dpars$mu
  info <- re_s2z_info(bframe)

  expect_false("Intercept" %in% info$fixef)
  expect_equal(info$qnames, c("Intercept[1]", "Intercept[2]"))
  expect_equal(dim(info$C), c(0L, 1L))
  expect_equal(info$a, 1)
  expect_equal(unname(info$H), matrix(-1, nrow = 2L, ncol = 1L))
  expect_true(is.na(info$match_q))
  expect_equal(info$active_q, info$threshold_q)
})

test_that("ordinal slope-only blocks carry centered means into thresholds", {
  form <- y_ord ~ w + z +
    (0 + w | gr(g, id = "slope-only", s2z = TRUE))
  bframe <- s2z_ordinal_frame(form, sratio())$dpars$mu
  info <- re_s2z_info(bframe)

  expect_equal(info$r$coef, "w")
  expect_equal(info$a, mean(s2z_validation_dat$w))
  expect_equal(unname(info$H[info$threshold_q, 1L]), rep(-info$a, 2L))
  expect_equal(unname(info$H[info$slope_q, 1L]), c(1, 0))
})

test_that("equidistant and grouped ordinal primitives are explicit", {
  equidistant_form <- y_ord ~ w +
    (1 + w | gr(g, id = "equidistant", s2z = TRUE))
  equidistant_frame <- s2z_ordinal_frame(
    equidistant_form, cumulative(threshold = "equidistant")
  )$dpars$mu
  equidistant <- re_s2z_info(equidistant_frame)
  expect_equal(equidistant$threshold_q, 1L)
  expect_equal(equidistant$threshold$kind, "first")
  expect_equal(equidistant$threshold$coef, "")
  expect_equal(equidistant$qnames, c("first_Intercept", "w"))
  expect_equal(unname(equidistant$H[1L, ]), -equidistant$a)

  grouped_form <- y_ord | thres(gr = h) ~ w +
    (1 + w | gr(g, id = "grouped", s2z = TRUE))
  grouped_frame <- s2z_ordinal_frame(
    grouped_form, cumulative()
  )$dpars$mu
  grouped <- re_s2z_info(grouped_frame)
  expect_equal(grouped$threshold_q, seq_len(8L))
  expect_equal(grouped$threshold$group_index, rep(seq_len(4L), each = 2L))
  expect_equal(grouped$threshold$threshold_index, rep(1:2, 4L))
  expect_equal(
    unname(grouped$H[grouped$threshold_q, , drop = FALSE]),
    matrix(rep(-grouped$a, 8L), nrow = 8L, byrow = TRUE)
  )
})

test_that("ordinal design errors distinguish names from affine identity", {
  missing_form <- y_ord ~ w +
    (1 + w + z | gr(g, id = "missing", s2z = TRUE))
  msg <- s2z_error_message(s2z_ordinal_frame(missing_form, cumulative()))
  expect_match(msg, "S2Z capability 'matching_name'", fixed = TRUE)
  expect_match(msg, "coefficient(s) 'z'", fixed = TRUE)
  expect_false(grepl("coefficient(s) 'Intercept'", msg, fixed = TRUE))

  form <- y_ord ~ w +
    (1 + w | gr(g, id = "affine", s2z = TRUE))
  bframe <- s2z_ordinal_frame(form, cumulative())$dpars$mu
  changed <- s2z_validation_dat
  changed$w[1L] <- changed$w[1L] + 1e-10
  msg <- s2z_error_message(validate_re_s2z_design(bframe, data = changed))
  expect_match(msg, "S2Z capability 'affine_identity'", fixed = TRUE)
  expect_match(msg, "coefficient(s) 'w'", fixed = TRUE)
  expect_match(msg, "Z = 1 a^T + Xc C", fixed = TRUE)
})

test_that("ordinal affine validation accepts exact general maps", {
  exact_form <- y_ord ~ I(2 * w + 3) +
    (0 + I(2 * w + 3) | gr(g, id = "exact-affine", s2z = TRUE))
  bframe <- s2z_ordinal_frame(exact_form, cumulative())$dpars$mu
  info <- re_s2z_info(bframe)

  expect_no_error(
    validate_re_s2z_design(bframe, data = s2z_validation_dat)
  )
  expect_true(info$affine_ok)
  expect_equal(unname(info$C), matrix(1, nrow = 1L, ncol = 1L))
  expect_equal(info$a, mean(2 * s2z_validation_dat$w + 3))

  renamed_form <- y_ord ~ I(2 * w + 3) +
    (0 + w | gr(g, id = "renamed-affine", s2z = TRUE))
  renamed_frame <- s2z_ordinal_frame(
    renamed_form, cumulative()
  )$dpars$mu
  renamed <- re_s2z_info(renamed_frame)
  expect_no_error(
    validate_re_s2z_design(renamed_frame, data = s2z_validation_dat)
  )
  expect_true(renamed$affine_ok)
  expect_true(is.na(renamed$match_slope))
  expect_equal(unname(renamed$C), matrix(0.5, nrow = 1L, ncol = 1L))
  expect_equal(renamed$a, mean(s2z_validation_dat$w))
})

test_that("ordinal threshold priors retain per-primitive metadata", {
  form <- y_ord4 ~ w +
    (1 + w | gr(g, id = "threshold-priors", s2z = TRUE))
  bframe <- s2z_ordinal_frame(form, cumulative())
  bprior <- prior(normal(0, 2), class = Intercept, coef = "1") +
    prior(student_t(5, 1, 3), class = Intercept, coef = "2") +
    prior(cauchy(-1, 4), class = Intercept, coef = "3")
  effective_prior <- validate_prior(
    bprior, formula = form, data = s2z_validation_dat,
    family = cumulative()
  )
  specs <- re_s2z_infos(
    bframe$dpars$mu, prior = effective_prior
  )[[1L]]$threshold_prior

  expect_equal(vapply(specs, `[[`, character(1), "dist"),
               c("normal", "student", "student"))
  expect_equal(vapply(specs, `[[`, character(1), "coef"), c("1", "2", "3"))
  expect_equal(vapply(specs, `[[`, character(1), "group"), rep("", 3L))
  expect_equal(vapply(specs, `[[`, numeric(1), "q"), 1:3)
  expect_equal(vapply(specs, `[[`, numeric(1), "df"), c(NA, 5, 1))
})

test_that("ordinal active logistic priors are rejected contextually", {
  form <- y_ord ~ w +
    (1 + w | gr(g, id = "ordinal-logistic", s2z = TRUE))
  msg <- s2z_error_message(stancode(
    form, data = s2z_validation_dat, family = cumulative(),
    prior = prior(logistic(0, 1), class = Intercept)
  ))
  expect_match(
    msg, "S2Z capability 'ordinal_active_prior_distribution'", fixed = TRUE
  )
  expect_match(msg, "threshold", fixed = TRUE)
  expect_match(msg, "prior 'logistic(0, 1)'", fixed = TRUE)

  msg <- s2z_error_message(stancode(
    form, data = s2z_validation_dat, family = cumulative(),
    prior = prior(logistic(0, 1), class = b, coef = w)
  ))
  expect_match(
    msg, "S2Z capability 'ordinal_active_prior_distribution'", fixed = TRUE
  )
  expect_match(msg, "coefficient(s) 'w'", fixed = TRUE)
})

test_that("bounded and tagged ordinal threshold priors stay gated", {
  form <- y_ord ~ w +
    (1 + w | gr(g, id = "threshold-prior-gate", s2z = TRUE))
  cases <- list(
    ordinal_threshold_prior_bounds = prior(
      normal(0, 1), class = Intercept, lb = -5
    ),
    ordinal_threshold_prior_tag = prior(
      normal(0, 1), class = Intercept, tag = "thresholds"
    )
  )

  for (capability in names(cases)) {
    msg <- s2z_error_message(stancode(
      form, data = s2z_validation_dat, family = cumulative(),
      prior = cases[[capability]]
    ))
    expect_match(
      msg, paste0("S2Z capability '", capability, "'"), fixed = TRUE
    )
    expect_match(msg, "ID 'threshold-prior-gate'", fixed = TRUE)
    expect_match(msg, "coefficient(s) 'threshold 1'", fixed = TRUE)
    expect_match(msg, "Remedy:", fixed = TRUE)
  }
})

test_that("unsupported ordinal threshold and category-specific maps stay gated", {
  sum_to_zero_form <- y_ord ~ w +
    (1 + w | gr(g, id = "threshold-stz", s2z = TRUE))
  msg <- s2z_error_message(s2z_ordinal_frame(
    sum_to_zero_form, cumulative(threshold = "sum_to_zero")
  ))
  expect_match(
    msg, "S2Z capability 'ordinal_sum_to_zero_thresholds'", fixed = TRUE
  )
  expect_match(msg, "common location coordinate", fixed = TRUE)

  cs_form <- y_ord ~ w +
    (cs(1) | gr(g, id = "category-specific", s2z = TRUE))
  msg <- s2z_error_message(s2z_ordinal_frame(cs_form, cumulative()))
  expect_match(
    msg, "S2Z capability 'ordinal_category_specific'", fixed = TRUE
  )
  expect_match(msg, "ID 'category-specific'", fixed = TRUE)
})

test_that("fixed mixture and custom ordinal thresholds stay gated", {
  mixture_form <- bf(
    y_ord ~ 1,
    mu1 ~ w +
      (1 + w | gr(g, id = "shared-thresholds", s2z = TRUE)),
    mu2 ~ 1
  )
  msg <- s2z_error_message(s2z_ordinal_frame(
    mixture_form, mixture(cumulative(), nmix = 2, order = "mu")
  ))
  expect_match(
    msg, "S2Z capability 'ordinal_fixed_thresholds'", fixed = TRUE
  )
  expect_match(msg, "Fixed or shared ordinal-mixture thresholds", fixed = TRUE)
  expect_match(msg, "ID 'shared-thresholds'", fixed = TRUE)
  expect_match(msg, "Remedy:", fixed = TRUE)

  custom_ordinal <- custom_family(
    "s2z_custom_ordinal", dpars = "mu", links = "identity",
    type = "int", specials = "ordinal", threshold = "flexible"
  )
  custom_form <- y_ord ~ w +
    (1 + w | gr(g, id = "custom-ordinal", s2z = TRUE))
  msg <- s2z_error_message(
    s2z_ordinal_frame(custom_form, custom_ordinal)
  )
  expect_match(msg, "S2Z capability 'ordinal_custom_family'", fixed = TRUE)
  expect_match(msg, "Custom ordinal likelihoods", fixed = TRUE)
  expect_match(msg, "ID 'custom-ordinal'", fixed = TRUE)
  expect_match(msg, "Remedy:", fixed = TRUE)
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

test_that("shared mvbind IDs use the monolithic cross-predictor system", {
  form <- bf(
    mvbind(y, z) ~ 1 + (1 | shared | gr(g, s2z = TRUE))
  ) + set_rescor(FALSE)
  code <- stancode(form, data = s2z_validation_dat)

  expect_s3_class(code, "brmsmodel")
  expect_match(code, "vector[1] theta_s2z_y;", fixed = TRUE)
  expect_match(code, "vector[1] theta_s2z_z;", fixed = TRUE)
  expect_match(code, "q_recovered_s2z_", fixed = TRUE)
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

test_that("known covariance composes with conventional S2Z blocks", {
  phylo_dat <- data.frame(
    phen = seq(-1, 1, length.out = 12),
    cofactor = rep(c(-0.5, 0.5), 6),
    phylo = factor(rep(letters[1:4], each = 3))
  )
  A <- diag(4)
  dimnames(A) <- list(levels(phylo_dat$phylo), levels(phylo_dat$phylo))
  code <- stancode(
    phen ~ cofactor + (1 | gr(phylo, cov = A, s2z = TRUE)),
    data = phylo_dat, data2 = list(A = A)
  )
  sdata <- standata(
    phen ~ cofactor + (1 | gr(phylo, cov = A, s2z = TRUE)),
    data = phylo_dat, data2 = list(A = A)
  )

  expect_s3_class(code, "brmsmodel")
  expect_match(code, "Lcov_1", fixed = TRUE)
  expect_equal(unname(sdata$Lcov_1), diag(4))
})
