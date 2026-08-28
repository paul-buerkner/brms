context("Stan-expression prior arguments for S2Z effects")

expect_match2 <- brms:::expect_match2

s2z_symbolic_dat <- local({
  n <- 30L
  dev <- rep(seq_len(5L), 6L) * 6
  data.frame(
    y = sin(seq_len(n) / 4),
    y_ord = ordered(rep(letters[1:3], length.out = n)),
    cum = 5000 * (1 - exp(-(dev / 45))),
    dev,
    x = rep(c(-1, 1), length.out = n),
    x2 = rep(c(-1, 0, 1), length.out = n),
    x3 = rep(c(-2, -1, 0, 1, 2), length.out = n),
    AY = factor(rep(seq_len(6L), each = 5L)),
    g = factor(rep(seq_len(6L), each = 5L)),
    h = factor(rep(seq_len(5L), each = 6L))
  )
})

test_that("S2Z prior parser retains scalar Stan expressions", {
  stanvars <-
    stanvar(8.51, "prior_location") +
    stanvar(0.25, "prior_scale") +
    stanvar(5, "prior_df")

  spec <- parse_re_s2z_prior(
    "student_t(prior_df, prior_location, prior_scale)",
    stanvars = stanvars
  )
  expect_identical(spec$dist, "student")
  expect_equal(
    unname(vapply(
      spec[c("df", "location", "scale")], `[[`, numeric(1), "value"
    )),
    c(5, 8.51, 0.25)
  )
  expect_identical(
    stan_s2z_arg_code(spec$df),
    "s2z_require_positive_brms(prior_df)"
  )
  expect_identical(
    stan_s2z_arg_code(spec$location),
    "s2z_require_finite_brms(prior_location)"
  )
  expect_identical(
    stan_s2z_arg_code(spec$scale),
    "s2z_require_positive_brms(prior_scale)"
  )

  logistic <- parse_re_s2z_prior(
    "logistic(prior_location, prior_scale)", stanvars = stanvars
  )
  expect_identical(
    stan_s2z_arg_code(logistic$location),
    "s2z_require_finite_brms(prior_location)"
  )
  expect_identical(
    stan_s2z_arg_code(logistic$scale),
    "s2z_require_positive_brms(prior_scale)"
  )

  constrained <- stanvar(
    0.25, "prior_scale", scode = "real<lower=0> prior_scale;"
  )
  expect_no_error(parse_re_s2z_prior(
    "normal(0, prior_scale)", stanvars = constrained
  ))
})

test_that("S2Z prior expressions enforce scope, shape, and support", {
  expect_error(
    parse_re_s2z_prior("normal(missing_location, 1)"),
    "refers to unknown variable"
  )
  arithmetic <- parse_re_s2z_prior(
    "normal(prior_location + 1, 1)",
    stanvars = stanvar(8.51, "prior_location")
  )
  expect_identical(
    stan_s2z_arg_code(arithmetic$location),
    "s2z_require_finite_brms(prior_location + 1)"
  )
  length_one <- parse_re_s2z_prior(
    "normal(prior_location, 1)",
    stanvars = stanvar(array(8.51), "prior_location")
  )
  expect_identical(
    stan_s2z_arg_code(length_one$location),
    "s2z_require_finite_brms(prior_location[1])"
  )
  expect_error(
    parse_re_s2z_prior(
      "normal(prior_location, 1)",
      stanvars = stanvar(Inf, "prior_location")
    ),
    "must resolve to one finite scalar value"
  )
  expect_error(
    parse_re_s2z_prior(
      "normal(prior_location, 1)",
      stanvars = stanvar(
        1 + 1i, "prior_location", scode = "real prior_location;"
      )
    ),
    "must resolve to one finite scalar value"
  )
  parameter_arg <- parse_re_s2z_prior(
    "normal(prior_location, 1)",
    stanvars = stanvar(
      name = "prior_location", scode = "real prior_location;",
      block = "parameters"
    )
  )
  expect_identical(
    stan_s2z_arg_code(parameter_arg$location),
    "s2z_require_finite_brms(prior_location)"
  )
  logical_arg <- parse_re_s2z_prior(
    "normal(prior_location, 1)",
    stanvars = stanvar(TRUE, "prior_location")
  )
  expect_identical(logical_arg$location$value, 1)
  expect_identical(
    stan_s2z_arg_code(logical_arg$location),
    "s2z_require_finite_brms(prior_location)"
  )
  array_arg <- parse_re_s2z_prior(
    "normal(prior_location, 1)",
    stanvars = stanvar(
      name = "prior_location", scode = "array[1] real prior_location;",
      block = "parameters"
    )
  )
  expect_identical(
    stan_s2z_arg_code(array_arg$location),
    "s2z_require_finite_brms(prior_location[1])"
  )
  expect_error(
    parse_re_s2z_prior(
      "normal(prior_location, 1)",
      stanvars = stanvar(
        name = "prior_location", scode = "real prior_location = 0;",
        block = "tparameters", position = "end"
      )
    ),
    "declared too late"
  )
  expect_error(
    parse_re_s2z_prior(
      "normal(0, prior_scale)",
      stanvars = stanvar(0, "prior_scale")
    ),
    "Scale and degrees-of-freedom arguments must be positive"
  )

  rng_named_data <- parse_re_s2z_prior(
    "normal(prior_rng, 1)", stanvars = stanvar(0, "prior_rng")
  )
  expect_identical(
    stan_s2z_arg_code(rng_named_data$location),
    "s2z_require_finite_brms(prior_rng)"
  )
  lp_named_data <- parse_re_s2z_prior(
    "normal(prior_lp, 1)", stanvars = stanvar(0, "prior_lp")
  )
  expect_identical(
    stan_s2z_arg_code(lp_named_data$location),
    "s2z_require_finite_brms(prior_lp)"
  )
  expect_error(
    parse_re_s2z_prior("normal(prior_rng(), 1)"),
    "contains an RNG or target-modifying function"
  )
  expect_error(
    parse_re_s2z_prior("normal(prior_lp(), 1)"),
    "contains an RNG or target-modifying function"
  )
  expect_error(
    parse_re_s2z_prior("normal(target(), 1)"),
    "contains an RNG or target-modifying function"
  )

  tdata_end <- parse_re_s2z_prior(
    "normal(prior_tdata_end, 1)",
    stanvars = stanvar(
      scode = "real prior_tdata_end = 0;", block = "tdata",
      position = "end"
    )
  )
  expect_identical(
    stan_s2z_arg_code(tdata_end$location),
    "s2z_require_finite_brms(prior_tdata_end)"
  )
  for (late_block in c("model", "likelihood", "genquant")) {
    expect_error(
      parse_re_s2z_prior(
        "normal(prior_late, 1)",
        stanvars = stanvar(
          scode = "real prior_late = 0;", block = late_block
        )
      ),
      "declared too late"
    )
  }
  expect_error(
    parse_re_s2z_prior(
      "normal(prior_hidden, 1)",
      stanvars = stanvar(
        scode = paste0(
          "real hidden_fun(real x) { ",
          "real prior_hidden = x; return prior_hidden; }"
        ),
        block = "functions"
      )
    ),
    "refers to unknown variable"
  )
  expect_error(
    parse_re_s2z_prior(
      "normal(prior_commented, 1)",
      stanvars = stanvar(
        scode = "// real prior_commented;\nreal prior_visible = 0;",
        block = "tdata"
      )
    ),
    "refers to unknown variable"
  )
  expect_error(
    parse_re_s2z_prior(
      "normal(prior_nested, 1)",
      stanvars = stanvar(
        name = "prior_nested", scode = "{ real prior_nested = 0; }",
        block = "tdata"
      )
    ),
    "refers to unknown variable"
  )
  expect_error(
    parse_re_s2z_prior(
      "normal(prior_commented, 1)",
      stanvars = stanvar(
        0, "prior_commented", scode = "// real prior_commented;"
      )
    ),
    "must resolve to one finite scalar value"
  )
  expect_error(
    parse_re_s2z_prior(
      "student_t(prior_df, 0, 1)",
      stanvars = stanvar(-1, "prior_df")
    ),
    "Scale and degrees-of-freedom arguments must be positive"
  )
  expect_error(
    validate_prior(
      prior(normal(prior_location, 1), class = Intercept),
      y ~ 1 + (1 | gr(g, s2z = TRUE)),
      data = s2z_symbolic_dat, stanvars = "not a stanvars object"
    ),
    "Argument 'stanvars' is invalid"
  )
})

test_that("data prior locations work for predictor-local nonlinear S2Z", {
  form <- bf(
    cum ~ exp(ult) *
      (1 - exp(-(dev / exp(theta))^exp(omega))),
    ult ~ 1 + (1 | gr(AY, s2z = TRUE, center = FALSE)),
    omega ~ 1 + (1 | gr(AY, s2z = TRUE, center = FALSE)),
    theta ~ 1 + (1 | gr(AY, s2z = TRUE, center = FALSE)),
    nl = TRUE
  )
  bprior <- c(
    prior(normal(prior_ult, 0.25), nlpar = "ult"),
    prior(normal(0, 0.5), nlpar = "omega"),
    prior(normal(3.80666248977032, 0.3), nlpar = "theta")
  )
  stanvars <- stanvar(array(8.51), "prior_ult")
  scode <- stancode(
    form, data = s2z_symbolic_dat, prior = bprior,
    stanvars = stanvars, parse = TRUE
  )

  expect_match2(scode, "vector[1] prior_ult;")
  expect_match2(
    scode,
    "prior_mean_s2z_1[1] = s2z_require_finite_brms(prior_ult[1]);"
  )
  expect_match2(
    scode, paste0(
      "normal_lpdf(qhat_s2z_1[1] | ",
      "s2z_require_finite_brms(prior_ult[1]), 0.25);"
    )
  )
  expect_no_error(validate_prior(
    bprior, form, data = s2z_symbolic_dat, stanvars = stanvars
  ))

  sdata <- standata(
    form, data = s2z_symbolic_dat, prior = bprior, stanvars = stanvars
  )
  expect_identical(sdata$prior_ult, array(8.51))
  expect_s3_class(brm(
    form, data = s2z_symbolic_dat, prior = bprior,
    stanvars = stanvars, empty = TRUE, chains = 1, cores = 1
  ), "brmsfit")

  changed_stanvars <- stanvar(array(9.25), "prior_ult")
  changed_code <- stancode(
    form, data = s2z_symbolic_dat, prior = bprior,
    stanvars = changed_stanvars
  )
  expect_identical(
    trimws(as.character(changed_code)), trimws(as.character(scode))
  )
  expect_identical(standata(
    form, data = s2z_symbolic_dat, prior = bprior,
    stanvars = changed_stanvars
  )$prior_ult, array(9.25))
})

test_that("scalar data priors reach all analytic S2Z systems", {
  stanvars <-
    stanvar(0.2, "prior_location") +
    stanvar(1.3, "prior_scale") +
    stanvar(5, "prior_df")

  logistic_code <- stancode(
    y ~ 1 + (1 | gr(g, s2z = TRUE, center = TRUE)),
    data = s2z_symbolic_dat,
    prior = prior(
      logistic(prior_location, prior_scale), class = Intercept
    ),
    stanvars = stanvars, parse = TRUE
  )
  expect_match2(
    logistic_code,
    paste0(
      "logistic_lpdf(q_explicit_s2z_1[1] | ",
      "s2z_require_finite_brms(prior_location), ",
      "s2z_require_positive_brms(prior_scale));"
    )
  )

  varying_code <- stancode(
    y ~ 1 + (1 | gr(g, s2z = TRUE, scale = "varying")),
    data = s2z_symbolic_dat,
    prior = prior(normal(prior_location, prior_scale), class = Intercept),
    stanvars = stanvars, parse = TRUE
  )
  expect_match2(
    varying_code, paste0(
      "prior_mean_s2z_1[1] = ",
      "s2z_require_finite_brms(prior_location);"
    )
  )
  expect_match2(
    varying_code, paste0(
      "prior_prec_s2z_1[1] = ",
      "inv_square(s2z_require_positive_brms(prior_scale));"
    )
  )

  matheron_code <- stancode(
    y ~ 1 +
      (1 | gr(g, s2z = TRUE)) +
      (1 | gr(h, s2z = TRUE)),
    data = s2z_symbolic_dat,
    prior = prior(
      student_t(prior_df, prior_location, prior_scale), class = Intercept
    ),
    stanvars = stanvars, parse = TRUE
  )
  expect_match2(matheron_code, "W_matheron_s2z_1")
  expect_match2(
    matheron_code, paste0(
      "inv_chi_square_lpdf(udf_b_s2z_1 | ",
      "s2z_require_positive_brms(prior_df));"
    )
  )
  expect_match2(
    matheron_code,
    paste0(
      "prior_scale_s2z_1[1] = ",
      "s2z_require_positive_brms(prior_scale) * ",
      "sqrt(s2z_require_positive_brms(prior_df) * udf_b_s2z_1);"
    )
  )

  ordinal_code <- stancode(
    y_ord ~ 1 + (1 | gr(g, s2z = TRUE)),
    data = s2z_symbolic_dat, family = cumulative(),
    prior = prior(
      normal(prior_location, prior_scale), class = Intercept
    ),
    stanvars = stanvars, parse = TRUE
  )
  expect_match2(ordinal_code, paste0(
    "prior_mean_s2z_1[1] = ",
    "s2z_require_finite_brms(prior_location);"
  ))
  expect_match2(
    ordinal_code, "inv_square(s2z_require_positive_brms(prior_scale))"
  )
})

test_that("data priors preserve nonlinear cross-ID mean noncentering", {
  form <- bf(
    cum ~ exp(ult) *
      (1 - exp(-(dev / exp(theta))^exp(omega))),
    ult ~ 1 + (1 | z | gr(AY, s2z = TRUE, center = FALSE)),
    omega ~ 1 + (1 | z | gr(AY, s2z = TRUE, center = FALSE)),
    theta ~ 1 + (1 | z | gr(AY, s2z = TRUE, center = FALSE)),
    nl = TRUE
  )
  bprior <- c(
    prior(normal(prior_ult, prior_ult_scale), nlpar = "ult"),
    prior(normal(0, 0.5), nlpar = "omega"),
    prior(normal(3.80666248977032, 0.3), nlpar = "theta")
  )
  stanvars <-
    stanvar(8.51, "prior_ult") +
    stanvar(0.25, "prior_ult_scale")
  scode <- stancode(
    form, data = s2z_symbolic_dat, prior = bprior,
    stanvars = stanvars, parse = TRUE
  )

  expect_match2(scode, "vector[1] z_theta_s2z_ult;")
  expect_match2(scode, paste0(
    "prior_mean_s2z_1[1] = s2z_require_finite_brms(prior_ult);"
  ))
  expect_match2(
    scode, paste0(
      "prior_prec_s2z_1[1] = ",
      "inv_square(s2z_require_positive_brms(prior_ult_scale));"
    )
  )
  expect_match2(
    scode,
    paste0(
      "normal_lpdf(theta_s2z_ult[1] | ",
      "s2z_require_finite_brms(prior_ult), ",
      "s2z_require_positive_brms(prior_ult_scale));"
    )
  )
})

test_that("early transformed-parameter priors use the physical cross chart", {
  form <- bf(
    cum ~ exp(ult) *
      (1 - exp(-(dev / exp(theta))^exp(omega))),
    ult ~ 1 + (1 | z | gr(AY, s2z = TRUE, center = FALSE)),
    omega ~ 1 + (1 | z | gr(AY, s2z = TRUE, center = FALSE)),
    theta ~ 1 + (1 | z | gr(AY, s2z = TRUE, center = FALSE)),
    nl = TRUE
  )
  bprior <- c(
    prior(normal(prior_tp, 0.25), nlpar = "ult"),
    prior(normal(0, 0.5), nlpar = "omega"),
    prior(normal(3.80666248977032, 0.3), nlpar = "theta")
  )
  stanvars <-
    stanvar(8.51, "prior_raw") +
    stanvar(
      scode = "real prior_tp = prior_raw + 0.1;",
      block = "tparameters", position = "start"
    )
  scode <- stancode(
    form, data = s2z_symbolic_dat, prior = bprior,
    stanvars = stanvars, parse = TRUE
  )

  expect_match2(scode, "vector[1] theta_s2z_ult;")
  expect_false(grepl("z_theta_s2z_ult", scode, fixed = TRUE))
  expect_false(grepl("L_mean_s2z_1", scode, fixed = TRUE))
  expect_match2(
    scode,
    "prior_mean_s2z_1[1] = s2z_require_finite_brms(prior_tp);"
  )
})

test_that("S2Z accepts ordinary scalar Stan prior expressions", {
  stanvars <-
    stanvar(2, "raw_mu") +
    stanvar(4, "raw_var") +
    stanvar(-2, "prior_alt") +
    stanvar(
      scode = "real prior_td = raw_mu + 0.5;",
      block = "tdata"
    ) +
    stanvar(
      scode = "real prior_tp = prior_td + 0.1;",
      block = "tparameters", position = "start"
    ) +
    stanvar(
      scode = "real<lower=0> prior_hyper;",
      block = "parameters"
    ) +
    stanvar(
      scode = paste0(
        "target += lognormal_lpdf(prior_hyper | 0, 1);"
      ),
      block = "model"
    ) +
    stanvar(
      scode = paste0(
        "real shift_prior_s2z(real x) { return x + 0.25; }"
      ),
      block = "functions"
    )
  bprior <- prior(
    normal(
      shift_prior_s2z(prior_tp) + log(raw_mu),
      sqrt(raw_var) * prior_hyper
    ),
    class = Intercept
  )

  normalized <- stancode(
    y ~ 1 + (1 | gr(g, s2z = TRUE)),
    data = s2z_symbolic_dat, prior = bprior,
    stanvars = stanvars, parse = TRUE
  )
  expect_match2(
    normalized,
    paste0(
      "s2z_require_finite_brms(shift_prior_s2z(prior_tp) + ",
      "log(raw_mu))"
    )
  )
  expect_match2(
    normalized,
    paste0(
      "s2z_require_positive_brms(sqrt(raw_var) * prior_hyper)"
    )
  )

  unnormalized <- stancode(
    y ~ 1 + (1 | gr(g, s2z = TRUE)),
    data = s2z_symbolic_dat, prior = bprior,
    stanvars = stanvars, normalize = FALSE, parse = TRUE
  )
  expect_match2(unnormalized, "normal_lupdf(qhat_s2z_1[1]")

  ternary <- stancode(
    y ~ 1 + (1 | gr(g, s2z = TRUE)),
    data = s2z_symbolic_dat,
    prior = set_prior(
      "normal(raw_mu > 0 ? raw_mu : prior_alt, 1)",
      class = "Intercept"
    ),
    stanvars = stanvars, parse = TRUE
  )
  expect_match2(
    ternary,
    "s2z_require_finite_brms(raw_mu > 0 ? raw_mu : prior_alt)"
  )

  nested_ternary <- stancode(
    y ~ 1 + (1 | gr(g, s2z = TRUE)),
    data = s2z_symbolic_dat,
    prior = set_prior(
      "normal(log(raw_mu > 0 ? raw_mu : prior_alt), 1)",
      class = "Intercept"
    ),
    stanvars = stanvars, parse = TRUE
  )
  expect_match2(
    nested_ternary,
    paste0(
      "s2z_require_finite_brms(log((raw_mu > 0 ? raw_mu : ",
      "prior_alt)))"
    )
  )

  spaced_signed <- stancode(
    y ~ 1 + (1 | gr(g, s2z = TRUE)),
    data = s2z_symbolic_dat,
    prior = set_prior(
      "normal(raw_mu < -1 ? prior_alt : raw_mu, 1)",
      class = "Intercept"
    ),
    stanvars = stanvars, parse = TRUE
  )
  expect_match2(
    spaced_signed,
    "s2z_require_finite_brms(raw_mu < -1 ? prior_alt : raw_mu)"
  )
})

test_that("parameter hyperarguments reach the Matheron system", {
  stanvars <-
    stanvar(
      scode = paste0(
        "real prior_mu; real<lower=0> prior_scale; ",
        "real<lower=0> prior_df;"
      ),
      block = "parameters"
    ) +
    stanvar(
      scode = paste0(
        "target += normal_lpdf(prior_mu | 0, 1); ",
        "target += lognormal_lpdf(prior_scale | 0, 1); ",
        "target += lognormal_lpdf(prior_df | 1, 0.5);"
      ),
      block = "model"
    )
  scode <- stancode(
    y ~ 1 +
      (1 | gr(g, s2z = TRUE)) +
      (1 | gr(h, s2z = TRUE)),
    data = s2z_symbolic_dat,
    prior = prior(
      student_t(prior_df, prior_mu, prior_scale), class = Intercept
    ),
    stanvars = stanvars, parse = TRUE
  )

  expect_match2(scode, "W_matheron_s2z_1")
  expect_match2(
    scode, paste0(
      "inv_chi_square_lpdf(udf_b_s2z_1 | ",
      "s2z_require_positive_brms(prior_df));"
    )
  )
  expect_match2(
    scode, "s2z_require_positive_brms(prior_scale) * sqrt("
  )
})

test_that("global vector priors retain original coefficient indices", {
  bprior <- prior(normal(prior_mu, prior_scale), class = b)
  stanvars <-
    stanvar(c(0.1, 0.2), "prior_mu") +
    stanvar(c(1.1, 1.2), "prior_scale")

  all_active <- stancode(
    y ~ 1 + x + x2 + (0 + x + x2 | gr(g, s2z = TRUE)),
    data = s2z_symbolic_dat, prior = bprior,
    stanvars = stanvars, parse = TRUE
  )
  expect_match2(
    all_active,
    "prior_mean_s2z_1[2] = s2z_require_finite_brms(prior_mu[1]);"
  )
  expect_match2(
    all_active,
    "prior_mean_s2z_1[3] = s2z_require_finite_brms(prior_mu[2]);"
  )
  expect_match2(
    all_active,
    "inv_square(s2z_require_positive_brms(prior_scale[2]))"
  )

  split <- stancode(
    y ~ 1 + x + x2 + (0 + x | gr(g, s2z = TRUE)),
    data = s2z_symbolic_dat, prior = bprior,
    stanvars = stanvars, parse = TRUE
  )
  expect_match2(
    split,
    paste0(
      "normal_lpdf(fixed_s2z[1] | ",
      "s2z_require_finite_brms(prior_mu[2]), ",
      "s2z_require_positive_brms(prior_scale[2]));"
    )
  )
  expect_match2(
    split,
    "prior_mean_s2z_1[2] = s2z_require_finite_brms(prior_mu[1]);"
  )

  active_override <- stancode(
    y ~ 1 + x + x2 + (0 + x | gr(g, s2z = TRUE)),
    data = s2z_symbolic_dat,
    prior = bprior + prior(normal(0, 1), class = b, coef = x),
    stanvars = stanvars, parse = TRUE
  )
  expect_match2(
    active_override,
    paste0(
      "normal_lpdf(fixed_s2z[1] | ",
      "s2z_require_finite_brms(prior_mu[2]), ",
      "s2z_require_positive_brms(prior_scale[2]));"
    )
  )

  parameter_code <- stancode(
    y ~ 1 + x + x2 + (0 + x + x2 | gr(g, s2z = TRUE)),
    data = s2z_symbolic_dat, prior = bprior,
    stanvars = stanvar(
      scode = paste0(
        "vector[2] prior_mu; ",
        "vector<lower=0>[2] prior_scale;"
      ),
      block = "parameters"
    ) + stanvar(
      scode = paste0(
        "target += normal_lpdf(prior_mu | 0, 1); ",
        "target += lognormal_lpdf(prior_scale | 0, 1);"
      ),
      block = "model"
    ),
    parse = TRUE
  )
  expect_match2(
    parameter_code,
    paste0(
      "prior_mean_s2z_1[3] = s2z_require_finite_brms(",
      "s2z_prior_coordinate_brms(prior_mu, 2, 2));"
    )
  )

  structured <- stancode(
    y ~ 1 + x + x2 + (0 + x + x2 | gr(g, s2z = TRUE)),
    data = s2z_symbolic_dat,
    prior = prior(
      normal(prior_mu + prior_shift, sqrt(prior_var)), class = b
    ),
    stanvars = stanvar(c(0.1, 0.2), "prior_mu") +
      stanvar(c(0.01, 0.02), "prior_shift") +
      stanvar(c(1.1, 1.2), "prior_var"),
    parse = TRUE
  )
  expect_match2(structured, paste0(
    "prior_mean_s2z_1[2] = s2z_require_finite_brms(",
    "s2z_prior_coordinate_brms(prior_mu + prior_shift, 1, 2));"
  ))
  expect_match2(structured, paste0(
    "prior_prec_s2z_1[3] = inv_square(s2z_require_positive_brms(",
    "s2z_prior_coordinate_brms(sqrt(prior_var), 2, 2)));"
  ))

  noncontiguous <- stancode(
    y ~ 1 + x + x2 + x3 +
      (0 + x + x3 | gr(g, s2z = TRUE)),
    data = s2z_symbolic_dat,
    prior = prior(normal(prior_mu, prior_scale), class = b),
    stanvars = stanvar(c(0.1, 0.2, 0.3), "prior_mu") +
      stanvar(c(1.1, 1.2, 1.3), "prior_scale"),
    parse = TRUE
  )
  expect_match2(
    noncontiguous,
    "prior_mean_s2z_1[2] = s2z_require_finite_brms(prior_mu[1]);"
  )
  expect_match2(
    noncontiguous,
    "prior_mean_s2z_1[4] = s2z_require_finite_brms(prior_mu[3]);"
  )
  expect_match2(noncontiguous, paste0(
    "normal_lpdf(fixed_s2z[1] | ",
    "s2z_require_finite_brms(prior_mu[2]), ",
    "s2z_require_positive_brms(prior_scale[2]));"
  ))

  expect_error(
    stancode(
      y ~ 1 + x + x2 + (0 + x + x2 | gr(g, s2z = TRUE)),
      data = s2z_symbolic_dat, prior = bprior,
      stanvars = stanvar(
        scode = paste0(
          "vector[3] prior_mu; ",
          "vector<lower=0>[3] prior_scale;"
        ),
        block = "parameters"
      )
    ),
    "has declared length '3'.*coefficient vector has length 2"
  )

  expect_error(
    stancode(
      y ~ 1 + x + x2 + (0 + x + x2 | gr(g, s2z = TRUE)),
      data = s2z_symbolic_dat, prior = bprior,
      stanvars = stanvar(c(0.1, 0.2, 0.3), "prior_mu") +
        stanvar(c(1.1, 1.2), "prior_scale")
    ),
    "has length 3 but the population-level coefficient vector has length 2"
  )

  ordinal <- stancode(
    y_ord ~ 1 + (1 | gr(g, s2z = TRUE)),
    data = s2z_symbolic_dat, family = cumulative(),
    prior = prior(normal(prior_mu, prior_scale), class = Intercept),
    stanvars = stanvars, parse = TRUE
  )
  expect_match2(
    ordinal,
    "prior_mean_s2z_1[1] = s2z_require_finite_brms(prior_mu[1]);"
  )
  expect_match2(
    ordinal,
    "prior_mean_s2z_1[2] = s2z_require_finite_brms(prior_mu[2]);"
  )

  local_vector <- stancode(
    y ~ 1 + x + x2 + (0 + x | gr(g, s2z = TRUE)),
    data = s2z_symbolic_dat,
    prior = c(
      prior(normal(prior_one, 1), class = b, coef = x),
      prior(logistic(0, 1), class = b, coef = x2)
    ),
    stanvars = stanvar(array(0.1), "prior_one"),
    parse = TRUE
  )
  expect_match2(
    local_vector,
    "logistic_lpdf(fixed_s2z[1] | 0, 1);"
  )
  expect_match2(
    local_vector,
    "prior_mean_s2z_1[2] = s2z_require_finite_brms(prior_one[1]);"
  )
})

test_that("S2Z prior data can be updated without recompilation", {
  form <- y ~ 1 + (1 | gr(g, s2z = TRUE))
  bprior <- prior(
    normal(prior_location, prior_scale), class = Intercept
  )
  stanvars <-
    stanvar(0, "prior_location") +
    stanvar(1, "prior_scale")
  fit <- brm(
    form, data = s2z_symbolic_dat, prior = bprior,
    stanvars = stanvars, empty = TRUE, chains = 1, cores = 1
  )

  updated <- update(
    fit, recompile = FALSE, algorithm = "meanfield", testmode = TRUE,
    stanvars = stanvar(2, "prior_location") +
      stanvar(1.5, "prior_scale")
  )
  expect_identical(updated$stanvars$prior_location$sdata, 2)
  expect_identical(updated$stanvars$prior_scale$sdata, 1.5)

  new_sdata <- standata(
    fit, newdata2 = list(prior_location = 3, prior_scale = 1.75)
  )
  expect_identical(new_sdata$prior_location, 3)
  expect_identical(new_sdata$prior_scale, 1.75)

  expect_error(
    update(
      fit, recompile = FALSE, algorithm = "meanfield", testmode = TRUE,
      stanvars = stanvar(2, "prior_location") +
        stanvar(0, "prior_scale")
    ),
    "Scale and degrees-of-freedom arguments must be positive"
  )
  expect_error(
    standata(fit, newdata2 = list(prior_scale = 0)),
    "Scale and degrees-of-freedom arguments must be positive"
  )
})

test_that("symbolic S2Z priors render unnormalized densities", {
  stanvars <-
    stanvar(0.2, "prior_location") +
    stanvar(1.3, "prior_scale") +
    stanvar(5, "prior_df")
  scode <- stancode(
    y ~ 1 + (1 | gr(g, s2z = TRUE)),
    data = s2z_symbolic_dat,
    prior = prior(
      student_t(prior_df, prior_location, prior_scale), class = Intercept
    ),
    stanvars = stanvars, normalize = FALSE, parse = TRUE
  )

  expect_match2(
    scode, paste0(
      "inv_chi_square_lupdf(udf_b_s2z_1 | ",
      "s2z_require_positive_brms(prior_df));"
    )
  )
  expect_match2(
    scode,
    paste0(
      "normal_lupdf(qhat_s2z_1[1] | ",
      "s2z_require_finite_brms(prior_location), ",
      "s2z_require_positive_brms(prior_scale) * ",
      "sqrt(s2z_require_positive_brms(prior_df) * udf_b_s2z_1));"
    )
  )
})

test_that("symbolic priors reach correlated and independent S2Z paths", {
  stanvars <-
    stanvar(0.2, "prior_location") +
    stanvar(1.3, "prior_scale")
  bprior <- prior(
    normal(prior_location, prior_scale), class = Intercept
  )

  correlated <- stancode(
    y ~ 1 + x + (1 + x | gr(g, s2z = TRUE)),
    data = s2z_symbolic_dat, prior = bprior,
    stanvars = stanvars, parse = TRUE
  )
  expect_match2(correlated, "matrix[M_1, M_1] P_s2z_1;")
  expect_match2(correlated, paste0(
    "prior_mean_s2z_1[1] = ",
    "s2z_require_finite_brms(prior_location);"
  ))
  expect_match2(correlated, paste0(
    "prior_prec_s2z_1[1] = ",
    "inv_square(s2z_require_positive_brms(prior_scale));"
  ))

  independent <- stancode(
    y ~ 1 + x + (1 + x || gr(g, s2z = TRUE)),
    data = s2z_symbolic_dat, prior = bprior,
    stanvars = stanvars, parse = TRUE
  )
  expect_match2(independent, "vector<lower=0>[M_1] D_diag_s2z_1;")
  expect_match2(independent, paste0(
    "prior_mean_s2z_1[1] = ",
    "s2z_require_finite_brms(prior_location);"
  ))
  expect_match2(independent, paste0(
    "prior_prec_s2z_1[1] = ",
    "inv_square(s2z_require_positive_brms(prior_scale));"
  ))
})

test_that("structured scalar priors reach automatic S2Z centering", {
  scode <- stancode(
    y ~ 1 + (1 | gr(g, s2z = TRUE, center = "auto")),
    data = s2z_symbolic_dat,
    prior = prior(
      normal(log(prior_raw_location), sqrt(prior_variance)),
      class = Intercept
    ),
    stanvars = stanvar(2, "prior_raw_location") +
      stanvar(4, "prior_variance"),
    parse = TRUE
  )
  expect_match2(
    scode, "s2z_require_finite_brms(log(prior_raw_location))"
  )
  expect_match2(
    scode, "s2z_require_positive_brms(sqrt(prior_variance))"
  )
  expect_match2(scode, "rho_s2z_1")
})

test_that("multi-element stanvar declarations get contextual S2Z errors", {
  stanvars <- stanvar(
    0.2, "prior_location",
    scode = c("real prior_location;", "real another_location;")
  )
  expect_error(
    parse_re_s2z_prior(
      "normal(prior_location, 1)", stanvars = stanvars
    ),
    "must resolve to one finite scalar value"
  )
})
