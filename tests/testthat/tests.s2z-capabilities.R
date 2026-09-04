context("Capability matrix for physical sum-to-zero group-level effects")

expect_match2 <- brms:::expect_match2

s2z_cap_count_fixed <- function(x, pattern) {
  unname(lengths(regmatches(
    x, gregexpr(pattern, x, fixed = TRUE)
  )))
}

s2z_cap_error <- function(expr) {
  error <- tryCatch(expr, error = identity)
  expect_s3_class(error, "error")
  conditionMessage(error)
}

s2z_cap_code <- function(case, normalize = TRUE) {
  args <- list(
    object = case$formula, data = s2z_cap_dat,
    normalize = normalize, parse = TRUE
  )
  if (!is.null(case$family)) {
    args$family <- case$family
  }
  if (!is.null(case$prior)) {
    args$prior <- case$prior
  }
  do.call(stancode, args)
}

s2z_cap_dat <- local({
  n <- 36L
  i <- seq_len(n)
  out <- data.frame(
    y = sin(i / 4) + cos(i / 9),
    y2 = cos(i / 5) - sin(i / 11),
    y_binary = as.integer(i %% 3L == 0L),
    y_count = i %% 5L,
    y_cat = factor(rep(1:3, length.out = n)),
    y_beta = rep(c(0, 1, 0.2, 0.65), length.out = n),
    rt = 0.4 + i / n,
    decision = i %% 2L,
    x = seq(-1.4, 1.6, length.out = n),
    z = cos(i / 6),
    w = sin(i / 7),
    g = factor(rep(seq_len(6L), each = 6L)),
    h = factor(rep(seq_len(4L), each = 9L))
  )

  ymulti <- cbind(
    ymulti1 = 1L + i %% 3L,
    ymulti2 = 1L + (i + 1L) %% 4L,
    ymulti3 = 1L + (i + 2L) %% 5L
  )
  out$ymulti <- ymulti
  out$multi_size <- rowSums(ymulti)

  ysimplex <- cbind(
    ysimplex1 = 1 + i %% 4L,
    ysimplex2 = 2 + (i + 1L) %% 5L,
    ysimplex3 = 3 + (i + 2L) %% 6L
  )
  out$ysimplex <- ysimplex / rowSums(ysimplex)
  out
})

test_that("S2Z Stan parses across core likelihood capabilities", {
  skip_if_not_installed("rstan")

  cases <- list(
    Gaussian = list(
      formula = y ~ x + (1 + x | gr(g, id = "gauss", s2z = TRUE)),
      family = gaussian(), prior = NULL, marker = "normal_id_glm"
    ),
    Bernoulli = list(
      formula = y_binary ~ x +
        (1 + x | gr(g, id = "bern", s2z = TRUE)),
      family = bernoulli(), prior = NULL, marker = "bernoulli_logit"
    ),
    Poisson = list(
      formula = y_count ~ x +
        (1 + x | gr(g, id = "pois", s2z = TRUE)),
      family = poisson(), prior = NULL, marker = "poisson_log"
    ),
    categorical = list(
      formula = bf(
        y_cat ~ 1,
        mu2 ~ x + (1 + x | gr(g, id = "cat-two", s2z = TRUE)),
        mu3 ~ x + (1 + x | gr(h, id = "cat-three", s2z = TRUE))
      ),
      family = categorical(), prior = NULL, marker = "categorical_logit"
    ),
    multinomial = list(
      formula = bf(
        ymulti | trials(multi_size) ~ 1,
        muymulti2 ~ x +
          (1 + x | gr(g, id = "multi-two", s2z = TRUE)),
        muymulti3 ~ x +
          (1 + x | gr(h, id = "multi-three", s2z = TRUE))
      ),
      family = multinomial(), prior = NULL, marker = "multinomial_logit2"
    ),
    simplex = list(
      formula = bf(
        ysimplex ~ 1,
        muysimplex2 ~ x +
          (1 + x | gr(g, id = "simplex-two", s2z = TRUE)),
        muysimplex3 ~ x +
          (1 + x | gr(h, id = "simplex-three", s2z = TRUE))
      ),
      family = dirichlet(), prior = NULL, marker = "dirichlet_logit"
    ),
    distributional = list(
      formula = bf(
        y ~ x,
        sigma ~ x + (1 + x | gr(g, id = "sigma", s2z = TRUE))
      ),
      family = gaussian(), prior = NULL, marker = "theta_s2z_sigma"
    ),
    nonlinear = list(
      formula = bf(
        y ~ a,
        a ~ 1 + x + (1 + x | gr(g, id = "nla", s2z = TRUE)),
        nl = TRUE
      ),
      family = gaussian(),
      prior = prior(normal(0, 1), nlpar = a), marker = "theta_s2z_a"
    ),
    multivariate = list(
      formula = bf(
        y ~ x + (1 + x | gr(g, id = "response-one", s2z = TRUE))
      ) +
        bf(
          y2 ~ x +
            (1 + x | gr(h, id = "response-two", s2z = TRUE))
        ) +
        set_rescor(FALSE),
      family = NULL, prior = NULL, marker = "theta_s2z_y2"
    )
  )

  for (name in names(cases)) {
    code <- s2z_cap_code(cases[[name]])
    expect_s3_class(code, "brmsmodel")
    expect_match2(code, "sum_to_zero_constrain_brms", info = name)
    expect_match2(code, "q_recovered_s2z_", info = name)
    expect_match2(code, cases[[name]]$marker, info = name)
  }
})

test_that("categorical shorthand allocates separate local S2Z IDs", {
  skip_if_not_installed("rstan")

  code <- stancode(
    y_cat ~ x + (1 + x | gr(g, s2z = TRUE)),
    data = s2z_cap_dat, family = categorical(), parse = TRUE
  )
  expect_equal(
    s2z_cap_count_fixed(code, "// physical sum-to-zero group-level effects"),
    2L
  )
  expect_match2(code, "q_recovered_s2z_1")
  expect_match2(code, "q_recovered_s2z_2")
  expect_match2(code, "b_mu2_Intercept")
  expect_match2(code, "b_mu3_Intercept")

  msg <- s2z_cap_error(stancode(
    y_cat ~ x + (1 + x | gr(g, id = "shared", s2z = TRUE)),
    data = s2z_cap_dat, family = categorical()
  ))
  expect_match(msg, "S2Z capability 'cross_predictor_id'", fixed = TRUE)
  expect_match(msg, "response 'y_cat'", fixed = TRUE)
  expect_match(msg, "family 'categorical'", fixed = TRUE)
  expect_match(msg, "ID 'shared'", fixed = TRUE)
  expect_match(msg, "Remedy:", fixed = TRUE)
})

test_that("default probability predictors use exact logistic S2Z means", {
  skip_if_not_installed("rstan")

  cases <- list(
    zi = list(
      formula = bf(
        y_count ~ 1,
        zi ~ 1 + (1 | gr(g, id = "zi", s2z = TRUE))
      ),
      family = zero_inflated_poisson(), prior = NULL,
      marker = "theta_s2z_zi", nlogistic = 1L
    ),
    hu = list(
      formula = bf(
        y_count ~ 1,
        hu ~ 1 + (1 | gr(g, id = "hu", s2z = TRUE))
      ),
      family = hurdle_poisson(), prior = NULL,
      marker = "theta_s2z_hu", nlogistic = 1L
    ),
    zoi_coi = list(
      formula = bf(
        y_beta ~ 1,
        zoi ~ 1 + (1 | gr(g, id = "zoi", s2z = TRUE)),
        coi ~ 1 + (1 | gr(h, id = "coi", s2z = TRUE))
      ),
      family = zero_one_inflated_beta(), prior = NULL,
      marker = "theta_s2z_coi", nlogistic = 2L
    ),
    mixture_theta = list(
      formula = bf(
        y ~ 1,
        theta1 ~ 1 + (1 | gr(g, id = "mixture", s2z = TRUE))
      ),
      family = mixture(gaussian, gaussian, order = "none"), prior = NULL,
      marker = "theta_s2z_theta1", nlogistic = 1L
    ),
    wiener_bias = list(
      formula = bf(
        rt | dec(decision) ~ 1,
        bias ~ 1 + (1 | gr(g, id = "bias", s2z = TRUE))
      ),
      family = wiener(), prior = NULL,
      marker = "theta_s2z_bias", nlogistic = 1L
    ),
    asym_laplace_quantile = list(
      formula = bf(
        y ~ 1,
        quantile ~ 1 + (1 | gr(g, id = "quantile", s2z = TRUE))
      ),
      family = asym_laplace(), prior = NULL,
      marker = "theta_s2z_quantile", nlogistic = 1L
    )
  )

  for (name in names(cases)) {
    code <- s2z_cap_code(cases[[name]])
    expect_s3_class(code, "brmsmodel")
    expect_match2(code, cases[[name]]$marker, info = name)
    expect_match2(code, "z_mean_s2z_", info = name)
    expect_match2(code, "q_explicit_s2z_", info = name)
    expect_equal(
      s2z_cap_count_fixed(
        code, "lprior += logistic_lpdf(q_explicit_s2z_"
      ),
      cases[[name]]$nlogistic,
      info = name
    )
  }
})

test_that("exact logistic means cover every physical S2Z kernel shape", {
  skip_if_not_installed("rstan")

  population_prior <- prior(logistic(-1, 2), class = Intercept) +
    prior(normal(0, 1), class = b, coef = x)
  cases <- list(
    scalar = list(
      formula = y ~ 1 + (1 | gr(g, id = "scalar", s2z = TRUE)),
      family = gaussian(),
      prior = prior(logistic(-1, 2), class = Intercept), nblocks = 1L
    ),
    independent = list(
      formula = y ~ x +
        (1 + x || gr(g, id = "independent", s2z = TRUE)),
      family = gaussian(), prior = population_prior, nblocks = 1L
    ),
    correlated = list(
      formula = y ~ x +
        (1 + x | gr(g, id = "correlated", s2z = TRUE)),
      family = gaussian(), prior = population_prior, nblocks = 1L
    ),
    Student = list(
      formula = y ~ x +
        (1 + x | gr(
          g, id = "student", dist = "student", s2z = TRUE
        )),
      family = gaussian(), prior = population_prior, nblocks = 1L
    ),
    multiblock = list(
      formula = y ~ x +
        (1 + x | gr(g, id = "first", s2z = TRUE)) +
        (1 | gr(
          h, id = "second", dist = "student", s2z = TRUE
        )),
      family = gaussian(), prior = population_prior, nblocks = 2L
    )
  )

  code <- lapply(cases, s2z_cap_code)
  for (name in names(code)) {
    expect_s3_class(code[[name]], "brmsmodel")
    expect_equal(
      s2z_cap_count_fixed(code[[name]], "// exact explicit-mean S2Z block"),
      cases[[name]]$nblocks,
      info = name
    )
    expect_equal(
      s2z_cap_count_fixed(code[[name]], "-0.5 * group_quad_s2z_"),
      cases[[name]]$nblocks,
      info = name
    )
    expect_equal(
      s2z_cap_count_fixed(code[[name]], "- 0.5 * N_"),
      cases[[name]]$nblocks,
      info = paste(name, "normalizing constants")
    )
    expect_equal(
      s2z_cap_count_fixed(
        code[[name]], "sum(log(diagonal(L_Sigma_s2z_"
      ),
      cases[[name]]$nblocks,
      info = paste(name, "restricted Jacobians")
    )
    expect_match2(
      code[[name]], "lprior += logistic_lpdf(q_explicit_s2z_", info = name
    )
    expect_match2(
      code[[name]], "z_mean_s2z_", info = paste(name, "omitted mean")
    )
    expect_match2(code[[name]], "q_recovered_s2z_", info = name)
  }

  expect_match2(
    code$independent, "L_Sigma_s2z_1 = diag_matrix(sd_1);"
  )
  expect_match2(
    code$correlated,
    "L_Sigma_s2z_1 = diag_pre_multiply(sd_1, L_1);"
  )
  expect_match2(code$Student, "group_scale_s2z_1 = dfm_1;")
  expect_match2(
    code$Student, "- M_1 * sum(log(group_scale_s2z_1))"
  )
  expect_match2(code$multiblock, "z_mean_s2z_1")
  expect_match2(code$multiblock, "z_mean_s2z_2")
  expect_equal(
    s2z_cap_count_fixed(
      code$multiblock, "q_explicit_s2z_1 -= H_s2z_"
    ),
    2L
  )

  unnormalized <- s2z_cap_code(cases$scalar, normalize = FALSE)
  expect_match2(
    unnormalized, "lprior += logistic_lupdf(q_explicit_s2z_1[1] | -1, 2);"
  )
  expect_false(grepl(
    "lprior += logistic_lpdf(q_explicit_s2z_1", unnormalized, fixed = TRUE
  ))
  expect_false(grepl(
    "- 0.5 * N_1 * M_1 * log(2 * pi())", unnormalized, fixed = TRUE
  ))
  expect_match2(
    unnormalized,
    "- (N_1 - 1) * sum(log(diagonal(L_Sigma_s2z_1)))"
  )
})

test_that("fixed-only coordinates keep arbitrary ordinary brms priors", {
  skip_if_not_installed("rstan")

  case <- list(
    formula = y ~ x + z + w +
      (1 + x | gr(g, id = "active", s2z = TRUE)),
    family = gaussian(),
    prior = prior(logistic(0, 1), class = Intercept) +
      prior(normal(0.2, 0.8), class = b, coef = x) +
      prior(
        double_exponential(0, 1.3), class = b, coef = z,
        tag = "fixedonly"
      ) +
      prior(constant(0.4), class = b, coef = w)
  )
  code <- s2z_cap_code(case)

  expect_match2(
    code,
    "vector[2] theta_s2z_active;  // S2Z-active finite-population coefficients"
  )
  expect_match2(
    code,
    "vector[2] fixed_s2z;  // S2Z-inactive regression coefficients"
  )
  expect_match2(code, "real par_fixed_s2z_1;")
  expect_match2(code, "fixed_s2z[2] = 0.4;")
  expect_match2(code, "theta_s2z[3] = fixed_s2z[1];")
  expect_match2(code, "theta_s2z[4] = fixed_s2z[2];")
  expect_match2(
    code,
    paste0(
      "lprior_fixedonly += double_exponential_lpdf(",
      "fixed_s2z[1] | 0, 1.3);"
    )
  )
  expect_match2(
    code, "lprior += logistic_lpdf(q_explicit_s2z_1[1] | 0, 1);"
  )
  expect_match2(
    code,
    paste0(
      "lprior += normal_lpdf(q_explicit_s2z_1[2] | ",
      "0.20000000000000001, 0.80000000000000004);"
    )
  )
  expect_false(grepl("q_explicit_s2z_1[3]", code, fixed = TRUE))
  expect_false(grepl("q_explicit_s2z_1[4]", code, fixed = TRUE))
  expect_false(grepl("b_s2z_inactive", code, fixed = TRUE))
})

test_that("additional S2Z gates reject in the intended phase with context", {
  structural <- list(
    sparse = bf(
      y ~ x + (1 + x + z | gr(g, id = "sparse", s2z = TRUE)),
      sparse = TRUE
    ),
    qr = bf(
      y ~ x + (1 + x + z | gr(g, id = "qr", s2z = TRUE)),
      decomp = "QR"
    )
  )
  for (capability in names(structural)) {
    msg <- s2z_cap_error(stancode(
      structural[[capability]], data = s2z_cap_dat
    ))
    expect_match(
      msg, paste0("S2Z capability '", capability, "'"), fixed = TRUE
    )
    expect_match(msg, "response 'y'", fixed = TRUE)
    expect_match(msg, paste0("ID '", capability, "'"), fixed = TRUE)
    expect_match(msg, "Remedy:", fixed = TRUE)
    expect_false(grepl("matching population-level", msg, fixed = TRUE))
  }

  active_priors <- list(
    active_prior_distribution = prior(
      double_exponential(0, 1), class = b, coef = x
    ),
    active_prior_tag = prior(
      normal(0, 1), class = b, coef = x, tag = "active"
    )
  )
  for (capability in names(active_priors)) {
    msg <- s2z_cap_error(stancode(
      y ~ x + (1 + x | gr(g, id = "prior", s2z = TRUE)),
      data = s2z_cap_dat, prior = active_priors[[capability]]
    ))
    expect_match(
      msg, paste0("S2Z capability '", capability, "'"), fixed = TRUE
    )
    expect_match(msg, "coefficient(s) 'x'", fixed = TRUE)
    expect_match(msg, "group 'g'", fixed = TRUE)
    expect_match(msg, "ID 'prior'", fixed = TRUE)
    expect_match(msg, "Remedy:", fixed = TRUE)
  }

  one_level <- transform(s2z_cap_dat, g = factor("only"))
  msg <- s2z_cap_error(stancode(
    y ~ x + (1 + x | gr(g, id = "one-level", s2z = TRUE)),
    data = one_level
  ))
  expect_match(msg, "S2Z capability 'minimum_levels'", fixed = TRUE)
  expect_match(msg, "ID 'one-level'", fixed = TRUE)
  expect_match(msg, "Remedy:", fixed = TRUE)

  retained_one_level <- transform(
    s2z_cap_dat, g = factor("used", levels = c("used", "unused"))
  )
  msg <- s2z_cap_error(standata(
    y ~ x + (1 + x | gr(g, id = "retained-one", s2z = TRUE)),
    data = retained_one_level, drop_unused_levels = FALSE
  ))
  expect_match(msg, "S2Z capability 'minimum_levels'", fixed = TRUE)
  expect_match(msg, "at least two observed grouping levels", fixed = TRUE)

  retained_unused <- transform(
    s2z_cap_dat,
    g = factor(
      rep(c("first", "second"), length.out = nrow(s2z_cap_dat)),
      levels = c("first", "second", "unused")
    )
  )
  msg <- s2z_cap_error(standata(
    y ~ x + (1 + x | gr(g, id = "retained-unused", s2z = TRUE)),
    data = retained_unused, drop_unused_levels = FALSE
  ))
  expect_match(msg, "S2Z capability 'unused_levels'", fixed = TRUE)
  expect_match(msg, "unused", fixed = TRUE)
  expect_match(msg, "drop_unused_levels = TRUE", fixed = TRUE)
})

test_that("foundation code contains no later S2Z APIs or state", {
  public_arguments <- names(formals(gr))
  expect_false("center" %in% public_arguments)
  expect_false("scale" %in% public_arguments)
  expect_error(
    gr(g, s2z = TRUE, center = "auto"),
    "expects only a single grouping term"
  )
  expect_error(
    gr(g, s2z = TRUE, scale = "varying"),
    "expects only a single grouping term"
  )

  code <- stancode(
    y ~ x + (1 + x | gr(g, s2z = TRUE)),
    data = s2z_cap_dat, parse = FALSE
  )
  forbidden <- c(
    "s2z_center", "s2z_center_auto", "rho_s2z", "mean_rho_s2z",
    "sdlog_s2z", "sd_level", "z_sd_s2z"
  )
  for (term in forbidden) {
    expect_false(grepl(term, code, fixed = TRUE), info = term)
  }

  # Audit the production and user-facing files that define this feature, not
  # only one representative Stan program. Construct the symbol names in
  # pieces so this test does not satisfy its own forbidden patterns.
  test_root <- normalizePath(
    file.path(testthat::test_path(), "..", ".."), mustWork = TRUE
  )
  root_candidates <- c(
    test_root, file.path(test_root, "00_pkg_src", "brms")
  )
  is_source_root <- file.exists(file.path(root_candidates, "R/re-s2z.R"))
  expect_true(any(is_source_root))
  package_root <- root_candidates[which(is_source_root)[1L]]
  audit_files <- file.path(package_root, c(
    "R/brmsframe.R", "R/exclude_pars.R", "R/formula-re.R",
    "R/priors.R", "R/re-s2z.R", "R/stan-likelihood.R",
    "R/stan-predictor.R", "R/stancode.R", "R/standata.R",
    "inst/chunks/fun_sum_to_zero.stan", "NEWS.md", "man/gr.Rd",
    "man/mm.Rd", "tests/testthat/tests.brmsterms.R",
    "tests/testthat/tests.s2z-code.R", "tests/testthat/tests.s2z.R",
    "tests/testthat/tests.s2z-optimizations.R",
    "tests/testthat/tests.s2z-validation.R"
  ))
  audit_text <- unlist(lapply(audit_files, readLines, warn = FALSE))
  forbidden_patterns <- c(
    paste0("s2z", "_center"),
    paste0("s2z", "_center_auto"),
    paste0("rho", "_s2z"),
    paste0("mean_rho", "_s2z"),
    paste0("sdlog", "_s2z"),
    paste0("sd", "_level"),
    paste0("z_sd", "_s2z"),
    paste0("relative_sd", "_s2z"),
    paste0("reference_sd", "_s2z"),
    "center[[:space:]]*=[[:space:]]*['\"](auto|fisher)['\"]",
    "scale[[:space:]]*=[[:space:]]*['\"]varying['\"]"
  )
  for (pattern in forbidden_patterns) {
    expect_false(any(grepl(pattern, audit_text)), info = pattern)
  }
})
