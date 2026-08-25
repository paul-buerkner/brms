context("Tests for nonlinear sampled-loading S2Z Fisher analysis")

s2z_fisher_factor_data <- local({
  out <- expand.grid(
    person = factor(seq_len(4L)), item = factor(letters[1:3]),
    KEEP.OUT.ATTRS = FALSE
  )
  out$batch <- factor(rep(seq_len(3L), length.out = nrow(out)))
  out$x <- seq_len(nrow(out)) / nrow(out)
  out$y <- sin(seq_len(nrow(out)))
  out
})

s2z_fisher_factor_frame <- function(outer = alpha + loading * eta,
                                    loading = loading ~ 0 + item,
                                    family = gaussian()) {
  outer <- substitute(outer)
  loading <- substitute(loading)
  outer <- as.formula(call("~", quote(y), outer), env = parent.frame())
  loading <- as.formula(loading, env = parent.frame())
  form <- bf(
    outer,
    alpha ~ 0 + item,
    loading,
    eta ~ 1 + (1 | gr(
      person, s2z = TRUE, center = "auto", latent = TRUE
    )),
    nl = TRUE,
    family = family
  )
  brmsframe(brmsterms(form), data = s2z_fisher_factor_data)
}

test_that("Gaussian nonlinear loading derivatives are analyzed for codegen", {
  bframe <- s2z_fisher_factor_frame()
  got <- re_s2z_fisher_nl_info(bframe, bframe$nlpars$eta)

  expect_true(got$supported)
  expect_equal(got$id, 1)
  expect_identical(got$target_nlpar, "eta")
  expect_identical(got$derivative, quote(loading))
  expect_identical(got$dependencies, "loading")
  expect_identical(got$nlpar_dependencies, "loading")
  expect_length(got$covariate_dependencies, 0L)
  expect_true(got$is_gaussian_identity)
  expect_false(got$is_student_identity)
  expect_true(got$is_location_identity)
  expect_identical(got$response_family, "gaussian")
  expect_true(got$is_observation_local)
  expect_true(got$score_independent)
  expect_false(got$has_response_addition_terms)
  expect_length(got$response_addition_terms, 0L)
  expect_identical(got$sigma_kind, "scalar_parameter")
  expect_true(got$sigma_is_scalar)
  expect_null(got$sigma_fixed_value)
  expect_identical(got$sigma_parameter, "sigma")
  expect_identical(got$obs_derivative_stan, "fisher_s2z_1_nlp_loading[n]")
  expect_identical(got$dmu_deta_stan, "fisher_s2z_1_nlp_loading[n]")
  expect_identical(got$cn, 1L)
  expect_identical(got$frame_rows, "1")
  expect_identical(got$row_keys, "1:eta:1")
  expect_identical(
    got$row_metadata[c("row_key", "nlpar", "coef", "cn")],
    data.frame(
      row_key = "1:eta:1", nlpar = "eta", coef = "Intercept", cn = 1L,
      stringsAsFactors = FALSE
    )
  )

  loading <- got$dependency_info$loading
  expect_identical(loading$x_name, "X_loading")
  expect_identical(loading$coefficient_name, "b_loading")
  expect_identical(loading$vector_name, "fisher_s2z_1_nlp_loading")
  expect_identical(
    loading$vector_expression, "X_loading * b_loading"
  )
})

test_that("dependency aliases are isolated across Fisher factor IDs", {
  form <- bf(
    y ~ alpha + loading * eta,
    alpha ~ 0 + item,
    loading ~ 0 + item,
    eta ~ 1 +
      (1 | gr(
        person, id = "person_factor", s2z = TRUE,
        center = "auto", latent = TRUE
      )) +
      (1 | gr(
        batch, id = "batch_factor", s2z = TRUE,
        center = "auto", latent = TRUE
      )),
    nl = TRUE
  )
  bframe <- brmsframe(brmsterms(form), data = s2z_fisher_factor_data)
  ids <- unique(bframe$nlpars$eta$frame$re$id)
  expect_length(ids, 2L)
  expect_error(
    re_s2z_fisher_nl_info(bframe, bframe$nlpars$eta),
    "requires an explicit group-level ID"
  )

  got <- lapply(ids, function(id) {
    re_s2z_fisher_nl_info(bframe, bframe$nlpars$eta, id = id)
  })
  aliases <- vapply(
    got, function(x) x$dependency_info$loading$vector_name, character(1)
  )
  expect_identical(
    aliases, paste0("fisher_s2z_", ids, "_nlp_loading")
  )
  expect_length(unique(aliases), 2L)
})

test_that("row metadata aligns dimensions shared across nonlinear targets", {
  form <- bf(
    y ~ alpha + loading1 * eta1 + loading2 * eta2,
    alpha ~ 0 + item,
    loading1 ~ 0 + item,
    loading2 ~ 0 + item,
    eta1 ~ 1 + (1 | gr(
      person, id = "factor", s2z = TRUE,
      center = "auto", latent = TRUE
    )),
    eta2 ~ 1 + (1 | gr(
      person, id = "factor", s2z = TRUE,
      center = "auto", latent = TRUE
    )),
    nl = TRUE
  )
  bframe <- brmsframe(brmsterms(form), data = s2z_fisher_factor_data)
  r <- bframe$frame$re
  id <- unique(r$id)
  expect_length(id, 1L)

  got <- lapply(c("eta1", "eta2"), function(nlpar) {
    re_s2z_fisher_nl_info(
      bframe, bframe$nlpars[[nlpar]], id = id
    )
  })
  metadata <- do.call(rbind, lapply(got, `[[`, "row_metadata"))
  global_keys <- paste(r$id, r$nlpar, r$cn, sep = ":")
  aligned <- metadata[match(global_keys, metadata$row_key), , drop = FALSE]

  expect_false(anyNA(match(global_keys, metadata$row_key)))
  expect_identical(aligned$row_key, global_keys)
  expect_identical(aligned$cn, unname(r$cn))
  expect_identical(aligned$nlpar, unname(r$nlpar))
  expect_identical(vapply(got, `[[`, character(1), "target_nlpar"),
                   c("eta1", "eta2"))
  expect_identical(got[[1]]$dmu_deta_stan,
                   "fisher_s2z_1_nlp_loading1[n]")
  expect_identical(got[[2]]$dmu_deta_stan,
                   "fisher_s2z_1_nlp_loading2[n]")
})

test_that("nonlinear Fisher analysis maps safe observation covariates", {
  bframe <- s2z_fisher_factor_frame(
    outer = alpha + x * loading * eta
  )
  got <- re_s2z_fisher_nl_info(bframe, bframe$nlpars$eta)

  expect_setequal(got$dependencies, c("x", "loading"))
  expect_identical(got$covariate_dependencies, "x")
  expect_match(got$obs_derivative_stan, "C_1[n]", fixed = TRUE)
  expect_match(
    got$obs_derivative_stan, "fisher_s2z_1_nlp_loading[n]", fixed = TRUE
  )
})

test_that("nonlinear Fisher analysis rejects score-dependent derivatives", {
  bframe <- s2z_fisher_factor_frame(
    outer = alpha + loading * eta^2
  )
  expect_error(
    re_s2z_fisher_nl_info(bframe, bframe$nlpars$eta),
    "depends on target nonlinear parameter 'eta'"
  )

  soft <- re_s2z_fisher_nl_info(
    bframe, bframe$nlpars$eta, strict = FALSE
  )
  expect_false(soft$supported)
  expect_match(soft$reason_unsupported, "target nonlinear parameter 'eta'")
})

test_that("nonlinear Fisher analysis rejects latent loading predictors", {
  bframe <- s2z_fisher_factor_frame(
    loading = loading ~ 0 + item + (1 | batch)
  )
  expect_error(
    re_s2z_fisher_nl_info(bframe, bframe$nlpars$eta),
    "latent or group-varying nonlinear parameter 'loading'"
  )
})

test_that("nonlinear Fisher analysis rejects centered loading predictors", {
  form <- bf(
    y ~ alpha + loading * eta,
    alpha ~ 0 + item,
    lf(loading ~ item, center = TRUE),
    eta ~ 1 + (1 | gr(
      person, s2z = TRUE, center = "auto", latent = TRUE
    )),
    nl = TRUE
  )
  bframe <- brmsframe(
    brmsterms(form), data = s2z_fisher_factor_data
  )
  expect_error(
    re_s2z_fisher_nl_info(bframe, bframe$nlpars$eta),
    "centered population-level predictors.*'loading'"
  )
})

test_that("nonlinear Fisher analysis reports likelihood capability metadata", {
  weighted <- bf(
    y | weights(x) ~ alpha + loading * eta,
    alpha ~ 0 + item,
    loading ~ 0 + item,
    eta ~ 1 + (1 | gr(
      person, s2z = TRUE, center = "auto", latent = TRUE
    )),
    nl = TRUE
  )
  bframe <- brmsframe(
    brmsterms(weighted), data = s2z_fisher_factor_data
  )
  got <- re_s2z_fisher_nl_info(bframe, bframe$nlpars$eta)
  expect_true(got$has_response_addition_terms)
  expect_identical(got$response_addition_terms, "weights")

  predicted_sigma <- bf(
    y ~ alpha + loading * eta,
    alpha ~ 0 + item,
    loading ~ 0 + item,
    eta ~ 1 + (1 | gr(
      person, s2z = TRUE, center = "auto", latent = TRUE
    )),
    sigma ~ 0 + item,
    nl = TRUE
  )
  bframe <- brmsframe(
    brmsterms(predicted_sigma), data = s2z_fisher_factor_data
  )
  got <- re_s2z_fisher_nl_info(bframe, bframe$nlpars$eta)
  expect_identical(got$sigma_kind, "predictor")
  expect_false(got$sigma_is_scalar)

  student_frame <- s2z_fisher_factor_frame(family = student())
  got <- re_s2z_fisher_nl_info(
    student_frame, student_frame$nlpars$eta
  )
  expect_identical(got$response_family, "student")
  expect_false(got$is_gaussian_identity)
  expect_true(got$is_student_identity)
  expect_true(got$is_location_identity)
})

test_that("constant sampled loadings retain the ordinary b symbol", {
  form <- bf(
    y ~ alpha + loading * eta,
    alpha ~ 0 + item,
    loading ~ 0 + item,
    eta ~ 1 + (1 | gr(
      person, s2z = TRUE, center = FALSE, latent = TRUE
    )),
    nl = TRUE
  )
  scode <- stancode(
    form, data = s2z_fisher_factor_data,
    prior = prior(constant(1), class = "b", nlpar = "loading")
  )
  expect_match2(scode, "vector[K_loading] b_loading;")
  expect_match2(scode, "b_loading = rep_vector(1, rows(b_loading));")
  expect_match2(scode, "nlp_loading += X_loading * b_loading;")
})

test_that("nonlinear Fisher analysis rejects unsupported outer models", {
  bernoulli_data <- s2z_fisher_factor_data
  bernoulli_data$y <- rep(0:1, length.out = nrow(bernoulli_data))
  form <- bf(
    y ~ alpha + loading * eta,
    alpha ~ 0 + item,
    loading ~ 0 + item,
    eta ~ 1 + (1 | gr(
      person, s2z = TRUE, center = "auto", latent = TRUE
    )),
    nl = TRUE,
    family = bernoulli()
  )
  bframe <- brmsframe(brmsterms(form), data = bernoulli_data)
  expect_error(
    re_s2z_fisher_nl_info(bframe, bframe$nlpars$eta),
    "requires a Gaussian or Student-t response with identity link"
  )

  unknown <- bf(
    y ~ alpha + mystery(eta),
    alpha ~ 0 + item,
    eta ~ 1 + (1 | gr(
      person, s2z = TRUE, center = "auto", latent = TRUE
    )),
    nl = TRUE
  )
  bframe <- brmsframe(
    brmsterms(unknown), data = s2z_fisher_factor_data
  )
  expect_error(
    re_s2z_fisher_nl_info(bframe, bframe$nlpars$eta),
    "could not symbolically differentiate"
  )
})
