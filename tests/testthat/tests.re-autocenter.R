context("Tests for automatic group-effect centering precursors")

test_that("autocenter controls validate the shared precursor contract", {
  control <- autocenter_control()
  expect_s3_class(control, "brms_re_autocenter_control")
  expect_identical(control$method, "pathfinder")
  expect_identical(control$aggregate, "median")
  expect_identical(control$fallback, "error")
  expect_identical(control$pilot_args, list())

  control <- autocenter_control(
    method = "hmc", aggregate = "mean", fallback = 0.35,
    pilot_args = list(iter = 300L, warmup = 250L)
  )
  expect_identical(control$method, "hmc")
  expect_identical(control$aggregate, "mean")
  expect_equal(control$fallback, 0.35)
  expect_identical(control$pilot_args$iter, 300L)

  expect_error(
    brms:::validate_re_autocenter_control("pathfinder"),
    "must be a list"
  )
  expect_error(
    brms:::validate_re_autocenter_control(list(method = "laplace")),
    "pathfinder.*hmc"
  )
  expect_error(
    brms:::validate_re_autocenter_control(list(aggregate = "mode")),
    "median.*mean"
  )
  expect_error(
    brms:::validate_re_autocenter_control(list(fallback = 1.1)),
    "in \\[0, 1\\]"
  )
  expect_error(
    brms:::validate_re_autocenter_control(list(pilot_args = list(300L))),
    "named list"
  )
  expect_error(
    brms:::validate_re_autocenter_control(list(unknown = TRUE)),
    "Invalid.*unknown"
  )
})

test_that("method-specific precursor arguments fail before dispatch", {
  call_pilot <- function(control) {
    brms:::run_re_autocenter_pilot(
      model = list(), sdata = list(),
      specs = list(`1` = list(id = 1L, G = 1L, M = 1L)),
      control = control, backend = "cmdstanr", chains = 1L, cores = 1L,
      threads = NULL, opencl = NULL, init = "random", seed = NA,
      silent = 2L
    )
  }

  expect_error(
    call_pilot(autocenter_control(
      method = "hmc", pilot_args = list(num_paths = 2L)
    )),
    "HMC.*cannot contain 'num_paths'"
  )
  expect_error(
    call_pilot(autocenter_control(
      method = "pathfinder", pilot_args = list(control = list(adapt_delta = 0.9))
    )),
    "available only for method.*hmc"
  )
  expect_error(
    call_pilot(autocenter_control(
      method = "hmc", pilot_args = list(iter = 50L, warmup = 50L)
    )),
    "warmup.*smaller than.*iter"
  )
})

test_that("Pathfinder dispatch flips only precursor data and reads raw draws", {
  capture <- new.env(parent = emptyenv())
  pathfinder <- function(data, seed = NULL, init = NULL, num_paths,
                         show_messages, show_exceptions, draws,
                         single_path_draws) {
    capture$data <- data
    capture$num_paths <- num_paths
    capture$draws <- draws
    capture$single_path_draws <- single_path_draws
    list(draws = function(variables, format) {
      out <- matrix(c(0.1, 0.3, 0.2, 0.8, 0.6, 0.7), 3L, 2L)
      colnames(out) <- c(
        "rho_center_candidate_1[1,1]",
        "rho_center_candidate_1[2,1]"
      )
      out
    })
  }
  specs <- list(`1` = list(
    id = 1L, G = 2L, M = 1L, levels = c("a", "b"),
    coefficients = "Intercept"
  ))
  final_data <- list(
    rho_s2z_1 = matrix(0.4, 2L, 1L),
    compute_rho_center_candidate_1 = 0L
  )
  summary <- brms:::run_re_autocenter_pilot(
    model = list(pathfinder = pathfinder), sdata = final_data,
    specs = specs, control = autocenter_control(), backend = "cmdstanr",
    chains = 4L, cores = 1L, threads = NULL, opencl = NULL,
    init = "random", seed = NA, silent = 2L
  )

  expect_equal(capture$data$rho_s2z_1, matrix(0, 2L, 1L))
  expect_identical(capture$data$compute_rho_center_candidate_1, 1L)
  expect_identical(final_data$compute_rho_center_candidate_1, 0L)
  expect_identical(capture$num_paths, 4L)
  expect_identical(capture$draws, 200L)
  expect_identical(capture$single_path_draws, 200L)
  expect_equal(
    unname(unclass(summary$rho[["1"]])), matrix(c(0.2, 0.7), 2L, 1L)
  )
  expect_identical(
    dimnames(summary$rho[["1"]]), list(c("a", "b"), "Intercept")
  )
})

test_that("candidate draws resolve to bounded named matrices", {
  draws <- cbind(
    unrelated = seq_len(4),
    a22 = c(NA, Inf, 1.2, -0.2),
    a11 = c(0.1, 0.3, 0.2, 0.8),
    b11 = c(0.8, 0.6, 0.7, 0.9),
    a12 = c(-1e-10, 0, 0.1, 1.4),
    a21 = c(0, 1, NA, Inf)
  )
  colnames(draws) <- c(
    "lp__",
    "rho_center_candidate_2[2,2]",
    "rho_center_candidate_2[1,1]",
    "rho_center_candidate_9[1,1]",
    "rho_center_candidate_2[1,2]",
    "rho_center_candidate_2[2,1]"
  )
  summary <- brms:::aggregate_re_autocenter_draws(
    draws, autocenter_control(fallback = 0.35)
  )

  expect_s3_class(summary, "brms_re_autocenter_summary")
  expect_setequal(names(summary$rho), c("2", "9"))
  expect_s3_class(summary$rho[["2"]], "brmsautocenter_resolved")
  expect_identical(dim(summary$rho[["2"]]), c(2L, 2L))
  expect_identical(
    dimnames(summary$rho[["2"]]),
    list(c("1", "2"), c("1", "2"))
  )
  expect_equal(
    unclass(summary$rho[["2"]]),
    matrix(c(0.25, 0.5, 0, 0.35), 2L, 2L,
           dimnames = list(c("1", "2"), c("1", "2")))
  )
  expect_equal(
    unname(unclass(summary$rho[["9"]])), matrix(0.75, 1L, 1L)
  )
  expect_identical(centering_weights(summary), summary$rho)

  fallback_row <- subset(
    summary$diagnostics,
    variable == "rho_center_candidate_2[2,2]"
  )
  expect_true(fallback_row$fallback_used)
  expect_equal(fallback_row$n_valid, 0L)
  expect_equal(fallback_row$n_out_of_bounds, 2L)

  mean_summary <- brms:::aggregate_re_autocenter_draws(
    draws, autocenter_control(aggregate = "mean", fallback = 0.35)
  )
  expect_equal(mean_summary$rho[["2"]][1, 1], 0.35)

  fit <- structure(
    list(autocenter = list(weights = summary$rho)), class = "brmsfit"
  )
  expect_identical(centering_weights(fit), summary$rho)
  multiple <- structure(
    list(autocenter = list(
      weights = summary$rho,
      by_fit = list(
        list(weights = summary$rho),
        list(weights = mean_summary$rho)
      )
    )),
    class = c("brmsfit_multiple", "brmsfit")
  )
  expect_identical(
    centering_weights(multiple),
    list(summary$rho, mean_summary$rho)
  )
  expect_error(
    centering_weights(structure(list(), class = "brmsfit")),
    "does not contain precursor centering weights"
  )
})

test_that("candidate draw aggregation diagnoses malformed output", {
  no_valid <- matrix(
    c(NA_real_, Inf), ncol = 1L,
    dimnames = list(NULL, "rho_center_candidate_1[1,1]")
  )
  expect_error(
    brms:::aggregate_re_autocenter_draws(no_valid),
    "No valid precursor draws"
  )
  expect_equal(
    brms:::aggregate_re_autocenter_draws(
      no_valid, autocenter_control(fallback = "centered")
    )$rho[["1"]][1, 1],
    1
  )

  incomplete <- matrix(
    runif(6), nrow = 2L,
    dimnames = list(NULL, c(
      "rho_center_candidate_1[1,1]",
      "rho_center_candidate_1[2,1]",
      "rho_center_candidate_1[1,2]"
    ))
  )
  expect_error(
    brms:::aggregate_re_autocenter_draws(incomplete),
    "complete G by K matrix"
  )
  expect_error(
    brms:::aggregate_re_autocenter_draws(
      matrix(1, 1L, 1L, dimnames = list(NULL, "lp__"))
    ),
    "No variables named"
  )
})

test_that("a resolved 1 by 1 matrix is not collapsed to a scalar", {
  resolved <- brms:::.as_brmsautocenter_resolved(matrix(
    0.42, 1L, 1L, dimnames = list("only-level", "Intercept")
  ))
  call <- gr(g, center = resolved)

  expect_true(call$s2z_center_auto)
  expect_equal(call$s2z_center, 0.5)
  expect_s3_class(call$s2z_center_data, "brmsautocenter_resolved")
  expect_identical(dim(call$s2z_center_data), c(1L, 1L))
  expect_equal(
    unname(unclass(call$s2z_center_data)), matrix(0.42, 1L, 1L)
  )
})

test_that("formula replacement freezes resolved matrices by framed occurrence", {
  source_env <- new.env(parent = globalenv())
  source_env$center_mode <- "auto"
  source_env$environment_sentinel <- 17L
  form <- as.formula(
    paste0(
      "y ~ (1 | gr(g, id = 'shared', center = center_mode)) + ",
      "(0 + x | brms::gr(g, id = 'shared', center = 'auto')) + ",
      "(1 | gr(h, center = 0.25))"
    ),
    env = source_env
  )
  first <- matrix(
    c(0.1, 0.2), 2L, 1L,
    dimnames = list(c("a", "b"), "Intercept")
  )
  second <- matrix(
    c(0.7, 0.8), 2L, 1L,
    dimnames = list(c("a", "b"), "x")
  )
  resolved <- brms:::replace_re_autocenter(
    form, center = list(first, second), id = c(3L, 3L),
    occurrence = c(1L, 2L)
  )
  resolved_text <- paste(deparse(resolved), collapse = " ")
  first_symbol <- ".brms_autocenter_id_3_occurrence_1"
  second_symbol <- ".brms_autocenter_id_3_occurrence_2"
  expect_match(resolved_text, first_symbol, fixed = TRUE)
  expect_match(resolved_text, second_symbol, fixed = TRUE)
  expect_match(resolved_text, "center = 0.25", fixed = TRUE)
  expect_identical(
    get("environment_sentinel", envir = environment(resolved)), 17L
  )
  expect_s3_class(
    get(first_symbol, envir = environment(resolved), inherits = FALSE),
    "brmsautocenter_resolved"
  )
  expect_equal(
    unname(unclass(get(second_symbol, envir = environment(resolved),
                       inherits = FALSE))),
    unname(second)
  )
  expect_match(paste(deparse(form), collapse = " "), "center_mode")

  expect_error(
    brms:::replace_re_autocenter(
      form, center = list(first), id = 3L, occurrence = 1L
    ),
    "More gr.*occurrences"
  )
  expect_error(
    brms:::replace_re_autocenter(
      form, center = list(first, second, first), id = c(3L, 3L, 4L),
      occurrence = c(1L, 2L, 1L)
    ),
    "Fewer gr.*occurrences"
  )
})

test_that("deduplicated source group terms reuse one resolved matrix", {
  dat <- data.frame(
    y = seq_len(6) / 10,
    g = factor(rep(c("a", "b"), each = 3))
  )
  form <- y ~
    (1 | gr(g, center = "auto")) +
    (1 | gr(g, center = "auto"))
  bframe <- brms:::brmsframe(brmsterms(form), dat)
  expect_equal(nrow(bframe$frame$re), 1L)
  rho <- matrix(
    c(0.25, 0.75), 2L, 1L,
    dimnames = list(levels(dat$g), "Intercept")
  )
  resolved <- brms:::replace_re_autocenter(
    form, center = list(rho), id = 1L, occurrence = 1L
  )
  resolved_text <- paste(deparse(resolved), collapse = " ")
  expect_equal(
    lengths(regmatches(
      resolved_text,
      gregexpr(".brms_autocenter_id_1_occurrence_1",
               resolved_text, fixed = TRUE)
    )),
    2L
  )
  expect_equal(
    unname(standata(resolved, data = dat)$rho_s2z_1),
    unname(rho)
  )
})

test_that("refits do not shadow resolved symbols when framed IDs shift", {
  dat <- data.frame(
    y = seq_len(12) / 10,
    z = factor(rep(c("z1", "z2"), 6)),
    a = factor(rep(c("a1", "a2", "a3"), each = 4))
  )
  rho_z <- matrix(
    c(0.2, 0.4), 2L, 1L,
    dimnames = list(levels(dat$z), "Intercept")
  )
  frozen <- brms:::replace_re_autocenter(
    y ~ (1 | gr(z, center = "auto")),
    center = list(rho_z), id = 1L, occurrence = 1L
  )
  expanded <- update(
    frozen, . ~ . + (1 | gr(a, center = "auto"))
  )
  bframe <- brms:::brmsframe(brmsterms(expanded), dat)
  expect_identical(names(brms:::re_autocenter_specs(bframe)), "1")
  rho_a <- brms:::.as_brmsautocenter_resolved(matrix(
    c(0.6, 0.7, 0.8), 3L, 1L,
    dimnames = list(levels(dat$a), "Intercept")
  ))
  occurrences <- brms:::re_autocenter_occurrence_weights(
    bframe, weights = setNames(list(rho_a), "1"), formula = expanded
  )
  resolved <- brms:::replace_re_autocenter(
    expanded, center = occurrences$center, id = occurrences$id,
    occurrence = occurrences$occurrence
  )
  resolved_text <- paste(deparse(resolved), collapse = " ")
  expect_match(resolved_text, "_refit_2", fixed = TRUE)
  sdata <- standata(resolved, data = dat)
  expect_equal(unname(sdata$rho_s2z_1), unname(unclass(rho_a)))
  expect_equal(unname(sdata$rho_s2z_2), unname(rho_z))
})

test_that("automatic frame specs preserve covariance and occurrence maps", {
  dat <- data.frame(
    y = seq_len(12) / 10,
    x = rep(c(-1, 0, 1), 4),
    g = factor(rep(letters[1:3], each = 4))
  )
  form <- y ~ x +
    (1 | gr(g, id = "shared", center = "auto")) +
    (0 + x | gr(g, id = "shared", center = "auto"))
  bframe <- brms:::brmsframe(brmsterms(form), dat)
  specs <- brms:::re_autocenter_specs(bframe)
  expect_identical(names(specs), "1")
  expect_identical(specs[[1]]$G, 3L)
  expect_identical(specs[[1]]$M, 2L)
  expect_identical(specs[[1]]$coefficients, c("Intercept", "x"))
  expect_identical(specs[[1]]$dimension, c(1L, 2L))

  rho <- matrix(
    seq(0.1, 0.6, length.out = 6), 3L, 2L,
    dimnames = list(levels(dat$g), c("Intercept", "x"))
  )
  occurrences <- brms:::re_autocenter_occurrence_weights(
    bframe, list(`1` = brms:::.as_brmsautocenter_resolved(rho)),
    formula = form
  )
  expect_identical(occurrences$id, c(1L, 1L))
  expect_identical(occurrences$occurrence, c(1L, 2L))
  expect_identical(colnames(occurrences$center[[1]]), "Intercept")
  expect_identical(colnames(occurrences$center[[2]]), "x")
  expect_equal(
    unclass(occurrences$center[[1]]), rho[, 1, drop = FALSE]
  )
  expect_equal(
    unclass(occurrences$center[[2]]), rho[, 2, drop = FALSE]
  )
})

test_that("occurrence weights follow source rather than framed group order", {
  dat <- data.frame(
    y = seq_len(12) / 10,
    z = factor(rep(c("z1", "z2"), 6)),
    a = factor(rep(c("a1", "a2", "a3"), each = 4))
  )
  form <- y ~
    (1 | gr(z, center = "auto")) +
    (1 | gr(a, center = "auto"))
  bframe <- brms:::brmsframe(brmsterms(form), dat)
  expect_identical(bframe$frame$re$group, c("a", "z"))
  expect_equal(bframe$frame$re$id, c(1L, 2L))

  rho_a <- brms:::.as_brmsautocenter_resolved(matrix(
    c(0.11, 0.22, 0.33), 3L, 1L,
    dimnames = list(levels(dat$a), "Intercept")
  ))
  rho_z <- brms:::.as_brmsautocenter_resolved(matrix(
    c(0.77, 0.88), 2L, 1L,
    dimnames = list(levels(dat$z), "Intercept")
  ))
  occurrences <- brms:::re_autocenter_occurrence_weights(
    bframe, weights = list(`1` = rho_a, `2` = rho_z), formula = form
  )
  expect_identical(occurrences$id, c(2L, 1L))
  expect_equal(unclass(occurrences$center[[1L]]), unclass(rho_z))
  expect_equal(unclass(occurrences$center[[2L]]), unclass(rho_a))

  resolved <- brms:::replace_re_autocenter(
    form, center = occurrences$center, id = occurrences$id,
    occurrence = occurrences$occurrence
  )
  sdata <- standata(resolved, data = dat)
  expect_equal(unname(sdata$rho_s2z_1), unname(unclass(rho_a)))
  expect_equal(unname(sdata$rho_s2z_2), unname(unclass(rho_z)))
})

test_that("explicit mu formulas retain source component order", {
  dat <- data.frame(
    y = seq_len(12) / 10,
    g_mu = factor(rep(c("m1", "m2"), 6)),
    g_sigma = factor(rep(c("s1", "s2", "s3"), each = 4))
  )
  form <- bf(
    y ~ 1,
    sigma ~ 1 + (1 | gr(g_sigma, center = "auto")),
    mu ~ 1 + (1 | gr(g_mu, center = "auto")),
    family = gaussian()
  )
  bframe <- brms:::brmsframe(brmsterms(form), dat)
  expect_identical(bframe$frame$re$group, c("g_mu", "g_sigma"))

  rho_mu <- brms:::.as_brmsautocenter_resolved(matrix(
    c(0.2, 0.4), 2L, 1L,
    dimnames = list(levels(dat$g_mu), "Intercept")
  ))
  rho_sigma <- brms:::.as_brmsautocenter_resolved(matrix(
    c(0.6, 0.7, 0.8), 3L, 1L,
    dimnames = list(levels(dat$g_sigma), "Intercept")
  ))
  occurrences <- brms:::re_autocenter_occurrence_weights(
    bframe,
    weights = setNames(list(rho_mu, rho_sigma), c("1", "2")),
    formula = form
  )
  expect_identical(occurrences$id, c(2L, 1L))
  expect_equal(unclass(occurrences$center[[1L]]), unclass(rho_sigma))
  expect_equal(unclass(occurrences$center[[2L]]), unclass(rho_mu))
})

test_that("multicategory shorthand shares aggregated automatic weights", {
  dat <- data.frame(
    y = factor(rep(c("a", "b", "c"), 4)),
    g = factor(rep(letters[1:4], each = 3))
  )
  form <- brms:::validate_formula(
    bf(
      y ~ 1 + (1 | gr(g, s2z = TRUE, center = "auto")),
      family = categorical()
    ),
    data = dat
  )
  bframe <- brms:::brmsframe(brmsterms(form), dat)
  specs <- brms:::re_autocenter_specs(bframe)
  expect_length(specs, 2L)
  rho_1 <- brms:::.as_brmsautocenter_resolved(matrix(
    c(0.1, 0.2, 0.3, 0.4), 4L, 1L,
    dimnames = list(levels(dat$g), "Intercept")
  ))
  rho_2 <- brms:::.as_brmsautocenter_resolved(matrix(
    c(0.5, 0.6, 0.7, 0.8), 4L, 1L,
    dimnames = list(levels(dat$g), "Intercept")
  ))
  occurrences <- brms:::re_autocenter_occurrence_weights(
    bframe,
    weights = setNames(list(rho_1, rho_2), names(specs)),
    formula = form,
    aggregate = "median"
  )
  shared <- (unclass(rho_1) + unclass(rho_2)) / 2
  expect_length(occurrences$center, 1L)
  expect_equal(unclass(occurrences$center[[1L]]), shared)

  resolved <- brms:::replace_re_autocenter(
    form, center = occurrences$center, id = occurrences$id,
    occurrence = occurrences$occurrence
  )
  sdata <- standata(resolved, data = dat)
  expect_equal(unname(sdata$rho_s2z_1), unname(shared))
  expect_equal(unname(sdata$rho_s2z_2), unname(shared))
})

test_that("partly resolved automatic IDs fail instead of filling zero weights", {
  dat <- data.frame(
    y = seq_len(12) / 10,
    x = rep(c(-1, 0, 1), 4),
    g = factor(rep(letters[1:3], each = 4))
  )
  fixed <- brms:::.as_brmsautocenter_resolved(matrix(
    c(0.2, 0.3, 0.4), 3L, 1L,
    dimnames = list(levels(dat$g), "Intercept")
  ))
  form <- y ~
    (1 | gr(g, id = "shared", center = fixed)) +
    (0 + x | gr(g, id = "shared", center = "auto"))
  bframe <- brms:::brmsframe(brmsterms(form), dat)
  expect_error(
    brms:::re_autocenter_specs(bframe),
    "mixes resolved and unresolved"
  )
})

test_that("resolved automatic weights are recoverable after reframing", {
  dat <- data.frame(
    y = seq_len(9) / 10,
    g = factor(rep(letters[1:3], each = 3))
  )
  fixed <- brms:::.as_brmsautocenter_resolved(matrix(
    c(0.15, 0.45, 0.75), 3L, 1L,
    dimnames = list(levels(dat$g), "Intercept")
  ))
  form <- y ~ (1 | gr(g, center = fixed))
  bframe <- brms:::brmsframe(brmsterms(form), dat)
  recovered <- brms:::resolved_re_autocenter_weights(bframe)
  expect_identical(names(recovered), "1")
  expect_s3_class(recovered[[1L]], "brmsautocenter_resolved")
  expect_equal(unclass(recovered[[1L]]), unclass(fixed))

  subset_dat <- droplevels(dat[dat$g != "c", , drop = FALSE])
  subset_frame <- brms:::brmsframe(brmsterms(form), subset_dat)
  subset_weights <- brms:::resolved_re_autocenter_weights(subset_frame)
  expect_equal(
    unclass(subset_weights[[1L]]),
    unclass(fixed[c("a", "b"), , drop = FALSE])
  )
})

test_that("strict latent specs deduplicate response aliases", {
  dat <- data.frame(
    person = factor(seq_len(6)),
    y1 = sin(seq_len(6)),
    y2 = cos(seq_len(6))
  )
  response_factor <- function(response) {
    bf(
      as.formula(paste0(response, " ~ loading * eta")),
      lf(loading ~ 1, center = FALSE),
      eta ~ 0 + (1 | gr(
        person, id = "score", s2z = TRUE, latent = TRUE,
        center = "auto"
      )),
      nl = TRUE
    )
  }
  form <- response_factor("y1") + response_factor("y2") +
    set_rescor(FALSE)
  bframe <- brms:::brmsframe(brmsterms(form), dat)
  specs <- brms:::re_autocenter_specs(bframe)
  expect_length(specs, 1L)
  expect_true(specs[[1]]$latent)
  expect_identical(specs[[1]]$M, 1L)
  expect_identical(specs[[1]]$dimension, c(1L, 1L))
  expect_identical(specs[[1]]$coefficients, "eta:Intercept")

  rho <- matrix(
    seq(0.1, 0.6, length.out = 6), 6L, 1L,
    dimnames = list(levels(dat$person), "eta:Intercept")
  )
  occurrences <- brms:::re_autocenter_occurrence_weights(
    bframe, list(`1` = brms:::.as_brmsautocenter_resolved(rho)),
    formula = form
  )
  expect_length(occurrences$center, 2L)
  expect_true(all(vapply(
    occurrences$center, inherits, logical(1), "brmsautocenter_resolved"
  )))
  expect_true(all(vapply(
    occurrences$center, function(x) {
      isTRUE(all.equal(unname(x), unname(rho), check.attributes = FALSE))
    },
    logical(1)
  )))
})
