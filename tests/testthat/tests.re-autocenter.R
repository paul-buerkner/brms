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
  expect_error(
    call_pilot(autocenter_control(pilot_args = list(calculate_lp = FALSE))),
    "calculate_lp"
  )
  expect_error(
    call_pilot(autocenter_control(pilot_args = list(psis_resample = TRUE))),
    "psis_resample"
  )
  for (draws in c(0, -1, 1.5, NA_real_, Inf)) {
    expect_error(
      call_pilot(autocenter_control(pilot_args = list(draws = draws))),
      "draws.*integer|single numeric value"
    )
  }
})

test_that("Pathfinder centering diagnoses and weights raw precursor draws", {
  capture <- new.env(parent = emptyenv())
  pathfinder <- function(data, seed = NULL, init = NULL, num_paths,
                         show_messages, show_exceptions, draws,
                         single_path_draws, calculate_lp, psis_resample) {
    capture$data <- data
    capture$num_paths <- num_paths
    capture$draws <- draws
    capture$single_path_draws <- single_path_draws
    capture$calculate_lp <- calculate_lp
    capture$psis_resample <- psis_resample
    list(draws = function(variables, format) {
      capture$variables <- variables
      out <- cbind(
        rep(c(0.1, 0.3, 0.2), 200L),
        rep(c(0.8, 0.6, 0.7), 200L),
        0.1 * qnorm((seq_len(600L) - 0.5) / 600L),
        0
      )
      colnames(out) <- c(
        "rho_center_candidate_1[1,1]",
        "rho_center_candidate_1[2,1]", "lp__", "lp_approx__"
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
  expect_equal(final_data$rho_s2z_1, matrix(0.4, 2L, 1L))
  expect_identical(final_data$compute_rho_center_candidate_1, 0L)
  expect_identical(capture$num_paths, 4L)
  expect_identical(capture$draws, 200L)
  expect_identical(capture$single_path_draws, 200L)
  expect_true(capture$calculate_lp)
  expect_false(capture$psis_resample)
  expect_identical(
    capture$variables,
    c("rho_center_candidate_1", "lp__", "lp_approx__")
  )
  expect_true(all(summary$diagnostics$pareto_k < 0.7))
  expect_true(all(summary$diagnostics$psis_ess > 0))
  expect_equal(
    unname(unclass(summary$rho[["1"]])), matrix(c(0.2, 0.7), 2L, 1L)
  )
  expect_identical(
    dimnames(summary$rho[["1"]]), list(c("a", "b"), "Intercept")
  )
})

test_that("HMC centering uses a centered precursor without changing final data", {
  capture <- new.env(parent = emptyenv())
  sample <- function(data, seed, init, iter_sampling, iter_warmup, chains,
                     thin, parallel_chains, show_messages, show_exceptions) {
    capture$data <- data
    list(draws = function(variables, format) {
      matrix(c(0.2, 0.4, 0.6), ncol = 1L, dimnames = list(
        NULL, "rho_center_candidate_1[1,1]"
      ))
    })
  }
  final_data <- list(
    rho_s2z_1 = matrix(0.5, 1L, 1L),
    compute_rho_center_candidate_1 = 0L,
    rho_s2z_2 = matrix(0.25, 2L, 1L)
  )
  specs <- list(`1` = list(
    id = 1L, G = 1L, M = 1L, levels = "a", coefficients = "Intercept"
  ))
  summary <- brms:::run_re_autocenter_pilot(
    model = list(sample = sample),
    sdata = final_data, specs = specs,
    control = autocenter_control(method = "hmc"),
    backend = "cmdstanr", chains = 4L, cores = 1L, threads = NULL,
    opencl = NULL, init = "random", seed = NA, silent = 2L
  )
  expect_equal(capture$data$rho_s2z_1, matrix(1, 1L, 1L))
  expect_identical(capture$data$compute_rho_center_candidate_1, 1L)
  expect_identical(capture$data$rho_s2z_2, final_data$rho_s2z_2)
  expect_equal(final_data$rho_s2z_1, matrix(0.5, 1L, 1L))
  expect_identical(final_data$compute_rho_center_candidate_1, 0L)
  expect_equal(summary$rho[["1"]][1, 1], 0.4)
  expect_true(all(is.na(summary$diagnostics[c("pareto_k", "psis_ess")])))
})

test_that("Pathfinder posterior weighting is deterministic and preserves RNG", {
  x <- qnorm((seq_len(1000L) - 0.5) / 1000L)
  draws <- cbind(
    `rho_center_candidate_1[1,1]` = as.numeric(x > 0),
    lp__ = dnorm(x, mean = 0.6, log = TRUE),
    lp_approx__ = dnorm(x, log = TRUE)
  )
  set.seed(192)
  rng_before <- .Random.seed
  weighted <- brms:::.re_autocenter_pathfinder_draws(draws, ndraws = 2000L)
  expect_identical(.Random.seed, rng_before)
  expect_identical(
    brms:::.re_autocenter_pathfinder_draws(draws, ndraws = 2000L),
    weighted
  )
  expect_equal(nrow(weighted$draws), 2000L)
  # Tilting a standard normal by these log-density ratios yields N(0.6, 1).
  # Resampled centering candidates should reflect that posterior probability.
  expect_equal(mean(draws[, 1]), 0.5)
  expect_equal(mean(weighted$draws[, 1]), pnorm(0.6), tolerance = 0.01)
})

test_that("Pathfinder centering requires a finite Pareto k below 0.7", {
  x <- qnorm((seq_len(200L) - 0.5) / 200L)
  draws <- cbind(
    `rho_center_candidate_1[1,1]` = plogis(x),
    lp__ = 0.1 * x, lp_approx__ = 0
  )
  psis <- loo::psis(draws[, "lp__"] - draws[, "lp_approx__"], r_eff = 1)
  check_k <- function(k, accept = FALSE) {
    controlled_psis <- psis
    controlled_psis$diagnostics$pareto_k <- k
    local_mocked_bindings(
      psis = function(log_ratios, r_eff, ...) controlled_psis,
      .package = "loo"
    )
    if (accept) {
      out <- brms:::.re_autocenter_pathfinder_draws(draws, ndraws = 200L)
      expect_equal(out$pareto_k, k)
    } else {
      expect_error(
        brms:::.re_autocenter_pathfinder_draws(draws, ndraws = 200L),
        "center_control.*autocenter_control.*hmc"
      )
    }
  }
  check_k(0.699, accept = TRUE)
  for (k in list(0.7, 0.8, NA_real_, Inf, -Inf, numeric(), c(0.1, 0.2))) {
    check_k(k)
  }
})

test_that("Pathfinder zero-weight draws do not affect centering candidates", {
  x <- qnorm((seq_len(200L) - 0.5) / 200L)
  draws <- cbind(
    `rho_center_candidate_1[1,1]` = c(1, rep(0, 199L)),
    lp__ = c(-Inf, 0.1 * x[-1L]), lp_approx__ = 0
  )
  weighted <- brms:::.re_autocenter_pathfinder_draws(draws, ndraws = 200L)
  expect_true(all(weighted$draws[, 1] == 0))

  draws[, "lp__"] <- -Inf
  expect_error(
    brms:::.re_autocenter_pathfinder_draws(draws, ndraws = 200L),
    "center_control.*autocenter_control.*hmc"
  )
})

test_that("Pathfinder centering rejects unavailable or malformed densities", {
  draws <- cbind(
    `rho_center_candidate_1[1,1]` = seq(0.1, 0.9, length.out = 100L),
    lp__ = seq(-1, 1, length.out = 100L), lp_approx__ = 0
  )
  invalid <- list(
    draws[, c(1, 2), drop = FALSE],
    draws[, c(1, 3), drop = FALSE],
    draws[FALSE, , drop = FALSE]
  )
  for (value in c(NA_real_, Inf)) {
    malformed <- draws
    malformed[1L, "lp__"] <- value
    invalid[[length(invalid) + 1L]] <- malformed
  }
  for (malformed in invalid) {
    expect_error(
      brms:::.re_autocenter_pathfinder_draws(malformed, ndraws = 200L),
      "center_control.*autocenter_control.*hmc"
    )
  }
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
