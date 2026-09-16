# Sampling validation for physical sum-to-zero group-level effects.
#
# This is intentionally a local test: it compiles and samples twelve Stan
# models, including literal-scalar, independent-component, and correlated
# paths.
# Run it against an installed development build of brms, for example with
#
#   Rscript tests/local/tests.models-s2z.R
#
# Useful controls (shown with their defaults) are
#
#   BRMS_S2Z_CHAINS=2
#   BRMS_S2Z_WARMUP=500
#   BRMS_S2Z_SAMPLING=500
#   BRMS_S2Z_CORES=min(chains, parallel::detectCores())
#   BRMS_S2Z_BACKEND=auto
#   BRMS_S2Z_ADAPT_DELTA=0.97
#   BRMS_S2Z_MAX_TREEDEPTH=12
#   BRMS_S2Z_EQUIV_Z=7
#   BRMS_S2Z_MAX_RHAT=1.10
#   BRMS_S2Z_MAX_DIVERGENCES=0
#   BRMS_S2Z_INVARIANT_TOL=1e-6
#   BRMS_S2Z_CACHE_DIR=
#
# Equivalence is assessed on both posterior means and posterior standard
# deviations. The Monte Carlo standard error from each independent fit is
# propagated as sqrt(mcse_conventional^2 + mcse_s2z^2).

suppressPackageStartupMessages({
  library(brms)
  library(posterior)
})

env_integer <- function(name, default, minimum = 1L) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(as.integer(default))
  }
  out <- suppressWarnings(as.integer(value))
  if (length(out) != 1L || is.na(out) || out < minimum) {
    stop(name, " must be an integer >= ", minimum, call. = FALSE)
  }
  out
}

env_number <- function(name, default, minimum = -Inf) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(as.numeric(default))
  }
  out <- suppressWarnings(as.numeric(value))
  if (length(out) != 1L || !is.finite(out) || out < minimum) {
    stop(name, " must be a finite number >= ", minimum, call. = FALSE)
  }
  out
}

choose_backend <- function() {
  requested <- tolower(Sys.getenv("BRMS_S2Z_BACKEND", unset = "auto"))
  if (requested %in% c("cmdstanr", "rstan")) {
    return(requested)
  }
  if (!identical(requested, "auto")) {
    stop("BRMS_S2Z_BACKEND must be auto, cmdstanr, or rstan.", call. = FALSE)
  }
  cmdstan_ok <- requireNamespace("cmdstanr", quietly = TRUE) &&
    isTRUE(tryCatch(
      nzchar(suppressWarnings(cmdstanr::cmdstan_path())),
      error = function(e) FALSE
    ))
  if (cmdstan_ok) {
    return("cmdstanr")
  }
  if (requireNamespace("rstan", quietly = TRUE)) {
    return("rstan")
  }
  stop("Neither a working CmdStan installation nor rstan was found.",
       call. = FALSE)
}

chains <- env_integer("BRMS_S2Z_CHAINS", 2L, minimum = 2L)
warmup <- env_integer("BRMS_S2Z_WARMUP", 500L)
sampling <- env_integer("BRMS_S2Z_SAMPLING", 500L)
detected_cores <- parallel::detectCores(logical = FALSE)
if (is.na(detected_cores)) {
  detected_cores <- 1L
}
cores <- env_integer(
  "BRMS_S2Z_CORES", min(chains, detected_cores), minimum = 1L
)
backend <- choose_backend()
adapt_delta <- env_number("BRMS_S2Z_ADAPT_DELTA", 0.97, minimum = 0.5)
max_treedepth <- env_integer("BRMS_S2Z_MAX_TREEDEPTH", 12L)
equiv_z <- env_number("BRMS_S2Z_EQUIV_Z", 7, minimum = 1)
max_rhat <- env_number("BRMS_S2Z_MAX_RHAT", 1.10, minimum = 1)
max_divergences <- env_integer(
  "BRMS_S2Z_MAX_DIVERGENCES", 0L, minimum = 0L
)
invariant_tolerance <- env_number(
  # CmdStan's default CSV precision introduces errors around 1e-7 when the
  # exactly equivalent generated quantities are read back into R.
  "BRMS_S2Z_INVARIANT_TOL", 1e-6, minimum = 0
)
cache_dir <- Sys.getenv("BRMS_S2Z_CACHE_DIR", unset = "")

options(mc.cores = cores, width = 140)

configuration <- data.frame(
  brms_version = as.character(utils::packageVersion("brms")),
  backend = backend,
  chains = chains,
  warmup = warmup,
  sampling = sampling,
  cores = cores,
  adapt_delta = adapt_delta,
  max_treedepth = max_treedepth,
  equiv_z = equiv_z,
  stringsAsFactors = FALSE
)
cat("S2Z sampling-equivalence configuration\n")
print(configuration, row.names = FALSE)

rmvn_rows <- function(n, sigma) {
  matrix(stats::rnorm(n * nrow(sigma)), nrow = n) %*% chol(sigma)
}

simulate_gaussian_interaction <- function(seed = 4101L, groups = 10L,
                                          per_group = 12L) {
  set.seed(seed)
  g_levels <- sprintf("g%02d", seq_len(groups))
  factorial_grid <- expand.grid(x = c(-1, 1), z = c(-1, 1))
  dat <- do.call(rbind, lapply(seq_len(groups), function(j) {
    grid <- factorial_grid[rep(seq_len(nrow(factorial_grid)),
                              length.out = per_group), ]
    data.frame(
      g = g_levels[j],
      x = grid$x + stats::rnorm(per_group, sd = 0.12),
      z = grid$z + stats::rnorm(per_group, sd = 0.12)
    )
  }))
  rownames(dat) <- NULL
  dat$g <- factor(dat$g, levels = g_levels)
  X <- stats::model.matrix(~ x * z, dat)
  beta <- c(0.45, 0.80, -0.55, 0.40)
  re_sd <- c(0.70, 0.45, 0.40, 0.30)
  re_cor <- outer(seq_along(re_sd), seq_along(re_sd), function(i, j) {
    0.30^abs(i - j)
  })
  sigma_re <- diag(re_sd) %*% re_cor %*% diag(re_sd)
  u <- rmvn_rows(groups, sigma_re)
  eta <- drop(X %*% beta) + rowSums(X * u[as.integer(dat$g), ])
  dat$y <- stats::rnorm(nrow(dat), eta, 0.45)
  dat
}

simulate_student_group <- function(seed = 4102L, groups = 12L,
                                   per_group = 10L, group_df = 5) {
  set.seed(seed)
  g_levels <- sprintf("g%02d", seq_len(groups))
  dat <- do.call(rbind, lapply(seq_len(groups), function(j) {
    data.frame(
      g = g_levels[j],
      x = seq(-1.25, 1.25, length.out = per_group) +
        stats::rnorm(per_group, sd = 0.08)
    )
  }))
  rownames(dat) <- NULL
  dat$g <- factor(dat$g, levels = g_levels)
  X <- stats::model.matrix(~ x, dat)
  beta <- c(-0.25, 0.75)
  re_sd <- c(0.80, 0.50)
  re_cor <- matrix(c(1, -0.25, -0.25, 1), nrow = 2L)
  sigma_re <- diag(re_sd) %*% re_cor %*% diag(re_sd)
  # One mixing scale per group gives a multivariate Student-t effect.
  u <- rmvn_rows(groups, sigma_re)
  u <- u * sqrt(group_df / stats::rchisq(groups, df = group_df))
  eta <- drop(X %*% beta) + rowSums(X * u[as.integer(dat$g), ])
  dat$y <- stats::rnorm(nrow(dat), eta, 0.55)
  dat
}

simulate_scalar_intercept <- function(seed = 4103L, groups = 14L,
                                      per_group = 10L) {
  set.seed(seed)
  g_levels <- sprintf("g%02d", seq_len(groups))
  dat <- data.frame(
    g = factor(
      rep(g_levels, each = per_group), levels = g_levels
    )
  )
  group_effect <- stats::rnorm(groups, sd = 0.75)
  eta <- 0.4 + group_effect[as.integer(dat$g)]
  dat$y <- stats::rnorm(nrow(dat), eta, 0.45)
  dat
}

simulate_scalar_student_slope <- function(seed = 4104L, groups = 14L,
                                          per_group = 10L,
                                          group_df = 5) {
  set.seed(seed)
  g_levels <- sprintf("g%02d", seq_len(groups))
  dat <- do.call(rbind, lapply(seq_len(groups), function(j) {
    data.frame(
      g = g_levels[j],
      x = seq(-1.5, 1.5, length.out = per_group) +
        stats::rnorm(per_group, sd = 0.08)
    )
  }))
  rownames(dat) <- NULL
  dat$g <- factor(dat$g, levels = g_levels)
  group_effect <- stats::rnorm(groups, sd = 0.55) *
    sqrt(group_df / stats::rchisq(groups, df = group_df))
  eta <- dat$x * (0.8 + group_effect[as.integer(dat$g)])
  dat$y <- stats::rnorm(nrow(dat), eta, 0.5)
  dat
}

model_prior <- function(student_group = FALSE, correlated = TRUE) {
  out <- prior(normal(0, 1.5), class = "b") +
    prior(normal(0, 1.5), class = "Intercept") +
    prior(exponential(1), class = "sd") +
    prior(exponential(1), class = "sigma")
  if (correlated) {
    out <- out + prior(lkj(2), class = "cor")
  }
  if (student_group) {
    out <- out + prior(gamma(2, 0.2), class = "df", group = "g")
  }
  out
}

scalar_prior <- function(intercept = TRUE, student_group = FALSE) {
  out <- if (intercept) {
    prior(normal(0, 1.5), class = "Intercept")
  } else {
    prior(normal(0, 1.5), class = "b", coef = "x")
  }
  out <- out + prior(exponential(1), class = "sd") +
    prior(exponential(1), class = "sigma")
  if (student_group) {
    out <- out + prior(gamma(2, 0.2), class = "df", group = "g")
  }
  out
}

cache_file <- function(label) {
  if (!nzchar(cache_dir)) {
    return(NULL)
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  suffix <- paste(backend, chains, warmup, sampling, sep = "-")
  file.path(cache_dir, paste0(label, "-", suffix, ".rds"))
}

sample_model <- function(formula, data, prior, seed, label) {
  cat("Sampling ", label, " ...\n", sep = "")
  args <- list(
    formula = formula,
    data = data,
    family = gaussian(),
    prior = prior,
    backend = backend,
    chains = chains,
    cores = cores,
    iter = warmup + sampling,
    warmup = warmup,
    seed = seed,
    refresh = 0,
    silent = 1,
    save_pars = save_pars(all = TRUE),
    control = list(
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth
    )
  )
  file <- cache_file(label)
  if (!is.null(file)) {
    args$file <- file
    args$file_refit <- "on_change"
  }
  do.call(brm, args)
}

exact_draw_matrix <- function(fit, variable) {
  missing <- setdiff(variable, variables(fit))
  if (length(missing)) {
    stop("Required saved variables are missing: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  as.matrix(as_draws_matrix(fit, variable = variable))
}

public_fixef_draws <- function(fit) {
  out <- as.matrix(fixef(fit, summary = FALSE))
  # Keep downstream comparisons interpretable even when this test detects that
  # an internal finite-population vector has accidentally leaked through
  # fixef().
  out[, !grepl("^s2z(\\[|$)", colnames(out)), drop = FALSE]
}

public_samples <- function(fit, data, group = "g") {
  fe <- public_fixef_draws(fit)
  re <- ranef(fit, summary = FALSE)[[group]]
  if (is.null(re)) {
    stop("No public group-level draws found for ", group, call. = FALSE)
  }
  out <- setNames(lapply(seq_len(ncol(fe)), function(k) fe[, k]),
                  paste0("b:", colnames(fe)))

  public_hyper <- grep(
    paste0("^(sd_", group, "__|cor_", group, "__|df_", group,
           "$|sigma$)"),
    variables(fit), value = TRUE
  )
  if (length(public_hyper)) {
    hyper <- exact_draw_matrix(fit, public_hyper)
    for (name in colnames(hyper)) {
      out[[paste0("parameter:", name)]] <- hyper[, name]
    }
  }

  selected_levels <- unique(c(1L, dim(re)[2L]))
  selected_coef <- unique(c(1L, dim(re)[3L]))
  for (j in selected_levels) {
    for (k in selected_coef) {
      label <- paste(dimnames(re)[[2L]][j], dimnames(re)[[3L]][k], sep = ",")
      out[[paste0("r:[", label, "]")]] <- re[, j, k]
      out[[paste0("b+r:[", label, "]")]] <- fe[, k] + re[, j, k]
    }
  }

  eta <- posterior_linpred(fit, transform = FALSE)
  selected_obs <- unique(round(seq(1, ncol(eta), length.out = 5L)))
  for (i in selected_obs) {
    out[[sprintf("eta:observation[%d]", i)]] <- eta[, i]
  }
  out[["eta:observation-average"]] <- rowMeans(eta)
  out
}

safe_mcse <- function(x, metric) {
  ans <- switch(
    metric,
    mean = posterior::mcse_mean(x),
    sd = posterior::mcse_sd(x),
    stop("Unknown comparison metric.", call. = FALSE)
  )
  if (!is.finite(ans)) {
    stop("Could not estimate MCSE for a sampled public quantity.",
         call. = FALSE)
  }
  ans
}

compare_public_posteriors <- function(conventional, s2z, case) {
  if (!setequal(names(conventional), names(s2z))) {
    stop("Public quantities differ between conventional and S2Z fits.\n",
         "Only conventional: ",
         paste(setdiff(names(conventional), names(s2z)), collapse = ", "),
         "\nOnly S2Z: ",
         paste(setdiff(names(s2z), names(conventional)), collapse = ", "),
         call. = FALSE)
  }
  quantity <- sort(names(conventional))
  rows <- vector("list", 2L * length(quantity))
  row <- 0L
  for (name in quantity) {
    x <- conventional[[name]]
    y <- s2z[[name]]
    for (metric in c("mean", "sd")) {
      row <- row + 1L
      estimate_x <- switch(metric, mean = mean(x), sd = stats::sd(x))
      estimate_y <- switch(metric, mean = mean(y), sd = stats::sd(y))
      mcse_x <- safe_mcse(x, metric)
      mcse_y <- safe_mcse(y, metric)
      mcse_difference <- sqrt(mcse_x^2 + mcse_y^2)
      difference <- estimate_y - estimate_x
      z_score <- if (mcse_difference > 0) {
        abs(difference) / mcse_difference
      } else if (identical(estimate_x, estimate_y)) {
        0
      } else {
        Inf
      }
      rows[[row]] <- data.frame(
        case = case,
        quantity = name,
        metric = metric,
        conventional = estimate_x,
        s2z = estimate_y,
        difference = difference,
        mcse_difference = mcse_difference,
        z_score = z_score,
        pass = z_score <= equiv_z,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  out[order(out$z_score, decreasing = TRUE), ]
}

s2z_coordinate_invariants <- function(fit, formula, data, prior,
                                      fixed_formula, case) {
  sdata <- standata(
    formula, data = data, family = gaussian(), prior = prior
  )
  X <- stats::model.matrix(fixed_formula, data)
  coef_names <- sub("^\\(Intercept\\)$", "Intercept", colnames(X))
  n_coef <- ncol(X)
  n_group <- sdata$N_1
  group_index <- as.integer(sdata$J_1)

  b_names <- sprintf("theta_s2z[%d]", seq_len(n_coef))
  b_s2z <- exact_draw_matrix(fit, b_names)
  n_draw <- nrow(b_s2z)
  r_s2z <- array(NA_real_, dim = c(n_draw, n_group, n_coef))
  for (k in seq_len(n_coef)) {
    r_names <- if (n_coef == 1L) {
      sprintf("r_s2z_1_1[%d]", seq_len(n_group))
    } else if (sprintf("r_s2z_1[1,%d]", k) %in% variables(fit)) {
      sprintf("r_s2z_1[%d,%d]", seq_len(n_group), k)
    } else {
      sprintf("r_s2z_1_%d[%d]", k, seq_len(n_group))
    }
    r_s2z[, , k] <- exact_draw_matrix(fit, r_names)
  }

  # Convert Stan's centered intercept coordinate to the raw model-matrix
  # intercept. Slopes and interaction coefficients are already raw.
  b_s2z_raw <- b_s2z
  if (n_coef > 1L) {
    means_X <- colMeans(X[, -1L, drop = FALSE])
    b_s2z_raw[, 1L] <- b_s2z[, 1L] -
      drop(b_s2z[, -1L, drop = FALSE] %*% means_X)
  }

  all_fe <- as.matrix(fixef(fit, summary = FALSE))
  fe <- public_fixef_draws(fit)
  re <- ranef(fit, summary = FALSE)$g
  if (!all(coef_names %in% colnames(fe)) ||
      !all(coef_names %in% dimnames(re)[[3L]])) {
    stop("Could not align public fixed/group coefficient names.",
         call. = FALSE)
  }
  fe <- fe[, coef_names, drop = FALSE]
  expected_levels <- levels(data$g)
  if (!all(expected_levels %in% dimnames(re)[[2L]])) {
    stop("Could not align public group levels with Stan group indices.",
         call. = FALSE)
  }
  re <- re[, expected_levels, coef_names, drop = FALSE]
  if (nrow(fe) != n_draw) {
    stop("Internal and public draw counts differ.", call. = FALSE)
  }

  eta_internal <- matrix(0, nrow = n_draw, ncol = nrow(X))
  eta_public <- matrix(0, nrow = n_draw, ncol = nrow(X))
  max_finite_error <- 0
  max_sum_to_zero_error <- 0
  for (k in seq_len(n_coef)) {
    internal_re <- matrix(r_s2z[, , k], nrow = n_draw, ncol = n_group)
    public_re <- matrix(re[, , k], nrow = n_draw, ncol = n_group)
    internal_finite <- sweep(internal_re, 1L, b_s2z_raw[, k], "+")
    public_finite <- sweep(public_re, 1L, fe[, k], "+")
    max_finite_error <- max(
      max_finite_error, max(abs(internal_finite - public_finite))
    )
    max_sum_to_zero_error <- max(
      max_sum_to_zero_error, max(abs(rowSums(internal_re)))
    )
    eta_internal <- eta_internal +
      sweep(internal_finite[, group_index, drop = FALSE],
            2L, X[, k], "*")
    eta_public <- eta_public +
      sweep(public_finite[, group_index, drop = FALSE],
            2L, X[, k], "*")
  }
  eta_brms <- posterior_linpred(fit, transform = FALSE)

  expected_public <- c(
    paste0("b_", coef_names),
    as.vector(outer(
      expected_levels, coef_names,
      function(level, coef) sprintf("r_g[%s,%s]", level, coef)
    ))
  )
  public_names_ok <- all(expected_public %in% variables(fit))
  fixef_names_ok <- setequal(colnames(all_fe), coef_names)
  out <- data.frame(
    case = case,
    check = c(
      "saved public b/r names use conventional API",
      "fixef exposes only public coefficient names",
      "internal group effects sum to zero",
      "public b+r equals internal finite coefficient",
      "posterior_linpred equals public b+r predictor",
      "posterior_linpred equals internal S2Z predictor"
    ),
    max_abs_error = c(
      if (public_names_ok) 0 else Inf,
      if (fixef_names_ok) 0 else Inf,
      max_sum_to_zero_error,
      max_finite_error,
      max(abs(eta_brms - eta_public)),
      max(abs(eta_brms - eta_internal))
    ),
    tolerance = invariant_tolerance,
    stringsAsFactors = FALSE
  )
  out$pass <- out$max_abs_error <= out$tolerance
  out
}

fit_quality <- function(fit, case, parameterization) {
  public <- grep(
    "^(b_|sd_g__|cor_g__|df_g$|sigma$|r_g\\[)",
    variables(fit), value = TRUE
  )
  public <- public[!grepl("^theta_s2z(\\[|$)", public)]
  draws <- as_draws_array(fit, variable = public)
  summary <- posterior::summarise_draws(
    draws,
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  )
  np <- nuts_params(fit)
  divergences <- sum(np$Value[np$Parameter == "divergent__"])
  data.frame(
    case = case,
    parameterization = parameterization,
    max_rhat = max(summary$rhat, na.rm = TRUE),
    min_bulk_ess = min(summary$ess_bulk, na.rm = TRUE),
    min_tail_ess = min(summary$ess_tail, na.rm = TRUE),
    divergences = divergences,
    pass = max(summary$rhat, na.rm = TRUE) <= max_rhat &&
      divergences <= max_divergences,
    stringsAsFactors = FALSE
  )
}

run_case <- function(case, data, conventional_formula, s2z_formula, prior,
                     fixed_formula, seed) {
  conventional <- sample_model(
    conventional_formula, data, prior, seed,
    paste0(case, "-conventional")
  )
  s2z <- sample_model(
    s2z_formula, data, prior, seed + 1L,
    paste0(case, "-s2z")
  )
  comparison <- compare_public_posteriors(
    public_samples(conventional, data),
    public_samples(s2z, data),
    case = case
  )
  invariants <- s2z_coordinate_invariants(
    s2z, s2z_formula, data, prior, fixed_formula, case
  )
  quality <- rbind(
    fit_quality(conventional, case, "conventional"),
    fit_quality(s2z, case, "s2z")
  )
  list(
    conventional = conventional,
    s2z = s2z,
    comparison = comparison,
    invariants = invariants,
    quality = quality
  )
}

gaussian_data <- simulate_gaussian_interaction()
gaussian_prior <- model_prior()
gaussian_result <- run_case(
  case = "gaussian-interaction",
  data = gaussian_data,
  conventional_formula = bf(y ~ x * z + (1 + x * z | g)),
  s2z_formula = bf(y ~ x * z + (1 + x * z | gr(g, s2z = TRUE))),
  prior = gaussian_prior,
  fixed_formula = ~ x * z,
  seed = 5101L
)

student_data <- simulate_student_group()
student_prior <- model_prior(student_group = TRUE)
student_result <- run_case(
  case = "student-group-effects",
  data = student_data,
  conventional_formula = bf(
    y ~ x + (1 + x | gr(g, dist = "student"))
  ),
  s2z_formula = bf(
    y ~ x + (1 + x | gr(g, dist = "student", s2z = TRUE))
  ),
  prior = student_prior,
  fixed_formula = ~ x,
  seed = 5102L
)

diagonal_gaussian_result <- run_case(
  case = "independent-gaussian-interaction",
  data = gaussian_data,
  conventional_formula = bf(y ~ x * z + (1 + x * z || g)),
  s2z_formula = bf(
    y ~ x * z + (1 + x * z || gr(g, s2z = TRUE))
  ),
  prior = model_prior(correlated = FALSE),
  fixed_formula = ~ x * z,
  seed = 5105L
)

diagonal_student_result <- run_case(
  case = "independent-student-intercept-slope",
  data = student_data,
  conventional_formula = bf(
    y ~ x + (1 + x || gr(g, dist = "student"))
  ),
  s2z_formula = bf(
    y ~ x +
      (1 + x || gr(g, dist = "student", s2z = TRUE))
  ),
  prior = model_prior(student_group = TRUE, correlated = FALSE),
  fixed_formula = ~ x,
  seed = 5106L
)

scalar_intercept_data <- simulate_scalar_intercept()
scalar_intercept_result <- run_case(
  case = "scalar-gaussian-intercept",
  data = scalar_intercept_data,
  conventional_formula = bf(y ~ 1 + (1 | g)),
  s2z_formula = bf(y ~ 1 + (1 | gr(g, s2z = TRUE))),
  prior = scalar_prior(intercept = TRUE),
  fixed_formula = ~ 1,
  seed = 5103L
)

scalar_student_data <- simulate_scalar_student_slope()
scalar_student_result <- run_case(
  case = "scalar-student-slope-no-intercept",
  data = scalar_student_data,
  conventional_formula = bf(
    y ~ 0 + x + (0 + x | gr(g, dist = "student"))
  ),
  s2z_formula = bf(
    y ~ 0 + x +
      (0 + x | gr(g, dist = "student", s2z = TRUE))
  ),
  prior = scalar_prior(intercept = FALSE, student_group = TRUE),
  fixed_formula = ~ 0 + x,
  seed = 5104L
)

comparison_results <- rbind(
  gaussian_result$comparison,
  student_result$comparison,
  diagonal_gaussian_result$comparison,
  diagonal_student_result$comparison,
  scalar_intercept_result$comparison,
  scalar_student_result$comparison
)
invariant_results <- rbind(
  gaussian_result$invariants,
  student_result$invariants,
  diagonal_gaussian_result$invariants,
  diagonal_student_result$invariants,
  scalar_intercept_result$invariants,
  scalar_student_result$invariants
)
quality_results <- rbind(
  gaussian_result$quality,
  student_result$quality,
  diagonal_gaussian_result$quality,
  diagonal_student_result$quality,
  scalar_intercept_result$quality,
  scalar_student_result$quality
)

cat("\nMCSE-aware conventional-vs-S2Z comparisons (largest z-scores first)\n")
print(comparison_results, row.names = FALSE, digits = 4)
cat("\nExact public/internal coordinate invariants\n")
print(invariant_results, row.names = FALSE, digits = 4)
cat("\nSampling quality\n")
print(quality_results, row.names = FALSE, digits = 4)

failed_comparisons <- subset(comparison_results, !pass)
failed_invariants <- subset(invariant_results, !pass)
failed_quality <- subset(quality_results, !pass)
if (nrow(failed_comparisons) || nrow(failed_invariants) ||
    nrow(failed_quality)) {
  stop(
    "S2Z sampling validation failed: ",
    nrow(failed_comparisons), " MCSE comparisons, ",
    nrow(failed_invariants), " coordinate invariants, and ",
    nrow(failed_quality), " sampling-quality checks failed.",
    call. = FALSE
  )
}

cat("\nPASS: conventional and S2Z posterior draws agree within propagated ",
    "MCSE, public b/r/predictor identities hold draw by draw, and all fits ",
    "meet the configured sampler checks.\n", sep = "")
