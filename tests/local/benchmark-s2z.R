# Reproducible benchmark for physical sum-to-zero group-level effects.
#
# This script compiles the conventional and S2Z programs once, then times only
# refits that reuse those compiled programs. It does not run as part of the
# package test suite. A modest default run is
#
#   Rscript tests/local/benchmark-s2z.R
#
# Scale it up with environment variables, for example
#
#   BRMS_S2Z_BENCH_CHAINS=4 \
#   BRMS_S2Z_BENCH_WARMUP=1000 \
#   BRMS_S2Z_BENCH_SAMPLING=1000 \
#   BRMS_S2Z_BENCH_GROUPS=80 \
#   BRMS_S2Z_BENCH_MODEL=orthogonal \
#   BRMS_S2Z_BENCH_K=10 \
#   BRMS_S2Z_BENCH_COVARIANCE=independent \
#   BRMS_S2Z_BENCH_REPS=3 \
#   BRMS_S2Z_BENCH_OUT=s2z-benchmark.csv \
#   Rscript tests/local/benchmark-s2z.R
#
# Defaults are deliberately suitable for a local smoke benchmark:
#
#   BRMS_S2Z_BENCH_CHAINS=2
#   BRMS_S2Z_BENCH_WARMUP=400
#   BRMS_S2Z_BENCH_SAMPLING=400
#   BRMS_S2Z_BENCH_GROUPS=24
#   BRMS_S2Z_BENCH_MODEL=interaction  # also scalar-*, orthogonal
#   BRMS_S2Z_BENCH_K=10               # orthogonal model only
#   BRMS_S2Z_BENCH_COVARIANCE=independent  # or correlated; orthogonal only
#   BRMS_S2Z_BENCH_GROUP_DIST=gaussian     # or student
#   BRMS_S2Z_BENCH_IDENTIFICATION=both     # or strong/weak
#   BRMS_S2Z_BENCH_SIGMA=0.08         # orthogonal residual SD
#   BRMS_S2Z_BENCH_STRONG_N=24       # observations per group
#   BRMS_S2Z_BENCH_WEAK_N=6
#   BRMS_S2Z_BENCH_REPS=1
#   BRMS_S2Z_BENCH_BACKEND=auto
#   BRMS_S2Z_BENCH_ADAPT_DELTA=0.95
#   BRMS_S2Z_BENCH_MAX_TREEDEPTH=12
#
# The orthogonal/independent setting exercises the component-wise (scalar)
# S2Z algebra with any K; "scalar" does not imply K = 1.
#
# `gradient_evals` counts post-warmup leapfrog gradients plus one initial
# gradient per transition. ESS/gradient is calculated on the same post-warmup
# transitions. Wall time is end-to-end refit time: it includes warmup,
# sampling, process/CSV overhead, and brms fit construction, but excludes the
# separately reported compile-and-smoke time.

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

env_choice <- function(name, default, choices) {
  value <- Sys.getenv(name, unset = default)
  if (!value %in% choices) {
    stop(name, " must be one of: ", paste(choices, collapse = ", "), ".",
         call. = FALSE)
  }
  value
}

choose_backend <- function() {
  requested <- tolower(Sys.getenv(
    "BRMS_S2Z_BENCH_BACKEND",
    unset = Sys.getenv("BRMS_S2Z_BACKEND", unset = "auto")
  ))
  if (requested %in% c("cmdstanr", "rstan")) {
    return(requested)
  }
  if (!identical(requested, "auto")) {
    stop("BRMS_S2Z_BENCH_BACKEND must be auto, cmdstanr, or rstan.",
         call. = FALSE)
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

chains <- env_integer("BRMS_S2Z_BENCH_CHAINS", 2L, minimum = 2L)
warmup <- env_integer("BRMS_S2Z_BENCH_WARMUP", 400L)
sampling <- env_integer("BRMS_S2Z_BENCH_SAMPLING", 400L)
groups <- env_integer("BRMS_S2Z_BENCH_GROUPS", 24L, minimum = 3L)
benchmark_model <- env_choice(
  "BRMS_S2Z_BENCH_MODEL", "interaction",
  c("interaction", "scalar-intercept", "scalar-slope", "orthogonal")
)
orthogonal_k <- env_integer("BRMS_S2Z_BENCH_K", 10L)
benchmark_covariance <- if (benchmark_model == "orthogonal") {
  env_choice(
    "BRMS_S2Z_BENCH_COVARIANCE", "independent",
    c("independent", "correlated")
  )
} else if (benchmark_model == "interaction") {
  "correlated"
} else {
  "independent"
}
group_dist <- env_choice(
  "BRMS_S2Z_BENCH_GROUP_DIST", "gaussian", c("gaussian", "student")
)
identification <- env_choice(
  "BRMS_S2Z_BENCH_IDENTIFICATION", "both",
  c("strong", "weak", "both")
)
identification_levels <- if (identification == "both") {
  c("strong", "weak")
} else {
  identification
}
orthogonal_sigma <- env_number(
  "BRMS_S2Z_BENCH_SIGMA", 0.08, minimum = .Machine$double.eps
)
strong_n <- env_integer("BRMS_S2Z_BENCH_STRONG_N", 24L, minimum = 4L)
weak_n <- env_integer("BRMS_S2Z_BENCH_WEAK_N", 6L, minimum = 4L)
reps <- env_integer("BRMS_S2Z_BENCH_REPS", 1L)
base_seed <- env_integer("BRMS_S2Z_BENCH_SEED", 6201L, minimum = 1L)
adapt_delta <- env_number(
  "BRMS_S2Z_BENCH_ADAPT_DELTA", 0.95, minimum = 0.5
)
max_treedepth <- env_integer("BRMS_S2Z_BENCH_MAX_TREEDEPTH", 12L)
detected_cores <- parallel::detectCores(logical = FALSE)
if (is.na(detected_cores)) {
  detected_cores <- 1L
}
cores <- env_integer(
  "BRMS_S2Z_BENCH_CORES", min(chains, detected_cores), minimum = 1L
)
backend <- choose_backend()
output_file <- Sys.getenv("BRMS_S2Z_BENCH_OUT", unset = "")
options(mc.cores = cores, width = 180)

benchmark_k <- switch(
  benchmark_model,
  interaction = 4L,
  `scalar-intercept` = 1L,
  `scalar-slope` = 1L,
  orthogonal = orthogonal_k
)
if (benchmark_model == "orthogonal" &&
    "strong" %in% identification_levels && strong_n < benchmark_k) {
  stop(
    "BRMS_S2Z_BENCH_STRONG_N must be at least BRMS_S2Z_BENCH_K ",
    "for the exact orthogonal design.", call. = FALSE
  )
}

backend_version <- if (backend == "cmdstanr") {
  paste(cmdstanr::cmdstan_version(), collapse = ".")
} else {
  as.character(utils::packageVersion("rstan"))
}
configuration <- data.frame(
  brms_version = as.character(utils::packageVersion("brms")),
  backend = backend,
  backend_version = backend_version,
  chains = chains,
  warmup = warmup,
  sampling = sampling,
  model = benchmark_model,
  k = benchmark_k,
  covariance = benchmark_covariance,
  group_dist = group_dist,
  identification = identification,
  groups = groups,
  strong_n = strong_n,
  weak_n = weak_n,
  orthogonal_sigma = orthogonal_sigma,
  reps = reps,
  seed = base_seed,
  cores = cores,
  adapt_delta = adapt_delta,
  max_treedepth = max_treedepth,
  stringsAsFactors = FALSE
)
cat("S2Z benchmark configuration\n")
print(configuration, row.names = FALSE)

rmvn_rows <- function(n, sigma) {
  matrix(stats::rnorm(n * nrow(sigma)), nrow = n) %*% chol(sigma)
}

draw_group_effects <- function(n, sigma, dist, df = 5) {
  out <- rmvn_rows(n, sigma)
  if (dist == "student") {
    # One scale per group gives an elliptically symmetric multivariate t.
    out <- out * sqrt(df / stats::rchisq(n, df = df))
  }
  out
}

dense_benchmark_design <- function(n, k, identification) {
  if (identification == "strong") {
    if (n < k) {
      stop("The strong orthogonal design requires n >= k.", call. = FALSE)
    }
    q <- qr.Q(qr(matrix(stats::rnorm(n * k), nrow = n)))
    out <- sqrt(n) * q[, seq_len(k), drop = FALSE]
    target <- diag(n, nrow = k, ncol = k)
    if (max(abs(crossprod(out) - target)) > 1e-8 * n) {
      stop("Failed to construct an orthogonal benchmark design.",
           call. = FALSE)
    }
  } else {
    weak_cor <- outer(seq_len(k), seq_len(k), function(i, j) {
      0.97^abs(i - j)
    })
    out <- matrix(stats::rnorm(n * k), nrow = n) %*% chol(weak_cor)
    out <- sweep(out, 2, sqrt(colMeans(out^2)), "/")
  }
  colnames(out) <- paste0("x", seq_len(k))
  out
}

simulate_benchmark_data <- function(identification, seed, groups,
                                    observations_per_group, model, k,
                                    covariance, group_dist,
                                    orthogonal_sigma) {
  identification <- match.arg(identification, c("strong", "weak"))
  set.seed(seed)
  g_levels <- sprintf("g%03d", seq_len(groups))
  pieces <- lapply(seq_len(groups), function(j) {
    n <- observations_per_group
    if (model == "orthogonal") {
      design <- dense_benchmark_design(n, k, identification)
      return(data.frame(g = g_levels[j], design, check.names = FALSE))
    } else if (identification == "strong") {
      grid <- expand.grid(x = c(-1, 1), z = c(-1, 1))
      grid <- grid[rep(seq_len(nrow(grid)), length.out = n), ]
      x <- grid$x + stats::rnorm(n, sd = 0.08)
      z <- grid$z + stats::rnorm(n, sd = 0.08)
    } else {
      latent <- stats::rnorm(n)
      x <- latent + stats::rnorm(n, sd = 0.08)
      z <- 0.95 * latent + stats::rnorm(n, sd = 0.16)
    }
    data.frame(g = g_levels[j], x = x, z = z)
  })
  dat <- do.call(rbind, pieces)
  rownames(dat) <- NULL
  dat$g <- factor(dat$g, levels = g_levels)

  if (model == "scalar-intercept") {
    beta <- 0.35
    u <- draw_group_effects(groups, matrix(0.75^2), group_dist)
    eta <- beta + u[as.integer(dat$g), 1]
  } else if (model == "scalar-slope") {
    beta <- 0.70
    u <- draw_group_effects(groups, matrix(0.50^2), group_dist)
    eta <- dat$x * (beta + u[as.integer(dat$g), 1])
  } else if (model == "orthogonal") {
    x_names <- paste0("x", seq_len(k))
    X <- as.matrix(dat[x_names])
    beta <- seq(0.20, 0.80, length.out = k) * (-1)^(seq_len(k) - 1L)
    re_sd <- seq(0.35, 0.80, length.out = k)
    re_cor <- if (covariance == "independent") {
      diag(k)
    } else {
      outer(seq_len(k), seq_len(k), function(i, j) 0.25^abs(i - j))
    }
    sigma_re <- diag(re_sd) %*% re_cor %*% diag(re_sd)
    u <- draw_group_effects(groups, sigma_re, group_dist)
    eta <- drop(X %*% beta) + rowSums(X * u[as.integer(dat$g), ])
  } else {
    X <- stats::model.matrix(~ x * z, dat)
    beta <- c(0.35, 0.70, -0.45, 0.30)
    re_sd <- c(0.75, 0.50, 0.45, 0.35)
    re_cor <- outer(seq_along(re_sd), seq_along(re_sd), function(i, j) {
      0.25^abs(i - j)
    })
    sigma_re <- diag(re_sd) %*% re_cor %*% diag(re_sd)
    u <- draw_group_effects(groups, sigma_re, group_dist)
    eta <- drop(X %*% beta) + rowSums(X * u[as.integer(dat$g), ])
  }
  residual_sd <- if (model == "orthogonal") {
    orthogonal_sigma
  } else if (identification == "strong") {
    0.35
  } else {
    1.40
  }
  dat$y <- stats::rnorm(nrow(dat), eta, residual_sd)
  dat
}

benchmark_data <- lapply(identification_levels, function(level) {
  simulate_benchmark_data(
    level,
    seed = base_seed + if (level == "strong") 10L else 20L,
    groups = groups,
    observations_per_group = if (level == "strong") strong_n else weak_n,
    model = benchmark_model,
    k = benchmark_k,
    covariance = benchmark_covariance,
    group_dist = group_dist,
    orthogonal_sigma = orthogonal_sigma
  )
})
names(benchmark_data) <- identification_levels

benchmark_prior <- if (benchmark_model == "scalar-intercept") {
  prior(normal(0, 2), class = "Intercept")
} else if (benchmark_model %in% c("scalar-slope", "orthogonal")) {
  prior(normal(0, 2), class = "b")
} else {
  prior(normal(0, 2), class = "b") +
    prior(normal(0, 2), class = "Intercept")
}
benchmark_prior <- benchmark_prior +
  prior(exponential(1), class = "sd") +
  prior(exponential(1), class = "sigma")
if (benchmark_covariance == "correlated" && benchmark_k > 1L) {
  benchmark_prior <- benchmark_prior + prior(lkj(2), class = "cor")
}
if (group_dist == "student") {
  benchmark_prior <- benchmark_prior +
    prior(gamma(2, 0.2), class = "df", group = "g")
}

group_expression <- function(s2z, explicit_cor = NULL) {
  args <- "g"
  if (!is.null(explicit_cor)) {
    args <- c(args, paste0("cor = ", as.character(explicit_cor)))
  }
  if (group_dist == "student") {
    args <- c(args, 'dist = "student"')
  }
  if (s2z) {
    args <- c(args, "s2z = TRUE")
  }
  if (length(args) == 1L) {
    return(args)
  }
  paste0("gr(", paste(args, collapse = ", "), ")")
}

orthogonal_cor <- identical(benchmark_covariance, "correlated")
conventional_group <- group_expression(
  FALSE,
  explicit_cor = if (benchmark_model == "orthogonal") orthogonal_cor else NULL
)
s2z_group <- group_expression(
  TRUE,
  explicit_cor = if (benchmark_model == "orthogonal") orthogonal_cor else NULL
)

benchmark_formula <- if (benchmark_model == "scalar-intercept") {
  list(
    conventional = bf(stats::as.formula(paste0(
      "y ~ 1 + (1 | ", conventional_group, ")"
    ))),
    s2z = bf(stats::as.formula(paste0(
      "y ~ 1 + (1 | ", s2z_group, ")"
    )))
  )
} else if (benchmark_model == "scalar-slope") {
  list(
    conventional = bf(stats::as.formula(paste0(
      "y ~ 0 + x + (0 + x | ", conventional_group, ")"
    ))),
    s2z = bf(stats::as.formula(paste0(
      "y ~ 0 + x + (0 + x | ", s2z_group, ")"
    )))
  )
} else if (benchmark_model == "orthogonal") {
  x_terms <- paste0("x", seq_len(benchmark_k), collapse = " + ")
  list(
    conventional = bf(stats::as.formula(paste0(
      "y ~ 0 + ", x_terms, " + (0 + ", x_terms, " | ",
      conventional_group, ")"
    ))),
    s2z = bf(stats::as.formula(paste0(
      "y ~ 0 + ", x_terms, " + (0 + ", x_terms, " | ",
      s2z_group, ")"
    )))
  )
} else {
  list(
    conventional = bf(stats::as.formula(paste0(
      "y ~ x * z + (1 + x * z | ", conventional_group, ")"
    ))),
    s2z = bf(stats::as.formula(paste0(
      "y ~ x * z + (1 + x * z | ", s2z_group, ")"
    )))
  )
}
coefficients_per_group <- benchmark_k

compile_template <- function(formula, data, parameterization) {
  cat("Compiling and smoke-running ", parameterization, " model ...\n",
      sep = "")
  timing <- system.time({
    fit <- brm(
      formula = formula,
      data = data,
      family = gaussian(),
      prior = benchmark_prior,
      backend = backend,
      chains = 1,
      cores = 1,
      iter = 2,
      warmup = 1,
      seed = base_seed + if (parameterization == "s2z") 2L else 1L,
      refresh = 0,
      silent = 1,
      save_pars = save_pars(all = TRUE),
      control = list(adapt_delta = 0.8, max_treedepth = 5)
    )
  })
  list(fit = fit, seconds = unname(timing[["elapsed"]]))
}

templates <- lapply(names(benchmark_formula), function(parameterization) {
  compile_template(
    benchmark_formula[[parameterization]], benchmark_data[[1]],
    parameterization
  )
})
names(templates) <- names(benchmark_formula)

public_benchmark_variables <- function(fit) {
  out <- grep(
    "^(b_|sd_g__|cor_g__|df_g($|__)|sigma$|r_g\\[)",
    variables(fit), value = TRUE
  )
  out[!grepl("^theta_s2z(\\[|$)", out)]
}

ebfmi_by_chain <- function(nuts) {
  energy <- nuts[nuts$Parameter == "energy__", , drop = FALSE]
  split_energy <- split(energy, energy$Chain)
  vapply(split_energy, function(x) {
    x <- x[order(x$Iteration), , drop = FALSE]
    value <- x$Value
    variance <- stats::var(value)
    if (length(value) < 3L || !is.finite(variance) || variance <= 0) {
      return(NA_real_)
    }
    mean(diff(value)^2) / variance
  }, numeric(1))
}

sampler_report <- function(fit) {
  monitored <- public_benchmark_variables(fit)
  summaries <- posterior::summarise_draws(
    as_draws_array(fit, variable = monitored),
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  )
  nuts <- nuts_params(fit)
  leapfrog <- nuts$Value[nuts$Parameter == "n_leapfrog__"]
  if (!length(leapfrog)) {
    stop("n_leapfrog__ was not available from nuts_params().", call. = FALSE)
  }
  gradient_evals <- sum(leapfrog + 1)
  divergences <- sum(nuts$Value[nuts$Parameter == "divergent__"])
  treedepth <- nuts$Value[nuts$Parameter == "treedepth__"]
  ebfmi <- ebfmi_by_chain(nuts)
  finite_rhat <- summaries$rhat[is.finite(summaries$rhat)]
  finite_bulk <- summaries$ess_bulk[is.finite(summaries$ess_bulk)]
  finite_tail <- summaries$ess_tail[is.finite(summaries$ess_tail)]
  if (!length(finite_rhat) || !length(finite_bulk) || !length(finite_tail)) {
    stop("Posterior diagnostics were not finite.", call. = FALSE)
  }
  data.frame(
    monitored_parameters = nrow(summaries),
    postwarmup_transitions = length(leapfrog),
    gradient_evals = gradient_evals,
    mean_leapfrog_steps = mean(leapfrog),
    median_bulk_ess = stats::median(finite_bulk),
    min_bulk_ess = min(finite_bulk),
    median_tail_ess = stats::median(finite_tail),
    min_tail_ess = min(finite_tail),
    median_bulk_ess_per_gradient =
      stats::median(finite_bulk) / gradient_evals,
    median_tail_ess_per_gradient =
      stats::median(finite_tail) / gradient_evals,
    median_bulk_ess_per_1000_gradients =
      1000 * stats::median(finite_bulk) / gradient_evals,
    median_tail_ess_per_1000_gradients =
      1000 * stats::median(finite_tail) / gradient_evals,
    divergences = divergences,
    treedepth_hits = sum(treedepth >= max_treedepth),
    min_ebfmi = if (all(is.na(ebfmi))) NA_real_ else min(ebfmi, na.rm = TRUE),
    max_rhat = max(finite_rhat),
    stringsAsFactors = FALSE
  )
}

sample_job <- function(template, data, seed) {
  timing <- system.time({
    fit <- update(
      template,
      newdata = data,
      recompile = FALSE,
      algorithm = "sampling",
      chains = chains,
      cores = cores,
      iter = warmup + sampling,
      warmup = warmup,
      seed = seed,
      refresh = 0,
      silent = 1,
      control = list(
        adapt_delta = adapt_delta,
        max_treedepth = max_treedepth
      )
    )
  })
  report <- sampler_report(fit)
  report$wall_seconds <- unname(timing[["elapsed"]])
  report$median_bulk_ess_per_second <-
    report$median_bulk_ess / report$wall_seconds
  report$median_tail_ess_per_second <-
    report$median_tail_ess / report$wall_seconds
  report
}

jobs <- expand.grid(
  identification = names(benchmark_data),
  parameterization = names(benchmark_formula),
  replicate = seq_len(reps),
  stringsAsFactors = FALSE
)
# Deterministic randomization reduces systematic thermal/order effects while
# keeping every run exactly reproducible from BRMS_S2Z_BENCH_SEED.
set.seed(base_seed + 100L)
jobs <- jobs[sample(seq_len(nrow(jobs))), , drop = FALSE]

rows <- vector("list", nrow(jobs))
for (i in seq_len(nrow(jobs))) {
  job <- jobs[i, ]
  identification_id <- match(job$identification, c("strong", "weak"))
  parameterization_id <- match(
    job$parameterization, c("conventional", "s2z")
  )
  job_seed <- base_seed + 1000L * job$replicate +
    100L * identification_id + 10L * parameterization_id
  cat(
    sprintf(
      "Benchmarking %s identification, %s parameterization, replicate %d (seed %d) ...\n",
      job$identification, job$parameterization, job$replicate, job_seed
    )
  )
  report <- sample_job(
    templates[[job$parameterization]]$fit,
    benchmark_data[[job$identification]],
    seed = job_seed
  )
  rows[[i]] <- cbind(
    data.frame(
      identification = job$identification,
      model = benchmark_model,
      k = benchmark_k,
      covariance = benchmark_covariance,
      group_dist = group_dist,
      parameterization = job$parameterization,
      replicate = job$replicate,
      seed = job_seed,
      observations = nrow(benchmark_data[[job$identification]]),
      groups = groups,
      coefficients_per_group = coefficients_per_group,
      compile_and_smoke_seconds =
        templates[[job$parameterization]]$seconds,
      stringsAsFactors = FALSE
    ),
    report
  )
  print(rows[[i]], row.names = FALSE, digits = 4)
  invisible(gc())
}

benchmark_results <- do.call(rbind, rows)
benchmark_results <- benchmark_results[
  order(
    match(benchmark_results$identification, c("strong", "weak")),
    benchmark_results$replicate,
    match(benchmark_results$parameterization, c("conventional", "s2z"))
  ),
]
rownames(benchmark_results) <- NULL

conventional <- subset(
  benchmark_results, parameterization == "conventional",
  select = -parameterization
)
s2z <- subset(
  benchmark_results, parameterization == "s2z",
  select = -parameterization
)
paired <- merge(
  conventional, s2z,
  by = c(
    "identification", "model", "k", "covariance", "group_dist",
    "replicate", "observations", "groups", "coefficients_per_group"
  ),
  suffixes = c("_conventional", "_s2z"), sort = FALSE
)
paired_results <- data.frame(
  identification = paired$identification,
  model = paired$model,
  k = paired$k,
  covariance = paired$covariance,
  group_dist = paired$group_dist,
  replicate = paired$replicate,
  wall_time_speedup_conventional_over_s2z =
    paired$wall_seconds_conventional / paired$wall_seconds_s2z,
  gradient_reduction_conventional_over_s2z =
    paired$gradient_evals_conventional / paired$gradient_evals_s2z,
  bulk_ess_per_gradient_gain_s2z_over_conventional =
    paired$median_bulk_ess_per_gradient_s2z /
      paired$median_bulk_ess_per_gradient_conventional,
  tail_ess_per_gradient_gain_s2z_over_conventional =
    paired$median_tail_ess_per_gradient_s2z /
      paired$median_tail_ess_per_gradient_conventional,
  bulk_ess_per_second_gain_s2z_over_conventional =
    paired$median_bulk_ess_per_second_s2z /
      paired$median_bulk_ess_per_second_conventional,
  stringsAsFactors = FALSE
)
paired_results <- paired_results[
  order(match(paired_results$identification, c("strong", "weak")),
        paired_results$replicate),
]
rownames(paired_results) <- NULL

cat("\nFull benchmark results\n")
print(benchmark_results, row.names = FALSE, digits = 5)
cat("\nPaired conventional/S2Z efficiency ratios (> 1 favors S2Z)\n")
print(paired_results, row.names = FALSE, digits = 5)

if (nzchar(output_file)) {
  output_dir <- dirname(output_file)
  if (!identical(output_dir, ".")) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  companion_csv <- function(path, suffix) {
    if (grepl("\\.csv$", path, ignore.case = TRUE)) {
      sub("\\.csv$", paste0(suffix, ".csv"), path, ignore.case = TRUE)
    } else {
      paste0(path, suffix, ".csv")
    }
  }
  utils::write.csv(benchmark_results, output_file, row.names = FALSE)
  paired_file <- companion_csv(output_file, "-paired")
  settings_file <- companion_csv(output_file, "-settings")
  utils::write.csv(paired_results, paired_file, row.names = FALSE)
  utils::write.csv(configuration, settings_file, row.names = FALSE)
  cat(
    "\nWrote ", normalizePath(output_file, mustWork = FALSE),
    ", ", normalizePath(paired_file, mustWork = FALSE),
    ", and ", normalizePath(settings_file, mustWork = FALSE), ".\n",
    sep = ""
  )
}
