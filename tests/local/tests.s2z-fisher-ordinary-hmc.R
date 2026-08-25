# Preliminary HMC benchmark for ordinary (latent = FALSE) Fisher-centered S2Z
# group effects. The three parameterizations below have the same posterior;
# only their coordinates differ.
# ESS per gradient below uses post-warmup transitions. Wall time includes warmup,
# and compile timings may include a cached executable, so compile times should not
# be compared across modes.
#
# Default run:
#   Rscript tests/local/tests.s2z-fisher-ordinary-hmc.R
#
# Useful controls:
#   BRMS_S2Z_FISHER_CHAINS=4
#   BRMS_S2Z_FISHER_WARMUP=400
#   BRMS_S2Z_FISHER_SAMPLING=500
#   BRMS_S2Z_FISHER_CORES=4
#   BRMS_S2Z_FISHER_CASES=gaussian_intercept,gaussian_slope,student_slope
#   BRMS_S2Z_FISHER_MODES=noncentered,centered,auto
#   BRMS_S2Z_FISHER_OUT=/tmp/brms-s2z-fisher-ordinary

suppressPackageStartupMessages({
  library(brms)
  library(cmdstanr)
  library(posterior)
})

env_integer <- function(name, default, minimum = 1L) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(as.integer(default))
  }
  value <- suppressWarnings(as.integer(value))
  if (length(value) != 1L || is.na(value) || value < minimum) {
    stop(name, " must be an integer >= ", minimum, call. = FALSE)
  }
  value
}

env_number <- function(name, default, minimum = -Inf, maximum = Inf) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(as.numeric(default))
  }
  value <- suppressWarnings(as.numeric(value))
  if (length(value) != 1L || !is.finite(value) ||
      value < minimum || value > maximum) {
    stop(name, " must be in [", minimum, ", ", maximum, "]", call. = FALSE)
  }
  value
}

env_subset <- function(name, default, choices) {
  value <- strsplit(Sys.getenv(name, unset = default), ",", fixed = TRUE)[[1L]]
  value <- trimws(value)
  if (!length(value) || any(!value %in% choices)) {
    stop(name, " must be a comma-separated subset of ",
         paste(choices, collapse = ", "), call. = FALSE)
  }
  unique(value)
}

chains <- env_integer("BRMS_S2Z_FISHER_CHAINS", 4L)
warmup <- env_integer("BRMS_S2Z_FISHER_WARMUP", 400L)
sampling <- env_integer("BRMS_S2Z_FISHER_SAMPLING", 500L)
cores <- env_integer("BRMS_S2Z_FISHER_CORES", min(chains, 4L))
adapt_delta <- env_number(
  "BRMS_S2Z_FISHER_ADAPT_DELTA", 0.9, minimum = 0.5, maximum = 0.9999
)
max_treedepth <- env_integer("BRMS_S2Z_FISHER_MAX_TREEDEPTH", 12L)
all_cases <- c("gaussian_intercept", "gaussian_slope", "student_slope")
all_modes <- c("noncentered", "centered", "auto")
cases <- env_subset(
  "BRMS_S2Z_FISHER_CASES", paste(all_cases, collapse = ","), all_cases
)
modes <- env_subset(
  "BRMS_S2Z_FISHER_MODES", paste(all_modes, collapse = ","), all_modes
)
output_dir <- Sys.getenv(
  "BRMS_S2Z_FISHER_OUT", unset = "/tmp/brms-s2z-fisher-ordinary"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "csv"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "stan"), recursive = TRUE, showWarnings = FALSE)

if (!nzchar(tryCatch(cmdstanr::cmdstan_path(), error = function(e) ""))) {
  stop("A working CmdStan installation is required.", call. = FALSE)
}

center_code <- c(
  noncentered = "FALSE", centered = "TRUE", auto = "\"auto\""
)

correlated_effects <- function(groups, sd, correlation) {
  stopifnot(length(sd) == 2L, abs(correlation) < 1)
  R <- matrix(c(1, correlation, correlation, 1), 2L, 2L)
  L <- diag(sd) %*% t(chol(R))
  t(L %*% matrix(rnorm(2L * groups), 2L, groups))
}

simulate_gaussian_intercept <- function(seed = 81401L) {
  set.seed(seed)
  groups <- 36L
  observations <- rep(c(1L, 8L, 60L), each = groups / 3L)
  g <- factor(rep(seq_len(groups), observations))
  group_effect <- rnorm(groups, 0, 0.35)
  data <- data.frame(
    y = 0.4 + group_effect[as.integer(g)] + rnorm(length(g), 0, 1),
    g = g
  )
  list(
    data = data,
    family = gaussian(),
    prior = prior(normal(0, 1), class = Intercept) +
      prior(exponential(1), class = sd) +
      prior(exponential(1), class = sigma),
    counts = observations
  )
}

simulate_gaussian_slope <- function(seed = 81402L) {
  set.seed(seed)
  groups <- 36L
  regime <- rep(c("weak", "medium", "strong"), each = groups / 3L)
  observations <- c(weak = 2L, medium = 8L, strong = 35L)[regime]
  b <- correlated_effects(groups, c(0.7, 0.32), 0.65)
  rows <- lapply(seq_len(groups), function(i) {
    x <- switch(
      regime[i],
      weak = seq(-0.025, 0.025, length.out = observations[i]),
      medium = seq(-0.6, 0.6, length.out = observations[i]),
      strong = seq(-1.5, 1.5, length.out = observations[i])
    )
    data.frame(
      g = factor(i, levels = seq_len(groups)),
      x = x,
      y = 0.3 + b[i, 1L] + (-0.45 + b[i, 2L]) * x +
        rnorm(length(x), 0, 0.8)
    )
  })
  list(
    data = do.call(rbind, rows),
    family = gaussian(),
    prior = prior(normal(0, 1), class = Intercept) +
      prior(normal(0, 1), class = b) +
      prior(exponential(1), class = sd) +
      prior(lkj(2), class = cor) +
      prior(exponential(1), class = sigma)
  )
}

simulate_student_slope <- function(seed = 81403L) {
  set.seed(seed)
  groups <- 36L
  regime <- rep(c("weak", "medium", "strong"), each = groups / 3L)
  observations <- c(weak = 2L, medium = 9L, strong = 30L)[regime]
  b <- correlated_effects(groups, c(0.6, 0.3), -0.5)
  nu <- 4
  rows <- lapply(seq_len(groups), function(i) {
    x <- switch(
      regime[i],
      weak = seq(-0.05, 0.05, length.out = observations[i]),
      medium = seq(-0.75, 0.75, length.out = observations[i]),
      strong = seq(-1.4, 1.4, length.out = observations[i])
    )
    data.frame(
      g = factor(i, levels = seq_len(groups)),
      x = x,
      y = -0.2 + b[i, 1L] + (0.55 + b[i, 2L]) * x +
        0.75 * rt(length(x), df = nu)
    )
  })
  list(
    data = do.call(rbind, rows),
    family = student(),
    prior = prior(normal(0, 1), class = Intercept) +
      prior(normal(0, 1), class = b) +
      prior(exponential(1), class = sd) +
      prior(lkj(2), class = cor) +
      prior(exponential(1), class = sigma) +
      prior(gamma(2, 0.5), class = nu)
  )
}

simulate_case <- function(case) {
  switch(
    case,
    gaussian_intercept = simulate_gaussian_intercept(),
    gaussian_slope = simulate_gaussian_slope(),
    student_slope = simulate_student_slope()
  )
}

case_formula <- function(case, mode) {
  center <- center_code[[mode]]
  rhs <- if (identical(case, "gaussian_intercept")) {
    paste0(
      "y ~ 1 + (1 | gr(g, s2z = TRUE, center = ", center, "))"
    )
  } else {
    paste0(
      "y ~ 1 + x + (1 + x | gr(g, s2z = TRUE, center = ",
      center, "))"
    )
  }
  as.formula(rhs, env = globalenv())
}

ebfmi <- function(energy) {
  variance <- stats::var(energy)
  if (length(energy) < 3L || !is.finite(variance) || variance <= 0) {
    return(NA_real_)
  }
  mean(diff(energy)^2) / variance
}

model_variables <- function(fit) {
  posterior::variables(fit$draws())
}

physical_variables <- function(fit) {
  variables <- model_variables(fit)
  hyper <- grep(
    "^(b_Intercept$|b\\[|sd_1\\[|cor_1\\[|sigma$|nu$)",
    variables, value = TRUE
  )
  group <- grep("^r_1_[[:digit:]]+\\[", variables, value = TRUE)
  list(hyper = hyper, group = group, all = c(hyper, group))
}

rho_summary <- function(fit, case, mode) {
  variables <- grep(
    "^rho_s2z_1\\[[[:digit:]]+,[[:digit:]]+\\]$",
    model_variables(fit), value = TRUE
  )
  if (!length(variables)) {
    endpoint <- if (identical(mode, "centered")) 1 else 0
    return(data.frame(
      case = case, mode = mode, group = NA_integer_, coefficient = NA_integer_,
      mean_rho = endpoint
    ))
  }
  index <- utils::strcapture(
    "^rho_s2z_1\\[([[:digit:]]+),([[:digit:]]+)\\]$",
    variables,
    data.frame(group = integer(), coefficient = integer())
  )
  draws <- posterior::as_draws_matrix(fit$draws(variables = variables))
  data.frame(
    case = case, mode = mode, index,
    mean_rho = colMeans(draws)
  )
}

scalar_rho_error <- function(fit, simulated, rho) {
  if (!length(simulated$counts) || all(is.na(rho$group))) {
    return(NA_real_)
  }
  variables <- sprintf(
    "rho_s2z_1[%d,1]", seq_along(simulated$counts)
  )
  rho_draws <- posterior::as_draws_matrix(
    fit$draws(variables = variables)
  )
  hyper <- posterior::as_draws_matrix(
    fit$draws(variables = c("sd_1[1]", "sigma"))
  )
  scale_ratio <- as.numeric(hyper[, "sd_1[1]"])^2 /
    as.numeric(hyper[, "sigma"])^2
  information <- outer(scale_ratio, simulated$counts)
  expected <- information / (1 + information)
  max(abs(as.matrix(rho_draws) - expected))
}

fit_metrics <- function(fit, case, mode, elapsed, compile_elapsed,
                        simulated, rho) {
  physical <- physical_variables(fit)
  if (!length(physical$all) || !length(physical$hyper)) {
    stop("No physical posterior variables were found.", call. = FALSE)
  }
  summary <- posterior::summarise_draws(
    fit$draws(variables = physical$all),
    mean, sd, mcse_mean, rhat, ess_bulk, ess_tail
  )
  hyper <- summary[summary$variable %in% physical$hyper, , drop = FALSE]
  group <- summary[summary$variable %in% physical$group, , drop = FALSE]
  diagnostics <- fit$sampler_diagnostics(format = "draws_array")
  diagnostic_variables <- posterior::variables(diagnostics)
  leapfrog_name <- intersect(
    c("n_leapfrog__", "num_steps__"), diagnostic_variables
  )[1L]
  if (is.na(leapfrog_name)) {
    stop("Leapfrog diagnostic was not found.", call. = FALSE)
  }
  leapfrog <- diagnostics[, , leapfrog_name]
  postwarmup_gradients <- sum(leapfrog + 1)
  energy <- diagnostics[, , "energy__", drop = FALSE]
  chain_ebfmi <- vapply(seq_len(dim(energy)[2L]), function(chain) {
    ebfmi(energy[, chain, 1L])
  }, numeric(1))
  rho_values <- rho$mean_rho[is.finite(rho$mean_rho)]
  data.frame(
    case = case,
    mode = mode,
    observations = nrow(simulated$data),
    parameters = nrow(summary),
    compile_seconds = compile_elapsed,
    sample_seconds = elapsed,
    postwarmup_gradients = postwarmup_gradients,
    mean_postwarmup_leapfrog = mean(leapfrog),
    divergences = sum(diagnostics[, , "divergent__"]),
    treedepth_hits = sum(diagnostics[, , "treedepth__"] >= max_treedepth),
    min_ebfmi = min(chain_ebfmi, na.rm = TRUE),
    max_physical_rhat = max(summary$rhat, na.rm = TRUE),
    min_physical_bulk_ess = min(summary$ess_bulk, na.rm = TRUE),
    hyper_bulk_ess = stats::median(hyper$ess_bulk, na.rm = TRUE),
    group_bulk_ess = stats::median(group$ess_bulk, na.rm = TRUE),
    hyper_ess_per_second =
      stats::median(hyper$ess_bulk, na.rm = TRUE) / elapsed,
    group_ess_per_second =
      stats::median(group$ess_bulk, na.rm = TRUE) / elapsed,
    hyper_ess_per_1000_postwarmup_gradients =
      1000 * stats::median(hyper$ess_bulk, na.rm = TRUE) /
      postwarmup_gradients,
    group_ess_per_1000_postwarmup_gradients =
      1000 * stats::median(group$ess_bulk, na.rm = TRUE) /
      postwarmup_gradients,
    rho_min = min(rho_values),
    rho_median = stats::median(rho_values),
    rho_max = max(rho_values),
    scalar_rho_error = scalar_rho_error(fit, simulated, rho),
    stringsAsFactors = FALSE
  )
}

physical_summary <- function(fit, case, mode) {
  physical <- physical_variables(fit)
  out <- posterior::summarise_draws(
    fit$draws(variables = physical$all), mean, sd, mcse_mean
  )
  out$parameter_class <- ifelse(
    out$variable %in% physical$hyper, "hyper", "group"
  )
  out$case <- case
  out$mode <- mode
  as.data.frame(out)
}

cat("Ordinary S2Z Fisher benchmark\n")
cat("  brms:", as.character(packageVersion("brms")), "\n")
cat("  CmdStan:", paste(cmdstanr::cmdstan_version(), collapse = "."), "\n")
cat("  chains/warmup/sampling:", chains, warmup, sampling, "\n")
cat("  cases:", paste(cases, collapse = ", "), "\n")
cat("  modes:", paste(modes, collapse = ", "), "\n\n")

metric_rows <- list()
rho_rows <- list()
summary_rows <- list()
job <- 0L
for (case in cases) {
  simulated <- simulate_case(case)
  for (mode in modes) {
    job <- job + 1L
    label <- paste(case, mode, sep = "__")
    cat("[", job, "/", length(cases) * length(modes), "] ", label, "\n", sep = "")
    formula <- case_formula(case, mode)
    code <- stancode(
      formula, data = simulated$data, family = simulated$family,
      prior = simulated$prior
    )
    sdata <- standata(
      formula, data = simulated$data, family = simulated$family,
      prior = simulated$prior
    )
    stan_file <- cmdstanr::write_stan_file(
      code, dir = file.path(output_dir, "stan")
    )
    compile_time <- system.time({
      model <- cmdstanr::cmdstan_model(stan_file, quiet = TRUE)
    })[["elapsed"]]
    csv_dir <- file.path(output_dir, "csv", label)
    dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
    sample_time <- system.time({
      fit <- model$sample(
        data = sdata,
        seed = 91000L + job,
        chains = chains,
        parallel_chains = min(chains, cores),
        iter_warmup = warmup,
        iter_sampling = sampling,
        adapt_delta = adapt_delta,
        max_treedepth = max_treedepth,
        refresh = 0,
        output_dir = csv_dir
      )
    })[["elapsed"]]
    rho <- rho_summary(fit, case, mode)
    metric_rows[[label]] <- fit_metrics(
      fit, case, mode, sample_time, compile_time, simulated, rho
    )
    rho_rows[[label]] <- rho
    summary_rows[[label]] <- physical_summary(fit, case, mode)
  }
}

metrics <- do.call(rbind, metric_rows)
rho <- do.call(rbind, rho_rows)
summaries <- do.call(rbind, summary_rows)
rownames(metrics) <- rownames(rho) <- rownames(summaries) <- NULL

agreement <- list()
for (case in unique(summaries$case)) {
  case_summary <- summaries[summaries$case == case, , drop = FALSE]
  if (length(unique(case_summary$mode)) < 2L) {
    next
  }
  mode_pairs <- combn(unique(case_summary$mode), 2L, simplify = FALSE)
  for (pair in mode_pairs) {
    left <- case_summary[case_summary$mode == pair[1L], , drop = FALSE]
    right <- case_summary[case_summary$mode == pair[2L], , drop = FALSE]
    merged <- merge(
      left, right, by = c("variable", "parameter_class"),
      suffixes = c("_left", "_right")
    )
    merged$case <- case
    merged$comparison <- paste(pair, collapse = "_vs_")
    merged$mean_mcse_z <- abs(merged$mean_left - merged$mean_right) /
      sqrt(merged$mcse_mean_left^2 + merged$mcse_mean_right^2)
    agreement[[length(agreement) + 1L]] <- merged[
      , c("case", "comparison", "parameter_class", "variable", "mean_mcse_z")
    ]
  }
}
if (length(agreement)) {
  agreement <- do.call(rbind, agreement)
} else {
  agreement <- data.frame(
    case = character(), comparison = character(), parameter_class = character(),
    variable = character(), mean_mcse_z = numeric()
  )
}
rownames(agreement) <- NULL

utils::write.csv(metrics, file.path(output_dir, "metrics.csv"), row.names = FALSE)
utils::write.csv(rho, file.path(output_dir, "rho.csv"), row.names = FALSE)
utils::write.csv(
  summaries, file.path(output_dir, "physical-summaries.csv"), row.names = FALSE
)
utils::write.csv(
  summaries[summaries$parameter_class == "hyper", , drop = FALSE],
  file.path(output_dir, "hyper-summaries.csv"), row.names = FALSE
)
utils::write.csv(
  agreement, file.path(output_dir, "posterior-agreement.csv"), row.names = FALSE
)

cat("\nGeometry summary\n")
print(metrics[, c(
  "case", "mode", "postwarmup_gradients", "divergences",
  "treedepth_hits", "min_ebfmi", "max_physical_rhat",
  "hyper_ess_per_second", "group_ess_per_second",
  "hyper_ess_per_1000_postwarmup_gradients",
  "group_ess_per_1000_postwarmup_gradients", "rho_min", "rho_median",
  "rho_max"
)], row.names = FALSE, digits = 4)

cat("\nMaximum posterior-mean discrepancy in combined MCSE units\n")
if (nrow(agreement)) {
  agreement_max <- aggregate(
    mean_mcse_z ~ case + comparison + parameter_class,
    agreement, max, na.rm = TRUE
  )
  print(agreement_max, row.names = FALSE, digits = 4)
} else {
  cat("  not computed (only one mode was requested)\n")
}
cat("\nResults written to", output_dir, "\n")
