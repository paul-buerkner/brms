# End-to-end sampling validation for multiple S2Z group-effect blocks.
#
# This is intentionally a local test because its default run compiles and
# samples seven CmdStan models. Run it against a freshly installed development
# build of brms, for example:
#
#   lib="$(mktemp -d)"
#   R CMD INSTALL --library="$lib" .
#   R_LIBS="$lib" Rscript tests/local/tests.models-s2z-multiblock.R
#
# Defaults and optional controls:
#
#   BRMS_S2Z_CHAINS=2
#   BRMS_S2Z_WARMUP=500
#   BRMS_S2Z_SAMPLING=500
#   BRMS_S2Z_CORES=min(chains, parallel::detectCores())
#   BRMS_S2Z_ADAPT_DELTA=0.98
#   BRMS_S2Z_MAX_TREEDEPTH=12
#   BRMS_S2Z_EQUIV_Z=7
#   BRMS_S2Z_MAX_RHAT=1.10
#   BRMS_S2Z_MAX_DIVERGENCES=0
#   BRMS_S2Z_MAX_TREEDEPTH_HITS=0
#   BRMS_S2Z_MIN_EBFMI=0.20
#   BRMS_S2Z_INVARIANT_TOL=1e-9
#   BRMS_S2Z_RUN_MIXED=true
#   BRMS_S2Z_CACHE_DIR=
#
# Posterior equivalence is checked for every public population coefficient,
# group-level SD/correlation/df, group effect, and selected linear predictor.
# Both means and SDs are compared using the propagated MCSE from the two
# independent fits. Exact draw-wise checks verify the joint recovery map:
#
#   theta_s2z = b + sum(block means),
#   r_s2z_b = r_b - mean(r_b),
#   theta_s2z + sum(r_s2z_b) = b + sum(r_b),
#
# along with the resulting observation-level linear predictor.

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

env_flag <- function(name, default = TRUE) {
  value <- tolower(Sys.getenv(name, unset = if (default) "true" else "false"))
  if (!value %in% c("true", "false", "1", "0", "yes", "no")) {
    stop(name, " must be true or false.", call. = FALSE)
  }
  value %in% c("true", "1", "yes")
}

if (!requireNamespace("cmdstanr", quietly = TRUE) ||
    !nzchar(tryCatch(cmdstanr::cmdstan_path(), error = function(e) ""))) {
  stop("This validation requires a working CmdStan installation.", call. = FALSE)
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
adapt_delta <- env_number("BRMS_S2Z_ADAPT_DELTA", 0.98, minimum = 0.5)
max_treedepth <- env_integer("BRMS_S2Z_MAX_TREEDEPTH", 12L)
equiv_z <- env_number("BRMS_S2Z_EQUIV_Z", 7, minimum = 1)
max_rhat <- env_number("BRMS_S2Z_MAX_RHAT", 1.10, minimum = 1)
max_divergences <- env_integer(
  "BRMS_S2Z_MAX_DIVERGENCES", 0L, minimum = 0L
)
max_treedepth_hits <- env_integer(
  "BRMS_S2Z_MAX_TREEDEPTH_HITS", 0L, minimum = 0L
)
min_ebfmi <- env_number("BRMS_S2Z_MIN_EBFMI", 0.20, minimum = 0)
invariant_tolerance <- env_number(
  "BRMS_S2Z_INVARIANT_TOL", 1e-9, minimum = 0
)
run_mixed <- env_flag("BRMS_S2Z_RUN_MIXED", TRUE)
cache_dir <- Sys.getenv("BRMS_S2Z_CACHE_DIR", unset = "")

options(mc.cores = cores, width = 160)

configuration <- data.frame(
  brms_version = as.character(utils::packageVersion("brms")),
  brms_library = normalizePath(find.package("brms")),
  cmdstan_version = as.character(cmdstanr::cmdstan_version()),
  chains = chains,
  warmup = warmup,
  sampling = sampling,
  cores = cores,
  sig_figs = 17L,
  adapt_delta = adapt_delta,
  max_treedepth = max_treedepth,
  run_mixed = run_mixed,
  stringsAsFactors = FALSE
)
cat("Multiblock S2Z sampling-equivalence configuration\n")
print(configuration, row.names = FALSE)

rmvn_rows <- function(n, sigma) {
  matrix(stats::rnorm(n * nrow(sigma)), nrow = n) %*% chol(sigma)
}

balanced_factor <- function(n, size, prefix) {
  values <- rep(sprintf("%s%02d", prefix, seq_len(size)), length.out = n)
  factor(sample(values), levels = sprintf("%s%02d", prefix, seq_len(size)))
}

simulate_three_intercepts <- function(seed = 7301L, n = 360L) {
  set.seed(seed)
  dat <- data.frame(
    g = balanced_factor(n, 8L, "g"),
    h = balanced_factor(n, 9L, "h"),
    site = balanced_factor(n, 10L, "s")
  )
  ug <- stats::rnorm(nlevels(dat$g), sd = 0.70)
  uh <- stats::rnorm(nlevels(dat$h), sd = 0.50)
  us <- stats::rnorm(nlevels(dat$site), sd = 0.35)
  eta <- 0.40 + ug[dat$g] + uh[dat$h] + us[dat$site]
  dat$y <- stats::rnorm(n, eta, 0.60)
  dat
}

simulate_overlapping_interactions <- function(seed = 7302L, n = 480L) {
  set.seed(seed)
  dat <- data.frame(
    g = balanced_factor(n, 10L, "g"),
    h = balanced_factor(n, 9L, "h"),
    x = stats::rnorm(n),
    z = sample(c(-1, 1), n, replace = TRUE)
  )
  X <- stats::model.matrix(~ x * z, dat)
  beta <- c(0.35, 0.75, -0.45, 0.40)
  g_sd <- c(0.65, 0.40, 0.35, 0.30)
  g_cor <- outer(seq_along(g_sd), seq_along(g_sd), function(i, j) {
    0.28^abs(i - j)
  })
  ug <- rmvn_rows(
    nlevels(dat$g), diag(g_sd) %*% g_cor %*% diag(g_sd)
  )
  h_sd <- c(0.45, 0.30, 0.28, 0.22)
  uh <- matrix(
    stats::rnorm(nlevels(dat$h) * length(h_sd)), ncol = length(h_sd)
  )
  uh <- sweep(uh, 2L, h_sd, "*")
  eta <- drop(X %*% beta) +
    rowSums(X * ug[as.integer(dat$g), ]) +
    rowSums(X * uh[as.integer(dat$h), ])
  dat$y <- stats::rnorm(n, eta, 0.65)
  dat
}

simulate_gaussian_student <- function(seed = 7303L, n = 360L,
                                      group_df = 5) {
  set.seed(seed)
  dat <- data.frame(
    g = balanced_factor(n, 10L, "g"),
    h = balanced_factor(n, 11L, "h"),
    x = stats::rnorm(n)
  )
  X <- stats::model.matrix(~ x, dat)
  beta <- c(-0.25, 0.70)
  g_sd <- c(0.60, 0.35)
  g_cor <- matrix(c(1, 0.25, 0.25, 1), nrow = 2L)
  ug <- rmvn_rows(
    nlevels(dat$g), diag(g_sd) %*% g_cor %*% diag(g_sd)
  )
  h_sd <- c(0.50, 0.30)
  uh <- matrix(
    stats::rnorm(nlevels(dat$h) * length(h_sd)), ncol = length(h_sd)
  )
  uh <- sweep(uh, 2L, h_sd, "*")
  uh <- uh * sqrt(group_df / stats::rchisq(nlevels(dat$h), df = group_df))
  eta <- drop(X %*% beta) +
    rowSums(X * ug[as.integer(dat$g), ]) +
    rowSums(X * uh[as.integer(dat$h), ])
  dat$y <- stats::rnorm(n, eta, 0.60)
  dat
}

base_prior <- function(correlated = TRUE) {
  out <- prior(normal(0, 1.5), class = "Intercept") +
    prior(normal(0, 1.0), class = "b") +
    prior(exponential(1), class = "sd") +
    prior(exponential(1), class = "sigma")
  if (correlated) {
    out <- out + prior(lkj(2), class = "cor")
  }
  out
}

cache_file <- function(label) {
  if (!nzchar(cache_dir)) {
    return(NULL)
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  suffix <- paste(chains, warmup, sampling, "sig17", sep = "-")
  file.path(cache_dir, paste0(label, "-", suffix, ".rds"))
}

sample_model <- function(formula, data, prior, seed, label) {
  cat("Sampling ", label, " ...\n", sep = "")
  args <- list(
    formula = formula,
    data = data,
    family = gaussian(),
    prior = prior,
    backend = "cmdstanr",
    chains = chains,
    cores = cores,
    iter = warmup + sampling,
    warmup = warmup,
    seed = seed,
    refresh = 0,
    silent = 1,
    sig_figs = 17,
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
  timing <- system.time(fit <- do.call(brm, args))
  list(fit = fit, elapsed = unname(timing[["elapsed"]]))
}

exact_draw_matrix <- function(fit, variable) {
  missing <- setdiff(variable, variables(fit))
  if (length(missing)) {
    stop(
      "Required saved variables are missing: ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  as.matrix(as_draws_matrix(fit, variable = variable))
}

as_draw_column <- function(x) {
  out <- as.numeric(x)
  if (!length(out) || anyNA(out)) {
    stop("A public posterior quantity has no finite saved draws.", call. = FALSE)
  }
  out
}

public_samples <- function(fit, groups, include_eta = TRUE) {
  fe <- as.matrix(fixef(fit, summary = FALSE))
  fe <- fe[, !grepl("^s2z(\\[|$)", colnames(fe)), drop = FALSE]
  out <- setNames(
    lapply(seq_len(ncol(fe)), function(k) as_draw_column(fe[, k])),
    paste0("b:", colnames(fe))
  )

  fit_variables <- variables(fit)
  for (group in groups) {
    hyper <- c(
      grep(paste0("^sd_", group, "__"), fit_variables, value = TRUE),
      grep(paste0("^cor_", group, "__"), fit_variables, value = TRUE),
      grep(paste0("^df_", group, "$"), fit_variables, value = TRUE)
    )
    if (length(hyper)) {
      draws <- exact_draw_matrix(fit, unique(hyper))
      for (name in colnames(draws)) {
        out[[paste0("parameter:", name)]] <- as_draw_column(draws[, name])
      }
    }
  }
  if ("sigma" %in% fit_variables) {
    out[["parameter:sigma"]] <- as_draw_column(
      exact_draw_matrix(fit, "sigma")[, 1L]
    )
  }

  re <- ranef(fit, summary = FALSE)
  for (group in groups) {
    block <- re[[group]]
    if (is.null(block)) {
      stop("No public group effects found for ", group, call. = FALSE)
    }
    for (j in seq_len(dim(block)[2L])) {
      for (k in seq_len(dim(block)[3L])) {
        label <- sprintf(
          "r:%s[%s,%s]", group,
          dimnames(block)[[2L]][j], dimnames(block)[[3L]][k]
        )
        out[[label]] <- as_draw_column(block[, j, k])
      }
    }
  }

  if (include_eta) {
    eta <- posterior_linpred(fit, transform = FALSE)
    selected <- unique(round(seq(1, ncol(eta), length.out = 7L)))
    for (i in selected) {
      out[[sprintf("eta:observation[%d]", i)]] <- eta[, i]
    }
    out[["eta:observation-average"]] <- rowMeans(eta)
  }
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
    stop("Could not estimate an MCSE for a public quantity.", call. = FALSE)
  }
  ans
}

compare_public_posteriors <- function(conventional, s2z, case) {
  if (!setequal(names(conventional), names(s2z))) {
    stop(
      "Public quantities differ between conventional and S2Z fits.\n",
      "Only conventional: ",
      paste(setdiff(names(conventional), names(s2z)), collapse = ", "),
      "\nOnly S2Z: ",
      paste(setdiff(names(s2z), names(conventional)), collapse = ", "),
      call. = FALSE
    )
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

internal_group_draws <- function(fit, id, n_group, n_coef) {
  fit_variables <- variables(fit)
  out <- NULL
  for (k in seq_len(n_coef)) {
    matrix_names <- sprintf(
      "r_s2z_%d[%d,%d]", id, seq_len(n_group), k
    )
    vector_names <- sprintf(
      "r_s2z_%d_%d[%d]", id, k, seq_len(n_group)
    )
    names <- if (all(matrix_names %in% fit_variables)) {
      matrix_names
    } else if (all(vector_names %in% fit_variables)) {
      vector_names
    } else {
      stop(
        "Could not find internal S2Z draws for group-level ID ", id,
        ", coefficient ", k, ".", call. = FALSE
      )
    }
    draws <- exact_draw_matrix(fit, names)
    if (is.null(out)) {
      out <- array(NA_real_, dim = c(nrow(draws), n_group, n_coef))
    }
    out[, , k] <- draws
  }
  out
}

draw_group_means <- function(x) {
  n_draw <- dim(x)[1L]
  n_coef <- dim(x)[3L]
  out <- matrix(NA_real_, nrow = n_draw, ncol = n_coef)
  for (k in seq_len(n_coef)) {
    out[, k] <- rowMeans(matrix(x[, , k], nrow = n_draw))
  }
  out
}

joint_coordinate_invariants <- function(fit, data, fixed_formula,
                                        groups, ids, case) {
  if (length(groups) != length(ids)) {
    stop("groups and ids must have equal length.", call. = FALSE)
  }
  X <- stats::model.matrix(fixed_formula, data)
  coef_names <- sub("^\\(Intercept\\)$", "Intercept", colnames(X))
  n_coef <- ncol(X)
  theta_names <- sprintf("theta_s2z[%d]", seq_len(n_coef))
  theta <- exact_draw_matrix(fit, theta_names)
  n_draw <- nrow(theta)

  # brms centers the population design matrix internally. Convert its
  # intercept coordinate to the raw model-matrix intercept used by fixef().
  theta_raw <- theta
  if ("(Intercept)" %in% colnames(X) && n_coef > 1L) {
    means_X <- colMeans(X[, -1L, drop = FALSE])
    theta_raw[, 1L] <- theta[, 1L] -
      drop(theta[, -1L, drop = FALSE] %*% means_X)
  }

  fe_all <- as.matrix(fixef(fit, summary = FALSE))
  fe <- fe_all[, !grepl("^s2z(\\[|$)", colnames(fe_all)), drop = FALSE]
  if (!all(coef_names %in% colnames(fe))) {
    stop("Could not align public fixed-effect names.", call. = FALSE)
  }
  fe <- fe[, coef_names, drop = FALSE]

  public_re <- ranef(fit, summary = FALSE)
  internal <- public <- means <- indices <- vector("list", length(groups))
  names(internal) <- names(public) <- names(means) <- names(indices) <- groups
  max_sum_to_zero <- 0
  max_centered_recovery <- 0

  for (b in seq_along(groups)) {
    group <- groups[b]
    block <- public_re[[group]]
    if (is.null(block) || !all(coef_names %in% dimnames(block)[[3L]])) {
      stop("Could not align public effects for group ", group, ".",
           call. = FALSE)
    }
    block <- block[, , coef_names, drop = FALSE]
    public[[b]] <- block
    internal[[b]] <- internal_group_draws(
      fit, ids[b], dim(block)[2L], n_coef
    )
    means[[b]] <- draw_group_means(block)
    indices[[b]] <- match(
      as.character(data[[group]]), dimnames(block)[[2L]]
    )
    if (anyNA(indices[[b]])) {
      stop("Could not align observations with group ", group, ".",
           call. = FALSE)
    }
    for (k in seq_len(n_coef)) {
      internal_k <- matrix(
        internal[[b]][, , k], nrow = n_draw, ncol = dim(block)[2L]
      )
      public_k <- matrix(
        block[, , k], nrow = n_draw, ncol = dim(block)[2L]
      )
      centered_public <- sweep(public_k, 1L, means[[b]][, k], "-")
      max_sum_to_zero <- max(
        max_sum_to_zero, max(abs(rowSums(internal_k)))
      )
      max_centered_recovery <- max(
        max_centered_recovery, max(abs(internal_k - centered_public))
      )
    }
  }

  sum_means <- Reduce(`+`, means)
  max_fixed_recovery <- max(abs(theta_raw - fe - sum_means))

  eta_internal <- theta_raw %*% t(X)
  eta_public <- fe %*% t(X)
  max_finite_recovery <- 0
  for (k in seq_len(n_coef)) {
    finite_internal <- matrix(theta_raw[, k], nrow = n_draw, ncol = nrow(X))
    finite_public <- matrix(fe[, k], nrow = n_draw, ncol = nrow(X))
    for (b in seq_along(groups)) {
      finite_internal <- finite_internal + matrix(
        internal[[b]][, indices[[b]], k],
        nrow = n_draw, ncol = nrow(X)
      )
      finite_public <- finite_public + matrix(
        public[[b]][, indices[[b]], k],
        nrow = n_draw, ncol = nrow(X)
      )
    }
    max_finite_recovery <- max(
      max_finite_recovery, max(abs(finite_internal - finite_public))
    )
    eta_internal <- eta_internal + sweep(
      finite_internal - matrix(theta_raw[, k], nrow = n_draw, ncol = nrow(X)),
      2L, X[, k], "*"
    )
    eta_public <- eta_public + sweep(
      finite_public - matrix(fe[, k], nrow = n_draw, ncol = nrow(X)),
      2L, X[, k], "*"
    )
  }
  eta_brms <- posterior_linpred(fit, transform = FALSE)

  out <- data.frame(
    case = case,
    check = c(
      "each internal group-effect column sums to zero",
      "public r minus its block mean equals internal r_s2z",
      "public b plus all block means equals theta_s2z",
      "public b + sum(r) equals internal theta + sum(r_s2z)",
      "posterior_linpred equals recovered public predictor",
      "posterior_linpred equals internal S2Z predictor"
    ),
    max_abs_error = c(
      max_sum_to_zero,
      max_centered_recovery,
      max_fixed_recovery,
      max_finite_recovery,
      max(abs(eta_brms - eta_public)),
      max(abs(eta_brms - eta_internal))
    ),
    tolerance = invariant_tolerance,
    stringsAsFactors = FALSE
  )
  out$pass <- out$max_abs_error <= out$tolerance
  out
}

public_variable_names <- function(fit, groups) {
  fit_variables <- variables(fit)
  out <- grep("^b_", fit_variables, value = TRUE)
  for (group in groups) {
    out <- c(
      out,
      grep(paste0("^sd_", group, "__"), fit_variables, value = TRUE),
      grep(paste0("^cor_", group, "__"), fit_variables, value = TRUE),
      grep(paste0("^df_", group, "$"), fit_variables, value = TRUE),
      grep(paste0("^r_", group, "\\["), fit_variables, value = TRUE)
    )
  }
  unique(c(out, intersect("sigma", fit_variables)))
}

mean_cmdstan_chain_time <- function(fit, phase = c("Total", "Sampling")) {
  phase <- match.arg(phase)
  stan_args <- tryCatch(
    methods::slot(fit$fit, "stan_args"),
    error = function(e) NULL
  )
  if (is.null(stan_args)) {
    return(NA_real_)
  }
  times <- vapply(stan_args, function(x) {
    line <- grep(
      paste0("\\(", phase, "\\)"), x$time_info,
      value = TRUE
    )
    if (length(line) != 1L) {
      return(NA_real_)
    }
    fields <- strsplit(
      trimws(sub("^#", "", line)), "[[:space:]]+"
    )[[1L]]
    suppressWarnings(as.numeric(fields[1L]))
  }, numeric(1))
  if (anyNA(times)) NA_real_ else mean(times)
}

fit_quality <- function(fit, groups, case, parameterization, elapsed) {
  public <- public_variable_names(fit, groups)
  if (!length(public)) {
    stop("No public variables found for diagnostics.", call. = FALSE)
  }
  draws <- as_draws_array(fit, variable = public)
  summary <- posterior::summarise_draws(
    draws,
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  )
  np <- nuts_params(fit)
  divergences <- sum(np$Value[np$Parameter == "divergent__"])
  treedepth_hits <- sum(
    np$Value[np$Parameter == "treedepth__"] >= max_treedepth
  )
  energy <- np[np$Parameter == "energy__", , drop = FALSE]
  ebfmi <- vapply(split(energy$Value, energy$Chain), function(x) {
    mean(diff(x)^2) / stats::var(x)
  }, numeric(1))
  min_chain_ebfmi <- min(ebfmi)
  observed_max_rhat <- max(summary$rhat, na.rm = TRUE)
  data.frame(
    case = case,
    parameterization = parameterization,
    elapsed_seconds = elapsed,
    mean_chain_seconds = mean_cmdstan_chain_time(fit, "Total"),
    mean_sampling_seconds = mean_cmdstan_chain_time(fit, "Sampling"),
    max_rhat = observed_max_rhat,
    min_bulk_ess = min(summary$ess_bulk, na.rm = TRUE),
    min_tail_ess = min(summary$ess_tail, na.rm = TRUE),
    divergences = divergences,
    treedepth_hits = treedepth_hits,
    min_chain_ebfmi = min_chain_ebfmi,
    pass = observed_max_rhat <= max_rhat &&
      divergences <= max_divergences &&
      treedepth_hits <= max_treedepth_hits &&
      min_chain_ebfmi >= min_ebfmi,
    stringsAsFactors = FALSE
  )
}

comparison_results <- list()
invariant_results <- list()
quality_results <- list()

three_data <- simulate_three_intercepts()
three_prior <- prior(normal(0, 1.5), class = "Intercept") +
  prior(exponential(1), class = "sd") +
  prior(exponential(1), class = "sigma")
three_conventional <- sample_model(
  y ~ 1 + (1 | g) + (1 | h) + (1 | site),
  three_data, three_prior, 8301L, "three-intercepts-conventional"
)
three_auto <- sample_model(
  y ~ 1 +
    (1 | gr(g, s2z = TRUE, center = "auto")) +
    (1 | gr(h, s2z = TRUE, center = "auto")) +
    (1 | gr(site, s2z = TRUE, center = "auto")),
  three_data, three_prior, 8302L, "three-intercepts-s2z-auto"
)
three_noncentered <- sample_model(
  y ~ 1 +
    (1 | gr(g, s2z = TRUE, center = FALSE)) +
    (1 | gr(h, s2z = TRUE, center = FALSE)) +
    (1 | gr(site, s2z = TRUE, center = FALSE)),
  three_data, three_prior, 8303L, "three-intercepts-s2z-noncentered"
)
three_groups <- c("g", "h", "site")
comparison_results[["three-auto"]] <- compare_public_posteriors(
  public_samples(three_conventional$fit, three_groups),
  public_samples(three_auto$fit, three_groups),
  "three-intercepts-s2z-auto"
)
comparison_results[["three-noncentered"]] <- compare_public_posteriors(
  public_samples(three_conventional$fit, three_groups),
  public_samples(three_noncentered$fit, three_groups),
  "three-intercepts-s2z-noncentered"
)
invariant_results[["three-auto"]] <- joint_coordinate_invariants(
  three_auto$fit, three_data, ~ 1, three_groups, 1:3,
  "three-intercepts-s2z-auto"
)
invariant_results[["three-noncentered"]] <- joint_coordinate_invariants(
  three_noncentered$fit, three_data, ~ 1, three_groups, 1:3,
  "three-intercepts-s2z-noncentered"
)
quality_results[["three-conventional"]] <- fit_quality(
  three_conventional$fit, three_groups, "three-intercepts", "conventional",
  three_conventional$elapsed
)
quality_results[["three-auto"]] <- fit_quality(
  three_auto$fit, three_groups, "three-intercepts", "s2z-auto",
  three_auto$elapsed
)
quality_results[["three-noncentered"]] <- fit_quality(
  three_noncentered$fit, three_groups, "three-intercepts", "s2z-noncentered",
  three_noncentered$elapsed
)

interaction_data <- simulate_overlapping_interactions()
interaction_prior <- base_prior(correlated = TRUE)
interaction_conventional <- sample_model(
  y ~ x * z + (1 + x * z || g) + (1 + x * z | h),
  interaction_data, interaction_prior, 8401L,
  "overlapping-interactions-conventional"
)
interaction_auto <- sample_model(
  y ~ x * z +
    (1 + x * z || gr(g, s2z = TRUE, center = "auto")) +
    (1 + x * z | gr(h, s2z = TRUE, center = "auto")),
  interaction_data, interaction_prior, 8402L,
  "overlapping-interactions-s2z-auto"
)
interaction_groups <- c("g", "h")
comparison_results[["interaction-auto"]] <- compare_public_posteriors(
  public_samples(interaction_conventional$fit, interaction_groups),
  public_samples(interaction_auto$fit, interaction_groups),
  "overlapping-interactions-s2z-auto"
)
invariant_results[["interaction-auto"]] <- joint_coordinate_invariants(
  interaction_auto$fit, interaction_data, ~ x * z,
  interaction_groups, 1:2, "overlapping-interactions-s2z-auto"
)
quality_results[["interaction-conventional"]] <- fit_quality(
  interaction_conventional$fit, interaction_groups,
  "overlapping-interactions", "conventional",
  interaction_conventional$elapsed
)
quality_results[["interaction-auto"]] <- fit_quality(
  interaction_auto$fit, interaction_groups,
  "overlapping-interactions", "s2z-auto", interaction_auto$elapsed
)

if (run_mixed) {
  mixed_data <- simulate_gaussian_student()
  mixed_prior <- base_prior(correlated = TRUE) +
    prior(gamma(2, 0.2), class = "df", group = "h")
  mixed_conventional <- sample_model(
    y ~ x +
      (1 + x | gr(g, dist = "gaussian")) +
      (1 + x || gr(h, dist = "student")),
    mixed_data, mixed_prior, 8501L, "gaussian-student-conventional"
  )
  mixed_auto <- sample_model(
    y ~ x +
      (1 + x | gr(
        g, dist = "gaussian", s2z = TRUE, center = "auto"
      )) +
      (1 + x || gr(
        h, dist = "student", s2z = TRUE, center = "auto"
      )),
    mixed_data, mixed_prior, 8502L, "gaussian-student-s2z-auto"
  )
  mixed_groups <- c("g", "h")
  comparison_results[["mixed-auto"]] <- compare_public_posteriors(
    public_samples(mixed_conventional$fit, mixed_groups),
    public_samples(mixed_auto$fit, mixed_groups),
    "gaussian-student-s2z-auto"
  )
  invariant_results[["mixed-auto"]] <- joint_coordinate_invariants(
    mixed_auto$fit, mixed_data, ~ x, mixed_groups, 1:2,
    "gaussian-student-s2z-auto"
  )
  quality_results[["mixed-conventional"]] <- fit_quality(
    mixed_conventional$fit, mixed_groups, "gaussian-student",
    "conventional", mixed_conventional$elapsed
  )
  quality_results[["mixed-auto"]] <- fit_quality(
    mixed_auto$fit, mixed_groups, "gaussian-student", "s2z-auto",
    mixed_auto$elapsed
  )
}

comparison_results <- do.call(rbind, comparison_results)
invariant_results <- do.call(rbind, invariant_results)
quality_results <- do.call(rbind, quality_results)

cat("\nMCSE-aware conventional-vs-S2Z comparisons (largest z-scores first)\n")
print(
  comparison_results[order(comparison_results$z_score, decreasing = TRUE), ],
  row.names = FALSE, digits = 5
)
cat("\nExact joint recovery invariants\n")
print(invariant_results, row.names = FALSE, digits = 5)
cat("\nSampling diagnostics and wall-clock timings\n")
print(quality_results, row.names = FALSE, digits = 5)

speed_results <- do.call(rbind, lapply(
  which(quality_results$parameterization != "conventional"),
  function(i) {
    conventional <- subset(
      quality_results,
      case == quality_results$case[i] & parameterization == "conventional"
    )
    data.frame(
      case = quality_results$case[i],
      parameterization = quality_results$parameterization[i],
      conventional_mean_chain_seconds = conventional$mean_chain_seconds,
      s2z_mean_chain_seconds = quality_results$mean_chain_seconds[i],
      chain_time_speedup = conventional$mean_chain_seconds /
        quality_results$mean_chain_seconds[i],
      stringsAsFactors = FALSE
    )
  }
))
cat("\nCmdStan mean-chain timing ratios (warmup plus sampling)\n")
print(speed_results, row.names = FALSE, digits = 5)

comparison_summary <- do.call(rbind, lapply(
  split(comparison_results, comparison_results$case),
  function(x) data.frame(
    case = x$case[1L],
    quantities = nrow(x) / 2L,
    mean_and_sd_checks = nrow(x),
    max_z_score = max(x$z_score),
    failed = sum(!x$pass),
    stringsAsFactors = FALSE
  )
))
cat("\nComparison summary\n")
print(comparison_summary, row.names = FALSE, digits = 5)

failed_comparisons <- subset(comparison_results, !pass)
failed_invariants <- subset(invariant_results, !pass)
failed_quality <- subset(quality_results, !pass)
if (nrow(failed_comparisons) || nrow(failed_invariants) ||
    nrow(failed_quality)) {
  stop(
    "Multiblock S2Z validation failed: ",
    nrow(failed_comparisons), " MCSE comparisons, ",
    nrow(failed_invariants), " exact invariants, and ",
    nrow(failed_quality), " sampling-quality checks failed.",
    call. = FALSE
  )
}

cat(
  "\nPASS: all public posterior means and SDs agree within propagated MCSE; ",
  "joint b/r recovery identities hold draw by draw; and all fits meet the ",
  "configured sampler checks.\n",
  sep = ""
)
