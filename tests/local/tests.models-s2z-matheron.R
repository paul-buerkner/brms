# End-to-end validation and benchmark for the Gaussian multiblock Matheron
# specialization. This is a local test because it compiles and samples four
# CmdStan models. Run it against a freshly installed source tree, for example:
#
#   lib="$(mktemp -d)"
#   R CMD INSTALL --library="$lib" .
#   R_LIBS="$lib" Rscript tests/local/tests.models-s2z-matheron.R
#
# The two cases deliberately exercise both Matheron implementations:
#
# * 12 crossed varying-intercept factors use the scalar W specialization;
# * two overlapping (intercept, slopes, interaction) blocks use a 4 by 4 W.
#
# Environment controls (defaults shown):
#
#   BRMS_S2Z_MATHERON_CHAINS=2
#   BRMS_S2Z_MATHERON_WARMUP=400
#   BRMS_S2Z_MATHERON_SAMPLING=400
#   BRMS_S2Z_MATHERON_SCALAR_FACTORS=12
#   BRMS_S2Z_MATHERON_EQUIV_Z=8
#   BRMS_S2Z_MATHERON_MAX_RHAT=1.10
#   BRMS_S2Z_MATHERON_MAX_DIVERGENCES=0
#   BRMS_S2Z_MATHERON_MAX_TREEDEPTH_HITS=0
#   BRMS_S2Z_MATHERON_MIN_EBFMI=0.20
#   BRMS_S2Z_MATHERON_INVARIANT_TOL=1e-9
#
# `gradient_evals` is the sum of n_leapfrog__ + 1 over post-warmup draws.
# Wall time includes compilation for each fit, so the paired wall ratio is a
# validation-run measurement rather than a compile-free production benchmark.

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

if (!requireNamespace("cmdstanr", quietly = TRUE) ||
    !nzchar(tryCatch(cmdstanr::cmdstan_path(), error = function(e) ""))) {
  stop("This validation requires a working CmdStan installation.",
       call. = FALSE)
}

chains <- env_integer("BRMS_S2Z_MATHERON_CHAINS", 2L, minimum = 2L)
warmup <- env_integer("BRMS_S2Z_MATHERON_WARMUP", 400L)
sampling <- env_integer("BRMS_S2Z_MATHERON_SAMPLING", 400L)
scalar_factors <- env_integer(
  "BRMS_S2Z_MATHERON_SCALAR_FACTORS", 12L, minimum = 10L
)
equiv_z <- env_number("BRMS_S2Z_MATHERON_EQUIV_Z", 8, minimum = 1)
max_rhat <- env_number("BRMS_S2Z_MATHERON_MAX_RHAT", 1.10, minimum = 1)
max_divergences <- env_integer(
  "BRMS_S2Z_MATHERON_MAX_DIVERGENCES", 0L, minimum = 0L
)
max_treedepth_hits <- env_integer(
  "BRMS_S2Z_MATHERON_MAX_TREEDEPTH_HITS", 0L, minimum = 0L
)
min_ebfmi <- env_number(
  "BRMS_S2Z_MATHERON_MIN_EBFMI", 0.20, minimum = 0
)
invariant_tolerance <- env_number(
  "BRMS_S2Z_MATHERON_INVARIANT_TOL", 1e-9, minimum = 0
)
adapt_delta <- env_number(
  "BRMS_S2Z_MATHERON_ADAPT_DELTA", 0.98, minimum = 0.5
)
max_treedepth <- env_integer(
  "BRMS_S2Z_MATHERON_MAX_TREEDEPTH", 12L, minimum = 1L
)
detected_cores <- parallel::detectCores(logical = FALSE)
if (is.na(detected_cores)) {
  detected_cores <- 1L
}
cores <- env_integer(
  "BRMS_S2Z_MATHERON_CORES", min(chains, detected_cores), minimum = 1L
)
options(mc.cores = cores, width = 180)

configuration <- data.frame(
  brms_version = as.character(utils::packageVersion("brms")),
  brms_library = normalizePath(find.package("brms")),
  cmdstan_version = paste(cmdstanr::cmdstan_version(), collapse = "."),
  chains = chains,
  warmup = warmup,
  sampling = sampling,
  scalar_factors = scalar_factors,
  cores = cores,
  adapt_delta = adapt_delta,
  max_treedepth = max_treedepth,
  stringsAsFactors = FALSE
)
cat("Gaussian Matheron S2Z validation configuration\n")
print(configuration, row.names = FALSE)

balanced_factor <- function(n, size, prefix) {
  levels <- sprintf("%s%02d", prefix, seq_len(size))
  factor(sample(rep(levels, length.out = n)), levels = levels)
}

rmvn_rows <- function(n, sigma) {
  matrix(stats::rnorm(n * nrow(sigma)), nrow = n) %*% chol(sigma)
}

simulate_crossed_scalar <- function(seed = 9101L, n = 360L,
                                    factors = scalar_factors,
                                    levels = 6L) {
  set.seed(seed)
  group_names <- sprintf("g%02d", seq_len(factors))
  dat <- setNames(lapply(group_names, function(group) {
    balanced_factor(n, levels, paste0(group, "_"))
  }), group_names)
  dat <- as.data.frame(dat)
  eta <- rep(0.25, n)
  for (group in group_names) {
    effect <- stats::rnorm(nlevels(dat[[group]]), sd = 0.25)
    eta <- eta + effect[as.integer(dat[[group]])]
  }
  dat$y <- stats::rnorm(n, eta, 0.55)
  dat
}

simulate_overlapping_vector <- function(seed = 9102L, n = 400L) {
  set.seed(seed)
  dat <- data.frame(
    g = balanced_factor(n, 10L, "g_"),
    h = balanced_factor(n, 11L, "h_"),
    x = stats::rnorm(n),
    z = sample(c(-1, 1), n, replace = TRUE)
  )
  X <- stats::model.matrix(~ x * z, dat)
  beta <- c(0.30, 0.70, -0.40, 0.35)
  g_sd <- c(0.55, 0.35, 0.30, 0.25)
  g_cor <- outer(seq_along(g_sd), seq_along(g_sd), function(i, j) {
    0.25^abs(i - j)
  })
  ug <- rmvn_rows(
    nlevels(dat$g), diag(g_sd) %*% g_cor %*% diag(g_sd)
  )
  h_sd <- c(0.45, 0.30, 0.25, 0.20)
  uh <- matrix(
    stats::rnorm(nlevels(dat$h) * length(h_sd)), ncol = length(h_sd)
  )
  uh <- sweep(uh, 2L, h_sd, "*")
  eta <- drop(X %*% beta) +
    rowSums(X * ug[as.integer(dat$g), ]) +
    rowSums(X * uh[as.integer(dat$h), ])
  dat$y <- stats::rnorm(n, eta, 0.60)
  dat
}

scalar_formulas <- function(data) {
  groups <- grep("^g[0-9]+$", names(data), value = TRUE)
  conventional_terms <- sprintf("(1 | %s)", groups)
  s2z_terms <- sprintf(
    "(1 | gr(%s, s2z = TRUE, center = FALSE))", groups
  )
  list(
    groups = groups,
    conventional = stats::as.formula(paste(
      "y ~ 1 +", paste(conventional_terms, collapse = " + ")
    )),
    s2z = stats::as.formula(paste(
      "y ~ 1 +", paste(s2z_terms, collapse = " + ")
    ))
  )
}

vector_formulas <- list(
  conventional = y ~ x * z + (1 + x * z | g) + (1 + x * z || h),
  s2z = y ~ x * z +
    (1 + x * z | gr(g, s2z = TRUE, center = 0.5)) +
    (1 + x * z || gr(h, s2z = TRUE, center = 0.5))
)

scalar_prior <- prior(normal(0, 1.5), class = "Intercept") +
  prior(exponential(2), class = "sd") +
  prior(exponential(1), class = "sigma")
vector_prior <- prior(normal(0, 1.5), class = "Intercept") +
  prior(normal(0, 1), class = "b") +
  prior(exponential(1.5), class = "sd") +
  prior(lkj(2), class = "cor") +
  prior(exponential(1), class = "sigma")

assert_matheron_code <- function(formula, data, prior, scalar) {
  code <- stancode(
    formula, data = data, family = gaussian(), prior = prior,
    backend = "cmdstanr"
  )
  if (!grepl("W_matheron_s2z", code, fixed = TRUE)) {
    stop("Generated code did not select the Gaussian Matheron path.",
         call. = FALSE)
  }
  if (grepl("P_group_s2z", code, fixed = TRUE) ||
      grepl("P_s2z_", code, fixed = TRUE)) {
    stop("Generated code retained the general omitted-mean precision path.",
         call. = FALSE)
  }
  if (scalar) {
    if (!grepl("real<lower=0> W_matheron_s2z", code, fixed = TRUE) ||
        grepl("L_W_matheron_s2z", code, fixed = TRUE)) {
      stop("The scalar case did not use scalar W without a Cholesky factor.",
           call. = FALSE)
    }
  } else if (!grepl("matrix[4, 4] W_matheron_s2z", code, fixed = TRUE) ||
             !grepl("L_W_matheron_s2z", code, fixed = TRUE)) {
    stop("The overlapping vector case did not use the 4 by 4 W path.",
         call. = FALSE)
  }
  invisible(code)
}

sample_model <- function(formula, data, prior, seed, label) {
  cat("Sampling ", label, " ...\n", sep = "")
  timing <- system.time(fit <- brm(
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
  ))
  list(fit = fit, elapsed = unname(timing[["elapsed"]]))
}

exact_draw_matrix <- function(fit, variable) {
  missing <- setdiff(variable, variables(fit))
  if (length(missing)) {
    stop("Required saved variables are missing: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  as.matrix(as_draws_matrix(fit, variable = variable))
}

as_draw_column <- function(x) {
  out <- as.numeric(x)
  if (!length(out) || any(!is.finite(out))) {
    stop("A public posterior quantity has no finite saved draws.",
         call. = FALSE)
  }
  out
}

public_samples <- function(fit, groups) {
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
      grep(paste0("^cor_", group, "__"), fit_variables, value = TRUE)
    )
    if (length(hyper)) {
      draws <- exact_draw_matrix(fit, unique(hyper))
      for (name in colnames(draws)) {
        out[[paste0("parameter:", name)]] <- as_draw_column(draws[, name])
      }
    }
  }
  out[["parameter:sigma"]] <- as_draw_column(
    exact_draw_matrix(fit, "sigma")[, 1L]
  )

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
  eta <- posterior_linpred(fit, transform = FALSE)
  selected <- unique(round(seq(1, ncol(eta), length.out = 9L)))
  for (i in selected) {
    out[[sprintf("eta:observation[%d]", i)]] <- eta[, i]
  }
  out[["eta:observation-average"]] <- rowMeans(eta)
  out
}

safe_mcse <- function(x, metric) {
  out <- switch(
    metric,
    mean = posterior::mcse_mean(x),
    sd = posterior::mcse_sd(x),
    stop("Unknown metric.", call. = FALSE)
  )
  if (!is.finite(out)) {
    stop("Could not estimate an MCSE.", call. = FALSE)
  }
  out
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
  rows <- list()
  row <- 0L
  for (name in sort(names(conventional))) {
    x <- conventional[[name]]
    y <- s2z[[name]]
    for (metric in c("mean", "sd")) {
      row <- row + 1L
      estimate_x <- if (metric == "mean") mean(x) else stats::sd(x)
      estimate_y <- if (metric == "mean") mean(y) else stats::sd(y)
      mcse_difference <- sqrt(
        safe_mcse(x, metric)^2 + safe_mcse(y, metric)^2
      )
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
    matrix_names <- sprintf("r_s2z_%d[%d,%d]", id, seq_len(n_group), k)
    vector_names <- sprintf(
      "r_s2z_%d_%d[%d]", id, k, seq_len(n_group)
    )
    names <- if (all(matrix_names %in% fit_variables)) {
      matrix_names
    } else if (all(vector_names %in% fit_variables)) {
      vector_names
    } else {
      stop("Could not find internal S2Z draws for ID ", id,
           ", coefficient ", k, ".", call. = FALSE)
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
  out <- matrix(NA_real_, nrow = dim(x)[1L], ncol = dim(x)[3L])
  for (k in seq_len(dim(x)[3L])) {
    out[, k] <- rowMeans(matrix(x[, , k], nrow = dim(x)[1L]))
  }
  out
}

joint_coordinate_invariants <- function(fit, data, fixed_formula,
                                        groups, ids, case) {
  X <- stats::model.matrix(fixed_formula, data)
  coef_names <- sub("^\\(Intercept\\)$", "Intercept", colnames(X))
  n_coef <- ncol(X)
  theta <- exact_draw_matrix(
    fit, sprintf("theta_s2z[%d]", seq_len(n_coef))
  )
  n_draw <- nrow(theta)

  theta_raw <- theta
  if ("(Intercept)" %in% colnames(X) && n_coef > 1L) {
    theta_raw[, 1L] <- theta[, 1L] -
      drop(theta[, -1L, drop = FALSE] %*% colMeans(X[, -1L, drop = FALSE]))
  }
  fe <- as.matrix(fixef(fit, summary = FALSE))
  fe <- fe[, !grepl("^s2z(\\[|$)", colnames(fe)), drop = FALSE]
  fe <- fe[, coef_names, drop = FALSE]

  public_re <- ranef(fit, summary = FALSE)
  internal <- public <- means <- indices <- vector("list", length(groups))
  max_sum_to_zero <- max_centered_recovery <- 0
  for (b in seq_along(groups)) {
    group <- groups[b]
    block <- public_re[[group]][, , coef_names, drop = FALSE]
    public[[b]] <- block
    internal[[b]] <- internal_group_draws(
      fit, ids[b], dim(block)[2L], n_coef
    )
    means[[b]] <- draw_group_means(block)
    indices[[b]] <- match(
      as.character(data[[group]]), dimnames(block)[[2L]]
    )
    if (anyNA(indices[[b]])) {
      stop("Could not align group levels for ", group, call. = FALSE)
    }
    for (k in seq_len(n_coef)) {
      internal_k <- matrix(
        internal[[b]][, , k], nrow = n_draw, ncol = dim(block)[2L]
      )
      public_k <- matrix(
        block[, , k], nrow = n_draw, ncol = dim(block)[2L]
      )
      centered_public <- sweep(public_k, 1L, means[[b]][, k], "-")
      max_sum_to_zero <- max(max_sum_to_zero, max(abs(rowSums(internal_k))))
      max_centered_recovery <- max(
        max_centered_recovery, max(abs(internal_k - centered_public))
      )
    }
  }

  max_fixed_recovery <- max(abs(theta_raw - fe - Reduce(`+`, means)))
  eta_internal <- theta_raw %*% t(X)
  eta_public <- fe %*% t(X)
  max_predictor_identity <- 0
  for (k in seq_len(n_coef)) {
    internal_term <- matrix(theta_raw[, k], nrow = n_draw, ncol = nrow(X))
    public_term <- matrix(fe[, k], nrow = n_draw, ncol = nrow(X))
    for (b in seq_along(groups)) {
      internal_term <- internal_term + matrix(
        internal[[b]][, indices[[b]], k],
        nrow = n_draw, ncol = nrow(X)
      )
      public_term <- public_term + matrix(
        public[[b]][, indices[[b]], k],
        nrow = n_draw, ncol = nrow(X)
      )
    }
    max_predictor_identity <- max(
      max_predictor_identity, max(abs(internal_term - public_term))
    )
    eta_internal <- eta_internal + sweep(
      internal_term - matrix(
        theta_raw[, k], nrow = n_draw, ncol = nrow(X)
      ), 2L, X[, k], "*"
    )
    eta_public <- eta_public + sweep(
      public_term - matrix(fe[, k], nrow = n_draw, ncol = nrow(X)),
      2L, X[, k], "*"
    )
  }
  eta_brms <- posterior_linpred(fit, transform = FALSE)

  out <- data.frame(
    case = case,
    check = c(
      "each internal group-effect column sums to zero",
      "public r minus block mean equals internal r_s2z",
      "public b plus all block means equals theta_s2z",
      "public b + sum(r) equals internal theta + sum(r_s2z)",
      "posterior_linpred equals recovered public predictor",
      "posterior_linpred equals internal S2Z predictor"
    ),
    max_abs_error = c(
      max_sum_to_zero,
      max_centered_recovery,
      max_fixed_recovery,
      max_predictor_identity,
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
      grep(paste0("^r_", group, "\\["), fit_variables, value = TRUE)
    )
  }
  unique(c(out, intersect("sigma", fit_variables)))
}

ebfmi_by_chain <- function(nuts) {
  energy <- nuts[nuts$Parameter == "energy__", , drop = FALSE]
  vapply(split(energy, energy$Chain), function(x) {
    x <- x[order(x$Iteration), , drop = FALSE]
    variance <- stats::var(x$Value)
    if (nrow(x) < 3L || !is.finite(variance) || variance <= 0) {
      return(NA_real_)
    }
    mean(diff(x$Value)^2) / variance
  }, numeric(1))
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
    line <- grep(paste0("\\(", phase, "\\)"), x$time_info, value = TRUE)
    if (length(line) != 1L) {
      return(NA_real_)
    }
    fields <- strsplit(trimws(sub("^#", "", line)), "[[:space:]]+")[[1L]]
    suppressWarnings(as.numeric(fields[1L]))
  }, numeric(1))
  if (anyNA(times)) NA_real_ else mean(times)
}

fit_quality <- function(fit, groups, case, parameterization, elapsed) {
  summary <- posterior::summarise_draws(
    as_draws_array(fit, variable = public_variable_names(fit, groups)),
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  )
  nuts <- nuts_params(fit)
  leapfrog <- nuts$Value[nuts$Parameter == "n_leapfrog__"]
  gradient_evals <- sum(leapfrog + 1)
  divergences <- sum(nuts$Value[nuts$Parameter == "divergent__"])
  treedepth_hits <- sum(
    nuts$Value[nuts$Parameter == "treedepth__"] >= max_treedepth
  )
  ebfmi <- ebfmi_by_chain(nuts)
  finite_rhat <- summary$rhat[is.finite(summary$rhat)]
  finite_bulk <- summary$ess_bulk[is.finite(summary$ess_bulk)]
  finite_tail <- summary$ess_tail[is.finite(summary$ess_tail)]
  row <- data.frame(
    case = case,
    parameterization = parameterization,
    wall_seconds = elapsed,
    mean_chain_total_seconds = mean_cmdstan_chain_time(fit, "Total"),
    mean_chain_sampling_seconds = mean_cmdstan_chain_time(fit, "Sampling"),
    monitored_parameters = nrow(summary),
    postwarmup_transitions = length(leapfrog),
    gradient_evals = gradient_evals,
    mean_leapfrog_steps = mean(leapfrog),
    median_bulk_ess = stats::median(finite_bulk),
    min_bulk_ess = min(finite_bulk),
    median_tail_ess = stats::median(finite_tail),
    min_tail_ess = min(finite_tail),
    median_bulk_ess_per_1000_gradients =
      1000 * stats::median(finite_bulk) / gradient_evals,
    median_tail_ess_per_1000_gradients =
      1000 * stats::median(finite_tail) / gradient_evals,
    max_rhat = max(finite_rhat),
    divergences = divergences,
    treedepth_hits = treedepth_hits,
    min_ebfmi = min(ebfmi, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  row$pass <- row$max_rhat <= max_rhat &&
    row$divergences <= max_divergences &&
    row$treedepth_hits <= max_treedepth_hits &&
    row$min_ebfmi >= min_ebfmi
  row
}

scalar_data <- simulate_crossed_scalar()
scalar <- scalar_formulas(scalar_data)
scalar_case <- paste0(scalar_factors, "-crossed-scalar-factors")
vector_data <- simulate_overlapping_vector()

cat("Checking generated-code routing ...\n")
invisible(assert_matheron_code(
  scalar$s2z, scalar_data, scalar_prior, scalar = TRUE
))
invisible(assert_matheron_code(
  vector_formulas$s2z, vector_data, vector_prior, scalar = FALSE
))
cat("PASS: scalar and vector models select their intended Matheron paths.\n")

scalar_conventional <- sample_model(
  scalar$conventional, scalar_data, scalar_prior, 9201L,
  paste0(scalar_factors, "-crossed-factors-conventional")
)
scalar_s2z <- sample_model(
  scalar$s2z, scalar_data, scalar_prior, 9202L,
  paste0(scalar_factors, "-crossed-factors-s2z-scalar-Matheron")
)
vector_conventional <- sample_model(
  vector_formulas$conventional, vector_data, vector_prior, 9301L,
  "overlapping-vector-conventional"
)
vector_s2z <- sample_model(
  vector_formulas$s2z, vector_data, vector_prior, 9302L,
  "overlapping-vector-s2z-Matheron"
)

comparison_results <- rbind(
  compare_public_posteriors(
    public_samples(scalar_conventional$fit, scalar$groups),
    public_samples(scalar_s2z$fit, scalar$groups),
    scalar_case
  ),
  compare_public_posteriors(
    public_samples(vector_conventional$fit, c("g", "h")),
    public_samples(vector_s2z$fit, c("g", "h")),
    "overlapping-vector-blocks"
  )
)
invariant_results <- rbind(
  joint_coordinate_invariants(
    scalar_s2z$fit, scalar_data, ~ 1, scalar$groups,
    seq_along(scalar$groups), scalar_case
  ),
  joint_coordinate_invariants(
    vector_s2z$fit, vector_data, ~ x * z, c("g", "h"), 1:2,
    "overlapping-vector-blocks"
  )
)
quality_results <- rbind(
  fit_quality(
    scalar_conventional$fit, scalar$groups,
    scalar_case, "conventional",
    scalar_conventional$elapsed
  ),
  fit_quality(
    scalar_s2z$fit, scalar$groups,
    scalar_case, "s2z-Matheron",
    scalar_s2z$elapsed
  ),
  fit_quality(
    vector_conventional$fit, c("g", "h"),
    "overlapping-vector-blocks", "conventional",
    vector_conventional$elapsed
  ),
  fit_quality(
    vector_s2z$fit, c("g", "h"),
    "overlapping-vector-blocks", "s2z-Matheron",
    vector_s2z$elapsed
  )
)

paired_efficiency <- do.call(rbind, lapply(
  unique(quality_results$case), function(case_name) {
    conventional <- subset(
      quality_results, quality_results$case == case_name &
        parameterization == "conventional"
    )
    s2z <- subset(
      quality_results, quality_results$case == case_name &
        parameterization == "s2z-Matheron"
    )
    data.frame(
      case = case_name,
      wall_speedup_conventional_over_s2z =
        conventional$wall_seconds / s2z$wall_seconds,
      mean_chain_total_speedup_conventional_over_s2z =
        conventional$mean_chain_total_seconds / s2z$mean_chain_total_seconds,
      mean_chain_sampling_speedup_conventional_over_s2z =
        conventional$mean_chain_sampling_seconds /
        s2z$mean_chain_sampling_seconds,
      gradient_reduction_conventional_over_s2z =
        conventional$gradient_evals / s2z$gradient_evals,
      bulk_ess_per_gradient_gain_s2z_over_conventional =
        (s2z$median_bulk_ess / s2z$gradient_evals) /
        (conventional$median_bulk_ess / conventional$gradient_evals),
      tail_ess_per_gradient_gain_s2z_over_conventional =
        (s2z$median_tail_ess / s2z$gradient_evals) /
        (conventional$median_tail_ess / conventional$gradient_evals),
      stringsAsFactors = FALSE
    )
  }
))

comparison_summary <- do.call(rbind, lapply(
  split(comparison_results, comparison_results$case), function(x) {
    data.frame(
      case = x$case[1L],
      quantities = nrow(x) / 2L,
      mean_and_sd_checks = nrow(x),
      max_z_score = max(x$z_score),
      failed = sum(!x$pass),
      stringsAsFactors = FALSE
    )
  }
))

cat("\nPosterior equivalence summary\n")
print(comparison_summary, row.names = FALSE, digits = 5)
cat("\nLargest MCSE-standardized differences\n")
print(utils::head(
  comparison_results[order(comparison_results$z_score, decreasing = TRUE), ],
  20L
), row.names = FALSE, digits = 5)
cat("\nExact draw-wise identities\n")
print(invariant_results, row.names = FALSE, digits = 5)
cat("\nSampling diagnostics and efficiency\n")
print(quality_results, row.names = FALSE, digits = 5)
cat("\nPaired efficiency ratios (> 1 favors S2Z Matheron)\n")
print(paired_efficiency, row.names = FALSE, digits = 5)

failed_comparisons <- subset(comparison_results, !pass)
failed_invariants <- subset(invariant_results, !pass)
failed_quality <- subset(quality_results, !pass)
if (nrow(failed_comparisons) || nrow(failed_invariants) ||
    nrow(failed_quality)) {
  stop(
    "Gaussian Matheron validation failed: ",
    nrow(failed_comparisons), " MCSE comparisons, ",
    nrow(failed_invariants), " exact invariants, and ",
    nrow(failed_quality), " sampling-quality checks failed.",
    call. = FALSE
  )
}

cat(
  "\nPASS: the scalar and vector Matheron paths reproduce the conventional ",
  "public posterior within propagated MCSE, satisfy every draw-wise recovery ",
  "identity, and meet the configured sampler checks.\n",
  sep = ""
)
