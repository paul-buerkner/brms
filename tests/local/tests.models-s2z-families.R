# Sampling validation for physical sum-to-zero effects in non-Gaussian
# families and distributional parameters.
#
# This local test complements tests.models-s2z.R. It fits conventional and
# exactly equivalent S2Z versions of five small, deterministically simulated
# models: Bernoulli, Poisson, zero-inflated Poisson with S2Z in zi, hurdle
# Poisson with S2Z in hu, and a two-component Poisson mixture with S2Z in
# theta1. Run it against an installed development build of brms, for example
# with
#
#   R_LIBS=/tmp/brms-foundation-lib \
#     Rscript tests/local/tests.models-s2z-families.R
#
# Useful controls (shown with their defaults) are
#
#   BRMS_S2Z_CHAINS=2
#   BRMS_S2Z_WARMUP=400
#   BRMS_S2Z_SAMPLING=400
#   BRMS_S2Z_CORES=min(chains, parallel::detectCores())
#   BRMS_S2Z_BACKEND=auto
#   BRMS_S2Z_ADAPT_DELTA=0.97
#   BRMS_S2Z_MAX_TREEDEPTH=12
#   BRMS_S2Z_EQUIV_Z=7
#   BRMS_S2Z_MAX_RHAT=1.10
#   BRMS_S2Z_MAX_DIVERGENCES=0
#   BRMS_S2Z_MAX_TREEDEPTH_HITS=0
#   BRMS_S2Z_MIN_EBFMI=0.20
#   BRMS_S2Z_INVARIANT_TOL=1e-6
#   BRMS_S2Z_REPORT_ROWS=30
#   BRMS_S2Z_CACHE_DIR=
#
# Posterior equivalence is tested for means and standard deviations. For each
# independently sampled pair, the acceptance threshold is fixed in advance at
# BRMS_S2Z_EQUIV_Z times
#
#   sqrt(mcse_conventional^2 + mcse_s2z^2).

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
    stop("BRMS_S2Z_BACKEND must be auto, cmdstanr, or rstan.",
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

chains <- env_integer("BRMS_S2Z_CHAINS", 2L, minimum = 2L)
warmup <- env_integer("BRMS_S2Z_WARMUP", 400L)
sampling <- env_integer("BRMS_S2Z_SAMPLING", 400L)
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
max_treedepth_hits <- env_integer(
  "BRMS_S2Z_MAX_TREEDEPTH_HITS", 0L, minimum = 0L
)
min_ebfmi <- env_number("BRMS_S2Z_MIN_EBFMI", 0.20, minimum = 0)
invariant_tolerance <- env_number(
  # CmdStan's default CSV precision introduces errors around 1e-7 when
  # generated quantities are read back into R.
  "BRMS_S2Z_INVARIANT_TOL", 1e-6, minimum = 0
)
report_rows <- env_integer("BRMS_S2Z_REPORT_ROWS", 30L)
cache_dir <- Sys.getenv("BRMS_S2Z_CACHE_DIR", unset = "")

options(mc.cores = cores, width = 150)

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
cat("S2Z non-Gaussian sampling-equivalence configuration\n")
print(configuration, row.names = FALSE)

group_frame <- function(groups, per_group) {
  levels <- sprintf("g%02d", seq_len(groups))
  data.frame(
    g = factor(rep(levels, each = per_group), levels = levels),
    within = rep(seq(-1, 1, length.out = per_group), groups)
  )
}

simulate_bernoulli <- function(seed = 6201L, groups = 10L,
                               per_group = 16L) {
  set.seed(seed)
  dat <- group_frame(groups, per_group)
  # Keep the group scale away from the weakly identified boundary so this
  # equivalence test diagnoses the transformation rather than a generic
  # near-zero hierarchical funnel.
  u <- 1.4 * sin(seq(0, 2 * pi, length.out = groups + 1L)[-1L])
  eta <- -0.15 + u[as.integer(dat$g)]
  dat$y <- stats::rbinom(nrow(dat), size = 1L, prob = plogis(eta))
  dat[, c("y", "g")]
}

simulate_poisson <- function(seed = 6202L, groups = 10L,
                             per_group = 10L) {
  set.seed(seed)
  dat <- group_frame(groups, per_group)
  u <- seq(-0.55, 0.55, length.out = groups)
  lambda <- exp(0.55 + u[as.integer(dat$g)])
  dat$y <- stats::rpois(nrow(dat), lambda)
  dat[, c("y", "g")]
}

simulate_zero_process <- function(type = c("zi", "hu"), seed = 6203L,
                                  groups = 10L, per_group = 14L) {
  type <- match.arg(type)
  set.seed(seed)
  dat <- group_frame(groups, per_group)
  u_limit <- if (type == "zi") 1.5 else 1.05
  u <- seq(-u_limit, u_limit, length.out = groups)
  pzero <- plogis(-0.65 + u[as.integer(dat$g)])
  # A higher ZIP rate separates structural from sampling zeros and makes the
  # zi group scale a well-identified test quantity.
  lambda <- rep(if (type == "zi") 4 else 2.4, nrow(dat))
  is_zero <- stats::rbinom(nrow(dat), size = 1L, prob = pzero) == 1L
  if (type == "zi") {
    dat$y <- ifelse(is_zero, 0L, stats::rpois(nrow(dat), lambda))
  } else {
    p0 <- exp(-lambda)
    positive <- stats::qpois(
      p0 + stats::runif(nrow(dat)) * (1 - p0), lambda = lambda
    )
    dat$y <- ifelse(is_zero, 0L, positive)
  }
  dat[, c("y", "g")]
}

simulate_poisson_mixture <- function(seed = 6205L, groups = 10L,
                                     per_group = 16L) {
  set.seed(seed)
  dat <- group_frame(groups, per_group)
  u <- seq(-1.2, 1.2, length.out = groups)
  theta1 <- plogis(0.15 + u[as.integer(dat$g)])
  component1 <- stats::rbinom(nrow(dat), size = 1L, prob = theta1) == 1L
  lambda <- ifelse(component1, 1.25, 6.0)
  dat$y <- stats::rpois(nrow(dat), lambda)
  dat[, c("y", "g")]
}

cache_file <- function(label) {
  if (!nzchar(cache_dir)) {
    return(NULL)
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  suffix <- paste(backend, chains, warmup, sampling, sep = "-")
  file.path(cache_dir, paste0(label, "-", suffix, ".rds"))
}

sample_model <- function(case, parameterization) {
  label <- paste(case$name, parameterization, sep = "-")
  cat("Sampling ", label, " ...\n", sep = "")
  formula <- case[[paste0(parameterization, "_formula")]]
  args <- list(
    formula = formula,
    data = case$data,
    family = case$family,
    prior = case$prior,
    backend = backend,
    chains = chains,
    cores = cores,
    iter = warmup + sampling,
    warmup = warmup,
    seed = case$seed + identical(parameterization, "s2z"),
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

selected_observations <- function(n) {
  unique(round(seq(1, n, length.out = min(5L, n))))
}

append_matrix_samples <- function(out, value, prefix) {
  selected <- selected_observations(ncol(value))
  for (i in selected) {
    out[[sprintf("%s:observation[%d]", prefix, i)]] <- value[, i]
  }
  out[[paste0(prefix, ":observation-average")]] <- rowMeans(value)
  out
}

public_samples <- function(fit, case) {
  fe <- as.matrix(fixef(fit, summary = FALSE))
  if (any(grepl("s2z", colnames(fe), fixed = TRUE))) {
    stop("An internal S2Z coefficient leaked through fixef().", call. = FALSE)
  }
  out <- setNames(
    lapply(seq_len(ncol(fe)), function(k) fe[, k]),
    paste0("b:", colnames(fe))
  )

  re <- ranef(fit, summary = FALSE)$g
  if (is.null(re) || !case$re_coef %in% dimnames(re)[[3L]]) {
    stop("No public group-level draws found for ", case$re_coef,
         call. = FALSE)
  }
  for (j in seq_len(dim(re)[2L])) {
    label <- paste(dimnames(re)[[2L]][j], case$re_coef, sep = ",")
    out[[paste0("r:[", label, "]")]] <- re[, j, case$re_coef]
  }

  hyper_name <- paste0("sd_g__", case$re_coef)
  hyper <- exact_draw_matrix(fit, hyper_name)
  out[[paste0("parameter:", hyper_name)]] <- hyper[, 1L]

  for (dpar in case$prediction_dpars) {
    eta <- posterior_linpred(fit, dpar = dpar, transform = FALSE)
    out <- append_matrix_samples(out, eta, paste0("eta:", dpar))
    response_dpar <- posterior_epred(fit, dpar = dpar)
    out <- append_matrix_samples(
      out, response_dpar, paste0("dpar-response:", dpar)
    )
  }
  epred <- posterior_epred(fit)
  append_matrix_samples(out, epred, "epred")
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
      threshold <- equiv_z * mcse_difference
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
        propagated_mcse = mcse_difference,
        threshold = threshold,
        z_score = z_score,
        pass = abs(difference) <= threshold,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  out[order(out$z_score, decreasing = TRUE), ]
}

inverse_link <- function(eta, dpar) {
  if (dpar %in% c("mu_bernoulli", "zi", "hu", "theta1")) {
    return(plogis(eta))
  }
  if (dpar %in% c("mu_poisson", "mu", "mu1", "mu2")) {
    return(exp(eta))
  }
  stop("No inverse-link identity registered for ", dpar, call. = FALSE)
}

response_mean_identity <- function(fit, case) {
  epred <- posterior_epred(fit)
  expected <- switch(
    case$mean_identity,
    bernoulli = posterior_epred(fit, dpar = "mu"),
    poisson = posterior_epred(fit, dpar = "mu"),
    zero_inflated_poisson = {
      mu <- posterior_epred(fit, dpar = "mu")
      zi <- posterior_epred(fit, dpar = "zi")
      (1 - zi) * mu
    },
    hurdle_poisson = {
      mu <- posterior_epred(fit, dpar = "mu")
      hu <- posterior_epred(fit, dpar = "hu")
      (1 - hu) * mu / (-expm1(-mu))
    },
    poisson_mixture = {
      mu1 <- posterior_epred(fit, dpar = "mu1")
      mu2 <- posterior_epred(fit, dpar = "mu2")
      theta1 <- posterior_epred(fit, dpar = "theta1")
      theta1 * mu1 + (1 - theta1) * mu2
    },
    stop("Unknown response-mean identity.", call. = FALSE)
  )
  max(abs(epred - expected))
}

public_prediction_invariants <- function(fit, case, parameterization) {
  fe <- as.matrix(fixef(fit, summary = FALSE))
  re <- ranef(fit, summary = FALSE)$g
  levels <- levels(case$data$g)
  group_index <- as.integer(case$data$g)
  eta <- posterior_linpred(
    fit, dpar = case$active_dpar, transform = FALSE
  )
  eta_public <- sweep(
    matrix(re[, levels, case$re_coef], nrow = nrow(fe)),
    1L, fe[, case$fe_coef], "+"
  )[, group_index, drop = FALSE]
  transformed <- posterior_epred(fit, dpar = case$active_dpar)
  dpar_link_name <- if (case$name == "bernoulli") {
    "mu_bernoulli"
  } else {
    case$active_dpar
  }
  expected_transformed <- inverse_link(eta, dpar_link_name)

  public_r_names <- if (identical(case$active_dpar, "mu")) {
    sprintf("r_g[%s,%s]", levels, case$re_coef)
  } else {
    sprintf("r_g__%s[%s,Intercept]", case$active_dpar, levels)
  }
  expected_public <- c(
    paste0("b_", case$fe_coef), public_r_names,
    paste0("sd_g__", case$re_coef)
  )
  public_names_ok <- all(expected_public %in% variables(fit))
  no_internal_fixef <- !any(grepl("s2z", colnames(fe), fixed = TRUE))
  rows <- data.frame(
    case = case$name,
    parameterization = parameterization,
    check = c(
      "saved public b/r/sd names use conventional API",
      "fixef excludes internal S2Z coordinates",
      "active dpar predictor equals public b+r",
      "active dpar inverse link is exact",
      "response-scale posterior mean identity is exact"
    ),
    max_abs_error = c(
      if (public_names_ok) 0 else Inf,
      if (no_internal_fixef) 0 else Inf,
      max(abs(eta - eta_public)),
      max(abs(transformed - expected_transformed)),
      response_mean_identity(fit, case)
    ),
    tolerance = invariant_tolerance,
    stringsAsFactors = FALSE
  )
  rows$pass <- rows$max_abs_error <= rows$tolerance
  rows
}

s2z_internal_invariants <- function(fit, case) {
  n_group <- nlevels(case$data$g)
  theta <- exact_draw_matrix(fit, case$internal_theta)[, 1L]
  internal_names <- sprintf(case$internal_re, seq_len(n_group))
  internal_re <- exact_draw_matrix(fit, internal_names)
  fe <- as.matrix(fixef(fit, summary = FALSE))[, case$fe_coef]
  re <- ranef(fit, summary = FALSE)$g[, levels(case$data$g), case$re_coef]
  public_finite <- sweep(matrix(re, nrow = length(fe)), 1L, fe, "+")
  internal_finite <- sweep(internal_re, 1L, theta, "+")
  rows <- data.frame(
    case = case$name,
    parameterization = "s2z",
    check = c(
      "internal physical group effects sum to zero",
      "public b+r equals internal finite coefficient"
    ),
    max_abs_error = c(
      max(abs(rowSums(internal_re))),
      max(abs(public_finite - internal_finite))
    ),
    tolerance = invariant_tolerance,
    stringsAsFactors = FALSE
  )
  rows$pass <- rows$max_abs_error <= rows$tolerance
  rows
}

chain_ebfmi <- function(nuts) {
  energy <- nuts[nuts$Parameter == "energy__", c("Chain", "Value")]
  split_energy <- split(energy$Value, energy$Chain)
  vapply(split_energy, function(x) {
    if (length(x) < 3L || !is.finite(stats::var(x)) || stats::var(x) == 0) {
      return(NA_real_)
    }
    mean(diff(x)^2) / stats::var(x)
  }, numeric(1L))
}

fit_quality <- function(fit, case, parameterization) {
  public <- grep(
    "^(b_|sd_g__|cor_g__|r_g\\[)", variables(fit), value = TRUE
  )
  draws <- as_draws_array(fit, variable = public)
  summary <- posterior::summarise_draws(
    draws,
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  )
  nuts <- nuts_params(fit)
  divergences <- sum(nuts$Value[nuts$Parameter == "divergent__"])
  treedepth_hits <- sum(
    nuts$Value[nuts$Parameter == "treedepth__"] >= max_treedepth
  )
  ebfmi <- chain_ebfmi(nuts)
  observed_max_rhat <- max(summary$rhat, na.rm = TRUE)
  observed_min_ebfmi <- min(ebfmi, na.rm = TRUE)
  data.frame(
    case = case,
    parameterization = parameterization,
    max_rhat = observed_max_rhat,
    min_bulk_ess = min(summary$ess_bulk, na.rm = TRUE),
    min_tail_ess = min(summary$ess_tail, na.rm = TRUE),
    divergences = divergences,
    treedepth_hits = treedepth_hits,
    min_ebfmi = observed_min_ebfmi,
    pass = observed_max_rhat <= max_rhat &&
      divergences <= max_divergences &&
      treedepth_hits <= max_treedepth_hits &&
      observed_min_ebfmi >= min_ebfmi,
    stringsAsFactors = FALSE
  )
}

representative_api_checks <- function(fit, case) {
  fe <- as.matrix(fixef(fit, summary = FALSE))
  re <- ranef(fit, summary = FALSE)$g
  co <- coef(fit, summary = FALSE)$g
  expected_coef <- sweep(
    matrix(re[, , case$re_coef], nrow = nrow(fe)),
    1L, fe[, case$fe_coef], "+"
  )
  coef_error <- max(abs(co[, , case$re_coef] - expected_coef))

  vc <- VarCorr(fit, summary = FALSE)$g
  sd_draws <- exact_draw_matrix(
    fit, paste0("sd_g__", case$re_coef)
  )[, 1L]
  varcorr_error <- max(abs(vc$sd[, case$re_coef] - sd_draws))

  observed <- posterior_epred(fit)
  first_level <- levels(case$data$g)[1L]
  newdata <- data.frame(g = c(first_level, "new-group"))
  set.seed(7211L)
  new_epred <- posterior_epred(
    fit, newdata = newdata, allow_new_levels = TRUE,
    sample_new_levels = "gaussian"
  )
  known_prediction_error <- max(abs(new_epred[, 1L] - observed[, 1L]))
  new_prediction_valid <- all(is.finite(new_epred)) &&
    all(new_epred >= 0) && all(new_epred <= 1)

  serialization_file <- tempfile(fileext = ".rds")
  saveRDS(fit, serialization_file)
  restored <- readRDS(serialization_file)
  unlink(serialization_file)
  serialization_error <- max(
    abs(as.matrix(fixef(restored, summary = FALSE)) - fe),
    abs(ranef(restored, summary = FALSE)$g - re),
    abs(posterior_epred(restored) - observed)
  )

  updated <- update(
    fit, newdata = case$data, recompile = FALSE,
    testmode = TRUE, silent = 2
  )
  update_ok <- inherits(updated, "brmsfit") &&
    isTRUE(all.equal(updated$formula, fit$formula)) &&
    identical(nrow(updated$data), nrow(case$data)) &&
    all(case$data$g == updated$data$g)

  rows <- data.frame(
    case = case$name,
    check = c(
      "fixef exposes the expected public coefficient",
      "ranef exposes the expected public group effect",
      "coef equals fixef plus ranef",
      "VarCorr sd equals the saved public sd",
      "observed-level newdata prediction is unchanged",
      "new-level predictions are finite probabilities",
      "serialized fit preserves public draws and predictions",
      "update validation path preserves S2Z formula and data"
    ),
    max_abs_error = c(
      if (case$fe_coef %in% colnames(fe)) 0 else Inf,
      if (case$re_coef %in% dimnames(re)[[3L]]) 0 else Inf,
      coef_error,
      varcorr_error,
      known_prediction_error,
      if (new_prediction_valid) 0 else Inf,
      serialization_error,
      if (update_ok) 0 else Inf
    ),
    tolerance = invariant_tolerance,
    stringsAsFactors = FALSE
  )
  rows$pass <- rows$max_abs_error <= rows$tolerance
  rows
}

run_case <- function(case) {
  conventional <- sample_model(case, "conventional")
  s2z <- sample_model(case, "s2z")
  comparison <- compare_public_posteriors(
    public_samples(conventional, case),
    public_samples(s2z, case),
    case = case$name
  )
  invariants <- rbind(
    public_prediction_invariants(conventional, case, "conventional"),
    public_prediction_invariants(s2z, case, "s2z"),
    s2z_internal_invariants(s2z, case)
  )
  quality <- rbind(
    fit_quality(conventional, case$name, "conventional"),
    fit_quality(s2z, case$name, "s2z")
  )
  list(
    conventional = conventional,
    s2z = s2z,
    comparison = comparison,
    invariants = invariants,
    quality = quality
  )
}

cases <- list(
  list(
    name = "bernoulli",
    data = simulate_bernoulli(),
    family = bernoulli(),
    conventional_formula = bf(y ~ 1 + (1 | gr(g, id = "binary"))),
    s2z_formula = bf(
      y ~ 1 + (1 | gr(g, id = "binary", s2z = TRUE))
    ),
    prior = prior(logistic(0, 1), class = "Intercept") +
      prior(exponential(1.5), class = "sd", group = "g"),
    seed = 7201L,
    active_dpar = "mu",
    prediction_dpars = "mu",
    fe_coef = "Intercept",
    re_coef = "Intercept",
    internal_theta = "theta_s2z[1]",
    internal_re = "r_s2z_1_1[%d]",
    mean_identity = "bernoulli"
  ),
  list(
    name = "poisson",
    data = simulate_poisson(),
    family = poisson(),
    conventional_formula = bf(y ~ 1 + (1 | gr(g, id = "count"))),
    s2z_formula = bf(
      y ~ 1 + (1 | gr(g, id = "count", s2z = TRUE))
    ),
    prior = prior(normal(0, 1.2), class = "Intercept") +
      prior(exponential(1.5), class = "sd", group = "g"),
    seed = 7202L,
    active_dpar = "mu",
    prediction_dpars = "mu",
    fe_coef = "Intercept",
    re_coef = "Intercept",
    internal_theta = "theta_s2z[1]",
    internal_re = "r_s2z_1_1[%d]",
    mean_identity = "poisson"
  ),
  list(
    name = "zero-inflated-poisson-zi",
    data = simulate_zero_process(
      "zi", seed = 6203L, per_group = 20L
    ),
    family = zero_inflated_poisson(),
    conventional_formula = bf(
      y ~ 1, zi ~ 1 + (1 | gr(g, id = "zero"))
    ),
    s2z_formula = bf(
      y ~ 1,
      zi ~ 1 + (1 | gr(g, id = "zero", s2z = TRUE))
    ),
    prior = prior(normal(log(4), 0.6), class = "Intercept") +
      prior(logistic(0, 1), class = "Intercept", dpar = "zi") +
      prior(
        exponential(1.5), class = "sd", group = "g", dpar = "zi"
      ),
    seed = 7203L,
    active_dpar = "zi",
    prediction_dpars = c("mu", "zi"),
    fe_coef = "zi_Intercept",
    re_coef = "zi_Intercept",
    internal_theta = "theta_s2z_zi[1]",
    internal_re = "r_s2z_1_zi_1[%d]",
    mean_identity = "zero_inflated_poisson"
  ),
  list(
    name = "hurdle-poisson-hu",
    data = simulate_zero_process("hu", seed = 6204L),
    family = hurdle_poisson(),
    conventional_formula = bf(
      y ~ 1, hu ~ 1 + (1 | gr(g, id = "hurdle"))
    ),
    s2z_formula = bf(
      y ~ 1,
      hu ~ 1 + (1 | gr(g, id = "hurdle", s2z = TRUE))
    ),
    prior = prior(normal(log(2.4), 0.6), class = "Intercept") +
      prior(logistic(0, 1), class = "Intercept", dpar = "hu") +
      prior(
        exponential(1.5), class = "sd", group = "g", dpar = "hu"
      ),
    seed = 7204L,
    active_dpar = "hu",
    prediction_dpars = c("mu", "hu"),
    fe_coef = "hu_Intercept",
    re_coef = "hu_Intercept",
    internal_theta = "theta_s2z_hu[1]",
    internal_re = "r_s2z_1_hu_1[%d]",
    mean_identity = "hurdle_poisson"
  ),
  list(
    name = "poisson-mixture-theta1",
    data = simulate_poisson_mixture(),
    # Ordering the component locations makes this a single-mode validation
    # fixture rather than a test of mixture label switching.
    family = mixture(poisson, poisson, order = "mu"),
    conventional_formula = bf(
      y ~ 1, mu2 ~ 1,
      theta1 ~ 1 + (1 | gr(g, id = "mixture"))
    ),
    s2z_formula = bf(
      y ~ 1, mu2 ~ 1,
      theta1 ~ 1 + (1 | gr(g, id = "mixture", s2z = TRUE))
    ),
    prior = prior(
      normal(log(1.25), 0.35), class = "Intercept", dpar = "mu1"
    ) +
      prior(
        normal(log(6), 0.35), class = "Intercept", dpar = "mu2"
      ) +
      prior(logistic(0, 1), class = "Intercept", dpar = "theta1") +
      prior(
        exponential(1.5), class = "sd", group = "g", dpar = "theta1"
      ),
    seed = 7205L,
    active_dpar = "theta1",
    prediction_dpars = c("mu1", "mu2", "theta1"),
    fe_coef = "theta1_Intercept",
    re_coef = "theta1_Intercept",
    internal_theta = "theta_s2z_theta1[1]",
    internal_re = "r_s2z_1_theta1_1[%d]",
    mean_identity = "poisson_mixture"
  )
)

results <- lapply(cases, run_case)
comparison_results <- do.call(rbind, lapply(results, `[[`, "comparison"))
invariant_results <- do.call(rbind, lapply(results, `[[`, "invariants"))
quality_results <- do.call(rbind, lapply(results, `[[`, "quality"))
api_results <- representative_api_checks(results[[1L]]$s2z, cases[[1L]])

cat("\nMCSE-aware conventional-vs-S2Z comparisons (largest z-scores first)\n")
comparison_results <- comparison_results[
  order(comparison_results$z_score, decreasing = TRUE),
]
print(
  utils::head(comparison_results, report_rows),
  row.names = FALSE, digits = 4
)
if (nrow(comparison_results) > report_rows) {
  cat("... ", nrow(comparison_results) - report_rows,
      " additional passing-or-failing comparisons retained for assertions.\n",
      sep = "")
}
cat("\nExact public prediction and S2Z coordinate invariants\n")
print(invariant_results, row.names = FALSE, digits = 4)
cat("\nSampling quality\n")
print(quality_results, row.names = FALSE, digits = 4)
cat("\nRepresentative public API checks\n")
print(api_results, row.names = FALSE, digits = 4)

failed_comparisons <- subset(comparison_results, !pass)
failed_invariants <- subset(invariant_results, !pass)
failed_quality <- subset(quality_results, !pass)
failed_api <- subset(api_results, !pass)
if (nrow(failed_comparisons) || nrow(failed_invariants) ||
    nrow(failed_quality) || nrow(failed_api)) {
  stop(
    "S2Z non-Gaussian sampling validation failed: ",
    nrow(failed_comparisons), " MCSE comparisons, ",
    nrow(failed_invariants), " exact invariants, ",
    nrow(failed_quality), " sampling-quality checks, and ",
    nrow(failed_api), " public API checks failed.",
    call. = FALSE
  )
}

cat(
  "\nPASS: all five conventional and S2Z posteriors agree within the ",
  "prespecified propagated-MCSE thresholds; dpar/public-coordinate ",
  "identities, sampler diagnostics, and representative public API ",
  "checks all pass.\n",
  sep = ""
)
