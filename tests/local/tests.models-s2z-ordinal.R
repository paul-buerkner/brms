# PLAN-02 sampling validation for physical sum-to-zero effects in ordinal
# location predictors.
#
# This is deliberately a local-only harness: it compiles and samples real
# Stan models and is not sourced by the package test suite. Run it against an
# installed development build of brms, for example with
#
#   R_LIBS=/tmp/brms-plan02-lib \
#     Rscript tests/local/tests.models-s2z-ordinal.R
#
# The default run fits conventional and exactly equivalent physical-S2Z
# versions of cumulative, cratio, sratio, acat, and hurdle_cumulative models.
# The primary models have varying intercepts and slopes; additional cumulative
# cases exercise slope-only, intercept-only, and Student group effects. Both
# weak and coefficient-specific informative threshold priors are represented.
# Useful controls (shown with their defaults) are
#
#   BRMS_S2Z_ORDINAL_CHAINS=4
#   BRMS_S2Z_ORDINAL_WARMUP=750
#   BRMS_S2Z_ORDINAL_SAMPLING=750
#   BRMS_S2Z_ORDINAL_CORES=min(chains, parallel::detectCores())
#   BRMS_S2Z_ORDINAL_BACKEND=auto
#   BRMS_S2Z_ORDINAL_CACHE_DIR=
#   BRMS_S2Z_ORDINAL_FILE_REFIT=on_change
#   BRMS_S2Z_ORDINAL_CASES=all
#   BRMS_S2Z_ORDINAL_STUDENT=true
#   BRMS_S2Z_ORDINAL_ADAPT_DELTA=0.97
#   BRMS_S2Z_ORDINAL_MAX_TREEDEPTH=12
#   BRMS_S2Z_ORDINAL_EQUIV_Z=7
#   BRMS_S2Z_ORDINAL_MAX_RHAT=1.01
#   BRMS_S2Z_ORDINAL_MIN_BULK_ESS=400
#   BRMS_S2Z_ORDINAL_MIN_TAIL_ESS=400
#   BRMS_S2Z_ORDINAL_MAX_DIVERGENCES=0
#   BRMS_S2Z_ORDINAL_MAX_TREEDEPTH_HITS=0
#   BRMS_S2Z_ORDINAL_MIN_EBFMI=0.20
#   BRMS_S2Z_ORDINAL_INVARIANT_TOL=1e-6
#   BRMS_S2Z_ORDINAL_REPORT_ROWS=40
#
# BRMS_S2Z_ORDINAL_CASES may be a comma-separated subset of the case names
# printed by this script. Posterior equivalence is assessed for both means and
# standard deviations. Each independently sampled difference must be no more
# than BRMS_S2Z_ORDINAL_EQUIV_Z times its propagated Monte Carlo standard
# error. Exact reconstruction identities have a separate numerical tolerance.

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

env_logical <- function(name, default) {
  value <- tolower(Sys.getenv(name, unset = ""))
  if (!nzchar(value)) {
    return(isTRUE(default))
  }
  if (value %in% c("true", "t", "yes", "y", "1")) {
    return(TRUE)
  }
  if (value %in% c("false", "f", "no", "n", "0")) {
    return(FALSE)
  }
  stop(name, " must be true or false.", call. = FALSE)
}

choose_backend <- function() {
  requested <- tolower(Sys.getenv(
    "BRMS_S2Z_ORDINAL_BACKEND", unset = "auto"
  ))
  if (requested %in% c("cmdstanr", "rstan")) {
    return(requested)
  }
  if (!identical(requested, "auto")) {
    stop(
      "BRMS_S2Z_ORDINAL_BACKEND must be auto, cmdstanr, or rstan.",
      call. = FALSE
    )
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
  stop(
    "Neither a working CmdStan installation nor rstan was found.",
    call. = FALSE
  )
}

chains <- env_integer("BRMS_S2Z_ORDINAL_CHAINS", 4L, minimum = 2L)
warmup <- env_integer("BRMS_S2Z_ORDINAL_WARMUP", 750L)
sampling <- env_integer("BRMS_S2Z_ORDINAL_SAMPLING", 750L)
detected_cores <- parallel::detectCores(logical = FALSE)
if (is.na(detected_cores)) {
  detected_cores <- 1L
}
cores <- env_integer(
  "BRMS_S2Z_ORDINAL_CORES", min(chains, detected_cores), minimum = 1L
)
backend <- choose_backend()
adapt_delta <- env_number(
  "BRMS_S2Z_ORDINAL_ADAPT_DELTA", 0.97, minimum = 0.5
)
max_treedepth <- env_integer("BRMS_S2Z_ORDINAL_MAX_TREEDEPTH", 12L)
equiv_z <- env_number("BRMS_S2Z_ORDINAL_EQUIV_Z", 7, minimum = 1)
max_rhat <- env_number("BRMS_S2Z_ORDINAL_MAX_RHAT", 1.01, minimum = 1)
min_bulk_ess <- env_number(
  "BRMS_S2Z_ORDINAL_MIN_BULK_ESS", 400, minimum = 0
)
min_tail_ess <- env_number(
  "BRMS_S2Z_ORDINAL_MIN_TAIL_ESS", 400, minimum = 0
)
max_divergences <- env_integer(
  "BRMS_S2Z_ORDINAL_MAX_DIVERGENCES", 0L, minimum = 0L
)
max_treedepth_hits <- env_integer(
  "BRMS_S2Z_ORDINAL_MAX_TREEDEPTH_HITS", 0L, minimum = 0L
)
min_ebfmi <- env_number(
  "BRMS_S2Z_ORDINAL_MIN_EBFMI", 0.20, minimum = 0
)
invariant_tolerance <- env_number(
  # CmdStan's default CSV precision can leave reconstruction errors near 1e-7.
  "BRMS_S2Z_ORDINAL_INVARIANT_TOL", 1e-6, minimum = 0
)
report_rows <- env_integer("BRMS_S2Z_ORDINAL_REPORT_ROWS", 40L)
include_student <- env_logical("BRMS_S2Z_ORDINAL_STUDENT", TRUE)
cache_dir <- Sys.getenv("BRMS_S2Z_ORDINAL_CACHE_DIR", unset = "")
file_refit <- Sys.getenv(
  "BRMS_S2Z_ORDINAL_FILE_REFIT", unset = "on_change"
)
if (!file_refit %in% c("never", "on_change", "always")) {
  stop(
    "BRMS_S2Z_ORDINAL_FILE_REFIT must be never, on_change, or always.",
    call. = FALSE
  )
}

options(mc.cores = cores, width = 160)

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
  max_rhat = max_rhat,
  min_bulk_ess = min_bulk_ess,
  min_tail_ess = min_tail_ess,
  stringsAsFactors = FALSE
)
cat("PLAN-02 ordinal S2Z sampling-equivalence configuration\n")
print(configuration, row.names = FALSE)

simulate_ordinal <- function(seed = 8201L, groups = 12L,
                             per_group = 14L) {
  set.seed(seed)
  n <- groups * per_group
  level_names <- sprintf("g%02d", seq_len(groups))
  g <- factor(rep(level_names, each = per_group), levels = level_names)
  within <- rep(seq(-1.2, 1.2, length.out = per_group), groups)
  x <- 0.65 + 0.80 * within + stats::rnorm(n, sd = 0.35)
  z <- stats::rnorm(n)

  raw0 <- stats::rnorm(groups)
  raw1 <- 0.30 * raw0 + sqrt(1 - 0.30^2) * stats::rnorm(groups)
  u0 <- 0.70 * (raw0 - mean(raw0)) / stats::sd(raw0)
  u1 <- 0.32 * (raw1 - mean(raw1)) / stats::sd(raw1)
  eta <- 0.65 * x + u0[as.integer(g)] + u1[as.integer(g)] * x
  latent <- eta + stats::rlogis(n)
  threshold <- c(-1.10, 0.15, 1.25)
  y <- ordered(
    cut(latent, breaks = c(-Inf, threshold, Inf), labels = FALSE),
    levels = seq_len(length(threshold) + 1L)
  )
  dat <- data.frame(y = y, x = x, z = z, g = g)
  if (any(table(dat$y) < 4L)) {
    stop("Ordinal simulation did not populate every category.", call. = FALSE)
  }
  dat
}

simulate_hurdle_ordinal <- function(seed = 8202L, groups = 12L,
                                    per_group = 14L) {
  dat <- simulate_ordinal(seed = seed, groups = groups, per_group = per_group)
  set.seed(seed + 1L)
  hurdle <- stats::rbinom(
    nrow(dat), size = 1L, prob = plogis(-1.15 + 0.30 * dat$z)
  ) == 1L
  dat$y <- ifelse(hurdle, 0L, as.integer(dat$y))
  if (!all(0:4 %in% dat$y)) {
    stop("Hurdle-ordinal simulation did not populate every category.",
         call. = FALSE)
  }
  dat
}

threshold_prior <- function(type = c("weak", "informative")) {
  type <- match.arg(type)
  if (type == "weak") {
    return(prior(normal(0, 2.5), class = Intercept))
  }
  # Coefficient-specific normal, Student-t, and Cauchy priors exercise the
  # conditional-Gaussian ordinal fallback as well as unequal prior locations.
  prior(normal(-1.10, 0.45), class = Intercept, coef = "1") +
    prior(student_t(6, 0.15, 0.50), class = Intercept, coef = "2") +
    prior(cauchy(1.25, 0.40), class = Intercept, coef = "3")
}

model_prior <- function(threshold = c("weak", "informative"),
                        correlated = TRUE, student_group = FALSE,
                        hurdle = FALSE) {
  threshold <- match.arg(threshold)
  out <- threshold_prior(threshold) +
    prior(normal(0, 1), class = b) +
    prior(exponential(1.25), class = sd, group = g)
  if (correlated) {
    out <- out + prior(lkj(2), class = cor, group = g)
  }
  if (student_group) {
    out <- out + prior(gamma(2, 0.2), class = df, group = g)
  }
  if (hurdle) {
    out <- out +
      prior(normal(-1, 0.8), class = Intercept, dpar = hu) +
      prior(normal(0, 0.8), class = b, dpar = hu)
  }
  out
}

ordinary_data <- simulate_ordinal()
hurdle_data <- simulate_hurdle_ordinal()

cases <- list(
  list(
    name = "cumulative-intercept-slope-weak",
    data = ordinary_data,
    family = cumulative(),
    conventional_formula = bf(
      y ~ x + (1 + x | gr(g, id = "cumulative"))
    ),
    s2z_formula = bf(
      y ~ x + (1 + x | gr(g, id = "cumulative", s2z = TRUE))
    ),
    prior = model_prior("weak"),
    group_coef = c("Intercept", "x"),
    response_values = 1:4,
    seed = 9201L
  ),
  list(
    name = "cratio-intercept-slope-informative",
    data = ordinary_data,
    family = cratio(),
    conventional_formula = bf(
      y ~ x + (1 + x | gr(g, id = "cratio"))
    ),
    s2z_formula = bf(
      y ~ x + (1 + x | gr(g, id = "cratio", s2z = TRUE))
    ),
    prior = model_prior("informative"),
    group_coef = c("Intercept", "x"),
    response_values = 1:4,
    seed = 9202L
  ),
  list(
    name = "sratio-intercept-slope-weak",
    data = ordinary_data,
    family = sratio(),
    conventional_formula = bf(
      y ~ x + (1 + x | gr(g, id = "sratio"))
    ),
    s2z_formula = bf(
      y ~ x + (1 + x | gr(g, id = "sratio", s2z = TRUE))
    ),
    prior = model_prior("weak"),
    group_coef = c("Intercept", "x"),
    response_values = 1:4,
    seed = 9203L
  ),
  list(
    name = "acat-intercept-slope-informative",
    data = ordinary_data,
    family = acat(),
    conventional_formula = bf(
      y ~ x + (1 + x | gr(g, id = "acat"))
    ),
    s2z_formula = bf(
      y ~ x + (1 + x | gr(g, id = "acat", s2z = TRUE))
    ),
    prior = model_prior("informative"),
    group_coef = c("Intercept", "x"),
    response_values = 1:4,
    seed = 9204L
  ),
  list(
    name = "hurdle-cumulative-intercept-slope-weak",
    data = hurdle_data,
    family = hurdle_cumulative(),
    conventional_formula = bf(
      y ~ x + (1 + x | gr(g, id = "hurdle-ordinal")),
      hu ~ 1 + z
    ),
    s2z_formula = bf(
      y ~ x + (1 + x | gr(
        g, id = "hurdle-ordinal", s2z = TRUE
      )),
      hu ~ 1 + z
    ),
    prior = model_prior("weak", hurdle = TRUE),
    group_coef = c("Intercept", "x"),
    response_values = 0:4,
    seed = 9205L
  ),
  list(
    name = "cumulative-slope-only-informative",
    data = ordinary_data,
    family = cumulative(),
    conventional_formula = bf(
      y ~ x + (0 + x | gr(g, id = "slope-only"))
    ),
    s2z_formula = bf(
      y ~ x + (0 + x | gr(g, id = "slope-only", s2z = TRUE))
    ),
    prior = model_prior("informative", correlated = FALSE),
    group_coef = "x",
    response_values = 1:4,
    seed = 9206L
  ),
  list(
    name = "cumulative-intercept-only-weak",
    data = ordinary_data,
    family = cumulative(),
    conventional_formula = bf(
      y ~ x + (1 | gr(g, id = "intercept-only"))
    ),
    s2z_formula = bf(
      y ~ x + (1 | gr(g, id = "intercept-only", s2z = TRUE))
    ),
    prior = model_prior("weak", correlated = FALSE),
    group_coef = "Intercept",
    response_values = 1:4,
    seed = 9207L
  ),
  list(
    name = "cumulative-student-intercept-slope-informative",
    data = ordinary_data,
    family = cumulative(),
    conventional_formula = bf(
      y ~ x + (1 + x | gr(g, id = "student", dist = "student"))
    ),
    s2z_formula = bf(
      y ~ x + (1 + x | gr(
        g, id = "student", dist = "student", s2z = TRUE
      ))
    ),
    prior = model_prior("informative", student_group = TRUE),
    group_coef = c("Intercept", "x"),
    response_values = 1:4,
    seed = 9208L,
    student = TRUE
  )
)

if (!include_student) {
  cases <- Filter(function(x) !isTRUE(x$student), cases)
}
requested_cases <- trimws(strsplit(Sys.getenv(
  "BRMS_S2Z_ORDINAL_CASES", unset = "all"
), ",", fixed = TRUE)[[1L]])
if (!identical(requested_cases, "all")) {
  unknown <- setdiff(requested_cases, vapply(cases, `[[`, character(1), "name"))
  if (length(unknown)) {
    stop("Unknown ordinal S2Z case(s): ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  cases <- Filter(function(x) x$name %in% requested_cases, cases)
}
if (!length(cases)) {
  stop("No ordinal S2Z sampling cases were selected.", call. = FALSE)
}
cat("Selected cases\n")
writeLines(paste0("  - ", vapply(cases, `[[`, character(1), "name")))

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
    args$file_refit <- file_refit
  }
  do.call(brm, args)
}

exact_draw_matrix <- function(fit, variable) {
  missing <- setdiff(variable, variables(fit))
  if (length(missing)) {
    stop(
      "Required saved variables are missing: ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  out <- as_draws_matrix(fit, variable = variable)
  # `[.draws_matrix` deliberately preserves a singleton variable dimension.
  # Strip the posterior subclass here so ordinary matrix subsetting below has
  # the base-R drop semantics used by the reconstruction algebra.
  matrix(
    as.numeric(out), nrow = nrow(out), ncol = ncol(out),
    dimnames = dimnames(out)
  )
}

selected_indices <- function(n, maximum = 5L) {
  unique(round(seq(1, n, length.out = min(maximum, n))))
}

threshold_columns <- function(fe) {
  out <- grep("^Intercept\\[[0-9]+\\]$", colnames(fe), value = TRUE)
  if (!length(out)) {
    stop("No public ordinal thresholds were returned by fixef().",
         call. = FALSE)
  }
  index <- as.integer(sub("^.*\\[([0-9]+)\\]$", "\\1", out))
  out[order(index)]
}

append_matrix_samples <- function(out, value, prefix) {
  value <- as.matrix(value)
  selected <- selected_indices(ncol(value))
  for (i in selected) {
    out[[sprintf("%s:observation[%d]", prefix, i)]] <- value[, i]
  }
  out[[paste0(prefix, ":observation-average")]] <- rowMeans(value)
  out
}

append_ordinal_probabilities <- function(out, value, prefix) {
  if (length(dim(value)) != 3L) {
    stop(prefix, " did not return draws x observations x categories.",
         call. = FALSE)
  }
  selected <- selected_indices(dim(value)[2L])
  for (i in selected) {
    for (k in seq_len(dim(value)[3L])) {
      out[[sprintf("%s:observation[%d],category[%d]", prefix, i, k)]] <-
        value[, i, k]
    }
  }
  for (k in seq_len(dim(value)[3L])) {
    out[[sprintf("%s:observation-average,category[%d]", prefix, k)]] <-
      rowMeans(value[, , k, drop = FALSE], dims = 1L)
  }
  out
}

append_group_samples <- function(out, value, prefix) {
  if (length(dim(value)) != 3L) {
    stop(prefix, " did not return draws x levels x coefficients.",
         call. = FALSE)
  }
  selected <- selected_indices(dim(value)[2L], maximum = 3L)
  for (j in selected) {
    for (k in seq_len(dim(value)[3L])) {
      label <- paste(dimnames(value)[[2L]][j],
                     dimnames(value)[[3L]][k], sep = ",")
      out[[paste0(prefix, ":[", label, "]")]] <- value[, j, k]
    }
  }
  out
}

public_samples <- function(fit, case) {
  fe <- as.matrix(fixef(fit, summary = FALSE))
  internal_name <- grepl(
    "s2z|theta_Intercept|finite_Intercept", colnames(fe)
  )
  if (any(internal_name)) {
    stop("An internal ordinal S2Z coordinate leaked through fixef().",
         call. = FALSE)
  }
  out <- setNames(
    lapply(seq_len(ncol(fe)), function(k) fe[, k]),
    paste0("fixef:", colnames(fe))
  )
  for (name in threshold_columns(fe)) {
    out[[paste0("threshold:", name)]] <- fe[, name]
  }

  re <- ranef(fit, summary = FALSE)$g
  co <- coef(fit, summary = FALSE)$g
  if (is.null(re) || is.null(co)) {
    stop("ranef() or coef() did not expose group g.", call. = FALSE)
  }
  out <- append_group_samples(out, re, "ranef")
  out <- append_group_samples(out, co, "coef")

  hyper_names <- grep(
    "^(sd_g__|cor_g__|df_g$)", variables(fit), value = TRUE
  )
  if (length(hyper_names)) {
    hyper <- exact_draw_matrix(fit, hyper_names)
    for (name in colnames(hyper)) {
      out[[paste0("parameter:", name)]] <- hyper[, name]
    }
  }

  eta <- posterior_linpred(fit, dpar = "mu", transform = FALSE)
  out <- append_matrix_samples(out, eta, "posterior_linpred")
  epred <- posterior_epred(fit)
  out <- append_ordinal_probabilities(out, epred, "posterior_epred")

  set.seed(case$seed + 101L)
  pp <- posterior_predict(fit)
  out <- append_matrix_samples(out, pp, "posterior_predict")
  set.seed(case$seed + 202L)
  pr <- predict(fit, summary = FALSE)
  append_matrix_samples(out, pr, "predict")
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
        propagated_mcse = mcse_difference,
        threshold = equiv_z * mcse_difference,
        z_score = z_score,
        pass = abs(difference) <= equiv_z * mcse_difference,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  out[order(out$z_score, decreasing = TRUE), ]
}

public_linear_predictor <- function(fit, case) {
  fe <- as.matrix(fixef(fit, summary = FALSE))
  re <- ranef(fit, summary = FALSE)$g
  group_index <- as.integer(case$data$g)
  eta <- outer(fe[, "x"], case$data$x)
  if ("Intercept" %in% case$group_coef) {
    eta <- eta + re[, group_index, "Intercept"]
  }
  if ("x" %in% case$group_coef) {
    eta <- eta + sweep(
      re[, group_index, "x"], 2L, case$data$x, "*"
    )
  }
  eta
}

public_prediction_invariants <- function(fit, case, parameterization) {
  fe <- as.matrix(fixef(fit, summary = FALSE))
  thresholds <- fe[, threshold_columns(fe), drop = FALSE]
  eta <- posterior_linpred(fit, dpar = "mu", transform = FALSE)
  eta_public <- public_linear_predictor(fit, case)
  epred <- posterior_epred(fit)
  probability_sum <- apply(epred, c(1L, 2L), sum)
  set.seed(case$seed + 303L)
  pp <- posterior_predict(fit)
  valid_prediction <- all(is.finite(pp)) &&
    all(pp %in% case$response_values)
  threshold_ordered <- all(apply(thresholds, 1L, function(x) all(diff(x) > 0)))
  public_names_ok <- all(
    c(sprintf("b_Intercept[%d]", seq_len(ncol(thresholds))), "b_x") %in%
      variables(fit)
  )
  no_internal_fixef <- !any(grepl(
    "s2z|theta_Intercept|finite_Intercept", colnames(fe)
  ))
  rows <- data.frame(
    case = case$name,
    parameterization = parameterization,
    check = c(
      "saved threshold and slope names use the conventional API",
      "fixef excludes internal finite ordinal coordinates",
      "public fixef and ranef reproduce posterior_linpred",
      "ordinal thresholds are strictly ordered",
      "posterior_epred category probabilities sum to one",
      "posterior_predict returns only response categories"
    ),
    max_abs_error = c(
      if (public_names_ok) 0 else Inf,
      if (no_internal_fixef) 0 else Inf,
      max(abs(eta - eta_public)),
      if (threshold_ordered) 0 else Inf,
      max(abs(probability_sum - 1)),
      if (valid_prediction) 0 else Inf
    ),
    tolerance = invariant_tolerance,
    stringsAsFactors = FALSE
  )
  rows$pass <- rows$max_abs_error <= rows$tolerance
  rows
}

s2z_internal_invariants <- function(fit, case) {
  n_group <- nlevels(case$data$g)
  n_threshold <- length(threshold_columns(
    as.matrix(fixef(fit, summary = FALSE))
  ))
  n_q <- n_threshold + 1L
  theta <- exact_draw_matrix(
    fit, sprintf("theta_s2z[%d]", seq_len(n_q))
  )
  recovered <- exact_draw_matrix(
    fit, sprintf("q_recovered_s2z_1[%d]", seq_len(n_q))
  )
  n_draw <- nrow(theta)
  internal <- array(
    NA_real_, dim = c(n_draw, n_group, length(case$group_coef)),
    dimnames = list(NULL, levels(case$data$g), case$group_coef)
  )
  for (k in seq_along(case$group_coef)) {
    internal[, , k] <- exact_draw_matrix(
      fit, sprintf("r_s2z_1_%d[%d]", k, seq_len(n_group))
    )
  }

  fe <- as.matrix(fixef(fit, summary = FALSE))
  threshold_name <- threshold_columns(fe)
  public_threshold <- fe[, threshold_name, drop = FALSE]
  public_re <- ranef(fit, summary = FALSE)$g[
    , levels(case$data$g), case$group_coef, drop = FALSE
  ]
  mean_x <- mean(case$data$x)
  recovered_public <- cbind(
    sweep(public_threshold, 1L, mean_x * fe[, "x"], "-"),
    fe[, "x"]
  )

  max_sum_error <- 0
  max_re_error <- 0
  for (k in seq_along(case$group_coef)) {
    internal_k <- matrix(internal[, , k], nrow = n_draw)
    public_k <- matrix(public_re[, , k], nrow = n_draw)
    omitted_mean <- rowMeans(public_k - internal_k)
    max_sum_error <- max(max_sum_error, max(abs(rowSums(internal_k))))
    max_re_error <- max(
      max_re_error,
      max(abs(sweep(public_k - internal_k, 1L, omitted_mean, "-")))
    )
  }

  group_index <- as.integer(case$data$g)
  finite_eta <- outer(theta[, n_q], case$data$x - mean_x)
  for (k in seq_along(case$group_coef)) {
    z <- if (case$group_coef[k] == "Intercept") {
      rep(1, nrow(case$data))
    } else if (case$group_coef[k] == "x") {
      case$data$x
    } else {
      stop("Unregistered group coefficient in ordinal invariant.",
           call. = FALSE)
    }
    finite_eta <- finite_eta + sweep(
      internal[, group_index, k], 2L, z, "*"
    )
  }
  public_eta <- public_linear_predictor(fit, case)
  max_gap_error <- 0
  for (k in seq_len(n_threshold)) {
    finite_gap <- sweep(-finite_eta, 1L, theta[, k], "+")
    public_gap <- sweep(-public_eta, 1L, public_threshold[, k], "+")
    max_gap_error <- max(max_gap_error, max(abs(finite_gap - public_gap)))
  }

  rows <- data.frame(
    case = case$name,
    parameterization = "s2z",
    check = c(
      "internal physical group effects sum exactly to zero",
      "public ranef equals physical ranef plus one omitted mean",
      "recovered thresholds and slope equal conventional fixef coordinates",
      "finite and conventional threshold-minus-predictor gaps agree"
    ),
    max_abs_error = c(
      max_sum_error,
      max_re_error,
      max(abs(recovered - recovered_public)),
      max_gap_error
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
    "^(b_|sd_g__|cor_g__|df_g$|r_g\\[)",
    variables(fit), value = TRUE
  )
  if (!length(public)) {
    stop("No public parameters were available for sampler diagnostics.",
         call. = FALSE)
  }
  summary <- posterior::summarise_draws(
    as_draws_array(fit, variable = public),
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  )
  finite_rhat <- summary$rhat[is.finite(summary$rhat)]
  finite_bulk <- summary$ess_bulk[is.finite(summary$ess_bulk)]
  finite_tail <- summary$ess_tail[is.finite(summary$ess_tail)]
  if (!length(finite_rhat) || !length(finite_bulk) || !length(finite_tail)) {
    stop("Rhat or ESS diagnostics were unavailable.", call. = FALSE)
  }
  nuts <- nuts_params(fit)
  divergences <- sum(nuts$Value[nuts$Parameter == "divergent__"])
  treedepth_hits <- sum(
    nuts$Value[nuts$Parameter == "treedepth__"] >= max_treedepth
  )
  ebfmi <- chain_ebfmi(nuts)
  observed_max_rhat <- max(finite_rhat)
  observed_min_bulk <- min(finite_bulk)
  observed_min_tail <- min(finite_tail)
  observed_min_ebfmi <- min(ebfmi, na.rm = TRUE)
  data.frame(
    case = case,
    parameterization = parameterization,
    max_rhat = observed_max_rhat,
    min_bulk_ess = observed_min_bulk,
    min_tail_ess = observed_min_tail,
    divergences = divergences,
    treedepth_hits = treedepth_hits,
    min_ebfmi = observed_min_ebfmi,
    pass = observed_max_rhat <= max_rhat &&
      observed_min_bulk >= min_bulk_ess &&
      observed_min_tail >= min_tail_ess &&
      divergences <= max_divergences &&
      treedepth_hits <= max_treedepth_hits &&
      observed_min_ebfmi >= min_ebfmi,
    stringsAsFactors = FALSE
  )
}

new_level_samples <- function(fit, case, seed) {
  first_level <- levels(case$data$g)[1L]
  newdata <- data.frame(
    x = c(-0.35, 0.85, 0.20),
    z = c(-0.50, 0.30, 1.10),
    g = c(first_level, "new-group-a", "new-group-b")
  )
  args <- list(
    object = fit,
    newdata = newdata,
    allow_new_levels = TRUE,
    sample_new_levels = "gaussian"
  )
  set.seed(seed)
  eta <- do.call(posterior_linpred, c(args, list(
    dpar = "mu", transform = FALSE
  )))
  set.seed(seed + 1L)
  epred <- do.call(posterior_epred, args)
  set.seed(seed + 2L)
  predicted <- do.call(predict, c(args, list(summary = FALSE)))
  out <- list()
  out <- append_matrix_samples(out, eta, "new-level:posterior_linpred")
  out <- append_ordinal_probabilities(
    out, epred, "new-level:posterior_epred"
  )
  append_matrix_samples(out, predicted, "new-level:predict")
}

max_abs_difference <- function(x, y) {
  if (!identical(dim(x), dim(y)) || !identical(dimnames(x), dimnames(y))) {
    return(Inf)
  }
  max(abs(x - y))
}

representative_api_checks <- function(fit, case) {
  fe <- as.matrix(fixef(fit, summary = FALSE))
  re <- ranef(fit, summary = FALSE)$g
  co <- coef(fit, summary = FALSE)$g
  eta <- posterior_linpred(fit, dpar = "mu", transform = FALSE)
  epred <- posterior_epred(fit)

  serialization_file <- tempfile(fileext = ".rds")
  saveRDS(fit, serialization_file)
  restored <- readRDS(serialization_file)
  unlink(serialization_file)
  set.seed(case$seed + 404L)
  prediction <- predict(fit, summary = FALSE)
  set.seed(case$seed + 404L)
  restored_prediction <- predict(restored, summary = FALSE)
  serialization_error <- max(
    max_abs_difference(as.matrix(fixef(restored, summary = FALSE)), fe),
    max_abs_difference(ranef(restored, summary = FALSE)$g, re),
    max_abs_difference(coef(restored, summary = FALSE)$g, co),
    max_abs_difference(
      posterior_linpred(restored, dpar = "mu", transform = FALSE), eta
    ),
    max_abs_difference(posterior_epred(restored), epred),
    max_abs_difference(restored_prediction, prediction)
  )

  updated <- update(
    fit, newdata = case$data, recompile = FALSE,
    testmode = TRUE, silent = 2
  )
  update_ok <- inherits(updated, "brmsfit") &&
    isTRUE(all.equal(updated$formula, fit$formula)) &&
    identical(nrow(updated$data), nrow(case$data)) &&
    all(updated$data$g == case$data$g)

  rows <- data.frame(
    case = case$name,
    check = c(
      "serialization preserves thresholds/fixef/ranef/coef and predictions",
      "update validation preserves the ordinal S2Z formula and data"
    ),
    max_abs_error = c(
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

results <- lapply(cases, run_case)
comparison_results <- do.call(rbind, lapply(results, `[[`, "comparison"))
invariant_results <- do.call(rbind, lapply(results, `[[`, "invariants"))
quality_results <- do.call(rbind, lapply(results, `[[`, "quality"))

# The first case is the representative API lifecycle fixture. Its newdata has
# one observed group and two unseen groups, so the same comparison covers both
# observed-level and Gaussian-sampled new-level behavior.
new_level_results <- compare_public_posteriors(
  new_level_samples(results[[1L]]$conventional, cases[[1L]], 10201L),
  new_level_samples(results[[1L]]$s2z, cases[[1L]], 10201L),
  case = paste0(cases[[1L]]$name, "-new-levels")
)
comparison_results <- rbind(comparison_results, new_level_results)
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
  cat(
    "... ", nrow(comparison_results) - report_rows,
    " additional comparisons retained for assertions.\n", sep = ""
  )
}
cat("\nExact public-prediction and ordinal S2Z reconstruction identities\n")
print(invariant_results, row.names = FALSE, digits = 4)
cat("\nSampling quality\n")
print(quality_results, row.names = FALSE, digits = 4)
cat("\nRepresentative serialization/update checks\n")
print(api_results, row.names = FALSE, digits = 4)

failed_comparisons <- subset(comparison_results, !pass)
failed_invariants <- subset(invariant_results, !pass)
failed_quality <- subset(quality_results, !pass)
failed_api <- subset(api_results, !pass)
if (nrow(failed_comparisons) || nrow(failed_invariants) ||
    nrow(failed_quality) || nrow(failed_api)) {
  stop(
    "PLAN-02 ordinal S2Z sampling validation failed: ",
    nrow(failed_comparisons), " MCSE comparisons, ",
    nrow(failed_invariants), " exact invariants, ",
    nrow(failed_quality), " sampling-quality checks, and ",
    nrow(failed_api), " public API lifecycle checks failed.",
    call. = FALSE
  )
}

cat(
  "\nPASS: ordinal conventional and physical-S2Z posteriors agree within ",
  "the prespecified propagated-MCSE thresholds; internal reconstruction, ",
  "observed/new-level prediction, serialization/update, and all sampler ",
  "diagnostic gates pass.\n",
  sep = ""
)
