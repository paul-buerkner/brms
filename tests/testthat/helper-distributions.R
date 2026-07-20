# Distribution registry, prep fixtures, and d/p/q/r expectations.
# posterior_predict expectations live in helper-posterior-predict.R
# (sourced after this file alphabetically).
#
# How to add a family
# -------------------
# 1. Register an entry with `dist_registry_add()` (or append in
#    `dist_registry_populate()`). Required fields:
#      name, backend, type, d/p/q/r (NULL if missing), params, support,
#      pp_fun, prep_builder, outputs, flags
# 2. Map `backend` to the name used by predict_*_helper / compute_*
#    (e.g. gaussian -> "norm", poisson -> "pois").
# 3. Set `type` to continuous|discrete|mixture|ordinal|circular.
# 4. List supported `outputs`: random|probability|pit|density|quantile.
# 5. Flags (logical list), common ones:
#      truncation          - exercise truncated CDF/PDF/quantile
#      numeric_q           - quantile solved numerically
#      no_random_truncation- random sampling ignores truncation (warns)
#      zi / hurdle         - mixture formulas vs a baseline d/p
#      skip_integrate      - skip integrate-to-one (e.g. circular, heavy tails)
#      skip_moments        - skip moment checks
#      skip_rng_cdf        - skip RNG-vs-CDF KS-style checks
#      stub                - registry placeholder without deep d/p/q checks
# 6. Implement `prep_builder()` returning a minimal
#    `structure(list(...), class = "brmsprep")` fixture.
#    Prefer constant dpars so PP vs d/p/q comparisons are exact.
#    entry$params MUST match args PP passes to the backend.
# 7. Prefer consolidating analytical / special cases into
#    tests.distributions_analytical.R.
#
# Families are registered when this helper is sourced via dist_registry_populate().

# ---------------------------------------------------------------------------
# Registry infrastructure
# ---------------------------------------------------------------------------

.dist_registry_env <- new.env(parent = emptyenv())
.dist_registry_env$entries <- list()

#' Create a family registry entry
#' @noRd
dist_entry <- function(
  name,
  backend,
  type = c("continuous", "discrete", "mixture", "ordinal", "circular"),
  d = NULL,
  p = NULL,
  q = NULL,
  r = NULL,
  params = list(),
  support = c(-Inf, Inf),
  ref = NULL,
  pp_fun = NULL,
  prep_builder = NULL,
  outputs = c("random", "probability", "pit", "density", "quantile"),
  flags = list(),
  q_ref = NULL,
  p_ref = 0.7,
  baseline = NULL
) {
  type <- match.arg(type)
  flags_default <- list(
    truncation = TRUE,
    numeric_q = FALSE,
    no_random_truncation = FALSE,
    zi = FALSE,
    hurdle = FALSE,
    skip_integrate = FALSE,
    skip_moments = FALSE,
    skip_rng_cdf = FALSE,
    skip_pdf_fd = FALSE,
    skip_p_tail_flags = FALSE,
    skip_d_sums = FALSE,
    pq_elementwise = FALSE,
    has_atoms = FALSE,
    stub = FALSE,
    discrete_support = type %in% c("discrete", "ordinal")
  )
  flags <- utils::modifyList(flags_default, flags)
  if (is.null(q_ref)) {
    q_ref <- if (flags$discrete_support) {
      max(0, ceiling(mean(support[is.finite(support)])))
    } else {
      0
    }
  }
  list(
    name = name,
    backend = backend,
    type = type,
    d = d,
    p = p,
    q = q,
    r = r,
    params = params,
    support = support,
    ref = ref,
    pp_fun = pp_fun,
    prep_builder = prep_builder,
    outputs = outputs,
    flags = flags,
    q_ref = q_ref,
    p_ref = p_ref,
    baseline = baseline
  )
}

dist_registry_clear <- function() {
  .dist_registry_env$entries <- list()
  invisible(.dist_registry_env$entries)
}

dist_registry_add <- function(entry) {
  stopifnot(is.list(entry), !is.null(entry$name))
  .dist_registry_env$entries[[entry$name]] <- entry
  invisible(entry)
}

dist_registry_get <- function(name = NULL, type = NULL, flag = NULL) {
  entries <- .dist_registry_env$entries
  if (!is.null(name)) {
    entries <- entries[intersect(name, names(entries))]
  }
  if (!is.null(type)) {
    entries <- Filter(function(e) e$type %in% type, entries)
  }
  if (!is.null(flag)) {
    entries <- Filter(function(e) isTRUE(e$flags[[flag]]), entries)
  }
  entries
}

dist_registry_names <- function(...) {
  names(dist_registry_get(...))
}

has_output <- function(entry, output) {
  output %in% entry$outputs
}

has_fun <- function(entry, which = c("d", "p", "q", "r", "pp_fun")) {
  which <- match.arg(which)
  !is.null(entry[[which]])
}

entry_requires <- function(entry) {
  req <- entry$ref$requires
  if (is.null(req)) character() else req
}

entry_deps_available <- function(entry) {
  req <- entry_requires(entry)
  if (!length(req)) return(TRUE)
  all(vapply(req, requireNamespace, logical(1), quietly = TRUE))
}

# Non-stub registry entries with dependencies available.
# Additional filters: type (passed to dist_registry_get), require_funs (d/p/q/r/pp_fun),
# truncation / discrete_support (match flags), has_baseline (non-NULL baseline).
dist_active_entries <- function(
    type = NULL,
    require_funs = NULL,
    truncation = NULL,
    discrete_support = NULL,
    has_baseline = NULL
) {
  Filter(function(e) {
    if (isTRUE(e$flags$stub)) return(FALSE)
    if (!entry_deps_available(e)) return(FALSE)
    if (!is.null(require_funs)) {
      ok_funs <- vapply(require_funs, function(f) has_fun(e, f), logical(1))
      if (!all(ok_funs)) return(FALSE)
    }
    if (!is.null(truncation) && !identical(isTRUE(e$flags$truncation), truncation)) {
      return(FALSE)
    }
    if (!is.null(discrete_support) &&
          !identical(isTRUE(e$flags$discrete_support), discrete_support)) {
      return(FALSE)
    }
    if (!is.null(has_baseline) && !identical(!is.null(e$baseline), has_baseline)) {
      return(FALSE)
    }
    TRUE
  }, dist_registry_get(type = type))
}

call_dist <- function(fun, x, params, ...) {
  do.call(fun, c(list(x), params, list(...)))
}

# ---------------------------------------------------------------------------
# Prep builders (minimal brmsprep fixtures)
# ---------------------------------------------------------------------------

.empty_bounds <- function(nobs) {
  list(lb = rep(NULL, nobs), ub = rep(NULL, nobs))
}

make_prep_location <- function(ns = 25, nobs = 4, mu = 0, sigma = 1,
                               extra = list(), Y = NULL, seed = NULL) {
  if (!is.null(seed)) withr::local_seed(seed)
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- c(
    list(
      mu = matrix(mu, nrow = ns, ncol = nobs),
      sigma = rep(sigma, ns)
    ),
    extra
  )
  if (is.null(Y)) Y <- rep(mu, nobs)
  prep$data <- c(list(Y = Y), .empty_bounds(nobs))
  prep
}

make_prep_positive <- function(ns = 25, nobs = 4, mu = 1, shape = 2,
                               extra = list(), Y = NULL, seed = NULL) {
  if (!is.null(seed)) withr::local_seed(seed)
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- c(
    list(
      mu = matrix(mu, nrow = ns, ncol = nobs),
      shape = rep(shape, ns),
      sigma = rep(1, ns),
      beta = rep(1, ns)
    ),
    extra
  )
  if (is.null(Y)) Y <- rep(mu, nobs)
  prep$data <- c(list(Y = Y), .empty_bounds(nobs))
  prep
}

make_prep_count <- function(ns = 25, nobs = 4, mu = 3, shape = 2,
                            trials = 10, zi = 0.2, hu = 0.2,
                            extra = list(), Y = NULL, seed = NULL) {
  if (!is.null(seed)) withr::local_seed(seed)
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- c(
    list(
      mu = matrix(mu, nrow = ns, ncol = nobs),
      shape = rep(shape, ns),
      phi = rep(5, ns),
      zi = rep(zi, ns),
      hu = rep(hu, ns)
    ),
    extra
  )
  if (is.null(Y)) Y <- rep(as.integer(mu), nobs)
  prep$data <- c(
    list(Y = Y, trials = rep(trials, nobs)),
    .empty_bounds(nobs)
  )
  prep
}

make_prep_beta_mix <- function(ns = 25, nobs = 4, mu = 0.4, phi = 5,
                               zoi = 0.2, coi = 0.4, Y = NULL, seed = NULL) {
  if (!is.null(seed)) withr::local_seed(seed)
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(mu, nrow = ns, ncol = nobs),
    phi = rep(phi, ns),
    zoi = rep(zoi, ns),
    coi = rep(coi, ns)
  )
  if (is.null(Y)) Y <- rep(mu, nobs)
  prep$data <- c(list(Y = Y), .empty_bounds(nobs))
  prep
}

make_prep_ordinal <- function(ns = 25, nobs = 4, nthres = 3,
                              family = "cumulative", link = "logit",
                              hu = NULL, Y = NULL, seed = NULL) {
  if (!is.null(seed)) withr::local_seed(seed)
  ncat <- nthres + 1L
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(
    mu = matrix(0, nrow = ns, ncol = nobs),
    disc = rep(1, ns)
  )
  if (!is.null(hu)) {
    prep$dpars$hu <- rep(hu, ns)
  }
  prep$thres$thres <- matrix(rep(c(-1, 0, 1)[seq_len(nthres)], each = ns),
                             nrow = ns)
  if (is.null(Y)) {
    Y <- if (is.null(hu)) {
      rep(seq_len(ncat), length.out = nobs)
    } else {
      rep(0:ncat, length.out = nobs)
    }
  }
  prep$data <- c(list(Y = Y, ncat = ncat), .empty_bounds(nobs))
  prep$family$family <- family
  prep$family$link <- link
  prep
}

make_prep_categorical <- function(ns = 25, nobs = 4, ncat = 3,
                                  Y = NULL, seed = NULL) {
  if (!is.null(seed)) withr::local_seed(seed)
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  # reference category inserted by posterior_predict_categorical
  mu_names <- paste0("mu", seq_len(ncat - 1L))
  prep$dpars <- setNames(
    lapply(mu_names, function(nm) matrix(0, nrow = ns, ncol = nobs)),
    mu_names
  )
  if (is.null(Y)) Y <- rep(seq_len(ncat), length.out = nobs)
  prep$data <- c(list(Y = Y, ncat = ncat), .empty_bounds(nobs))
  prep$family <- categorical()
  prep$refcat <- 1L
  prep
}

set_trunc_bounds <- function(prep, lb = NULL, ub = NULL, i = NULL) {
  nobs <- prep$nobs
  if (is.null(i)) i <- seq_len(nobs)
  # Match existing fixtures: numeric vectors when truncating; indexing
  # with prep$data$lb[i] must return a scalar (see predict_*_helper).
  if (!is.null(lb)) {
    lb_vec <- rep(lb, nobs)
    prep$data$lb <- lb_vec
  }
  if (!is.null(ub)) {
    ub_vec <- rep(ub, nobs)
    prep$data$ub <- ub_vec
  }
  prep
}

# ---------------------------------------------------------------------------
# Registry population
# ---------------------------------------------------------------------------

dist_registry_populate <- function(reset = TRUE) {
  if (reset) dist_registry_clear()

  dist_registry_add(dist_entry(
    name = "gaussian",
    backend = "norm",
    type = "continuous",
    d = stats::dnorm,
    p = stats::pnorm,
    q = stats::qnorm,
    r = stats::rnorm,
    params = list(mean = 0, sd = 1),
    support = c(-Inf, Inf),
    ref = list(d = stats::dnorm, p = stats::pnorm, q = stats::qnorm),
    pp_fun = brms:::posterior_predict_gaussian,
    prep_builder = function(...) make_prep_location(mu = 0, sigma = 1, ...),
    q_ref = 0.25,
    flags = list(truncation = TRUE)
  ))

  dist_registry_add(dist_entry(
    name = "student_t",
    backend = "student_t",
    type = "continuous",
    d = dstudent_t,
    p = pstudent_t,
    q = qstudent_t,
    r = rstudent_t,
    params = list(df = 7, mu = 0, sigma = 1),
    support = c(-Inf, Inf),
    pp_fun = brms:::posterior_predict_student,
    prep_builder = function(ns = 25, ...) {
      make_prep_location(
        ns = ns, mu = 0, sigma = 1,
        extra = list(nu = rep(7, ns)),
        ...
      )
    },
    q_ref = 0.25,
    flags = list(truncation = TRUE)
  ))

  dist_registry_add(dist_entry(
    name = "inv_gaussian",
    backend = "inv_gaussian",
    type = "continuous",
    d = dinv_gaussian,
    p = pinv_gaussian,
    q = qinv_gaussian,
    r = rinv_gaussian,
    params = list(mu = 1, shape = 2),
    support = c(0, Inf),
    pp_fun = brms:::posterior_predict_inverse.gaussian,
    prep_builder = function(ns = 25, nobs = 4, ...) {
      make_prep_positive(
        ns = ns, nobs = nobs, mu = 1, shape = 2,
        Y = rep(1, nobs), ...
      )
    },
    q_ref = 1.2,
    flags = list(truncation = TRUE, numeric_q = TRUE)
  ))

  dist_registry_add(dist_entry(
    name = "exgaussian",
    backend = "exgaussian",
    type = "continuous",
    d = dexgaussian,
    p = pexgaussian,
    q = qexgaussian,
    r = rexgaussian,
    params = list(mu = 0, sigma = 1, beta = 1),
    support = c(-Inf, Inf),
    pp_fun = brms:::posterior_predict_exgaussian,
    prep_builder = function(ns = 25, nobs = 4, ...) {
      make_prep_positive(
        ns = ns, nobs = nobs, mu = 0, shape = 1,
        extra = list(beta = rep(1, ns), sigma = rep(1, ns)),
        Y = rep(0, nobs),
        ...
      )
    },
    q_ref = 0.5,
    flags = list(truncation = TRUE, numeric_q = TRUE)
  ))

  dist_registry_add(dist_entry(
    name = "von_mises",
    backend = "von_mises",
    type = "circular",
    d = dvon_mises,
    p = pvon_mises,
    q = qvon_mises,
    r = rvon_mises,
    params = list(mu = 0, kappa = 2),
    support = c(-pi, pi),
    pp_fun = brms:::posterior_predict_von_mises,
    prep_builder = function(ns = 25, nobs = 4, ...) {
      make_prep_location(
        ns = ns, nobs = nobs, mu = 0, sigma = 1,
        extra = list(kappa = rep(2, ns)),
        Y = rep(0, nobs),
        ...
      )
    },
    q_ref = 0.2,
    flags = list(
      truncation = FALSE,
      numeric_q = TRUE,
      skip_integrate = FALSE,
      skip_moments = TRUE,
      skip_pdf_fd = TRUE
    )
  ))

  dist_registry_add(dist_entry(
    name = "pois",
    backend = "pois",
    type = "discrete",
    d = stats::dpois,
    p = stats::ppois,
    q = stats::qpois,
    r = stats::rpois,
    params = list(lambda = 3),
    support = c(0, Inf),
    pp_fun = brms:::posterior_predict_poisson,
    prep_builder = function(ns = 25, nobs = 4, ...) {
      make_prep_count(ns = ns, nobs = nobs, mu = 3, Y = rep(3L, nobs), ...)
    },
    q_ref = 3,
    flags = list(truncation = TRUE)
  ))

  dist_registry_add(dist_entry(
    name = "binom",
    backend = "binom",
    type = "discrete",
    d = stats::dbinom,
    p = stats::pbinom,
    q = stats::qbinom,
    r = stats::rbinom,
    params = list(size = 10, prob = 0.4),
    support = c(0, 10),
    pp_fun = brms:::posterior_predict_binomial,
    prep_builder = function(ns = 25, nobs = 4, ...) {
      make_prep_count(
        ns = ns, nobs = nobs, mu = 0.4, trials = 10,
        Y = rep(4L, nobs), ...
      )
    },
    q_ref = 4,
    flags = list(truncation = TRUE)
  ))

  dist_registry_add(dist_entry(
    name = "beta_binomial",
    backend = "beta_binomial",
    type = "discrete",
    d = brms:::dbeta_binomial,
    p = brms:::pbeta_binomial,
    q = brms:::qbeta_binomial,
    r = brms:::rbeta_binomial,
    params = list(size = 12, mu = 0.35, phi = 8),
    support = c(0, 12),
    pp_fun = brms:::posterior_predict_beta_binomial,
    prep_builder = function(ns = 25, nobs = 4, ...) {
      p <- make_prep_count(
        ns = ns, nobs = nobs, mu = 0.35, trials = 12,
        Y = rep(4L, nobs), ...
      )
      p$dpars$phi <- rep(8, p$ndraws)
      p
    },
    q_ref = 4,
    flags = list(truncation = TRUE),
    # beta_binomial d/p/q/r require extraDistr
    ref = list(requires = "extraDistr")
  ))

  dist_registry_add(dist_entry(
    name = "nbinom",
    backend = "nbinom",
    type = "discrete",
    d = stats::dnbinom,
    p = stats::pnbinom,
    q = stats::qnbinom,
    r = stats::rnbinom,
    params = list(mu = 4, size = 2.5),
    support = c(0, Inf),
    pp_fun = brms:::posterior_predict_negbinomial,
    prep_builder = function(ns = 25, nobs = 4, ...) {
      make_prep_count(
        ns = ns, nobs = nobs, mu = 4, shape = 2.5,
        Y = rep(4L, nobs), ...
      )
    },
    q_ref = 3,
    flags = list(truncation = TRUE)
  ))

  dist_registry_add(dist_entry(
    name = "com_poisson",
    backend = "com_poisson",
    type = "discrete",
    d = brms:::dcom_poisson,
    p = brms:::pcom_poisson,
    q = brms:::qcom_poisson,
    r = brms:::rcom_poisson,
    params = list(mu = 2, shape = 0.8),
    support = c(0, Inf),
    pp_fun = brms:::posterior_predict_com_poisson,
    prep_builder = function(ns = 25, nobs = 4, ...) {
      make_prep_count(
        ns = ns, nobs = nobs, mu = 2, shape = 0.8,
        Y = rep(2L, nobs), ...
      )
    },
    q_ref = 2,
    flags = list(
      truncation = TRUE,
      skip_moments = TRUE,
      skip_rng_cdf = TRUE,
      skip_d_sums = TRUE,
      pq_elementwise = TRUE
    )
  ))

  dist_registry_add(dist_entry(
    name = "hurdle_poisson",
    backend = "hurdle_poisson",
    type = "discrete",
    d = brms:::dhurdle_poisson,
    p = brms:::phurdle_poisson,
    q = brms:::qhurdle_poisson,
    r = NULL,
    params = list(lambda = 2, hu = 0.2),
    support = c(0, Inf),
    pp_fun = brms:::posterior_predict_hurdle_poisson,
    prep_builder = function(ns = 25, nobs = 4, ...) {
      make_prep_count(
        ns = ns, nobs = nobs, mu = 2, hu = 0.2,
        Y = rep(2L, nobs), ...
      )
    },
    q_ref = 2,
    baseline = list(
      d = stats::dpois,
      p = stats::ppois,
      params = list(lambda = 2),
      mix = "hu"
    ),
    flags = list(
      truncation = TRUE,
      hurdle = TRUE,
      no_random_truncation = TRUE,
      skip_moments = TRUE
    )
  ))

  dist_registry_add(dist_entry(
    name = "zero_inflated_negbinomial",
    backend = "zero_inflated_negbinomial",
    type = "discrete",
    d = brms:::dzero_inflated_negbinomial,
    p = brms:::pzero_inflated_negbinomial,
    q = brms:::qzero_inflated_negbinomial,
    r = NULL,
    params = list(mu = 4, shape = 2.5, zi = 0.3),
    support = c(0, Inf),
    pp_fun = brms:::posterior_predict_zero_inflated_negbinomial,
    prep_builder = function(ns = 25, nobs = 4, ...) {
      make_prep_count(
        ns = ns, nobs = nobs, mu = 4, shape = 2.5, zi = 0.3,
        Y = rep(3L, nobs), ...
      )
    },
    q_ref = 3,
    baseline = list(
      d = stats::dnbinom,
      p = stats::pnbinom,
      params = list(mu = 4, size = 2.5),
      mix = "zi"
    ),
    flags = list(
      truncation = TRUE,
      zi = TRUE,
      skip_moments = TRUE
    )
  ))

  dist_registry_add(dist_entry(
    name = "zero_one_inflated_beta",
    backend = "zero_one_inflated_beta",
    type = "mixture",
    d = brms:::dzero_one_inflated_beta,
    p = brms:::pzero_one_inflated_beta,
    q = brms:::qzero_one_inflated_beta,
    r = NULL,
    params = list(shape1 = 2, shape2 = 3, zoi = 0.25, coi = 0.4),
    support = c(0, 1),
    pp_fun = brms:::posterior_predict_zero_one_inflated_beta,
    prep_builder = function(...) {
      # mu * phi = shape1, (1-mu)*phi = shape2 => mu=2/5, phi=5
      make_prep_beta_mix(mu = 0.4, phi = 5, zoi = 0.25, coi = 0.4, ...)
    },
    q_ref = 0.4,
    flags = list(
      truncation = FALSE,
      skip_integrate = TRUE,
      skip_moments = TRUE,
      skip_pdf_fd = TRUE,
      has_atoms = TRUE
    )
  ))

  dist_registry_add(dist_entry(
    name = "zero_inflated_asym_laplace",
    backend = "zero_inflated_asym_laplace",
    type = "mixture",
    d = brms:::dzero_inflated_asym_laplace,
    p = brms:::pzero_inflated_asym_laplace,
    q = brms:::qzero_inflated_asym_laplace,
    r = NULL,
    params = list(mu = 0, sigma = 1, quantile = 0.5, zi = 0.3),
    support = c(-Inf, Inf),
    pp_fun = brms:::posterior_predict_zero_inflated_asym_laplace,
    prep_builder = function(ns = 25, ...) {
      make_prep_location(
        ns = ns, mu = 0, sigma = 1,
        extra = list(
          quantile = rep(0.5, ns),
          zi = rep(0.3, ns)
        ),
        ...
      )
    },
    q_ref = 0.25,
    flags = list(
      truncation = FALSE,
      zi = TRUE,
      skip_integrate = TRUE,
      skip_moments = TRUE,
      skip_pdf_fd = TRUE,
      has_atoms = TRUE
    )
  ))

  dist_registry_add(dist_entry(
    name = "ordinal_cumulative",
    backend = "ordinal",
    type = "ordinal",
    d = function(x, ...) {
      brms:::dordinal(x, family = "cumulative", link = "logit", ...)
    },
    p = function(q, ...) {
      brms:::pordinal(q, family = "cumulative", link = "logit", ...)
    },
    q = function(p, ...) {
      brms:::qordinal(p, family = "cumulative", link = "logit", ...)
    },
    r = function(n, ...) {
      brms:::rordinal(n, family = "cumulative", link = "logit", ...)
    },
    params = list(
      eta = 0,
      thres = matrix(c(-1, 0, 1), nrow = 1),
      disc = 1
    ),
    support = c(1, 4),
    pp_fun = brms:::posterior_predict_ordinal,
    prep_builder = function(...) {
      make_prep_ordinal(family = "cumulative", ...)
    },
    q_ref = 2,
    flags = list(
      truncation = FALSE,
      skip_integrate = TRUE,
      skip_moments = TRUE,
      skip_pdf_fd = TRUE,
      skip_rng_cdf = TRUE,
      skip_p_tail_flags = TRUE
    )
  ))

  dist_registry_add(dist_entry(
    name = "categorical",
    backend = "categorical",
    type = "ordinal",
    d = brms:::dcategorical,
    p = brms:::pcategorical,
    q = brms:::qcategorical,
    r = brms:::rcategorical,
    params = list(eta = matrix(c(0, 0, 0), nrow = 1)),
    support = c(1, 3),
    pp_fun = brms:::posterior_predict_categorical,
    prep_builder = function(...) make_prep_categorical(ncat = 3, ...),
    q_ref = 2,
    flags = list(
      truncation = FALSE,
      skip_integrate = TRUE,
      skip_moments = TRUE,
      skip_pdf_fd = TRUE,
      skip_rng_cdf = TRUE,
      skip_p_tail_flags = TRUE
    )
  ))

  dist_registry_add(dist_entry(
    name = "hurdle_cumulative",
    backend = "hurdle_cumulative",
    type = "ordinal",
    d = function(x, ...) brms:::dhurdle_cumulative(x, link = "logit", ...),
    p = function(q, ...) brms:::phurdle_cumulative(q, link = "logit", ...),
    q = function(p, ...) brms:::qhurdle_cumulative(p, link = "logit", ...),
    r = function(n, ...) brms:::rhurdle_cumulative(n, link = "logit", ...),
    params = list(
      eta = 0,
      thres = matrix(c(-1, 0, 1), nrow = 1),
      hu = 0.25,
      disc = 1
    ),
    support = c(0, 4),
    pp_fun = brms:::posterior_predict_hurdle_cumulative,
    prep_builder = function(...) {
      make_prep_ordinal(family = "hurdle_cumulative", hu = 0.25, ...)
    },
    q_ref = 2,
    flags = list(
      truncation = FALSE,
      hurdle = TRUE,
      skip_integrate = TRUE,
      skip_moments = TRUE,
      skip_pdf_fd = TRUE,
      skip_rng_cdf = TRUE,
      skip_p_tail_flags = TRUE
    )
  ))

  # ---- Continuous / positive ----

  dist_registry_add(dist_entry(
    name = "lognormal",
    backend = "lnorm",
    type = "continuous",
    d = stats::dlnorm,
    p = stats::plnorm,
    q = stats::qlnorm,
    r = stats::rlnorm,
    params = list(meanlog = 0, sdlog = 1),
    support = c(0, Inf),
    ref = list(d = stats::dlnorm, p = stats::plnorm, q = stats::qlnorm),
    pp_fun = brms:::posterior_predict_lognormal,
    prep_builder = function(...) make_prep_location(mu = 0, sigma = 1, ...),
    q_ref = 1,
    flags = list(truncation = TRUE)
  ))

  dist_registry_add(dist_entry(
    name = "shifted_lognormal",
    backend = "shifted_lnorm",
    type = "continuous",
    d = dshifted_lnorm,
    p = pshifted_lnorm,
    q = qshifted_lnorm,
    r = rshifted_lnorm,
    params = list(meanlog = 0, sdlog = 1, shift = 0.3),
    support = c(0.3, Inf),
    pp_fun = brms:::posterior_predict_shifted_lognormal,
    prep_builder = function(ns = 25, ...) {
      make_prep_location(
        ns = ns, mu = 0, sigma = 1,
        extra = list(ndt = rep(0.3, ns)),
        ...
      )
    },
    q_ref = 1.3,
    flags = list(truncation = TRUE)
  ))

  dist_registry_add(dist_entry(
    name = "skew_normal",
    backend = "skew_normal",
    type = "continuous",
    d = dskew_normal,
    p = pskew_normal,
    q = qskew_normal,
    r = rskew_normal,
    params = list(mu = 0, sigma = 1, alpha = 2),
    support = c(-Inf, Inf),
    ref = list(requires = "mnormt"),
    pp_fun = brms:::posterior_predict_skew_normal,
    prep_builder = function(ns = 25, ...) {
      make_prep_location(
        ns = ns, mu = 0, sigma = 1,
        extra = list(alpha = rep(2, ns)),
        ...
      )
    },
    q_ref = 0.5,
    flags = list(
      truncation = TRUE,
      skip_pdf_fd = TRUE,
      skip_moments = TRUE
    )
  ))

  dist_registry_add(dist_entry(
    name = "exponential",
    backend = "exp",
    type = "continuous",
    d = stats::dexp,
    p = stats::pexp,
    q = stats::qexp,
    r = stats::rexp,
    params = list(rate = 1),
    support = c(0, Inf),
    pp_fun = brms:::posterior_predict_exponential,
    prep_builder = function(ns = 25, nobs = 4, ...) {
      make_prep_positive(ns = ns, nobs = nobs, mu = 1, Y = rep(1, nobs), ...)
    },
    q_ref = 0.7,
    flags = list(truncation = TRUE)
  ))

  dist_registry_add(dist_entry(
    name = "gamma",
    backend = "gamma",
    type = "continuous",
    d = stats::dgamma,
    p = stats::pgamma,
    q = stats::qgamma,
    r = stats::rgamma,
    params = list(shape = 2, scale = 1),
    support = c(0, Inf),
    pp_fun = brms:::posterior_predict_gamma,
    prep_builder = function(ns = 25, nobs = 4, ...) {
      # PP: scale = mu / shape => mu = shape * scale = 2
      make_prep_positive(
        ns = ns, nobs = nobs, mu = 2, shape = 2,
        Y = rep(2, nobs), ...
      )
    },
    q_ref = 1.5,
    flags = list(truncation = TRUE)
  ))

  weibull_shape <- 2
  weibull_scale <- 1
  weibull_mu <- weibull_scale * gamma(1 + 1 / weibull_shape)
  dist_registry_add(dist_entry(
    name = "weibull",
    backend = "weibull",
    type = "continuous",
    d = stats::dweibull,
    p = stats::pweibull,
    q = stats::qweibull,
    r = stats::rweibull,
    params = list(shape = weibull_shape, scale = weibull_scale),
    support = c(0, Inf),
    pp_fun = brms:::posterior_predict_weibull,
    prep_builder = function(ns = 25, nobs = 4, ...) {
      make_prep_positive(
        ns = ns, nobs = nobs, mu = weibull_mu, shape = weibull_shape,
        Y = rep(weibull_mu, nobs), ...
      )
    },
    q_ref = 1,
    flags = list(truncation = TRUE)
  ))

  frechet_nu <- 3
  frechet_scale <- 1
  frechet_mu <- frechet_scale * gamma(1 - 1 / frechet_nu)
  dist_registry_add(dist_entry(
    name = "frechet",
    backend = "frechet",
    type = "continuous",
    d = dfrechet,
    p = pfrechet,
    q = qfrechet,
    r = rfrechet,
    params = list(scale = frechet_scale, shape = frechet_nu),
    support = c(0, Inf),
    pp_fun = brms:::posterior_predict_frechet,
    prep_builder = function(ns = 25, nobs = 4, ...) {
      make_prep_positive(
        ns = ns, nobs = nobs, mu = frechet_mu,
        extra = list(nu = rep(frechet_nu, ns)),
        Y = rep(frechet_mu, nobs), ...
      )
    },
    q_ref = 1.2,
    flags = list(truncation = TRUE, skip_moments = TRUE)
  ))

  dist_registry_add(dist_entry(
    name = "gen_extreme_value",
    backend = "gen_extreme_value",
    type = "continuous",
    d = dgen_extreme_value,
    p = pgen_extreme_value,
    q = qgen_extreme_value,
    r = rgen_extreme_value,
    params = list(mu = 0, sigma = 1, xi = 0.1),
    support = c(-Inf, Inf),
    pp_fun = brms:::posterior_predict_gen_extreme_value,
    prep_builder = function(ns = 25, ...) {
      make_prep_location(
        ns = ns, mu = 0, sigma = 1,
        extra = list(xi = rep(0.1, ns)),
        ...
      )
    },
    q_ref = 0.5,
    flags = list(truncation = TRUE, skip_moments = TRUE)
  ))

  dist_registry_add(dist_entry(
    name = "beta",
    backend = "beta",
    type = "continuous",
    d = stats::dbeta,
    p = stats::pbeta,
    q = stats::qbeta,
    r = stats::rbeta,
    # PP: shape1 = mu * phi, shape2 = (1 - mu) * phi
    params = list(shape1 = 2, shape2 = 3),
    support = c(0, 1),
    pp_fun = brms:::posterior_predict_beta,
    prep_builder = function(ns = 25, nobs = 4, ...) {
      make_prep_beta_mix(
        ns = ns, nobs = nobs, mu = 0.4, phi = 5,
        zoi = 0, coi = 0, Y = rep(0.4, nobs), ...
      )
    },
    q_ref = 0.4,
    flags = list(truncation = TRUE)
  ))

  dist_registry_add(dist_entry(
    name = "asym_laplace",
    backend = "asym_laplace",
    type = "continuous",
    d = dasym_laplace,
    p = pasym_laplace,
    q = qasym_laplace,
    r = rasym_laplace,
    params = list(mu = 0, sigma = 1, quantile = 0.5),
    support = c(-Inf, Inf),
    pp_fun = brms:::posterior_predict_asym_laplace,
    prep_builder = function(ns = 25, ...) {
      make_prep_location(
        ns = ns, mu = 0, sigma = 1,
        extra = list(quantile = rep(0.5, ns)),
        ...
      )
    },
    q_ref = 0.25,
    flags = list(truncation = TRUE)
  ))

  # ---- Discrete ----

  dist_registry_add(dist_entry(
    name = "bernoulli",
    backend = "binom",
    type = "discrete",
    d = stats::dbinom,
    p = stats::pbinom,
    q = stats::qbinom,
    r = stats::rbinom,
    params = list(size = 1, prob = 0.4),
    support = c(0, 1),
    pp_fun = brms:::posterior_predict_bernoulli,
    prep_builder = function(ns = 25, nobs = 4, ...) {
      make_prep_count(
        ns = ns, nobs = nobs, mu = 0.4, trials = 1,
        Y = rep(0L, nobs), ...
      )
    },
    q_ref = 0,
    flags = list(truncation = FALSE)
  ))

  dist_registry_add(dist_entry(
    name = "negbinomial2",
    backend = "nbinom",
    type = "discrete",
    d = stats::dnbinom,
    p = stats::pnbinom,
    q = stats::qnbinom,
    r = stats::rnbinom,
    # PP: size = 1 / sigma
    params = list(mu = 4, size = 2.5),
    support = c(0, Inf),
    pp_fun = brms:::posterior_predict_negbinomial2,
    prep_builder = function(ns = 25, nobs = 4, ...) {
      p <- make_prep_count(
        ns = ns, nobs = nobs, mu = 4,
        Y = rep(4L, nobs), ...
      )
      p$dpars$sigma <- rep(1 / 2.5, p$ndraws)
      p
    },
    q_ref = 3,
    flags = list(truncation = TRUE)
  ))

  dist_registry_add(dist_entry(
    name = "geometric",
    backend = "nbinom",
    type = "discrete",
    d = stats::dnbinom,
    p = stats::pnbinom,
    q = stats::qnbinom,
    r = stats::rnbinom,
    params = list(mu = 3, size = 1),
    support = c(0, Inf),
    pp_fun = brms:::posterior_predict_geometric,
    prep_builder = function(ns = 25, nobs = 4, ...) {
      make_prep_count(
        ns = ns, nobs = nobs, mu = 3, shape = 1,
        Y = rep(3L, nobs), ...
      )
    },
    q_ref = 2,
    flags = list(truncation = TRUE)
  ))

  dist_registry_add(dist_entry(
    name = "discrete_weibull",
    backend = "discrete_weibull",
    type = "discrete",
    d = ddiscrete_weibull,
    p = pdiscrete_weibull,
    q = qdiscrete_weibull,
    r = rdiscrete_weibull,
    params = list(mu = 0.7, shape = 1.5),
    support = c(0, Inf),
    pp_fun = brms:::posterior_predict_discrete_weibull,
    prep_builder = function(ns = 25, nobs = 4, ...) {
      make_prep_count(
        ns = ns, nobs = nobs, mu = 0.7, shape = 1.5,
        Y = rep(1L, nobs), ...
      )
    },
    q_ref = 1,
    flags = list(truncation = TRUE, skip_moments = TRUE)
  ))

  # ---- Hurdle / ZI ----

  dist_registry_add(dist_entry(
    name = "hurdle_negbinomial",
    backend = "hurdle_negbinomial",
    type = "discrete",
    d = dhurdle_negbinomial,
    p = phurdle_negbinomial,
    q = qhurdle_negbinomial,
    r = NULL,
    params = list(mu = 4, shape = 2.5, hu = 0.2),
    support = c(0, Inf),
    pp_fun = brms:::posterior_predict_hurdle_negbinomial,
    prep_builder = function(ns = 25, nobs = 4, ...) {
      make_prep_count(
        ns = ns, nobs = nobs, mu = 4, shape = 2.5, hu = 0.2,
        Y = rep(3L, nobs), ...
      )
    },
    q_ref = 3,
    baseline = list(
      d = stats::dnbinom,
      p = stats::pnbinom,
      params = list(mu = 4, size = 2.5),
      mix = "hu"
    ),
    flags = list(
      truncation = TRUE,
      hurdle = TRUE,
      no_random_truncation = TRUE,
      skip_moments = TRUE
    )
  ))

  dist_registry_add(dist_entry(
    name = "hurdle_gamma",
    backend = "hurdle_gamma",
    type = "mixture",
    d = dhurdle_gamma,
    p = phurdle_gamma,
    q = qhurdle_gamma,
    r = NULL,
    params = list(shape = 2, scale = 1, hu = 0.2),
    support = c(0, Inf),
    pp_fun = brms:::posterior_predict_hurdle_gamma,
    prep_builder = function(ns = 25, nobs = 4, ...) {
      make_prep_positive(
        ns = ns, nobs = nobs, mu = 2, shape = 2,
        extra = list(hu = rep(0.2, ns)),
        Y = rep(2, nobs), ...
      )
    },
    q_ref = 1.5,
    baseline = list(
      d = stats::dgamma,
      p = stats::pgamma,
      params = list(shape = 2, scale = 1),
      mix = "hu"
    ),
    flags = list(
      truncation = FALSE,
      hurdle = TRUE,
      no_random_truncation = TRUE,
      skip_integrate = TRUE,
      skip_moments = TRUE,
      skip_pdf_fd = TRUE,
      has_atoms = TRUE
    )
  ))

  dist_registry_add(dist_entry(
    name = "hurdle_lognormal",
    backend = "hurdle_lognormal",
    type = "mixture",
    d = dhurdle_lognormal,
    p = phurdle_lognormal,
    q = qhurdle_lognormal,
    r = NULL,
    params = list(mu = 0, sigma = 1, hu = 0.2),
    support = c(0, Inf),
    pp_fun = brms:::posterior_predict_hurdle_lognormal,
    prep_builder = function(ns = 25, nobs = 4, ...) {
      make_prep_location(
        ns = ns, nobs = nobs, mu = 0, sigma = 1,
        extra = list(hu = rep(0.2, ns)),
        Y = rep(1, nobs), ...
      )
    },
    q_ref = 1,
    baseline = list(
      d = stats::dlnorm,
      p = stats::plnorm,
      params = list(meanlog = 0, sdlog = 1),
      mix = "hu"
    ),
    flags = list(
      truncation = FALSE,
      hurdle = TRUE,
      no_random_truncation = TRUE,
      skip_integrate = TRUE,
      skip_moments = TRUE,
      skip_pdf_fd = TRUE,
      has_atoms = TRUE
    )
  ))

  dist_registry_add(dist_entry(
    name = "zero_inflated_poisson",
    backend = "zero_inflated_poisson",
    type = "discrete",
    d = dzero_inflated_poisson,
    p = pzero_inflated_poisson,
    q = qzero_inflated_poisson,
    r = NULL,
    params = list(lambda = 3, zi = 0.25),
    support = c(0, Inf),
    pp_fun = brms:::posterior_predict_zero_inflated_poisson,
    prep_builder = function(ns = 25, nobs = 4, ...) {
      make_prep_count(
        ns = ns, nobs = nobs, mu = 3, zi = 0.25,
        Y = rep(2L, nobs), ...
      )
    },
    q_ref = 2,
    baseline = list(
      d = stats::dpois,
      p = stats::ppois,
      params = list(lambda = 3),
      mix = "zi"
    ),
    flags = list(
      truncation = TRUE,
      zi = TRUE,
      skip_moments = TRUE
    )
  ))

  dist_registry_add(dist_entry(
    name = "zero_inflated_binomial",
    backend = "zero_inflated_binomial",
    type = "discrete",
    d = dzero_inflated_binomial,
    p = pzero_inflated_binomial,
    q = qzero_inflated_binomial,
    r = NULL,
    params = list(size = 10, prob = 0.4, zi = 0.2),
    support = c(0, 10),
    pp_fun = brms:::posterior_predict_zero_inflated_binomial,
    prep_builder = function(ns = 25, nobs = 4, ...) {
      make_prep_count(
        ns = ns, nobs = nobs, mu = 0.4, trials = 10, zi = 0.2,
        Y = rep(4L, nobs), ...
      )
    },
    q_ref = 3,
    baseline = list(
      d = stats::dbinom,
      p = stats::pbinom,
      params = list(size = 10, prob = 0.4),
      mix = "zi"
    ),
    flags = list(
      truncation = TRUE,
      zi = TRUE,
      skip_moments = TRUE
    )
  ))

  dist_registry_add(dist_entry(
    name = "zero_inflated_beta_binomial",
    backend = "zero_inflated_beta_binomial",
    type = "discrete",
    d = dzero_inflated_beta_binomial,
    p = pzero_inflated_beta_binomial,
    q = qzero_inflated_beta_binomial,
    r = NULL,
    params = list(size = 12, mu = 0.35, phi = 8, zi = 0.2),
    support = c(0, 12),
    ref = list(requires = "extraDistr"),
    pp_fun = brms:::posterior_predict_zero_inflated_beta_binomial,
    prep_builder = function(ns = 25, nobs = 4, ...) {
      p <- make_prep_count(
        ns = ns, nobs = nobs, mu = 0.35, trials = 12, zi = 0.2,
        Y = rep(4L, nobs), ...
      )
      p$dpars$phi <- rep(8, p$ndraws)
      p
    },
    q_ref = 3,
    baseline = list(
      d = brms:::dbeta_binomial,
      p = brms:::pbeta_binomial,
      params = list(size = 12, mu = 0.35, phi = 8),
      mix = "zi"
    ),
    flags = list(
      truncation = TRUE,
      zi = TRUE,
      skip_moments = TRUE
    )
  ))

  dist_registry_add(dist_entry(
    name = "zero_inflated_beta",
    backend = "zero_inflated_beta",
    type = "mixture",
    d = dzero_inflated_beta,
    p = pzero_inflated_beta,
    q = qzero_inflated_beta,
    r = NULL,
    params = list(shape1 = 2, shape2 = 3, zi = 0.25),
    support = c(0, 1),
    pp_fun = brms:::posterior_predict_zero_inflated_beta,
    prep_builder = function(ns = 25, nobs = 4, ...) {
      p <- make_prep_beta_mix(
        ns = ns, nobs = nobs, mu = 0.4, phi = 5,
        zoi = 0, coi = 0, Y = rep(0.4, nobs), ...
      )
      p$dpars$zi <- rep(0.25, p$ndraws)
      p
    },
    q_ref = 0.4,
    flags = list(
      truncation = FALSE,
      zi = TRUE,
      skip_integrate = TRUE,
      skip_moments = TRUE,
      skip_pdf_fd = TRUE,
      has_atoms = TRUE
    )
  ))

  # ---- Ordinal variants ----

  register_ordinal_variant <- function(fam) {
    dist_registry_add(dist_entry(
      name = paste0("ordinal_", fam),
      backend = "ordinal",
      type = "ordinal",
      d = function(x, ...) {
        brms:::dordinal(x, family = fam, link = "logit", ...)
      },
      p = function(q, ...) {
        brms:::pordinal(q, family = fam, link = "logit", ...)
      },
      q = function(p, ...) {
        brms:::qordinal(p, family = fam, link = "logit", ...)
      },
      r = function(n, ...) {
        brms:::rordinal(n, family = fam, link = "logit", ...)
      },
      params = list(
        eta = 0,
        thres = matrix(c(-1, 0, 1), nrow = 1),
        disc = 1
      ),
      support = c(1, 4),
      pp_fun = brms:::posterior_predict_ordinal,
      prep_builder = function(...) {
        make_prep_ordinal(family = fam, ...)
      },
      q_ref = 2,
      flags = list(
        truncation = FALSE,
        skip_integrate = TRUE,
        skip_moments = TRUE,
        skip_pdf_fd = TRUE,
        skip_rng_cdf = TRUE,
        skip_p_tail_flags = TRUE
      )
    ))
  }
  for (ord_fam in c("sratio", "cratio", "acat")) {
    register_ordinal_variant(ord_fam)
  }

  # ---- Stubs (no deep d/p/q) ----

  stub_names_random <- c(
    "gaussian_mv", "student_mv",
    "gaussian_time", "student_time",
    "gaussian_lagsar", "student_lagsar",
    "gaussian_errorsar", "student_errorsar",
    "gaussian_fcor", "student_fcor",
    "multinomial", "dirichlet_multinomial",
    "dirichlet", "dirichlet2", "logistic_normal",
    "mixture"
  )
  for (nm in stub_names_random) {
    dist_registry_add(dist_entry(
      name = nm,
      backend = nm,
      type = "continuous",
      d = NULL,
      p = NULL,
      q = NULL,
      r = NULL,
      params = list(),
      support = c(-Inf, Inf),
      pp_fun = NULL,
      prep_builder = NULL,
      outputs = c("random"),
      flags = list(
        stub = TRUE,
        truncation = FALSE,
        skip_integrate = TRUE,
        skip_moments = TRUE,
        skip_rng_cdf = TRUE,
        skip_pdf_fd = TRUE
      )
    ))
  }

  dist_registry_add(dist_entry(
    name = "cox",
    backend = "cox",
    type = "continuous",
    d = NULL,
    p = NULL,
    q = NULL,
    r = NULL,
    params = list(),
    support = c(0, Inf),
    pp_fun = NULL,
    prep_builder = NULL,
    outputs = character(0),
    flags = list(
      stub = TRUE,
      truncation = FALSE,
      skip_integrate = TRUE,
      skip_moments = TRUE,
      skip_rng_cdf = TRUE,
      skip_pdf_fd = TRUE
    )
  ))

  invisible(.dist_registry_env$entries)
}

# Ensure registry is available when helpers are sourced
dist_registry_populate()

# ---------------------------------------------------------------------------
# Custom expectations: d/p/q/r relations
# ---------------------------------------------------------------------------

expect_pq_inverse <- function(entry, p = c(0.01, 0.1, 0.5, 0.85, 0.99),
                              tol = 1e-6) {
  testthat::expect_true(has_fun(entry, "p") && has_fun(entry, "q"))
  # Ordinal/categorical q*() treat length(p) as ndraws when eta is scalar;
  # evaluate each probability separately. Mixture atoms make p(q(p))=p fail.
  elementwise <- entry$type %in% c("ordinal", "mixture") ||
    isTRUE(entry$flags$pq_elementwise)
  use_discrete_def <- isTRUE(entry$flags$discrete_support) ||
    entry$type == "mixture" ||
    isTRUE(entry$flags$has_atoms)

  if (elementwise) {
    q <- vapply(p, function(pj) {
      as.numeric(call_dist(entry$q, pj, entry$params))[1]
    }, numeric(1))
    F_q <- vapply(q, function(qj) {
      as.numeric(call_dist(entry$p, qj, entry$params))[1]
    }, numeric(1))
    F_qm1 <- vapply(q, function(qj) {
      as.numeric(call_dist(entry$p, qj - 1, entry$params))[1]
    }, numeric(1))
  } else {
    q <- as.numeric(call_dist(entry$q, p, entry$params))
    F_q <- as.numeric(call_dist(entry$p, q, entry$params))
    F_qm1 <- as.numeric(call_dist(entry$p, q - 1, entry$params))
  }

  if (use_discrete_def) {
    testthat::expect_true(all(F_q + tol >= p),
                          info = paste(entry$name, "F(q) >= p"))
    if (isTRUE(entry$flags$discrete_support)) {
      ok <- ifelse(q > entry$support[1], F_qm1 < p + tol, TRUE)
      testthat::expect_true(all(ok),
                            info = paste(entry$name, "F(q-1) < p"))
    }
  } else {
    testthat::expect_equal(F_q, p, tolerance = tol,
                           info = paste(entry$name, "p(q(p))"))
  }
  invisible(TRUE)
}

expect_qp_inverse <- function(entry, x = NULL, tol = 1e-5) {
  testthat::expect_true(has_fun(entry, "p") && has_fun(entry, "q"))
  if (isTRUE(entry$flags$discrete_support) ||
      entry$type == "mixture" ||
      isTRUE(entry$flags$has_atoms)) {
    return(invisible(TRUE))
  }
  if (is.null(x)) {
    lo <- entry$support[1]
    hi <- entry$support[2]
    if (!is.finite(lo)) lo <- call_dist(entry$q, 0.05, entry$params)
    if (!is.finite(hi)) hi <- call_dist(entry$q, 0.95, entry$params)
    x <- as.numeric(seq(lo, hi, length.out = 7))
  }
  # avoid endpoints where q(p) is unstable for numeric inverses
  p <- call_dist(entry$p, x, entry$params)
  keep <- is.finite(p) & p > 1e-4 & p < 1 - 1e-4
  if (!any(keep)) return(invisible(TRUE))
  x <- x[keep]
  p <- p[keep]
  x2 <- call_dist(entry$q, p, entry$params)
  testthat::expect_equal(as.numeric(x2), as.numeric(x), tolerance = tol,
                         info = paste(entry$name, "q(p(x))"))
  invisible(TRUE)
}

expect_d_sums_to_p <- function(entry, k = NULL, tol = 1e-8) {
  testthat::expect_true(has_fun(entry, "d") && has_fun(entry, "p"))
  testthat::expect_true(isTRUE(entry$flags$discrete_support))
  if (isTRUE(entry$flags$skip_d_sums)) return(invisible(TRUE))
  if (is.null(k)) {
    hi <- entry$support[2]
    if (!is.finite(hi)) {
      hi <- as.integer(call_dist(entry$q, 0.999, entry$params)) + 2L
    }
    lo <- max(0L, as.integer(entry$support[1]))
    k <- lo:hi
  }
  dens <- as.numeric(call_dist(entry$d, k, entry$params))
  if (isTRUE(entry$flags$pq_elementwise)) {
    cdf <- vapply(k, function(kj) {
      as.numeric(call_dist(entry$p, kj, entry$params))[1]
    }, numeric(1))
  } else {
    cdf <- as.numeric(call_dist(entry$p, k, entry$params))
  }
  # F(k) = sum_{j<=k} f(j); compare consecutive differences
  pmf_from_cdf <- c(cdf[1], diff(cdf))
  # for support starting above min(k), first mass may include left tail
  if (k[1] > entry$support[1] || k[1] > 0) {
    F_prev <- as.numeric(call_dist(entry$p, k[1] - 1, entry$params))[1]
    pmf_from_cdf[1] <- cdf[1] - F_prev
  }
  testthat::expect_equal(dens, pmf_from_cdf, tolerance = tol,
                         info = paste(entry$name, "d == diff(p)"))
  invisible(TRUE)
}

expect_pdf_matches_cdf_fd <- function(entry, x = NULL, eps = 1e-5,
                                      tol = 5e-3) {
  testthat::expect_true(has_fun(entry, "d") && has_fun(entry, "p"))
  if (isTRUE(entry$flags$discrete_support) || isTRUE(entry$flags$skip_pdf_fd)) {
    return(invisible(TRUE))
  }
  if (is.null(x)) {
    qs <- call_dist(entry$q, c(0.2, 0.4, 0.6, 0.8), entry$params)
    x <- as.numeric(qs)
  }
  dens <- as.numeric(call_dist(entry$d, x, entry$params))
  fd <- (as.numeric(call_dist(entry$p, x + eps, entry$params)) -
           as.numeric(call_dist(entry$p, x - eps, entry$params))) / (2 * eps)
  testthat::expect_equal(dens, fd, tolerance = tol,
                         info = paste(entry$name, "d ≈ F'"))
  invisible(TRUE)
}

expect_integrates_to_one <- function(entry, tol = 1e-3) {
  testthat::expect_true(has_fun(entry, "d"))
  if (isTRUE(entry$flags$skip_integrate)) return(invisible(TRUE))
  if (isTRUE(entry$flags$discrete_support)) {
    hi <- entry$support[2]
    if (!is.finite(hi)) {
      hi <- as.integer(call_dist(entry$q, 0.9999, entry$params)) + 20L
    }
    lo <- max(0L, as.integer(entry$support[1]))
    s <- sum(as.numeric(call_dist(entry$d, lo:hi, entry$params)))
    testthat::expect_equal(s, 1, tolerance = tol,
                           info = paste(entry$name, "sum d = 1"))
  } else {
    lo <- entry$support[1]
    hi <- entry$support[2]
    if (!is.finite(lo)) lo <- as.numeric(call_dist(entry$q, 1e-6, entry$params)) - 1
    if (!is.finite(hi)) hi <- as.numeric(call_dist(entry$q, 1 - 1e-6, entry$params)) + 1
    integrand <- function(x) {
      as.numeric(call_dist(entry$d, x, entry$params))
    }
    s <- stats::integrate(integrand, lo, hi, rel.tol = 1e-4)$value
    testthat::expect_equal(s, 1, tolerance = tol,
                           info = paste(entry$name, "∫ d = 1"))
  }
  invisible(TRUE)
}

expect_log_and_tail_flags <- function(entry, x = NULL, p = 0.2) {
  if (is.null(x)) {
    x <- if (isTRUE(entry$flags$discrete_support)) {
      entry$q_ref
    } else {
      as.numeric(call_dist(entry$q, 0.4, entry$params))
    }
  }
  if (has_fun(entry, "d")) {
    d0 <- as.numeric(call_dist(entry$d, x, entry$params, log = FALSE))
    d1 <- as.numeric(call_dist(entry$d, x, entry$params, log = TRUE))
    testthat::expect_equal(d1, log(d0), tolerance = 1e-8,
                           info = paste(entry$name, "d log"))
  }
  if (has_fun(entry, "p") && !isTRUE(entry$flags$skip_p_tail_flags)) {
    pl <- as.numeric(call_dist(entry$p, x, entry$params,
                               lower.tail = TRUE, log.p = FALSE))
    pu <- as.numeric(call_dist(entry$p, x, entry$params,
                               lower.tail = FALSE, log.p = FALSE))
    testthat::expect_equal(pl + pu, 1, tolerance = 1e-8,
                           info = paste(entry$name, "p tails"))
    pll <- as.numeric(call_dist(entry$p, x, entry$params,
                                lower.tail = TRUE, log.p = TRUE))
    testthat::expect_equal(pll, log(pl), tolerance = 1e-8,
                           info = paste(entry$name, "p log.p"))
  }
  if (has_fun(entry, "q")) {
    q_up <- as.numeric(call_dist(entry$q, p, entry$params,
                                 lower.tail = FALSE, log.p = FALSE))
    q_lo <- as.numeric(call_dist(entry$q, 1 - p, entry$params,
                                 lower.tail = TRUE, log.p = FALSE))
    testthat::expect_equal(q_up, q_lo, tolerance = 1e-6,
                           info = paste(entry$name, "q lower.tail"))
    q_log <- as.numeric(call_dist(entry$q, log(p), entry$params,
                                  lower.tail = TRUE, log.p = TRUE))
    q_ref <- as.numeric(call_dist(entry$q, p, entry$params,
                                  lower.tail = TRUE, log.p = FALSE))
    testthat::expect_equal(q_log, q_ref, tolerance = 1e-6,
                           info = paste(entry$name, "q log.p"))
  }
  invisible(TRUE)
}

expect_moments_match <- function(entry, n = 8000, tol = 0.15, seed = 1) {
  if (!has_fun(entry, "r") || isTRUE(entry$flags$skip_moments)) {
    return(invisible(TRUE))
  }
  withr::local_seed(seed)
  draws <- as.numeric(do.call(entry$r, c(list(n), entry$params)))
  # approximate mean via quantile grid when closed form unknown
  ps <- seq(0.001, 0.999, length.out = 500)
  qs <- as.numeric(call_dist(entry$q, ps, entry$params))
  mean_ref <- mean(qs)
  testthat::expect_equal(mean(draws), mean_ref, tolerance = tol,
                         info = paste(entry$name, "E[X]"))
  invisible(TRUE)
}

expect_rng_fits_cdf <- function(entry, n = 4000, tol = 0.08, seed = 2) {
  if (!has_fun(entry, "r") || !has_fun(entry, "p") ||
      isTRUE(entry$flags$skip_rng_cdf)) {
    return(invisible(TRUE))
  }
  withr::local_seed(seed)
  draws <- as.numeric(do.call(entry$r, c(list(n), entry$params)))
  u <- as.numeric(call_dist(entry$p, draws, entry$params))
  # continuous: PIT ~ U(0,1); discrete: not uniform — compare ECDF at probes
  if (isTRUE(entry$flags$discrete_support)) {
    probes <- as.numeric(call_dist(entry$q, c(0.25, 0.5, 0.75), entry$params))
    for (qk in unique(probes)) {
      emp <- mean(draws <= qk)
      th <- as.numeric(call_dist(entry$p, qk, entry$params))
      testthat::expect_equal(emp, th, tolerance = tol,
                             info = paste(entry$name, "ECDF at", qk))
    }
  } else {
    testthat::expect_equal(mean(u), 0.5, tolerance = tol,
                           info = paste(entry$name, "mean PIT"))
    testthat::expect_true(
      abs(mean(u < 0.5) - 0.5) < tol,
      info = paste(entry$name, "P(PIT<0.5)")
    )
  }
  invisible(TRUE)
}

expect_truncated_cdf_density_quantile <- function(entry, lb, ub,
                                                   q = NULL, p = 0.4,
                                                   tol = 1e-7) {
  if (!has_fun(entry, "p") || !isTRUE(entry$flags$truncation)) {
    return(invisible(TRUE))
  }
  if (is.null(q)) q <- entry$q_ref
  F <- function(x) as.numeric(call_dist(entry$p, x, entry$params))
  denom <- F(ub) - F(lb)
  testthat::expect_true(denom > 0, info = paste(entry$name, "trunc mass"))

  # truncated CDF (cheatsheet)
  F_tr <- (F(q) - F(lb)) / denom
  F_tr_pp <- do.call(
    brms:::compute_cdf,
    c(list(q = q, distribution = entry$backend, lb = lb, ub = ub,
           randomized = FALSE), entry$params)
  )
  testthat::expect_equal(as.numeric(F_tr_pp), F_tr, tolerance = tol,
                         info = paste(entry$name, "trunc CDF"))

  if (has_fun(entry, "d")) {
    dens <- as.numeric(call_dist(entry$d, q, entry$params))
    dens_tr <- if (q < lb || q > ub) 0 else dens / denom
    dens_pp <- do.call(
      brms:::compute_density,
      c(list(q = q, distribution = entry$backend, lb = lb, ub = ub),
        entry$params)
    )
    testthat::expect_equal(as.numeric(dens_pp), dens_tr, tolerance = tol,
                           info = paste(entry$name, "trunc density"))
  }

  if (has_fun(entry, "q")) {
    p_star <- p * denom + F(lb)
    q_tr <- as.numeric(call_dist(entry$q, p_star, entry$params))
    q_pp <- do.call(
      brms:::compute_quantile,
      c(list(p = p, distribution = entry$backend, lb = lb, ub = ub),
        entry$params)
    )
    testthat::expect_equal(as.numeric(q_pp), q_tr, tolerance = 1e-5,
                           info = paste(entry$name, "trunc quantile"))
  }
  invisible(TRUE)
}

expect_zi_hurdle_mixture <- function(entry, y = 0:8, tol = 1e-8) {
  if (is.null(entry$baseline)) return(invisible(TRUE))
  base <- entry$baseline
  pi <- entry$params[[base$mix]]
  g_d <- function(x) as.numeric(call_dist(base$d, x, base$params))
  g_p <- function(x) as.numeric(call_dist(base$p, x, base$params))

  if (isTRUE(entry$flags$zi)) {
    # P(Y=0) = pi + (1-pi) G(0); P(Y=y>0) = (1-pi) g(y)
    # For PMF comparison use d; for y=0 use point mass formula
    d0 <- as.numeric(call_dist(entry$d, 0, entry$params))
    testthat::expect_equal(
      d0, pi + (1 - pi) * g_d(0), tolerance = tol,
      info = paste(entry$name, "ZI d(0)")
    )
    ypos <- y[y > 0]
    if (length(ypos)) {
      d_y <- as.numeric(call_dist(entry$d, ypos, entry$params))
      testthat::expect_equal(
        d_y, (1 - pi) * g_d(ypos), tolerance = tol,
        info = paste(entry$name, "ZI d(y>0)")
      )
    }
  }
  if (isTRUE(entry$flags$hurdle)) {
    # P(Y=0)=h; P(Y=y>0)=(1-h) g(y)/(1-G(0))
    d0 <- as.numeric(call_dist(entry$d, 0, entry$params))
    testthat::expect_equal(d0, pi, tolerance = tol,
                           info = paste(entry$name, "hurdle d(0)"))
    ypos <- y[y > 0]
    if (length(ypos)) {
      d_y <- as.numeric(call_dist(entry$d, ypos, entry$params))
      testthat::expect_equal(
        d_y, (1 - pi) * g_d(ypos) / (1 - g_p(0)), tolerance = tol,
        info = paste(entry$name, "hurdle d(y>0)")
      )
    }
  }
  invisible(TRUE)
}
