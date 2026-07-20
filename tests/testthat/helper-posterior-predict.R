# posterior_predict contract helpers
# Depends on helper-distributions.R (sourced first alphabetically).

# Registry entries suitable for posterior_predict contract tests
pp_test_entries <- function(truncation = FALSE) {
  dist_active_entries(
    require_funs = "pp_fun",
    truncation = if (truncation) TRUE else NULL
  )
}

skip_if_entry_missing_deps <- function(entry) {
  for (pkg in entry_requires(entry)) {
    testthat::skip_if_not_installed(pkg)
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# posterior_predict expectations
# ---------------------------------------------------------------------------

expect_pp_output_matches_dist <- function(entry, i = 1L, tol = 1e-7,
                                          seed = 42) {
  testthat::expect_true(has_fun(entry, "pp_fun"))
  prep <- entry$prep_builder(ns = 15, nobs = max(3, i), seed = seed)
  q <- entry$q_ref
  p <- entry$p_ref
  pp <- entry$pp_fun

  if (has_output(entry, "probability") && has_fun(entry, "p")) {
    got <- pp(i, prep = prep, output = "probability", q = q)
    # constant dpars => all draws equal scalar CDF
    exp_val <- as.numeric(call_dist(entry$p, q, entry$params))
    # ordinal/categorical may return length-1 or matrix; coerce
    if (length(exp_val) > 1L) exp_val <- exp_val[1]
    testthat::expect_equal(
      as.numeric(got), rep(as.numeric(exp_val), prep$ndraws),
      tolerance = tol,
      info = paste(entry$name, "PP probability")
    )
  }

  if (has_output(entry, "density") && has_fun(entry, "d")) {
    got <- pp(i, prep = prep, output = "density", q = q)
    exp_val <- as.numeric(call_dist(entry$d, q, entry$params))
    if (length(exp_val) > 1L) exp_val <- exp_val[1]
    testthat::expect_equal(
      as.numeric(got), rep(as.numeric(exp_val), prep$ndraws),
      tolerance = tol,
      info = paste(entry$name, "PP density")
    )
  }

  if (has_output(entry, "quantile") && has_fun(entry, "q")) {
    got <- pp(i, prep = prep, output = "quantile", p = p)
    exp_val <- as.numeric(call_dist(entry$q, p, entry$params))
    if (length(exp_val) > 1L) exp_val <- exp_val[1]
    testthat::expect_equal(
      as.numeric(got), rep(as.numeric(exp_val), prep$ndraws),
      tolerance = 1e-5,
      info = paste(entry$name, "PP quantile")
    )
  }

  if (has_output(entry, "random")) {
    set.seed(seed)
    got <- pp(i, prep = prep, output = "random")
    testthat::expect_length(got, prep$ndraws)
    testthat::expect_true(all(is.finite(got)))
    lo <- entry$support[1]
    hi <- entry$support[2]
    if (is.finite(lo)) testthat::expect_true(all(got >= lo - 1e-8))
    if (is.finite(hi)) testthat::expect_true(all(got <= hi + 1e-8))
  }
  invisible(TRUE)
}

# Map prep dpars to compute_cdf / backend args for randomized PIT checks.
# Keyed by registry entry$name; ordinal_* share the "ordinal" builder.
.pp_dist_args_from_prep <- list(
  ordinal = function(prep, i) {
    list(
      eta = prep$dpars$mu[, i],
      thres = prep$thres$thres,
      disc = prep$dpars$disc,
      family = prep$family$family,
      link = prep$family$link
    )
  },
  categorical = function(prep, i) {
    Mu <- brms:::get_Mu(prep, i = i)
    list(eta = brms:::insert_refcat(Mu, refcat = prep$refcat))
  },
  hurdle_cumulative = function(prep, i) {
    list(
      eta = prep$dpars$mu[, i],
      thres = prep$thres$thres,
      hu = prep$dpars$hu,
      disc = prep$dpars$disc,
      link = "logit"
    )
  },
  pois = function(prep, i) {
    list(lambda = prep$dpars$mu[, i])
  },
  binom = function(prep, i) {
    list(size = prep$data$trials[i], prob = prep$dpars$mu[, i])
  },
  bernoulli = function(prep, i) {
    list(size = prep$data$trials[i], prob = prep$dpars$mu[, i])
  },
  beta_binomial = function(prep, i) {
    list(
      size = prep$data$trials[i],
      mu = prep$dpars$mu[, i],
      phi = prep$dpars$phi
    )
  },
  nbinom = function(prep, i) {
    list(mu = prep$dpars$mu[, i], size = prep$dpars$shape)
  },
  negbinomial2 = function(prep, i) {
    list(mu = prep$dpars$mu[, i], size = 1 / prep$dpars$sigma)
  },
  geometric = function(prep, i) {
    list(mu = prep$dpars$mu[, i], size = rep(1, prep$ndraws))
  },
  com_poisson = function(prep, i) {
    list(mu = prep$dpars$mu[, i], shape = prep$dpars$shape)
  },
  discrete_weibull = function(prep, i) {
    list(mu = prep$dpars$mu[, i], shape = prep$dpars$shape)
  },
  hurdle_poisson = function(prep, i) {
    list(lambda = prep$dpars$mu[, i], hu = prep$dpars$hu)
  },
  hurdle_negbinomial = function(prep, i) {
    list(
      mu = prep$dpars$mu[, i],
      shape = prep$dpars$shape,
      hu = prep$dpars$hu
    )
  },
  zero_inflated_negbinomial = function(prep, i) {
    list(
      mu = prep$dpars$mu[, i],
      shape = prep$dpars$shape,
      zi = prep$dpars$zi
    )
  },
  zero_inflated_poisson = function(prep, i) {
    list(lambda = prep$dpars$mu[, i], zi = prep$dpars$zi)
  },
  zero_inflated_binomial = function(prep, i) {
    list(
      size = prep$data$trials[i],
      prob = prep$dpars$mu[, i],
      zi = prep$dpars$zi
    )
  },
  zero_inflated_beta_binomial = function(prep, i) {
    list(
      size = prep$data$trials[i],
      mu = prep$dpars$mu[, i],
      phi = prep$dpars$phi,
      zi = prep$dpars$zi
    )
  }
)

expand_scalar_params <- function(params, ndraws) {
  lapply(params, function(a) {
    if (is.atomic(a) && !is.null(a) && length(a) == 1L) {
      rep(a, ndraws)
    } else {
      a
    }
  })
}

get_pp_dist_args <- function(entry, prep, i) {
  key <- entry$name
  if (entry$backend == "ordinal" || grepl("^ordinal_", key)) {
    key <- "ordinal"
  }
  fun <- .pp_dist_args_from_prep[[key]]
  if (is.null(fun)) {
    expand_scalar_params(entry$params, prep$ndraws)
  } else {
    fun(prep, i)
  }
}

expect_pp_pit_contract <- function(entry, i = 1L, seed = 99, tol = 1e-8) {
  testthat::expect_true(has_fun(entry, "pp_fun"))
  if (!has_output(entry, "pit") || !has_output(entry, "probability")) {
    return(invisible(TRUE))
  }
  prep <- entry$prep_builder(ns = 40, nobs = max(3, i), seed = seed)
  q <- entry$q_ref
  pp <- entry$pp_fun

  if (entry$type %in% c("continuous", "circular", "mixture") &&
      !isTRUE(entry$flags$discrete_support)) {
    prob <- pp(i, prep = prep, output = "probability", q = q)
    pit <- pp(i, prep = prep, output = "pit", q = q)
    testthat::expect_equal(pit, prob, tolerance = tol,
                           info = paste(entry$name, "PIT ≡ probability"))
  } else {
    # discrete / ordinal: randomized PIT = F(q-1) + V*(F(q)-F(q-1))
    set.seed(seed)
    pit <- pp(i, prep = prep, output = "pit", q = q)
    set.seed(seed)
    args <- get_pp_dist_args(entry, prep, i)
    expected <- do.call(
      brms:::compute_cdf,
      c(list(q = q, distribution = entry$backend, lb = NULL, ub = NULL,
             randomized = TRUE), args)
    )
    testthat::expect_equal(pit, expected, tolerance = tol,
                           info = paste(entry$name, "randomized PIT"))
    Fq <- do.call(
      brms:::compute_cdf,
      c(list(q = q, distribution = entry$backend, lb = NULL, ub = NULL,
             randomized = FALSE), args)
    )
    Fqm1 <- do.call(
      brms:::compute_cdf,
      c(list(q = q - 1, distribution = entry$backend, lb = NULL, ub = NULL,
             randomized = FALSE), args)
    )
    testthat::expect_true(all(pit >= pmin(Fqm1, Fq) - 1e-10 &
                                pit <= pmax(Fqm1, Fq) + 1e-10),
                          info = paste(entry$name, "PIT in [F(q-1), F(q)]"))
  }
  invisible(TRUE)
}

expect_pp_truncation <- function(entry, i = 1L, lb = NULL, ub = NULL,
                                 tol = 1e-7, seed = 7) {
  if (!isTRUE(entry$flags$truncation) || !has_fun(entry, "pp_fun")) {
    return(invisible(TRUE))
  }
  if (is.null(lb) || is.null(ub)) {
    if (isTRUE(entry$flags$discrete_support)) {
      lb <- 1
      ub <- 6
    } else {
      qs <- as.numeric(call_dist(entry$q, c(0.15, 0.85), entry$params))
      lb <- qs[1]
      ub <- qs[2]
    }
  }
  prep <- entry$prep_builder(ns = 12, nobs = max(3, i), seed = seed)
  prep <- set_trunc_bounds(prep, lb = lb, ub = ub, i = i)
  q <- entry$q_ref
  # clamp q into [lb, ub] for density/prob checks
  q_in <- min(max(q, lb), ub)

  got_p <- entry$pp_fun(i, prep = prep, output = "probability", q = q_in)
  exp_p <- do.call(
    brms:::compute_cdf,
    c(list(q = q_in, distribution = entry$backend, lb = lb, ub = ub,
           randomized = FALSE), entry$params)
  )
  testthat::expect_equal(got_p, rep(as.numeric(exp_p), prep$ndraws),
                         tolerance = tol,
                         info = paste(entry$name, "PP trunc probability"))

  if (has_output(entry, "density")) {
    got_d <- entry$pp_fun(i, prep = prep, output = "density", q = q_in)
    exp_d <- do.call(
      brms:::compute_density,
      c(list(q = q_in, distribution = entry$backend, lb = lb, ub = ub),
        entry$params)
    )
    testthat::expect_equal(got_d, rep(as.numeric(exp_d), prep$ndraws),
                           tolerance = tol,
                           info = paste(entry$name, "PP trunc density"))
  }

  if (has_output(entry, "quantile")) {
    got_q <- entry$pp_fun(i, prep = prep, output = "quantile", p = 0.4)
    exp_q <- do.call(
      brms:::compute_quantile,
      c(list(p = 0.4, distribution = entry$backend, lb = lb, ub = ub),
        entry$params)
    )
    testthat::expect_equal(got_q, rep(as.numeric(exp_q), prep$ndraws),
                           tolerance = 1e-5,
                           info = paste(entry$name, "PP trunc quantile"))
  }

  if (has_output(entry, "random") &&
      !isTRUE(entry$flags$no_random_truncation)) {
    set.seed(seed)
    got_r <- entry$pp_fun(i, prep = prep, output = "random", ntrys = 20)
    testthat::expect_true(
      all(got_r >= lb - 1e-8 & got_r <= ub + 1e-8, na.rm = TRUE),
      info = paste(entry$name, "PP trunc random in [lb, ub]")
    )
  }
  invisible(TRUE)
}

expect_pp_log_tail_flags <- function(entry, i = 1L, seed = 11) {
  if (!has_fun(entry, "pp_fun")) return(invisible(TRUE))
  prep <- entry$prep_builder(ns = 8, nobs = max(3, i), seed = seed)
  q <- entry$q_ref
  p <- entry$p_ref
  pp <- entry$pp_fun

  if (has_output(entry, "probability")) {
    pl <- pp(i, prep = prep, output = "probability", q = q,
             lower.tail = TRUE, log.p = FALSE)
    pu <- pp(i, prep = prep, output = "probability", q = q,
             lower.tail = FALSE, log.p = FALSE)
    testthat::expect_equal(pl + pu, rep(1, prep$ndraws), tolerance = 1e-7,
                           info = paste(entry$name, "PP p tails"))
    pll <- pp(i, prep = prep, output = "probability", q = q,
              lower.tail = TRUE, log.p = TRUE)
    testthat::expect_equal(pll, log(pl), tolerance = 1e-7,
                           info = paste(entry$name, "PP p log.p"))
  }
  if (has_output(entry, "density")) {
    d0 <- pp(i, prep = prep, output = "density", q = q, log = FALSE)
    d1 <- pp(i, prep = prep, output = "density", q = q, log = TRUE)
    testthat::expect_equal(d1, log(d0), tolerance = 1e-7,
                           info = paste(entry$name, "PP d log"))
  }
  if (has_output(entry, "quantile")) {
    q_up <- pp(i, prep = prep, output = "quantile", p = p,
               lower.tail = FALSE)
    q_lo <- pp(i, prep = prep, output = "quantile", p = 1 - p,
               lower.tail = TRUE)
    testthat::expect_equal(q_up, q_lo, tolerance = 1e-5,
                           info = paste(entry$name, "PP q lower.tail"))
  }
  invisible(TRUE)
}
