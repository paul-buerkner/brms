context("Distribution validation: d/p/q/r relations")

test_that("p/q inverses hold for families with p and q", {
  for (entry in dist_active_entries(require_funs = c("p", "q"))) {
    expect_pq_inverse(entry)
    expect_qp_inverse(entry)
  }
})

test_that("discrete PMF matches CDF differences", {
  for (entry in dist_active_entries(
    require_funs = c("d", "p"),
    discrete_support = TRUE
  )) {
    expect_d_sums_to_p(entry)
  }
})

test_that("continuous PDF matches finite-difference of CDF", {
  for (entry in dist_active_entries(type = c("continuous", "circular"))) {
    expect_pdf_matches_cdf_fd(entry)
  }
})

test_that("densities integrate or sum to one", {
  for (entry in dist_active_entries(require_funs = "d")) {
    expect_integrates_to_one(entry)
  }
})

test_that("log and lower.tail / log.p flags are consistent", {
  for (entry in dist_active_entries()) {
    expect_log_and_tail_flags(entry)
  }
})

test_that("RNG moments and CDF fit for families with r()", {
  skip_on_cran()
  for (entry in dist_active_entries()) {
    expect_moments_match(entry)
    expect_rng_fits_cdf(entry)
  }
})

test_that("truncation formulas match compute_* helpers", {
  for (entry in dist_active_entries(truncation = TRUE)) {
    if (isTRUE(entry$flags$discrete_support)) {
      expect_truncated_cdf_density_quantile(entry, lb = 1, ub = 6)
      # a non-integer bound admits the same support as ceiling(lb), which
      # separates ceiling(lb) - 1 from a plain lb - 1. See #1923
      expect_truncated_cdf_density_quantile(entry, lb = 1.5, ub = 6)
    } else {
      qs <- as.numeric(call_dist(entry$q, c(0.2, 0.8), entry$params))
      expect_truncated_cdf_density_quantile(entry, lb = qs[1], ub = qs[2])
    }
  }
})

test_that("zero-inflated and hurdle mixture formulas match baseline d/p", {
  for (entry in dist_active_entries(has_baseline = TRUE)) {
    expect_zi_hurdle_mixture(entry)
  }
})
