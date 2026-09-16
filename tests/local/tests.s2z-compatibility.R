# Compare independently installed brms revisions, using a fresh R process for
# each library. This is deliberately separate from implicit/explicit s2z=FALSE
# checks: both sides must come from the specified, separately installed builds.
#
# Usage:
#   Rscript tests/local/tests.s2z-compatibility.R \
#     MODE BASELINE_LIB CANDIDATE_LIB [OUTPUT_DIR]
#
# MODE is conventional or s2z. For conventional models, use the agreed upstream
# revision plus any separately reviewed prerequisite fixes as BASELINE_LIB. For
# s2z, use the pre-refactoring foundation revision. The runner requires identical
# Stan source and Stan data, including attributes, and exits nonzero on
# an unexpected generation error or difference. S2Z comparisons allow only the
# documented deletion of the exact, unused chol2inv_brms definition; its body is
# checked verbatim; any reference outside the definition fails the comparison.
# Unsupported S2Z sparse/QR designs must retain their exact diagnostic messages.
# OUTPUT_DIR retains source, data, package
# provenance, logs, and a comparison table (a temporary directory by default).
#
# No Stan compiler, sampling backend, testthat, or internet access is required.
# Use the sampling-equivalence scripts for statistical validation of intentional
# changes to generated mathematics; do not relax this comparison to conceal one.

compatibility_cases <- function(mode) {
  i <- seq_len(48L)
  dat <- data.frame(
    y = sin(i / 3) + i / 100,
    y2 = cos(i / 5) - i / 80,
    binary = as.integer(i %% 3L == 0L),
    count = i %% 5L,
    x = (i - 20) / 24,
    z = cos(i / 7),
    w = sin(i / 11),
    g = factor(rep(letters[1:6], each = 8L)),
    h = factor(rep(LETTERS[1:4], length.out = 48L)),
    time = i,
    se = rep(0.2, 48L)
  )
  dat$by <- factor(rep(c("first", "second"), each = 24L))
  dat$pw <- rep(seq(0.8, 1.3, length.out = 6L), each = 8L)
  A <- outer(seq_len(6L), seq_len(6L), function(a, b) 0.6^abs(a - b))
  dimnames(A) <- list(levels(dat$g), levels(dat$g))
  W <- outer(i, i, function(a, b) as.numeric(abs(a - b) == 1L))
  # SAR uses a row-standardized matrix; CAR uses the symmetric adjacency matrix.
  W_sar <- W / rowSums(W)
  normal_prior <- prior(normal(0, 1.7), class = Intercept) +
    prior(normal(0, 1.2), class = b)
  if (identical(mode, "conventional")) {
    cases <- list(
      gaussian_correlated = list(object = y ~ x + (1 + x | g)),
      gaussian_independent = list(object = y ~ x + (1 + x || g)),
      student_group = list(
        object = y ~ x + (1 + x | gr(g, dist = "student"))
      ),
      student_family = list(
        object = y ~ x + (1 + x | g), family = student()
      ),
      known_covariance = list(
        object = y ~ x + (1 + x | gr(g, cov = A)), data2 = list(A = A)
      ),
      by_group = list(object = y ~ x + (1 + x | gr(g, by = by))),
      weighted_group = list(object = y ~ x + (1 + x | gr(g, pw = pw))),
      multivariate = list(
        object = bf(y ~ x + (1 + x | p | g)) +
          bf(y2 ~ z + (1 + z | p | g)) + set_rescor(TRUE)
      ),
      distributional = list(object = bf(y ~ x + (1 | g), sigma ~ z)),
      sparse = list(object = bf(y ~ x + z + (1 | g), sparse = TRUE)),
      qr = list(object = bf(y ~ x + z + (1 + x | g), decomp = "QR")),
      gp = list(object = y ~ x + gp(z, w, iso = FALSE)),
      ar_covariance = list(object = y ~ x + ar(time, cov = TRUE)),
      ar_se_covariance = list(object = y | se(se, sigma = TRUE) ~
        x + ar(time, cov = TRUE)),
      car = list(object = y ~ x + car(W), data2 = list(W = W)),
      sar_lag = list(object = y ~ x + sar(W_sar, type = "lag"),
        data2 = list(W_sar = W_sar)),
      sar_error = list(object = y ~ x + sar(W_sar, type = "error"),
        data2 = list(W_sar = W_sar)),
      constant_unnormalized = list(
        object = y ~ x + (1 + x | g), normalize = FALSE,
        prior = prior(constant(1), class = sd, group = g, coef = Intercept)
      ),
      glm = list(object = binary ~ x + z + (1 | g), family = bernoulli()),
      threaded_glm = list(object = count ~ x + (1 | g), family = poisson(),
        threads = threading(2L, grainsize = 8L, static = TRUE)),
      exact_prior_lookup = list(
        object = y ~ x + z + (1 + x | g),
        prior = prior(normal(0, 1), class = b) +
          prior(normal(3, 0.7), class = b, coef = x) +
          prior(constant(0.2), class = sd, group = g, coef = x)
      )
    )
  } else {
    cases <- list(
      scalar_default = list(object = y ~ 1 + (1 | gr(g, s2z = TRUE))),
      scalar_normal = list(object = y ~ 1 + (1 | gr(g, s2z = TRUE)),
        prior = prior(normal(0, 2), class = Intercept)),
      scalar_logistic = list(object = y ~ 1 + (1 | gr(g, s2z = TRUE)),
        prior = prior(logistic(0, 1), class = Intercept)),
      scalar_student_slope = list(object = y ~ 0 + x +
        (0 + x | gr(g, s2z = TRUE, dist = "student")),
        prior = prior(student_t(5, 0, 1), class = b)),
      correlated = list(object = y ~ x + (1 + x | gr(g, s2z = TRUE))),
      independent = list(object = y ~ x + z +
        (1 + x + z || gr(g, s2z = TRUE))),
      correlated_student = list(object = y ~ x +
        (1 + x | gr(g, s2z = TRUE, dist = "student")),
        prior = prior(student_t(7, 0, 2), class = b)),
      independent_student = list(object = y ~ x +
        (1 + x || gr(g, s2z = TRUE, dist = "student"))),
      multiblock_matheron = list(object = y ~ x +
        (1 + x | gr(g, s2z = TRUE)) + (1 + x | gr(h, s2z = TRUE)),
        prior = normal_prior),
      multiblock_independent = list(object = y ~ x + z +
        (1 + x || gr(g, s2z = TRUE)) + (0 + z || gr(h, s2z = TRUE)),
        prior = normal_prior),
      multiblock_dense = list(object = y ~ x + z +
        (1 + x | gr(g, s2z = TRUE, dist = "student")) +
        (1 + z | gr(h, s2z = TRUE))),
      multiblock_logistic = list(object = y ~ x +
        (1 + x | gr(g, s2z = TRUE)) + (1 | gr(h, s2z = TRUE)),
        prior = prior(logistic(0.2, 1.3), class = b)),
      mixed_ordinary_s2z = list(object = y ~ x + z +
        (1 + z | h) + (1 + x | gr(g, s2z = TRUE)), prior = normal_prior),
      original_prior_indices = list(object = y ~ x + z + w +
        (0 + z | gr(g, s2z = TRUE)),
        prior = prior(normal(0, 2), class = Intercept) +
          prior(double_exponential(0.4, 0.8), class = b, coef = x) +
          prior(student_t(6, -0.3, 1.4), class = b, coef = z) +
          prior(normal(2, 0.6), class = b, coef = w)),
      inactive_horseshoe = list(object = y ~ x + z +
        (1 | gr(g, s2z = TRUE)),
        prior = prior(horseshoe(), class = b) +
          prior(normal(0, 2), class = Intercept)),
      inactive_r2d2 = list(object = y ~ x + z + (1 | gr(g, s2z = TRUE)),
        prior = prior(R2D2(), class = b) +
          prior(normal(0, 2), class = Intercept)),
      nonconsecutive_active = list(object = y ~ x * z +
        (1 + x:z | gr(g, s2z = TRUE)) + (1 + x:z | gr(h, s2z = TRUE)),
        prior = normal_prior),
      sparse = list(object = bf(y ~ x + z +
        (1 + x | gr(g, s2z = TRUE)), sparse = TRUE),
        expected_error = "S2Z capability 'sparse'"),
      qr = list(object = bf(y ~ x + z +
        (1 + x | gr(g, s2z = TRUE)), decomp = "QR"),
        expected_error = "S2Z capability 'qr'"),
      unnormalized = list(object = y ~ x + (1 + x | gr(g, s2z = TRUE)),
        prior = normal_prior, normalize = FALSE),
      glm = list(object = binary ~ x + (1 + x | gr(g, s2z = TRUE)),
        family = bernoulli(), prior = normal_prior),
      threaded_glm = list(object = count ~ x + (1 + x | gr(g, s2z = TRUE)),
        family = poisson(), prior = normal_prior,
        threads = threading(2L, grainsize = 8L, static = TRUE)),
      distributional = list(object = bf(y ~ x + (1 | gr(g, s2z = TRUE)),
        sigma ~ z + (1 | gr(h, s2z = TRUE)))),
      multivariate = list(object = bf(y ~ x + (1 | gr(g, s2z = TRUE))) +
        bf(y2 ~ z + (1 | gr(h, s2z = TRUE))) + set_rescor(FALSE))
    )
  }
  lapply(cases, function(case) c(case, list(data = dat)))
}

generate_compatibility <- function(mode, lib, output) {
  lib <- normalizePath(lib, mustWork = TRUE)
  if (!dir.exists(file.path(lib, "brms"))) {
    stop("No installed brms in requested library: ", lib)
  }
  .libPaths(c(lib, .libPaths()))
  suppressPackageStartupMessages(library(brms, lib.loc = lib))
  loaded <- normalizePath(find.package("brms"), mustWork = TRUE)
  stopifnot(identical(loaded, normalizePath(file.path(lib, "brms"))))
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  provenance <- list(
    mode = mode, library = lib, package_path = loaded,
    description = utils::packageDescription("brms"),
    session = utils::sessionInfo()
  )
  saveRDS(provenance, file.path(output, "provenance.rds"))
  writeLines(capture.output(str(provenance)),
    file.path(output, "provenance.txt"))
  cases <- compatibility_cases(mode)
  status <- lapply(names(cases), function(name) {
    # Fixed seeds also stabilize approximate GP or other data preparation that
    # might use RNG in later revisions. Each fixture starts from the same seed.
    set.seed(1925L)
    case <- cases[[name]]
    expected <- case$expected_error
    case$expected_error <- NULL
    code <- tryCatch(do.call(stancode, case), error = identity)
    set.seed(1925L)
    data_case <- case[setdiff(names(case), "normalize")]
    sdata <- tryCatch(do.call(standata, data_case), error = identity)
    code_error <- inherits(code, "error")
    data_error <- inherits(sdata, "error")
    if (code_error) {
      writeLines(conditionMessage(code),
        file.path(output, paste0(name, ".error")))
    } else {
      writeLines(as.character(code), file.path(output, paste0(name, ".stan")))
    }
    if (data_error) {
      writeLines(conditionMessage(sdata),
        file.path(output, paste0(name, ".data-error")))
    } else {
      saveRDS(sdata, file.path(output, paste0(name, ".rds")))
    }
    code_ok <- if (is.null(expected)) !code_error else {
      code_error && grepl(expected, conditionMessage(code), fixed = TRUE)
    }
    data_ok <- if (is.null(expected)) !data_error else {
      data_error && grepl(expected, conditionMessage(sdata), fixed = TRUE)
    }
    data.frame(case = name, code_ok = code_ok, data_ok = data_ok,
      expected_error = !is.null(expected))
  })
  status <- do.call(rbind, status)
  utils::write.csv(status, file.path(output, "generation.csv"),
    row.names = FALSE)
  print(status, row.names = FALSE)
  if (!all(status$code_ok & status$data_ok)) {
    stop("Generation failed; inspect *.error and *.data-error in ", output)
  }
}

remove_unused_legacy_helper <- function(code) {
  # A literal historical definition, not a whitespace or Stan-code normalizer.
  # Check every occurrence so this exception can never mask a new caller.
  helper <- c(
    "",
    "  matrix chol2inv_brms(matrix L) {",
    "    int K = rows(L);",
    "    if (cols(L) != K) {",
    "      return chol2inv(L);",
    "    }",
    "    if (K == 1 && L[1, 1] != 0.0) {",
    "      return rep_matrix(inv_square(L[1, 1]), 1, 1);",
    "    }",
    "    if (K == 2 && L[1, 2] == 0.0 &&",
    "        L[1, 1] != 0.0 && L[2, 2] != 0.0) {",
    "      matrix[2, 2] precision;",
    "      real a = inv(L[1, 1]);",
    "      real d = inv(L[2, 2]);",
    "      real c = -L[2, 1] * a / L[2, 2];",
    "      precision[1, 1] = square(a) + square(c);",
    "      precision[1, 2] = c * d;",
    "      precision[2, 1] = precision[1, 2];",
    "      precision[2, 2] = square(d);",
    "      return precision;",
    "    }",
    "    return chol2inv(L);",
    "  }"
  )
  occurrences <- gregexpr(
    "\\bchol2inv_brms\\b", paste(code, collapse = "\n")
  )[[1L]]
  if (identical(as.integer(occurrences), -1L)) return(code)
  if (length(occurrences) != 1L) {
    stop("Cannot ignore chol2inv_brms: generated code contains a caller.")
  }
  start <- which(code == helper[[2L]]) - 1L
  index <- if (length(start) == 1L) {
    seq.int(start, length.out = length(helper))
  } else integer()
  if (!length(index) || !identical(code[index], helper)) {
    stop("Cannot ignore chol2inv_brms: definition differs from audited helper.")
  }
  code[-index]
}

compare_compatibility <- function(mode, baseline, candidate, output, script) {
  output <- normalizePath(output, mustWork = FALSE)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  for (label in c("baseline", "candidate")) {
    lib <- if (label == "baseline") baseline else candidate
    target <- file.path(output, label)
    log <- file.path(output, paste0(label, ".log"))
    args <- c(script, "--generate", mode, lib, target)
    result <- system2(file.path(R.home("bin"), "Rscript"), shQuote(args),
      stdout = log, stderr = log)
    if (result != 0L) {
      cat(readLines(log, warn = FALSE), sep = "\n")
      stop(label, " generation failed; see ", log)
    }
  }
  before <- utils::read.csv(file.path(output, "baseline", "generation.csv"))
  after <- utils::read.csv(file.path(output, "candidate", "generation.csv"))
  stopifnot(identical(before, after))
  result <- lapply(seq_len(nrow(before)), function(i) {
    name <- before$case[[i]]
    files <- file.path(output, c("baseline", "candidate"), name)
    if (before$expected_error[[i]]) {
      code <- lapply(paste0(files, ".error"), readLines, warn = FALSE)
      sdata <- lapply(paste0(files, ".data-error"), readLines, warn = FALSE)
    } else {
      code <- lapply(paste0(files, ".stan"), readLines, warn = FALSE)
      sdata <- lapply(paste0(files, ".rds"), readRDS)
    }
    code_equal <- identical(code[[1L]], code[[2L]])
    if (mode == "s2z" && !before$expected_error[[i]]) {
      code <- lapply(code, remove_unused_legacy_helper)
    }
    code_compatible <- identical(code[[1L]], code[[2L]])
    data_equal <- identical(sdata[[1L]], sdata[[2L]])
    if (!data_equal) {
      writeLines(capture.output(all.equal(sdata[[1L]], sdata[[2L]])),
        file.path(output, paste0(name, ".data-diff.txt")))
    }
    data.frame(case = name, code_identical = code_equal,
      unused_helper_cleanup = !code_equal && code_compatible,
      code_compatible = code_compatible, data_identical = data_equal)
  })
  result <- do.call(rbind, result)
  utils::write.csv(result, file.path(output, "comparison.csv"),
    row.names = FALSE)
  options(width = max(120L, getOption("width")))
  print(result, row.names = FALSE)
  cat("Artifacts: ", output, "\n", sep = "")
  if (!all(result$code_compatible & result$data_identical)) {
    stop("Compatibility comparison failed; inspect retained source and data.")
  }
  cat(nrow(result), " ", mode,
    " fixtures agree in generated Stan source/data or expected diagnostics",
    if (any(result$unused_helper_cleanup)) {
      " after the unused-helper deletion"
    } else "",
    ".\n", sep = "")
}

args <- commandArgs(trailingOnly = TRUE)
worker <- length(args) && identical(args[[1L]], "--generate")
if (worker) args <- args[-1L]
if (length(args) < 3L || length(args) > 4L ||
    !args[[1L]] %in% c("conventional", "s2z")) {
  stop(paste(
    "Usage: Rscript tests/local/tests.s2z-compatibility.R",
    "conventional|s2z BASELINE_LIB CANDIDATE_LIB [OUTPUT_DIR]"
  ), call. = FALSE)
}
if (worker) {
  if (length(args) != 3L) stop("Invalid generation arguments")
  generate_compatibility(args[[1L]], args[[2L]], args[[3L]])
} else {
  invocation <- commandArgs(trailingOnly = FALSE)
  script <- sub("^--file=", "", invocation[grepl("^--file=", invocation)])
  stopifnot(length(script) == 1L)
  output <- if (length(args) == 4L) {
    args[[4L]]
  } else tempfile("s2z-compatibility-")
  compare_compatibility(args[[1L]], args[[2L]], args[[3L]], output,
    normalizePath(script, mustWork = TRUE))
}
