# Real-model validation of sum-to-zero group effects against PosteriorDB
# models.
#
# This is intentionally a local test. It downloads data and model provenance
# from a pinned PosteriorDB commit, compiles conventional and selected S2Z
# brms models, and then samples each compiled model without recompiling it.
# Four cases are exact brms representations of their PosteriorDB likelihoods
# and priors:
#
#   - eight schools: scalar Gaussian group effects with known standard errors;
#   - GLMM1: scalar Poisson varying site intercepts;
#   - radon: independent Gaussian varying intercepts and slopes;
#   - rats: independent Gaussian varying intercepts and slopes.
#
# A fifth case uses the PosteriorDB rats data with the conventional, correlated
# `weight ~ age + (age | rat)` model and brms default priors. It is a brms-native
# extension, not an exact representation of a PosteriorDB posterior. Its model
# is from Aki Vehtari's rats case study:
# https://github.com/avehtari/casestudies/blob/
# fc9e544f7133af78557f55c329507916e73d6b45/rats/rats.R#L41-L99
#
# Eight schools is also compared with its official 10,000-draw PosteriorDB
# reference posterior. PosteriorDB currently has no reference draws for the
# other three cases.
#
# Run against an installed development build of brms, for example:
#
#   R_LIBS=/path/to/development-library \
#   Rscript tests/local/tests.models-s2z-posteriordb.R
#
# Useful controls (shown with their defaults) are:
#
#   BRMS_S2Z_PDB_CHAINS=2
#   BRMS_S2Z_PDB_WARMUP=500
#   BRMS_S2Z_PDB_SAMPLING=500
#   BRMS_S2Z_PDB_CORES=min(chains, physical cores)
#   BRMS_S2Z_PDB_BACKEND=auto
#   BRMS_S2Z_PDB_ADAPT_DELTA=0.99
#   BRMS_S2Z_PDB_MAX_TREEDEPTH=14
#   BRMS_S2Z_PDB_REPS=1
#   BRMS_S2Z_PDB_CASES=eight-schools,glmm1,radon,rats,rats-brms-correlated
#   BRMS_S2Z_PDB_S2Z_MODES=centered
#   BRMS_S2Z_PDB_EQUIV_Z=7
#   BRMS_S2Z_PDB_EQUIV_ABS=0.06
#   BRMS_S2Z_PDB_EQUIV_REL=0.125
#   BRMS_S2Z_PDB_MAX_RHAT=1.10
#   BRMS_S2Z_PDB_MIN_BULK_ESS=40
#   BRMS_S2Z_PDB_MIN_TAIL_ESS=50
#   BRMS_S2Z_PDB_MIN_EBFMI=0.30
#   BRMS_S2Z_PDB_MAX_DIVERGENCES=0
#   BRMS_S2Z_PDB_MAX_TREEDEPTH_HITS=0
#   BRMS_S2Z_PDB_SIG_FIGS=17
#   BRMS_S2Z_PDB_INVARIANT_TOL=1e-8
#   BRMS_S2Z_PDB_INVARIANT_REL_TOL=1e-10
#   BRMS_S2Z_PDB_SEED=9100
#   BRMS_S2Z_PDB_OUTPUT_DIR=/private/tmp/brms-s2z-posteriordb-output
#   BRMS_S2Z_PDB_CACHE_DIR=/private/tmp/brms-s2z-posteriordb-cache
#   BRMS_S2Z_PDB_DOWNLOAD_DIR=/private/tmp/brms-s2z-posteriordb-downloads
#   BRMS_S2Z_PDB_ROOT=/path/to/posteriordb
#
# The default keeps the original centered-S2Z validation and runtime. To run
# the brms-only four-way Radon comparison (conventional noncentered, centered
# S2Z, noncentered S2Z, and automatic expected-Fisher S2Z), use:
#
#   BRMS_S2Z_PDB_CASES=radon \
#   BRMS_S2Z_PDB_S2Z_MODES=centered,noncentered,auto \
#   Rscript tests/local/tests.models-s2z-posteriordb.R
#
# Equivalence is checked for posterior means, standard deviations, and the
# 10th, 50th, and 90th percentiles. Every comparison must pass both a
# propagated-MCSE gate and a prespecified absolute-plus-relative margin.

suppressPackageStartupMessages({
  library(brms)
  library(posterior)
})

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required.", call. = FALSE)
}

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

env_number <- function(name, default, minimum = -Inf, maximum = Inf) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(as.numeric(default))
  }
  out <- suppressWarnings(as.numeric(value))
  if (length(out) != 1L || !is.finite(out) || out < minimum ||
      out > maximum) {
    stop(
      name, " must be a finite number in [", minimum, ", ", maximum, "].",
      call. = FALSE
    )
  }
  out
}

env_flag <- function(name, default = FALSE) {
  value <- tolower(Sys.getenv(name, unset = if (default) "true" else "false"))
  if (!value %in% c("true", "false", "1", "0", "yes", "no")) {
    stop(name, " must be true or false.", call. = FALSE)
  }
  value %in% c("true", "1", "yes")
}

choose_backend <- function() {
  requested <- tolower(Sys.getenv("BRMS_S2Z_PDB_BACKEND", unset = "auto"))
  if (requested %in% c("cmdstanr", "rstan")) {
    return(requested)
  }
  if (!identical(requested, "auto")) {
    stop(
      "BRMS_S2Z_PDB_BACKEND must be auto, cmdstanr, or rstan.",
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
  stop("Neither a working CmdStan installation nor rstan was found.",
       call. = FALSE)
}

chains <- env_integer("BRMS_S2Z_PDB_CHAINS", 2L, minimum = 2L)
warmup <- env_integer("BRMS_S2Z_PDB_WARMUP", 500L)
sampling <- env_integer("BRMS_S2Z_PDB_SAMPLING", 500L)
detected_cores <- parallel::detectCores(logical = FALSE)
if (is.na(detected_cores)) {
  detected_cores <- 1L
}
cores <- env_integer(
  "BRMS_S2Z_PDB_CORES", min(chains, detected_cores), minimum = 1L
)
backend <- choose_backend()
adapt_delta <- env_number(
  "BRMS_S2Z_PDB_ADAPT_DELTA", 0.99, minimum = 0.5, maximum = 1
)
max_treedepth <- env_integer("BRMS_S2Z_PDB_MAX_TREEDEPTH", 14L)
reps <- env_integer("BRMS_S2Z_PDB_REPS", 1L)
equiv_z <- env_number("BRMS_S2Z_PDB_EQUIV_Z", 7, minimum = 1)
equiv_abs <- env_number("BRMS_S2Z_PDB_EQUIV_ABS", 0.06, minimum = 0)
equiv_rel <- env_number("BRMS_S2Z_PDB_EQUIV_REL", 0.125, minimum = 0)
max_rhat <- env_number("BRMS_S2Z_PDB_MAX_RHAT", 1.10, minimum = 1)
min_bulk_ess <- env_number(
  "BRMS_S2Z_PDB_MIN_BULK_ESS", 40, minimum = 0
)
min_tail_ess <- env_number(
  "BRMS_S2Z_PDB_MIN_TAIL_ESS", 50, minimum = 0
)
min_ebfmi <- env_number(
  "BRMS_S2Z_PDB_MIN_EBFMI", 0.30, minimum = 0
)
max_divergences <- env_integer(
  "BRMS_S2Z_PDB_MAX_DIVERGENCES", 0L, minimum = 0L
)
max_treedepth_hits <- env_integer(
  "BRMS_S2Z_PDB_MAX_TREEDEPTH_HITS", 0L, minimum = 0L
)
invariant_tolerance <- env_number(
  "BRMS_S2Z_PDB_INVARIANT_TOL", 1e-8, minimum = 0
)
invariant_relative_tolerance <- env_number(
  "BRMS_S2Z_PDB_INVARIANT_REL_TOL", 1e-10, minimum = 0
)
output_sig_figs <- env_integer("BRMS_S2Z_PDB_SIG_FIGS", 17L, minimum = 1L)
if (output_sig_figs > 18L) {
  stop("BRMS_S2Z_PDB_SIG_FIGS must be at most 18.", call. = FALSE)
}
base_seed <- env_integer("BRMS_S2Z_PDB_SEED", 9100L, minimum = 1L)
output_dir <- Sys.getenv(
  "BRMS_S2Z_PDB_OUTPUT_DIR",
  unset = "/private/tmp/brms-s2z-posteriordb-output"
)
cache_dir <- Sys.getenv(
  "BRMS_S2Z_PDB_CACHE_DIR",
  unset = "/private/tmp/brms-s2z-posteriordb-cache"
)
download_dir <- Sys.getenv(
  "BRMS_S2Z_PDB_DOWNLOAD_DIR",
  unset = "/private/tmp/brms-s2z-posteriordb-downloads"
)
force_download <- env_flag("BRMS_S2Z_PDB_FORCE_DOWNLOAD", FALSE)
pdb_root <- Sys.getenv(
  "BRMS_S2Z_PDB_ROOT",
  unset = Sys.getenv("POSTERIOR_DB_PATH", unset = "")
)

# Pin every downloaded data, model, metadata, and reference-draw artifact to
# one immutable PosteriorDB tree. Override only to audit a deliberate update.
pdb_commit <- Sys.getenv(
  "BRMS_S2Z_PDB_COMMIT",
  unset = "28f8d3d6e975315f42aa274a8399f21e07a43b30"
)
if (!grepl("^[[:xdigit:]]{40}$", pdb_commit)) {
  stop("BRMS_S2Z_PDB_COMMIT must be a 40-character Git SHA.", call. = FALSE)
}
# Keep cached downloads from different immutable trees physically separate.
download_dir <- file.path(download_dir, pdb_commit)
pdb_base_url <- paste0(
  "https://raw.githubusercontent.com/stan-dev/posteriordb/", pdb_commit,
  "/posterior_database/"
)
if (nzchar(pdb_root)) {
  pdb_root <- normalizePath(pdb_root, mustWork = TRUE)
  if (basename(pdb_root) == "posterior_database") {
    pdb_database_root <- pdb_root
    pdb_root <- dirname(pdb_root)
  } else {
    pdb_database_root <- file.path(pdb_root, "posterior_database")
  }
  if (!dir.exists(pdb_database_root)) {
    stop(
      "BRMS_S2Z_PDB_ROOT must be a PosteriorDB repository root or its ",
      "posterior_database directory: ", pdb_root, call. = FALSE
    )
  }
  pdb_database_root <- normalizePath(pdb_database_root, mustWork = TRUE)
  local_commit <- suppressWarnings(system2(
    "git", c("-C", shQuote(pdb_root), "rev-parse", "HEAD"),
    stdout = TRUE, stderr = TRUE
  ))
  local_status <- attr(local_commit, "status")
  if (is.null(local_status)) {
    local_status <- 0L
  }
  if (local_status != 0L || length(local_commit) != 1L ||
      !identical(tolower(local_commit), tolower(pdb_commit))) {
    stop(
      "BRMS_S2Z_PDB_ROOT must be a Git checkout at the pinned commit ",
      pdb_commit, "; found ", paste(local_commit, collapse = " "), ".",
      call. = FALSE
    )
  }
  local_changes <- suppressWarnings(system2(
    "git",
    c(
      "-C", shQuote(pdb_root), "status", "--porcelain",
      "--untracked-files=no", "--", "posterior_database"
    ),
    stdout = TRUE, stderr = TRUE
  ))
  changes_status <- attr(local_changes, "status")
  if (is.null(changes_status)) {
    changes_status <- 0L
  }
  if (changes_status != 0L || length(local_changes)) {
    stop(
      "BRMS_S2Z_PDB_ROOT has tracked changes under posterior_database; ",
      "use a clean checkout of the pinned commit.", call. = FALSE
    )
  }
} else {
  pdb_database_root <- ""
}

all_case_names <- c(
  "eight-schools", "glmm1", "radon", "rats", "rats-brms-correlated"
)
requested_cases <- strsplit(
  Sys.getenv(
    "BRMS_S2Z_PDB_CASES", unset = paste(all_case_names, collapse = ",")
  ),
  ",", fixed = TRUE
)[[1]]
requested_cases <- trimws(requested_cases)
requested_cases <- requested_cases[nzchar(requested_cases)]
unknown_cases <- setdiff(requested_cases, all_case_names)
if (length(unknown_cases)) {
  stop(
    "Unknown BRMS_S2Z_PDB_CASES value(s): ",
    paste(unknown_cases, collapse = ", "), call. = FALSE
  )
}
if (!length(requested_cases)) {
  stop("BRMS_S2Z_PDB_CASES selected no cases.", call. = FALSE)
}

all_s2z_modes <- c("centered", "noncentered", "auto")
all_parameterizations <- c(
  "conventional", paste0("s2z-", all_s2z_modes)
)
requested_s2z_modes <- strsplit(
  Sys.getenv("BRMS_S2Z_PDB_S2Z_MODES", unset = "centered"),
  ",", fixed = TRUE
)[[1]]
requested_s2z_modes <- unique(trimws(requested_s2z_modes))
requested_s2z_modes <- requested_s2z_modes[nzchar(requested_s2z_modes)]
unknown_s2z_modes <- setdiff(requested_s2z_modes, all_s2z_modes)
if (length(unknown_s2z_modes)) {
  stop(
    "Unknown BRMS_S2Z_PDB_S2Z_MODES value(s): ",
    paste(unknown_s2z_modes, collapse = ", "), call. = FALSE
  )
}
if (!length(requested_s2z_modes)) {
  stop("BRMS_S2Z_PDB_S2Z_MODES selected no modes.", call. = FALSE)
}
parameterizations <- c(
  "conventional", paste0("s2z-", requested_s2z_modes)
)
s2z_parameterizations <- setdiff(parameterizations, "conventional")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)
options(mc.cores = cores, timeout = max(300, getOption("timeout")), width = 150)
if (backend == "cmdstanr") {
  # Keep the Stan sources and executables beside the serialized templates so
  # cached CmdStan models remain usable after the current R session ends.
  options(cmdstanr_write_stan_file_dir = cache_dir)
}

git_output <- function(args) {
  tryCatch(
    paste(system2("git", args, stdout = TRUE, stderr = TRUE), collapse = "\n"),
    error = function(e) NA_character_
  )
}

harness_path <- normalizePath(
  file.path("tests", "local", "tests.models-s2z-posteriordb.R"),
  mustWork = TRUE
)
configuration <- data.frame(
  brms_version = as.character(utils::packageVersion("brms")),
  posterior_version = as.character(utils::packageVersion("posterior")),
  jsonlite_version = as.character(utils::packageVersion("jsonlite")),
  backend = backend,
  chains = chains,
  warmup = warmup,
  sampling = sampling,
  cores = cores,
  adapt_delta = adapt_delta,
  max_treedepth = max_treedepth,
  reps = reps,
  equiv_z = equiv_z,
  equiv_abs = equiv_abs,
  equiv_rel = equiv_rel,
  max_rhat = max_rhat,
  min_bulk_ess = min_bulk_ess,
  min_tail_ess = min_tail_ess,
  min_ebfmi = min_ebfmi,
  max_divergences = max_divergences,
  max_treedepth_hits = max_treedepth_hits,
  invariant_tolerance = invariant_tolerance,
  invariant_relative_tolerance = invariant_relative_tolerance,
  output_sig_figs = if (backend == "cmdstanr") output_sig_figs else NA_integer_,
  base_seed = base_seed,
  requested_cases = paste(requested_cases, collapse = ","),
  requested_s2z_modes = paste(requested_s2z_modes, collapse = ","),
  parameterizations = paste(parameterizations, collapse = ","),
  pdb_commit = pdb_commit,
  pdb_root = if (nzchar(pdb_root)) pdb_root else NA_character_,
  git_head = git_output(c("rev-parse", "HEAD")),
  git_status = git_output(c("status", "--short")),
  harness_md5 = unname(tools::md5sum(harness_path)),
  output_dir = normalizePath(output_dir, mustWork = TRUE),
  cache_dir = normalizePath(cache_dir, mustWork = TRUE),
  download_dir = normalizePath(download_dir, mustWork = TRUE),
  stringsAsFactors = FALSE
)
utils::write.csv(
  configuration, file.path(output_dir, "s2z-posteriordb-settings.csv"),
  row.names = FALSE
)
cat("PosteriorDB S2Z validation configuration\n")
print(configuration, row.names = FALSE)

artifact_paths <- list(
  `eight-schools` = c(
    data = "data/data/eight_schools.json.zip",
    model = "models/stan/eight_schools_noncentered.stan",
    posterior = paste0(
      "posteriors/eight_schools-eight_schools_noncentered.json"
    ),
    reference = paste0(
      "reference_posteriors/draws/draws/",
      "eight_schools-eight_schools_noncentered.json.zip"
    )
  ),
  glmm1 = c(
    data = "data/data/GLMM_data.json.zip",
    model = "models/stan/GLMM1_model.stan",
    posterior = "posteriors/GLMM_data-GLMM1_model.json"
  ),
  radon = c(
    data = "data/data/radon_mn.json.zip",
    model = paste0(
      "models/stan/radon_variable_intercept_slope_noncentered.stan"
    ),
    posterior = paste0(
      "posteriors/",
      "radon_mn-radon_variable_intercept_slope_noncentered.json"
    )
  ),
  rats = c(
    data = "data/data/rats_data.json.zip",
    model = "models/stan/rats_model.stan",
    posterior = "posteriors/rats_data-rats_model.json"
  ),
  `rats-brms-correlated` = c(
    data = "data/data/rats_data.json.zip"
  )
)

download_artifact <- function(relative_path) {
  if (nzchar(pdb_database_root)) {
    local_path <- file.path(pdb_database_root, relative_path)
    if (!file.exists(local_path)) {
      stop(
        "PosteriorDB artifact is missing under BRMS_S2Z_PDB_ROOT: ",
        local_path, call. = FALSE
      )
    }
    return(normalizePath(local_path, mustWork = TRUE))
  }
  destination <- file.path(download_dir, relative_path)
  if (file.exists(destination) && !force_download) {
    return(normalizePath(destination))
  }
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  url <- paste0(pdb_base_url, relative_path)
  temporary <- tempfile(
    pattern = paste0(basename(destination), "-"), tmpdir = dirname(destination)
  )
  on.exit(unlink(temporary), add = TRUE)
  cat("Downloading ", url, " ...\n", sep = "")
  status <- tryCatch(
    utils::download.file(url, temporary, mode = "wb", quiet = TRUE),
    error = function(e) e
  )
  if (inherits(status, "error") || !identical(status, 0L) ||
      !file.exists(temporary) || file.info(temporary)$size <= 0) {
    detail <- if (inherits(status, "error")) conditionMessage(status) else status
    stop(
      "Could not download pinned PosteriorDB artifact ", url, ": ", detail,
      call. = FALSE
    )
  }
  if (!file.rename(temporary, destination)) {
    stop("Could not move downloaded artifact to ", destination, call. = FALSE)
  }
  normalizePath(destination)
}

artifact_rows <- list()
artifact_files <- list()
row <- 0L
for (case_name in requested_cases) {
  artifact_files[[case_name]] <- list()
  for (kind in names(artifact_paths[[case_name]])) {
    relative_path <- artifact_paths[[case_name]][[kind]]
    path <- download_artifact(relative_path)
    artifact_files[[case_name]][[kind]] <- path
    row <- row + 1L
    artifact_rows[[row]] <- data.frame(
      case = case_name,
      kind = kind,
      relative_path = relative_path,
      source = if (nzchar(pdb_database_root)) {
        paste0("local:", path)
      } else {
        paste0(pdb_base_url, relative_path)
      },
      local_path = path,
      bytes = file.info(path)$size,
      md5 = unname(tools::md5sum(path)),
      stringsAsFactors = FALSE
    )
  }
}
artifact_manifest <- do.call(rbind, artifact_rows)
utils::write.csv(
  artifact_manifest,
  file.path(output_dir, "s2z-posteriordb-artifacts.csv"), row.names = FALSE
)

read_zipped_json <- function(path, simplify = TRUE) {
  members <- utils::unzip(path, list = TRUE)$Name
  json_members <- members[grepl("[.]json$", members)]
  if (length(json_members) != 1L) {
    stop(
      "Expected exactly one JSON member in ", path, "; found ",
      length(json_members), ".", call. = FALSE
    )
  }
  connection <- unz(path, json_members, open = "rb")
  on.exit(close(connection), add = TRUE)
  jsonlite::parse_json(
    paste(readLines(connection, warn = FALSE), collapse = "\n"),
    simplifyVector = simplify
  )
}

assert_index <- function(index, n_level, label) {
  index <- as.integer(index)
  if (!length(index) || anyNA(index) || any(index < 1L) ||
      any(index > n_level)) {
    stop("Invalid grouping index in ", label, ".", call. = FALSE)
  }
  missing <- setdiff(seq_len(n_level), unique(index))
  if (length(missing)) {
    stop(
      label, " contains grouping levels with no likelihood observation: ",
      paste(missing, collapse = ", "), ".", call. = FALSE
    )
  }
  index
}

make_eight_schools_data <- function(path) {
  x <- read_zipped_json(path)
  stopifnot(length(x$y) == x$J, length(x$sigma) == x$J)
  data.frame(
    y = as.numeric(x$y),
    sei = as.numeric(x$sigma),
    one = 1,
    school = factor(seq_len(x$J), levels = seq_len(x$J))
  )
}

make_glmm1_data <- function(path) {
  x <- read_zipped_json(path)
  site <- assert_index(x$obssite, x$nsite, "GLMM_data")
  stopifnot(length(x$obs) == x$nobs, length(site) == x$nobs)
  data.frame(
    obs = as.integer(x$obs),
    one = 1,
    site = factor(site, levels = seq_len(x$nsite))
  )
}

make_radon_data <- function(path) {
  x <- read_zipped_json(path)
  county <- assert_index(x$county_idx, x$J, "radon_mn")
  stopifnot(
    length(x$log_radon) == x$N,
    length(x$floor_measure) == x$N,
    length(county) == x$N
  )
  data.frame(
    y = as.numeric(x$log_radon),
    one = 1,
    floor = as.numeric(x$floor_measure),
    county = factor(county, levels = seq_len(x$J))
  )
}

make_rats_data <- function(path) {
  x <- read_zipped_json(path)
  rat <- assert_index(x$rat, x$N, "rats_data")
  stopifnot(
    length(x$y) == x$Npts, length(x$x) == x$Npts,
    length(rat) == x$Npts, length(x$xbar) == 1L
  )
  data.frame(
    y = as.numeric(x$y),
    weight = as.numeric(x$y),
    age = as.numeric(x$x),
    one = 1,
    xc = as.numeric(x$x) - as.numeric(x$xbar),
    rat = factor(rat, levels = seq_len(x$N))
  )
}

make_case <- function(name) {
  files <- artifact_files[[name]]
  if (name == "eight-schools") {
    data <- make_eight_schools_data(files$data)
    return(list(
      name = name,
      posterior_name = "eight_schools-eight_schools_noncentered",
      data = data,
      family = gaussian(),
      conventional_formula = bf(
        y | se(sei, sigma = FALSE) ~ 0 + one + (0 + one | school)
      ),
      s2z_centered_formula = bf(
        y | se(sei, sigma = FALSE) ~ 0 + one +
          (0 + one | gr(school, s2z = TRUE))
      ),
      s2z_noncentered_formula = bf(
        y | se(sei, sigma = FALSE) ~ 0 + one +
          (0 + one | gr(school, s2z = TRUE, center = FALSE))
      ),
      s2z_auto_formula = NULL,
      auto_supported = FALSE,
      auto_skip_reason = paste0(
        "the known-standard-error response addition is not yet included ",
        "in the automatic expected-Fisher target"
      ),
      prior = prior(normal(0, 5), class = "b", coef = "one") +
        prior(cauchy(0, 5), class = "sd", group = "school"),
      fixed_formula = ~ 0 + one,
      group = "school",
      group_coef = "one",
      include_sigma = FALSE,
      reference_path = files$reference
    ))
  }
  if (name == "glmm1") {
    data <- make_glmm1_data(files$data)
    return(list(
      name = name,
      posterior_name = "GLMM_data-GLMM1_model",
      data = data,
      family = poisson(),
      conventional_formula = bf(
        obs ~ 0 + one + (0 + one | site)
      ),
      s2z_centered_formula = bf(
        obs ~ 0 + one + (0 + one | gr(site, s2z = TRUE))
      ),
      s2z_noncentered_formula = bf(
        obs ~ 0 + one +
          (0 + one | gr(site, s2z = TRUE, center = FALSE))
      ),
      s2z_auto_formula = NULL,
      auto_supported = FALSE,
      auto_skip_reason = paste0(
        "automatic expected-Fisher centering currently requires a ",
        "Gaussian or Student-t identity likelihood"
      ),
      prior = prior(normal(0, 10), class = "b", coef = "one") +
        prior(uniform(0, 5), class = "sd", group = "site", ub = 5),
      fixed_formula = ~ 0 + one,
      group = "site",
      group_coef = "one",
      include_sigma = FALSE,
      reference_path = NULL
    ))
  }
  if (name == "radon") {
    data <- make_radon_data(files$data)
    return(list(
      name = name,
      posterior_name = paste0(
        "radon_mn-radon_variable_intercept_slope_noncentered"
      ),
      data = data,
      family = gaussian(),
      conventional_formula = bf(
        y ~ 0 + one + floor + (0 + one + floor || county)
      ),
      s2z_centered_formula = bf(
        y ~ 0 + one + floor +
          (0 + one + floor || gr(county, s2z = TRUE))
      ),
      s2z_noncentered_formula = bf(
        y ~ 0 + one + floor +
          (0 + one + floor ||
            gr(county, s2z = TRUE, center = FALSE))
      ),
      s2z_auto_formula = bf(
        y ~ 0 + one + floor +
          (0 + one + floor ||
            gr(county, s2z = TRUE, center = "auto"))
      ),
      auto_supported = TRUE,
      auto_skip_reason = NA_character_,
      prior = prior(normal(0, 10), class = "b") +
        prior(normal(0, 1), class = "sd", group = "county") +
        prior(normal(0, 1), class = "sigma"),
      fixed_formula = ~ 0 + one + floor,
      group = "county",
      group_coef = c("one", "floor"),
      include_sigma = TRUE,
      reference_path = NULL
    ))
  }
  if (name == "rats") {
    data <- make_rats_data(files$data)
    return(list(
      name = name,
      posterior_name = "rats_data-rats_model",
      data = data,
      family = gaussian(),
      conventional_formula = bf(
        y ~ 0 + one + xc + (0 + one + xc || rat)
      ),
      s2z_centered_formula = bf(
        y ~ 0 + one + xc +
          (0 + one + xc || gr(rat, s2z = TRUE))
      ),
      s2z_noncentered_formula = bf(
        y ~ 0 + one + xc +
          (0 + one + xc || gr(rat, s2z = TRUE, center = FALSE))
      ),
      s2z_auto_formula = bf(
        y ~ 0 + one + xc +
          (0 + one + xc || gr(rat, s2z = TRUE, center = "auto"))
      ),
      auto_supported = TRUE,
      auto_skip_reason = NA_character_,
      prior = prior(normal(0, 100), class = "b") +
        # A global empty sd prior suppresses the brms default and reproduces
        # PosteriorDB's improper flat priors on the two positive group scales.
        prior("", class = "sd") +
        prior("", class = "sigma"),
      fixed_formula = ~ 0 + one + xc,
      group = "rat",
      group_coef = c("one", "xc"),
      centered_intercept = FALSE,
      include_sigma = TRUE,
      reference_path = NULL
    ))
  }
  if (name == "rats-brms-correlated") {
    data <- make_rats_data(files$data)
    return(list(
      name = name,
      posterior_name = NA_character_,
      case_scope = paste0(
        "brms-native correlated extension on PosteriorDB rats data; ",
        "not an exact PosteriorDB posterior"
      ),
      source_url = paste0(
        "https://github.com/avehtari/casestudies/blob/",
        "fc9e544f7133af78557f55c329507916e73d6b45/",
        "rats/rats.R#L41-L99"
      ),
      data = data,
      family = gaussian(),
      conventional_formula = bf(weight ~ age + (age | rat)),
      s2z_centered_formula = bf(
        weight ~ age + (age | gr(rat, s2z = TRUE))
      ),
      s2z_noncentered_formula = bf(
        weight ~ age +
          (age | gr(rat, s2z = TRUE, center = FALSE))
      ),
      s2z_auto_formula = bf(
        weight ~ age +
          (age | gr(rat, s2z = TRUE, center = "auto"))
      ),
      auto_supported = TRUE,
      auto_skip_reason = NA_character_,
      prior = NULL,
      fixed_formula = ~ age,
      group = "rat",
      group_coef = c("Intercept", "age"),
      centered_intercept = TRUE,
      include_sigma = TRUE,
      reference_path = NULL
    ))
  }
  stop("Unknown case: ", name, call. = FALSE)
}

cases <- setNames(lapply(requested_cases, make_case), requested_cases)

case_formula <- function(case, parameterization) {
  if (identical(parameterization, "s2z-auto") &&
      !isTRUE(case$auto_supported)) {
    stop(
      "Automatic expected-Fisher S2Z is unavailable for ", case$name,
      ": ", case$auto_skip_reason, ".", call. = FALSE
    )
  }
  switch(
    parameterization,
    conventional = case$conventional_formula,
    `s2z-centered` = case$s2z_centered_formula,
    `s2z-noncentered` = case$s2z_noncentered_formula,
    `s2z-auto` = case$s2z_auto_formula,
    stop("Unknown parameterization: ", parameterization, call. = FALSE)
  )
}

selected_parameterizations_for_case <- function(case) {
  if (isTRUE(case$auto_supported)) {
    return(parameterizations)
  }
  setdiff(parameterizations, "s2z-auto")
}

parameterizations_by_case <- lapply(
  cases, selected_parameterizations_for_case
)
selected_s2z_by_case <- lapply(
  parameterizations_by_case, setdiff, y = "conventional"
)
if (!length(unique(unlist(selected_s2z_by_case, use.names = FALSE)))) {
  stop(
    "No requested S2Z parameterization is supported by the selected cases. ",
    "In particular, center = \"auto\" is skipped for eight-schools and ",
    "glmm1; select radon or a rats case to test automatic expected-Fisher ",
    "centering.", call. = FALSE
  )
}

auto_eligibility <- do.call(rbind, lapply(cases, function(case) {
  requested <- "s2z-auto" %in% parameterizations
  supported <- isTRUE(case$auto_supported)
  data.frame(
    case = case$name,
    auto_requested = requested,
    auto_supported = supported,
    action = if (!requested) {
      "not requested"
    } else if (supported) {
      "run"
    } else {
      "skip"
    },
    reason = if (supported) NA_character_ else case$auto_skip_reason,
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(
  auto_eligibility,
  file.path(output_dir, "s2z-posteriordb-auto-eligibility.csv"),
  row.names = FALSE
)
if (any(auto_eligibility$action == "skip")) {
  cat("\nSkipping unsupported automatic expected-Fisher jobs\n")
  print(
    auto_eligibility[auto_eligibility$action == "skip", ],
    row.names = FALSE
  )
}

case_manifest <- do.call(rbind, lapply(cases, function(case) {
  data.frame(
    case = case$name,
    case_scope = if (!is.null(case$case_scope)) {
      case$case_scope
    } else {
      "exact PosteriorDB likelihood and prior representation"
    },
    posterior_name = if (!is.null(case$posterior_name)) {
      case$posterior_name
    } else {
      NA_character_
    },
    auto_supported = isTRUE(case$auto_supported),
    auto_skip_reason = case$auto_skip_reason,
    selected_parameterizations = paste(
      parameterizations_by_case[[case$name]], collapse = ","
    ),
    pdb_data_url = paste0(
      "https://github.com/stan-dev/posteriordb/blob/", pdb_commit,
      "/posterior_database/", artifact_paths[[case$name]][["data"]]
    ),
    source_url = if (!is.null(case$source_url)) {
      case$source_url
    } else {
      paste0(
        "https://github.com/stan-dev/posteriordb/tree/", pdb_commit,
        "/posterior_database"
      )
    },
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(
  case_manifest, file.path(output_dir, "s2z-posteriordb-cases.csv"),
  row.names = FALSE
)

validate_case_design <- function(case) {
  X <- stats::model.matrix(case$fixed_formula, case$data)
  fixed_names <- sub("^\\(Intercept\\)$", "Intercept", colnames(X))
  if (!identical(fixed_names, case$group_coef)) {
    stop(
      "Fixed/group coefficient mismatch for ", case$name, ": fixed = ",
      paste(fixed_names, collapse = ", "), "; group = ",
      paste(case$group_coef, collapse = ", "), ".", call. = FALSE
    )
  }
  # Each selected S2Z code-generation call runs the production validation,
  # including exact equality of the fixed and group design columns.
  case_parameterizations <- parameterizations_by_case[[case$name]]
  codes <- setNames(lapply(case_parameterizations, function(parameterization) {
    formula <- case_formula(case, parameterization)
    invisible(standata(
      formula, data = case$data, family = case$family, prior = case$prior
    ))
    stancode(
      formula, data = case$data, family = case$family, prior = case$prior
    )
  }), case_parameterizations)
  if (case$name == "rats") {
    scale_prior_pattern <- "_lpdf(sd_1 |"
    sigma_prior_pattern <- "_lpdf(sigma |"
    if (any(vapply(
      codes,
      function(code) {
        grepl(scale_prior_pattern, code, fixed = TRUE) ||
          grepl(sigma_prior_pattern, code, fixed = TRUE)
      },
      logical(1)
    ))) {
      stop(
        "The exact rats target unexpectedly contains an sd or sigma prior.",
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}
invisible(lapply(cases, validate_case_design))

if ("s2z-auto" %in% parameterizations) {
  for (case_name in names(cases)[vapply(
    cases, function(case) isTRUE(case$auto_supported), logical(1)
  )]) {
    case <- cases[[case_name]]
    sdata <- standata(
      case$s2z_auto_formula, data = case$data, family = case$family,
      prior = case$prior
    )
    if (any(grepl("^rho_s2z_", names(sdata)))) {
      stop(
        "Automatic expected-Fisher weights for ", case_name,
        " were unexpectedly supplied as data.", call. = FALSE
      )
    }
    code <- stancode(
      case$s2z_auto_formula, data = case$data, family = case$family,
      prior = case$prior
    )
    if (!grepl("marginal Fisher centering fractions", code, fixed = TRUE) ||
        !grepl("rho_s2z_", code, fixed = TRUE)) {
      stop(
        "Generated code for ", case_name,
        " did not contain Stan-side expected-Fisher weights.",
        call. = FALSE
      )
    }
  }
  cat(
    "\nAutomatic S2Z weights are computed in Stan from expected Fisher ",
    "information; no rho_s2z_* values are supplied as data.\n",
    sep = ""
  )
}

template_path <- function(case_name, parameterization) {
  file.path(
    cache_dir,
    paste0("template-", case_name, "-", parameterization, "-", backend)
  )
}

compile_template <- function(case, parameterization) {
  formula <- case_formula(case, parameterization)
  label <- paste(case$name, parameterization, sep = "-")
  file <- template_path(case$name, parameterization)
  existed <- file.exists(paste0(file, ".rds")) || file.exists(file)
  cat("Compiling and smoke-running ", label, " ...\n", sep = "")
  timing <- system.time({
    fit <- brm(
      formula = formula,
      data = case$data,
      family = case$family,
      prior = case$prior,
      backend = backend,
      chains = 1,
      cores = 1,
      iter = 2,
      warmup = 1,
      seed = base_seed + 10L * match(case$name, all_case_names) +
        match(parameterization, all_parameterizations),
      refresh = 0,
      silent = 1,
      normalize = TRUE,
      save_pars = save_pars(all = TRUE),
      file = file,
      file_refit = "on_change",
      control = list(adapt_delta = 0.8, max_treedepth = 5)
    )
  })
  list(
    fit = fit,
    report = data.frame(
      case = case$name,
      parameterization = parameterization,
      template_preexisted = existed,
      prepare_and_smoke_seconds = unname(timing[["elapsed"]]),
      stringsAsFactors = FALSE
    )
  )
}

templates <- list()
compile_rows <- list()
row <- 0L
for (case_name in requested_cases) {
  templates[[case_name]] <- list()
  for (parameterization in parameterizations_by_case[[case_name]]) {
    result <- compile_template(cases[[case_name]], parameterization)
    templates[[case_name]][[parameterization]] <- result$fit
    row <- row + 1L
    compile_rows[[row]] <- result$report
  }
}
compile_results <- do.call(rbind, compile_rows)
utils::write.csv(
  compile_results,
  file.path(output_dir, "s2z-posteriordb-compile-smoke.csv"),
  row.names = FALSE
)

draw_variable <- function(fit, variable) {
  if (!variable %in% variables(fit)) {
    stop("Required saved variable is missing: ", variable, call. = FALSE)
  }
  out <- as_draws_array(fit, variable = variable)
  matrix(
    out[, , 1L], nrow = dim(out)[1L], ncol = dim(out)[2L],
    dimnames = dimnames(out)[1:2]
  )
}

reference_variable <- function(draws, variable) {
  if (!variable %in% variables(draws)) {
    stop("Required reference variable is missing: ", variable,
         call. = FALSE)
  }
  out <- as_draws_array(subset_draws(draws, variable = variable))
  matrix(
    out[, , 1L], nrow = dim(out)[1L], ncol = dim(out)[2L],
    dimnames = dimnames(out)[1:2]
  )
}

add_quantity <- function(out, name, value) {
  if (!is.matrix(value) || !is.numeric(value) || anyNA(value)) {
    stop("Quantity '", name, "' must be a numeric iteration-by-chain matrix.",
         call. = FALSE)
  }
  out[[name]] <- value
  out
}

public_quantities <- function(fit, case) {
  out <- list()
  fixed <- setNames(
    lapply(case$group_coef, function(coef) draw_variable(fit, paste0("b_", coef))),
    case$group_coef
  )
  for (coef in case$group_coef) {
    out <- add_quantity(out, paste0("b:", coef), fixed[[coef]])
  }

  group_levels <- levels(case$data[[case$group]])
  totals <- vector("list", length(case$group_coef))
  names(totals) <- case$group_coef
  for (coef in case$group_coef) {
    sd_name <- paste0("sd_", case$group, "__", coef)
    out <- add_quantity(
      out, paste0("parameter:", sd_name), draw_variable(fit, sd_name)
    )
    totals[[coef]] <- vector("list", length(group_levels))
    names(totals[[coef]]) <- group_levels
    for (level in group_levels) {
      r_name <- paste0("r_", case$group, "[", level, ",", coef, "]")
      r <- draw_variable(fit, r_name)
      total <- fixed[[coef]] + r
      out <- add_quantity(
        out,
        paste0("r:", case$group, "[", level, ",", coef, "]"), r
      )
      out <- add_quantity(
        out,
        paste0("coef:", case$group, "[", level, ",", coef, "]"), total
      )
      totals[[coef]][[level]] <- total
    }
  }

  cor_names <- grep(
    paste0("^cor_", case$group, "__"), variables(fit), value = TRUE
  )
  for (cor_name in cor_names) {
    out <- add_quantity(
      out, paste0("parameter:", cor_name), draw_variable(fit, cor_name)
    )
  }

  if (isTRUE(case$include_sigma) && "sigma" %in% variables(fit)) {
    out <- add_quantity(out, "parameter:sigma", draw_variable(fit, "sigma"))
  }

  X <- stats::model.matrix(case$fixed_formula, case$data)
  group_index <- as.integer(case$data[[case$group]])
  selected_obs <- unique(round(seq(1, nrow(X), length.out = 7L)))
  eta_average <- matrix(
    0, nrow = nrow(fixed[[1L]]), ncol = ncol(fixed[[1L]])
  )
  for (n in seq_len(nrow(X))) {
    eta <- matrix(0, nrow = nrow(eta_average), ncol = ncol(eta_average))
    level <- group_levels[group_index[n]]
    for (k in seq_along(case$group_coef)) {
      coef <- case$group_coef[k]
      eta <- eta + X[n, k] * totals[[coef]][[level]]
    }
    eta_average <- eta_average + eta / nrow(X)
    if (n %in% selected_obs) {
      out <- add_quantity(out, sprintf("eta:observation[%d]", n), eta)
    }
  }
  add_quantity(out, "eta:observation-average", eta_average)
}

read_reference_draws <- function(path) {
  raw <- read_zipped_json(path, simplify = FALSE)
  raw <- lapply(raw, function(chain) lapply(chain, as.numeric))
  posterior::as_draws_list(raw)
}

eight_schools_reference_quantities <- function(path, case) {
  draws <- read_reference_draws(path)
  mu <- reference_variable(draws, "mu")
  tau <- reference_variable(draws, "tau")
  out <- list()
  out <- add_quantity(out, "b:one", mu)
  out <- add_quantity(out, "parameter:sd_school__one", tau)
  eta_average <- matrix(0, nrow = nrow(mu), ncol = ncol(mu))
  group_levels <- levels(case$data$school)
  selected_obs <- unique(round(seq(
    1, length(group_levels), length.out = 7L
  )))
  for (j in seq_along(group_levels)) {
    theta <- reference_variable(draws, sprintf("theta[%d]", j))
    r <- theta - mu
    out <- add_quantity(out, paste0("r:school[", j, ",one]"), r)
    out <- add_quantity(out, paste0("coef:school[", j, ",one]"), theta)
    if (j %in% selected_obs) {
      out <- add_quantity(out, sprintf("eta:observation[%d]", j), theta)
    }
    eta_average <- eta_average + theta / length(group_levels)
  }
  add_quantity(out, "eta:observation-average", eta_average)
}

one_variable_draws <- function(x, name = "x") {
  out <- array(
    x,
    dim = c(nrow(x), ncol(x), 1L),
    dimnames = list(
      iteration = NULL,
      chain = paste0("chain", seq_len(ncol(x))),
      variable = name
    )
  )
  as_draws_array(out)
}

quantity_summary <- function(x) {
  draws <- one_variable_draws(x)
  ans <- summarise_draws(
    draws,
    mean = base::mean,
    sd = stats::sd,
    mcse_mean = mcse_mean,
    mcse_sd = mcse_sd,
    rhat = rhat,
    ess_bulk = ess_bulk,
    ess_tail = ess_tail
  )
  out <- as.list(ans[1L, c(
    "mean", "sd", "mcse_mean", "mcse_sd", "rhat", "ess_bulk", "ess_tail"
  )])
  quantiles <- stats::quantile(c(x), probs = c(0.1, 0.5, 0.9))
  names(quantiles) <- c("q10", "q50", "q90")
  quantile_mcse <- mcse_quantile(draws, probs = c(0.1, 0.5, 0.9))
  c(out, as.list(quantiles), as.list(quantile_mcse))
}

compare_quantity_sets <- function(lhs, rhs, case_name, comparison, replicate) {
  if (!setequal(names(lhs), names(rhs))) {
    stop(
      "Canonical quantities differ for ", case_name, " (", comparison,
      "). Only lhs: ", paste(setdiff(names(lhs), names(rhs)), collapse = ", "),
      "; only rhs: ", paste(setdiff(names(rhs), names(lhs)), collapse = ", "),
      call. = FALSE
    )
  }
  rows <- list()
  row <- 0L
  for (name in sort(names(lhs))) {
    x <- quantity_summary(lhs[[name]])
    y <- quantity_summary(rhs[[name]])
    for (metric in c("mean", "sd", "q10", "q50", "q90")) {
      row <- row + 1L
      mcse_name <- paste0("mcse_", metric)
      difference <- y[[metric]] - x[[metric]]
      mcse_difference <- sqrt(x[[mcse_name]]^2 + y[[mcse_name]]^2)
      z_score <- if (is.finite(mcse_difference) && mcse_difference > 0) {
        abs(difference) / mcse_difference
      } else if (identical(x[[metric]], y[[metric]])) {
        0
      } else {
        Inf
      }
      # Scale the practical margin by posterior uncertainty as well as the
      # statistic itself. This avoids making near-zero means and quantiles
      # arbitrarily harder to match than the rest of the same distribution.
      equiv_margin <- equiv_abs + equiv_rel * max(
        abs(x[[metric]]), abs(y[[metric]]), x$sd, y$sd
      )
      mcse_pass <- is.finite(z_score) && z_score <= equiv_z
      margin_pass <- is.finite(difference) &&
        abs(difference) <= equiv_margin
      rows[[row]] <- data.frame(
        case = case_name,
        replicate = replicate,
        comparison = comparison,
        quantity = name,
        metric = metric,
        lhs = x[[metric]],
        rhs = y[[metric]],
        difference = difference,
        mcse_difference = mcse_difference,
        z_score = z_score,
        equiv_margin = equiv_margin,
        mcse_pass = mcse_pass,
        margin_pass = margin_pass,
        pass = mcse_pass && margin_pass,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  out[order(out$z_score, decreasing = TRUE), ]
}

internal_r_variable <- function(fit, k, n_group, n_coef) {
  candidates <- if (n_coef == 1L) {
    sprintf("r_s2z_1_1[%d]", seq_len(n_group))
  } else if ("r_s2z_1[1,1]" %in% variables(fit)) {
    sprintf("r_s2z_1[%d,%d]", seq_len(n_group), k)
  } else {
    sprintf("r_s2z_1_%d[%d]", k, seq_len(n_group))
  }
  lapply(candidates, function(variable) draw_variable(fit, variable))
}

s2z_invariants <- function(fit, case, parameterization, replicate) {
  X <- stats::model.matrix(case$fixed_formula, case$data)
  group_levels <- levels(case$data[[case$group]])
  group_index <- as.integer(case$data[[case$group]])
  n_group <- length(group_levels)
  n_coef <- length(case$group_coef)
  theta <- lapply(
    seq_len(n_coef),
    function(k) draw_variable(fit, sprintf("theta_s2z[%d]", k))
  )
  theta_raw <- theta
  if (isTRUE(case$centered_intercept) && n_coef > 1L) {
    means_X <- colMeans(X[, -1L, drop = FALSE])
    for (k in 2:n_coef) {
      theta_raw[[1L]] <- theta_raw[[1L]] - means_X[k - 1L] * theta[[k]]
    }
  }
  internal_r <- lapply(
    seq_len(n_coef),
    function(k) internal_r_variable(fit, k, n_group, n_coef)
  )

  max_sum_to_zero_error <- 0
  max_sum_to_zero_scale <- 1
  max_finite_error <- 0
  max_finite_scale <- 1
  n_iteration <- nrow(theta[[1L]])
  n_chain <- ncol(theta[[1L]])
  eta_internal <- array(0, dim = c(n_iteration, n_chain, nrow(X)))
  for (k in seq_len(n_coef)) {
    r_sum <- Reduce(`+`, internal_r[[k]])
    max_sum_to_zero_error <- max(
      max_sum_to_zero_error, max(abs(r_sum))
    )
    max_sum_to_zero_scale <- max(
      max_sum_to_zero_scale,
      max(vapply(internal_r[[k]], function(x) max(abs(x)), numeric(1)))
    )
  }
  # Construct public total coefficients from separately saved b and r draws.
  eta_internal[] <- 0
  for (k in seq_len(n_coef)) {
    public_b <- draw_variable(fit, paste0("b_", case$group_coef[k]))
    for (j in seq_len(n_group)) {
      public_r <- draw_variable(
        fit,
        paste0(
          "r_", case$group, "[", group_levels[j], ",",
          case$group_coef[k], "]"
        )
      )
      public_total <- public_b + public_r
      internal_total <- theta_raw[[k]] + internal_r[[k]][[j]]
      max_finite_error <- max(
        max_finite_error, max(abs(public_total - internal_total))
      )
      max_finite_scale <- max(
        max_finite_scale, max(abs(public_total)), max(abs(internal_total))
      )
    }
    for (n in seq_len(nrow(X))) {
      eta_internal[, , n] <- eta_internal[, , n] + X[n, k] *
        (theta_raw[[k]] + internal_r[[k]][[group_index[n]]])
    }
  }

  eta_brms <- posterior_linpred(fit, transform = FALSE)
  eta_flat <- matrix(NA_real_, nrow = n_iteration * n_chain, ncol = nrow(X))
  for (n in seq_len(nrow(X))) {
    eta_flat[, n] <- c(eta_internal[, , n])
  }
  if (!identical(dim(eta_brms), dim(eta_flat))) {
    stop("posterior_linpred shape mismatch for ", case$name, call. = FALSE)
  }
  max_predictor_error <- max(abs(eta_brms - eta_flat))
  max_predictor_scale <- max(1, abs(eta_brms), abs(eta_flat))

  # Verify that the standard brms extractors expose the same conventional
  # draws as the saved public parameters. Exact algebraic identities use a
  # prespecified absolute-plus-relative tolerance below; CmdStan runs request
  # round-trip precision so serialization error remains well below this gate.
  extracted_fixef <- fixef(fit, summary = FALSE)
  extracted_ranef <- ranef(
    fit, summary = FALSE, groups = case$group
  )[[case$group]]
  extracted_coef <- coef(
    fit, summary = FALSE, groups = case$group
  )[[case$group]]
  extractor_shape_ok <-
    !is.null(extracted_ranef) && !is.null(extracted_coef) &&
    setequal(colnames(extracted_fixef), case$group_coef) &&
    setequal(dimnames(extracted_ranef)[[2L]], group_levels) &&
    setequal(dimnames(extracted_ranef)[[3L]], case$group_coef) &&
    setequal(dimnames(extracted_coef)[[2L]], group_levels) &&
    setequal(dimnames(extracted_coef)[[3L]], case$group_coef)
  max_fixef_extractor_error <- Inf
  max_ranef_extractor_error <- Inf
  max_coef_extractor_error <- Inf
  max_fixef_extractor_scale <- 1
  max_ranef_extractor_scale <- 1
  max_coef_extractor_scale <- 1
  if (extractor_shape_ok) {
    max_fixef_extractor_error <- 0
    max_ranef_extractor_error <- 0
    max_coef_extractor_error <- 0
    for (k in seq_len(n_coef)) {
      coef_name <- case$group_coef[k]
      raw_b <- c(draw_variable(fit, paste0("b_", coef_name)))
      extractor_b <- extracted_fixef[, coef_name]
      max_fixef_extractor_error <- max(
        max_fixef_extractor_error, max(abs(raw_b - extractor_b))
      )
      max_fixef_extractor_scale <- max(
        max_fixef_extractor_scale, abs(raw_b), abs(extractor_b)
      )
      for (j in seq_len(n_group)) {
        level <- group_levels[j]
        raw_r <- c(draw_variable(
          fit,
          paste0("r_", case$group, "[", level, ",", coef_name, "]")
        ))
        extractor_r <- extracted_ranef[, level, coef_name]
        extractor_total <- extracted_coef[, level, coef_name]
        raw_total <- raw_b + raw_r
        max_ranef_extractor_error <- max(
          max_ranef_extractor_error, max(abs(raw_r - extractor_r))
        )
        max_coef_extractor_error <- max(
          max_coef_extractor_error, max(abs(raw_total - extractor_total))
        )
        max_ranef_extractor_scale <- max(
          max_ranef_extractor_scale, abs(raw_r), abs(extractor_r)
        )
        max_coef_extractor_scale <- max(
          max_coef_extractor_scale, abs(raw_total), abs(extractor_total)
        )
      }
    }
  }

  public_names <- c(
    paste0("b_", case$group_coef),
    as.vector(outer(
      group_levels, case$group_coef,
      function(level, coef) {
        paste0("r_", case$group, "[", level, ",", coef, "]")
      }
    ))
  )
  public_names_ok <- all(public_names %in% variables(fit))
  fixef_names_ok <- setequal(
    colnames(fixef(fit, summary = FALSE)), case$group_coef
  )
  out <- data.frame(
    case = case$name,
    replicate = replicate,
    parameterization = parameterization,
    check = c(
      "saved public b/r names use conventional API",
      "fixef exposes only public coefficient names",
      "fixef draws equal saved public b",
      "ranef draws equal saved public r",
      "coef draws equal saved public b+r",
      "internal group effects sum to zero",
      "public b+r equals internal finite coefficient",
      "posterior_linpred equals internal S2Z predictor"
    ),
    max_abs_error = c(
      if (public_names_ok) 0 else Inf,
      if (fixef_names_ok) 0 else Inf,
      max_fixef_extractor_error,
      max_ranef_extractor_error,
      max_coef_extractor_error,
      max_sum_to_zero_error,
      max_finite_error,
      max_predictor_error
    ),
    comparison_scale = c(
      1,
      1,
      max_fixef_extractor_scale,
      max_ranef_extractor_scale,
      max_coef_extractor_scale,
      max_sum_to_zero_scale,
      max_finite_scale,
      max_predictor_scale
    ),
    absolute_tolerance = invariant_tolerance,
    relative_tolerance = invariant_relative_tolerance,
    stringsAsFactors = FALSE
  )
  out$tolerance <- out$absolute_tolerance +
    out$relative_tolerance * out$comparison_scale
  out$pass <- out$max_abs_error <= out$tolerance
  out
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

sampling_report <- function(fit, case, parameterization, replicate,
                            wall_seconds, seed) {
  public <- grep(
    paste0(
      "^(b_|sd_", case$group, "__|cor_", case$group, "__|df_",
      case$group, "$|sigma$|r_", case$group, "\\[)"
    ),
    variables(fit), value = TRUE
  )
  public <- public[!grepl("^theta_s2z(\\[|$)", public)]
  summary <- summarise_draws(
    as_draws_array(fit, variable = public),
    rhat = rhat,
    ess_bulk = ess_bulk,
    ess_tail = ess_tail
  )
  nuts <- nuts_params(fit)
  leapfrog <- nuts$Value[nuts$Parameter == "n_leapfrog__"]
  if (!length(leapfrog)) {
    stop("n_leapfrog__ was unavailable for ", case$name, call. = FALSE)
  }
  gradient_evals <- sum(leapfrog + 1)
  divergences <- sum(nuts$Value[nuts$Parameter == "divergent__"])
  treedepth <- nuts$Value[nuts$Parameter == "treedepth__"]
  ebfmi <- ebfmi_by_chain(nuts)
  finite_rhat <- summary$rhat[is.finite(summary$rhat)]
  finite_bulk <- summary$ess_bulk[is.finite(summary$ess_bulk)]
  finite_tail <- summary$ess_tail[is.finite(summary$ess_tail)]
  if (!length(finite_rhat) || !length(finite_bulk) || !length(finite_tail)) {
    stop("Posterior diagnostics were not finite for ", case$name,
         call. = FALSE)
  }
  out <- data.frame(
    case = case$name,
    replicate = replicate,
    parameterization = parameterization,
    seed = seed,
    observations = nrow(case$data),
    groups = nlevels(case$data[[case$group]]),
    coefficients_per_group = length(case$group_coef),
    monitored_parameters = nrow(summary),
    postwarmup_transitions = length(leapfrog),
    wall_seconds = wall_seconds,
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
    median_bulk_ess_per_second =
      stats::median(finite_bulk) / wall_seconds,
    median_tail_ess_per_second =
      stats::median(finite_tail) / wall_seconds,
    divergences = divergences,
    treedepth_hits = sum(treedepth >= max_treedepth),
    min_ebfmi = if (all(is.na(ebfmi))) NA_real_ else min(ebfmi, na.rm = TRUE),
    max_rhat = max(finite_rhat),
    stringsAsFactors = FALSE
  )
  out$pass <- out$max_rhat <= max_rhat &&
    out$min_bulk_ess >= min_bulk_ess &&
    out$min_tail_ess >= min_tail_ess &&
    is.finite(out$min_ebfmi) && out$min_ebfmi >= min_ebfmi &&
    out$divergences <= max_divergences &&
    out$treedepth_hits <= max_treedepth_hits
  out
}

sample_model <- function(template, case, parameterization, replicate, seed) {
  cat(
    "Sampling ", case$name, " / ", parameterization, " / replicate ",
    replicate, " (seed ", seed, ") ...\n", sep = ""
  )
  update_args <- list(
    object = template,
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
  if (backend == "cmdstanr") {
    update_args$sig_figs <- output_sig_figs
  }
  timing <- system.time({
    fit <- do.call(update, update_args)
  })
  wall_seconds <- unname(timing[["elapsed"]])
  list(
    fit = fit,
    report = sampling_report(
      fit, case, parameterization, replicate, wall_seconds, seed
    )
  )
}

jobs <- expand.grid(
  case = requested_cases,
  replicate = seq_len(reps),
  stringsAsFactors = FALSE
)
set.seed(base_seed + 90000L)
jobs <- jobs[sample(seq_len(nrow(jobs))), , drop = FALSE]

comparison_rows <- list()
invariant_rows <- list()
quality_rows <- list()
comparison_row <- 0L
invariant_row <- 0L
quality_row <- 0L

reference_quantities <- NULL
if ("eight-schools" %in% requested_cases) {
  reference_quantities <- eight_schools_reference_quantities(
    cases[["eight-schools"]]$reference_path, cases[["eight-schools"]]
  )
}

for (i in seq_len(nrow(jobs))) {
  case <- cases[[jobs$case[i]]]
  replicate <- jobs$replicate[i]
  case_parameterizations <- parameterizations_by_case[[case$name]]
  parameterization_order <- sample(case_parameterizations)
  results <- list()
  for (parameterization in parameterization_order) {
    seed <- base_seed + 1000L * replicate +
      100L * match(case$name, all_case_names) +
      match(parameterization, all_parameterizations)
    results[[parameterization]] <- sample_model(
      templates[[case$name]][[parameterization]], case,
      parameterization, replicate, seed
    )
    quality_row <- quality_row + 1L
    quality_rows[[quality_row]] <- results[[parameterization]]$report
  }

  public <- setNames(lapply(case_parameterizations, function(parameterization) {
    public_quantities(results[[parameterization]]$fit, case)
  }), case_parameterizations)
  comparison_pairs <- if (length(case_parameterizations) >= 2L) {
    utils::combn(case_parameterizations, 2L, simplify = FALSE)
  } else {
    list()
  }
  for (pair in comparison_pairs) {
    comparison_row <- comparison_row + 1L
    comparison_rows[[comparison_row]] <- compare_quantity_sets(
      public[[pair[1L]]], public[[pair[2L]]], case$name,
      paste(pair, collapse = "-vs-"), replicate
    )
  }
  for (parameterization in intersect(
    s2z_parameterizations, case_parameterizations
  )) {
    invariant_row <- invariant_row + 1L
    invariant_rows[[invariant_row]] <- s2z_invariants(
      results[[parameterization]]$fit, case, parameterization, replicate
    )
  }

  if (case$name == "eight-schools") {
    for (parameterization in case_parameterizations) {
      comparison_row <- comparison_row + 1L
      comparison_rows[[comparison_row]] <- compare_quantity_sets(
        reference_quantities, public[[parameterization]], case$name,
        paste0("reference-vs-", parameterization), replicate
      )
    }
  }
  rm(results, public)
  invisible(gc())
}

comparison_results <- do.call(rbind, comparison_rows)
comparison_results <- comparison_results[
  order(comparison_results$z_score, decreasing = TRUE),
]
invariant_results <- do.call(rbind, invariant_rows)
quality_results <- do.call(rbind, quality_rows)
quality_results <- quality_results[
  order(
    match(quality_results$case, requested_cases), quality_results$replicate,
    match(quality_results$parameterization, parameterizations)
  ),
]

conventional_quality <- subset(
  quality_results, parameterization == "conventional",
  select = -parameterization
)
efficiency_results <- do.call(rbind, lapply(
  s2z_parameterizations, function(s2z_parameterization) {
    s2z_quality <- subset(
      quality_results, parameterization == s2z_parameterization,
      select = -parameterization
    )
    paired <- merge(
      conventional_quality, s2z_quality,
      by = c(
        "case", "replicate", "observations", "groups",
        "coefficients_per_group"
      ),
      suffixes = c("_conventional", "_s2z"), sort = FALSE
    )
    data.frame(
      case = paired$case,
      replicate = paired$replicate,
      s2z_parameterization = s2z_parameterization,
      observations = paired$observations,
      groups = paired$groups,
      coefficients_per_group = paired$coefficients_per_group,
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
      tail_ess_per_second_gain_s2z_over_conventional =
        paired$median_tail_ess_per_second_s2z /
          paired$median_tail_ess_per_second_conventional,
      stringsAsFactors = FALSE
    )
  }
))

utils::write.csv(
  comparison_results,
  file.path(output_dir, "s2z-posteriordb-equivalence.csv"), row.names = FALSE
)
utils::write.csv(
  invariant_results,
  file.path(output_dir, "s2z-posteriordb-invariants.csv"), row.names = FALSE
)
utils::write.csv(
  quality_results,
  file.path(output_dir, "s2z-posteriordb-quality.csv"), row.names = FALSE
)
utils::write.csv(
  efficiency_results,
  file.path(output_dir, "s2z-posteriordb-efficiency.csv"), row.names = FALSE
)

cat("\nMCSE- and margin-aware posterior comparisons (largest 40 z-scores)\n")
print(utils::head(comparison_results, 40L), row.names = FALSE, digits = 5)
cat("\nExact S2Z coordinate invariants\n")
print(invariant_results, row.names = FALSE, digits = 5)
cat("\nSampling quality and raw efficiency\n")
print(quality_results, row.names = FALSE, digits = 5)
cat("\nPaired efficiency ratios by S2Z mode (> 1 favors S2Z)\n")
print(efficiency_results, row.names = FALSE, digits = 5)
cat("\nFull outcomes written to ", normalizePath(output_dir), "\n", sep = "")

failed_comparisons <- subset(comparison_results, !pass)
failed_invariants <- subset(invariant_results, !pass)
failed_quality <- subset(quality_results, !pass)
if (nrow(failed_comparisons) || nrow(failed_invariants) ||
    nrow(failed_quality)) {
  stop(
    "PosteriorDB S2Z validation failed: ",
    nrow(failed_comparisons), " posterior comparisons, ",
    nrow(failed_invariants), " exact invariants, and ",
    nrow(failed_quality), " sampling-quality checks failed.",
    call. = FALSE
  )
}

cat(
  "\nPASS: all selected conventional and S2Z real-model posteriors agree ",
  "pairwise for means, SDs, and quantiles within propagated MCSE and ",
  "prespecified margins; ",
  if ("eight-schools" %in% requested_cases) {
    "eight schools also agrees with its official PosteriorDB reference; "
  } else {
    ""
  },
  "all coordinate and sampler checks pass.\n",
  sep = ""
)
