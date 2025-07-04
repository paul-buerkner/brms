# Store original option value
original_cache_folder <- getOption("brms.cache_folder")

# Set a temporary value for tests
test_cache_dir <- tempdir()
options(brms.cache_folder = test_cache_dir)

# Optional: only show message in interactive mode
if (interactive()) {
  message(glue::glue("Temporarily setting option(brms.cache_folder = '{test_cache_dir}')"))
}

# Restore the original option after tests
withr::defer(
  {
    if (is.null(original_cache_folder)) {
      options(brms.cache_folder = NULL)
    } else {
      options(brms.cache_folder = original_cache_folder)
    }
  },
  teardown_env()
)
