#' @details By default, this method uses \code{sample_new_levels = "gaussian"}
#'   to sample parameter values for new grouping-factor levels (see also
#'   \code{\link{prepare_predictions}}). This default will fail for models with
#'   non-Gaussian group-level effects. In this case, we recommend setting
#'   \code{sample_new_levels = "uncertainty"}.
