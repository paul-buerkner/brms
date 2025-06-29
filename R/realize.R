#' Realize a brm Call from Preview
#'
#' Converts a `brm_call_preview` object created by `brm(..., preview = TRUE)` into a full model fit.
#' This function evaluates the original call with `preview = FALSE`, allowing the user to "realize" the model
#' when they are ready to run it for real.
#'
#' @param x An object of class `brm_call_preview`, typically returned by calling `brm(..., preview = TRUE)`.
#' @param ... Optional arguments to override or supplement the original call. These arguments
#'   will be merged into the original call before evaluation.
#'
#' @details
#' This function is useful in interactive workflows or debugging sessions where you first want
#' to preview what a model call would look like, and then later decide to run it. It ensures
#' that `preview = FALSE` is set in the final call. Additional parameters can be supplied
#' to change or extend the original model fitting request.
#'
#' @return The result of the evaluated `brm()` call (typically a `brmsfit` object).
#'
#' @seealso [brms::brm()], [print.brm_call_preview()]
#'
#' @examples
#' \dontrun{
#' # Preview model setup
#' preview_obj <- brm(
#'   count ~ zAge + zBase * Trt + (1|patient),
#'   data = epilepsy, family = poisson(), preview = TRUE
#' )
#'
#' # Realize the full model fit
#' fit <- realize(preview_obj)
#'
#' # Override any parameter from the original call if needed
#' fit2 <- realize(preview_obj, family = gaussian())
#' }
#'
#' @export

realize <- function(x, ...) {
  if (!inherits(x, "brm_call_preview")) {
    stop("Object must be of class 'brm_call_preview'")
  }

  # Retrieve the original call
  call <- x$call
  eval_env <- parent.frame()

  # Turn the call into a brm call and expand dots
  call <- match.call(definition = brm, call = call, expand.dots = TRUE)

  # Convert to list, update preview to FALSE, merge additional args
  call_list <- as.list(call)
  call_list$preview <- FALSE

  # Let user override or add parameters via ...
  call_list <- modifyList(call_list, list(...))

  # Convert back to call
  call <- as.call(call_list)

  # Evaluate
  eval(call, envir = eval_env)
}
