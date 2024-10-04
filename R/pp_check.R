#' Posterior Predictive Checks for \code{brmsfit} Objects
#'
#' Perform posterior predictive checks with the help
#' of the \pkg{bayesplot} package.
#'
#' @aliases pp_check
#'
#' @param object An object of class \code{brmsfit}.
#' @param type Type of the ppc plot as given by a character string.
#'   See \code{\link[bayesplot:PPC-overview]{PPC}} for an overview
#'   of currently supported types. You may also use an invalid
#'   type (e.g. \code{type = "xyz"}) to get a list of supported
#'   types in the resulting error message.
#' @param ndraws Positive integer indicating how many
#'  posterior draws should be used.
#'  If \code{NULL} all draws are used. If not specified,
#'  the number of posterior draws is chosen automatically.
#'  Ignored if \code{draw_ids} is not \code{NULL}.
#' @param prefix The prefix of the \pkg{bayesplot} function to be applied.
#'  Either `"ppc"` (posterior predictive check; the default)
#'  or `"ppd"` (posterior predictive distribution), the latter being the same
#'  as the former except that the observed data is not shown for `"ppd"`.
#' @param group Optional name of a factor variable in the model
#'  by which to stratify the ppc plot. This argument is required for
#'  ppc \code{*_grouped} types and ignored otherwise.
#' @param x Optional name of a variable in the model.
#'  Only used for ppc types having an \code{x} argument
#'  and ignored otherwise.
#' @param ... Further arguments passed to \code{\link{predict.brmsfit}}
#'   as well as to the PPC function specified in \code{type}.
#' @inheritParams prepare_predictions.brmsfit
#'
#' @return A ggplot object that can be further
#'  customized using the \pkg{ggplot2} package.
#'
#' @details For a detailed explanation of each of the ppc functions,
#' see the \code{\link[bayesplot:PPC-overview]{PPC}}
#' documentation of the \pkg{\link[bayesplot:bayesplot-package]{bayesplot}}
#' package.
#'
#' @examples
#' \dontrun{
#' fit <-  brm(count ~ zAge + zBase * Trt
#'             + (1|patient) + (1|obs),
#'             data = epilepsy, family = poisson())
#'
#' pp_check(fit)  # shows dens_overlay plot by default
#' pp_check(fit, type = "error_hist", ndraws = 11)
#' pp_check(fit, type = "scatter_avg", ndraws = 100)
#' pp_check(fit, type = "stat_2d")
#' pp_check(fit, type = "rootogram")
#' pp_check(fit, type = "loo_pit")
#'
#' ## get an overview of all valid types
#' pp_check(fit, type = "xyz")
#'
#' ## get a plot without the observed data
#' pp_check(fit, prefix = "ppd")
#' }
#'
#' @importFrom bayesplot pp_check
#' @export pp_check
#' @export
pp_check.brmsfit <- function(object, type, ndraws = NULL, prefix = c("ppc", "ppd"),
                             group = NULL, x = NULL, newdata = NULL, resp = NULL,
                             draw_ids = NULL, nsamples = NULL, subset = NULL, ...) {
  dots <- list(...)
  if (missing(type)) {
    type <- "dens_overlay"
  }
  type <- as_one_character(type)
  prefix <- match.arg(prefix)
  if (!is.null(group)) {
    group <- as_one_character(group)
  }
  if (!is.null(x)) {
    x <- as_one_character(x)
  }
  ndraws_given <- any(c("ndraws", "nsamples") %in% names(match.call()))
  ndraws <- use_alias(ndraws, nsamples)
  draw_ids <- use_alias(draw_ids, subset)
  resp <- validate_resp(resp, object, multiple = FALSE)
  if (prefix == "ppc") {
    # no type checking for prefix 'ppd' yet
    valid_types <- as.character(bayesplot::available_ppc(""))
    valid_types <- sub("^ppc_", "", valid_types)
    if (!type %in% valid_types) {
      stop2("Type '", type, "' is not a valid ppc type. ",
            "Valid types are:\n", collapse_comma(valid_types))
    }
  }
  ppc_fun <- get(paste0(prefix, "_", type), asNamespace("bayesplot"))

  object <- restructure(object)
  stopifnot_resp(object, resp)
  family <- family(object, resp = resp)
  if (has_multicol(family)) {
    stop2("'pp_check' is not implemented for this family.")
  }
  valid_vars <- names(model.frame(object))
  if ("group" %in% names(formals(ppc_fun))) {
    if (is.null(group)) {
      stop2("Argument 'group' is required for ppc type '", type, "'.")
    }
    if (!group %in% valid_vars) {
      stop2("Variable '", group, "' could not be found in the data.")
    }
  }
  if ("x" %in% names(formals(ppc_fun))) {
    if (!is.null(x) && !x %in% valid_vars) {
      stop2("Variable '", x, "' could not be found in the data.")
    }
  }
  if (type == "error_binned") {
    if (is_polytomous(family)) {
      stop2("Type '", type, "' is not available for polytomous models.")
    }
    method <- "posterior_epred"
  } else {
    method <- "posterior_predict"
  }
  if (!ndraws_given) {
    aps_types <- c(
      "error_scatter_avg", "error_scatter_avg_vs_x",
      "intervals", "intervals_grouped",
      "loo_intervals", "loo_pit", "loo_pit_overlay",
      "loo_pit_qq", "loo_ribbon",
      'pit_ecdf', 'pit_ecdf_grouped',
      "ribbon", "ribbon_grouped",
      "rootogram", "scatter_avg", "scatter_avg_grouped",
      "stat", "stat_2d", "stat_freqpoly_grouped", "stat_grouped",
      "violin_grouped"
    )
    if (!is.null(draw_ids)) {
      ndraws <- NULL
    } else if (type %in% aps_types) {
      ndraws <- NULL
      message("Using all posterior draws for ppc type '",
              type, "' by default.")
    } else {
      ndraws <- 10
      message("Using 10 posterior draws for ppc type '",
              type, "' by default.")
    }
  }

  y <- NULL
  if (prefix == "ppc") {
    # y is ignored in prefix 'ppd' plots
    y <- get_y(object, resp = resp, newdata = newdata, ...)
  }
  draw_ids <- validate_draw_ids(object, draw_ids, ndraws)
  pred_args <- list(
    object, newdata = newdata, resp = resp,
    draw_ids = draw_ids, ...
  )
  yrep <- do_call(method, pred_args)

  if (anyNA(y)) {
    warning2("NA responses are not shown in 'pp_check'.")
    take <- !is.na(y)
    y <- y[take]
    yrep <- yrep[, take, drop = FALSE]
  }

  data <- current_data(
    object, newdata = newdata, resp = resp,
    re_formula = NA, check_response = TRUE, ...
  )

  # prepare plotting arguments
  ppc_args <- list()
  if (prefix == "ppc") {
    ppc_args$y <- y
    ppc_args$yrep <- yrep
  } else if (prefix == "ppd") {
    ppc_args$ypred <- yrep
  }
  if (!is.null(group)) {
    ppc_args$group <- data[[group]]
  }
  if (!is.null(x)) {
    ppc_args$x <- data[[x]]
    if (!is_like_factor(ppc_args$x)) {
      ppc_args$x <- as.numeric(ppc_args$x)
    }
  }
  if ("lw" %in% setdiff(names(formals(ppc_fun)), names(ppc_args))) {
    # run loo instead of psis to allow for moment matching
    loo_object <- do_call(loo, c(pred_args, save_psis = TRUE))
    ppc_args$lw <- weights(loo_object$psis_object, log = TRUE)
  } else if ("psis_object" %in% setdiff(names(formals(ppc_fun)), names(ppc_args))) {
    # some PPCs may only support 'psis_object' but not 'lw' for whatever reason
    loo_object <- do_call(loo, c(pred_args, save_psis = TRUE))
    ppc_args$psis_object <- loo_object$psis_object
  }

  # censored responses are misleading when displayed in pp_check
  bterms <- brmsterms(object$formula)
  cens <- get_cens(bterms, data, resp = resp)
  is_censoring_type <- type %in% c("km_overlay", "km_overlay_grouped")
  if (!is.null(cens)) {
    if (is_censoring_type) {
      if (any(cens %in% c(-1, 2))) {
        warning2("Left and interval censored responses are not included.")
      }
      take <- cens %in% c(0, 1)
      ppc_args$status_y <- 1 - cens[take]
    } else {
      warning2("Censored responses are not included.")
      take <- !cens
    }
    if (!any(take)) {
      stop2("No valid responses found to include.")
    }
    ppc_args$y <- ppc_args$y[take]
    ppc_args$yrep <- ppc_args$yrep[, take, drop = FALSE]
    if (!is.null(ppc_args$group)) {
      ppc_args$group <- ppc_args$group[take]
    }
    if (!is.null(ppc_args$x)) {
      ppc_args$x <- ppc_args$x[take]
    }
    if (!is.null(ppc_args$lw)) {
      ppc_args$lw <- ppc_args$lw[, take]
    } else if (!is.null(ppc_args$psis_object)) {
      ppc_args$psis_object <- subset(ppc_args$psis_object, take)
    }
  } else if (is_censoring_type) {
    # status_y is mandatory for some ppc types
    ppc_args$status_y <- rep(1, length(ppc_args$y))
  }

  # most ... arguments are meant for the prediction function
  for_pred <- names(dots) %in% names(formals(prepare_predictions.brmsfit))
  ppc_args <- c(ppc_args, dots[!for_pred])

  do_call(ppc_fun, ppc_args)
}
