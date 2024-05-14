#' @title Stan data for Bayesian models
#'
#' @description \code{standata} is a generic function that can be used to
#'   generate data for Bayesian models to be passed to Stan. Its original use is
#'   within the \pkg{brms} package, but new methods for use
#'   with objects from other packages can be registered to the same generic.
#'
#' @param object A formula object whose class will determine which method will
#'   be used. A symbolic description of the model to be fitted.
#' @param formula Synonym of \code{object} for use in \code{make_standata}.
#' @param ... Further arguments passed to the specific method.
#'
#' @return A named list of objects containing the required data to fit a
#'   Bayesian model with \pkg{Stan}.
#'
#' @details
#' See \code{\link{standata.default}} for the default method applied for
#' \pkg{brms} models. You can view the available methods by typing
#' \code{methods(standata)}. The \code{make_standata} function is an alias
#' of \code{standata}.
#'
#' @examples
#' sdata1 <- standata(rating ~ treat + period + carry + (1|subject),
#'                    data = inhaler, family = "cumulative")
#' str(sdata1)
#'
#' @seealso
#'   \code{\link{standata.default}}, \code{\link{standata.brmsfit}}
#'
#' @export
standata <- function(object, ...) {
  UseMethod("standata")
}

#' @rdname standata
#' @export
make_standata <- function(formula, ...) {
  # became an alias of standata in 2.20.14.
  standata(formula, ...)
}

#' Data for \pkg{brms} Models
#'
#' Generate data for \pkg{brms} models to be passed to \pkg{Stan}.
#'
#' @inheritParams brm
#' @param object An object of class \code{\link[stats:formula]{formula}},
#'   \code{\link{brmsformula}}, or \code{\link{mvbrmsformula}} (or one that can
#'   be coerced to that classes): A symbolic description of the model to be
#'   fitted. The details of model specification are explained in
#'   \code{\link{brmsformula}}.
#' @param ... Other arguments for internal use.
#'
#' @return A named list of objects containing the required data
#'   to fit a \pkg{brms} model with \pkg{Stan}.
#'
#' @examples
#' sdata1 <- standata(rating ~ treat + period + carry + (1|subject),
#'                    data = inhaler, family = "cumulative")
#' str(sdata1)
#'
#' sdata2 <- standata(count ~ zAge + zBase * Trt + (1|patient),
#'                    data = epilepsy, family = "poisson")
#' str(sdata2)
#'
#' @export
standata.default <- function(object, data, family = gaussian(), prior = NULL,
                             autocor = NULL, data2 = NULL, cov_ranef = NULL,
                             sample_prior = "no", stanvars = NULL,
                             threads = getOption("brms.threads", NULL),
                             knots = NULL, drop_unused_levels = TRUE, ...) {

  object <- validate_formula(
    object, data = data, family = family,
    autocor = autocor, cov_ranef = cov_ranef
  )
  bterms <- brmsterms(object)
  data2 <- validate_data2(
    data2, bterms = bterms,
    get_data2_autocor(object),
    get_data2_cov_ranef(object)
  )
  data <- validate_data(
    data, bterms = bterms,
    knots = knots, data2 = data2,
    drop_unused_levels = drop_unused_levels
  )
  bframe <- brmsframe(bterms, data)
  prior <- .validate_prior(
    prior, bframe = bframe,
    sample_prior = sample_prior
  )
  stanvars <- validate_stanvars(stanvars)
  threads <- validate_threads(threads)
  .standata(
    bframe, data = data, prior = prior,
    data2 = data2, stanvars = stanvars,
    threads = threads, ...
  )
}

# internal work function of 'standata'
# @param check_response check validity of the response?
# @param only_response extract data related to the response only?
# @param internal prepare Stan data for use in post-processing methods?
# @param basis original Stan data as prepared by 'frame_basis'
# @param ... currently ignored
# @return names list of data passed to Stan
.standata <- function(bframe, data, prior, stanvars, data2,
                      threads = threading(), check_response = TRUE,
                      only_response = FALSE, internal = FALSE, ...) {

  stopifnot(is.anybrmsframe(bframe))
  check_response <- as_one_logical(check_response)
  only_response <- as_one_logical(only_response)
  internal <- as_one_logical(internal)
  # order data for use in autocorrelation models
  data <- order_data(data, bterms = bframe)
  out <- data_response(
    bframe, data, check_response = check_response,
    internal = internal
  )
  if (!only_response) {
    # pass as sdata so that data_special_prior knows about data_gr_global
    # TODO: compute sdata_gr_global in brmsframe in brms 3.0
    # this would require passing data2 to brmsframe
    sdata_gr_global <- data_gr_global(bframe, data2 = data2)
    c(out) <- data_predictor(
      bframe, data = data, prior = prior, data2 = data2,
      sdata = sdata_gr_global
    )
    c(out) <- sdata_gr_global
    c(out) <- data_Xme(bframe, data = data)
  }
  out$prior_only <- as.integer(is_prior_only(prior))
  if (use_threading(threads)) {
    out$grainsize <- threads$grainsize
    if (is.null(out$grainsize)) {
      out$grainsize <- ceiling(out$N / (2 * threads$threads))
      out$grainsize <- max(100, out$grainsize)
    }
  }
  if (is.stanvars(stanvars)) {
    stanvars <- subset_stanvars(stanvars, block = "data")
    inv_names <- intersect(names(stanvars), names(out))
    if (length(inv_names)) {
      stop2("Cannot overwrite existing variables: ",
            collapse_comma(inv_names))
    }
    out[names(stanvars)] <- from_list(stanvars, "sdata")
  }
  if (internal) {
    # allows to recover the original order of the data
    attr(out, "old_order") <- attr(data, "old_order")
    # ensures currently used grouping levels are known in post-processing
    set_levels(out, "used") <- get_levels(bframe, prefix = "used")
  }
  structure(out, class = c("standata", "list"))
}

#' Extract data passed to Stan from \code{brmsfit} objects
#'
#' Extract all data that was used by Stan to fit a \pkg{brms} model.
#'
#' @param object An object of class \code{brmsfit}.
#' @param ... More arguments passed to
#'   \code{\link[brms:standata.default]{standata.default}}.
#'   and \code{\link{validate_newdata}}.
#' @inheritParams prepare_predictions
#'
#' @return A named list containing the data passed to Stan.
#'
#' @export
standata.brmsfit <- function(object, newdata = NULL, re_formula = NULL,
                             newdata2 = NULL, new_objects = NULL,
                             incl_autocor = TRUE, ...) {

  # allows functions to fall back to old default behavior
  # which was used when originally fitting the model
  options(.brmsfit_version = object$version$brms)
  on.exit(options(.brmsfit_version = NULL))

  object <- exclude_terms(object, incl_autocor = incl_autocor)
  formula <- update_re_terms(object$formula, re_formula)
  bterms <- brmsterms(formula)

  newdata2 <- use_alias(newdata2, new_objects)
  data2 <- current_data2(object, newdata2)
  data <- current_data(
    object, newdata, newdata2 = data2,
    re_formula = re_formula, ...
  )
  stanvars <- add_newdata_stanvars(object$stanvars, data2)

  basis <- object$basis
  if (is.null(basis)) {
    # this case should not happen actually, perhaps when people use
    # the 'empty' feature. But computing it here will be fine
    # for almost all models, only causing potential problems for processing
    # of splines on new machines (#1465)
    basis <- frame_basis(bterms, data = object$data)
  }
  bframe <- brmsframe(bterms, data = data, basis = basis)
  .standata(
    bframe, data = data, prior = object$prior,
    data2 = data2, stanvars = stanvars,
    threads = object$threads, ...
  )
}
