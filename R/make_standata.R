#' Data for \pkg{brms} Models
#' 
#' Generate data for \pkg{brms} models to be passed to \pkg{Stan}
#'
#' @inheritParams brm
#' @param check_response Logical; check validity of the response?
#' @param only_response Logical; extract data related to the response only?
#' @param control A named list currently for internal usage only
#' @param ... Other potential arguments
#' 
#' @return A named list of objects containing the required data 
#'   to fit a \pkg{brms} model with \pkg{Stan}. 
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' data1 <- make_standata(rating ~ treat + period + carry + (1|subject), 
#'                        data = inhaler, family = "cumulative")
#' names(data1)
#' 
#' data2 <- make_standata(count ~ zAge + zBase * Trt + (1|patient),
#'                        data = epilepsy, family = "poisson")
#' names(data2)
#'          
#' @export
make_standata <- function(formula, data, family = gaussian(), 
                          prior = NULL, data2 = NULL, 
                          autocor = NULL, cov_ranef = NULL, 
                          sample_prior = c("no", "yes", "only"), 
                          stanvars = NULL, knots = NULL, 
                          check_response = TRUE, only_response = FALSE, 
                          control = list(), ...) {
  # control arguments:
  #   internal: is make_standata called for internal use in S3 methods?
  #   new: is make_standata is called with new data?
  #   save_order: should the initial order of the data be saved?
  #   old_sdata: list of stan data computed from the orginal data
  #   terms_attr: list of attributes of the original model.frame
  dots <- list(...)
  # some input checks
  if (is.brmsfit(formula)) {
    stop2("Use 'standata' to extract Stan data from 'brmsfit' objects.")
  }
  check_response <- as_one_logical(check_response)
  only_response <- as_one_logical(only_response)
  internal <- isTRUE(control$internal)
  new <- isTRUE(control$new)
  formula <- validate_formula(
    formula, data = data, family = family, autocor = autocor
  )
  bterms <- brmsterms(formula)
  sample_prior <- check_sample_prior(sample_prior)
  check_prior_content(prior, warn = FALSE)
  prior <- check_prior_special(
    prior, bterms = bterms, data = data, 
    check_nlpar_prior = FALSE
  )
  na_action <- if (new) na.pass else na.omit2
  data <- validate_data(
    data, bterms = bterms, na.action = na_action, 
    drop.unused.levels = !new, knots = knots,
    terms_attr = control$terms_attr
  )
  # order data for use in autocorrelation models
  data <- order_data(data, bterms = bterms)
  data2 <- validate_data2(
    data2, bterms = bterms, 
    get_data2_autocor(formula)
  )
  
  out <- data_response(
    bterms, data, check_response = check_response,
    internal = internal, old_sdata = control$old_sdata
  )
  if (!only_response) {
    ranef <- tidy_ranef(bterms, data, old_levels = control$old_levels)
    c(out) <- data_predictor(
      bterms, data = data, prior = prior, data2 = data2,
      ranef = ranef, old_sdata = control$old_sdata
    )
    c(out) <- data_gr_global(ranef, cov_ranef = cov_ranef, internal = internal)
    meef <- tidy_meef(bterms, data, old_levels = control$old_levels)
    c(out) <- data_Xme(meef, data = data)
  }
  out$prior_only <- as.integer(identical(sample_prior, "only"))
  stanvars <- validate_stanvars(stanvars)
  if (is.stanvars(stanvars)) {
    stanvars <- subset_stanvars(stanvars, block = "data")
    inv_names <- intersect(names(stanvars), names(out))
    if (length(inv_names)) {
      stop2("Cannot overwrite existing variables: ", 
            collapse_comma(inv_names))
    }
    out[names(stanvars)] <- lapply(stanvars, "[[", "sdata")
  }
  if (isTRUE(control$save_order)) {
    attr(out, "old_order") <- attr(data, "old_order")
  }
  structure(out, class = "standata")
}

#' Extract data passed to Stan
#' 
#' Extract all data that was used by Stan to fit the model.
#' 
#' @aliases standata.brmsfit
#' 
#' @param object An object of class \code{brmsfit}.
#' @param internal Logical, indicates if the data should be prepared 
#'   for internal use in other post-processing methods.
#' @param control A named list currently for internal usage only.
#' @param ... More arguments passed to \code{\link{make_standata}}.
#' @inheritParams extract_draws
#' 
#' @return A named list containing the data originally passed to Stan.
#' 
#' @export
standata.brmsfit <- function(object, newdata = NULL, re_formula = NULL, 
                             newdata2 = NULL, new_objects = NULL,
                             incl_autocor = TRUE, internal = FALSE,
                             control = list(), ...) {
  object <- restructure(object)
  object <- exclude_terms(object, incl_autocor = incl_autocor)
  newdata2 <- use_alias(newdata2, new_objects)
  internal <- as_one_logical(internal)
  is_old_data <- isTRUE(attr(newdata, "old"))
  if (is.null(newdata)) {
    newdata <- object$data
    is_old_data <- TRUE
  }
  if (is.null(newdata2)) {
    newdata2 <- object$data2
  }
  new_formula <- update_re_terms(object$formula, re_formula)
  bterms <- brmsterms(new_formula)
  newdata2 <- validate_data2(newdata2, bterms = bterms)
  version <- object$version$brms
  if (is_old_data) {
    if (version <= "2.8.6" && has_smooths(bterms)) {
      # the spline penality has changed in 2.8.7 (#646)
      control$old_sdata <- extract_old_standata(
        bterms, data = object$data, version = version
      )
    }
  } else {
    if (!isTRUE(attr(newdata, "valid"))) {
      newdata <- validate_newdata(
        newdata, object, re_formula = re_formula, ...
      )
    }
    object <- add_new_stanvars(object, newdata2)
    control$new <- TRUE
    # ensure correct handling of functions like poly or scale
    old_terms <- attr(object$data, "terms")
    terms_attr <- c("variables", "predvars")
    control$terms_attr <- attributes(old_terms)[terms_attr]
    control$old_sdata <- extract_old_standata(
      bterms, data = object$data, version = version
    )
    control$old_levels <- get_levels(
      tidy_ranef(bterms, object$data),
      tidy_meef(bterms, object$data)
    )
  }
  if (internal) {
    control$internal <- TRUE
    control$save_order <- TRUE
  }
  sample_prior <- attr(object$prior, "sample_prior")
  make_standata(
    formula = new_formula, data = newdata, 
    prior = object$prior, data2 = newdata2,
    cov_ranef = object$cov_ranef, 
    sample_prior = sample_prior, 
    stanvars = object$stanvars, 
    control = control, ...
  )
}

#' @rdname standata.brmsfit
#' @export
standata <- function(object, ...) {
  UseMethod("standata")
}
