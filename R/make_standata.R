#' Data for \pkg{brms} Models
#' 
#' Generate data for \pkg{brms} models to be passed to \pkg{Stan}
#'
#' @inheritParams brm
#' @param ... Other arguments for internal use.
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
make_standata <- function(formula, data, family = gaussian(), prior = NULL, 
                          autocor = NULL, data2 = NULL, cov_ranef = NULL, 
                          sample_prior = "no", stanvars = NULL, knots = NULL, 
                          ...) {

  if (is.brmsfit(formula)) {
    stop2("Use 'standata' to extract Stan data from 'brmsfit' objects.")
  }
  formula <- validate_formula(
    formula, data = data, family = family, 
    autocor = autocor, cov_ranef = cov_ranef
  )
  bterms <- brmsterms(formula)
  data <- validate_data(data, bterms = bterms, knots = knots)
  prior <- validate_prior(
    prior, bterms = bterms, data = data,
    sample_prior = sample_prior,
    require_nlpar_prior = FALSE
  )
  data2 <- validate_data2(
    data2, bterms = bterms, 
    get_data2_autocor(formula),
    get_data2_cov_ranef(formula)
  )
  stanvars <- validate_stanvars(stanvars)
  .make_standata(
    bterms, data = data, prior = prior,
    data2 = data2, stanvars = stanvars,
    ...
  )
}

# internal work function of 'make_stancode'
# @param check_response check validity of the response?
# @param only_response extract data related to the response only?
# @param internal prepare Stan data for use in post-processing methods?
# @param old_sdata original Stan data as prepared by 'extract_old_standata'
# @param old_levels named list of original group levels
# @param ... currently ignored
# @return names list of data passed to Stan
.make_standata <- function(bterms, data, prior, stanvars, data2, 
                           check_response = TRUE, only_response = FALSE, 
                           internal = FALSE, old_sdata = NULL,
                           old_levels = NULL, ...) {
  
  check_response <- as_one_logical(check_response)
  only_response <- as_one_logical(only_response)
  internal <- as_one_logical(internal)
  # order data for use in autocorrelation models
  data <- order_data(data, bterms = bterms)
  out <- data_response(
    bterms, data, check_response = check_response,
    internal = internal, old_sdata = old_sdata
  )
  if (!only_response) {
    ranef <- tidy_ranef(bterms, data, old_levels = old_levels)
    meef <- tidy_meef(bterms, data, old_levels = old_levels)
    c(out) <- data_predictor(
      bterms, data = data, prior = prior, data2 = data2,
      ranef = ranef, old_sdata = old_sdata
    )
    c(out) <- data_gr_global(ranef, data2 = data2)
    c(out) <- data_Xme(meef, data = data)
  }
  out$prior_only <- as.integer(is_equal(get_sample_prior(prior), "only"))
  if (is.stanvars(stanvars)) {
    stanvars <- subset_stanvars(stanvars, block = "data")
    inv_names <- intersect(names(stanvars), names(out))
    if (length(inv_names)) {
      stop2("Cannot overwrite existing variables: ", 
            collapse_comma(inv_names))
    }
    out[names(stanvars)] <- lapply(stanvars, "[[", "sdata")
  }
  if (internal) {
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
#' @param ... More arguments passed to \code{\link{make_standata}}
#'   and \code{\link{validate_newdata}}.
#' @inheritParams extract_draws
#' 
#' @return A named list containing the data originally passed to Stan.
#' 
#' @export
standata.brmsfit <- function(object, newdata = NULL, re_formula = NULL, 
                             newdata2 = NULL, new_objects = NULL,
                             incl_autocor = TRUE, ...) {
  object <- restructure(object)
  object <- exclude_terms(object, incl_autocor = incl_autocor)
  version <- object$version$brms
  newdata2 <- use_alias(newdata2, new_objects)
  formula <- update_re_terms(object$formula, re_formula)
  bterms <- brmsterms(formula)
  data <- current_data(object, newdata, re_formula = re_formula, ...)
  
  old_sdata <- old_levels <- NULL
  if (is.null(newdata)) {
    if (version <= "2.8.6" && has_smooths(bterms)) {
      # the spline penality has changed in 2.8.7 (#646)
      old_sdata <- extract_old_standata(
        bterms, data = object$data, version = version
      )
    }
  } else {
    old_levels <- get_levels(
      tidy_ranef(bterms, object$data),
      tidy_meef(bterms, object$data)
    )
    old_sdata <- extract_old_standata(
      bterms, data = object$data, version = version
    )
  }
  stanvars <- object$stanvars
  if (is.null(newdata2)) {
    data2 <- object$data2
  } else {
    data2 <- validate_data2(newdata2, bterms = bterms)
    stanvars <- add_newdata_stanvars(stanvars, data2)
  }
  
  .make_standata(
    bterms, data = data, prior = object$prior,
    data2 = data2, stanvars = stanvars, 
    old_sdata = old_sdata, old_levels = old_levels,
    ...
  )
}

#' @rdname standata.brmsfit
#' @export
standata <- function(object, ...) {
  UseMethod("standata")
}
