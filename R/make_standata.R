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
    # allows to recover the original order of the data
    attr(out, "old_order") <- attr(data, "old_order")
    # ensures current grouping levels are known in post-processing
    attr(out, "levels") <- get_levels(
      tidy_meef(bterms, data), tidy_ranef(bterms, data)
    )
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
    # TODO: move old_levels inside old_sdata
    old_levels <- get_levels(
      tidy_meef(bterms, object$data),
      tidy_ranef(bterms, object$data)
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

# TODO: refactor preparation and storage of old standata
# helper function for validate_newdata to extract
# old standata required for the computation of new standata
extract_old_standata <- function(x, data, ...) {
  UseMethod("extract_old_standata")
}

#' @export
extract_old_standata.default <- function(x, data, ...) {
  NULL
}

#' @export
extract_old_standata.mvbrmsterms <- function(x, data, ...) {
  out <- named_list(names(x$responses))
  for (i in seq_along(out)) {
    out[[i]] <- extract_old_standata(x$terms[[i]], data, ...)
  }
  out
}

#' @export
extract_old_standata.brmsterms <- function(x, data, ...) {
  out <- named_list(c(names(x$dpars), names(x$nlpars)))
  data <- subset_data(data, x)
  for (dp in names(x$dpars)) {
    out[[dp]] <- extract_old_standata(x$dpars[[dp]], data, ...)
  }
  for (nlp in names(x$nlpars)) {
    out[[nlp]] <- extract_old_standata(x$nlpars[[nlp]], data, ...)
  }
  if (has_trials(x$family)) {
    # trials should not be computed based on new data
    datr <- data_response(x, data, check_response = FALSE, internal = TRUE)
    # partially match via $ to be independent of the response suffix
    out$trials <- datr$trials
  }
  if (is_binary(x$family) || is_categorical(x$family)) {
    Y <- model.response(model.frame(x$respform, data, na.action = na.pass))
    out$resp_levels <- levels(as.factor(Y))
  }
  if (is_cox(x$family)) {
    # compute basis matrix of the baseline hazard for the Cox model
    datr <- data_response(x, data, check_response = FALSE, internal = TRUE)
    out$bhaz_basis <- bhaz_basis_matrix(datr$Y, args = x$family$bhaz)
  }
  out
}

#' @export
extract_old_standata.btnl <- function(x, data, ...) {
  NULL
}

#' @export
extract_old_standata.btl <- function(x, data, ...) {
  out <- list()
  out$smooths <- make_sm_list(x, data, ...)
  out$gps <- make_gp_list(x, data, ...)
  out$Jmo <- make_Jmo_list(x, data, ...)
  if (has_ac_class(x, "car")) {
    gr <- get_ac_vars(x, "gr", class = "car")
    stopifnot(length(gr) <= 1L)
    if (isTRUE(nzchar(gr))) {
      out$locations <- levels(factor(get(gr, data)))
    } else {
      out$locations <- NA
    }
  }
  out
}

# extract data related to smooth terms
# for use in extract_old_standata
# @param version optional brms version number
make_sm_list <- function(x, data, version = NULL, ...) {
  stopifnot(is.btl(x))
  smterms <- all_terms(x[["sm"]])
  out <- named_list(smterms)
  if (length(smterms)) {
    knots <- get_knots(data)
    data <- rm_attr(data, "terms")
    # the spline penality has changed in 2.8.7 (#646)
    diagonal.penalty <- !isTRUE(version <= "2.8.6")
    gam_args <- list(
      data = data, knots = knots, 
      absorb.cons = TRUE, modCon = 3,
      diagonal.penalty = diagonal.penalty
    )
    for (i in seq_along(smterms)) {
      sc_args <- c(list(eval2(smterms[i])), gam_args)
      out[[i]] <- do_call(smoothCon, sc_args)
    }
  }
  out
}

# extract data related to gaussian processes
# for use in extract_old_standata
make_gp_list <- function(x, data, ...) {
  stopifnot(is.btl(x))
  out <- data_gp(x, data, raw = TRUE)
  out <- out[grepl("^(dmax)|(cmeans)", names(out))]
  out
}

# extract data related to monotonic effects
# for use in extract_old_standata
make_Jmo_list <- function(x, data, ...) {
  stopifnot(is.btl(x))
  out <- NULL
  if (length(attr(x$sp, "uni_mo"))) {
    # do it like data_sp()
    spef <- tidy_spef(x, data)
    Xmo <- lapply(unlist(spef$calls_mo), get_mo_values, data = data)
    out <- as.array(ulapply(Xmo, max))
  }
  out
}
