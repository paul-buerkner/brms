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
    prior, bterms = bframe,
    sample_prior = sample_prior
  )
  stanvars <- validate_stanvars(stanvars)
  threads <- validate_threads(threads)
  .standata(
    bterms, data = data, prior = prior,
    data2 = data2, stanvars = stanvars,
    threads = threads, ...
  )
}

# internal work function of 'standata'
# @param check_response check validity of the response?
# @param only_response extract data related to the response only?
# @param internal prepare Stan data for use in post-processing methods?
# @param basis original Stan data as prepared by 'standata_basis'
# @param ... currently ignored
# @return names list of data passed to Stan
.standata <- function(bterms, data, prior, stanvars, data2,
                           threads = threading(), check_response = TRUE,
                           only_response = FALSE, internal = FALSE,
                           basis = NULL, ...) {

  check_response <- as_one_logical(check_response)
  only_response <- as_one_logical(only_response)
  internal <- as_one_logical(internal)
  # order data for use in autocorrelation models
  data <- order_data(data, bterms = bterms)
  out <- data_response(
    bterms, data, check_response = check_response,
    internal = internal, basis = basis
  )
  if (!only_response) {
    ranef <- tidy_ranef(bterms, data, old_levels = basis$levels)
    meef <- tidy_meef(bterms, data, old_levels = basis$levels)
    index <- tidy_index(bterms, data)
    # pass as sdata so that data_special_prior knows about data_gr_global
    sdata_gr_global <- data_gr_global(ranef, data2 = data2)
    c(out) <- data_predictor(
      bterms, data = data, prior = prior, data2 = data2, ranef = ranef,
      sdata = sdata_gr_global, index = index, basis = basis
    )
    c(out) <- sdata_gr_global
    c(out) <- data_Xme(meef, data = data)
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
    # ensures current grouping levels are known in post-processing
    ranef_new <- tidy_ranef(bterms, data)
    meef_new <- tidy_meef(bterms, data)
    attr(out, "levels") <- get_levels(ranef_new, meef_new)
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
    basis <- standata_basis(bterms, data = object$data)
  }
  .standata(
    bterms, data = data, prior = object$prior,
    data2 = data2, stanvars = stanvars,
    threads = object$threads, basis = basis, ...
  )
}

# prepare basis data required for correct predictions from new data
# TODO: eventually export this function if we want to ensure full compatibility
#   with the 'empty' feature. see ?rename_pars for an example
standata_basis <- function(x, data, ...) {
  UseMethod("standata_basis")
}

#' @export
standata_basis.default <- function(x, data, ...) {
  list()
}

#' @export
standata_basis.mvbrmsterms <- function(x, data, ...) {
  out <- list()
  for (r in names(x$terms)) {
    out$resps[[r]] <- standata_basis(x$terms[[r]], data, ...)
  }
  out$levels <- get_levels(tidy_meef(x, data), tidy_ranef(x, data))
  out
}

#' @export
standata_basis.brmsterms <- function(x, data, ...) {
  out <- list()
  data <- subset_data(data, x)
  for (dp in names(x$dpars)) {
    out$dpars[[dp]] <- standata_basis(x$dpars[[dp]], data, ...)
  }
  for (nlp in names(x$nlpars)) {
    out$nlpars[[nlp]] <- standata_basis(x$nlpars[[nlp]], data, ...)
  }
  # old levels are required to select the right indices for new levels
  out$levels <- get_levels(tidy_meef(x, data), tidy_ranef(x, data))
  if (is_binary(x$family) || is_categorical(x$family)) {
    y <- model.response(model.frame(x$respform, data, na.action = na.pass))
    out$resp_levels <- levels(as.factor(y))
  }
  out
}

#' @export
standata_basis.btnl <- function(x, data, ...) {
  list()
}

#' @export
standata_basis.btl <- function(x, data, ...) {
  out <- list()
  out$sm <- standata_basis_sm(x, data, ...)
  out$gp <- standata_basis_gp(x, data, ...)
  out$sp <- standata_basis_sp(x, data, ...)
  out$ac <- standata_basis_ac(x, data, ...)
  out$bhaz <- standata_basis_bhaz(x, data, ...)
  out
}

# prepare basis data related to smooth terms
standata_basis_sm <- function(x, data, ...) {
  stopifnot(is.btl(x))
  smterms <- all_terms(x[["sm"]])
  out <- named_list(smterms)
  if (length(smterms)) {
    knots <- get_knots(data)
    data <- rm_attr(data, "terms")
    # the spline penalty has changed in 2.8.7 (#646)
    diagonal.penalty <- !require_old_default("2.8.7")
    gam_args <- list(
      data = data, knots = knots,
      absorb.cons = TRUE, modCon = 3,
      diagonal.penalty = diagonal.penalty
    )
    for (i in seq_along(smterms)) {
      sc_args <- c(list(eval2(smterms[i])), gam_args)
      sm <- do_call(smoothCon, sc_args)
      re <- vector("list", length(sm))
      for (j in seq_along(sm)) {
        re[[j]] <- mgcv::smooth2random(sm[[j]], names(data), type = 2)
      }
      out[[i]]$sm <- sm
      out[[i]]$re <- re
    }
  }
  out
}

# prepare basis data related to gaussian processes
standata_basis_gp <- function(x, data, ...) {
  stopifnot(is.btl(x))
  out <- data_gp(x, data, internal = TRUE)
  out <- out[grepl("^((Xgp)|(dmax)|(cmeans))", names(out))]
  out
}

# prepare basis data related to special terms
standata_basis_sp <- function(x, data, ...) {
  stopifnot(is.btl(x))
  out <- list()
  if (length(attr(x$sp, "uni_mo"))) {
    # do it like data_sp()
    spef <- tidy_spef(x, data)
    Xmo <- lapply(unlist(spef$calls_mo), get_mo_values, data = data)
    out$Jmo <- as.array(ulapply(Xmo, attr, "max"))
  }
  out
}

# prepare basis data related to autocorrelation structures
standata_basis_ac <- function(x, data, ...) {
  out <- list()
  if (has_ac_class(x, "car")) {
    gr <- get_ac_vars(x, "gr", class = "car")
    if (isTRUE(nzchar(gr))) {
      out$locations <- extract_levels(get(gr, data))
    } else {
      out$locations <- NA
    }
  }
  if (has_ac_class(x, "unstr")) {
    time <- get_ac_vars(x, "time", dim = "time")
    out$times <- extract_levels(get(time, data))
  }
  out
}

# prepare basis data for baseline hazards of the cox model
standata_basis_bhaz <- function(x, data, ...) {
  out <- list()
  if (is_cox(x$family)) {
    # compute basis matrix of the baseline hazard for the Cox model
    y <- model.response(model.frame(x$respform, data, na.action = na.pass))
    out$basis_matrix <- bhaz_basis_matrix(y, args = x$family$bhaz)
  }
  out
}
