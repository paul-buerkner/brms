# update data for use in brms functions
# @param data the data passed by the user
# @param bterms object of class brmsterms
# @param na_action function defining how to treat NAs
# @param drop_unused_levels should unused factor levels be removed?
# @param attr_terms a list of attributes of the terms object of
#   the original model.frame; only used with newdata;
#   this ensures that (1) calls to 'poly' work correctly
#   and (2) that the number of variables matches the number
#   of variable names; fixes issue #73
# @param knots: a list of knot values for GAMMs
# @return model.frame for use in brms functions
validate_data <- function(data, bterms, data2 = list(), knots = NULL,
                          na_action = na_omit, drop_unused_levels = TRUE,
                          attr_terms = NULL) {
  if (missing(data)) {
    stop2("Data must be specified using the 'data' argument.")
  }
  if (is.null(knots)) {
    knots <- get_knots(data)
  }
  data <- try(as.data.frame(data), silent = TRUE)
  if (is(data, "try-error")) {
    stop2("Argument 'data' must be coercible to a data.frame.")
  }
  if (!isTRUE(nrow(data) > 0L)) {
    stop2("Argument 'data' does not contain observations.")
  }
  data <- data_rsv_intercept(data, bterms = bterms)
  all_vars_formula <- bterms$allvars
  missing_vars <- setdiff(all_vars(all_vars_formula), names(data))
  if (length(missing_vars)) {
    missing_vars2 <- setdiff(missing_vars, names(data2))
    if (length(missing_vars2)) {
      stop2("The following variables can neither be found in ",
            "'data' nor in 'data2':\n", collapse_comma(missing_vars2))
    }
    # all initially missing variables can be found in 'data2'
    # they are not necessarily of the length required for 'data'
    # so need to be excluded from the evaluation of 'model.frame'
    missing_vars_formula <- paste0(". ~ . ", collapse(" - ", missing_vars))
    all_vars_formula <- update(all_vars_formula, missing_vars_formula)
  }
  all_vars_terms <- terms(all_vars_formula)
  # ensure that 'data2' comes first in the search path
  # during the evaluation of model.frame
  terms_env <- environment(all_vars_terms)
  environment(all_vars_terms) <- as.environment(as.list(data2))
  parent.env(environment(all_vars_terms)) <- terms_env
  attributes(all_vars_terms)[names(attr_terms)] <- attr_terms
  # 'terms' prevents correct validation in 'model.frame'
  attr(data, "terms") <- NULL
  data <- model.frame(
    all_vars_terms, data, na.action = na.pass,
    drop.unused.levels = drop_unused_levels
  )
  data <- na_action(data, bterms = bterms)
  if (any(grepl("__|_$", colnames(data)))) {
    stop2("Variable names may not contain double underscores ",
          "or underscores at the end.")
  }
  if (!isTRUE(nrow(data) > 0L)) {
    stop2("All observations in the data were removed ",
          "presumably because of NA values.")
  }
  groups <- get_group_vars(bterms)
  data <- combine_groups(data, groups)
  data <- fix_factor_contrasts(data, ignore = groups)
  attr(data, "knots") <- knots
  attr(data, "drop_unused_levels") <- drop_unused_levels
  data
}

# validate the 'data2' argument
# @param data2 a named list of data objects
# @param bterms object returned by 'brmsterms'
# @param ... more named list to pass objects to data2 from other sources
#   only required for backwards compatibility with deprecated arguments
# @return a validated named list of data objects
validate_data2 <- function(data2, bterms, ...) {
  # TODO: specify spline-related matrices in 'data2'
  # this requires adding another parser layer with bterms and data as input
  if (is.null(data2)) {
    data2 <- list()
  }
  if (!is.list(data2)) {
    stop2("'data2' must be a list.")
  }
  if (length(data2) && !is_named(data2)) {
    stop2("All elements of 'data2' must be named.")
  }
  dots <- list(...)
  for (i in seq_along(dots)) {
    if (length(dots[[i]])) {
      stopifnot(is.list(dots[[i]]), is_named(dots[[i]]))
      data2[names(dots[[i]])] <- dots[[i]]
    }
  }
  # validate autocorrelation matrices
  acef <- tidy_acef(bterms)
  sar_M_names <- get_ac_vars(acef, "M", class = "sar")
  for (M in sar_M_names) {
    data2[[M]] <- validate_sar_matrix(get_from_data2(M, data2))
    attr(data2[[M]], "obs_based_matrix") <- TRUE
  }
  car_M_names <- get_ac_vars(acef, "M", class = "car")
  for (M in car_M_names) {
    data2[[M]] <- validate_car_matrix(get_from_data2(M, data2))
    # observation based CAR matrices are deprecated and
    # there is no need to label them as observation based
  }
  fcor_M_names <- get_ac_vars(acef, "M", class = "fcor")
  for (M in fcor_M_names) {
    data2[[M]] <- validate_fcor_matrix(get_from_data2(M, data2))
    attr(data2[[M]], "obs_based_matrix") <- TRUE
  }
  # validate within-group covariance matrices
  cov_names <- ufrom_list(get_re(bterms)$gcall, "cov")
  cov_names <- cov_names[nzchar(cov_names)]
  for (cov in cov_names) {
    data2[[cov]] <- validate_recov_matrix(get_from_data2(cov, data2))
  }
  data2
}

# get an object from the 'data2' argument
get_from_data2 <- function(x, data2) {
  if (!x %in% names(data2)) {
    stop2("Object '", x, "' was not found in 'data2'.")
  }
  get(x, data2)
}

# index observation based elements in 'data2'
# @param data2 a named list of objects
# @param i observation based indices
# @return data2 with potentially indexed elements
subset_data2 <- function(data2, i) {
  if (!length(data2)) {
    return(data2)
  }
  stopifnot(is.list(data2), is_named(data2))
  for (var in names(data2)) {
    if (isTRUE(attr(data2[[var]], "obs_based_matrix"))) {
      # matrices with dimensions equal to the number of observations
      data2[[var]] <- data2[[var]][i, i, drop = FALSE]
      attr(data2[[var]], "obs_based_matrix") <- TRUE
    }
  }
  data2
}

# add the reserved intercept variables to the data
data_rsv_intercept <- function(data, bterms) {
  fe_forms <- get_effect(bterms, "fe")
  if (any(ulapply(fe_forms, no_int))) {
    if ("intercept" %in% ulapply(fe_forms, all_vars)) {
      warning2("Reserved variable name 'intercept' is deprecated. ",
               "Please use 'Intercept' instead.")
    }
    if (any(data[["intercept"]] != 1)) {
      stop2("Variable name 'intercept' is reserved in models ",
            "without a population-level intercept.")
    }
    if (any(data[["Intercept"]] != 1)) {
      stop2("Variable name 'Intercept' is reserved in models ",
            "without a population-level intercept.")
    }
    data$intercept <- data$Intercept <- rep(1, length(data[[1]]))
  }
  data
}

# combine grouping factors to form new variables
# @param data data.frame to be updated
# @param ... the grouping factors to be combined
# @return 'data' including the new combined grouping factors
combine_groups <- function(data, ...) {
  group <- c(...)
  for (i in seq_along(group)) {
    sgroup <- unlist(strsplit(group[[i]], ":"))
    if (length(sgroup) > 1L && !group[[i]] %in% names(data)) {
      new_var <- get(sgroup[1], data)
      for (j in 2:length(sgroup)) {
        new_var <- paste0(new_var, "_", get(sgroup[j], data))
      }
      data[[group[[i]]]] <- new_var
    }
  }
  data
}

# hard code factor contrasts to be independent of the global "contrasts" option
# @param data data.frame to be updated
# @param olddata: optional data.frame from which contrasts are taken if present
# @param ignore: names of variables for which not to fix contrasts
# @return 'data' with amended contrasts attributes
fix_factor_contrasts <- function(data, olddata = NULL, ignore = NULL) {
  stopifnot(is(data, "data.frame"))
  stopifnot(is.null(olddata) || is.list(olddata))
  olddata <- as.data.frame(olddata)  # fixes issue #105
  for (i in seq_along(data)) {
    needs_contrast <- is.factor(data[[i]]) && !names(data)[i] %in% ignore
    if (needs_contrast && is.null(attr(data[[i]], "contrasts"))) {
      old_contrasts <- attr(olddata[[names(data)[i]]], "contrasts")
      if (!is.null(old_contrasts)) {
        # take contrasts from olddata
        contrasts(data[[i]]) <- old_contrasts
      } else if (length(unique(data[[i]])) > 1L) {
        # avoid error when supplying only a single level
        # hard code current global "contrasts" option
        contrasts(data[[i]]) <- contrasts(data[[i]])
      }
    }
  }
  data
}

# order data for use in time-series models
# @param data data.frame to be ordered
# @param bterms brmsterms of mvbrmsterms object
# @return 'data' potentially ordered differently
order_data <- function(data, bterms) {
  # ordering does only matter for time-series models
  time <- get_ac_vars(bterms, "time", dim = "time")
  gr <- get_ac_vars(bterms, "gr", dim = "time")
  if (length(time) > 1L || length(gr) > 1L) {
    stop2("All time-series structures must have the same ",
          "'time' and 'gr' variables.")
  }
  if (length(time) || length(gr)) {
    if (length(gr)) {
      gv <- data[[gr]]
    } else {
      gv <- rep(1L, nrow(data))
    }
    if (length(time)) {
      tv <- data[[time]]
    } else {
      tv <- seq_rows(data)
    }
    if (any(duplicated(data.frame(gv, tv)))) {
      stop2("Time points within groups must be unique.")
    }
    new_order <- do_call(order, list(gv, tv))
    data <- data[new_order, , drop = FALSE]
    # old_order will allow to retrieve the initial order of the data
    attr(data, "old_order") <- order(new_order)
  }
  data
}

# subset data according to addition argument 'subset'
subset_data <- function(data, bterms) {
  if (has_subset(bterms)) {
    # only evaluate a subset of the data
    subset <- as.logical(get_ad_values(bterms, "subset", "subset", data))
    if (length(subset) != nrow(data)) {
      stop2("Length of 'subset' does not match the rows of 'data'.")
    }
    if (anyNA(subset)) {
      stop2("Subset variables may not contain NAs.")
    }
    # cross-formula indexing is no longer trivial for subsetted models
    check_cross_formula_indexing(bterms)
    data <- data[subset, , drop = FALSE]
    attr(data, "subset") <- subset
  }
  if (!NROW(data)) {
    stop2(
      "All rows of 'data' were removed via 'subset'. ",
      "Please make sure that variables do not contain NAs ",
      "for observations in which they are supposed to be used. ",
      "Please also make sure that each subset variable is ",
      "TRUE for at least one observation."
    )
  }
  data
}

# like stats:::na.omit.data.frame but allows to certain NA values
na_omit <- function(object, bterms, ...) {
  stopifnot(is.data.frame(object))
  nobs <- nrow(object)
  if (is.mvbrmsterms(bterms)) {
    responses <- names(bterms$terms)
    subsets <- lapply(bterms$terms, get_ad_values, "subset", "subset", object)
    vars_sub <- lapply(bterms$terms, function(x) all_vars(x$allvars))
  }
  vars_keep_na <- vars_keep_na(bterms)
  omit <- logical(nobs)
  for (v in names(object)) {
    x <- object[[v]]
    vars_v <- all_vars(v)
    keep_all_na <- all(vars_v %in% vars_keep_na)
    if (!is.atomic(x) || keep_all_na) {
      next
    }
    if (!is.mvbrmsterms(bterms)) {
      # remove all NAs in this variable
      keep_na <- rep(FALSE, nobs)
    } else {
      # allow to retain NAs in subsetted variables
      keep_na <- rep(TRUE, nobs)
      for (r in responses) {
        if (any(vars_v %in% vars_sub[[r]])) {
          if (!is.null(subsets[[r]])) {
            # keep NAs ignored because of 'subset'
            keep_na <- keep_na & !subsets[[r]]
          } else {
            # remove all NAs in this variable
            keep_na <- keep_na & FALSE
          }
        }
      }
    }
    is_na <- is.na(x)
    d <- dim(is_na)
    if (is.null(d) || length(d) != 2L) {
      omit <- omit | (is_na & !keep_na)
    } else {
      for (ii in seq_len(d[2L])) {
        omit <- omit | (is_na[, ii] & !keep_na)
      }
    }
  }
  if (any(omit > 0L)) {
    out <- object[!omit, , drop = FALSE]
    temp <- setNames(seq(omit)[omit], attr(object, "row.names")[omit])
    attr(temp, "class") <- "omit"
    attr(out, "na.action") <- temp
    warning2("Rows containing NAs were excluded from the model.")
  } else {
    out <- object
  }
  out
}

# get a single value per group
# @param x vector of values to extract one value per group
# @param gr vector of grouping values
# @return a vector of the same length as unique(group)
get_one_value_per_group <- function(x, gr) {
  stopifnot(length(x) == length(gr))
  not_dupl_gr <- !duplicated(gr)
  gr_unique <- gr[not_dupl_gr]
  to_order <- order(gr_unique)
  gr_unique <- gr_unique[to_order]
  out <- x[not_dupl_gr][to_order]
  names(out) <- gr_unique
  out
}

# extract knots values for use in spline terms
get_knots <- function(data) {
  attr(data, "knots", TRUE)
}

get_drop_unused_levels <- function(data) {
  out <- attr(data, "drop_unused_levels", TRUE) %||% TRUE
}

# extract name of the data as originally passed by the user
get_data_name <- function(data) {
  out <- attr(data, "data_name", TRUE)
  if (is.null(out)) {
    out <- "NULL"
  }
  out
}

#' Validate New Data
#'
#' Validate new data passed to post-processing methods of \pkg{brms}. Unless you
#' are a package developer, you will rarely need to call \code{validate_newdata}
#' directly.
#'
#' @inheritParams prepare_predictions
#' @param newdata A \code{data.frame} containing new data to be validated.
#' @param object A \code{brmsfit} object.
#' @param check_response Logical; Indicates if response variables should
#'   be checked as well. Defaults to \code{TRUE}.
#' @param group_vars Optional names of grouping variables to be validated.
#'   Defaults to all grouping variables in the model.
#' @param req_vars Optional names of variables required in \code{newdata}.
#'   If \code{NULL} (the default), all variables in the original data
#'   are required (unless ignored for some other reason).
#' @param ... Currently ignored.
#'
#' @return A validated \code{'data.frame'} based on \code{newdata}.
#'
#' @export
validate_newdata <- function(
  newdata, object, re_formula = NULL, allow_new_levels = FALSE,
  newdata2 = NULL, resp = NULL, check_response = TRUE,
  incl_autocor = TRUE, group_vars = NULL, req_vars = NULL, ...
) {
  newdata <- try(as.data.frame(newdata), silent = TRUE)
  if (is(newdata, "try-error")) {
    stop2("Argument 'newdata' must be coercible to a data.frame.")
  }
  object <- restructure(object)
  object <- exclude_terms(object, incl_autocor = incl_autocor)
  resp <- validate_resp(resp, object)
  new_formula <- update_re_terms(formula(object), re_formula)
  bterms <- brmsterms(new_formula, resp_rhs_all = FALSE)

  # fill values of not required variables
  all_vars <- all.vars(bterms$allvars)
  if (is.null(req_vars)) {
    req_vars <- all_vars
  } else {
    req_vars <- as.character(req_vars)
    req_vars <- intersect(req_vars, all_vars)
  }
  if (is.mvbrmsterms(bterms) && !is.null(resp)) {
    # variables not used in the included model parts
    # do not need to be specified in newdata
    resp <- validate_resp(resp, bterms$responses)
    form_req_vars <- from_list(bterms$terms[resp], "allvars")
    form_req_vars <- allvars_formula(form_req_vars)
    req_vars <- intersect(req_vars, all.vars(form_req_vars))
  }
  not_req_vars <- setdiff(all_vars, req_vars)
  not_req_vars <- setdiff(not_req_vars, names(newdata))
  newdata <- fill_newdata(newdata, not_req_vars, object$data)
  # check response and addition variables
  only_resp <- all.vars(bterms$respform)
  only_resp <- setdiff(only_resp, all.vars(rhs(bterms$allvars)))
  # always require 'dec' variables to be specified
  dec_vars <- get_ad_vars(bterms, "dec")
  missing_resp <- setdiff(c(only_resp, dec_vars), names(newdata))
  if (length(missing_resp)) {
    if (check_response) {
      stop2("Response variables must be specified in 'newdata'.\n",
            "Missing variables: ", collapse_comma(missing_resp))
    } else {
      newdata <- fill_newdata(newdata, missing_resp)
    }
  }
  # censoring and weighting vars are unused in post-processing methods
  cens_vars <- get_ad_vars(bterms, "cens")
  for (v in setdiff(cens_vars, names(newdata))) {
    newdata[[v]] <- 0
  }
  weights_vars <- get_ad_vars(bterms, "weights")
  for (v in setdiff(weights_vars, names(newdata))) {
    newdata[[v]] <- 1
  }
  mf <- model.frame(object)
  for (i in seq_along(mf)) {
    if (is_like_factor(mf[[i]])) {
      mf[[i]] <- as.factor(mf[[i]])
    }
  }
  # fixes issue #279
  newdata <- data_rsv_intercept(newdata, bterms)
  new_group_vars <- get_group_vars(bterms)
  if (allow_new_levels && length(new_group_vars)) {
    # grouping factors do not need to be specified
    # by the user if new levels are allowed
    mis_group_vars <- new_group_vars[!grepl(":", new_group_vars)]
    mis_group_vars <- setdiff(mis_group_vars, names(newdata))
    newdata <- fill_newdata(newdata, mis_group_vars)
  }
  newdata <- combine_groups(newdata, new_group_vars)
  # validate factor levels in newdata
  if (is.null(group_vars)) {
    group_vars <- get_group_vars(object)
  }
  do_check <- union(get_pred_vars(bterms), get_int_vars(bterms))
  # do not check variables from the 'unused' argument #1238
  unused_arg_vars <- get_unused_arg_vars(bterms)
  dont_check <- unique(c(group_vars, cens_vars, unused_arg_vars))
  dont_check <- setdiff(dont_check, do_check)
  dont_check <- names(mf) %in% dont_check
  is_factor <- ulapply(mf, is.factor)
  factors <- mf[is_factor & !dont_check]
  if (length(factors)) {
    factor_names <- names(factors)
    for (i in seq_along(factors)) {
      new_factor <- newdata[[factor_names[i]]]
      if (!is.null(new_factor)) {
        if (!is.factor(new_factor)) {
          new_factor <- factor(new_factor)
        }
        old_levels <- levels(factors[[i]])
        if (length(old_levels) <= 1L) {
          # contrasts are not defined for factors with 1 or fewer levels
          next
        }
        new_levels <- levels(new_factor)
        old_contrasts <- contrasts(factors[[i]])
        old_ordered <- is.ordered(factors[[i]])
        to_zero <- is.na(new_factor) | new_factor %in% "zero__"
        # don't add the 'zero__' level to response variables
        is_resp <- factor_names[i] %in% all.vars(bterms$respform)
        if (!is_resp && any(to_zero)) {
          levels(new_factor) <- c(new_levels, "zero__")
          new_factor[to_zero] <- "zero__"
          old_levels <- c(old_levels, "zero__")
          old_contrasts <- rbind(old_contrasts, zero__ = 0)
        }
        if (any(!new_levels %in% old_levels)) {
          stop2(
            "New factor levels are not allowed.",
            "\nLevels allowed: ", collapse_comma(old_levels),
            "\nLevels found: ", collapse_comma(new_levels)
          )
        }
        newdata[[factor_names[i]]] <-
          factor(new_factor, old_levels, ordered = old_ordered)
        # don't use contrasts(.) here to avoid dimension checks
        attr(newdata[[factor_names[i]]], "contrasts") <- old_contrasts
      }
    }
  }
  # check if originally numeric variables are still numeric
  num_names <- names(mf)[!is_factor]
  num_names <- setdiff(num_names, group_vars)
  for (nm in intersect(num_names, names(newdata))) {
    if (!anyNA(newdata[[nm]]) && !is.numeric(newdata[[nm]])) {
      stop2("Variable '", nm, "' was originally ",
            "numeric but is not in 'newdata'.")
    }
  }
  # validate monotonic variables
  mo_vars <- get_sp_vars(bterms, "mo")
  if (length(mo_vars)) {
    # factors have already been checked
    num_mo_vars <- names(mf)[!is_factor & names(mf) %in% mo_vars]
    for (v in num_mo_vars) {
      new_values <- get(v, newdata)
      min_value <- min(mf[[v]])
      invalid <- new_values < min_value | new_values > max(mf[[v]])
      invalid <- invalid | !is_wholenumber(new_values)
      if (sum(invalid)) {
        stop2("Invalid values in variable '", v, "': ",
              collapse_comma(new_values[invalid]))
      }
      attr(newdata[[v]], "min") <- min_value
    }
  }
  # update_data expects all original variables to be present
  used_vars <- c(names(newdata), all.vars(bterms$allvars))
  used_vars <- union(used_vars, rsv_vars(bterms))
  all_vars <- all.vars(str2formula(names(mf)))
  unused_vars <- setdiff(all_vars, used_vars)
  newdata <- fill_newdata(newdata, unused_vars)
  # validate grouping factors
  new_ranef <- tidy_ranef(bterms, data = mf)
  new_meef <- tidy_meef(bterms, data = mf)
  old_levels <- get_levels(new_ranef, new_meef)
  if (!allow_new_levels) {
    new_levels <- get_levels(
      tidy_ranef(bterms, data = newdata),
      tidy_meef(bterms, data = newdata)
    )
    for (g in names(old_levels)) {
      unknown_levels <- setdiff(new_levels[[g]], old_levels[[g]])
      if (length(unknown_levels)) {
        unknown_levels <- collapse_comma(unknown_levels)
        stop2(
          "Levels ", unknown_levels, " of grouping factor '", g, "' ",
          "cannot be found in the fitted model. ",
          "Consider setting argument 'allow_new_levels' to TRUE."
        )
      }
    }
  }
  # ensure correct handling of functions like 'poly' or 'scale'
  old_terms <- attr(object$data, "terms")
  attr_terms <- c("variables", "predvars")
  attr_terms <- attributes(old_terms)[attr_terms]
  newdata <- validate_data(
    newdata, bterms = bterms, na_action = na.pass,
    drop_unused_levels = FALSE, attr_terms = attr_terms,
    data2 = current_data2(object, newdata2),
    knots = get_knots(object$data)
  )
  newdata
}

# fill newdata with values for not required variables
# @param newdata data.frame to be filled
# @param vars character vector of not required variables
# @param olddata optional data.frame to take values from
# @param n row number of olddata to extract values from
fill_newdata <- function(newdata, vars, olddata = NULL, n = 1L) {
  stopifnot(is.data.frame(newdata), is.character(vars))
  vars <- setdiff(vars, names(newdata))
  if (is.null(olddata)) {
    if (length(vars)) {
      newdata[, vars] <- NA
    }
    return(newdata)
  }
  stopifnot(is.data.frame(olddata), length(n) == 1L)
  for (v in vars) {
    # using NA for variables is not safe in all cases
    # for example when processing splines using mgcv
    # hence it is safer to use existing data values
    cval <- olddata[n, v] %||% NA
    if (length(dim(cval)) == 2L) {
      # matrix columns don't have automatic broadcasting apparently
      cval <- matrix(cval, nrow(newdata), ncol(cval), byrow = TRUE)
    }
    newdata[[v]] <- cval
  }
  newdata
}

# validate new data2
validate_newdata2 <- function(newdata2, object, ...) {
  stopifnot(is.brmsfit(object))
  bterms <- brmsterms(object$formula)
  validate_data2(newdata2, bterms = bterms, ...)
}

# extract the current data
current_data <- function(object, newdata = NULL, skip_validate = FALSE, ...) {
  stopifnot(is.brmsfit(object))
  if (is.null(newdata)) {
    data <- object$data
  } else if(skip_validate) {
    data <- newdata
  } else {
    data <- validate_newdata(newdata, object = object, ...)
  }
  data
}

# extract the current data2
current_data2 <- function(object, newdata2 = NULL, ...) {
  stopifnot(is.brmsfit(object))
  if (is.null(newdata2)) {
    data2 <- object$data2
  } else {
    data2 <- validate_newdata2(newdata2, object = object, ...)
  }
  data2
}
