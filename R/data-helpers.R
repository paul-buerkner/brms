# update data for use in brms functions
# @param data the data passed by the user
# @param bterms object of class brmsterms
# @param na.action function defining how to treat NAs
# @param drop.unused.levels should unused factor levels be removed?
# @param attr_terms a list of attributes of the terms object of 
#   the original model.frame; only used with newdata;
#   this ensures that (1) calls to 'poly' work correctly
#   and (2) that the number of variables matches the number 
#   of variable names; fixes issue #73
# @param knots: a list of knot values for GAMMs
# @return model.frame for use in brms functions
validate_data <- function(data, bterms, na.action = na.omit2,
                          drop.unused.levels = TRUE, knots = NULL,
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
  terms_all <- terms(bterms$allvars)
  attributes(terms_all)[names(attr_terms)] <- attr_terms
  data <- data_rsv_intercept(data, bterms = bterms)
  missing_vars <- setdiff(all.vars(terms_all), names(data))
  if (length(missing_vars)) {
    stop2("The following variables are missing in 'data':\n",
          collapse_comma(missing_vars))
  }
  for (v in intersect(vars_keep_na(bterms), names(data))) {
    attr(data[[v]], "keep_na") <- TRUE
  }
  # 'terms' prevents correct validation in 'model.frame'
  attr(data, "terms") <- NULL
  data <- model.frame(
    terms_all, data, na.action = na.action,
    drop.unused.levels = drop.unused.levels
  )
  if (any(grepl("__|_$", colnames(data)))) {
    stop2("Variable names may not contain double underscores ",
          "or underscores at the end.")
  }
  if (!isTRUE(nrow(data) > 0L)) {
    stop2("All observations in the data were removed presumably because of ",
          "NA values. If you are using the 'subset' addition term, please ",
          "make sure that variables are not NA even if locally unused.") 
  }
  groups <- get_group_vars(bterms)
  data <- combine_groups(data, groups)
  data <- fix_factor_contrasts(data, ignore = groups)
  attr(data, "knots") <- knots
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
  cov_names <- ulapply(get_re(bterms)$gcall, "[[", "cov")
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
      new.var <- get(sgroup[1], data)
      for (j in 2:length(sgroup)) {
        new.var <- paste0(new.var, "_", get(sgroup[j], data))
      }
      data[[group[[i]]]] <- new.var
    }
  } 
  data
}

# hard code factor contrasts to be independent of the global "contrasts" option
# @param data data.frame to be updated
# @param optdata: optional data.frame from which contrasts are taken if present
# @param ignore: names of variables for which not to fix contrasts
# @return 'data' with amended contrasts attributes
fix_factor_contrasts <- function(data, optdata = NULL, ignore = NULL) {
  stopifnot(is(data, "data.frame"))
  stopifnot(is.null(optdata) || is.list(optdata))
  optdata <- as.data.frame(optdata)  # fixes issue #105
  for (i in seq_along(data)) {
    needs_contrast <- is.factor(data[[i]]) && !names(data)[i] %in% ignore
    if (needs_contrast && is.null(attr(data[[i]], "contrasts"))) {
      old_contrasts <- attr(optdata[[names(data)[i]]], "contrasts")
      if (!is.null(old_contrasts)) {
        # take contrasts from optdata
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
  if (is.formula(bterms$adforms$subset)) {
    # only evaluate a subset of the data
    subset <- eval_rhs(bterms$adforms$subset)
    subset <- as.logical(eval2(subset$vars$subset, data))
    if (length(subset) != nrow(data)) {
      stop2("Length of 'subset' does not match the rows of 'data'.")
    }
    if (anyNA(subset)) {
      stop2("Subset variables may not contain NAs.")
    }
    data <- data[subset, , drop = FALSE]
  }
  if (!NROW(data)) {
    stop2(
      "All rows of 'data' were removed via 'subset'. ",
      "Please make sure that variables do not contain NAs ",
      "even in rows unused by the subsetted model. ",
      "Please also make sure that each subset variable is ",
      "TRUE for at least one observation."
    )
  }
  data
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
#' @param all_group_vars Optional names of grouping variables to be validated.
#'   Defaults to all grouping variables in the model.
#' @param ... Currently ignored.
#' 
#' @return A validated \code{'data.frame'} based on \code{newdata}.
#' 
#' @export
validate_newdata <- function(
  newdata, object, re_formula = NULL, allow_new_levels = FALSE,
  resp = NULL, check_response = TRUE, incl_autocor = TRUE,
  all_group_vars = NULL, ...
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
  if (is.mvbrmsterms(bterms) && !is.null(resp)) {
    # variables not used in the included model parts
    # do not need to be specified in newdata
    resp <- validate_resp(resp, bterms$responses)
    reqvars <- allvars_formula(lapply(bterms$terms[resp], "[[", "allvars"))
    not_reqvars <- setdiff(all.vars(bterms$allvars), all.vars(reqvars))
    not_reqvars <- setdiff(not_reqvars, names(newdata))
    if (length(not_reqvars)) {
      newdata[, not_reqvars] <- NA
    }
  }
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
      newdata[, missing_resp] <- NA
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
    newdata[, mis_group_vars] <- NA
  }
  newdata <- combine_groups(newdata, new_group_vars)
  # validate factor levels in newdata
  if (is.null(all_group_vars)) {
    all_group_vars <- get_group_vars(object) 
  }
  dont_check <- c(all_group_vars, cens_vars)
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
        new_levels <- levels(new_factor)
        old_levels <- levels(factors[[i]])
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
  num_names <- setdiff(num_names, all_group_vars)
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
  if (length(unused_vars)) {
    newdata[, unused_vars] <- NA
  }
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
    newdata, bterms = bterms, na.action = na.pass, 
    drop.unused.levels = FALSE, attr_terms = attr_terms,
    knots = get_knots(object$data)
  )
  newdata
}

# extract the current data set
current_data <- function(object, newdata = NULL, ...) {
  stopifnot(is.brmsfit(object))
  if (is.null(newdata)) {
    data <- object$data
  } else {
    data <- validate_newdata(newdata, object = object, ...)
  }
  data
}
