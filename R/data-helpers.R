update_data <- function(data, bterms, na.action = na.omit2,
                        drop.unused.levels = TRUE,
                        terms_attr = NULL, knots = NULL) {
  # Update data for use in brms functions
  # Args:
  #   data: the original data.frame
  #   bterms: object of class brmsterms
  #   na.action: function defining how to treat NAs
  #   drop.unused.levels: indicates if unused factor levels
  #     should be removed
  #   terms_attr: a list of attributes of the terms object of 
  #     the original model.frame; only used with newdata;
  #     this ensures that (1) calls to 'poly' work correctly
  #     and (2) that the number of variables matches the number 
  #     of variable names; fixes issue #73
  #   knots: a list of knot values for GAMMs
  # Returns:
  #   model.frame for use in brms functions
  if (missing(data)) {
    stop2("Argument 'data' is missing.")
  }
  if (isTRUE(attr(data, "brmsframe"))) {
    return(data)
  }
  if (is.null(knots)) {
    knots <- attr(data, "knots", TRUE)
  }
  data <- try(as.data.frame(data), silent = TRUE)
  if (is(data, "try-error")) {
    stop2("Argument 'data' must be coercible to a data.frame.")
  }
  if (!isTRUE(nrow(data) > 0L)) {
    stop2("Argument 'data' does not contain observations.")
  }
  bterms$allvars <- terms(bterms$allvars)
  attributes(bterms$allvars)[names(terms_attr)] <- terms_attr
  data <- data_rsv_intercept(data, bterms = bterms)
  missing_vars <- setdiff(all.vars(bterms$allvars), names(data))
  if (length(missing_vars)) {
    message(
      "The following variables are missing in 'data' and ",
      "will be searched in the enclosing environments: ",
      collapse_comma(missing_vars)
    )
  }
  data <- model.frame(
    bterms$allvars, data, na.action = na.pass,
    drop.unused.levels = drop.unused.levels
  )
  nrow_with_na <- nrow(data)
  data <- na.action(data, ignore = vars_keep_na(bterms))
  if (nrow(data) != nrow_with_na) {
    warning2("Rows containing NAs were excluded from the model")
  }
  if (any(grepl("__|_$", colnames(data)))) {
    stop2("Variable names may not contain double underscores ",
          "or underscores at the end.")
  }
  groups <- get_groups(bterms)
  data <- combine_groups(data, groups)
  data <- fix_factor_contrasts(data, ignore = groups)
  attr(data, "knots") <- knots
  attr(data, "brmsframe") <- TRUE
  data
}

data_rsv_intercept <- function(data, bterms) {
  # add the resevered variable 'intercept' to the data
  # Args:
  #   data: data.frame or list
  #   bterms: object of class brmsterms
  fe_forms <- get_effect(bterms, "fe")
  rsv_int <- any(ulapply(fe_forms, attr, "rsv_intercept"))
  if (rsv_int) {
    if (any(data[["intercept"]] != 1)) {
      stop2("Variable name 'intercept' is resevered in models ",
            "without a population-level intercept.")
    }
    data$intercept <- rep(1, length(data[[1]]))
  }
  data
}

combine_groups <- function(data, ...) {
  # combine grouping factors
  # Args:
  #   data: a data.frame
  #   ...: the grouping factors to be combined. 
  # Returns:
  #   a data.frame containing all old variables and 
  #   the new combined grouping factors
  group <- c(...)
  for (i in seq_along(group)) {
    sgroup <- unlist(strsplit(group[[i]], ":"))
    if (length(sgroup) > 1L) {
      new.var <- get(sgroup[1], data)
      for (j in 2:length(sgroup)) {
        new.var <- paste0(new.var, "_", get(sgroup[j], data))
      }
      data[[group[[i]]]] <- new.var
    }
  } 
  data
}

fix_factor_contrasts <- function(data, optdata = NULL, ignore = NULL) {
  # hard code factor contrasts to be independent
  # of the global "contrasts" option
  # Args:
  #   data: a data.frame
  #   optdata: optional data.frame from which contrasts
  #     are taken if present
  #   ignore: names of variables for which not to fix contrasts
  # Returns:
  #   a data.frame with amended contrasts attributes
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

order_data <- function(data, bterms) {
  # order data for use in time-series models
  # Args:
  #   data: data.frame to be ordered
  #   bterms: brmsterms of mvbrmsterms object
  # Returns:
  #   potentially ordered data
  time <- get_autocor_vars(bterms, "time")
  # ordering does not matter for the CAR structure
  group <- get_autocor_vars(bterms, "group", incl_car = FALSE)
  if (length(time) > 1L || length(group) > 1L) {
    stop2("All autocorrelation structures must have the same ",
          "time and group variables.")
  }
  time_var <- if (length(time)) data[[time]]
  group_var <- if (length(group)) data[[group]]
  to_order <- rmNULL(list(group_var, time_var))
  if (length(to_order)) {
    new_order <- do.call(order, to_order)
    data <- data[new_order, , drop = FALSE]
    # old_order will allow to retrieve the initial order of the data
    attr(data, "old_order") <- order(new_order)
  }
  data
}

validate_newdata <- function(
  newdata, fit, re_formula = NULL, allow_new_levels = FALSE,
  resp = NULL, check_response = FALSE, only_response = FALSE,
  incl_autocor = TRUE, return_standata = TRUE,
  all_group_vars = NULL, new_objects = list()
) {
  # amend newdata passed to predict and fitted methods
  # Args:
  #   newdata: a data.frame containing new data for prediction 
  #   fit: an object of class brmsfit
  #   re_formula: a group-level formula
  #   allow_new_levels: Are new group-levels allowed?
  #   resp: optional name of response variables whose 
  #     variables should be checked
  #   check_response: Should response variables be checked
  #     for existence and validity?
  #   only_response: compute only response related stuff
  #     in make_standata?
  #   incl_autocor: Check data of autocorrelation terms?
  #   return_standata: Compute the data to be passed to Stan
  #     or just return the updated newdata?
  #   all_group_vars: optional names of all grouping 
  #     variables in the model
  #   new_objects: see function 'add_new_objects'
  # Returns:
  #   updated data.frame being compatible with formula(fit)
  if (!incl_autocor) {
    fit <- remove_autocor(fit) 
  }
  if (is.null(newdata)) {
    # to shorten expressions in S3 methods such as predict.brmsfit
    if (return_standata) {
      control <- list(not4stan = TRUE, save_order = TRUE)
      newdata <- standata(
        fit, re_formula = re_formula, check_response = check_response, 
        only_response = only_response, control = control
      )
    } else {
      newdata <- model.frame(fit)
    }
    return(newdata)
  }
  newdata <- try(as.data.frame(newdata), silent = TRUE)
  if (is(newdata, "try-error")) {
    stop2("Argument 'newdata' must be coercible to a data.frame.")
  }
  new_formula <- update_re_terms(formula(fit), re_formula = re_formula)
  bterms <- parse_bf(new_formula, resp_rhs_all = FALSE)
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
  dec_vars <- get_advars(bterms, "dec")
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
  cens_vars <- get_advars(bterms, "cens")
  for (v in setdiff(cens_vars, names(newdata))) {
    newdata[[v]] <- 0
  }
  weights_vars <- get_advars(bterms, "weights")
  for (v in setdiff(weights_vars, names(newdata))) {
    newdata[[v]] <- 1
  }
  mf <- model.frame(fit)
  for (i in seq_along(mf)) {
    if (is_like_factor(mf[[i]])) {
      mf[[i]] <- as.factor(mf[[i]])
    }
  }
  # fixes issue #279
  newdata <- data_rsv_intercept(newdata, bterms)
  new_ranef <- tidy_ranef(bterms, data = mf)
  group_vars <- get_group_vars(new_ranef)
  group_vars <- union(group_vars, get_autocor_vars(bterms, "group"))
  if (allow_new_levels && length(group_vars)) {
    # grouping factors do not need to be specified 
    # by the user if new levels are allowed
    new_gf <- unique(unlist(strsplit(group_vars, split = ":")))
    missing_gf <- setdiff(new_gf, names(newdata))
    newdata[, missing_gf] <- NA
  }
  newdata <- combine_groups(newdata, group_vars)
  # validate factor levels in newdata
  if (is.null(all_group_vars)) {
    all_group_vars <- get_group_vars(fit) 
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
        newdata[[factor_names[i]]] <- factor(new_factor, old_levels)
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
  old_levels <- attr(new_ranef, "levels")
  if (!allow_new_levels) {
    new_levels <- attr(tidy_ranef(bterms, data = newdata), "levels")
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
  if (return_standata) {
    fit <- add_new_objects(fit, newdata, new_objects)
    control <- list(
      new = TRUE, not4stan = TRUE, 
      old_levels = old_levels, save_order = TRUE
    )
    # ensure correct handling of functions like poly or scale
    old_terms <- attr(mf, "terms")
    terms_attr <- c("variables", "predvars")
    control$terms_attr <- attributes(old_terms)[terms_attr]
    control$old_standata <- extract_old_standata(bterms, data = mf)
    knots <- attr(mf, "knots")
    newdata <- make_standata(
      new_formula, data = newdata, 
      knots = knots, check_response = check_response,
      only_response = only_response, control = control
    )
  }
  newdata
}

add_new_objects <- function(x, newdata, new_objects = list()) {
  # allows for updating of objects containing new data
  # which cannot be passed via argument 'newdata'
  # Args:
  #   x: object of class 'brmsfit'
  #   new_objects: optional list of new objects
  # Return:
  #   a possibly updated 'brmsfit' object
  stopifnot(is.brmsfit(x), is.data.frame(newdata))
  if (is.cor_sar(x$autocor)) {
    if ("W" %in% names(new_objects)) {
      x$autocor <- cor_sar(new_objects$W, type = x$autocor$type)
    } else {
      message("Using the identity matrix as weighting matrix by default")
      x$autocor$W <- diag(nrow(newdata))
    }
  }
  # do not include cor_car here as the adjacency matrix
  # (or subsets of it) should be the same for newdata 
  if (is.cor_fixed(x$autocor)) {
    if ("V" %in% names(new_objects)) {
      x$autocor <- cor_fixed(new_objects$V)
    } else {
      message("Using the median variance by default")
      median_V <- median(diag(x$autocor$V), na.rm = TRUE)
      x$autocor$V <- diag(median_V, nrow(newdata)) 
    }
  }
  x
}

get_model_matrix <- function(formula, data = environment(formula),
                             cols2remove = NULL, rename = TRUE, ...) {
  # Construct Design Matrices for \code{brms} models
  # Args:
  #   formula: An object of class formula
  #   data: A data frame created with model.frame. 
  #         If another sort of object, model.frame is called first.
  #   cols2remove: names of the columns to remove from 
  #                the model matrix (mainly used for intercepts)
  #   rename: rename column names via brms:::rename()?
  #   ...: currently ignored
  # Returns:
  #   The design matrix for a regression-like model 
  #   with the specified formula and data. 
  #   For details see the documentation of \code{model.matrix}.
  stopifnot(is.atomic(cols2remove))
  terms <- validate_terms(formula)
  if (is.null(terms)) {
    return(NULL)
  }
  if (isTRUE(attr(terms, "rm_intercept"))) {
    cols2remove <- union(cols2remove, "(Intercept)")
  }
  X <- stats::model.matrix(terms, data)
  cols2remove <- which(colnames(X) %in% cols2remove)
  if (rename) {
    colnames(X) <- rename(colnames(X), check_dup = TRUE) 
  }
  if (length(cols2remove)) {
    X <- X[, -cols2remove, drop = FALSE]
  }
  X
}

arr_design_matrix <- function(Y, r, group)  { 
  # compute the design matrix for ARR effects
  # Args:
  #   Y: a vector containing the response variable
  #   r: ARR order
  #   group: vector containing the grouping variable
  # Notes: 
  #   expects Y to be sorted after group already
  # Returns:
  #   the design matrix for ARR effects
  stopifnot(length(Y) == length(group))
  if (r > 0) {
    U_group <- unique(group)
    N_group <- length(U_group)
    out <- matrix(0, nrow = length(Y), ncol = r)
    ptsum <- rep(0, N_group + 1)
    for (j in seq_len(N_group)) {
      ptsum[j + 1] <- ptsum[j] + sum(group == U_group[j])
      for (i in seq_len(r)) {
        if (ptsum[j] + i + 1 <= ptsum[j + 1]) {
          out[(ptsum[j] + i + 1):ptsum[j + 1], i] <- 
            Y[(ptsum[j] + 1):(ptsum[j + 1] - i)]
        }
      }
    }
  } else {
    out <- NULL
  } 
  out
}

extract_old_standata <- function(x, data, ...) {
  # helper function for validate_newdata to extract
  # old standata required for the computation of new standata
  UseMethod("extract_old_standata")
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
  out <- named_list(names(x$dpars))
  for (i in seq_along(out)) {
    out[[i]] <- extract_old_standata(x$dpars[[i]], data, ...)
  }
  if (has_trials(x$family) || has_cat(x$family)) {
    # trials and ncat should not be computed based on new data
    data_response <- data_response(x, data)
    # partially match via $ to be independent of the response suffix
    out$trials <- data_response$trials
    out$ncat <- data_response$ncat
  }
  if (is.cor_car(x$autocor)) {
    if (isTRUE(nzchar(x$time$group))) {
      out$locations <- levels(factor(get(x$time$group, data)))
    }
  }
  out
}

#' @export
extract_old_standata.btnl <- function(x, data, ...) {
  out <- named_list(names(x$nlpars))
  for (i in seq_along(out)) {
    out[[i]] <- extract_old_standata(x$nlpars[[i]], data, ...)
  }
  out
}

#' @export
extract_old_standata.btl <- function(x, data, ...) {
  list(
    smooths = make_smooth_list(x, data, ...),
    gps = make_gp_list(x, data, ...),
    Jmo = make_Jmo_list(x, data, ...)
  )
}

make_smooth_list <- function(x, data, ...) {
  # extract data related to smooth terms
  # for use in extract_old_standata
  stopifnot(is.btl(x))
  sm_labels <- get_sm_labels(x)
  out <- named_list(sm_labels)
  if (length(sm_labels)) {
    knots <- attr(data, "knots")
    data <- rm_attr(data, "terms")
    gam_args <- list(
      data = data, knots = knots, 
      absorb.cons = TRUE, modCon = 3
    )
    for (i in seq_along(sm_labels)) {
      sc_args <- c(list(eval2(sm_labels[i])), gam_args)
      out[[i]] <- do.call(mgcv::smoothCon, sc_args)
    }
  }
  out
}

make_gp_list <- function(x, data, ...) {
  # extract data related to gaussian processes
  # for use in extract_old_standata
  stopifnot(is.btl(x))
  gpef <- get_gp_labels(x)
  out <- named_list(gpef)
  for (i in seq_along(gpef)) {
    gp <- eval2(gpef[i])
    Xgp <- do.call(cbind, lapply(gp$term, eval2, data))
    out[[i]] <- list(dmax = sqrt(max(diff_quad(Xgp))))
  }
  out
}

make_Jmo_list <- function(x, data, ...) {
  # extract data related to monotonic effects
  # for use in extract_old_standata
  stopifnot(is.btl(x))
  out <- NULL
  if (length(attr(x$sp, "uni_mo"))) {
    # do it like data_sp()
    spef <- tidy_spef(x, data)
    Xmo_fun <- function(x) attr(eval2(x, data), "var")
    Xmo <- lapply(unlist(spef$call_mo), Xmo_fun)
    out <- as.array(ulapply(Xmo, max))
  }
  out
}
