update_data <- function(data, family, bterms,
                        na.action = na.omit,
                        drop.unused.levels = TRUE,
                        terms_attr = NULL, knots = NULL) {
  # Update data for use in brms functions
  # Args:
  #   data: the original data.frame
  #   family: the model family
  #   bterms: object of class brmsterms
  #   na.action: function defining how to treat NAs
  #   drop.unused.levels: indicates if unused factor levels
  #                       should be removed
  #   terms_attr: a list of attributes of the terms object of 
  #     the original model.frame; only used with newdata;
  #     this ensures that (1) calls to 'poly' work correctly
  #     and (2) that the number of variables matches the number 
  #     of variable names; fixes issue #73
  #   knots: a list of knot values for GAMMS
  # Returns:
  #   model.frame in long format with combined grouping variables if present
  if (missing(data)) {
    stop2("Argument 'data' is missing.")
  }
  if (is.null(attr(data, "terms")) && "brms.frame" %in% class(data)) {
    # to avoid error described in #30
    # brms.frame class is deprecated as of brms > 0.7.0
    data <- as.data.frame(data)
  }
  if (!(isTRUE(attr(data, "brmsframe")) || "brms.frame" %in% class(data))) {
    data <- try(as.data.frame(data, silent = TRUE))
    if (is(data, "try-error")) {
      stop2("Argument 'data' must be coercible to a data.frame.")
    }
    if (!isTRUE(nrow(data) > 0L)) {
      stop2("Argument 'data' does not contain observations.")
    }
    bterms$allvars <- terms(bterms$allvars)
    attributes(bterms$allvars)[names(terms_attr)] <- terms_attr
    if (isTRUE(attr(bterms$formula, "old_mv"))) {
      data <- melt_data(data, family = family, bterms = bterms)
    } else {
      check_data_old_mv(data, family = family, bterms = bterms)
    }
    data <- data_rsv_intercept(data, bterms = bterms)
    missing_vars <- setdiff(all.vars(bterms$allvars), names(data))
    if (length(missing_vars)) {
      stop2("The following variables are missing in 'data':\n",
            collapse_comma(missing_vars))
    }
    data <- model.frame(bterms$allvars, data, na.action = na.pass,
                        drop.unused.levels = drop.unused.levels)
    nrow_with_NA <- nrow(data)
    data <- na.action(data)
    if (nrow(data) != nrow_with_NA) {
      warning2("Rows containing NAs were excluded from the model")
    }
    if (any(grepl("__|_$", colnames(data)))) {
      stop2("Variable names may not contain double underscores ",
            "or underscores at the end.")
    }
    groups <- c(get_re(bterms)$group, bterms$time$group)
    data <- combine_groups(data, groups)
    data <- fix_factor_contrasts(data, ignore = groups)
    attr(data, "knots") <- knots
    attr(data, "brmsframe") <- TRUE
  }
  data
}

melt_data <- function(data, family, bterms) {
  # add reserved variables to the data
  # and transform it into long format for mv models
  # DEPRECATED as of brms 1.0.0
  # Args:
  #   data: a data.frame
  #   family: the model family
  #   bterms: object of class brmsterms
  response <- bterms$response
  nresp <- length(response)
  if (is_mv(family, response = response)) {
    if (!is(data, "data.frame")) {
      stop2("Argument 'data' must be a data.frame for this model.")
    }
    # only keep variables that are relevant for the model
    rel_vars <- c(all.vars(attr(terms(bterms$allvars), "variables")), 
                  all.vars(bterms$respform))
    data <- data[, which(names(data) %in% rel_vars), drop = FALSE]
    rsv_vars <- intersect(c("trait", "response"), names(data))
    if (length(rsv_vars)) {
      rsv_vars <- paste0("'", rsv_vars, "'", collapse = ", ")
      stop2(paste(rsv_vars, "is a reserved variable name."))
    }
    if (is_categorical(family)) {
      # no parameters are modeled for the reference category
      response <- response[-1]
    }
    # prepare the response variable
    # use na.pass as otherwise cbind will complain
    # when data contains NAs in the response
    nobs <- nrow(data)
    trait <- factor(rep(response, each = nobs), levels = response)
    new_cols <- data.frame(trait = trait)
    model_response <- model.response(model.frame(
      bterms$respform, data = data, na.action = na.pass))
    # allow to remove NA responses later on
    rows2remove <- which(!complete.cases(model_response))
    if (is_linear(family)) {
      model_response[rows2remove, ] <- NA
      model_response <- as.vector(model_response)
    } else if (is_categorical(family)) {
      model_response[rows2remove] <- NA
    } else if (is_forked(family)) {
      model_response[rows2remove] <- NA
      rsv_vars <- intersect(c(response[2], "main", "spec"), names(data))
      if (length(rsv_vars)) {
        rsv_vars <- paste0("'", rsv_vars, "'", collapse = ", ")
        stop2(paste(rsv_vars, "is a reserved variable name."))
      }
      one <- rep(1, nobs)
      zero <- rep(0, nobs)
      new_cols$main <- c(one, zero)
      new_cols$spec <- c(zero, one)
      # dummy responses not actually used in Stan
      model_response <- rep(model_response, 2)
    }
    new_cols$response <- model_response
    old_data <- data
    data <- replicate(length(response), old_data, simplify = FALSE)
    data <- cbind(do.call(rbind, data), new_cols)
    data <- fix_factor_contrasts(data, optdata = old_data)
  }
  data
}

check_data_old_mv <- function(data, family, bterms) {
  # check if the deprecated MV syntax was used in a new model
  # Args:
  #   see update_data
  rsv_vars <- rsv_vars(bterms, incl_intercept = FALSE)
  rsv_vars <- setdiff(rsv_vars, names(data))
  used_rsv_vars <- intersect(rsv_vars, all.vars(bterms$allvars))
  if (length(used_rsv_vars)) {
    stop2("It is no longer possible (and necessary) to specify models ", 
          "using the multivariate 'trait' syntax. See help(brmsformula) ",
          "for details on the new syntax.")
  }
  invisible(NULL)
}

data_rsv_intercept <- function(data, bterms) {
  # introduces the resevered variable "intercept" to the data
  # Args:
  #   data: data.frame or list
  #   bterms: object of class brmsterms
  rsv_int <- ulapply(bterms$auxpars, function(x) attr(x$fe, "rsv_intercept"))
  if (any(rsv_int)) {
    if ("intercept" %in% names(data)) {
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
    if (length(sgroup) > 1) {
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
  #            are taken if present
  #   ignore: names of variables for which not to fix contrasts
  # Returns:
  #   a data.frame with amended contrasts attributes
  stopifnot(is(data, "data.frame"))
  stopifnot(is.null(optdata) || is.list(optdata))
  optdata <- as.data.frame(optdata)  # fixes issue #105
  for (i in seq_along(data)) {
    needs_contrast <- is.factor(data[[i]]) && !names(data)[i] %in% ignore
    if (needs_contrast && is.null(attr(data[[i]], "contrasts"))) {
      if (!is.null(attr(optdata[[names(data)[i]]], "contrasts"))) {
        # take contrasts from optdata
        contrasts(data[[i]]) <- attr(optdata[[names(data)[i]]], "contrasts")
      } else if (length(unique(data[[i]])) > 1L) {
        # avoid error when supplying only a single level
        # hard code current global "contrasts" option
        contrasts(data[[i]]) <- contrasts(data[[i]])
      }
    }
  }
  data
}

amend_newdata <- function(newdata, fit, re_formula = NULL, 
                          allow_new_levels = FALSE,
                          check_response = FALSE,
                          only_response = FALSE,
                          incl_autocor = TRUE,
                          return_standata = TRUE) {
  # amend newdata passed to predict and fitted methods
  # Args:
  #   newdata: a data.frame containing new data for prediction 
  #   fit: an object of class brmsfit
  #   re_formula: a group-level formula
  #   allow_new_levels: Are new group-leveld allowed?
  #   check_response: Should response variables be checked
  #                   for existence and validity?
  #   incl_autocor: Check data of autocorrelation terms?
  #   return_standata: Compute the data to be passed to Stan
  #                    or just return the updated newdata?
  # Returns:
  #   updated data.frame being compatible with formula(fit)
  if (!incl_autocor) {
    fit$autocor <- NULL
  }
  if (is.null(newdata)) {
    # to shorten expressions in S3 methods such as predict.brmsfit
    if (return_standata) {
      control <- list(
        not4stan = TRUE, save_order = TRUE,
        omit_response = !check_response,
        only_response = only_response
      )
      newdata <- standata(fit, re_formula = re_formula, control = control)
    } else {
      newdata <- model.frame(fit)
    }
    return(newdata)
  }
  newdata <- try(as.data.frame(newdata, silent = TRUE))
  if (is(newdata, "try-error")) {
    stop2("Argument 'newdata' must be coercible to a data.frame.")
  }
  # standata will be based on an updated formula if re_formula is specified
  new_formula <- update_re_terms(formula(fit), re_formula = re_formula)
  bterms <- parse_bf(new_formula, family = family(fit),
                     autocor = fit$autocor, resp_rhs_all = FALSE)
  resp_only_vars <- setdiff(all.vars(bterms$respform), 
                            all.vars(rhs(bterms$allvars)))
  # fixes #162
  resp_only_vars <- c(resp_only_vars, all.vars(bterms$adforms$dec))
  missing_resp <- setdiff(resp_only_vars, names(newdata))
  if (check_response && length(missing_resp)) {
    stop2("Response variables must be specified in 'newdata'.\n",
          "Missing variables: ", collapse_comma(missing_resp))
  } else {
    for (resp in missing_resp) {
      # add irrelevant response variables but make sure they pass all checks
      newdata[[resp]] <- NA 
    }
  }
  cens_vars <- all.vars(bterms$adforms$cens)
  for (v in setdiff(cens_vars, names(newdata))) {
    # censoring vars are unused in predict and related methods
    newdata[[v]] <- 0
  }
  weights_vars <- all.vars(bterms$adforms$weights)
  for (v in setdiff(weights_vars, names(newdata))) {
    # weighting vars are unused in predict and related methods
    newdata[[v]] <- 1
  }
  new_ranef <- tidy_ranef(bterms, data = model.frame(fit))
  group_vars <- unique(ulapply(new_ranef$gcall, "[[", "groups"))
  if (nrow(fit$ranef)) {
    if (nrow(new_ranef) && allow_new_levels) {
      # grouping factors do not need to be specified 
      # by the user if new levels are allowed
      new_gf <- unique(unlist(strsplit(group_vars, split = ":")))
      missing_gf <- setdiff(new_gf, names(newdata))
      newdata[, missing_gf] <- NA
    }
  }
  newdata <- combine_groups(newdata, group_vars)
  # try to validate factor levels in newdata
  if (is.data.frame(fit$data)) {
    # validating is possible (implies brms > 0.5.0)
    list_data <- lapply(as.list(fit$data), function(x)
      if (is.numeric(x)) x else as.factor(x))
    is_factor <- sapply(list_data, is.factor)
    all_group_vars <- ulapply(fit$ranef$gcall, "[[", "groups")
    dont_check <- c(all_group_vars, cens_vars)
    dont_check <- names(list_data) %in% dont_check
    factors <- list_data[is_factor & !dont_check]
    if (length(factors)) {
      factor_names <- names(factors)
      factor_levels <- lapply(factors, levels) 
      for (i in seq_along(factors)) {
        new_factor <- newdata[[factor_names[i]]]
        if (!is.null(new_factor)) {
          if (!is.factor(new_factor)) {
            new_factor <- factor(new_factor)
          }
          new_levels <- levels(new_factor)
          old_levels <- factor_levels[[i]]
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
            stop2("New factor levels are not allowed.",
                  "\nLevels found: ", paste(new_levels, collapse = ", ") ,
                  "\nLevels allowed: ", paste(old_levels, collapse = ", "))
          }
          newdata[[factor_names[i]]] <- factor(new_factor, old_levels)
          # don't use contrasts(.) here to avoid dimension checks
          attr(newdata[[factor_names[i]]], "contrasts") <- old_contrasts
        }
      }
    }
    # validate monotonic variables
    if (is.formula(bterms$mo)) {
      take_num <- !is_factor & names(list_data) %in% all.vars(bterms$mo)
      # factors have already been checked
      num_mo_vars <- names(list_data)[take_num]
      for (v in num_mo_vars) {
        # use 'get' to check whether v is defined in newdata
        new_values <- get(v, newdata)
        min_value <- min(list_data[[v]])
        invalid <- new_values < min_value | 
                   new_values > max(list_data[[v]]) |
                   !is_wholenumber(new_values)
        if (sum(invalid)) {
          stop2("Invalid values in variable '", v, "': ",
                paste(new_values[invalid], collapse = ","))
        }
        attr(newdata[[v]], "min") <- min_value
      }
    }
    # brms:::update_data expects all original variables to be present
    # even if not actually used later on
    rsv_vars <- rsv_vars(bterms)
    used_vars <- unique(c(names(newdata), all.vars(bterms$allvars), rsv_vars))
    all_vars <- all.vars(str2formula(names(model.frame(fit))))
    unused_vars <- setdiff(all_vars, used_vars)
    if (length(unused_vars)) {
      newdata[, unused_vars] <- NA
    }
  } else {
    warning2("Validity of factors cannot be checked for ", 
             "fitted model objects created with brms <= 0.5.0.")
  }
  # validate grouping factors
  gnames <- unique(new_ranef$group)
  old_levels <- attr(new_ranef, "levels")
  new_levels <- attr(tidy_ranef(bterms, data = newdata), "levels")
  for (g in gnames) {
    unknown_levels <- setdiff(new_levels[[g]], old_levels[[g]])
    if (!allow_new_levels && length(unknown_levels)) {
      unknown_levels <- paste0("'", unknown_levels, "'", collapse = ", ")
      stop2("Levels ", unknown_levels, " of grouping factor '", g, "' ",
            "cannot be found in the fitted model. ",
            "Consider setting argument 'allow_new_levels' to TRUE.")
    }
  }
  if (return_standata) {
    control <- list(
      is_newdata = TRUE, not4stan = TRUE, 
      old_levels = old_levels, save_order = TRUE, 
      omit_response = !check_response,
      only_response = only_response,
      old_cat <- is_old_categorical(fit)
    )
    old_terms <- attr(model.frame(fit), "terms")
    terms_attr <- c("variables", "predvars")
    control$terms_attr <- attributes(old_terms)[terms_attr]
    has_mo <- length(get_effect(bterms, "mo")) > 0L
    if (has_trials(fit$family) || has_cat(fit$family) || has_mo) {
      # some components should not be computed based on newdata
      pars <- c(names(bterms$nlpars), intersect(auxpars(), names(bterms)))
      comp <- c("trials", "ncat", paste0("Jmo", c("", paste0("_", pars))))
      old_standata <- rmNULL(standata(fit)[comp])
      control[c("trials", "ncat")] <- old_standata[c("trials", "ncat")]
      Jmo <- old_standata[grepl("^Jmo", names(old_standata))]
      names(Jmo) <- sub("^Jmo$", "mu", sub("^Jmo_", "", names(Jmo)))
      control[["Jmo"]] <- Jmo 
    }
    control$smooth <- make_smooth_list(bterms, model.frame(fit))
    if (is(fit$autocor, "cov_fixed")) {
      median_V <- median(diag(fit$autocor$V), na.rm = TRUE)
      fit$autocor$V <- diag(median_V, nrow(newdata))
    }
    knots <- attr(model.frame(fit), "knots")
    newdata <- make_standata(
      new_formula, data = newdata, family = fit$family, 
      autocor = fit$autocor, knots = knots, control = control
    )
  }
  newdata
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
  terms <- amend_terms(formula)
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

prepare_mo_vars <- function(formula, data, check = TRUE) {
  # prepare monotonic variables for use in Stan
  # Args:
  #   formula: formula containing mononotic terms
  #   data: a data.frame or named list
  #   check: check the number of levels? 
  # Returns:
  #   'data' with amended monotonic variables
  data <- model.frame(formula, data)
  vars <- names(data)
  for (i in seq_along(vars)) {
    # validate predictors to be modeled as monotonic
    if (is.ordered(data[[vars[i]]])) {
      # counting starts at zero
      data[[vars[i]]] <- as.numeric(data[[vars[i]]]) - 1 
    } else if (all(is_wholenumber(data[[vars[i]]]))) {
      min_value <- attr(data[[vars[i]]], "min")
      if (is.null(min_value)) {
        min_value <- min(data[[vars[i]]])
      }
      data[[vars[i]]] <- data[[vars[i]]] - min_value
    } else {
      stop2("Monotonic predictors must be either integers or ",
            "ordered factors. Error occured for variable '", 
            vars[i], "'.")
    }
    if (check && max(data[[vars[i]]]) < 2L) {
      stop2("Monotonic predictors must have at least 3 different ", 
            "values. Error occured for variable '", vars[i], "'.")
    }
  }
  out <- get_model_matrix(formula, data, cols2remove = "(Intercept)")
  if (any(grepl(":", colnames(out), fixed = TRUE))) {
    stop2("Modeling interactions as monotonic is not meaningful.")
  }
  out
}

#' @export
make_smooth_list.btl <- function(x, data, ...) {
  if (has_smooths(x)) {
    knots <- attr(data, "knots")
    data <- rm_attr(data, "terms")
    gam_args <- list(data = data, knots = knots, 
                     absorb.cons = TRUE, modCon = 3)
    sm_labels <- get_sm_labels(x)
    smooth <- named_list(sm_labels)
    for (i in seq_along(sm_labels)) {
      sc_args <- c(list(eval_smooth(sm_labels[i])), gam_args)
      smooth[[i]] <- do.call(mgcv::smoothCon, sc_args)
    }
  } else {
    smooth <- list()
  }
  smooth
}

#' @export
make_smooth_list.btnl <- function(x, data, ...) {
  smooth <- named_list(names(x$nlpars))
  for (i in seq_along(smooth)) {
    smooth[[i]] <- make_smooth_list(x$nlpars[[i]], data, ...)
  }
  smooth
}

#' @export
make_smooth_list.brmsterms <- function(x, data, ...) {
  smooth <- named_list(names(x$auxpars))
  for (i in seq_along(smooth)) {
    smooth[[i]] <- make_smooth_list(x$auxpars[[i]], data, ...)
  }
  smooth
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
