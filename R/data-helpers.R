melt_data <- function(data, family, effects, na.action = na.omit) {
  # melt data frame for multinormal models
  #
  # Args:
  #   data: a data.frame
  #   family: the model family
  #   effects: a named list as returned by extract_effects
  #
  # Returns:
  #   data in long format 
  response <- effects$response
  nresp <- length(response)
  if (nresp == 2 && is.forked(family) || nresp > 1 && is.linear(family)) {
    if (!is(data, "data.frame")) {
      stop("data must be a data.frame for multivariate models", 
           call. = FALSE)
    }
    # only keep variables that are relevant for the model
    rel_vars <- c(all.vars(effects$all), all.vars(effects$respform))
    data <- data[, which(names(data) %in% rel_vars), drop = FALSE]
    if ("trait" %in% names(data)) {
      stop("trait is a resevered variable name in multivariate models",
           call. = FALSE)
    }
    if ("response" %in% names(data)) {
      stop("response is a resevered variable name in multivariate models",
           call. = FALSE)
    }
    nobs <- nrow(data)
    trait <- factor(rep(response, each = nobs), levels = response)
    new_cols <- data.frame(trait = trait)
    # prepare the response variable
    # use na.pass as otherwise cbind will complain
    # when data contains NAs in the response
    temp_mf <- model.frame(effects$respform, data = data, 
                           na.action = na.pass)
    model_response <- model.response(temp_mf)
    # allow to remove NA responses later on
    rows2remove <- which(!complete.cases(model_response))
    if (is.linear(family)) {
      model_response[rows2remove, ] <- NA
      model_response <- as.vector(model_response)
    } else if (is.forked(family)) {
      model_response[rows2remove] <- NA
      reserved <- c(response[2], "main", "spec")
      reserved <- reserved[reserved %in% names(data)]
      if (length(reserved)) {
        stop(paste(paste(reserved, collapse = ", "), 
                   "is a resevered variable name"),
             call. = FALSE)
      }
      one <- rep(1, nobs)
      zero <- rep(0, nobs)
      new_cols$main <- c(one, zero)
      new_cols$spec <- c(zero, one)
      # dummy responses not actually used in Stan
      model_response <- rep(model_response, 2)
    }
    new_cols$response <- model_response
    data <- replicate(length(response), data, simplify = FALSE)
    data <- do.call(na.action, list(cbind(do.call(rbind, data), new_cols)))
  } else if (nresp > 1) {
    stop("invalid multivariate model", call. = FALSE)
  }
  if (isTRUE(attr(effects$fixed, "rsv_intercept"))) {
    if ("intercept" %in% names(data)) {
      stop(paste("intercept is a reserved variable name in models",
                 "without a fixed effects intercept"), call. = FALSE)
    }
    data$intercept <- 1
  }
  data
}  

combine_groups <- function(data, ...) {
  # combine grouping factors
  #
  # Args:
  #   data: a data.frame
  #   ...: the grouping factors to be combined. 
  #
  # Returns:
  #   a data.frame containing all old variables and 
  #   the new combined grouping factors
  group <- c(...)
  if (length(group)) {
    for (i in 1:length(group)) {
      sgroup <- unlist(strsplit(group[[i]], ":"))
      if (length(sgroup) > 1) {
        new.var <- get(sgroup[1], data)
        for (j in 2:length(sgroup)) {
          new.var <- paste0(new.var, "_", get(sgroup[j], data))
        }
        data[[group[[i]]]] <- new.var
      }
    } 
  }
  data
}

update_data <- function(data, family, effects, ..., 
                        na.action = na.omit,
                        drop.unused.levels = TRUE) {
  # update data for use in brm
  #
  # Args:
  #   data: the original data.frame
  #   family: the model family
  #   effects: output of extract_effects (see validate.R)
  #   ...: More formulae passed to combine_groups
  #        Currently only used for autocorrelation structures
  #   na.action: function defining how to treat NAs
  #   drop.unused.levels: indicates if unused factor levels
  #                       should be removed
  #
  # Returns:
  #   model.frame in long format with combined grouping variables if present
  if (is.null(attr(data, "terms")) && "brms.frame" %in% class(data)) {
    # to avoid error described in #30
    # brms.frame class is deprecated as of brms > 0.7.0
    data <- as.data.frame(data)
  }
  if (!(isTRUE(attr(data, "brmsframe")) || "brms.frame" %in% class(data))) {
    data <- melt_data(data, family = family, effects = effects,
                      na.action = na.action)
    data <- model.frame(effects$all, data = data, na.action = na.action,
                        drop.unused.levels = drop.unused.levels)
    if (any(grepl("__", colnames(data))))
      stop("variable names may not contain double underscores '__'",
           call. = FALSE)
    data <- combine_groups(data, effects$random$group, ...)
    attr(data, "brmsframe") <- TRUE
  }
  data
}

amend_newdata <- function(newdata, fit, re_formula = NULL, 
                          allow_new_levels = FALSE,
                          return_standata = TRUE,
                          check_response = FALSE) {
  # amend newdata passed to predict and fitted methods
  # 
  # Args:
  #   newdata: a data.frame containing new data for prediction 
  #   fit: an object of class brmsfit
  #   re_formula: a random effects formula
  #   allow_new_levels: are new random effects levels allowed?
  #   return_standata: logical; compute the data to be passed 
  #                    to Stan, or just return the updated newdata?
  #   check_response: Should response variables be checked
  #                   for existence and validity?
  #
  # Notes:
  #   used in predict.brmsfit, fitted.brmsfit and linear_predictor.brmsfit
  #
  # Returns:
  #   updated data.frame being compatible with fit$formula
  if (is.null(newdata) || is(newdata, "list")) {
    # to shorten expressions in S3 methods such as predict.brmsfit
    if (return_standata && is.null(newdata)) {
      control <- list(keep_intercept = TRUE, save_order = TRUE)
      newdata <- standata(fit, re_formula = re_formula, control = control)
    }
    return(newdata)
  } else if (!"data.frame" %in% class(newdata)) {
    stop("newdata must be a data.frame")
  }
  # standata will be based on an updated formula if re_formula is specified
  new_ranef <- check_re_formula(re_formula, old_ranef = fit$ranef,
                                data = fit$data)
  new_formula <- update_re_terms(fit$formula, re_formula = re_formula)
  new_nonlinear <- lapply(fit$nonlinear, update_re_terms, 
                          re_formula = re_formula)
  et <- extract_time(fit$autocor$formula)
  ee <- extract_effects(new_formula, et$all, family = fit$family,
                        nonlinear = new_nonlinear, resp_rhs_all = FALSE)
  resp_only_vars <- setdiff(all.vars(ee$respform), all.vars(rhs(ee$all)))
  missing_resp <- setdiff(resp_only_vars, names(newdata))
  check_response <- check_response || 
                    (has_arma(fit$autocor) && !use_cov(fit$autocor))
  if (check_response && length(missing_resp)) {
    stop("Response variables must be specified in newdata for this model.",
         call. = FALSE)
  } else {
    for (resp in missing_resp) {
      # add irrelevant response variables but make sure they pass all checks
      newdata[[resp]] <- NA 
    }
  }
  if (is.formula(ee$cens)) {
    for (cens in setdiff(all.vars(ee$cens), names(newdata))) { 
      newdata[[cens]] <- 0 # add irrelevant censor variables
    }
  }
  if (allow_new_levels) {
    # random effects grouping factors do not need to be specified 
    # by the user if new_levels are allowed
    if (length(new_ranef)) {
      all_gf <- unique(unlist(strsplit(names(new_ranef), split = ":")))
      missing_gf <- all_gf[!all_gf %in% names(newdata)]
      newdata[, missing_gf] <- NA
    }
  }
  newdata <- combine_groups(newdata, get_random(ee)$group, et$group)
  # try to validate factor levels in newdata
  if (is.data.frame(fit$data)) {
    # validating is possible (implies brms > 0.5.0)
    list_data <- lapply(as.list(fit$data), function(x)
      if (is.numeric(x)) x else as.factor(x))
    is_factor <- sapply(list_data, is.factor)
    is_group <- names(list_data) %in% names(new_ranef)
    factors <- list_data[is_factor & !is_group]
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
          if (any(!new_levels %in% factor_levels[[i]])) {
            stop(paste("New factor levels are not allowed. \n",
                 "Levels found:", paste(new_levels, collapse = ", ") , "\n",
                 "Levels allowed:", paste(factor_levels[[i]], collapse = ", ")),
                 call. = FALSE)
          }
          newdata[[factor_names[i]]] <- factor(new_factor, factor_levels[[i]])
        }
      }
    }
  } else {
    warning(paste("Validity of factors cannot be checked for", 
                  "fitted model objects created with brms <= 0.5.0"),
            call. = FALSE)
  }
  # validate grouping factors
  if (length(new_ranef)) {
    gnames <- names(new_ranef)
    for (i in seq_along(gnames)) {
      gf <- as.character(get(gnames[i], newdata))
      new_levels <- unique(gf)
      old_levels <- attr(new_ranef[[i]], "levels")
      unknown_levels <- setdiff(new_levels, old_levels)
      if (!allow_new_levels && length(unknown_levels)) {
        stop(paste("levels", paste0(unknown_levels, collapse = ", "), 
                   "of grouping factor", gnames[i], 
                   "not found in the fitted model"), call. = FALSE)
      } 
      if (return_standata) {
        # transform grouping factor levels into their corresponding integers
        # to match the output of make_standata
        newdata[[gnames[i]]] <- sapply(gf, match, table = old_levels)
      }
    }
  }
  if (return_standata) {
    control <- list(is_newdata = TRUE, keep_intercept = TRUE,
                    save_order = TRUE, omit_response = !check_response)
    if (has_trials(fit$family) || has_cat(fit$family)) {
      # if trials or cat are not explicitly part of the formula
      # the will be computed based on the response variable,
      # which should be avoided for newdata
      control[c("trials", "ncat")] <- standata(fit)[c("trials", "ncat")]
    }
    newdata <- make_standata(new_formula, data = newdata, family = fit$family, 
                             autocor = fit$autocor, nonlinear = new_nonlinear,
                             partial = fit$partial, control = control)
  }
  newdata
}

get_model_matrix <- function(formula, data = environment(formula), ...) {
  # Construct Design Matrices for \code{brms} models
  # 
  # Args:
  #   formula: An object of class formula
  #   data: A data frame created with model.frame. 
  #         If another sort of object, model.frame is called first.
  #   ...: Further arguments passed to amend_terms
  # 
  # Returns:
  #   The design matrix for a regression-like model 
  #   with the specified formula and data. 
  #   For details see the documentation of \code{model.matrix}.
  terms <- amend_terms(formula, ...)
  if (is.null(terms)) return(NULL)
  X <- stats::model.matrix(terms, data)
  new_colnames <- rename(colnames(X), check_dup = TRUE)
  if (isTRUE(attr(terms, "rm_intercept")) && "Intercept" %in% new_colnames) {
    X <- X[, -1, drop = FALSE]
    if (ncol(X)) {
      colnames(X) <- new_colnames[2:length(new_colnames)]
    }
  } else colnames(X) <- new_colnames
  X   
}

arr_design_matrix <- function(Y, r, group)  { 
  # calculate design matrix for autoregressive effects of the response
  #
  # Args:
  #   Y: a vector containing the response variable
  #   r: ARR order
  #   group: vector containing the grouping variable for each observation
  #
  # Notes: 
  #   expects Y to be sorted after group already
  # 
  # Returns:
  #   the design matrix for ARR effects
  if (length(Y) != length(group)) 
    stop("Y and group must have the same length")
  if (r > 0) {
    U_group <- unique(group)
    N_group <- length(U_group)
    out <- matrix(0, nrow = length(Y), ncol = r)
    ptsum <- rep(0, N_group + 1)
    for (j in 1:N_group) {
      ptsum[j+1] <- ptsum[j] + sum(group == U_group[j])
      for (i in 1:r) {
        if (ptsum[j] + i + 1 <= ptsum[j + 1])
          out[(ptsum[j] + i + 1):ptsum[j + 1], i] <- 
            Y[(ptsum[j] + 1):(ptsum[j + 1] - i)]
      }
    }
  }
  else out <- NULL
  out
}

.addition <- function(formula, data = NULL) {
  # computes data for addition arguments
  if (!is.formula(formula))
    formula <- as.formula(formula)
  eval(formula[[2]], data, environment(formula))
}

.se <- function(x) {
  # standard errors for meta-analysis
  if (!is.numeric(x)) 
    stop("SEs must be numeric")
  if (min(x) < 0) 
    stop("standard errors must be non-negative", call. = FALSE)
  x  
}

.weights <- function(x) {
  # weights to be applied on any model
  if (!is.numeric(x)) 
    stop("weights must be numeric")
  if (min(x) < 0) 
    stop("weights must be non-negative", call. = FALSE)
  x
}

.disp <- function(x) {
  # dispersion factors
  if (!is.numeric(x)) 
    stop("dispersion factors must be numeric")
  if (min(x) < 0) 
    stop("dispersion factors must be non-negative", call. = FALSE)
  x  
}

.trials <- function(x) {
  # trials for binomial models
  if (any(!is.wholenumber(x) || x < 1))
    stop("number of trials must be positive integers", call. = FALSE)
  x
}

.cat <- function(x) {
  # number of categories for categorical and ordinal models
  if (any(!is.wholenumber(x) || x < 1))
    stop("number of categories must be positive integers", call. = FALSE)
  x
}

.cens <- function(x) {
  # indicator for censoring
  if (is.factor(x)) x <- as.character(x)
  cens <- unname(sapply(x, function(x) {
    if (grepl(paste0("^", x), "right") || isTRUE(x)) x <- 1
    else if (grepl(paste0("^", x), "none") || isFALSE(x)) x <- 0
    else if (grepl(paste0("^", x), "left")) x <- -1
    else x
  }))
  if (!all(unique(cens) %in% c(-1:1)))
    stop (paste0("Invalid censoring data. Accepted values are ", 
                 "'left', 'none', and 'right' \n(abbreviations are allowed) ", 
                 "or -1, 0, and 1. TRUE and FALSE are also accepted \n",
                 "and refer to 'right' and 'none' respectively."),
          call. = FALSE)
  cens
}

.trunc <- function(lb = -Inf, ub = Inf) {
  lb <- as.numeric(lb)
  ub <- as.numeric(ub)
  if (length(lb) != 1 || length(ub) != 1) {
    stop("invalid truncation values", call. = FALSE)
  }
  nlist(lb, ub)
}
  