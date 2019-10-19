# update data for use in brms functions
# @param data the data passed by the user
# @param bterms object of class brmsterms
# @param na.action function defining how to treat NAs
# @param drop.unused.levels should unused factor levels be removed?
# @param terms_attr a list of attributes of the terms object of 
#   the original model.frame; only used with newdata;
#   this ensures that (1) calls to 'poly' work correctly
#   and (2) that the number of variables matches the number 
#   of variable names; fixes issue #73
# @param knots: a list of knot values for GAMMs
# @return model.frame for use in brms functions
update_data <- function(data, bterms, na.action = na.omit2,
                        drop.unused.levels = TRUE,
                        terms_attr = NULL, knots = NULL) {
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
    stop2("The following variables are missing in 'data':\n",
          collapse_comma(missing_vars))
  }
  for (v in intersect(vars_keep_na(bterms), names(data))) {
    attr(data[[v]], "keep_na") <- TRUE
  }
  data <- model.frame(
    bterms$allvars, data, na.action = na.action,
    drop.unused.levels = drop.unused.levels
  )
  if (any(grepl("__|_$", colnames(data)))) {
    stop2("Variable names may not contain double underscores ",
          "or underscores at the end.")
  }
  groups <- get_group_vars(bterms)
  data <- combine_groups(data, groups)
  data <- fix_factor_contrasts(data, ignore = groups)
  attr(data, "knots") <- knots
  attr(data, "brmsframe") <- TRUE
  data
}

# add the resevered intercept variables to the data
data_rsv_intercept <- function(data, bterms) {
  fe_forms <- get_effect(bterms, "fe")
  if (any(ulapply(fe_forms, no_int))) {
    if (any(data[["intercept"]] != 1)) {
      stop2("Variable name 'intercept' is resevered in models ",
            "without a population-level intercept.")
    }
    if (any(data[["Intercept"]] != 1)) {
      stop2("Variable name 'Intercept' is resevered in models ",
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
  time <- get_autocor_vars(bterms, "time")
  # ordering does not matter for the CAR structure
  group <- get_autocor_vars(bterms, "group", incl_car = FALSE)
  if (length(time) > 1L || length(group) > 1L) {
    stop2("All autocorrelation structures must have the same ",
          "time and group variables.")
  }
  if (length(time) || length(group)) {
    if (length(group)) {
      gv <- data[[group]]
    } else {
      gv <- rep(1, nrow(data))
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

#' Validate New Data
#' 
#' Validate new data passed to post-processing methods of \pkg{brms}. Unless you
#' are a package developer, you will rarely need to call \code{validate_newdata}
#' directly.
#' 
#' @inheritParams extract_draws
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
  if (is.null(newdata)) {
    newdata <- structure(object$data, valid = TRUE, old = TRUE)
  }
  if (isTRUE(attr(newdata, "valid"))) {
    return(newdata)
  }
  newdata <- try(as.data.frame(newdata), silent = TRUE)
  if (is(newdata, "try-error")) {
    stop2("Argument 'newdata' must be coercible to a data.frame.")
  }
  newdata <- rm_attr(newdata, c("terms", "brmsframe"))
  stopifnot(is.brmsfit(object))
  resp <- validate_resp(resp, object)
  if (!incl_autocor) {
    object <- remove_autocor(object) 
  }
  new_formula <- update_re_terms(formula(object), re_formula)
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
  structure(newdata, valid = TRUE)
}

# allows for updating of objects containing new data
# which cannot be passed via argument 'newdata'
# @param x object of class 'brmsfit'
# @param new_objects: optional named list of new objects
# @return a possibly updated 'brmsfit' object
add_new_objects <- function(x, newdata, new_objects = list()) {
  stopifnot(is.brmsfit(x), is.data.frame(newdata))
  # update autocor variables with new objects
  # do not include 'cor_car' here as the adjacency matrix
  # (or subsets of it) should be the same for newdata
  .update_autocor <- function(autocor, resp = "") {
    resp <- usc(resp)
    if (is.cor_sar(autocor)) {
      W_name <- paste0("W", resp)
      if (W_name %in% names(new_objects)) {
        autocor$W <- cor_sar(new_objects[[W_name]])$W
      } else {
        message("Using the identity matrix as weighting matrix by default")
        autocor$W <- diag(nrow(newdata))
      }
    } else if (is.cor_fixed(autocor)) {
      V_name <- paste0("V", resp)
      if (V_name %in% names(new_objects)) {
        autocor$V <- cor_fixed(new_objects[[V_name]])$V
      } else {
        message("Using the median variance by default")
        median_V <- median(diag(autocor$V), na.rm = TRUE)
        autocor$V <- diag(median_V, nrow(newdata)) 
      }
    }
    return(autocor)
  }
  if (!isTRUE(attr(x, "autocor_updated"))) {
    # attribute is set by subset_autocor() to prevent double updating
    if (is_mv(x)) {
      resps <- names(x$formula$forms)
      for (i in seq_along(resps)) {
        new_autocor <- autocor(x, resp = resps[i])
        new_autocor <- .update_autocor(new_autocor, resps[i])
        x$formula$forms[[i]]$autocor <- x$autocor[[i]] <- new_autocor
      }
    } else {
      x$formula$autocor <- x$autocor <- .update_autocor(autocor(x))
    }
  }
  stanvars_data <- subset_stanvars(x$stanvars, block = "data")
  for (name in names(stanvars_data)) {
    if (name %in% names(new_objects)) {
      x$stanvars[[name]]$sdata <- new_objects[[name]]
    }
  }
  x
}

# Construct design matrices for brms models
# @param formula a formula object
# @param data A data frame created with model.frame.
#   If another sort of object, model.frame is called first.
# @param cols2remove names of the columns to remove from 
#   the model matrix; mainly used for intercepts
# @param rename rename column names via rename()?
# @param ... currently ignored
# @return
#   The design matrix for the given formula and data.
#   For details see ?stats::model.matrix
get_model_matrix <- function(formula, data = environment(formula),
                             cols2remove = NULL, rename = TRUE, ...) {
  stopifnot(is.atomic(cols2remove))
  terms <- validate_terms(formula)
  if (is.null(terms)) {
    return(NULL)
  }
  if (no_int(terms)) {
    cols2remove <- union(cols2remove, "(Intercept)")
  }
  X <- stats::model.matrix(terms, data)
  cols2remove <- which(colnames(X) %in% cols2remove)
  if (length(cols2remove)) {
    X <- X[, -cols2remove, drop = FALSE]
  }
  if (rename) {
    colnames(X) <- rename(colnames(X), check_dup = TRUE) 
  }
  X
}

# convenient wrapper around mgcv::PredictMat
PredictMat <- function(object, data, ...) {
  data <- rm_attr(data, "terms")
  out <- mgcv::PredictMat(object, data = data, ...)
  if (length(dim(out)) < 2L) {
    # fixes issue #494
    out <- matrix(out, nrow = 1)
  }
  out
}

# convenient wrapper around mgcv::smoothCon
smoothCon <- function(object, data, ...) {
  data <- rm_attr(data, "terms")
  vars <- setdiff(c(object$term, object$by), "NA")
  for (v in vars) {
    # allow factor-like variables #562
    if (is_like_factor(data[[v]])) {
      data[[v]] <- as.factor(data[[v]])
    }
  }
  mgcv::smoothCon(object, data = data, ...)
}

# Aid prediction from smooths represented as 'type = 2'
# originally provided by Simon Wood 
# @param sm output of mgcv::smoothCon
# @param data new data supplied for prediction
# @return A list of the same structure as returned by mgcv::smoothCon
s2rPred <- function(sm, data) {
  re <- mgcv::smooth2random(sm, names(data), type = 2)
  # prediction matrix for new data
  X <- PredictMat(sm, data)   
  # transform to RE parameterization
  if (!is.null(re$trans.U)) {
    X <- X %*% re$trans.U
  }
  X <- t(t(X) * re$trans.D)
  # re-order columns according to random effect re-ordering
  X[, re$rind] <- X[, re$pen.ind != 0] 
  # re-order penalization index in same way  
  pen.ind <- re$pen.ind
  pen.ind[re$rind] <- pen.ind[pen.ind > 0]
  # start returning the object
  Xf <- X[, which(re$pen.ind == 0), drop = FALSE]
  out <- list(rand = list(), Xf = Xf)
  for (i in seq_along(re$rand)) { 
    # loop over random effect matrices
    out$rand[[i]] <- X[, which(pen.ind == i), drop = FALSE]
    attr(out$rand[[i]], "s.label") <- attr(re$rand[[i]], "s.label")
  }
  names(out$rand) <- names(re$rand)
  out
}

# Basis matrices for baseline hazard functions of the Cox model
# @param y vector of response values
# @param args arguments passed to the spline generating functions
# @param integrate compute the I-spline instead of the M-spline basis?
# @param basis optional precomputed basis matrix
# @return the design matrix of the baseline hazard function
bhaz_basis_matrix <- function(y, args = list(), integrate = FALSE, 
                              basis = NULL) {
  require_package("splines2")
  if (!is.null(basis)) {
    # perform predictions based on an existing basis matrix
    stopifnot(inherits(basis, "mSpline"))
    if (integrate) {
      # for predictions just the attibutes are required
      # which are the same of M-Splines and I-Splines
      class(basis) <- c("matrix", "iSpline")
    }
    return(predict(basis, y))
  }
  stopifnot(is.list(args))
  args$x <- y
  if (!is.null(args$intercept)) {
    args$intercept <- as_one_logical(args$intercept) 
  }
  if (is.null(args$Boundary.knots)) {
    if (isTRUE(args$intercept)) {
      lower_knot <- min(y)
      upper_knot <- max(y)
    } else {
      # we need a smaller lower boundary knot to avoid lp = -Inf 
      # the below choices are ad-hoc and may need further thought
      lower_knot <- max(min(y) - mad(y, na.rm = TRUE) / 10, 0)
      upper_knot <- max(y) + mad(y, na.rm = TRUE) / 10
    }
    args$Boundary.knots <- c(lower_knot, upper_knot)
  }
  if (integrate) {
    out <- do_call(splines2::iSpline, args)
  } else {
    out <- do_call(splines2::mSpline, args)
  }
  out
}

#' Extract response values
#' 
#' Extract response values from a \code{\link{brmsfit}} object.
#' 
#' @param x A \code{\link{brmsfit}} object.
#' @param resp Optional names of response variables for which to extract values.
#' @param warn For internal use only.
#' @param ... Further arguments passed to \code{\link{standata}}.
#' 
#' @return Returns a vector of response values for univariate models and a
#'   matrix of response values with one column per response variable for
#'   multivariate models.
#' 
#' @keywords internal
#' @export
get_y <- function(x, resp = NULL, warn = FALSE, ...) {
  stopifnot(is.brmsfit(x))
  resp <- validate_resp(resp, x)
  warn <- as_one_logical(warn)
  args <- list(x, resp = resp, ...)
  args$re_formula <- NA
  args$check_response <- TRUE
  args$only_response <- TRUE
  args$internal <- TRUE
  sdata <- do_call(standata, args)
  if (warn) {
    if (any(paste0("cens", usc(resp)) %in% names(sdata))) {
      warning2("Results may not be meaningful for censored models.")
    }
  }
  Ynames <- paste0("Y", usc(resp))
  if (length(Ynames) > 1L) {
    out <- do_call(cbind, sdata[Ynames])
    colnames(out) <- resp
  } else {
    out <- sdata[[Ynames]]
  }
  structure(out, old_order = attr(sdata, "old_order"))
}

# coerce censored values into the right format
# @param x vector of censoring indicators
# @return transformed vector of censoring indicators
prepare_cens <- function(x) {
  .prepare_cens <- function(x) {  
    stopifnot(length(x) == 1L)
    regx <- paste0("^", x)
    if (grepl(regx, "left")) {
      x <- -1
    } else if (grepl(regx, "none") || isFALSE(x)) {
      x <- 0
    } else if (grepl(regx, "right") || isTRUE(x)) {
      x <- 1
    } else if (grepl(regx, "interval")) {
      x <- 2
    }
    return(x)
  }
  x <- unname(x)
  if (is.factor(x)) {
    x <- as.character(x)
  }
  ulapply(x, .prepare_cens)
}

# extract information on censoring of the response variable
# @param x a brmsfit object
# @param resp optional names of response variables for which to extract values
# @return vector of censoring indicators or NULL in case of no censoring
get_cens <- function(x, resp = NULL, newdata = NULL) {
  stopifnot(is.brmsfit(x))
  resp <- validate_resp(resp, x, multiple = FALSE)
  bterms <- parse_bf(x$formula)
  if (!is.null(resp)) {
    bterms <- bterms$terms[[resp]]
  }
  if (is.null(newdata)) {
    newdata <- model.frame(x)
  }
  out <- NULL
  if (is.formula(bterms$adforms$cens)) {
    cens <- eval_rhs(bterms$adforms$cens)
    out <- eval2(cens$vars$cens, newdata)
    out <- prepare_cens(out)
  }
  out
}

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
  if (has_trials(x$family) || has_cat(x$family)) {
    # trials and ncat should not be computed based on new data
    data_response <- data_response(
      x, data, check_response = FALSE, not4stan = TRUE
    )
    # partially match via $ to be independent of the response suffix
    out$trials <- data_response$trials
    out$ncat <- data_response$ncat
  }
  if (is_binary(x$family) || is_categorical(x$family)) {
    Y <- model.response(model.frame(x$respform, data, na.action = na.pass))
    out$resp_levels <- levels(as.factor(Y))
  }
  if (is.cor_car(x$autocor)) {
    if (isTRUE(nzchar(x$time$group))) {
      out$locations <- levels(factor(get(x$time$group, data)))
    }
  }
  if (is_cox(x$family)) {
    # compute basis matrix of the baseline hazard for the Cox model
    data_response <- data_response(
      x, data, check_response = FALSE, not4stan = TRUE
    )
    out$bhaz_basis <- bhaz_basis_matrix(
      data_response$y, args = x$family$baseline
    )
  }
  out
}

#' @export
extract_old_standata.btnl <- function(x, data, ...) {
  NULL
}

#' @export
extract_old_standata.btl <- function(x, data, ...) {
  list(
    smooths = make_sm_list(x, data, ...),
    gps = make_gp_list(x, data, ...),
    Jmo = make_Jmo_list(x, data, ...)
  )
}

# extract data related to smooth terms
# for use in extract_old_standata
# @param version optional brms version number
make_sm_list <- function(x, data, version = NULL, ...) {
  stopifnot(is.btl(x))
  smterms <- all_terms(x[["sm"]])
  out <- named_list(smterms)
  if (length(smterms)) {
    knots <- attr(data, "knots")
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
