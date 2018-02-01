# marginal_effects requires quit a few special helper functions
# which are currently only used in this particular method

get_var_combs <- function(..., alist = list()) {
  # get all variable combinations occuring in elements of ...
  # Args:
  #   ...: character vectors or formulas
  #   alist: a list of character vectors or formulas
  dots <- c(list(...), alist)
  for (i in seq_along(dots)) {
    if (is(dots[[i]], "formula")) {
      dots[[i]] <- attr(terms(dots[[i]]), "term.labels")
    }
    dots[[i]] <- lapply(dots[[i]], function(y) all.vars(parse(text = y)))
  }
  unique(unlist(dots, recursive = FALSE))
}

get_all_effects <- function(x, ...) {
  # extract combinations of predictor variables
  UseMethod("get_all_effects")
}

#' @export
get_all_effects.mvbrmsterms <- function(x, ...) {
  out <- lapply(x$terms, get_all_effects, ...)
  unique(unlist(out, recursive = FALSE))
}

#' @export
get_all_effects.brmsterms <- function(x, rsv_vars = NULL, 
                                      comb_all = FALSE) {
  # get all effects for use in marginal_effects
  # Args:
  #   bterms: object of class brmsterms
  #   rsv_vars: character vector of reserved variables
  #   comb_all: include all main effects and two-way interactions?
  # Returns:
  #   a list with one element per valid effect / effects combination
  #   excludes all 3-way or higher interactions
  stopifnot(is.atomic(rsv_vars))
  out <- list()
  for (dp in names(x$dpars)) {
    out <- c(out, get_all_effects(x$dpars[[dp]]))
  }
  out <- rmNULL(lapply(out, setdiff, y = rsv_vars))
  if (comb_all) {
    out <- unique(unlist(out))
    int <- expand.grid(out, out, stringsAsFactors = FALSE)
    int <- int[int[, 1] != int[, 2], ]
    int <- as.list(as.data.frame(t(int), stringsAsFactors = FALSE))
    int <- unique(unname(lapply(int, sort)))
    out <- c(as.list(out), int)
  }
  unique(out[lengths(out) <= 2L])
}

#' @export
get_all_effects.btl <- function(x, ...) {
  int_formula <- function(x) {
    formula(paste("~", paste(x, collapse = "*")))
  }
  covars <- attr(x[["sm"]], "covars")
  byvars <- attr(x[["sm"]], "byvars")
  svars <- mapply(c, covars, byvars, SIMPLIFY = FALSE)
  alist <- lapply(svars, int_formula)
  get_var_combs(x[["fe"]], x[["sp"]], x[["cs"]], x[["gp"]], alist = alist)
}

#' @export
get_all_effects.btnl <- function(x, ...) {
  covars <- all.vars(rhs(x$covars))
  covars_comb <- as.list(covars)
  if (length(covars) > 1L) {
    covars_comb <- c(covars_comb, 
      utils::combn(covars, 2, simplify = FALSE)
    )
  }
  nl_effects <- lapply(x$nlpars, get_all_effects)
  nl_effects <- unlist(nl_effects, recursive = FALSE)
  unique(c(covars_comb, nl_effects))
}

get_int_vars <- function(x, ...) {
  # extract names of variables treated as integers
  UseMethod("get_int_vars")
}

#' @export
get_int_vars.mvbrmsterms <- function(x, ...) {
  unique(ulapply(x$terms, get_int_vars))
}

#' @export
get_int_vars.brmsterms <- function(x, ...) {
  advars <- ulapply(rmNULL(x$adforms[c("trials", "cat")]), all.vars)
  unique(c(advars, get_mo_vars(x)))
}

prepare_conditions <- function(fit, conditions = NULL, effects = NULL,
                               re_formula = NA, rsv_vars = NULL) {
  # prepare conditions for use in marginal_effects
  # Args:
  #   fit: an object of class 'brmsfit'
  #   conditions: optional data.frame containing user defined conditions
  #   effects: see marginal_effects
  #   re_formula: see marginal_effects
  #   rsv_vars: names of reserved variables
  # Returns:
  #   A data.frame with (possibly updated) conditions
  mf <- model.frame(fit)
  new_formula <- update_re_terms(fit$formula, re_formula = re_formula)
  bterms <- parse_bf(new_formula)
  if (any(grepl_expr("^(as\\.)?factor(.+)$", bterms$allvars))) {
    # conditions are chosen based the variables stored in the data
    # this approach cannot take into account possible transformations
    # to factors happening inside the model formula
    warning2(
      "Using 'factor' or 'as.factor' in the model formula ",
      "might lead to problems in 'marginal_effects'.",
      "Please convert your variables to factors beforehand."
    )
  }
  req_vars <- all.vars(rhs(bterms$allvars))
  req_vars <- setdiff(req_vars, rsv_vars)
  if (is.null(conditions)) {
    conditions <- as.data.frame(as.list(rep(NA, length(req_vars))))
    names(conditions) <- req_vars
  } else {
    conditions <- as.data.frame(conditions)
    if (!nrow(conditions)) {
      stop2("Argument 'conditions' must have a least one row.")
    }
    if (any(duplicated(rownames(conditions)))) {
      stop2("Row names of 'conditions' should be unique.")
    }
    conditions <- unique(conditions)
    req_vars <- setdiff(req_vars, names(conditions))
  }
  trial_vars <- all.vars(bterms$adforms$trials)
  if (length(trial_vars)) {
    write_msg <- any(ulapply(trial_vars, function(x) 
      !isTRUE(x %in% names(conditions)) || anyNA(conditions[[x]])
    ))
    if (write_msg) {
      message("Using the median number of trials by ", 
              "default if not specified otherwise.")
    }
  }
  # use default values for unspecified variables
  int_vars <- get_int_vars(bterms)
  for (v in req_vars) {
    if (!is_like_factor(mf[[v]])) {
      # treat variable as numeric
      if (v %in% int_vars) {
        conditions[[v]] <- round(median(mf[[v]]))
      } else {
        conditions[[v]] <- mean(mf[[v]], na.rm = TRUE)
      }
    } else {
      # use reference category for factors
      levels <- attr(as.factor(mf[[v]]), "levels")
      conditions[[v]] <- factor(
        levels[1], levels = levels, ordered = is.ordered(mf[[v]])
      )
    }
  }
  unused_vars <- setdiff(names(conditions), all.vars(bterms$allvars))
  if (length(unused_vars)) {
    warning2(
      "The following variables in 'conditions' are not ", 
      "part of the model:\n", collapse_comma(unused_vars)
    )
  }
  validate_newdata(
    conditions, fit = fit, re_formula = re_formula,
    allow_new_levels = TRUE, incl_autocor = FALSE, 
    return_standata = FALSE
  )
}

prepare_marg_data <- function(data, conditions, int_conditions = NULL,
                              int_vars = NULL, surface = FALSE, 
                              resolution = 100, reorder = TRUE) {
  # prepare data to be used in marginal_effects
  # Args:
  #  data: data.frame containing only data of the predictors of interest
  #  conditions: see argument 'conditions' of marginal_effects
  #  int_conditions: see argument 'int_conditions' of marginal_effects
  #  int_vars: names of variables being treated as integers
  #  surface: generate surface plots later on?
  #  resolution: number of distinct points at which to evaluate
  #              the predictors of interest
  #  reorder: reorder predictors so that numeric ones come first?
  effects <- names(data)
  stopifnot(length(effects) %in% c(1L, 2L))
  pred_types <- ifelse(ulapply(data, is_like_factor), "factor", "numeric")
  # numeric effects should come first
  if (reorder) {
    new_order <- order(pred_types, decreasing = TRUE)
    effects <- effects[new_order]
    pred_types <- pred_types[new_order]
  }
  mono <- effects %in% int_vars
  if (pred_types[1] == "numeric") {
    min1 <- min(data[, effects[1]], na.rm = TRUE)
    max1 <- max(data[, effects[1]], na.rm = TRUE)
    if (mono[1]) {
      values <- seq(min1, max1, by = 1)
    } else {
      values <- seq(min1, max1, length.out = resolution)
    }
  } else {
    values <- unique(data[, effects[1]])
  }
  if (length(effects) == 2L) {
    values <- setNames(list(values, NA), effects)
    if (pred_types[2] == "numeric") {
      if (surface) {
        min2 <- min(data[, effects[2]], na.rm = TRUE)
        max2 <- max(data[, effects[2]], na.rm = TRUE)
        if (mono[2]) {
          values[[2]] <- seq(min2, max2, by = 1)
        } else {
          values[[2]] <- seq(min2, max2, length.out = resolution)
        }
      } else {
        if (effects[2] %in% names(int_conditions)) {
          int_cond <- int_conditions[[effects[2]]]
          if (is.function(int_cond)) {
            int_cond <- int_cond(data[, effects[2]])
          }
          values[[2]] <- int_cond
        } else if (mono[2]) {
          median2 <- median(data[, effects[2]])
          mad2 <- mad(data[, effects[2]])
          values[[2]] <- round((-1:1) * mad2 + median2)
        } else {
          mean2 <- mean(data[, effects[2]], na.rm = TRUE)
          sd2 <- sd(data[, effects[2]], na.rm = TRUE)
          values[[2]] <- (-1:1) * sd2 + mean2
        }
      }
    } else {
      values[[2]] <- unique(data[, effects[2]])
    }
    data <- do.call(expand.grid, values)
  } else {
    data <- structure(data.frame(values), names = effects)
  }
  # no need to have the same value combination more than once
  data <- unique(data)
  data <- data[do.call(order, as.list(data)), , drop = FALSE]
  data <- replicate(nrow(conditions), data, simplify = FALSE)
  marg_vars <- setdiff(names(conditions), effects)
  for (j in seq_len(nrow(conditions))) {
    data[[j]][, marg_vars] <- conditions[j, marg_vars]
    data[[j]]$cond__ <- rownames(conditions)[j]
  }
  data <- do.call(rbind, data)
  data$cond__ <- factor(data$cond__, rownames(conditions))
  structure(data, effects = effects, types = pred_types, mono = mono)
}

marginal_effects_internal <- function(x, ...) {
  # compute fitted values for use in marginal_effects
  UseMethod("marginal_effects_internal")
}

#' @export
marginal_effects_internal.mvbrmsterms <- function(x, resp = NULL, ...) {
  # Returns: a list of summarized prediction matrices
  resp <- validate_resp(resp, x$responses)
  x$terms <- x$terms[resp]
  out <- lapply(x$terms, marginal_effects_internal, ...)
  unlist(out, recursive = FALSE)
}

#' @export
marginal_effects_internal.brmsterms <- function(
  x, fit, marg_data, int_conditions, method, 
  surface, spaghetti, ordinal, probs, robust, ...
) {
  # Returns: a list with the summarized prediction matrix as the only element
  stopifnot(is.brmsfit(fit))
  if (is_categorical(x$family)) {
    stop2("'marginal_effects' is not implemented for categorical models.")
  } else if (is_ordinal(x$family) && !ordinal) {
    warning2(
      "Predictions are treated as continuous variables ",
      "in 'marginal_effects' by default, which is likely ",
      "invalid for family '", x$family$family, "'. ",
      "Consider setting 'ordinal' to TRUE."
    )
  }
  ordinal <- ordinal && is_ordinal(x$family)
  pred_args <- list(
    fit, newdata = marg_data, allow_new_levels = TRUE, 
    incl_autocor = FALSE, summary = FALSE, 
    resp = if (nzchar(x$resp)) x$resp, ...
  )
  out <- do.call(method, pred_args)
  if (is_ordinal(x$family) && !ordinal && method == "fitted") {
    for (k in seq_len(dim(out)[3])) {
      out[, , k] <- out[, , k] * k
    }
    out <- lapply(seq_len(dim(out)[2]), function(s) rowSums(out[, s, ]))
    out <- do.call(cbind, out)
  }
  rownames(marg_data) <- NULL
  eff <- attr(marg_data, "effects")
  types <- attr(marg_data, "types")
  first_numeric <- types[1] %in% "numeric"
  second_numeric <- types[2] %in% "numeric"
  both_numeric <- first_numeric && second_numeric
  if (second_numeric && !surface) {
    # can only be converted to factor after having called method
    mde2 <- round(marg_data[[eff[2]]], 2)
    levels2 <- sort(unique(mde2), TRUE)
    marg_data[[eff[2]]] <- factor(mde2, levels = levels2)
    labels2 <- names(int_conditions[[eff[2]]])
    if (length(labels2) == length(levels2)) {
      levels(marg_data[[eff[2]]]) <- labels2
    }
  }
  spaghetti_data <- NULL
  if (first_numeric && spaghetti) {
    if (surface) {
      stop2("Cannot use 'spaghetti' and 'surface' at the same time.")
    }
    sample <- rep(seq_len(nrow(out)), each = ncol(out))
    if (length(types) == 2L) {
      # samples should be unique across plotting groups
      sample <- paste0(sample, "_", marg_data[[eff[2]]])
    }
    spaghetti_data <- data.frame(as.numeric(t(out)), factor(sample))
    colnames(spaghetti_data) <- c("estimate__", "sample__")
    spaghetti_data <- cbind(marg_data, spaghetti_data)
  }
  if (ordinal) {
    if (method == "fitted") {
      out <- posterior_summary(out, probs = probs, robust = robust)[, 1, ]
    } else if (method == "predict") {
      out <- posterior_table(out)
    }
    cats <- rep(seq_len(ncol(out)), each = nrow(out))
    out <- cbind(estimate__ = as.vector(out), cats__ = cats)
  } else {
    out <- posterior_summary(out, probs = probs, robust = robust)
    colnames(out) <- c("estimate__", "se__", "lower__", "upper__")
  }
  out <- cbind(marg_data, out)
  attr(out, "effects") <- c(eff, if (ordinal) "cats__")
  attr(out, "response") <- as.character(x$formula[2])
  attr(out, "surface") <- unname(both_numeric && surface)
  attr(out, "ordinal") <- ordinal
  attr(out, "spaghetti") <- spaghetti_data
  attr(out, "points") <- make_point_frame(x, fit$data, eff, ...)
  name <- paste0(usc(x$resp, "suffix"), paste0(eff, collapse = ":"))
  setNames(list(out), name)
}

make_point_frame <- function(bterms, mf, effects, conditions, 
                             select_points = 0, transform = NULL, ...) {
  # helper function for marginal_effects
  # allowing add data points to the marginal effects plots
  # Returns:
  #   a data.frame containing the data points to be plotted
  stopifnot(is.brmsterms(bterms), is.data.frame(mf))
  points <- mf[, effects, drop = FALSE]
  points$resp__ <- model.response(
    model.frame(bterms$respform, mf, na.action = na.pass)
  )
  req_vars <- names(mf)
  groups <- get_re(bterms)$group
  if (length(groups)) {
    req_vars <- c(req_vars, unlist(strsplit(groups, ":")))
  }
  req_vars <- unique(setdiff(req_vars, effects))
  req_vars <- intersect(req_vars, names(conditions))
  if (length(req_vars)) {
    # find out which data point is valid for which condition
    mf <- mf[, req_vars, drop = FALSE]
    conditions <- conditions[, req_vars, drop = FALSE]
    points$cond__ <- NA
    points <- replicate(nrow(conditions), points, simplify = FALSE)
    for (i in seq_along(points)) {
      cond <- conditions[i, , drop = FALSE]
      not_na <- !(is.na(cond) | cond == "zero__")
      cond <- cond[, not_na, drop = FALSE]
      mf_tmp <- mf[, not_na, drop = FALSE]
      if (ncol(mf_tmp)) {
        is_num <- sapply(mf_tmp, is.numeric) 
        is_num <- is_num & !names(mf_tmp) %in% groups
        if (sum(is_num)) {
          # handle numeric variables
          stopifnot(select_points >= 0)
          if (select_points > 0) {
            for (v in names(mf_tmp)[is_num]) {
              min <- min(mf_tmp[, v], na.rm = TRUE)
              max <- max(mf_tmp[, v], na.rm = TRUE)
              unit <- scale_unit(mf_tmp[, v], min, max)
              unit_cond <- scale_unit(cond[, v], min, max)
              unit_diff <- abs(unit - unit_cond)
              close_enough <- unit_diff <= select_points
              mf_tmp[close_enough, v] <- cond[, v]
              mf_tmp[!close_enough, v] <- NA
            }
          } else {
            # take all numeric values if select_points is zero
            cond <- cond[, !is_num, drop = FALSE]
            mf_tmp <- mf_tmp[, !is_num, drop = FALSE]
          }
        }
      }
      if (ncol(mf_tmp)) {
        # handle factors and grouping variables
        # do it like base::duplicated
        K <- do.call("paste", c(mf_tmp, sep = "\r")) %in%
          do.call("paste", c(cond, sep = "\r"))
      } else {
        K <- seq_len(nrow(mf))
      }
      # cond__ allows to assign points to conditions
      points[[i]]$cond__[K] <- rownames(conditions)[i]
    }
    points <- do.call(rbind, points)
    points <- points[!is.na(points$cond__), , drop = FALSE]
    points$cond__ <- factor(points$cond__, rownames(conditions))
  }
  if (!is.numeric(points$resp__)) {
    points$resp__ <- as.numeric(as.factor(points$resp__))
    if (is_binary(bterms$family)) {
      points$resp__ <- points$resp__ - 1
    }
  }
  if (!is.null(transform)) {
    points$resp__ <- do.call(transform, list(points$resp__))
  }
  points
}
