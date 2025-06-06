#' Display Conditional Effects of Predictors
#'
#' Display conditional effects of one or more numeric and/or categorical
#' predictors including two-way interaction effects.
#'
#' @aliases marginal_effects marginal_effects.brmsfit
#'
#' @param x An object of class \code{brmsfit}.
#' @param effects An optional character vector naming effects (main effects or
#'   interactions) for which to compute conditional plots. Interactions are
#'   specified by a \code{:} between variable names. If \code{NULL} (the
#'   default), plots are generated for all main effects and two-way interactions
#'   estimated in the model. When specifying \code{effects} manually, \emph{all}
#'   two-way interactions (including grouping variables) may be plotted
#'   even if not originally modeled.
#' @param conditions An optional \code{data.frame} containing variable values
#'   to condition on. Each effect defined in \code{effects} will
#'   be plotted separately for each row of \code{conditions}. Values in the
#'   \code{cond__} column will be used as titles of the subplots. If \code{cond__}
#'   is not given, the row names will be used for this purpose instead.
#'   It is recommended to only define a few rows in order to keep the plots clear.
#'   See \code{\link{make_conditions}} for an easy way to define conditions.
#'   If \code{NULL} (the default), numeric variables will be conditionalized by
#'   using their means and factors will get their first level assigned.
#'   \code{NA} values within factors are interpreted as if all dummy
#'   variables of this factor are zero. This allows, for instance, to make
#'   predictions of the grand mean when using sum coding.
#' @param int_conditions An optional named \code{list} whose elements are
#'   vectors of values of the variables specified in \code{effects}.
#'   At these values, predictions are evaluated. The names of
#'   \code{int_conditions} have to match the variable names exactly.
#'   Additionally, the elements of the vectors may be named themselves,
#'   in which case their names appear as labels for the conditions in the plots.
#'   Instead of vectors, functions returning vectors may be passed and are
#'   applied on the original values of the corresponding variable.
#'   If \code{NULL} (the default), predictions are evaluated at the
#'   \eqn{mean} and at \eqn{mean +/- sd} for numeric predictors and at
#'   all categories for factor-like predictors.
#' @param re_formula A formula containing group-level effects to be considered
#'   in the conditional predictions. If \code{NULL}, include all group-level
#'   effects; if \code{NA} (default), include no group-level effects.
#' @param robust If \code{TRUE} (the default) the median is used as the
#'   measure of central tendency. If \code{FALSE} the mean is used instead.
#' @param prob A value between 0 and 1 indicating the desired probability
#'   to be covered by the uncertainty intervals. The default is 0.95.
#' @param probs (Deprecated) The quantiles to be used in the computation of
#'   uncertainty intervals. Please use argument \code{prob} instead.
#' @param method Method used to obtain predictions. Can be set to
#'   \code{"posterior_epred"} (the default), \code{"posterior_predict"},
#'   or \code{"posterior_linpred"}. For more details, see the respective
#'   function documentations.
#' @param spaghetti Logical. Indicates if predictions should
#'   be visualized via spaghetti plots. Only applied for numeric
#'   predictors. If \code{TRUE}, it is recommended
#'   to set argument \code{ndraws} to a relatively small value
#'   (e.g., \code{100}) in order to reduce computation time.
#' @param surface Logical. Indicates if interactions or
#'   two-dimensional smooths should be visualized as a surface.
#'   Defaults to \code{FALSE}. The surface type can be controlled
#'   via argument \code{stype} of the related plotting method.
#' @param categorical Logical. Indicates if effects of categorical
#'   or ordinal models should be shown in terms of probabilities
#'   of response categories. Defaults to \code{FALSE}.
#' @param ordinal (Deprecated) Please use argument \code{categorical}.
#'   Logical. Indicates if effects in ordinal models
#'   should be visualized as a raster with the response categories
#'   on the y-axis. Defaults to \code{FALSE}.
#' @param transform A function or a character string naming
#'   a function to be applied on the predicted responses
#'   before summary statistics are computed. Only allowed
#'   if \code{method = "posterior_predict"}.
#' @param resolution Number of support points used to generate
#'   the plots. Higher resolution leads to smoother plots.
#'   Defaults to \code{100}. If \code{surface} is \code{TRUE},
#'   this implies \code{10000} support points for interaction terms,
#'   so it might be necessary to reduce \code{resolution}
#'   when only few RAM is available.
#' @param too_far Positive number.
#'   For surface plots only: Grid points that are too
#'   far away from the actual data points can be excluded from the plot.
#'   \code{too_far} determines what is too far. The grid is scaled into
#'   the unit square and then grid points more than \code{too_far}
#'   from the predictor variables are excluded. By default, all
#'   grid points are used. Ignored for non-surface plots.
#' @param select_points Positive number.
#'   Only relevant if \code{points} or \code{rug} are set to \code{TRUE}:
#'   Actual data points of numeric variables that
#'   are too far away from the values specified in \code{conditions}
#'   can be excluded from the plot. Values are scaled into
#'   the unit interval and then points more than \code{select_points}
#'   from the values in \code{conditions} are excluded.
#'   By default, all points are used.
#' @param ... Further arguments such as \code{draw_ids} or \code{ndraws}
#'   passed to \code{\link{posterior_predict}} or \code{\link{posterior_epred}}.
#' @inheritParams plot.brmsfit
#' @param ncol Number of plots to display per column for each effect.
#'   If \code{NULL} (default), \code{ncol} is computed internally based
#'   on the number of rows of \code{conditions}.
#' @param points Logical. Indicates if the original data points should be added
#'   via \code{\link[ggplot2:geom_jitter]{geom_jitter}}. Default is
#'   \code{FALSE}. Can be controlled globally via the \code{brms.plot_points}
#'   option. Note that only those data points will be added that match the
#'   specified conditions defined in \code{conditions}. For categorical
#'   predictors, the conditions have to match exactly. For numeric predictors,
#'   argument \code{select_points} is used to determine, which points do match a
#'   condition.
#' @param rug Logical. Indicates if a rug representation of predictor values
#'   should be added via \code{\link[ggplot2:geom_rug]{geom_rug}}. Default is
#'   \code{FALSE}. Depends on \code{select_points} in the same way as
#'   \code{points} does. Can be controlled globally via the \code{brms.plot_rug}
#'   option.
#' @param mean Logical. Only relevant for spaghetti plots.
#'   If \code{TRUE} (the default), display the mean regression
#'   line on top of the regression lines for each sample.
#' @param jitter_width Only used if \code{points = TRUE}:
#'   Amount of horizontal jittering of the data points.
#'   Mainly useful for ordinal models. Defaults to \code{0} that
#'   is no jittering.
#' @param stype Indicates how surface plots should be displayed.
#'   Either \code{"contour"} or \code{"raster"}.
#' @param line_args Only used in plots of continuous predictors:
#'   A named list of arguments passed to
#'   \code{\link[ggplot2:geom_smooth]{geom_smooth}}.
#' @param cat_args Only used in plots of categorical predictors:
#'   A named list of arguments passed to
#'   \code{\link[ggplot2:geom_point]{geom_point}}.
#' @param errorbar_args Only used in plots of categorical predictors:
#'   A named list of arguments passed to
#'   \code{\link[ggplot2:geom_errorbar]{geom_errorbar}}.
#' @param surface_args Only used in surface plots:
#'   A named list of arguments passed to
#'   \code{\link[ggplot2:geom_contour]{geom_contour}} or
#'   \code{\link[ggplot2:geom_raster]{geom_raster}}
#'   (depending on argument \code{stype}).
#' @param spaghetti_args Only used in spaghetti plots:
#'   A named list of arguments passed to
#'   \code{\link[ggplot2:geom_smooth]{geom_smooth}}.
#' @param point_args Only used if \code{points = TRUE}:
#'   A named list of arguments passed to
#'   \code{\link[ggplot2:geom_jitter]{geom_jitter}}.
#' @param rug_args Only used if \code{rug = TRUE}:
#'   A named list of arguments passed to
#'   \code{\link[ggplot2:geom_rug]{geom_rug}}.
#' @param facet_args Only used if if multiple conditions are provided:
#'   A named list of arguments passed to
#'   \code{\link[ggplot2:facet_wrap]{facet_wrap}}.
#'
#' @return An object of class \code{'brms_conditional_effects'} which is a
#'   named list with one data.frame per effect containing all information
#'   required to generate conditional effects plots. Among others, these
#'   data.frames contain some special variables, namely \code{estimate__}
#'   (predicted values of the response), \code{se__} (standard error of the
#'   predicted response), \code{lower__} and \code{upper__} (lower and upper
#'   bounds of the uncertainty interval of the response), as well as
#'   \code{cond__} (used in faceting when \code{conditions} contains multiple
#'   rows).
#'
#'   The corresponding \code{plot} method returns a named
#'   list of \code{\link[ggplot2:ggplot]{ggplot}} objects, which can be further
#'   customized using the \pkg{ggplot2} package.
#'
#' @details When creating \code{conditional_effects} for a particular predictor
#'   (or interaction of two predictors), one has to choose the values of all
#'   other predictors to condition on. By default, the mean is used for
#'   continuous variables and the reference category is used for factors, but
#'   you may change these values via argument \code{conditions}. This also has
#'   an implication for the \code{points} argument: In the created plots, only
#'   those points will be shown that correspond to the factor levels actually
#'   used in the conditioning, in order not to create the false impression of
#'   bad model fit, where it is just due to conditioning on certain factor
#'   levels.
#'
#'   To fully change colors of the created plots, one has to amend both
#'   \code{scale_colour} and \code{scale_fill}. See
#'   \code{\link[ggplot2:scale_color_grey]{scale_colour_grey}} or
#'   \code{\link[ggplot2:scale_color_gradient]{scale_colour_gradient}} for
#'   more details.
#'
#' @examples
#' \dontrun{
#' fit <- brm(count ~ zAge + zBase * Trt + (1 | patient),
#'   data = epilepsy, family = poisson()
#' )
#'
#' ## plot all conditional effects
#' plot(conditional_effects(fit), ask = FALSE)
#'
#' ## change colours to grey scale
#' library(ggplot2)
#' ce <- conditional_effects(fit, "zBase:Trt")
#' plot(ce, plot = FALSE)[[1]] +
#'   scale_color_grey() +
#'   scale_fill_grey()
#'
#' ## only plot the conditional interaction effect of 'zBase:Trt'
#' ## for different values for 'zAge'
#' conditions <- data.frame(zAge = c(-1, 0, 1))
#' plot(conditional_effects(fit,
#'   effects = "zBase:Trt",
#'   conditions = conditions
#' ))
#'
#' ## also incorporate group-level effects variance over patients
#' ## also add data points and a rug representation of predictor values
#' plot(
#'   conditional_effects(fit,
#'     effects = "zBase:Trt",
#'     conditions = conditions, re_formula = NULL
#'   ),
#'   points = TRUE, rug = TRUE
#' )
#'
#' ## change handling of two-way interactions
#' int_conditions <- list(
#'   zBase = setNames(c(-2, 1, 0), c("b", "c", "a"))
#' )
#' conditional_effects(fit,
#'   effects = "Trt:zBase",
#'   int_conditions = int_conditions
#' )
#' conditional_effects(fit,
#'   effects = "Trt:zBase",
#'   int_conditions = list(zBase = quantile)
#' )
#'
#' ## fit a model to illustrate how to plot 3-way interactions
#' fit3way <- brm(count ~ zAge * zBase * Trt, data = epilepsy)
#' conditions <- make_conditions(fit3way, "zAge")
#' conditional_effects(fit3way, "zBase:Trt", conditions = conditions)
#' ## only include points close to the specified values of zAge
#' ce <- conditional_effects(
#'   fit3way, "zBase:Trt",
#'   conditions = conditions,
#'   select_points = 0.1
#' )
#' plot(ce, points = TRUE)
#' }
#'
#' @export
conditional_effects.brmsfit <- function(x, effects = NULL, conditions = NULL,
                                        int_conditions = NULL, re_formula = NA,
                                        prob = 0.95, robust = TRUE,
                                        method = "posterior_epred",
                                        spaghetti = FALSE, surface = FALSE,
                                        categorical = FALSE, ordinal = FALSE,
                                        transform = NULL, resolution = 100,
                                        select_points = 0, too_far = 0,
                                        probs = NULL, ...) {
  probs <- validate_ci_bounds(prob, probs = probs)
  method <- validate_pp_method(method)
  spaghetti <- as_one_logical(spaghetti)
  surface <- as_one_logical(surface)
  categorical <- as_one_logical(categorical)
  ordinal <- as_one_logical(ordinal)
  contains_draws(x)
  x <- restructure(x)
  new_formula <- update_re_terms(x$formula, re_formula = re_formula)
  bterms <- brmsterms(new_formula)

  if (!is.null(transform) && method != "posterior_predict") {
    stop2("'transform' is only allowed if 'method = posterior_predict'.")
  }
  if (ordinal) {
    warning2(
      "Argument 'ordinal' is deprecated. ",
      "Please use 'categorical' instead."
    )
  }
  rsv_vars <- rsv_vars(bterms)
  use_def_effects <- is.null(effects)
  if (use_def_effects) {
    effects <- get_all_effects(bterms, rsv_vars = rsv_vars)
  } else {
    # allow to define interactions in any order
    effects <- strsplit(as.character(effects), split = ":")
    if (any(unique(unlist(effects)) %in% rsv_vars)) {
      stop2(
        "Variables ", collapse_comma(rsv_vars),
        " should not be used as effects for this model"
      )
    }
    if (any(lengths(effects) > 2L)) {
      stop2(
        "To display interactions of order higher than 2 ",
        "please use the 'conditions' argument."
      )
    }
    all_effects <- get_all_effects(
      bterms,
      rsv_vars = rsv_vars, comb_all = TRUE
    )
    ae_coll <- all_effects[lengths(all_effects) == 1L]
    ae_coll <- ulapply(ae_coll, paste, collapse = ":")
    matches <- match(lapply(all_effects, sort), lapply(effects, sort), 0L)
    if (sum(matches) > 0 && sum(matches > 0) < length(effects)) {
      invalid <- effects[setdiff(seq_along(effects), sort(matches))]
      invalid <- ulapply(invalid, paste, collapse = ":")
      warning2(
        "Some specified effects are invalid for this model: ",
        collapse_comma(invalid), "\nValid effects are ",
        "(combinations of): ", collapse_comma(ae_coll)
      )
    }
    effects <- unique(effects[sort(matches)])
    if (!length(effects)) {
      stop2(
        "All specified effects are invalid for this model.\n",
        "Valid effects are (combinations of): ",
        collapse_comma(ae_coll)
      )
    }
  }
  if (categorical || ordinal) {
    int_effs <- lengths(effects) == 2L
    if (any(int_effs)) {
      effects <- effects[!int_effs]
      warning2(
        "Interactions cannot be plotted directly if 'categorical' ",
        "is TRUE. Please use argument 'conditions' instead."
      )
    }
  }
  if (!length(effects)) {
    stop2("No valid effects detected.")
  }
  mf <- model.frame(x)
  conditions <- prepare_conditions(
    x,
    conditions = conditions, effects = effects,
    re_formula = re_formula, rsv_vars = rsv_vars
  )
  int_conditions <- lapply(
    int_conditions,
    function(x) if (is.numeric(x)) sort(x, TRUE) else x
  )
  int_vars <- get_int_vars(bterms)
  group_vars <- get_group_vars(bterms)
  out <- list()
  for (i in seq_along(effects)) {
    eff <- effects[[i]]
    cond_data <- prepare_cond_data(
      mf[, eff, drop = FALSE],
      conditions = conditions,
      int_conditions = int_conditions, int_vars = int_vars,
      group_vars = group_vars, surface = surface,
      resolution = resolution, reorder = use_def_effects
    )
    if (surface && length(eff) == 2L && too_far > 0) {
      # exclude prediction grid points too far from data
      ex_too_far <- mgcv::exclude.too.far(
        g1 = cond_data[[eff[1]]],
        g2 = cond_data[[eff[2]]],
        d1 = mf[, eff[1]],
        d2 = mf[, eff[2]],
        dist = too_far
      )
      cond_data <- cond_data[!ex_too_far, ]
    }
    c(out) <- conditional_effects(
      bterms,
      fit = x, cond_data = cond_data, method = method,
      surface = surface, spaghetti = spaghetti, categorical = categorical,
      ordinal = ordinal, re_formula = re_formula, transform = transform,
      conditions = conditions, int_conditions = int_conditions,
      select_points = select_points, probs = probs, robust = robust,
      ...
    )
  }
  structure(out, class = "brms_conditional_effects")
}

#' @rdname conditional_effects.brmsfit
#' @export
conditional_effects <- function(x, ...) {
  UseMethod("conditional_effects")
}

# compute expected values of MV models for use in conditional_effects
# @return a list of summarized prediction matrices
#' @export
conditional_effects.mvbrmsterms <- function(x, resp = NULL, ...) {
  resp <- validate_resp(resp, x$responses)
  x$terms <- x$terms[resp]
  out <- lapply(x$terms, conditional_effects, ...)
  unlist(out, recursive = FALSE)
}

# conditional_effects for univariate model
# @return a list with the summarized prediction matrix as the only element
# @note argument 'resp' exists only to be excluded from '...' (#589)
#' @export
conditional_effects.brmsterms <- function(
    x, fit, cond_data, int_conditions, method, surface,
    spaghetti, categorical, ordinal, probs, robust,
    dpar = NULL, nlpar = NULL, resp = NULL, ...) {
  stopifnot(is.brmsfit(fit))
  effects <- attr(cond_data, "effects")
  types <- attr(cond_data, "types")
  catscale <- NULL
  pred_args <- list(
    fit,
    newdata = cond_data, allow_new_levels = TRUE,
    dpar = dpar, nlpar = nlpar, resp = if (nzchar(x$resp)) x$resp,
    incl_autocor = FALSE, ...
  )
  if (method != "posterior_predict") {
    # 'transform' creates problems in 'posterior_linpred'
    pred_args$transform <- NULL
  }
  out <- do_call(method, pred_args)
  rownames(cond_data) <- NULL

  if (categorical || ordinal) {
    if (method != "posterior_epred") {
      stop2("Can only use 'categorical' with method = 'posterior_epred'.")
    }
    if (!is_polytomous(x)) {
      stop2(
        "Argument 'categorical' may only be used ",
        "for categorical or ordinal models."
      )
    }
    if (categorical && ordinal) {
      stop2("Please use argument 'categorical' instead of 'ordinal'.")
    }
    catscale <- str_if(is_multinomial(x), "Count", "Probability")
    cats <- dimnames(out)[[3]]
    if (is.null(cats)) cats <- seq_dim(out, 3)
    cond_data <- repl(cond_data, length(cats))
    cond_data <- do_call(rbind, cond_data)
    cond_data$cats__ <- factor(rep(cats, each = ncol(out)), levels = cats)
    effects[2] <- "cats__"
    types[2] <- "factor"
  } else {
    if (conv_cats_dpars(x$family) && is.null(dpar)) {
      stop2("Please set 'categorical' to TRUE.")
    }
    if (is_ordinal(x$family) && is.null(dpar) && method != "posterior_linpred") {
      warning2(
        "Predictions are treated as continuous variables in ",
        "'conditional_effects' by default which is likely invalid ",
        "for ordinal families. Please set 'categorical' to TRUE."
      )
      if (method == "posterior_epred") {
        out <- ordinal_probs_continuous(out)
      }
    }
  }

  cond_data <- add_effects__(cond_data, effects)
  first_numeric <- types[1] %in% "numeric"
  second_numeric <- types[2] %in% "numeric"
  both_numeric <- first_numeric && second_numeric
  if (second_numeric && !surface) {
    # only convert 'effect2__' to factor so that the original
    # second effect variable remains unchanged in the data
    mde2 <- round(cond_data[[effects[2]]], 2)
    levels2 <- sort(unique(mde2), TRUE)
    cond_data$effect2__ <- factor(mde2, levels = levels2)
    labels2 <- names(int_conditions[[effects[2]]])
    if (length(labels2) == length(levels2)) {
      levels(cond_data$effect2__) <- labels2
    }
  }

  spag <- NULL
  if (first_numeric && spaghetti) {
    if (surface) {
      stop2("Cannot use 'spaghetti' and 'surface' at the same time.")
    }
    spag <- out
    if (categorical) {
      spag <- do_call(cbind, array2list(spag))
    }
    sample <- rep(seq_rows(spag), each = ncol(spag))
    if (length(types) == 2L) {
      # draws should be unique across plotting groups
      sample <- paste0(sample, "_", cond_data[[effects[2]]])
    }
    spag <- data.frame(as.numeric(t(spag)), factor(sample))
    colnames(spag) <- c("estimate__", "sample__")
    # ensures that 'cbind' works even in the presence of matrix columns
    cond_data_spag <- repl(cond_data, nrow(spag) / nrow(cond_data))
    cond_data_spag <- Reduce(rbind, cond_data_spag)
    spag <- cbind(cond_data_spag, spag)
  }

  out <- posterior_summary(out, probs = probs, robust = robust)
  if (categorical || ordinal) {
    out <- do_call(rbind, array2list(out))
  }
  colnames(out) <- c("estimate__", "se__", "lower__", "upper__")
  out <- cbind(cond_data, out)
  if (!is.null(dpar)) {
    response <- dpar
  } else if (!is.null(nlpar)) {
    response <- nlpar
  } else {
    response <- as.character(x$formula[2])
  }
  attr(out, "effects") <- effects
  attr(out, "response") <- response
  attr(out, "surface") <- unname(both_numeric && surface)
  attr(out, "categorical") <- categorical
  attr(out, "catscale") <- catscale
  attr(out, "ordinal") <- ordinal
  attr(out, "spaghetti") <- spag
  attr(out, "points") <- make_point_frame(x, fit$data, effects, ...)
  name <- paste0(usc(x$resp, "suffix"), paste0(effects, collapse = ":"))
  setNames(list(out), name)
}

# get combinations of variables used in predictor terms
# @param ... character vectors or formulas
# @param alist a list of character vectors or formulas
get_var_combs <- function(..., alist = list()) {
  dots <- c(list(...), alist)
  for (i in seq_along(dots)) {
    if (is.formula(dots[[i]])) {
      dots[[i]] <- attr(terms(dots[[i]]), "term.labels")
    }
    dots[[i]] <- lapply(dots[[i]], all_vars)
  }
  unique(unlist(dots, recursive = FALSE))
}

# extract combinations of predictor variables
get_all_effects <- function(x, ...) {
  UseMethod("get_all_effects")
}

#' @export
get_all_effects.default <- function(x, ...) {
  NULL
}

#' @export
get_all_effects.mvbrmsterms <- function(x, ...) {
  out <- lapply(x$terms, get_all_effects, ...)
  unique(unlist(out, recursive = FALSE))
}

# get all effects for use in conditional_effects
# @param bterms object of class brmsterms
# @param rsv_vars character vector of reserved variables
# @param comb_all include all main effects and two-way interactions?
# @return a list with one element per valid effect / effects combination
#   excludes all 3-way or higher interactions
#' @export
get_all_effects.brmsterms <- function(x, rsv_vars = NULL, comb_all = FALSE, ...) {
  stopifnot(is_atomic_or_null(rsv_vars))
  out <- list()
  for (dp in names(x$dpars)) {
    out <- c(out, get_all_effects(x$dpars[[dp]]))
  }
  for (nlp in names(x$nlpars)) {
    out <- c(out, get_all_effects(x$nlpars[[nlp]]))
  }
  out <- rmNULL(lapply(out, setdiff, y = rsv_vars))
  if (comb_all) {
    # allow to combine all variables with each other
    out <- unique(unlist(out))
    out <- c(out, get_group_vars(x))
    if (length(out)) {
      int <- expand.grid(out, out, stringsAsFactors = FALSE)
      int <- int[int[, 1] != int[, 2], ]
      int <- as.list(as.data.frame(t(int), stringsAsFactors = FALSE))
      int <- unique(unname(lapply(int, sort)))
      out <- c(as.list(out), int)
    }
  }
  unique(out[lengths(out) <= 2L])
}

#' @export
get_all_effects.btl <- function(x, ...) {
  c(
    get_var_combs(x[["fe"]], x[["cs"]]),
    get_all_effects_type(x, "sp"),
    get_all_effects_type(x, "sm"),
    get_all_effects_type(x, "gp")
  )
}

# extract variable combinations from special terms
get_all_effects_type <- function(x, type) {
  stopifnot(is.btl(x))
  type <- as_one_character(type)
  regex_type <- regex_sp(type)
  terms <- all_terms(x[[type]])
  out <- named_list(terms)
  for (i in seq_along(terms)) {
    # some special terms can appear within interactions
    # we did not allow ":" within these terms so we can use it for splitting
    term_parts <- unlist(strsplit(terms[i], split = ":"))
    vars <- vector("list", length(term_parts))
    for (j in seq_along(term_parts)) {
      matches <- get_matches_expr(regex_type, term_parts[j])
      for (k in seq_along(matches)) {
        # evaluate special terms to extract variables
        tmp <- eval2(matches[[k]])
        c(vars[[j]]) <- setdiff(unique(c(tmp$term, tmp$by)), "NA")
      }
      # extract all variables not part of any special term
      c(vars[[j]]) <- setdiff(all_vars(term_parts[j]), all_vars(matches))
    }
    vars <- unique(unlist(vars))
    out[[i]] <- str2formula(vars, collapse = "*")
  }
  get_var_combs(alist = out)
}

#' @export
get_all_effects.btnl <- function(x, ...) {
  covars <- all_vars(rhs(x$covars))
  out <- as.list(covars)
  if (length(covars) > 1L) {
    c(out) <- utils::combn(covars, 2, simplify = FALSE)
  }
  unique(out)
}

# extract names of predictor variables
get_pred_vars <- function(x) {
  unique(unlist(get_all_effects(x)))
}

# extract names of variables treated as integers
get_int_vars <- function(x, ...) {
  UseMethod("get_int_vars")
}

#' @export
get_int_vars.mvbrmsterms <- function(x, ...) {
  unique(ulapply(x$terms, get_int_vars))
}

#' @export
get_int_vars.brmsterms <- function(x, ...) {
  adterms <- c("trials", "thres", "vint")
  advars <- ulapply(rmNULL(x$adforms[adterms]), all_vars)
  unique(c(advars, get_sp_vars(x, "mo")))
}

# transform posterior draws of ordinal probabilities to a
# continuous scale assuming equidistance between adjacent categories
# @param x an ndraws x nobs x ncat array of posterior draws
# @return an ndraws x nobs matrix of posterior draws
ordinal_probs_continuous <- function(x) {
  stopifnot(length(dim(x)) == 3)
  for (k in seq_dim(x, 3)) {
    x[, , k] <- x[, , k] * k
  }
  x <- lapply(seq_dim(x, 2), function(s) rowSums(x[, s, ]))
  do_call(cbind, x)
}

#' Prepare Fully Crossed Conditions
#'
#' This is a helper function to prepare fully crossed conditions primarily
#' for use with the \code{conditions} argument of \code{\link{conditional_effects}}.
#' Automatically creates labels for each row in the \code{cond__} column.
#'
#' @param x An \R object from which to extract the variables
#'   that should be part of the conditions.
#' @param vars Names of the variables that should be part of the conditions.
#' @param ... Arguments passed to \code{\link{rows2labels}}.
#'
#' @return A \code{data.frame} where each row indicates a condition.
#'
#' @details For factor like variables, all levels are used as conditions.
#'   For numeric variables, \code{mean + (-1:1) * SD} are used as conditions.
#'
#' @seealso \code{\link{conditional_effects}}, \code{\link{rows2labels}}
#'
#' @examples
#' df <- data.frame(x = c("a", "b"), y = rnorm(10))
#' make_conditions(df, vars = c("x", "y"))
#'
#' @export
make_conditions <- function(x, vars, ...) {
  # rev ensures that the last variable varies fastest in expand.grid
  vars <- rev(as.character(vars))
  if (!is.data.frame(x) && "data" %in% names(x)) {
    x <- x$data
  }
  x <- as.data.frame(x)
  out <- named_list(vars)
  for (v in vars) {
    tmp <- get(v, x)
    if (is_like_factor(tmp)) {
      tmp <- levels(as.factor(tmp))
    } else {
      tmp <- mean(tmp, na.rm = TRUE) + (-1:1) * sd(tmp, na.rm = TRUE)
    }
    out[[v]] <- tmp
  }
  out <- rev(expand.grid(out))
  out$cond__ <- rows2labels(out, ...)
  out
}

# extract the cond__ variable used for faceting
get_cond__ <- function(x) {
  out <- x[["cond__"]]
  if (is.null(out)) {
    out <- rownames(x)
  }
  as.character(out)
}

#' Convert Rows to Labels
#'
#' Convert information in rows to labels for each row.
#'
#' @param x A \code{data.frame} for which to extract labels.
#' @param digits Minimal number of decimal places shown in
#'   the labels of numeric variables.
#' @param sep A single character string defining the separator
#'   between variables used in the labels.
#' @param incl_vars Indicates if variable names should
#'   be part of the labels. Defaults to \code{TRUE}.
#' @param ... Currently unused.
#'
#' @return A character vector of the same length as the number
#'   of rows of \code{x}.
#'
#' @seealso \code{\link{make_conditions}}, \code{\link{conditional_effects}}
#'
#' @export
rows2labels <- function(x, digits = 2, sep = " & ", incl_vars = TRUE, ...) {
  x <- as.data.frame(x)
  incl_vars <- as_one_logical(incl_vars)
  out <- x
  for (i in seq_along(out)) {
    if (!is_like_factor(out[[i]])) {
      out[[i]] <- round(out[[i]], digits)
    }
    if (incl_vars) {
      out[[i]] <- paste0(names(out)[i], " = ", out[[i]])
    }
  }
  paste_sep <- function(..., sep__ = sep) {
    paste(..., sep = sep__)
  }
  Reduce(paste_sep, out)
}

# prepare conditions for use in conditional_effects
# @param fit an object of class 'brmsfit'
# @param conditions optional data.frame containing user defined conditions
# @param effects see conditional_effects
# @param re_formula see conditional_effects
# @param rsv_vars names of reserved variables
# @return a data.frame with (possibly updated) conditions
prepare_conditions <- function(fit, conditions = NULL, effects = NULL,
                               re_formula = NA, rsv_vars = NULL) {
  mf <- model.frame(fit)
  new_formula <- update_re_terms(fit$formula, re_formula = re_formula)
  bterms <- brmsterms(new_formula)
  if (any(grepl_expr("^(as\\.)?factor(.+)$", bterms$allvars))) {
    # conditions are chosen based the variables stored in the data
    # this approach cannot take into account possible transformations
    # to factors happening inside the model formula
    warning2(
      "Using 'factor' or 'as.factor' in the model formula ",
      "might lead to problems in 'conditional_effects'.",
      "Please convert your variables to factors beforehand."
    )
  }
  req_vars <- all_vars(rhs(bterms$allvars))
  req_vars <- setdiff(req_vars, rsv_vars)
  req_vars <- setdiff(req_vars, names(fit$data2))
  if (is.null(conditions)) {
    conditions <- as.data.frame(as.list(rep(NA, length(req_vars))))
    names(conditions) <- req_vars
  } else {
    conditions <- as.data.frame(conditions)
    if (!nrow(conditions)) {
      stop2("Argument 'conditions' must have a least one row.")
    }
    conditions <- unique(conditions)
    if (any(duplicated(get_cond__(conditions)))) {
      stop2("Condition labels should be unique.")
    }
    req_vars <- setdiff(req_vars, names(conditions))
  }
  # special treatment for 'trials' addition variables
  trial_vars <- all_vars(bterms$adforms$trials)
  trial_vars <- trial_vars[!vars_specified(trial_vars, conditions)]
  if (length(trial_vars)) {
    message(
      "Setting all 'trials' variables to 1 by ",
      "default if not specified otherwise."
    )
    req_vars <- setdiff(req_vars, trial_vars)
    for (v in trial_vars) {
      conditions[[v]] <- 1L
    }
  }
  # use sensible default values for unspecified variables
  subset_vars <- get_ad_vars(bterms, "subset")
  int_vars <- get_int_vars(bterms)
  group_vars <- get_group_vars(bterms)
  req_vars <- setdiff(req_vars, group_vars)
  for (v in req_vars) {
    if (is_like_factor(mf[[v]])) {
      # factor-like variable
      if (v %in% subset_vars) {
        # avoid unintentional subsetting of newdata (#755)
        conditions[[v]] <- TRUE
      } else {
        # use reference category for factors
        levels <- levels(as.factor(mf[[v]]))
        ordered <- is.ordered(mf[[v]])
        conditions[[v]] <- factor(levels[1], levels, ordered = ordered)
      }
    } else {
      # numeric-like variable
      if (v %in% subset_vars) {
        # avoid unintentional subsetting of newdata (#755)
        conditions[[v]] <- 1
      } else if (v %in% int_vars) {
        # ensure valid integer values
        conditions[[v]] <- round(median(mf[[v]], na.rm = TRUE))
      } else {
        conditions[[v]] <- mean(mf[[v]], na.rm = TRUE)
      }
    }
  }
  all_vars <- c(all_vars(bterms$allvars), "cond__")
  unused_vars <- setdiff(names(conditions), all_vars)
  if (length(unused_vars)) {
    warning2(
      "The following variables in 'conditions' are not ",
      "part of the model:\n", collapse_comma(unused_vars)
    )
  }
  cond__ <- conditions$cond__
  conditions <- validate_newdata(
    conditions, fit,
    re_formula = re_formula,
    allow_new_levels = TRUE, check_response = FALSE,
    incl_autocor = FALSE
  )
  conditions$cond__ <- cond__
  conditions
}

# prepare data to be used in conditional_effects
# @param data data.frame containing only data of the predictors of interest
# @param conditions see argument 'conditions' of conditional_effects
# @param int_conditions see argument 'int_conditions' of conditional_effects
# @param int_vars names of variables being treated as integers
# @param group_vars names of grouping variables
# @param surface generate surface plots later on?
# @param resolution number of distinct points at which to evaluate
#   the predictors of interest
# @param reorder reorder predictors so that numeric ones come first?
prepare_cond_data <- function(data, conditions, int_conditions = NULL,
                              int_vars = NULL, group_vars = NULL,
                              surface = FALSE, resolution = 100,
                              reorder = TRUE) {
  effects <- names(data)
  stopifnot(length(effects) %in% c(1L, 2L))
  is_factor <- ulapply(data, is_like_factor) | names(data) %in% group_vars
  types <- ifelse(is_factor, "factor", "numeric")
  # numeric effects should come first
  if (reorder) {
    new_order <- order(types, decreasing = TRUE)
    effects <- effects[new_order]
    types <- types[new_order]
  }
  # handle first predictor
  if (effects[1] %in% names(int_conditions)) {
    # first predictor has pre-specified conditions
    int_cond <- int_conditions[[effects[1]]]
    if (is.function(int_cond)) {
      int_cond <- int_cond(data[[effects[1]]])
    }
    values <- int_cond
  } else if (types[1] == "factor") {
    # first predictor is factor-like
    values <- factor(unique(data[[effects[1]]]))
  } else {
    # first predictor is numeric
    min1 <- min(data[[effects[1]]], na.rm = TRUE)
    max1 <- max(data[[effects[1]]], na.rm = TRUE)
    if (effects[1] %in% int_vars) {
      values <- seq(min1, max1, by = 1)
    } else {
      values <- seq(min1, max1, length.out = resolution)
    }
  }
  if (length(effects) == 2L) {
    # handle second predictor
    values <- setNames(list(values, NA), effects)
    if (effects[2] %in% names(int_conditions)) {
      # second predictor has pre-specified conditions
      int_cond <- int_conditions[[effects[2]]]
      if (is.function(int_cond)) {
        int_cond <- int_cond(data[[effects[2]]])
      }
      values[[2]] <- int_cond
    } else if (types[2] == "factor") {
      # second predictor is factor-like
      values[[2]] <- factor(unique(data[[effects[2]]]))
    } else {
      # second predictor is numeric
      if (surface) {
        min2 <- min(data[[effects[2]]], na.rm = TRUE)
        max2 <- max(data[[effects[2]]], na.rm = TRUE)
        if (effects[2] %in% int_vars) {
          values[[2]] <- seq(min2, max2, by = 1)
        } else {
          values[[2]] <- seq(min2, max2, length.out = resolution)
        }
      } else {
        if (effects[2] %in% int_vars) {
          median2 <- median(data[[effects[2]]])
          mad2 <- mad(data[[effects[2]]])
          values[[2]] <- round((-1:1) * mad2 + median2)
        } else {
          mean2 <- mean(data[[effects[2]]], na.rm = TRUE)
          sd2 <- sd(data[[effects[2]]], na.rm = TRUE)
          values[[2]] <- (-1:1) * sd2 + mean2
        }
      }
    }
    data <- do_call(expand.grid, values)
  } else {
    stopifnot(length(effects) == 1L)
    data <- structure(data.frame(values), names = effects)
  }
  # no need to have the same value combination more than once
  data <- unique(data)
  data <- data[do_call(order, unname(as.list(data))), , drop = FALSE]
  data <- replicate(nrow(conditions), data, simplify = FALSE)
  cond_vars <- setdiff(names(conditions), effects)
  cond__ <- get_cond__(conditions)
  for (j in seq_rows(conditions)) {
    data[[j]] <- fill_newdata(data[[j]], cond_vars, conditions, n = j)
    data[[j]]$cond__ <- cond__[j]
  }
  data <- do_call(rbind, data)
  data$cond__ <- factor(data$cond__, cond__)
  structure(data, effects = effects, types = types)
}

# which variables in 'vars' are specified in 'data'?
vars_specified <- function(vars, data) {
  .fun <- function(v) isTRUE(v %in% names(data)) && any(!is.na(data[[v]]))
  as.logical(ulapply(vars, .fun))
}

# prepare data points based on the provided conditions
# allows to add data points to conditional effects plots
# @return a data.frame containing the data points to be plotted
make_point_frame <- function(bterms, mf, effects, conditions,
                             select_points = 0, transform = NULL, ...) {
  stopifnot(is.brmsterms(bterms), is.data.frame(mf))
  effects <- intersect(effects, names(mf))
  points <- mf[, effects, drop = FALSE]
  points$resp__ <- model.response(
    model.frame(bterms$respform, mf, na.action = na.pass)
  )
  req_vars <- names(mf)
  groups <- get_re_group_vars(bterms)
  if (length(groups)) {
    c(req_vars) <- unlist(strsplit(groups, ":"))
  }
  req_vars <- unique(setdiff(req_vars, effects))
  req_vars <- intersect(req_vars, names(conditions))
  if (length(req_vars)) {
    # find out which data point is valid for which condition
    cond__ <- get_cond__(conditions)
    mf <- mf[, req_vars, drop = FALSE]
    conditions <- conditions[, req_vars, drop = FALSE]
    points$cond__ <- NA
    points <- replicate(nrow(conditions), points, simplify = FALSE)
    for (i in seq_along(points)) {
      cond <- conditions[i, , drop = FALSE]
      # ensures correct handling of matrix columns
      not_na <- function(x) !any(is.na(x) | x %in% "zero__")
      not_na <- ulapply(cond, not_na)
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
              mf_tmp[[v]][close_enough] <- cond[, v]
              mf_tmp[[v]][!close_enough] <- NA
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
        K <- do_call("paste", c(mf_tmp, sep = "\r")) %in%
          do_call("paste", c(cond, sep = "\r"))
      } else {
        K <- seq_rows(mf)
      }
      # cond__ allows to assign points to conditions
      points[[i]]$cond__[K] <- cond__[i]
    }
    points <- do_call(rbind, points)
    points <- points[!is.na(points$cond__), , drop = FALSE]
    points$cond__ <- factor(points$cond__, cond__)
  }
  points <- add_effects__(points, effects)
  if (!is.numeric(points$resp__)) {
    points$resp__ <- as.numeric(as.factor(points$resp__))
    if (is_binary(bterms$family)) {
      points$resp__ <- points$resp__ - 1
    }
  }
  if (!is.null(transform)) {
    points$resp__ <- do_call(transform, list(points$resp__))
  }
  points
}

# add effect<i>__ variables to the data
add_effects__ <- function(data, effects) {
  for (i in seq_along(effects)) {
    data[[paste0("effect", i, "__")]] <- eval2(effects[i], data)
  }
  data
}

#' @export
print.brms_conditional_effects <- function(x, ...) {
  plot(x, ...)
}

#' @rdname conditional_effects.brmsfit
#' @method plot brms_conditional_effects
#' @importFrom rlang .data
#' @export
plot.brms_conditional_effects <- function(
    x, ncol = NULL, points = getOption("brms.plot_points", FALSE),
    rug = getOption("brms.plot_rug", FALSE), mean = TRUE,
    jitter_width = 0, stype = c("contour", "raster"),
    line_args = list(), cat_args = list(), errorbar_args = list(),
    surface_args = list(), spaghetti_args = list(), point_args = list(),
    rug_args = list(), facet_args = list(), theme = NULL, ask = TRUE,
    plot = TRUE, ...) {
  dots <- list(...)
  plot <- use_alias(plot, dots$do_plot)
  stype <- match.arg(stype)
  smooths_only <- isTRUE(attr(x, "smooths_only"))
  if (points && smooths_only) {
    stop2(
      "Argument 'points' is invalid for objects ",
      "returned by 'conditional_smooths'."
    )
  }
  if (!is_equal(jitter_width, 0)) {
    warning2(
      "'jitter_width' is deprecated. Please use ",
      "'point_args = list(width = <width>)' instead."
    )
  }
  if (!is.null(theme) && !is.theme(theme)) {
    stop2("Argument 'theme' should be a 'theme' object.")
  }
  if (plot) {
    default_ask <- devAskNewPage()
    on.exit(devAskNewPage(default_ask))
    devAskNewPage(ask = FALSE)
  }
  dont_replace <- c("mapping", "data", "inherit.aes")
  plots <- named_list(names(x))
  for (i in seq_along(x)) {
    response <- attr(x[[i]], "response")
    effects <- attr(x[[i]], "effects")
    ncond <- length(unique(x[[i]]$cond__))
    df_points <- attr(x[[i]], "points")
    categorical <- isTRUE(attr(x[[i]], "categorical"))
    catscale <- attr(x[[i]], "catscale")
    surface <- isTRUE(attr(x[[i]], "surface"))
    # deprecated as of brms 2.4.3
    ordinal <- isTRUE(attr(x[[i]], "ordinal"))
    if (surface || ordinal) {
      # surface plots for two dimensional interactions or ordinal plots
      plots[[i]] <- ggplot(x[[i]]) +
        aes(.data[["effect1__"]], .data[["effect2__"]]) +
        labs(x = effects[1], y = effects[2])
      if (ordinal) {
        width <- ifelse(is_like_factor(x[[i]]$effect1__), 0.9, 1)
        .surface_args <- nlist(
          mapping = aes(fill = .data[["estimate__"]]),
          height = 0.9, width = width
        )
        replace_args(.surface_args, dont_replace) <- surface_args
        plots[[i]] <- plots[[i]] +
          do_call(geom_tile, .surface_args) +
          scale_fill_gradientn(colors = viridis6(), name = catscale) +
          ylab(response)
      } else if (stype == "contour") {
        .surface_args <- nlist(
          mapping = aes(
            z = .data[["estimate__"]],
            colour = after_stat(.data[["level"]])
          ),
          bins = 30, linewidth = 1.3
        )
        replace_args(.surface_args, dont_replace) <- surface_args
        plots[[i]] <- plots[[i]] +
          do_call(geom_contour, .surface_args) +
          scale_color_gradientn(colors = viridis6(), name = response)
      } else if (stype == "raster") {
        .surface_args <- nlist(mapping = aes(fill = .data[["estimate__"]]))
        replace_args(.surface_args, dont_replace) <- surface_args
        plots[[i]] <- plots[[i]] +
          do_call(geom_raster, .surface_args) +
          scale_fill_gradientn(colors = viridis6(), name = response)
      }
    } else {
      # plot effects of single predictors or two-way interactions
      gvar <- if (length(effects) == 2L) "effect2__"
      spaghetti <- attr(x[[i]], "spaghetti")
      aes_tmp <- aes(x = .data[["effect1__"]], y = .data[["estimate__"]])
      if (!is.null(gvar)) {
        aes_tmp$colour <- aes(colour = .data[[gvar]])$colour
      }
      plots[[i]] <- ggplot(x[[i]]) +
        aes_tmp +
        labs(x = effects[1], y = response, colour = effects[2])
      if (is.null(spaghetti)) {
        aes_tmp <- aes(ymin = .data[["lower__"]], ymax = .data[["upper__"]])
        if (!is.null(gvar)) {
          aes_tmp$fill <- aes(fill = .data[[gvar]])$fill
        }
        plots[[i]] <- plots[[i]] + aes_tmp +
          labs(fill = effects[2])
      }
      # extract suggested colors for later use
      colors <- ggplot_build(plots[[i]])
      colors <- unique(colors$data[[1]][["colour"]])
      if (points && !categorical && !surface) {
        # add points first so that they appear behind the predictions
        .point_args <- list(
          mapping = aes(x = .data[["effect1__"]], y = .data[["resp__"]]),
          data = df_points, inherit.aes = FALSE,
          size = 2 / ncond^0.25, height = 0, width = jitter_width
        )
        if (is_like_factor(df_points[, gvar])) {
          .point_args$mapping[c("colour", "fill")] <-
            aes(colour = .data[[gvar]], fill = .data[[gvar]])
        }
        replace_args(.point_args, dont_replace) <- point_args
        plots[[i]] <- plots[[i]] +
          do_call(geom_jitter, .point_args)
      }
      if (!is.null(spaghetti)) {
        # add a regression line for each sample separately
        .spaghetti_args <- list(
          aes(group = .data[["sample__"]]),
          data = spaghetti, stat = "identity", linewidth = 0.5
        )
        if (!is.null(gvar)) {
          .spaghetti_args[[1]]$colour <- aes(colour = .data[[gvar]])$colour
        }
        if (length(effects) == 1L) {
          .spaghetti_args$colour <- alpha("blue", 0.1)
        } else {
          # workaround to get transparent lines
          plots[[i]] <- plots[[i]] +
            scale_color_manual(values = alpha(colors, 0.1))
        }
        replace_args(.spaghetti_args, dont_replace) <- spaghetti_args
        plots[[i]] <- plots[[i]] +
          do_call(geom_smooth, .spaghetti_args)
      }
      if (is.numeric(x[[i]]$effect1__)) {
        # line plots for numeric predictors
        .line_args <- list(stat = "identity")
        if (!is.null(spaghetti)) {
          # display a white mean regression line
          if (!is.null(gvar)) {
            .line_args$mapping <- aes(group = .data[[gvar]])
          }
          .line_args$colour <- alpha("white", 0.8)
        }
        replace_args(.line_args, dont_replace) <- line_args
        if (mean || is.null(spaghetti)) {
          plots[[i]] <- plots[[i]] +
            do_call(geom_smooth, .line_args)
        }
        if (rug) {
          .rug_args <- list(
            aes(x = .data[["effect1__"]]),
            sides = "b",
            data = df_points, inherit.aes = FALSE
          )
          if (is_like_factor(df_points[, gvar])) {
            .point_args$mapping[c("colour", "fill")] <-
              aes(colour = .data[[gvar]], fill = .data[[gvar]])
          }
          replace_args(.rug_args, dont_replace) <- rug_args
          plots[[i]] <- plots[[i]] +
            do_call(geom_rug, .rug_args)
        }
      } else {
        # points and errorbars for factors
        .cat_args <- list(
          position = position_dodge(width = 0.4),
          size = 4 / ncond^0.25
        )
        .errorbar_args <- list(
          position = position_dodge(width = 0.4),
          width = 0.3
        )
        replace_args(.cat_args, dont_replace) <- cat_args
        replace_args(.errorbar_args, dont_replace) <- errorbar_args
        plots[[i]] <- plots[[i]] +
          do_call(geom_point, .cat_args) +
          do_call(geom_errorbar, .errorbar_args)
      }
      if (categorical) {
        plots[[i]] <- plots[[i]] + ylab(catscale) +
          labs(fill = response, color = response)
      }
    }
    if (ncond > 1L) {
      # one plot per row of conditions
      if (is.null(ncol)) {
        ncol <- max(floor(sqrt(ncond)), 3)
      }
      .facet_args <- nlist(facets = "cond__", ncol)
      replace_args(.facet_args, dont_replace) <- facet_args
      plots[[i]] <- plots[[i]] +
        do_call(facet_wrap, .facet_args)
    }
    plots[[i]] <- plots[[i]] + theme
    if (plot) {
      plot(plots[[i]])
      if (i == 1) {
        devAskNewPage(ask = ask)
      }
    }
  }
  invisible(plots)
}

# the name 'marginal_effects' is deprecated as of brms 2.10.3
# do not remove it eventually as it has been used in the brms papers
#' @export
marginal_effects <- function(x, ...) {
  UseMethod("marginal_effects")
}

#' @export
marginal_effects.brmsfit <- function(x, ...) {
  warning2(
    "Method 'marginal_effects' is deprecated. ",
    "Please use 'conditional_effects' instead."
  )
  conditional_effects.brmsfit(x, ...)
}

#' @export
print.brmsMarginalEffects <- function(x, ...) {
  class(x) <- "brms_conditional_effects"
  print(x, ...)
}

#' @export
plot.brmsMarginalEffects <- function(x, ...) {
  class(x) <- "brms_conditional_effects"
  plot(x, ...)
}
