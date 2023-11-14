#' Efficient approximate leave-one-out cross-validation (LOO)
#'
#' Perform approximate leave-one-out cross-validation based
#' on the posterior likelihood using the \pkg{loo} package.
#' For more details see \code{\link[loo:loo]{loo}}.
#'
#' @aliases loo LOO LOO.brmsfit
#'
#' @param x A \code{brmsfit} object.
#' @param ... More \code{brmsfit} objects or further arguments
#'   passed to the underlying post-processing functions.
#'   In particular, see \code{\link{prepare_predictions}} for further
#'   supported arguments.
#' @param compare A flag indicating if the information criteria
#'  of the models should be compared to each other
#'  via \code{\link{loo_compare}}.
#' @param pointwise A flag indicating whether to compute the full
#'  log-likelihood matrix at once or separately for each observation.
#'  The latter approach is usually considerably slower but
#'  requires much less working memory. Accordingly, if one runs
#'  into memory issues, \code{pointwise = TRUE} is the way to go.
#' @param moment_match Logical; Indicate whether \code{\link{loo_moment_match}}
#'  should be applied on problematic observations. Defaults to \code{FALSE}.
#'  For most models, moment matching will only work if you have set
#'  \code{save_pars = save_pars(all = TRUE)} when fitting the model with
#'  \code{\link{brm}}. See \code{\link{loo_moment_match.brmsfit}} for more
#'  details.
#' @param reloo Logical; Indicate whether \code{\link{reloo}}
#'  should be applied on problematic observations. Defaults to \code{FALSE}.
#' @param k_threshold The threshold at which pareto \eqn{k}
#'   estimates are treated as problematic. Defaults to \code{0.7}.
#'   Only used if argument \code{reloo} is \code{TRUE}.
#'   See \code{\link[loo:pareto-k-diagnostic]{pareto_k_ids}} for more details.
#' @param save_psis Should the \code{"psis"} object created internally be saved
#'   in the returned object? For more details see \code{\link[loo:loo]{loo}}.
#' @param moment_match_args Optional \code{list} of additional arguments passed to
#'   \code{\link{loo_moment_match}}.
#' @param reloo_args Optional \code{list} of additional arguments passed to
#'   \code{\link{reloo}}.
#' @param model_names If \code{NULL} (the default) will use model names
#'   derived from deparsing the call. Otherwise will use the passed
#'   values as model names.
#' @inheritParams predict.brmsfit
#'
#' @details See \code{\link{loo_compare}} for details on model comparisons.
#'  For \code{brmsfit} objects, \code{LOO} is an alias of \code{loo}.
#'  Use method \code{\link{add_criterion}} to store
#'  information criteria in the fitted model object for later usage.
#'
#' @return If just one object is provided, an object of class \code{loo}.
#'  If multiple objects are provided, an object of class \code{loolist}.
#'
#' @examples
#' \dontrun{
#' # model with population-level effects only
#' fit1 <- brm(rating ~ treat + period + carry,
#'             data = inhaler)
#' (loo1 <- loo(fit1))
#'
#' # model with an additional varying intercept for subjects
#' fit2 <- brm(rating ~ treat + period + carry + (1|subject),
#'             data = inhaler)
#' (loo2 <- loo(fit2))
#'
#' # compare both models
#' loo_compare(loo1, loo2)
#' }
#'
#' @references
#' Vehtari, A., Gelman, A., & Gabry J. (2016). Practical Bayesian model
#' evaluation using leave-one-out cross-validation and WAIC. In Statistics
#' and Computing, doi:10.1007/s11222-016-9696-4. arXiv preprint arXiv:1507.04544.
#'
#' Gelman, A., Hwang, J., & Vehtari, A. (2014).
#' Understanding predictive information criteria for Bayesian models.
#' Statistics and Computing, 24, 997-1016.
#'
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation
#' and widely applicable information criterion in singular learning theory.
#' The Journal of Machine Learning Research, 11, 3571-3594.
#'
#' @importFrom loo loo is.loo
#' @export loo
#' @export
loo.brmsfit <-  function(x, ..., compare = TRUE, resp = NULL,
                         pointwise = FALSE, moment_match = FALSE,
                         reloo = FALSE, k_threshold = 0.7, save_psis = FALSE,
                         moment_match_args = list(), reloo_args = list(),
                         model_names = NULL) {
  args <- split_dots(x, ..., model_names = model_names)
  c(args) <- nlist(
    criterion = "loo", pointwise, compare,
    resp, k_threshold, save_psis, moment_match,
    reloo, moment_match_args, reloo_args
  )
  do_call(compute_loolist, args)
}

#' @export
LOO.brmsfit <- function(x, ..., compare = TRUE, resp = NULL,
                        pointwise = FALSE, moment_match = FALSE,
                        reloo = FALSE, k_threshold = 0.7, save_psis = FALSE,
                        moment_match_args = list(), reloo_args = list(),
                        model_names = NULL) {
  cl <- match.call()
  cl[[1]] <- quote(loo)
  eval(cl, parent.frame())
}

#' @export
LOO <- function(x, ...) {
  UseMethod("LOO")
}

#' Widely Applicable Information Criterion (WAIC)
#'
#' Compute the widely applicable information criterion (WAIC)
#' based on the posterior likelihood using the \pkg{loo} package.
#' For more details see \code{\link[loo:waic]{waic}}.
#'
#' @aliases waic WAIC WAIC.brmsfit
#'
#' @inheritParams loo.brmsfit
#'
#' @details See \code{\link{loo_compare}} for details on model comparisons.
#'  For \code{brmsfit} objects, \code{WAIC} is an alias of \code{waic}.
#'  Use method \code{\link[brms:add_criterion]{add_criterion}} to store
#'  information criteria in the fitted model object for later usage.
#'
#' @return If just one object is provided, an object of class \code{loo}.
#'  If multiple objects are provided, an object of class \code{loolist}.
#'
#' @examples
#' \dontrun{
#' # model with population-level effects only
#' fit1 <- brm(rating ~ treat + period + carry,
#'             data = inhaler)
#' (waic1 <- waic(fit1))
#'
#' # model with an additional varying intercept for subjects
#' fit2 <- brm(rating ~ treat + period + carry + (1|subject),
#'             data = inhaler)
#' (waic2 <- waic(fit2))
#'
#' # compare both models
#' loo_compare(waic1, waic2)
#' }
#'
#' @references
#' Vehtari, A., Gelman, A., & Gabry J. (2016). Practical Bayesian model
#' evaluation using leave-one-out cross-validation and WAIC. In Statistics
#' and Computing, doi:10.1007/s11222-016-9696-4. arXiv preprint arXiv:1507.04544.
#'
#' Gelman, A., Hwang, J., & Vehtari, A. (2014).
#' Understanding predictive information criteria for Bayesian models.
#' Statistics and Computing, 24, 997-1016.
#'
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation
#' and widely applicable information criterion in singular learning theory.
#' The Journal of Machine Learning Research, 11, 3571-3594.
#'
#' @importFrom loo waic
#' @export waic
#' @export
waic.brmsfit <- function(x, ..., compare = TRUE, resp = NULL,
                         pointwise = FALSE, model_names = NULL) {
  args <- split_dots(x, ..., model_names = model_names)
  c(args) <- nlist(criterion = "waic", pointwise, compare, resp)
  do_call(compute_loolist, args)
}

#' @export
WAIC.brmsfit <- function(x, ..., compare = TRUE, resp = NULL,
                         pointwise = FALSE, model_names = NULL) {
  cl <- match.call()
  cl[[1]] <- quote(waic)
  eval(cl, parent.frame())
}

#' @export
WAIC <- function(x, ...) {
  UseMethod("WAIC")
}

# helper function used to create (lists of) 'loo' objects
# @param models list of brmsfit objects
# @param criterion name of the criterion to compute
# @param use_stored use precomputed criterion objects if possible?
# @param compare compare models using 'loo_compare'?
# @param ... more arguments passed to compute_loo
# @return If length(models) > 1 an object of class 'loolist'
#   If length(models) == 1 an object of class 'loo'
compute_loolist <- function(models, criterion, use_stored = TRUE,
                            compare = TRUE, ...) {
  criterion <- match.arg(criterion, loo_criteria())
  args <- nlist(criterion, ...)
  for (i in seq_along(models)) {
    models[[i]] <- restructure(models[[i]])
  }
  if (length(models) > 1L) {
    if (!match_nobs(models)) {
      stop2("Models have different number of observations.")
    }
    if (length(use_stored) == 1L) {
      use_stored <- rep(use_stored, length(models))
    }
    out <- list(loos = named_list(names(models)))
    for (i in seq_along(models)) {
      args$x <- models[[i]]
      args$model_name <- names(models)[i]
      args$use_stored <- use_stored[i]
      out$loos[[i]] <- do_call(compute_loo, args)
    }
    compare <- as_one_logical(compare)
    if (compare) {
      out$diffs <- loo_compare(out$loos)
      # for backwards compatibility; remove in brms 3.0
      out$ic_diffs__ <- SW(compare_ic(x = out$loos))$ic_diffs__
    }
    class(out) <- "loolist"
  } else {
    args$x <- models[[1]]
    args$model_name <- names(models)
    args$use_stored <- use_stored
    out <- do_call(compute_loo, args)
  }
  out
}

# compute model fit criteria using the 'loo' package
# @param x an object of class brmsfit
# @param criterion the criterion to be computed
# @param newdata optional data.frame of new data
# @param resp optional names of the predicted response variables
# @param model_name original variable name of object 'x'
# @param use_stored use precomputed criterion objects if possible?
# @param ... passed to the individual methods
# @return an object of class 'loo'
compute_loo <- function(x, criterion, newdata = NULL, resp = NULL,
                        model_name = "", use_stored = TRUE, ...) {
  criterion <- match.arg(criterion, loo_criteria())
  model_name <- as_one_character(model_name)
  use_stored <- as_one_logical(use_stored)
  out <- get_criterion(x, criterion)
  if (!(use_stored && is.loo(out))) {
    args <- nlist(x, newdata, resp, model_name, ...)
    out <- do_call(paste0(".", criterion), args)
    attr(out, "yhash") <- hash_response(x, newdata = newdata, resp = resp)
  }
  attr(out, "model_name") <- model_name
  out
}

# possible criteria to evaluate via the loo package
loo_criteria <- function() {
  c("loo", "waic", "psis", "kfold", "loo_subsample")
}

# compute 'loo' criterion using the 'loo' package
.loo <- function(x, pointwise, k_threshold, moment_match, reloo,
                 moment_match_args, reloo_args, newdata,
                 resp, model_name, save_psis, ...) {
  loo_args <- prepare_loo_args(
    x, newdata = newdata, resp = resp,
    pointwise = pointwise, save_psis = save_psis,
    ...
  )
  out <- SW(do_call("loo", loo_args, pkg = "loo"))
  if (moment_match) {
    c(moment_match_args) <- nlist(
      x, loo = out, newdata, resp,
      k_threshold, check = FALSE, ...
    )
    out <- do_call("loo_moment_match", moment_match_args)
  }
  if (reloo) {
    c(reloo_args) <- nlist(
      x, loo = out, newdata, resp,
      k_threshold, check = FALSE, ...
    )
    out <- do_call("reloo", reloo_args)
  }
  recommend_loo_options(out, k_threshold, moment_match, model_name)
  out
}

# compute 'waic' criterion using the 'loo' package
# @param model_name ignored but included to avoid being passed to '...'
.waic <- function(x, pointwise, newdata, resp, model_name, ...) {
  loo_args <- prepare_loo_args(
    x, newdata = newdata, resp = resp,
    pointwise = pointwise, ...
  )
  do_call("waic", loo_args, pkg = "loo")
}

# alias of psis for convenient use in compute_loo()
.psis <- function(x, newdata, resp, model_name, ...) {
  psis(x, newdata = newdata, resp = resp, model_name = model_name, ...)
}

#' @inherit loo::psis return title description details seealso references
#'
#' @aliases psis psis.brmsfit
#'
#' @param log_ratios A fitted model object of class \code{brmsfit}.
#'   Argument is named "log_ratios" to match the argument name of the
#'   \code{\link[loo:psis]{loo::psis}} generic function.
#' @param model_name Currently ignored.
#' @param ... Further arguments passed to \code{\link{log_lik}} and
#'   \code{\link[loo:psis]{loo::psis}}.
#' @inheritParams log_lik.brmsfit
#'
#' @examples
#' \dontrun{
#' fit <- brm(rating ~ treat + period + carry, data = inhaler)
#' psis(fit)
#'}
#' @importFrom loo psis
#' @export psis
#' @export
psis.brmsfit <- function(log_ratios, newdata = NULL, resp = NULL,
                         model_name, ...) {
  loo_args <- prepare_loo_args(
    log_ratios, newdata = newdata, resp = resp,
    pointwise = FALSE, ...
  )
  loo_args$log_ratios <- -loo_args$x
  loo_args$x <- NULL
  do_call("psis", loo_args, pkg = "loo")
}

# prepare arguments passed to the methods of the `loo` package
prepare_loo_args <- function(x, newdata, resp, pointwise, ...) {
  pointwise <- as_one_logical(pointwise)
  loo_args <- list(...)
  ll_args <- nlist(object = x, newdata, resp, pointwise, ...)
  loo_args$x <- do_call(log_lik, ll_args)
  if (pointwise) {
    loo_args$draws <- attr(loo_args$x, "draws")
    loo_args$data <- attr(loo_args$x, "data")
  }
  # compute pointwise relative efficiencies
  r_eff_args <- loo_args
  r_eff_args$fit <- x
  loo_args$r_eff <- do_call(r_eff_log_lik, r_eff_args)
  loo_args
}

#' Model comparison with the \pkg{loo} package
#'
#' For more details see \code{\link[loo:loo_compare]{loo_compare}}.
#'
#' @aliases loo_compare
#'
#' @inheritParams loo.brmsfit
#' @param ... More \code{brmsfit} objects.
#' @param criterion The name of the criterion to be extracted
#'   from \code{brmsfit} objects.
#'
#' @details All \code{brmsfit} objects should contain precomputed
#'   criterion objects. See \code{\link{add_criterion}} for more help.
#'
#' @return An object of class "\code{compare.loo}".
#'
#' @examples
#' \dontrun{
#' # model with population-level effects only
#' fit1 <- brm(rating ~ treat + period + carry,
#'             data = inhaler)
#' fit1 <- add_criterion(fit1, "waic")
#'
#' # model with an additional varying intercept for subjects
#' fit2 <- brm(rating ~ treat + period + carry + (1|subject),
#'             data = inhaler)
#' fit2 <- add_criterion(fit2, "waic")
#'
#' # compare both models
#' loo_compare(fit1, fit2, criterion = "waic")
#' }
#'
#' @importFrom loo loo_compare
#' @export loo_compare
#' @export
loo_compare.brmsfit <- function(x, ..., criterion = c("loo", "waic", "kfold"),
                                model_names = NULL) {
  criterion <- match.arg(criterion)
  models <- split_dots(x, ..., model_names = model_names, other = FALSE)
  loos <- named_list(names(models))
  for (i in seq_along(models)) {
    models[[i]] <- restructure(models[[i]])
    loos[[i]] <- get_criterion(models[[i]], criterion)
    if (is.null(loos[[i]])) {
      stop2(
        "Model '", names(models)[i], "' does not contain a precomputed '",
        criterion, "' criterion. See ?loo_compare.brmsfit for help."
      )
    }
  }
  loo_compare(loos)
}

#' Model averaging via stacking or pseudo-BMA weighting.
#'
#' Compute model weights for \code{brmsfit} objects via stacking
#' or pseudo-BMA weighting. For more details, see
#' \code{\link[loo:loo_model_weights]{loo::loo_model_weights}}.
#'
#' @aliases loo_model_weights
#'
#' @inheritParams loo.brmsfit
#'
#' @return A named vector of model weights.
#'
#' @examples
#' \dontrun{
#' # model with population-level effects only
#' fit1 <- brm(rating ~ treat + period + carry,
#'             data = inhaler, family = "gaussian")
#' # model with an additional varying intercept for subjects
#' fit2 <- brm(rating ~ treat + period + carry + (1|subject),
#'             data = inhaler, family = "gaussian")
#' loo_model_weights(fit1, fit2)
#' }
#'
#' @method loo_model_weights brmsfit
#' @importFrom loo loo_model_weights
#' @export loo_model_weights
#' @export
loo_model_weights.brmsfit <- function(x, ..., model_names = NULL) {
  args <- split_dots(x, ..., model_names = model_names)
  models <- args$models
  args$models <- NULL
  log_lik_list <- lapply(models, function(x)
    do_call(log_lik, c(list(x), args))
  )
  args$x <- log_lik_list
  args$r_eff_list <- mapply(
    r_eff_log_lik, log_lik_list,
    fit = models, SIMPLIFY = FALSE
  )
  out <- do_call(loo::loo_model_weights, args)
  names(out) <- names(models)
  out
}

#' Add model fit criteria to model objects
#'
#' @param x An \R object typically of class \code{brmsfit}.
#' @param criterion Names of model fit criteria
#'   to compute. Currently supported are \code{"loo"},
#'   \code{"waic"}, \code{"kfold"}, \code{"loo_subsample"},
#'   \code{"bayes_R2"} (Bayesian R-squared),
#'   \code{"loo_R2"} (LOO-adjusted R-squared), and
#'   \code{"marglik"} (log marginal likelihood).
#' @param model_name Optional name of the model. If \code{NULL}
#'   (the default) the name is taken from the call to \code{x}.
#' @param overwrite Logical; Indicates if already stored fit
#'   indices should be overwritten. Defaults to \code{FALSE}.
#' @param file Either \code{NULL} or a character string. In the latter case, the
#'   fitted model object including the newly added criterion values is saved via
#'   \code{\link{saveRDS}} in a file named after the string supplied in
#'   \code{file}. The \code{.rds} extension is added automatically. If \code{x}
#'   was already stored in a file before, the file name will be reused
#'   automatically (with a message) unless overwritten by \code{file}. In any
#'   case, \code{file} only applies if new criteria were actually added via
#'   \code{add_criterion} or if \code{force_save} was set to \code{TRUE}.
#' @param force_save Logical; only relevant if \code{file} is specified and
#'   ignored otherwise. If \code{TRUE}, the fitted model object will be saved
#'   regardless of whether new criteria were added via \code{add_criterion}.
#' @param ... Further arguments passed to the underlying
#'   functions computing the model fit criteria.
#'
#' @return An object of the same class as \code{x}, but
#'   with model fit criteria added for later usage.
#'
#' @details Functions \code{add_loo} and \code{add_waic} are aliases of
#'   \code{add_criterion} with fixed values for the \code{criterion} argument.
#'
#' @examples
#' \dontrun{
#' fit <- brm(count ~ Trt, data = epilepsy)
#' # add both LOO and WAIC at once
#' fit <- add_criterion(fit, c("loo", "waic"))
#' print(fit$criteria$loo)
#' print(fit$criteria$waic)
#' }
#'
#' @export
add_criterion <- function(x, ...) {
  UseMethod("add_criterion")
}

#' @rdname add_criterion
#' @export
add_criterion.brmsfit <- function(x, criterion, model_name = NULL,
                                  overwrite = FALSE, file = NULL,
                                  force_save = FALSE, ...) {
  if (!is.null(model_name)) {
    model_name <- as_one_character(model_name)
  } else {
    model_name <- deparse0(substitute(x))
  }
  criterion <- unique(as.character(criterion))
  if (any(criterion == "R2")) {
    # deprecated as of version 2.10.4
    warning2("Criterion 'R2' is deprecated. Please use 'bayes_R2' instead.")
    criterion[criterion == "R2"] <- "bayes_R2"
  }
  loo_options <- c("loo", "waic", "kfold", "loo_subsample")
  options <- c(loo_options, "bayes_R2", "loo_R2", "marglik")
  if (!length(criterion) || !all(criterion %in% options)) {
    stop2("Argument 'criterion' should be a subset of ",
          collapse_comma(options))
  }
  auto_save <- FALSE
  if (!is.null(file)) {
    file <- paste0(as_one_character(file), ".rds")
  } else {
    file <- x$file
    if (!is.null(file)) auto_save <- TRUE
  }
  force_save <- as_one_logical(force_save)
  overwrite <- as_one_logical(overwrite)
  if (overwrite) {
    # recompute all criteria
    new_criteria <- criterion
  } else {
    # only computed criteria not already stored
    new_criteria <- criterion[ulapply(x$criteria[criterion], is.null)]
  }
  # remove all criteria that are to be recomputed
  x$criteria[new_criteria] <- NULL
  args <- list(x, ...)
  for (fun in intersect(new_criteria, loo_options)) {
    args$model_names <- model_name
    x$criteria[[fun]] <- do_call(fun, args)
  }
  if ("bayes_R2" %in% new_criteria) {
    args$summary <- FALSE
    x$criteria$bayes_R2 <- do_call(bayes_R2, args)
  }
  if ("loo_R2" %in% new_criteria) {
    args$summary <- FALSE
    x$criteria$loo_R2 <- do_call(loo_R2, args)
  }
  if ("marglik" %in% new_criteria) {
    x$criteria$marglik <- do_call(bridge_sampler, args)
  }
  if (!is.null(file) && (force_save || length(new_criteria))) {
    if (auto_save) {
      message("Automatically saving the model object in '", file, "'")
    }
    x$file <- file
    saveRDS(x, file = file)
  }
  x
}

# extract a recomputed model fit criterion
get_criterion <- function(x, criterion) {
  stopifnot(is.brmsfit(x))
  criterion <- as_one_character(criterion)
  x$criteria[[criterion]]
}

# create a hash based on the response of a model
hash_response <- function(x, newdata = NULL, resp = NULL, ...) {
  require_package("digest")
  stopifnot(is.brmsfit(x))
  sdata <- standata(
    x, newdata = newdata, re_formula = NA, internal = TRUE,
    check_response = TRUE, only_response = TRUE
  )
  add_funs <- lsp("brms", what = "exports", pattern = "^resp_")
  regex <- c("Y", sub("^resp_", "", add_funs))
  regex <- outer(regex, escape_all(usc(resp)), FUN = paste0)
  regex <- paste0("(", as.vector(regex), ")", collapse = "|")
  regex <- paste0("^(", regex, ")(_|$)")
  out <- sdata[grepl(regex, names(sdata))]
  out <- as.matrix(as.data.frame(rmNULL(out)))
  out <- p(out, attr(sdata, "old_order"))
  # see issue #642
  attributes(out) <- NULL
  digest::sha1(x = out, ...)
}

# compare the response parts of multiple brmsfit objects
# @param models A list of brmsfit objects
# @param ... passed to hash_response
# @return TRUE if the response parts of all models match and FALSE otherwise
match_response <- function(models, ...) {
  if (length(models) <= 1L) {
    out <- TRUE
  } else {
    yhash <- lapply(models, hash_response, ...)
    yhash_check <- ulapply(yhash, is_equal, yhash[[1]])
    if (all(yhash_check)) {
      out <- TRUE
    } else {
      out <- FALSE
    }
  }
  out
}

# compare number of observations of multipe models
# @param models A list of brmsfit objects
# @param ... currently ignored
# @return TRUE if the number of rows match
match_nobs <- function(models, ...) {
  if (length(models) <= 1L) {
    out <- TRUE
  } else {
    nobs <- lapply(models, nobs)
    nobs_check <- ulapply(nobs, is_equal, nobs[[1]])
    if (all(nobs_check)) {
      out <- TRUE
    } else {
      out <- FALSE
    }
  }
  out
}

# validate models passed to loo and related methods
# @param models list of fitted model objects
# @param model_names names specified by the user
# @param sub_names names inferred by substitute()
validate_models <- function(models, model_names, sub_names) {
  stopifnot(is.list(models))
  model_names <- as.character(model_names)
  if (!length(model_names)) {
    model_names <- as.character(sub_names)
  }
  if (length(model_names) != length(models)) {
    stop2("Number of model names is not equal to the number of models.")
  }
  names(models) <- model_names
  for (i in seq_along(models)) {
    if (!is.brmsfit(models[[i]])) {
      stop2("Object '", names(models)[i], "' is not of class 'brmsfit'.")
    }
  }
  models
}

# recommend options if approximate loo fails for some observations
# @param moment_match has moment matching already been performed?
recommend_loo_options <- function(loo, k_threshold, moment_match = FALSE,
                                  model_name = "") {
  if (isTRUE(nzchar(model_name))) {
    model_name <- paste0(" in model '", model_name, "'")
  } else {
    model_name <- ""
  }
  n <- length(loo::pareto_k_ids(loo, threshold = k_threshold))
  if (!moment_match && n > 0) {
    warning2(
      "Found ", n, " observations with a pareto_k > ", k_threshold,
      model_name, ". It is recommended to set 'moment_match = TRUE' in order ",
      "to perform moment matching for problematic observations. "
    )
    out <- "loo_moment_match"
  } else if (n > 0 && n <= 10) {
    warning2(
      "Found ", n, " observations with a pareto_k > ", k_threshold,
      model_name, ". It is recommended to set 'reloo = TRUE' in order to ",
      "calculate the ELPD without the assumption that these observations " ,
      "are negligible. This will refit the model ", n, " times to compute ",
      "the ELPDs for the problematic observations directly."
    )
    out <- "reloo"
  } else if (n > 10) {
    warning2(
      "Found ", n, " observations with a pareto_k > ", k_threshold,
      model_name, ". With this many problematic observations, it may be more ",
      "appropriate to use 'kfold' with argument 'K = 10' to perform ",
      "10-fold cross-validation rather than LOO."
    )
    out <- "kfold"
  } else {
    out <- "loo"
  }
  invisible(out)
}

# helper function to compute relative efficiences
# @param x matrix of posterior draws
# @param fit a brmsfit object to extract metadata from
# @param allow_na allow NA values in the output?
# @return a numeric vector of length NCOL(x)
r_eff_helper <- function(x, chain_id, allow_na = TRUE, ...) {
  out <- loo::relative_eff(x, chain_id = chain_id, ...)
  if (!allow_na && anyNA(out)) {
    # avoid error in loo if some but not all r_effs are NA
    out <- rep(1, length(out))
    warning2(
      "Ignoring relative efficiencies as some were NA. ",
      "See argument 'r_eff' in ?loo::loo for more details."
    )
  }
  out
}

# wrapper around r_eff_helper to compute efficiency
# of likelihood draws based on log-likelihood draws
r_eff_log_lik <- function(x, ...) {
  UseMethod("r_eff_log_lik")
}

#' @export
r_eff_log_lik.matrix <- function(x, fit, allow_na = FALSE, ...) {
  if (is.brmsfit_multiple(fit)) {
    # due to stacking of chains from multiple models
    # efficiency computations will likely be incorrect
    # assume relative efficiency of 1 for now
    return(rep(1, ncol(x)))
  }
  chain_id <- get_chain_id(nrow(x), fit)
  r_eff_helper(exp(x), chain_id = chain_id, allow_na = allow_na, ...)
}

#' @export
r_eff_log_lik.function <- function(x, fit, draws, allow_na = FALSE, ...) {
  if (is.brmsfit_multiple(fit)) {
    # due to stacking of chains from multiple models
    # efficiency computations will likely be incorrect
    # assume relative efficiency of 1 for now
    return(rep(1, draws$nobs))
  }
  lik_fun <- function(data_i, draws, ...) {
    exp(x(data_i, draws, ...))
  }
  chain_id <- get_chain_id(draws$ndraws, fit)
  r_eff_helper(
    lik_fun, chain_id = chain_id, draws = draws,
    allow_na = allow_na, ...
  )
}

# get chain IDs per posterior draw
get_chain_id <- function(ndraws, fit) {
  if (ndraws != ndraws(fit)) {
    # don't know the chain IDs of a subset of draws
    chain_id <- rep(1L, ndraws)
  } else {
    nchains <- fit$fit@sim$chains
    chain_id <- rep(seq_len(nchains), each = ndraws / nchains)
  }
  chain_id
}

# print the output of a list of loo objects
#' @export
print.loolist <- function(x, digits = 1, ...) {
  model_names <- loo::find_model_names(x$loos)
  for (i in seq_along(x$loos)) {
    cat(paste0("Output of model '", model_names[i], "':\n"))
    print(x$loos[[i]], digits = digits, ...)
    cat("\n")
  }
  if (!is.null(x$diffs)) {
    cat("Model comparisons:\n")
    print(x$diffs, digits = digits, ...)
  }
  invisible(x)
}

# ---------- deprecated functions ----------
#' @rdname add_ic
#' @export
add_loo <- function(x, model_name = NULL, ...) {
  warning2("'add_loo' is deprecated. Please use 'add_criterion' instead.")
  if (!is.null(model_name)) {
    model_name <- as_one_character(model_name)
  } else {
    model_name <- deparse0(substitute(x))
  }
  add_criterion(x, criterion = "loo", model_name = model_name, ...)
}

#' @rdname add_ic
#' @export
add_waic <- function(x, model_name = NULL, ...) {
  warning2("'add_waic' is deprecated. Please use 'add_criterion' instead.")
  if (!is.null(model_name)) {
    model_name <- as_one_character(model_name)
  } else {
    model_name <- deparse0(substitute(x))
  }
  add_criterion(x, criterion = "waic", model_name = model_name, ...)
}

#' Add model fit criteria to model objects
#'
#' Deprecated aliases of \code{\link{add_criterion}}.
#'
#' @inheritParams add_criterion
#' @param ic,value Names of model fit criteria
#'   to compute. Currently supported are \code{"loo"},
#'   \code{"waic"}, \code{"kfold"}, \code{"R2"} (R-squared), and
#'   \code{"marglik"} (log marginal likelihood).
#'
#' @return An object of the same class as \code{x}, but
#'   with model fit criteria added for later usage.
#'   Previously computed criterion objects will be overwritten.
#'
#' @export
add_ic <- function(x, ...) {
  UseMethod("add_ic")
}

#' @rdname add_ic
#' @export
add_ic.brmsfit <- function(x, ic = "loo", model_name = NULL, ...) {
  warning2("'add_ic' is deprecated. Please use 'add_criterion' instead.")
  if (!is.null(model_name)) {
    model_name <- as_one_character(model_name)
  } else {
    model_name <- deparse0(substitute(x))
  }
  add_criterion(x, criterion = ic, model_name = model_name, ...)
}

#' @rdname add_ic
#' @export
'add_ic<-' <- function(x, ..., value) {
  add_ic(x, ic = value, ...)
}

#' Compare Information Criteria of Different Models
#'
#' Compare information criteria of different models fitted
#' with \code{\link{waic}} or \code{\link{loo}}.
#' Deprecated and will be removed in the future. Please use
#' \code{\link{loo_compare}} instead.
#'
#' @param ... At least two objects returned by
#'   \code{\link{waic}} or \code{\link{loo}}.
#'   Alternatively, \code{brmsfit} objects with information
#'   criteria precomputed via \code{\link{add_ic}}
#'   may be passed, as well.
#' @param x A \code{list} containing the same types of objects as
#'   can be passed via \code{...}.
#' @param ic The name of the information criterion to be extracted
#'   from \code{brmsfit} objects. Ignored if information
#'   criterion objects are only passed directly.
#'
#' @return An object of class \code{iclist}.
#'
#' @details See \code{\link{loo_compare}} for the recommended way
#'   of comparing models with the \pkg{loo} package.
#'
#' @seealso
#'   \code{\link{loo}},
#'   \code{\link{loo_compare}}
#'   \code{\link{add_criterion}}
#'
#' @examples
#' \dontrun{
#' # model with population-level effects only
#' fit1 <- brm(rating ~ treat + period + carry,
#'             data = inhaler)
#' waic1 <- waic(fit1)
#'
#' # model with an additional varying intercept for subjects
#' fit2 <- brm(rating ~ treat + period + carry + (1|subject),
#'             data = inhaler)
#' waic2 <- waic(fit2)
#'
#' # compare both models
#' compare_ic(waic1, waic2)
#' }
#'
#' @export
compare_ic <- function(..., x = NULL, ic = c("loo", "waic", "kfold")) {
  # will be removed in brms 3.0
  warning2(
    "'compare_ic' is deprecated and will be removed ",
    "in the future. Please use 'loo_compare' instead."
  )
  ic <- match.arg(ic)
  if (!(is.null(x) || is.list(x))) {
    stop2("Argument 'x' should be a list.")
  }
  x$ic_diffs__ <- NULL
  x <- c(list(...), x)
  for (i in seq_along(x)) {
    # extract precomputed values from brmsfit objects
    if (is.brmsfit(x[[i]]) && !is.null(x[[i]][[ic]])) {
      x[[i]] <- x[[i]][[ic]]
    }
  }
  if (!all(sapply(x, inherits, "loo"))) {
    stop2("All inputs should have class 'loo' ",
          "or contain precomputed 'loo' objects.")
  }
  if (length(x) < 2L) {
    stop2("Expecting at least two objects.")
  }
  ics <- unname(sapply(x, function(y) rownames(y$estimates)[3]))
  if (!all(ics %in% ics[1])) {
    stop2("All inputs should be from the same criterion.")
  }
  yhash <- lapply(x, attr, which = "yhash")
  yhash_check <- ulapply(yhash, is_equal, yhash[[1]])
  if (!all(yhash_check)) {
    warning2(
      "Model comparisons are likely invalid as the response ",
      "values of at least two models do not match."
    )
  }
  names(x) <- loo::find_model_names(x)
  n_models <- length(x)
  ic_diffs <- matrix(0, nrow = n_models * (n_models - 1) / 2, ncol = 2)
  rnames <- rep("", nrow(ic_diffs))
  # pairwise comparision to get differences in ICs and their SEs
  n <- 1
  for (i in seq_len(n_models - 1)) {
    for (j in (i + 1):n_models) {
      tmp <- SW(loo::compare(x[[j]], x[[i]]))
      ic_diffs[n, ] <- c(-2 * tmp[["elpd_diff"]], 2 * tmp[["se"]])
      rnames[n] <- paste(names(x)[i], "-", names(x)[j])
      n <- n + 1
    }
  }
  rownames(ic_diffs) <- rnames
  colnames(ic_diffs) <- c(toupper(ics[1]), "SE")
  x$ic_diffs__ <- ic_diffs
  class(x) <- "iclist"
  x
}

# print the output of LOO and WAIC with multiple models
# deprecated as of brms > 2.5.0 and will be removed in brms 3.0
#' @export
print.iclist <- function(x, digits = 2, ...) {
  m <- x
  m$ic_diffs__ <- NULL
  if (length(m)) {
    ic <- rownames(m[[1]]$estimates)[3]
    mat <- matrix(0, nrow = length(m), ncol = 2)
    dimnames(mat) <- list(names(m), c(toupper(ic), "SE"))
    for (i in seq_along(m)) {
      mat[i, ] <- m[[i]]$estimates[3, ]
    }
  } else {
    mat <- ic <- NULL
  }
  ic_diffs <- x$ic_diffs__
  if (is.matrix(attr(x, "compare"))) {
    # deprecated as of brms 1.4.0
    ic_diffs <- attr(x, "compare")
  }
  if (is.matrix(ic_diffs)) {
    # models were compared using the compare_ic function
    mat <- rbind(mat, ic_diffs)
  }
  print(round(mat, digits = digits), na.print = "")
  invisible(x)
}
