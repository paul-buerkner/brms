#' @export
LOO <- function(x, ...) {
  UseMethod("LOO")
}

#' @export
WAIC <- function(x, ...) {
  UseMethod("WAIC")
}

#' @export
loo_R2 <- function(object, ...) {
  # temporary generic until available in loo
  UseMethod("loo_R2")
}

compute_loos <- function(
  models, criterion = c("loo", "waic", "psis", "psislw", "kfold"),
  use_stored = TRUE, compare = TRUE, ...
) {
  # helper function used to create (lists of) 'loo' objects
  # Args:
  #   models: list of brmsfit objects
  #   criterion: name of the criterion to compute
  #   use_stored: use precomputed criterion objects if possible?
  #   ...: more arguments passed to compute_loo
  # Returns:
  #   If length(models) > 1 an object of class 'loolist'
  #   If length(models) == 1 an object of class 'loo'
  criterion <- match.arg(criterion)
  args <- nlist(criterion, ...)
  if (length(models) > 1L) {
    warning2(
      "Passing multiple brmsfit objects to 'loo' and related methods is ",
      "deprecated. Please see ?loo.brmsfit for the recommended workflow."
    )
    if (!match_nobs(models)) {
      stop2("Models have different number of observations.")
    }
    if (length(use_stored) == 1L) {
      use_stored <- rep(use_stored, length(models))
    }
    out <- named_list(names(models))
    for (i in seq_along(models)) {
      args$x <- models[[i]]
      args$model_name <- names(models)[i]
      args$use_stored <- use_stored[i]
      out[[i]] <- run(compute_loo, args) 
    }
    compare <- as_one_logical(compare)
    if (compare) {
      diffs__ <- loo_compare(out)
      out <- SW(compare_ic(x = out))
      out$diffs__ <- diffs__
    }
    class(out) <- "loolist"
  } else {
    args$x <- models[[1]]
    args$model_name <- names(models)
    args$use_stored <- use_stored
    out <- run(compute_loo, args) 
  }
  out
}

compute_loo <- function(x, criterion = c("loo", "waic", "psis", "kfold"),
                        reloo = FALSE, k_threshold = 0.7, reloo_args = list(),
                        pointwise = FALSE, newdata = NULL, model_name = "", 
                        use_stored = TRUE, ...) {
  # compute information criteria using the 'loo' package
  # Args:
  #   x: an object of class brmsfit
  #   criterion: the criterion to be computed
  #   model_name: original variable name of object 'x'
  #   use_stored: use precomputed criterion objects if possible?
  #   newdata: optional data.frame of new data
  #   reloo: call 'reloo' after computing 'loo'?
  #   reloo_args: list of arguments passed to 'reloo'
  #   pointwise: compute log-likelihood point-by-point?
  #   ...: passed to other post-processing methods
  # Returns:
  #   an object of class 'loo'
  criterion <- match.arg(criterion)
  model_name <- as_one_character(model_name)
  use_stored <- as_one_logical(use_stored)
  if (use_stored && is.loo(x[[criterion]])) {
    # extract the stored criterion
    out <- x[[criterion]]
  } else {
    # compute the criterion
    if (criterion == "kfold") {
      kfold_args <- nlist(x, newdata, ...)
      out <- run(kfold_internal, kfold_args)
    } else {
      contains_samples(x)
      pointwise <- as_one_logical(pointwise)
      loo_args <- list(...)
      ll_args <- nlist(object = x, newdata, pointwise, ...)
      loo_args$x <- run(log_lik, ll_args)
      if (pointwise) {
        loo_args$draws <- attr(loo_args$x, "draws")
        loo_args$data <- attr(loo_args$x, "data")
      }
      if (criterion == "psis") {
        if (pointwise) {
          stop2("Cannot use pointwise evaluation for 'psis'.")
        }
        loo_args$log_ratios <- -loo_args$x
        loo_args$x <- NULL
      }
      out <- SW(run(criterion, loo_args, pkg = "loo"))
    }
    attr(out, "yhash") <- hash_response(x, newdata = newdata)
  }
  attr(out, "model_name") <- model_name
  if (criterion == "loo") {
    if (reloo) {
      c(reloo_args) <- nlist(
        x = out, fit = x, newdata, k_threshold, check = FALSE, ...
      )
      out <- run("reloo", reloo_args)
    } else {
      n_bad_obs <- length(loo::pareto_k_ids(out, threshold = k_threshold))
      recommend_loo_options(n_bad_obs, k_threshold, model_name) 
    }
  }
  out
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
loo_compare.brmsfit <- function(
  x, ..., criterion = c("loo", "waic", "kfold"),
  model_names = NULL
) {
  criterion <- match.arg(criterion)
  models <- split_dots(x, ..., model_names = model_names, other = FALSE)
  loos <- named_list(names(models))
  for (i in seq_along(models)) {
    if (is.null(models[[i]][[criterion]])) {
      stop2(
        "Model '", names(models)[i], "' does not contain a precomputed '",
        criterion, "' criterion. See ?loo_compare.brmsfit for help."
      )
    }
    loos[[i]] <- models[[i]][[criterion]]
  }
  loo_compare(loos)
}

#' Add model fit criteria to model objects
#' 
#' @param x An \R object typically of class \code{brmsfit}.
#' @param criterion Names of model fit criteria
#'   to compute. Currently supported are \code{"loo"}, 
#'   \code{"waic"}, \code{"kfold"}, \code{"R2"} (R-squared), and 
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
#' print(fit$loo)
#' print(fit$waic)
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
    model_name <- deparse_combine(substitute(x)) 
  }
  criterion <- unique(as.character(criterion))
  options <- c("loo", "waic", "kfold", "R2", "marglik")
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
    # remove previously stored criterion objects
    x[criterion] <- list(NULL)
  }
  new_criteria <- criterion[ulapply(x[criterion], is.null)]
  args <- list(x, ...)
  for (fun in intersect(criterion, c("loo", "waic", "kfold"))) {
    args$model_names <- model_name
    x[[fun]] <- run(fun, args)
  }
  if ("R2" %in% criterion) {
    args$summary <- FALSE
    x$R2 <- run(bayes_R2, args)
  }
  if ("marglik" %in% criterion) {
    x$marglik <- run(bridge_sampler, args)
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

#' @rdname add_criterion
#' @export
add_loo <- function(x, model_name = NULL, ...) {
  if (!is.null(model_name)) {
    model_name <- as_one_character(model_name)
  } else {
    model_name <- deparse_combine(substitute(x)) 
  }
  add_criterion(x, criterion = "loo", model_name = model_name, ...)
}

#' @rdname add_criterion
#' @export
add_waic <- function(x, model_name = NULL, ...) {
  if (!is.null(model_name)) {
    model_name <- as_one_character(model_name)
  } else {
    model_name <- deparse_combine(substitute(x)) 
  }
  add_criterion(x, criterion = "waic", model_name = model_name, ...)
}

hash_response <- function(x, newdata = NULL, ...) {
  # create a hash based on the response of a model
  require_package("digest")
  stopifnot(is.brmsfit(x))
  sdata <- standata(
    x, newdata = newdata, re_formula = NA, internal = TRUE, 
    check_response = TRUE, only_response = TRUE
  )
  add_funs <- lsp("brms", what = "exports", pattern = "^resp_")
  regex <- c("Y", sub("^resp_", "", add_funs))
  regex <- paste0("(", regex, ")", collapse = "|")
  regex <- paste0("^(", regex, ")(_|$)")
  out <- sdata[grepl(regex, names(sdata))]
  out <- as.matrix(as.data.frame(rmNULL(out)))
  out <- p(out, attr(sdata, "old_order"))
  digest::sha1(x = out, ...)
}

match_response <- function(models, ...) {
  # compare the response parts of multiple brmsfit objects
  # Args:
  #   models: A list of brmsfit objects
  #   ...: passed to hash_response
  # Returns:
  #   TRUE if the response parts of all models match and FALSE otherwise
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

match_nobs <- function(models, ...) {
  # compare number of observations of multipe models
  # Args:
  #   models: A list of brmsfit objects
  # Returns:
  #   TRUE if the number of rows match
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

validate_models <- function(models, model_names, sub_names) {
  # validate models passed to loo and related methods
  # Args:
  #   models: list of fitted model objects
  #   model_names: names specified by the user
  #   sub_names: names inferred by substitute()
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

#' Compute exact cross-validation for problematic observations
#' 
#' Compute exact cross-validation for problematic observations
#' for which approximate leave-one-out cross-validation may
#' return incorrect results.
#' 
#' @inheritParams predict.brmsfit
#' @param x An \R object of class \code{loo}.
#' @param fit An \R object of class \code{brmsfit}.
#' @param k_threshold The threshold at which pareto \eqn{k} 
#'   estimates are treated as problematic. Defaults to \code{0.7}. 
#'   See \code{\link[loo:pareto_k_ids]{pareto_k_ids}}
#'   for more details.
#' @param check Logical; If \code{TRUE} (the default), some checks
#'   check are performed if the \code{loo} object was generated
#'   from the \code{brmsfit} object passed to argument \code{fit}.
#' @param ... Further arguments passed to 
#'   \code{\link{update.brmsfit}} and \code{\link{log_lik.brmsfit}}.
#'   
#' @return An object of the class \code{loo}.
#' 
#' @details 
#' Warnings about Pareto \eqn{k} estimates indicate observations
#' for which the approximation to LOO is problematic (this is described in
#' detail in Vehtari, Gelman, and Gabry (2017) and the 
#' \pkg{\link[loo:loo-package]{loo}} package documentation).
#' If there are \eqn{J} observations with \eqn{k} estimates above
#' \code{k_threshold}, then \code{reloo} will refit the original model 
#' \eqn{J} times, each time leaving out one of the \eqn{J} 
#' problematic observations. The pointwise contributions of these observations
#' to the total ELPD are then computed directly and substituted for the
#' previous estimates from these \eqn{J} observations that are stored in the
#' original \code{loo} object.
#' 
#' @seealso \code{\link{loo}}, \code{\link{kfold}}
#' 
#' @examples 
#' \dontrun{
#' fit1 <- brm(count ~ log_Age_c + log_Base4_c * Trt + (1|patient),
#'             data = epilepsy, family = poisson())
#' # throws warning about some pareto k estimates being too high
#' (loo1 <- loo(fit1))
#' (reloo1 <- reloo(loo1, fit1, chains = 1))
#' }
#' 
#' @export
reloo <- function(x, fit, k_threshold = 0.7, newdata = NULL, 
                  resp = NULL, check = TRUE, ...) {
  stopifnot(is.loo(x), is.brmsfit(fit))
  if (is.null(newdata)) {
    mf <- model.frame(fit) 
  } else {
    mf <- as.data.frame(newdata)
  }
  if (NROW(mf) != NROW(x$pointwise)) {
    stop2("Number of observations in 'x' and 'fit' do not match.")
  }
  if (check) {
    yhash_loo <- attr(x, "yhash")
    yhash_fit <- hash_response(fit, newdata = newdata)
    if (!is_equal(yhash_loo, yhash_fit)) {
      stop2(
        "Response values used in 'x' and 'fit' do not match. ",
        "If this is a false positive, please set 'check' to FALSE."
      )
    }
  }
  if (is.null(x$diagnostics$pareto_k)) {
    stop2("No Pareto k estimates found in the 'loo' object.")
  }
  obs <- loo::pareto_k_ids(x, k_threshold)
  J <- length(obs)
  if (J == 0L) {
    message(
      "No problematic observations found. ",
      "Returning the original 'loo' object."
    )
    return(x)
  }
  lls <- vector("list", J)
  message(
    J, " problematic observation(s) found.", 
    "\nThe model will be refit ", J, " times."
  )
  for (j in seq_len(J)) {
    message(
      "\nFitting model ", j, " out of ", J,
      " (leaving out observation ", obs[j], ")"
    )
    omitted <- obs[j]
    mf_omitted <- mf[-omitted, , drop = FALSE]
    fit_j <- subset_autocor(fit, -omitted)
    fit_j <- SW(update(fit_j, newdata = mf_omitted, refresh = 0, ...))
    fit_j <- subset_autocor(fit_j, omitted, autocor = x$autocor)
    lls[[j]] <- log_lik(
      fit_j, newdata = mf[omitted, , drop = FALSE],
      allow_new_levels = TRUE, resp = resp, ...
    )
  }
  # most of the following code is taken from rstanarm:::reloo
  # compute elpd_{loo,j} for each of the held out observations
  elpd_loo <- ulapply(lls, log_mean_exp)
  # compute \hat{lpd}_j for each of the held out observations (using log-lik
  # matrix from full posterior, not the leave-one-out posteriors)
  fit <- subset_autocor(fit, obs)
  ll_x <- log_lik(fit, newdata = mf[obs, , drop = FALSE])
  hat_lpd <- apply(ll_x, 2, log_mean_exp)
  # compute effective number of parameters
  p_loo <- hat_lpd - elpd_loo
  # replace parts of the loo object with these computed quantities
  sel <- c("elpd_loo", "p_loo", "looic")
  x$pointwise[obs, sel] <- cbind(elpd_loo, p_loo, -2 * elpd_loo)
  new_pw <- x$pointwise[, sel, drop = FALSE]
  x$estimates[, 1] <- colSums(new_pw)
  x$estimates[, 2] <- sqrt(nrow(x$pointwise) * apply(new_pw, 2, var))
  # what should we do about pareto-k? for now setting them to 0
  x$diagnostics$pareto_k[obs] <- 0
  x
}

kfold_internal <- function(x, K = 10, Ksub = NULL, folds = NULL, 
                           group = NULL, newdata = NULL, resp = NULL,
                           save_fits = FALSE, ...) {
  # helper function to perform k-fold cross-validation
  # Args: see kfold.brmsfit
  stopifnot(is.brmsfit(x))
  if (is.null(newdata)) {
    mf <- model.frame(x) 
  } else {
    mf <- as.data.frame(newdata)
  }
  N <- nrow(mf)
  # validate argument 'group'
  if (!is.null(group)) {
    valid_groups <- get_cat_vars(x)
    if (length(group) != 1L || !group %in% valid_groups) {
      stop2("Group '", group, "' is not a valid grouping factor. ",
            "Valid groups are: \n", collapse_comma(valid_groups))
    }
    gvar <- factor(get(group, mf))
  }
  # validate argument 'folds'
  if (is.null(folds)) {
    if (is.null(group)) {
      fold_type <- "random"
      folds <- loo::kfold_split_random(K, N)
    } else {
      fold_type <- "group"
      folds <- as.numeric(gvar)
      K <- length(levels(gvar))
      message("Setting 'K' to the number of levels of '", group, "' (", K, ")")
    }
  } else if (is.character(folds) && length(folds) == 1L) {
    opts <- c("stratified", "balanced", "loo")
    fold_type <- match.arg(folds, opts)
    if (fold_type == "loo") {
      folds <- seq_len(N)
      K <- N
      message("Setting 'K' to the number of observations (", K, ")")
    } else if (fold_type == "stratified") {
      if (is.null(group)) {
        stop2("Argument 'group' is required for stratified folds.")
      }
      folds <- loo::kfold_split_stratified(K, gvar)
    } else if (fold_type == "balanced") {
      if (is.null(group)) {
        stop2("Argument 'group' is required for balanced folds.")
      }
      folds <- loo::kfold_split_balanced(K, gvar)
    }
  } else {
    fold_type <- "custom"
    folds <- as.numeric(factor(folds))
    if (length(folds) != N) {
      stop2("If 'folds' is a vector, it must be of length N.")
    }
    K <- max(folds)
    message("Setting 'K' to the number of folds (", K, ")")
  }
  # validate argument 'Ksub'
  if (is.null(Ksub)) {
    Ksub <- seq_len(K)
  } else {
    # see issue #441 for reasons to check for arrays
    is_array_Ksub <- is.array(Ksub)
    Ksub <- as.integer(Ksub)
    if (any(Ksub <= 0 | Ksub > K)) {
      stop2("'Ksub' must contain positive integers not larger than 'K'.")
    }
    if (length(Ksub) == 1L && !is_array_Ksub) {
      Ksub <- sample(seq_len(K), Ksub)
    } else {
      Ksub <- unique(Ksub)
    }
    Ksub <- sort(Ksub)
  }
  
  if (save_fits) {
    fits <- array(list(), dim = c(length(Ksub), 3))
    dimnames(fits) <- list(NULL, c("fit", "omitted", "predicted"))
  }
  lppds <- obs_order <- vector("list", length(Ksub))
  for (k in Ksub) {
    message("Fitting model ", k, " out of ", K)
    ks <- match(k, Ksub)
    if (fold_type == "loo" && !is.null(group)) {
      omitted <- which(folds == folds[k])
      predicted <- k
    } else {
      omitted <- predicted <- which(folds == k)
    }
    mf_omitted <- mf[-omitted, , drop = FALSE]
    fit_k <- subset_autocor(x, -omitted)
    fit_k <- SW(update(fit_k, newdata = mf_omitted, refresh = 0, ...))
    if (save_fits) {
      fits[ks, ] <- nlist(fit = fit_k, omitted, predicted) 
    }
    fit_k <- subset_autocor(fit_k, predicted, autocor = x$autocor)
    obs_order[[ks]] <- predicted
    lppds[[ks]] <- log_lik(
      fit_k, newdata = mf[predicted, , drop = FALSE], 
      allow_new_levels = TRUE, resp = resp
    )
  }
  elpds <- ulapply(lppds, function(x) apply(x, 2, log_mean_exp))
  # make sure elpds are put back in the right order
  elpds <- elpds[order(unlist(obs_order))]
  elpd_kfold <- sum(elpds)
  se_elpd_kfold <- sqrt(length(elpds) * var(elpds))
  rnames <- c("elpd_kfold", "p_kfold", "kfoldic")
  cnames <- c("Estimate", "SE")
  estimates <- matrix(nrow = 3, ncol = 2, dimnames = list(rnames, cnames))
  estimates[1, ] <- c(elpd_kfold, se_elpd_kfold)
  estimates[3, ] <- c(-2 * elpd_kfold, 2 * se_elpd_kfold)
  out <- nlist(estimates, pointwise = cbind(elpd_kfold = elpds))
  atts <- nlist(K, Ksub, group, folds, fold_type)
  attributes(out)[names(atts)] <- atts
  if (save_fits) {
    out$fits <- fits
    out$data <- mf
  }
  structure(out, class = c("kfold", "loo"))
}

#' Predictions from K-Fold Cross-Validation
#' 
#' Compute and evaluate predictions after performing K-fold 
#' cross-validation via \code{\link{kfold}}. 
#' 
#' @param x Object of class \code{'kfold'} computed by \code{\link{kfold}}.
#'   For \code{kfold_predict} to work, the fitted model objects need to have
#'   been stored via argument \code{save_fits} of \code{\link{kfold}}.
#' @param method The method used to make predictions. Either \code{"predict"}
#'   or \code{"fitted"}. See \code{\link{predict.brmsfit}} for details.
#' @inheritParams predict.brmsfit
#' 
#' @return A \code{list} with two slots named \code{'y'} and \code{'yrep'}.
#'   Slot \code{y} contains the vector of observed responses.
#'   Slot \code{yrep} contains the matrix of predicted responses,
#'   with rows being posterior draws and columns being observations.
#'   
#' @seealso \code{\link{kfold}}
#'   
#' @examples 
#' \dontrun{
#' fit <- brm(count ~ log_Base4_c * Trt + (1|patient),
#'            data = epilepsy, family = poisson())
#'             
#' # perform k-fold cross validation
#' (kf <- kfold(fit, save_fits = TRUE, chains = 1))
#' 
#' # define a loss function
#' rmse <- function(y, yrep) {
#'   yrep_mean <- colMeans(yrep)
#'   sqrt(mean((yrep_mean - y)^2))
#' }
#' 
#' # predict responses and evaluate the loss
#' kfp <- kfold_predict(kf)
#' rmse(y = kfp$y, yrep = kfp$yrep)
#' }
#'   
#' @export
kfold_predict <- function(x, method = c("predict", "fitted"), 
                          resp = NULL, ...) {
  if (!inherits(x, "kfold")) {
    stop2("'x' must be a 'kfold' object.")
  }
  if (!all(c("fits", "data") %in% names(x))) {
    stop2(
      "Slots 'fits' and 'data' are required. ", 
      "Please run kfold with 'save_fits = TRUE'."
    )
  }
  method <- get(match.arg(method), mode = "function")
  resp <- validate_resp(resp, x$fits[[1, "fit"]], multiple = FALSE)
  all_predicted <- as.character(sort(unlist(x$fits[, "predicted"])))
  npredicted <- length(all_predicted)
  nsamples <- nsamples(x$fits[[1, "fit"]])
  y <- rep(NA, npredicted)
  yrep <- matrix(NA, nrow = nsamples, ncol = npredicted)
  names(y) <- colnames(yrep) <- all_predicted
  for (k in seq_rows(x$fits)) {
    fit_k <- x$fits[[k, "fit"]]
    predicted_k <- x$fits[[k, "predicted"]]
    obs_names <- as.character(predicted_k)
    newdata <- x$data[predicted_k, , drop = FALSE]
    y[obs_names] <- get_y(fit_k, resp, newdata = newdata, ...)
    yrep[, obs_names] <- method(
      fit_k, newdata = newdata, resp = resp, 
      allow_new_levels = TRUE, summary = FALSE, ...
    )
  }
  nlist(y, yrep)
}

recommend_loo_options <- function(n, k_threshold, model_name = "") {
  model_name <- if (isTRUE(nzchar(model_name))) {
    paste0(" in model '", model_name, "'")
  }
  if (n > 0 && n <= 10) {
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

r_eff_helper <- function(log_lik, fit) {
  # helper function to compute relative efficiences
  stopifnot(is.matrix(log_lik), is.brmsfit(fit))
  chains <- fit$fit@sim$chains
  chain_id <- rep(seq_len(chains), each = nrow(log_lik) / chains)
  loo::relative_eff(log_lik, chain_id = chain_id)
}

#' @export
print.loolist <- function(x, digits = 1, ...) {
  # print the output of LOO and WAIC with multiple models
  loos <- x
  loos$ic_diffs__ <- loos$diffs__ <- NULL
  for (i in seq_along(loos)) {
    name <- attributes(loos[[i]])[["model_name"]]
    if (is.null(name)) {
      name <- loos[[i]][["model_name"]]
    }
    if (is.null(name)) {
      name <- paste0("model", i)
    }
    cat(paste0("Output of model '", name, "':\n"))
    print(loos[[i]], digits = digits, ...)
    cat("\n")
  }
  if (!is.null(x$diffs__)) {
    cat("Model comparisons:\n")
    print(x$diffs__, digits = digits, ...)
  }
  invisible(x)
}

# ---------- deprecated functions ----------
#' Add model fit criteria to model objects
#' 
#' Deprecated alias of \code{\link{add_criterion}}.
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
    model_name <- deparse_combine(substitute(x)) 
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

#' @export
print.iclist <- function(x, digits = 2, ...) {
  # deprecated as of brms > 2.5.0 and will be removed in brms 3.0
  # print the output of LOO and WAIC with multiple models
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
