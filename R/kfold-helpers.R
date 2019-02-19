#' K-Fold Cross-Validation
#' 
#' Perform exact K-fold cross-validation by refitting the model \eqn{K}
#' times each leaving out one-\eqn{K}th of the original data.
#' 
#' @inheritParams loo.brmsfit
#' @param K The number of subsets of equal (if possible) size
#'   into which the data will be partitioned for performing
#'   \eqn{K}-fold cross-validation. The model is refit \code{K} times, each time
#'   leaving out one of the \code{K} subsets. If \code{K} is equal to the total
#'   number of observations in the data then \eqn{K}-fold cross-validation is
#'   equivalent to exact leave-one-out cross-validation.
#' @param Ksub Optional number of subsets (of those subsets defined by \code{K}) 
#'   to be evaluated. If \code{NULL} (the default), \eqn{K}-fold cross-validation 
#'   will be performed on all subsets. If \code{Ksub} is a single integer, 
#'   \code{Ksub} subsets (out of all \code{K}) subsets will be randomly chosen.
#'   If \code{Ksub} consists of multiple integers or a one-dimensional array
#'   (created via \code{as.array}) potentially of length one, the corresponding 
#'   subsets will be used. This argument is primarily useful, if evaluation of 
#'   all subsets is infeasible for some reason.
#' @param folds Determines how the subsets are being constructed.
#'   Possible values are \code{NULL} (the default), \code{"stratified"},
#'   \code{"balanced"}, or \code{"loo"}. May also be a vector of length
#'   equal to the number of observations in the data. Alters the way
#'   \code{group} is handled. More information is provided in the 'Details'
#'   section.
#' @param group Optional name of a grouping variable or factor in the model.
#'   What exactly is done with this variable depends on argument \code{folds}.
#'   More information is provided in the 'Details' section.
#' @param exact_loo Deprecated! Please use \code{folds = "loo"} instead.
#' @param save_fits If \code{TRUE}, components \code{fits} and \code{data} are 
#'   added to the returned object. \code{fits} stores the cross-validated 
#'   \code{brmsfit} objects and the indices of the omitted and predicted 
#'   observations for each fold. \code{data} returns the full data frame
#'   used in the cross-validation. Defaults to \code{FALSE}.
#'   
#' @return \code{kfold} returns an object that has a similar structure as the 
#'   objects returned by the \code{loo} and \code{waic} methods.
#'    
#' @details The \code{kfold} function performs exact \eqn{K}-fold
#'   cross-validation. First the data are partitioned into \eqn{K} folds 
#'   (i.e. subsets) of equal (or as close to equal as possible) size by default. 
#'   Then the model is refit \eqn{K} times, each time leaving out one of the 
#'   \code{K} subsets. If \eqn{K} is equal to the total number of observations 
#'   in the data then \eqn{K}-fold cross-validation is equivalent to exact 
#'   leave-one-out cross-validation (to which \code{loo} is an efficient 
#'   approximation). The \code{compare_ic} function is also compatible with 
#'   the objects returned by \code{kfold}.
#'   
#'   The subsets can be constructed in multiple different ways: 
#'   \itemize{
#'   \item If both \code{folds} and \code{group} are \code{NULL}, the subsets 
#'   are randomly chosen so that they have equal (or as close to equal as 
#'   possible) size. 
#'   \item If \code{folds} is \code{NULL} but \code{group} is specified, the 
#'   data is split up into subsets, each time omitting all observations of one 
#'   of the factor levels, while ignoring argument \code{K}. 
#'   \item If \code{folds = "stratified"} the subsets are stratified after 
#'   \code{group} using \code{\link[loo:kfold-helpers]{loo::kfold_split_stratified}}.
#'   \item If \code{folds = "balanced"} the subsets are balanced by
#'   \code{group} using \code{\link[loo:kfold-helpers]{loo::kfold_split_balanced}}.
#'   \item If \code{folds = "loo"} exact leave-one-out cross-validation
#'   will be performed and \code{K} will be ignored. Further, if \code{group}
#'   is specified, all observations corresponding to the factor level of the 
#'   currently predicted single value are omitted. Thus, in this case, the 
#'   predicted values are only a subset of the omitted ones.
#'   \item If \code{folds} is a numeric vector, it must contain one element per 
#'   observation in the data. Each element of the vector is an integer in 
#'   \code{1:K} indicating to which of the \code{K} folds the corresponding 
#'   observation belongs. There are some convenience functions available in 
#'   the \pkg{loo} package that create integer vectors to use for this purpose 
#'   (see the Examples section below and also the 
#'   \link[loo:kfold-helpers]{kfold-helpers} page).
#' }
#'   
#' @seealso \code{\link{loo}}, \code{\link{reloo}}, \code{\link{kfold_predict}}
#'   
#' @examples 
#' \dontrun{
#' fit1 <- brm(count ~ zAge + zBase * Trt + 
#'               (1|patient) + (1|obs),
#'            data = epilepsy, family = poisson())
#' # throws warning about some pareto k estimates being too high
#' (loo1 <- loo(fit1))
#' # perform 10-fold cross validation
#' (kfold1 <- kfold(fit1, chains = 1)
#' }   
#'  
#' @export
kfold <- function(x, ...) {
  UseMethod("kfold")
}

kfold_internal <- function(x, K = 10, Ksub = NULL, folds = NULL, 
                           group = NULL, newdata = NULL, resp = NULL,
                           save_fits = FALSE, ...) {
  # most of the code is taken from rstanarm::kfold
  # Args:
  #   group: character string of length one naming 
  #     a variable to group excluded chunks
  stopifnot(is.brmsfit(x))
  if (is.null(newdata)) {
    mf <- model.frame(x) 
  } else {
    mf <- as.data.frame(newdata)
  }
  mf <- rm_attr(mf, c("terms", "brmsframe"))
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
  out <- nlist(
    estimates, pointwise = cbind(elpd_kfold = elpds),
    model_name = deparse(substitute(x)), K, Ksub, 
    group, folds, fold_type
  )
  if (save_fits) {
    out$fits <- fits 
    out$data <- mf
  }
  structure(out, class = c("kfold", "loo"))
}

#' Predictions of K-Fold Cross-Validation
#' 
#' Compute and evaluate predictions after performing from K-fold 
#' cross-validation via \code{\link{kfold}}. 
#' 
#' @param x Object of class \code{'kfold'} computed by \code{\link{kfold}}.
#'   For \code{kfold_predict} to work, the fitted model objects need to be
#'   stored in \code{x} via argument \code{save_fits} of \code{\link{kfold}}.
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
#' fit <- brm(count ~ zBase * Trt + (1|patient),
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
  stopifnot(is.kfold(x))
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

is.kfold <- function(x) {
  inherits(x, "kfold")
}

