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
