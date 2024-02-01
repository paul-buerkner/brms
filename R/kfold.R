#' K-Fold Cross-Validation
#'
#' Perform exact K-fold cross-validation by refitting the model \eqn{K}
#' times each leaving out one-\eqn{K}th of the original data.
#' Folds can be run in parallel using the \pkg{future} package.
#'
#' @aliases kfold
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
#'   \code{"grouped"}, or \code{"loo"}. May also be a vector of length
#'   equal to the number of observations in the data. Alters the way
#'   \code{group} is handled. More information is provided in the 'Details'
#'   section.
#' @param group Optional name of a grouping variable or factor in the model.
#'   What exactly is done with this variable depends on argument \code{folds}.
#'   More information is provided in the 'Details' section.
#' @param joint Logical; If \code{TRUE}, the joint log likelihoods per
#'   fold are used for ELPD computation instead of the pointwise
#'   (per observations) likelihoods. The latter is the default
#'   for comparability of \code{kfold} with \code{loo}.
#' @param save_fits If \code{TRUE}, a component \code{fits} is added to
#'   the returned object to store the cross-validated \code{brmsfit}
#'   objects and the indices of the omitted observations for each fold.
#'   Defaults to \code{FALSE}.
#' @param recompile Logical, indicating whether the Stan model should be
#'   recompiled. This may be necessary if you are running \code{reloo} on
#'   another machine than the one used to fit the model.
#' @param future_args A list of further arguments passed to
#'   \code{\link[future:future]{future}} for additional control over parallel
#'   execution if activated.
#' @param ... Further arguments passed to \code{\link{brm}}.
#'
#' @return \code{kfold} returns an object that has a similar structure as the
#'   objects returned by the \code{loo} and \code{waic} methods and
#'   can be used with the same post-processing functions.
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
#'   \item If \code{folds = "grouped"} the subsets are split by
#'   \code{group} using \code{\link[loo:kfold-helpers]{loo::kfold_split_grouped}}.
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
#'   }
#'
#'   When running \code{kfold} on a \code{brmsfit} created with the
#'   \pkg{cmdstanr} backend in a different \R session, several recompilations
#'   will be triggered because by default, \pkg{cmdstanr} writes the model
#'   executable to a temporary directory. To avoid that, set option
#'   \code{"cmdstanr_write_stan_file_dir"} to a nontemporary path of your choice
#'   before creating the original \code{brmsfit} (see section 'Examples' below).
#'
#' @examples
#' \dontrun{
#' fit1 <- brm(count ~ zAge + zBase * Trt + (1|patient) + (1|obs),
#'            data = epilepsy, family = poisson())
#' # throws warning about some pareto k estimates being too high
#' (loo1 <- loo(fit1))
#' # perform 10-fold cross validation
#' (kfold1 <- kfold(fit1, chains = 1))
#'
#' # use joint likelihoods per fold for ELPD evaluation
#' kfold(fit1, chains = 1, joint = TRUE)
#'
#' # use the future package for parallelization of models
#' # that is to fit models belonging to different folds in parallel
#' library(future)
#' plan(multisession, workers = 4)
#' kfold(fit1, chains = 1)
#' plan(sequential)
#'
#' ## to avoid recompilations when running kfold() on a 'cmdstanr'-backend fit
#' ## in a fresh R session, set option 'cmdstanr_write_stan_file_dir' before
#' ## creating the initial 'brmsfit'
#' ## CAUTION: the following code creates some files in the current working
#' ## directory: two 'model_<hash>.stan' files, one 'model_<hash>(.exe)'
#' ## executable, and one 'fit_cmdstanr_<some_number>.rds' file
#' set.seed(7)
#' fname <- paste0("fit_cmdstanr_", sample.int(.Machine$integer.max, 1))
#' options(cmdstanr_write_stan_file_dir = getwd())
#' fit_cmdstanr <- brm(rate ~ conc + state, data = Puromycin,
#'                     backend = "cmdstanr", file = fname)
#'
#' # now restart the R session and run the following (after attaching 'brms')
#' set.seed(7)
#' fname <- paste0("fit_cmdstanr_", sample.int(.Machine$integer.max, 1))
#' fit_cmdstanr <- brm(rate ~ conc + state,
#'                     data = Puromycin,
#'                     backend = "cmdstanr",
#'                     file = fname)
#' kfold_cmdstanr <- kfold(fit_cmdstanr, K = 2)
#' }
#'
#' @seealso \code{\link{loo}}, \code{\link{reloo}}
#'
#' @importFrom loo kfold
#' @export kfold
#' @export
kfold.brmsfit <- function(x, ..., K = 10, Ksub = NULL, folds = NULL,
                          group = NULL, joint = FALSE, compare = TRUE,
                          resp = NULL, model_names = NULL, save_fits = FALSE,
                          recompile = NULL, future_args = list()) {
  args <- split_dots(x, ..., model_names = model_names)
  if (!"use_stored" %in% names(args)) {
    further_arg_names <- c(
      "K", "Ksub", "folds", "group", "joint", "resp", "save_fits"
    )
    args$use_stored <- all(names(args) %in% "models") &&
      !any(further_arg_names %in% names(match.call()))
  }
  c(args) <- nlist(
    criterion = "kfold", K, Ksub, folds, group, joint,
    compare, resp, save_fits, recompile, future_args
  )
  do_call(compute_loolist, args)
}

# helper function to perform k-fold cross-validation
# @inheritParams kfold.brmsfit
# @param model_name ignored but included to avoid being passed to '...'
.kfold <- function(x, K, Ksub, folds, group, joint, save_fits,
                   newdata, resp, model_name, recompile = NULL,
                   future_args = list(), newdata2 = NULL, ...) {
  stopifnot(is.brmsfit(x), is.list(future_args))
  if (is.brmsfit_multiple(x)) {
    warn_brmsfit_multiple(x)
    class(x) <- "brmsfit"
  }
  if (is.null(newdata)) {
    newdata <- x$data
  } else {
    newdata <- as.data.frame(newdata)
  }
  if (is.null(newdata2)) {
    newdata2 <- x$data2
  } else {
    bterms <- brmsterms(x$formula)
    newdata2 <- validate_data2(newdata2, bterms)
  }
  N <- nrow(newdata)
  joint <- as_one_logical(joint)
  # validate argument 'group'
  if (!is.null(group)) {
    valid_groups <- get_cat_vars(x)
    if (length(group) != 1L || !group %in% valid_groups) {
      stop2("Group '", group, "' is not a valid grouping factor. ",
            "Valid groups are: \n", collapse_comma(valid_groups))
    }
    gvar <- factor(get(group, newdata))
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
    opts <- c("loo", "stratified", "grouped")
    fold_type <- match.arg(folds, opts)
    req_group_opts <- c("stratified", "grouped")
    if (fold_type %in% req_group_opts && is.null(group)) {
      stop2("Argument 'group' is required for fold type '", fold_type, "'.")
    }
    if (fold_type == "loo") {
      folds <- seq_len(N)
      K <- N
      message("Setting 'K' to the number of observations (", K, ")")
    } else if (fold_type == "stratified") {
      folds <- loo::kfold_split_stratified(K, gvar)
    } else if (fold_type == "grouped") {
      folds <- loo::kfold_split_grouped(K, gvar)
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

  # ensure that the model can be run in the current R session
  x <- recompile_model(x, recompile = recompile)

  # split dots for use in log_lik and update
  dots <- list(...)
  ll_arg_names <- arg_names("log_lik")
  ll_args <- dots[intersect(names(dots), ll_arg_names)]
  ll_args$allow_new_levels <- TRUE
  ll_args$sample_new_levels <-
    first_not_null(ll_args$sample_new_levels, "gaussian")
  ll_args$resp <- resp
  ll_args$combine <- TRUE
  up_args <- dots[setdiff(names(dots), ll_arg_names)]
  up_args$object <- x
  up_args$refresh <- 0

  # function to be run inside future::future
  .kfold_k <- function(k) {
    message("Fitting model ", k, " out of ", K)
    if (fold_type == "loo" && !is.null(group)) {
      omitted <- which(folds == folds[k])
      predicted <- k
    } else {
      omitted <- predicted <- which(folds == k)
    }
    newdata_omitted <- newdata[-omitted, , drop = FALSE]
    up_args$newdata <- newdata_omitted
    up_args$data2 <- subset_data2(newdata2, -omitted)
    fit <- SW(do_call(update, up_args))

    ll_args$object <- fit
    ll_args$newdata <- newdata[predicted, , drop = FALSE]
    ll_args$newdata2 <- subset_data2(newdata2, predicted)
    lppds <- do_call(log_lik, ll_args)
    if (joint) {
      # compute the joint log score over all observations within a fold
      lppds <- rowSums(lppds)
    }

    out <- nlist(lppds, omitted, predicted)
    if (save_fits) {
      out$fit <- fit
    }
    return(out)
  }

  # TODO: separate parallel and non-parallel code to enable better printing?
  future_args$X <- Ksub
  future_args$FUN <- .kfold_k
  future_args$future.seed <- TRUE
  res <- do_call("future_lapply", future_args, pkg = "future.apply")

  lppds <- pred_obs_list <- vector("list", length(Ksub))
  if (save_fits) {
    fits <- array(list(), dim = c(length(Ksub), 3))
    dimnames(fits) <- list(NULL, c("fit", "omitted", "predicted"))
  }
  for (i in seq_along(Ksub)) {
    if (save_fits) {
      fits[i, ] <- res[[i]][c("fit", "omitted", "predicted")]
    }
    pred_obs_list[[i]] <- res[[i]]$predicted
    lppds[[i]] <- res[[i]]$lppds
  }

  lppds <- do_call(cbind, lppds)
  elpds <- apply(lppds, 2, log_mean_exp)
  pred_obs <- unlist(pred_obs_list)
  if (!joint) {
    # bring back elpds into the original observation order
    elpds <- elpds[order(pred_obs)]
  }

  # compute effective number of parameters
  ll_args$object <- x
  ll_args$newdata <- newdata
  ll_args$newdata2 <- newdata2
  pred_obs_sorted <- sort(pred_obs)
  if (length(Ksub) < K) {
    # select the correct subset of predicted observations in the original order
    ll_args$newdata <- ll_args$newdata[pred_obs_sorted, , drop = FALSE]
    ll_args$newdata2 <- subset_data2(ll_args$newdata2, pred_obs_sorted)
  }
  ll_full <- do_call(log_lik, ll_args)
  if (joint) {
    # compute the joint log score over all observations within a fold
    ll_full_marg <- matrix(nrow = nrow(ll_full), ncol = length(Ksub))
    for (i in seq_along(Ksub)) {
      sel_obs <- match(pred_obs_list[[i]], pred_obs_sorted)
      ll_full_marg[, i] <- rowSums(ll_full[, sel_obs, drop = FALSE])
    }
    ll_full <- ll_full_marg
  }
  lpds <- apply(ll_full, 2, log_mean_exp)
  ps <- lpds - elpds

  # put everything together in a loo object
  pointwise <- cbind(elpd_kfold = elpds, p_kfold = ps, kfoldic = -2 * elpds)
  est <- colSums(pointwise)
  se_est <- sqrt(nrow(pointwise) * apply(pointwise, 2, var))
  estimates <- cbind(Estimate = est, SE = se_est)
  rownames(estimates) <- colnames(pointwise)
  out <- nlist(estimates, pointwise)
  atts <- nlist(K, Ksub, group, folds, fold_type, joint)
  attributes(out)[names(atts)] <- atts
  if (save_fits) {
    out$fits <- fits
    out$data <- newdata
    out$data2 <- newdata2
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
  ndraws <- ndraws(x$fits[[1, "fit"]])
  y <- rep(NA, npredicted)
  yrep <- matrix(NA, nrow = ndraws, ncol = npredicted)
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
