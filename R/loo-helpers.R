compute_ics <- function(models, model_names, 
                        ic = c("loo", "waic", "psislw", "kfold"),
                        use_stored_ic = FALSE, 
                        compare = TRUE, ...) {
  # helper function used to create (lists of) 'ic' objects
  # Args:
  #   models: list of brmsfit objects
  #   model_names: character vector containing model names
  #   ic: name of the information criterion to compute
  #   use_stored_ic: use recomputed ic objects if possible?
  #   ...: more arguments passed to compute_ic
  # Returns:
  #   If length(models) > 1 an object of class 'iclist'
  #   If length(models) == 1 an object of class 'ic'
  ic <- match.arg(ic)
  args <- nlist(ic, ...)
  if (length(models) > 1L) {
    stopifnot(length(models) == length(model_names))
    if (length(use_stored_ic) == 1L) {
      use_stored_ic <- rep(use_stored_ic, length(models))
    }
    out <- named_list(model_names)
    for (i in seq_along(models)) {
      # verify all model objects before computing any ICs
      if (!is.brmsfit(models[[i]])) {
        stop2("Object '", model_names[i], "' is not of class 'brmsfit'.")
      }
    }
    for (i in seq_along(models)) {
      ic_obj <- models[[i]][[ic]]
      if (use_stored_ic[i] && is.ic(ic_obj)) {
        out[[i]] <- ic_obj
      } else {
        args$x <- models[[i]]
        args$model_name <- model_names[i]
        out[[i]] <- do.call(compute_ic, args) 
      }
    }
    if (compare) {
      match_response(models)
      out <- compare_ic(x = out)
    }
    class(out) <- "iclist"
  } else {
    ic_obj <- models[[1]][[ic]]
    stopifnot(length(use_stored_ic) == 1L)
    if (use_stored_ic && is.ic(ic_obj)) {
      out <- ic_obj
    } else {
      args$x <- models[[1]]
      args$model_name <- model_names
      out <- do.call(compute_ic, args) 
    }
  }
  out
}

compute_ic <- function(x, ic = c("loo", "waic", "psislw", "kfold"),
                       model_name = "", reloo = FALSE, k_threshold = 0.7,
                       loo_args = list(), update_args = list(), ...) {
  # compute information criteria using the 'loo' package
  # Args:
  #   x: an object of class brmsfit
  #   ic: the information criterion to be computed
  #   model_name: original variable name of object 'x'
  #   reloo: call 'reloo' after computing 'loo'?
  #   loo_args: passed to functions of the loo package
  #   update_args: passed to update.brmsfit
  #   ...: passed to log_lik.brmsfit
  # Returns:
  #   an object of class 'ic' which inherits from class 'loo'
  stopifnot(is.list(loo_args))
  ic <- match.arg(ic)
  contains_samples(x)
  if (ic == "kfold") {
    IC <- do.call(kfold_internal, c(list(x), update_args))
  } else {
    loo_args$x <- log_lik(x, ...)
    pointwise <- is.function(loo_args$x)
    if (pointwise) {
      loo_args$args <- attr(loo_args$x, "args")
      attr(loo_args$x, "args") <- NULL
    }
    if (ic == "psislw") {
      if (pointwise) {
        loo_args[["llfun"]] <- loo_args[["x"]]
        loo_args[["llargs"]] <- loo_args[["args"]]
        loo_args[["x"]] <- loo_args[["args"]] <- NULL
      } else {
        loo_args[["lw"]] <- -loo_args[["x"]]
        loo_args[["x"]] <- NULL
      }
    }
    IC <- SW(do.call(eval(parse(text = paste0("loo::", ic))), loo_args)) 
  }
  IC$model_name <- model_name
  class(IC) <- c("ic", class(IC))
  if (ic == "loo") {
    if (reloo) {
      reloo_args <- nlist(x = IC, fit = x, k_threshold, check = FALSE)
      IC <- do.call(reloo.loo, c(reloo_args, update_args))
    } else {
      n_bad_obs <- length(loo::pareto_k_ids(IC, threshold = k_threshold))
      recommend_loo_options(n_bad_obs, model_name) 
    }
  }
  IC
}

#' Compare Information Criteria of Different Models
#'
#' Compare information criteria of different models fitted
#' with \code{\link[brms:WAIC]{WAIC}} or \code{\link[brms:LOO]{LOO}}.
#' 
#' @param ... At least two objects returned by 
#'   \code{\link[brms:WAIC]{WAIC}} or \code{\link[brms:LOO]{LOO}}.
#'   Alternatively, \code{brmsfit} objects with information 
#'   criteria precomputed via \code{\link[brms:add_ic]{add_ic}}
#'   may be passed, as well.
#' @param x A \code{list} containing the same types of objects as
#'   can be passed via \code{...}.
#' @param ic The name of the information criterion to be extracted 
#'   from \code{brmsfit} objects. Ignored if information 
#'   criterion objects are only passed directly.
#'   
#' @return An object of class \code{iclist}.
#' 
#' @details For more details see \code{\link[loo:compare]{compare}}.
#' 
#' @seealso 
#'   \code{\link[brms:WAIC]{WAIC}}, 
#'   \code{\link[brms:LOO]{LOO}},
#'   \code{\link[brms:add_ic]{add_ic}},
#'   \code{\link[loo:compare]{compare}}
#'   
#' @examples 
#' \dontrun{
#' # model with population-level effects only
#' fit1 <- brm(rating ~ treat + period + carry,
#'             data = inhaler, family = "gaussian")
#' w1 <- WAIC(fit1)
#' 
#' # model with an additional varying intercept for subjects
#' fit2 <- brm(rating ~ treat + period + carry + (1|subject),
#'             data = inhaler, family = "gaussian")
#' w2 <- WAIC(fit2)
#' 
#' # compare both models
#' compare_ic(w1, w2)
#' }
#' 
#' @export
compare_ic <- function(..., x = NULL, ic = c("loo", "waic")) {
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
  if (!all(sapply(x, inherits, "ic"))) {
    stop2("All inputs should have class 'ic' ", 
          "or contain precomputed 'ic' objects.")
  }
  if (length(x) < 2L) {
    stop2("Expecting at least two objects.")
  }
  ics <- unname(sapply(x, function(y) names(y)[3]))
  if (!all(ics %in% ics[1])) {
    stop2("All inputs should be from the same criterion.")
  }
  if (ics[1] == "kfoldic") {
    Ks <- sapply(x, "[[", "K")
    if (!all(Ks %in% Ks[1])) {
      stop2("'K' differs across kfold objects.")
    }
  }
  names(x) <- ulapply(x, "[[", "model_name")
  n_models <- length(x)
  ic_diffs <- matrix(0, nrow = n_models * (n_models - 1) / 2, ncol = 2)
  rnames <- rep("", nrow(ic_diffs))
  # pairwise comparision to get differences in ICs and their SEs
  n <- 1
  for (i in seq_len(n_models - 1)) {
    for (j in (i + 1):n_models) {
      temp <- loo::compare(x[[j]], x[[i]])
      ic_diffs[n, ] <- c(-2 * temp[["elpd_diff"]], 2 * temp[["se"]]) 
      rnames[n] <- paste(names(x)[i], "-", names(x)[j])
      n <- n + 1
    }
  }
  rownames(ic_diffs) <- rnames
  colnames(ic_diffs) <- c(toupper(ics[1]), "SE")
  x$ic_diffs__ <- ic_diffs
  class(x) <- c("iclist", "list")
  x
}

#' @rdname add_ic
#' @export
add_ic.brmsfit <- function(x, ic = "loo", ...) {
  dots <- list(...)
  unused_args <- intersect(names(dots), args_not_for_reloo())
  if (length(unused_args)) {
    unused_args <- collapse_comma(unused_args)
    stop2("Cannot use arguments ", unused_args," in calls to 'add_ic'.")
  }
  model_name <- deparse(substitute(x))
  ic <- unique(tolower(as.character(ic)))
  valid_ics <- c("loo", "waic", "kfold")
  if (!length(ic) || !all(ic %in% valid_ics)) {
    stop2("Argument 'ic' should be a subset of ",
          collapse_comma(valid_ics))
  }
  for (i in seq_along(ic)) {
    x[[ic[i]]] <- do.call(ic[i], c(list(x), dots))
    x[[ic[i]]]$model_name <- model_name
  }
  x
}

#' @rdname add_loo
#' @export
add_loo.brmsfit <- function(x, ...) {
  add_ic(x, ic = "loo", ...)
}

#' @rdname add_waic
#' @export
add_waic.brmsfit <- function(x, ...) {
  add_ic(x, ic = "waic", ...)
}

loo_weights <- function(object, lw = NULL, log = FALSE, 
                        loo_args = list(), ...) {
  # compute loo weights for use in loo_predict and related methods
  # Args:
  #   object: a brmsfit object
  #   lw: precomputed log weights matrix
  #   log: return log weights?
  #   loo_args: further arguments passed to functions of loo
  #   ...: further arguments passed to compute_ic
  # Returns:
  #   an S x N matrix
  if (!is.null(lw)) {
    stopifnot(is.matrix(lw))
  } else {
    message("Running PSIS to compute weights")
    psis <- compute_ic(object, ic = "psislw", loo_args = loo_args, ...)
    lw <- psis[["lw_smooth"]]
  }
  if (!log) {
    lw <- exp(lw) 
  } 
  lw
}

set_pointwise <- function(x, pointwise = NULL, newdata = NULL, 
                          subset = NULL, thres = 1e+08) {
  # set the pointwise argument based on the model size
  # Args:
  #   x: a brmsfit object
  #   newdata: optional data.frame containing new data
  #   subset: a vector to indicate a subset of the posterior samples
  #   thres: threshold above which pointwise is set to TRUE
  # Returns:
  #   TRUE or FALSE
  if (!is.null(pointwise)) {
    pointwise <- as.logical(pointwise)
    if (length(pointwise) != 1L || anyNA(pointwise)) {
      stop2("Argument 'pointwise' must be either TRUE or FALSE.")
    }
  } else {
    nsamples <- nsamples(x, subset = subset)
    if (is.data.frame(newdata)) {
      nobs <- nrow(newdata)
    } else {
      nobs <- nobs(x)
    }
    pointwise <- nsamples * nobs > thres
    if (pointwise) {
      message(
        "Switching to pointwise evaluation to reduce ",  
        "RAM requirements.\nThis will likely increase ",
        "computation time."
      )
    }
  }
  pointwise
}

args_not_for_reloo <- function() {
  # arguments not usable with 'reloo'
  # the same arguments cannot be used in add_ic
  c("newdata", "re_formula", "subset", "nsamples",
    "allow_new_levels", "sample_new_levels", "new_objects")
}

match_response <- function(models) {
  # compare the response parts of multiple brmsfit objects
  # Args:
  #   models: A list of brmsfit objects
  # Returns:
  #   TRUE if the response parts of all models match and FALSE else
  if (length(models) <= 1L) {
    out <- TRUE  
  } else {
    add_funs <- lsp("brms", what = "exports", pattern = "^resp_")
    match_vars <- c("Y", sub("^resp_", "", add_funs))
    .match_fun <- function(x, y) {
      # checks if all relevant parts of the response are the same 
      # Args:
      #   x, y: named lists as returned by standata
      old_order_x <- attr(x, "old_order")
      old_order_y <- attr(y, "old_order")
      all(ulapply(match_vars, function(v) {
        a <- p(as.vector(x[[v]]), old_order_x)
        b <- p(as.vector(y[[v]]), old_order_y)
        is_equal(a, b)
      }))
    }
    sdatas <- lapply(models, standata, control = list(save_order = TRUE))
    matches <- ulapply(sdatas[-1], .match_fun, y = sdatas[[1]]) 
    if (all(matches)) {
      out <- TRUE
    } else {
      out <- FALSE
      warning2("Model comparisons are likely invalid as the response ", 
               "parts of at least two models do not match.")
    }
  }
  invisible(out)
}

#' @rdname reloo
#' @export
reloo.loo <- function(x, fit, k_threshold = 0.7, check = TRUE, ...) {
  # most of the code is taken from rstanarm:::reloo
  stopifnot(is.brmsfit(fit))
  model_name <- deparse(substitute(fit))
  if (check && !is_equal(model_name, x$model_name)) {
    loo_name <- deparse(substitute(x))
    stop2(
      "Object '", loo_name, "' appears to be generated from ",
      "a brmsfit object other than '", model_name, "'. ",
      "If this is a false positive, please set 'check' to FALSE."
    )
  }
  if (is.null(x$pareto_k)) {
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
  mf <- model.frame(fit)
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
    fit_j <- SW(update(fit, newdata = mf_omitted, refresh = 0, ...))
    lls[[j]] <- log_lik(
      fit_j, newdata = mf[omitted, , drop = FALSE],
      allow_new_levels = TRUE
    )
  }
  # compute elpd_{loo,j} for each of the held out observations
  elpd_loo <- unlist(lapply(lls, log_mean_exp))
  # compute \hat{lpd}_j for each of the held out observations (using log-lik
  # matrix from full posterior, not the leave-one-out posteriors)
  ll_x <- log_lik(fit, newdata = mf[obs,, drop=FALSE])
  hat_lpd <- apply(ll_x, 2, log_mean_exp)
  # compute effective number of parameters
  p_loo <- hat_lpd - elpd_loo
  # replace parts of the loo object with these computed quantities
  sel <- c("elpd_loo", "p_loo", "looic")
  x$pointwise[obs, sel] <- cbind(elpd_loo, p_loo, -2 * elpd_loo)
  x[sel] <- colSums(x$pointwise[, sel])
  x[paste0("se_", sel)] <- 
    sqrt(nrow(x$pointwise) * apply(x$pointwise[, sel], 2, var))
  # what should we do about pareto k? for now setting them to 0
  x$pareto_k[obs] <- 0
  x
}

kfold_internal <- function(x, K = 10, save_fits = FALSE, ...) {
  # most of the code is taken from rstanarm::kfold
  stopifnot(is.brmsfit(x))
  N <- nobs(x)
  if (K < 1 || K > N) {
    stop2("'K' must be greater than one and smaller or ", 
          "equal to the number of observations in the model.")
  }
  mf <- model.frame(x)
  perm <- sample.int(N)
  idx <- ceiling(seq(from = 1, to = N, length.out = K + 1))
  bin <- .bincode(perm, breaks = idx, right = FALSE, include.lowest = TRUE)
  
  lppds <- list()
  if (save_fits) {
    fits <- array(list(), c(K, 2), list(NULL, c("fit", "omitted")))    
  }
  for (k in seq_len(K)) {
    message("Fitting model ", k, " out of ", K)
    omitted <- which(bin == k)
    mf_omitted <- mf[-omitted, , drop = FALSE]
    fit_k <- SW(update(x, newdata = mf_omitted, refresh = 0, ...))
    lppds[[k]] <- log_lik(
      fit_k, newdata = mf[omitted, , drop = FALSE], 
      allow_new_levels = TRUE
    )
    if (save_fits) {
      fits[k, ] <- list(fit = fit_k, omitted = omitted) 
    }
  }
  elpds <- ulapply(lppds, function(x) apply(x, 2, log_mean_exp))
  elpd_kfold <- sum(elpds)
  se_elpd_kfold <- sqrt(N * var(elpds))
  out <- nlist(
    elpd_kfold, p_kfold = NA, kfoldic = - 2 * elpd_kfold,
    se_elpd_kfold, se_p_kfold = NA, se_kfoldic = 2 * se_elpd_kfold,
    pointwise = cbind(elpd_kfold = elpds),
    K = K, model_name = deparse(substitute(x))
  )
  if (save_fits) {
    out$fits <- fits 
  }
  structure(out, class = "loo")
}

recommend_loo_options <- function(n, model_name = "") {
  model_name <- if (isTRUE(nzchar(model_name))) {
    paste0(" in model '", model_name, "'")
  }
  if (n > 0 && n <= 10) {
    warning2(
      "Found ", n, " observations with a pareto_k > 0.7", model_name, ". ",
      "It is recommended to set 'reloo = TRUE' in order to calculate ",
      "the ELPD without the assumption that these observations are ", 
      "negligible. This will refit the model ", n, " times to compute ", 
      "the ELPDs for the problematic observations directly."
    )
    out <- "reloo"
  } else if (n > 10) {
    warning2(
      "Found ", n, " observations with a pareto_k > 0.7", model_name, ". ",
      "With this many problematic observations, it may be more ", 
      "appropriate to use 'kfold' with argument 'K = 10' to perform ", 
      "10-fold cross-validation rather than LOO."
    )
    out <- "kfold"
  } else {
    out <- "loo"
  }
  invisible(out)
}

#' @export
print.ic <- function(x, digits = 2, ...) {
  # print the output of LOO(x) and WAIC(x)
  ic <- names(x)[3]
  mat <- matrix(
    c(x[[ic]], x[[paste0("se_",ic)]]), ncol = 2, 
      dimnames = list("", c(toupper(ic), "SE"))
  )
  print(round(mat, digits = digits))
  if (is_equal(ic, "kfoldic")) {
    cat(paste0("\nBased on ", x$K, "-fold cross-validation\n"))
  }
  invisible(x)
}

#' @export
print.iclist <- function(x, digits = 2, ...) {
  # print the output of LOO and WAIC with multiple models
  m <- x
  m$ic_diffs__ <- NULL
  if (length(m)) {
    ic <- names(m[[1]])[3]
    mat <- matrix(0, nrow = length(m), ncol = 2)
    dimnames(mat) <- list(names(m), c(toupper(ic), "SE"))
    for (i in seq_along(m)) { 
      mat[i, ] <- c(m[[i]][[ic]], m[[i]][[paste0("se_", ic)]])
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
  if (is_equal(ic, "kfoldic")) {
    cat(paste0("\nBased on ", x[[1]]$K, "-fold cross-validation\n"))
  }
  invisible(x)
}

is.loo <- function(x) {
  inherits(x, "loo")
}

is.ic <- function(x) {
  # objects of class 'ic' are returned by LOO and WAIC
  inherits(x, "ic")
}
