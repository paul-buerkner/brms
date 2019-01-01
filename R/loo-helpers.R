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

compute_ics <- function(models, ic = c("loo", "waic", "psis", "psislw", "kfold"),
                        use_stored_ic = FALSE, compare = TRUE, ...) {
  # helper function used to create (lists of) 'ic' objects
  # Args:
  #   models: list of brmsfit objects
  #   ic: name of the information criterion to compute
  #   use_stored_ic: use recomputed ic objects if possible?
  #   ...: more arguments passed to compute_ic
  # Returns:
  #   If length(models) > 1 an object of class 'iclist'
  #   If length(models) == 1 an object of class 'ic'
  ic <- match.arg(ic)
  args <- nlist(ic, ...)
  if (length(models) > 1L) {
    if (!match_nobs(models)) {
      stop2("Models have different number of observations.")
    }
    if (length(use_stored_ic) == 1L) {
      use_stored_ic <- rep(use_stored_ic, length(models))
    }
    out <- named_list(names(models))
    for (i in seq_along(models)) {
      ic_obj <- models[[i]][[ic]]
      if (use_stored_ic[i] && is.ic(ic_obj)) {
        out[[i]] <- ic_obj
        out[[i]]$model_name <- names(models)[i]
      } else {
        args$x <- models[[i]]
        args$model_name <- names(models)[i]
        out[[i]] <- do_call(compute_ic, args) 
      }
    }
    compare <- as_one_logical(compare)
    if (compare) {
      out <- compare_ic(x = out)
    }
    class(out) <- "iclist"
  } else {
    ic_obj <- models[[1]][[ic]]
    use_stored_ic <- as_one_logical(use_stored_ic)
    if (use_stored_ic && is.ic(ic_obj)) {
      out <- ic_obj
      out$model_name <- names(models)
    } else {
      args$x <- models[[1]]
      args$model_name <- names(models)
      out <- do_call(compute_ic, args) 
    }
  }
  out
}

compute_ic <- function(x, ic = c("loo", "waic", "psis", "kfold"),
                       reloo = FALSE, k_threshold = 0.7, reloo_args = list(),
                       pointwise = FALSE, model_name = "", ...) {
  # compute information criteria using the 'loo' package
  # Args:
  #   x: an object of class brmsfit
  #   ic: the information criterion to be computed
  #   model_name: original variable name of object 'x'
  #   reloo: call 'reloo' after computing 'loo'?
  #   reloo_args: list of arguments passed to 'reloo'
  #   pointwise: compute log-likelihood point-by-point?
  #   ...: passed to other post-processing methods
  # Returns:
  #   an object of class 'ic' which inherits from class 'loo'
  ic <- match.arg(ic)
  if (ic == "kfold") {
    out <- do_call(kfold_internal, list(x, ...))
  } else {
    contains_samples(x)
    pointwise <- as_one_logical(pointwise)
    loo_args <- list(...)
    loo_args$x <- log_lik(x, pointwise = pointwise, ...)
    if (pointwise) {
      loo_args$draws <- attr(loo_args$x, "draws")
      loo_args$data <- attr(loo_args$x, "data")
    }
    if (ic == "psis") {
      if (pointwise) {
        stop2("Cannot use pointwise evaluation for 'psis'.")
      }
      loo_args$log_ratios <- -loo_args$x
      loo_args$x <- NULL
    }
    out <- SW(do_call(ic, loo_args, pkg = "loo"))
  }
  out$model_name <- model_name
  class(out) <- c("ic", class(out))
  # TODO: fix hashing when new data is passed
  attr(out, "yhash") <- hash_response(x)
  if (ic == "loo") {
    if (reloo) {
      c(reloo_args) <- nlist(x = out, fit = x, k_threshold, check = FALSE)
      out <- do_call(reloo.loo, reloo_args)
    } else {
      n_bad_obs <- length(loo::pareto_k_ids(out, threshold = k_threshold))
      recommend_loo_options(n_bad_obs, k_threshold, model_name) 
    }
  }
  out
}

#' Compare Information Criteria of Different Models
#'
#' Compare information criteria of different models fitted
#' with \code{\link{WAIC}} or \code{\link{loo}}.
#' 
#' @param ... At least two objects returned by 
#'   \code{\link{WAIC}} or \code{\link{loo}}.
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
#' @details For more details see \code{\link[loo:compare]{compare}}.
#' 
#' @seealso 
#'   \code{\link{WAIC}}, 
#'   \code{\link{loo}},
#'   \code{\link{add_ic}},
#'   \code{\link[loo:compare]{compare}}
#'   
#' @examples 
#' \dontrun{
#' # model with population-level effects only
#' fit1 <- brm(rating ~ treat + period + carry,
#'             data = inhaler, family = "gaussian")
#' waic1 <- WAIC(fit1)
#' 
#' # model with an additional varying intercept for subjects
#' fit2 <- brm(rating ~ treat + period + carry + (1|subject),
#'             data = inhaler, family = "gaussian")
#' waic2 <- WAIC(fit2)
#' 
#' # compare both models
#' compare_ic(waic1, waic2)
#' }
#' 
#' @export
compare_ic <- function(..., x = NULL, ic = c("loo", "waic", "kfold")) {
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
  ics <- unname(sapply(x, function(y) rownames(y$estimates)[3]))
  if (!all(ics %in% ics[1])) {
    stop2("All inputs should be from the same criterion.")
  }
  if (ics[1] == "kfoldic") {
    Ks <- sapply(x, "[[", "K")
    if (!all(Ks %in% Ks[1])) {
      stop2("'K' differs across kfold objects.")
    }
    subs <- lengths(lapply(x, "[[", "Ksub"))
    subs <- ifelse(subs %in% 0, Ks, subs)
    if (!all(subs %in% subs[1])) {
      stop2("The number of subsets differs across kfold objects.")
    }
  }
  yhash <- lapply(x, attr, which = "yhash")
  yhash_check <- ulapply(yhash, is_equal, yhash[[1]])
  if (!all(yhash_check)) {
    warning2(
      "Model comparisons are likely invalid as the response ", 
      "values of at least two models do not match."
    )
  }
  names(x) <- ulapply(x, "[[", "model_name")
  n_models <- length(x)
  ic_diffs <- matrix(0, nrow = n_models * (n_models - 1) / 2, ncol = 2)
  rnames <- rep("", nrow(ic_diffs))
  # pairwise comparision to get differences in ICs and their SEs
  n <- 1
  for (i in seq_len(n_models - 1)) {
    for (j in (i + 1):n_models) {
      tmp <- loo::compare(x[[j]], x[[i]])
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

#' Add information criteria and fit indices to fitted model objects
#' 
#' @param x An \R object typically of class \code{brmsfit}.
#' @param ic,value Names of the information criteria / fit indices 
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
#'   \code{add_ic} or if \code{force_save} was set to \code{TRUE}.
#' @param force_save Logical; only relevant if \code{file} is specified and
#'   ignored otherwise. If \code{TRUE}, the fitted model object will be saved
#'   regardless of whether new criteria were added via \code{add_ic}.
#' @param ... Further arguments passed to the underlying 
#'   functions computing the information criteria or fit indices.
#'   
#' @return An object of the same class as \code{x}, but
#'   with information criteria added for later usage.
#'   
#' @details The methods \code{add_loo} and add \code{add_waic}
#'   are just convenient wrappers around \code{add_ic}.
#'   
#' @examples
#' \dontrun{
#' fit <- brm(count ~ Trt, epilepsy, poisson())
#' # add both LOO and WAIC at once
#' fit <- add_ic(fit, ic = c("loo", "waic"))
#' print(fit$loo)
#' print(fit$waic)
#' }
#' 
#' @export
add_ic <- function(x, ...) {
  UseMethod("add_ic")
}

#' @rdname add_ic
#' @export
add_ic.brmsfit <- function(x, ic = "loo", model_name = NULL, 
                           overwrite = FALSE, file = NULL,
                           force_save = FALSE, ...) {
  unused_args <- intersect(names(list(...)), args_not_for_reloo())
  if (length(unused_args)) {
    unused_args <- collapse_comma(unused_args)
    stop2("Cannot use arguments ", unused_args," in calls to 'add_ic'.")
  }
  if (!is.null(model_name)) {
    model_name <- as_one_character(model_name)
  } else {
    model_name <- deparse_combine(substitute(x)) 
  }
  ic <- unique(as.character(ic))
  valid_ics <- c("loo", "waic", "kfold", "R2", "marglik")
  if (!length(ic) || !all(ic %in% valid_ics)) {
    stop2("Argument 'ic' should be a subset of ",
          collapse_comma(valid_ics))
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
    x[ic] <- list(NULL)
  }
  new_ic <- ic[ulapply(x[ic], is.null)]
  args <- list(x, ...)
  for (fun in intersect(ic, c("loo", "waic", "kfold"))) {
    x[[fun]] <- do_call(fun, args)
    x[[fun]]$model_name <- model_name
  }
  if ("R2" %in% ic) {
    args$summary <- FALSE
    x$R2 <- do_call(bayes_R2, args)
  }
  if ("marglik" %in% ic) {
    x$marglik <- do_call(bridge_sampler, args)
  }
  if (!is.null(file) && (force_save || length(new_ic))) {
    if (auto_save) {
      message("Automatically saving the model object in '", file, "'")
    }
    x$file <- file
    saveRDS(x, file = file)
  } 
  x
}

#' @rdname add_ic 
#' @export
'add_ic<-' <- function(x, ..., value) {
  add_ic(x, ic = value, ...)
}

#' @rdname add_ic
#' @export
add_loo <- function(x, ...) {
  add_ic(x, ic = "loo", ...)
}

#' @rdname add_ic
#' @export
add_waic <- function(x, ...) {
  add_ic(x, ic = "waic", ...)
}

args_not_for_reloo <- function() {
  # arguments not usable with 'reloo'
  # the same arguments cannot be used in add_ic
  c("newdata", "re_formula", "allow_new_levels", 
    "sample_new_levels", "new_objects")
}

hash_response <- function(x, ...) {
  # create a hash based on the response of a model
  require_package("digest")
  stopifnot(is.brmsfit(x))
  sdata <- standata(x, internal = TRUE)
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
  #   TRUE if the response parts of all models match and FALSE else
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
#' @inheritParams loo.brmsfit
#' @param x An \R object typically of class \code{loo}.
#' @param fit An \R object typically of class \code{brmsfit}.
#' @param k_threshold The threshold at which pareto \eqn{k} 
#'   estimates are treated as problematic. Defaults to \code{0.7}. 
#'   See \code{\link[loo:pareto_k_ids]{pareto_k_ids}}
#'   for more details.
#' @param check Logical; If \code{TRUE} (the default), a crude 
#'   check is performed if the \code{loo} object was generated
#'   from the \code{brmsfit} object passed to argument \code{fit}.
#' @param ... Further arguments passed to 
#'   \code{\link{update.brmsfit}} such
#'   as \code{iter}, \code{chains}, or \code{cores}.
#'   
#' @return An object of the class as \code{x}.
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
#'            data = epilepsy, family = poisson())
#' # throws warning about some pareto k estimates being too high
#' (loo1 <- loo(fit1))
#' (loo1 <- reloo(loo1, fit1))
#' }
#' 
#' @export
reloo <- function(x, ...) {
  UseMethod("reloo")
}

#' @rdname reloo
#' @export
reloo.loo <- function(x, fit, k_threshold = 0.7, check = TRUE,
                      resp = NULL, ...) {
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
    fit_j <- subset_autocor(fit, -omitted)
    fit_j <- SW(update(fit_j, newdata = mf_omitted, refresh = 0, ...))
    fit_j <- subset_autocor(fit_j, omitted, autocor = x$autocor)
    lls[[j]] <- log_lik(
      fit_j, newdata = mf[omitted, , drop = FALSE],
      allow_new_levels = TRUE, resp = resp
    )
  }
  # compute elpd_{loo,j} for each of the held out observations
  elpd_loo <- unlist(lapply(lls, log_mean_exp))
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
print.iclist <- function(x, digits = 2, ...) {
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

#' @importFrom loo print_dims
#' @export
print_dims.kfold <- function(x, ...) {
  sub <- length(x$Ksub)
  sub <- ifelse(sub > 0 & sub < x$K, paste0(sub, " subsets of "), "")
  cat(paste0("Based on ", sub, x$K, "-fold cross-validation\n"))
}

is.loo <- function(x) {
  # class from the loo package
  inherits(x, "loo")
}

is.ic <- function(x) {
  # class 'ic' inherits from class 'loo'
  inherits(x, "ic")
}
