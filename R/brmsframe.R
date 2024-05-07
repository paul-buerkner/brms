brmsframe <- function(x, ...) {
  # TODO: store output of standata_basis in brmsframe?
  UseMethod("brmsframe")
}

#' @export
brmsframe.mvbrmsterms <- function(x, data, old_levels = NULL, ...) {
  # this is a univariate model so brmsterms is at the top level
  x$frame <- list(
    re = tidy_ranef(x, data = data, old_levels = old_levels),
    me = tidy_meef(x, data = data, old_levels = old_levels)
  )
  for (r in names(x$terms)) {
    x$terms[[r]] <- brmsframe(
      x$terms[[r]], data = data, frame = x$frame,
      old_levels = old_levels, ...
    )
  }
  class(x) <- c("mvbrmsframe", class(x))
  x
}

#' @export
brmsframe.brmsterms <- function(x, data, frame = NULL,
                                old_levels = NULL, ...) {
  if (is.null(frame)) {
    # this is a univariate model so brmsterms is at the top level
    x$frame <- list(
      re = tidy_ranef(x, data = data, old_levels = old_levels),
      me = tidy_meef(x, data = data, old_levels = old_levels)
    )
  } else {
    x$frame <- frame
    x$frame$re <- subset(x$frame$re, resp = x$resp)
  }
  data <- subset_data(data, x)
  x$frame$resp <- frame_resp(x, data)
  for (dp in names(x$dpars)) {
    x$dpars[[dp]] <- brmsframe(x$dpars[[dp]], data, frame = x$frame, ...)
  }
  for (nlp in names(x$nlpars)) {
    x$nlpars[[nlp]] <- brmsframe(x$nlpars[[nlp]], data, frame = x$frame, ...)
  }
  class(x) <- c("brmsframe", class(x))
  x
}

#' @export
brmsframe.btl <- function(x, data, frame = NULL, ...) {
  # TODO: rename the tidy_ functions to frame_ functions and move them here?
  # TODO: create a proper fixef data.frame as for the other terms?
  # TODO: store more the data_ functions?
  x$sdata <- list(
    fe = data_fe(x, data),
    sm = data_sm(x, data),
    gp = data_gp(x, data, internal = TRUE),
    offset = data_offset(x, data)
  )
  px <- check_prefix(x)
  x$frame <- list(
    fe = frame_fe(x),
    re = subset2(frame$re, ls = px),
    sp = tidy_spef(x, data),
    me = frame$me,
    cs = colnames(get_model_matrix(x$cs, data = data)),
    gp = tidy_gpef(x, data),
    sm = tidy_smef(x),
    ac = tidy_acef(x)
  )
  class(x) <- c("bfrl", class(x))
  x
}

#' @export
brmsframe.btnl <- function(x, data, ...) {
  x$sdata <- list(
    cnl = data_cnl(x, data)
  )
  x$frame <- list(
    cnl = frame_cnl(x),
    ac = tidy_acef(x)
  )
  class(x) <- c("bfrnl", class(x))
  x
}

#' @export
brmsframe.default <- function(x, ...) {
  brmsframe(brmsterms(x), ...)
}

frame_resp <- function(x, data, ....) {
  stopifnot(is.brmsterms(x))
  out <- list(
    values = model.response(model.frame(x$respform, data, na.action = na.pass)),
    bounds = trunc_bounds(x, data),
    Ybounds = trunc_bounds(x, data, incl_family = TRUE, stan = TRUE)
  )
  out
}

frame_fe <- function(x) {
  stopifnot(is.btl(x), !is.null(x$sdata))
  out <- list(
    vars = colnames(x$sdata$fe$X),
    center = stan_center_X(x),
    sparse = is_sparse(x$fe),
    decomp = get_decomp(x$fe)
  )
  out$vars_stan <- out$vars
  if (out$center) {
    out$vars_stan <- setdiff(out$vars_stan, "Intercept")
  }
  out
}

frame_cnl <- function(x, ...) {
  stopifnot(is.btnl(x), !is.null(x$sdata))
  covars <- all.vars(x$covars)
  if (!length(covars)) {
    return(empty_data_frame())
  }
  out <- data.frame(
    covar = covars, integer = FALSE,
    matrix = FALSE, dim2 = 0
  )
  p <- usc(combine_prefix(x))
  for (i in seq_along(covars)) {
    cname <- glue("C{p}_{i}")
    cvalues <- x$sdata$cnl[[cname]]
    out$integer[i] <- is.integer(cvalues)
    out$matrix[i] <- is.matrix(cvalues)
    if (out$matrix[i]) {
      out$dim2[i] <- dim(cvalues)[2]
    }
  }
  out
}
