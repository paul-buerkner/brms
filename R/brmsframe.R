brmsframe <- function(x, ...) {
  # TODO: store output of standata_basis in brmsframe?
  UseMethod("brmsframe")
}

#' @export
brmsframe.mvbrmsterms <- function(x, data, basis = NULL, ...) {
  # this is a univariate model so brmsterms is at the top level
  x$frame <- list(
    re = tidy_ranef(x, data = data, old_levels = basis$levels),
    me = tidy_meef(x, data = data, old_levels = basis$levels),
    index = tidy_index(x, data = data)
  )
  for (r in names(x$terms)) {
    x$terms[[r]] <- brmsframe(
      x$terms[[r]], data = data, frame = x$frame,
      basis = basis$resps[[r]], ...
    )
  }
  class(x) <- c("mvbrmsframe", class(x))
  x
}

# @param basis information from original Stan data used to correctly
#   predict from new data. See 'standata_basis' for details.
#' @export
brmsframe.brmsterms <- function(x, data, frame = NULL, basis = NULL, ...) {
  if (is.null(frame)) {
    # this is a univariate model so brmsterms is at the top level
    x$frame <- list(
      re = tidy_ranef(x, data = data, old_levels = basis$levels),
      me = tidy_meef(x, data = data, old_levels = basis$levels),
      index = tidy_index(x, data = data)
    )
  } else {
    # this must be a multivariate model
    stopifnot(is.list(frame))
    x$frame <- frame
    x$frame$re <- subset(x$frame$re, resp = x$resp)
  }
  data <- subset_data(data, x)
  x$frame$resp <- frame_resp(x, data)
  x$frame$ac <- tidy_acef(x, data)
  for (dp in names(x$dpars)) {
    x$dpars[[dp]] <- brmsframe(
      x$dpars[[dp]], data, frame = x$frame,
      basis = basis$dpars[[dp]], ...
    )
  }
  for (nlp in names(x$nlpars)) {
    x$nlpars[[nlp]] <- brmsframe(
      x$nlpars[[nlp]], data, frame = x$frame,
      basis = basis$nlpars[[nlp]], ...
    )
  }
  class(x) <- c("brmsframe", class(x))
  x
}

#' @export
brmsframe.btl <- function(x, data, frame = list(), basis = NULL, ...) {
  # TODO: rename the tidy_ functions to frame_ functions?
  stopifnot(is.list(frame))
  # TODO: add more comments on the relation of data_* and frame_* functions
  # the outputs of these data_* functions are required in the corresponding
  # frame_* functions (but not vice versa) and are thus evaluated first
  x$frame <- frame
  x$basis <- basis
  x$sdata <- list(
    fe = data_fe(x, data),
    cs = data_cs(x, data),
    sm = data_sm(x, data)
  )
  # this enables overwriting of frames if necessary
  x$frame$fe <- frame_fe(x)
  x$frame$cs <- frame_cs(x)
  x$frame$sm <- tidy_smef(x)
  x$frame$sp <- tidy_spef(x, data = data)
  x$frame$gp <- tidy_gpef(x, data = data)
  x$frame$ac <- tidy_acef(x, data = data)
  # only store the ranefs of this specific linear formula
  px <- check_prefix(x)
  x$frame$re <- subset2(frame$re, ls = px)
  # these data_* functions require the outputs of the corresponding
  # frame_* functions (but not vice versa) and are thus evaluated last
  x$sdata$gp <- data_gp(x, data, internal = TRUE)
  x$sdata$offset <- data_offset(x, data)
  class(x) <- c("bfrl", class(x))
  x
}

#' @export
brmsframe.btnl <- function(x, data, frame = list(), basis = NULL, ...) {
  stopifnot(is.list(frame))
  x$frame <- frame
  x$basis <- basis
  x$sdata <- list(
    cnl = data_cnl(x, data)
  )
  x$frame$cnl <- frame_cnl(x)
  x$frame$ac <- tidy_acef(x)
  class(x) <- c("bfrnl", class(x))
  x
}

#' @export
brmsframe.default <- function(x, ...) {
  brmsframe(brmsterms(x), ...)
}

frame_resp <- function(x, data, ....) {
  stopifnot(is.brmsterms(x))
  y <- model.response(model.frame(x$respform, data, na.action = na.pass))
  out <- list(
    values = y,
    bounds = trunc_bounds(x, data),
    Ybounds = trunc_bounds(x, data, incl_family = TRUE, stan = TRUE),
    Jmi = as.array(which(is.na(y))),
    subset = attr(data, "subset")
  )
  out
}

frame_fe <- function(x, data = NULL, ...) {
  stopifnot(is.btl(x))
  sdata <- x$sdata$fe
  if (is.null(sdata)) {
    sdata <- data_fe(x, data)
  }
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

frame_cs <- function(x, data = NULL, ...) {
  stopifnot(is.btl(x))
  sdata <- x$sdata$cs
  if (is.null(sdata)) {
    sdata <- data_cs(x, data)
  }
  out <- list(vars = colnames(x$sdata$cs$Xcs))
  out
}

frame_cnl <- function(x, data, ...) {
  stopifnot(is.btnl(x))
  covars <- all.vars(x$covars)
  if (!length(covars)) {
    return(empty_data_frame())
  }
  sdata <- x$sdata$cnl
  if (is.null(sdata)) {
    sdata <- data_cnl(x, data)
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

is.bfrl <- function(x) {
  inherits(x, "bfrl")
}

is.bfrnl <- function(x) {
  inherits(x, "bfrnl")
}
