brmsframe <- function(x, ...) {
  # TODO: store output of standata_basis in brmsframe?
  UseMethod("brmsframe")
}

#' @export
brmsframe.mvbrmsterms <- function(x, data, basis = NULL, ...) {
  x$frame <- initialize_frame(x, data = data, basis = basis, ...)
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
    x$frame <- initialize_frame(x, data = data, basis = basis, ...)
  } else {
    # this must be a multivariate model
    stopifnot(is.list(frame))
    x$frame <- frame
    x$frame$re <- subset(x$frame$re, resp = x$resp)
  }
  data <- subset_data(data, x)
  x$frame$resp <- frame_resp(x, data = data)
  x$frame$ac <- frame_ac(x, data = data)
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
  x$frame$sm <- frame_sm(x)
  x$frame$sp <- frame_sp(x, data = data)
  x$frame$gp <- frame_gp(x, data = data)
  x$frame$ac <- frame_ac(x, data = data)
  # only store the ranefs of this specific linear formula
  x$frame$re <- subset2(frame$re, ls = check_prefix(x))
  class(x) <- c("bframel", class(x))
  # these data_* functions require the outputs of the corresponding
  # frame_* functions (but not vice versa) and are thus evaluated last
  x$sdata$gp <- data_gp(x, data, internal = TRUE)
  x$sdata$offset <- data_offset(x, data)
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
  x$frame$ac <- frame_ac(x, data = data)
  class(x) <- c("bframenl", class(x))
  x
}

#' @export
brmsframe.default <- function(x, ...) {
  brmsframe(brmsterms(x), ...)
}

# initialize the frame list with general information
initialize_frame <- function(x, data, basis = NULL, ...) {
  out <- list(
    re = frame_re(x, data = data, old_levels = basis$levels),
    me = frame_me(x, data = data, old_levels = basis$levels),
    index = frame_index(x, data = data)
  )
  set_levels(out) <- get_levels(ls = out)
  # TODO: store both old and current (new) levels if basis$levels is provided?
  # this will enable to avoid repeated calls of frame_re and frame_me
  # in different places just to get the right levels
  out
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

is.brmsframe <- function(x) {
  inherits(x, "brmsframe")
}

is.mvbrmsframe <- function(x) {
  inherits(x, "mvbrmsframe")
}

# useful for functions that require either of the two objects
is.anybrmsframe <- function(x) {
  is.brmsframe(x) || is.mvbrmsframe(x)
}

is.bframel <- function(x) {
  inherits(x, "bframel")
}

is.bframenl <- function(x) {
  inherits(x, "bframenl")
}
