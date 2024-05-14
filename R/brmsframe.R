# The brmsframe methods are combining formula with data information in
# such a way that it can be used across the full range of primary brms
# functions, including brm, stancode, standata, prepare_predictions etc.
# Before brmsframe was introduced, a lot of of the frame_ functions had to
# be run many times at different places of the pre- and post-processing functions,
# which was uncessarily wasteful. To avoid this, brmsframe also computes
# some parts of the Stan data already (see brmsframe.btl), which is automatically
# reused in standata to avoid redundant function evaluations.
brmsframe <- function(x, ...) {
  UseMethod("brmsframe")
}

# @param basis information from original Stan data used to correctly
#   predict from newdata. See 'frame_basis' for details.
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

# This methods handles the intricate relationship of frame_ and data_ functions.
# In some cases, frame_ functions are computed from their data_ function's
# output, but in some other cases, it is the other way around. This is slightly
# inconsistent but avoids code duplication as much as possible, reflecting
# the different ways formula terms are evaluated in brms
#' @export
brmsframe.btl <- function(x, data, frame = list(), basis = NULL, ...) {
  stopifnot(is.list(frame))
  # the outputs of these data_ functions are required in the corresponding
  # frame_ functions (but not vice versa) and are thus evaluated first
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
  # these data_ functions may require the outputs of the corresponding
  # frame_ functions (but not vice versa) and are thus evaluated last
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

# initialize the $frame list with general information
initialize_frame <- function(x, data, basis = NULL, ...) {
  old_levels <- basis$group_levels
  out <- list(
    re = frame_re(x, data = data, old_levels = old_levels),
    me = frame_me(x, data = data, old_levels = old_levels),
    index = frame_index(x, data = data)
  )
  if (!is.null(old_levels)) {
    # this can only happen in post-processing potentially with newdata
    # knowing both new and old indices in important in prepare_predictions
    set_levels(out) <- old_levels
    set_levels(out, "used") <- get_levels(out, prefix = "used")
  } else {
    set_levels(out) <- get_levels(out)
  }
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

# assignment function to store levels as an attribute
'set_levels<-' <- function(x, prefix = "", value) {
  prefix_ <- usc(prefix, "suffix")
  attr_name <- paste0(prefix_, "levels")
  attr(x, attr_name) <- value
  x
}

# extract list of levels with one element per grouping factor
# assumes that levels have been stored as a 'levels' attribute
get_levels <- function(x, ...) {
  UseMethod("get_levels")
}

#' @export
get_levels.default <- function(x, prefix = "", ...) {
  prefix_ <- usc(prefix, "suffix")
  attr_name <- paste0(prefix_, "levels")
  attr(x, attr_name, exact = TRUE)
}

#' @export
get_levels.list <- function(x, ...) {
  out <- get_levels.default(x, ...)
  if (!is.null(out)) {
    return(out)
  }
  out <- vector("list", length(x))
  for (i in seq_along(out)) {
    levels <- get_levels(x[[i]], ...)
    if (is.list(levels)) {
      stopifnot(!is.null(names(levels)))
      out[[i]] <- as.list(levels)
    } else if (!is.null(levels)) {
      stopifnot(isTRUE(nzchar(names(x)[i])))
      out[[i]] <- setNames(list(levels), names(x)[[i]])
    }
  }
  out <- unlist(out, recursive = FALSE)
  out[!duplicated(names(out))]
}

#' @export
get_levels.brmsterms <- function(x, data = NULL, ...) {
  # if available, precomputed levels are stored in x$frame
  out <- get_levels(x$frame, ...)
  if (!is.null(out)) {
    return(out)
  }
  if (!is.null(data)) {
    ls <- list(frame_re(x, data), frame_me(x, data))
    out <- get_levels(ls)
  }
  out
}

#' @export
get_levels.mvbrmsterms <- function(x, data = NULL, ...) {
  get_levels.brmsterms(x, data = data, ...)
}

# prepare basis data required for correct predictions from new data
# TODO: eventually export this function if we want to ensure full compatibility
#   with the 'empty' feature. see ?rename_pars for an example
frame_basis <- function(x, data, ...) {
  UseMethod("frame_basis")
}

#' @export
frame_basis.default <- function(x, data, ...) {
  list()
}

#' @export
frame_basis.mvbrmsterms <- function(x, data, ...) {
  out <- list()
  # old levels are required to select the right indices for new levels
  levels <- get_levels(x, data = data)
  for (r in names(x$terms)) {
    out$resps[[r]] <- frame_basis(x$terms[[r]], data, levels = levels, ...)
  }
  # store levels as list element rather than as attribute (via set_levels)
  # to differentiate more easily whether or not old levels were provided
  out$group_levels <- levels
  out
}

#' @export
frame_basis.brmsterms <- function(x, data, levels = NULL, ...) {
  out <- list()
  data <- subset_data(data, x)
  for (dp in names(x$dpars)) {
    out$dpars[[dp]] <- frame_basis(x$dpars[[dp]], data, ...)
  }
  for (nlp in names(x$nlpars)) {
    out$nlpars[[nlp]] <- frame_basis(x$nlpars[[nlp]], data, ...)
  }
  # old levels are required to select the right indices for new levels
  if (is.null(levels)) {
    levels <- get_levels(x, data = data)
  }
  # store levels as list element rather than as attribute (via set_levels)
  # to differentiate more easily whether or not old levels were provided
  out$group_levels <- levels
  if (is_binary(x$family) || is_categorical(x$family)) {
    y <- model.response(model.frame(x$respform, data, na.action = na.pass))
    out$resp_levels <- levels(as.factor(y))
  }
  out
}

#' @export
frame_basis.btnl <- function(x, data, ...) {
  list()
}

#' @export
frame_basis.btl <- function(x, data, ...) {
  out <- list()
  out$sm <- frame_basis_sm(x, data, ...)
  out$gp <- frame_basis_gp(x, data, ...)
  out$sp <- frame_basis_sp(x, data, ...)
  out$ac <- frame_basis_ac(x, data, ...)
  out$bhaz <- frame_basis_bhaz(x, data, ...)
  out
}

# prepare basis data related to smooth terms
frame_basis_sm <- function(x, data, ...) {
  stopifnot(is.btl(x))
  smterms <- all_terms(x[["sm"]])
  out <- named_list(smterms)
  if (length(smterms)) {
    knots <- get_knots(data)
    data <- rm_attr(data, "terms")
    # the spline penalty has changed in 2.8.7 (#646)
    diagonal.penalty <- !require_old_default("2.8.7")
    gam_args <- list(
      data = data, knots = knots,
      absorb.cons = TRUE, modCon = 3,
      diagonal.penalty = diagonal.penalty
    )
    for (i in seq_along(smterms)) {
      sc_args <- c(list(eval2(smterms[i])), gam_args)
      sm <- do_call(smoothCon, sc_args)
      re <- vector("list", length(sm))
      for (j in seq_along(sm)) {
        re[[j]] <- mgcv::smooth2random(sm[[j]], names(data), type = 2)
      }
      out[[i]]$sm <- sm
      out[[i]]$re <- re
    }
  }
  out
}

# prepare basis data related to gaussian processes
frame_basis_gp <- function(x, data, ...) {
  stopifnot(is.btl(x))
  out <- data_gp(x, data, internal = TRUE)
  out <- out[grepl("^((Xgp)|(dmax)|(cmeans))", names(out))]
  out
}

# prepare basis data related to special terms
frame_basis_sp <- function(x, data, ...) {
  stopifnot(is.btl(x))
  out <- list()
  if (length(attr(x$sp, "uni_mo"))) {
    # do it like data_sp()
    spframe <- frame_sp(x, data)
    Xmo <- lapply(unlist(spframe$calls_mo), get_mo_values, data = data)
    out$Jmo <- as.array(ulapply(Xmo, attr, "max"))
  }
  out
}

# prepare basis data related to autocorrelation structures
frame_basis_ac <- function(x, data, ...) {
  out <- list()
  if (has_ac_class(x, "car")) {
    gr <- get_ac_vars(x, "gr", class = "car")
    if (isTRUE(nzchar(gr))) {
      out$locations <- extract_levels(get(gr, data))
    } else {
      out$locations <- NA
    }
  }
  if (has_ac_class(x, "unstr")) {
    time <- get_ac_vars(x, "time", dim = "time")
    out$times <- extract_levels(get(time, data))
  }
  out
}

# prepare basis data for baseline hazards of the cox model
frame_basis_bhaz <- function(x, data, ...) {
  out <- list()
  if (is_cox(x$family)) {
    # compute basis matrix of the baseline hazard for the Cox model
    y <- model.response(model.frame(x$respform, data, na.action = na.pass))
    out$basis_matrix <- bhaz_basis_matrix(y, args = x$family$bhaz)
  }
  out
}

