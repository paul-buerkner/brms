# type-stable indexing of vector and matrix type objects
# @param x an R object typically a vector or matrix
# @param i optional index; if NULL, x is returned unchanged
# @param row indicating if rows or cols should be indexed
#   only relevant if x has two or three dimensions
p <- function(x, i = NULL, row = TRUE) {
  # TODO: replace by "slice"
  if (isTRUE(length(dim(x)) > 3L)) {
    stop2("'p' can only handle objects up to 3 dimensions.")
  }
  if (!length(i)) {
    out <- x
  } else if (length(dim(x)) == 2L) {
    if (row) {
      out <- x[i, , drop = FALSE]
    } else {
      out <- x[, i, drop = FALSE]
    }
  } else if (length(dim(x)) == 3L) {
    if (row) {
      out <- x[i, , , drop = FALSE]
    } else {
      out <- x[, i, , drop = FALSE]
    }
  } else {
    out <- x[i]
  }
  out
}

# extract parts of an object with selective dropping of dimensions
# @param x,...,drop same as in x[..., drop]
# @param drop_dim Optional numeric or logical vector controlling
#   which dimensions to drop. Will overwrite argument 'drop'.
extract <- function(x, ..., drop = FALSE, drop_dim = NULL) {
  if (!length(dim(x))) {
    return(x[...])
  }
  if (length(drop_dim)) {
    drop <- FALSE
  } else {
    drop <- as_one_logical(drop)
  }
  out <- x[..., drop = drop]
  if (drop || !length(drop_dim) || any(dim(out) == 0L)) {
    return(out)
  }
  if (is.numeric(drop_dim)) {
    drop_dim <- seq_along(dim(x)) %in% drop_dim
  }
  if (!is.logical(drop_dim)) {
    stop2("'drop_dim' needs to be logical or numeric.")
  }
  keep <- dim(out) > 1L | !drop_dim
  new_dim <- dim(out)[keep]
  if (length(new_dim) <= 1L) {
    # use vectors instead of 1D arrays
    new_dim <- NULL
  }
  dim(out) <- new_dim
  out
}

# extract slices of one array dimension without dropping other dimensions
# @param x an array
# @param dim dimension from which to take the slice
# @param i slice index
# @param drop Logical (length 1) indicating whether to drop dimension `dim`.
slice <- function(x, dim, i, drop = TRUE) {
  ndim <- length(dim(x))
  commas1 <- collapse(rep(", ", dim - 1))
  commas2 <- collapse(rep(", ", ndim - dim))
  drop_dim <- ifelse(drop, ", drop_dim = dim", "")
  expr <- paste0("extract(x, ", commas1, "i", commas2, drop_dim, ")")
  eval2(expr)
}

# slice out columns without dropping other dimensions
# @param x an array; a vector or 1D array is treated as already sliced
# @param i column index
slice_col <- function(x, i) {
  if (length(dim(x)) < 2L) {
    # a vector or 1D array is treated as already sliced
    return(x)
  }
  slice(x, 2, i)
}

seq_rows <- function(x) {
  seq_len(NROW(x))
}

seq_cols <- function(x) {
  seq_len(NCOL(x))
}

seq_dim <- function(x, dim) {
  dim <- as_one_numeric(dim)
  if (dim == 1) {
    len <- NROW(x)
  } else if (dim == 2) {
    len <- NCOL(x)
  } else {
    len <- dim(x)[dim]
  }
  if (length(len) == 1L && !isNA(len)) {
    out <- seq_len(len)
  } else {
    out <- integer(0)
  }
  out
}

# match rows in x with rows in y
match_rows <- function(x, y, ...) {
  x <- as.data.frame(x)
  y <- as.data.frame(y)
  x <- do.call("paste", c(x, sep = "\r"))
  y <- do.call("paste", c(y, sep = "\r"))
  match(x, y, ...)
}

# find elements of 'x' matching sub-elements passed via 'ls' and '...'
find_elements <- function(x, ..., ls = list(), fun = '%in%') {
  x <- as.list(x)
  if (!length(x)) {
    return(logical(0))
  }
  out <- rep(TRUE, length(x))
  ls <- c(ls, list(...))
  if (!length(ls)) {
    return(out)
  }
  if (is.null(names(ls))) {
    stop("Argument 'ls' must be named.")
  }
  for (name in names(ls)) {
    tmp <- from_list(x, name)
    out <- out & do_call(fun, list(tmp, ls[[name]]))
  }
  out
}

# find rows of 'x' matching columns passed via 'ls' and '...'
# similar to 'find_elements' but for matrix like objects
find_rows <- function(x, ..., ls = list(), fun = '%in%') {
  x <- as.data.frame(x)
  if (!nrow(x)) {
    return(logical(0))
  }
  out <- rep(TRUE, nrow(x))
  ls <- c(ls, list(...))
  if (!length(ls)) {
    return(out)
  }
  if (is.null(names(ls))) {
    stop("Argument 'ls' must be named.")
  }
  for (name in names(ls)) {
    out <- out & do_call(fun, list(x[[name]], ls[[name]]))
  }
  out
}

# subset 'x' using arguments passed via 'ls' and '...'
subset2 <- function(x, ..., ls = list(), fun = '%in%') {
  x[find_rows(x, ..., ls = ls, fun = fun), , drop = FALSE]
}

# not-in operator
"%notin%" <- Negate("%in%")

# convert array to list of elements with reduced dimension
# @param x an arrary of dimension d
# @return a list of arrays of dimension d-1
array2list <- function(x) {
  if (is.null(dim(x))) {
    return(as.list(x))
  }
  ndim <- length(dim(x))
  out <- list(length = dim(x)[ndim])
  ind <- collapse(rep(",", ndim - 1))
  for (i in seq_len(dim(x)[ndim])) {
    out[[i]] <- eval2(paste0("x[", ind, i, "]"))
    if (length(dim(x)) > 2) {
      # avoid accidental dropping of other dimensions
      dim(out[[i]]) <- dim(x)[-ndim]
    }
  }
  names(out) <- dimnames(x)[[ndim]]
  out
}

# move elements to the start of a named object
move2start <- function(x, first) {
  x[c(first, setdiff(names(x), first))]
}

# move elements to the end of a named object
move2end <- function(x, last) {
  x[c(setdiff(names(x), last), last)]
}

# wrapper around replicate but without simplifying
repl <- function(expr, n) {
  replicate(n, expr, simplify = FALSE)
}

# find the first element in A that is greater than target
# @param A a matrix
# @param target a vector of length nrow(A)
# @param i column of A being checked first
# @return a vector of the same length as target containing the
#   column ids where A[,i] was first greater than target
first_greater <- function(A, target, i = 1) {
  ifelse(target <= A[, i] | ncol(A) == i, i, first_greater(A, target, i + 1))
}

# check if an object is NULL
isNULL <- function(x) {
  is.null(x) || ifelse(is.vector(x), all(sapply(x, is.null)), FALSE)
}

# recursively removes NULL entries from an object
rmNULL <- function(x, recursive = TRUE) {
  x <- Filter(Negate(isNULL), x)
  if (recursive) {
    x <- lapply(x, function(x) if (is.list(x)) rmNULL(x) else x)
  }
  x
}

# find the first argument that is not NULL
first_not_null <- function(...) {
  dots <- list(...)
  out <- NULL
  i <- 1L
  while (isNULL(out) && i <= length(dots)) {
    if (!isNULL(dots[[i]])) {
      out <- dots[[i]]
    }
    i <- i + 1L
  }
  out
}

is_atomic_or_null <- function(x) {
  is.atomic(x) || is.null(x)
}

isNA <- function(x) {
  length(x) == 1L && is.na(x)
}

is_equal <- function(x, y, check.attributes = FALSE, ...) {
  isTRUE(all.equal(x, y, check.attributes = check.attributes, ...))
}

# extract factor levels from an arbitrary variable
extract_levels <- function(x) {
  # do not check for NAs according to #1355
  if (!is.factor(x)) {
    x <- factor(x)
  }
  levels(x)
}

# check if 'x' will behave like a factor in design matrices
is_like_factor <- function(x) {
  is.factor(x) || is.character(x) || is.logical(x)
}

# as.factor but allows to pass levels
as_factor <- function(x, levels = NULL) {
  if (is.null(levels)) {
    out <- as.factor(x)
  } else {
    out <- factor(x, levels = levels)
  }
  out
}

# coerce 'x' to a single logical value
as_one_logical <- function(x, allow_na = FALSE) {
  s <- substitute(x)
  x <- as.logical(x)
  if (length(x) != 1L || anyNA(x) && !allow_na) {
    s <- deparse0(s, max_char = 100L)
    stop2("Cannot coerce '", s, "' to a single logical value.")
  }
  x
}

# coerce 'x' to a single integer value
as_one_integer <- function(x, allow_na = FALSE) {
  s <- substitute(x)
  x <- SW(as.integer(x))
  if (length(x) != 1L || anyNA(x) && !allow_na) {
    s <- deparse0(s, max_char = 100L)
    stop2("Cannot coerce '", s, "' to a single integer value.")
  }
  x
}

# coerce 'x' to a single numeric value
as_one_numeric <- function(x, allow_na = FALSE) {
  s <- substitute(x)
  x <- SW(as.numeric(x))
  if (length(x) != 1L || anyNA(x) && !allow_na) {
    s <- deparse0(s, max_char = 100L)
    stop2("Cannot coerce '", s, "' to a single numeric value.")
  }
  x
}

# coerce 'x' to a single character string
as_one_character <- function(x, allow_na = FALSE) {
  s <- substitute(x)
  x <- as.character(x)
  if (length(x) != 1L || anyNA(x) && !allow_na) {
    s <- deparse0(s, max_char = 100L)
    stop2("Cannot coerce '", s, "' to a single character value.")
  }
  x
}

# coerce 'x' to a single character variable name
as_one_variable <- function(x, allow_na = TRUE) {
  x <- as_one_character(x)
  if (x == "NA" && allow_na) {
    return(x)
  }
  if (!nzchar(x) || !is_equal(x, all_vars(x))) {
    stop2("Cannot coerce '", x, "' to a single variable name.")
  }
  x
}

has_rows <- function(x) {
  isTRUE(nrow(x) > 0L)
}

has_cols <- function(x) {
  isTRUE(ncol(x) > 0L)
}

# expand arguments to the same length
# @param ... arguments to expand
# @param length optional expansion length
#   otherwise taken to be the largest supplied length
# @return a data.frame with one variable per element in '...'
expand <- function(..., dots = list(), length = NULL) {
  dots <- c(dots, list(...))
  max_dim <- NULL
  if (is.null(length)) {
    lengths <- lengths(dots)
    length <- max(lengths)
    max_dim <- dim(dots[[match(length, lengths)]])
  }
  out <- as.data.frame(lapply(dots, rep, length.out = length))
  structure(out, max_dim = max_dim)
}

# structure but ignore NULL
structure_not_null <- function(.Data, ...) {
  if (!is.null(.Data)) {
    .Data <- structure(.Data, ...)
  }
  .Data
}

# remove specified attributes
rm_attr <- function(x, attr) {
  attributes(x)[attr] <- NULL
  x
}

# unidimensional subsetting while keeping attributes
subset_keep_attr <- function(x, y) {
  att <- attributes(x)
  x <- x[y]
  att$names <- names(x)
  attributes(x) <- att
  x
}

'%||%' <- function(x, y) {
  if (is.null(x)) x <- y
  x
}

# check if 'x' is a whole number (integer)
is_wholenumber <- function(x, tol = .Machine$double.eps) {
  if (is.numeric(x)) {
    out <- abs(x - round(x)) < tol
  } else {
    out <- rep(FALSE, length(x))
  }
  dim(out) <- dim(x)
  out
}

# helper function to check symmetry of a matrix
is_symmetric <- function(x, tol = sqrt(.Machine$double.eps)) {
  isSymmetric(x, tol = tol, check.attributes = FALSE)
}

# unlist lapply output
ulapply <- function(X, FUN, ..., recursive = TRUE, use.names = TRUE) {
  unlist(lapply(X, FUN, ...), recursive, use.names)
}

# rbind lapply output
rblapply <- function(X, FUN, ...) {
  do.call(rbind, lapply(X, FUN, ...))
}

# cbind lapply output
cblapply <- function(X, FUN, ...) {
  do.call(cbind, lapply(X, FUN, ...))
}

# parallel lapply sensitive to the operating system
# args:
#  .psock: use a PSOCK cluster? Default is TRUE until
#.    the zombie worker issue #1658 has been fully resolved
plapply <- function(X, FUN, .cores = 1, .psock = TRUE, ...) {
  if (.cores == 1) {
    out <- lapply(X, FUN, ...)
  } else {
    if (!os_is_windows() && !.psock) {
      out <- parallel::mclapply(X = X, FUN = FUN, mc.cores = .cores, ...)
    } else {
      cl <- parallel::makePSOCKcluster(.cores)
      on.exit(parallel::stopCluster(cl))
      out <- parallel::parLapply(cl = cl, X = X, fun = FUN, ...)
    }
    # The version below was suggested to prevent the spawning of zombies
    # but it does not always succeed in that. It also seems to cause
    # other issues as discussed in #1658, so commented out for now.
    # cl_type <- ifelse(os_is_windows(), "PSOCK", "FORK")
    # cl <- parallel::makeCluster(.cores, type = cl_type)
    # # Register a cleanup for the cluster in case the function fails
    # # Need to wrap in a tryCatch to avoid error if cluster is already stopped
    # on.exit(tryCatch(
    #   { parallel::stopCluster(cl) },
    #   error = function(e) invisible(NULL)
    # ))
    # out <- parallel::parLapply(cl = cl, X = X, fun = FUN, ...)
    # parallel::stopCluster(cl)
  }
  out
}

# extract objects stored in each element of a list
# @param x a list-like object
# @param name name of the object to extract
from_list <- function(x, name, ...) {
  lapply(x, "[[", name, ...)
}

# unlist from_list output
ufrom_list <- function(x, name, ..., recursive = TRUE, use.names = TRUE) {
  unlist(from_list(x, name, ...), recursive, use.names)
}

# check if the operating system is Windows
os_is_windows <- function() {
  isTRUE(Sys.info()[['sysname']] == "Windows")
}

# find variables in a character string or expression
all_vars <- function(expr, ...) {
  if (is.character(expr)) {
    expr <- str2expression(expr)
  }
  all.vars(expr, ...)
}

# reimplemented for older R versions
# see ?parse in R 3.6 or higher
str2expression <- function(x) {
  parse(text = x, keep.source = FALSE)
}

# reimplemented for older R versions
# see ?parse in R 3.6 or higher
str2lang <- function(x) {
  str2expression(x)[[1]]
}

# append list(...) to x
lc <- function(x, ...) {
  dots <- rmNULL(list(...), recursive = FALSE)
  c(x, dots)
}

'c<-' <- function(x, value) {
  c(x, value)
}

'lc<-' <- function(x, value) {
  lc(x, value)
}

collapse <- function(..., sep = "") {
  paste(..., sep = sep, collapse = "")
}

collapse_comma <- function(...) {
  paste0("'", ..., "'", collapse = ", ")
}

# add characters to an existing string
'str_add<-' <- function(x, start = FALSE, value) {
  if (start) paste0(value, x) else paste0(x, value)
}

# add list of characters to an existing list
'str_add_list<-' <- function(x, start = FALSE, value) {
  stopifnot(is.list(x), is.list(value))
  out <- if (start) list(value, x) else list(x, value)
  collapse_lists(ls = out)
}

# type-stable if clause for strings with default else output
str_if <- function(cond, yes, no = "") {
  cond <- as_one_logical(cond)
  if (cond) as.character(yes) else as.character(no)
}

# select elements which match a regex pattern
str_subset <- function(x, pattern, ...) {
  x[grepl(pattern, x, ...)]
}

# similar to glue::glue but specialized for generating Stan code
glue <- function(..., sep = "", collapse = NULL, envir = parent.frame(),
                 open = "{", close = "}", na = "NA") {
  dots <- list(...)
  dots <- dots[lengths(dots) > 0L]
  args <- list(
    .x = NULL, .sep = sep, .envir = envir, .open = open,
    .close = close, .na = na, .trim = FALSE,
    .transformer = zero_length_transformer
  )
  out <- do.call(glue::glue_data, c(dots, args))
  if (!is.null(collapse)) {
    collapse <- as_one_character(collapse)
    out <- paste0(out, collapse = collapse)
  }
  out
}

# used in 'glue' to handle zero-length inputs
zero_length_transformer <- function(text, envir) {
  out <- glue::identity_transformer(text, envir)
  if (!length(out)) {
    out <- ""
  }
  out
}

# collapse strings evaluated with glue
cglue <- function(..., envir = parent.frame()) {
  glue(..., envir = envir, collapse = "")
}

# check if a certain package is installed
# @param package package name
# @param version optional minimal version number to require
require_package <- function(package, version = NULL) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop2("Please install the '", package, "' package.")
  }
  if (!is.null(version)) {
    version <- as.package_version(version)
    if (utils::packageVersion(package) < version) {
      stop2("Please install package '", package,
            "' version ", version, " or higher.")
    }
  }
  invisible(TRUE)
}

# rename specified patterns in a character vector
# @param x a character vector to be renamed
# @param pattern the regular expressions in x to be replaced
# @param replacement the replacements
# @param fixed same as for 'gsub'
# @param check_dup: logical; check for duplications in x after renaming
# @param ... passed to 'gsub'
# @return renamed character vector of the same length as x
rename <- function(x, pattern = NULL, replacement = NULL,
                   fixed = TRUE, check_dup = FALSE, ...) {
  pattern <- as.character(pattern)
  replacement <- as.character(replacement)
  if (!length(pattern) && !length(replacement)) {
    # default renaming to avoid special characters in coeffcient names
    pattern <- c(
      " ", "(", ")", "[", "]", ",", "\"", "'",
      "?", "+", "-", "*", "/", "^", "=", "$"
    )
    replacement <- c(rep("", 9), "P", "M", "MU", "D", "E", "EQ", "USD")
  }
  if (length(replacement) == 1L) {
    replacement <- rep(replacement, length(pattern))
  }
  stopifnot(length(pattern) == length(replacement))
  # avoid zero-length pattern error
  has_chars <- nzchar(pattern)
  pattern <- pattern[has_chars]
  replacement <- replacement[has_chars]
  out <- x
  for (i in seq_along(pattern)) {
    out <- gsub(pattern[i], replacement[i], out, fixed = fixed, ...)
  }
  dup <- duplicated(out)
  if (check_dup && any(dup)) {
    dup <- x[out %in% out[dup]]
    stop2("Internal renaming led to duplicated names. ",
          "Consider renaming your variables to have different suffixes.\n",
          "Occured for: ", collapse_comma(dup))
  }
  out
}

# collapse strings having the same name in different lists
# @param ... named lists
# @param ls a list of named lists
# @param a named list containing the collapsed strings
collapse_lists <- function(..., ls = list()) {
  ls <- c(list(...), ls)
  elements <- unique(unlist(lapply(ls, names)))
  args <- c(FUN = collapse, lapply(ls, "[", elements), SIMPLIFY = FALSE)
  out <- do.call(mapply, args)
  names(out) <- elements
  out
}

# create a named list using object names
nlist <- function(...) {
  m <- match.call()
  dots <- list(...)
  no_names <- is.null(names(dots))
  has_name <- if (no_names) FALSE else nzchar(names(dots))
  if (all(has_name)) return(dots)
  nms <- as.character(m)[-1]
  if (no_names) {
    names(dots) <- nms
  } else {
    names(dots)[!has_name] <- nms[!has_name]
  }
  dots
}

# initialize a named list
# @param names names of the elements
# @param values optional values of the elements
named_list <- function(names, values = NULL) {
  if (!is.null(values)) {
    if (length(values) <= 1L) {
      values <- replicate(length(names), values)
    }
    values <- as.list(values)
    stopifnot(length(values) == length(names))
  } else {
    values <- vector("list", length(names))
  }
  setNames(values, names)
}

# is an object named?
is_named <- function(x) {
  names <- names(x)
  if (is.null(names)) {
    return(FALSE)
  }
  if (any(!nzchar(names) | is.na(names))) {
    return(FALSE)
  }
  TRUE
}

#' Execute a Function Call
#'
#' Execute a function call similar to \code{\link{do.call}}, but without
#' deparsing function arguments. For large number of arguments (i.e., more
#' than a few thousand) this function currently is somewhat inefficient
#' and should be used with care in this case.
#'
#' @param what Either a function or a non-empty character string naming the
#'   function to be called.
#' @param args A list of arguments to the function call. The names attribute of
#'   \code{args} gives the argument names.
#' @param pkg Optional name of the package in which to search for the
#'   function if \code{what} is a character string.
#' @param envir An environment within which to evaluate the call.
#'
#' @return The result of the (evaluated) function call.
#'
#' @keywords internal
#' @export
do_call <- function(what, args, pkg = NULL, envir = parent.frame()) {
  call <- ""
  if (length(args)) {
    if (!is.list(args)) {
      stop2("'args' must be a list.")
    }
    fun_args <- names(args)
    if (is.null(fun_args)) {
      fun_args <- rep("", length(args))
    } else {
      nzc <- nzchar(fun_args)
      fun_args[nzc] <- paste0("`", fun_args[nzc], "` = ")
    }
    names(args) <- paste0(".x", seq_along(args))
    call <- paste0(fun_args, names(args), collapse = ",")
  } else {
    args <- list()
  }
  if (is.function(what)) {
    args$.fun <- what
    what <- ".fun"
  } else {
    what <- paste0("`", as_one_character(what), "`")
    if (!is.null(pkg)) {
      what <- paste0(as_one_character(pkg), "::", what)
    }
  }
  call <- paste0(what, "(", call, ")")
  eval2(call, envir = args, enclos = envir)
}

# create an empty data frame
empty_data_frame <- function() {
  as.data.frame(matrix(nrow = 0, ncol = 0))
}

# replace elements in x with elements in value
# @param x named list-like object
# @param value another named list-like object
# @param dont_replace names of elements that cannot be replaced
'replace_args<-' <- function(x, dont_replace = NULL, value) {
  value_name <- deparse0(substitute(value), max_char = 100L)
  value <- as.list(value)
  if (length(value) && is.null(names(value))) {
    stop2("Argument '", value_name, "' must be named.")
  }
  invalid <- names(value)[names(value) %in% dont_replace]
  if (length(invalid)) {
    invalid <- collapse_comma(invalid)
    stop2("Argument(s) ", invalid, " cannot be replaced.")
  }
  x[names(value)] <- value
  x
}

# deparse0 'x' if it is no string
deparse_no_string <- function(x) {
  if (!is.character(x)) {
    x <- deparse0(x)
  }
  x
}

# combine deparse lines into one string
# since R 4.0 we also have base::deparse1 for this purpose
deparse0 <- function(x, max_char = NULL, ...) {
  out <- collapse(deparse(x, ...))
  if (isTRUE(max_char > 0)) {
    out <- substr(out, 1L, max_char)
  }
  out
}

# like 'eval' but parses characters before evaluation
eval2 <- function(expr, envir = parent.frame(), ...) {
  if (is.character(expr)) {
    expr <- str2expression(expr)
  }
  eval(expr, envir, ...)
}

# evaluate an expression without printing output or messages
# @param expr expression to be evaluated
# @param type type of output to be suppressed (see ?sink)
# @param try wrap evaluation of expr in 'try' and
#   not suppress outputs if evaluation fails?
# @param silent actually evaluate silently?
eval_silent <- function(expr, type = "output", try = FALSE,
                        silent = TRUE, ...) {
  try <- as_one_logical(try)
  silent <- as_one_logical(silent)
  type <- match.arg(type, c("output", "message"))
  expr <- substitute(expr)
  envir <- parent.frame()
  if (silent) {
    if (try && type == "message") {
      try_out <- try(utils::capture.output(
        out <- eval(expr, envir), type = type, ...
      ))
      if (is_try_error(try_out)) {
        # try again without suppressing error messages
        out <- eval(expr, envir)
      }
    } else {
      utils::capture.output(out <- eval(expr, envir), type = type, ...)
    }
  } else {
    out <- eval(expr, envir)
  }
  out
}

# find the name that 'x' had in a specific environment
substitute_name <- function(x, envir = parent.frame(), nchar = 50) {
  out <- substitute(x)
  out <- eval2(paste0("substitute(", out, ")"), envir = envir)
  if (missing(out)) {
    return(NULL)
  }
  substr(collapse(deparse(out)), 1, nchar)
}

# recursive sorting of dependencies
# @param x named list of dependencies per element
# @param sorted already sorted element names
# @return a vector of sorted element names
sort_dependencies <- function(x, sorted = NULL) {
  if (!length(x)) {
    return(NULL)
  }
  if (length(names(x)) != length(x)) {
    stop2("Argument 'x' must be named.")
  }
  take <- !ulapply(x, function(dep) any(!dep %in% sorted))
  new <- setdiff(names(x)[take], sorted)
  out <- union(sorted, new)
  if (length(new)) {
    out <- union(out, sort_dependencies(x, sorted = out))
  } else if (!all(names(x) %in% out)) {
    stop2("Cannot handle circular dependency structures.")
  }
  out
}

# older version of stop2 was preserved here to compare
# in case of any error
# stop2 <- function(...) {
#   stop(..., call. = FALSE)
# }

stop2 <- function(message = "", ..., .subclass = NULL, call = rlang::caller_call()) {
  rlang::abort(
    message = glue::glue(message, ...),
    .subclass = c(.subclass, "brms_error"),
    call = call
  )
}

warning2 <- function(...) {
  warning(..., call. = FALSE)
}

# get first occurrence of 'x' in '...' objects
# @param x The name of the required element
# @param ... named R objects that may contain 'x'
get_arg <- function(x, ...) {
  dots <- list(...)
  i <- 1
  out <- NULL
  while (i <= length(dots) && is.null(out)) {
    if (!is.null(dots[[i]][[x]])) {
      out <- dots[[i]][[x]]
    } else {
      i <- i + 1
    }
  }
  out
}

SW <- function(expr) {
  base::suppressWarnings(expr)
}

# get pattern matches in text as vector
# @param simplify return an atomic vector of matches?
# @param first only return the first match in each string?
# @return character vector containing matches
get_matches <- function(pattern, text, simplify = TRUE,
                        first = FALSE, ...) {
  x <- regmatches(text, gregexpr(pattern, text, ...))
  if (first) {
    x <- lapply(x, function(t) if (length(t)) t[1] else t)
  }
  if (simplify) {
    if (first) {
      x <- lapply(x, function(t) if (length(t)) t else "")
    }
    x <- unlist(x)
  }
  x
}

# find matches in the parse tree of an expression
# @param pattern pattern to be matched
# @param expr expression to be searched in
# @return character vector containing matches
get_matches_expr <- function(pattern, expr, ...) {
  if (is.character(expr)) {
    expr <- str2expression(expr)
  }
  out <- NULL
  for (i in seq_along(expr)) {
    sexpr <- try(expr[[i]], silent = TRUE)
    if (!is_try_error(sexpr)) {
      sexpr_char <- deparse0(sexpr)
      out <- c(out, get_matches(pattern, sexpr_char, ...))
    }
    if (is.call(sexpr) || is.expression(sexpr)) {
      out <- c(out, get_matches_expr(pattern, sexpr, ...))
    }
  }
  trim_wsp(unique(out))
}

# like 'grepl' but handles (parse trees of) expressions
grepl_expr <- function(pattern, expr, ...) {
  as.logical(ulapply(expr, function(e)
    length(get_matches_expr(pattern, e, ...)) > 0L))
}

# combine character vectors into a joint regular 'or' expression
# @param x a character vector
# @param escape escape all special characters in 'x'?
regex_or <- function(x, escape = FALSE) {
  if (escape) {
    x <- escape_all(x)
  }
  paste0("(", paste0("(", x, ")", collapse = "|"), ")")
}

# escape dots in character strings
escape_dot <- function(x) {
  gsub(".", "\\.", x, fixed = TRUE)
}

# escape all special characters in character strings
escape_all <- function(x) {
  specials <- c(".", "*", "+", "?", "^", "$", "(", ")", "[", "]", "|")
  for (s in specials) {
    x <- gsub(s, paste0("\\", s), x, fixed = TRUE)
  }
  x
}

# add an underscore to non-empty character strings
# @param x a character vector
# @param pos position of the underscore
usc <- function(x, pos = c("prefix", "suffix")) {
  pos <- match.arg(pos)
  x <- as.character(x)
  if (!length(x)) x <- ""
  if (pos == "prefix") {
    x <- ifelse(nzchar(x), paste0("_", x), "")
  } else {
    x <- ifelse(nzchar(x), paste0(x, "_"), "")
  }
  x
}

# round using the largest remainder method
round_largest_remainder <- function(x) {
  x <- as.numeric(x)
  total <- round(sum(x))
  out <- floor(x)
  diff <- x - out
  J <- order(diff, decreasing = TRUE)
  I <- seq_len(total - floor(sum(out)))
  out[J[I]] <- out[J[I]] + 1
  out
}

# add leading and trailing white spaces
# @param x object accepted by paste
# @param nsp number of white spaces to add
wsp <- function(x = "", nsp = 1) {
  sp <- collapse(rep(" ", nsp))
  if (length(x)) {
    out <- ifelse(nzchar(x), paste0(sp, x, sp), sp)
  } else {
    out <- NULL
  }
  out
}

# add white space per line the the strings
# @param x object accepted by paste
# @param nsp number of white spaces to add
wsp_per_line <- function(x, nsp) {
  sp <- collapse(rep(" ", nsp))
  x <- paste0(sp, x)
  x <- gsub("\\n(?=.+)", paste0("\n", sp), x, perl = TRUE)
  x
}

# remove whitespaces in character strings
rm_wsp <- function(x) {
  out <- gsub("[ \t\r\n]+", "", x, perl = TRUE)
  dim(out) <- dim(x)
  out
}

# trim whitespaces in character strings
trim_wsp <- function(x) {
  out <- gsub("[ \t\r\n]+", " ", x, perl = TRUE)
  dim(out) <- dim(x)
  out
}

# limit the number of characters of a vector
# @param x a character vector
# @param chars maximum number of characters to show
# @param lsuffix number of characters to keep at the end of the strings
# @return possible truncated character vector
limit_chars <- function(x, chars = NULL, lsuffix = 4) {
  stopifnot(is.character(x))
  if (!is.null(chars)) {
    chars_x <- nchar(x) - lsuffix
    suffix <- substr(x, chars_x + 1, chars_x + lsuffix)
    x <- substr(x, 1, chars_x)
    x <- ifelse(chars_x <= chars, x, paste0(substr(x, 1, chars - 3), "..."))
    x <- paste0(x, suffix)
  }
  x
}

# ensure that deprecated arguments still work
# @param arg input to the new argument
# @param alias input to the deprecated argument
# @param default the default value of alias
# @param warn should a warning be printed if alias is specified?
use_alias <- function(arg, alias = NULL, default = NULL, warn = TRUE) {
  arg_name <- Reduce(paste, deparse(substitute(arg)))
  alias_name <- Reduce(paste, deparse(substitute(alias)))
  if (!is_equal(alias, default)) {
    arg <- alias
    if (grepl("^dots\\$", alias_name)) {
      alias_name <- gsub("^dots\\$", "", alias_name)
    } else if (grepl("^dots\\[\\[", alias_name)) {
      alias_name <- gsub("^dots\\[\\[\"|\"\\]\\]$", "", alias_name)
    }
    if (warn) {
      warning2("Argument '", alias_name, "' is deprecated. ",
               "Please use argument '", arg_name, "' instead.")
    }
  }
  arg
}

warn_deprecated <- function(new, old = as.character(sys.call(sys.parent()))[1]) {
  msg <- paste0("Function '", old, "' is deprecated.")
  if (!missing(new)) {
    msg <- paste0(msg, " Please use '", new, "' instead.")
  }
  warning2(msg)
  invisible(NULL)
}

# check if x is a try-error resulting from try()
is_try_error <- function(x) {
  inherits(x, "try-error")
}

# check if verbose mode is activated
is_verbose <- function() {
  as_one_logical(getOption("brms.verbose", FALSE))
}

viridis6 <- function() {
  c("#440154", "#414487", "#2A788E", "#22A884", "#7AD151", "#FDE725")
}

expect_match2 <- function(object, regexp, ..., all = TRUE) {
  testthat::expect_match(object, regexp, fixed = TRUE, ..., all = all)
}

# startup messages for brms
.onAttach <- function(libname, pkgname) {
  version <- utils::packageVersion("brms")
  packageStartupMessage(
    "Loading 'brms' package (version ", version, "). Useful instructions\n",
    "can be found by typing help('brms'). A more detailed introduction\n",
    "to the package is available through vignette('brms_overview')."
  )
  invisible(NULL)
}

# code to execute when loading brms
.onLoad <- function(libname, pkgname) {
  # ensure compatibility with older R versions
  backports::import(pkgname)
  # dynamically register the 'recover_data' and 'emm_basis'
  # methods needed by 'emmeans', if that package is installed
  if (requireNamespace("emmeans", quietly = TRUE) &&
      utils::packageVersion("emmeans") >= "1.4.0") {
    emmeans::.emm_register("brmsfit", pkgname)
  }
  invisible(NULL)
}
