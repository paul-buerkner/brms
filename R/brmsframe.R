brmsframe <- function(x, ...) {
  # TODO: store output of standata_basis in brmsframe?
  UseMethod("brmsframe")
}

#' @export
brmsframe.mvbrmsterms <- function(x, data, old_levels = NULL, ...) {
  # this is a univariate model so brmsterms is at the top level
  x$ranef <- tidy_ranef(x, data = data, old_levels = old_levels)
  x$meef <- tidy_meef(x, data = data, old_levels = old_levels)
  for (r in names(x$terms)) {
    x$terms[[r]] <- brmsframe(
      x$terms[[r]], data = data, mv = TRUE,
      old_levels = old_levels, ...
    )
  }
  class(x) <- c("mvbrmsframe", class(x))
  x
}

#' @export
brmsframe.brmsterms <- function(x, data, old_levels = NULL,
                                mv = FALSE, ...) {
  mv <- as_one_logical(mv)
  if (!mv) {
    # this is a univariate model so brmsterms is at the top level
    x$ranef <- tidy_ranef(x, data = data, old_levels = old_levels)
    x$meef <- tidy_meef(x, data = data, old_levels = old_levels)
  }
  data <- subset_data(data, x)
  for (dp in names(x$dpars)) {
    x$dpars[[dp]] <- brmsframe(x$dpars[[dp]], data = data, ...)
  }
  for (nlp in names(x$nlpars)) {
    x$nlpars[[nlp]] <- brmsframe(x$nlpars[[nlp]], data = data, ...)
  }
  x$resp_values <- model.response(model.frame(x$respform, data, na.action = na.pass))
  x$data <- data
  class(x) <- c("brmsframe", class(x))
  x
}

#' @export
brmsframe.btl <- function(x, data, ...) {
  # TODO: rename the tidy_ functions to brmsframe_ functions and move them here?
  # TODO: create a proper fixef data.frame as for the other terms?
  # TODO: store all the data_ functions?
  x$sdata <- list(
    fe = data_fe(x, data),
    sm = data_sm(x, data),
    gp = data_gp(x, data, internal = TRUE),
    offset = data_offset(x, data)
  )
  x$effects <- list(
    fe = colnames(x$sdata$fe$X),
    sp = tidy_spef(x, data),
    cs = colnames(get_model_matrix(x$cs, data = data)),
    gp = tidy_gpef(x, data),
    sm = tidy_smef(x),
    ac = tidy_acef(x)
  )
  class(x) <- c("bfrl", class(x))
  x
}

#' @export
brmsframe.btnl <- function(x, ...) {
  x$acef <- tidy_acef(x)
  class(x) <- c("bfrnl", class(x))
  x
}

#' @export
brmsframe.default <- function(x, ...) {
  brmsframe(brmsterms(x), ...)
}
