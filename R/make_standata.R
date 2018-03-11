#' Data for \pkg{brms} Models
#' 
#' Generate data for \pkg{brms} models to be passed to \pkg{Stan}
#'
#' @inheritParams brm
#' @param check_response Logical; check validity of the response?
#' @param only_response Logical; extract data related to the response only?
#' @param control A named list currently for internal usage only
#' @param ... Other potential arguments
#' 
#' @return A named list of objects containing the required data 
#'   to fit a \pkg{brms} model with \pkg{Stan}. 
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' data1 <- make_standata(rating ~ treat + period + carry + (1|subject), 
#'                        data = inhaler, family = "cumulative")
#' names(data1)
#' 
#' data2 <- make_standata(count ~ log_Age_c + log_Base4_c * Trt_c 
#'                        + (1|patient) + (1|visit), 
#'                        data = epilepsy, family = "poisson")
#' names(data2)
#'          
#' @export
make_standata <- function(formula, data, family = gaussian(), 
                          prior = NULL, autocor = NULL, cov_ranef = NULL,
                          sample_prior = c("no", "yes", "only"), 
                          stan_vars = NULL, knots = NULL, 
                          check_response = TRUE, only_response = FALSE, 
                          control = list(), ...) {
  # internal control arguments:
  #   new: is make_standata is called with new data?
  #   not4stan: is make_standata called for use in S3 methods?
  #   save_order: should the initial order of the data be saved?
  #   old_sdata: list of stan data computed from the orginal data
  #   terms_attr: list of attributes of the original model.frame
  dots <- list(...)
  # some input checks
  if (is.brmsfit(formula)) {
    stop2("Use 'standata' to extract Stan data from 'brmsfit' objects.")
  }
  check_response <- as_one_logical(check_response)
  only_response <- as_one_logical(only_response)
  not4stan <- isTRUE(control$not4stan)
  new <- isTRUE(control$new)
  formula <- validate_formula(
    formula, data = data, family = family, autocor = autocor
  )
  bterms <- parse_bf(formula)
  sample_prior <- check_sample_prior(sample_prior)
  check_prior_content(prior, warn = FALSE)
  prior <- check_prior_special(
    prior, bterms = bterms, data = data, 
    check_nlpar_prior = FALSE
  )
  na_action <- if (new) na.pass else na.omit2
  data <- update_data(
    data, bterms = bterms, na.action = na_action, 
    drop.unused.levels = !new, knots = knots,
    terms_attr = control$terms_attr
  )
  if (has_arma(autocor) || is.cor_bsts(autocor)) {
    # order data in case of autocorrelation models
    data <- order_data(data, bterms = bterms)
  }
  
  out <- c(
    list(N = nrow(data)), 
    data_response(
      bterms, data, check_response = check_response,
      not4stan = not4stan, new = new, 
      old_sdata = control$old_sdata
    )
  )
  if (!only_response) {
    ranef <- tidy_ranef(
      bterms, data, old_levels = control$old_levels,
      old_sdata = control$old_sdata  
    )
    meef <- tidy_meef(bterms, data)
    args_eff <- nlist(
      x = bterms, data, prior, ranef, meef, cov_ranef, 
      knots, not4stan, old_sdata = control$old_sdata
    )
    out <- c(out, do.call(data_effects, args_eff))
  }
  out$prior_only <- as.integer(identical(sample_prior, "only"))
  stan_vars <- validate_stanvars(stan_vars)
  if (is.stanvars(stan_vars)) {
    inv_names <- intersect(names(stan_vars), names(out))
    if (length(inv_names)) {
      stop2("Cannot overwrite existing variables: ", 
            collapse_comma(inv_names))
    }
    out[names(stan_vars)] <- lapply(stan_vars, "[[", "sdata")
  }
  if (isTRUE(control$save_order)) {
    attr(out, "old_order") <- attr(data, "old_order")
  }
  structure(out, class = "standata")
}
