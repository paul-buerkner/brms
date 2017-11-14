#' Data for \pkg{brms} Models
#' 
#' Generate data for \pkg{brms} models to be passed to \pkg{Stan}
#'
#' @inheritParams brm
#' @param control A named list currently for internal usage only
#' @param ... Other potential arguments
#' 
#' @aliases brmdata
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
                          prior = NULL, autocor = NULL, 
                          nonlinear = NULL, cov_ranef = NULL, 
                          sample_prior = c("no", "yes", "only"), 
                          knots = NULL, control = list(), ...) {
  # internal control arguments:
  #   is_newdata: is make_standata is called with new data?
  #   not4stan: is make_standata called for use in S3 methods?
  #   save_order: should the initial order of the data be saved?
  #   omit_response: omit checking of the response?
  #   only_response: only compute data related to the response?
  #   ntrials, ncat, Jmo: standata based on the original data
  dots <- list(...)
  not4stan <- isTRUE(control$not4stan)
  is_newdata <- isTRUE(control$is_newdata)
  # use deprecated arguments if specified
  cov_ranef <- use_alias(cov_ranef, dots$cov.ranef, warn = FALSE)
  
  formula <- amend_formula(
    formula, data = data, family = family, 
    autocor = autocor, nonlinear = nonlinear
  )
  bterms <- parse_bf(formula)
  sample_prior <- check_sample_prior(sample_prior)
  check_prior_content(prior, warn = FALSE)
  prior <- check_prior_special(bterms, prior = prior)
  na_action <- if (is_newdata) na.pass else na.omit
  data <- update_data(
    data, bterms = bterms, na.action = na_action, 
    drop.unused.levels = !is_newdata, knots = knots,
    terms_attr = control$terms_attr
  )
  if (has_arma(autocor) || is.cor_bsts(autocor)) {
    # order data in case of autocorrelation models
    data <- order_data(data, bterms = bterms)
  }
  
  out <- list(N = nrow(data))
  check_response <- !isTRUE(control$omit_response)
  # TODO: pass control$trials and control$ncat
  out <- c(out, data_response(bterms, data, check_response))
  only_response <- isTRUE(control$only_response)
  if (!only_response) {
    ranef <- tidy_ranef(
      bterms, data, ncat = control$ncat, 
      old_levels = control$old_levels
    )
    args_eff <- nlist(
      x = bterms, data, prior, ranef, cov_ranef, knots, 
      not4stan, smooths = control$smooths, 
      gps = control$gps, Jmo = control$Jmo, 
      old_locations = control$old_locations
    )
    out <- c(out, do.call(data_effects, args_eff))
  }
  out$prior_only <- as.integer(identical(sample_prior, "only"))
  if (isTRUE(control$save_order)) {
    attr(out, "old_order") <- attr(data, "old_order")
  }
  structure(out, class = "standata")
}
