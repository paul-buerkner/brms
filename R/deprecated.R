#' Parameters of interest for \code{brms} models (deprecated)
#' 
#' @aliases brm.pars
#' 
#' @inheritParams brm
#' @param ranef logical; indicating if random effects estimates should be returned
#' 
#' @return A vector of character strings specifying parameters of interest for models produced by the \code{brms} package.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @details This function is depricated. Parameters of interest are now chosen by exclusion not by inclusion.
#' 
#' @examples 
#' brmpars(rating ~ treat + period + carry + (1|subject),
#'          data = inhaler, family = "cumulative")
#
#  brmpars(count ~ log_Age_c + log_Base4_c * Trt_c + (1|patient) + (1|visit),
#         data = epilepsy, family = c("poisson", "log"))
#'               
#' @export
brmpars <- function(formula, data = NULL, family = "gaussian", autocor = NULL, partial = NULL,
                    threshold = "flexible", ranef = TRUE) {
  family <- check_family(family[1])
  if (is.null(autocor)) autocor <- cor_arma()
  if (!is(autocor,"cor_brms")) stop("cor must be of class cor_brms")
  ee <- extract_effects(formula = formula, family = family, partial = partial)
  data <- update_data(data, family = family, effects = ee)
  
  is_linear <- family %in% c("gaussian", "student", "cauchy")
  is_ordinal <- family  %in% c("cumulative","cratio","sratio","acat")
  is_skew <- family %in% c("gamma", "weibull", "exponential")
  if (!(is_linear || is_ordinal || is_skew || family %in% 
        c("poisson", "negbinomial", "geometric", "binomial","bernoulli", "categorical")))
    stop(paste(family,"is not a valid family"))
  
  f <- colnames(get_model_matrix(ee$fixed, data, rm_intercept = is_ordinal))
  r <- lapply(lapply(ee$random, get_model_matrix, data = data), colnames)
  p <- colnames(get_model_matrix(partial, data, rm_intercept = TRUE))
  out <- NULL
  if (is_ordinal && threshold == "flexible") out <- c(out, "b_Intercept")
  if (is_ordinal && threshold == "equidistant") out <- c(out, "b_Intercept1", "delta")
  if (length(f) && family != "categorical") out <- c(out, "b")
  if (is_ordinal && length(p) || family == "categorical") out <- c(out, "bp")
  if (is_linear && !is(ee$se,"formula") && length(ee$response) == 1) out <- c(out, "sigma")
  if (family == "gaussian" && length(ee$response) > 1) out <- c(out, "sigma", "rescor")
  if (family == "student") out <- c(out,"nu")
  if (family %in% c("gamma", "weibull", "negbinomial")) out <- c(out, "shape")
  if (autocor$p > 0) out <- c(out, "ar")
  if (autocor$q > 0) out <- c(out, "ma")
  if (length(ee$group)) {
    out <- c(out, paste0("sd_",ee$group))
    out <- c(out, unlist(lapply(1:length(ee$group), function(i)
      if (length(r[[i]])>1 && ee$cor[[i]]) paste0("cor_",ee$group[[i]]))))
    if (ranef) out <- c(out, paste0("r_",ee$group))
  }
  return(out)
}

#' @export
brm.pars <- function(formula, data = NULL, family = "gaussian", autocor = NULL, partial = NULL,
                     threshold = "flexible", ranef = TRUE) {
  # deprecated alias of brm.pars
  brmpars(formula = formula, data = data, family = family, autocor = autocor, partial = partial,
          threshold = threshold, ranef = ranef)
}