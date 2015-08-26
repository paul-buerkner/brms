#' Parameters of interest for \code{brms} models (deprecated)
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
#' brm.pars(rating ~ treat + period + carry + (1|subject),
#'          data = inhaler, family = "cumulative")
#
#  brm.pars(count ~ log_Age_c + log_Base4_c * Trt_c + (1|patient) + (1|visit),
#         data = epilepsy, family = c("poisson", "log"))
#'               
#' @export
brm.pars = function(formula, data = NULL, family = "gaussian", autocor = NULL, partial = NULL,
                    threshold = "flexible", ranef = TRUE) {
  family <- check_family(family[1])
  if (is.null(autocor)) autocor <- cor.arma()
  if (!is(autocor,"cor.brms")) stop("cor must be of class cor.brms")
  ee <- extract.effects(formula = formula, family = family, partial = partial)
  data <- update_data(data, family = family, effects = ee)
  
  is.linear <- family %in% c("gaussian", "student", "cauchy")
  is.ordinal <- family  %in% c("cumulative","cratio","sratio","acat")
  is.skew <- family %in% c("gamma", "weibull", "exponential")
  if (!(is.linear || is.ordinal || is.skew || family %in% 
        c("poisson", "negbinomial", "geometric", "binomial","bernoulli", "categorical")))
    stop(paste(family,"is not a valid family"))
  
  f <- colnames(brm.model.matrix(ee$fixed, data, rm.int = is.ordinal))
  r <- lapply(lapply(ee$random, brm.model.matrix, data = data), colnames)
  p <- colnames(brm.model.matrix(partial, data, rm.int = TRUE))
  out <- NULL
  if (is.ordinal && threshold == "flexible") out <- c(out, "b_Intercept")
  if (is.ordinal && threshold == "equidistant") out <- c(out, "b_Intercept1", "delta")
  if (length(f) && family != "categorical") out <- c(out, "b")
  if (is.ordinal && length(p) || family == "categorical") out <- c(out, "bp")
  if (is.linear && !is(ee$se,"formula") && length(ee$response) == 1) out <- c(out, "sigma")
  if (family == "gaussian" && length(ee$response) > 1) out <- c(out, "sigma", "rescor")
  if (family == "student") out <- c(out,"nu")
  if (family %in% c("gamma", "weibull", "negbinomial")) out <- c(out,"shape")
  if (autocor$p > 0) out <- c(out,"ar")
  if (autocor$q > 0) out <- c(out,"ma")
  if (length(ee$group)) {
    out <- c(out, paste0("sd_",ee$group))
    out <- c(out, unlist(lapply(1:length(ee$group), function(i)
      if (length(r[[i]])>1 && ee$cor[[i]]) paste0("cor_",ee$group[[i]]))))
    if (ranef) out <- c(out, paste0("r_",ee$group))
  }
  return(out)
}