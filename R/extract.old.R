#' Trace and density plots for MCMC samples (deprecated)
#' 
#' @param fit An object of class \code{stanfit} or an object that can be coerced to an 
#'  \code{mcmc.list} via \code{as.mcmc} of package \code{coda}.
#' @param family Name of the family of parameters to plot, as given by a character vector or a regular expression. A family of parameters is considered to 
#'   be any group of parameters with the same name but different numerical value between square brackets (as beta[1], beta[2], etc). 
#'   By default, all parameters except for random effects of each level are plotted. 
#'   For details see the documentation of \code{ggmcmc}.
#'   
#' @details This function is deprecated and will soon be defunct. For models created with brms 0.2.0 or higher, 
#' we recommend using the S3 plot method for class \code{brmsfit} (see also \code{\link[brms:plot.brmsfit]{plot.brmsfit}}).
#' 
#' @return NULL
#' 
#' @import ggplot2
#' @export
brm.plot <- function(fit, family = "^[^(r_)]") {
  pfit <- ggmcmc::ggs(ifelse(is(fit,"stanfit"), list(fit), list(coda::as.mcmc(fit)))[[1]])   
  gridExtra::grid.arrange(ggmcmc::ggs_traceplot(pfit, family = family) + 
                            ggplot2::theme(legend.position = "none"), 
                          ggmcmc::ggs_density(pfit, family = family), ncol = 2, nrow = 1)
}