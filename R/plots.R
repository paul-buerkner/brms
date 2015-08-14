# make a trace and density plot for one parameter
td_plot <- function(par, x) {
  if (!is.character(par) || length(par) != 1)
    stop("par must be a character string")
  if (!is.data.frame(x))
    stop("x must be a data.frame")
  names(x)[match(par, names(x))] <- "value" 
  trace <- ggplot(x, aes_string(x = "iter", y = "value", group = "chains", colour = "chains")) +
    geom_line(alpha = 0.7) + xlab("") + ggtitle(par) + 
    theme(legend.position = "none",
          plot.title = element_text(size=15, vjust=1),
          plot.margin = grid::unit(c(0.2,0,-0.8,0), "lines"))
  density <- ggplot(x, aes_string(x = "value")) + 
    geom_density(aes_string(fill = "chains"), alpha = 0.5) + xlab("") + ggtitle(par) + 
    theme(plot.title = element_text(size=15, vjust=1),
          plot.margin = grid::unit(c(0.2,0,-0.8,0), "lines"))
  return(gridExtra::arrangeGrob(trace, density, ncol = 2, nrow = 1))
}

#' Trace and density plots for MCMC samples
#' 
#' @param x An object of class \code{brmsfit}.
#' @param parameters Name of the parameters to plot, as given by a character vector or a regular expression.
#'   By default, all parameters except for random effects, posterior predictives, and log likelihood values are plotted. 
#' @param N The number of parameters plotted per page.
#' @param ask logical; Indicates if the user is prompted before a new page is plotted.   
#' @param ... Currently ignored.
#' 
#' @return NULL
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' \dontrun{ 
#' fit_e <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1|patient) + (1|visit), 
#'              data = epilepsy, family = "poisson")
#' ## plot fixed effects as well as standard devations of the random effects
#' plot(fit_e)
#' ## plot fixed effects only and combine the chains into one posterior
#' plot(fit_e, parameters = "^b_", combine = TRUE) 
#' }
#' 
#' @method plot brmsfit
#' @import ggplot2
#' @export
plot.brmsfit <- function(x, parameters = NA, N = 5, ask = TRUE, ...) {
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  if (is.na(parameters)) 
    parameters <- c("^b_", "^sd_", "^cor_", "^sigma", "^rescor", "^nu$", 
                    "^shape$", "^delta$", "^ar", "^ma")
  samples <- posterior.samples(x, parameters = parameters, add.chains = TRUE)
  pars <- names(samples)[which(!names(samples) %in% c("chains", "iter"))] 
  
  default.ask <- grDevices::devAskNewPage()
  grDevices::devAskNewPage(ask = FALSE)
  for (i in 1:ceiling(length(pars)/N)) {
    plots <- lapply(pars[((i-1)*N+1):min(i*N,length(pars))], td_plot, x = samples)
    plot(gridExtra::arrangeGrob(grobs = plots, nrow = length(plots), ncol = 1))
    if (i == 1) grDevices::devAskNewPage(ask = ask)
  }
  grDevices::devAskNewPage(default.ask)
}