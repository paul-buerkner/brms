#' Trace and Density Plots for MCMC Samples
#' 
#' @param x An object of class \code{brmsfit}.
#' @param pars Names of the parameters to plot, as given by a character vector 
#'   or a regular expression. By default, all parameters except 
#'   for group-level and smooth effects are plotted. 
#' @param combo A character vector with at least two elements. 
#'   Each element of \code{combo} corresponds to a column in the resulting 
#'   graphic and should be the name of one of the available 
#'   \code{\link[bayesplot:MCMC-overview]{MCMC}} functions 
#'   (omitting the \code{mcmc_} prefix).
#' @param N The number of parameters plotted per page.
#' @param theme A \code{\link[ggplot2:theme]{theme}} object 
#'   modifying the appearance of the plots. 
#'   For some basic themes see \code{\link[ggplot2:ggtheme]{ggtheme}}
#'   and \code{\link[bayesplot:theme_default]{theme_default}}.
#' @param fixed Indicates whether parameter names 
#'   should be matched exactly (\code{TRUE}) or treated as
#'   regular expressions (\code{FALSE}). Default is \code{FALSE}.
#' @param exact_match Deprecated alias of argument \code{fixed}.
#' @param plot Logical; indicates if plots should be
#'   plotted directly in the active graphic device.
#'   Defaults to \code{TRUE}.
#' @param ask Logical; indicates if the user is prompted 
#'   before a new page is plotted. 
#'   Only used if \code{plot} is \code{TRUE}.
#' @param newpage Logical; indicates if the first set of plots
#'   should be plotted to a new page. 
#'   Only used if \code{plot} is \code{TRUE}.
#' @param ... Further arguments passed to 
#'   \code{\link[bayesplot:MCMC-combos]{mcmc_combo}}.
#' 
#' @return An invisible list of 
#'   \code{\link[gtable:gtable]{gtable}} objects.
#' 
#' @examples
#' \dontrun{ 
#' fit <- brm(count ~ zAge + zBase * Trt 
#'            + (1|patient) + (1|visit), 
#'            data = epilepsy, family = "poisson")
#' plot(fit)
#' ## plot population-level effects only
#' plot(fit, pars = "^b_") 
#' }
#' 
#' @method plot brmsfit
#' @import ggplot2
#' @importFrom graphics plot
#' @importFrom grDevices devAskNewPage
#' @export
plot.brmsfit <- function(x, pars = NA, combo = c("dens", "trace"), 
                         N = 5, fixed = FALSE, exact_match = FALSE,
                         theme = NULL, plot = TRUE, ask = TRUE, 
                         newpage = TRUE, ...) {
  contains_samples(x)
  if (!is_wholenumber(N) || N < 1) {
    stop2("Argument 'N' must be a positive integer.")
  }
  fixed <- check_deprecated_fixed(fixed, exact_match)
  if (!is.character(pars)) {
    pars <- default_plot_pars(x)
    fixed <- FALSE
  }
  samples <- as.data.frame(x, pars = pars, add_chain = TRUE, fixed = fixed)
  pars <- names(samples)[!names(samples) %in% c("chain", "iter")] 
  if (!length(pars)) {
    stop2("No valid parameters selected.")
  }
  
  if (plot) {
    default_ask <- devAskNewPage()
    on.exit(devAskNewPage(default_ask))
    devAskNewPage(ask = FALSE)
  }
  n_plots <- ceiling(length(pars) / N)
  plots <- vector(mode = "list", length = n_plots)
  for (i in seq_len(n_plots)) {
    sub_pars <- pars[((i - 1) * N + 1):min(i * N, length(pars))]
    sub_samples <- samples[, c(sub_pars, "chain"), drop = FALSE]
    plots[[i]] <- bayesplot::mcmc_combo(
      sub_samples, combo = combo, gg_theme = theme, ...
    )
    if (plot) {
      plot(plots[[i]], newpage = newpage || i > 1)
      if (i == 1) {
        devAskNewPage(ask = ask)
      }
    }
  }
  invisible(plots) 
}

# list all parameter classes to be included in plots by default
default_plot_pars <- function(family) {
  c(fixef_pars(), "^sd_", "^cor_", "^sigma_", "^rescor_", 
    paste0("^", valid_dpars(family), "$"), "^delta$",
    "^theta", "^ar", "^ma", "^arr", "^sderr", "^lagsar", "^errorsar", 
    "^car", "^sdcar", "^sds_", "^sdgp_", "^lscale_")
}

#' MCMC Plots Implemented in \pkg{bayesplot} 
#' 
#' Convenient way to call MCMC plotting functions 
#' implemented in the \pkg{bayesplot} package.
#' 
#' @aliases stanplot stanplot.brmsfit 
#' 
#' @inheritParams plot.brmsfit
#' @param object An \R object typically of class \code{brmsfit}
#' @param pars Names of parameters to be plotted, 
#'   as given by a character vector or regular expressions. 
#'   By default, all parameters except for group-level and 
#'   smooth effects are plotted. May be ignored for some plots.
#' @param type The type of the plot. 
#'   Supported types are (as names) \code{hist}, \code{dens}, 
#'   \code{hist_by_chain}, \code{dens_overlay}, 
#'   \code{violin}, \code{intervals}, \code{areas}, \code{acf}, 
#'   \code{acf_bar},\code{trace}, \code{trace_highlight}, \code{scatter},
#'   \code{rhat}, \code{rhat_hist}, \code{neff}, \code{neff_hist}
#'   \code{nuts_acceptance}, \code{nuts_divergence},
#'   \code{nuts_stepsize}, \code{nuts_treedepth}, and \code{nuts_energy}. 
#'   For an overview on the various plot types see
#'   \code{\link[bayesplot:MCMC-overview]{MCMC-overview}}.
#' @param ... Additional arguments passed to the plotting functions.
#'   See \code{\link[bayesplot:MCMC-overview]{MCMC-overview}} for
#'   more details.
#' 
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object 
#'   that can be further customized using the \pkg{ggplot2} package.
#' 
#' @details 
#'   Also consider using the \pkg{shinystan} package available via 
#'   method \code{\link{launch_shinystan}} in \pkg{brms} for flexible 
#'   and interactive visual analysis. 
#' 
#' @examples
#' \dontrun{
#' model <- brm(count ~ zAge + zBase * Trt + (1|patient),
#'              data = epilepsy, family = "poisson")
#'              
#' # plot posterior intervals
#' mcmc_plot(model)
#' 
#' # only show population-level effects in the plots
#' mcmc_plot(model, pars = "^b_")
#' 
#' # show histograms of the posterior distributions
#' mcmc_plot(model, type = "hist")
#' 
#' # plot some diagnostics of the sampler
#' mcmc_plot(model, type = "neff")
#' mcmc_plot(model, type = "rhat")
#' 
#' # plot some diagnostics specific to the NUTS sampler
#' mcmc_plot(model, type = "nuts_acceptance")
#' mcmc_plot(model, type = "nuts_divergence")
#' }
#' 
#' @export
mcmc_plot.brmsfit <- function(object, pars = NA, type = "intervals", 
                              fixed = FALSE, exact_match = FALSE, ...) {
  contains_samples(object)
  object <- restructure(object)
  type <- as_one_character(type)
  fixed <- check_deprecated_fixed(fixed, exact_match)
  if (!is.character(pars)) {
    pars <- default_plot_pars(object)
    fixed <- FALSE
  }
  valid_types <- as.character(bayesplot::available_mcmc(""))
  valid_types <- sub("^mcmc_", "", valid_types)
  if (!type %in% valid_types) {
    stop2("Invalid plot type. Valid plot types are: \n",
          collapse_comma(valid_types))
  }
  mcmc_fun <- get(paste0("mcmc_", type), asNamespace("bayesplot"))
  mcmc_arg_names <- names(formals(mcmc_fun))
  mcmc_args <- list(...)
  if ("x" %in% mcmc_arg_names) {
    if (grepl("^nuts_", type)) {
      # x refers to a molten data.frame of NUTS parameters
      mcmc_args$x <- nuts_params(object)
    } else {
      # x refers to a data.frame of samples
      # TODO: replace once bayesplot supports posterior
      samples <- as.data.frame(
        object, pars = pars, add_chain = TRUE, fixed = fixed
      )
      if (!length(samples)) {
        stop2("No valid parameters selected.")
      }
      samples$iter <- NULL
      sel_pars <- names(samples)[!names(samples) %in% "chain"]
      if (type %in% c("scatter", "hex") && length(sel_pars) != 2L) {
        stop2("Exactly 2 parameters must be selected for this type.",
              "\nParameters selected: ", collapse_comma(sel_pars))
      }
      mcmc_args$x <- samples
    }
  }
  if ("lp" %in% mcmc_arg_names) {
    mcmc_args$lp <- log_posterior(object)
  }
  use_nuts <- isTRUE(object$algorithm == "sampling")
  if ("np" %in% mcmc_arg_names && use_nuts) {
    mcmc_args$np <- nuts_params(object)
  }
  interval_type <- type %in% c("intervals", "areas")
  if ("rhat" %in% mcmc_arg_names && !interval_type) {
    mcmc_args$rhat <- rhat(object)
  }
  if ("ratio" %in% mcmc_arg_names) {
    mcmc_args$ratio <- neff_ratio(object)
  }
  do_call(mcmc_fun, mcmc_args)
}

#' @rdname mcmc_plot.brmsfit
#' @export
mcmc_plot <- function(object, ...) {
  UseMethod("mcmc_plot")
}

# 'stanplot' has been deprecated in brms 2.10.6; remove in brms 3.0
#' @export
stanplot <- function(object, ...) {
  UseMethod("stanplot")
}

#' @export
stanplot.brmsfit <- function(object, ...) {
  warning2("Method 'stanplot' is deprecated. Please use 'mcmc_plot' instead.")
  mcmc_plot.brmsfit(object, ...)
}

#' Create a matrix of output plots from a \code{brmsfit} object
#'
#' A \code{\link[graphics:pairs]{pairs}} 
#' method that is customized for MCMC output.
#' 
#' @param x An object of class \code{brmsfit}
#' @inheritParams plot.brmsfit
#' @param ... Further arguments to be passed to 
#'   \code{\link[bayesplot:MCMC-scatterplots]{mcmc_pairs}}.
#'  
#' @details For a detailed description see  
#'   \code{\link[bayesplot:MCMC-scatterplots]{mcmc_pairs}}.
#'  
#' @examples 
#' \dontrun{
#' fit <- brm(count ~ zAge + zBase * Trt 
#'            + (1|patient) + (1|visit), 
#'            data = epilepsy, family = "poisson")  
#' pairs(fit, pars = parnames(fit)[1:3], fixed = TRUE)
#' pairs(fit, pars = "^sd_")
#' }
#'
#' @export
pairs.brmsfit <- function(x, pars = NA, fixed = FALSE, exact_match = FALSE, ...) {
  fixed <- check_deprecated_fixed(fixed, exact_match)
  if (!is.character(pars)) {
    pars <- default_plot_pars(x)
    fixed <- FALSE
  }
  samples <- as.data.frame(x, pars = pars, add_chain = TRUE, fixed = fixed)
  samples$iter <- NULL
  bayesplot::mcmc_pairs(samples, ...)
}

#' Default \pkg{bayesplot} Theme for \pkg{ggplot2} Graphics
#' 
#' This theme is imported from the \pkg{bayesplot} package.
#' See \code{\link[bayesplot:theme_default]{theme_default}}
#' for a complete documentation.
#' 
#' @name theme_default
#' 
#' @param base_size base font size
#' @param base_family base font family
#' 
#' @return A \code{theme} object used in \pkg{ggplot2} graphics.
#' 
#' @importFrom bayesplot theme_default
#' @export theme_default
NULL

