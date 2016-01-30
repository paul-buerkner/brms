trace_density_plot <- function(par, x, theme = "classic") {
  # trace and density plots for one parameter
  #
  # Args:
  #   par: a single character string 
  #   x: a data.frame containing the samples
  #
  # Returns:
  #   a list containing the trace and density plot of one parameter
  if (!is.character(par) || length(par) != 1)
    stop("par must be a character string")
  if (!is.data.frame(x))
    stop("x must be a data.frame")
  names(x)[match(par, names(x))] <- "value" 
  # trace plot
  trace <- ggplot(x, aes_string(x = "iter", y = "value", group = "chain", 
                                colour = "chain")) +
    geom_line(alpha = 0.7) + 
    xlab("") + ylab("") + ggtitle(paste("Trace of", par)) + 
    do.call(paste0("theme_", theme), args = list()) + 
    theme(legend.position = "none",
          plot.title = element_text(size = 15, vjust = 1),
          plot.margin = grid::unit(c(0.2, 0, -0.5, -0.5), "lines"))
  # density plot
  density <- ggplot(x, aes_string(x = "value")) + 
    geom_density(aes_string(fill = "chain"), alpha = 0.5) + 
    xlab("") + ylab("") + ggtitle(paste("Density of", par)) + 
    do.call(paste0("theme_", theme), args = list()) +
    theme(plot.title = element_text(size = 15, vjust = 1),
          plot.margin = grid::unit(c(0.2, 0, -0.5, -0.5), "lines"))
  list(trace, density)
}

#' @rdname marginal_effects
#' @method plot brmsMarginalEffects
#' @export 
plot.brmsMarginalEffects <- function(x, ncol = NULL, rug = FALSE,
                                     theme = "gray", ask = TRUE,
                                     do_plot = TRUE, ...) {
  # Compute marginal effects plots using ggplot2
  # Returns:
  #   A list of ggplot objects
  if (do_plot) {
    default_ask <- devAskNewPage()
    on.exit(devAskNewPage(default_ask))
    devAskNewPage(ask = FALSE)
  }
  plots <- list()
  for (i in seq_along(x)) {
    response <- attributes(x[[i]])$response
    effects <- attributes(x[[i]])$effects
    plots[[i]] <- ggplot(data = x[[i]]) + 
      aes_string(x = effects, y = "Estimate", ymin = "`2.5%ile`",
                 ymax = "`97.5%ile`") + ylab(response) +
      do.call(paste0("theme_", theme), args = list())
    nMargins <- length(unique(x[[i]]$MargRow))
    if (nMargins > 1) {
      # one plot per row of marginal_data
      if (is.null(ncol)) ncol <- max(floor(sqrt(nMargins)), 3) 
      plots[[i]] <- plots[[i]] + 
        facet_wrap("MargRow", ncol = ncol)
    }
    if (length(effects) == 2) {
      # differentiate by colour in case of interaction effects
      plots[[i]] <- plots[[i]] + 
        aes_string(colour = effects[2], fill = effects[2])
    }
    if (is.numeric(x[[i]][, effects[1]])) {
      # smooth plots for numeric predictors
      plots[[i]] <- plots[[i]] + geom_smooth(stat = "identity")
      if (rug) {
        plots[[i]] <- plots[[i]] + 
          geom_rug(aes_string(x = effects[1]), sides = "b", 
                   data = attributes(x[[i]])$rug, inherit.aes = FALSE)
      }
    } else {
      # barplots for factors
      plots[[i]] <- plots[[i]] + 
        geom_bar(stat = "identity", position = "dodge") +
        geom_errorbar(position = position_dodge(0.9), colour = "black")
    }
    if (do_plot) {
      plot(plots[[i]])
      if (i == 1) devAskNewPage(ask = ask)
    }
  }
  invisible(plots)
}

#' @rdname hypothesis
#' @method plot brmshypothesis
#' @export
plot.brmshypothesis <- function(x, N = 5, ignore_prior = FALSE, 
                                theme = "classic", ask = TRUE, 
                                do_plot = TRUE, newpage = TRUE, ...) {
  if (!is.data.frame(x$samples)) {
    stop("No posterior samples found")
  }
  .plot_fun <- function(i) {
    # create the ggplot object for each hypothesis
    # Args: i: index variable
    ggplot(x$samples, aes_string(x = hypnames[i])) + 
      geom_density(aes_string(fill = "Type"), alpha = 0.5, na.rm = TRUE) + 
      scale_fill_manual(values = c("red", "blue")) + 
      xlab("") + ylab("") + ggtitle(hyps[i]) + 
      do.call(paste0("theme_", theme), args = list())
  }
  if (ignore_prior) {
    x$samples <- subset(x$samples, x$samples$Type == "posterior")
  }
  if (do_plot) {
    default_ask <- devAskNewPage()
    on.exit(devAskNewPage(default_ask))
    devAskNewPage(ask = FALSE)
  }
  hyps <- rownames(x$hypothesis)
  hypnames <- names(x$samples)[seq_along(hyps)]
  n_plots <- ceiling(length(hyps) / N)
  plots <- vector(mode = "list", length = n_plots)
  for (i in 1:n_plots) {
    I <- ((i - 1) * N + 1):min(i * N, length(hyps))
    temp_plot <- lapply(I, .plot_fun)
    plots[[i]] <- arrangeGrob(grobs = temp_plot, ncol = 1, 
                              nrow = length(temp_plot), ...)
    if (do_plot) {
      if (newpage || i > 1) grid.newpage()
      grid.draw(plots[[i]])
      if (i == 1) devAskNewPage(ask = ask)
    }
  }
  invisible(plots) 
}