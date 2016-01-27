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

margins_plot_internal <- function(x, ncol = NULL, theme = "gray",
                                  do_plot = TRUE) {
  # Compute marginal plots using ggplot2
  # Args:
  #   x: A list of marginal predictions computed by margins_plot
  #   ncol: Number of columns in each plot
  #   theme: The ggplot theme
  #   do_plot: logical; Should plots be printed on the active device?
  # Returns:
  #   A list of ggplot objects
  plot_list = list()
  for (i in seq_along(x)) {
    response = attributes(x[[i]])$response
    effects = attributes(x[[i]])$effects
    plot_list[[i]] <- ggplot(data = x[[i]]) + 
      aes_string(x = effects, y = "Estimate", ymin = "`2.5%ile`",
                 ymax = "`97.5%ile`") + ylab(response) +
      do.call(paste0("theme_", theme), args = list())
    nMargins <- length(unique(x[[i]]$MargRow))
    if (nMargins > 1) {
      # one plot per row of marginal_data
      if (is.null(ncol)) ncol <- max(floor(sqrt(nMargins)), 3) 
      plot_list[[i]] <- plot_list[[i]] + 
        facet_wrap("MargRow", ncol = ncol)
    }
    if (length(effects) == 2) {
      # differentiate by colour in case of interaction effects
      plot_list[[i]] <- plot_list[[i]] + 
        aes_string(colour = effects[2], fill = effects[2])
    }
    if (is.numeric(x[[i]][, effects[1]])) {
      # smooth plots for numeric predictors
      plot_list[[i]] <- plot_list[[i]] + geom_smooth(stat = "identity")
      if (!is.null(attributes(x[[i]])$rug)) {
        plot_list[[i]] <- plot_list[[i]] + 
          geom_rug(aes_string(x = effects[1]), sides = "b", 
                   data = attributes(x[[i]])$rug, inherit.aes = FALSE)
      }
    } else {
      # barplots for factors
      plot_list[[i]] = plot_list[[i]] + 
        geom_bar(stat = "identity", position = "dodge") +
        geom_errorbar(position = position_dodge(0.9), colour = "black")
    }
    if (do_plot) {
      plot(plot_list[[i]])
    }
  }
  invisible(plot_list)
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
  if (do_plot) {
    invisible(plots) 
  } else {
    plots
  }
}