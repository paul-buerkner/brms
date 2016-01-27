marginal_plot_internal <- function(x, ncol = NULL, theme = "gray",
                                   do_plot = TRUE) {
  # Compute marginal plots using ggplot2
  # Args:
  #   x: A list of marginal predictions computed by marginal_plot
  #   do_plot: logical; Should plots be printed on the active device?
  plot_list = list()
  for (i in seq_along(x)) {
    response = attributes(x[[i]])$response
    preds = attributes(x[[i]])$predictors
    plot_list[[i]] <- ggplot(data = x[[i]]) + 
      aes_string(x = preds, y = "Estimate", ymin = "`2.5%ile`",
                 ymax = "`97.5%ile`") + ylab(response) +
      do.call(paste0("theme_", theme), args = list())
    nMargins <- length(unique(x[[i]]$MarginalRow))
    if (nMargins > 1) {
      # one plot per row of marginal_data
      if (is.null(ncol)) ncol <- max(floor(sqrt(nMargins)), 3) 
      plot_list[[i]] <- plot_list[[i]] + 
        facet_wrap("MarginalRow", ncol = ncol)
    }
    if (length(preds) == 2) {
      # differentiate by colour in case of interaction effects
      plot_list[[i]] <- plot_list[[i]] + 
        aes_string(colour = preds[2], fill = preds[2])
    }
    if (is.numeric(x[[i]][, preds[1]])) {
      # smooth plots for numeric predictors
      plot_list[[i]] <- plot_list[[i]] + geom_smooth(stat = "identity")
      if (!is.null(attributes(x[[i]])$rug)) {
        plot_list[[i]] <- plot_list[[i]] + 
          geom_rug(aes_string(x = preds[1]), sides = "b", 
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