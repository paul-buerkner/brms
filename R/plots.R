trace_density_plot <- function(x, theme = ggplot2::theme()) {
  # trace and density plots for one parameter
  #
  # Args:
  #   x: a data.frame containing the samples
  #   theme: a ggplot2 theme object
  #
  # Returns:
  #   a list containing trace and density plots 
  if (!is.data.frame(x))
    stop("x must be a data.frame")
  out <- list()
  out$trace <- ggplot(x, aes_string(x = "iter", y = "values", 
                                    colour = "chain")) +
    geom_line(alpha = 0.7) +
    facet_wrap("ind", ncol = 1, scales = "free") +
    xlab("") + ylab("") + ggtitle(paste("Trace")) + theme + 
    theme(legend.position = "none",
          plot.title = element_text(size = 16),
          strip.text.x = element_text(size = 15),
          plot.margin = grid::unit(c(0.2, 0.2, -0.5, -0.5), "lines"))
  
  out$density <- ggplot(x, aes_string(x = "values")) + 
    geom_density(aes_string(fill = "chain"), alpha = 0.5) +
    facet_wrap("ind", ncol = 1, scales = "free") + 
    xlab("") + ylab("") + ggtitle(paste("Density")) + theme +
    theme(plot.title = element_text(size = 16),
          strip.text.x = element_text(size = 15),
          plot.margin = grid::unit(c(0.2, 0, -0.5, -0.5), "lines"))
  out
}

#' @rdname marginal_effects
#' @method plot brmsMarginalEffects
#' @export 
plot.brmsMarginalEffects <- function(x, ncol = NULL, rug = FALSE,
                                     theme = ggplot2::theme(), 
                                     ask = TRUE, do_plot = TRUE, ...) {
  # Compute marginal effects plots using ggplot2
  # Returns:
  #   A list of ggplot objects
  if (do_plot) {
    default_ask <- devAskNewPage()
    on.exit(devAskNewPage(default_ask))
    devAskNewPage(ask = FALSE)
  }
  plots <- setNames(vector("list", length(x)), names(x))
  for (i in seq_along(x)) {
    response <- attributes(x[[i]])$response
    effects <- attributes(x[[i]])$effects
    plots[[i]] <- ggplot(data = x[[i]]) + 
      aes_string(x = effects, y = "Estimate", ymin = "lowerCI",
                 ymax = "upperCI") + ylab(response) + theme
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
      # pointrange for factors
      plots[[i]] <- plots[[i]] + 
        geom_pointrange(position = position_dodge(width = 0.4), fatten = 7)
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
                                theme = ggplot2::theme(), ask = TRUE, 
                                do_plot = TRUE, ...) {
  if (!is.data.frame(x$samples)) {
    stop("No posterior samples found")
  }
  .plot_fun <- function(samples) {
    ggplot(samples, aes_string(x = "values")) + 
      facet_wrap("ind", ncol = 1, scales = "free") +
      geom_density(aes_string(fill = "Type"), alpha = 0.5, na.rm = TRUE) + 
      scale_fill_manual(values = c("red", "blue")) + 
      ggtitle(paste("Hypothesis for class", x$class)) + 
      xlab("") + ylab("") + theme
  }
  if (ignore_prior) {
    x$samples[x$samples$Type == "prior", ] <- NA
  }
  if (do_plot) {
    default_ask <- devAskNewPage()
    on.exit(devAskNewPage(default_ask))
    devAskNewPage(ask = FALSE)
  }
  hyps <- rownames(x$hypothesis)
  names(x$samples)[seq_along(hyps)] <- hyps
  n_plots <- ceiling(length(hyps) / N)
  plots <- vector(mode = "list", length = n_plots)
  for (i in 1:n_plots) {
    rel_hyps <- hyps[((i - 1) * N + 1):min(i * N, length(hyps))]
    sub_samples <- cbind(utils::stack(x$samples[, rel_hyps, drop = FALSE]),
                         x$samples[, "Type", drop = FALSE])
    # make sure that parameters appear in the original order
    sub_samples$ind <- with(sub_samples, factor(ind, levels = unique(ind)))
    plots[[i]] <- .plot_fun(sub_samples)
    if (do_plot) {
      plot(plots[[i]])
      if (i == 1) devAskNewPage(ask = ask)
    }
  }
  invisible(plots) 
}