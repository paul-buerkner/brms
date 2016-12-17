#' @rdname marginal_effects
#' @method plot brmsMarginalEffects
#' @export 
plot.brmsMarginalEffects <- function(x, ncol = NULL, 
                                     points = FALSE, rug = FALSE,
                                     theme = bayesplot::theme_default(), 
                                     ask = TRUE, plot = TRUE, ...) {
  # Compute marginal effects plots using ggplot2
  # Returns:
  #   A list of ggplot objects
  dots <- list(...)
  plot <- use_alias(plot, dots$do_plot)
  smooths_only <- isTRUE(attr(x, "smooths_only"))
  if (points && smooths_only) {
    stop2("Argument 'points' is invalid for objects ", 
          "returned by 'marginal_smooths'.")
  }
  if (plot) {
    default_ask <- devAskNewPage()
    on.exit(devAskNewPage(default_ask))
    devAskNewPage(ask = FALSE)
  }
  plots <- named_list(names(x))
  for (i in seq_along(x)) {
    response <- attributes(x[[i]])$response
    effects <- attributes(x[[i]])$effects
    if (smooths_only && length(effects) == 2L) {
      # contour plot for two dimensional smooths
      # colours taken from the viridis package
      magma5 <- c("#000004FF", "#51127CFF", "#B63679FF", 
                  "#FB8861FF", "#FCFDBFFF")
      plots[[i]] <- ggplot(x[[i]], aes_string(effects[1], effects[2])) + 
        geom_contour(aes_string(z = "Estimate", colour = "..level.."), 
                     bins = 30, size = 1.3) +
        scale_color_gradientn(colors = magma5, name = response) +
        theme
    } else {
      # plot effects of single predictors / smooths
      # as well as two-way interactions
      gvar <- if (length(effects) == 2L) effects[2]
      plots[[i]] <- ggplot(x[[i]]) + 
        aes_string(x = effects, y = "Estimate", ymin = "lowerCI",
                   ymax = "upperCI", colour = gvar, fill = gvar) + 
        ylab(response) + theme
      nCond <- length(unique(x[[i]]$MargCond))
      if (points) {
        # show the data as points in the plot
        # add points first so that they appear behind the regression lines
        aes_points <- aes_string(x = effects[1], y = ".RESP")
        if (is.factor(attr(x[[i]], "points")[, gvar])) {
          aes_points$colour <- parse(text = gvar)[[1]]
        }
        plots[[i]] <- plots[[i]] + 
          geom_point(aes_points, shape = 1, size = 4 / nCond^0.25,
                     data = attr(x[[i]], "points"), inherit.aes = FALSE)
      }
      if (nCond > 1L) {
        # one plot per row of marginal_data
        if (is.null(ncol)) ncol <- max(floor(sqrt(nCond)), 3) 
        plots[[i]] <- plots[[i]] + 
          facet_wrap("MargCond", ncol = ncol)
      }
      if (is.numeric(x[[i]][, effects[1]])) {
        # smooth plots for numeric predictors
        plots[[i]] <- plots[[i]] + geom_smooth(stat = "identity")
        if (rug) {
          plots[[i]] <- plots[[i]] +
            geom_rug(aes_string(x = effects[1]),
                     sides = "b", data = attr(x[[i]], "points"), 
                     inherit.aes = FALSE)
        }
      } else {
        # points and errorbars for factors
        plots[[i]] <- plots[[i]] + 
          geom_point(position = position_dodge(width = 0.4),
                     size = 4 / nCond^0.25) + 
          geom_errorbar(position = position_dodge(width = 0.4),
                        width = 0.3)
      }
    }
    if (plot) {
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
                                chars = 40, colors = c("red", "blue"),
                                theme = bayesplot::theme_default(),
                                ask = TRUE, plot = TRUE,  ...) {
  dots <- list(...)
  if (!is.data.frame(x$samples)) {
    stop("No posterior samples found", call. = FALSE)
  }
  plot <- use_alias(plot, dots$do_plot)
  
  .plot_fun <- function(samples) {
    ggplot(samples, aes_string(x = "values")) + 
      facet_wrap("ind", ncol = 1, scales = "free") +
      geom_density(aes_string(fill = "Type"), 
                   alpha = 0.5, na.rm = TRUE) + 
      scale_fill_manual(values = colors) + 
      xlab("") + ylab("") + theme
  }
  
  if (ignore_prior) {
    x$samples[x$samples$Type == "prior", ] <- NA
  }
  if (plot) {
    default_ask <- devAskNewPage()
    on.exit(devAskNewPage(default_ask))
    devAskNewPage(ask = FALSE)
  }
  hyps <- limit_chars(rownames(x$hypothesis), chars = chars)
  names(x$samples)[seq_along(hyps)] <- hyps
  n_plots <- ceiling(length(hyps) / N)
  plots <- vector(mode = "list", length = n_plots)
  for (i in seq_len(n_plots)) {
    rel_hyps <- hyps[((i - 1) * N + 1):min(i * N, length(hyps))]
    sub_samples <- cbind(utils::stack(x$samples[, rel_hyps, drop = FALSE]),
                         x$samples[, "Type", drop = FALSE])
    # make sure that parameters appear in the original order
    sub_samples$ind <- with(sub_samples, factor(ind, levels = unique(ind)))
    plots[[i]] <- .plot_fun(sub_samples)
    if (plot) {
      plot(plots[[i]])
      if (i == 1) devAskNewPage(ask = ask)
    }
  }
  invisible(plots) 
}
