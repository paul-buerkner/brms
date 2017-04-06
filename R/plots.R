#' @rdname marginal_effects
#' @method plot brmsMarginalEffects
#' @export 
plot.brmsMarginalEffects <- function(x, ncol = NULL, points = FALSE, 
                                     rug = FALSE, jitter_width = 0,
                                     stype = c("contour", "raster"),
                                     theme = bayesplot::theme_default(), 
                                     ask = TRUE, plot = TRUE, ...) {
  # Compute marginal effects plots using ggplot2
  # Returns:
  #   A list of ggplot objects
  dots <- list(...)
  plot <- use_alias(plot, dots$do_plot)
  stype <- match.arg(stype)
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
    response <- attr(x[[i]], "response")
    effects <- attr(x[[i]], "effects")
    ncond <- length(unique(x[[i]]$cond__))
    surface <- isTRUE(attr(x[[i]], "surface"))
    # for backwards compatibility with brms < 1.4.0
    surface <- surface || isTRUE(attr(x[[i]], "contour"))
    if (surface) {
      # surface plots for two dimensional smooths
      plots[[i]] <- ggplot(x[[i]]) +
        aes_string(effects[1], effects[2])
      if (stype == "contour") {
        plots[[i]] <- plots[[i]] + 
          geom_contour(aes_(z = ~ estimate__, colour = ~ ..level..),
                       bins = 30, size = 1.3) +
          scale_color_gradientn(colors = viridis6(), name = response)
      } else if (stype == "raster") {
        plots[[i]] <- plots[[i]] + 
          geom_raster(aes_(fill = ~ estimate__)) + 
          scale_fill_gradientn(colors = viridis6(), name = response)
      }
    } else {
      # plot effects of single predictors / smooths
      # as well as two-way interactions
      gvar <- if (length(effects) == 2L) effects[2]
      plots[[i]] <- ggplot(x[[i]]) + 
        aes_string(x = effects[1], y = "estimate__",
                   ymin = "lower__", ymax = "upper__", 
                   colour = gvar, fill = gvar) + 
        ylab(response)
      if (points) {
        # show the data as points in the plot
        # add points first so that they appear behind the regression lines
        aes_points <- aes_string(x = effects[1], y = "resp__")
        if (is.factor(attr(x[[i]], "points")[, gvar])) {
          aes_points$colour <- parse(text = gvar)[[1]]
        }
        plots[[i]] <- plots[[i]] + 
          geom_jitter(aes_points, shape = 1, size = 4 / ncond^0.25,
                      data = attr(x[[i]], "points"), inherit.aes = FALSE,
                      height = 0, width = jitter_width)
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
                     size = 4 / ncond^0.25) + 
          geom_errorbar(position = position_dodge(width = 0.4),
                        width = 0.3)
      }
    }
    if (ncond > 1L) {
      # one plot per row of conditions
      if (is.null(ncol)) {
        ncol <- max(floor(sqrt(ncond)), 3)
      }
      plots[[i]] <- plots[[i]] + 
        facet_wrap("cond__", ncol = ncol)
    }
    plots[[i]] <- plots[[i]] + theme
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
                                chars = 40, colors = NULL,
                                theme = bayesplot::theme_default(),
                                ask = TRUE, plot = TRUE,  ...) {
  dots <- list(...)
  if (!is.data.frame(x$samples)) {
    stop("No posterior samples found", call. = FALSE)
  }
  plot <- use_alias(plot, dots$do_plot)
  if (is.null(colors)) {
    colors <- bayesplot::color_scheme_get()[c(6, 2)]
    colors <- unname(unlist(colors))
  }
  if (length(colors) != 2L) {
    stop2("Argument 'colors' must be of length 2.")
  }
  
  .plot_fun <- function(samples) {
    ggplot(samples, aes_string(x = "values")) + 
      facet_wrap("ind", ncol = 1, scales = "free") +
      geom_density(aes_string(fill = "Type"), 
                   alpha = 0.7, na.rm = TRUE) + 
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

viridis6 <- function() {
  # colours taken from the viridis package
  c("#440154", "#414487", "#2A788E", "#22A884", "#7AD151", "#FDE725")
}
