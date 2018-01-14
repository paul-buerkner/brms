#' @rdname marginal_effects
#' @method plot brmsMarginalEffects
#' @export 
plot.brmsMarginalEffects <- function(
  x, ncol = NULL, points = FALSE, rug = FALSE, mean = TRUE, 
  jitter_width = 0, stype = c("contour", "raster"),
  line_args = list(), cat_args = list(), errorbar_args = list(), 
  surface_args = list(), spaghetti_args = list(), point_args = list(), 
  rug_args = list(), theme = NULL, ask = TRUE, plot = TRUE, ...
) {
  dots <- list(...)
  plot <- use_alias(plot, dots$do_plot)
  stype <- match.arg(stype)
  smooths_only <- isTRUE(attr(x, "smooths_only"))
  if (points && smooths_only) {
    stop2("Argument 'points' is invalid for objects ", 
          "returned by 'marginal_smooths'.")
  }
  if (!is_equal(jitter_width, 0)) {
    warning2("'jitter_width' is deprecated. Please use ",
             "'point_args = list(width = <width>)' instead.")
  }
  if (!is.null(theme)) {
    if (!is.theme(theme)) {
      stop2("Argument 'theme' should be a 'theme' object.")
    }
    pb_colour <- theme$plot.background$colour
  } else {
    pb_colour <- theme_get()$plot.background$colour
  }
  is_theme_black <- isTRUE(pb_colour == "black")
  if (plot) {
    default_ask <- devAskNewPage()
    on.exit(devAskNewPage(default_ask))
    devAskNewPage(ask = FALSE)
  }
  dont_replace <- c("mapping", "data", "inherit.aes") 
  plots <- named_list(names(x))
  for (i in seq_along(x)) {
    response <- attr(x[[i]], "response")
    effects <- attr(x[[i]], "effects")
    ncond <- length(unique(x[[i]]$cond__))
    surface <- isTRUE(attr(x[[i]], "surface"))
    # for backwards compatibility with brms < 1.4.0
    surface <- surface || isTRUE(attr(x[[i]], "contour"))
    ordinal <- isTRUE(attr(x[[i]], "ordinal"))
    if (surface || ordinal) {
      # surface plots for two dimensional interactions or ordinal plots
      plots[[i]] <- ggplot(x[[i]]) + aes_string(effects[1], effects[2])
      if (ordinal) {
        width <- ifelse(is_like_factor(x[[i]][[effects[1]]]), 0.9, 1)
        .surface_args <- nlist(
          mapping = aes_(fill = ~ estimate__), 
          height = 0.9, width = width
        )
        replace_args(.surface_args, dont_replace) <- surface_args
        plots[[i]] <- plots[[i]] + 
          do.call(geom_tile, .surface_args) + 
          scale_fill_gradientn(colors = viridis6(), name = "Probability") + 
          ylab(response)
      } else if (stype == "contour") {
        .surface_args <- nlist(
          mapping = aes_(z = ~ estimate__, colour = ~ ..level..),
          bins = 30, size = 1.3
        )
        replace_args(.surface_args, dont_replace) <- surface_args
        plots[[i]] <- plots[[i]] + 
          do.call(geom_contour, .surface_args) +
          scale_color_gradientn(colors = viridis6(), name = response)
      } else if (stype == "raster") {
        .surface_args <- nlist(mapping = aes_(fill = ~ estimate__))
        replace_args(.surface_args, dont_replace) <- surface_args
        plots[[i]] <- plots[[i]] + 
          do.call(geom_raster, .surface_args) + 
          scale_fill_gradientn(colors = viridis6(), name = response)
      }
    } else {
      # plot effects of single predictors or two-way interactions
      gvar <- if (length(effects) == 2L) effects[2]
      spaghetti <- attr(x[[i]], "spaghetti")
      plots[[i]] <- ggplot(x[[i]]) + 
        aes_string(x = effects[1], y = "estimate__", colour = gvar) +
        ylab(response)
      if (is.null(spaghetti)) {
        plots[[i]] <- plots[[i]] +
          aes_string(ymin = "lower__", ymax = "upper__", fill = gvar)  
      }
      # extract suggested colors for later use
      colors <- ggplot_build(plots[[i]])
      colors <- unique(colors$data[[1]][["colour"]])
      if (points) {
        # add points first so that they appear behind the predictions
        .point_args <- list(
          mapping = aes_string(x = effects[1], y = "resp__"),
          data = attr(x[[i]], "points"), inherit.aes = FALSE,
          size = 2 / ncond^0.25, height = 0, width = jitter_width
        )
        is_factor_gvar <- is_like_factor(attr(x[[i]], "points")[, gvar])
        if (is_factor_gvar) {
          expr_gvar <- parse(text = gvar)[[1]]
          .point_args$mapping$colour <- expr_gvar
          .point_args$mapping$fill <- expr_gvar
        } else if (is_theme_black) {
          .point_args$colour <- "white"
        }
        replace_args(.point_args, dont_replace) <- point_args
        plots[[i]] <- plots[[i]] + 
          do.call(geom_jitter, .point_args)
      }
      if (!is.null(spaghetti)) {
        # add a regression line for each sample separately
        .spaghetti_args <- list(
          aes_string(group = "sample__", colour = gvar),
          data = spaghetti, stat = "identity", size = 0.5
        )
        if (length(effects) == 1L) {
          .spaghetti_args$colour <- alpha("blue", 0.1)
        } else {
          # workaround to get transparent lines
          plots[[i]] <- plots[[i]] +
            scale_color_manual(values = alpha(colors, 0.1))
        }
        replace_args(.spaghetti_args, dont_replace) <- spaghetti_args
        plots[[i]] <- plots[[i]] +
          do.call(geom_smooth, .spaghetti_args)
      }
      if (is.numeric(x[[i]][, effects[1]])) {
        # line plots for numeric predictors
        .line_args <- list(stat = "identity")
        if (!is.null(spaghetti)) {
          # display a white mean regression line
          .line_args$mapping <- aes_string(group = gvar)
          .line_args$colour <- alpha("white", 0.8)
        }
        replace_args(.line_args, dont_replace) <- line_args
        if (mean || is.null(spaghetti)) {
          plots[[i]] <- plots[[i]] + 
            do.call(geom_smooth, .line_args)
        }
        if (rug) {
          .rug_args <- list(
            aes_string(x = effects[1]), sides = "b", 
            data = attr(x[[i]], "points"), inherit.aes = FALSE
          )
          if (is.null(gvar) && is_theme_black) {
            .rug_args$colour <- "white"
          }
          replace_args(.rug_args, dont_replace) <- rug_args
          plots[[i]] <- plots[[i]] + 
            do.call(geom_rug, .rug_args)
        }
      } else {
        # points and errorbars for factors
        .cat_args <- list(
          position = position_dodge(width = 0.4),
          size = 4 / ncond^0.25
        )
        .errorbar_args <- list(
          position = position_dodge(width = 0.4), 
          width = 0.3
        )
        if (is.null(gvar) && is_theme_black) {
          .cat_args$colour <- .errorbar_args$colour <- "white"
        }
        replace_args(.cat_args, dont_replace) <- cat_args
        replace_args(.errorbar_args, dont_replace) <- errorbar_args
        plots[[i]] <- plots[[i]] + 
          do.call(geom_point, .cat_args) +
          do.call(geom_errorbar, .errorbar_args)
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
      if (i == 1) {
        devAskNewPage(ask = ask)
      }
    }
  }
  invisible(plots)
}

#' @rdname brmshypothesis
#' @method plot brmshypothesis
#' @export
plot.brmshypothesis <- function(x, N = 5, ignore_prior = FALSE,
                                chars = 40, colors = NULL,
                                theme = NULL, ask = TRUE, 
                                plot = TRUE,  ...) {
  dots <- list(...)
  if (!is.data.frame(x$samples)) {
    stop2("No posterior samples found")
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
  
  samples <- cbind(x$samples, Type = "Posterior")
  if (!ignore_prior) {
    samples <- rbind(samples, cbind(x$prior_samples, Type = "Prior"))
  }
  if (plot) {
    default_ask <- devAskNewPage()
    on.exit(devAskNewPage(default_ask))
    devAskNewPage(ask = FALSE)
  }
  hyps <- limit_chars(rownames(x$hypothesis), chars = chars)
  names(samples)[seq_along(hyps)] <- hyps
  n_plots <- ceiling(length(hyps) / N)
  plots <- vector(mode = "list", length = n_plots)
  for (i in seq_len(n_plots)) {
    rel_hyps <- hyps[((i - 1) * N + 1):min(i * N, length(hyps))]
    sub_samples <- cbind(
      utils::stack(samples[, rel_hyps, drop = FALSE]),
      samples[, "Type", drop = FALSE]
    )
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
