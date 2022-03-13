#' (Deprecated) Black Theme for \pkg{ggplot2} Graphics
#'
#' A black theme for ggplot graphics inspired by a blog post of Jon Lefcheck
#' (\url{https://jonlefcheck.net/2013/03/11/black-theme-for-ggplot2-2/}).
#'
#' @param base_size base font size
#' @param base_family base font family
#'
#' @return A \code{theme} object used in \pkg{ggplot2} graphics.
#'
#' @details When using \code{theme_black} in plots powered by the
#' \pkg{bayesplot} package such as \code{pp_check} or \code{stanplot},
#' I recommend using the \code{"viridisC"} color scheme (see examples).
#'
#' @examples
#' \dontrun{
#' # change default ggplot theme
#' ggplot2::theme_set(theme_black())
#'
#' # change default bayesplot color scheme
#' bayesplot::color_scheme_set("viridisC")
#'
#' # fit a simple model
#' fit <- brm(count ~ zAge + zBase * Trt + (1|patient),
#'            data = epilepsy, family = poisson(), chains = 2)
#' summary(fit)
#'
#' # create various plots
#' plot(marginal_effects(fit), ask = FALSE)
#' pp_check(fit)
#' mcmc_plot(fit, type = "hex", variable = c("b_Intercept", "b_Trt1"))
#' }
#'
#' @export
theme_black = function(base_size = 12, base_family = "") {
  warning2("'theme_black' is deprecated. Please use the 'ggdark' package ",
           "for dark ggplot themes.")
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # axis options
      axis.line = element_blank(),
      axis.text.x = element_text(
        size = base_size * 0.8, color = "white", lineheight = 0.9
      ),
      axis.text.y = element_text(
        size = base_size * 0.8, color = "white", lineheight = 0.9
      ),
      axis.ticks = element_line(color = "white", size  =  0.2),
      axis.title.x = element_text(
        size = base_size, color = "white", margin = margin(10, 0, 0, 0)
      ),
      axis.title.y = element_text(
        size = base_size, color = "white", angle = 90,
        margin = margin(0, 10, 0, 0)
      ),
      axis.ticks.length = unit(0.3, "lines"),
      # legend options
      legend.background = element_rect(color = NA, fill = "black"),
      legend.key = element_rect(color = "white", fill = "black"),
      legend.key.size = unit(1.2, "lines"),
      legend.key.height = NULL,
      legend.key.width = NULL,
      legend.text = element_text(size = base_size * 0.8, color = "white"),
      legend.title = element_text(
        size = base_size * 0.8, face = "bold", hjust = 0, color = "white"
      ),
      legend.position = "right",
      legend.text.align = NULL,
      legend.title.align = NULL,
      legend.direction = "vertical",
      legend.box = NULL,
      # panel options
      panel.background = element_rect(fill = "black", color  =  NA),
      panel.border = element_rect(fill = NA, color = "white"),
      panel.grid.major = element_line(color = "grey35"),
      panel.grid.minor = element_line(color = "grey20"),
      panel.spacing = unit(0.5, "lines"),
      # facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),
      strip.text.x = element_text(
        size = base_size * 0.8, color = "white", margin = margin(3, 0, 4, 0)
      ),
      strip.text.y = element_text(
        size = base_size * 0.8, color = "white", angle = -90
      ),
      # plot options
      plot.background = element_rect(color = "black", fill = "black"),
      plot.title = element_text(size = base_size * 1.2, color = "white"),
      plot.margin = unit(rep(1, 4), "lines")
    )
}
