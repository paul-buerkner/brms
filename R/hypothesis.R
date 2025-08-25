#' Non-Linear Hypothesis Testing
#'
#' Perform non-linear hypothesis testing for all model parameters.
#'
#' @param x An \code{R} object. If it is no \code{brmsfit} object,
#'  it must be coercible to a \code{data.frame}.
#'  In the latter case, the variables used in the \code{hypothesis} argument
#'  need to correspond to column names of \code{x}, while the rows
#'  are treated as representing posterior draws of the variables.
#' @param hypothesis A character vector specifying one or more
#'  non-linear hypothesis concerning parameters of the model.
#' @param class A string specifying the class of parameters being tested.
#'  Default is "b" for population-level effects.
#'  Other typical options are "sd" or "cor".
#'  If \code{class = NULL}, all parameters can be tested
#'  against each other, but have to be specified with their full name
#'  (see also \code{\link[brms:draws-index-brms]{variables}})
#' @param group Name of a grouping factor to evaluate only
#'  group-level effects parameters related to this grouping factor.
#' @param alpha The alpha-level of the tests (default is 0.05;
#'  see 'Details' for more information).
#' @param robust If \code{FALSE} (the default) the mean is used as
#'  the measure of central tendency and the standard deviation as
#'  the measure of variability. If \code{TRUE}, the median and the
#'  median absolute deviation (MAD) are applied instead.
#' @param scope Indicates where to look for the variables specified in
#'  \code{hypothesis}. If \code{"standard"}, use the full parameter names
#'  (subject to the restriction given by \code{class} and \code{group}).
#'  If \code{"coef"} or \code{"ranef"}, compute the hypothesis for all levels
#'  of the grouping factor given in \code{"group"}, based on the
#'  output of \code{\link{coef.brmsfit}} and \code{\link{ranef.brmsfit}},
#'  respectively.
#' @param seed A single numeric value passed to \code{\link{set.seed}}
#'  to make results reproducible.
#' @param ... Currently ignored.
#'
#' @details Among others, \code{hypothesis} computes an evidence ratio
#'   (\code{Evid.Ratio}) for each hypothesis. For a one-sided hypothesis, this
#'   is just the posterior probability (\code{Post.Prob}) under the hypothesis
#'   against its alternative. That is, when the hypothesis is of the form
#'   \code{a > b}, the evidence ratio is the ratio of the posterior probability
#'   of \code{a > b} and the posterior probability of \code{a < b}. In this
#'   example, values greater than one indicate that the evidence in favor of
#'   \code{a > b} is larger than evidence in favor of \code{a < b}. For an
#'   two-sided (point) hypothesis, the evidence ratio is a Bayes factor between
#'   the hypothesis and its alternative computed via the Savage-Dickey density
#'   ratio method. That is the posterior density at the point of interest
#'   divided by the prior density at that point. Values greater than one
#'   indicate that evidence in favor of the point hypothesis has increased after
#'   seeing the data. In order to calculate this Bayes factor, all parameters
#'   related to the hypothesis must have proper priors and argument
#'   \code{sample_prior} of function \code{brm} must be set to \code{"yes"}.
#'   Otherwise \code{Evid.Ratio} (and \code{Post.Prob}) will be \code{NA}.
#'   Please note that, for technical reasons, we cannot sample from priors of
#'   certain parameters classes. Most notably, these include overall intercept
#'   parameters (prior class \code{"Intercept"}) as well as group-level
#'   coefficients. When interpreting Bayes factors, make sure that your priors
#'   are reasonable and carefully chosen, as the result will depend heavily on
#'   the priors. In particular, avoid using default priors.
#'
#'   The \code{Evid.Ratio} may sometimes be \code{0} or \code{Inf} implying very
#'   small or large evidence, respectively, in favor of the tested hypothesis.
#'   For one-sided hypotheses pairs, this basically means that all posterior
#'   draws are on the same side of the value dividing the two hypotheses. In
#'   that sense, instead of \code{0} or \code{Inf,} you may rather read it as
#'   \code{Evid.Ratio} smaller \code{1 / S} or greater \code{S}, respectively,
#'   where \code{S} denotes the number of posterior draws used in the
#'   computations.
#'
#'   The argument \code{alpha} specifies the size of the credible interval
#'   (i.e., Bayesian confidence interval). For instance, if we tested a
#'   two-sided hypothesis and set \code{alpha = 0.05} (5\%) an, the credible
#'   interval will contain \code{1 - alpha = 0.95} (95\%) of the posterior
#'   values. Hence, \code{alpha * 100}\% of the posterior values will
#'   lie outside of the credible interval. Although this allows testing of
#'   hypotheses in a similar manner as in the frequentist null-hypothesis
#'   testing framework, we strongly argue against using arbitrary cutoffs (e.g.,
#'   \code{p < .05}) to determine the 'existence' of an effect.
#'
#' @return A \code{\link{brmshypothesis}} object.
#'
#' @seealso \code{\link{brmshypothesis}}
#'
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#'
#' @examples
#' \dontrun{
#' ## define priors
#' prior <- c(set_prior("normal(0,2)", class = "b"),
#'            set_prior("student_t(10,0,1)", class = "sigma"),
#'            set_prior("student_t(10,0,1)", class = "sd"))
#'
#' ## fit a linear mixed effects models
#' fit <- brm(time ~ age + sex + disease + (1 + age|patient),
#'            data = kidney, family = lognormal(),
#'            prior = prior, sample_prior = "yes",
#'            control = list(adapt_delta = 0.95))
#'
#' ## perform two-sided hypothesis testing
#' (hyp1 <- hypothesis(fit, "sexfemale = age + diseasePKD"))
#' plot(hyp1)
#' hypothesis(fit, "exp(age) - 3 = 0", alpha = 0.01)
#'
#' ## perform one-sided hypothesis testing
#' hypothesis(fit, "diseasePKD + diseaseGN - 3 < 0")
#'
#' hypothesis(fit, "age < Intercept",
#'            class = "sd", group  = "patient")
#'
#' ## test the amount of random intercept variance on all variance
#' h <- paste("sd_patient__Intercept^2 / (sd_patient__Intercept^2 +",
#'            "sd_patient__age^2 + sigma^2) = 0")
#' (hyp2 <- hypothesis(fit, h, class = NULL))
#' plot(hyp2)
#'
#' ## test more than one hypothesis at once
#' h <- c("diseaseGN = diseaseAN", "2 * diseaseGN - diseasePKD = 0")
#' (hyp3 <- hypothesis(fit, h))
#' plot(hyp3, ignore_prior = TRUE)
#'
#' ## compute hypotheses for all levels of a grouping factor
#' hypothesis(fit, "age = 0", scope = "coef", group = "patient")
#'
#' ## use the default method
#' dat <- as.data.frame(fit)
#' str(dat)
#' hypothesis(dat, "b_age > 0")
#' }
#'
#' @export
hypothesis.brmsfit <- function(x, hypothesis, class = "b", group = "",
                               scope = c("standard", "ranef", "coef"),
                               alpha = 0.05, robust = FALSE, seed = NULL,
                               ...) {
  # use a seed as prior_draws.brmsfit randomly permutes draws
  if (!is.null(seed)) {
    set.seed(seed)
  }
  contains_draws(x)
  x <- restructure(x)
  group <- as_one_character(group)
  scope <- match.arg(scope)
  if (scope == "standard") {
    if (!length(class)) {
      class <- ""
    }
    class <- as_one_character(class)
    if (nzchar(group)) {
      class <- paste0(class, "_", group, "__")
    } else if (nzchar(class)) {
      class <- paste0(class, "_")
    }
    out <- .hypothesis(
      x, hypothesis, class = class, alpha = alpha,
      robust = robust, ...
    )
  } else {
    co <- do_call(scope, list(x, summary = FALSE))
    if (!group %in% names(co)) {
      stop2("'group' should be one of ", collapse_comma(names(co)))
    }
    out <- hypothesis_coef(
      co[[group]], hypothesis, alpha = alpha,
      robust = robust, ...
    )
  }
  out
}

#' @rdname hypothesis.brmsfit
#' @export
hypothesis <- function(x, ...) {
  UseMethod("hypothesis")
}

#' @rdname hypothesis.brmsfit
#' @export
hypothesis.default <- function(x, hypothesis, alpha = 0.05,
                               robust = FALSE, ...) {
  x <- as.data.frame(x)
  .hypothesis(
    x, hypothesis, class = "", alpha = alpha,
    robust = robust, ...
  )
}

#' Descriptions of \code{brmshypothesis} Objects
#'
#' A \code{brmshypothesis} object contains posterior draws
#' as well as summary statistics of non-linear hypotheses as
#' returned by \code{\link{hypothesis}}.
#'
#' @name brmshypothesis
#'
#' @param ignore_prior A flag indicating if prior distributions
#'  should also be plotted. Only used if priors were specified on
#'  the relevant parameters.
#' @param digits Minimal number of significant digits,
#'   see \code{\link[base:print.default]{print.default}}.
#' @param chars Maximum number of characters of each hypothesis
#'  to print or plot. If \code{NULL}, print the full hypotheses.
#'  Defaults to \code{20}.
#' @param colors Two values specifying the colors of the posterior
#'  and prior density respectively. If \code{NULL} (the default)
#'  colors are taken from the current color scheme of
#'  the \pkg{bayesplot} package.
#' @param ... Currently ignored.
#' @inheritParams plot.brmsfit
#'
#' @details
#' The two most important elements of a \code{brmshypothesis} object are
#' \code{hypothesis}, which is a data.frame containing the summary estimates
#' of the hypotheses, and \code{samples}, which is a data.frame containing
#' the corresponding posterior draws.
#'
#' @seealso \code{\link{hypothesis}}
NULL

# internal function to evaluate hypotheses
# @param x the primary object passed to the hypothesis method;
#   needs to be a brmsfit object or coercible to a data.frame
# @param hypothesis vector of character strings containing the hypotheses
# @param class prefix of the parameters in the hypotheses
# @param alpha the 'alpha-level' as understood by frequentist statistics
# @return a 'brmshypothesis' object
.hypothesis <- function(x, hypothesis, class, alpha, robust,
                        combine = TRUE, ...) {
  if (!is.character(hypothesis) || !length(hypothesis)) {
    stop2("Argument 'hypothesis' must be a character vector.")
  }
  if (length(alpha) != 1L || alpha < 0 || alpha > 1) {
    stop2("Argument 'alpha' must be a single value in [0,1].")
  }
  class <- as_one_character(class)
  robust <- as_one_logical(robust)
  out <- vector("list", length(hypothesis))
  for (i in seq_along(out)) {
    out[[i]] <- eval_hypothesis(
      hypothesis[i], x = x, class = class,
      alpha = alpha, robust = robust,
      name = names(hypothesis)[i]
    )
  }
  if (combine) {
    out <- combine_hlist(out, class = class, alpha = alpha)
  }
  out
}

# evaluate hypotheses for an arrary of ranefs or coefs
# seperaly for each grouping-factor level
hypothesis_coef <- function(x, hypothesis, alpha, ...) {
  stopifnot(is.array(x), length(dim(x)) == 3L)
  levels <- dimnames(x)[[2]]
  coefs <- dimnames(x)[[3]]
  x <- lapply(seq_along(levels), function(l)
    structure(as.data.frame(x[, l, ]), names = coefs)
  )
  out <- vector("list", length(levels))
  for (l in seq_along(levels)) {
    out[[l]] <- .hypothesis(
      x[[l]], hypothesis, class = "",
      alpha = alpha, combine = FALSE, ...
    )
    for (i in seq_along(out[[l]])) {
      out[[l]][[i]]$summary$Group <- levels[l]
    }
  }
  out <- unlist(out, recursive = FALSE)
  out <- as.list(matrix(out, ncol = length(hypothesis), byrow = TRUE))
  out <- combine_hlist(out, class = "", alpha = alpha)
  out$hypothesis$Group <- factor(out$hypothesis$Group, levels)
  out$hypothesis <- move2start(out$hypothesis, "Group")
  out
}

# combine list of outputs of eval_hypothesis
# @param hlist list of evaluate hypothesis
# @return a 'brmshypothesis' object
combine_hlist <- function(hlist, class, alpha) {
  stopifnot(is.list(hlist))
  hs <- do_call(rbind, lapply(hlist, function(h) h$summary))
  rownames(hs) <- NULL
  samples <- lapply(hlist, function(h) h$samples)
  samples <- as.data.frame(do_call(cbind, samples))
  prior_samples <- lapply(hlist, function(h) h$prior_samples)
  prior_samples <- as.data.frame(do_call(cbind, prior_samples))
  names(samples) <- names(prior_samples) <- paste0("H", seq_along(hlist))
  class <- sub("_+$", "", class)
  # TODO: rename 'samples' to 'draws' in brms 3.0
  out <- nlist(hypothesis = hs, samples, prior_samples, class, alpha)
  structure(out, class = "brmshypothesis")
}

# evaluate a single hypothesis based on the posterior draws
eval_hypothesis <- function(h, x, class, alpha, robust, name = NULL) {
  stopifnot(length(h) == 1L && is.character(h))
  pars <- variables(x)[grepl(paste0("^", class), variables(x))]
  # parse hypothesis string
  h <- gsub("[ \t\r\n]", "", h)
  sign <- get_matches("=|<|>", h)
  lr <- get_matches("[^=<>]+", h)
  if (length(sign) != 1L || length(lr) != 2L) {
    stop2("Every hypothesis must be of the form 'left (= OR < OR >) right'.")
  }
  h <- paste0("(", lr[1], ")")
  h <- paste0(h, ifelse(lr[2] != "0", paste0("-(", lr[2], ")"), ""))
  varsH <- find_vars(h)
  parsH <- paste0(class, varsH)
  miss_pars <- setdiff(parsH, pars)
  if (length(miss_pars)) {
    miss_pars <- collapse_comma(miss_pars)
    stop2("Some parameters cannot be found in the model: \n", miss_pars)
  }
  # rename hypothesis for correct evaluation
  h_renamed <- rename(h, c(":", "[", "]", ","),  c("___", ".", ".", ".."))
  # get posterior and prior draws
  pattern <- c(paste0("^", class), ":", "\\[", "\\]", ",")
  repl <- c("", "___", ".", ".", "..")
  samples <- as.data.frame(x, variable = parsH)
  names(samples) <- rename(names(samples), pattern, repl, fixed = FALSE)
  samples <- as.matrix(eval2(h_renamed, samples))
  prior_samples <- prior_draws(x, variable = parsH)
  if (!is.null(prior_samples) && ncol(prior_samples) == length(varsH)) {
    names(prior_samples) <- rename(
      names(prior_samples), pattern, repl, fixed = FALSE
    )
    prior_samples <- as.matrix(eval2(h_renamed, prior_samples))
  } else {
    prior_samples <- NULL
  }
  # summarize hypothesis
  wsign <- switch(sign, "=" = "equal", "<" = "less", ">" = "greater")
  probs <- switch(sign,
                  "=" = c(alpha / 2, 1 - alpha / 2),
                  "<" = c(alpha, 1 - alpha), ">" = c(alpha, 1 - alpha)
  )
  if (robust) {
    measures <- c("median", "mad")
  } else {
    measures <- c("mean", "sd")
  }
  measures <- c(measures, "quantile", "evidence_ratio")
  sm <- lapply(
    measures, get_estimate, draws = samples, probs = probs,
    wsign = wsign, prior_samples = prior_samples
  )
  sm <- as.data.frame(matrix(unlist(sm), nrow = 1))
  names(sm) <- c("Estimate", "Est.Error", "CI.Lower", "CI.Upper", "Evid.Ratio")
  sm$Post.Prob <- sm$Evid.Ratio / (1 + sm$Evid.Ratio)
  if (is.infinite(sm$Evid.Ratio)) {
    sm$Post.Prob <- 1
  }
  if (sign == "=") {
    sm$Star <- str_if(!(sm$CI.Lower <= 0 && 0 <= sm$CI.Upper), "*")
  } else {
    sm$Star <- str_if(sm$Post.Prob > 1 - alpha, "*")
  }
  if (!length(name) || !nzchar(name)) {
    name <- paste(h, sign, "0")
  }
  sm$Hypothesis <- as_one_character(name)
  sm <- move2start(sm, "Hypothesis")
  if (is.null(prior_samples)) {
    prior_samples <- as.matrix(rep(NA, nrow(samples)))
  }
  nlist(summary = sm, samples, prior_samples)
}

# find all valid variable names in a string
# @param x a character string
# @param dot are dots allowed in variable names?
# @param brackets allow brackets at the end of variable names?
# @return all valid variable names within the string
# @note does not use the R parser itself to allow for double points,
#   square brackets, and commas at the end of names
find_vars <- function(x, dot = TRUE, brackets = TRUE) {
  x <- gsub("[[:space:]]", "", as_one_character(x))
  dot <- as_one_logical(dot)
  brackets <- as_one_logical(brackets)
  regex_all <- paste0(
    "([^([:digit:]|[:punct:])]", if (dot) "|\\.", ")",
    "[[:alnum:]_\\:", if (dot) "\\.", "]*",
    if (brackets) "(\\[[^],]+(,[^],]+)*\\])?"
  )
  pos_all <- gregexpr(regex_all, x)[[1]]
  regex_fun <- paste0(
    "([^([:digit:]|[:punct:])]", if (dot) "|\\.", ")",
    "[[:alnum:]_", if (dot) "\\.", "]*\\("
  )
  pos_fun <- gregexpr(regex_fun, x)[[1]]
  pos_decnum <- gregexpr("\\.[[:digit:]]+", x)[[1]]
  keep <- !pos_all %in% c(pos_fun, pos_decnum)
  pos_var <- pos_all[keep]
  attr(pos_var, "match.length") <- attributes(pos_all)$match.length[keep]
  if (length(pos_var)) {
    out <- unique(unlist(regmatches(x, list(pos_var))))
  } else {
    out <- character(0)
  }
  out
}

#' Compute Density Ratios
#'
#' Compute the ratio of two densities at given points based on draws of the
#' corresponding distributions.
#'
#' @param x Vector of draws from the first distribution, usually the posterior
#'   distribution of the quantity of interest.
#' @param y Optional vector of draws from the second distribution, usually the
#'   prior distribution of the quantity of interest. If \code{NULL} (the
#'   default), only the density of \code{x} will be evaluated.
#' @param point Numeric values at which to evaluate and compare the densities.
#'   Defaults to \code{0}.
#' @param n Single numeric value. Influences the accuracy of the density
#'   estimation. See \code{\link[stats:density]{density}} for details.
#' @param ... Further arguments passed to \code{\link[stats:density]{density}}.
#'
#' @return A vector of length equal to \code{length(point)}. If \code{y} is
#'   provided, the density ratio of \code{x} against \code{y} is returned. Else,
#'   only the density of \code{x} is returned.
#'
#' @details In order to achieve sufficient accuracy in the density estimation,
#'   more draws than usual are required. That is you may need an effective
#'   sample size of 10,000 or more to reliably estimate the densities.
#'
#' @examples
#' x <- rnorm(10000)
#' y <- rnorm(10000, mean = 1)
#' density_ratio(x, y, point = c(0, 1))
#'
#' @export
density_ratio <- function(x, y = NULL, point = 0, n = 4096, ...) {
  x <- as.numeric(x)
  point <- as.numeric(point)
  dots <- list(...)
  dots <- dots[names(dots) %in% names(formals("density.default"))]
  dots$n <- n

  eval_density <- function(x, point) {
    # evaluate density of x at point
    from <- min(x)
    to <- max(x)
    if (from > point) {
      from <- point - sd(x) / 4
    } else if (to < point) {
      to <- point + sd(x) / 4
    }
    density <- KernSmooth::bkde(x, range.x = c(from,to), gridsize = n)
    stats::approx(density$x, density$y, xout = point)$y
  }

  out <- ulapply(point, eval_density, x = x)
  if (!is.null(y)) {
    y <- as.numeric(y)
    out <- out / ulapply(point, eval_density, x = y)
  }
  out
}

# compute the evidence ratio between two disjunct hypotheses
# @param x posterior draws
# @param cut the cut point between the two hypotheses
# @param wsign direction of the hypothesis
# @param prior_samples optional prior draws for two-sided hypothesis
# @param ... optional arguments passed to density_ratio
# @return the evidence ratio of the two hypothesis
evidence_ratio <- function(x, cut = 0, wsign = c("equal", "less", "greater"),
                           prior_samples = NULL, ...) {
  wsign <- match.arg(wsign)
  if (wsign == "equal") {
    if (is.null(prior_samples)) {
      out <- NA
    } else {
      out <- density_ratio(x, prior_samples, point = cut, ...)
    }
  } else if (wsign == "less") {
    out <- length(which(x < cut))
    out <- out / (length(x) - out)
  } else if (wsign == "greater") {
    out <- length(which(x > cut))
    out <- out / (length(x) - out)
  }
  out
}

# round all numeric elements of a list-like object
round_numeric <- function(x, digits = 2) {
  stopifnot(is.list(x))
  for (i in seq_along(x)) {
    if (is.numeric(x[[i]])) {
      x[[i]] <- round(x[[i]], digits = digits)
    }
  }
  x
}

#' @rdname brmshypothesis
#' @export
print.brmshypothesis <- function(x, digits = 2, chars = 20, ...) {
  # make sure hypothesis names are not too long
  x$hypothesis$Hypothesis <- limit_chars(
    x$hypothesis$Hypothesis, chars = chars
  )
  cat(paste0("Hypothesis Tests for class ", x$class, ":\n"))
  x$hypothesis <- round_numeric(x$hypothesis, digits = digits)
  print(x$hypothesis, quote = FALSE)
  pone <- (1 - x$alpha * 2) * 100
  ptwo <- (1 - x$alpha) * 100
  cat(glue(
    "---\n'CI': {pone}%-CI for one-sided and {ptwo}%-CI for two-sided hypotheses.\n",
    "'*': For one-sided hypotheses, the posterior probability exceeds {ptwo}%;\n",
    "for two-sided hypotheses, the value tested against lies outside the {ptwo}%-CI.\n",
    "Posterior probabilities of point hypotheses assume equal prior probabilities.\n"
  ))
  invisible(x)
}

#' @rdname brmshypothesis
#' @method plot brmshypothesis
#' @export
plot.brmshypothesis <- function(x, nvariables = 5, N = NULL,
                                ignore_prior = FALSE, chars = 40, colors = NULL,
                                theme = NULL, ask = TRUE, plot = TRUE,  ...) {
  dots <- list(...)
  nvariables <- use_alias(nvariables, N)
  if (!is.data.frame(x$samples)) {
    stop2("No posterior draws found.")
  }
  plot <- use_alias(plot, dots$do_plot)
  if (is.null(colors)) {
    colors <- bayesplot::color_scheme_get()[c(4, 2)]
    colors <- unname(unlist(colors))
  }
  if (length(colors) != 2L) {
    stop2("Argument 'colors' must be of length 2.")
  }

  .plot_fun <- function(samples) {
    samples <- na.omit(samples)
    # if no prior draws are present, there is no need to plot a legend
    ignore_prior <- ignore_prior || length(unique(samples$Type)) == 1L
    gg <- ggplot(samples, aes(x = .data[["values"]])) +
      facet_wrap("ind", ncol = 1, scales = "free") +
      xlab("") + ylab("") + theme +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    if (ignore_prior) {
      gg <- gg +
        geom_density(alpha = 0.7, fill = colors[1], na.rm = TRUE)
    } else {
      gg <- gg +
        geom_density(aes(fill = .data[["Type"]]), alpha = 0.7, na.rm = TRUE) +
        scale_fill_manual(values = colors)
    }
    return(gg)
  }

  samples <- cbind(x$samples, Type = "Posterior")
  if (!ignore_prior) {
    prior_samples <- cbind(x$prior_samples, Type = "Prior")
    samples <- rbind(samples, prior_samples)
  }
  if (plot) {
    default_ask <- devAskNewPage()
    on.exit(devAskNewPage(default_ask))
    devAskNewPage(ask = FALSE)
  }
  hyps <- limit_chars(x$hypothesis$Hypothesis, chars = chars)
  if (!is.null(x$hypothesis$Group)) {
    hyps <- paste0(x$hypothesis$Group, ":  ", hyps)
  }
  names(samples)[seq_along(hyps)] <- hyps
  nplots <- ceiling(length(hyps) / nvariables)
  plots <- vector(mode = "list", length = nplots)
  for (i in seq_len(nplots)) {
    sub <- ((i - 1) * nvariables + 1):min(i * nvariables, length(hyps))
    sub_hyps <- hyps[sub]
    sub_samples <- cbind(
      utils::stack(samples[, sub_hyps, drop = FALSE]),
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
