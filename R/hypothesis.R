#' Non-Linear Hypothesis Testing
#' 
#' Perform non-linear hypothesis testing for all model parameters. 
#' 
#' @param x An \code{R} object. If it is no \code{brmsfit} object,
#'  it must be coercible to a \code{data.frame}.
#' @param hypothesis A character vector specifying one or more 
#'  non-linear hypothesis concerning parameters of the model.
#' @param class A string specifying the class of parameters being tested. 
#'  Default is "b" for population-level effects. 
#'  Other typical options are "sd" or "cor". 
#'  If \code{class = NULL}, all parameters can be tested
#'  against each other, but have to be specified with their full name 
#'  (see also \code{\link[brms:parnames]{parnames}}) 
#' @param group Name of a grouping factor to evaluate only 
#'  group-level effects parameters related to this grouping factor.
#'  Ignored if \code{class} is not \code{"sd"} or \code{"cor"}.
#' @param alpha The alpha-level of the tests (default is 0.05;
#'  see 'Details' for more information).
#' @param scope Indicates where to look for the variables specified in
#'  \code{hypothesis}. If \code{"standard"}, use the full parameter names
#'  (subject to the restriction given by \code{class}). If \code{"coef"}
#'  or \code{"ranef"} compute the hypothesis for all levels of 
#'  the grouping factor given in \code{"group"}, based on the 
#'  output of \code{\link{coef.brmsfit}} and \code{\link{ranef.brmsfit}},
#'  respectively.
#' @param seed A single numeric value passed to \code{\link{set.seed}}
#'  to make results reproducible.
#' @param ... Currently ignored.
#' 
#' @details Among others, \code{hypothesis} computes an 
#'  evidence ratio (\code{Evid.Ratio}) for each hypothesis. 
#'  For a directed hypothesis, this is just the posterior probability 
#'  under the hypothesis against its alternative.
#'  That is, when the hypothesis if of the form \code{a > b}, 
#'  the evidence ratio is the ratio of the posterior probability 
#'  of \code{a > b} and the posterior probability of \code{a < b}.
#'  In this example, values greater than one indicate that the evidence in
#'  favor of \code{a > b} is larger than evidence in favor of \code{a < b}.
#'  For an undirected (point) hypothesis, the evidence ratio 
#'  is a Bayes factor between the hypothesis and its alternative
#'  computed via the Savage-Dickey density ratio method.
#'  That is the posterior density at the point of interest divided
#'  by the prior density at that point.
#'  Values greater than one indicate that evidence in favor of the point
#'  hypothesis has increased after seeing the data.
#'  In order to calculate this Bayes factor, all parameters related 
#'  to the hypothesis must have proper priors
#'  and argument \code{sample_prior} of function \code{brm} 
#'  must be set to \code{TRUE}. 
#'  When interpreting Bayes factors, make sure 
#'  that your priors are reasonable and carefully chosen,
#'  as the result will depend heavily on the priors. 
#'  In particular, avoid using default priors.
#'  
#'  The argument \code{alpha} specifies the size of the credible interval
#'  (i.e., Bayesian confidence interval).
#'  For instance, if \code{alpha = 0.05} (5\%), the credible interval
#'  will contain \code{1 - alpha = 0.95} (95\%) of the posterior values.
#'  Hence, \code{alpha * 100}\% of the posterior values will lie
#'  outside of the credible interval. Although this allows testing of
#'  hypotheses in a similar manner as in the frequentist null-hypothesis
#'  testing framework, we strongly argue against using arbitrary cutoffs 
#'  (e.g., \code{p < .05}) to determine the 'existence' of an effect.
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
#'            prior = prior, sample_prior = TRUE, 
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
#' hypothesis(dat, "b_age > 0")
#' }
#' 
#' @export
hypothesis <- function(x, ...) {
  UseMethod("hypothesis")
}

#' @rdname hypothesis
#' @export
hypothesis.default <- function(x, hypothesis, alpha = 0.05, ...) {
  x <- as.data.frame(x)
  hypothesis_internal(x, hypothesis, class = "", alpha = alpha, ...)
}

hypothesis_internal <- function(x, hypothesis, class, alpha,
                                combine = TRUE, ...) {
  # internal function to evaluate hypotheses
  # Args:
  #   x: the primary object passed to the hypothesis method;
  #     Needs to be a brmsfit object or coercible to a data.frame
  #   hypothesis: Vector of character strings containing the hypotheses
  #   class: prefix of the parameters in the hypotheses
  #   alpha: alpha-level
  # Returns:
  #   an object of class 'brmshypothesis'
  if (!is.character(hypothesis)) {
    stop2("Argument 'hypothesis' must be a character vector.")
  }
  if (length(alpha) != 1L || alpha < 0 || alpha > 1) {
    stop2("Argument 'alpha' must be a single value in [0,1].")
  }
  class <- as_one_character(class)
  out <- vector("list", length(hypothesis))
  for (i in seq_along(out)) {
    out[[i]] <- eval_hypothesis(
      hypothesis[i], x = x, class = class, 
      alpha = alpha, name = names(hypothesis)[i]
    )
  }
  if (combine) {
    out <- combine_hlist(out, class = class, alpha = alpha)
  }
  out
}

hypothesis_coef <- function(x, hypothesis, alpha, ...) {
  # evaluate hypotheses for an arrary of ranefs or coefs
  # seperaly for each grouping-factor level
  stopifnot(is.array(x), length(dim(x)) == 3L)
  levels <- dimnames(x)[[2]]
  coefs <- dimnames(x)[[3]]
  x <- lapply(seq_along(levels), function(l)
    structure(as.data.frame(x[, l, ]), names = coefs)
  )
  out <- vector("list", length(levels))
  for (l in seq_along(levels)) {
    out[[l]] <- hypothesis_internal(
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
  out$hypothesis <- move2start(out$hypothesis, "Group")
  out
}

combine_hlist <- function(hlist, class, alpha) {
  # combine list of outputs of eval_hypothesis
  # Returns: a brmshypothesis object
  stopifnot(is.list(hlist))
  hs <- do.call(rbind, lapply(hlist, function(h) h$summary))
  rownames(hs) <- NULL
  samples <- lapply(hlist, function(h) h$samples)
  samples <- as.data.frame(do.call(cbind, samples))
  prior_samples <- lapply(hlist, function(h) h$prior_samples)
  prior_samples <- as.data.frame(do.call(cbind, prior_samples))
  names(samples) <- names(prior_samples) <- paste0("H", seq_along(hlist))
  class <- sub("_+$", "", class)
  out <- nlist(hypothesis = hs, samples, prior_samples, class, alpha)
  structure(out, class = "brmshypothesis")
}

eval_hypothesis <- function(h, x, class, alpha, name = NULL) {
  stopifnot(length(h) == 1L && is.character(h))
  pars <- parnames(x)[grepl(paste0("^", class), parnames(x))]
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
  # get posterior and prior samples
  pattern <- c(paste0("^", class), ":", "\\[", "\\]", ",")
  repl <- c("", "___", ".", ".", "..")
  samples <- posterior_samples(x, pars = parsH, exact_match = TRUE)
  names(samples) <- rename(names(samples), pattern, repl, fixed = FALSE)
  samples <- as.matrix(eval2(h_renamed, samples))
  prior_samples <- prior_samples(x, pars = parsH, exact_match = TRUE)
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
    "<" = c(0, 1 - alpha), ">" = c(alpha, 1)
  )
  sm <- lapply(
    c("mean", "sd", "quantile", "evidence_ratio"), 
    get_estimate, samples = samples, probs = probs, 
    wsign = wsign, prior_samples = prior_samples
  )
  sm <- as.data.frame(matrix(unlist(sm), nrow = 1))
  names(sm) <- c("Estimate", "Est.Error", "CI.Lower", "CI.Upper", "Evid.Ratio")
  if (sign == "<") {
    sm[1, 3] <- -Inf
  } else if (sign == ">") {
    sm[1, 4] <- Inf
  }
  sm$Star <- ifelse(!(sm[1, 3] <= 0 && 0 <= sm[1, 4]), '*', '')
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

find_vars <- function(x) {
  # find all valid variable names in a string 
  # Args:
  #   x: a character string
  # Notes:
  #   Does not use the R parser itself to allow for double points, 
  #   square brackets and commas at the end of names.
  #   currently only used in 'hypothesis_internal'
  # Returns:
  #   all valid variable names within the string
  x <- gsub("[[:space:]]", "", as_one_character(x))
  regex_all <- "([^([:digit:]|[:punct:])]|\\.)[[:alnum:]_\\.\\:]*"
  regex_all <- paste0(regex_all, "(\\[[^],]+(,[^],]+)*\\])?")
  pos_all <- gregexpr(regex_all, x)[[1]]
  regex_fun <- "([^([:digit:]|[:punct:])]|\\.)[[:alnum:]_\\.]*\\("
  pos_fun <- gregexpr(regex_fun, x)[[1]]
  pos_decnum <- gregexpr("\\.[[:digit:]]+", x)[[1]]
  pos_var <- list(rmMatch(pos_all, pos_fun, pos_decnum))
  unique(unlist(regmatches(x, pos_var)))
}

evidence_ratio <- function(x, cut = 0, wsign = c("equal", "less", "greater"), 
                           prior_samples = NULL, pow = 12, ...) {
  # compute the evidence ratio between two disjunct hypotheses
  # Args:
  #   x: posterior samples 
  #   cut: the cut point between the two hypotheses
  #   wsign: direction of the hypothesis
  #   prior_samples: optional prior samples for undirected hypothesis
  #   pow: influences the accuracy of the density
  #   ...: optional arguments passed to density.default
  # Returns:
  #   the evidence ratio of the two hypothesis
  wsign <- match.arg(wsign)
  if (wsign == "equal") {
    if (is.null(prior_samples)) {
      out <- NA
    } else {
      dots <- list(...)
      dots <- dots[names(dots) %in% names(formals("density.default"))]
      args <- c(list(n = 2^pow), dots)
      eval_dens <- function(x) {
        # evaluate density of x at cut
        from <- min(x)
        to <- max(x)
        if (from > cut) {
          from <- cut - sd(x) / 4
        } else if (to < cut) {
          to <- cut + sd(x) / 4
        }
        dens <- do.call(density, c(nlist(x, from, to), args))
        spline(dens$x, dens$y, xout = cut)$y
      }
      out <- eval_dens(x) / eval_dens(prior_samples)
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

round_numeric <- function(x, digits = 2) {
  # round all numeric elements of a list like object
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
  cat(paste0(
    "---\n'*': The expected value under the hypothesis ", 
    "lies outside the ", (1 - x$alpha) * 100, "%-CI.\n"
  ))
  invisible(x)
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
  hyps <- limit_chars(x$hypothesis$Hypothesis, chars = chars)
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
