#' Bayesian Regression Models using 'Stan'
#'
#' @name brms-package
#' @aliases brms
#'
#' @description
#' \if{html}{
#'    \figure{stanlogo.png}{options: width="50" alt="https://mc-stan.org/about/logo/"}
#'    \emph{Stan Development Team}
#' }
#'
#' The \pkg{brms} package provides an interface to fit Bayesian generalized
#' multivariate (non-)linear multilevel models using \pkg{Stan}, which is a C++
#' package for obtaining full Bayesian inference (see
#' \url{https://mc-stan.org/}). The formula syntax is an extended version of the
#' syntax applied in the \pkg{lme4} package to provide a familiar and simple
#' interface for performing regression analyses.
#'
#' @details
#' The main function of \pkg{brms} is \code{\link{brm}}, which uses
#' formula syntax to specify a wide range of complex Bayesian models
#' (see \code{\link{brmsformula}} for details). Based on the supplied
#' formulas, data, and additional information, it writes the Stan code
#' on the fly via \code{\link[brms:stancode.default]{stancode}}, prepares the data via
#' \code{\link[brms:standata.default]{standata}} and fits the model using
#' \pkg{\link[rstan:rstan]{Stan}}.
#'
#' Subsequently, a large number of post-processing methods can be applied:
#' To get an overview on the estimated parameters,
#' \code{\link[brms:summary.brmsfit]{summary}} or
#' \code{\link[brms:conditional_effects.brmsfit]{conditional_effects}}
#' are perfectly suited. Detailed visual analyses can be performed by applying
#' the \code{\link{pp_check}} and \code{\link{stanplot}} methods, which both
#' rely on the \pkg{\link[bayesplot:bayesplot-package]{bayesplot}} package.
#' Model comparisons can be done via \code{\link{loo}} and \code{\link{waic}},
#' which make use of the \pkg{\link[loo:loo-package]{loo}} package as well as
#' via \code{\link{bayes_factor}} which relies on the \pkg{bridgesampling} package.
#' For a full list of methods to apply, type \code{methods(class = "brmsfit")}.
#'
#' Because \pkg{brms} is based on \pkg{Stan}, a C++ compiler is required. The
#' program Rtools (available on
#' \url{https://cran.r-project.org/bin/windows/Rtools/}) comes with a C++
#' compiler for Windows. On Mac, you should use Xcode. For further instructions
#' on how to get the compilers running, see the prerequisites section at the
#' \href{https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started}{RStan-Getting-Started}
#' page.
#'
#' When comparing other packages fitting multilevel models to \pkg{brms}, keep
#' in mind that the latter needs to compile models before actually fitting them,
#' which will require between 20 and 40 seconds depending on your machine,
#' operating system and overall model complexity.
#'
#' Thus, fitting smaller models may be relatively slow as compilation time makes
#' up the majority of the whole running time. For larger / more complex
#' models however, fitting my take several minutes or even hours, so that the
#' compilation time won't make much of a difference for these models.
#'
#' See \code{vignette("brms_overview")} and \code{vignette("brms_multilevel")}
#' for a general introduction and overview of \pkg{brms}. For a full list of
#' available vignettes, type \code{vignette(package = "brms")}.
#'
#' @references
#' Paul-Christian Buerkner (2017). brms: An R Package for Bayesian Multilevel
#' Models Using Stan. \emph{Journal of Statistical Software}, 80(1), 1-28.
#' \code{doi:10.18637/jss.v080.i01}
#'
#' Paul-Christian Buerkner (2018). Advanced Bayesian Multilevel Modeling
#' with the R Package brms. \emph{The R Journal}. 10(1), 395–411.
#' \code{doi:10.32614/RJ-2018-017}
#'
#' The Stan Development Team. \emph{Stan Modeling Language User's Guide and
#' Reference Manual}. \url{https://mc-stan.org/users/documentation/}.
#'
#' Stan Development Team (2020). RStan: the R interface to Stan. R package
#' version 2.21.2. \url{https://mc-stan.org/}
#'
#' @seealso
#' \code{\link{brm}},
#' \code{\link{brmsformula}},
#' \code{\link{brmsfamily}},
#' \code{\link{brmsfit}}
#'
"_PACKAGE"
