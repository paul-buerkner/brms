#' Bayesian Regression Models using Stan
#'
#' @docType package
#' @name brms-package
#' @aliases brms
#'
#' @description
#' \if{html}{
#'    \figure{stanlogo.png}{options: width="50px" alt="http://mc-stan.org/about/logo/"}
#'    \emph{Stan Development Team}
#' }
#' 
#' The \pkg{brms} package provides an interface to fit Bayesian generalized
#' (non-)linear multilevel models using \pkg{Stan}, which is a C++ package 
#' for obtaining full Bayesian inference (see \url{http://mc-stan.org/}). 
#' The formula syntax is an extended version of the syntax applied in the 
#' \pkg{lme4} package to provide a familiar and simple interface for 
#' performing regression analyses.
#'  
#' @details 
#' The main function of the brms package is \code{\link{brm}}, which creates the
#' model in Stan language and fits it using \pkg{\link[rstan:rstan]{Stan}}.
#' Subsequently, a large number of methods can be applied:
#' To get an overview on the estimated parameters, 
#' \code{\link[brms:summary.brmsfit]{summary}} or 
#' \code{\link[brms:marginal_effects.brmsfit]{marginal_effects}} 
#' are perfectly suited. 
#' Detailed visual analyses can be performed by applying the
#' \pkg{\link[shinystan:shinystan-package]{shinystan}} package, which can be
#' called directly within \pkg{brms} using
#' \code{\link[brms:launch_shinystan]{launch_shinystan}}. Information Criteria
#' are also readily available via \code{\link[brms:LOO]{LOO}} and
#' \code{\link[brms:WAIC]{WAIC}}, both relying on the
#' \pkg{\link[loo:loo-package]{loo}} package. For a full list of methods to
#' apply, type \code{methods(class = "brmsfit")}.
#' 
#' Because \pkg{brms} is based on \pkg{Stan}, a C++ compiler is required. 
#' The program Rtools (available on \url{https://cran.r-project.org/bin/windows/Rtools/}) 
#' comes with a C++ compiler for Windows. On Mac, you should use Xcode.
#' For further instructions on how to get the compilers running, see the
#' prerequisites section at the
#' \href{https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started}{RStan-Getting-Started}
#' page.
#' 
#' When comparing other packages fitting GLMMs to \pkg{brms}, keep in mind that
#' the latter needs to compile models before actually fitting them, which will
#' require between 20 and 40 seconds depending on your machine, operating system
#' and overall model complexity.
#' 
#' Thus, fitting smaller models may be relatively slow as compilation time makes
#' up the majority of the whole running time. For larger / more complicated
#' models however, fitting my take several minutes or even hours, so that the
#' compilation time won't make much of a difference here.
#' 
#' See \code{vignette("brms_overview")} for a general introduction and overview
#' of \pkg{brms}. For a full list of available vignettes, type
#' \code{vignette(package = "brms")}.
#' 
#' @references 
#' Paul-Christian Buerkner (2017). brms: An R Package for Bayesian Multilevel 
#' Models Using Stan. \emph{Journal of Statistical Software}, 80(1), 1-28. 
#' \code{doi:10.18637/jss.v080.i01}
#' 
#' The Stan Development Team. \emph{Stan Modeling Language User's Guide and
#' Reference Manual}. \url{http://mc-stan.org/users/documentation}.
#' 
#' Stan Development Team (2018). RStan: the R interface to Stan. R package
#' version 2.18.1. \url{http://mc-stan.org}
#' 
#' @seealso 
#' \code{\link[brms:brm]{brm}}, 
#' \code{\link[brms:brmsformula]{brmsformula}}, 
#' \code{\link[brms:brmsfamily]{brmsfamily}},
#' \code{\link[brms:brmsfit-class]{brmsfit}}
#'
NULL
