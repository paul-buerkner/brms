\name{brms-package}
\alias{brms-package}
\alias{brms}
\docType{package}
\title{
Bayesian Regression Models using Stan
}
\description{
The \pkg{brms} package provides an interface to fit Bayesian generalized linear mixed models using Stan,
which is a C++ package for obtaining Bayesian inference using the No-U-turn sampler (see \url{http://mc-stan.org/}). The formula syntax is very similar to that of the package \pkg{lme4} to provide a familiar and simple interface for performing regression analyses.
}

\details{
%\tabular{ll}{
%Package: \tab brms\cr
%Type: \tab Package\cr
%Version: \tab 0.0.0.9000\cr
%Date: \tab 2015-04-02\cr
%License: \tab GPL >= 2 \cr
%}

The main function of the brms package is \code{\link[brms:brm]{brm}}, which creates the model in Stan language and fits it using \pkg{Stan}. Subsequently, a large number of methods can be applied: To get an overview on the estimated parameters, \code{\link[brms:summary.brmsfit]{summary}} or \code{\link[brms:plot.brmsfit]{plot}} are perfectly suited. Detailed visual analyses can be performed by applying the \pkg{shinystan} package, which can be called directly within \pkg{brms} using \code{\link[brms:launch_shiny]{launch_shiny}}. Information Criteria are also readily available via \code{\link[brms:WAIC]{WAIC}} and \code{\link[brms:LOO]{LOO}} both relying on the \pkg{loo} package.

Because \pkg{brms} is based on \pkg{Stan}, a C++ compiler is required. The program Rtools (available on https://cran.r-project.org/bin/windows/Rtools/) comes with a C++ compiler for Windows. On Mac, you should use Xcode. We recommend to install rstan (available on CRAN) before installing \pkg{brms}. For further instructions see https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started.
}

\author{
Paul-Christian Buerkner

Maintainer: Paul-Christian Buerkner <paul.buerkner@gmail.com>
}
\references{
The Stan Development Team Stan Modeling Language User's Guide and Reference Manual. \url{http://mc-stan.org/}.
}
\keyword{package}
\seealso{
\code{\link[brms:brm]{brm}}, \code{\link[brms:brmsfit-class]{brmsfit}}
}
