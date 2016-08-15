\name{brms-package}
\alias{brms-package}
\alias{brms}
\docType{package}
\title{
Bayesian Regression Models using Stan
}
\description{
The \pkg{brms} package provides an interface to fit Bayesian generalized (non)-linear mixed models using Stan,
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

The main function of the brms package is \code{\link[brms:brm]{brm}}, 
which creates the model in Stan language and fits it using \pkg{Stan}. 
Subsequently, a large number of methods can be applied: 
To get an overview on the estimated parameters, 
\code{\link[brms:summary.brmsfit]{summary}} or \code{\link[brms:plot.brmsfit]{plot}} 
are perfectly suited. 
Detailed visual analyses can be performed by applying the \pkg{shinystan} package, 
which can be called directly within \pkg{brms} using 
\code{\link[brms:launch_shiny]{launch_shiny}}. 
Information Criteria are also readily available via \code{\link[brms:WAIC]{WAIC}}
and \code{\link[brms:LOO]{LOO}} both relying on the \pkg{loo} package.

Because \pkg{brms} is based on \pkg{Stan}, a C++ compiler is required. 
The program Rtools (available on https://cran.r-project.org/bin/windows/Rtools/) 
comes with a C++ compiler for Windows. On Mac, you should use Xcode.
For further instructions on how to get the compilers running, 
see the prerequisites section on 
https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started.

When comparing other packages fitting GLMMs to \pkg{brms}, keep in mind 
that the latter needs to compile models before actually fitting them, 
which will require between 20 and 40 seconds depending on your machine, 
operating system and overall model complexity.
Thus, fitting smaller models may be relatively slow as compilation time makes
up the majority of the whole running time. For larger / more complicated models 
however, fitting my take several minutes or even hours, 
so that the compilation time won't make much of a difference here.
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
\code{\link[brms:brm]{brm}}, \code{\link[brms:brmsformula]{brmsformula}}, \code{\link[brms:brmsfit-class]{brmsfit}}
}
