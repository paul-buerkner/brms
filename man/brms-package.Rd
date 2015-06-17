\name{brms-package}
\alias{brms-package}
\alias{brms}
\docType{package}
\title{
Bayesian Regression Models using Stan
}
\description{
The brms package provides an interface to fit Bayesian generalized linear mixed models using Stan,
which is a C++ package for obtaining Bayesian inference using the No-U-turn sampler (see \url{http://mc-stan.org/}). The formula syntax is very similar to that of the package lme4 to provide a familiar and simple interface for performing regression analyses.

}
\details{
%\tabular{ll}{
%Package: \tab brms\cr
%Type: \tab Package\cr
%Version: \tab 0.0.0.9000\cr
%Date: \tab 2015-04-02\cr
%License: \tab GPL >= 2 \cr
%}

The main function of the brms package is \code{\link[brms:brm]{brm}}, which creates the model in Stan language and fits it using the package rstan, if installed. Subsequently, methods suchs as \code{\link[brms:summary.brmsfit]{summary}} or \code{\link[brms:plot.brmsfit]{plot}} can be applied to investigate the obtained results

Without having rstan installed, the function \code{brm} will not return posterior samples, but all prerequisites to fit the model in Stan, such as the model code and the relevant data. If rstan is installed, the model is fitted automatically. Unfortunately, rstan has to be installed manually, as it is not on CRAN, yet. However, the developers of Stan and rstan are currently working on a version to be uploaded on CRAN. In the meantime, instructions on how to install rstan can be found at \url{https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started}.  

%First, you need a C++ compiler. See \url{https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started#prerequisites} for instructions on how to get one. Second, install rstan by running the following R code (the number behind 'j' in the first line corresponds to the number of cores to use for the installation). This may take a few minutes and you should restart R after the installation.

%\code{Sys.setenv(MAKEFLAGS = "-j1")} \cr
%\code{source('http://mc-stan.org/rstan/install.R', echo = TRUE, max.deparse.length = 2000)} \cr
%\code{install_rstan()}
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

%\examples{
%## Not run: 
%fit <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1|patient), data = epilepsy, 
%           family = c("poisson", "log"))
%summary(fit)
%plot(fit)       
%
%## End(Not run)
%}
