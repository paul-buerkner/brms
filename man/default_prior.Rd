% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/priors.R
\name{default_prior}
\alias{default_prior}
\alias{get_prior}
\title{Default priors for Bayesian models}
\usage{
default_prior(object, ...)

get_prior(formula, ...)
}
\arguments{
\item{object}{An object whose class will determine which method will
be used. A symbolic description of the model to be fitted.}

\item{...}{Further arguments passed to the specific method.}

\item{formula}{Synonym of \code{object} for use in \code{get_prior}.}
}
\value{
Usually, a \code{brmsprior} object. See
  \code{\link{default_prior.default}} for more details.
}
\description{
\code{default_prior} is a generic function that can be used to
  get default priors for Bayesian models. Its original use is
  within the \pkg{brms} package, but new methods for use
  with objects from other packages can be registered to the same generic.
}
\details{
See \code{\link{default_prior.default}} for the default method applied for
\pkg{brms} models. You can view the available methods by typing
\code{methods(default_prior)}.
}
\examples{
## get all parameters and parameters classes to define priors on
(prior <- default_prior(count ~ zAge + zBase * Trt + (1|patient) + (1|obs),
                        data = epilepsy, family = poisson()))

}
\seealso{
\code{\link{set_prior}}, \code{\link{default_prior.default}}
}
