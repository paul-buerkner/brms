% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distributions.R
\name{SkewNormal}
\alias{SkewNormal}
\alias{dskew_normal}
\alias{pskew_normal}
\alias{qskew_normal}
\alias{rskew_normal}
\title{The Skew-Normal Distribution}
\usage{
dskew_normal(
  x,
  mu = 0,
  sigma = 1,
  alpha = 0,
  xi = NULL,
  omega = NULL,
  log = FALSE
)

pskew_normal(
  q,
  mu = 0,
  sigma = 1,
  alpha = 0,
  xi = NULL,
  omega = NULL,
  lower.tail = TRUE,
  log.p = FALSE
)

qskew_normal(
  p,
  mu = 0,
  sigma = 1,
  alpha = 0,
  xi = NULL,
  omega = NULL,
  lower.tail = TRUE,
  log.p = FALSE,
  tol = 1e-08
)

rskew_normal(n, mu = 0, sigma = 1, alpha = 0, xi = NULL, omega = NULL)
}
\arguments{
\item{x, q}{Vector of quantiles.}

\item{mu}{Vector of mean values.}

\item{sigma}{Vector of standard deviation values.}

\item{alpha}{Vector of skewness values.}

\item{xi}{Optional vector of location values.
If \code{NULL} (the default), will be computed internally.}

\item{omega}{Optional vector of scale values.
If \code{NULL} (the default), will be computed internally.}

\item{log}{Logical; If \code{TRUE}, values are returned on the log scale.}

\item{lower.tail}{Logical; If \code{TRUE} (default), return P(X <= x).
Else, return P(X > x) .}

\item{log.p}{Logical; If \code{TRUE}, values are returned on the log scale.}

\item{p}{Vector of probabilities.}

\item{tol}{Tolerance of the approximation used in the
computation of quantiles.}

\item{n}{Number of draws to sample from the distribution.}
}
\description{
Density, distribution function, and random generation for the
skew-normal distribution with mean \code{mu},
standard deviation \code{sigma}, and skewness \code{alpha}.
}
\details{
See \code{vignette("brms_families")} for details
on the parameterization.
}
