\name{brmsfit-class}
\alias{brmsfit-class}
\alias{brmsfit}
\docType{class}
\title{
Class \code{brmsfit} of fitted mixed-effect models
}
\description{
Models fitted with the \code{\link[brms:brms]{brms}} package are represented as a \code{brmsfit} object, which contains, among others, the posterior samples (if package \code{rstan} is installed), the model code, the relevant data, the model formula, etc.
}

\details{
  Currently, the following S3 methods are available for class \code{brmsfit}: 
  
  \code{\link[brms:summary.brmsfit]{summary}}, \code{\link[brms:print.brmsfit]{print}}, 
  \code{\link[brms:plot.brmsfit]{plot}}, \code{\link[brms:fixef.brmsfit]{fixef}}
  \code{\link[brms:ranef.brmsfit]{ranef}},   \code{\link[brms:VarCorr.brmsfit]{VarCorr}}, 
  \code{\link[stats:predict]{predict}}, \code{\link[stats:formula]{formula}},
  \code{\link[stats:nobs]{nobs}}, \code{\link[brms:ngrps.brmsfit]{ngrps}},
  \code{\link[brms:residuals.brmsfit]{residuals}}, \code{\link[brms:loglik.brmsfit]{loglik}},
  \code{\link[brms:hypothesis]{hypothesis}}, \code{\link[brms:parnames]{parnames}},
  \code{\link[brms:posterior_samples]{posterior_samples}}, \code{\link[brms:prior_samples]{prior_samples}},
  \code{\link[brms:linear_predictor.brmsfit]{linear_predictor}}, \code{\link[brms:launch_shiny.brmsfit]{launch_shiny}},
  \code{\link[brms:WAIC.brmsfit]{WAIC}}, \code{\link[brms:LOO.brmsfit]{LOO}},
}

\section{Slots}{
  \describe{
    \item{\code{formula}:}{the model formula; an object of class \code{formula}} 
    \item{\code{family}:}{the model family; a character string}
    \item{\code{link}:}{the link function; a character string}
    \item{\code{data.name}:}{the name of the data frame; a character string} 
    \item{\code{data}:}{the data passed to Stan; a list of vectors and matrices}
    \item{\code{model}:}{the model in Stan language; a character string}
    \item{\code{exclude}:}{parameters for which samples are not saved; a character vector}
    \item{\code{autcor}:}{an object of class 'cor.brms' containing the autocorrelation structure}
    \item{\code{partial}:}{the formula of the partial effects applied in categorical and ordinal models}
    \item{\code{fit}:}{the fitted results including the posterior samples; (usually) an object of class \code{stanfit}}
  }
}

\seealso{
\code{\link[brms:brms]{brms}}, \code{\link[brms:brm]{brm}}   
}