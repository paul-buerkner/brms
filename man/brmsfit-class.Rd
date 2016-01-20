\name{brmsfit-class}
\alias{brmsfit-class}
\alias{brmsfit}
\docType{class}
\title{
Class \code{brmsfit} of fitted mixed-effect models
}
\description{
Models fitted with the \code{\link[brms:brms]{brms}} package are represented as a \code{brmsfit} object, which contains the posterior samples, the model code, the relevant data, the model formula, and some other useful information.
}

\details{
  Currently, the following S3 methods are available for class \code{brmsfit}: 
  
  \code{\link[brms:coef.brmsfit]{coef}},
  \code{\link[stats:family]{family}},
  \code{\link[brms:fitted.brmsfit]{fitted}},
  \code{\link[brms:fixef.brmsfit]{fixef}},
  \code{\link[stats:formula]{formula}},
  \code{\link[brms:hypothesis]{hypothesis}},
  \code{\link[brms:launch_shiny.brmsfit]{launch_shiny}},
  \code{\link[brms:logLik.brmsfit]{logLik}},
  \code{\link[brms:LOO.brmsfit]{LOO}},
  \code{\link[stats:model.frame]{model.frame}},
  \code{\link[brms:ngrps.brmsfit]{ngrps}},
  \code{\link[stats:nobs]{nobs}},
  \code{\link[brms:pairs.brmsfit]{pairs}},
  \code{\link[brms:parnames]{parnames}},
  \code{\link[brms:plot.brmsfit]{plot}},
  \code{\link[brms:posterior_samples]{posterior_samples}}, 
  \code{\link[brms:predict.brmsfit]{predict}}, 
  \code{\link[brms:print.brmsfit]{print}},
  \code{\link[brms:prior_samples]{prior_samples}},
  \code{\link[brms:ranef.brmsfit]{ranef}}, 
  \code{\link[brms:residuals.brmsfit]{residuals}}, 
  \code{\link[brms:stancode.brmsfit]{stancode}},
  \code{\link[brms:standata.brmsfit]{standata}},
  \code{\link[brms:stanplot.brmsfit]{stanplot}},
  \code{\link[brms:summary.brmsfit]{summary}}, 
  \code{\link[brms:update.brmsfit]{update}},   
  \code{\link[brms:VarCorr.brmsfit]{VarCorr}}, 
  \code{\link[brms:vcov.brmsfit]{vcov}}, 
  \code{\link[brms:WAIC.brmsfit]{WAIC}}
}

\section{Slots}{
  \describe{
    \item{\code{formula}:}{model formula; an object of class \code{formula}} 
    \item{\code{family}:}{model family; an object of class \code{family}}
    \item{\code{link}:}{link function; a character string (deprecated)}
    \item{\code{data.name}:}{name of the data frame; a character string} 
    \item{\code{data}:}{model.frame containing all variables used in the model}
    \item{\code{model}:}{model in Stan language; a character string}
    \item{\code{exclude}:}{parameters for which samples are not saved; 
    a character vector}
    \item{\code{prior}:}{priors applied in the model; a matrix}
    \item{\code{ranef}:}{random effects structure; a named list}
    \item{\code{autocor}:}{an object of class \code{cor_brms} 
    containing the autocorrelation structure}
    \item{\code{partial}:}{formula of the category specific effects 
    applied in ordinal models}
    \item{\code{cov_ranef}:}{a list of customized random effects 
    covariance matrices}
    \item{\code{fit}:}{fitted results including the posterior samples; 
    an object of class \code{stanfit}}
    \item{\code{algorithm}:}{the name of the algorithm used to fit the model}
  }
}

\seealso{
\code{\link[brms:brms]{brms}}, \code{\link[brms:brm]{brm}}   
}