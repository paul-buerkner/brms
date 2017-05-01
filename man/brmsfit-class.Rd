\name{brmsfit-class}
\alias{brmsfit-class}
\alias{brmsfit}
\docType{class}
\title{
Class \code{brmsfit} of models fitted with the \pkg{brms} package
}
\description{
Models fitted with the \code{\link[brms:brms]{brms}} package are represented as a \code{brmsfit} object, which contains the posterior samples, model formula, \pkg{Stan} code, relevant data, and other useful information.
}

\details{
  Currently, the following S3 methods are available for class \code{brmsfit}: 
  
  \code{\link[brms:posterior_samples]{as.data.frame}},
  \code{\link[brms:posterior_samples]{as.matrix}},
  \code{\link[brms:as.mcmc.brmsfit]{as.mcmc}},
  \code{\link[brms:coef.brmsfit]{coef}},
  \code{\link[brms:expose_functions]{expose_functions}},
  \code{\link[stats:family]{family}},
  \code{\link[brms:fitted.brmsfit]{fitted}},
  \code{\link[brms:fixef.brmsfit]{fixef}},
  \code{\link[stats:formula]{formula}},
  \code{\link[brms:hypothesis]{hypothesis}},
  \code{\link[brms:launch_shiny.brmsfit]{launch_shiny}},
  \code{\link[brms:logLik.brmsfit]{logLik}},
  \code{\link[brms:LOO.brmsfit]{LOO}},
  \code{\link[brms:marginal_effects.brmsfit]{marginal_effects}},
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
  \code{\link[brms:WAIC.brmsfit]{WAIC}}.
  
  See also \code{methods(class = "brmsfit")}.
}

\section{Slots}{
  \describe{
    \item{\code{formula}:}{model formula; an object of class \code{brmsformula}} 
    \item{\code{family}:}{model family; an object of class \code{brmsfamily}}
    \item{\code{data}:}{model.frame containing all variables used in the model}
    \item{\code{data.name}:}{name of the data frame; a character string} 
    \item{\code{model}:}{model in Stan language; a character string}
    \item{\code{prior}:}{priors applied in the model; a tidy data.frame}
    \item{\code{autocor}:}{an object of class \code{cor_brms} 
    containing the autocorrelation structure}
    \item{\code{threshold}:}{character string defining the threshold
    type used in ordinal models}
        \item{\code{ranef}:}{structure of group-level terms; a tidy data.frame}
    \item{\code{cov_ranef}:}{a list of customized group-level 
    covariance matrices}
    \item{\code{loo}:}{empty slot for adding the \code{loo} information 
    criterion after model fitting}
    \item{\code{waic}:}{empty slot for adding the \code{waic} information 
    criterion after model fitting}
    \item{\code{fit}:}{fitted results including the posterior samples; 
    an object of class \code{stanfit}}
    \item{\code{exclude}:}{names of parameters for which samples 
    are not saved; a character vector}
    \item{\code{algorithm}:}{the name of the algorithm used to fit the model}
     \item{\code{version}:}{the version of \pkg{brms} and \pkg{rstan} with 
     which the model was fitted}
  }
}

\seealso{
\code{\link[brms:brms]{brms}}, 
\code{\link[brms:brm]{brm}},
\code{\link[brms:brmsformula]{brmsformula}}, 
\code{\link[brms:brmsfamily]{brmsfamily}}
}