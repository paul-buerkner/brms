\name{NEWS}
\title{brms News}
\encoding{UTF-8}

\section{CHANGES IN VERSION 1.6.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Fit finite mixture models using family
      function \code{mixture}.
      \item Implement different ways to sample new levels
      of grouping factors in \code{predict} and related 
      methods through argument \code{sample_new_levels}.
      Thanks to Tom Wallis and Jonah Gabry for a detailed 
      discussion about this feature.
      \item Add methods \code{loo_predict}, \code{loo_linpred},
      and \code{loo_predictive_interval} for computing
      LOO predictions thanks to Aki Vehtari and Jonah Gabry.
      \item Allow using \code{offset} in formulas
      of non-linear and auxiliary parameters.
      \item Allow sparse matrix multiplication in
      non-linear and distributional models.
      \item Allow using the \code{identity} link for
      all auxiliary parameters.
      \item Introduce method \code{control_params} to
      conveniently extract control parameters of the 
      NUTS sampler.
      \item Improve flexibility of the \code{conditions} 
      argument of \code{marginal_effects}.
      \item Extend method \code{stanplot} to correctly 
      handle some new \code{mcmc_} plots of the 
      \pkg{bayesplot} package.
    }
  }
  \subsection{OTHER CHANGES}{
    \itemize{
      \item Improve the \code{update} method to
      only recompile models when the \pkg{Stan} code
      changes.
      \item Warn about divergent transitions when calling
      \code{summary} or \code{print} on \code{brmsfit} objects.
      \item Warn about unused variables in argument 
      \code{conditions} when calling \code{marginal_effects}.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fix problems with the inclusion of offsets
      occuring for more complicated formulas thanks to 
      Christian Stock.
      \item Fix a bug that led to invalid Stan code when 
      sampling from priors in intercept only models thanks 
      to Tom Wallis.
      \item Correctly check for category specific
      group-level effects in non-ordinal models thanks to
      Wayne Folta.
      \item Fix problems in \code{pp_check} when specifying
      argument \code{newdata} together with arguments 
      \code{x} or \code{group}. 
      \item Rename the last column in the output of 
      \code{hypothesis} to \code{"star"} in order to avoid 
      problems with zero length column names thanks to 
      the GitHub user puterleat.
      \item Add a missing new line statement at the end 
      of the \code{summary} output thanks to Thomas Kluth.
    }
  }
}

\section{CHANGES IN VERSION 1.5.1}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Allow \code{horseshoe} and \code{lasso} 
      priors to be applied on population-level effects 
      of non-linear and auxiliary parameters.
      \item Force recompiling \pkg{Stan} models
      in \code{update.brmsfit} via argument
      \code{recompile}.
    }
  }
  \subsection{OTHER CHANGES}{
    \itemize{
      \item Avoid indexing of matrices in non-linear
      models to slightly improve sampling speed.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fix a severe problem (introduced in version 1.5.0), 
      when predicting \code{Beta} models thanks to Vivian Lam.
      \item Fix problems when summarizing some models
      fitted with older version of \pkg{brms} thanks
      to Vivian Lam.
      \item Fix checks of argument \code{group} in
      method \code{pp_check} thanks to Thomas K.
      \item Get arguments \code{subset} and \code{nsamples}
      working correctly in \code{marginal_smooths}.
    }
  }
}

\section{CHANGES IN VERSION 1.5.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Implement the generalized extreme value
      distribution via family \code{gen_extreme_value}.
      \item Improve flexibility of the \code{horseshoe}
      prior thanks to Juho Piironen.
      \item Introduce auxiliary parameter \code{mu}
      as an alternative to specifying effects within
      the \code{formula} argument in function
      \code{brmsformula}.
      \item Return fitted values of auxiliary parameters 
      via argument \code{auxpar} of method \code{fitted}.
      \item Add vignette \code{"brms_multilevel"}, in which
      the advanced formula syntax of \pkg{brms} is explained
      in detail using several examples.
    }
  }
  \subsection{OTHER CHANGES}{
    \itemize{
      \item Refactor various parts of the package
      to ease implementation of mixture and multivariate
      models in future updates. This should not have
      any user visible effects.
      \item Save the version number of \pkg{rstan} in 
      element \code{version} of \code{brmsfit} objects.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fix a rare error when predicting \code{von_mises}
      models thanks to John Kirwan.
    }
  }
}

\section{CHANGES IN VERSION 1.4.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Fit quantile regression models via family
      \code{asym_laplace} (asymmetric Laplace distribution).
      \item Specify non-linear models in a (hopefully) more
      intuitive way using \code{brmsformula}.
      \item Fix auxiliary parameters to certain values
      through \code{brmsformula}.
      \item Allow \code{family} to be specified in
      \code{brmsformula}.
      \item Introduce family \code{frechet} for modelling
      strictly positive responses.
      \item Allow truncation and censoring at the same time.
      \item Introduce function \code{prior_} allowing
      to specify priors using one-sided formulas or \code{quote}.
      \item Pass priors to \pkg{Stan} directly without
      performing any checks by setting \code{check = FALSE} 
      in \code{set_prior}.
      \item Introduce method \code{nsamples} to extract
      the number of posterior samples.
      \item Export the main formula parsing function 
      \code{parse_bf}.
      \item Add more options to customize two-dimensional surface
      plots created by \code{marginal_effects} or \code{marginal_smooths}.
    }
  }
  \subsection{OTHER CHANGES}{
    \itemize{
      \item Change structure of \code{brmsformula}
      objects to be more reliable and easier to extend.
      \item Make sure that parameter \code{nu} never 
      falls below \code{1} to reduce convergence problems
      when using family \code{student}.
      \item Deprecate argument \code{nonlinear}.
      \item Deprecate family \code{geometric}.
      \item Rename \code{cov_fixed} to \code{cor_fixed}.
      \item Make handling of addition terms more transparent
      by exporting and documenting related functions.
      \item Refactor helper functions of the \code{fitted}
      method to be easier to extend in the future.
      \item Remove many units tests of internal functions
      and add tests of user-facing functions instead.
      \item Import some generics from \pkg{nlme} instead
      of \pkg{lme4} to remove dependency on the latter one.
      \item Do not apply \code{structure} to \code{NULL}
      anymore to get rid of warnings in \R-devel.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fix problems when fitting smoothing terms
      with factors as \code{by} variables thanks to
      Milani Chaloupka.
      \item Fix a bug that could cause some monotonic 
      effects to be ignored in the \pkg{Stan} code thanks 
      to the GitHub user bschneider.
      \item Make sure that the data of models with
      only a single observation are compatible with 
      the generated \pkg{Stan} code.
      \item Handle argument \code{algorithm} 
      correctly in \code{update.brmsfit}.
      \item Fix a bug sometimes causing an error in
      \code{marginal_effects} when using family
      \code{wiener} thanks to Andrew Ellis.
      \item Fix problems in \code{fitted} when applied
      to \code{zero_inflated_beta} models thanks to
      Milani Chaloupka.
      \item Fix minor problems related to the prediction
      of autocorrelated models.
      \item Fix a few minor bugs related to the backwards 
      compatibility of multivariate and related models 
      fitted with \pkg{brms} < 1.0.0.
    }
  }
}

\section{CHANGES IN VERSION 1.3.1}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Introduce the auxiliary parameter \code{disc}
      ('discrimination') to be used in ordinal models.
      By default it is not estimated but fixed to one.
      \item Create \code{marginal_effects} plots of 
      two-way interactions of variables that were 
      not explicitely modeled as interacting.
    }
  }
  \subsection{OTHER CHANGES}{
    \itemize{
      \item Move \pkg{rstan} to 'Imports' and \pkg{Rcpp} 
      to 'Depends' in order to avoid loading \pkg{rstan} 
      into the global environment automatically.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fix a bug leading to unexpected errors
      in some S3 methods when applied to ordinal models.
    }
  }
}

\section{CHANGES IN VERSION 1.3.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Fit error-in-variables models
      using function \code{me} in the model formulae.
      \item Fit multi-membership models using function
      \code{mm} in grouping terms.
      \item Add families \code{exgaussian}
      (exponentially modified Gaussian distribution) 
      and \code{wiener} (Wiener diffusion model distribution)
      specifically suited to handle for response times.
      \item Add the \code{lasso} prior as an alternative
      to the \code{horseshoe} prior for sparse models.
      \item Add the methods \code{log_posterior},
      \code{nuts_params}, \code{rhat}, and \code{neff_ratio}
      for \code{brmsfit} objects to conveniently access
      quantities used to diagnose sampling behavior.
      \item Combine chains in method \code{as.mcmc} using
      argument \code{combine_chains}.
      \item Estimate the auxiliary parameter
      \code{sigma} in models with known standard errors of
      the response by setting argument \code{sigma} to 
      \code{TRUE} in addition function \code{se}.
      \item Allow visualizing two-dimensional smooths
      with the \code{marginal_smooths} method.
    }
  }
  \subsection{OTHER CHANGES}{
    \itemize{
      \item Require argument \code{data} to be explicitely 
      specified in all user facing functions.
      \item Refactor the \code{stanplot} method
      to use \pkg{bayesplot} on the backend.
      \item Use the \pkg{bayesplot} theme as the default
      in all plotting functions.
      \item Add the abbreviations \code{mo} and \code{cs}
      to specify monotonic and category specific effects
      respectively.
      \item Rename generated variables in the data.frames 
      returned by \code{marginal_effects} to avoid potential
      naming conflicts.
      \item Deprecate argument \code{cluster} and use 
      the native \code{cores} argument of \pkg{rstan} instead.
      \item Remove argument \code{cluster_type} as it is
      no longer required to apply forking.
      \item Remove the deprecated \code{partial} argument.
    }
  }
}

\section{CHANGES IN VERSION 1.2.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Add the new family \code{hurdle_lognormal}
      specifically suited for zero-inflated continuous responses.
      \item Introduce the \code{pp_check} method to perform 
      various posterior predictive checks
      using the \pkg{bayesplot} package.
      \item Introduce the \code{marginal_smooths} method to
      better visualize smooth terms.
      \item Allow varying the scale of global shrinkage 
      parameter of the \code{horseshoe} prior.
      \item Add functions \code{prior} and \code{prior_string}
      as aliases of \code{set_prior}, the former 
      allowing to pass arguments without quotes \code{""}
      using non-standard evaluation.
      \item Introduce four new vignettes explaining how to fit
      non-linear models, distributional models, phylogenetic models, 
      and monotonic effects respectively.
      \item Extend the \code{coef} method to better 
      handle category specific group-level effects.
      \item Introduce the \code{prior_summary} method
      for \code{brmsfit} objects to obtain a summary
      of prior distributions applied.
      \item Sample from the prior of the original population-level 
      intercept when \code{sample_prior = TRUE} even in models
      with an internal temporary intercept used to improve
      sampling efficiency.
      \item Introduce methods \code{posterior_predict}, 
      \code{predictive_error} and \code{log_lik} as
      (partial) aliases of \code{predict}, \code{residuals},
      and \code{logLik} respectively.
    }
  }
  \subsection{OTHER CHANGES}{
    \itemize{
      \item Improve computation of Bayes factors
      in the \code{hypothesis} method to be less 
      influenced by MCMC error.
      \item Improve documentation of default priors.
      \item Refactor internal structure of some 
      formula and prior evaluating functions. 
      This should not have any user visible effects.
      \item Use the \pkg{bayesplot} package as the
      new backend of \code{plot.brmsfit}.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{ 
      \item Better mimic \pkg{mgcv} when parsing smooth terms
      to make sure all arguments are correctly handled.
      \item Avoid an error occuring during the prediction 
      of new data when grouping factors with only a single 
      factor level were supplied thanks to Tom Wallis.
      \item Fix \code{marginal_effects} to consistently 
      produce plots for all covariates in non-linear models 
      thanks to David Auty.
      \item Improve the \code{update} method to better recognize
      situations where recompliation of the \pkg{Stan} code
      is necessary thanks to Raphael P.H.
      \item Allow to correctly \code{update} the \code{sample_prior}
      argument to value \code{"only"}.
      \item Fix an unexpected error occuring in many S3 methods
      when the thinning rate is not a divisor of the total
      number of posterior samples thanks to Paul Zerr.
    }
  }
}

\section{CHANGES IN VERSION 1.1.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Estimate monotonic group-level effects.
      \item Estimate category specific group-level effects.
      \item Allow \code{t2} smooth terms based on
      multiple covariates.
      \item Estimate interval censored data via the
      addition argument \code{cens} in the model formula.
      \item Allow to compute \code{residuals} also based 
      on predicted values instead of fitted values.
    }
  }
  \subsection{OTHER CHANGES}{
    \itemize{
      \item Use the prefix \code{bcs} in parameter names
      of category specific effects and the prefix \code{bm} 
      in parameter names of monotonic effects (instead
      of the prefix \code{b}) to simplify their identifaction.
      \item Ensure full compatibility with 
      \pkg{ggplot2} version 2.2.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{ 
      \item Fix a bug that could result in incorrect
      threshold estimates for \code{cumulative} and
      \code{sratio} models thanks to Peter Congdon.
      \item Fix a bug that sometimes kept distributional 
      \code{gamma} models from being compiled 
      thanks to Tim Beechey.
      \item Fix a bug causing an error in \code{predict} 
      and related methods when two-level factors or 
      logical variables were used as covariates in
      non-linear models thanks to Martin Schmettow.
      \item Fix a bug causing an error when passing 
      lists to additional arguments of smoothing 
      functions thanks to Wayne Folta.
      \item Fix a bug causing an error in the 
      \code{prior_samples} method for models with
      multiple group-level terms that refer to the same 
      grouping factor thanks to Marco Tullio Liuzza.
      \item Fix a bug sometimes causing an error when 
      calling \code{marginal_effects} for weighted models.
    }
  }
}

\section{CHANGES IN VERSION 1.0.1}{
  \subsection{MINOR CHANGES}{
    \itemize{
      \item Center design matrices inside the Stan code
      instead of inside \code{make_standata}.
      \item Get rid of several warning messages
      occuring on CRAN.
    }
  }
}

\section{CHANGES IN VERSION 1.0.0}{
  This is one of the largest updates of \code{brms} since its 
  initial release. In addition to many new features, 
  the multivariate \code{'trait'} syntax has been removed 
  from the package as it was confusing for users, required 
  much special case coding, and was hard to maintain.
  See \code{help(brmsformula)} for details of the formula
  syntax applied in \pkg{brms}.
  \subsection{NEW FEATURES}{
    \itemize{ 
      \item Allow estimating correlations between 
      group-level effects defined across multiple formulae 
      (e.g., in non-linear models) by specifying IDs in
      each grouping term via an extended \pkg{lme4} syntax.
      \item Implement distributional regression models
      allowing to fully predict auxiliary parameters 
      of the response distribution. Among many other 
      possibilities, this can be used to model 
      heterogeneity of variances.
      \item Zero-inflated and hurdle models do not use
      multivariate syntax anymore but instead have
      special auxiliary parameters named \code{zi} and 
      \code{hu} defining zero-inflation / hurdle probabilities.
      \item Implement the \code{von_mises} family to model
      circular responses.
      \item Introduce the \code{brmsfamily} function for
      convenient specification of \code{family} objects.
      \item Allow predictions of \code{t2} smoothing
      terms for new data.
      \item Feature vectors as arguments for the addition 
      argument \code{trunc} in order to model varying
      truncation points.
    }
  }
  \subsection{OTHER CHANGES}{
    \itemize{
      \item Remove the \code{cauchy} family
      after several months of deprecation.
      \item Make sure that group-level parameter names
      are unambiguous by adding double underscores 
      thanks to the idea of the GitHub user schmettow.
      \item The \code{predict} method now returns predicted 
      probabilities instead of absolute frequencies of 
      samples for ordinal and categorical models.
      \item Compute the linear predictor in the model
      block of the Stan program instead of in the
      transformed parameters block. This avoids saving
      samples of unnecessary parameters to disk.
      Thanks goes to Rick Arrano for pointing me 
      to this issue.
      \item Colour points in \code{marginal_effects} plots if 
      sensible.
      \item Set the default of the \code{robust} argument 
      to \code{TRUE} in \code{marginal_effects.brmsfit}.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{ 
      \item Fix a bug that could occur when predicting             
      factorial response variables for new data.
      Only affects categorical and ordinal models.
      \item Fix a bug that could lead to duplicated  
      variable names in the Stan code when sampling 
      from priors in non-linear models thanks to Tom Wallis.
      \item Fix problems when trying to pointwise 
      evaluate non-linear formulae in 
      \code{logLik.brmsfit} thanks to Tom Wallis.
      \item Ensure full compatibility of the \code{ranef}
      and \code{coef} methods with non-linear models.
      \item Fix problems that occasionally occured when
      handling \pkg{dplyr} datasets thanks to the
      GitHub user Atan1988.
    }
  }
}

\section{CHANGES IN VERSION 0.10.0}{
  \subsection{NEW FEATURES}{
    \itemize{ 
      \item Add support for generalized additive mixed models
      (GAMMs). Smoothing terms can be specified using
      the \code{s} and \code{t2} functions in the model formula. 
      \item Introduce \code{as.data.frame} and \code{as.matrix}
      methods for \code{brmsfit} objects.
    }
  }
  \subsection{OTHER CHANGES}{
    \itemize{
      \item The \code{gaussian("log")} family no longer implies
      a log-normal distribution, but a normal distribution with 
      log-link to match the behavior of \code{glm}. 
      The log-normal distribution can now be specified via
      family \code{lognormal}.
      \item Update syntax of \pkg{Stan} models to match the 
      recommended syntax of \pkg{Stan} 2.10.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{ 
      \item The \code{ngrps} method should now always
      return the correct result for non-linear models.
      \item Fix problems in \code{marginal_effects} for
      models using the reserved variable \code{intercept}
      thanks to Frederik Aust.
      \item Fix a bug in the \code{print} method of
      \code{brmshypothesis} objects that could lead to 
      duplicated and thus invalid row names.
      \item Residual standard deviation parameters of 
      multivariate models are again correctly displayed 
      in the output of the \code{summary} method.
      \item Fix problems when using variational Bayes
      algorithms with \pkg{brms} while having 
      \pkg{rstan} >= 2.10.0 installed thanks to the 
      GitHub user cwerner87.
    }
  }
}  

\section{CHANGES IN VERSION 0.9.1}{
  \subsection{NEW FEATURES}{
    \itemize{ 
      \item Allow the '/' symbol in group-level terms in 
      the \code{formula} argument to indicate nested 
      grouping structures.
      \item Allow to compute \code{WAIC} and \code{LOO}
      based on the pointwise log-likelihood using argument 
      \code{pointwise} to substantially reduce memory requirements.
    }
  }
  \subsection{OTHER CHANGES}{
    \itemize{ 
      \item Add horizontal lines to the errorbars in 
      \code{marginal_effects} plots for factors.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fix a bug that could lead to a cryptic error
      message when changing some parts of the model 
      \code{formula} using the \code{update} method.
      \item Fix a bug that could lead to an error when
      calling \code{marginal_effects} for predictors
      that were generated with the \code{base::scale} 
      function thanks to Tom Wallis.
      \item Allow interactions of numeric and categorical
      predictors in \code{marginal_effects} to be passed 
      to the \code{effects} argument in any order.
      \item Fix a bug that could lead to incorrect results
      of \code{predict} and related methods when called with
      \code{newdata} in models using the \code{poly} function
      thanks to Brock Ferguson.
      \item Make sure that user-specified factor contrasts 
      are always applied in multivariate models.
    }
  }
}

\section{CHANGES IN VERSION 0.9.0}{
  \subsection{NEW FEATURES}{
    \itemize{ 
      \item Add support for \code{monotonic} effects
      allowing to use ordinal predictors without 
      assuming their categories to be equidistant.
      \item Apply multivariate formula syntax in categorical 
      models to considerably increase modeling flexibility. 
      \item Add the addition argument \code{disp} to define 
      multiplicative factors on dispersion parameters. 
      For linear models, \code{disp} applies to the residual
      standard deviation \code{sigma} so that it can be
      used to weight observations.
      \item Treat the fixed effects design matrix as sparse
      by using the \code{sparse} argument of \code{brm}. 
      This can considerably reduce working memory 
      requirements if the predictors contain many zeros.
      \item Add the \code{cor_fixed} correlation structure to
      allow for fixed user-defined covariance matrices of the 
      response variable.
      \item Allow to pass self-defined \pkg{Stan} functions
      via argument \code{stan_funs} of \code{brm}. 
      \item Add the \code{expose_functions} method allowing to
      expose self-defined \pkg{Stan} functions in \code{R}.
      \item Extend the functionality of the \code{update} 
      method to allow all model parts to be updated.
      \item Center the fixed effects design matrix also 
      in multivariate models. This may lead to increased
      sampling speed in models with many predictors.
    }
  }
  \subsection{OTHER CHANGES}{
    \itemize{
      \item Refactor \pkg{Stan} code and data generating 
      functions to be more consistent and easier to extent.
      \item Improve checks of user-define prior specifications.
      \item Warn about models that have not converged.
      \item Make sure that regression curves computed by
      the \code{marginal_effects} method are always smooth.
      \item Allow to define category specific effects in
      ordinal models directly within the \code{formula}
      argument.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fix problems in the generated \pkg{Stan} code
      when using very long non-linear model formulas
      thanks to Emmanuel Charpentier.
      \item Fix a bug that prohibited to change priors
      on single standard deviation parameters
      in non-linear models thanks to Emmanuel Charpentier.
      \item Fix a bug that prohibited to use nested
      grouping factors in non-linear models thanks to
      Tom Wallis.
      \item Fix a bug in the linear predictor computation
      within \code{R}, occuring for ordinal models 
      with multiple category specific effects. This
      could lead to incorrect outputs of \code{predict},
      \code{fitted}, and \code{logLik} for these models.
      \item Make sure that the global \code{"contrasts"} option 
      is not used when post-processing a model.
    }
  }
}

\section{CHANGES IN VERSION 0.8.0}{
  \subsection{NEW FEATURES}{
    \itemize{ 
      \item Implement generalized non-linear models, which
      can be specified with the help of the \code{nonlinear} 
      argument in \code{brm}.
      \item Compute and plot marginal effects using the
      \code{marginal_effects} method thanks to the help 
      of Ruben Arslan.
      \item Implement zero-inflated beta models through
      family \code{zero_inflated_beta} thanks to the 
      idea of Ali Roshan Ghias.
      \item Allow to restrict domain of fixed effects and
      autocorrelation parameters using new arguments 
      \code{lb} and \code{ub} in function \code{set_prior}
      thanks to the idea of Joel Gombin.
      \item Add an \code{as.mcmc} method for compatibility 
      with the \pkg{coda} package.
      \item Allow to call the \code{WAIC}, \code{LOO}, 
      and \code{logLik} methods with new data.
    }
  }
  \subsection{OTHER CHANGES}{
    \itemize{
      \item Make sure that \pkg{brms} is fully compatible 
      with \pkg{loo} version 0.1.5.
      \item Optionally define the intercept
      as an ordinary fixed effect to avoid the
      reparametrization via centering of the
      fixed effects design matrix.
      \item Do not compute the WAIC in \code{summary}
      by default anymore to reduce computation time
      of the method for larger models.
      \item The \code{cauchy} family is now deprecated 
      and will be removed soon as it often has convergence 
      issues and not much practical application anyway.
      \item Change the default settings of the number of 
      chains and warmup samples to the defaults of \pkg{rstan}
      (i.e., \code{chains = 4} and \code{warmup = iter / 2}).
      \item Do not remove bad behaving chains anymore as
      they may point to general convergence problems that 
      are dangerous to ignore.
      \item Improve flexibility of the \code{theme} argument
      in all plotting functions.
      \item Only show the legend once per page, when computing
      trace and density plots with the \code{plot} method.
      \item Move code of self-defined \pkg{Stan} functions
      to \code{inst/chunks} and incorporate them into the
      models using \code{rstan::stanc_builder}.
      Also, add unit tests for these functions.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fix problems when predicting with \code{newdata}
      for zero-inflated and hurdle models thanks to Ruben Arslan.
      \item Fix problems when predicting with \code{newdata}
      if it is a subset of the data stored in a 
      \code{brmsfit} object thanks to Ruben Arslan.
      \item Fix data preparation for multivariate models
      if some responses are \code{NA} thanks to Raphael Royaute.
      \item Fix a bug in the \code{predict} method occurring 
      for some multivariate models so that it now always 
      returns the predictions of all response variables, 
      not just the first one.
      \item Fix a bug in the log-likelihood computation of
      \code{hurdle_poisson} and \code{hurdle_negbinomial} models.
      This may lead to minor changes in the values obtained by
      \code{WAIC} and \code{LOO} for these models.
      \item Fix some backwards compatibility issues of models fitted 
      with version <= 0.5.0 thanks to Ulf Koether.
    }
  }
}

\section{CHANGES IN VERSION 0.7.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Use variational inference algorithms
      as alternative to the NUTS sampler by specifying
      argument \code{algorithm} in the \code{brm} function.
      \item Implement beta regression models through family
      \code{Beta}.
      \item Implement zero-inflated binomial models through family
      \code{zero_inflated_binomial}.
      \item Implement multiplicative effects for family
      \code{bernoulli} to fit (among others) 2PL IRT models.
      \item Generalize the \code{formula} argument for zero-inflated
      and hurdle models so that predictors can be included in only 
      one of the two model parts thanks to the idea of Wade Blanchard.
      \item Combine fixed and random effects estimates using
      the new \code{coef} method.
      \item Call the \code{residuals} method with \code{newdata} thanks
      to the idea of Friederike Holz-Ebeling.
      \item Allow new levels of random effects grouping
      factors in the \code{predict}, \code{fitted}, and \code{residuals}
      methods using argument \code{allow_new_levels}.
      \item Selectively exclude random effects in the \code{predict}, 
      \code{fitted}, and \code{residuals} methods using 
      argument \code{re_formula}.
      \item Add a \code{plot} method for objects returned by
      method \code{hypothesis} to visualize prior and posterior
      distributions of the hypotheses being tested.
    }
  }
  \subsection{OTHER CHANGES}{
    \itemize{
      \item Improve evaluation of the response
      part of the \code{formula} argument to
      reliably allow terms with more than one variable
      (e.g., \code{y/x ~ 1}).
      \item Improve sampling efficiency of models containing
      many fixed effects through centering the fixed effects
      design matrix thanks to Wayne Folta.
      \item Improve sampling efficiency of models containing
      uncorrelated random effects specified by means
      of \code{(random || group)} terms in \code{formula}
      thanks to Ali Roshan Ghias.
      \item Utilize user-defined functions in the \pkg{Stan} code 
      of ordinal models to improve readability as well as
      sampling efficiency.
      \item Make sure that model comparisons using \code{LOO} 
      or \code{WAIC} are only performed when models are 
      based on the same responses.
      \item Use some generic functions of the \pkg{lme4}
      package to avoid unnecessary function masking. This 
      leads to a change in the argument order of
      method \code{VarCorr}.
      \item Change the \code{ggplot} theme in the
      \code{plot} method through argument \code{theme}.
      \item Remove the \code{n.} prefix in arguments
      \code{n.iter}, \code{n.warmup}, \code{n.thin}, \code{n.chains}, 
      and \code{n.cluster} of the \code{brm} function.
      The old argument names remain usable as deprecated aliases.
      \item Amend names of random effects parameters to simplify
      matching with their respective grouping factor levels.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fix a bug in the \code{hypothesis} method 
      that could cause valid model parameters to be falsely 
      reported as invalid.
      \item Fix a bug in the \code{prior_samples} method 
      that could cause prior samples of parameters 
      of the same class to be artifically correlated.
      \item Fix \pkg{Stan} code of linear models with 
      moving-average effects and non-identity link functions 
      so that they no longer contain code related solely 
      to autoregressive effects.
      \item Fix a bug in the evaluation of \code{formula} that
      could cause complicated random effects terms to be 
      falsely treated as fixed effects.
      \item Fix several bugs when calling the \code{fitted} 
      and \code{predict} methods with \code{newdata} thanks
      to Ali Roshan Ghias.
    }
  }
}

\section{CHANGES IN VERSION 0.6.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Add support for zero-inflated and hurdle models thanks
      to the idea of Scott Baldwin.
      \item Implement inverse gaussian models through family 
      \code{inverse.gaussian}.
      \item Allow to specify truncation boundaries of the response 
      variable thanks to the idea of Maciej Beresewicz.
      \item Add support for autoregressive (AR) effects of residuals, which can
      be modeled using the \code{cor_ar} and \code{cor_arma} functions.
      \item Stationary autoregressive-moving-average (ARMA) effects 
      of order one can now also be fitted using special covariance matrices.
      \item Implement multivariate student-t models.
      \item Binomial and ordinal families now support 
      the \code{cauchit} link function.
      \item Allow family functions to be used in the 
      \code{family} argument.
      \item Easy access to various \code{rstan} plotting functions 
      using the \code{stanplot} method.
      \item Implement horseshoe priors to model sparsity in 
      fixed effects coefficients thanks to the idea of Josh Chang.
      \item Automatically scale default standard deviation priors so that
      they remain only weakly informative independent on the response scale.
      \item Report model weights computed by the \pkg{loo} package 
      when comparing multiple fitted models.
    }
  }
  \subsection{OTHER CHANGES}{
    \itemize{
      \item Separate the fixed effects Intercept from other fixed effects
      in the \pkg{Stan} code to slightly improve sampling efficiency.
      \item Move autoregressive (AR) effects of the response from 
      the \code{cor_ar} to the \code{cor_arr} function as the result 
      of implementing AR effects of residuals.
      \item Improve checks on argument \code{newdata} used in the
      \code{fitted} and \code{predict} method.
      \item Method \code{standata} is now the only way to extract data
      that was passed to \pkg{Stan} from a \code{brmsfit} object.
      \item Slightly improve \pkg{Stan} code for models 
      containing no random effects.
      \item Change the default prior of the degrees of freedom of the 
      \code{student} family to \code{gamma(2,0.1)}.
      \item Improve readability of the output of method \code{VarCorr}.
      \item Export the \code{make_stancode} function to give users 
      direct access to \pkg{Stan} code generated by \pkg{brms}.
      \item Rename the \code{brmdata} function to \code{make_standata}.
      The former remains usable as a deprecated alias.
      \item Improve documenation to better explain differences in 
      autoregressive effects across R packages.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fix a bug that could cause an unexpected error
      when the \code{predict} method was called with \code{newdata}.
      \item Avoid side effects of the \code{rstan} compilation 
      routines that could occasionally cause R to crash.
      \item Make \pkg{brms} work correctly with \pkg{loo} version 0.1.3
      thanks to Mauricio Garnier Villarreal and Jonah Gabry.
      \item Fix a bug that could cause WAIC and LOO estimates to be slightly 
      incorrect for \code{gaussian} models with \code{log} link.
    }
  }
}

\section{CHANGES IN VERSION 0.5.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Compute the Watanabe-Akaike information criterion (WAIC) and 
      leave-one-out cross-validation (LOO) using the \pkg{loo} package.
      \item Provide an interface to \pkg{shinystan} with 
      S3 method \code{launch_shiny}.
      \item New functions \code{get_prior} and \code{set_prior} to make 
      prior specifications easier.
      \item Log-likelihood values and posterior predictive samples can now 
      be calculated within R after the model has been fitted. 
      \item Make predictions based on new data using S3 method \code{predict}.
      \item Allow for customized covariance structures of 
      grouping factors with multiple random effects.
      \item New S3 methods \code{fitted} and \code{residuals} 
      to compute fitted values and residuals, respectively.
    }
  }
  \subsection{OTHER CHANGES}{
    \itemize{
      \item Arguments \code{WAIC} and \code{predict} are removed 
      from the \code{brm} function, as they are no longer necessary.
      \item New argument \code{cluster_type} in function \code{brm} 
      allowing to choose the cluster type created by the parallel package.
      \item Remove chains that fail to initialize while sampling 
      in parallel leaving the other chains untouched.
      \item Redesign trace and density plots to be faster and more stable.
      \item S3 method \code{VarCorr} now always returns covariance 
      matrices regardless of whether correlations were estimated.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
       \item Fix a bug in S3 method \code{hypothesis} related 
       to the calculation of Bayes-factors for point hypotheses.
       \item User-defined covariance matrices that are not 
       strictly positive definite for numerical reasons should 
       now be handled correctly.
       \item Fix problems when a factor is used as fixed effect
       and as random effects grouping variable at the same time
       thanks to Ulf Koether.
       \item Fix minor issues with internal parameter naming.
       \item Perform additional checking on user defined priors.
    }
  }
}

\section{CHANGES IN VERSION 0.4.1}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Allow for sampling from all specified proper 
      priors in the model.
      \item Compute Bayes-factors for point hypotheses 
      in S3 method \code{hypothesis}.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fix a bug that could cause an error for models with 
      multiple grouping factors thanks to Jonathan Williams.
      \item Fix a bug that could cause an error for
      weighted poisson and exponential models. 
    }
  }
}

\section{CHANGES IN VERSION 0.4.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Implement the Watanabe-Akaike Information Criterion (WAIC).
      \item Implement the \code{||}-syntax for random effects allowing 
      for the estimation of random effects standard deviations 
      without the estimation of correlations.
      \item Allow to combine multiple grouping factors within one random 
      effects argument using the interaction symbol \code{:}.
      \item Generalize S3 method \code{hypothesis} to be used with 
      all parameter classes not just fixed effects. 
      In addition, one-sided hypothesis testing is now possible.
      \item Introduce new family \code{multigaussian} allowing 
      for multivariate normal regression.
      \item Introduce new family \code{bernoulli} for dichotomous 
      response variables as a more efficient alternative to families
      \code{binomial} or \code{categorical} in this special case.
    }
  }
  \subsection{OTHER CHANGES}{
    \itemize{
      \item Slightly change the internal structure of brms to reflect that 
      \pkg{rstan} is finally on CRAN.
      \item Thoroughly check validity of the response variable 
      before the data is passed to \pkg{Stan}.
      \item Prohibit variable names containing double underscores 
      \code{__} to avoid naming conflicts.
      \item Allow function calls with several arguments 
      (e.g. \code{poly(x,3)}) in the formula argument of 
      function \code{brm}.
      \item Always center random effects estimates returned 
      by S3 method \code{ranef} around zero.
      \item Prevent the use of customized covariance matrices for 
      grouping factors with multiple random effects for now. 
      \item Remove any experimental \pkg{JAGS} code from the package. 
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fix a bug in S3 method \code{hypothesis} 
      leading to an error when numbers with decimal places 
      were used in the formulation of the hypotheses. 
      \item Fix a bug in S3 method \code{ranef} that caused 
      an error for grouping factors with only one random effect.
      \item Fix a bug that could cause the fixed intercept to be 
      wrongly estimated in the presence of multiple random intercepts
      thanks to Jarrod Hadfield.
    }
  }
}

\section{CHANGES IN VERSION 0.3.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Introduce new methods \code{parnames} and 
      \code{posterior_samples} for class 'brmsfit' to extract 
      parameter names and posterior samples for given parameters, 
      respectively.
      \item Introduce new method \code{hypothesis} for class \code{brmsfit} 
      allowing to test non-linear hypotheses concerning fixed effects.
      \item Introduce new argument \code{addition} in function brm 
      to get a more flexible approach in specifying additional information 
      on the response variable (e.g., standard errors for meta-analysis). 
      Alternatively, this information can also be passed to the 
      \code{formula} argument directly.
      \item Introduce weighted and censored regressions through argument 
      \code{addition} of function brm.
      \item Introduce new argument \code{cov.ranef} in the \code{brm}
      function allowing for customized covariance structures of 
      random effects thanks to the idea of Boby Mathew.
      \item Introduce new argument \code{autocor} in function brm 
      allowing for autocorrelation of the response variable.
      \item Introduce new functions \code{cor.ar}, \code{cor.ma}, 
      and \code{cor.arma}, to be used with argument \code{autocor} 
      for modeling autoregressive, moving-average, 
      and autoregressive-moving-average models. 
    }
  }
  \subsection{OTHER CHANGES}{
    \itemize{
      \item Amend parametrization of random effects to 
      increase efficiency of the sampling algorithms.
      \item Improve vectorization of sampling statements.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fix a bug that could cause an error when fitting 
      poisson models while \code{predict = TRUE}.
      \item Fix a bug that caused an error when sampling only one chain 
      while \code{silent = TRUE}.  
    }
  }
}

\section{CHANGES IN VERSION 0.2.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item New S3 class \code{brmsfit} to be returned by 
      the \code{brm} function.
      \item New methods for class \code{brmsfit}: 
      \code{summary}, \code{print}, \code{plot}, \code{predict}, 
      \code{fixef}, \code{ranef}, \code{VarCorr}, 
      \code{nobs}, \code{ngrps}, and \code{formula}.
      \item Introduce new argument \code{silent} in the 
      \code{brm} function, allowing to suppress most of 
      \pkg{Stan}'s intermediate output.
      \item Introduce new families \code{negbinomial} (negative binomial) 
      and \code{geometric} to allow for more flexibility 
      in modeling count data.
    }
  }
  \subsection{OTHER CHANGES}{
    \itemize{
      \item Amend warning and error messages to make 
      them more informative.
      \item Correct examples in the documentation.
      \item Extend the README file. 
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fix a bug that caused problems when formulas contained 
      more complicated function calls.
      \item Fix a bug that caused an error when posterior predictives were 
      sampled for family \code{cumulative}
      \item Fix a bug that prohibited to use of improper flat priors 
      for parameters that have proper priors by default.
    }
  }
}

\section{CHANGES IN VERSION 0.1.0}{
  \itemize{
    \item Initial release version
  }
}