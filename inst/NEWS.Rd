\name{NEWS}
\title{brms News}
\encoding{UTF-8}

\section{CHANGES IN VERSION 0.7.0.9000}{
  \subsection{NEW FEATURES}{
    \itemize{  
      \item Restrict domain of fixed effects parameters
      using new arguments \code{lb} and \code{ub}
      in function \code{set_prior} thanks to the idea
      of Joel Gombin.
    }
  }
  \subsection{OTHER CHANGES}{
    \itemize{
      \item Optionally define the intercept
      as an ordinary fixed effect to avoid the
      reparametrization via centering of the
      fixed effects design matrix.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fix problems when predicting with \code{newdata}
      for zero-inflated and hurdle models thanks to
      Ruben Arslan.
      \item Fix problems when predicting with \code{newdata}
      if it is a subset of the data stored in a 
      \code{brmsfit} object thanks to Ruben Arslan.
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