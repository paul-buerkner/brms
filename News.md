---
title: "News"
output: html_document
---

# 0.4.0
----------------------------------------------------------------

## new features
* implement the Wakanabe-Akaike Information Criterion (WAIC)
* implement the ||-syntax for random effects allowing for the estimation of random effects standard deviations without the estimation of correlations.
* allow to combine multiple grouping factors within one random effects argument using the interaction symbol ':'
* generalize S3 method 'hypothesis' to be used with all parameter classes not just fixed effects. In addition, one-sided hypothesis testing is now possible.
* introduce new family 'multigaussian' allowing for multivariate normal regression.
* introduce new family 'bernoulli' for dichotomous response variables as a more efficient alternative to families 'binomial' or 'categorical' in this special case.

## minor changes
* slightly change the internal structure of brms to reflect that rstan is finally on CRAN.
* thoroughly check validity of the response variable before the data is passed to Stan.
* prohibit variable names containing double underscores '__' to avoid naming conflicts.
* allow function calls with several arguments (e.g. poly(x,3)) in the formula argument of function 'brm.
* random effects estimates returned by S3 method 'ranef' are now always centered around zero.
* prevent the use of customized covariance matrices for grouping factors with multiple random effects for now. 

## bug fixes
* fix a bug in S3 method 'hypothesis' leading to an error when numbers with decimal places were used in the formulation of the hypotheses. 
* fix a bug in S3 method 'ranef' that caused an error for grouping factors with only one random effect.
* fix a bug that could cause the fixed intercept to be wrongly estimated in the presence of multiple random intercepts.

# brms 0.3.0
----------------------------------------------------------------

* introduced new methods 'par.names' and 'posterior.samples' for class 'brmsfit' to extract parameter names and posterior samples for given parameters, respectively.
* introduced new method 'hypothesis' for class 'brmsfit' allowing to test non-linear hypotheses concerning fixed effects
* introduced new argument 'addition' in function brm to get a more flexible approach in specifying additional information on the response variable (e.g., standard errors for meta-analysis). Alternatively, this information can also be passed to the formula argument directly.
* introduced weighted and censored regressions through argument 'addition' of function brm
* introduced new argument 'cov.ranef' in function brm allowing for customized covariance structures of random effects
* introduced new argument 'autocor' in function brm allowing for autocorrelation of the response variable.
* introduced new functions 'cor.ar', 'cor.ma', and 'cor.arma', to be used with argument 'autocor' for modeling autoregressive, moving-average, and autoregressive-moving-average models. 
* amended parametrization of random effects to increase efficiency of the sampling algorithms
* improved vectorization of sampling statements
* fixed a bug that could cause an error when fitting poisson models while predict = TRUE
* fixed a bug that caused an error when sampling only one chain while silent = TRUE 

# brms 0.2.0
----------------------------------------------------------------

* new S3 class 'brmsfit' to be returned by function brm
* new methods for class 'brmsfit': 
  summary, print, plot, predict, fixef, ranef, and VarCorr, nobs, ngrps, formula
* introduced new argument 'silent' in function brm, allowing to suppress most 
  of stan's intermediate output
* introduced new families 'negbinomial' (negative binomial) and 'geometric' to allow for more flexibility in modeling count data
* fixed a bug that caused an error when formulas contained 
  more complicated function calls
* fixed a bug that caused an error when posterior predictives were sampled for family 'cumulative'
* fixed a bug that prohibited to use improper flat priors for parameters that have proper priors by default
* amended warning and error messages to make them more informative
* corrected examples in the documentation
* extended the README file

# brms 0.1.0 
----------------------------------------------------------------

* inital release version