---
title: "News"
output: html_document
---

# under development
----------------------------------------------------------------

* variable names containing double underscores '__' are now prohibited to avoid naming conflicts
* multiple grouping factors can now be combined within one random effects argument using the interaction symbol ':'
* generalized S3 method 'hypothesis' to be used with all parameter classes not just fixed effects
* introduced new family 'multigaussian' allowing for multivariate normal regression
* function calls with several arguments (e.g. poly(x,3)) are now allowed in the formula argument of function 'brm'
* fixed a bug in S3 method 'hypothesis' leading to an error when numbers with decimal places were used in the formulation of the hypotheses. 

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