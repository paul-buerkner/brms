---
title: "News"
output: html_document
---

# under development
----------------------------------------------------------------

## new features
* compute the Watanabe-Akaike information criterion (WAIC) and the leave-one-out cross-validation (LOO) using the loo package.
* provide an interface to shinystan with S3 method 'launch_shiny'.
* log-likelihood values and posterior predictive samples can now be calculated after the model has been fitted. 
* make predictions based on new data using S3 method 'predict'. Currently only available for fixed effects models. 
* allow for customized covariance structures of grouping factors with multiple random effects
* new S3 method 'parnames' for the formula class returning all names of parameters for which priors can be specified.
* new S3 methods 'fitted' and 'residuals' to compute fitted values and residuals, respectively.

## other changes
* arguments WAIC and predict are removed from function brm as they are no longer necessary
* remove chains that fail to initialize while sampling in parallel leaving the other chains untouched
* redesign trace and density plots to be faster and more stable
* S3 method 'VarCorr' now always returns covariance matrices regardless of whether correlations were estimated
* perform additional checking on user specified priors

## bug fixes
* fix a bug in S3 method 'hypothesis' related to the calculation of Bayes factors for point hypotheses
* user defined covariance matrices that are not strictly positive definite for numerical reasons should now be handled correctly.
* fix minor issues with internal parameter naming

# brms 0.4.1
----------------------------------------------------------------

* allow for sampling from all specified proper priors in the model
* calculate Bayes factors for point hypotheses in S3 method 'hypothesis'
* fix a bug that could cause an error for models with multiple grouping factors
* fix a bug that could cause an error for weighted poisson and exponential models

# brms 0.4.0
----------------------------------------------------------------

## new features
* implement the Watanabe-Akaike Information Criterion (WAIC)
* implement the ||-syntax for random effects allowing for the estimation of random effects standard deviations without the estimation of correlations.
* allow to combine multiple grouping factors within one random effects argument using the interaction symbol ':'
* generalize S3 method 'hypothesis' to be used with all parameter classes not just fixed effects. In addition, one-sided hypothesis testing is now possible.
* introduce new family 'multigaussian' allowing for multivariate normal regression.
* introduce new family 'bernoulli' for dichotomous response variables as a more efficient alternative to families 'binomial' or 'categorical' in this special case.

## other changes
* slightly change the internal structure of brms to reflect that rstan is finally on CRAN.
* thoroughly check validity of the response variable before the data is passed to Stan.
* prohibit variable names containing double underscores '__' to avoid naming conflicts.
* allow function calls with several arguments (e.g. poly(x,3)) in the formula argument of function 'brm.
* always center random effects estimates returned by S3 method 'ranef' around zero.
* prevent the use of customized covariance matrices for grouping factors with multiple random effects for now. 
* remove any experimental JAGS code from the package. 

## bug fixes
* fix a bug in S3 method 'hypothesis' leading to an error when numbers with decimal places were used in the formulation of the hypotheses. 
* fix a bug in S3 method 'ranef' that caused an error for grouping factors with only one random effect.
* fix a bug that could cause the fixed intercept to be wrongly estimated in the presence of multiple random intercepts.

# brms 0.3.0
----------------------------------------------------------------

* introduce new methods 'parnames' and 'posterior_samples' for class 'brmsfit' to extract parameter names and posterior samples for given parameters, respectively.
* introduce new method 'hypothesis' for class 'brmsfit' allowing to test non-linear hypotheses concerning fixed effects
* introduce new argument 'addition' in function brm to get a more flexible approach in specifying additional information on the response variable (e.g., standard errors for meta-analysis). Alternatively, this information can also be passed to the formula argument directly.
* introduce weighted and censored regressions through argument 'addition' of function brm
* introduce new argument 'cov.ranef' in function brm allowing for customized covariance structures of random effects
* introduce new argument 'autocor' in function brm allowing for autocorrelation of the response variable.
* introduce new functions 'cor.ar', 'cor.ma', and 'cor.arma', to be used with argument 'autocor' for modeling autoregressive, moving-average, and autoregressive-moving-average models. 
* amend parametrization of random effects to increase efficiency of the sampling algorithms
* improve vectorization of sampling statements
* fix a bug that could cause an error when fitting poisson models while predict = TRUE
* fix a bug that caused an error when sampling only one chain while silent = TRUE 

# brms 0.2.0
----------------------------------------------------------------

* new S3 class 'brmsfit' to be returned by function brm
* new methods for class 'brmsfit': 
  summary, print, plot, predict, fixef, ranef, and VarCorr, nobs, ngrps, formula
* introduce new argument 'silent' in function brm, allowing to suppress most 
  of stan's intermediate output
* introduce new families 'negbinomial' (negative binomial) and 'geometric' to allow for more flexibility in modeling count data
* fix a bug that caused an error when formulas contained more complicated function calls
* fix a bug that caused an error when posterior predictives were sampled for family 'cumulative'
* fix a bug that prohibited to use of improper flat priors for parameters that have proper priors by default
* amend warning and error messages to make them more informative
* correct examples in the documentation
* extend the README file

# brms 0.1.0 
----------------------------------------------------------------

* inital release version