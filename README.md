<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Build Status](https://travis-ci.org/paul-buerkner/brms.svg?branch=master)](https://travis-ci.org/paul-buerkner/brms) [![Coverage Status](https://codecov.io/github/paul-buerkner/brms/coverage.svg?branch=master)](https://codecov.io/github/paul-buerkner/brms?branch=master) [![CRAN Version](http://www.r-pkg.org/badges/version/brms)](http://cran.r-project.org/package=brms)

brms
====

The <b>brms</b> package provides an interface to fit Bayesian generalized (non-)linear mixed models using Stan, which is a C++ package for obtaining Bayesian inference using the No-U-turn sampler (see <http://mc-stan.org/>). The formula syntax is very similar to that of the package lme4 to provide a familiar and simple interface for performing regression analyses.

<!--

-->
How to use brms
===============

``` r
library(brms)
```

As a simple example, we use poisson regression to model the seizure counts in epileptic patients to investigate whether the treatment (represented by variable Trt\_c) can reduce the seizure counts. Three random intercepts are incorporated to account for the variance between patients and visits, as well as for the residual variance.

``` r
fit <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1|patient) + (1|visit) + (1|obs), 
           data = epilepsy, family = "poisson")
#> Compiling the C++ model
```

The results (i.e. posterior samples) can be investigated using

``` r
summary(fit, waic = TRUE) 
#>  Family: poisson (log) 
#> Formula: count ~ log_Age_c + log_Base4_c * Trt_c + (1 | patient) + (1 | visit) + (1 | obs) 
#>    Data: epilepsy (Number of observations: 236) 
#> Samples: 4 chains, each with iter = 2000; warmup = 1000; thin = 1; 
#>          total post-warmup samples = 4000
#>    WAIC: 1149.34
#>  
#> Random Effects: 
#> ~obs (Number of levels: 236) 
#>               Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
#> sd(Intercept)     0.37      0.04     0.29     0.45       1433    1
#> 
#> ~patient (Number of levels: 59) 
#>               Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
#> sd(Intercept)     0.51      0.07     0.38     0.66       1154    1
#> 
#> ~visit (Number of levels: 4) 
#>               Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
#> sd(Intercept)     0.11      0.11        0     0.41       1129    1
#> 
#> Fixed Effects: 
#>                   Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
#> Intercept             1.56      0.10     1.37     1.76       1790    1
#> log_Age_c             0.48      0.36    -0.24     1.20       1520    1
#> log_Base4_c           1.06      0.11     0.85     1.28       1640    1
#> Trt_c                -0.34      0.16    -0.65    -0.03       1724    1
#> log_Base4_c:Trt_c     0.35      0.22    -0.09     0.77       1600    1
#> 
#> Samples were drawn using sampling(NUTS). For each parameter, Eff.Sample 
#> is a crude measure of effective sample size, and Rhat is the potential 
#> scale reduction factor on split chains (at convergence, Rhat = 1).
```

On the top of the output, some general information on the model is given, such as family, formula, number of iterations and chains, as well as the WAIC, which is an information criterion for Bayesian models. Next, random effects are displayed seperately for each grouping factor in terms of standard deviations and (in case of more than one random effect per grouping factor; not displayed here) correlations between random effects. On the bottom of the output, fixed effects are displayed. If incorporated, autocorrelation effects and family specific parameters (e.g., the residual standard deviation 'sigma' in normal models) are also given.

In general, every parameter is summarized using the mean ('Estimate') and the standard deviation ('Est.Error') of the posterior distribution as well as two-sided 95% Credible intervals ('l-95% CI' and 'u-95% CI') based on quantiles. The last two values ('Eff.Sample' and 'Rhat') provide information on how well the algorithm could estimate the posterior distribution of this parameter. If 'Rhat' is considerably greater than 1, the algorithm has not yet converged and it is necessary to run more iterations and / or set stronger priors.

To visually investigate the chains as well as the posterior distributions, you can use

``` r
plot(fit) 
```

An even more detailed investigation can be achieved by applying the shinystan package:

``` r
launch_shiny(fit) 
```

There are several methods to compute and visualize model predictions. Suppose that we want to predict responses (i.e. seizure counts) of a person in the treatment group (`Trt_c = 0.5`) and in the control group (`Trt_c = -0.5`) with average age and average number of previous seizures at the fourth visit. Than we can use

``` r
newdata <- data.frame(Trt_c = c(0.5, -0.5), log_Age_c = 0, 
                      log_Base4_c = 0, visit = 4)
predict(fit, newdata = newdata, allow_new_levels = TRUE, probs = c(0.05, 0.95))
#>   Estimate Est.Error 5%ile 95%ile
#> 1  4.75325  3.938558     0     12
#> 2  6.66450  5.346362     1     17
```

We need to set `allow_new_levels = TRUE` because we want to predict responses of a person that was not present in the data used to fit the model. While the `predict` method returns predictions of the responses, the `fitted` method returns predictions of the regression line.

``` r
fitted(fit, newdata = newdata, allow_new_levels = TRUE, probs = c(0.05, 0.95))
#>   Estimate Est.Error    5%ile   95%ile
#> 1 4.770850  3.246336 1.363647 11.01045
#> 2 6.678924  4.605785 1.939766 15.40695
```

Both methods return the same etimate (up to random error), while the latter has smaller variance, because the uncertainty in the regression line is smaller than the uncertainty in each response. If we want to predict values of the original data, we can just leave the `newdata` argument empty.

A related feature is the computation and visualization of marginal effects, which can help in better understanding the influence of the predictors on the response.

``` r
plot(marginal_effects(fit, probs = c(0.05, 0.95)))
```

For a complete list of methods to apply on <b>brms</b> models see

``` r
methods(class = "brmsfit") 
#>  [1] as.mcmc           coef              family            fitted           
#>  [5] fixef             formula           hypothesis        launch_shiny     
#>  [9] logLik            LOO               marginal_effects  model.frame      
#> [13] ngrps             nobs              pairs             parnames         
#> [17] plot              posterior_samples predict           print            
#> [21] prior_samples     ranef             residuals         stancode         
#> [25] standata          stanplot          summary           update           
#> [29] VarCorr           vcov              WAIC             
#> see '?methods' for accessing help and source code
```

Details on formula syntax, families and link functions, as well as prior distributions can be found on the help page of the brm function:

``` r
help(brm) 
```

More instructions on how to use <b>brms</b> are given in the package's vignette.

``` r
vignette("brms") 
```

How to install brms
===================

To install the latest release version from CRAN use

``` r
install.packages("brms")
```

The current developmental version can be downloaded from github via

``` r
library(devtools)
install_github("paul-buerkner/brms")
```

Because <b>brms</b> is based on Stan, a C++ compiler is required. The program Rtools (available on <https://cran.r-project.org/bin/windows/Rtools/>) comes with a C++ compiler for Windows. On Mac, you should use Xcode. For further instructions on how to get the compilers running, see the prerequisites section on <https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started>.
