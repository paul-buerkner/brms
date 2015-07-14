<!-- README.md is generated from README.Rmd. Please edit that file -->
brms
====

The <b>brms</b> package provides an interface to fit bayesian generalized linear mixed models using Stan, which is a C++ package for obtaining Bayesian inference using the No-U-turn sampler (see <http://mc-stan.org/>). The formula syntax is very similar to that of the package lme4 to provide a familiar and simple interface for performing regression analyses.

<!--

-->
How to use brms
===============

``` r
library(brms)
```

As a simple example, we use poisson regression to model the seizure counts in epileptic patients to investigate whether the treatment (represented by variable Trt\_c) can reduce the seizure counts. Two random intercepts are incorporated to account for the variance between patients and visits.

``` r
fit <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1|patient) + (1|visit), 
           data = epilepsy, family = "poisson")
```

If rstan is not installed, brm will return the Stan model, the required data, and the parameters of interest, which are the important prerequisites to fit the model in Stan. If rstan is installed, the model is fitted automatically and the results (i.e. posterior samples) can be investigated using

``` r
summary(fit) 
#>  Family: poisson (log) 
#> Formula: count ~ log_Age_c + log_Base4_c + Trt_c + (1 | patient) + (1 | visit) + log_Base4_c:Trt_c 
#>    Data: epilepsy (Number of observations: 236) 
#> Samples: 2 chains, each with n.iter = 2000; n.warmup = 500; n.thin = 1; 
#>          total post-warmup samples = 3000
#>  
#> Random Effects: 
#> ~patient (Number of levels: 59) 
#>               Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
#> sd(Intercept)     0.55      0.07     0.43     0.71        587    1
#> 
#> ~visit (Number of levels: 4) 
#>               Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
#> sd(Intercept)     0.13      0.12     0.02     0.47        504    1
#> 
#> Fixed Effects: 
#>                   Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
#> Intercept             0.81      0.06     0.69     0.92        663    1
#> log_Age_c             0.48      0.38    -0.26     1.25        690    1
#> log_Base4_c           1.06      0.11     0.85     1.28        780    1
#> Trt_c                -0.34      0.16    -0.67    -0.02        582    1
#> log_Base4_c:Trt_c     0.34      0.23    -0.10     0.77        607    1
#> 
#> Samples were drawn using NUTS(diag_e). For each parameter, Eff.Sample is a 
#> crude measure of effective sample size, and Rhat is the potential scale 
#> reduction factor on split chains (at convergence, Rhat = 1).
```

On the top of the output, some general information on the model is given, such as family, formula, as well as number of iterations and chains. Next, random effects are displayed seperately for each grouping factor in terms of standard deviations and (in case of more than one random effect per grouping factor; not displayed here) correlations between random effects. On the bottom of the output, fixed effects are displayed. If incorporated, autocorrelation effects and family specific parameters (e.g., the residual standard deviation 'sigma' in normal models) are also given.

In general, every parameter is summarized using the mean ('Estimate') and the standard deviation ('Est.Error') of the posterior distribution as well as two-sided 95% Credible intervals ('l-95% CI' and 'u-95% CI') based on quantiles. The last two values ('Eff.Sample' and 'Rhat') provide information on how well the algorithm could estimate the posterior distribution of this parameter. If 'Rhat' is considerably greater than 1, the algorithm has not yet converged and it is necessary to run more iterations and / or set stronger priors.

To visually investigate the chains as well as the posterior, you can use

``` r
plot(fit) 
```

For a complete list of methods to apply on <b>brms</b> models see

``` r
methods(class = "brmsfit") 
#>  [1] fixef             formula           hypothesis        ngrps             nobs              par.names        
#>  [7] plot              posterior.samples predict           print             ranef             summary          
#> [13] VarCorr          
#> see '?methods' for accessing help and source code
```

For details on formula syntax, families and link functions, as well as prior distributions see the help page of the <b>brm</b> function:

``` r
help(brm) 
```

How to install brms
===================

``` r
install.packages("brms")
```

Without having rstan installed, the function brm will return the Stan model, the required data, and the parameters of interest. To allow brm to fit the model automatically, the package rstan has to be installed manually, as it is not on CRAN, yet. However, the developers of Stan and rstan are currently working on a version to be uploaded on CRAN. In the meantime, instructions on how to install rstan can be found at <https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started>.

<!-- Before you will be able to actually fit bayesian models with brms, the package rstan has to be installed manually, as it is not on CRAN, yet. First, you need a C++ compiler. See https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started#prerequisites for instructions on how to get one. Second, install rstan by running the following R code (the number behind 'j' in the first line corresponds to the number of cores to use for the installation). This may take a few minutes and you should restart R after the installation.


```r
Sys.setenv(MAKEFLAGS = "-j1") 
source('http://mc-stan.org/rstan/install.R', echo = TRUE, max.deparse.length = 2000)
install_rstan()
```
-->
