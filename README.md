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
#> Formula: count ~ log_Age_c + log_Base4_c * Trt_c + (1 | patient) + (1 | visit) 
#>    Data: epilepsy (Number of observations: 236) 
#> Samples: 2 chains, each with n.iter = 2000; n.warmup = 500; n.thin = 1; 
#>          total post-warmup samples = 3000
#>  
#> Random Effects: 
#> ~patient (Number of levels: 59) 
#>               Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
#> sd(Intercept)     0.55      0.07     0.43      0.7       1868    1
#> 
#> ~visit (Number of levels: 4) 
#>               Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
#> sd(Intercept)     0.14      0.14     0.02     0.49        959    1
#> 
#> Fixed Effects: 
#>                   Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
#> Intercept             0.81      0.06     0.70     0.92       1280 1.00
#> log_Age_c             0.47      0.38    -0.27     1.21        488 1.01
#> log_Base4_c           1.06      0.11     0.86     1.27        549 1.00
#> Trt_c                -0.34      0.16    -0.66    -0.03        515 1.01
#> log_Base4_c:Trt_c     0.34      0.22    -0.09     0.78        611 1.00
#> 
#> Samples were drawn using NUTS(diag_e). For each parameter, Eff.Sample is a 
#> crude measure of effective sample size, and Rhat is the potential scale 
#> reduction factor on split chains (at convergence, Rhat = 1).
```

To visually investigate the chains as well as the posterior, you can use

``` r
plot(fit) 
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
