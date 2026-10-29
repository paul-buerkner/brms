# Finite Mixture Families in brms

Set up a finite mixture family for use in brms.

## Usage

``` r
mixture(..., flist = NULL, nmix = 1, order = NULL, refcat = NULL)
```

## Arguments

- ...:

  One or more objects providing a description of the response
  distributions to be combined in the mixture model. These can be family
  functions, calls to family functions or character strings naming the
  families. For details of supported families see
  [`brmsfamily`](https://paulbuerkner.com/brms/reference/brmsfamily.md).

- flist:

  Optional list of objects, which are treated in the same way as objects
  passed via the `...` argument.

- nmix:

  Optional numeric vector specifying the number of times each family is
  repeated. If specified, it must have the same length as the number of
  families passed via `...` and `flist`.

- order:

  Ordering constraint to identify mixture components. If `'mu'` or
  `TRUE`, population-level intercepts of the mean parameters are ordered
  in non-ordinal models and fixed to the same value in ordinal models
  (see details). If `'none'` or `FALSE`, no ordering constraint is
  applied. If `NULL` (the default), `order` is set to `'mu'` if all
  families are the same and `'none'` otherwise. Other ordering
  constraints may be implemented in the future.

- refcat:

  Optional reference category for the mixing proportions. By default
  (`NULL`), when the mixing proportions are predicted, all but one of
  them have to be modeled while the remaining one serves as the
  reference category (its linear predictor is fixed to zero) to identify
  the model. If `refcat = NA`, no reference category is used and all
  mixing proportions are predicted via a 'softmax' transformation. This
  is analogous to `refcat = NA` in
  [`categorical`](https://paulbuerkner.com/brms/reference/brmsfamily.md)
  models. As the resulting model is only weakly identified, informative
  priors on the `theta` parameters are strongly recommended (see the
  examples in
  [`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md)).

## Value

An object of class `mixfamily`.

## Details

Most families supported by brms can be used to form mixtures. The
response variable has to be valid for all components of the mixture
family. Currently, the number of mixture components has to be specified
by the user. It is not yet possible to estimate the number of mixture
components from the data.

Ordering intercepts in mixtures of ordinal families is not possible as
each family has itself a set of vector of intercepts (i.e. ordinal
thresholds). Instead, brms will fix the vector of intercepts across
components in ordinal mixtures, if desired, so that users can try to
identify the mixture model via selective inclusion of predictors.

For most mixture models, you may want to specify priors on the
population-level intercepts via
[`set_prior`](https://paulbuerkner.com/brms/reference/set_prior.md) to
improve posterior inference. In addition, it is sometimes necessary to
set `init = 0` in the call to
[`brm`](https://paulbuerkner.com/brms/reference/brm.md) to allow chains
to initialize properly.

For more details on the specification of mixture models, see
[`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md).

## Examples

``` r
# \dontrun{
## simulate some data
set.seed(1234)
dat <- data.frame(
  y = c(rnorm(200), rnorm(100, 6)),
  x = rnorm(300),
  z = sample(0:1, 300, TRUE)
)

## fit a simple normal mixture model
mix <- mixture(gaussian, gaussian)
#> Setting order = 'mu' for mixtures of the same family.
prior <- c(
  prior(normal(0, 7), Intercept, dpar = mu1),
  prior(normal(5, 7), Intercept, dpar = mu2)
)
fit1 <- brm(bf(y ~ x + z), dat, family = mix,
            prior = prior, chains = 2)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.00018 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 1.8 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 1.291 seconds (Warm-up)
#> Chain 1:                0.954 seconds (Sampling)
#> Chain 1:                2.245 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.000124 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 1.24 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 25.671 seconds (Warm-up)
#> Chain 2:                28.67 seconds (Sampling)
#> Chain 2:                54.341 seconds (Total)
#> Chain 2: 
#> Warning: There were 8 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: There were 167 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
#> https://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: The largest R-hat is 1.9, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
summary(fit1)
#> Warning: Inference for the model posterior has not converged (some Rhats are > 1.05). Be careful when analysing the results! We recommend running more iterations or setting stronger priors.
#> Warning: There were 8 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#>  Family: mixture(gaussian, gaussian) 
#>   Links: mu1 = identity; mu2 = identity 
#> Formula: y ~ x + z 
#>    Data: dat (Number of observations: 300) 
#>   Draws: 2 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 2000
#> 
#> Regression Coefficients:
#>                         Estimate          Est.Error            l-95% CI
#> mu1_Intercept               1.28               1.24               -0.13
#> mu2_Intercept -62080756601433.02 212486694935603.22 -498758314800336.06
#> mu1_x                      -0.04               0.16               -0.43
#> mu1_z                      -0.54               0.46               -1.50
#> mu2_x         -41296644183408.45  97568956732815.98 -282145670377849.75
#> mu2_z         109362994232676.08 390519451599940.69 -798256814877106.75
#>                         u-95% CI Rhat Bulk_ESS Tail_ESS
#> mu1_Intercept               2.94 1.83        3       88
#> mu2_Intercept 427193169573408.06 1.87        4       26
#> mu1_x                       0.21 1.23        6       98
#> mu1_z                       0.08 1.66        3       82
#> mu2_x         124548944645316.80 1.90       13       27
#> mu2_z         907569476893578.50 1.87        4       26
#> 
#> Further Distributional Parameters:
#>        Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma1     2.06      1.03     0.95     3.29 1.83        3       58
#> sigma2     1.73      2.00     0.15     7.41 1.60       18       61
#> theta1     0.83      0.17     0.62     1.00 1.83        3       67
#> theta2     0.17      0.17     0.00     0.38 1.83        3       67
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
pp_check(fit1)
#> Using 10 posterior draws for ppc type 'dens_overlay' by default.


## use different predictors for the components
fit2 <- brm(bf(y ~ 1, mu1 ~ x, mu2 ~ z), dat, family = mix,
            prior = prior, chains = 2)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000154 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 1.54 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 1.138 seconds (Warm-up)
#> Chain 1:                0.914 seconds (Sampling)
#> Chain 1:                2.052 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.000153 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 1.53 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 1.289 seconds (Warm-up)
#> Chain 2:                0.942 seconds (Sampling)
#> Chain 2:                2.231 seconds (Total)
#> Chain 2: 
#> Warning: The largest R-hat is 1.83, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
summary(fit2)
#> Warning: Inference for the model posterior has not converged (some Rhats are > 1.05). Be careful when analysing the results! We recommend running more iterations or setting stronger priors.
#>  Family: mixture(gaussian, gaussian) 
#>   Links: mu1 = identity; mu2 = identity 
#> Formula: y ~ 1 
#>          mu1 ~ x
#>          mu2 ~ z
#>    Data: dat (Number of observations: 300) 
#>   Draws: 2 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 2000
#> 
#> Regression Coefficients:
#>               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> mu1_Intercept    -0.04      0.09    -0.22     0.14 1.02     2253     1163
#> mu2_Intercept     3.19      3.04    -0.06     6.45 1.83        3       75
#> mu1_x            -0.10      0.33    -0.91     0.40 1.36       12       67
#> mu2_z            -0.18      0.20    -0.54     0.22 1.10       13      139
#> 
#> Further Distributional Parameters:
#>        Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma1     3.03      2.00     0.96     5.59 1.83        3       84
#> sigma2     0.85      0.11     0.65     1.06 1.54        4       70
#> theta1     0.60      0.08     0.47     0.71 1.80        3       81
#> theta2     0.40      0.08     0.29     0.53 1.80        3       81
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).

## fix the mixing proportions
fit3 <- brm(bf(y ~ x + z, theta1 = 1, theta2 = 2),
            dat, family = mix, prior = prior,
            init = 0, chains = 2)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000121 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 1.21 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.937 seconds (Warm-up)
#> Chain 1:                0.858 seconds (Sampling)
#> Chain 1:                1.795 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.00011 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 1.1 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 0.916 seconds (Warm-up)
#> Chain 2:                0.79 seconds (Sampling)
#> Chain 2:                1.706 seconds (Total)
#> Chain 2: 
summary(fit3)
#>  Family: mixture(gaussian, gaussian) 
#>   Links: mu1 = identity; mu2 = identity 
#> Formula: y ~ x + z 
#>          theta1 = 0.33
#>          theta2 = 0.67
#>    Data: dat (Number of observations: 300) 
#>   Draws: 2 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 2000
#> 
#> Regression Coefficients:
#>               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> mu1_Intercept    -0.19      0.11    -0.41     0.04 1.00     2974     1510
#> mu2_Intercept     4.49      0.42     3.67     5.30 1.00     3074     1837
#> mu1_x            -0.01      0.08    -0.15     0.14 1.00     3494     1428
#> mu1_z            -0.23      0.15    -0.52     0.07 1.00     2820     1408
#> mu2_x            -0.22      0.28    -0.76     0.34 1.00     2652     1640
#> mu2_z            -1.28      0.57    -2.39    -0.14 1.00     3327     1683
#> 
#> Further Distributional Parameters:
#>        Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma1     0.68      0.06     0.57     0.81 1.00     3682     1787
#> sigma2     2.98      0.19     2.62     3.37 1.00     3157     1662
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
pp_check(fit3)
#> Using 10 posterior draws for ppc type 'dens_overlay' by default.


## predict the mixing proportions
fit4 <- brm(bf(y ~ x + z, theta2 ~ x),
            dat, family = mix, prior = prior,
            init = 0, chains = 2)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000219 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 2.19 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 1.87 seconds (Warm-up)
#> Chain 1:                1.418 seconds (Sampling)
#> Chain 1:                3.288 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.000179 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 1.79 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 1.936 seconds (Warm-up)
#> Chain 2:                1.41 seconds (Sampling)
#> Chain 2:                3.346 seconds (Total)
#> Chain 2: 
summary(fit4)
#>  Family: mixture(gaussian, gaussian) 
#>   Links: mu1 = identity; mu2 = identity; theta2 = identity 
#> Formula: y ~ x + z 
#>          theta2 ~ x
#>    Data: dat (Number of observations: 300) 
#>   Draws: 2 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 2000
#> 
#> Regression Coefficients:
#>                  Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> mu1_Intercept        0.05      0.12    -0.19     0.29 1.00     2883     1629
#> mu2_Intercept        6.21      0.13     5.96     6.46 1.00     3036     1501
#> theta2_Intercept    -0.71      0.13    -0.95    -0.47 1.00     3983     1474
#> mu1_x                0.05      0.07    -0.08     0.18 1.00     2929     1557
#> mu1_z               -0.17      0.15    -0.46     0.13 1.00     3956     1735
#> mu2_x               -0.08      0.10    -0.29     0.13 1.00     4463     1517
#> mu2_z               -0.10      0.19    -0.50     0.26 1.00     4268     1658
#> theta2_x            -0.11      0.12    -0.35     0.12 1.00     4406     1552
#> 
#> Further Distributional Parameters:
#>        Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma1     1.04      0.05     0.94     1.15 1.00     4016     1390
#> sigma2     0.93      0.07     0.80     1.09 1.01     3609     1783
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
pp_check(fit4)
#> Using 10 posterior draws for ppc type 'dens_overlay' by default.


## predict all mixing proportions without a reference category
## (requires informative priors on the theta parameters)
mix_na <- mixture(gaussian, gaussian, refcat = NA)
#> Setting order = 'mu' for mixtures of the same family.
theta_prior <- c(prior, prior(normal(0, 1), Intercept, dpar = theta1),
                  prior(normal(0, 1), Intercept, dpar = theta2))
fit5 <- brm(bf(y ~ x + z, theta1 ~ x, theta2 ~ x),
            dat, family = mix_na, prior = theta_prior,
            init = 0, chains = 2)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000217 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 2.17 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 152.635 seconds (Warm-up)
#> Chain 1:                187.432 seconds (Sampling)
#> Chain 1:                340.067 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.000195 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 1.95 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 148.671 seconds (Warm-up)
#> Chain 2:                183.174 seconds (Sampling)
#> Chain 2:                331.845 seconds (Total)
#> Chain 2: 
#> Warning: There were 196 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: There were 1803 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
#> https://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: The largest R-hat is 1.06, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
summary(fit5)
#> Warning: Inference for the model posterior has not converged (some Rhats are > 1.05). Be careful when analysing the results! We recommend running more iterations or setting stronger priors.
#> Warning: There were 196 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#>  Family: mixture(gaussian, gaussian) 
#>   Links: mu1 = identity; mu2 = identity; theta1 = identity; theta2 = identity 
#> Formula: y ~ x + z 
#>          theta1 ~ x
#>          theta2 ~ x
#>    Data: dat (Number of observations: 300) 
#>   Draws: 2 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 2000
#> 
#> Regression Coefficients:
#>                  Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> mu1_Intercept        0.05      0.11    -0.18     0.28 1.03       61      153
#> mu2_Intercept        6.21      0.14     5.90     6.50 1.02       70      107
#> theta1_Intercept     1.03      6.89   -10.81    12.11 1.05       23      160
#> theta2_Intercept     0.32      6.89   -11.59    11.44 1.05       23      161
#> mu1_x                0.07      0.07    -0.06     0.21 1.01      102      181
#> mu1_z               -0.17      0.14    -0.42     0.11 1.02       85      206
#> mu2_x               -0.09      0.09    -0.27     0.08 1.02       76      203
#> mu2_z               -0.09      0.21    -0.54     0.31 1.02       49       71
#> theta1_x            12.52    119.70  -197.71   198.21 1.05       23      122
#> theta2_x            12.40    119.70  -197.79   197.99 1.05       23      122
#> 
#> Further Distributional Parameters:
#>        Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma1     1.04      0.05     0.94     1.14 1.02       78      117
#> sigma2     0.93      0.07     0.81     1.11 1.02       39       77
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
pp_check(fit5)
#> Using 10 posterior draws for ppc type 'dens_overlay' by default.


## compare model fit
loo(fit1, fit2, fit3, fit4)
#> Output of model 'fit1':
#> 
#> Computed from 2000 by 300 log-likelihood matrix.
#> 
#>          Estimate   SE
#> elpd_loo   -721.4 10.0
#> p_loo        62.1  3.8
#> looic      1442.7 20.0
#> ------
#> MCSE of elpd_loo is 5.0.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.0, 1.5]).
#> 
#> All Pareto k estimates are good (k < 0.7).
#> See help('pareto-k-diagnostic') for details.
#> 
#> Output of model 'fit2':
#> 
#> Computed from 2000 by 300 log-likelihood matrix.
#> 
#>          Estimate   SE
#> elpd_loo   -721.4 17.0
#> p_loo        71.9  5.8
#> looic      1442.8 34.0
#> ------
#> MCSE of elpd_loo is 4.7.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.0, 1.4]).
#> 
#> All Pareto k estimates are good (k < 0.7).
#> See help('pareto-k-diagnostic') for details.
#> 
#> Output of model 'fit3':
#> 
#> Computed from 2000 by 300 log-likelihood matrix.
#> 
#>          Estimate   SE
#> elpd_loo   -683.6 11.1
#> p_loo         7.6  0.6
#> looic      1367.3 22.2
#> ------
#> MCSE of elpd_loo is 0.1.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.8, 1.9]).
#> 
#> All Pareto k estimates are good (k < 0.7).
#> See help('pareto-k-diagnostic') for details.
#> 
#> Output of model 'fit4':
#> 
#> Computed from 2000 by 300 log-likelihood matrix.
#> 
#>          Estimate   SE
#> elpd_loo   -620.5 13.9
#> p_loo        10.3  1.0
#> looic      1241.1 27.7
#> ------
#> MCSE of elpd_loo is 0.1.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [1.2, 2.3]).
#> 
#> All Pareto k estimates are good (k < 0.7).
#> See help('pareto-k-diagnostic') for details.
#> 
#> Model comparisons:
#>  model elpd_diff se_diff p_worse diag_diff diag_elpd
#>   fit4       0.0     0.0      NA                    
#>   fit3     -63.1     9.8    1.00                    
#>   fit1    -100.8     5.9    1.00                    
#>   fit2    -100.8    10.2    1.00                    
# }
```
