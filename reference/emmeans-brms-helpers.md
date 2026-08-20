# Support Functions for emmeans

Functions required for compatibility of brms with emmeans. Users are not
required to call these functions themselves. Instead, they will be
called automatically by the `emmeans` function of the emmeans package.

## Usage

``` r
recover_data.brmsfit(
  object,
  data,
  resp = NULL,
  dpar = NULL,
  nlpar = NULL,
  re_formula = NA,
  epred = FALSE,
  ...
)

emm_basis.brmsfit(
  object,
  trms,
  xlev,
  grid,
  vcov.,
  resp = NULL,
  dpar = NULL,
  nlpar = NULL,
  re_formula = NA,
  epred = FALSE,
  ...
)
```

## Arguments

- object:

  An object of class `brmsfit`.

- data, trms, xlev, grid, vcov.:

  Arguments required by emmeans.

- resp:

  Optional names of response variables. If specified, predictions are
  performed only for the specified response variables.

- dpar:

  Optional name of a predicted distributional parameter. If specified,
  expected predictions of this parameters are returned.

- nlpar:

  Optional name of a predicted non-linear parameter. If specified,
  expected predictions of this parameters are returned.

- re_formula:

  Optional formula containing group-level effects to be considered in
  the prediction. If `NULL`, include all group-level effects; if `NA`
  (default), include no group-level effects.

- epred:

  Logical. If `TRUE` compute predictions of the posterior predictive
  distribution's mean (see
  [`posterior_epred.brmsfit`](https://paulbuerkner.com/brms/reference/posterior_epred.brmsfit.md))
  while ignoring arguments `dpar` and `nlpar`. Defaults to `FALSE`. If
  you have specified a response transformation within the formula, you
  need to set `epred` to `TRUE` for emmeans to detect this
  transformation.

- ...:

  Additional arguments passed to emmeans.

## Details

In order to ensure compatibility of most brms models with emmeans,
predictions are not generated 'manually' via a design matrix and
coefficient vector, but rather via
[`posterior_linpred.brmsfit`](https://paulbuerkner.com/brms/reference/posterior_linpred.brmsfit.md).
This appears to generally work well, but note that it produces an
\`.@linfct\` slot that contains the computed predictions as columns
instead of the coefficients.

## Examples

``` r
# \dontrun{
fit1 <- brm(time | cens(censored) ~ age * sex + disease + (1|patient),
            data = kidney, family = lognormal())
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: Rejecting initial value:
#> Chain 1:   Log probability evaluates to log(0), i.e. negative infinity.
#> Chain 1:   Stan can't start sampling from this initial value.
#> Chain 1: Rejecting initial value:
#> Chain 1:   Log probability evaluates to log(0), i.e. negative infinity.
#> Chain 1:   Stan can't start sampling from this initial value.
#> Chain 1: 
#> Chain 1: Gradient evaluation took 2.8e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.28 seconds.
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
#> Chain 1:  Elapsed Time: 1.062 seconds (Warm-up)
#> Chain 1:                0.384 seconds (Sampling)
#> Chain 1:                1.446 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: Rejecting initial value:
#> Chain 2:   Log probability evaluates to log(0), i.e. negative infinity.
#> Chain 2:   Stan can't start sampling from this initial value.
#> Chain 2: Rejecting initial value:
#> Chain 2:   Log probability evaluates to log(0), i.e. negative infinity.
#> Chain 2:   Stan can't start sampling from this initial value.
#> Chain 2: Rejecting initial value:
#> Chain 2:   Log probability evaluates to log(0), i.e. negative infinity.
#> Chain 2:   Stan can't start sampling from this initial value.
#> Chain 2: Rejecting initial value:
#> Chain 2:   Log probability evaluates to log(0), i.e. negative infinity.
#> Chain 2:   Stan can't start sampling from this initial value.
#> Chain 2: 
#> Chain 2: Gradient evaluation took 2e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.2 seconds.
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
#> Chain 2:  Elapsed Time: 0.947 seconds (Warm-up)
#> Chain 2:                0.373 seconds (Sampling)
#> Chain 2:                1.32 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 2e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.2 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 0.96 seconds (Warm-up)
#> Chain 3:                0.385 seconds (Sampling)
#> Chain 3:                1.345 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: Rejecting initial value:
#> Chain 4:   Log probability evaluates to log(0), i.e. negative infinity.
#> Chain 4:   Stan can't start sampling from this initial value.
#> Chain 4: 
#> Chain 4: Gradient evaluation took 2.6e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.26 seconds.
#> Chain 4: Adjust your expectations accordingly!
#> Chain 4: 
#> Chain 4: 
#> Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 4: 
#> Chain 4:  Elapsed Time: 0.949 seconds (Warm-up)
#> Chain 4:                0.257 seconds (Sampling)
#> Chain 4:                1.206 seconds (Total)
#> Chain 4: 
summary(fit1)
#>  Family: lognormal 
#>   Links: mu = identity 
#> Formula: time | cens(censored) ~ age * sex + disease + (1 | patient) 
#>    Data: kidney (Number of observations: 76) 
#>   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 4000
#> 
#> Multilevel Hyperparameters:
#> ~patient (Number of levels: 38) 
#>               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sd(Intercept)     0.41      0.25     0.02     0.94 1.01      845     1607
#> 
#> Regression Coefficients:
#>               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept         2.61      0.98     0.70     4.56 1.00     2415     2877
#> age               0.02      0.02    -0.03     0.07 1.00     1885     2280
#> sexfemale         2.59      1.16     0.33     4.80 1.00     2261     2620
#> diseaseGN        -0.45      0.53    -1.50     0.61 1.00     2367     2761
#> diseaseAN        -0.54      0.51    -1.55     0.46 1.00     2599     2653
#> diseasePKD        0.57      0.70    -0.83     1.97 1.00     2496     2352
#> age:sexfemale    -0.03      0.03    -0.08     0.03 1.00     2244     2393
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma     1.16      0.13     0.93     1.44 1.00     2225     2357
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).

# summarize via 'emmeans'
library(emmeans)
#> Welcome to emmeans.
#> Caution: You lose important information if you filter this package's results.
#> See '? untidy'
rg <- ref_grid(fit1)
em <- emmeans(rg, "disease")
summary(em, point.est = mean)
#>  disease emmean lower.HPD upper.HPD
#>  other     4.13      3.41      4.83
#>  GN        3.68      2.99      4.40
#>  AN        3.58      2.93      4.24
#>  PKD       4.70      3.64      5.80
#> 
#> Results are averaged over the levels of: sex 
#> Point estimate displayed: mean 
#> HPD interval probability: 0.95 

# obtain estimates for the posterior predictive distribution's mean
epred <- emmeans(fit1, "disease", epred = TRUE)
summary(epred, point.est = mean)
#>  disease emmean lower.HPD upper.HPD
#>  other    172.5      63.5       322
#>  GN       110.6      38.7       206
#>  AN        98.7      41.3       175
#>  PKD      336.2      53.9       752
#> 
#> Results are averaged over the levels of: censored, sex 
#> Point estimate displayed: mean 
#> HPD interval probability: 0.95 


# model with transformed response variable
fit2 <- brm(log(mpg) ~ factor(cyl), data = mtcars)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 8e-06 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.08 seconds.
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
#> Chain 1:  Elapsed Time: 0.013 seconds (Warm-up)
#> Chain 1:                0.011 seconds (Sampling)
#> Chain 1:                0.024 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 4e-06 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.04 seconds.
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
#> Chain 2:  Elapsed Time: 0.013 seconds (Warm-up)
#> Chain 2:                0.011 seconds (Sampling)
#> Chain 2:                0.024 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 4e-06 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.04 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 0.013 seconds (Warm-up)
#> Chain 3:                0.011 seconds (Sampling)
#> Chain 3:                0.024 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 4e-06 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.04 seconds.
#> Chain 4: Adjust your expectations accordingly!
#> Chain 4: 
#> Chain 4: 
#> Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 4: 
#> Chain 4:  Elapsed Time: 0.014 seconds (Warm-up)
#> Chain 4:                0.011 seconds (Sampling)
#> Chain 4:                0.025 seconds (Total)
#> Chain 4: 
summary(fit2)
#>  Family: gaussian 
#>   Links: mu = identity 
#> Formula: log(mpg) ~ factor(cyl) 
#>    Data: mtcars (Number of observations: 32) 
#>   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 4000
#> 
#> Regression Coefficients:
#>            Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept      3.27      0.05     3.17     3.37 1.00     2785     1847
#> factorcyl6    -0.29      0.08    -0.45    -0.13 1.00     2885     2749
#> factorcyl8    -0.57      0.07    -0.70    -0.44 1.00     3295     2639
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma     0.17      0.02     0.13     0.22 1.00     3005     2483
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).

# results will be on the log scale by default
emmeans(fit2, ~ cyl)
#>  cyl emmean lower.HPD upper.HPD
#>    4   3.27      3.16      3.37
#>    6   2.98      2.86      3.11
#>    8   2.70      2.61      2.79
#> 
#> Point estimate displayed: median 
#> HPD interval probability: 0.95 
# log transform is detected and can be adjusted automatically
emmeans(fit2, ~ cyl, epred = TRUE, type = "response")
#>  cyl response lower.HPD upper.HPD
#>    4     26.4      23.6      28.9
#>    6     19.7      17.4      22.2
#>    8     14.9      13.6      16.2
#> 
#> Point estimate displayed: median 
#> Results are back-transformed from the log scale 
#> HPD interval probability: 0.95 
# }
```
