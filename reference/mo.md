# Monotonic Predictors in brms Models

Specify a monotonic predictor term in brms. The function does not
evaluate its arguments – it exists purely to help set up a model.

## Usage

``` r
mo(x, id = NA)
```

## Arguments

- x:

  An integer variable or an ordered factor to be modeled as monotonic.

- id:

  Optional character string. All monotonic terms with the same `id`
  within one formula will be modeled as having the same simplex (shape)
  parameter vector. If all monotonic terms of the same predictor have
  the same `id`, the resulting predictions will be conditionally
  monotonic for all values of interacting covariates (Bürkner &
  Charpentier, 2020).

## Details

See Bürkner and Charpentier (2020) for the underlying theory. For
detailed documentation of the formula syntax used for monotonic terms,
see
[`help(brmsformula)`](https://paulbuerkner.com/brms/reference/brmsformula.md)
as well as
[`vignette("brms_monotonic")`](https://paulbuerkner.com/brms/articles/brms_monotonic.md).

## References

Bürkner P. C. & Charpentier E. (2020). Modeling Monotonic Effects of
Ordinal Predictors in Regression Models. British Journal of Mathematical
and Statistical Psychology. doi:10.1111/bmsp.12195

## See also

[`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md)

## Examples

``` r
# \dontrun{
# generate some data
income_options <- c("below_20", "20_to_40", "40_to_100", "greater_100")
income <- factor(sample(income_options, 100, TRUE),
                 levels = income_options, ordered = TRUE)
mean_ls <- c(30, 60, 70, 75)
ls <- mean_ls[income] + rnorm(100, sd = 7)
dat <- data.frame(income, ls)

# fit a simple monotonic model
fit1 <- brm(ls ~ mo(income), data = dat)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 2.6e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.26 seconds.
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
#> Chain 1:  Elapsed Time: 0.201 seconds (Warm-up)
#> Chain 1:                0.173 seconds (Sampling)
#> Chain 1:                0.374 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 1.9e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.19 seconds.
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
#> Chain 2:  Elapsed Time: 0.202 seconds (Warm-up)
#> Chain 2:                0.176 seconds (Sampling)
#> Chain 2:                0.378 seconds (Total)
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
#> Chain 3:  Elapsed Time: 0.196 seconds (Warm-up)
#> Chain 3:                0.187 seconds (Sampling)
#> Chain 3:                0.383 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 2.1e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.21 seconds.
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
#> Chain 4:  Elapsed Time: 0.197 seconds (Warm-up)
#> Chain 4:                0.178 seconds (Sampling)
#> Chain 4:                0.375 seconds (Total)
#> Chain 4: 
summary(fit1)
#>  Family: gaussian 
#>   Links: mu = identity 
#> Formula: ls ~ mo(income) 
#>    Data: dat (Number of observations: 100) 
#>   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 4000
#> 
#> Regression Coefficients:
#>           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept    29.41      1.47    26.50    32.31 1.00     2609     2083
#> moincome     15.09      0.68    13.78    16.40 1.00     2420     2137
#> 
#> Monotonic Simplex Parameters:
#>              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> moincome1[1]     0.64      0.04     0.57     0.71 1.00     3150     2489
#> moincome1[2]     0.28      0.04     0.20     0.37 1.00     3851     2229
#> moincome1[3]     0.08      0.04     0.01     0.15 1.00     2605     1540
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma     7.04      0.51     6.12     8.14 1.00     2876     2192
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
plot(fit1, N = 6)
#> Warning: Argument 'N' is deprecated. Please use argument 'nvariables' instead.

plot(conditional_effects(fit1), points = TRUE)


# model interaction with other variables
dat$x <- sample(c("a", "b", "c"), 100, TRUE)
fit2 <- brm(ls ~ mo(income)*x, data = dat)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 5.3e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.53 seconds.
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
#> Chain 1:  Elapsed Time: 0.917 seconds (Warm-up)
#> Chain 1:                0.654 seconds (Sampling)
#> Chain 1:                1.571 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 9.2e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.92 seconds.
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
#> Chain 2:  Elapsed Time: 0.874 seconds (Warm-up)
#> Chain 2:                0.561 seconds (Sampling)
#> Chain 2:                1.435 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 0.000105 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 1.05 seconds.
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
#> Chain 3:  Elapsed Time: 0.897 seconds (Warm-up)
#> Chain 3:                0.626 seconds (Sampling)
#> Chain 3:                1.523 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 4.8e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.48 seconds.
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
#> Chain 4:  Elapsed Time: 0.966 seconds (Warm-up)
#> Chain 4:                0.558 seconds (Sampling)
#> Chain 4:                1.524 seconds (Total)
#> Chain 4: 
summary(fit2)
#>  Family: gaussian 
#>   Links: mu = identity 
#> Formula: ls ~ mo(income) * x 
#>    Data: dat (Number of observations: 100) 
#>   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 4000
#> 
#> Regression Coefficients:
#>             Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept      29.57      1.95    25.90    33.43 1.00     2421     2337
#> xb              1.47      2.87    -4.40     6.92 1.00     1839     2140
#> xc             -2.26      2.84    -7.79     3.54 1.00     2093     2429
#> moincome       15.61      1.00    13.66    17.65 1.00     1738     2383
#> moincome:xb    -1.21      1.55    -4.26     1.96 1.00     1630     2151
#> moincome:xc    -0.45      1.52    -3.52     2.44 1.00     1981     2518
#> 
#> Monotonic Simplex Parameters:
#>                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> moincome1[1]        0.63      0.04     0.55     0.72 1.00     2861     2692
#> moincome1[2]        0.29      0.04     0.20     0.37 1.00     5405     3029
#> moincome1[3]        0.08      0.04     0.01     0.17 1.00     2450     1444
#> moincome:xb1[1]     0.29      0.22     0.01     0.79 1.00     3393     1977
#> moincome:xb1[2]     0.32      0.23     0.01     0.82 1.00     4500     2794
#> moincome:xb1[3]     0.39      0.25     0.02     0.89 1.00     3621     2213
#> moincome:xc1[1]     0.33      0.23     0.02     0.83 1.00     3992     2011
#> moincome:xc1[2]     0.32      0.23     0.01     0.82 1.00     4411     2297
#> moincome:xc1[3]     0.35      0.23     0.01     0.84 1.00     4811     2943
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma     7.03      0.52     6.11     8.12 1.00     4805     3196
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
plot(conditional_effects(fit2), points = TRUE)




# ensure conditional monotonicity
fit3 <- brm(ls ~ mo(income, id = "i")*x, data = dat)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 5.2e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.52 seconds.
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
#> Chain 1:  Elapsed Time: 0.962 seconds (Warm-up)
#> Chain 1:                0.797 seconds (Sampling)
#> Chain 1:                1.759 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 4.7e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.47 seconds.
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
#> Chain 2:  Elapsed Time: 0.967 seconds (Warm-up)
#> Chain 2:                0.766 seconds (Sampling)
#> Chain 2:                1.733 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 4.6e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.46 seconds.
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
#> Chain 3:  Elapsed Time: 0.968 seconds (Warm-up)
#> Chain 3:                0.66 seconds (Sampling)
#> Chain 3:                1.628 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 4.6e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.46 seconds.
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
#> Chain 4:  Elapsed Time: 0.937 seconds (Warm-up)
#> Chain 4:                0.7 seconds (Sampling)
#> Chain 4:                1.637 seconds (Total)
#> Chain 4: 
summary(fit3)
#>  Family: gaussian 
#>   Links: mu = identity 
#> Formula: ls ~ mo(income, id = "i") * x 
#>    Data: dat (Number of observations: 100) 
#>   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 4000
#> 
#> Regression Coefficients:
#>                  Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> Intercept           30.12      2.06    26.17    34.29 1.00     2178     2741
#> xb                   0.19      3.42    -6.59     6.93 1.00     1969     2571
#> xc                  -3.11      3.32    -9.74     3.44 1.00     1966     2331
#> moincomeidEQi       15.14      0.98    13.18    17.03 1.00     1719     2195
#> moincomeidEQi:xb    -0.16      1.51    -3.17     2.79 1.00     1830     2389
#> moincomeidEQi:xc     0.11      1.46    -2.73     2.93 1.00     1825     2117
#> 
#> Monotonic Simplex Parameters:
#>                      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS
#> moincomeidEQi1[1]        0.65      0.04     0.57     0.72 1.00     3483
#> moincomeidEQi1[2]        0.28      0.04     0.20     0.36 1.00     4371
#> moincomeidEQi1[3]        0.07      0.04     0.01     0.15 1.00     2525
#> moincomeidEQi:xb1[1]     0.65      0.04     0.57     0.72 1.00     3483
#> moincomeidEQi:xb1[2]     0.28      0.04     0.20     0.36 1.00     4371
#> moincomeidEQi:xb1[3]     0.07      0.04     0.01     0.15 1.00     2525
#> moincomeidEQi:xc1[1]     0.65      0.04     0.57     0.72 1.00     3483
#> moincomeidEQi:xc1[2]     0.28      0.04     0.20     0.36 1.00     4371
#> moincomeidEQi:xc1[3]     0.07      0.04     0.01     0.15 1.00     2525
#>                      Tail_ESS
#> moincomeidEQi1[1]        2729
#> moincomeidEQi1[2]        2994
#> moincomeidEQi1[3]        1552
#> moincomeidEQi:xb1[1]     2729
#> moincomeidEQi:xb1[2]     2994
#> moincomeidEQi:xb1[3]     1552
#> moincomeidEQi:xc1[1]     2729
#> moincomeidEQi:xc1[2]     2994
#> moincomeidEQi:xc1[3]     1552
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma     7.07      0.53     6.13     8.16 1.00     3767     3093
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
plot(conditional_effects(fit3), points = TRUE)



# }
```
