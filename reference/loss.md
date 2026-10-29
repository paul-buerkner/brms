# Cumulative Insurance Loss Payments

This dataset, discussed in Gesmann & Morris (2020), contains cumulative
insurance loss payments over the course of ten years.

## Usage

``` r
loss
```

## Format

A data frame of 55 observations containing information on the following
4 variables.

- AY:

  Origin year of the insurance (1991 to 2000)

- dev:

  Deviation from the origin year in months

- cum:

  Cumulative loss payments

- premium:

  Achieved premiums for the given origin year

## Source

Gesmann M. & Morris J. (2020). Hierarchical Compartmental Reserving
Models. *CAS Research Papers*.

## Examples

``` r
# \dontrun{
# non-linear model to predict cumulative loss payments
fit_loss <- brm(
  bf(cum ~ ult * (1 - exp(-(dev/theta)^omega)),
     ult ~ 1 + (1|AY), omega ~ 1, theta ~ 1,
     nl = TRUE),
  data = loss, family = gaussian(),
  prior = c(
    prior(normal(5000, 1000), nlpar = "ult"),
    prior(normal(1, 2), nlpar = "omega"),
    prior(normal(45, 10), nlpar = "theta")
  ),
  control = list(adapt_delta = 0.9)
)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 3.3e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.33 seconds.
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
#> Chain 1:  Elapsed Time: 2.031 seconds (Warm-up)
#> Chain 1:                0.979 seconds (Sampling)
#> Chain 1:                3.01 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: Rejecting initial value:
#> Chain 2:   Error evaluating the log probability at the initial value.
#> Chain 2: Exception: normal_lpdf: Location parameter[1] is nan, but must be finite! (in 'anon_model', line 68, column 4 to column 41)
#> Chain 2: 
#> Chain 2: Gradient evaluation took 2.4e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.24 seconds.
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
#> Chain 2:  Elapsed Time: 2.118 seconds (Warm-up)
#> Chain 2:                1.132 seconds (Sampling)
#> Chain 2:                3.25 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: Rejecting initial value:
#> Chain 3:   Error evaluating the log probability at the initial value.
#> Chain 3: Exception: normal_lpdf: Location parameter[1] is nan, but must be finite! (in 'anon_model', line 68, column 4 to column 41)
#> Chain 3: Rejecting initial value:
#> Chain 3:   Error evaluating the log probability at the initial value.
#> Chain 3: Exception: normal_lpdf: Location parameter[1] is nan, but must be finite! (in 'anon_model', line 68, column 4 to column 41)
#> Chain 3: 
#> Chain 3: Gradient evaluation took 2.3e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.23 seconds.
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
#> Chain 3:  Elapsed Time: 1.938 seconds (Warm-up)
#> Chain 3:                0.928 seconds (Sampling)
#> Chain 3:                2.866 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: Rejecting initial value:
#> Chain 4:   Error evaluating the log probability at the initial value.
#> Chain 4: Exception: normal_lpdf: Location parameter[1] is nan, but must be finite! (in 'anon_model', line 68, column 4 to column 41)
#> Chain 4: Rejecting initial value:
#> Chain 4:   Error evaluating the log probability at the initial value.
#> Chain 4: Exception: normal_lpdf: Location parameter[1] is nan, but must be finite! (in 'anon_model', line 68, column 4 to column 41)
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
#> Chain 4:  Elapsed Time: 2.167 seconds (Warm-up)
#> Chain 4:                1.177 seconds (Sampling)
#> Chain 4:                3.344 seconds (Total)
#> Chain 4: 

# basic summaries
summary(fit_loss)
#>  Family: gaussian 
#>   Links: mu = identity 
#> Formula: cum ~ ult * (1 - exp(-(dev/theta)^omega)) 
#>          ult ~ 1 + (1 | AY)
#>          omega ~ 1
#>          theta ~ 1
#>    Data: loss (Number of observations: 55) 
#>   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
#>          total post-warmup draws = 4000
#> 
#> Multilevel Hyperparameters:
#> ~AY (Number of levels: 10) 
#>                   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sd(ult_Intercept)   732.58    222.78   433.87  1284.77 1.00     1162     2052
#> 
#> Regression Coefficients:
#>                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> ult_Intercept    5300.74    294.32  4742.62  5905.04 1.00     1062     1757
#> omega_Intercept     1.34      0.05     1.24     1.44 1.00     2324     2280
#> theta_Intercept    46.24      2.18    42.39    50.92 1.00     2083     2139
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma   140.41     15.94   113.81   176.04 1.00     2569     2446
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
conditional_effects(fit_loss)


# plot predictions per origin year
conditions <- data.frame(AY = unique(loss$AY))
rownames(conditions) <- unique(loss$AY)
me_loss <- conditional_effects(
  fit_loss, conditions = conditions,
  re_formula = NULL, method = "predict"
)
plot(me_loss, ncol = 5, points = TRUE)

# }
```
