# Priors of `brms` models

Extract priors of models fitted with brms.

## Usage

``` r
# S3 method for class 'brmsfit'
prior_summary(object, all = TRUE, ...)
```

## Arguments

- object:

  An object of class `brmsfit`.

- all:

  Logical; Show all parameters in the model which may have priors
  (`TRUE`) or only those with proper priors (`FALSE`)?

- ...:

  Further arguments passed to or from other methods.

## Value

An `brmsprior` object.

## Examples

``` r
# \dontrun{
fit <- brm(
  count ~ zAge + zBase * Trt + (1|patient) + (1|obs),
  data = epilepsy, family = poisson(),
  prior = prior(student_t(5,0,10), class = b) +
    prior(cauchy(0,2), class = sd)
)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 4.5e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.45 seconds.
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
#> Chain 1:  Elapsed Time: 3.913 seconds (Warm-up)
#> Chain 1:                3.264 seconds (Sampling)
#> Chain 1:                7.177 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 3.8e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.38 seconds.
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
#> Chain 2:  Elapsed Time: 3.916 seconds (Warm-up)
#> Chain 2:                2.331 seconds (Sampling)
#> Chain 2:                6.247 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 4e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.4 seconds.
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
#> Chain 3:  Elapsed Time: 3.803 seconds (Warm-up)
#> Chain 3:                2.333 seconds (Sampling)
#> Chain 3:                6.136 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 3.9e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.39 seconds.
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
#> Chain 4:  Elapsed Time: 3.843 seconds (Warm-up)
#> Chain 4:                2.318 seconds (Sampling)
#> Chain 4:                6.161 seconds (Total)
#> Chain 4: 

prior_summary(fit)
#>                   prior     class       coef   group resp dpar nlpar lb ub tag
#>  student_t(3, 1.4, 2.5) Intercept                                             
#>     student_t(5, 0, 10)         b                                             
#>     student_t(5, 0, 10)         b       Trt1                                  
#>     student_t(5, 0, 10)         b       zAge                                  
#>     student_t(5, 0, 10)         b      zBase                                  
#>     student_t(5, 0, 10)         b zBase:Trt1                                  
#>            cauchy(0, 2)        sd                                     0       
#>            cauchy(0, 2)        sd                obs                  0       
#>            cauchy(0, 2)        sd  Intercept     obs                  0       
#>            cauchy(0, 2)        sd            patient                  0       
#>            cauchy(0, 2)        sd  Intercept patient                  0       
#>        source
#>       default
#>          user
#>  (vectorized)
#>  (vectorized)
#>  (vectorized)
#>  (vectorized)
#>          user
#>  (vectorized)
#>  (vectorized)
#>  (vectorized)
#>  (vectorized)
prior_summary(fit, all = FALSE)
#>                   prior     class coef group resp dpar nlpar lb ub tag  source
#>  student_t(3, 1.4, 2.5) Intercept                                      default
#>     student_t(5, 0, 10)         b                                         user
#>            cauchy(0, 2)        sd                             0           user
print(prior_summary(fit, all = FALSE), show_df = FALSE)
#> Intercept ~ student_t(3, 1.4, 2.5)
#> b ~ student_t(5, 0, 10)
#> <lower=0> sd ~ cauchy(0, 2)
# }
```
