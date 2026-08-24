# Trace and Density Plots for MCMC Draws

Trace and Density Plots for MCMC Draws

## Usage

``` r
# S3 method for class 'brmsfit'
plot(
  x,
  pars = NA,
  combo = c("hist", "trace"),
  nvariables = 5,
  N = NULL,
  variable = NULL,
  regex = FALSE,
  fixed = FALSE,
  bins = 30,
  theme = NULL,
  plot = TRUE,
  ask = TRUE,
  newpage = TRUE,
  ...
)
```

## Arguments

- x:

  An object of class `brmsfit`.

- pars:

  Deprecated alias of `variable`. Names of the parameters to plot, as
  given by a character vector or a regular expression.

- combo:

  A character vector with at least two elements. Each element of `combo`
  corresponds to a column in the resulting graphic and should be the
  name of one of the available
  [`MCMC`](https://mc-stan.org/bayesplot/reference/MCMC-overview.html)
  functions (omitting the `mcmc_` prefix).

- nvariables:

  The number of variables (parameters) plotted per page.

- N:

  Deprecated alias of `nvariables`.

- variable:

  Names of the variables (parameters) to plot, as given by a character
  vector or a regular expression (if `regex = TRUE`). By default, a
  hopefully not too large selection of variables is plotted.

- regex:

  Logical; Indicates whether `variable` should be treated as regular
  expressions. Defaults to `FALSE`.

- fixed:

  (Deprecated) Indicates whether parameter names should be matched
  exactly (`TRUE`) or treated as regular expressions (`FALSE`). Default
  is `FALSE` and only works with argument `pars`.

- bins:

  Number of bins used for posterior histograms (defaults to 30).

- theme:

  A [`theme`](https://ggplot2.tidyverse.org/reference/theme.html) object
  modifying the appearance of the plots. For some basic themes see
  `ggtheme` and
  [`theme_default`](https://mc-stan.org/bayesplot/reference/theme_default.html).

- plot:

  Logical; indicates if plots should be plotted directly in the active
  graphic device. Defaults to `TRUE`.

- ask:

  Logical; indicates if the user is prompted before a new page is
  plotted. Only used if `plot` is `TRUE`.

- newpage:

  Logical; indicates if the first set of plots should be plotted to a
  new page. Only used if `plot` is `TRUE`.

- ...:

  Further arguments passed to
  [`mcmc_combo`](https://mc-stan.org/bayesplot/reference/MCMC-combos.html).

## Value

An invisible list of
[`gtable`](https://gtable.r-lib.org/reference/gtable.html) objects.

## Examples

``` r
# \dontrun{
fit <- brm(count ~ zAge + zBase * Trt
           + (1|patient) + (1|visit),
           data = epilepsy, family = "poisson")
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 4.2e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.42 seconds.
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
#> Chain 1:  Elapsed Time: 2.766 seconds (Warm-up)
#> Chain 1:                1.861 seconds (Sampling)
#> Chain 1:                4.627 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 3.4e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.34 seconds.
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
#> Chain 2:  Elapsed Time: 2.711 seconds (Warm-up)
#> Chain 2:                2.107 seconds (Sampling)
#> Chain 2:                4.818 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 3.5e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.35 seconds.
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
#> Chain 3:  Elapsed Time: 2.475 seconds (Warm-up)
#> Chain 3:                1.789 seconds (Sampling)
#> Chain 3:                4.264 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 3.5e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.35 seconds.
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
#> Chain 4:  Elapsed Time: 2.526 seconds (Warm-up)
#> Chain 4:                1.884 seconds (Sampling)
#> Chain 4:                4.41 seconds (Total)
#> Chain 4: 
#> Warning: There were 1 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
plot(fit)


## plot population-level effects only
plot(fit, variable = "^b_", regex = TRUE)

# }
```
