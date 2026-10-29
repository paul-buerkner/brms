# Non-Linear Hypothesis Testing

Perform non-linear hypothesis testing for all model parameters.

## Usage

``` r
# S3 method for class 'brmsfit'
hypothesis(
  x,
  hypothesis,
  class = "b",
  group = "",
  scope = c("standard", "ranef", "coef"),
  alpha = 0.05,
  robust = FALSE,
  seed = NULL,
  ...
)

hypothesis(x, ...)

# Default S3 method
hypothesis(x, hypothesis, alpha = 0.05, robust = FALSE, ...)
```

## Arguments

- x:

  An `R` object. If it is no `brmsfit` object, it must be coercible to a
  `data.frame`. In the latter case, the variables used in the
  `hypothesis` argument need to correspond to column names of `x`, while
  the rows are treated as representing posterior draws of the variables.

- hypothesis:

  A character vector specifying one or more non-linear hypothesis
  concerning parameters of the model.

- class:

  A string specifying the class of parameters being tested. Default is
  "b" for population-level effects. Other typical options are "sd" or
  "cor". If `class = NULL`, all parameters can be tested against each
  other, but have to be specified with their full name (see also
  [`variables`](https://paulbuerkner.com/brms/reference/draws-index-brms.md))

- group:

  Name of a grouping factor to evaluate only group-level effects
  parameters related to this grouping factor.

- scope:

  Indicates where to look for the variables specified in `hypothesis`.
  If `"standard"`, use the full parameter names (subject to the
  restriction given by `class` and `group`). If `"coef"` or `"ranef"`,
  compute the hypothesis for all levels of the grouping factor given in
  `"group"`, based on the output of
  [`coef.brmsfit`](https://paulbuerkner.com/brms/reference/coef.brmsfit.md)
  and
  [`ranef.brmsfit`](https://paulbuerkner.com/brms/reference/ranef.brmsfit.md),
  respectively.

- alpha:

  The alpha-level of the tests (default is 0.05; see 'Details' for more
  information).

- robust:

  If `FALSE` (the default) the mean is used as the measure of central
  tendency and the standard deviation as the measure of variability. If
  `TRUE`, the median and the median absolute deviation (MAD) are applied
  instead.

- seed:

  A single numeric value passed to
  [`set.seed`](https://rdrr.io/r/base/Random.html) to make results
  reproducible. This is currently only relevant for point hypotheses
  that scope over at least two parameters (see Details).

- ...:

  Currently ignored.

## Value

A
[`brmshypothesis`](https://paulbuerkner.com/brms/reference/brmshypothesis.md)
object.

## Details

Among others, `hypothesis` computes an evidence ratio (`Evid.Ratio`) for
each hypothesis. For a one-sided hypothesis, this is just the posterior
probability (`Post.Prob`) under the hypothesis against its alternative.
That is, when the hypothesis is of the form `a > b`, the evidence ratio
is the ratio of the posterior probability of `a > b` and the posterior
probability of `a < b`. In this example, values greater than one
indicate that the evidence in favor of `a > b` is larger than evidence
in favor of `a < b`. For a two-sided (point) hypothesis, the evidence
ratio is a Bayes factor between the hypothesis and its alternative
computed via the Savage-Dickey density ratio method. That is the
posterior density at the point of interest divided by the prior density
at that point. Values greater than one indicate that evidence in favor
of the point hypothesis has increased after seeing the data. In order to
calculate this Bayes factor, all parameters related to the hypothesis
must have proper priors and argument `sample_prior` of function `brm`
must be set to `"yes"`. Otherwise `Evid.Ratio` (and `Post.Prob`) will be
`NA`.

Please note that the Savage-Dickey density ratio as implemented here
provides only a very basic test of point hypotheses. It is recommended
that you use bridge sampling instead (via
[`bayes_factor`](https://paulbuerkner.com/brms/reference/bayes_factor.brmsfit.md)
which relies on the bridgesampling package). When interpreting Bayes
factors for point hypotheses, make sure that your priors are reasonable
and carefully chosen, as the result will depend heavily on the priors.
In particular, avoid using default priors. Additionally, note that point
hypotheses that scope over more than one parameter (e.g., when testing
equality between two parameters) involve random sampling of the priors
over those parameters (to accommodate the the assumption that priors for
different parameters are independent of each other). This introduces an
element of randomness into such hypothesis tests. Consider repeating the
test to assure results are sufficiently stable, and use the argument
`seed` for reproducibility. Finally, note that, for technical reasons,
we cannot sample from priors of certain parameters classes. Most
notably, these include overall intercept parameters (prior class
`"Intercept"`) as well as group-level coefficients.

For one-sided hypotheses, the `Evid.Ratio` may sometimes be `0` or `Inf`
implying very small or large evidence, respectively, in favor of the
tested hypothesis. For one-sided hypotheses pairs, this basically means
that all posterior draws are on the same side of the value dividing the
two hypotheses. In that sense, instead of `0` or `Inf,` you may rather
read it as `Evid.Ratio` smaller `1 / S` or greater `S`, respectively,
where `S` denotes the number of posterior draws used in the
computations.

The argument `alpha` specifies the size of the credible interval (i.e.,
Bayesian confidence interval). For instance, if we tested a two-sided
hypothesis and set `alpha = 0.05` (5%) an, the credible interval will
contain `1 - alpha = 0.95` (95%) of the posterior values. Hence,
`alpha * 100`% of the posterior values will lie outside of the credible
interval. Although this allows testing of hypotheses in a similar manner
as in the frequentist null-hypothesis testing framework, we strongly
argue against using arbitrary cutoffs (e.g., `p < .05`) to determine the
'existence' of an effect.

## See also

[`brmshypothesis`](https://paulbuerkner.com/brms/reference/brmshypothesis.md)

## Examples

``` r
# \dontrun{
## define priors
prior <- c(set_prior("normal(0,2)", class = "b"),
           set_prior("student_t(10,0,1)", class = "sigma"),
           set_prior("student_t(10,0,1)", class = "sd"))

## fit a linear mixed effects models
fit <- brm(time ~ age + sex + disease + (1 + age|patient),
           data = kidney, family = lognormal(),
           prior = prior, sample_prior = "yes",
           control = list(adapt_delta = 0.95))
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 3.2e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.32 seconds.
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
#> Chain 1:  Elapsed Time: 1.467 seconds (Warm-up)
#> Chain 1:                0.578 seconds (Sampling)
#> Chain 1:                2.045 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 2.3e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.23 seconds.
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
#> Chain 2:  Elapsed Time: 1.488 seconds (Warm-up)
#> Chain 2:                0.579 seconds (Sampling)
#> Chain 2:                2.067 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
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
#> Chain 3:  Elapsed Time: 1.294 seconds (Warm-up)
#> Chain 3:                0.583 seconds (Sampling)
#> Chain 3:                1.877 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 3.1e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.31 seconds.
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
#> Chain 4:  Elapsed Time: 1.454 seconds (Warm-up)
#> Chain 4:                1.076 seconds (Sampling)
#> Chain 4:                2.53 seconds (Total)
#> Chain 4: 
#> Warning: There were 2 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems

## perform two-sided hypothesis testing
(hyp1 <- hypothesis(fit, "sexfemale = age + diseasePKD"))
#> Hypothesis Tests for class b:
#>                 Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio
#> 1 (sexfemale)-(age+... = 0     0.69      0.66    -0.63        2        2.8
#>   Post.Prob Star
#> 1      0.74     
#> ---
#> 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
#> '*': For one-sided hypotheses, the posterior probability exceeds 95%;
#> for two-sided hypotheses, the value tested against lies outside the 95%-CI.
#> Posterior probabilities of point hypotheses assume equal prior probabilities.
plot(hyp1)

hypothesis(fit, "exp(age) - 3 = 0", alpha = 0.01)
#> Hypothesis Tests for class b:
#>         Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob
#> 1 (exp(age)-3) = 0    -2.01      0.01    -2.04    -1.97          0         0
#>   Star
#> 1    *
#> ---
#> 'CI': 98%-CI for one-sided and 99%-CI for two-sided hypotheses.
#> '*': For one-sided hypotheses, the posterior probability exceeds 99%;
#> for two-sided hypotheses, the value tested against lies outside the 99%-CI.
#> Posterior probabilities of point hypotheses assume equal prior probabilities.

## perform one-sided hypothesis testing
hypothesis(fit, "diseasePKD + diseaseGN - 3 < 0")
#> Hypothesis Tests for class b:
#>                 Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio
#> 1 (diseasePKD+disea... < 0    -2.84      0.91    -4.35    -1.32        999
#>   Post.Prob Star
#> 1         1    *
#> ---
#> 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
#> '*': For one-sided hypotheses, the posterior probability exceeds 95%;
#> for two-sided hypotheses, the value tested against lies outside the 95%-CI.
#> Posterior probabilities of point hypotheses assume equal prior probabilities.

hypothesis(fit, "age < Intercept",
           class = "sd", group  = "patient")
#> Hypothesis Tests for class sd_patient:
#>              Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio
#> 1 (age)-(Intercept) < 0    -0.33      0.25    -0.81    -0.03        124
#>   Post.Prob Star
#> 1      0.99    *
#> ---
#> 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
#> '*': For one-sided hypotheses, the posterior probability exceeds 95%;
#> for two-sided hypotheses, the value tested against lies outside the 95%-CI.
#> Posterior probabilities of point hypotheses assume equal prior probabilities.

## test the amount of random intercept variance on all variance
h <- paste("sd_patient__Intercept^2 / (sd_patient__Intercept^2 +",
           "sd_patient__age^2 + sigma^2) = 0")
(hyp2 <- hypothesis(fit, h, class = NULL))
#> Hypothesis Tests for class :
#>                 Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio
#> 1 (sd_patient__Inte... = 0     0.09      0.11        0     0.39       4.23
#>   Post.Prob Star
#> 1      0.81    *
#> ---
#> 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
#> '*': For one-sided hypotheses, the posterior probability exceeds 95%;
#> for two-sided hypotheses, the value tested against lies outside the 95%-CI.
#> Posterior probabilities of point hypotheses assume equal prior probabilities.
plot(hyp2)


## test more than one hypothesis at once
h <- c("diseaseGN = diseaseAN", "2 * diseaseGN - diseasePKD = 0")
(hyp3 <- hypothesis(fit, h))
#> Hypothesis Tests for class b:
#>                 Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio
#> 1 (diseaseGN)-(dise... = 0     0.03      0.46    -0.86     0.93       6.35
#> 2 (2*diseaseGN-dise... = 0    -0.90      0.88    -2.60     0.78       2.87
#>   Post.Prob Star
#> 1      0.86     
#> 2      0.74     
#> ---
#> 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
#> '*': For one-sided hypotheses, the posterior probability exceeds 95%;
#> for two-sided hypotheses, the value tested against lies outside the 95%-CI.
#> Posterior probabilities of point hypotheses assume equal prior probabilities.
plot(hyp3, ignore_prior = TRUE)


## compute hypotheses for all levels of a grouping factor
hypothesis(fit, "age = 0", scope = "coef", group = "patient")
#> Hypothesis Tests for class :
#>    Group Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob
#> 1      1  (age) = 0    -0.01      0.02    -0.04     0.03         NA        NA
#> 2      2  (age) = 0    -0.01      0.02    -0.04     0.02         NA        NA
#> 3      3  (age) = 0    -0.01      0.02    -0.04     0.03         NA        NA
#> 4      4  (age) = 0     0.00      0.02    -0.04     0.03         NA        NA
#> 5      5  (age) = 0    -0.01      0.02    -0.04     0.03         NA        NA
#> 6      6  (age) = 0    -0.01      0.02    -0.04     0.02         NA        NA
#> 7      7  (age) = 0    -0.01      0.02    -0.04     0.02         NA        NA
#> 8      8  (age) = 0     0.00      0.02    -0.03     0.03         NA        NA
#> 9      9  (age) = 0     0.00      0.02    -0.03     0.03         NA        NA
#> 10    10  (age) = 0     0.00      0.02    -0.03     0.03         NA        NA
#> 11    11  (age) = 0    -0.01      0.02    -0.04     0.03         NA        NA
#> 12    12  (age) = 0    -0.01      0.02    -0.04     0.02         NA        NA
#> 13    13  (age) = 0    -0.01      0.02    -0.04     0.03         NA        NA
#> 14    14  (age) = 0     0.00      0.02    -0.04     0.03         NA        NA
#> 15    15  (age) = 0    -0.01      0.02    -0.04     0.03         NA        NA
#> 16    16  (age) = 0    -0.01      0.02    -0.04     0.02         NA        NA
#> 17    17  (age) = 0     0.00      0.02    -0.03     0.03         NA        NA
#> 18    18  (age) = 0     0.00      0.02    -0.04     0.03         NA        NA
#> 19    19  (age) = 0    -0.01      0.02    -0.04     0.03         NA        NA
#> 20    20  (age) = 0    -0.01      0.02    -0.04     0.02         NA        NA
#> 21    21  (age) = 0     0.00      0.02    -0.03     0.04         NA        NA
#> 22    22  (age) = 0    -0.01      0.02    -0.04     0.03         NA        NA
#> 23    23  (age) = 0    -0.01      0.02    -0.04     0.02         NA        NA
#> 24    24  (age) = 0    -0.01      0.02    -0.04     0.03         NA        NA
#> 25    25  (age) = 0    -0.01      0.02    -0.04     0.03         NA        NA
#> 26    26  (age) = 0     0.00      0.02    -0.03     0.03         NA        NA
#> 27    27  (age) = 0    -0.01      0.02    -0.04     0.03         NA        NA
#> 28    28  (age) = 0    -0.01      0.02    -0.04     0.02         NA        NA
#> 29    29  (age) = 0    -0.01      0.02    -0.04     0.02         NA        NA
#> 30    30  (age) = 0    -0.01      0.02    -0.04     0.03         NA        NA
#> 31    31  (age) = 0    -0.01      0.02    -0.04     0.02         NA        NA
#> 32    32  (age) = 0    -0.01      0.02    -0.05     0.02         NA        NA
#> 33    33  (age) = 0    -0.01      0.02    -0.04     0.02         NA        NA
#> 34    34  (age) = 0    -0.01      0.02    -0.04     0.02         NA        NA
#> 35    35  (age) = 0    -0.01      0.02    -0.04     0.03         NA        NA
#> 36    36  (age) = 0    -0.01      0.02    -0.04     0.02         NA        NA
#> 37    37  (age) = 0    -0.01      0.02    -0.05     0.02         NA        NA
#> 38    38  (age) = 0    -0.01      0.02    -0.04     0.02         NA        NA
#>    Star
#> 1      
#> 2      
#> 3      
#> 4      
#> 5      
#> 6      
#> 7      
#> 8      
#> 9      
#> 10     
#> 11     
#> 12     
#> 13     
#> 14     
#> 15     
#> 16     
#> 17     
#> 18     
#> 19     
#> 20     
#> 21     
#> 22     
#> 23     
#> 24     
#> 25     
#> 26     
#> 27     
#> 28     
#> 29     
#> 30     
#> 31     
#> 32     
#> 33     
#> 34     
#> 35     
#> 36     
#> 37     
#> 38     
#> ---
#> 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
#> '*': For one-sided hypotheses, the posterior probability exceeds 95%;
#> for two-sided hypotheses, the value tested against lies outside the 95%-CI.
#> Posterior probabilities of point hypotheses assume equal prior probabilities.

## use the default method
dat <- as.data.frame(fit)
str(dat)
#> 'data.frame':    4000 obs. of  94 variables:
#>  $ b_Intercept                : num  3.74 3.74 2.78 3.81 3.82 ...
#>  $ b_age                      : num  -0.000802 -0.002204 0.011277 -0.007372 -0.021224 ...
#>  $ b_sexfemale                : num  1.097 0.834 1.045 0.953 1.367 ...
#>  $ b_diseaseGN                : num  -0.885 -0.531 -0.696 -0.519 -0.127 ...
#>  $ b_diseaseAN                : num  -0.4118 -1.1497 -0.4139 -0.6296 -0.0936 ...
#>  $ b_diseasePKD               : num  -0.0884 -0.4201 -0.2354 0.1743 1.703 ...
#>  $ sd_patient__Intercept      : num  0.284 0.154 0.244 0.489 0.271 ...
#>  $ sd_patient__age            : num  0.00219 0.00519 0.00306 0.01669 0.00821 ...
#>  $ cor_patient__Intercept__age: num  -0.5404 -0.9144 0.8005 -0.0607 -0.7528 ...
#>  $ sigma                      : num  1.25 1.16 1.22 1.12 1.64 ...
#>  $ Intercept                  : num  4.17 3.73 3.73 3.88 4.02 ...
#>  $ r_patient[1,Intercept]     : num  0.0628 -0.0472 -0.0766 0.5533 -0.3003 ...
#>  $ r_patient[2,Intercept]     : num  0.00539 -0.1831 -0.26451 0.60243 -0.34227 ...
#>  $ r_patient[3,Intercept]     : num  -0.152 -0.107 -0.477 0.468 -0.283 ...
#>  $ r_patient[4,Intercept]     : num  0.175 -0.129 -0.107 0.318 -0.144 ...
#>  $ r_patient[5,Intercept]     : num  -0.169 0.1103 0.1626 0.0785 -0.0649 ...
#>  $ r_patient[6,Intercept]     : num  -0.22328 -0.17357 -0.00207 -0.62118 0.48207 ...
#>  $ r_patient[7,Intercept]     : num  -0.2547 -0.0217 -0.0622 0.6282 -0.1007 ...
#>  $ r_patient[8,Intercept]     : num  0.5927 0.1066 0.0864 -0.5588 0.1316 ...
#>  $ r_patient[9,Intercept]     : num  -0.2666 0.0177 0.0308 -0.5143 0.029 ...
#>  $ r_patient[10,Intercept]    : num  0.01679 0.02087 -0.17534 0.00349 -0.2909 ...
#>  $ r_patient[11,Intercept]    : num  0.0426 -0.3315 0.0204 0.6106 -0.303 ...
#>  $ r_patient[12,Intercept]    : num  0.16 -0.127 -0.379 0.464 -0.123 ...
#>  $ r_patient[13,Intercept]    : num  0.2235 0.0928 0.1322 -0.5379 0.2948 ...
#>  $ r_patient[14,Intercept]    : num  -0.033 0.0779 -0.0752 0.0954 -0.075 ...
#>  $ r_patient[15,Intercept]    : num  0.1641 -0.2129 -0.1444 0.2154 -0.0907 ...
#>  $ r_patient[16,Intercept]    : num  -0.0131 0.2838 -0.1326 0.127 -0.0115 ...
#>  $ r_patient[17,Intercept]    : num  -0.3671 0.0857 0.5296 -0.4918 0.0973 ...
#>  $ r_patient[18,Intercept]    : num  0.08089 0.15734 0.00574 0.37839 -0.33488 ...
#>  $ r_patient[19,Intercept]    : num  -0.0952 -0.1257 0.2462 0.0977 -0.1249 ...
#>  $ r_patient[20,Intercept]    : num  0.1364 0.0761 0.4431 -0.0506 0.0373 ...
#>  $ r_patient[21,Intercept]    : num  0.456 0.252 0.326 -1.306 0.163 ...
#>  $ r_patient[22,Intercept]    : num  0.1638 0.0617 -0.0364 -0.2897 0.1527 ...
#>  $ r_patient[23,Intercept]    : num  0.3498 0.0649 0.00707 -0.88331 0.57557 ...
#>  $ r_patient[24,Intercept]    : num  0.3349 0.1369 -0.0876 -0.2936 0.2247 ...
#>  $ r_patient[25,Intercept]    : num  -0.0803 0.0615 0.1833 0.3415 -0.1028 ...
#>  $ r_patient[26,Intercept]    : num  0.317 -0.135 -0.319 0.276 -0.362 ...
#>  $ r_patient[27,Intercept]    : num  0.169 -0.205 -0.574 0.706 -0.208 ...
#>  $ r_patient[28,Intercept]    : num  -0.0209 -0.0598 0.2585 -0.2529 0.2309 ...
#>  $ r_patient[29,Intercept]    : num  0.1506 0.0639 -0.2208 0.1074 0.0801 ...
#>  $ r_patient[30,Intercept]    : num  0.5416 0.1508 0.0235 -0.3952 0.1948 ...
#>  $ r_patient[31,Intercept]    : num  0.0835 -0.0609 0.1463 0.0341 0.0222 ...
#>  $ r_patient[32,Intercept]    : num  -0.29 0.226 -0.302 -0.141 0.175 ...
#>  $ r_patient[33,Intercept]    : num  0.155 -0.138 -0.196 0.374 -0.184 ...
#>  $ r_patient[34,Intercept]    : num  -0.024 -0.394 -0.245 -0.23 0.188 ...
#>  $ r_patient[35,Intercept]    : num  -0.4489 -0.2108 -0.0961 -0.7732 0.435 ...
#>  $ r_patient[36,Intercept]    : num  0.1217 -0.2292 -0.0508 -0.7978 0.6087 ...
#>  $ r_patient[37,Intercept]    : num  -0.4585 0.2441 0.0657 0.0293 -0.023 ...
#>  $ r_patient[38,Intercept]    : num  -0.4501 0.0639 0.2947 -0.3232 0.3829 ...
#>  $ r_patient[1,age]           : num  -0.003007 0.000236 -0.001798 -0.023301 0.013269 ...
#>  $ r_patient[2,age]           : num  0.000539 0.008098 0.000899 -0.031718 0.016198 ...
#>  $ r_patient[3,age]           : num  0.00141 0.00279 -0.00419 0.00327 0.00522 ...
#>  $ r_patient[4,age]           : num  -4.16e-03 6.25e-03 8.52e-05 1.53e-02 -1.37e-03 ...
#>  $ r_patient[5,age]           : num  0.000479 -0.004057 0.000498 -0.021921 0.008381 ...
#>  $ r_patient[6,age]           : num  -0.001187 0.007545 -0.000998 0.02559 -0.017242 ...
#>  $ r_patient[7,age]           : num  -3.04e-03 -1.20e-03 -4.43e-03 2.75e-04 -8.31e-05 ...
#>  $ r_patient[8,age]           : num  -3.61e-05 -5.35e-03 -7.47e-04 2.81e-02 -8.70e-03 ...
#>  $ r_patient[9,age]           : num  3.60e-03 2.99e-03 -6.52e-04 9.77e-03 -4.09e-05 ...
#>  $ r_patient[10,age]          : num  -0.000329 -0.001759 0.00156 0.007668 0.005307 ...
#>  $ r_patient[11,age]          : num  -0.001285 0.00885 0.000166 0.008395 0.003361 ...
#>  $ r_patient[12,age]          : num  0.000185 0.001239 -0.007 0.012701 -0.003081 ...
#>  $ r_patient[13,age]          : num  -0.002959 -0.004307 -0.000363 -0.001697 -0.005594 ...
#>  $ r_patient[14,age]          : num  0.00317 -0.00139 0.00137 0.02271 -0.0051 ...
#>  $ r_patient[15,age]          : num  -0.000832 0.005951 -0.003799 0.005096 0.000253 ...
#>  $ r_patient[16,age]          : num  -0.00211 -0.0089 -0.00181 -0.01522 0.00451 ...
#>  $ r_patient[17,age]          : num  0.00125 -0.00396 0.00368 0.01703 -0.00428 ...
#>  $ r_patient[18,age]          : num  -0.000629 -0.004414 -0.001227 0.013283 0.005808 ...
#>  $ r_patient[19,age]          : num  0.003435 0.004819 0.000462 0.021763 -0.004557 ...
#>  $ r_patient[20,age]          : num  0.00107 -0.00108 0.00605 -0.02045 0.00639 ...
#>  $ r_patient[21,age]          : num  -0.000551 -0.010716 0.003788 0.018187 -0.00273 ...
#>  $ r_patient[22,age]          : num  0.003314 -0.000575 0.000491 -0.017794 0.003343 ...
#>  $ r_patient[23,age]          : num  -0.00176 -0.00244 -0.00115 0.00264 -0.01457 ...
#>  $ r_patient[24,age]          : num  -0.00221 -0.00464 -0.00312 0.02439 -0.01359 ...
#>  $ r_patient[25,age]          : num  -0.002776 -0.000575 0.000215 0.016771 -0.001138 ...
#>  $ r_patient[26,age]          : num  -0.00173 0.00402 -0.00141 0.00273 0.01017 ...
#>  $ r_patient[27,age]          : num  0.00172 0.00501 -0.00331 0.00407 0.00335 ...
#>  $ r_patient[28,age]          : num  0.00078 0.005116 0.003779 -0.024375 0.000172 ...
#>  $ r_patient[29,age]          : num  -0.00206 -0.0019 -0.00222 -0.02794 0.00578 ...
#>  $ r_patient[30,age]          : num  -0.00288 -0.007 -0.00158 0.01894 -0.01025 ...
#>  $ r_patient[31,age]          : num  0.000162 0.000953 -0.000242 0.004564 -0.001694 ...
#>  $ r_patient[32,age]          : num  0.0035 -0.00793 -0.00277 -0.02342 0.00104 ...
#>  $ r_patient[33,age]          : num  -0.00072 0.00454 -0.0023 0.00485 0.00156 ...
#>  $ r_patient[34,age]          : num  0.000443 0.007381 -0.001151 -0.008265 -0.004437 ...
#>  $ r_patient[35,age]          : num  0.00469 0.00785 -0.00282 -0.00018 -0.01122 ...
#>  $ r_patient[36,age]          : num  -0.00378 0.00832 -0.00193 -0.01035 -0.012 ...
#>  $ r_patient[37,age]          : num  0.000217 -0.009887 -0.000324 -0.00624 0.00048 ...
#>  $ r_patient[38,age]          : num  0.00207 -0.00298 0.00273 0.00326 -0.01002 ...
#>  $ prior_Intercept            : num  2.48 5.04 5.37 8.36 5 ...
#>  $ prior_b                    : num  -0.245 1.239 -0.888 -1.436 -0.452 ...
#>  $ prior_sigma                : num  1.348 0.196 1.062 0.16 0.702 ...
#>  $ prior_sd_patient           : num  0.481 0.4477 3.1253 0.8301 0.0892 ...
#>  $ prior_cor_patient          : num  -0.856 -0.439 0.105 -0.289 0.114 ...
#>  $ lprior                     : num  -12.6 -12.4 -12.4 -12.4 -13.4 ...
#>  $ lp__                       : num  -544 -546 -543 -540 -544 ...
hypothesis(dat, "b_age > 0")
#> Hypothesis Tests for class :
#>    Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob Star
#> 1 (b_age) > 0    -0.01      0.01    -0.03     0.02       0.47      0.32     
#> ---
#> 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
#> '*': For one-sided hypotheses, the posterior probability exceeds 95%;
#> for two-sided hypotheses, the value tested against lies outside the 95%-CI.
#> Posterior probabilities of point hypotheses assume equal prior probabilities.
# }
```
