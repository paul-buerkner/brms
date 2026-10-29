# Handle Missing Values with brms

## Introduction

Many real world data sets contain missing values for various reasons.
Generally, we have quite a few options to handle those missing values.
The easiest solution is to remove all rows from the data set, where one
or more variables are missing. However, if values are not missing
completely at random, this will likely lead to bias in our analysis.
Accordingly, we usually want to impute missing values in one way or the
other. Here, we will consider two very general approaches using
**brms**: (1) Impute missing values *before* the model fitting with
multiple imputation, and (2) impute missing values on the fly *during*
model fitting[^1]. As a simple example, we will use the `nhanes` data
set, which contains information on participants’ `age`, `bmi` (body mass
index), `hyp` (hypertensive), and `chl` (total serum cholesterol). For
the purpose of the present vignette, we are primarily interested in
predicting `bmi` by `age` and `chl`.

``` r

data("nhanes", package = "mice")
head(nhanes)
```

      age  bmi hyp chl
    1   1   NA  NA  NA
    2   2 22.7   1 187
    3   1   NA   1 187
    4   3   NA  NA  NA
    5   1 20.4   1 113
    6   3   NA  NA 184

## Imputation before model fitting

There are many approaches allowing us to impute missing data before the
actual model fitting takes place. From a statistical perspective,
multiple imputation is one of the best solutions. Each missing value is
not imputed once but `m` times leading to a total of `m` fully imputed
data sets. The model can then be fitted to each of those data sets
separately and results are pooled across models, afterwards. One widely
applied package for multiple imputation is **mice** (Buuren &
Groothuis-Oudshoorn, 2010) and we will use it in the following in
combination with **brms**. Here, we apply the default settings of
**mice**, which means that all variables will be used to impute missing
values in all other variables and imputation functions automatically
chosen based on the variables’ characteristics.

``` r

library(mice)
m <- 5
imp <- mice(nhanes, m = m, print = FALSE)
```

Now, we have `m = 5` imputed data sets stored within the `imp` object.
In practice, we will likely need more than `5` of those to accurately
account for the uncertainty induced by the missingness, perhaps even in
the area of `100` imputed data sets (Zhou & Reiter, 2010). Of course,
this increases the computational burden by a lot and so we stick to
`m = 5` for the purpose of this vignette. Regardless of the value of
`m`, we can either extract those data sets and then pass them to the
actual model fitting function as a list of data frames, or pass `imp`
directly. The latter works because **brms** offers special support for
data imputed by **mice**. We will go with the latter approach, since it
is less typing. Fitting our model of interest with **brms** to the
multiple imputed data sets is straightforward.

``` r

fit_imp1 <- brm_multiple(bmi ~ age*chl, data = imp, chains = 2)
```

The returned fitted model is an ordinary `brmsfit` object containing the
posterior draws of all `m` submodels. While pooling across models is not
necessarily straightforward in classical statistics, it is trivial in a
Bayesian framework. Here, pooling results of multiple imputed data sets
is simply achieved by combining the posterior draws of the submodels.
Accordingly, all post-processing methods can be used out of the box
without having to worry about pooling at all.

``` r

summary(fit_imp1)
```

     Family: gaussian 
      Links: mu = identity 
    Formula: bmi ~ age * chl 
       Data: imp (Number of observations: 25) 
      Draws: 10 chains, each with iter = 2000; warmup = 1000; thin = 1;
             total post-warmup draws = 10000

    Regression Coefficients:
              Estimate Est.Error l-95% CI u-95% CI
    Intercept    14.71      8.19    -0.80    30.41
    age           1.75      4.86    -7.37    12.25
    chl           0.09      0.04     0.00     0.17
    age:chl      -0.02      0.02    -0.07     0.02

    Further Distributional Parameters:
          Estimate Est.Error l-95% CI u-95% CI
    sigma     3.33      0.76     2.13     5.06

    Draws were sampled using sampling(NUTS). Overall Rhat and ESS estimates
    are not informative for brm_multiple models and are hence not displayed.
    Please see ?brm_multiple for how to assess convergence in case of such models.

In the summary output, we notice that some `Rhat` values are higher than
\\1.1\\ indicating possible convergence problems. For models based on
multiple imputed data sets, this is often a **false positive**: Chains
of different submodels may not overlay each other exactly, since there
were fitted to different data. We can see the chains on the right-hand
side of

``` r

plot(fit_imp1, variable = "^b", regex = TRUE)
```

![](brms_missings_files/figure-html/unnamed-chunk-5-1.png)

Such non-overlaying chains imply high `Rhat` values without there
actually being any convergence issue. Accordingly, we have to
investigate the convergence of the submodels separately, which we can do
for example via:

``` r

library(posterior)
draws <- as_draws_array(fit_imp1)
# every dataset has nc = 2 chains in this example
nc <- nchains(fit_imp1) / m
draws_per_dat <- lapply(1:m, 
  \(i) subset_draws(draws, chain = ((i-1)*nc+1):(i*nc))
)
lapply(draws_per_dat, summarise_draws, default_convergence_measures())
```

    [[1]]
     [38;5;246m# A tibble: 8 × 4 [39m
      variable     rhat ess_bulk ess_tail
       [3m [38;5;246m<chr> [39m [23m        [3m [38;5;246m<dbl> [39m [23m     [3m [38;5;246m<dbl> [39m [23m     [3m [38;5;246m<dbl> [39m [23m
     [38;5;250m1 [39m b_Intercept  1.00     575.     555.
     [38;5;250m2 [39m b_age        1.00     586.     576.
     [38;5;250m3 [39m b_chl        1.00     600.     637.
     [38;5;250m4 [39m b_age:chl    1.01     537.     490.
     [38;5;250m5 [39m sigma        1.00     918.     [4m1 [24m115.
     [38;5;250m6 [39m Intercept    1.00     [4m1 [24m260.     [4m1 [24m160.
     [38;5;250m7 [39m lprior       1.00     856.     [4m1 [24m064.
     [38;5;250m8 [39m lp__         1.00     602.     813.

    [[2]]
     [38;5;246m# A tibble: 8 × 4 [39m
      variable     rhat ess_bulk ess_tail
       [3m [38;5;246m<chr> [39m [23m        [3m [38;5;246m<dbl> [39m [23m     [3m [38;5;246m<dbl> [39m [23m     [3m [38;5;246m<dbl> [39m [23m
     [38;5;250m1 [39m b_Intercept  1.01     412.     539.
     [38;5;250m2 [39m b_age        1.01     413.     418.
     [38;5;250m3 [39m b_chl        1.01     451.     549.
     [38;5;250m4 [39m b_age:chl    1.01     396.     433.
     [38;5;250m5 [39m sigma        1.00     950.     966.
     [38;5;250m6 [39m Intercept    1.00     [4m1 [24m117.     [4m1 [24m130.
     [38;5;250m7 [39m lprior       1.00     954.     902.
     [38;5;250m8 [39m lp__         1.00     618.     708.

    [[3]]
     [38;5;246m# A tibble: 8 × 4 [39m
      variable     rhat ess_bulk ess_tail
       [3m [38;5;246m<chr> [39m [23m        [3m [38;5;246m<dbl> [39m [23m     [3m [38;5;246m<dbl> [39m [23m     [3m [38;5;246m<dbl> [39m [23m
     [38;5;250m1 [39m b_Intercept  1.00     890.     [4m1 [24m070.
     [38;5;250m2 [39m b_age        1.00     827.     937.
     [38;5;250m3 [39m b_chl        1.00     973.     [4m1 [24m103.
     [38;5;250m4 [39m b_age:chl    1.00     818.     949.
     [38;5;250m5 [39m sigma        1.01     [4m1 [24m139.     [4m1 [24m234.
     [38;5;250m6 [39m Intercept    1.00     [4m1 [24m246.     [4m1 [24m034.
     [38;5;250m7 [39m lprior       1.00     892.     [4m1 [24m091.
     [38;5;250m8 [39m lp__         1.00     594.     [4m1 [24m062.

    [[4]]
     [38;5;246m# A tibble: 8 × 4 [39m
      variable     rhat ess_bulk ess_tail
       [3m [38;5;246m<chr> [39m [23m        [3m [38;5;246m<dbl> [39m [23m     [3m [38;5;246m<dbl> [39m [23m     [3m [38;5;246m<dbl> [39m [23m
     [38;5;250m1 [39m b_Intercept  1.01     691.     972.
     [38;5;250m2 [39m b_age        1.01     671.     985.
     [38;5;250m3 [39m b_chl        1.01     733.     973.
     [38;5;250m4 [39m b_age:chl    1.01     643.     890.
     [38;5;250m5 [39m sigma        1.00     957.     931.
     [38;5;250m6 [39m Intercept    1.00     [4m1 [24m238.     [4m1 [24m061.
     [38;5;250m7 [39m lprior       1.00     879.     931.
     [38;5;250m8 [39m lp__         1.00     705.     [4m1 [24m043.

    [[5]]
     [38;5;246m# A tibble: 8 × 4 [39m
      variable     rhat ess_bulk ess_tail
       [3m [38;5;246m<chr> [39m [23m        [3m [38;5;246m<dbl> [39m [23m     [3m [38;5;246m<dbl> [39m [23m     [3m [38;5;246m<dbl> [39m [23m
     [38;5;250m1 [39m b_Intercept  1.00     487.     836.
     [38;5;250m2 [39m b_age        1.01     462.     805.
     [38;5;250m3 [39m b_chl        1.00     514.     756.
     [38;5;250m4 [39m b_age:chl    1.01     447.     763.
     [38;5;250m5 [39m sigma        1.00     865.     836.
     [38;5;250m6 [39m Intercept    1.00     [4m1 [24m028.     905.
     [38;5;250m7 [39m lprior       1.00     809.     789.
     [38;5;250m8 [39m lp__         1.00     677.     928.

The convergence of each of the submodels looks good. Accordingly, we can
proceed with further post-processing and interpretation of the results.
For instance, we could investigate the combined effect of `age` and
`chl`.

``` r

conditional_effects(fit_imp1, "age:chl")
```

![](brms_missings_files/figure-html/unnamed-chunk-7-1.png)

To summarize, the advantages of multiple imputation are obvious: One can
apply it to all kinds of models, since model fitting functions do not
need to know that the data sets were imputed, beforehand. Also, we do
not need to worry about pooling across submodels when using fully
Bayesian methods. The only drawback is the amount of time required for
model fitting. Estimating Bayesian models is already quite slow with
just a single data set and it only gets worse when working with multiple
imputation.

### Compatibility with other multiple imputation packages

**brms** offers built-in support for **mice** mainly because I use the
latter in some of my own research projects. Nevertheless, `brm_multiple`
supports all kinds of multiple imputation packages as it also accepts a
*list* of data frames as input for its `data` argument. Thus, you just
need to extract the imputed data frames in the form of a list, which can
then be passed to `brm_multiple`. Most multiple imputation packages have
some built-in functionality for this task. When using the **mi**
package, for instance, you simply need to call the `mi::complete`
function to get the desired output.

## Imputation during model fitting

Imputation during model fitting is generally thought to be more complex
than imputation before model fitting, because one has to take care of
everything within one step. This remains true when imputing missing
values with **brms**, but possibly to a somewhat smaller degree.
Consider again the `nhanes` data with the goal to predict `bmi` by
`age`, and `chl`. Since `age` contains no missing values, we only have
to take special care of `bmi` and `chl`. We need to tell the model two
things. (1) Which variables contain missing values and how they should
be predicted, as well as (2) which of these imputed variables should be
used as predictors. In **brms** we can do this as follows:

``` r

bform <- bf(bmi | mi() ~ age * mi(chl)) +
  bf(chl | mi() ~ age) + set_rescor(FALSE)
fit_imp2 <- brm(bform, data = nhanes)
```

The model has become multivariate, as we no longer only predict `bmi`
but also `chl` (see
[`vignette("brms_multivariate")`](https://paulbuerkner.com/brms/articles/brms_multivariate.md)
for details about the multivariate syntax of **brms**). We ensure that
missings in both variables will be modeled rather than excluded by
adding `| mi()` on the left-hand side of the formulas[^2]. We write
`mi(chl)` on the right-hand side of the formula for `bmi` to ensure that
the estimated missing values of `chl` will be used in the prediction of
`bmi`. The summary is a bit more cluttered as we get coefficients for
both response variables, but apart from that we can interpret
coefficients in the usual way.

``` r

summary(fit_imp2)
```

     Family: MV(gaussian, gaussian) 
      Links: mu = identity
             mu = identity 
    Formula: bmi | mi() ~ age * mi(chl) 
             chl | mi() ~ age 
       Data: nhanes (Number of observations: 25) 
      Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
             total post-warmup draws = 4000

    Regression Coefficients:
                  Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    bmi_Intercept    13.94      8.95    -3.19    31.91 1.00     1593     1589
    chl_Intercept   136.79     22.70    90.85   180.71 1.00     3662     3195
    bmi_age          -4.86      5.80   -16.81     5.78 1.00     1326     1584
    chl_age          31.60     11.55     8.88    54.50 1.00     3658     3111
    bmi_michl         0.12      0.05     0.03     0.21 1.00     1662     1656
    bmi_michl:age    -0.01      0.03    -0.05     0.05 1.00     1355     1434

    Further Distributional Parameters:
              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    sigma_bmi     3.27      0.72     2.16     4.97 1.00     1469     2358
    sigma_chl    36.54      6.29    26.31    50.85 1.00     2414     3163

    Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    and Tail_ESS are effective sample size measures, and Rhat is the potential
    scale reduction factor on split chains (at convergence, Rhat = 1).

``` r

conditional_effects(fit_imp2, "age:chl", resp = "bmi")
```

![](brms_missings_files/figure-html/unnamed-chunk-9-1.png)

The results look pretty similar to those obtained from multiple
imputation, but be aware that this may not be generally the case. In
multiple imputation, the default is to impute all variables based on all
other variables, while in the ‘one-step’ approach, we have to explicitly
specify the variables used in the imputation. Thus, arguably, multiple
imputation is easier to apply. An obvious advantage of the ‘one-step’
approach is that the model needs to be fitted only once instead of `m`
times. Also, within the **brms** framework, we can use multilevel
structure and complex non-linear relationships for the imputation of
missing values, which is not achieved as easily in standard multiple
imputation software. On the downside, it is currently not possible to
impute discrete variables, because **Stan** (the engine behind **brms**)
does not allow estimating discrete parameters.

### Combining measurement error and missing values

Missing value terms in **brms** cannot only handle missing values but
also measurement error, or arbitrary combinations of the two. In fact,
we can think of a missing value as a value with infinite measurement
error. Thus, `mi` terms are a natural (and somewhat more verbose)
generalization of the now soft deprecated `me` terms. Suppose we had
measured the variable `chl` with some known error:

``` r

nhanes$se <- rexp(nrow(nhanes), 2)
```

Then we can go ahead an include this information into the model as
follows:

``` r

bform <- bf(bmi | mi() ~ age * mi(chl)) +
  bf(chl | mi(se) ~ age) + set_rescor(FALSE)
fit_imp3 <- brm(bform, data = nhanes)
```

Summarizing and post-processing the model continues to work as usual.

## References

Buuren, S. V. & Groothuis-Oudshoorn, K. (2010). mice: Multivariate
imputation by chained equations in R. *Journal of Statistical Software*,
1-68. doi.org/10.18637/jss.v045.i03

Zhou, X. & Reiter, J. P. (2010). A Note on Bayesian Inference After
Multiple Imputation. *The American Statistician*, 64(2), 159-163.
doi.org/10.1198/tast.2010.09109

[^1]: Actually, there is a third approach that only applies to missings
    in response variables. If we want to impute missing responses, we
    just fit the model using the observed responses and than impute the
    missings *after* fitting the model by means of posterior prediction.
    That is, we supply the predictor values corresponding to missing
    responses to the `predict` method.

[^2]: We don’t really need this for `bmi`, since `bmi` is not used as a
    predictor for another variable. Accordingly, we could also – and
    equivalently – impute missing values of `bmi` *after* model fitting
    by means of posterior prediction.
