# Stan Code for Bayesian models

`stancode` is a generic function that can be used to generate Stan code
for Bayesian models. Its original use is within the brms package, but
new methods for use with objects from other packages can be registered
to the same generic.

## Usage

``` r
stancode(object, ...)

make_stancode(formula, ...)
```

## Arguments

- object:

  An object whose class will determine which method to apply. Usually,
  it will be some kind of symbolic description of the model form which
  Stan code should be generated.

- ...:

  Further arguments passed to the specific method.

- formula:

  Synonym of `object` for use in `make_stancode`.

## Value

Usually, a character string containing the generated Stan code. For
pretty printing, we recommend the returned object to be of class
`c("character", "brmsmodel")`.

## Details

See
[`stancode.default`](https://paulbuerkner.com/brms/reference/stancode.default.md)
for the default method applied for brms models. You can view the
available methods by typing: `methods(stancode)` The `make_stancode`
function is an alias of `stancode`.

## See also

[`stancode.default`](https://paulbuerkner.com/brms/reference/stancode.default.md),
[`stancode.brmsfit`](https://paulbuerkner.com/brms/reference/stancode.brmsfit.md)

## Examples

``` r
stancode(rating ~ treat + period + carry + (1|subject),
         data = inhaler, family = "cumulative")
#> // generated with brms 2.23.1
#> functions {
#>   /* cumulative-logit log-PDF for a single response
#>    * Args:
#>    *   y: response category
#>    *   mu: latent mean parameter
#>    *   disc: discrimination parameter
#>    *   thres: ordinal thresholds
#>    * Returns:
#>    *   a scalar to be added to the log posterior
#>    */
#>    real cumulative_logit_lpmf(int y, real mu, real disc, vector thres) {
#>      int nthres = num_elements(thres);
#>      if (y == 1) {
#>        return log_inv_logit(disc * (thres[1] - mu));
#>      } else if (y == nthres + 1) {
#>        return log1m_inv_logit(disc * (thres[nthres] - mu));
#>      } else {
#>        return log_inv_logit_diff(disc * (thres[y] - mu), disc * (thres[y - 1] - mu));
#>      }
#>    }
#>   /* cumulative-logit log-PDF for a single response and merged thresholds
#>    * Args:
#>    *   y: response category
#>    *   mu: latent mean parameter
#>    *   disc: discrimination parameter
#>    *   thres: vector of merged ordinal thresholds
#>    *   j: start and end index for the applid threshold within 'thres'
#>    * Returns:
#>    *   a scalar to be added to the log posterior
#>    */
#>    real cumulative_logit_merged_lpmf(int y, real mu, real disc, vector thres, array[] int j) {
#>      return cumulative_logit_lpmf(y | mu, disc, thres[j[1]:j[2]]);
#>    }
#>   /* ordered-logistic log-PDF for a single response and merged thresholds
#>    * Args:
#>    *   y: response category
#>    *   mu: latent mean parameter
#>    *   thres: vector of merged ordinal thresholds
#>    *   j: start and end index for the applid threshold within 'thres'
#>    * Returns:
#>    *   a scalar to be added to the log posterior
#>    */
#>    real ordered_logistic_merged_lpmf(int y, real mu, vector thres, array[] int j) {
#>      return ordered_logistic_lpmf(y | mu, thres[j[1]:j[2]]);
#>    }
#> 
#> }
#> data {
#>   int<lower=1> N;  // total number of observations
#>   array[N] int Y;  // response variable
#>   int<lower=2> nthres;  // number of thresholds
#>   int<lower=1> K;  // number of population-level effects
#>   matrix[N, K] X;  // population-level design matrix
#>   int<lower=1> Kc;  // number of population-level effects after centering
#>   // data for group-level effects of ID 1
#>   int<lower=1> N_1;  // number of grouping levels
#>   int<lower=1> M_1;  // number of coefficients per level
#>   array[N] int<lower=1> J_1;  // grouping indicator per observation
#>   // group-level predictor values
#>   vector[N] Z_1_1;
#>   int prior_only;  // should the likelihood be ignored?
#> }
#> transformed data {
#>   matrix[N, Kc] Xc;  // centered version of X
#>   vector[Kc] means_X;  // column means of X before centering
#>   for (i in 1:K) {
#>     means_X[i] = mean(X[, i]);
#>     Xc[, i] = X[, i] - means_X[i];
#>   }
#> }
#> parameters {
#>   vector[Kc] b;  // regression coefficients
#>   ordered[nthres] Intercept;  // temporary thresholds for centered predictors
#>   vector<lower=0>[M_1] sd_1;  // group-level standard deviations
#>   array[M_1] vector[N_1] z_1;  // standardized group-level effects
#> }
#> transformed parameters {
#>   real disc = 1;  // discrimination parameters
#>   vector[N_1] r_1_1;  // actual group-level effects
#>   // prior contributions to the log posterior
#>   real lprior = 0;
#>   r_1_1 = (sd_1[1] * (z_1[1]));
#>   lprior += student_t_lpdf(Intercept | 3, 0, 2.5);
#>   lprior += student_t_lpdf(sd_1 | 3, 0, 2.5)
#>     - 1 * student_t_lccdf(0 | 3, 0, 2.5);
#> }
#> model {
#>   // likelihood including constants
#>   if (!prior_only) {
#>     // initialize linear predictor term
#>     vector[N] mu = rep_vector(0.0, N);
#>     mu += Xc * b;
#>     for (n in 1:N) {
#>       // add more terms to the linear predictor
#>       mu[n] += r_1_1[J_1[n]] * Z_1_1[n];
#>     }
#>     for (n in 1:N) {
#>       target += ordered_logistic_lpmf(Y[n] | mu[n], Intercept);
#>     }
#>   }
#>   // priors including constants
#>   target += lprior;
#>   target += std_normal_lpdf(z_1[1]);
#> }
#> generated quantities {
#>   // compute actual thresholds
#>   vector[nthres] b_Intercept = Intercept + dot_product(means_X, b);
#> }
```
