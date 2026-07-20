# Constant priors in brms

Function used to set up constant priors in brms. The function does not
evaluate its arguments – it exists purely to help set up the model.

## Usage

``` r
constant(const, broadcast = TRUE)
```

## Arguments

- const:

  Numeric value, vector, matrix of values to which the parameters should
  be fixed to. Can also be a valid Stan variable in the model.

- broadcast:

  Should `const` be automatically broadcasted to the correct size of the
  parameter? Defaults to `TRUE`. If you supply vectors or matrices in
  `const` or vector/matrix valued Stan variables, you need to set
  `broadcast` to `TRUE` (see Examples).

## Value

A named list with elements `const` and `broadcast`.

## See also

[`set_prior`](https://paulbuerkner.com/brms/reference/set_prior.md)

## Examples

``` r
stancode(count ~ Base + Age, data = epilepsy,
         prior = prior(constant(1), class = "b"))
#> // generated with brms 2.23.1
#> functions {
#> }
#> data {
#>   int<lower=1> N;  // total number of observations
#>   vector[N] Y;  // response variable
#>   int<lower=1> K;  // number of population-level effects
#>   matrix[N, K] X;  // population-level design matrix
#>   int<lower=1> Kc;  // number of population-level effects after centering
#>   int prior_only;  // should the likelihood be ignored?
#> }
#> transformed data {
#>   matrix[N, Kc] Xc;  // centered version of X without an intercept
#>   vector[Kc] means_X;  // column means of X before centering
#>   for (i in 2:K) {
#>     means_X[i - 1] = mean(X[, i]);
#>     Xc[, i - 1] = X[, i] - means_X[i - 1];
#>   }
#> }
#> parameters {
#>   real Intercept;  // temporary intercept for centered predictors
#>   real<lower=0> sigma;  // dispersion parameter
#> }
#> transformed parameters {
#>   vector[Kc] b;  // regression coefficients
#>   // prior contributions to the log posterior
#>   real lprior = 0;
#>   b = rep_vector(1, rows(b));
#>   lprior += student_t_lpdf(Intercept | 3, 4, 4.4);
#>   lprior += student_t_lpdf(sigma | 3, 0, 4.4)
#>     - 1 * student_t_lccdf(0 | 3, 0, 4.4);
#> }
#> model {
#>   // likelihood including constants
#>   if (!prior_only) {
#>     target += normal_id_glm_lpdf(Y | Xc, Intercept, b, sigma);
#>   }
#>   // priors including constants
#>   target += lprior;
#> }
#> generated quantities {
#>   // actual population-level intercept
#>   real b_Intercept = Intercept - dot_product(means_X, b);
#> }

# will fail parsing because brms will try to broadcast a vector into a vector
stancode(count ~ Base + Age, data = epilepsy,
         prior = prior(constant(alpha), class = "b"),
         stanvars = stanvar(c(1, 0), name = "alpha"))
#> // generated with brms 2.23.1
#> functions {
#> }
#> data {
#>   int<lower=1> N;  // total number of observations
#>   vector[N] Y;  // response variable
#>   int<lower=1> K;  // number of population-level effects
#>   matrix[N, K] X;  // population-level design matrix
#>   int<lower=1> Kc;  // number of population-level effects after centering
#>   int prior_only;  // should the likelihood be ignored?
#>   vector[2] alpha;
#> }
#> transformed data {
#>   matrix[N, Kc] Xc;  // centered version of X without an intercept
#>   vector[Kc] means_X;  // column means of X before centering
#>   for (i in 2:K) {
#>     means_X[i - 1] = mean(X[, i]);
#>     Xc[, i - 1] = X[, i] - means_X[i - 1];
#>   }
#> }
#> parameters {
#>   real Intercept;  // temporary intercept for centered predictors
#>   real<lower=0> sigma;  // dispersion parameter
#> }
#> transformed parameters {
#>   vector[Kc] b;  // regression coefficients
#>   // prior contributions to the log posterior
#>   real lprior = 0;
#>   b = rep_vector(alpha, rows(b));
#>   lprior += student_t_lpdf(Intercept | 3, 4, 4.4);
#>   lprior += student_t_lpdf(sigma | 3, 0, 4.4)
#>     - 1 * student_t_lccdf(0 | 3, 0, 4.4);
#> }
#> model {
#>   // likelihood including constants
#>   if (!prior_only) {
#>     target += normal_id_glm_lpdf(Y | Xc, Intercept, b, sigma);
#>   }
#>   // priors including constants
#>   target += lprior;
#> }
#> generated quantities {
#>   // actual population-level intercept
#>   real b_Intercept = Intercept - dot_product(means_X, b);
#> }

stancode(count ~ Base + Age, data = epilepsy,
         prior = prior(constant(alpha, broadcast = FALSE), class = "b"),
         stanvars = stanvar(c(1, 0), name = "alpha"))
#> // generated with brms 2.23.1
#> functions {
#> }
#> data {
#>   int<lower=1> N;  // total number of observations
#>   vector[N] Y;  // response variable
#>   int<lower=1> K;  // number of population-level effects
#>   matrix[N, K] X;  // population-level design matrix
#>   int<lower=1> Kc;  // number of population-level effects after centering
#>   int prior_only;  // should the likelihood be ignored?
#>   vector[2] alpha;
#> }
#> transformed data {
#>   matrix[N, Kc] Xc;  // centered version of X without an intercept
#>   vector[Kc] means_X;  // column means of X before centering
#>   for (i in 2:K) {
#>     means_X[i - 1] = mean(X[, i]);
#>     Xc[, i - 1] = X[, i] - means_X[i - 1];
#>   }
#> }
#> parameters {
#>   real Intercept;  // temporary intercept for centered predictors
#>   real<lower=0> sigma;  // dispersion parameter
#> }
#> transformed parameters {
#>   vector[Kc] b;  // regression coefficients
#>   // prior contributions to the log posterior
#>   real lprior = 0;
#>   b = alpha;
#>   lprior += student_t_lpdf(Intercept | 3, 4, 4.4);
#>   lprior += student_t_lpdf(sigma | 3, 0, 4.4)
#>     - 1 * student_t_lccdf(0 | 3, 0, 4.4);
#> }
#> model {
#>   // likelihood including constants
#>   if (!prior_only) {
#>     target += normal_id_glm_lpdf(Y | Xc, Intercept, b, sigma);
#>   }
#>   // priors including constants
#>   target += lprior;
#> }
#> generated quantities {
#>   // actual population-level intercept
#>   real b_Intercept = Intercept - dot_product(means_X, b);
#> }
```
