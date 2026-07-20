# User-defined variables passed to Stan

Prepare user-defined variables to be passed to one of Stan's program
blocks. This is primarily useful for defining more complex priors, for
refitting models without recompilation despite changing priors, or for
defining custom Stan functions.

## Usage

``` r
stanvar(
  x = NULL,
  name = NULL,
  scode = NULL,
  block = "data",
  position = "start",
  pll_args = NULL
)
```

## Arguments

- x:

  An R object containing data to be passed to Stan. Only required if
  `block = 'data'` and ignored otherwise.

- name:

  Optional character string providing the desired variable name of the
  object in `x`. If `NULL` (the default) the variable name is directly
  inferred from `x`.

- scode:

  Line of Stan code to define the variable in Stan language. If
  `block = 'data'`, the Stan code is inferred based on the class of `x`
  by default.

- block:

  Name of one of Stan's program blocks in which the variable should be
  defined. Can be `'data'`, `'tdata'` (transformed data),
  `'parameters'`, `'tparameters'` (transformed parameters), `'model'`,
  `'likelihood'` (part of the model block where the likelihood is
  given), `'genquant'` (generated quantities) or `'functions'`.

- position:

  Name of the position within the block where the Stan code should be
  placed. Currently allowed are `'start'` (the default) and `'end'` of
  the block.

- pll_args:

  Optional Stan code to be put into the header of `partial_log_lik`
  functions. This ensures that the variables specified in `scode` can be
  used in the likelihood even when within-chain parallelization is
  activated via
  [`threading`](https://paulbuerkner.com/brms/reference/threading.md).

## Value

An object of class `stanvars`.

## Details

The `stanvar` function is not vectorized. Instead, multiple `stanvars`
objects can be added together via `+` (see Examples).

Special attention is necessary when using `stanvars` to inject code into
the `'likelihood'` block while having
[`threading`](https://paulbuerkner.com/brms/reference/threading.md)
activated. In this case, your custom Stan code may need adjustments to
ensure correct observation indexing. Please investigate the generated
Stan code via
[`stancode`](https://paulbuerkner.com/brms/reference/stancode.default.md)
to see which adjustments are necessary in your case.

## Examples

``` r
bprior <- prior(normal(mean_intercept, 10), class = "Intercept")
stanvars <- stanvar(5, name = "mean_intercept")
stancode(count ~ Trt, epilepsy, prior = bprior,
         stanvars = stanvars)
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
#>   real mean_intercept;
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
#>   vector[Kc] b;  // regression coefficients
#>   real Intercept;  // temporary intercept for centered predictors
#>   real<lower=0> sigma;  // dispersion parameter
#> }
#> transformed parameters {
#>   // prior contributions to the log posterior
#>   real lprior = 0;
#>   lprior += normal_lpdf(Intercept | mean_intercept, 10);
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

# define a multi-normal prior with known covariance matrix
bprior <- prior(multi_normal(M, V), class = "b")
stanvars <- stanvar(rep(0, 2), "M", scode = "  vector[K] M;") +
  stanvar(diag(2), "V", scode = "  matrix[K, K] V;")
stancode(count ~ Trt + zBase, epilepsy,
         prior = bprior, stanvars = stanvars)
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
#>     vector[K] M;
#>     matrix[K, K] V;
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
#>   vector[Kc] b;  // regression coefficients
#>   real Intercept;  // temporary intercept for centered predictors
#>   real<lower=0> sigma;  // dispersion parameter
#> }
#> transformed parameters {
#>   // prior contributions to the log posterior
#>   real lprior = 0;
#>   lprior += multi_normal_lpdf(b | M, V);
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

# define a hierachical prior on the regression coefficients
bprior <- set_prior("normal(0, tau)", class = "b") +
  set_prior("target += normal_lpdf(tau | 0, 10)", check = FALSE)
stanvars <- stanvar(scode = "real<lower=0> tau;",
                    block = "parameters")
stancode(count ~ Trt + zBase, epilepsy,
         prior = bprior, stanvars = stanvars)
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
#>   vector[Kc] b;  // regression coefficients
#>   real Intercept;  // temporary intercept for centered predictors
#>   real<lower=0> sigma;  // dispersion parameter
#>   real<lower=0> tau;
#> }
#> transformed parameters {
#>   // prior contributions to the log posterior
#>   real lprior = 0;
#>   lprior += normal_lpdf(b | 0, tau);
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
#>   target += normal_lpdf(tau | 0, 10);
#> }
#> generated quantities {
#>   // actual population-level intercept
#>   real b_Intercept = Intercept - dot_product(means_X, b);
#> }

# ensure that 'tau' is passed to the likelihood of a threaded model
# not necessary for this example but may be necessary in other cases
stanvars <- stanvar(scode = "real<lower=0> tau;",
                    block = "parameters", pll_args = "real tau")
stancode(count ~ Trt + zBase, epilepsy,
         stanvars = stanvars, threads = threading(2))
#> // generated with brms 2.23.1
#> functions {
#>   /* integer sequence of values
#>    * Args:
#>    *   start: starting integer
#>    *   end: ending integer
#>    * Returns:
#>    *   an integer sequence from start to end
#>    */
#>   array[] int sequence(int start, int end) {
#>     array[end - start + 1] int seq;
#>     for (n in 1:num_elements(seq)) {
#>       seq[n] = n + start - 1;
#>     }
#>     return seq;
#>   }
#>   // compute partial sums of the log-likelihood
#>   real partial_log_lik_lpmf(array[] int seq, int start, int end, data vector Y, data matrix Xc, vector b, real Intercept, real sigma, real tau) {
#>     real ptarget = 0;
#>     int N = end - start + 1;
#>     ptarget += normal_id_glm_lpdf(Y[start:end] | Xc[start:end], Intercept, b, sigma);
#>     return ptarget;
#>   }
#> }
#> data {
#>   int<lower=1> N;  // total number of observations
#>   vector[N] Y;  // response variable
#>   int<lower=1> K;  // number of population-level effects
#>   matrix[N, K] X;  // population-level design matrix
#>   int<lower=1> Kc;  // number of population-level effects after centering
#>   int grainsize;  // grainsize for threading
#>   int prior_only;  // should the likelihood be ignored?
#> }
#> transformed data {
#>   matrix[N, Kc] Xc;  // centered version of X without an intercept
#>   vector[Kc] means_X;  // column means of X before centering
#>   array[N] int seq = sequence(1, N);
#>   for (i in 2:K) {
#>     means_X[i - 1] = mean(X[, i]);
#>     Xc[, i - 1] = X[, i] - means_X[i - 1];
#>   }
#> }
#> parameters {
#>   vector[Kc] b;  // regression coefficients
#>   real Intercept;  // temporary intercept for centered predictors
#>   real<lower=0> sigma;  // dispersion parameter
#>   real<lower=0> tau;
#> }
#> transformed parameters {
#>   // prior contributions to the log posterior
#>   real lprior = 0;
#>   lprior += student_t_lpdf(Intercept | 3, 4, 4.4);
#>   lprior += student_t_lpdf(sigma | 3, 0, 4.4)
#>     - 1 * student_t_lccdf(0 | 3, 0, 4.4);
#> }
#> model {
#>   // likelihood including constants
#>   if (!prior_only) {
#>     target += reduce_sum(partial_log_lik_lpmf, seq, grainsize, Y, Xc, b, Intercept, sigma, tau);
#>   }
#>   // priors including constants
#>   target += lprior;
#> }
#> generated quantities {
#>   // actual population-level intercept
#>   real b_Intercept = Intercept - dot_product(means_X, b);
#> }
```
