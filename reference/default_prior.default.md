# Default Priors for brms Models

Get information on all parameters (and parameter classes) for which
priors may be specified including default priors.

## Usage

``` r
# Default S3 method
default_prior(
  object,
  data,
  family = gaussian(),
  autocor = NULL,
  data2 = NULL,
  knots = NULL,
  drop_unused_levels = TRUE,
  sparse = NULL,
  ...
)
```

## Arguments

- object:

  An object of class [`formula`](https://rdrr.io/r/stats/formula.html),
  [`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md),
  or
  [`mvbrmsformula`](https://paulbuerkner.com/brms/reference/mvbrmsformula.md)
  (or one that can be coerced to that classes): A symbolic description
  of the model to be fitted. The details of model specification are
  explained in
  [`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md).

- data:

  An object of class `data.frame` (or one that can be coerced to that
  class) containing data of all variables used in the model.

- family:

  A description of the response distribution and link function to be
  used in the model. This can be a family function, a call to a family
  function or a character string naming the family. Every family
  function has a `link` argument allowing to specify the link function
  to be applied on the response variable. If not specified, default
  links are used. For details of supported families see
  [`brmsfamily`](https://paulbuerkner.com/brms/reference/brmsfamily.md).
  By default, a linear `gaussian` model is applied. In multivariate
  models, `family` might also be a list of families.

- autocor:

  (Deprecated) An optional
  [`cor_brms`](https://paulbuerkner.com/brms/reference/cor_brms.md)
  object describing the correlation structure within the response
  variable (i.e., the 'autocorrelation'). See the documentation of
  [`cor_brms`](https://paulbuerkner.com/brms/reference/cor_brms.md) for
  a description of the available correlation structures. Defaults to
  `NULL`, corresponding to no correlations. In multivariate models,
  `autocor` might also be a list of autocorrelation structures. It is
  now recommend to specify autocorrelation terms directly within
  `formula`. See
  [`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md)
  for more details.

- data2:

  A named `list` of objects containing data, which cannot be passed via
  argument `data`. Required for some objects used in autocorrelation
  structures to specify dependency structures as well as for
  within-group covariance matrices.

- knots:

  Optional list containing user specified knot values to be used for
  basis construction of smoothing terms. See
  [`gamm`](https://rdrr.io/pkg/mgcv/man/gamm.html) for more details.

- drop_unused_levels:

  Should unused factors levels in the data be dropped? Defaults to
  `TRUE`.

- sparse:

  (Deprecated) Logical; indicates whether the population-level design
  matrices should be treated as sparse (defaults to `FALSE`). For design
  matrices with many zeros, this can considerably reduce required
  memory. Sampling speed is currently not improved or even slightly
  decreased. It is now recommended to use the `sparse` argument of
  [`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md)
  and related functions.

- ...:

  Other arguments for internal usage only.

## Value

A `brmsprior` object. That is, a data.frame with specific columns
including `prior`, `class`, `coef`, and `group` and several rows, each
providing information on a parameter (or parameter class) on which
priors can be specified. The prior column is empty except for internal
default priors.

## See also

[`default_prior`](https://paulbuerkner.com/brms/reference/default_prior.md),
[`set_prior`](https://paulbuerkner.com/brms/reference/set_prior.md)

## Examples

``` r
# get all parameters and parameters classes to define priors on
(prior <- default_prior(count ~ zAge + zBase * Trt + (1|patient) + (1|obs),
                        data = epilepsy, family = poisson()))
#>                   prior     class       coef   group resp dpar nlpar lb ub tag
#>  student_t(3, 1.4, 2.5) Intercept                                             
#>                  (flat)         b                                             
#>                  (flat)         b       Trt1                                  
#>                  (flat)         b       zAge                                  
#>                  (flat)         b      zBase                                  
#>                  (flat)         b zBase:Trt1                                  
#>    student_t(3, 0, 2.5)        sd                                     0       
#>    student_t(3, 0, 2.5)        sd                obs                  0       
#>    student_t(3, 0, 2.5)        sd  Intercept     obs                  0       
#>    student_t(3, 0, 2.5)        sd            patient                  0       
#>    student_t(3, 0, 2.5)        sd  Intercept patient                  0       
#>        source
#>       default
#>       default
#>  (vectorized)
#>  (vectorized)
#>  (vectorized)
#>  (vectorized)
#>       default
#>  (vectorized)
#>  (vectorized)
#>  (vectorized)
#>  (vectorized)

# define a prior on all population-level effects a once
prior$prior[1] <- "normal(0,10)"

# define a specific prior on the population-level effect of Trt
prior$prior[5] <- "student_t(10, 0, 5)"

# verify that the priors indeed found their way into Stan's model code
stancode(count ~ zAge + zBase * Trt + (1|patient) + (1|obs),
         data = epilepsy, family = poisson(),
         prior = prior)
#> // generated with brms 2.23.1
#> functions {
#> }
#> data {
#>   int<lower=1> N;  // total number of observations
#>   array[N] int Y;  // response variable
#>   int<lower=1> K;  // number of population-level effects
#>   matrix[N, K] X;  // population-level design matrix
#>   int<lower=1> Kc;  // number of population-level effects after centering
#>   // data for group-level effects of ID 1
#>   int<lower=1> N_1;  // number of grouping levels
#>   int<lower=1> M_1;  // number of coefficients per level
#>   array[N] int<lower=1> J_1;  // grouping indicator per observation
#>   // group-level predictor values
#>   vector[N] Z_1_1;
#>   // data for group-level effects of ID 2
#>   int<lower=1> N_2;  // number of grouping levels
#>   int<lower=1> M_2;  // number of coefficients per level
#>   array[N] int<lower=1> J_2;  // grouping indicator per observation
#>   // group-level predictor values
#>   vector[N] Z_2_1;
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
#>   vector<lower=0>[M_1] sd_1;  // group-level standard deviations
#>   array[M_1] vector[N_1] z_1;  // standardized group-level effects
#>   vector<lower=0>[M_2] sd_2;  // group-level standard deviations
#>   array[M_2] vector[N_2] z_2;  // standardized group-level effects
#> }
#> transformed parameters {
#>   vector[N_1] r_1_1;  // actual group-level effects
#>   vector[N_2] r_2_1;  // actual group-level effects
#>   // prior contributions to the log posterior
#>   real lprior = 0;
#>   r_1_1 = (sd_1[1] * (z_1[1]));
#>   r_2_1 = (sd_2[1] * (z_2[1]));
#>   lprior += student_t_lpdf(b[2] | 10, 0, 5);
#>   lprior += normal_lpdf(Intercept | 0,10);
#>   lprior += student_t_lpdf(sd_1 | 3, 0, 2.5)
#>     - 1 * student_t_lccdf(0 | 3, 0, 2.5);
#>   lprior += student_t_lpdf(sd_2 | 3, 0, 2.5)
#>     - 1 * student_t_lccdf(0 | 3, 0, 2.5);
#> }
#> model {
#>   // likelihood including constants
#>   if (!prior_only) {
#>     // initialize linear predictor term
#>     vector[N] mu = rep_vector(0.0, N);
#>     mu += Intercept;
#>     for (n in 1:N) {
#>       // add more terms to the linear predictor
#>       mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_2_1[J_2[n]] * Z_2_1[n];
#>     }
#>     target += poisson_log_glm_lpmf(Y | Xc, mu, b);
#>   }
#>   // priors including constants
#>   target += lprior;
#>   target += std_normal_lpdf(z_1[1]);
#>   target += std_normal_lpdf(z_2[1]);
#> }
#> generated quantities {
#>   // actual population-level intercept
#>   real b_Intercept = Intercept - dot_product(means_X, b);
#> }
```
