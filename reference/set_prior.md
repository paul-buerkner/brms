# Prior Definitions for brms Models

Define priors for specific parameters or classes of parameters.

## Usage

``` r
set_prior(
  prior,
  class = "b",
  coef = "",
  group = "",
  resp = "",
  dpar = "",
  nlpar = "",
  lb = NA,
  ub = NA,
  tag = "",
  check = TRUE
)

prior(prior, ...)

prior_(prior, ...)

prior_string(prior, ...)

empty_prior()
```

## Arguments

- prior:

  A character string defining a distribution in Stan language

- class:

  The parameter class. Defaults to `"b"` (i.e. population-level
  effects). See 'Details' for other valid parameter classes.

- coef:

  Name of the coefficient within the parameter class.

- group:

  Grouping factor for group-level parameters.

- resp:

  Name of the response variable. Only used in multivariate models.

- dpar:

  Name of a distributional parameter. Only used in distributional
  models.

- nlpar:

  Name of a non-linear parameter. Only used in non-linear models.

- lb:

  Lower bound for parameter restriction. Currently only allowed for
  classes `"b"`. Defaults to `NA`, that is no restriction.

- ub:

  Upper bound for parameter restriction. Currently only allowed for
  classes `"b"`. Defaults to `NA`, that is no restriction.

- tag:

  Character to append to the `lprior` variable in the Stan code. Used
  for selectively checking sensitivity of priors in `priorsense`.

- check:

  Logical; Indicates whether priors should be checked for validity (as
  far as possible). Defaults to `TRUE`. If `FALSE`, `prior` is passed to
  the Stan code as is, and all other arguments are ignored.

- ...:

  Arguments passed to `set_prior`.

## Value

An object of class `brmsprior` to be used in the `prior` argument of
[`brm`](https://paulbuerkner.com/brms/reference/brm.md).

## Details

`set_prior` is used to define prior distributions for parameters in brms
models. The functions `prior`, `prior_`, and `prior_string` are aliases
of `set_prior` each allowing for a different kind of argument
specification. `prior` allows specifying arguments as expression without
quotation marks using non-standard evaluation. `prior_` allows
specifying arguments as one-sided formulas or wrapped in `quote`.
`prior_string` allows specifying arguments as strings just as
`set_prior` itself.

Below, we explain its usage and list some common prior distributions for
parameters. A complete overview on possible prior distributions is given
in the Stan Reference Manual available at <https://mc-stan.org/>.

To combine multiple priors, use `c(...)` or the `+` operator (see
'Examples'). brms does not check if the priors are written in correct
Stan language. Instead, Stan will check their syntactical correctness
when the model is parsed to `C++` and returns an error if they are not.
This, however, does not imply that priors are always meaningful if they
are accepted by Stan. Although brms trys to find common problems (e.g.,
setting bounded priors on unbounded parameters), there is no guarantee
that the defined priors are reasonable for the model. Below, we list the
types of parameters in brms models, for which the user can specify prior
distributions.

Below, we provide details for the individual parameter classes that you
can set priors on. Often, it may not be immediately clear, which
parameters are present in the model. To get a full list of parameters
and parameter classes for which priors can be specified (depending on
the model) use function
[`default_prior`](https://paulbuerkner.com/brms/reference/default_prior.default.md).

1\. Population-level ('fixed') effects

Every Population-level effect has its own regression parameter
represents the name of the corresponding population-level effect.
Suppose, for instance, that `y` is predicted by `x1` and `x2` (i.e.,
`y ~ x1 + x2` in formula syntax). Then, `x1` and `x2` have regression
parameters `b_x1` and `b_x2` respectively. The default prior for
population-level effects (including monotonic and category specific
effects) is an improper flat prior over the reals. Other common options
are normal priors or student-t priors. If we want to have a normal prior
with mean 0 and standard deviation 5 for `x1`, and a unit student-t
prior with 10 degrees of freedom for `x2`, we can specify this via
`set_prior("normal(0,5)", class = "b", coef = "x1")` and  
`set_prior("student_t(10, 0, 1)", class = "b", coef = "x2")`. To put the
same prior on all population-level effects at once, we may write as a
shortcut `set_prior("<prior>", class = "b")`. This also leads to faster
sampling, because priors can be vectorized in this case. Both ways of
defining priors can be combined using for instance
`set_prior("normal(0, 2)", class = "b")` and  
`set_prior("normal(0, 10)", class = "b", coef = "x1")` at the same time.
This will set a `normal(0, 10)` prior on the effect of `x1` and a
`normal(0, 2)` prior on all other population-level effects. However,
this will break vectorization and may slow down the sampling procedure a
bit.

In case of the default intercept parameterization (discussed in the
'Details' section of
[`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md)),
general priors on class `"b"` will *not* affect the intercept. Instead,
the intercept has its own parameter class named `"Intercept"` and priors
can thus be specified via `set_prior("<prior>", class = "Intercept")`.
Setting a prior on the intercept will not break vectorization of the
other population-level effects. Note that technically, this prior is set
on an intercept that results when internally centering all
population-level predictors around zero to improve sampling efficiency.
On this centered intercept, specifying a prior is actually much easier
and intuitive than on the original intercept, since the former
represents the expected response value when all predictors are at their
means. To treat the intercept as an ordinary population-level effect and
avoid the centering parameterization, use `0 + Intercept` on the
right-hand side of the model formula.

In non-linear models, population-level effects are defined separately
for each non-linear parameter. Accordingly, it is necessary to specify
the non-linear parameter in `set_prior` so that priors we can be
assigned correctly. If, for instance, `alpha` is the parameter and `x`
the predictor for which we want to define the prior, we can write
`set_prior("<prior>", coef = "x", nlpar = "alpha")`. As a shortcut we
can use `set_prior("<prior>", nlpar = "alpha")` to set the same prior on
all population-level effects of `alpha` at once.

The same goes for specifying priors for specific distributional
parameters in the context of distributional regression, for example,
`set_prior("<prior>", coef = "x", dpar = "sigma")`. For most other
parameter classes (see below), you need to indicate non-linear and
distributional parameters in the same way as shown here.

If desired, population-level effects can be restricted to fall only
within a certain interval using the `lb` and `ub` arguments of
`set_prior`. This is often required when defining priors that are not
defined everywhere on the real line, such as uniform or gamma priors.
When defining a `uniform(2,4)` prior, you should write
`set_prior("uniform(2,4)", lb = 2, ub = 4)`. When using a prior that is
defined on the positive reals only (such as a gamma prior) set `lb = 0`.
In most situations, it is not useful to restrict population-level
parameters through bounded priors (non-linear models are an important
exception), but if you really want to this is the way to go.

2\. Group-level ('random') effects

Each group-level effect of each grouping factor has a standard deviation
named `sd_<group>_<coef>`. Consider, for instance, the formula
`y ~ x1 + x2 + (1 + x1 | g)`. We see that the intercept as well as `x1`
are group-level effects nested in the grouping factor `g`. The
corresponding standard deviation parameters are named as
`sd_g_Intercept` and `sd_g_x1` respectively. These parameters are
restricted to be non-negative and, by default, have a half student-t
prior with 3 degrees of freedom and a scale parameter that depends on
the standard deviation of the response after applying the link function.
Minimally, the scale parameter is 2.5. This prior is used (a) to be only
weakly informative in order to influence results as few as possible,
while (b) providing at least some regularization to considerably improve
convergence and sampling efficiency. To define a prior distribution only
for standard deviations of a specific grouping factor, use  
`set_prior("<prior>", class = "sd", group = "<group>")`. To define a
prior distribution only for a specific standard deviation of a specific
grouping factor, you may write  
`set_prior("<prior>", class = "sd", group = "<group>", coef = "<coef>")`.

If there is more than one group-level effect per grouping factor, the
correlations between those effects have to be estimated. The prior
`lkj_corr_cholesky(eta)` or in short `lkj(eta)` with `eta > 0` is
essentially the only prior for (Cholesky factors) of correlation
matrices. If `eta = 1` (the default) all correlations matrices are
equally likely a priori. If `eta > 1`, extreme correlations become less
likely, whereas `0 < eta < 1` results in higher probabilities for
extreme correlations. Correlation matrix parameters in `brms` models are
named as `cor_<group>`, (e.g., `cor_g` if `g` is the grouping factor).
To set the same prior on every correlation matrix, use for instance
`set_prior("lkj(2)", class = "cor")`. Internally, the priors are
transformed to be put on the Cholesky factors of the correlation
matrices to improve efficiency and numerical stability. The
corresponding parameter class of the Cholesky factors is `L`, but it is
not recommended to specify priors for this parameter class directly.

4\. Smoothing Splines

Smoothing splines are implemented in brms using the 'random effects'
formulation as explained in
[`gamm`](https://rdrr.io/pkg/mgcv/man/gamm.html)). Thus, each spline has
its corresponding standard deviations modeling the variability within
this term. In brms, this parameter class is called `sds` and priors can
be specified via
`set_prior("<prior>", class = "sds", coef = "<term label>")`. The
default prior is the same as for standard deviations of group-level
effects.

5\. Gaussian processes

Gaussian processes as currently implemented in brms have two parameters,
the standard deviation parameter `sdgp`, and characteristic length-scale
parameter `lscale` (see
[`gp`](https://paulbuerkner.com/brms/reference/gp.md) for more details).
The default prior of `sdgp` is the same as for standard deviations of
group-level effects. The default prior of `lscale` is an informative
inverse-gamma prior specifically tuned to the covariates of the Gaussian
process (for more details see
<https://betanalpha.github.io/assets/case_studies/gp_part3/part3.html>).
This tuned prior may be overly informative in some cases, so please
consider other priors as well to make sure inference is robust to the
prior specification. If tuning fails, a half-normal prior is used
instead.

6\. Autocorrelation parameters

The autocorrelation parameters currently implemented are named `ar`
(autoregression), `ma` (moving average), `sderr` (standard deviation of
latent residuals in latent ARMA models), `cosy` (compound symmetry
correlation), `car` (spatial conditional autoregression), as well as
`lagsar` and `errorsar` (spatial simultaneous autoregression).

Priors can be defined by `set_prior("<prior>", class = "ar")` for `ar`
and similar for other autocorrelation parameters. By default, `ar` and
`ma` are bounded between `-1` and `1`; `cosy`, `car`, `lagsar`, and
`errorsar` are bounded between `0` and `1`. The default priors are flat
over the respective definition areas.

7\. Parameters of measurement error terms

Latent variables induced via measurement error
[`me`](https://paulbuerkner.com/brms/reference/me.md) terms require both
mean and standard deviation parameters, whose prior classes are named
`"meanme"` and `"sdme"`, respectively. If multiple latent variables are
induced this way, their correlation matrix will be modeled as well and
corresponding priors can be specified via the `"corme"` class. All of
the above parameters have flat priors over their respective definition
spaces by default.

8\. Distance parameters of monotonic effects

As explained in the details section of
[`brm`](https://paulbuerkner.com/brms/reference/brm.md), monotonic
effects make use of a special parameter vector to estimate the
'normalized distances' between consecutive predictor categories. This is
realized in Stan using the `simplex` parameter type. This class is named
`"simo"` (short for simplex monotonic) in brms. The only valid prior for
simplex parameters is the dirichlet prior, which accepts a vector of
length `K - 1` (K = number of predictor categories) as input defining
the 'concentration' of the distribution. Explaining the dirichlet prior
is beyond the scope of this documentation, but we want to describe how
to define this prior syntactically correct. If a predictor `x` with `K`
categories is modeled as monotonic, we can define a prior on its
corresponding simplex via  
`prior(dirichlet(<vector>), class = simo, coef = mox1)`. The `1` in the
end of `coef` indicates that this is the first simplex in this term. If
interactions between multiple monotonic variables are modeled, multiple
simplexes per term are required. For `<vector>`, we can put in any `R`
expression defining a vector of length `K - 1`. The default is a uniform
prior (i.e. `<vector> = rep(1, K-1)`) over all simplexes of the
respective dimension.

9\. Parameters for specific families

Some families need additional parameters to be estimated. Families
`gaussian`, `student`, `skew_normal`, `lognormal`, and
`gen_extreme_value` need the parameter `sigma` to account for the
residual standard deviation. By default, `sigma` has a half student-t
prior that scales in the same way as the group-level standard
deviations. Further, family `student` needs the parameter `nu`
representing the degrees of freedom of Student-t distribution. By
default, `nu` has prior `gamma(2, 0.1)`, which is close to a penalized
complexity prior (see Stan prior choice Wiki), and a fixed lower bound
of `1`. Family `negbinomial` needs a `shape` parameter that has by
default `inv_gamma(0.4, 0.3)` prior which is close to a penalized
complexity prior (see Stan prior choice Wiki). Families `gamma`,
`weibull`, and `inverse.gaussian`, need a `shape` parameter that has a
`gamma(0.01, 0.01)` prior by default. For families `cumulative`,
`cratio`, `sratio`, and `acat`, and only if `threshold = "equidistant"`,
the parameter `delta` is used to model the distance between two adjacent
thresholds. By default, `delta` has an improper flat prior over the
reals. The `von_mises` family needs the parameter `kappa`, representing
the concentration parameter. By default, `kappa` has prior
`gamma(2, 0.01)`.

Every family specific parameter has its own prior class, so that
`set_prior("<prior>", class = "<parameter>")` is the right way to go.
All of these priors are chosen to be weakly informative, having only
minimal influence on the estimations, while improving convergence and
sampling efficiency.

10\. Shrinkage priors

To reduce the danger of overfitting in models with many predictor terms
fit on comparably sparse data, brms supports special shrinkage priors,
namely the (regularized)
[`horseshoe`](https://paulbuerkner.com/brms/reference/horseshoe.md) and
the [`R2D2`](https://paulbuerkner.com/brms/reference/R2D2.md) prior.
These priors can be applied on many parameter classes, either directly
on the coefficient classes (e.g., class `b`), if directly setting priors
on them is supported, or on the corresponding standard deviation
hyperparameters (e.g., class `sd`) otherwise. Currently, the following
classes support shrinkage priors: `b` (overall regression coefficients),
`sds` (SDs of smoothing splines), `sdgp` (SDs of Gaussian processes),
`ar` (autoregressive coefficients), `ma` (moving average coefficients),
`sderr` (SD of latent residuals), `sdcar` (SD of spatial CAR
structures), `sd` (SD of varying coefficients).

11\. Fixing parameters to constants

Fixing parameters to constants is possible by using the `constant`
function, for example, `constant(1)` to fix a parameter to 1.
Broadcasting to vectors and matrices is done automatically.

## Functions

- `prior()`: Alias of `set_prior` allowing to specify arguments as
  expressions without quotation marks.

- `prior_()`: Alias of `set_prior` allowing to specify arguments as as
  one-sided formulas or wrapped in `quote`.

- `prior_string()`: Alias of `set_prior` allowing to specify arguments
  as strings.

- `empty_prior()`: Create an empty `brmsprior` object.

## See also

[`default_prior`](https://paulbuerkner.com/brms/reference/default_prior.default.md)

## Examples

``` r
## use alias functions
(prior1 <- prior(cauchy(0, 1), class = sd))
#> sd ~ cauchy(0, 1)
(prior2 <- prior_(~cauchy(0, 1), class = ~sd))
#> sd ~ cauchy(0, 1)
(prior3 <- prior_string("cauchy(0, 1)", class = "sd"))
#> sd ~ cauchy(0, 1)
identical(prior1, prior2)
#> [1] TRUE
identical(prior1, prior3)
#> [1] TRUE

# check which parameters can have priors
default_prior(rating ~ treat + period + carry + (1|subject),
             data = inhaler, family = cumulative())
#>                 prior     class      coef   group resp dpar nlpar lb ub tag
#>  student_t(3, 0, 2.5) Intercept                                            
#>  student_t(3, 0, 2.5) Intercept         1                                  
#>  student_t(3, 0, 2.5) Intercept         2                                  
#>  student_t(3, 0, 2.5) Intercept         3                                  
#>                (flat)         b                                            
#>                (flat)         b     carry                                  
#>                (flat)         b    period                                  
#>                (flat)         b     treat                                  
#>  student_t(3, 0, 2.5)        sd                                    0       
#>  student_t(3, 0, 2.5)        sd           subject                  0       
#>  student_t(3, 0, 2.5)        sd Intercept subject                  0       
#>        source
#>       default
#>  (vectorized)
#>  (vectorized)
#>  (vectorized)
#>       default
#>  (vectorized)
#>  (vectorized)
#>  (vectorized)
#>       default
#>  (vectorized)
#>  (vectorized)

# define some priors
bprior <- c(prior_string("normal(0,10)", class = "b"),
            prior(normal(1,2), class = b, coef = treat),
            prior_(~cauchy(0,2), class = ~sd,
                   group = ~subject, coef = ~Intercept))

# verify that the priors indeed found their way into Stan's model code
stancode(rating ~ treat + period + carry + (1|subject),
         data = inhaler, family = cumulative(),
         prior = bprior)
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
#>   lprior += normal_lpdf(b[1] | 1, 2);
#>   lprior += normal_lpdf(b[2] | 0,10);
#>   lprior += normal_lpdf(b[3] | 0,10);
#>   lprior += student_t_lpdf(Intercept | 3, 0, 2.5);
#>   lprior += cauchy_lpdf(sd_1[1] | 0, 2)
#>     - 1 * cauchy_lccdf(0 | 0, 2);
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

# use the horseshoe prior to model sparsity in regression coefficients
stancode(count ~ zAge + zBase * Trt,
         data = epilepsy, family = poisson(),
         prior = set_prior("horseshoe(3)"))
#> // generated with brms 2.23.1
#> functions {
#>   /* Efficient computation of the horseshoe scale parameters
#>    * see Appendix C.1 in https://projecteuclid.org/euclid.ejs/1513306866
#>    * Args:
#>    *   lambda: local shrinkage parameters
#>    *   tau: global shrinkage parameter
#>    *   c2: slap regularization parameter
#>    * Returns:
#>    *   scale parameter vector of the horseshoe prior
#>    */
#>   vector scales_horseshoe(vector lambda, real tau, real c2) {
#>     int K = rows(lambda);
#>     vector[K] lambda2 = square(lambda);
#>     vector[K] lambda_tilde = sqrt(c2 * lambda2 ./ (c2 + tau^2 * lambda2));
#>     return lambda_tilde * tau;
#>   }
#> }
#> data {
#>   int<lower=1> N;  // total number of observations
#>   array[N] int Y;  // response variable
#>   int<lower=1> K;  // number of population-level effects
#>   matrix[N, K] X;  // population-level design matrix
#>   int<lower=1> Kc;  // number of population-level effects after centering
#>   int<lower=1> Kscales;  // number of local scale parameters
#>   // data for the horseshoe prior
#>   real<lower=0> hs_df;  // local degrees of freedom
#>   real<lower=0> hs_df_global;  // global degrees of freedom
#>   real<lower=0> hs_df_slab;  // slab degrees of freedom
#>   real<lower=0> hs_scale_global;  // global prior scale
#>   real<lower=0> hs_scale_slab;  // slab prior scale
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
#>   vector[Kc] zb;  // unscaled coefficients
#>   real Intercept;  // temporary intercept for centered predictors
#>   // horseshoe shrinkage parameters
#>   real<lower=0> hs_global;  // global shrinkage parameter
#>   real<lower=0> hs_slab;  // slab regularization parameter
#>   vector<lower=0>[Kscales] hs_local;  // local parameters for the horseshoe prior
#> }
#> transformed parameters {
#>   vector[Kc] b;  // scaled coefficients
#>   vector<lower=0>[Kc] sdb;  // SDs of the coefficients
#>   vector<lower=0>[Kscales] scales;  // local horseshoe scale parameters
#>   // prior contributions to the log posterior
#>   real lprior = 0;
#>   // compute horseshoe scale parameters
#>   scales = scales_horseshoe(hs_local, hs_global, hs_scale_slab^2 * hs_slab);
#>   sdb = scales[(1):(Kc)];
#>   b = zb .* sdb;  // scale coefficients
#>   lprior += student_t_lpdf(Intercept | 3, 1.4, 2.5);
#>   lprior += student_t_lpdf(hs_global | hs_df_global, 0, hs_scale_global)
#>     - 1 * log(0.5);
#>   lprior += inv_gamma_lpdf(hs_slab | 0.5 * hs_df_slab, 0.5 * hs_df_slab);
#> }
#> model {
#>   // likelihood including constants
#>   if (!prior_only) {
#>     target += poisson_log_glm_lpmf(Y | Xc, Intercept, b);
#>   }
#>   // priors including constants
#>   target += lprior;
#>   target += std_normal_lpdf(zb);
#>   target += student_t_lpdf(hs_local | hs_df, 0, 1)
#>     - rows(hs_local) * log(0.5);
#> }
#> generated quantities {
#>   // actual population-level intercept
#>   real b_Intercept = Intercept - dot_product(means_X, b);
#> }

# fix certain priors to constants
bprior <- prior(constant(1), class = "b") +
  prior(constant(2), class = "b", coef = "zBase") +
  prior(constant(0.5), class = "sd")
stancode(count ~ zAge + zBase + (1 | patient),
              data = epilepsy, prior = bprior)
#> // generated with brms 2.23.1
#> functions {
#> }
#> data {
#>   int<lower=1> N;  // total number of observations
#>   vector[N] Y;  // response variable
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
#>   array[M_1] vector[N_1] z_1;  // standardized group-level effects
#> }
#> transformed parameters {
#>   vector[Kc] b;  // regression coefficients
#>   vector<lower=0>[M_1] sd_1;  // group-level standard deviations
#>   vector[N_1] r_1_1;  // actual group-level effects
#>   // prior contributions to the log posterior
#>   real lprior = 0;
#>   b[1] = 1;
#>   b[2] = 2;
#>   sd_1 = rep_vector(0.5, rows(sd_1));
#>   r_1_1 = (sd_1[1] * (z_1[1]));
#>   lprior += student_t_lpdf(Intercept | 3, 4, 4.4);
#>   lprior += student_t_lpdf(sigma | 3, 0, 4.4)
#>     - 1 * student_t_lccdf(0 | 3, 0, 4.4);
#> }
#> model {
#>   // likelihood including constants
#>   if (!prior_only) {
#>     // initialize linear predictor term
#>     vector[N] mu = rep_vector(0.0, N);
#>     mu += Intercept;
#>     for (n in 1:N) {
#>       // add more terms to the linear predictor
#>       mu[n] += r_1_1[J_1[n]] * Z_1_1[n];
#>     }
#>     target += normal_id_glm_lpdf(Y | Xc, mu, b, sigma);
#>   }
#>   // priors including constants
#>   target += lprior;
#>   target += std_normal_lpdf(z_1[1]);
#> }
#> generated quantities {
#>   // actual population-level intercept
#>   real b_Intercept = Intercept - dot_product(means_X, b);
#> }

# pass priors to Stan without checking
prior <- prior_string("target += normal_lpdf(b[1] | 0, 1)", check = FALSE)
stancode(count ~ Trt, data = epilepsy, prior = prior)
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
#>     target += normal_id_glm_lpdf(Y | Xc, Intercept, b, sigma);
#>   }
#>   // priors including constants
#>   target += lprior;
#>   target += normal_lpdf(b[1] | 0, 1);
#> }
#> generated quantities {
#>   // actual population-level intercept
#>   real b_Intercept = Intercept - dot_product(means_X, b);
#> }

# define priors in a vectorized manner
# useful in particular for categorical or multivariate models
set_prior("normal(0, 2)", dpar = c("muX", "muY", "muZ"))
#>         prior class coef group resp dpar nlpar   lb   ub tag source
#>  normal(0, 2)     b                  muX       <NA> <NA>       user
#>  normal(0, 2)     b                  muY       <NA> <NA>       user
#>  normal(0, 2)     b                  muZ       <NA> <NA>       user

# specify tags for different priors for sensitivity analysis using priorsense
# It is then possible to check the sensitivity when changing the priors with
# the same tag while leaving others the same
prior_cov <- prior(normal(0, 10), class = "b", tag = "covariates")
prior_trt <- prior(normal(0, 1), class = "b", coef = "Trt1", tag = "treatment")
stancode(count ~ Trt + zAge + zBase + (1 | patient),
         data = epilepsy, prior = c(prior_cov, prior_trt))
#> // generated with brms 2.23.1
#> functions {
#> }
#> data {
#>   int<lower=1> N;  // total number of observations
#>   vector[N] Y;  // response variable
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
#>   vector<lower=0>[M_1] sd_1;  // group-level standard deviations
#>   array[M_1] vector[N_1] z_1;  // standardized group-level effects
#> }
#> transformed parameters {
#>   vector[N_1] r_1_1;  // actual group-level effects
#>   // prior contributions to the log posterior
#>   real lprior = 0;
#>   real lprior_covariates = 0;
#>   real lprior_treatment = 0;
#>   r_1_1 = (sd_1[1] * (z_1[1]));
#>   lprior_treatment += normal_lpdf(b[1] | 0, 1);
#>   lprior_covariates += normal_lpdf(b[2] | 0, 10);
#>   lprior_covariates += normal_lpdf(b[3] | 0, 10);
#>   lprior += student_t_lpdf(Intercept | 3, 4, 4.4);
#>   lprior += student_t_lpdf(sigma | 3, 0, 4.4)
#>     - 1 * student_t_lccdf(0 | 3, 0, 4.4);
#>   lprior += student_t_lpdf(sd_1 | 3, 0, 4.4)
#>     - 1 * student_t_lccdf(0 | 3, 0, 4.4);
#>   lprior += lprior_covariates;
#>   lprior += lprior_treatment;
#> }
#> model {
#>   // likelihood including constants
#>   if (!prior_only) {
#>     // initialize linear predictor term
#>     vector[N] mu = rep_vector(0.0, N);
#>     mu += Intercept;
#>     for (n in 1:N) {
#>       // add more terms to the linear predictor
#>       mu[n] += r_1_1[J_1[n]] * Z_1_1[n];
#>     }
#>     target += normal_id_glm_lpdf(Y | Xc, mu, b, sigma);
#>   }
#>   // priors including constants
#>   target += lprior;
#>   target += std_normal_lpdf(z_1[1]);
#> }
#> generated quantities {
#>   // actual population-level intercept
#>   real b_Intercept = Intercept - dot_product(means_X, b);
#> }
```
