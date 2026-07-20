# Class `brmsfit` of models fitted with the brms package

Models fitted with the
[`brms`](https://paulbuerkner.com/brms/reference/brms-package.md)
package are represented as a `brmsfit` object, which contains the
posterior draws (samples), model formula, Stan code, relevant data, and
other information.

## Details

See `methods(class = "brmsfit")` for an overview of available methods.

## Slots

- `formula`:

  A
  [`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md)
  object.

- `data`:

  A `data.frame` containing all variables used in the model.

- `data2`:

  A `list` of data objects which cannot be passed via `data`.

- `prior`:

  A [`brmsprior`](https://paulbuerkner.com/brms/reference/set_prior.md)
  object containing information on the priors used in the model.

- `stanvars`:

  A [`stanvars`](https://paulbuerkner.com/brms/reference/stanvar.md)
  object.

- `model`:

  The model code in Stan language.

- `exclude`:

  The names of the parameters for which draws are not saved.

- `algorithm`:

  The name of the algorithm used to fit the model.

- `backend`:

  The name of the backend used to fit the model.

- `threads`:

  An object of class \`brmsthreads\` created by
  [`threading`](https://paulbuerkner.com/brms/reference/threading.md).

- `opencl`:

  An object of class \`brmsopencl\` created by
  [`opencl`](https://paulbuerkner.com/brms/reference/opencl.md).

- `stan_args`:

  Named list of additional control arguments that were passed to the
  Stan backend directly.

- `fit`:

  An object of class
  [`stanfit`](https://mc-stan.org/rstan/reference/stanfit-class.html)
  among others containing the posterior draws.

- `basis`:

  An object that contains a small subset of the Stan data created at
  fitting time, which is needed to process new data correctly.

- `criteria`:

  An empty `list` for adding model fit criteria after estimation of the
  model.

- `file`:

  Optional name of a file in which the model object was stored in or
  loaded from.

- `version`:

  The versions of brms and rstan with which the model was fitted.

- `family`:

  (Deprecated) A
  [`brmsfamily`](https://paulbuerkner.com/brms/reference/brmsfamily.md)
  object.

- `autocor`:

  (Deprecated) An
  [`cor_brms`](https://paulbuerkner.com/brms/reference/cor_brms.md)
  object containing the autocorrelation structure if specified.

- `ranef`:

  (Deprecated) A `data.frame` containing the group-level structure.

- `cov_ranef`:

  (Deprecated) A `list` of customized group-level covariance matrices.

- `stan_funs`:

  (Deprecated) A character string of length one or `NULL`.

- `data.name`:

  (Deprecated) The name of `data` as specified by the user.

## See also

[`brms`](https://paulbuerkner.com/brms/reference/brms-package.md),
[`brm`](https://paulbuerkner.com/brms/reference/brm.md),
[`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md),
[`brmsfamily`](https://paulbuerkner.com/brms/reference/brmsfamily.md)
