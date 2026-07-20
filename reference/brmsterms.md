# Parse Formulas of brms Models

Parse formulas objects for use in brms.

## Usage

``` r
brmsterms(formula, ...)

# Default S3 method
brmsterms(formula, ...)

# S3 method for class 'brmsformula'
brmsterms(formula, check_response = TRUE, resp_rhs_all = TRUE, ...)

# S3 method for class 'mvbrmsformula'
brmsterms(formula, ...)
```

## Arguments

- formula:

  An object of class [`formula`](https://rdrr.io/r/stats/formula.html),
  [`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md),
  or
  [`mvbrmsformula`](https://paulbuerkner.com/brms/reference/mvbrmsformula.md)
  (or one that can be coerced to that classes): A symbolic description
  of the model to be fitted. The details of model specification are
  explained in
  [`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md).

- ...:

  Further arguments passed to or from other methods.

- check_response:

  Logical; Indicates whether the left-hand side of `formula` (i.e.
  response variables and addition arguments) should be parsed. If
  `FALSE`, `formula` may also be one-sided.

- resp_rhs_all:

  Logical; Indicates whether to also include response variables on the
  right-hand side of formula `.$allvars`, where `.` represents the
  output of `brmsterms`.

## Value

An object of class `brmsterms` or `mvbrmsterms` (for multivariate
models), which is a `list` containing all required information initially
stored in `formula` in an easier to use format, basically a list of
formulas (not an abstract syntax tree).

## Details

This is the main formula parsing function of brms. It should usually not
be called directly, but is exported to allow package developers making
use of the formula syntax implemented in brms. As long as no other
packages depend on this functions, it may be changed without deprecation
warnings, when new features make this necessary.

## See also

[`brm`](https://paulbuerkner.com/brms/reference/brm.md),
[`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md),
[`mvbrmsformula`](https://paulbuerkner.com/brms/reference/mvbrmsformula.md)
