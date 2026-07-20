# Check if cached fit can be used.

Checks whether a given cached fit can be used without refitting when
`file_refit = "on_change"` is used. This function is internal and
exposed only to facilitate debugging problems with cached fits. The
function may change or be removed in future versions and scripts should
not use it.

## Usage

``` r
brmsfit_needs_refit(
  fit,
  sdata = NULL,
  scode = NULL,
  data = NULL,
  algorithm = NULL,
  silent = FALSE,
  verbose = FALSE
)
```

## Arguments

- fit:

  Old `brmsfit` object (e.g., loaded from file).

- sdata:

  New Stan data (result of a call to
  [`standata`](https://paulbuerkner.com/brms/reference/standata.default.md)).
  Pass `NULL` to avoid this data check.

- scode:

  New Stan code (result of a call to
  [`stancode`](https://paulbuerkner.com/brms/reference/stancode.default.md)).
  Pass `NULL` to avoid this code check.

- data:

  New data to check consistency of factor level names. Pass `NULL` to
  avoid this data check.

- algorithm:

  New algorithm. Pass `NULL` to avoid algorithm check.

- silent:

  Logical. If `TRUE`, no messages will be given.

- verbose:

  Logical. If `TRUE` detailed report of the differences is printed to
  the console.

## Value

A boolean indicating whether a refit is needed.

## Details

Use with `verbose = TRUE` to get additional info on how the stored fit
differs from the given data and code.
