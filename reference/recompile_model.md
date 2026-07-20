# Recompile Stan models in `brmsfit` objects

Recompile the Stan model inside a `brmsfit` object, if necessary. This
does not change the model, it simply recreates the executable so that
sampling is possible again.

## Usage

``` r
recompile_model(x, recompile = NULL)
```

## Arguments

- x:

  An object of class `brmsfit`.

- recompile:

  Logical, indicating whether the Stan model should be recompiled. If
  `NULL` (the default), `recompile_model` tries to figure out
  internally, if recompilation is necessary. Setting it to `FALSE` will
  cause `recompile_model` to always return the `brmsfit` object
  unchanged.

## Value

A (possibly updated) `brmsfit` object.
