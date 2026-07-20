# Scaled inverse logit-link

Computes `inv_logit(x) * (ub - lb) + lb`

## Usage

``` r
inv_logit_scaled(x, lb = 0, ub = 1)
```

## Arguments

- x:

  A numeric or complex vector.

- lb:

  Lower bound defaulting to `0`.

- ub:

  Upper bound defaulting to `1`.

## Value

A numeric or complex vector between `lb` and `ub`.
