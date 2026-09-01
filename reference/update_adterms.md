# Update Formula Addition Terms

Update additions terms used in formulas of brms. See
[`addition-terms`](https://paulbuerkner.com/brms/reference/addition-terms.md)
for details.

## Usage

``` r
update_adterms(formula, adform, action = c("update", "replace"))
```

## Arguments

- formula:

  Two-sided formula to be updated.

- adform:

  One-sided formula containing addition terms to update `formula` with.

- action:

  Indicates what should happen to the existing addition terms in
  `formula`. If `"update"` (the default), old addition terms that have
  no corresponding term in `adform` will be kept. If `"replace"`, all
  old addition terms will be removed.

## Value

An object of class `formula`.

## Examples

``` r
form <- y | trials(size) ~ x
update_adterms(form, ~ trials(10))
#> y | trials(10) ~ x
#> <environment: 0x564056a15238>
update_adterms(form, ~ weights(w))
#> y | trials(size) + weights(w) ~ x
#> <environment: 0x564056a15238>
update_adterms(form, ~ weights(w), action = "replace")
#> y | weights(w) ~ x
#> <environment: 0x564056a15238>
update_adterms(y ~ x, ~ trials(10))
#> y | trials(10) ~ x
#> <environment: 0x564056a15238>
```
