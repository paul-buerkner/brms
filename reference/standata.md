# Stan data for Bayesian models

`standata` is a generic function that can be used to generate data for
Bayesian models to be passed to Stan. Its original use is within the
brms package, but new methods for use with objects from other packages
can be registered to the same generic.

## Usage

``` r
standata(object, ...)

make_standata(formula, ...)
```

## Arguments

- object:

  A formula object whose class will determine which method will be used.
  A symbolic description of the model to be fitted.

- ...:

  Further arguments passed to the specific method.

- formula:

  Synonym of `object` for use in `make_standata`.

## Value

A named list of objects containing the required data to fit a Bayesian
model with Stan.

## Details

See
[`standata.default`](https://paulbuerkner.com/brms/reference/standata.default.md)
for the default method applied for brms models. You can view the
available methods by typing `methods(standata)`. The `make_standata`
function is an alias of `standata`.

## See also

[`standata.default`](https://paulbuerkner.com/brms/reference/standata.default.md),
[`standata.brmsfit`](https://paulbuerkner.com/brms/reference/standata.brmsfit.md)

## Examples

``` r
sdata1 <- standata(rating ~ treat + period + carry + (1|subject),
                   data = inhaler, family = "cumulative")
str(sdata1)
#> List of 13
#>  $ N         : int 572
#>  $ Y         : num [1:572(1d)] 1 1 1 1 1 1 1 1 1 1 ...
#>  $ nthres    : int 3
#>  $ K         : int 3
#>  $ Kc        : num 3
#>  $ X         : num [1:572, 1:3] 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:572] "1" "2" "3" "4" ...
#>   .. ..$ : chr [1:3] "treat" "period" "carry"
#>  $ Z_1_1     : num [1:572(1d)] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..- attr(*, "dimnames")=List of 1
#>   .. ..$ : chr [1:572] "1" "2" "3" "4" ...
#>  $ disc      : num 1
#>  $ J_1       : int [1:572(1d)] 1 2 3 4 5 6 7 8 9 10 ...
#>  $ N_1       : int 286
#>  $ M_1       : int 1
#>  $ NC_1      : int 0
#>  $ prior_only: int 0
#>  - attr(*, "class")= chr [1:2] "standata" "list"
```
