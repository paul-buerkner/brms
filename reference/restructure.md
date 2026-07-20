# Restructure Old R Objects

`restructure` is a generic function used to restructure old R objects to
work with newer versions of the package that generated them. Its
original use is within the brms package, but new methods for use with
objects from other packages can be registered to the same generic.

## Usage

``` r
restructure(x, ...)
```

## Arguments

- x:

  An object to be restructured. The object's class will determine which
  method to apply

- ...:

  Additional arguments to pass to the specific methods

## Value

An object of the same class as `x` compatible with the latest version of
the package that generated it.

## Details

Usually the version of the package that generated the object will be
stored somewhere in the object and this information will be used by the
specific method to determine what transformations to apply. See
[`restructure.brmsfit`](https://paulbuerkner.com/brms/reference/restructure.brmsfit.md)
for the default method applied for brms models. You can view the
available methods by typing: `methods(restructure)`

## See also

[`restructure.brmsfit`](https://paulbuerkner.com/brms/reference/restructure.brmsfit.md)
