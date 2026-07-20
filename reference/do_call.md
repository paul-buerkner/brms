# Execute a Function Call

Execute a function call similar to
[`do.call`](https://rdrr.io/r/base/do.call.html), but without deparsing
function arguments. For large number of arguments (i.e., more than a few
thousand) this function currently is somewhat inefficient and should be
used with care in this case.

## Usage

``` r
do_call(what, args, pkg = NULL, envir = parent.frame())
```

## Arguments

- what:

  Either a function or a non-empty character string naming the function
  to be called.

- args:

  A list of arguments to the function call. The names attribute of
  `args` gives the argument names.

- pkg:

  Optional name of the package in which to search for the function if
  `what` is a character string.

- envir:

  An environment within which to evaluate the call.

## Value

The result of the (evaluated) function call.
