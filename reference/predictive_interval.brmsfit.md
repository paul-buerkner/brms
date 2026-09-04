# Predictive Intervals

Compute intervals from the posterior predictive distribution.

## Usage

``` r
# S3 method for class 'brmsfit'
predictive_interval(object, prob = 0.9, ...)
```

## Arguments

- object:

  An R object of class `brmsfit`.

- prob:

  A number p (0 \< p \< 1) indicating the desired probability mass to
  include in the intervals. Defaults to `0.9`.

- ...:

  Further arguments passed to
  [`posterior_predict`](https://paulbuerkner.com/brms/reference/posterior_predict.brmsfit.md).

## Value

A matrix with 2 columns for the lower and upper bounds of the intervals,
respectively, and as many rows as observations being predicted.

## Examples

``` r
# \dontrun{
fit <- brm(count ~ zBase, data = epilepsy, family = poisson())
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 2.4e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.24 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.046 seconds (Warm-up)
#> Chain 1:                0.044 seconds (Sampling)
#> Chain 1:                0.09 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 1.1e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.11 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 0.046 seconds (Warm-up)
#> Chain 2:                0.045 seconds (Sampling)
#> Chain 2:                0.091 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 1.1e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.11 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 0.045 seconds (Warm-up)
#> Chain 3:                0.046 seconds (Sampling)
#> Chain 3:                0.091 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 1e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.1 seconds.
#> Chain 4: Adjust your expectations accordingly!
#> Chain 4: 
#> Chain 4: 
#> Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 4: 
#> Chain 4:  Elapsed Time: 0.044 seconds (Warm-up)
#> Chain 4:                0.041 seconds (Sampling)
#> Chain 4:                0.085 seconds (Total)
#> Chain 4: 
predictive_interval(fit)
#>        5%   95%
#>   [1,]  1  8.00
#>   [2,]  1  8.00
#>   [3,]  1  7.00
#>   [4,]  1  7.00
#>   [5,]  8 19.00
#>   [6,]  2 10.00
#>   [7,]  1  8.00
#>   [8,]  5 15.00
#>   [9,]  2  9.00
#>  [10,]  1  8.00
#>  [11,]  5 16.00
#>  [12,]  3 11.00
#>  [13,]  1  9.00
#>  [14,]  4 13.00
#>  [15,] 14 29.00
#>  [16,]  5 15.00
#>  [17,]  1  9.00
#>  [18,] 25 45.00
#>  [19,]  2  9.00
#>  [20,]  2  9.00
#>  [21,]  1  8.00
#>  [22,]  1  8.00
#>  [23,]  2  9.00
#>  [24,]  2 10.00
#>  [25,]  6 16.00
#>  [26,]  1  8.00
#>  [27,]  1  8.00
#>  [28,]  4 14.00
#>  [29,] 10 23.00
#>  [30,]  3 12.00
#>  [31,]  2  9.00
#>  [32,]  1  8.00
#>  [33,]  2  9.00
#>  [34,]  2  9.00
#>  [35,]  3 11.00
#>  [36,]  1  8.00
#>  [37,]  1  8.00
#>  [38,]  8 20.00
#>  [39,]  3 13.00
#>  [40,]  1  7.00
#>  [41,]  2  9.00
#>  [42,]  1  8.00
#>  [43,]  4 14.00
#>  [44,]  3 12.00
#>  [45,]  3 12.00
#>  [46,]  1  7.00
#>  [47,]  3 12.00
#>  [48,]  1  8.00
#>  [49,] 66 98.00
#>  [50,]  2 10.00
#>  [51,]  3 13.00
#>  [52,]  3 11.00
#>  [53,]  6 17.00
#>  [54,]  2 10.00
#>  [55,]  1  8.00
#>  [56,]  2  9.00
#>  [57,]  2 10.00
#>  [58,]  1  8.00
#>  [59,]  1  8.00
#>  [60,]  1  8.00
#>  [61,]  1  8.00
#>  [62,]  1  7.00
#>  [63,]  1  7.00
#>  [64,]  8 20.00
#>  [65,]  2 10.00
#>  [66,]  1  8.00
#>  [67,]  5 15.00
#>  [68,]  2  9.00
#>  [69,]  1  8.00
#>  [70,]  5 15.00
#>  [71,]  3 11.00
#>  [72,]  1  9.00
#>  [73,]  4 13.00
#>  [74,] 14 29.00
#>  [75,]  5 15.00
#>  [76,]  2  9.00
#>  [77,] 25 45.00
#>  [78,]  1  9.00
#>  [79,]  2  9.00
#>  [80,]  1  8.00
#>  [81,]  1  7.00
#>  [82,]  1  9.00
#>  [83,]  2 10.00
#>  [84,]  5 16.00
#>  [85,]  1  7.00
#>  [86,]  1  8.00
#>  [87,]  4 14.00
#>  [88,] 10 24.00
#>  [89,]  3 12.00
#>  [90,]  2  9.00
#>  [91,]  1  8.00
#>  [92,]  2  9.00
#>  [93,]  2 10.00
#>  [94,]  3 11.00
#>  [95,]  1  8.00
#>  [96,]  1  8.00
#>  [97,]  8 20.00
#>  [98,]  4 13.00
#>  [99,]  1  7.00
#> [100,]  2  9.00
#> [101,]  1  8.00
#> [102,]  4 14.00
#> [103,]  3 12.00
#> [104,]  3 12.00
#> [105,]  1  7.00
#> [106,]  3 12.00
#> [107,]  1  8.00
#> [108,] 67 98.00
#> [109,]  2  9.00
#> [110,]  4 13.00
#> [111,]  3 11.00
#> [112,]  6 16.00
#> [113,]  2 10.00
#> [114,]  1  8.00
#> [115,]  2  9.00
#> [116,]  2 10.00
#> [117,]  1  8.00
#> [118,]  1  8.00
#> [119,]  1  8.00
#> [120,]  1  8.00
#> [121,]  1  7.00
#> [122,]  1  7.00
#> [123,]  8 20.00
#> [124,]  2 10.00
#> [125,]  1  8.00
#> [126,]  5 15.00
#> [127,]  2  9.00
#> [128,]  1  8.00
#> [129,]  5 16.00
#> [130,]  3 11.00
#> [131,]  2  9.00
#> [132,]  4 13.00
#> [133,] 14 29.00
#> [134,]  5 15.00
#> [135,]  2  9.00
#> [136,] 25 45.00
#> [137,]  2  9.00
#> [138,]  2  9.00
#> [139,]  1  8.00
#> [140,]  1  7.00
#> [141,]  1  8.00
#> [142,]  2 10.00
#> [143,]  6 16.00
#> [144,]  1  8.00
#> [145,]  1  7.00
#> [146,]  4 14.00
#> [147,] 10 23.00
#> [148,]  3 12.00
#> [149,]  2  9.00
#> [150,]  1  8.00
#> [151,]  2  9.00
#> [152,]  2  9.00
#> [153,]  3 11.00
#> [154,]  1  8.00
#> [155,]  1  8.00
#> [156,]  8 20.00
#> [157,]  4 13.00
#> [158,]  1  7.00
#> [159,]  2  9.00
#> [160,]  1  8.00
#> [161,]  4 14.00
#> [162,]  3 12.00
#> [163,]  3 12.00
#> [164,]  1  7.00
#> [165,]  3 12.00
#> [166,]  1  8.00
#> [167,] 66 99.00
#> [168,]  2  9.00
#> [169,]  4 13.00
#> [170,]  3 11.00
#> [171,]  6 17.00
#> [172,]  2 10.00
#> [173,]  1  8.00
#> [174,]  2  9.00
#> [175,]  2 10.00
#> [176,]  1  8.00
#> [177,]  1  8.00
#> [178,]  1  8.00
#> [179,]  1  8.00
#> [180,]  1  7.00
#> [181,]  1  7.00
#> [182,]  8 20.00
#> [183,]  2 10.00
#> [184,]  1  8.00
#> [185,]  5 15.00
#> [186,]  2  9.00
#> [187,]  1  8.00
#> [188,]  5 15.00
#> [189,]  3 11.00
#> [190,]  1  9.00
#> [191,]  4 13.00
#> [192,] 14 29.00
#> [193,]  5 15.00
#> [194,]  2  8.05
#> [195,] 25 45.00
#> [196,]  2  9.00
#> [197,]  2  9.00
#> [198,]  1  8.00
#> [199,]  1  8.00
#> [200,]  1  8.00
#> [201,]  2 10.00
#> [202,]  6 16.00
#> [203,]  1  8.00
#> [204,]  1  8.00
#> [205,]  4 14.00
#> [206,] 10 24.00
#> [207,]  3 12.00
#> [208,]  2  9.00
#> [209,]  1  8.00
#> [210,]  2  9.00
#> [211,]  2 10.00
#> [212,]  3 11.00
#> [213,]  1  8.00
#> [214,]  1  8.00
#> [215,]  8 20.00
#> [216,]  4 13.00
#> [217,]  1  7.00
#> [218,]  2  9.00
#> [219,]  1  8.00
#> [220,]  4 14.00
#> [221,]  3 12.00
#> [222,]  3 12.00
#> [223,]  1  7.00
#> [224,]  3 12.00
#> [225,]  1  8.00
#> [226,] 66 98.00
#> [227,]  2  9.00
#> [228,]  4 13.00
#> [229,]  3 11.00
#> [230,]  6 17.00
#> [231,]  2 10.00
#> [232,]  1  8.00
#> [233,]  2  9.00
#> [234,]  2 10.00
#> [235,]  1  8.00
#> [236,]  1  8.00
# }
```
