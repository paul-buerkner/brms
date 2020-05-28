[![Build Status](https://travis-ci.org/paul-buerkner/brms.svg?branch=master)](https://travis-ci.org/paul-buerkner/brms)[![CRAN Version](http://www.r-pkg.org/badges/version/brms)](https://cran.r-project.org/package=brms)[![Coverage Status](https://codecov.io/github/paul-buerkner/brms/coverage.svg?branch=master)](https://codecov.io/github/paul-buerkner/brms?branch=master)


<br><br><br>

<div style="text-align:left; padding:-40px;">
<a href="http://mc-stan.org">
<img src="https://raw.githubusercontent.com/stan-dev/logos/master/logo_tm.png" align="right" width=100 alt="Stan Logo"/>
<img src="https://raw.githubusercontent.com/paul-buerkner/brms/master/man/figures/brms.png" width = 100 alt="brms Logo"/>
</a><h2><strong>brms</strong></h2>
<h4>Bayesian regression models using Stan</h4>
</div>

****

The **brms** package provides an interface to fit Bayesian generalized
(non-)linear multivariate multilevel models using Stan. The formula syntax is
very similar to that of the package lme4 to provide a familiar and simple
interface for performing regression analyses.

A wide range of distributions and link functions are supported, allowing users
to fit -- among others -- linear, robust linear, count data, survival, response
times, ordinal, zero-inflated, hurdle, and even self-defined mixture models all
in a multilevel context. Further modeling options include non-linear and smooth
terms, auto-correlation structures, censored data, meta-analytic standard
errors, and quite a few more. In addition, all parameters of the response
distribution can be predicted in order to perform distributional regression.
Prior specifications are flexible and explicitly encourage users to apply prior
distributions that actually reflect their beliefs. Model fit can easily be
assessed and compared with posterior predictive checks and leave-one-out
cross-validation.

## Getting Started

If you are new to **brms** we recommend starting with the [vignettes](https://paul-buerkner.github.io/brms/articles/) and these 
other resources:

* [Introduction to brms](https://www.jstatsoft.org/article/view/v080i01) 
(Journal of Statistical Software)
* [Advanced multilevel modeling with brms](https://journal.r-project.org/archive/2018/RJ-2018-017/index.html) 
(The R Journal)
* [Blog posts](https://paul-buerkner.github.io/blog/brms-blogposts/) 
(List of blog posts about brms)
* [Ask a question](http://discourse.mc-stan.org/) (Stan Forums on Discourse)


## Installation


Install the latest release from **CRAN**

```r
install.packages("brms")
```

Install the latest development version from **GitHub**

```r
if (!require(devtools)) {
  install.packages("devtools")
}
devtools::install_github("paul-buerkner/brms", build_vignettes = FALSE)
```

You can also set `build_vignettes=TRUE` but this will slow down the installation
drastically (the vignettes can always be accessed online anytime at
[paul-buerkner.github.io/brms/articles](https://paul-buerkner.github.io/brms/articles/)).


