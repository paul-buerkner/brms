# Bayesian Regression Models using 'Stan'

*Stan Development Team*

The brms package provides an interface to fit Bayesian generalized
multivariate (non-)linear multilevel models using Stan, which is a C++
package for obtaining full Bayesian inference (see
<https://mc-stan.org/>). The formula syntax is an extended version of
the syntax applied in the lme4 package to provide a familiar and simple
interface for performing regression analyses.

## Details

The main function of brms is
[`brm`](https://paulbuerkner.com/brms/reference/brm.md), which uses
formula syntax to specify a wide range of complex Bayesian models (see
[`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md)
for details). Based on the supplied formulas, data, and additional
information, it writes the Stan code on the fly via
[`stancode`](https://paulbuerkner.com/brms/reference/stancode.default.md),
prepares the data via
[`standata`](https://paulbuerkner.com/brms/reference/standata.default.md)
and fits the model using
[Stan](https://mc-stan.org/rstan/reference/rstan.html).

Subsequently, a large number of post-processing methods can be applied:
To get an overview on the estimated parameters,
[`summary`](https://paulbuerkner.com/brms/reference/summary.brmsfit.md)
or
[`conditional_effects`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.md)
are perfectly suited. Detailed visual analyses can be performed by
applying the
[`pp_check`](https://paulbuerkner.com/brms/reference/pp_check.brmsfit.md)
and
[`stanplot`](https://paulbuerkner.com/brms/reference/mcmc_plot.brmsfit.md)
methods, which both rely on the
[bayesplot](https://mc-stan.org/bayesplot/reference/bayesplot-package.html)
package. Model comparisons can be done via
[`loo`](https://paulbuerkner.com/brms/reference/loo.brmsfit.md) and
[`waic`](https://paulbuerkner.com/brms/reference/waic.brmsfit.md), which
make use of the
[loo](https://mc-stan.org/loo/reference/loo-package.html) package as
well as via
[`bayes_factor`](https://paulbuerkner.com/brms/reference/bayes_factor.brmsfit.md)
which relies on the bridgesampling package. For a full list of methods
to apply, type `methods(class = "brmsfit")`.

Because brms is based on Stan, a C++ compiler is required. The program
Rtools (available on <https://cran.r-project.org/bin/windows/Rtools/>)
comes with a C++ compiler for Windows. On Mac, you should use Xcode. For
further instructions on how to get the compilers running, see the
prerequisites section at the
[RStan-Getting-Started](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started)
page.

When comparing other packages fitting multilevel models to brms, keep in
mind that the latter needs to compile models before actually fitting
them, which will require between 20 and 40 seconds depending on your
machine, operating system and overall model complexity.

Thus, fitting smaller models may be relatively slow as compilation time
makes up the majority of the whole running time. For larger / more
complex models however, fitting my take several minutes or even hours,
so that the compilation time won't make much of a difference for these
models.

See `vignette("brms_overview")` and `vignette("brms_multilevel")` for a
general introduction and overview of brms. For a full list of available
vignettes, type `vignette(package = "brms")`.

## References

Paul-Christian Buerkner (2017). brms: An R Package for Bayesian
Multilevel Models Using Stan. *Journal of Statistical Software*, 80(1),
1-28. `doi:10.18637/jss.v080.i01`

Paul-Christian Buerkner (2018). Advanced Bayesian Multilevel Modeling
with the R Package brms. *The R Journal*. 10(1), 395–411.
`doi:10.32614/RJ-2018-017`

The Stan Development Team. *Stan Modeling Language User's Guide and
Reference Manual*. <https://mc-stan.org/users/documentation/>.

Stan Development Team (2020). RStan: the R interface to Stan. R package
version 2.21.2. <https://mc-stan.org/>

## See also

[`brm`](https://paulbuerkner.com/brms/reference/brm.md),
[`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.md),
[`brmsfamily`](https://paulbuerkner.com/brms/reference/brmsfamily.md),
[`brmsfit`](https://paulbuerkner.com/brms/reference/brmsfit-class.md)

## Author

**Maintainer**: Paul-Christian Bürkner <paul.buerkner@gmail.com>

Other contributors:

- Jonah Gabry \[contributor\]

- Sebastian Weber \[contributor\]

- Andrew Johnson \[contributor\]

- Martin Modrak \[contributor\]

- Hamada S. Badr \[contributor\]

- Frank Weber \[contributor\]

- Aki Vehtari \[contributor\]

- Mattan S. Ben-Shachar \[contributor\]

- Hayden Rabel \[contributor\]

- Simon C. Mills \[contributor\]

- Stephen Wild \[contributor\]

- Ven Popov \[contributor\]

- Ioannis Kosmidis \[contributor\]

- Ben Schneider \[contributor\]

- Noa Kallioinen \[contributor\]

- Wellington J. Silva \[contributor\]

- Luiz Carvalho \[contributor\]

- Sermet Pekin \[contributor\]

- Florence Bockting \[contributor\]
