set.seed(1234)
ndraws <- 5
nobs <- 4
ncat <- 3
x_test <- array(rnorm(ndraws * nobs * (ncat - 1)),
                dim = c(ndraws, nobs, ncat - 1))
nx_test <- -x_test
exp_nx_cumprod <- aperm(apply(exp(nx_test), c(1, 2), cumprod),
                        perm = c(2, 3, 1))
