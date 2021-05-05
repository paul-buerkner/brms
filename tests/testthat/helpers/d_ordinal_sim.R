# This test corresponds to a single observation.
set.seed(1234)
ndraws <- 5
ncat <- 3
thres_test <- matrix(rnorm(ndraws * (ncat - 1)), nrow = ndraws)
# Emulate no category-specific effects (i.e., only a single vector of linear
# predictors) as well as category-specific effects (i.e., a matrix of linear
# predictors):
eta_test_list <- list(rnorm(ndraws),
                      matrix(rnorm(ndraws * (ncat - 1)), nrow = ndraws))
