context("Tests for data helper functions")

test_that("validate_newdata handles factors correctly", {
  fit <- brms:::rename_pars(brms:::brmsfit_example1)
  fit$data$fac <- factor(sample(1:3, nrow(fit$data), TRUE))
  newdata <- fit$data[1:5, ]
  expect_silent(brms:::validate_newdata(newdata, fit))
  newdata$visit <- 1:5
  expect_error(brms:::validate_newdata(newdata, fit),
               "Levels '5' of grouping factor 'visit' cannot")
  newdata$fac <- 1:5
  expect_error(brms:::validate_newdata(newdata, fit),
               "New factor levels are not allowed")
})

test_that("validate_newdata handles factor matrix columns", {
  skip_if_not_installed("mgcv")
  G <- 4L
  nodes <- paste0("n", seq_len(G))
  # each observation carries a factor matrix of nodes and matching weights
  make_data <- function(weights) {
    out <- data.frame(y = as.numeric(weights %*% seq(0.1, 0.4, by = 0.1)))
    node <- factor(rep(nodes, each = nrow(weights)), levels = nodes)
    dim(node) <- dim(weights)
    out$node <- node
    out$weight <- weights
    out
  }
  dat <- make_data(rbind(diag(G), diag(G)))
  # smooth terms are evaluated in the environment of the brms namespace, so
  # the path-graph Laplacian serving as MRF penalty is specified inline
  # (its row/column order matches the order of the node levels)
  fit <- brm(
    y ~ 0 + s(node, by = weight, bs = "mrf", k = 4,
              xt = list(penalty = crossprod(diff(diag(4))))),
    data = dat, empty = TRUE
  )

  # the matrix structure of the factor must survive newdata validation
  newdata <- brms:::validate_newdata(dat, fit)
  expect_equal(dim(newdata$node), dim(dat$node))
  expect_equal(levels(newdata$node), nodes)
  sdata <- standata(fit)
  expect_equal(standata(fit, newdata = dat), sdata)

  # genuinely new data must be evaluated with the basis of the fitted smooth
  new_weights <- rbind(c(0.5, 0.5, 0, 0), c(0, 0, 1, 0), rep(0.25, G))
  sdata_new <- standata(fit, newdata = make_data(new_weights))
  # mgcv sums the basis over the columns of the matrix predictors, so the
  # basis of new rows is the weighted average of the single-node basis rows
  expect_equal(sdata_new$Zs_1_1, new_weights %*% sdata$Zs_1_1[seq_len(G), ],
               check.attributes = FALSE)

  # subsets of the original data reproduce the original basis rows
  sdata_sub <- standata(fit, newdata = dat[1:3, ])
  expect_equal(sdata_sub$Zs_1_1, sdata$Zs_1_1[1:3, ], check.attributes = FALSE)

  # a single observation using only one node is valid as well
  sdata_one <- standata(fit, newdata = make_data(rbind(c(1, 0, 0, 0))))
  expect_equal(sdata_one$Zs_1_1, sdata$Zs_1_1[1, , drop = FALSE],
               check.attributes = FALSE)

  # node levels may also be passed as a plain character matrix
  chardata <- make_data(new_weights)
  chardata$node <- matrix(as.character(chardata$node), nrow(new_weights), G)
  expect_equal(standata(fit, newdata = chardata)$Zs_1_1, sdata_new$Zs_1_1,
               check.attributes = FALSE)
})

test_that("validate_data returns correct model.frames", {
  dat <- data.frame(y = 1:5, x = 1:5, z = 6:10, g = 5:1)

  bterms <- brmsterms(y ~ as.numeric(x) + (as.factor(z) | g),
                      family = gaussian())
  mf <- brms:::validate_data(dat, bterms = bterms)
  expect_true(all(c("x", "z") %in% names(mf)))

  bterms <- brmsterms(y ~ 1 + (1|g/x/z), family = gaussian())
  mf <- brms:::validate_data(dat, bterms = bterms)
  expect_equal(mf[["g:x"]], paste0(dat$g, "_", dat$x))
  expect_equal(mf[["g:x:z"]], paste0(dat$g, "_", dat$x, "_", dat$z))
})

