# toy objects reused across tests
f1 <- mpg ~ wt
f2 <- mpg ~ wt + cyl            # same data, different formula
d1 <- mtcars
d2 <- mtcars[sample(nrow(mtcars)), ]   # same rows, shuffled order
fam1 <- gaussian()
fam2 <- student()
base <- nlist(formula = f1, data = d1, family = fam1)
h_base <- hash_dots(base, data_policy = "hash")

test_that("hash function check" , {
  aa1 =   hash_dots( count ~  zAge + zBase * Trt + (1|patient),
                         data = epilepsy, family = poisson() ,    file_auto = TRUE )
  aa2 =   hash_dots( count ~  zAge + zBase * Trt + (1|patient),
                         data = epilepsy, family = poisson() ,    file_auto = TRUE )
  expect_equal( aa1 , aa2 )
  g1 =  gaussian()
  g2 = gaussian()
  expect_true(hash_dots(g1) == hash_dots(g2))
  expect_true(hash_dots(gaussian()) == hash_dots(gaussian()))
  expect_false(hash_dots(poisson()) == hash_dots(gaussian()))
  expect_true(hash_dots(epilepsy) == hash_dots(epilepsy))
  expect_false(hash_dots(epilepsy[-c(2) , ]) == hash_dots(epilepsy))
})

test_that("identical arguments give identical hash", {
  a <- nlist(formula = f1, data = d1, family = fam1)
  h1 <- hash_dots(a, data_policy = "hash")
  h2 <- hash_dots(a, data_policy = "hash")
  expect_identical(h1, h2)
})

test_that("ordering of list elements is irrelevant", {
  a <- nlist(formula = f1, data = d1, family = fam1)
  b <- a[c("family", "data", "formula")]   # reorder
  h1 <- hash_dots(a, data_policy = "hash")
  h2 <- hash_dots(b, data_policy = "hash")
  expect_identical(h1, h2)
})

test_that("formula environment ignored", {
  env <- new.env()
  environment(f1) <- env
  a <- nlist(formula = f1, data = d1, family = fam1)
  b <- nlist(formula = mpg ~ wt, data = d1, family = fam1)
  h1 <- hash_dots(a, data_policy = "hash")
  h2 <- hash_dots(b, data_policy = "hash")
  expect_identical(h1, h2)
})

test_that("changing *one* statistical input changes the hash", {

  # 1. Different formula
  h_formula <- hash_dots(
    modifyList(base, list(formula = f2)), data_policy = "hash"
  )
  expect_false(identical(h_base, h_formula))
  # 2. Different family
  h_family <- hash_dots(
    modifyList(base, list(family = fam2)), data_policy = "hash"
  )
  expect_false(identical(h_base, h_family))
  # 3. Different data (row order → different digest with data_policy = "hash")
  h_data <- hash_dots(
    modifyList(base, list(data = d2)), data_policy = "hash"
  )
  expect_false(identical(h_base, h_data))
})


# -------------------------------------------------------------------------
# Data: just something to pass; content is irrelevant for hashing logic
d <- mtcars

# Helper: a small convenience to shorten calls
h <- function(...) hash_dots(..., data = d, family = gaussian(), data_policy = "names")

context("hash_dots() – complex formula structures")

# ─────────────────────────────────────────────────────────────────────────
test_that("Hidden environments in nested brmsformula are ignored", {
  # Non-linear single-response model
  bf1 <- bf(y ~ mu, mu ~ wt + (1|cyl), nl = TRUE)
  bf2 <- bf1
  # Attach a random environment to *internal* formula
  environment(bf2$formula) <- new.env(parent = emptyenv())

  expect_identical(h(formula = bf1),
                   h(formula = bf2))
})



# ─────────────────────────────────────────────────────────────────────────
test_that("Changing an internal sub-formula *does* change the hash", {
  base <- bf(y ~ mu, mu ~ wt, nl = TRUE)
  changed <- bf(y ~ mu, mu ~ wt + hp, nl = TRUE)  # extra predictor

  expect_false(identical(h(formula = base),
                         h(formula = changed)))
})

# ─────────────────────────────────────────────────────────────────────────
test_that("Changing only the prior alters the hash", {
  f <- bf(mpg ~ wt)

  h_base <- h(formula = f,
              prior   = NULL)

  h_prior <- h(formula = f,
               prior   = prior(normal(0, 5), class = "b"))

  expect_false(identical(h_base, h_prior))
})

# ─────────────────────────────────────────────────────────────────────────
test_that("List argument order is normalised", {
  f <- bf(mpg ~ wt)

  h1 <- hash_dots(formula = f, data = d, family = gaussian(), iter = 2000)
  h2 <- hash_dots(iter = 2000, family = gaussian(), data = d, formula = f)

  expect_identical(h1, h2)
})


