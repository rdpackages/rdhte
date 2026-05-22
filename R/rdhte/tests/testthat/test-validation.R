## Input validation hardening (ported from checkup_rdhte_validation.R,
## 2026-05-08). Each block mirrors the rdrobust/nprobust validation
## pattern (level<100 strict, vector-length checks, c/p/q validation).

test_that("c is a single finite numeric", {
  fx <- make_fixture()
  expect_error(rdhte(fx$y, fx$x, c = NA,     covs.hte = fx$W), "Cutoff 'c'")
  expect_error(rdhte(fx$y, fx$x, c = Inf,    covs.hte = fx$W), "Cutoff 'c'")
  expect_error(rdhte(fx$y, fx$x, c = c(0,0), covs.hte = fx$W), "Cutoff 'c'")
})

test_that("p must be a non-negative integer", {
  fx <- make_fixture()
  expect_error(rdhte(fx$y, fx$x, p = -1,  covs.hte = fx$W), "Polynomial order 'p'")
  expect_error(rdhte(fx$y, fx$x, p = 1.5, covs.hte = fx$W), "Polynomial order 'p'")
})

test_that("q must satisfy q > p", {
  fx <- make_fixture()
  expect_error(rdhte(fx$y, fx$x, p = 1, q = 1, covs.hte = fx$W), "q > p")
  expect_error(rdhte(fx$y, fx$x, p = 2, q = 1, covs.hte = fx$W), "q > p")
})

test_that("level must be in (0, 100) and finite", {
  fx <- make_fixture()
  expect_error(rdhte(fx$y, fx$x, covs.hte = fx$W, level = 100), "level")
  expect_error(rdhte(fx$y, fx$x, covs.hte = fx$W, level = 0),   "level")
  expect_error(rdhte(fx$y, fx$x, covs.hte = fx$W, level = NA),  "level")
})

test_that("auxiliary-vector lengths must match length(x)", {
  fx <- make_fixture()
  expect_error(rdhte(fx$y[1:10], fx$x, covs.hte = fx$W),
               "must have equal length")
  expect_error(rdhte(fx$y, fx$x, covs.hte = fx$W, weights = rep(1, 10)),
               "weights")
  expect_error(rdhte(fx$y, fx$x, covs.hte = fx$W, cluster = fx$cl[1:10]),
               "cluster")
  expect_error(rdhte(fx$y, fx$x, covs.hte = fx$W, covs.eff = fx$zw[1:10]),
               "covs.eff")
  expect_error(rdhte(fx$y, fx$x, covs.hte = fx$W[1:10]),
               "covs.hte")
})

test_that("subset bad shapes/values are rejected", {
  fx <- make_fixture()
  expect_error(rdhte(fx$y, fx$x, covs.hte = fx$W, subset = rep(TRUE, 10)),
               "subset")
  expect_error(rdhte(fx$y, fx$x, covs.hte = fx$W, subset = c(1, 999999)),
               "subset")
})

test_that("h / h.l / h.r must be positive and finite", {
  fx <- make_fixture()
  expect_error(rdhte(fx$y, fx$x, covs.hte = fx$W, h   = 0),  "'h' must be a positive")
  expect_error(rdhte(fx$y, fx$x, covs.hte = fx$W, h.l = -1), "h.l")
  expect_error(rdhte(fx$y, fx$x, covs.hte = fx$W, h.r = NA), "h.r")
})

test_that("rdbwhte mirrors the cluster / cutoff validation", {
  fx <- make_fixture()
  expect_error(rdbwhte(fx$y, fx$x, covs.hte = fx$W, cluster = fx$cl[1:10]),
               "cluster")
  expect_error(rdbwhte(fx$y, fx$x, covs.hte = fx$W, c = NA),
               "Cutoff 'c'")
})
