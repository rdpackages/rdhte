## vce surface for cluster:
##   - default with cluster -> CR1 (was: silent HC1 alias).
##   - cr1/cr2/cr3 are canonical cluster-vce names.
##   - hc{0,1,2,3} + cluster: warn + remap to cr1/cr1/cr2/cr3.
##   - cr1/cr2/cr3 without cluster: warn + fall back to hc1/hc2/hc3.

make_fixture <- function(n = 600, seed = 20260509) {
  set.seed(seed)
  x  <- runif(n, -1, 1)
  W  <- factor(rbinom(n, 1, 0.5))
  cl <- ceiling(20 * runif(n))
  y  <- 0.3 + 1.0 * (x >= 0) + 0.5 * x + 1.5 * (W == "1") * (x >= 0) + rnorm(n, sd = 0.4)
  list(y = y, x = x, W = W, cl = cl)
}

test_that("default with cluster is CR1 (silent)", {
  d <- make_fixture()
  expect_silent(
    m_default <- rdhte(y = d$y, x = d$x, covs.hte = d$W, cluster = d$cl)
  )
  m_cr1 <- rdhte(y = d$y, x = d$x, covs.hte = d$W, cluster = d$cl, vce = "cr1")
  expect_equal(m_default$se.rb, m_cr1$se.rb, tolerance = 1e-12)
  expect_identical(m_default$vce, "CR1")
  expect_identical(m_cr1$vce, "CR1")
})

test_that("hc{0,1,2,3} + cluster warns and remaps to cr1/cr1/cr2/cr3", {
  d <- make_fixture()
  expect_warning(
    m_hc0 <- rdhte(y = d$y, x = d$x, covs.hte = d$W, cluster = d$cl, vce = "hc0"),
    "Switching to vce='cr1'"
  )
  expect_warning(
    m_hc1 <- rdhte(y = d$y, x = d$x, covs.hte = d$W, cluster = d$cl, vce = "hc1"),
    "Switching to vce='cr1'"
  )
  expect_warning(
    m_hc2 <- rdhte(y = d$y, x = d$x, covs.hte = d$W, cluster = d$cl, vce = "hc2"),
    "Switching to vce='cr2'"
  )
  expect_warning(
    m_hc3 <- rdhte(y = d$y, x = d$x, covs.hte = d$W, cluster = d$cl, vce = "hc3"),
    "Switching to vce='cr3'"
  )
  m_cr1 <- rdhte(y = d$y, x = d$x, covs.hte = d$W, cluster = d$cl, vce = "cr1")
  m_cr2 <- rdhte(y = d$y, x = d$x, covs.hte = d$W, cluster = d$cl, vce = "cr2")
  m_cr3 <- rdhte(y = d$y, x = d$x, covs.hte = d$W, cluster = d$cl, vce = "cr3")
  expect_equal(m_hc0$se.rb, m_cr1$se.rb, tolerance = 1e-12)
  expect_equal(m_hc1$se.rb, m_cr1$se.rb, tolerance = 1e-12)
  expect_equal(m_hc2$se.rb, m_cr2$se.rb, tolerance = 1e-12)
  expect_equal(m_hc3$se.rb, m_cr3$se.rb, tolerance = 1e-12)
  expect_identical(m_hc0$vce, "CR1")
  expect_identical(m_hc1$vce, "CR1")
  expect_identical(m_hc2$vce, "CR2")
  expect_identical(m_hc3$vce, "CR3")
})

test_that("CR1, CR2, CR3 produce different SEs (sanity)", {
  d <- make_fixture()
  m1 <- rdhte(y = d$y, x = d$x, covs.hte = d$W, cluster = d$cl, vce = "cr1")
  m2 <- rdhte(y = d$y, x = d$x, covs.hte = d$W, cluster = d$cl, vce = "cr2")
  m3 <- rdhte(y = d$y, x = d$x, covs.hte = d$W, cluster = d$cl, vce = "cr3")
  expect_false(isTRUE(all.equal(as.numeric(m1$se.rb), as.numeric(m2$se.rb))))
  expect_false(isTRUE(all.equal(as.numeric(m1$se.rb), as.numeric(m3$se.rb))))
  ## CR3 should be >= CR1 elementwise (leverage-adjusted residuals).
  expect_true(all(as.numeric(m3$se.rb) >= as.numeric(m1$se.rb)))
})

test_that("cr1/cr2/cr3 without cluster warn + fall back to hc1/hc2/hc3", {
  d <- make_fixture()
  expect_warning(
    m_cr1 <- rdhte(y = d$y, x = d$x, covs.hte = d$W, vce = "cr1"),
    "requires a cluster variable"
  )
  expect_warning(
    m_cr2 <- rdhte(y = d$y, x = d$x, covs.hte = d$W, vce = "cr2"),
    "requires a cluster variable"
  )
  expect_warning(
    m_cr3 <- rdhte(y = d$y, x = d$x, covs.hte = d$W, vce = "cr3"),
    "requires a cluster variable"
  )
  m_hc1 <- rdhte(y = d$y, x = d$x, covs.hte = d$W, vce = "hc1")
  m_hc2 <- rdhte(y = d$y, x = d$x, covs.hte = d$W, vce = "hc2")
  m_hc3 <- rdhte(y = d$y, x = d$x, covs.hte = d$W, vce = "hc3")
  expect_equal(m_cr1$se.rb, m_hc1$se.rb, tolerance = 1e-12)
  expect_equal(m_cr2$se.rb, m_hc2$se.rb, tolerance = 1e-12)
  expect_equal(m_cr3$se.rb, m_hc3$se.rb, tolerance = 1e-12)
  expect_identical(m_cr1$vce, "HC1")
  expect_identical(m_cr2$vce, "HC2")
  expect_identical(m_cr3$vce, "HC3")
})

test_that("non-cluster default still HC3", {
  d <- make_fixture()
  m <- rdhte(y = d$y, x = d$x, covs.hte = d$W)
  expect_identical(m$vce, "HC3")
})

test_that("invalid vce strings are rejected", {
  d <- make_fixture(n = 200)
  expect_error(rdhte(y = d$y, x = d$x, covs.hte = d$W, vce = "robust"), "vce")
  expect_error(rdhte(y = d$y, x = d$x, covs.hte = d$W, vce = c("hc1", "hc3")), "vce")
  expect_error(rdhte(y = d$y, x = d$x, covs.hte = d$W, vce = 1L), "vce")
})

test_that("rdbwhte mirrors the cluster-aware default + cr aliases", {
  d <- make_fixture()
  expect_silent(
    bw_default <- rdbwhte(y = d$y, x = d$x, covs.hte = d$W, cluster = d$cl)
  )
  bw_cr1 <- rdbwhte(y = d$y, x = d$x, covs.hte = d$W, cluster = d$cl, vce = "cr1")
  expect_equal(bw_default$h, bw_cr1$h, tolerance = 1e-12)
  expect_identical(bw_default$vce, "CR1")
  expect_identical(bw_cr1$vce, "CR1")

  expect_warning(
    bw_hc1 <- rdbwhte(y = d$y, x = d$x, covs.hte = d$W, cluster = d$cl, vce = "hc1"),
    "Switching to vce='cr1'"
  )
  expect_equal(bw_hc1$h, bw_cr1$h, tolerance = 1e-12)
})
