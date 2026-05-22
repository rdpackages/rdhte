## plot.rdhte tests

skip_if_no_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    testthat::skip("ggplot2 not installed")
  }
}

make_plot_fixture <- function(n = 600, seed = 20260510) {
  set.seed(seed)
  x <- runif(n, -1, 1)
  W <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  y <- 0.5 +
       (x >= 0) * (W == "A") * 1.0 +
       (x >= 0) * (W == "B") * 2.0 +
       (x >= 0) * (W == "C") * 0.3 +
       rnorm(n, sd = 0.4)
  list(y = y, x = x, W = W)
}

test_that("plot.rdhte returns a ggplot object for categorical W", {
  skip_if_no_ggplot2()
  d <- make_plot_fixture()
  m <- rdhte(y = d$y, x = d$x, covs.hte = d$W)
  ## suppress the side-effect print() so the test runner stays quiet.
  p <- suppressMessages(plot(m))
  expect_s3_class(p, "ggplot")
})

test_that("plot data has one row per group with the right columns", {
  skip_if_no_ggplot2()
  d <- make_plot_fixture()
  m <- rdhte(y = d$y, x = d$x, covs.hte = d$W)
  p <- suppressMessages(plot(m))
  expect_equal(nrow(p$data), length(m$W.lev))
  expect_true(all(c("group", "estimate", "ci_low", "ci_high") %in% names(p$data)))
  ## point estimates come from m$Estimate, CIs from m$ci.rb
  ord <- order(as.character(p$data$group))
  expect_equal(p$data$estimate[ord],
               as.numeric(m$Estimate)[match(p$data$group[ord], m$W.lev)],
               tolerance = 1e-12)
  expect_equal(p$data$ci_low[ord],
               as.numeric(m$ci.rb[, 1])[match(p$data$group[ord], m$W.lev)],
               tolerance = 1e-12)
  expect_equal(p$data$ci_high[ord],
               as.numeric(m$ci.rb[, 2])[match(p$data$group[ord], m$W.lev)],
               tolerance = 1e-12)
})

test_that("sort = TRUE reorders the x-axis factor levels by effect", {
  skip_if_no_ggplot2()
  d <- make_plot_fixture()
  m <- rdhte(y = d$y, x = d$x, covs.hte = d$W)
  p_sorted <- suppressMessages(plot(m, sort = TRUE))
  est <- p_sorted$data$estimate
  expect_equal(est, sort(est))
  expect_equal(levels(p_sorted$data$group), as.character(p_sorted$data$group))
})

test_that("plot.rdhte errors on continuous covs.hte", {
  skip_if_no_ggplot2()
  set.seed(11); n <- 400
  x <- runif(n, -1, 1)
  z <- runif(n)
  y <- (x >= 0) * (1 + 0.5 * z) + rnorm(n)
  m <- rdhte(y = y, x = x, covs.hte = z)
  expect_error(plot(m), "categorical")
})

test_that("plot.rdhte errors on rdhte fit with no covs.hte", {
  skip_if_no_ggplot2()
  set.seed(11); n <- 400
  x <- runif(n, -1, 1)
  y <- (x >= 0) + rnorm(n)
  m <- rdhte(y = y, x = x)
  expect_error(plot(m), "categorical")
})
