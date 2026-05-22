## When covs.hte is a binary factor, rdhte's per-subgroup tau must match
## a plain rdrobust() fit on the same subgroup to machine precision.
## Mirrors the long-standing checkup_rdhte_vs_rdrobust.R shape but pinned
## here so R CMD check exercises it.

test_that("rdhte per-subgroup tau matches rdrobust on the subgroup", {
  testthat::skip_if_not_installed("rdrobust")

  set.seed(20260509)
  n  <- 800
  x  <- runif(n, -1, 1)
  W  <- rbinom(n, 1, 0.5)
  y  <- 0.5 + 1.0 * (x >= 0) + 0.5 * x + 2 * W * (x >= 0) + rnorm(n, sd = 0.3)
  Wf <- factor(W)

  m <- rdhte(y = y, x = x, covs.hte = Wf)

  for (lev in c(0, 1)) {
    keep <- W == lev
    rd   <- rdrobust::rdrobust(y = y[keep], x = x[keep],
                               h = m$h[lev + 1, 1],
                               b = m$h[lev + 1, 1])
    ## Conventional point estimate
    expect_equal(as.numeric(m$Estimate[lev + 1]),
                 as.numeric(rd$coef["Conventional", 1]),
                 tolerance = 1e-10,
                 label = sprintf("subgroup %d Conventional", lev))
    ## Bias-corrected
    expect_equal(as.numeric(m$Estimate.bc[lev + 1]),
                 as.numeric(rd$coef["Bias-Corrected", 1]),
                 tolerance = 1e-10,
                 label = sprintf("subgroup %d Bias-Corrected", lev))
  }
})
