## Shared fixture for rdhte test suite.
make_fixture <- function(n = 600, seed = 20260508) {
  set.seed(seed)
  x  <- runif(n, -1, 1)
  y  <- 0.5 * x + 1 * (x >= 0) + rnorm(n)
  W  <- factor(rbinom(n, 1, 0.5))
  cl <- ceiling(20 * runif(n))
  zw <- rnorm(n)
  list(x = x, y = y, W = W, cl = cl, zw = zw, n = n)
}

## Helper that captures error messages so we can assert on a substring.
## rdhte's validators use stop() with informative text; testthat's
## expect_error(regexp = ...) is the natural fit, but we also want to
## flag warnings as misses. capture_error_message returns the error text
## or NA if the expression succeeded.
capture_error_message <- function(expr) {
  tryCatch({ force(expr); NA_character_ },
           error = function(e) conditionMessage(e))
}
