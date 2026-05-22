## print/summary contracts:
##   - print(x) returns x invisibly
##   - summary(x) returns an S3 `summary.<class>` object
##   - print(summary(x)) returns its argument invisibly
##   - broom::tidy/glance dispatch correctly when broom is loaded

fit_for <- function() {
  fx <- make_fixture(n = 400)
  rdhte(y = fx$y, x = fx$x, covs.hte = fx$W)
}

bw_for <- function() {
  fx <- make_fixture(n = 400)
  rdbwhte(y = fx$y, x = fx$x, covs.hte = fx$W)
}

test_that("print methods return their input invisibly", {
  m  <- fit_for()
  bw <- bw_for()

  silence <- file(tempfile(), open = "wt")
  sink(silence)
  on.exit({ sink(); close(silence) }, add = TRUE)

  for (obj in list(m, bw)) {
    res <- withVisible(print(obj))
    expect_false(res$visible)
    expect_identical(res$value, obj)
  }
})

test_that("summary returns an S3 summary.<class> object", {
  m  <- fit_for()
  bw <- bw_for()

  capture.output({
    sm  <- summary(m)
    sbw <- summary(bw)
  })

  expect_s3_class(sm,  "summary.rdhte")
  expect_s3_class(sbw, "summary.rdbwhte")
  ## Inheritance preserved so other generics still dispatch
  expect_s3_class(sm,  "rdhte")
  expect_s3_class(sbw, "rdbwhte")
})

test_that("broom tidy/glance dispatch for rdhte and rdbwhte", {
  testthat::skip_if_not_installed("broom")
  m  <- fit_for()
  bw <- bw_for()

  td_m <- broom::tidy(m)
  expect_s3_class(td_m, "data.frame")
  expect_true(all(c("term","estimate","std.error","statistic","p.value",
                    "conf.low","conf.high") %in% names(td_m)))
  expect_equal(nrow(td_m), length(unique(make_fixture(n = 400)$W)))

  gl_m <- broom::glance(m)
  expect_equal(nrow(gl_m), 1L)
  expect_true(all(c("n","kernel","vce","bwselect") %in% names(gl_m)))

  td_bw <- broom::tidy(bw)
  expect_s3_class(td_bw, "data.frame")
  expect_true(all(c("term","h.left","h.right") %in% names(td_bw)))

  gl_bw <- broom::glance(bw)
  expect_equal(nrow(gl_bw), 1L)
})
