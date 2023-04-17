test_that("Test various functions", {
  v <- rnorm(100, mean=2, sd=5)
  v <- stdVar(v)
  expect_equal(0, mean(v))
  expect_equal(1, var(v))
})
