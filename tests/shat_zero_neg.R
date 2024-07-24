library(testthat)

test_that("handle singular Gb matrix", {
  # Creating a singular matrix B
  B <- matrix(c(1, 0, 0, 0), nrow = 2)
  P <- matrix(rnorm(4), nrow = 2)
  result <- shat_r(P, B)
  expect_equal(result, 1e-20)
})
