test_that("performance with large matrices", {
  dim_size <- 100  # Testing with a 100x100 matrix
  P <- matrix(rnorm(dim_size^2), nrow = dim_size)
  B <- matrix(rnorm(dim_size^2), nrow = dim_size)
  
  start_time <- Sys.time()
  result <- shat_r(P, B)
  end_time <- Sys.time()
  
  print(end_time - start_time)  # Print the time taken to execute
  expect_true(result > 0)  # Basic check to ensure function executes
})
