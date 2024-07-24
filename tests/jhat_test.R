rm(list=ls())
library(Rcpp)
library(RcppArmadillo)

# Source R implementation
source("~/Documents/uni/master-dissertation/code-cont/chen/jhat.R")

# Source C++ implementation 
Rcpp::sourceCpp("chen/jhat.cpp")


# Test function
test_jhat <- function(n_tests = 100, max_dim = 100, max_nL = 10) {
  
  set.seed(42)  # For reproducibility
  d = 0
  # browser()
  
  for (i in 1:n_tests) {
    # cat("Test", i, "\n")
    
    # Generate random inputs
    n <- sample(50:1000, 1)
    nL <- sample(2:max_nL, 1)
    M <- runif(1, 0.5, 2)
    
    # Generate PP and BB matrices
    dim_p <- c(sample(10:max_dim, 1), sample(10:max_dim, 1))
    dim_b <- c(sample(10:max_dim, 1), dim_p[2])
    PP <- matrix(rnorm(prod(dim_p)), nrow = dim_p[1])
    BB <- matrix(rnorm(prod(dim_b)), nrow = dim_b[1])
    
    # Generate CJ, CK, and TJ
    CJ <- sort(c(0, sample(1:(dim_p[2]-1), nL, replace = FALSE), dim_p[2]))
    CK <- sort(c(0, sample(1:(dim_b[2]-1), nL, replace = FALSE), dim_b[2]))
    TJ <- sort(runif(nL+1, 1, 10))
    
    
    # Run R implementation
    # browser()
    if (i == 20){
      # browser()
      r_result <- tryCatch({
        # print("r:")
        jhat_r(PP, BB, CJ, CK, TJ, M, n, nL)
      }, error = function(e) {
        cat("R implementation error:", conditionMessage(e), "\n")
        return(NULL)
      })
    
      # Run C++ implementation
      cpp_result <- tryCatch({
        # print("cpp:")
        jhat(PP, BB, CJ, CK, TJ, M, n, nL)  # Adjust for 0-based indexing in C++
      }, error = function(e) {
        cat("C++ implementation error:", conditionMessage(e), "\n")
        return(NULL)
      })
    
    
    # Compare results
    if (!is.null(r_result) && !is.null(cpp_result)) {
      if (r_result$LL != cpp_result$LL || r_result$flag != cpp_result$flag) {
        cat("Results differ:\n")
        print(i)
        cat("R result:  LL =", r_result$LL, ", flag =", r_result$flag, "\n")
        cat("C++ result: LL =", cpp_result$LL, ", flag =", cpp_result$flag, "\n")
      } else {
        
        d = d + 1
      }
    }
    
    # cat("\n")
    }
  }
  if (d >0) {
    # print(d)
  }
}

# Run tests
test_jhat()

# # Additional edge case tests
# test_edge_cases <- function() {
#   cat("Testing edge cases:\n")
#   
#   # Test case 1: All zeros
#   PP <- matrix(0, 10, 10)
#   BB <- matrix(0, 10, 10)
#   CJ <- c(0, 5, 10)
#   CK <- c(0, 5, 10)
#   TJ <- c(1, 5, 10)
#   M <- 1
#   n <- 100
#   nL <- 2
#   
#   r_result <- Jhat(PP, BB, CJ, CK, TJ, M, n, nL)
#   cpp_result <- jhat(PP, BB, CJ-1, CK-1, TJ, M, n, nL)
#   
#   cat("All zeros case:\n")
#   cat("R result:  ", r_result, "\n")
#   cat("C++ result:", cpp_result, "\n\n")
#   
#   # Test case 2: Very large values
#   PP <- matrix(1e10, 10, 10)
#   BB <- matrix(1e10, 10, 10)
#   
#   r_result <- Jhat(PP, BB, CJ, CK, TJ, M, n, nL)
#   cpp_result <- jhat(PP, BB, CJ-1, CK-1, TJ, M, n, nL)
#   
#   cat("Very large values case:\n")
#   cat("R result:  ", r_result, "\n")
#   cat("C++ result:", cpp_result, "\n\n")
#   
#   # Add more edge cases as needed
# }
# 
# # Run edge case tests
# test_edge_cases()