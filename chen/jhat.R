# Source the required file and load necessary libraries
source("~/Documents/uni/master-dissertation/code-cont/chen/shat.R")
library(parallel)
library(doParallel)

# Jhat function
Jhat <- function(PP, BB, CJ, CK, TJ, M, n, nL) {
  lb <- numeric(nL + 1)
  ub <- numeric(nL + 1)
  
  # Create a cluster
  cl <- makeCluster(detectCores() - 2)
  
  # Export necessary variables and functions to the cluster
  clusterExport(cl, varlist = c("PP", "PP", "CJ", "CJ", "TJ", "M", "n", "shat"))
  clusterEvalQ(cl, {
    library(expm)
  })
  
  results <- parLapply(cl, 1:(nL + 1), function(ll) {
    s <- tryCatch({
      shat(PP[, (CJ[ll] + 1):(CJ[ll + 1])], PP[, (CJ[ll] + 1):(CJ[ll + 1])])
    }, error = function(e) {
      1e-20
    })
    
    J <- TJ[ll]
    lb_value <- J * sqrt(log(J)) * max(0 * (log(n))^4, 1 / s)
    return(lb_value)
  })
  
  stopCluster(cl) # Stop the cluster
  
  lb <- unlist(results)
  ub[1:nL] <- lb[2:(nL + 1)]
  
  L <- which(lb <= 2 * M * sqrt(n) & 2 * M * sqrt(n) <= ub) - 1
  f <- 0
  
  if (length(L) == 0) {
    L <- which(lb <= 2 * M * sqrt(n), arr.ind = TRUE)[length(which(lb <= 2 * M * sqrt(n)))] - 1
    f <- 1
  }
  
  if (length(L) == 0) {
    L <- 1
    f <- 2
  }
  
  LL <- max(L, 1)
  flag <- f
  
  return(list(LL = LL, flag = flag))
}
