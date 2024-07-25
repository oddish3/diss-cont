source("~/Documents/uni/master-dissertation/code-cont/chen/shat.R")

# Compute \hat{J}_{\max}
# For computational reasons will return the resolution level, rather than J

jhat_r <- function(PP, BB, CJ, CK, TJ, M, n, nL) {
  # 
  # cat("M ", M)
  # cat("\n")
  # cat("n", n)
  # cat("\n")
  # cat("nL", nL)
  # cat("\n")
  # cat("CJ:", CJ)
  # cat("\n")
  # cat("CK:", CK)
  # cat("\n")
  # cat("TJ:", TJ)
  
 
  
  lb <- numeric(nL + 1)
  ub <- numeric(nL + 1)
  
  for (ll in 1:(nL + 1)) {
    s <- tryCatch({
      shat(PP[, (CJ[ll] + 1):(CJ[ll + 1])], BB[, (CK[ll] + 1):(CK[ll + 1])])
    }, error = function(e) {
      1e-20
    })
    
    # cat("ll:", ll)
    # cat("\n")
    # cat("s:", s)
    # cat("\n")
    # 
    J <- TJ[ll]
    # cat("J:",J)
    lb[ll] <- J * sqrt(log(J)) * max(0 * (log(n))^4, 1 / s)
    # cat("\n")
    # cat("lb:", lb[ll])
  }
  # browser()
  ub[1:nL] <- lb[2:(nL + 1)]
  
  # cat("ub:", ub)
  # cat("\n")
  # cat("lb:", lb)
  # cat("\n")
  threshold <- 2 * M * sqrt(n)
  L <- which(lb <= 2 * M * sqrt(n) & 2 * M * sqrt(n) <= ub) - 1
  f <- 0
  # cat("thresh", threshold)
  # cat("\n")
  # cat("l:", L)
  # cat("\n")
  
  if (length(L) == 0) {
    L <- which(lb <= 2 * M * sqrt(n), arr.ind = TRUE)[length(which(lb <= 2 * M * sqrt(n)))] - 1
    f <- 1
  }
  # cat("l1:", L)
  # cat("\n")
  
  if (length(L) == 0) {
    L <- 1
    f <- 2
  }
  # cat("l2:", L)
  # cat("\n")
  
  LL <- max(L, 1)
  # cat("LL:", LL)
  flag <- f
  
  return(list(LL = LL, flag = flag))
}
