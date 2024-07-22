source("~/Documents/uni/master-dissertation/code-cont/chen/shat.R")

# Compute \hat{J}_{\max}
# For computational reasons will return the resolution level, rather than J

Jhat <- function(PP, BB, CJ, CK, TJ, M, n, nL) {
  lb <- numeric(nL + 1)
  ub <- numeric(nL + 1)
  
  for (ll in 1:(nL + 1)) {
    s <- tryCatch({
      # debugonce(shat)
      shat(PP[, (CJ[ll] + 1):(CJ[ll + 1])], BB[, (CK[ll] + 1):(CK[ll + 1])])
    }, error = function(e) {
      1e-20
    })
    
    J <- TJ[ll]
    lb[ll] <- J * sqrt(log(J)) * max(0 * (log(n))^4, 1 / s)
  }
  # browser()
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
