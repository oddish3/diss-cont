ucb_cv <- function(Ltil, Lhat, Px, PP, BB, CJ, CK, y, n, nb, type, alpha) {
  Lmax <- Lhat + 1
  omega <- matrix(rnorm(n * nb), n, nb)
  ZZ <- matrix(0, max(Ltil, Lmax - 2), nb)
  
  for (i in 1:max(Ltil, Lmax - 2)) {
    
    # Precompute that which can be pre-computed
    Px1 <- Px[, (CJ[i] + 1):(CJ[i + 1])]
    PP1 <- PP[, (CJ[i] + 1):(CJ[i + 1])]
    BB1 <- BB[, (CK[i] + 1):(CK[i + 1])]
    
    # Note: multiply Px1 by Q1 to get Lx1
    npiv_result <- npiv(PP1, BB1, y)
    u1 <- npiv_result$u
    Q1 <- npiv_result$Q
    
    Bu1 <- sweep(BB1, 1, u1, `*`)
    
    OL1 <- Q1 %*% (t(Bu1) %*% Bu1 / n) %*% t(Q1)
    
    # Variance term
    tden <- numeric(nrow(Px))
    for (x in 1:nrow(Px)) {
      s1 <- Px1[x, ] %*% OL1 %*% Px1[x, ]
      tden[x] <- sqrt(s1)
    }
    
    # Bootstrap
    for (b in 1:nb) {
      Buw1 <- t(Bu1) %*% omega[, b] / sqrt(n)
      
      # Compute bootstrapped sup-t-stat at (J, J2)
      tnum <- Px1 %*% Q1 %*% Buw1
      if (type == 0) {
        ZZ[i, b] <- max(abs(tnum / tden))
      } else if (type == -1) {
        ZZ[i, b] <- max(tnum / tden)
      } else if (type == 1) {
        ZZ[i, b] <- min(tnum / tden)
      }
    }
  }
  
  # Critical value
  if (type == 0 || type == -1) {
    z <- apply(ZZ, 2, max)
    cv <- quantile(z, 1 - alpha)
  } else if (type == 1) {
    z <- apply(ZZ, 2, min)
    cv <- -quantile(z, alpha)
  }
  
  return(cv)
}
