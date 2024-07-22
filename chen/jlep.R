library(foreach)
library(doParallel)

# Define the Jlep function with optimizations
Jlep <- function(Lhat, Px, PP, BB, CJ, CK, TJ, y, n, nb) {
  Lmax <- Lhat + 1
  Jmax <- TJ[Lhat + 1]
  omega <- matrix(rnorm(n * nb), n, nb)
  
  # Step 1: compute critical value
  ZZ <- array(0, dim = c(Lmax, Lmax, nb))
  HH <- matrix(0, nrow = Lmax, ncol = Lmax)
  
  # Register parallel backend
  registerDoParallel(cores = detectCores() - 1)
  
  # Use parallel processing for the main computation
  foreach(i = 1:(Lmax - 1), .combine = 'c', .packages = c("foreach", "parallel", "expm", "MASS")) %dopar% {
    foreach(j = (i + 1):Lmax, .combine = 'c') %dopar% {
      Px1 <- Px[, (CJ[i] + 1):(CJ[i + 1])]
      PP1 <- PP[, (CJ[i] + 1):(CJ[i + 1])]
      BB1 <- BB[, (CK[i] + 1):(CK[i + 1])]
      Px2 <- Px[, (CJ[j] + 1):(CJ[j + 1])]
      PP2 <- PP[, (CJ[j] + 1):(CJ[j + 1])]
      BB2 <- BB[, (CK[j] + 1):(CK[j + 1])]
      
      npiv1 <- npiv(PP1, BB1, y)
      npiv2 <- npiv(PP2, BB2, y)
      c1 <- npiv1$c
      u1 <- npiv1$u
      Q1 <- npiv1$Q
      c2 <- npiv2$c
      u2 <- npiv2$u
      Q2 <- npiv2$Q
      
      Bu1 <- sweep(BB1, 1, u1, `*`)
      Bu2 <- sweep(BB2, 1, u2, `*`)
      
      OL1 <- Q1 %*% (t(Bu1) %*% Bu1 / n) %*% t(Q1)
      OL2 <- Q2 %*% (t(Bu2) %*% Bu2 / n) %*% t(Q2)
      OL12 <- Q1 %*% (t(Bu1) %*% Bu2 / n) %*% t(Q2)
      
      tden <- sqrt(rowSums((Px1 %*% OL1) * Px1) + rowSums((Px2 %*% OL2) * Px2) - 2 * rowSums((Px1 %*% OL12) * Px2))
      tnum <- sqrt(n) * (Px1 %*% c1 - Px2 %*% c2)
      HH[i, j] <- max(abs(tnum / tden))
      
      bootstrapped_stats <- foreach(b = 1:nb, .combine = 'c') %dopar% {
        Buw1 <- t(Bu1) %*% omega[, b] / sqrt(n)
        Buw2 <- t(Bu2) %*% omega[, b] / sqrt(n)
        tnum_boot <- Px1 %*% Q1 %*% Buw1 - Px2 %*% Q2 %*% Buw2
        max(abs(tnum_boot / tden))
      }
      
      ZZ[i, j, ] <- bootstrapped_stats
    }
  }
  
  # Critical value
  z <- apply(ZZ, 3, max)
  theta <- quantile(z, max(0.5, 1 - sqrt(log(Jmax) / Jmax)))
  
  # Step 2: cut-off rule
  LL <- which.max(apply(HH, 1, max) <= 1.1 * theta) - 1
  
  # Stop parallel backend
  stopImplicitCluster()
  
  return(list(LL = LL, theta = theta))
}