library(foreach)
library(doParallel)
library(future)
library(future.apply)

Jlep <- function(Lhat, Px, PP, BB, CJ, CK, TJ, y, n, nb) {
  Lmax <- Lhat + 1
  Jmax <- TJ[Lhat + 1]
  omega <- matrix(rnorm(n * nb), n, nb)
  
  # Step 1: compute critical value
  ZZ <- array(0, dim = c(Lmax, Lmax, nb))
  HH <- matrix(0, nrow = Lmax, ncol = Lmax)
  
  # Set up parallel backend
  cores <- detectCores() - 1
  plan(multisession, workers = cores)
  
  # Pre-compute npiv results
  npiv_results <- future_lapply(1:Lmax, function(idx) {
    PP_sub <- PP[, (CJ[idx] + 1):(CJ[idx + 1])]
    BB_sub <- BB[, (CK[idx] + 1):(CK[idx + 1])]
    npiv(PP_sub, BB_sub, y)
  })
  
  # Main computation
  result <- future_lapply(1:(Lmax-1), function(i) {
    lapply((i+1):Lmax, function(j) {
      Px1 <- Px[, (CJ[i] + 1):(CJ[i + 1])]
      Px2 <- Px[, (CJ[j] + 1):(CJ[j + 1])]
      
      npiv1 <- npiv_results[[i]]
      npiv2 <- npiv_results[[j]]
      
      c1 <- npiv1$c
      u1 <- npiv1$u
      Q1 <- npiv1$Q
      c2 <- npiv2$c
      u2 <- npiv2$u
      Q2 <- npiv2$Q
      
      BB1 <- BB[, (CK[i] + 1):(CK[i + 1])]
      BB2 <- BB[, (CK[j] + 1):(CK[j + 1])]
      
      Bu1 <- sweep(BB1, 1, u1, `*`)
      Bu2 <- sweep(BB2, 1, u2, `*`)
      
      OL1 <- Q1 %*% (t(Bu1) %*% Bu1 / n) %*% t(Q1)
      OL2 <- Q2 %*% (t(Bu2) %*% Bu2 / n) %*% t(Q2)
      OL12 <- Q1 %*% (t(Bu1) %*% Bu2 / n) %*% t(Q2)
      
      tden <- sqrt(rowSums((Px1 %*% OL1) * Px1) + rowSums((Px2 %*% OL2) * Px2) - 2 * rowSums((Px1 %*% OL12) * Px2))
      tnum <- sqrt(n) * (Px1 %*% c1 - Px2 %*% c2)
      HH_val <- max(abs(tnum / tden))
      
      ZZ_vals <- sapply(1:nb, function(b) {
        Buw1 <- t(Bu1) %*% omega[, b] / sqrt(n)
        Buw2 <- t(Bu2) %*% omega[, b] / sqrt(n)
        tnum_boot <- Px1 %*% Q1 %*% Buw1 - Px2 %*% Q2 %*% Buw2
        max(abs(tnum_boot / tden))
      })
      
      list(HH = HH_val, ZZ = ZZ_vals, i = i, j = j)
    })
  })
  
  # Unpack results
  for (res in unlist(result, recursive = FALSE)) {
    HH[res$i, res$j] <- res$HH
    ZZ[res$i, res$j, ] <- res$ZZ
  }
  
  # Critical value
  z <- apply(ZZ, 3, max)
  theta <- quantile(z, max(0.5, 1 - sqrt(log(Jmax) / Jmax)))
  
  # Step 2: cut-off rule
  LL <- which.max(apply(HH, 1, max) <= 1.1 * theta) - 1
  
  # Stop parallel backend
  plan(sequential)
  
  return(list(LL = LL, theta = theta))
}