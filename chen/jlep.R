# Lepski method
# For computational reasons will return the resolution level, rather than J

Jlep <- function(Lhat, Px, PP, BB, CJ, CK, TJ, y, n, nb) {
 
  Lmax <- Lhat + 1
  Jmax <- TJ[Lhat + 1]
  omega <- matrix(rnorm(n * nb), n, nb)
  # browser()
  
  # Step 1: compute critical value
  ZZ <- array(0, dim = c(Lmax, Lmax, nb))
  HH <- matrix(0, Lmax, Lmax)
  
  for (i in 1:Lmax) {
    if (i == 8) {
      # browser()
    }
    # cat(sprintf("i = %d \n", i))
    
    if (i < Lmax) { #change due to R indexing
      for (j in (i + 1):Lmax) {
        # cat(sprintf("j = %d \n", j))
        
        # Precompute that which can be pre-computed
        Px1 <- Px[, (CJ[i] + 1):(CJ[i + 1])]
        PP1 <- PP[, (CJ[i] + 1):(CJ[i + 1])]
        BB1 <- BB[, (CK[i] + 1):(CK[i + 1])]
        Px2 <- Px[, (CJ[j] + 1):(CJ[j + 1])]
        PP2 <- PP[, (CJ[j] + 1):(CJ[j + 1])]
        BB2 <- BB[, (CK[j] + 1):(CK[j + 1])]
        
        # Note: multiply Px1 by Q1 to get Lx1
        # 
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
        
        # Variance term
        
        tden <- numeric(nrow(Px))
        for (x in 1:nrow(Px)) {
          s1 <- Px1[x, ] %*% OL1 %*% Px1[x, ]
          s2 <- Px2[x, ] %*% OL2 %*% Px2[x, ]
          s12 <- Px1[x, ] %*% OL12 %*% Px2[x, ]
          tden[x] <- sqrt(s1 + s2 - 2 * s12)
        }
        # browser()
        # Compute sup-t-stat at (J, J2)
        tnum <- sqrt(n) * (Px1 %*% c1 - Px2 %*% c2)
        HH[i, j] <- max(abs(tnum / tden))
        
        # Bootstrap
        for (b in 1:nb) {
          
          Buw1 <- t(Bu1) %*% omega[, b] / sqrt(n)
          Buw2 <- t(Bu2) %*% omega[, b] / sqrt(n)
          
          # Compute bootstrapped sup-t-stat at (J, J2)
          tnum <- Px1 %*% Q1 %*% Buw1 - Px2 %*% Q2 %*% Buw2
          ZZ[i, j, b] <- max(abs(tnum / tden))
          
        }
      }    
    }
  
    
  }
  

  # Critical value
  z <- numeric(nb)
  for (b in 1:nb) {
    z[b] <- max(ZZ[, , b])
  }

  # browser()
  theta <- quantile(z, max(0.5, 1 - sqrt(log(Jmax) / Jmax)))
  
  # Step 2: cut-off rule
  LL <- which.max(apply(HH, 1, max) <= 1.1 * theta) - 1
  
  return(list(LL = LL, theta = theta))
}
