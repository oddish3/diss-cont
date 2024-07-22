# Estimate and compute sieve variance

npiv_estimate <- function(Ltil, Px, PP, BB, CJ, CK, y, n) {
  Ltil <- Ltil + 1
  
  Px1 <- Px[, (CJ[Ltil] + 1):(CJ[Ltil + 1])]
  PP1 <- PP[, (CJ[Ltil] + 1):(CJ[Ltil + 1])]
  BB1 <- BB[, (CK[Ltil] + 1):(CK[Ltil + 1])]
  
  # Note: multiply Px1 by Q1 to get Lx1
  npiv_result <- npiv(PP1, BB1, y)
  c1 <- npiv_result$c
  u1 <- npiv_result$uhat
  Q1 <- npiv_result$QQ
  
  Bu1 <- sweep(BB1, 1, u1, `*`)
  OL1 <- Q1 %*% (t(Bu1) %*% Bu1 / n) %*% t(Q1)
  
  # Variance term
  tden <- numeric(nrow(Px))
  for (x in 1:nrow(Px)) {
    s1 <- Px1[x, ] %*% OL1 %*% Px1[x, ]
    tden[x] <- sqrt(s1)
  }
  
  hhat <- Px1 %*% c1
  sigh <- tden / sqrt(n)
  
  return(list(hhat = hhat, sigh = sigh))
}
