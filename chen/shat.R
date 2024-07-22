library(expm)
shat <- function(P, B) {
  Gp <- t(P) %*% P
  Gb <- t(B) %*% B
  S  <- t(B) %*% P
  
  if (min(eigen(Gb)$values) > 0) {
    ss <- svd(solve(sqrtm(Gb)) %*% S %*% solve(sqrtm(Gp)))$d
    s <- min(ss)
  } else {
    s <- 1e-20
  }
  return(s)
}