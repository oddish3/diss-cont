# Estimate function coefficient, return residuals
library(MASS)
npiv <- function(P, B, y) {
  Q <- ginv(t(P) %*% B %*% ginv(t(B) %*% B) %*% t(B) %*% P) %*% t(P) %*% B %*% ginv(t(B) %*% B)
  c <- Q %*% t(B) %*% y
  uhat <- y - P %*% c
  QQ <- NULL
  
  if (exists("QQ")) {
    QQ <- Q * length(y)
  }
  
  return(list(c = c, uhat = uhat, QQ = QQ))
}
