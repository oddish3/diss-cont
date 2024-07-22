bspline <- function(x, l, r, kts = NULL) {
  N <- length(x)
  m <- 2^l - 1
  r <- r + 1
  
  # Define the augmented knot set
  if (is.null(kts)) {
    if (l == 0) {
      kts <- c(rep(0, r-1), rep(1, r-1))
    } else if (l >= 1) {
      kts <- c(rep(0, r-2), seq(0, 1, length.out = 2^l + 1), rep(1, r-2))
    }
  }
  
  # Initialize for recursion
  BB <- array(0, dim = c(N, m + 2*r - 2, r-1))
  for (i in 1:N) {
    ix <- which(x[i] >= kts[(r-1):(r+m-1)] & x[i] <= kts[r:(r+m)])[1]
    BB[i, ix + r - 2, 1] <- 1
  }
  
  # Recursion
  for (j in 2:(r-1)) {
    for (i in 1:(m + 2*r - 2 - j)) {
      if (i + j + 1 <= m + 2*r) {
        if (kts[i + j - 1] - kts[i] != 0) {
          a1 <- (x - kts[i]) / (kts[i + j - 1] - kts[i])
        } else {
          a1 <- rep(0, N)
        }
        if (kts[i + j] - kts[i + 1] != 0) {
          a2 <- (x - kts[i + 1]) / (kts[i + j] - kts[i + 1])
        } else {
          a2 <- rep(0, N)
        }
        BB[, i, j] <- a1 * BB[, i, j - 1] + (1 - a2) * BB[, i + 1, j - 1]
      } else if(i + j <= m + 2*r) {
        if (kts[i + j] - kts[i] != 0) {
          a1 <- (x - kts[i]) / (kts[i + j] - kts[i])
        } else {
          a1 <- rep(0, N)
        }
        BB[, i, j] <- a1 * BB[, i, j - 1]
      }
    }
  }
  # browser()
  XX <- BB[, 1:(2^l + r - 2), r - 1]
  
  return(XX)
}
