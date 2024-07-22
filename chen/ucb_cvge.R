ucb_cvge <- function(h0, hhat, sigh, zast, theta, A) {
  # browser()
  loss <- max(abs(hhat - h0))
  tmax <- max(abs((hhat - h0) / sigh))
  
  check <- matrix(0, length(zast), length(A))
  
  for (i in 1:length(zast)) {
    for (j in 1:length(A)) {
      check[i, j] <- as.numeric(tmax <= zast[i] + A[j] * theta)
    }
  }
  
  return(list(check = check, loss = loss))
}
