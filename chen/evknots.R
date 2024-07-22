# Define the evknots function
evknots <- function(k, x) {
  # If k is 1, set Xknots to x and return
  if (k == 1) {
    Xknots <- matrix(x, ncol = 1)
    Deriv <- rep(1, length(x))
    return(list(Xknots = Xknots, k = k, Deriv = Deriv))
  }
  
  # Define knots based on quantiles of x
  Xknots <- matrix(0, nrow = length(x), ncol = k)
  for (i in 1:k) {
    # Create a knot at the i-th quantile of x
    Xknots[, i] <- quantile(x, i / k)
  }
  
  # Set derivatives to 1
  Deriv <- matrix(1, nrow = length(x), ncol = k)
  return(list(Xknots = Xknots, k = k, Deriv = Deriv))
}