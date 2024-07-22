source("~/Documents/uni/master-dissertation/code-cont/chen/evknots.R")

make_data <- function(M, knots) {
  
  ## Make Data File from Matrix M
  filter <- which(rep(1, nrow(M)) == 1)
  d <- list()
  d$M <- M
  d$N <- length(filter)
  d$R_I <- model.matrix(~ M[filter, 1] - 1)
  d$R_J <- model.matrix(~ M[filter, 2] - 1)
  d$R_N_ij <- M[filter, 3]
  d$R_xbar_ij <- M[filter, 4]
  d$R_zij <- M[filter, 5]
  d$R_nE_ij <- M[filter, 10]
  
  ## Create Fixed Effects Dummies
  d$FEA <- cbind(d$R_I, matrix(0, nrow(d$R_I), ncol(d$R_I)), d$R_J[, -ncol(d$R_J)], matrix(0, nrow(d$R_J), ncol(d$R_J) - 1))
  d$FEB <- cbind(matrix(0, nrow(d$R_I), ncol(d$R_I)), d$R_I, matrix(0, nrow(d$R_J), ncol(d$R_J) - 1), d$R_J[, -ncol(d$R_J)])
  
  ## MAKE KNOTS
  k <- knots
  if (knots == 1) {
    Xknots <- d$R_nE_ij
    Zknots <- d$R_zij
    d$ZA_K <- cbind(Zknots, d$FEA)
    d$ZB_K <- cbind(Zknots, d$FEB)
    d$Deriv <- rep(1, length(Zknots))
    k <- NULL
  } else {
    x <- d$R_nE_ij
    evknots_result <- evknots(k, x)
    Xknots <- evknots_result$Xknots
    k <- evknots_result$k
    Deriv <- evknots_result$Deriv
    Zknots <- (abs(Xknots) > 0) * matrix(rep(d$R_zij, ncol(Xknots)), ncol = ncol(Xknots))
    d$Deriv <- Deriv
    d$ZA_K <- cbind(Zknots, d$FEA)
    d$ZB_K <- cbind(Zknots, d$FEB)
  }
  
  d$k <- k
  d$X_K <- Xknots
  
  ## Pre-Compute Projection Matrices (to save on computing cost down the road)
  d$Z_K <- rbind(d$ZA_K, d$ZB_K)
  d$PZ_K <- d$Z_K %*% solve(t(d$Z_K) %*% d$Z_K) %*% t(d$Z_K)
  d$FE <- rbind(d$FEA, d$FEB)
  
  return(d)
}


