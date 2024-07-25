# Main file to run nonparametric regression on real data
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)

library(R.matlab)
omega_data <- readMat('/home/oddish3/Documents/M_folder/CCK-OG/CCK/omega_values.mat')
omega <- omega_data$omega

# Source your C++ functions
sourceCpp("~/Documents/uni/master-dissertation/code-cont/chen/jhat.cpp")
sourceCpp("~/Documents/uni/master-dissertation/code-cont/chen/jlep.cpp")
sourceCpp("~/Documents/uni/master-dissertation/code-cont/chen/npiv_estimate.cpp")
sourceCpp("~/Documents/uni/master-dissertation/code-cont/chen/ucb_cc.cpp")
source("~/Documents/uni/master-dissertation/code-cont/chen/ucb_cvge.R")
source("~/Documents/uni/master-dissertation/code-cont/chen/ucb_cv.R")
source("~/Documents/uni/master-dissertation/code-cont/chen/bspline.R")

# Number of iterations
num_iterations <- 10

for (iter in 1:num_iterations) {
  cat("Iteration", iter, "\n")
  set.seed(1234567)
  # Load and prepare data
  data <- read.csv('~/Downloads/medicare1.csv')
  data <- data[data$medicare_share_1983 != 0, ]
  x <- data$medicare_share_1983
  y <- data$d_capital_labor_ratio
  n <- length(x)
  
  # Inputs
  nx <- 1000
  nL <- 9
  r <- 4
  M <- 5
  alpha <- c(0.10, 0.05, 0.01)
  
  # Pre-compute
  TJ <- 2^((0:nL) + 0) + r - 1
  CJ <- c(0, cumsum(TJ))
  
  # Create grid for x
  x_min <- min(x)
  x_max <- max(x)
  Xx <- seq(x_min, x_max, length.out = nx)
  Xx_sub <- Xx
  
  # Compute basis functions
  Px <- matrix(0, nrow = length(Xx_sub), ncol = CJ[length(CJ)])
  for (ll in 0:nL) {
    Px[, (CJ[ll + 1] + 1):CJ[ll + 2]] <- bspline(Xx_sub, ll, r)
  }
  PP <- matrix(0, nrow = n, ncol = CJ[length(CJ)])
  for (ll in 0:nL) {
    PP[, (CJ[ll + 1] + 1):CJ[ll + 2]] <- bspline(x, ll, r)
  }
  
  # Compute \hat{J}_{\max} resolution level
  Jhat_result <- jhat(PP, PP, CJ, CJ, TJ, M, n, nL)
  Lhat <- Jhat_result[["LL"]]
  flag <- Jhat_result[["flag"]]
  
  # Compute Lepski method resolution level
  Jlep_result <- jlep(Lhat, Px, PP, PP, as.integer(CJ), as.integer(CJ), as.integer(TJ), 
                      y, as.integer(n), 1000, omega)
  Llep <- Jlep_result[["LL"]]
  thet <- Jlep_result[["theta"]]
  
  # Compute \tilde{J} resolution level
  Ltil <- max(min(Llep, Lhat - 1), 0)
  
  # Compute estimator and pre-asymptotic standard error
  npiv_result <- npiv_estimate_cpp(Ltil, Px, PP, PP, CJ, CJ, y, n)
  hhat <- npiv_result$hhat
  sigh <- npiv_result$sigh
  
  # Print debug information
  cat("Number of columns in Px1:", ncol(npiv_result$basis_functions), "\n")
  cat("Lhat:", Lhat, "\n")
  cat("Llep:", Llep, "\n")
  cat("Ltil:", Ltil, "\n")
  # cat("CJ[Ltil + 1]:", CJ[Ltil + 2], "\n")
  # cat("CJ[Ltil + 2]:", CJ[Ltil + 3], "\n")
  cat("-------------------\n")
}