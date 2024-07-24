# Main file to run regression simulations
rm(list = ls())
cat("\014")

library(splines2)
library(tictoc)
library(Rcpp)
library(RcppArmadillo)
library(dplyr)
library(tibble)

library(R.matlab)

dat <- readMat("/home/oddish3/Documents/M_folder/CCK2/data_iteration_1.mat")
u <- dat$u
x <- dat$x
y <- dat$y


debugging <- T

sourceCpp("~/Documents/uni/master-dissertation/code-cont/chen/jhat.cpp")
# sourceCpp("~/Documents/uni/master-dissertation/code-cont/chen/npiv.cpp")
sourceCpp("~/Documents/uni/master-dissertation/code-cont/chen/jlep.cpp")
sourceCpp("~/Documents/uni/master-dissertation/code-cont/chen/npiv_estimate.cpp")
# sourceCpp("~/Documents/uni/master-dissertation/code-cont/chen/ucb_cc.cpp")
# sourceCpp("~/Documents/uni/master-dissertation/code-cont/chen/bspline.cpp")
# sourceCpp("~/Documents/uni/master-dissertation/code-cont/chen/ucb_cv.cpp")

source("~/Documents/uni/master-dissertation/code-cont/chen/ucb_cvge.R")
source("~/Documents/uni/master-dissertation/code-cont/chen/ucb_cv.R")
source("~/Documents/uni/master-dissertation/code-cont/chen/bspline.R")

array_value <- 1 # as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (array_value <= 4) {
  trimming <- TRUE
  n_index <- array_value
} else {
  trimming <- FALSE
  n_index <- array_value - 4
}
print(trimming)
print(n_index)

# Inputs
nn <- c(1250, 2500, 5000, 10000) # sample sizes
nm <- 1000                      # number of replications
nb <- 1000                      # number of bootstrap draws per replication
nx <- 1000                      # number of points for grid of x values
nL <- 9                         # maximum resolution level for J
r <- 4                          # B-spline order
M <- 5                          # ex ante upper bound on \sup_x h_0(x)
alpha <- c(0.10, 0.05, 0.01)    # level of significance
nj <- 4


# Pre-compute
TJ <- 2^((0:nL) + 0) + r - 1
CJ <- c(0, cumsum(TJ))
Lhat <- rep(0, nm)
Llep <- rep(0, nm)
Ltil <- rep(0, nm)
flag <- rep(0, nm)
thet <- rep(0, nm)
zast <- matrix(0, nm, length(alpha))
zdet <- array(0, dim = c(nm, length(alpha), nj))
cvge <- array(0, dim = c(nm, length(alpha), nj + 1))
loss <- array(0, dim = c(nm, nj + 1))
rati <- array(0, dim = c(nm, nj))

Xx <- seq(0, 1, by = 1/nx)

if (trimming == TRUE) {
  Xx_sub <- Xx[Xx > 0.01 & Xx <= 0.99]
} else {
  Xx_sub <- Xx
}
# browser()
Px <- matrix(0, nrow = length(Xx_sub), ncol = CJ[length(CJ)])
for (ll in 0:nL) {
  Px[, (CJ[ll + 1] + 1):CJ[ll + 2]] <- bspline(Xx_sub, ll, r)  # Note: bspline function needs to be defined in R
}

h0 <- sin(15 * pi * Xx) * cos(Xx)

# Simulations
set.seed(1234567)

n <- nn[n_index]


for (j in if(debugging) 1 else 1:nm) { #  j=1
  
  if (j %% 25 == 0) {
    cat(sprintf("j = %d \n", j))
  }
  if (j == nm) {
    cat("\n")
  }
  # browser()
  # Simulate data
  # x <- runif(n)
  # u <- rnorm(n)
  # y <- sin(15 * pi * x) * cos(x) + u
  
  # Pre-compute basis functions and store in arrays PP and BB
  PP <- matrix(0, nrow = n, ncol = CJ[length(CJ)])
  for (ll in 0:nL) {
    start_col <- CJ[ll + 1] +1
    end_col <- CJ[ll + 2]
    # cat("Iteration:", ll, "\n")
    # cat("Start column:", start_col, "\n")
    # cat("End column:", end_col, "\n")
    # cat("Dimensions of bspline output:", dim(bspline(x, ll, r)), "\n")
    PP[, start_col:end_col] <- bspline(x, ll, r)
    # cat("Dimensions of PP after assignment:", dim(PP), "\n\n")
  }
  
  # Compute \hat{J}_{\max} resolution level
  tic()
  result1 <- jhat(PP, PP, as.integer(CJ), as.integer(CJ), as.integer(TJ), M, as.integer(n), as.integer(nL))
  toc()
  Lhat[j] <- result1$LL
  flag[j] <- result1$flag
  
  # Compute Lepski method resolution level
  tic()
  result2 <- jlep(as.integer(Lhat[j]), Px, PP, PP, as.integer(CJ), as.integer(CJ), as.integer(TJ), y, as.integer(n), as.integer(nb))
  toc()
  Llep[j] <- result2$LL
  thet[j] <- result2$theta
  
  # Compute \tilde{J} resolution level
  Ltil[j] <- max(min(Llep[j], Lhat[j] - 1), 0)
  
  # Compute estimator and pre-asymptotic standard error
  # debugonce(npiv_estimate)
  npiv_result <- npiv_estimate_cpp(Ltil[j], Px, PP, PP, CJ, CJ, y, n)
  hhat <- npiv_result$hhat
  sigh <- npiv_result$sigh
  
  # Compute critical value for UCB
  zast[j, ] <- ucb_cv(Ltil[j], Lhat[j], Px, PP, PP, CJ, CJ, y, n, nb, 0, alpha)
  
  # Compute sup-norm loss and excess width
  loss[j, 1] <- max(abs(h0[which(Xx %in% Xx_sub)] - hhat))
  
  # Compute coverage
  results3 <-  ucb_cvge(h0[which(Xx %in% Xx_sub)], hhat, sigh, zast[j, ], thet[j], log(log(TJ[Llep[j] + 1])))
  cvge[j, , 1] <- results3$check
  
  for (k in if(debugging) 1 else 1:nj) {
    # Compute undersmoothed estimator and pre-asymptotic standard error
    result <- npiv_estimate_cpp(k + 2, Px, PP, PP, CJ, CJ, y, n)
    hha1 <- result$hha
    sig1 <- result$sig
    
    # Compute deterministic J critical value for undersmoothed UCB
    zdet[j, , k] <- ucb_cc(k + 2, Px, PP, PP, CJ, CJ, y, n, nb, 0, alpha)
    # browser()
    # Compute sup-norm loss and excess width
    loss[j, 1 + k] <- max(abs(h0[Xx %in% Xx_sub] - hha1))
    rati[j, k] <- max(sig1) / max(sigh)
    
    # Compute coverage
    results4 <- ucb_cvge(h0[Xx %in% Xx_sub], hha1, sig1, zdet[j, , k], 0, 0)
    cvge[j, , 1 + k] <- results4$check
  }
}

# Calculate mean and median for the loss matrix
calculate_mean_median <- function(x) {
  x_non_zero <- x[x != 0]
  c(mean = mean(x_non_zero), median = median(x_non_zero))
}

# Calculate mean and median for the loss matrix
loss_summary <- apply(loss, 2, calculate_mean_median)

# Calculate mean and median for the cvge matrix
cvge_summary <- apply(cvge, 2, calculate_mean_median)

# Calculate mean and median for the rati matrix
rati_summary <- apply(rati, 2, calculate_mean_median)

# Print the summaries in a formatted table
print(loss_summary)
print(cvge_summary)
print(rati_summary)


