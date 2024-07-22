# Main file to run regression simulations
rm(list = ls())
cat("\014")

library(splines2)
library(tictoc)

source("~/Documents/uni/master-dissertation/code-cont/chen/jhat.R")
source("~/Documents/uni/master-dissertation/code-cont/chen/jlep.R")
source("~/Documents/uni/master-dissertation/code-cont/chen/npiv.R")
source("~/Documents/uni/master-dissertation/code-cont/chen/npiv_estimate.R")
source("~/Documents/uni/master-dissertation/code-cont/chen/bspline.R")
source("~/Documents/uni/master-dissertation/code-cont/chen/npiv_estimate.R")
source("~/Documents/uni/master-dissertation/code-cont/chen/ucb_cv.R")
source("~/Documents/uni/master-dissertation/code-cont/chen/ucb_cvge.R")
source("~/Documents/uni/master-dissertation/code-cont/ucb_cc.R")

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
bwd <- array(0, dim = c(980, 4, 10))

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

for (j in 1:nm) { # j=1
  
  if (j %% 25 == 0) {
    cat(sprintf("j = %d \n", j))
  }
  if (j == nm) {
    cat("\n")
  }
  # browser()
  # Simulate data
  x <- runif(n)
  u <- rnorm(n)
  y <- sin(15 * pi * x) * cos(x) + u
  
  # Pre-compute basis functions and store in arrays PP and BB
  PP <- matrix(0, nrow = n, ncol = CJ[length(CJ)])
  for (ll in 0:nL) {
    PP[, (CJ[ll + 1] + 1):CJ[ll + 2]] <- bspline(x, ll, r)
  }
  
  # Compute \hat{J}_{\max} resolution level
  tic()
  result1 <- Jhat(PP, PP, CJ, CJ, TJ, M, n, nL) 
  toc()
  Lhat[j] <- result1$LL
  flag[j] <- result1$flag
  
  # Compute Lepski method resolution level
  tic()
  result2 <- Jlep(Lhat[j], Px, PP, PP, CJ, CJ, TJ, y, n, nb)
  toc()
  Llep[j] <- result2$LL
  thet[j] <- result2$theta
  
  # Compute \tilde{J} resolution level
  Ltil[j] <- max(min(Llep[j], Lhat[j] - 1), 0)
  
  # Compute estimator and pre-asymptotic standard error
  # debugonce(npiv_estimate)
  npiv_result <- npiv_estimate(Ltil[j], Px, PP, PP, CJ, CJ, y, n)
  hhat <- npiv_result$hhat
  sigh <- npiv_result$sigh
  
  # Compute critical value for UCB
  zast[j, ] <- ucb_cv(Ltil[j], Lhat[j], Px, PP, PP, CJ, CJ, y, n, nb, 0, alpha)
  
  # Compute sup-norm loss and excess width
  loss[j, 1] <- max(abs(h0[which(Xx %in% Xx_sub)] - hhat))
  
  # Compute coverage
  result3 <- ucb_cvge(h0[which(Xx %in% Xx_sub)], hhat, sigh, zast[j, ], thet[j], log(log(TJ[Llep[j] + 1])))
  cvge[j, , 1] <- result3$check
  
  for (k in 1:nj) {
    # Compute undersmoothed estimator and pre-asymptotic standard error
    result <- npiv_estimate(k + 2, Px, PP, PP, CJ, CJ, y, n)
    hha1 <- result$hha
    sig1 <- result$sig
    
    # Compute deterministic J critical value for undersmoothed UCB
    zdet[j, , k] <- ucb_cc(k + 2, Px, PP, PP, CJ, CJ, y, n, nb, 0, alpha)
    
    # Compute sup-norm loss and excess width
    loss[j, 1 + k] <- max(abs(h0[Xx %in% Xx_sub] - hha1))
    rati[j, k] <- max(sig1) / max(sigh)
    
    # Compute coverage
    result4 <- ucb_cvge(h0[Xx %in% Xx_sub], hha1, sig1, zdet[j, , k], 0, 0)
    cvge[j, , 1 + k] <- result4$check
  }
}

colMeans(loss[1:12,])
colMeans(cvge[1:12, , ])
colMeans(rati[1:12,])

save_filename <- paste0("./results/regression_", array_value, ".RData")
save(loss, cvge, rati, zdet, zast, Llep, Ltil, Lhat, thet, TJ, 
     file = save_filename)