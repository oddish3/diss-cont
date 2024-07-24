# Main file to run nonparametric regression on real data
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)


sourceCpp("~/Documents/uni/master-dissertation/code-cont/chen/jhat.cpp")
# sourceCpp("~/Documents/uni/master-dissertation/code-cont/chen/npiv.cpp")
sourceCpp("~/Documents/uni/master-dissertation/code-cont/chen/jlep.cpp")
sourceCpp("~/Documents/uni/master-dissertation/code-cont/chen/npiv_estimate.cpp")
sourceCpp("~/Documents/uni/master-dissertation/code-cont/chen/ucb_cc.cpp")
# sourceCpp("~/Documents/uni/master-dissertation/code-cont/chen/bspline.cpp")
# sourceCpp("~/Documents/uni/master-dissertation/code-cont/chen/ucb_cv.cpp")

source("~/Documents/uni/master-dissertation/code-cont/chen/ucb_cvge.R")
source("~/Documents/uni/master-dissertation/code-cont/chen/ucb_cv.R")
source("~/Documents/uni/master-dissertation/code-cont/chen/bspline.R")


# Load and prepare data
data <- read.csv('~/Downloads/medicare1.csv')

# Extract the columns (x is Medicare share, y is capital-labor ratio)
x <- data$medicare_share_1983
y <- data$d_capital_labor_ratio

# Remove values where x (medicare_share_1983) is 0
nonzero_indices <- x != 0
x <- x[nonzero_indices]
y <- y[nonzero_indices]

n <- length(x)

# Inputs
nx <- 1000 # number of points for grid of x values
nL <- 9 # maximum resolution level for J
r <- 4 # B-spline order
M <- 5 # ex ante upper bound on sup_x h_0(x)
alpha <- c(0.10, 0.05, 0.01) # level of significance


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
                    y, as.integer(n), 1000)
Llep <- Jlep_result[["LL"]]
thet <- Jlep_result[["theta"]]

# Compute \tilde{J} resolution level
Ltil <- max(min(Llep, Lhat - 1), 0)

# Compute estimator and pre-asymptotic standard error
npiv_result <- npiv_estimate_cpp(Ltil, Px, PP, PP, CJ, CJ, y, n)
hhat <- npiv_result$hhat
sigh <- npiv_result$sigh

# Compute critical value for UCB
zast <- ucb_cv(Ltil, Lhat, Px, PP, PP, CJ, CJ, y, n, 1000, 0, alpha)

# Plot results
plot(x, y, pch = 20, col = rgb(0.75, 0.75, 0.75, 0.5), 
     xlab = 'Medicare Share 1983', ylab = 'Capital-Labor Ratio',
     main = 'Nonparametric Regression: Capital-Labor Ratio vs Medicare Share')
lines(Xx_sub, hhat, col = 'black', lwd = 2)
lines(Xx_sub, hhat + (zast[2] + thet * log(log(TJ[Llep + 1]))) * sigh, col = 'black', lty = 2, lwd = 2)
lines(Xx_sub, hhat - (zast[2] + thet * log(log(TJ[Llep + 1]))) * sigh, col = 'black', lty = 2, lwd = 2)
legend('topright', legend = c('Data', 'Estimated Function', 'Confidence Bands'),
       pch = c(20, NA, NA), lty = c(NA, 1, 2), col = c('grey', 'black', 'black'))
