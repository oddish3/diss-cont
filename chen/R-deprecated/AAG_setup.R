# Main file to run empirical application

# Clear all variables and console
rm(list=ls())
cat("\014")

# Load necessary package
library(splines2)
library(readr)

source("~/Documents/uni/master-dissertation/code-cont/chen/make_data.R")
# Set seed for reproducibility
set.seed(1234567)

trimming <- TRUE

nL <- 10                      # maximum resolution level for J
r <- 4                        # cubic spline as in AAG
l <- 2                        # link between resolution levels for J and K
M <- 10
nx <- 1000                    # number of points for grid of x values
nb <- 1000                    # number of bootstrap draws per replication
alpha <- 0.05                 # level of significance

# Adjustment factors
scale <- 10
center <- 1

TJ <- 2^((0:nL) + 0) + r - 1
TK <- 2^((0:nL) + l) + r
# TK <- l * TJ
CJ <- c(0, cumsum(TJ))
CK <- c(0, cumsum(TK))

# Constants from AAG
kappa_tau <- 0.359
sigma <- 3.9
kappa_tilde_tau <- (1 - sigma) * kappa_tau
kappa_tilde_rho <- 1 / kappa_tilde_tau

# Import data using AAG's command
# Read the CSV file
data <- read_csv("~/Documents/M_folder/CCK1/data/2012.csv")

# Call the make_data function
# debug(make_data)
d <- make_data(data, 3)

# Display the structure of the resulting list
str(d)

# Extract relevant variables and change notations
d$y <- d$R_xbar_ij - kappa_tilde_tau * d$R_zij
d$x <- d$R_nE_ij
d$w <- d$R_zij
d$F <- d$FEB
n <- length(d$x)

# Rescale to (0, 1)
xtemp <- pmin(pmax(d$x / scale + center, 0), 1)
wtemp <- sapply(d$w, function(wi) sum(d$w <= wi) / n)

# Pre-compute basis functions
Xx <- seq(0, 1, length.out = nx + 1)
if (trimming) {
  Xx_sub <- Xx[Xx >= (log(0.001) / scale + center) & Xx <= (log(0.5) / scale + center)]
} else {
  Xx_sub <- Xx
}

Px <- matrix(0, nrow = length(Xx_sub), ncol = CJ[length(CJ)])
Dx <- matrix(0, nrow = length(Xx_sub), ncol = CJ[length(CJ)])
BB <- matrix(0, nrow = n, ncol = CK[length(CK)])
PP <- matrix(0, nrow = n, ncol = CJ[length(CJ)])
DP <- matrix(0, nrow = n, ncol = CJ[length(CJ)])

for (ll in 0:nL) {
  # Compute knots using quantiles
  if (ll == 0) {
    kts <- c(rep(0, r-1), rep(1, r-1))
  } else if (ll == 1) {
    kts <- c(rep(0, r-1), quantile(xtemp, 1 / 2), rep(1, r-1))
  } else {
    kts <- c(rep(0, r-1), quantile(xtemp, seq(1 / 2^ll, 1 - 1 / 2^ll, by = 1 / 2^ll)), rep(1, r-1))
  }
  
  # Compute B-spline basis and derivatives
  PP[, (CJ[ll + 1] + 1):(CJ[ll + 2])] <- bSpline(xtemp, knots = kts, degree = r - 1, intercept = TRUE)
  Px[, (CJ[ll + 1] + 1):(CJ[ll + 2])] <- bSpline(Xx_sub, knots = kts, degree = r - 1, intercept = TRUE)
  Dx[, (CJ[ll + 1] + 1):(CJ[ll + 2])] <- deriv(bSpline(Xx_sub, knots = kts, degree = r - 1, intercept = TRUE))
  BB[, (CK[ll + 1] + 1):(CK[ll + 2])] <- bSpline(wtemp, degree = r + l, intercept = TRUE)
}