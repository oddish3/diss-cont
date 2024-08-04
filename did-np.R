# Main file to run nonparametric regression on real data
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(ggplot2)
library(fixest)

set.seed(1234567)


sourceCpp("~/Documents/uni/master-dissertation/code-cont/chen/jhat.cpp")
# sourceCpp("~/Documents/uni/master-dissertation/code-cont/chen/npiv.cpp")
sourceCpp("~/Documents/uni/master-dissertation/code-cont/chen/jlep.cpp")
sourceCpp("~/Documents/uni/master-dissertation/code-cont/chen/npiv_estimate.cpp")
sourceCpp("~/Documents/uni/master-dissertation/code-cont/chen/ucb_cc.cpp")
# sourceCpp("~/Documents/uni/master-dissertation/code-cont/chen/bspline.cpp") doesnt work yet
# sourceCpp("~/Documents/uni/master-dissertation/code-cont/chen/ucb_cv.cpp")

source("~/Documents/uni/master-dissertation/code-cont/chen/ucb_cvge.R")
source("~/Documents/uni/master-dissertation/code-cont/chen/ucb_cv.R")
source("~/Documents/uni/master-dissertation/code-cont/chen/bspline.R")


# Load and prepare data
dat <- read.csv('~/Downloads/medicare1.csv')

# Filter out rows where medicare_share_1983 is 0
data <- dat[dat$medicare_share_1983 != 0, ]

# Extract x and y and sort them in increasing order
sorted_indices <- order(data$medicare_share_1983)
x <- data$medicare_share_1983[sorted_indices]
y <- data$d_capital_labor_ratio[sorted_indices]
hospital_id <- data$hospital_id[sorted_indices]

n <- length(x)

# Inputs
nx <- 880 # number of points for grid of x values
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
Xx <- seq(x_min, x_max, length.out = nx)  # Note: increasing order to match sorted x
Xx_sub <- Xx

# Compute basis functions
PP <- matrix(0, nrow = n, ncol = CJ[length(CJ)])
Px <- matrix(0, nrow = length(Xx_sub), ncol = CJ[length(CJ)])
Dx <- matrix(0, nrow = length(Xx_sub), ncol = CJ[length(CJ)])
for (ll in 0:nL) {
  PP[, (CJ[ll + 1] + 1):CJ[ll + 2]] <- bspline(x, ll, r)$XX
  Px[, (CJ[ll + 1] + 1):CJ[ll + 2]] <- bspline(Xx_sub, ll, r)$XX
  Dx[, (CJ[ll + 1] + 1):CJ[ll + 2]] <- bspline(Xx_sub, ll, r)$DX
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
spline_dosage_SR <- npiv_result$basis_function

npiv_result <- npiv_estimate_cpp(Ltil, Dx, PP, PP, CJ, CJ, y, n)
dhat <- npiv_result$hhat
sigd <- npiv_result$sigh
spline_dosage_SR1 <- npiv_result$basis_function

# Compute critical value for UCB
zast <- ucb_cv(Ltil, Lhat, Px, PP, PP, CJ, CJ, y, n, 1000, 0, alpha)
zast <- ucb_cv(Ltil, Lhat, Dx, PP, PP, CJ, CJ, y, n, 1000, 0, alpha)

# Plot results
plot(x, y, pch = 20, col = rgb(0.75, 0.75, 0.75, 0.5), 
     xlab = 'Medicare Share 1983', ylab = 'Capital-Labor Ratio',
     main = 'Nonparametric Regression: Capital-Labor Ratio vs Medicare Share')
lines(Xx_sub, hhat, col = 'black', lwd = 2)
lines(Xx_sub, hhat + (zast[2] + thet * log(log(TJ[Llep + 1]))) * sigh, col = 'black', lty = 2, lwd = 2)
lines(Xx_sub, hhat - (zast[2] + thet * log(log(TJ[Llep + 1]))) * sigh, col = 'black', lty = 2, lwd = 2)
legend('topright', legend = c('Data', 'Estimated Function', 'Confidence Bands'),
       pch = c(20, NA, NA), lty = c(NA, 1, 2), col = c('grey', 'black', 'black'))
treated_hospitals <- data[data$medicare_share_1983 > 0, ]

# Calculate ATT^o ----------------------------------------------
dat$binary <- ifelse(dat$medicare_share_1983>0, 1, 0)
binarised <- feols(d_capital_labor_ratio ~ binary, data=dat)
summary(binarised)

find_closest <- function(x, grid) {
  which.min(abs(grid - x))
}

# Plot the data points
plot(Xx_sub, dhat, pch = 20, col = rgb(0.75, 0.75, 0.75, 0.5), 
     xlab = 'Medicare Share 1983', ylab = 'Capital-Labor Ratio (Derivative)',
     main = 'Nonparametric Regression: Derivative of Capital-Labor Ratio vs Medicare Share')
# Plot the estimated derivative function dhat
lines(Xx_sub, dhat, col = 'black', lwd = 2)
# Plot the upper confidence band
lines(Xx_sub, dhat + (zast[2] + thet * log(log(TJ[Llep + 1]))) * sigd, col = 'black', lty = 2, lwd = 2)
# Plot the lower confidence band
lines(Xx_sub, dhat - (zast[2] + thet * log(log(TJ[Llep + 1]))) * sigd, col = 'black', lty = 2, lwd = 2)
# Add legend
legend('topright', legend = c('Data', 'Estimated Derivative', 'Confidence Bands'),
       pch = c(20, NA, NA), lty = c(NA, 1, 2), col = c('grey', 'black', 'black'))



################################################################################

# Filter for positive doses
positive_doses <- dat$medicare_share_1983[dat$medicare_share_1983 > 0]
n_positive <- length(positive_doses)

# Calculate ACR^o using the actual doses
ACR_o <- mean(dhat[findInterval(positive_doses, Xx_sub)])

# Calculate E_n[ΔY|D=0]
E_DY_D0 <- mean(dat$d_capital_labor_ratio[dat$medicare_share_1983 == 0])

calculate_eta <- function(d_i) {
  idx <- findInterval(d_i, Xx_sub)
  ACR_K_D_i <- dhat[idx]
  E_ACR_K_D <- mean(dhat[findInterval(positive_doses, Xx_sub)])
  
  d_psi_K <- Dx[idx, ]
  psi_K <- Px[idx, ]
  
  # Calculate û_i
  u_i <- y[which(x == d_i)] - E_DY_D0 - hhat[idx]
  
  # Add small constant for numerical stability
  PP_inv <- solve(t(PP) %*% PP + diag(1e-8, ncol(PP)))
  
  eta <- ACR_K_D_i - E_ACR_K_D +
    t(d_psi_K) %*% PP_inv %*% PP[which(x == d_i), ] * u_i
  
  return(eta)
}

# Calculate σ^2_ACR^o
eta_values <- sapply(positive_doses, calculate_eta)
sigma_sq_ACR_o <- mean(eta_values^2)

# Standard error for ACR^o
se_ACR_o <- sqrt(sigma_sq_ACR_o / n_positive)

# Calculate confidence intervals
alpha <- 0.05  # for 95% CI
z_score <- qnorm(1 - alpha/2)
ci_lower <- ACR_o - z_score * se_ACR_o
ci_upper <- ACR_o + z_score * se_ACR_o

# Print results
print(paste("ACR^o estimate:", ACR_o))
print(paste("Standard Error:", se_ACR_o))
print(paste("95% CI: [", ci_lower, ",", ci_upper, "]"))






# Create a data frame with sorted x, y, and basis functions
sorted_data <- data.frame(
  medicare_share_1983 = x,
  d_capital_labor_ratio = y,
  hospital_id = hospital_id
)

# Add the spline basis functions to the sorted data frame
for (i in 1:ncol(spline_dosage_SR)) {
  sorted_data[paste0("spline_", i)] <- spline_dosage_SR[, i]
}

# Create a formula string for all spline columns
spline_formula <- paste(paste0("spline_", 1:ncol(spline_dosage_SR)), collapse = " + ")

# Fit the model using the spline columns in the sorted data
splines_SR <- fixest::feols(as.formula(paste("d_capital_labor_ratio ~", spline_formula, "- 1")),
                            data = sorted_data,
                            cluster = ~ hospital_id)

# Compute ATT(d|d)
att_SR <- predict(splines_SR)

# Compute influence functions
n_treated_SR <- nrow(sorted_data)
spline_matrix <- as.matrix(sorted_data[, paste0("spline_", 1:ncol(spline_dosage_SR))])
infl_reg_SR <- splines_SR$residuals * spline_matrix %*% 
  (MASS::ginv(t(spline_matrix) %*% spline_matrix / n_treated_SR))
infl_att_SR <- infl_reg_SR %*% t(spline_matrix)

# Compute standard error
se_att_SR <- sqrt(colMeans(infl_att_SR^2) / n_treated_SR)

# Create results dataframe
results_cdid_SR <- data.frame(
  d = sorted_data$medicare_share_1983,
  att = att_SR,
  p_uci_att = (att_SR + 1.96 * se_att_SR),
  p_lci_att = (att_SR - 1.96 * se_att_SR)
)

# Plot ATT(d|d)
p_att_SR <- ggplot(data = results_cdid_SR,
                   mapping = aes(x = d, y = att)) +
  geom_ribbon(aes(ymin = p_lci_att, ymax = p_uci_att),
              alpha = 0.2,
              size = 1,
              fill = "#E69F00") +
  geom_line(linewidth = 1.2,
            alpha = 2,
            colour = "#E69F00") +
  geom_hline(yintercept = 0,
             colour = "black",
             linewidth = 0.5,
             linetype = "dotted") +
  xlab("County prospectivity score") +
  ylab("Average differences in log(capital-labor ratio)") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(color = "black", size = 12),
        plot.title = element_text(size = 14, hjust = 0))

print(p_att_SR)

data$weight <- ifelse(data$medicare_share_1983 > 0, data$medicare_share_1983, NA_real_)
weighted_att_SR <- sum(results_cdid_SR$att * data$weight / n_treated_SR, na.rm = TRUE)
print(paste("Weighted ATT:", weighted_att_SR))
print(paste("Mean ATT:", mean(results_cdid_SR$att)))

# Assuming npiv_result contains the derivative of the basis functions
derivative_spline_dosage_SR <- npiv_result$derivative_basis_function

# If it doesn't, you might need to compute it numerically
# This is a simple numerical differentiation and might need to be adjusted
# h <- 1e-5  # small step size
# derivative_spline_dosage_SR <- (npiv_result$basis_function[2:nrow(spline_dosage_SR),] - 
#                                 npiv_result$basis_function[1:(nrow(spline_dosage_SR)-1),]) / h

# Compute ACRT using the same coefficients from splines_SR
coefficients_SR <- coef(splines_SR)
acrt_SR <- as.vector(spline_dosage_SR1 %*% coefficients_SR)

# Create a results data frame for ACRT
results_acrt_SR <- data.frame(d = data$medicare_share_1983,
                              acrt = acrt_SR)

# Plot ACRT(d|d)
p_acrt_SR <- ggplot(data = results_acrt_SR,
                    mapping = aes(x = d, y = acrt)) +
  geom_line(linewidth = 1.2,
            alpha = 2,
            colour = "#E69F00") +
  geom_hline(yintercept = 0,
             colour = "black",
             linewidth = 0.5,
             linetype = "dotted") +
  xlab("County prospectivity score") +
  ylab("Average Causal Response to Treatment") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(color = "black", size = 12),
        plot.title = element_text(size = 14, hjust = 0))

print(p_acrt_SR)



# Compute ACRT for positive doses
acrt_SR_filtered <- acrt_SR[data$medicare_share_1983 > 0]
ACR_hat_0 <- mean(acrt_SR_filtered)
print(paste("Mean ACRT for positive doses:", ACR_hat_0))

