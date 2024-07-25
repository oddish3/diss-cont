rm(list=ls())

data1 <- read.csv('~/Downloads/medicare1.csv')
data <- data1
# Filter out rows where medicare_share_1983 is 0
# data <- data[data$medicare_share_1983 != 0, ]

# Extract x and y
x <- data$medicare_share_1983
y <- data$d_capital_labor_ratio
set.seed(123)
u <- rnorm(1000, 0, 1) # Generate 1000 random normal disturbance terms

# d_capital_labor_ratio <- seq(0, 1, length.out = 1000) # Generate linearly increasing x from 0 to 1
# medicare_share_1983 <- sin(15 * pi * d_capital_labor_ratio) + cos(d_capital_labor_ratio) + u # Calculate y as the requested function of x and u
# hospital_id <- 0:(length(d_capital_labor_ratio) - 1) # Generate an incremental hospital_id
# medicare1 <- data.frame(hospital_id, d_capital_labor_ratio, medicare_share_1983)
# ggplot(medicare1, aes(x = d_capital_labor_ratio, y = medicare_share_1983)) +
#   geom_point() +
#   # geom_smooth(method = "loess", se = FALSE) +
#   labs(x = "d_capital_labor_ratio", y = "medicare_share_1983", title = "Medicare Share 1983")




# load("/home/oddish3/Downloads/medicare2.RData")
# medicare1 <- data
library(tidyr)
library(dplyr)
library(fixest)
library(broom)
library(stringr)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(latex2exp)
library(patchwork)
library(ggtext)
library(splines2)

########################################################################################
# Compute ATT(d|d) for all values of the dose d ---------------------------
########################################################################################
# compute cubic splines
splines_SR <- fixest::feols(d_capital_labor_ratio ~ bSpline(medicare_share_1983, 
                                          knots = quantile(medicare_share_1983,probs = c(0.25, 0.5, 0.75)),
                                          degree = 3,
                                          intercept = TRUE) -1,
                            data = data,
                            cluster = ~ hospital_id)

# Compute ATT(d|d) for all values of the dose d
att_SR <- predict(splines_SR)
mean(att_SR)



# Compute influence functions ----------------------------------------------
# Get splines for dosage
spline_dosage_SR <- bSpline(data$medicare_share_1983,
                            knots = quantile(data$medicare_share_1983,probs = c(0.25, 0.5, 0.75)),
                            degree = 3,
                            intercept = TRUE)
# Sample Size
n_treated_SR <- length(data$hospital_id)

# compute influence function of spline beta
infl_reg_SR <- splines_SR$residuals * spline_dosage_SR %*% (MASS::ginv(t(spline_dosage_SR)%*%spline_dosage_SR/n_treated_SR))
infl_att_SR <-  infl_reg_SR %*% t(spline_dosage_SR)

# Compute standard error
se_att_SR <- sqrt(colMeans(infl_att_SR^2)/(n_treated_SR))
results_cdid_SR <- data.frame(d = data$medicare_share_1983,
                              att = att_SR,
                              p_uci_att = (att_SR + 1.96*se_att_SR),
                              p_lci_att = (att_SR - 1.96*se_att_SR))
# Plot ATT(d|d) for Short-Run
p_att_SR <- ggplot(data = results_cdid_SR,
                   mapping = aes(x = d, y = att)) +
  geom_ribbon(aes(ymin= p_lci_att, ymax=  p_uci_att),
              alpha = 0.2,
              size = 1,
              fill = "#E69F00")+
  geom_line(linewidth = 1.2,
            alpha = 2,
            colour = "#E69F00") +
  geom_hline(yintercept = 0,
             colour="black",
             linewidth = 0.5,
             linetype = "dotted")+
  #ylim(range(-0.04, 0.12))+
  #scale_y_continuous(breaks = seq(-0.04, 0.12, 0.02), limits = c(-0.04,0.12))+
  xlab("County prospectivity score") +
  ylab("Average differences in log(employment)") +
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title = element_text(color="black",
                                  size = 12))+
  theme(plot.title=ggtext::element_markdown(size=14,
                                            #face = "bold",
                                            color="black",
                                            hjust=0,
                                            lineheight=1.2)
  )
p_att_SR

# Print results
mean(results_cdid_SR$att)
mean(results_cdid_SR$att[data$medicare_share_1983 != 0])

########################################################################################
# Compute ACRT(d|d) for all values of the dose d ---------------------------
########################################################################################
# Compute the derivative of the spline basis functions
derivative_spline_dosage_SR <- dbs(data$medicare_share_1983,
                                   knots = quantile(data$medicare_share_1983, probs = c(0.25, 0.5, 0.75)),
                                   degree = 3,
                                   intercept = TRUE)

# Compute ACRT using the same coefficients from splines_SR
coefficients_SR <- coef(splines_SR)
acrt_SR <- as.vector(derivative_spline_dosage_SR %*% coefficients_SR)

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
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title = element_text(color = "black", size = 12)) +
  theme(plot.title = ggtext::element_markdown(size = 14, color = "black", hjust = 0, lineheight = 1.2))

print(p_acrt_SR)

# Filter data to exclude dose = 0
medicare2 <- data[data$medicare_share_1983 > 0,]

# Ensure ACRT matches the filtered data length
acrt_SR_filtered <- acrt_SR[data$medicare_share_1983 > 0]

# Sample size with positive dose
n_D_gt_0 <- nrow(medicare2)

# Compute ACRT for positive doses
ACR_hat_0 <- mean(acrt_SR_filtered)

# Print results
print(ACR_hat_0)








