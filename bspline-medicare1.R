rm(list=ls())
load("/home/oddish3/Downloads/medicare1.RData")

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
                            data = medicare1,
                            cluster = ~ hospital_id)

# Compute ATT(d|d) for all values of the dose d
att_SR <- predict(splines_SR)

# Compute influence functions ----------------------------------------------
# Get splines for dosage
spline_dosage_SR <- bSpline(medicare1$medicare_share_1983,
                            knots = quantile(medicare1$medicare_share_1983,probs = c(0.25, 0.5, 0.75)),
                            degree = 3,
                            intercept = TRUE)
# Sample Size
n_treated_SR <- length(medicare1$hospital_id)

# compute influence function of spline beta
infl_reg_SR <- splines_SR$residuals * spline_dosage_SR %*% (MASS::ginv(t(spline_dosage_SR)%*%spline_dosage_SR/n_treated_SR))
infl_att_SR <-  infl_reg_SR %*% t(spline_dosage_SR)

# Compute standard error
se_att_SR <- sqrt(colMeans(infl_att_SR^2)/(n_treated_SR))
results_cdid_SR <- data.frame(d = medicare1$medicare_share_1983,
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
# weight results_cdid_SR$att by dose distribution itself ignoring dose = 0
medicare1$weight <- ifelse(medicare1$medicare_share_1983>0 ,medicare1$medicare_share_1983, NA_real_)

weighted_att_SR <- results_cdid_SR$att * medicare1$medicare_share_1983 / n_treated_SR

weighted_att_SR <- sum(weighted_att_SR, na.rm = TRUE)

mean(results_cdid_SR$att)

########################################################################################
# Compute ACRT(d|d) for all values of the dose d ---------------------------
########################################################################################
# Compute the derivative of the spline basis functions
derivative_spline_dosage_SR <- dbs(medicare1$medicare_share_1983,
                                   knots = quantile(medicare1$medicare_share_1983, probs = c(0.25, 0.5, 0.75)),
                                   degree = 3,
                                   intercept = TRUE)

# Compute ACRT using the same coefficients from splines_SR
coefficients_SR <- coef(splines_SR)
acrt_SR <- as.vector(derivative_spline_dosage_SR %*% coefficients_SR)

# Create a results data frame for ACRT
results_acrt_SR <- data.frame(d = medicare1$medicare_share_1983,
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
medicare2 <- medicare1[medicare1$medicare_share_1983 > 0,]

# Ensure ACRT matches the filtered data length
acrt_SR_filtered <- acrt_SR[medicare1$medicare_share_1983 > 0]

# Sample size with positive dose
n_D_gt_0 <- nrow(medicare2)

# Compute ACRT for positive doses
ACR_hat_0 <- mean(acrt_SR_filtered)

# Print results
print(ACR_hat_0)








