set.seed(1234)
library(truncnorm)
library(binscatteR)
generate_dataset <- function() {
  # Define the number of hospitals
  num_schools <- 1000
  
  # Generate a binary treatment indicator where 10% are untreated (0) and 90% are treated (1)
  treatment <- rbinom(num_schools, 1, 0.9)
  
  # Initialize dose with zeros for all hospitals
  dose <- rep(0, num_schools)
  
  # For the treated hospitals, generate dose from a normal distribution centered around 0.5 with sd of 0.16
  dose[treatment == 1] <- rtruncnorm(sum(treatment), a = 0, b = 1, mean = 0.5, sd = 0.16)
  
  # Calculate the dependent variable dy using the quadratic equation
  dy <- -4 * (dose-.5)^2 + 1 
  dy[treatment == 0] <- rnorm(sum(treatment == 0), mean = 0.004687925, sd = 0.1013251)
  # Create DataFrame
  data.frame(
    uni_ident = 1:num_schools,
    dy = dy, 
    dose = dose
  )
}

#debugonce(generate_dataset)
# Generate the dataset
data <- generate_dataset()
colnames(data) <- c("hospital_id", "d_capital_labor_ratio", "medicare_share_1983")
save(data, file="/home/oddish3/Downloads/medicare2.RData")

dose <- data$medicare_share_1983
dy <- data$d_capital_labor_ratio

ggplot(data.frame(dose=dose), aes(x=dose)) + 
  geom_histogram()
# a non trivial fraction of units are untreated while common values of traetment
# are around 0.5 and this decreases as we move towards 0 and 1

binnedout <- binscatter(data=data, x="medicare_share_1983", y="d_capital_labor_ratio")
binnedout

twfe <- lm(dy ~ dose)
# summary(twfe)$coefficients
# if we were hoping to estimate ACRT^o, then we did not do it well here, 
# recall it is 0, but here the coefficient on the dose is 1.08 and strongly stat
# diff from 0

#' nonparametric estimates of att(d|d) and acrt(d|d)
#' @param dy the change in the outcome over time
#' @param dose the amount of the treatment
#' @return list( 
#'            local_effects - data frame containing the dose and estimates of 
#'              att(dose) and acrt(dose)
#'            att.overall - an estimate of the overall att
#'            acrt.overall - an estimate of the overall acrt
#'          )
cont_did <- function(dy, dose) {
  # choose bandwidth
  bw <- np::npregbw(formula=dy ~ dose,
                    regtype="ll",
                    bws=1.06,
                    bwscaling=TRUE,
                    bandwidth.compute=FALSE)
  # estimate att and acrt nonparametrically
  out <- np::npreg(bws=bw, gradients=TRUE, exdat=dose)
  
  # order from smallest to largest dose and drop untreated
  this_order <- order(dose)
  dose <- dose[this_order]
  dy <- dy[this_order]
  att.d <- out$mean[this_order]
  acrt.d <- out$grad[,1][this_order]
  att.d <- att.d[dose>0]
  acrt.d <- acrt.d[dose>0]
  att.overall <- mean(att.d)
  acrt.overall <- mean(acrt.d)
  
  return(list(local_effects=data.frame(dose=dose[dose>0],
                                       att.d=att.d,
                                       acrt.d=acrt.d),
              att.overall=att.overall,
              acrt.overall=acrt.overall))
}
#' @param l a particular value of the treatment for which to compute weights
#' @param D an nx1 vector containing doses for all units
cont_twfe_weights <- function(l, D) {
  wt <- ( ( mean(D[D>=l]) - mean(D) ) * mean(1*(D>=l)) ) / var(D)
  wt
}
# debug(np::npreg)
cont_res <- cont_did(dy, dose)

plot_df <- cont_res$local_effects

colnames(plot_df) <- c("dose", "att", "acrt")
ggplot(plot_df, aes(x=dose, att)) +
  geom_hline(yintercept=0, color="red", linetype="dashed") +
  geom_line() +
  theme_bw()

ggplot(plot_df, aes(x=dose, acrt)) +
  geom_hline(yintercept=0, color="red", linetype="dashed") +
  geom_line() +
  theme_bw()

#-----------------------------------------------------------------------------
dL <- min(dose[dose>0])
dU <- max(dose)
# density of the dose
dose_grid <- seq(dL, dU, length.out=100)
frq_weights_plot <- ggplot(data.frame(dose=dose[dose>0]), aes(x=dose)) +
  geom_density(colour = "darkblue", linewidth = 1.2) +
  xlim(c(min(dose_grid), max(dose_grid)))+
  ylab("Density weights") +
  xlab("Dose") +
  ylim(c(0,3)) + 
  labs(title="Density of dose")
frq_weights_plot

twfe_weights <- sapply(dose_grid, cont_twfe_weights, D=dose)
cont_res <- cont_did(dy, dose)

# cont_res$att.overall
# cont_res$acrt.overall


plot_df <- cbind.data.frame(twfe_weights, dose_grid)

twfe_weights_plot <- ggplot(data=plot_df,
                            mapping=aes(x = dose_grid,
                                        y = twfe_weights)) +
  geom_line(colour = "darkblue", linewidth = 1.2) +
  xlim(c(min(dose_grid),
         max(dose_grid)))+
  ylab("TWFE weights") +
  xlab("Dose") +
  geom_vline(xintercept = mean(dose),
             colour="black",
             linewidth = 0.5,
             linetype = "dotted") +
  ylim(c(0,3)) +
  labs(title="TWFE weights")

twfe_weights_plot

library(fixest)
# split data into medicareshare greater than 1
data$binary <- ifelse(data$medicare_share_1983>0, 1, 0)
binarised <- feols(d_capital_labor_ratio ~ binary, data=data)
twfe <- lm(dy ~ dose)
summary(binarised)$coefficients[2]
cont_res$att.overall
summary(twfe)$coefficients[2,1]
print(cont_res$acrt.overall)
