rm(list=ls())
library(stats)

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


n = 1000
Xsi.ps = 0.75


# Define the DGP function
generate_dgp <- function(n, Xsi.ps) {
  # Generate covariates
  x1 <- stats::rnorm(n, mean = 0, sd = 1)
  x2 <- stats::rnorm(n, mean = 0, sd = 1)
  x3 <- stats::rnorm(n, mean = 0, sd = 1)
  x4 <- stats::rnorm(n, mean = 0, sd = 1)
  
  z1 <- exp(x1 / 2)
  z2 <- x2 / (1 + exp(x1)) + 10
  z3 <- (x1 * x3 / 25 + 0.6)^3
  z4 <- (x1 + x4 + 20)^2
  
  mean.z1 <- mean(z1)
  mean.z2 <- mean(z2)
  mean.z3 <- mean(z3)
  mean.z4 <- mean(z4)
  
  sd.z1 <- sd(z1)
  sd.z2 <- sd(z2)
  sd.z3 <- sd(z3)
  sd.z4 <- sd(z4)
  
  z1 <- (z1 - mean.z1) / sd.z1
  z2 <- (z2 - mean.z2) / sd.z2
  z3 <- (z3 - mean.z3) / sd.z3
  z4 <- (z4 - mean.z4) / sd.z4
  
  x <- cbind(x1, x2, x3, x4)
  z <- cbind(z1, z2, z3, z4)
  
  # Generate groups
  # Propensity score
  pi <- stats::plogis(Xsi.ps * (-z1 + 0.5 * z2 - 0.25 * z3 - 0.1 * z4))
  
  # Use the propensity score as the continuous treatment variable
  d <- pi
  # d <- ifelse(pi<0.1, 0, pi)
  
  # Generate aux indexes for the potential outcomes
  index.lin <- 210 + 27.4 * z1 + 13.7 * (z2 + z3 + z4)
  index.unobs.het <- d * (index.lin)
  index.att <- 0  # 10 + 27.4 * x1 + 13.7 * (x2 + x3 + x4)
  
  # This is the key for consistency of outcome regression
  index.trend <- 210 + 27.4 * z1 + 13.7 * (z2 + z3 + z4)
  
  # v is the unobserved heterogeneity
  v <- stats::rnorm(n, mean = index.unobs.het, sd = 1)
  
  # Generate realized outcome at time 0
  y0 <- index.lin + v + stats::rnorm(n)
  
  # Generate outcomes at time 1
  # First let's generate potential outcomes: y_1_potential
  y10 <- index.lin + v + stats::rnorm(n, mean = 0, sd = 1) +  # This is the baseline
    index.trend  # This is for the trend based on X
  
  y11 <- index.lin + v + stats::rnorm(n, mean = 0, sd = 1) +  # This is the baseline
    index.trend +  # This is for the trend based on X
    index.att  # This is the treatment effect
  
  # Get unfeasible att
  att.unf <- (mean(d * y11) - mean(d * y10)) / mean(d)
  
  # Generate realized outcome at time 1
  y1 <- d * y11 + (1 - d) * y10
  
  dy = y1 - y0
  dy =  (dy - min(dy)) / (max(dy) - min(dy))
  return(data.frame(i = 1:n, dy = dy, dose = d #x = z, y1 = y1, y0 = y0, pi = pi
              #att = index.att, att.unf = att.unf, eff = 11.10
              ))
}

data <- generate_dgp(n, Xsi.ps)
# plot(data$dose)

# summary(data$dose)
# summary(data$dy)

ggplot(data.frame(dose=data$dose), aes(x=data$dose)) + 
  geom_histogram()

binnedout <- binscatter(data=data, x="dose", y="dy")
binnedout

twfe <- lm(dy ~ dose, data = data)
# summary(twfe)$coefficients[2,1]

cont_res <- cont_did(data$dy, data$dose)
cont_res$att.overall #0.35
cont_res$acrt.overall #-0.25

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
dL <- min(data$dose[data$dose>0])
dU <- max(data$dose)
# density of the dose
dose_grid <- seq(dL, dU, length.out=100)
frq_weights_plot <- ggplot(data.frame(dose=data$dose[data$dose > 0]), aes(x=dose)) +
  geom_density(colour = "darkblue", linewidth = 1.2) +
  xlim(c(min(data$dose), max(data$dose))) +  # Use data$dose instead of dose_grid
  ylab("Density weights") +
  xlab("Dose") +
  ylim(c(0, 3)) + 
  labs(title="Density of dose")
frq_weights_plot

twfe_weights <- sapply(dose_grid, cont_twfe_weights, D=data$dose)
cont_res <- cont_did(data$dy, data$dose)

plot_df <- cbind.data.frame(twfe_weights, dose_grid)

twfe_weights_plot <- ggplot(data=plot_df,
                            mapping=aes(x = dose_grid,
                                        y = twfe_weights)) +
  geom_line(colour = "darkblue", linewidth = 1.2) +
  xlim(c(min(dose_grid),
         max(dose_grid)))+
  ylab("TWFE weights") +
  xlab("Dose") +
  geom_vline(xintercept = mean(data$dose),
             colour="black",
             linewidth = 0.5,
             linetype = "dotted") +
  ylim(c(0,3)) +
  labs(title="TWFE weights")

twfe_weights_plot


cont_res$att.overall

cont_res$acrt.overall
summary(twfe)$coefficients[2,1]




