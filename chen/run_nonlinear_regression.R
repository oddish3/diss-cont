# Main file to run regression simulations
rm(list = ls())
cat("\014")
library(splines2)
library(tictoc)
library(parallel)
library(future.apply)
library(data.table)

# Source functions (consider combining these into a single file)
source_files <- c(
  "~/Documents/uni/master-dissertation/code-cont/chen/jhat.R",
  "~/Documents/uni/master-dissertation/code-cont/chen/jlep.R",
  "~/Documents/uni/master-dissertation/code-cont/chen/npiv.R",
  "~/Documents/uni/master-dissertation/code-cont/chen/npiv_estimate.R",
  "~/Documents/uni/master-dissertation/code-cont/chen/bspline.R",
  "~/Documents/uni/master-dissertation/code-cont/chen/ucb_cv.R",
  "~/Documents/uni/master-dissertation/code-cont/chen/ucb_cvge.R"
)
invisible(lapply(source_files, source))

array_value <- 1 # as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
trimming <- array_value <= 4
n_index <- if (trimming) array_value else array_value - 4

# Inputs
nn <- c(1250, 2500, 5000, 10000) # sample sizes
nm <- 1                         # number of replications
nb <- 1000                       # number of bootstrap draws per replication
nx <- 1000                       # number of points for grid of x values
nL <- 9                          # maximum resolution level for J
r <- 4                           # B-spline order
M <- 5                           # ex ante upper bound on \sup_x h_0(x)
alpha <- c(0.10, 0.05, 0.01)     # level of significance
nj <- 4

# Pre-compute
TJ <- 2^((0:nL) + 0) + r - 1
CJ <- c(0, cumsum(TJ))
Xx <- seq(0, 1, length.out = nx)
Xx_sub <- if (trimming) Xx[Xx > 0.01 & Xx <= 0.99] else Xx
Px <- matrix(0, nrow = length(Xx_sub), ncol = CJ[length(CJ)])
for (ll in 0:nL) {
  Px[, (CJ[ll + 1] + 1):CJ[ll + 2]] <- bspline(Xx_sub, ll, r)
}
h0 <- sin(15 * pi * Xx) * cos(Xx)

# Simulations
set.seed(1234567)
n <- nn[n_index]

# Function to run a single simulation
run_simulation <- function(j, n, nL, r, CJ, TJ, M, Px, nb, alpha, Xx, Xx_sub, h0) {
  x <- runif(n)
  u <- rnorm(n)
  y <- sin(15 * pi * x) * cos(x) + u
  
  PP <- matrix(0, nrow = n, ncol = CJ[length(CJ)])
  for (ll in 0:nL) {
    PP[, (CJ[ll + 1] + 1):CJ[ll + 2]] <- bspline(x, ll, r)
  }
  
  result1 <- Jhat(PP, PP, CJ, CJ, TJ, M, n, nL)
  Lhat <- result1$LL
  flag <- result1$flag
  
  result2 <- Jlep(Lhat, Px, PP, PP, CJ, CJ, TJ, y, n, nb)
  Llep <- result2$LL
  thet <- result2$theta
  
  Ltil <- max(min(Llep, Lhat - 1), 0)
  
  npiv_result <- npiv_estimate(Ltil, Px, PP, PP, CJ, CJ, y, n)
  hhat <- npiv_result$hhat
  sigh <- npiv_result$sigh
  
  zast <- ucb_cv(Ltil, Lhat, Px, PP, PP, CJ, CJ, y, n, nb, 0, alpha)
  
  loss <- max(abs(h0[which(Xx %in% Xx_sub)] - hhat))
  
  results3 <- ucb_cvge(h0[which(Xx %in% Xx_sub)], hhat, sigh, zast, thet, log(log(TJ[Llep + 1])))
  cvge <- results3$check
  
  list(Lhat = Lhat, flag = flag, Llep = Llep, thet = thet, Ltil = Ltil, 
       zast = zast, loss = loss, cvge = cvge)
}

# Set up parallel backend
num_cores <- detectCores() - 1
plan(multisession, workers = num_cores)

# Run simulations in parallel
tic()
results <- future_lapply(1:nm, function(j) {
  run_simulation(j, n, nL, r, CJ, TJ, M, Px, nb, alpha, Xx, Xx_sub, h0)
}, future.seed = TRUE)
toc()

# Process results using data.table for faster operations
results_dt <- rbindlist(lapply(results, as.data.table))

# Further analysis or saving results can be done here
# For example:
# fwrite(results_dt, "simulation_results.csv")