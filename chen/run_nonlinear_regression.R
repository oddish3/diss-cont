# Main file to run regression simulations
rm(list = ls())
cat("\014")

library(splines2)
library(tictoc)
library(parallel)
library(data.table)

# Define source files
source_files <- c(
  "~/Documents/uni/master-dissertation/code-cont/chen/jhat.R",
  "~/Documents/uni/master-dissertation/code-cont/chen/jlep.R",
  "~/Documents/uni/master-dissertation/code-cont/chen/npiv.R",
  "~/Documents/uni/master-dissertation/code-cont/chen/npiv_estimate.R",
  "~/Documents/uni/master-dissertation/code-cont/chen/bspline.R",
  "~/Documents/uni/master-dissertation/code-cont/chen/ucb_cv.R",
  "~/Documents/uni/master-dissertation/code-cont/chen/ucb_cvge.R",
  "~/Documents/uni/master-dissertation/code-cont/ucb_cc.R"
)

# Source all required functions
invisible(lapply(source_files, source))

# Setup parameters
array_value <- 1 # as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
trimming <- if (array_value <= 4) TRUE else FALSE
n_index <- if (trimming) array_value else array_value - 4

# Inputs
nn <- c(1250, 2500, 5000, 10000)
nm <- 1000
nb <- 1000
nx <- 1000
nL <- 9
r <- 4
M <- 5
alpha <- c(0.10, 0.05, 0.01)
nj <- 4

# Pre-compute
TJ <- 2^((0:nL) + 0) + r - 1
CJ <- c(0, cumsum(TJ))
Xx <- seq(0, 1, length.out = nx + 1)
Xx_sub <- if (trimming) Xx[Xx > 0.01 & Xx <= 0.99] else Xx
h0 <- sin(15 * pi * Xx) * cos(Xx)

# Precompute Px
Px <- do.call(cbind, lapply(0:nL, function(ll) bspline(Xx_sub, ll, r)))

# Main simulation function
run_simulation <- function(j, nn, n_index, nb, nx, nL, r, M, alpha, nj, TJ, CJ, Xx, Xx_sub, h0, Px) {
  set.seed(1234567 + j)  # Ensure reproducibility in parallel
  
  log_file <- paste0("debug/debug_log_", j, ".txt")
  
  cat(paste("Running simulation for job:", j, "\n"), file = log_file)
  
  n <- nn[n_index]
  
  # Simulate data
  x <- runif(n)
  u <- rnorm(n)
  y <- sin(15 * pi * x) * cos(x) + u
  
  cat("Data simulation complete\n", file = log_file, append = TRUE)
  
  # Pre-compute basis functions
  PP <- tryCatch({
    do.call(cbind, lapply(0:nL, function(ll) bspline(x, ll, r)))
  }, error = function(e) {
    cat(paste("Error in creating PP:", e, "\n"), file = log_file, append = TRUE)
    NULL
  })
  
  if (is.null(PP)) {
    cat("PP object not created\n", file = log_file, append = TRUE)
    stop("PP object not found")
  } else {
    cat(paste("PP object created with dimensions:", dim(PP), "\n"), file = log_file, append = TRUE)
  }
  
  cat("Basis functions (PP) computed for job:", j, "\n", file = log_file, append = TRUE)
  
  # Compute resolution levels and estimates
  result1 <- Jhat(PP, PP, CJ, CJ, TJ, M, n, nL)
  result2 <- Jlep(result1$LL, Px, PP, PP, CJ, CJ, TJ, y, n, nb)
  Ltil <- max(min(result2$LL, result1$LL - 1), 0)
  
  cat("Resolution levels computed for job:", j, "\n", file = log_file, append = TRUE)
  
  npiv_result <- npiv_estimate(Ltil, Px, PP, PP, CJ, CJ, y, n)
  
  cat("NPIV estimates computed for job:", j, "\n", file = log_file, append = TRUE)
  
  # Compute critical values and coverage
  zast_j <- ucb_cv(Ltil, result1$LL, Px, PP, PP, CJ, CJ, y, n, nb, 0, alpha)
  cvge_result <- ucb_cvge(h0[which(Xx %in% Xx_sub)], npiv_result$hhat, npiv_result$sigh, zast_j, result2$theta, log(log(TJ[result2$LL + 1])))
  
  cat("Critical values and coverage computed for job:", j, "\n", file = log_file, append = TRUE)
  
  # Compute undersmoothed estimates
  us_results <- lapply(1:nj, function(k) {
    res <- npiv_estimate(k + 2, Px, PP, PP, CJ, CJ, y, n)
    zdet_k <- ucb_cc(k + 2, Px, PP, PP, CJ, CJ, y, n, nb, 0, alpha)
    cvge_k <- ucb_cvge(h0[Xx %in% Xx_sub], res$hha, res$sig, zdet_k, 0, 0)
    list(hha1 = res$hha, sig1 = res$sig, zdet_k = zdet_k, cvge_k = cvge_k$check)
  })
  
  cat("Undersmoothed estimates computed for job:", j, "\n", file = log_file, append = TRUE)
  
  # Compile results
  result <- list(
    Lhat = result1$LL, Llep = result2$LL, Ltil = Ltil,
    flag = result1$flag, thet = result2$theta,
    zast = zast_j, hhat = npiv_result$hhat, sigh = npiv_result$sigh,
    cvge = c(cvge_result$check, sapply(us_results, function(x) x$cvge_k)),
    loss = c(max(abs(h0[which(Xx %in% Xx_sub)] - npiv_result$hhat)),
             sapply(us_results, function(x) max(abs(h0[Xx %in% Xx_sub] - x$hha1)))),
    rati = sapply(us_results, function(x) max(x$sig1) / max(npiv_result$sigh)),
    zdet = sapply(us_results, function(x) x$zdet_k)
  )
  
  # Save intermediate results for debugging
  save(result, file = paste0("result_debug_", j, ".RData"))
  
  cat("Results saved for job:", j, "\n", file = log_file, append = TRUE)
  
  return(result)
}

# debugonce(run_simulation)
# run_simulation(1, nn, n_index, nb, nx, nL, r, M, alpha, nj, TJ, CJ, Xx, Xx_sub, h0, Px)

# Run simulations in parallel
cl <- makeCluster(detectCores() - 1)

# Export all necessary objects and functions to the cluster
clusterExport(cl, c("nn", "n_index", "nb", "nx", "nL", "r", "M", "alpha", "nj", "TJ", "CJ", "Xx", "Xx_sub", "h0", "Px", 
                    "bspline", "Jhat", "Jlep", "npiv_estimate", "ucb_cv", "ucb_cvge", "ucb_cc"))


# Load necessary libraries and source files on each cluster node
clusterEvalQ(cl, {
  library(splines2)
  # source the necessary files directly
  source("~/Documents/uni/master-dissertation/code-cont/chen/jhat.R")
  source("~/Documents/uni/master-dissertation/code-cont/chen/jlep.R")
  source("~/Documents/uni/master-dissertation/code-cont/chen/npiv.R")
  source("~/Documents/uni/master-dissertation/code-cont/chen/npiv_estimate.R")
  source("~/Documents/uni/master-dissertation/code-cont/chen/bspline.R")
  source("~/Documents/uni/master-dissertation/code-cont/chen/ucb_cv.R")
  source("~/Documents/uni/master-dissertation/code-cont/chen/ucb_cvge.R")
  source("~/Documents/uni/master-dissertation/code-cont/ucb_cc.R")
})


tic()
test_js <- c(1, 3, 4, 5, 7, 8, 9)
sim_results <- parLapply(cl, test_js, run_simulation, 
                         nn, n_index, nb, nx, nL, r, M, alpha, nj, TJ, CJ, Xx, Xx_sub, h0, Px)
toc()

stopCluster(cl)

# Process results
process_results <- function(sim_results) {
  results <- data.table(
    Lhat = sapply(sim_results, function(x) x$Lhat),
    Llep = sapply(sim_results, function(x) x$Llep),
    Ltil = sapply(sim_results, function(x) x$Ltil),
    flag = sapply(sim_results, function(x) x$flag),
    thet = sapply(sim_results, function(x) x$thet)
  )
  zast <- do.call(rbind, lapply(sim_results, function(x) x$zast))
  zdet <- array(unlist(lapply(sim_results, function(x) x$zdet)), dim = c(nm, length(alpha), nj))
  cvge <- do.call(rbind, lapply(sim_results, function(x) x$cvge))
  loss <- do.call(rbind, lapply(sim_results, function(x) x$loss))
  rati <- do.call(rbind, lapply(sim_results, function(x) x$rati))
  
  list(results = results, zast = zast, zdet = zdet, cvge = cvge, loss = loss, rati = rati)
}

final_results <- process_results(sim_results)

# Output results
print(colMeans(final_results$loss[1:12,]))
print(colMeans(final_results$cvge[1:12,]))
print(colMeans(final_results$rati[1:12,]))

# Save results
save_filename <- paste0("./results/regression_", array_value, ".RData")
save(list = c("final_results", "TJ"), file = save_filename)
