library(mvtnorm)
library(coda)
library(MASS)
library(ggplot2)
# Set seed for reproducibility
set.seed(123)

# Parameters
d <- 100
s <- 0.01 * (1:d)
Sigma <- diag(s^2)
Sigma_inv <- diag(1/s^2)
n_iter <- 1000
burn_in <- 100
M <- 100  # Number of independent MCMC runs
n_values <- c(0.1, 1, 10)  # Values of n to test for control variates

# Optimal gamma for GI-MALA
gamma_gi_mala <- 1

# Log-density and gradient functions
log_pi <- function(x) -0.5 * t(x) %*% Sigma_inv %*% x - 0.5 * sum(log(s^2)) - 0.5 * d * log(2 * pi)
grad_log_pi <- function(x) as.vector(-Sigma_inv %*% x)

# Proposal function for GI-MALA
gi_mala_proposal <- function(x, gamma) {
  as.vector(rmvnorm(1, mean = x + gamma * Sigma %*% grad_log_pi(x), sigma = (2 * gamma - gamma^2) * Sigma))
}

# Log-proposal density for GI-MALA
log_proposal_density <- function(y, x, gamma) {
  mu <- x + gamma * Sigma %*% grad_log_pi(x)
  dmvnorm(y, mean = mu, sigma = (2 * gamma - gamma^2) * Sigma, log = TRUE)
}

# Control variate functions
h_x <- function(x, n) {
  x * exp(-x^2 / n)  # h_i(x) = x_i * exp(-x_i^2/n) for F(x) = x
}

E_x2_exp <- function(s_i, n) {
  n * s_i^2 / (2 * s_i^2 + n)  # E[x_i^2 * exp(-x_i^2/n)]
}

h_xxT_diag <- function(x_i, s_i, n) {
  x_i^2 * exp(-x_i^2 / n) - E_x2_exp(s_i, n)  # For xx^T diagonal
}

h_xxT_offdiag <- function(x_i, x_j, n) {
  x_i * x_j * exp(-(x_i^2 + x_j^2) / n)  # For xx^T off-diagonal
}

# Modified MH algorithm to compute estimators
mh_algorithm <- function(gamma, type, n) {
  x <- rep(0, d)
  chain <- matrix(NA, n_iter, d)
  accept <- 0
  
  for (i in 1:n_iter) {
    y <- gi_mala_proposal(x, gamma)
    log_alpha <- log_pi(y) - log_pi(x) + 
                 log_proposal_density(x, y, gamma) - 
                 log_proposal_density(y, x, gamma)
    if (runif(1) < min(1, exp(log_alpha))) {
      x <- y
      accept <- accept + 1
    }
    chain[i, ] <- x
  }
  
  # Use post-burn-in samples
  chain_post_burn <- chain[-(1:burn_in), ]
  N_eff <- nrow(chain_post_burn)
  
  # Standard Monte Carlo estimates
  est_x <- colMeans(chain_post_burn)  # For F(x) = x
  est_xxT_11 <- mean(chain_post_burn[,1]^2)  # xx^T (1,1)
  est_xxT_100100 <- mean(chain_post_burn[,100]^2)  # xx^T (100,100)
  est_xxT_12 <- mean(chain_post_burn[,1] * chain_post_burn[,2])  # xx^T (1,2)
  
  # Control variate for F(x) = x
  h_x_values <- t(apply(chain_post_burn, 1, h_x, n = n))  # N_eff x d
  cov_x_h <- colMeans(chain_post_burn * h_x_values - 
                      colMeans(chain_post_burn) * colMeans(h_x_values))
  var_h <- apply(h_x_values, 2, var) * (N_eff - 1) / N_eff
  c_opt <- cov_x_h / var_h
  est_x_cv <- est_x - c_opt * colMeans(h_x_values)
  
  # Control variate for xx^T (1,1)
  h_xxT_11 <- sapply(chain_post_burn[,1], h_xxT_diag, s_i = s[1], n = n)
  cov_xxT_11 <- mean((chain_post_burn[,1]^2 - mean(chain_post_burn[,1]^2)) * 
                     (h_xxT_11 - mean(h_xxT_11)))
  var_h_xxT_11 <- var(h_xxT_11) * (N_eff - 1) / N_eff
  c_opt_11 <- cov_xxT_11 / var_h_xxT_11
  est_xxT_11_cv <- est_xxT_11 - c_opt_11 * mean(h_xxT_11)
  
  # Control variate for xx^T (100,100)
  h_xxT_100 <- sapply(chain_post_burn[,100], h_xxT_diag, s_i = s[100], n = n)
  cov_xxT_100 <- mean((chain_post_burn[,100]^2 - mean(chain_post_burn[,100]^2)) * 
                      (h_xxT_100 - mean(h_xxT_100)))
  var_h_xxT_100 <- var(h_xxT_100) * (N_eff - 1) / N_eff
  c_opt_100 <- cov_xxT_100 / var_h_xxT_100
  est_xxT_100100_cv <- est_xxT_100100 - c_opt_100 * mean(h_xxT_100)
  
  # Control variate for xx^T (1,2)
  h_xxT_12 <- mapply(h_xxT_offdiag, chain_post_burn[,1], chain_post_burn[,2], 
                     MoreArgs = list(n = n))
  cov_xxT_12 <- mean((chain_post_burn[,1] * chain_post_burn[,2] - 
                      mean(chain_post_burn[,1] * chain_post_burn[,2])) * 
                     (h_xxT_12 - mean(h_xxT_12)))
  var_h_xxT_12 <- var(h_xxT_12) * (N_eff - 1) / N_eff
  c_opt_12 <- cov_xxT_12 / var_h_xxT_12
  est_xxT_12_cv <- est_xxT_12 - c_opt_12 * mean(h_xxT_12)
  
  list(
    ar = accept/n_iter,
    est_x = est_x, est_x_cv = est_x_cv,
    est_xxT_11 = est_xxT_11, est_xxT_11_cv = est_xxT_11_cv,
    est_xxT_100100 = est_xxT_100100, est_xxT_100100_cv = est_xxT_100100_cv,
    est_xxT_12 = est_xxT_12, est_xxT_12_cv = est_xxT_12_cv
  )
}

# Run M simulations for each n and compute variances
results <- lapply(n_values, function(n) {
  run_results <- replicate(M, mh_algorithm(gamma_gi_mala, "gi_mala", n), 
                          simplify = FALSE)
  
  # Extract estimates
  est_x <- t(sapply(run_results, function(r) r$est_x))
  est_x_cv <- t(sapply(run_results, function(r) r$est_x_cv))
  est_xxT_11 <- sapply(run_results, function(r) r$est_xxT_11)
  est_xxT_11_cv <- sapply(run_results, function(r) r$est_xxT_11_cv)
  est_xxT_100100 <- sapply(run_results, function(r) r$est_xxT_100100)
  est_xxT_100100_cv <- sapply(run_results, function(r) r$est_xxT_100100_cv)
  est_xxT_12 <- sapply(run_results, function(r) r$est_xxT_12)
  est_xxT_12_cv <- sapply(run_results, function(r) r$est_xxT_12_cv)
  
  # Compute variances
  var_x <- apply(est_x, 2, var)
  var_x_cv <- apply(est_x_cv, 2, var)
  var_xxT_11 <- var(est_xxT_11)
  var_xxT_11_cv <- var(est_xxT_11_cv)
  var_xxT_100100 <- var(est_xxT_100100)
  var_xxT_100100_cv <- var(est_xxT_100100_cv)
  var_xxT_12 <- var(est_xxT_12)
  var_xxT_12_cv <- var(est_xxT_12_cv)
  
  # Average acceptance rate
  ar <- mean(sapply(run_results, function(r) r$ar))
  
  return(list(
    n = n,
    ar = ar,
    var_x = var_x, var_x_cv = var_x_cv,
    var_xxT_11 = var_xxT_11, var_xxT_11_cv = var_xxT_11_cv,
    var_xxT_100100 = var_xxT_100100, var_xxT_100100_cv = var_xxT_100100_cv,
    var_xxT_12 = var_xxT_12, var_xxT_12_cv = var_xxT_12_cv
  ))
})

# Summarize results
cat("Variance Reduction Results for GI-MALA (n_iter = 1000, M = 100):\n")
for (res in results) {
  cat(sprintf("\nn = %g\n", res$n))
  cat(sprintf("Average Acceptance Rate: %g\n", res$ar))
  cat("F(x) = x (mean variance across dimensions):\n")
  cat(sprintf("Standard MC: %g\n", mean(res$var_x)))
  cat(sprintf("Control Variate: %g\n", mean(res$var_x_cv)))
  cat(sprintf("Reduction: %g%%\n", 100 * (mean(res$var_x) - mean(res$var_x_cv)) / mean(res$var_x)))
  
  cat("F(x) = xx^T (1,1):\n")
  cat(sprintf("Standard MC: %g\n", res$var_xxT_11))
  cat(sprintf("Control Variate: %g\n", res$var_xxT_11_cv))
  cat(sprintf("Reduction: %g%%\n", 100 * (res$var_xxT_11 - res$var_xxT_11_cv) / res$var_xxT_11))
  
  cat("F(x) = xx^T (100,100):\n")
  cat(sprintf("Standard MC: %g\n", res$var_xxT_100100))
  cat(sprintf("Control Variate: %g\n", res$var_xxT_100100_cv))
  cat(sprintf("Reduction: %g%%\n", 100 * (res$var_xxT_100100 - res$var_xxT_100100_cv) / res$var_xxT_100100))
  
  cat("F(x) = xx^T (1,2):\n")
  cat(sprintf("Standard MC: %g\n", res$var_xxT_12))
  cat(sprintf("Control Variate: %g\n", res$var_xxT_12_cv))
  cat(sprintf("Reduction: %g%%\n", 100 * (res$var_xxT_12 - res$var_xxT_12_cv) / res$var_xxT_12))
}


library(ggplot2)

# Data from provided results
results <- list(
  list(
    n = 0.1,
    var_x = 0.000376327, var_x_cv = 0.000342864,
    var_xxT_11 = 1.70798e-11, var_xxT_11_cv = 2.75506e-16,
    var_xxT_100100 = 0.00205851, var_xxT_100100_cv = 0.0031882,
    var_xxT_12 = 5.05443e-11, var_xxT_12_cv = 5.23725e-15
  ),
  list(
    n = 1,
    var_x = 0.000377875, var_x_cv = 0.000148285,
    var_xxT_11 = 2.32793e-11, var_xxT_11_cv = 2.59525e-18,
    var_xxT_100100 = 0.00159223, var_xxT_100100_cv = 0.00102077,
    var_xxT_12 = 4.4299e-11, var_xxT_12_cv = 4.52837e-17
  ),
  list(
    n = 10,
    var_x = 0.000370942, var_x_cv = 7.57708e-06,
    var_xxT_11 = 1.67181e-11, var_xxT_11_cv = 2.85198e-20,
    var_xxT_100100 = 0.00275132, var_xxT_100100_cv = 0.000363311,
    var_xxT_12 = 5.08204e-11, var_xxT_12_cv = 5.31915e-19
  )
)



