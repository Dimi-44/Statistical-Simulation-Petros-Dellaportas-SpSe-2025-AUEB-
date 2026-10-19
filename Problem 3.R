library(mvtnorm)
library(coda)
library(ggplot2)

set.seed(123)

# =============================================================================
# GLOBAL PARAMETERS & TARGET DISTRIBUTION (d = 100 Anisotropic Gaussian)
# =============================================================================
d <- 100
s <- 0.01 * (1:d)          # Standard deviations: range from 0.01 to 1.0
Sigma <- diag(s^2)
Sigma_inv <- diag(1 / s^2)

# Target log-density and gradient
log_pi <- function(x) {
  -0.5 * sum((x / s)^2) - 0.5 * sum(log(s^2)) - 0.5 * d * log(2 * pi)
}

grad_log_pi <- function(x) {
  -x / (s^2)
}

# =============================================================================
# PART (i): METROPOLIS-HASTINGS ALGORITHM SAMPLERS
# =============================================================================

# (A) Random Walk Metropolis (RWM)
run_rwm <- function(n_iter, sigma_step) {
  x <- rep(0, d)
  chain <- matrix(NA, n_iter, d)
  accept <- 0
  
  for (i in 1:n_iter) {
    y <- x + rnorm(d, mean = 0, sd = sigma_step)
    log_alpha <- log_pi(y) - log_pi(x)
    
    if (log(runif(1)) < log_alpha) {
      x <- y
      accept <- accept + 1
    }
    chain[i, ] <- x
  }
  list(chain = chain, ar = accept / n_iter)
}

# (B) Standard Metropolis Adjusted Langevin Algorithm (MALA)
run_mala <- function(n_iter, dt) {
  x <- rep(0, d)
  chain <- matrix(NA, n_iter, d)
  accept <- 0
  
  for (i in 1:n_iter) {
    grad_x <- grad_log_pi(x)
    mu_y <- x + 0.5 * (dt^2) * grad_x
    y <- mu_y + rnorm(d, mean = 0, sd = dt)
    
    grad_y <- grad_log_pi(y)
    mu_x <- y + 0.5 * (dt^2) * grad_y
    
    log_q_yx <- -0.5 * sum(((y - mu_y) / dt)^2)
    log_q_xy <- -0.5 * sum(((x - mu_x) / dt)^2)
    
    log_alpha <- log_pi(y) - log_pi(x) + log_q_xy - log_q_yx
    
    if (log(runif(1)) < log_alpha) {
      x <- y
      accept <- accept + 1
    }
    chain[i, ] <- x
  }
  list(chain = chain, ar = accept / n_iter)
}

# (C) Gaussian Invariant MALA (GI-MALA)
run_gi_mala <- function(n_iter, gamma) {
  x <- rep(0, d)
  chain <- matrix(NA, n_iter, d)
  accept <- 0
  
  prop_sd <- sqrt(2 * gamma - gamma^2) * s
  
  for (i in 1:n_iter) {
    grad_x <- grad_log_pi(x)
    mu_y <- x + gamma * (s^2) * grad_x
    y <- mu_y + rnorm(d, mean = 0, sd = prop_sd)
    
    grad_y <- grad_log_pi(y)
    mu_x <- y + gamma * (s^2) * grad_y
    
    log_q_yx <- -0.5 * sum(((y - mu_y) / prop_sd)^2)
    log_q_xy <- -0.5 * sum(((x - mu_x) / prop_sd)^2)
    
    log_alpha <- log_pi(y) - log_pi(x) + log_q_xy - log_q_yx
    
    if (log(runif(1)) < log_alpha) {
      x <- y
      accept <- accept + 1
    }
    chain[i, ] <- x
  }
  list(chain = chain, ar = accept / n_iter)
}

# =============================================================================
# PART (i) BENCHMARK EXECUTION & COMPARISON
# =============================================================================
cat("\n=================================================================\n")
cat("  PART (i): BENCHMARKING THREE METROPOLIS-HASTINGS ALGORITHMS    \n")
cat("=================================================================\n")

n_iter_bench <- 5000
burn_in_bench <- 1000

# Step sizes tuned for stability/performance
sigma_rwm <- 0.002 
dt_mala   <- 0.003 
gamma_gi  <- 1.0   

res_rwm     <- run_rwm(n_iter_bench, sigma_rwm)
res_mala    <- run_mala(n_iter_bench, dt_mala)
res_gi_mala <- run_gi_mala(n_iter_bench, gamma_gi)

# Extract post-burn-in chains
chain_rwm     <- mcmc(res_rwm$chain[-(1:burn_in_bench), ])
chain_mala    <- mcmc(res_mala$chain[-(1:burn_in_bench), ])
chain_gi_mala <- mcmc(res_gi_mala$chain[-(1:burn_in_bench), ])

# Effective Sample Size (ESS)
ess_rwm     <- effectiveSize(chain_rwm)
ess_mala    <- effectiveSize(chain_mala)
ess_gi_mala <- effectiveSize(chain_gi_mala)

summary_df <- data.frame(
  Algorithm       = c("RWM", "Standard MALA", "GI-MALA"),
  Acceptance_Rate = c(res_rwm$ar, res_mala$ar, res_gi_mala$ar),
  Mean_ESS        = c(mean(ess_rwm), mean(ess_mala), mean(ess_gi_mala)),
  Min_ESS_Dim1    = c(ess_rwm[1], ess_mala[1], ess_gi_mala[1]),
  Max_ESS_Dim100  = c(ess_rwm[100], ess_mala[100], ess_gi_mala[100])
)

print(summary_df, row.names = FALSE)


# =============================================================================
# PART (ii): CONTROL VARIATE VARIANCE REDUCTION IN GI-MALA
# =============================================================================
cat("\n=================================================================\n")
cat("  PART (ii): CONTROL VARIATE VARIANCE REDUCTION IN GI-MALA       \n")
cat("=================================================================\n")

# Control Variate Auxiliary Functions
h_x <- function(x, n) {
  x * exp(-x^2 / n)
}

E_x2_exp <- function(s_i, n) {
  n * s_i^2 / (2 * s_i^2 + n)
}

h_xxT_diag <- function(x_i, s_i, n) {
  x_i^2 * exp(-x_i^2 / n) - E_x2_exp(s_i, n)
}

h_xxT_offdiag <- function(x_i, x_j, n) {
  x_i * x_j * exp(-(x_i^2 + x_j^2) / n)
}

# Monte Carlo estimator with Control Variate calculation
run_gi_mala_cv <- function(n_iter, burn_in, gamma, n) {
  res <- run_gi_mala(n_iter, gamma)
  chain_post <- res$chain[-(1:burn_in), ]
  N_eff <- nrow(chain_post)
  
  # 1. Standard Monte Carlo Estimates
  est_x          <- colMeans(chain_post)
  est_xxT_11     <- mean(chain_post[, 1]^2)
  est_xxT_100100 <- mean(chain_post[, 100]^2)
  est_xxT_12     <- mean(chain_post[, 1] * chain_post[, 2])
  
  # 2. Control Variate Estimates for F(x) = x
  h_x_vals <- t(apply(chain_post, 1, h_x, n = n))
  cov_x_h  <- colMeans(chain_post * h_x_vals - colMeans(chain_post) * colMeans(h_x_vals))
  var_h    <- apply(h_x_vals, 2, var) * (N_eff - 1) / N_eff
  c_opt    <- ifelse(var_h == 0, 0, cov_x_h / var_h)
  est_x_cv <- est_x - c_opt * colMeans(h_x_vals)
  
  # 3. Control Variate for xx^T (1,1)
  h_11     <- sapply(chain_post[, 1], h_xxT_diag, s_i = s[1], n = n)
  cov_11   <- mean((chain_post[, 1]^2 - est_xxT_11) * (h_11 - mean(h_11)))
  var_h_11 <- var(h_11) * (N_eff - 1) / N_eff
  c_opt_11 <- ifelse(var_h_11 == 0, 0, cov_11 / var_h_11)
  est_xxT_11_cv <- est_xxT_11 - c_opt_11 * mean(h_11)
  
  # 4. Control Variate for xx^T (100,100)
  h_100     <- sapply(chain_post[, 100], h_xxT_diag, s_i = s[100], n = n)
  cov_100   <- mean((chain_post[, 100]^2 - est_xxT_100100) * (h_100 - mean(h_100)))
  var_h_100 <- var(h_100) * (N_eff - 1) / N_eff
  c_opt_100 <- ifelse(var_h_100 == 0, 0, cov_100 / var_h_100)
  est_xxT_100100_cv <- est_xxT_100100 - c_opt_100 * mean(h_100)
  
  # 5. Control Variate for xx^T (1,2)
  h_12     <- mapply(h_xxT_offdiag, chain_post[, 1], chain_post[, 2], MoreArgs = list(n = n))
  cov_12   <- mean((chain_post[, 1] * chain_post[, 2] - est_xxT_12) * (h_12 - mean(h_12)))
  var_h_12 <- var(h_12) * (N_eff - 1) / N_eff
  c_opt_12 <- ifelse(var_h_12 == 0, 0, cov_12 / var_h_12)
  est_xxT_12_cv <- est_xxT_12 - c_opt_12 * mean(h_12)
  
  list(
    ar = res$ar,
    est_x = est_x, est_x_cv = est_x_cv,
    est_xxT_11 = est_xxT_11, est_xxT_11_cv = est_xxT_11_cv,
    est_xxT_100100 = est_xxT_100100, est_xxT_100100_cv = est_xxT_100100_cv,
    est_xxT_12 = est_xxT_12, est_xxT_12_cv = est_xxT_12_cv
  )
}

# Run M independent runs across different n values
M <- 100
n_iter_cv <- 1000
burn_in_cv <- 100
n_values <- c(0.1, 1, 10)

cv_results <- lapply(n_values, function(n_val) {
  runs <- replicate(M, run_gi_mala_cv(n_iter_cv, burn_in_cv, gamma_gi, n_val), simplify = FALSE)
  
  # Extract estimates across M runs
  est_x          <- t(sapply(runs, function(r) r$est_x))
  est_x_cv       <- t(sapply(runs, function(r) r$est_x_cv))
  est_xxT_11     <- sapply(runs, function(r) r$est_xxT_11)
  est_xxT_11_cv  <- sapply(runs, function(r) r$est_xxT_11_cv)
  est_xxT_100100 <- sapply(runs, function(r) r$est_xxT_100100)
  est_xxT_100100_cv <- sapply(runs, function(r) r$est_xxT_100100_cv)
  est_xxT_12     <- sapply(runs, function(r) r$est_xxT_12)
  est_xxT_12_cv  <- sapply(runs, function(r) r$est_xxT_12_cv)
  
  # Compute variances
  var_x          <- apply(est_x, 2, var)
  var_x_cv       <- apply(est_x_cv, 2, var)
  var_xxT_11     <- var(est_xxT_11)
  var_xxT_11_cv  <- var(est_xxT_11_cv)
  var_xxT_100100 <- var(est_xxT_100100)
  var_xxT_100100_cv <- var(est_xxT_100100_cv)
  var_xxT_12     <- var(est_xxT_12)
  var_xxT_12_cv  <- var(est_xxT_12_cv)
  
  list(
    n = n_val,
    mean_var_x = mean(var_x), mean_var_x_cv = mean(var_x_cv),
    var_11 = var_xxT_11, var_11_cv = var_xxT_11_cv,
    var_100 = var_xxT_100100, var_100_cv = var_xxT_100100_cv,
    var_12 = var_xxT_12, var_12_cv = var_xxT_12_cv
  )
})

# Print Results
for (res in cv_results) {
  cat(sprintf("\n--- Results for Control Variate Tuning Scale n = %g ---\n", res$n))
  
  pct_x <- 100 * (res$mean_var_x - res$mean_var_x_cv) / res$mean_var_x
  cat(sprintf("F(x) = x (Mean Var Across Dimensions):\n  Standard MC: %.6e | CV: %.6e | Reduction: %.2f%%\n",
              res$mean_var_x, res$mean_var_x_cv, pct_x))
  
  pct_11 <- 100 * (res$var_11 - res$var_11_cv) / res$var_11
  cat(sprintf("F(x) = xx^T (1,1):\n  Standard MC: %.6e | CV: %.6e | Reduction: %.2f%%\n",
              res$var_11, res$var_11_cv, pct_11))
  
  pct_100 <- 100 * (res$var_100 - res$var_100_cv) / res$var_100
  cat(sprintf("F(x) = xx^T (100,100):\n  Standard MC: %.6e | CV: %.6e | Reduction: %.2f%%\n",
              res$var_100, res$var_100_cv, pct_100))
  
  pct_12 <- 100 * (res$var_12 - res$var_12_cv) / res$var_12
  cat(sprintf("F(x) = xx^T (1,2):\n  Standard MC: %.6e | CV: %.6e | Reduction: %.2f%%\n",
              res$var_12, res$var_12_cv, pct_12))
}
