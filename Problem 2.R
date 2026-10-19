##############################################################################
############## The Delaporte Distribution: A great risk story ################
##############################################################################

#No MCMC approach

set.seed(123)

# 1. True parameters and theoretical moments
r_true <- 2
p_true <- 0.5
lambda_true <- 3
n <- 100

theoretical_mean <- lambda_true +
  r_true * (1 - p_true) / p_true

theoretical_variance <- lambda_true +
  r_true * (1 - p_true) / p_true^2

cat("Theoretical mean:", theoretical_mean, "\n")
cat("Theoretical variance:", theoretical_variance, "\n")

# 2. Simulate Delaporte observations
Y <- rnbinom(n, size = r_true, prob = p_true)
Z <- rpois(n, lambda = lambda_true)
X <- Y + Z

cat("Observed sample mean:", mean(X), "\n")
cat("Observed sample variance:", var(X), "\n")

# 3. Numerically stable log-sum-exp
log_sum_exp <- function(a) {
  m <- max(a)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(a - m)))
}

# 4. Negative log-likelihood
# Optimize unconstrained parameters:
# eta[1] = log(r), eta[2] = logit(p), eta[3] = log(lambda)

negative_loglik <- function(eta, x) {

  r <- exp(eta[1])
  p <- plogis(eta[2])
  lambda <- exp(eta[3])

  loglik_i <- vapply(x, function(xi) {

    y <- 0:xi

    log_terms <-
      dnbinom(y, size = r, prob = p, log = TRUE) +
      dpois(xi - y, lambda = lambda, log = TRUE)

    log_sum_exp(log_terms)

  }, numeric(1))

  if (any(!is.finite(loglik_i))) return(1e100)

  -sum(loglik_i)
}

# 5. Estimate parameters by MLE
initial <- c(
  log(2),
  qlogis(0.5),
  log(3)
)

fit <- optim(
  par = initial,
  fn = negative_loglik,
  x = X,
  method = "BFGS",
  hessian = TRUE,
  control = list(maxit = 5000, reltol = 1e-10)
)

estimates <- c(
  r = exp(fit$par[1]),
  p = plogis(fit$par[2]),
  lambda = exp(fit$par[3])
)

cat("\nOptimization convergence code:", fit$convergence, "\n")
cat("0 indicates that optim reports successful convergence.\n")
cat("\nTrue parameters:\n")
print(c(r = r_true, p = p_true, lambda = lambda_true))

cat("\nMLE estimates:\n")
print(estimates)

cat("\nMaximized log-likelihood:", -fit$value, "\n")

# 6. Compare fitted and theoretical moments
fitted_mean <- estimates["lambda"] +
  estimates["r"] * (1 - estimates["p"]) / estimates["p"]

fitted_variance <- estimates["lambda"] +
  estimates["r"] * (1 - estimates["p"]) / estimates["p"]^2

cat("\nFitted theoretical mean:", fitted_mean, "\n")
cat("Fitted theoretical variance:", fitted_variance, "\n")

# 7. Compare observed and fitted frequencies
max_count <- max(X)
counts <- 0:max_count

observed_freq <- tabulate(X + 1, nbins = max_count + 1)
fitted_prob <- vapply(counts, function(k) {
  log_terms <- dnbinom(
    0:k,
    size = estimates["r"],
    prob = estimates["p"],
    log = TRUE
  ) + dpois(
    k - (0:k),
    lambda = estimates["lambda"],
    log = TRUE
  )
  exp(log_sum_exp(log_terms))
}, numeric(1))

comparison <- data.frame(
  Count = counts,
  Observed = observed_freq,
  Expected = n * fitted_prob
)

print(comparison)

# Visual comparison
plot(
  counts, observed_freq,
  type = "h", lwd = 3,
  xlab = "Simulated count",
  ylab = "Frequency",
  main = "Observed versus fitted Delaporte distribution"
)

points(
  counts, n * fitted_prob,
  pch = 19
)

legend(
  "topright",
  legend = c("Observed", "Fitted expected frequency"),
  lty = c(1, NA),
  pch = c(NA, 19),
  bty = "n"
)


#Gibbs Sampler approach 
library(coda)

set.seed(123)

# 1. Simulate Toy Delaporte Data (Grecontina)
n <- 100
r_true <- 2
p_true <- 0.5
lambda_true <- 3

Y1_true <- rnbinom(n, size = r_true, prob = p_true)
Y2_true <- rpois(n, lambda = lambda_true)
X <- Y1_true + Y2_true

# 2. Log-binomial coefficient helper
log_binom <- function(y, r) {
  lgamma(y + r) - lgamma(r) - lgamma(y + 1)
}

# 3. Griddy Gibbs Sampler Function
run_delaporte_gibbs <- function(X, n_iter = 15000, burn_in = 2000, thin = 1,
                                init_r = 2, init_p = 0.5, init_lambda = 3) {
  n <- length(X)
  r <- init_r
  p <- init_p
  lambda <- init_lambda
  Y1 <- pmin(X, rpois(n, lambda))
  
  n_store <- floor((n_iter - burn_in) / thin)
  samples <- matrix(NA, nrow = n_store, ncol = 3)
  colnames(samples) <- c("r", "p", "lambda")
  
  r_grid <- seq(0.05, 12.0, length.out = 300)
  store_idx <- 1
  
  for (iter in 1:n_iter) {
    # Step A: Sample Latent Y1_i | X_i, r, p, lambda
    for (i in 1:n) {
      y_vals <- 0:X[i]
      log_p1 <- log_binom(y_vals, r) + r * log(p) + y_vals * log(1 - p)
      log_p2 <- (X[i] - y_vals) * log(lambda) - lfactorial(X[i] - y_vals) - lambda
      
      log_probs <- log_p1 + log_p2
      probs <- exp(log_probs - max(log_probs))
      probs <- probs / sum(probs)
      
      Y1[i] <- sample(y_vals, size = 1, prob = probs)
    }
    
    # Step B: Sample p | Y1, r ~ Beta(r*n + 1, sum(Y1) + 1)
    p <- rbeta(1, r * n + 1, sum(Y1) + 1)
    p <- pmin(pmax(p, 0.001), 0.999)
    
    # Step C: Sample lambda | X, Y1 ~ Gamma(2 + sum(X - Y1), 1 + n)
    lambda <- rgamma(1, shape = 2 + sum(X - Y1), rate = 1 + n)
    lambda <- pmax(lambda, 0.01)
    
    # Step D: Griddy Gibbs for r with Gamma(2, 1) Prior
    log_dens_r <- sapply(r_grid, function(rg) {
      log_prior <- (2 - 1) * log(rg) - 1 * rg
      log_like  <- n * rg * log(p) + sum(log_binom(Y1, rg))
      return(log_prior + log_like)
    })
    
    probs_r <- exp(log_dens_r - max(log_dens_r))
    probs_r <- probs_r / sum(probs_r)
    r <- sample(r_grid, size = 1, prob = probs_r)
    
    # Storage
    if (iter > burn_in && (iter - burn_in) %% thin == 0) {
      samples[store_idx, ] <- c(r, p, lambda)
      store_idx <- store_idx + 1
    }
  }
  
  return(mcmc(samples, start = burn_in + 1, end = n_iter, thin = thin))
}

# 4. Run 4 Chains with Diverse Starting Values
cat("Running 4 Gibbs Chains...\n")
chain1 <- run_delaporte_gibbs(X, init_r = 1.0, init_p = 0.2, init_lambda = 1.5)
chain2 <- run_delaporte_gibbs(X, init_r = 4.0, init_p = 0.7, init_lambda = 4.5)
chain3 <- run_delaporte_gibbs(X, init_r = 2.0, init_p = 0.5, init_lambda = 3.0)
chain4 <- run_delaporte_gibbs(X, init_r = 6.0, init_p = 0.3, init_lambda = 2.0)

mcmc_list <- mcmc.list(chain1, chain2, chain3, chain4)

# 5. Convergence Verification Metrics
cat("\n--- Gelman-Rubin Diagnostic (R-hat) ---\n")
print(gelman.diag(mcmc_list)$psrf)

cat("\n--- Effective Sample Size (ESS) ---\n")
print(effectiveSize(mcmc_list))

# 6. Plot Exact Diagnostic Panel (2x3 Layout)
par(mfrow = c(3, 2), mar = c(4.5, 4.5, 2.5, 1.5))

# Colors matching multi-chain trace overlay (Black, Red, Green, Blue)
chain_cols <- c("black", "red", "green3", "deepskyblue")

# Combine all post-burn-in iterations for ACF computation
combined_samples <- as.matrix(mcmc_list)
n_post <- nrow(chain1)

# --- Row 1: r Parameter ---
plot(1:n_post, chain1[, "r"], type = "l", col = chain_cols[1],
     ylim = c(0, 12), xlab = "Post-burn-in iteration", ylab = "r",
     main = "Trace: r")
lines(1:n_post, chain2[, "r"], col = chain_cols[2])
lines(1:n_post, chain3[, "r"], col = chain_cols[3])
lines(1:n_post, chain4[, "r"], col = chain_cols[4])

acf(combined_samples[, "r"], lag.max = 50, main = "ACF: r", ylab = "ACF", lwd = 1.5)

# --- Row 2: p Parameter ---
plot(1:n_post, chain1[, "p"], type = "l", col = chain_cols[1],
     ylim = c(0, 1.0), xlab = "Post-burn-in iteration", ylab = "p",
     main = "Trace: p")
lines(1:n_post, chain2[, "p"], col = chain_cols[2])
lines(1:n_post, chain3[, "p"], col = chain_cols[3])
lines(1:n_post, chain4[, "p"], col = chain_cols[4])

acf(combined_samples[, "p"], lag.max = 50, main = "ACF: p", ylab = "ACF", lwd = 1.5)

# --- Row 3: lambda Parameter ---
plot(1:n_post, chain1[, "lambda"], type = "l", col = chain_cols[1],
     ylim = c(0.5, 5.5), xlab = "Post-burn-in iteration", ylab = "lambda",
     main = "Trace: lambda")
lines(1:n_post, chain2[, "lambda"], col = chain_cols[2])
lines(1:n_post, chain3[, "lambda"], col = chain_cols[3])
lines(1:n_post, chain4[, "lambda"], col = chain_cols[4])

acf(combined_samples[, "lambda"], lag.max = 50, main = "ACF: lambda", ylab = "ACF", lwd = 1.5)

# Reset plotting grid
par(mfrow = c(1, 1))
