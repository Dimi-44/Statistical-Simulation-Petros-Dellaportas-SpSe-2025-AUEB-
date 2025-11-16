##############################################################################
############## The Delaporte Distribution: A great risk story ################
##############################################################################

######Task 1: simple simulation process ######
# Set seed for reproducibility
set.seed(123)

# Parameters for Delaporte distribution
r <- 2
p <- 0.5
lambda <- 3
n <- 100
theoretical_mean <- lambda + r * (1 - p) / p
theoretical_variance <- lambda + r * (1 - p) / (p^2)

sample_mean<-rep(0,10000)
sample_variance<-rep(0,10000)
for (i in 1:10000) {
# Simulate Negative Binomial and Poisson components
Y1 <- rnbinom(n, size = r, prob = p)  # Negative Binomial(r, p)
Y2 <- rpois(n, lambda = lambda)       # Poisson(lambda)
Y1;Y2
# Sum to obtain Delaporte samples
X <- Y1 + Y2
sample_mean[i] <- mean(X)
sample_variance[i] <- var(X)
}
mean(sample_mean);mean(sample_variance)

######Task 2.1: Gibbs with rejection!!!######

# Revised Gibbs Sampler for Delaporte Parameter Estimation
set.seed(123)

# Simulated data (from Task 1)
n <- 100
r_true <- 2
p_true <- 0.5
lambda_true <- 3
Y1 <- rnbinom(n, size = r_true, prob = p_true)
Y2 <- rpois(n, lambda = lambda_true)
X <- Y1 + Y2

# Initialize parameters closer to true values
r <- 1.7
p <- 0.8
lambda <- 2.6
Y1 <- pmin(X, rpois(n, lambda))  # Initial latent variables
n_iter <- 10000
burn_in <- 3000
thin <- 40

# Storage for samples (adjusted for thinning)
n_store <- floor((n_iter - burn_in) / thin)
r_samples <- numeric(n_store)
p_samples <- numeric(n_store)
lambda_samples <- numeric(n_store)
store_idx <- 1

# Function to compute log-binomial coefficient safely
log_binom <- function(y, r) {
  if (y < 0 || r <= 0) return(-Inf)
  lgamma(y + r) - lgamma(r) - lgamma(y + 1)
}

# Rejection sampling for r with Gamma(2, 1) proposal
sample_r <- function(p, Y1, n) {
  target <- function(r) {
    if (r <= 0.01 || r >= 8) return(-Inf)
    log_p_term <- n * r * log(p)
    log_prior <- (0.01 - 1) * log(r) - 0.01 * r
    log_binom_sum <- sum(sapply(Y1, function(y) log_binom(y, r)))
    log_p_term + log_prior + log_binom_sum
  }
  
  # Proposal: Gamma(2, 1)
  proposal <- function() rgamma(1, shape = 2, rate = 1)
  log_proposal <- function(r) dgamma(r, shape = 2, rate = 1, log = TRUE)
  
  # Numerically find M
  log_f <- function(r) target(r) - log_proposal(r)
  opt <- optimize(log_f, interval = c(0.01, 8), maximum = TRUE)
  log_M <- opt$objective
  
  # Rejection sampling loop
  while (TRUE) {
    r_cand <- proposal()
    if (r_cand < 0.01 || r_cand > 8) next
    log_ratio <- log_f(r_cand) - log_M
    if (is.na(log_ratio) || is.nan(log_ratio)) next
    if (log(runif(1)) <= log_ratio) return(r_cand)
  }
}

# Gibbs Sampler
for (iter in 1:n_iter) {
  # Sample Y_{i1}
  for (i in 1:n) {
    max_y <- X[i]
    y_vals <- 0:max_y
    log_probs <- sapply(y_vals, function(y) {
      if (y > X[i]) return(-Inf)
      log_binom(y, r) + y * log(1 - p) + (X[i] - y) * log(lambda) - lfactorial(X[i] - y)
    })
    log_probs[!is.finite(log_probs)] <- -Inf
    if (all(log_probs == -Inf)) log_probs[1] <- 0  # Avoid all -Inf
    probs <- exp(log_probs - max(log_probs))
    probs <- pmax(probs, 1e-10)  # Prevent zero probabilities
    probs <- probs / sum(probs)
    Y1[i] <- sample(y_vals, 1, prob = probs)
  }
  
  # Sample p ~ Beta(r * n + 1, sum(Y1) + 1)
  p <- rbeta(1, r * n + 1, sum(Y1) + 1)
  
  # Sample lambda ~ Gamma(0.01 + sum(X - Y1), 0.01 + n)
  lambda <- rgamma(1, shape = 0.01 + sum(X - Y1), rate = 0.01 + n)
  lambda <- pmin(pmax(lambda, 0.01), 10)  # Bound lambda
  
  # Sample r using rejection sampling
  r <- sample_r(p, Y1, n)
  
  # Store samples with thinning
  if (iter > burn_in && (iter - burn_in) %% thin == 0) {
    r_samples[store_idx] <- r
    p_samples[store_idx] <- p
    lambda_samples[store_idx] <- lambda
    store_idx <- store_idx + 1
  }

 if (iter %% 100 == 0) {cat(paste0('  iteration: ', iter, '    '), '\r')}
}

# Summary statistics
cat("Mean r:", mean(r_samples), "(True:", r_true, ")\n")
cat("Mean p:", mean(p_samples), "(True:", p_true, ")\n")
cat("Mean lambda:", mean(lambda_samples), "(True:", lambda_true, ")\n")


# Posterior diagnostics and visualization
par(mar = c(4, 5, 2, 1))
layout(matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE), widths = c(1.5, 1, 1), heights = c(1, 1, 1))

# r diagnostics
plot(r_samples, col = "blue", type = 'l', lwd = 1.5, 
     cex.lab = 1.5, cex.axis = 1.5, xlab = 'Iteration', ylab = 'r')
acf(r_samples, lag.max = 100, main = 'r ACF', lwd = 1.5, 
    cex.lab = 1.5, cex.axis = 1.5)
hist(r_samples, cex.lab = 1.5, cex.axis = 1.5, xlab = 'r', freq = FALSE, 
     ylab = 'Posterior Density', col = 'lightgray', main = '')

# p diagnostics
plot(p_samples, col = "darkgreen", type = 'l', lwd = 1.5, 
     cex.lab = 1.5, cex.axis = 1.5, xlab = 'Iteration', ylab = 'p')
acf(p_samples, lag.max = 100, main = 'p ACF', lwd = 1.5, 
    cex.lab = 1.5, cex.axis = 1.5)
hist(p_samples, cex.lab = 1.5, cex.axis = 1.5, xlab = 'p', freq = FALSE, 
     ylab = 'Posterior Density', col = 'lightgray', main = '')

# lambda diagnostics
plot(lambda_samples, col = "purple", type = 'l', lwd = 1.5, 
     cex.lab = 1.5, cex.axis = 1.5, xlab = 'Iteration', ylab = 'lambda')
acf(lambda_samples, lag.max = 100, main = 'lambda ACF', lwd = 1.5, 
    cex.lab = 1.5, cex.axis = 1.5)
hist(lambda_samples, cex.lab = 1.5, cex.axis = 1.5, xlab = 'lambda', freq = FALSE, 
     ylab = 'Posterior Density', col = 'lightgray', main = '')

###### "Importance Sampling" ######
set.seed(123)

# Simulated data
n <- 100
r_true <- 2
p_true <- 0.5
lambda_true <- 3
Y1 <- rnbinom(n, size = r_true, prob = p_true)
Y2 <- rpois(n, lambda = lambda_true)
X <- Y1 + Y2

# Initialize parameters
r <- 2
p <- 0.8
lambda <- 2.2
Y1 <- pmin(X, rpois(n, lambda))  # Initial latent variables
n_iter <- 10000
burn_in <- 3000
thin <- 10

# Storage for samples
n_store <- floor((n_iter - burn_in) / thin)
r_samples <- numeric(n_store)
p_samples <- numeric(n_store)
lambda_samples <- numeric(n_store)
store_idx <- 1

# Function to compute log-binomial coefficient vectorized
log_binom <- function(y, r) {
  ifelse(y < 0 | r <= 0, -Inf, lgamma(y + r) - lgamma(r) - lgamma(y + 1))
}

# Importance sampling for r with optimized Gamma proposal
sample_r <- function(p, Y1, n, r_prev, max_r = 10) {
  K <- 1000  # Reduced from 10,000
  shape <- 2  # Tuned proposal
  rate <- 1
  
  # Draw candidates
  r_cands <- rgamma(K, shape = shape, rate = rate)
  r_cands <- r_cands[r_cands > 0 & r_cands <= max_r]
  if (length(r_cands) < K) r_cands <- c(r_cands, rep(r_prev, K - length(r_cands)))
  r_cands <- r_cands[1:K]
  
  # Compute log unnormalized target density (vectorized)
  log_target <- function(r) {
    if (r <= 0 || r > max_r) return(-Inf)
    n * r * log(p) + (0.01 - 1) * log(r) - 0.01 * r + sum(log_binom(Y1, r))
  }
  log_targets <- vapply(r_cands, log_target, numeric(1))
  
  # Compute log proposal density
  log_proposals <- dgamma(r_cands, shape = shape, rate = rate, log = TRUE)
  
  # Compute importance weights
  log_weights <- log_targets - log_proposals
  max_log_weight <- max(log_weights, na.rm = TRUE)
  if (is.infinite(max_log_weight) || all(is.na(log_weights))) return(r_prev)
  weights <- exp(log_weights - max_log_weight)
  norm_weights <- weights / sum(weights)
  
  # Calculate Effective Sample Size (ESS)
  ess <- 1 / sum(norm_weights^2)
  if (ess < 10 || any(is.na(norm_weights))) return(r_prev)
  
  # Sample new r
  r_new <- sample(r_cands, 1, prob = norm_weights)
  return(r_new)
}

# Gibbs Sampler
for (iter in 1:n_iter) {
  # Sample Y_{i1}
  for (i in 1:n) {
    max_y <- X[i]
    y_vals <- 0:max_y  # Use full range for simplicity
    log_probs <- vapply(y_vals, function(y) {
      if (y > X[i]) return(-Inf)
      log_binom(y, r) +
        y * log(1 - p) + (X[i] - y) * log(lambda) - lfactorial(X[i] - y)
    }, numeric(1))
    log_probs[!is.finite(log_probs)] <- -Inf
    if (all(log_probs == -Inf)) {
      log_probs <- rep(-Inf, length(y_vals))  # Ensure same length
      log_probs[1] <- 0  # Default to y = 0
    }
    probs <- exp(log_probs - max(log_probs))
    probs <- pmax(probs, 1e-10)
    probs <- probs / sum(probs)
    # Ensure lengths match
    if (length(y_vals) != length(probs)) {
      y_vals <- 0:max_y
      log_probs <- vapply(y_vals, function(y) {
        if (y > X[i]) return(-Inf)
        log_binom(y, r) +
          y * log(1 - p) + (X[i] - y) * log(lambda) - lfactorial(X[i] - y)
      }, numeric(1))
      log_probs[!is.finite(log_probs)] <- -Inf
      if (all(log_probs == -Inf)) {
        log_probs <- rep(-Inf, length(y_vals))
        log_probs[1] <- 0
      }
      probs <- exp(log_probs - max(log_probs))
      probs <- pmax(probs, 1e-10)
      probs <- probs / sum(probs)
    }
    Y1[i] <- sample(y_vals, 1, prob = probs)
  }
  
  # Sample p ~ Beta(r * n + 1, sum(Y1) + 1)
  p <- rbeta(1, r * n + 1, sum(Y1) + 1)
  
  # Sample lambda ~ Gamma(0.01 + sum(X - Y1), 0.01 + n)
  lambda <- rgamma(1, shape = 0.01 + sum(X - Y1), rate = 0.01 + n)
  lambda <- pmin(pmax(lambda, 0.01), 10)
  
  # Sample r using importance sampling
  r <- sample_r(p, Y1, n, r)
  
  # Store samples with thinning
  if (iter > burn_in && (iter - burn_in) %% thin == 0) {
    r_samples[store_idx] <- r
    p_samples[store_idx] <- p
    lambda_samples[store_idx] <- lambda
    store_idx <- store_idx + 1
  }
if (iter %% 100 == 0) {cat(paste0('  iteration: ', iter, '    '), '\r')}
}

# Summary statistics
cat("Mean r:", mean(r_samples), "(True:", r_true, ")\n")
cat("Mean p:", mean(p_samples), "(True:", p_true, ")\n")
cat("Mean lambda:", mean(lambda_samples), "(True:", lambda_true, ")\n")

# Posterior diagnostics and visualization
par(mar = c(4, 5, 2, 1))
layout(matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE), widths = c(1.5, 1, 1), heights = c(1, 1, 1))

# r diagnostics
plot(r_samples, col = "blue", type = 'l', lwd = 1.5, 
     cex.lab = 1.5, cex.axis = 1.5, xlab = 'Iteration', ylab = 'r')
acf(r_samples, lag.max = 100, main = 'r ACF', lwd = 1.5, 
    cex.lab = 1.5, cex.axis = 1.5)
hist(r_samples, cex.lab = 1.5, cex.axis = 1.5, xlab = 'r', freq = FALSE, 
     ylab = 'Posterior Density', col = 'lightgray', main = '')

# p diagnostics
plot(p_samples, col = "darkgreen", type = 'l', lwd = 1.5, 
     cex.lab = 1.5, cex.axis = 1.5, xlab = 'Iteration', ylab = 'p')
acf(p_samples, lag.max = 100, main = 'p ACF', lwd = 1.5, 
    cex.lab = 1.5, cex.axis = 1.5)
hist(p_samples, cex.lab = 1.5, cex.axis = 1.5, xlab = 'p', freq = FALSE, 
     ylab = 'Posterior Density', col = 'lightgray', main = '')

# lambda diagnostics
plot(lambda_samples, col = "purple", type = 'l', lwd = 1.5, 
     cex.lab = 1.5, cex.axis = 1.5, xlab = 'Iteration', ylab = 'lambda')
acf(lambda_samples, lag.max = 100, main = 'lambda ACF', lwd = 1.5, 
    cex.lab = 1.5, cex.axis = 1.5)
hist(lambda_samples, cex.lab = 1.5, cex.axis = 1.5, xlab = 'lambda', freq = FALSE, 
     ylab = 'Posterior Density', col = 'lightgray', main = '')


